/*
 * 文件名: modelsolver07.cpp
 * 功能描述:
 * 1. 夹层型非均质裂缝（多表皮独立分配）计算核心。
 * 2. 精准映射: 距离水平井跟端的绝对距离 xw，经无因次化和偏移成为积分中心的 xwD。
 * 3. 严格使用 16.118095 进行有因次压力转换，使用合并项 phi_ct 防止微小标度除法截断误差。
 * 4. PWD_composite 全面使用各段 S_vec 进行 2nf+1 维系统的等效流场重构。
 */

#include "modelsolver07.h"
#include "pressurederivativecalculator.h"
#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <QtConcurrent>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ================= 底层安全贝塞尔函数 =================
double ModelSolver07::safe_bessel_k(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    try { return boost::math::cyl_bessel_k(v, x); } catch (...) { return 0.0; }
}

double ModelSolver07::safe_bessel_k_sc(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    if (x > 600.0) return std::sqrt(M_PI / (2.0 * x));
    try { return boost::math::cyl_bessel_k(v, x) * std::exp(x); } catch (...) { return 0.0; }
}

double ModelSolver07::safe_bessel_i_sc(int v, double x) {
    if (x < 0) x = -x;
    if (x > 600.0) return 1.0 / std::sqrt(2.0 * M_PI * x);
    try { return boost::math::cyl_bessel_i(v, x) * std::exp(-x); } catch (...) { return 0.0; }
}

ModelSolver07::ModelSolver07(ModelType type) :
    m_type(type), m_highPrecision(true), m_currentN(0) {
    precomputeStehfestCoeffs(8);
}

ModelSolver07::~ModelSolver07() {}

void ModelSolver07::setHighPrecision(bool h) { m_highPrecision = h; }

QString ModelSolver07::getModelName(ModelType t, bool v) {
    int id = (int)t + 1;
    return QString("独立表皮非均质夹层型模型%1").arg(id);
}

QVector<double> ModelSolver07::generateLogTimeSteps(int c, double s, double e) {
    QVector<double> t;
    if (c <= 0) return t;
    t.reserve(c);
    for (int i = 0; i < c; ++i) {
        t.append(std::pow(10.0, s + (e - s) * i / (c - 1)));
    }
    return t;
}

ModelCurveData ModelSolver07::calculateTheoreticalCurve(const QMap<QString, double>& p, const QVector<double>& pt) {
    QVector<double> tP = pt.isEmpty() ? generateLogTimeSteps(100, -3.0, 4.0) : pt;

    double phi = p.value("phi", 0.05);
    double mu = p.value("mu", 0.5);
    double B = p.value("B", 1.2);
    double Ct = p.value("Ct", 5e-4);
    double q = p.value("q", 50.0);
    double h = p.value("h", 20.0);
    double kf = p.value("kf", 50.0);
    double L = p.value("L", 1000.0) / 2.0;

    if (L < 1e-9) L = 500.0;
    if (phi < 1e-12 || mu < 1e-12 || Ct < 1e-12 || kf < 1e-12) {
        return std::make_tuple(tP, QVector<double>(tP.size(), 0.0), QVector<double>(tP.size(), 0.0));
    }

    // 采用合并因子，防除法导致截断漂移
    double phi_ct = phi * Ct;
    double td_c = 3.6 * kf / (phi_ct * mu * std::pow(L, 2.0));
    // 应用校准的高精度压力换算乘子
    double p_c = 16.118095 * q * mu * B / (kf * h);

    QVector<double> tD_v;
    tD_v.reserve(tP.size());
    for(double t : tP) tD_v.append(td_c * t);

    QMap<QString, double> cP = p;
    int N = (int)cP.value("N", 8);
    cP["N"] = N;
    precomputeStehfestCoeffs(N);

    if (!cP.contains("nf") || cP["nf"] < 1) cP["nf"] = 9;
    cP["L_half"] = L;

    QVector<double> PD_v, Der_v;
    auto func = std::bind(&ModelSolver07::flaplace_composite, this, std::placeholders::_1, std::placeholders::_2);
    calculatePDandDeriv(tD_v, cP, func, PD_v, Der_v);

    QVector<double> fP(tP.size()), fDP(tP.size());
    for(int i = 0; i < tP.size(); ++i) fP[i] = p_c * PD_v[i];

    if (tP.size() > 2) fDP = PressureDerivativeCalculator::calculateBourdetDerivative(tP, fP, 0.2);
    else fDP.fill(0.0);

    return std::make_tuple(tP, fP, fDP);
}

void ModelSolver07::calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& p,
                                        std::function<double(double, const QMap<QString, double>&)> f,
                                        QVector<double>& oPD, QVector<double>& oDer)
{
    int nP = tD.size();
    oPD.resize(nP); oDer.resize(nP);
    int N = (int)p.value("N", 8);
    double ln2 = 0.6931471805599453;
    double gD = p.value("gamaD", 0.02);

    QVector<int> idx(nP);
    std::iota(idx.begin(), idx.end(), 0);

    auto calcP = [&](int k) {
        double t = tD[k];
        if (t <= 1e-10) { oPD[k] = 0.0; return; }

        double pd_val = 0.0;
        for (int m = 1; m <= N; ++m) {
            double pf = f(m * ln2 / t, p);
            if (std::isnan(pf) || std::isinf(pf)) pf = 0.0;
            pd_val += getStehfestCoeff(m, N) * pf;
        }
        double pd_r = pd_val * ln2 / t;

        // 压敏切线外推，处理极小储层能量时的对数越界崩溃
        if (gD > 1e-9) {
            double arg = 1.0 - gD * pd_r;
            double a_min = 1e-3;
            if (arg >= a_min) pd_r = -1.0 / gD * std::log(arg);
            else {
                double val_at_min = -1.0 / gD * std::log(a_min);
                double slope_at_min = -1.0 / (gD * a_min);
                pd_r = val_at_min + slope_at_min * (arg - a_min);
            }
        }
        if (pd_r <= 1e-15) pd_r = 1e-15;
        oPD[k] = pd_r;
    };
    QtConcurrent::blockingMap(idx, calcP);
}

double ModelSolver07::flaplace_composite(double z, const QMap<QString, double>& p) {
    int r_t = ((int)m_type / 12) + 1;

    double kf = p.value("kf", 50.0);
    double k2 = p.value("k2", 10.0);
    if (k2 < 1e-12) k2 = 1e-12;
    // 导压系数比直接由物性推导得出，屏蔽人工主观干扰
    double M12 = kf / k2;
    double eta12 = M12;

    double L = p.value("L_half", 500);
    double rm = p.value("rm", 1500.0);
    double re = p.value("re", 20000.0);
    double rmD = (L > 1e-9) ? rm / L : 1.25;
    double reD = (L > 1e-9) ? re / L : 25.0;

    int nF = (int)p.value("nf", 9);

    // 【重要重构】：解析独立向量 Lf_vec, xw_vec 和 表皮 S_vec
    QVector<double> LfD_vec(nF, 0.1);
    QVector<double> xwD(nF, 0.0);
    QVector<double> S_vec(nF, 0.0);

    double Lf_def = p.value("Lf", 50.0);
    double S_def = p.value("S", 1.0);

    for(int i = 0; i < nF; ++i) {
        double Lfi = p.value(QString("Lf_%1").arg(i), Lf_def);
        LfD_vec[i] = (L > 1e-9) ? (Lfi / L) : 0.1;

        // 基于跟端距离 (xw) 计算积分网格绝对坐标 (偏移 L，除以 L)
        double fallback_pos = (0.1 + 0.8 * (double)i / (double)std::max(1, nF - 1)) * (2.0 * L);
        double xwi = p.value(QString("xw_%1").arg(i), fallback_pos);
        xwD[i] = (xwi - L) / L;

        // 读取该段专属的表皮系数
        S_vec[i] = p.value(QString("S_%1").arg(i), S_def);
        if (S_vec[i] < 0.0) S_vec[i] = 0.0;
    }

    // 内外区流场组合分配
    double fs1 = 1.0, fs2 = 1.0;
    int id = (int)m_type + 1;

    if (id <= 24) {
        double om1 = p.value("omega1", 0.4);
        double rm1 = p.value("lambda1", 1e-3);
        double o_m1 = 1.0 - om1;
        fs1 = (om1 * o_m1 * z + rm1) / (o_m1 * z + rm1);
    }
    if (id <= 12) {
        double om2 = p.value("omega2", 0.08);
        double rm2 = p.value("lambda2", 1e-4);
        double o_m2 = 1.0 - om2;
        fs2 = eta12 * (om2 * o_m2 * eta12 * z + rm2) / (o_m2 * eta12 * z + rm2);
    } else {
        fs2 = eta12;
    }

    double CD = p.value("cD", 0.1);
    double alpha = p.value("alpha", 1e-1);
    double C_p = p.value("C_phi", 1e-4);

    // 将独立的 S_vec 透传进干涉矩阵中进行计算
    return PWD_composite(z, fs1, fs2, M12, LfD_vec, xwD, rmD, reD, nF, m_type, S_vec, CD, C_p, alpha);
}

// 【核心模块】包含独立表皮向量的 2nf+1 维封闭耦合求解
double ModelSolver07::PWD_composite(double z, double fs1, double fs2, double M12,
                                    const QVector<double>& LfD_vec, const QVector<double>& xwD,
                                    double rmD, double reD, int n_fracs, ModelType type,
                                    const QVector<double>& S_vec, double CD, double C_phi, double alpha) {
    int gI = ((int)type) % 12;
    bool isI = (gI < 4);
    bool isC = (gI >= 4 && gI < 8);
    bool isP = (gI >= 8);

    double ga1 = std::sqrt(z * fs1);
    double ga2 = std::sqrt(z * fs2);
    double a_g1_r = ga1 * rmD;
    double a_g2_r = ga2 * rmD;

    double k0g1 = safe_bessel_k_sc(0, a_g1_r);
    double k1g1 = safe_bessel_k_sc(1, a_g1_r);
    double i0g1 = safe_bessel_i_sc(0, a_g1_r);
    double i1g1 = safe_bessel_i_sc(1, a_g1_r);

    double k0g2 = safe_bessel_k_sc(0, a_g2_r);
    double k1g2 = safe_bessel_k_sc(1, a_g2_r);
    double i0g2 = safe_bessel_i_sc(0, a_g2_r);
    double i1g2 = safe_bessel_i_sc(1, a_g2_r);

    double T1p = k0g2;
    double T2p = -k1g2;

    if (!isI && reD > 1e-5) {
        double a_re = ga2 * reD;
        double k0re = safe_bessel_k_sc(0, a_re);
        double k1re = safe_bessel_k_sc(1, a_re);
        double i0re = safe_bessel_i_sc(0, a_re);
        double i1re = safe_bessel_i_sc(1, a_re);
        double ef = std::exp(2.0 * ga2 * (rmD - reD));

        if (isC) {
            double r = k1re / std::max(i1re, 1e-100);
            T1p = r * i0g2 * ef + k0g2;
            T2p = r * i1g2 * ef - k1g2;
        } else if (isP) {
            double r = -k0re / std::max(i0re, 1e-100);
            T1p = r * i0g2 * ef + k0g2;
            T2p = r * i1g2 * ef - k1g2;
        }
    }

    double Acu = M12 * ga1 * k1g1 * T1p + ga2 * k0g1 * T2p;
    double Acd = M12 * ga1 * i1g1 * T1p - ga2 * i0g1 * T2p;
    if (std::abs(Acd) < 1e-100) Acd = (Acd >= 0) ? 1e-100 : -1e-100;
    double Ac_c = Acu / Acd;

    // ========== 1. 纯地层缝间空间流场积分干涉矩阵 ==========
    Eigen::MatrixXd pfD_mat = Eigen::MatrixXd::Zero(n_fracs, n_fracs);
    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) {
            double dx = std::abs(xwD[i] - xwD[j]);
            double LfDj = LfD_vec[j];
            double val = 0.0;

            if (i == j) {
                double LfDi = LfD_vec[i];
                double x_o = 0.732 * LfDi;
                auto fn = [&](double a) -> double {
                    double dt = std::abs(x_o - a);
                    if (dt < 1e-15) return 0.0;
                    double a_d = ga1 * dt;
                    return safe_bessel_k(0, a_d) + Ac_c * safe_bessel_i_sc(0, a_d) * std::exp(a_d - 2.0 * a_g1_r);
                };
                double err;
                double I_l = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(fn, -LfDi, x_o, 15, 1e-6, &err);
                double I_r = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(fn, x_o, LfDi, 15, 1e-6, &err);
                val = (I_l + I_r) / (z * 2.0 * LfDi);
            } else {
                auto fn = [&](double a) -> double {
                    double a_d = ga1 * std::sqrt(dx * dx + a * a);
                    return safe_bessel_k(0, a_d) + Ac_c * safe_bessel_i_sc(0, a_d) * std::exp(a_d - 2.0 * a_g1_r);
                };
                double err;
                double I_f = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(fn, -LfDj, LfDj, 15, 1e-6, &err);
                val = I_f / (z * 2.0 * LfDj);
            }
            pfD_mat(i, j) = val;
        }
    }

    // ========== 2. 井筒与独立表皮流场条件映射 ==========
    int storageType = ((int)type) % 4;
    double Ceff = CD;
    QVector<double> S_eff_vec(n_fracs, 0.0);

    // 精确实现 MATLAB "get_effective_wellbore_terms_v8_1"
    if (storageType == 1) {
        Ceff = CD;
        S_eff_vec = S_vec;
    } else if (storageType == 2) {
        Ceff = 0.0;
        S_eff_vec.fill(0.0); // 线源解：退化所有井筒污染效应
    } else if (storageType == 3) {
        double P_phiD = (std::abs(alpha) > 1e-12) ? (C_phi / (z * (1.0 + alpha * z))) : 0.0;
        Ceff = CD * (1.0 + z * z * P_phiD);
        S_eff_vec = S_vec;
    } else if (storageType == 4) {
        double P_phiD = 0.0;
        if (std::abs(alpha) > 1e-12) {
            double xx = z * alpha / 2.0;
            double exp_erfc = (xx < 10.0) ? (std::exp(xx * xx) * std::erfc(xx)) : (1.0 / (std::sqrt(M_PI) * xx));
            P_phiD = (C_phi / z) * exp_erfc;
        }
        Ceff = CD * (1.0 + z * z * P_phiD);
        S_eff_vec = S_vec;
    }

    // ========== 3. 构造非对称阻力 2nf+1 大矩阵系统 ==========
    int sys_size = 2 * n_fracs + 1;
    Eigen::MatrixXd A_sys = Eigen::MatrixXd::Zero(sys_size, sys_size);
    Eigen::VectorXd b_sys = Eigen::VectorXd::Zero(sys_size);

    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) A_sys(i, j) = z * pfD_mat(i, j);
        A_sys(i, n_fracs + i) = -1.0;
    }

    // 【关键突破点】此处的流动压降施加，提取自 S_eff_vec[i]
    for (int i = 0; i < n_fracs; ++i) {
        int row = n_fracs + i;
        A_sys(row, i) = -S_eff_vec[i];
        A_sys(row, n_fracs + i) = -1.0;
        A_sys(row, 2 * n_fracs) = 1.0;
    }

    for (int j = 0; j < n_fracs; ++j) A_sys(2 * n_fracs, j) = 1.0;
    A_sys(2 * n_fracs, 2 * n_fracs) = Ceff * z;
    b_sys(2 * n_fracs) = 1.0 / z;

    // 自适应吉洪诺夫正则化
    double max_norm = A_sys.cwiseAbs().colwise().sum().maxCoeff();
    double lam_reg = std::max(1e-12, 1e-10 * max_norm);
    A_sys += lam_reg * Eigen::MatrixXd::Identity(sys_size, sys_size);

    Eigen::VectorXd x_sol = A_sys.fullPivLu().solve(b_sys);
    return x_sol(2 * n_fracs);
}

void ModelSolver07::precomputeStehfestCoeffs(int N) {
    if (m_currentN == N && !m_stehfestCoeffs.isEmpty()) return;
    m_currentN = N;
    m_stehfestCoeffs.resize(N + 1);
    for (int i = 1; i <= N; ++i) {
        double s = 0.0;
        for (int k = (i + 1) / 2; k <= std::min(i, N / 2); ++k) {
            double den = factorial(N / 2 - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i);
            if (den != 0) s += (std::pow((double)k, N / 2.0) * factorial(2 * k)) / den;
        }
        m_stehfestCoeffs[i] = (((i + N / 2) % 2 == 0) ? 1.0 : -1.0) * s;
    }
}

double ModelSolver07::getStehfestCoeff(int i, int N) {
    if (m_currentN != N || i < 1 || i > N) return 0.0;
    return m_stehfestCoeffs[i];
}

double ModelSolver07::factorial(int n) {
    double r = 1.0;
    for(int i = 2; i <= n; ++i) r *= i;
    return r;
}
