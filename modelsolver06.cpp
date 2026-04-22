/*
 * 文件名: modelsolver06.cpp
 * 功能描述:
 * 1. 模型 6 混积型非均质裂缝的物理闭合计算核心。
 * 2. 积分域严格执行 `[-LfDj, LfDj]` 源裂缝长度限制，观察点精准绑定为 `0.732 * LfDi`。
 * 3. 严格使用 `phi_ct` 合并规避尺度漂移，并引入 16.118095 压力换算高精度标度。
 * 4. 自带矩阵条件数检测与自适应吉洪诺夫微小对角线正则化。
 */

#include "modelsolver06.h"
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
double ModelSolver06::safe_bessel_k(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    try { return boost::math::cyl_bessel_k(v, x); } catch (...) { return 0.0; }
}

double ModelSolver06::safe_bessel_k_sc(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    if (x > 600.0) return std::sqrt(M_PI / (2.0 * x));
    try { return boost::math::cyl_bessel_k(v, x) * std::exp(x); } catch (...) { return 0.0; }
}

double ModelSolver06::safe_bessel_i_sc(int v, double x) {
    if (x < 0) x = -x;
    if (x > 600.0) return 1.0 / std::sqrt(2.0 * M_PI * x);
    try { return boost::math::cyl_bessel_i(v, x) * std::exp(-x); } catch (...) { return 0.0; }
}

ModelSolver06::ModelSolver06(ModelType type) :
    m_type(type), m_highPrecision(true), m_currentN(0) {
    precomputeStehfestCoeffs(8);
}

ModelSolver06::~ModelSolver06() {}

void ModelSolver06::setHighPrecision(bool h) { m_highPrecision = h; }

QString ModelSolver06::getModelName(ModelType t, bool v) {
    int id = (int)t + 1;
    QString bN = QString("不等长非均匀三重孔隙混积型模型%1").arg(id);
    return bN;
}

QVector<double> ModelSolver06::generateLogTimeSteps(int c, double s, double e) {
    QVector<double> t;
    if (c <= 0) return t;
    t.reserve(c);
    for (int i = 0; i < c; ++i) {
        t.append(std::pow(10.0, s + (e - s) * i / (c - 1)));
    }
    return t;
}

// 主时间循环入口，融合高精度无因次参数
ModelCurveData ModelSolver06::calculateTheoreticalCurve(const QMap<QString, double>& p, const QVector<double>& pt) {
    QVector<double> tP = pt.isEmpty() ? generateLogTimeSteps(100, -3.0, 4.0) : pt;

    double phi = p.value("phi", 0.05);
    double mu = p.value("mu", 5.0);
    double B = p.value("B", 1.2);
    double Ct = p.value("Ct", 5e-3);
    double q = p.value("q", 10.0);
    double h = p.value("h", 10.0);
    double kf = p.value("kf", 10.0);
    double L = p.value("L", 1000.0) / 2.0;

    if (L < 1e-9) L = 500.0;
    if (phi < 1e-12 || mu < 1e-12 || Ct < 1e-12 || kf < 1e-12) {
        return std::make_tuple(tP, QVector<double>(tP.size(), 0.0), QVector<double>(tP.size(), 0.0));
    }

    // 使用合并因子控制精度衰减
    double phi_ct = phi * Ct;
    double td_c = 3.6 * kf / (phi_ct * mu * std::pow(L, 2.0));
    double p_c = 16.118095 * q * mu * B / (kf * h);

    QVector<double> tD_v;
    tD_v.reserve(tP.size());
    for(double t : tP) tD_v.append(td_c * t);

    QMap<QString, double> cP = p;
    int N = (int)cP.value("N", 8);
    cP["N"] = N;
    precomputeStehfestCoeffs(N);

    if (!cP.contains("nf") || cP["nf"] < 1) cP["nf"] = 10;
    cP["L_half"] = L;

    QVector<double> PD_v, Der_v;
    auto func = std::bind(&ModelSolver06::flaplace_composite, this, std::placeholders::_1, std::placeholders::_2);
    calculatePDandDeriv(tD_v, cP, func, PD_v, Der_v);

    QVector<double> fP(tP.size()), fDP(tP.size());
    for(int i = 0; i < tP.size(); ++i) fP[i] = p_c * PD_v[i];

    if (tP.size() > 2) fDP = PressureDerivativeCalculator::calculateBourdetDerivative(tP, fP, 0.2);
    else fDP.fill(0.0);

    return std::make_tuple(tP, fP, fDP);
}

// 并发 Stehfest 反演与压敏保护
void ModelSolver06::calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& p,
                                        std::function<double(double, const QMap<QString, double>&)> f,
                                        QVector<double>& oPD, QVector<double>& oDer)
{
    int nP = tD.size();
    oPD.resize(nP); oDer.resize(nP);
    int N = (int)p.value("N", 8);
    double ln2 = 0.6931471805599453;
    double gD = p.value("gamaD", 0.006);

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

        // 压敏非线性扰动：一阶泰勒展开外推，防止对数为负溢出
        if (gD > 1e-9) {
            double arg = 1.0 - gD * pd_r;
            double a_min = 1e-3;
            if (arg >= a_min) {
                pd_r = -1.0 / gD * std::log(arg);
            } else {
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

// 页岩 Tanh 特征函数
double ModelSolver06::calc_fs_shale(double u, double o, double l) {
    if (u < 1e-15) return 1.0;
    double om = std::max(1.0 - o, 1e-9);
    double lbd = std::max(l, 1e-15);
    double isq = 3.0 * om * u / lbd;
    return o + std::sqrt((lbd * om) / (3.0 * u)) * std::tanh(std::sqrt(isq));
}

// 复合模型特征函数构建
double ModelSolver06::flaplace_composite(double z, const QMap<QString, double>& p) {
    int r_t = ((int)m_type / 12) + 1;

    double kf = p.value("kf", 50.0);
    double k2 = p.value("k2", 10.0);
    if (k2 < 1e-12) k2 = 1e-12;
    double M12 = kf / k2;
    double eta12 = M12;

    double L = p.value("L_half", 500);
    double rm = p.value("rm", 5000);
    double re = p.value("re", 5e5);
    double rmD = (L > 1e-9) ? rm / L : 16.0;
    double reD = (L > 1e-9) ? re / L : 1000.0;

    int nF = (int)p.value("nf", 10);

    // 【向量化提取】解析各缝长度与非均匀中心位置
    QVector<double> LfD_vec(nF, 0.1);
    QVector<double> xwD(nF, 0.0);
    double Lf_default = p.value("Lf", 50.0);

    for(int i = 0; i < nF; ++i) {
        double Lfi = p.value(QString("Lf_%1").arg(i), Lf_default);
        LfD_vec[i] = (L > 1e-9) ? (Lfi / L) : 0.1;

        double fallback_pos = (0.1 + 0.8 * (double)i / (double)std::max(1, nF - 1)) * (2.0 * L);
        double xwi = p.value(QString("xw_%1").arg(i), fallback_pos);
        xwD[i] = (xwi - L) / L;
    }

    // 三重孔隙内区特征函数组装
    double of1 = p.value("omega_f1", 0.02);
    double ov1 = p.value("omega_v1", 0.01);
    double rm1 = p.value("lambda_m1", 1e-4);
    double rv1 = p.value("lambda_v1", 1e-1);

    double om1 = 1.0 - of1 - ov1;
    double dm1 = om1 * z + rm1;
    double dv1 = ov1 * z + rv1;
    double fs1 = of1 + (dm1 == 0 ? 0 : (rm1 * om1) / dm1) + (dv1 == 0 ? 0 : (rv1 * ov1) / dv1);

    // 外区特征函数
    double of2 = p.value("omega_f2", 0.008);
    double rm2 = p.value("lambda_m2", 1e-7);
    double fs2 = eta12;

    if (r_t == 1) { // 夹层
        double om2 = 1.0 - of2;
        double dm2 = om2 * eta12 * z + rm2;
        fs2 = eta12 * (of2 + (dm2 == 0 ? 0 : (rm2 * om2) / dm2));
    } else if (r_t == 2) { // 页岩
        fs2 = eta12 * calc_fs_shale(eta12 * z, of2, rm2);
    }

    double CD = p.value("cD", 0.1);
    double S = p.value("S", 10.0);
    if (S < 0.0) S = 0.0;
    double alpha = p.value("alpha", 0.1);
    double C_p = p.value("C_phi", 1e-4);

    return PWD_composite(z, fs1, fs2, M12, LfD_vec, xwD, rmD, reD, nF, m_type, S, CD, C_p, alpha);
}

// 【核心模块】2nf+1 维干扰大矩阵耦合装配与求逆
double ModelSolver06::PWD_composite(double z, double fs1, double fs2, double M12,
                                    const QVector<double>& LfD_vec, const QVector<double>& xwD,
                                    double rmD, double reD, int n_fracs, ModelType type,
                                    double S, double CD, double C_phi, double alpha) {
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

    // ========== 1. 构建流场干扰基底矩阵 ==========
    Eigen::MatrixXd pfD_mat = Eigen::MatrixXd::Zero(n_fracs, n_fracs);
    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) {
            double dx = std::abs(xwD[i] - xwD[j]);
            double LfDj = LfD_vec[j]; // 积分微元界限限定为作用源裂缝半长
            double val = 0.0;

            if (i == j) {
                double LfDi = LfD_vec[i];
                double x_o = 0.732 * LfDi; // 自身干涉等效观察点基于自身长度
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

    // ========== 2. 井筒闭合属性映射 ==========
    int storageType = ((int)type) % 4;
    double Ceff = CD;
    double S_eff = S;

    if (storageType == 1) {
        Ceff = 0.0; S_eff = 0.0;
    } else if (storageType == 2) {
        double P_phiD = (std::abs(alpha) > 1e-12) ? (C_phi / (z * (1.0 + alpha * z))) : 0.0;
        Ceff = CD * (1.0 + z * z * P_phiD);
    } else if (storageType == 3) {
        double P_phiD = 0.0;
        if (std::abs(alpha) > 1e-12) {
            double xx = z * alpha / 2.0;
            double exp_erfc = (xx < 10.0) ? (std::exp(xx * xx) * std::erfc(xx)) : (1.0 / (std::sqrt(M_PI) * xx));
            P_phiD = (C_phi / z) * exp_erfc;
        }
        Ceff = CD * (1.0 + z * z * P_phiD);
    }

    // ========== 3. 构建 2nf+1 耦合大矩阵 ==========
    int sys_size = 2 * n_fracs + 1;
    Eigen::MatrixXd A_sys = Eigen::MatrixXd::Zero(sys_size, sys_size);
    Eigen::VectorXd b_sys = Eigen::VectorXd::Zero(sys_size);

    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) A_sys(i, j) = z * pfD_mat(i, j);
        A_sys(i, n_fracs + i) = -1.0;
    }

    for (int i = 0; i < n_fracs; ++i) {
        int row = n_fracs + i;
        A_sys(row, i) = -S_eff;
        A_sys(row, n_fracs + i) = -1.0;
        A_sys(row, 2 * n_fracs) = 1.0;
    }

    for (int j = 0; j < n_fracs; ++j) A_sys(2 * n_fracs, j) = 1.0;
    A_sys(2 * n_fracs, 2 * n_fracs) = Ceff * z;
    b_sys(2 * n_fracs) = 1.0 / z;

    double max_norm = A_sys.cwiseAbs().colwise().sum().maxCoeff();
    double lam_reg = std::max(1e-12, 1e-10 * max_norm);
    A_sys += lam_reg * Eigen::MatrixXd::Identity(sys_size, sys_size);

    Eigen::VectorXd x_sol = A_sys.fullPivLu().solve(b_sys);
    return x_sol(2 * n_fracs);
}

void ModelSolver06::precomputeStehfestCoeffs(int N) {
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

double ModelSolver06::getStehfestCoeff(int i, int N) {
    if (m_currentN != N || i < 1 || i > N) return 0.0;
    return m_stehfestCoeffs[i];
}

double ModelSolver06::factorial(int n) {
    double r = 1.0;
    for(int i = 2; i <= n; ++i) r *= i;
    return r;
}
