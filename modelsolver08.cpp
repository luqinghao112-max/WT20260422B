/*
 * 文件名: modelsolver08.cpp
 * 功能描述:
 * 1. 模型 8 (页岩型非均质独立表皮模型) 的核心物理计算实现。
 * 2. 精准映射: 每条裂缝具有独立长度 Lf_vec 与独立位置 xw，且具有专属局部表皮 S_vec。
 * 3. 严格使用 16.118095 进行真实压力换算，利用 phi_ct 合并储容控制量纲转换精度。
 * 4. PWD_composite 全面使用各段独立表皮 S_vec 进行局部流场压降的闭合计算。
 */

#include "modelsolver08.h"
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

// ================= 1. 底层安全贝塞尔函数实现 =================
double ModelSolver08::safe_bessel_k(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    try { return boost::math::cyl_bessel_k(v, x); } catch (...) { return 0.0; }
}

double ModelSolver08::safe_bessel_k_scaled(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    if (x > 600.0) return std::sqrt(M_PI / (2.0 * x)); // 渐近展开防止大自变量崩溃
    try { return boost::math::cyl_bessel_k(v, x) * std::exp(x); } catch (...) { return 0.0; }
}

double ModelSolver08::safe_bessel_i_scaled(int v, double x) {
    if (x < 0) x = -x;
    if (x > 600.0) return 1.0 / std::sqrt(2.0 * M_PI * x); // 渐近展开防止大自变量崩溃
    try { return boost::math::cyl_bessel_i(v, x) * std::exp(-x); } catch (...) { return 0.0; }
}

// 构造函数：初始化 N=8 阶 Stehfest 系数
ModelSolver08::ModelSolver08(ModelType type)
    : m_type(type), m_highPrecision(true), m_currentN(0) {
    precomputeStehfestCoeffs(8);
}

ModelSolver08::~ModelSolver08() {}

void ModelSolver08::setHighPrecision(bool high) { m_highPrecision = high; }

QString ModelSolver08::getModelName(ModelType type, bool verbose) {
    int id = (int)type + 1;
    return QString("独立表皮非均质页岩型模型%1").arg(id);
}

QVector<double> ModelSolver08::generateLogTimeSteps(int count, double startExp, double endExp) {
    QVector<double> t;
    if (count <= 0) return t;
    t.reserve(count);
    for (int i = 0; i < count; ++i) {
        t.append(std::pow(10.0, startExp + (endExp - startExp) * i / (count - 1)));
    }
    return t;
}

// ================= 2. 理论曲线计算主入口 =================
ModelCurveData ModelSolver08::calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime) {
    QVector<double> tPoints = providedTime;
    if (tPoints.isEmpty()) tPoints = generateLogTimeSteps(100, -3.0, 4.0);

    // 提取宏观地层物性参数
    double phi = params.value("phi", 0.05);
    double mu = params.value("mu", 0.5);
    double B = params.value("B", 1.2);
    double Ct = params.value("Ct", 5e-4);
    double q = params.value("q", 50.0);
    double h = params.value("h", 20.0);
    double kf = params.value("kf", 50.0);
    double L = params.value("L", 1000.0) / 2.0;

    if (L < 1e-9) L = 500.0;
    if (phi < 1e-12 || mu < 1e-12 || Ct < 1e-12 || kf < 1e-12) {
        return std::make_tuple(tPoints, QVector<double>(tPoints.size(), 0.0), QVector<double>(tPoints.size(), 0.0));
    }

    // 计算高精度无因次转换系数（综合储容系数化简处理）
    double phi_ct = phi * Ct;
    double td_coeff = 3.6 * kf / (phi_ct * mu * std::pow(L, 2.0));
    double p_coeff = 16.118095 * q * mu * B / (kf * h);

    QVector<double> tD_vec;
    tD_vec.reserve(tPoints.size());
    for(double t : tPoints) tD_vec.append(td_coeff * t);

    QMap<QString, double> calcParams = params;
    int N = (int)calcParams.value("N", 8);
    calcParams["N"] = N;
    precomputeStehfestCoeffs(N);

    if (!calcParams.contains("nf") || calcParams["nf"] < 1) calcParams["nf"] = 9;
    calcParams["L_half"] = L;

    // 并发调度拉氏逆变换
    QVector<double> PD_vec, Deriv_vec;
    auto func = std::bind(&ModelSolver08::flaplace_composite, this, std::placeholders::_1, std::placeholders::_2);
    calculatePDandDeriv(tD_vec, calcParams, func, PD_vec, Deriv_vec);

    // 转换至真实压力并计算 Bourdet 导数
    QVector<double> finalP(tPoints.size()), finalDP(tPoints.size());
    for(int i = 0; i < tPoints.size(); ++i) finalP[i] = p_coeff * PD_vec[i];

    if (tPoints.size() > 2) finalDP = PressureDerivativeCalculator::calculateBourdetDerivative(tPoints, finalP, 0.2);
    else finalDP.fill(0.0);

    return std::make_tuple(tPoints, finalP, finalDP);
}

// ================= 3. 时域并发数值转换与压敏切线外推 =================
void ModelSolver08::calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                                        std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                                        QVector<double>& outPD, QVector<double>& outDeriv)
{
    int numPoints = tD.size();
    outPD.resize(numPoints); outDeriv.resize(numPoints);
    int N = (int)params.value("N", 8);
    double ln2 = 0.6931471805599453;
    double gamaD = params.value("gamaD", 0.02); // 渗透率模量

    QVector<int> indexes(numPoints);
    std::iota(indexes.begin(), indexes.end(), 0);

    auto calculateSinglePoint = [&](int k) {
        double t = tD[k];
        if (t <= 1e-10) { outPD[k] = 0.0; return; }

        double pd_val = 0.0;
        for (int m = 1; m <= N; ++m) {
            double z = m * ln2 / t;
            double pf = laplaceFunc(z, params);
            if (std::isnan(pf) || std::isinf(pf)) pf = 0.0;
            pd_val += getStehfestCoeff(m, N) * pf;
        }
        double pd_real = pd_val * ln2 / t;

        // 引入压敏非线性摄动的泰勒保护机制
        if (gamaD > 1e-9) {
            double arg = 1.0 - gamaD * pd_real;
            double arg_min = 1e-3;
            if (arg >= arg_min) pd_real = -1.0 / gamaD * std::log(arg);
            else {
                // 当储层能量极低导致对数自变量崩溃时，切换至切线平滑外推
                double val_at_min = -1.0 / gamaD * std::log(arg_min);
                double slope_at_min = -1.0 / (gamaD * arg_min);
                pd_real = val_at_min + slope_at_min * (arg - arg_min);
            }
        }

        if (pd_real <= 1e-15) pd_real = 1e-15;
        outPD[k] = pd_real;
    };
    QtConcurrent::blockingMap(indexes, calculateSinglePoint);
}

// 外区双孔模型伪稳态窜流特征函数
double ModelSolver08::calc_fs_dual(double u, double omega, double lambda) {
    double one_minus = 1.0 - omega;
    double den = one_minus * u + lambda;
    if (std::abs(den) < 1e-20) return 0.0;
    return (omega * one_minus * u + lambda) / den;
}

// 页岩型基质瞬态扩散特征函数 (包含 Tanh 非线性核)
double ModelSolver08::calc_fs_shale(double u, double omega, double lambda) {
    if (u < 1e-15) return 1.0;
    double one_minus = 1.0 - omega;
    if (one_minus < 1e-9) one_minus = 1e-9;
    if (lambda < 1e-15) lambda = 1e-15;
    double inside_sqrt = 3.0 * one_minus * u / lambda;
    double arg_tanh = std::sqrt(inside_sqrt);
    double front_sqrt = std::sqrt( (lambda * one_minus) / (3.0 * u) );
    return omega + front_sqrt * std::tanh(arg_tanh);
}

// ================= 4. 拉普拉斯空间分配与独立特征提取 =================
double ModelSolver08::flaplace_composite(double z, const QMap<QString, double>& p) {
    double kf = p.value("kf", 50.0);
    double k2 = p.value("k2", 10.0);
    if (k2 < 1e-12) k2 = 1e-12;

    double M12 = kf / k2; // 导压系数比直接计算，保证理论正确性
    double eta12 = M12;

    double L = p.value("L_half", 500.0);
    double rm = p.value("rm", 1500.0);
    double re = p.value("re", 20000.0);

    double rmD = (L > 1e-9) ? rm / L : 1.25;
    double reD = (L > 1e-9) ? re / L : 25.0;

    int n_fracs = (int)p.value("nf", 9);

    // 建立存储多缝信息的物理属性向量
    QVector<double> LfD_vec(n_fracs, 0.1);
    QVector<double> xwD(n_fracs, 0.0);
    QVector<double> S_vec(n_fracs, 0.0); // 独立表皮系数向量

    double Lf_default = p.value("Lf", 50.0);
    double S_default = p.value("S", 1.0);

    // 逐段提取非均质裂缝参数
    for(int i = 0; i < n_fracs; ++i) {
        // 读取半长 (缺省使用统一的 Lf 托底)
        double Lfi = p.value(QString("Lf_%1").arg(i), Lf_default);
        LfD_vec[i] = (L > 1e-9) ? (Lfi / L) : 0.1;

        // 读取裂缝跟端距离，转换为中心原点的绝对无因次坐标
        double fallback_pos = (0.1 + 0.8 * (double)i / (double)std::max(1, n_fracs - 1)) * (2.0 * L);
        double xwi = p.value(QString("xw_%1").arg(i), fallback_pos);
        xwD[i] = (xwi - L) / L;

        // 提取每段对应的局部表皮损伤
        S_vec[i] = p.value(QString("S_%1").arg(i), S_default);
        if (S_vec[i] < 0.0) S_vec[i] = 0.0;
    }

    // 内区强制应用页岩 Tanh 窜流特征
    double omga1 = p.value("omega1", 0.4);
    double remda1 = p.value("lambda1", 1e-3);
    double fs1 = calc_fs_shale(z, omga1, remda1);

    double z_outer = eta12 * z;
    int id = (int)m_type + 1;
    double fs2;

    // 外区根据模型类型决定特征函数组合
    if (id <= 12) {
        double omga2 = p.value("omega2", 0.08);
        double remda2 = p.value("lambda2", 1e-4);
        fs2 = eta12 * calc_fs_shale(z_outer, omga2, remda2);
    } else if (id <= 24) {
        fs2 = eta12;
    } else {
        double omga2 = p.value("omega2", 0.08);
        double remda2 = p.value("lambda2", 1e-4);
        fs2 = eta12 * calc_fs_dual(z_outer, omga2, remda2);
    }

    double CD = p.value("cD", 0.1);
    double alpha = p.value("alpha", 1e-1);
    double C_phi = p.value("C_phi", 1e-4);

    // 将独立的 S_vec透传进干扰大矩阵
    return PWD_composite(z, fs1, fs2, M12, LfD_vec, xwD, rmD, reD, n_fracs, m_type, S_vec, CD, C_phi, alpha);
}

// ================= 5. 【核心模块】包含独立表皮向量的 2nf+1 维封闭耦合求解 =================
double ModelSolver08::PWD_composite(double z, double fs1, double fs2, double M12,
                                    const QVector<double>& LfD_vec, const QVector<double>& xwD,
                                    double rmD, double reD, int n_fracs, ModelType type,
                                    const QVector<double>& S_vec, double CD, double C_phi, double alpha) {
    int id = (int)type + 1;
    int groupIdx = (id - 1) % 12;
    bool isInfinite = (groupIdx < 4);
    bool isClosed = (groupIdx >= 4 && groupIdx < 8);
    bool isConstP = (groupIdx >= 8);

    double gama1 = std::sqrt(z * fs1);
    double gama2 = std::sqrt(z * fs2);
    double arg_g1_rm = gama1 * rmD;
    double arg_g2_rm = gama2 * rmD;

    double k0_g1 = safe_bessel_k_scaled(0, arg_g1_rm);
    double k1_g1 = safe_bessel_k_scaled(1, arg_g1_rm);
    double i0_g1 = safe_bessel_i_scaled(0, arg_g1_rm);
    double i1_g1 = safe_bessel_i_scaled(1, arg_g1_rm);

    double k0_g2 = safe_bessel_k_scaled(0, arg_g2_rm);
    double k1_g2 = safe_bessel_k_scaled(1, arg_g2_rm);
    double i0_g2 = safe_bessel_i_scaled(0, arg_g2_rm);
    double i1_g2 = safe_bessel_i_scaled(1, arg_g2_rm);

    double T1_prime = k0_g2;
    double T2_prime = -k1_g2;

    if (!isInfinite && reD > 1e-5) {
        double arg_re = gama2 * reD;
        double k0_re = safe_bessel_k_scaled(0, arg_re);
        double k1_re = safe_bessel_k_scaled(1, arg_re);
        double i0_re = safe_bessel_i_scaled(0, arg_re);
        double i1_re = safe_bessel_i_scaled(1, arg_re);

        double exp_factor = std::exp(2.0 * gama2 * (rmD - reD));

        if (isClosed) {
            double ratio = k1_re / std::max(i1_re, 1e-100);
            T1_prime = ratio * i0_g2 * exp_factor + k0_g2;
            T2_prime = ratio * i1_g2 * exp_factor - k1_g2;
        } else if (isConstP) {
            double ratio = -k0_re / std::max(i0_re, 1e-100);
            T1_prime = ratio * i0_g2 * exp_factor + k0_g2;
            T2_prime = ratio * i1_g2 * exp_factor - k1_g2;
        }
    }

    double Acup_prime   = M12 * gama1 * k1_g1 * T1_prime + gama2 * k0_g1 * T2_prime;
    double Acdown_prime = M12 * gama1 * i1_g1 * T1_prime - gama2 * i0_g1 * T2_prime;
    if (std::abs(Acdown_prime) < 1e-100) Acdown_prime = (Acdown_prime >= 0) ? 1e-100 : -1e-100;

    double Ac_core = Acup_prime / Acdown_prime;

    // ----------------- (1) 装配空间流场干涉基底矩阵 -----------------
    Eigen::MatrixXd pfD_mat = Eigen::MatrixXd::Zero(n_fracs, n_fracs);
    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) {
            double dx = std::abs(xwD[i] - xwD[j]);
            double LfDj = LfD_vec[j];
            double val = 0.0;

            if (i == j) {
                double LfDi = LfD_vec[i];
                double x_obs = 0.732 * LfDi;
                auto integrand_self = [&](double a) -> double {
                    double dist = std::abs(x_obs - a);
                    if (dist < 1e-15) return 0.0;
                    double arg_dist = gama1 * dist;
                    return safe_bessel_k(0, arg_dist) + Ac_core * safe_bessel_i_scaled(0, arg_dist) * std::exp(arg_dist - 2.0 * arg_g1_rm);
                };

                unsigned max_depth = 15; double tol = 1e-6; double err_est;
                double I_left = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_self, -LfDi, x_obs, max_depth, tol, &err_est);
                double I_right = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_self, x_obs, LfDi, max_depth, tol, &err_est);
                val = (I_left + I_right) / (z * 2.0 * LfDi);
            } else {
                auto integrand = [&](double a) -> double {
                    double dist_val = std::sqrt(dx * dx + a * a);
                    double arg_dist = gama1 * dist_val;
                    return safe_bessel_k(0, arg_dist) + Ac_core * safe_bessel_i_scaled(0, arg_dist) * std::exp(arg_dist - 2.0 * arg_g1_rm);
                };

                unsigned max_depth = 15; double tol = 1e-6; double err_est;
                double I_full = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, -LfDj, LfDj, max_depth, tol, &err_est);
                val = I_full / (z * 2.0 * LfDj);
            }
            pfD_mat(i, j) = val;
        }
    }

    // ----------------- (2) 处理井筒与多段独立表皮向量的边界约束 -----------------
    int storageType = ((int)type) % 4;
    double Ceff = CD;
    QVector<double> S_eff_vec(n_fracs, 0.0);

    // 等效井筒效应重构：将多段表皮与井储效应分类处理
    if (storageType == 1) {
        Ceff = CD;
        S_eff_vec = S_vec;
    } else if (storageType == 2) {
        Ceff = 0.0;
        S_eff_vec.fill(0.0); // 线源解完全忽略表皮
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

    // ----------------- (3) 构建 2nf+1 耦合流场大矩阵 -----------------
    int sys_size = 2 * n_fracs + 1;
    Eigen::MatrixXd A_sys = Eigen::MatrixXd::Zero(sys_size, sys_size);
    Eigen::VectorXd b_sys = Eigen::VectorXd::Zero(sys_size);

    // 第一块：无阻力裂缝面压降映射
    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) {
            A_sys(i, j) = z * pfD_mat(i, j);
        }
        A_sys(i, n_fracs + i) = -1.0;
    }

    // 第二块：缝筒耦合（局部流动阻力施加）
    // 【关键改动】使用提取出的各段独立 S_eff_vec 施加到矩阵对角线上
    for (int i = 0; i < n_fracs; ++i) {
        int row = n_fracs + i;
        A_sys(row, i) = -S_eff_vec[i];
        A_sys(row, n_fracs + i) = -1.0;
        A_sys(row, 2 * n_fracs) = 1.0;
    }

    // 第三块：宏观流量守恒分配约束
    for (int j = 0; j < n_fracs; ++j) {
        A_sys(2 * n_fracs, j) = 1.0;
    }
    A_sys(2 * n_fracs, 2 * n_fracs) = Ceff * z;
    b_sys(2 * n_fracs) = 1.0 / z;

    // 应用正则化防护，针对密集裂缝分布防奇性
    double max_norm = A_sys.cwiseAbs().colwise().sum().maxCoeff();
    double lam_reg = std::max(1e-12, 1e-10 * max_norm);
    A_sys += lam_reg * Eigen::MatrixXd::Identity(sys_size, sys_size);

    Eigen::VectorXd x_sol = A_sys.fullPivLu().solve(b_sys);
    return x_sol(2 * n_fracs);
}

// 预计算 N=8 阶 Stehfest 拉氏逆变换常系数
void ModelSolver08::precomputeStehfestCoeffs(int N) {
    if (m_currentN == N && !m_stehfestCoeffs.isEmpty()) return;
    m_currentN = N; m_stehfestCoeffs.resize(N + 1);
    for (int i = 1; i <= N; ++i) {
        double s = 0.0;
        int k1 = (i + 1) / 2;
        int k2 = std::min(i, N / 2);
        for (int k = k1; k <= k2; ++k) {
            double num = std::pow((double)k, N / 2.0) * factorial(2 * k);
            double den = factorial(N / 2 - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i);
            if (den != 0) s += num / den;
        }
        double sign = ((i + N / 2) % 2 == 0) ? 1.0 : -1.0;
        m_stehfestCoeffs[i] = sign * s;
    }
}

// 提取 Stehfest 系数
double ModelSolver08::getStehfestCoeff(int i, int N) {
    if (m_currentN != N || i < 1 || i > N) return 0.0;
    return m_stehfestCoeffs[i];
}

// 阶乘辅助运算
double ModelSolver08::factorial(int n) {
    if(n <= 1) return 1.0;
    double r = 1.0;
    for(int i = 2; i <= n; ++i) r *= i;
    return r;
}
