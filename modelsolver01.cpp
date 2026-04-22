/*
 * 文件名: modelsolver01.cpp
 * 作用与功能:
 * 1. 夹层型储层 + 等长裂缝 算法核心。
 * 2. 【核心更新】采用 2nf+1 维线性方程组，将表皮、井储、各缝干涉流量与井底流压一并联立闭合求解。
 * 3. 强制使用 phi_ct 合并参数控制无因次时间缩放，确保数值计算在任意物理尺度下不丢精度。
 * 4. N=8 阶的高精度 Stehfest 变换，及完整的压敏效应泰勒切线外推。
 */

#include "modelsolver01.h"
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
double ModelSolver01::safe_bessel_k(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    try { return boost::math::cyl_bessel_k(v, x); } catch (...) { return 0.0; }
}

double ModelSolver01::safe_bessel_k_scaled(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    if (x > 600.0) return std::sqrt(M_PI / (2.0 * x)); // 渐进展开防溢出
    try { return boost::math::cyl_bessel_k(v, x) * std::exp(x); } catch (...) { return 0.0; }
}

double ModelSolver01::safe_bessel_i_scaled(int v, double x) {
    if (x < 0) x = -x;
    if (x > 600.0) return 1.0 / std::sqrt(2.0 * M_PI * x); // 渐进展开防溢出
    try { return boost::math::cyl_bessel_i(v, x) * std::exp(-x); } catch (...) { return 0.0; }
}

ModelSolver01::ModelSolver01(ModelType type)
    : m_type(type), m_highPrecision(true), m_currentN(0) {
    precomputeStehfestCoeffs(8); // N锁定为8阶，平衡精度与计算速度
}

ModelSolver01::~ModelSolver01() {}

void ModelSolver01::setHighPrecision(bool high) { m_highPrecision = high; }

QString ModelSolver01::getModelName(ModelType type, bool verbose) {
    int id = (int)type + 1;
    QString baseName;
    QString subType;

    if (id <= 12) {
        baseName = QString("夹层型储层试井解释模型%1").arg(id);
        subType = "夹层型+夹层型";
    } else if (id <= 24) {
        baseName = QString("夹层型储层试井解释模型%1").arg(id);
        subType = "夹层型+均质";
    } else {
        baseName = QString("径向复合模型%1").arg(id - 24);
        subType = "均质+均质";
    }

    if (!verbose) return baseName;

    int rem4 = (id - 1) % 4;
    QString strStorage;
    if (rem4 == 0) strStorage = "定井储";
    else if (rem4 == 1) strStorage = "线源解";
    else if (rem4 == 2) strStorage = "Fair模型";
    else strStorage = "Hegeman模型";

    int groupIdx = (id - 1) % 12;
    QString strBoundary;
    if (groupIdx < 4) strBoundary = "无限大外边界";
    else if (groupIdx < 8) strBoundary = "封闭边界";
    else strBoundary = "定压边界";

    return QString("%1\n(%2、%3、%4)").arg(baseName).arg(strStorage).arg(strBoundary).arg(subType);
}

QVector<double> ModelSolver01::generateLogTimeSteps(int count, double startExp, double endExp) {
    QVector<double> t;
    if (count <= 0) return t;
    t.reserve(count);
    for (int i = 0; i < count; ++i) {
        t.append(std::pow(10.0, startExp + (endExp - startExp) * i / (count - 1)));
    }
    return t;
}

ModelCurveData ModelSolver01::calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime) {
    QVector<double> tPoints = providedTime;
    if (tPoints.isEmpty()) tPoints = generateLogTimeSteps(100, -3.0, 4.0);

    // 获取并校验物理参数
    double phi = params.value("phi", 0.05);
    double mu = params.value("mu", 0.5);
    double B = params.value("B", 1.2);
    double Ct = params.value("Ct", 5e-4);
    double q = params.value("q", 50.0);
    double h = params.value("h", 20.0);
    double kf = params.value("kf", 50.0);
    double L_total = params.value("L", 1000.0);
    double L = L_total / 2.0;

    if (L < 1e-9) L = 500.0;
    if (phi < 1e-12 || mu < 1e-12 || Ct < 1e-12 || kf < 1e-12) {
        return std::make_tuple(tPoints, QVector<double>(tPoints.size(), 0.0), QVector<double>(tPoints.size(), 0.0));
    }

    // 【重要约束】：使用合并项 phi_ct 参与计算，维持无因次时间的尺度精度
    double phi_ct = phi * Ct;
    double td_coeff = 3.6 * kf / (phi_ct * mu * std::pow(L, 2.0));

    QVector<double> tD_vec;
    tD_vec.reserve(tPoints.size());
    for(double t : tPoints) tD_vec.append(td_coeff * t);

    QMap<QString, double> calcParams = params;
    int N = (int)calcParams.value("N", 8);
    if (N < 4 || N > 18 || N % 2 != 0) N = 8;
    calcParams["N"] = N;
    precomputeStehfestCoeffs(N);

    if (!calcParams.contains("nf") || calcParams["nf"] < 1) calcParams["nf"] = 9;
    calcParams["L_half"] = L;

    QVector<double> PD_vec, Deriv_vec;
    auto func = std::bind(&ModelSolver01::flaplace_composite, this, std::placeholders::_1, std::placeholders::_2);

    // QtConcurrent 并发执行拉氏逆变换
    calculatePDandDeriv(tD_vec, calcParams, func, PD_vec, Deriv_vec);

    double p_coeff = 1.842 * q * mu * B / (kf * h);
    QVector<double> finalP(tPoints.size()), finalDP(tPoints.size());

    for(int i = 0; i < tPoints.size(); ++i) finalP[i] = p_coeff * PD_vec[i];

    if (tPoints.size() > 2) finalDP = PressureDerivativeCalculator::calculateBourdetDerivative(tPoints, finalP, 0.2);
    else finalDP.fill(0.0);

    return std::make_tuple(tPoints, finalP, finalDP);
}

void ModelSolver01::calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                                        std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                                        QVector<double>& outPD, QVector<double>& outDeriv)
{
    int numPoints = tD.size();
    outPD.resize(numPoints); outDeriv.resize(numPoints);
    int N = (int)params.value("N", 8);
    double ln2 = 0.6931471805599453;
    double gamaD = params.value("gamaD", 0.02);

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

        // 【算法升级】：完整的压敏泰勒相切外推，防止对数运算在极端情况崩溃
        if (gamaD > 1e-9) {
            double arg = 1.0 - gamaD * pd_real;
            double arg_min = 1e-3;
            if (arg >= arg_min) {
                pd_real = -1.0 / gamaD * std::log(arg);
            } else {
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

// 复合模型在拉氏空间的组织分配入口
double ModelSolver01::flaplace_composite(double z, const QMap<QString, double>& p) {
    double kf = p.value("kf", 50.0);
    double k2 = p.value("k2", 10.0);
    if (k2 < 1e-12) k2 = 1e-12;

    double M12 = kf / k2;
    double eta12 = M12;

    double L = p.value("L_half", 500.0);
    double Lf = p.value("Lf", 50.0);
    double rm = p.value("rm", 1500.0);
    double re = p.value("re", 20000.0);

    double LfD = (L > 1e-9) ? Lf / L : 0.05;
    double rmD = (L > 1e-9) ? rm / L : 1.25;
    double reD = (L > 1e-9) ? re / L : 25.0;

    int n_fracs = (int)p.value("nf", 9);

    double fs1 = 1.0, fs2 = 1.0;
    int id = (int)m_type + 1;

    // 内外区介质特征函数分离
    if (id <= 24) {
        double omga1 = p.value("omega1", 0.4);
        double remda1 = p.value("lambda1", 1e-3);
        double one_minus = 1.0 - omga1;
        fs1 = (omga1 * one_minus * z + remda1) / (one_minus * z + remda1);
    }
    if (id <= 12) {
        double omga2 = p.value("omega2", 0.08);
        double remda2 = p.value("lambda2", 1e-4);
        double one_minus = 1.0 - omga2;
        fs2 = eta12 * (omga2 * one_minus * eta12 * z + remda2) / (one_minus * eta12 * z + remda2);
    } else {
        fs2 = eta12;
    }

    // 提取井筒附加参数，向下透传进行闭合联立
    double CD = p.value("cD", 0.1);
    double S = p.value("S", 1.0);
    if (S < 0.0) S = 0.0;
    double alpha = p.value("alpha", 1e-1);
    double C_phi = p.value("C_phi", 1e-4);

    return PWD_composite(z, fs1, fs2, M12, LfD, rmD, reD, n_fracs, m_type, S, CD, C_phi, alpha);
}

// 【重构核心】PWD_composite (核心数值求解器；2nf+1维井储+表皮闭合联立版)
double ModelSolver01::PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD,
                                    int n_fracs, ModelType type, double S, double CD, double C_phi, double alpha) {
    int id = (int)type + 1;
    int groupIdx = (id - 1) % 12;
    bool isInfinite = (groupIdx < 4);
    bool isClosed = (groupIdx >= 4 && groupIdx < 8);
    bool isConstP = (groupIdx >= 8);

    double gama1 = std::sqrt(z * fs1);
    double gama2 = std::sqrt(z * fs2);
    double arg_g1_rm = gama1 * rmD;
    double arg_g2_rm = gama2 * rmD;

    // 使用 Boost 库提供的修正贝塞尔函数，在 C++ 层面实现 MATLAB 的缩放贝塞尔以规避数值溢出
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

        double exp_factor = std::exp(2.0 * gama2 * (rmD - reD)); // 安全负指数

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

    // 裂缝在水平井段分布位置 (与 V8.2 严格等长等间距设定吻合)
    QVector<double> xwD(n_fracs);
    if (n_fracs == 1) {
        xwD[0] = 0.0; // 单缝居中
    } else {
        for (int k = 0; k < n_fracs; ++k) {
            xwD[k] = -0.9 + 1.8 * (double)k / (double)(n_fracs - 1);
        }
    }

    // ========== 1. 计算纯地层缝间干涉积分矩阵 ==========
    Eigen::MatrixXd pfD_mat = Eigen::MatrixXd::Zero(n_fracs, n_fracs);
    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) {
            double dx = std::abs(xwD[i] - xwD[j]);
            double val = 0.0;

            if (i == j) {
                // 自干涉计算，引入 Gringarten 理论观察点 0.732
                double x_obs = 0.732 * LfD;
                auto integrand_self = [&](double a) -> double {
                    double dist = std::abs(x_obs - a);
                    if (dist < 1e-15) return 0.0;
                    double arg_dist = gama1 * dist;
                    return safe_bessel_k(0, arg_dist) + Ac_core * safe_bessel_i_scaled(0, arg_dist) * std::exp(arg_dist - 2.0 * arg_g1_rm);
                };

                unsigned max_depth = 15; double tol = 1e-6; double err_est;
                double I_left = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_self, -LfD, x_obs, max_depth, tol, &err_est);
                double I_right = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_self, x_obs, LfD, max_depth, tol, &err_est);
                val = (I_left + I_right) / (z * 2.0 * LfD);
            } else {
                // 互干涉积分
                auto integrand = [&](double a) -> double {
                    double dist_val = std::sqrt(dx * dx + a * a);
                    double arg_dist = gama1 * dist_val;
                    return safe_bessel_k(0, arg_dist) + Ac_core * safe_bessel_i_scaled(0, arg_dist) * std::exp(arg_dist - 2.0 * arg_g1_rm);
                };

                unsigned max_depth = 15; double tol = 1e-6; double err_est;
                double I_full = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, -LfD, LfD, max_depth, tol, &err_est);
                val = I_full / (z * 2.0 * LfD);
            }
            pfD_mat(i, j) = val;
        }
    }

    // ========== 2. 计算井筒等效闭合参数 (V8.2 架构) ==========
    int storageType = ((int)type) % 4; // 0:定井储, 1:线源解, 2:Fair, 3:Hegeman
    double Ceff = CD;
    double S_eff = S;

    if (storageType == 1) {
        // 线源解：退化井筒效应
        Ceff = 0.0;
        S_eff = 0.0;
    } else if (storageType == 2) {
        // Fair 相重分布
        double P_phiD = (std::abs(alpha) > 1e-12) ? (C_phi / (z * (1.0 + alpha * z))) : 0.0;
        Ceff = CD * (1.0 + z * z * P_phiD);
    } else if (storageType == 3) {
        // Hegeman 相重分布
        double P_phiD = 0.0;
        if (std::abs(alpha) > 1e-12) {
            double xx = z * alpha / 2.0;
            double exp_erfc = (xx < 10.0) ? (std::exp(xx * xx) * std::erfc(xx)) : (1.0 / (std::sqrt(M_PI) * xx));
            P_phiD = (C_phi / z) * exp_erfc;
        }
        Ceff = CD * (1.0 + z * z * P_phiD);
    }

    // ========== 3. 组装 2nf+1 维闭合联立系统 ==========
    int sys_size = 2 * n_fracs + 1;
    Eigen::MatrixXd A_sys = Eigen::MatrixXd::Zero(sys_size, sys_size);
    Eigen::VectorXd b_sys = Eigen::VectorXd::Zero(sys_size);

    // 约束1：各缝裂缝面物理响应方程 [z * pfD * q_jD - p_sf,iD = 0]
    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_fracs; ++j) {
            A_sys(i, j) = z * pfD_mat(i, j);
        }
        A_sys(i, n_fracs + i) = -1.0;
    }

    // 约束2：局部表皮流动附加压降 [ -S_eff * q_iD - p_sf,iD + p_wD = 0 ]
    for (int i = 0; i < n_fracs; ++i) {
        int row = n_fracs + i;
        A_sys(row, i) = -S_eff;
        A_sys(row, n_fracs + i) = -1.0;
        A_sys(row, 2 * n_fracs) = 1.0;
    }

    // 约束3：地面定产与井储效应能量分配守恒 [ sum(q_jD) + Ceff * z * p_wD = 1/z ]
    for (int j = 0; j < n_fracs; ++j) {
        A_sys(2 * n_fracs, j) = 1.0;
    }
    A_sys(2 * n_fracs, 2 * n_fracs) = Ceff * z;
    b_sys(2 * n_fracs) = 1.0 / z;

    // ========== 4. 矩阵正则化防崩溃解法 ==========
    // 使用基于矩阵元素极值的自适应微小对角线正则化参数 (吉洪诺夫正则化)
    double max_norm = A_sys.cwiseAbs().colwise().sum().maxCoeff();
    double lam_reg = std::max(1e-12, 1e-10 * max_norm);
    A_sys += lam_reg * Eigen::MatrixXd::Identity(sys_size, sys_size);

    Eigen::VectorXd x_sol = A_sys.fullPivLu().solve(b_sys);

    // 提取全局统一井底流压 p_wD
    return x_sol(2 * n_fracs);
}

void ModelSolver01::precomputeStehfestCoeffs(int N) {
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

double ModelSolver01::getStehfestCoeff(int i, int N) {
    if (m_currentN != N || i < 1 || i > N) return 0.0;
    return m_stehfestCoeffs[i];
}

double ModelSolver01::factorial(int n) {
    if(n <= 1) return 1.0;
    double r = 1.0;
    for(int i = 2; i <= n; ++i) r *= i;
    return r;
}
