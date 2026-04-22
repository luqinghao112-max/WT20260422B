/*
 * 文件名: modelsolver09.h
 * 功能描述:
 * 1. 压裂水平井混积型（三重孔隙）深化试井模型 (对应 V8.2 Number 9)。
 * 2. 【终极版架构】支持：非等长裂缝分布 (Lf_vec)、非均匀间距 (xw_vec) 以及 各裂缝独立表皮系数 (S_vec)。
 * 3. 采用 2nf+1 维线性闭合方程组，实现不同裂缝污染程度的精确模拟与缝间流量动态耦合。
 */

#ifndef MODELSOLVER09_H
#define MODELSOLVER09_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 定义返回值元组：<真实时间序列, 无因次/有因次压力, 压力导数>
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver09
{
public:
    enum ModelType {
        Model_1 = 0, Model_2, Model_3, Model_4,
        Model_5, Model_6, Model_7, Model_8,
        Model_9, Model_10, Model_11, Model_12,
        Model_13, Model_14, Model_15, Model_16,
        Model_17, Model_18, Model_19, Model_20,
        Model_21, Model_22, Model_23, Model_24,
        Model_25, Model_26, Model_27, Model_28,
        Model_29, Model_30, Model_31, Model_32,
        Model_33, Model_34, Model_35, Model_36
    };

    explicit ModelSolver09(ModelType type);
    virtual ~ModelSolver09();

    // 设置高精度积分标定
    void setHighPrecision(bool h);

    // 理论曲线计算主控接口：提取物性、降维时间、调度多线程反演
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& p, const QVector<double>& t = QVector<double>());

    static QString getModelName(ModelType type, bool v = true);
    static QVector<double> generateLogTimeSteps(int c, double s, double e);

private:
    // 时域转换核心：Stehfest N=8 数值反演，内置压敏非线性扰动防护
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& p,
                             std::function<double(double, const QMap<QString, double>&)> f,
                             QVector<double>& oPD, QVector<double>& oDeriv);

    // 复合拉普拉斯空间参数解析：组装三重孔隙特征函数并提取多维独立裂缝向量
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 【重构终极核心算子】2nf+1 维闭合干扰大矩阵求解器，接收各段独立表皮向量
    double PWD_composite(double z, double fs1, double fs2, double M12,
                         const QVector<double>& LfD_vec, const QVector<double>& xwD,
                         double rmD, double reD, int nF, ModelType type,
                         const QVector<double>& S_vec, double CD, double C_phi, double alpha);

    // 基础数学安全贝塞尔函数群 (缩放保护)
    static double safe_bessel_i_sc(int v, double x);
    static double safe_bessel_k_sc(int v, double x);
    static double safe_bessel_k(int v, double x);

    // 页岩型外区 Tanh 特征函数计算
    double calc_fs_shale(double u, double o, double l);

    // 数值反演系数发生器
    double getStehfestCoeff(int i, int N);
    void precomputeStehfestCoeffs(int N);
    double factorial(int n);

private:
    ModelType m_type;
    bool m_highPrecision;
    QVector<double> m_stehfestCoeffs;
    int m_currentN;
};

#endif // MODELSOLVER09_H
