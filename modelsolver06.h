/*
 * 文件名: modelsolver06.h
 * 功能描述:
 * 1. 混积型（三重孔隙）试井解释非均质深化模型 (对应 V8.2 Number 6)。
 * 2. 【核心架构】引入了支持非等长 (Lf_vec)、非均匀间距 (xw_vec) 的裂缝向量传递机制。
 * 3. 使用 2nf+1 维线性方程组闭合联立表皮压降、井筒储集、以及各缝地层干涉流量。
 * 4. 全面采用 Boost 标准库处理复杂拉氏空间的贝塞尔运算与高斯-克罗内德数值积分。
 */

#ifndef MODELSOLVER06_H
#define MODELSOLVER06_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 返回格式：<时间序列, 无因次/有因次压力, 压力导数>
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver06
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

    explicit ModelSolver06(ModelType type);
    virtual ~ModelSolver06();

    // 设置高精度模式标定
    void setHighPrecision(bool h);

    // 主计算接口：解析物性参数，生成时间步，并执行拉氏空间到真实时域的转换
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& p, const QVector<double>& t = QVector<double>());

    static QString getModelName(ModelType type, bool v = true);
    static QVector<double> generateLogTimeSteps(int c, double s, double e);

private:
    // 时域转换核心：Stehfest N=8 数值反演，内置 gamaD 压敏摄动泰勒外推
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& p,
                             std::function<double(double, const QMap<QString, double>&)> f,
                             QVector<double>& oPD, QVector<double>& oDeriv);

    // 复合拉普拉斯空间特征函数装配 (处理三重孔隙内区与多类型外区)
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 【核心算子】2nf+1 维闭合干扰矩阵求解器，动态接收各缝的长度与位置向量
    double PWD_composite(double z, double fs1, double fs2, double M12,
                         const QVector<double>& LfD_vec, const QVector<double>& xwD,
                         double rmD, double reD, int nF, ModelType type,
                         double S, double CD, double C_phi, double alpha);

    // 基础数学安全防御函数
    static double safe_bessel_i_sc(int v, double x);
    static double safe_bessel_k_sc(int v, double x);
    static double safe_bessel_k(int v, double x);

    // 页岩型外区 Tanh 特征函数
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

#endif // MODELSOLVER06_H
