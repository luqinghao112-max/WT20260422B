/*
 * 文件名: modelsolver04.h
 * 功能描述:
 * 1. 压裂水平井夹层型深化模型 (Model 4对应模型，涵盖双重孔隙等效等)。
 * 2. 【核心升级】支持多段不等长裂缝 (Lf_vec) 与 非均匀间距布置 (xw_vec)。
 * 3. 井筒与裂缝系统采用 2nf+1 维闭合联立矩阵求解，统一单个表皮系数 S。
 */

#ifndef MODELSOLVER04_H
#define MODELSOLVER04_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver04
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

    explicit ModelSolver04(ModelType type);
    virtual ~ModelSolver04();

    void setHighPrecision(bool high);

    // 计算理论曲线主入口
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime = QVector<double>());

    static QString getModelName(ModelType type, bool verbose = true);
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

private:
    // 时域转换数值核心 (Stehfest反演 + 压敏切线外推)
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    // 拉普拉斯空间分配与组装
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 【核心重构】底层解析器，支持传入独立长度 LfD_vec 和位置 xwD 向量
    double PWD_composite(double z, double fs1, double fs2, double M12,
                         const QVector<double>& LfD_vec, const QVector<double>& xwD,
                         double rmD, double reD, int n_fracs, ModelType type,
                         double S, double CD, double C_phi, double alpha);

    // Boost标准库安全贝塞尔函数
    static double safe_bessel_i_scaled(int v, double x);
    static double safe_bessel_k_scaled(int v, double x);
    static double safe_bessel_k(int v, double x);

    double getStehfestCoeff(int i, int N);
    void precomputeStehfestCoeffs(int N);
    double factorial(int n);

private:
    ModelType m_type;
    bool m_highPrecision;
    QVector<double> m_stehfestCoeffs;
    int m_currentN;
};

#endif // MODELSOLVER04_H
