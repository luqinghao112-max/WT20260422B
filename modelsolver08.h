/*
 * 文件名: modelsolver08.h
 * 功能描述:
 * 1. 压裂水平井页岩型深化试井模型 (对应 V8.2 Number 8)。
 * 2. 【核心架构】全面支持页岩瞬态窜流特征，并引入独立表皮向量 S_vec。
 * 3. 采用 2nf+1 维线性方程组同步联立，每一条缝具有独立的物理属性（长度、位置、表皮污染）。
 */

#ifndef MODELSOLVER08_H
#define MODELSOLVER08_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 定义返回值元组：<真实时间序列, 真实压力序列, 压力导数序列>
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver08
{
public:
    // 枚举 36 种子模型 (根据内外边界和井筒存储条件组合)
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

    explicit ModelSolver08(ModelType type);
    virtual ~ModelSolver08();

    // 设置是否启用高精度数值防溢出模式
    void setHighPrecision(bool high);

    // 理论曲线计算主入口：处理输入参数，执行量纲转换与逆变换调度
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime = QVector<double>());

    // 获取当前模型名称的文本描述
    static QString getModelName(ModelType type, bool verbose = true);

    // 生成对数等比分布的时间步向量
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

private:
    // 时域转换核心：Stehfest N=8 数值反演，内置压敏非线性摄动及泰勒外推防护
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    // 拉普拉斯空间分配中心：组装页岩特征函数，提取独立的 LfD_vec, xwD 以及 S_vec
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 【重构核心算子】2nf+1 维闭合干扰大矩阵求解器，接收多段独立表皮
    double PWD_composite(double z, double fs1, double fs2, double M12,
                         const QVector<double>& LfD_vec, const QVector<double>& xwD,
                         double rmD, double reD, int n_fracs, ModelType type,
                         const QVector<double>& S_vec, double CD, double C_phi, double alpha);

    // 外部双孔伪稳态特征函数计算
    double calc_fs_dual(double u, double omega, double lambda);

    // 页岩型基质瞬态扩散特征函数 (包含 Tanh 非线性核)
    double calc_fs_shale(double u, double omega, double lambda);

    // 带有缩放保护机制的底层安全贝塞尔函数 (针对极端物理量纲)
    static double safe_bessel_i_scaled(int v, double x);
    static double safe_bessel_k_scaled(int v, double x);
    static double safe_bessel_k(int v, double x);

    // Stehfest 逆变换系数计算模块
    double getStehfestCoeff(int i, int N);
    void precomputeStehfestCoeffs(int N);
    double factorial(int n);

private:
    ModelType m_type;
    bool m_highPrecision;
    QVector<double> m_stehfestCoeffs;
    int m_currentN;
};

#endif // MODELSOLVER08_H
