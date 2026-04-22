/*
 * 文件名: modelsolver05.h
 * 功能描述:
 * 1. 压裂水平井页岩型复合深化模型 (Model 5对应模型，涵盖页岩瞬态窜流等效)。
 * 2. 【核心升级】全面支持多段不等长裂缝 (Lf_vec) 与 非均匀间距布置 (xw_vec)。
 * 3. 井筒与裂缝系统采用 2nf+1 维闭合联立矩阵求解，统一单个表皮系数 S。
 * 4. 包含 Boost 库支持的高精度数值积分及 Tanh 特征函数解析。
 */

#ifndef MODELSOLVER05_H
#define MODELSOLVER05_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 定义返回值类型：包含时间、压力、压力导数三个向量
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver05
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

    explicit ModelSolver05(ModelType type);
    virtual ~ModelSolver05();

    // 设置是否启用高精度计算
    void setHighPrecision(bool high);

    // 计算理论曲线主入口：解析参数，时间无因次化，调用拉氏逆变换，并还原为真实物理量
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime = QVector<double>());

    // 获取当前模型名称的描述字符串
    static QString getModelName(ModelType type, bool verbose = true);

    // 生成对数分布的时间步向量
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

private:
    // 时域转换数值核心 (Stehfest反演 + 压敏切线外推)
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    // 拉普拉斯空间分配与组装入口：计算页岩特征函数 fs1/fs2 并提取裂缝数组
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 【核心重构】底层解析器，支持传入独立长度 LfD_vec 和位置 xwD 向量，联立求解井底压力
    double PWD_composite(double z, double fs1, double fs2, double M12,
                         const QVector<double>& LfD_vec, const QVector<double>& xwD,
                         double rmD, double reD, int n_fracs, ModelType type,
                         double S, double CD, double C_phi, double alpha);

    // 计算外区双孔介质的特征函数 (伪稳态)
    double calc_fs_dual(double u, double omega, double lambda);

    // 计算页岩特有的瞬态窜流特征函数 (Tanh扩散核)
    double calc_fs_shale(double u, double omega, double lambda);

    // Boost标准库安全贝塞尔函数，加入 scaled 机制防止大自变量溢出
    static double safe_bessel_i_scaled(int v, double x);
    static double safe_bessel_k_scaled(int v, double x);
    static double safe_bessel_k(int v, double x);

    // Stehfest 逆变换系数计算相关
    double getStehfestCoeff(int i, int N);
    void precomputeStehfestCoeffs(int N);
    double factorial(int n);

private:
    ModelType m_type;
    bool m_highPrecision;
    QVector<double> m_stehfestCoeffs;
    int m_currentN;
};

#endif // MODELSOLVER05_H
