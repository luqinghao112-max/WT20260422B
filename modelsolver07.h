/*
 * 文件名: modelsolver07.h
 * 功能描述:
 * 1. 压裂水平井夹层型深化模型 (对应 V8.2 Number 7)。
 * 2. 【终极架构升级】全面支持：非等长裂缝 (Lf_vec)、非均匀位置 (xw_vec)、以及各裂缝独立表皮 (S_vec)。
 * 3. 底层解析器使用 2nf+1 维线性方程组闭合联立，对角线上的表皮压降由独立向量提供。
 */

#ifndef MODELSOLVER07_H
#define MODELSOLVER07_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 返回数据组：时间 t, 压力 P, 压力导数 P'
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver07
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

    explicit ModelSolver07(ModelType type);
    virtual ~ModelSolver07();

    // 设置高精度积分标志
    void setHighPrecision(bool h);

    // 主计算接口：负责量纲转换、参数拆包及多线程拉氏逆变换调度
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& p, const QVector<double>& t = QVector<double>());

    static QString getModelName(ModelType type, bool v = true);
    static QVector<double> generateLogTimeSteps(int c, double s, double e);

private:
    // 时域转换核心：Stehfest N=8 反演与压敏外推防护
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& p,
                             std::function<double(double, const QMap<QString, double>&)> f,
                             QVector<double>& oPD, QVector<double>& oDeriv);

    // 复合拉普拉斯空间参数解析：组装夹层特征函数并提取多维裂缝数组（包含表皮）
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 【重构核心算子】2nf+1 维干扰矩阵构建与求解
    // 引入 S_vec 替换标量 S，实现各缝流动污染独立分配
    double PWD_composite(double z, double fs1, double fs2, double M12,
                         const QVector<double>& LfD_vec, const QVector<double>& xwD,
                         double rmD, double reD, int nF, ModelType type,
                         const QVector<double>& S_vec, double CD, double C_phi, double alpha);

    // 缩放型安全贝塞尔函数群 (规避数值溢出)
    static double safe_bessel_i_sc(int v, double x);
    static double safe_bessel_k_sc(int v, double x);
    static double safe_bessel_k(int v, double x);

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

#endif // MODELSOLVER07_H
