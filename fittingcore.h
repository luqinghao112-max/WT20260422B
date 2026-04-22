/*
 * 文件名: fittingcore.h
 * 文件作用: 试井拟合核心算法（反演引擎）类头文件
 * * 功能与架构描述:
 * 1. 【核心算法】封装基于 Levenberg-Marquardt (LM) 的信赖域非线性优化算法。
 * 2. 【数据抽样】处理庞大实测数据点集的对数等距抽样 (getLogSampledData)，加速雅可比矩阵运算。
 * 3. 【多表皮追踪】全面兼容 Model 07~09 引入的独立表皮向量 (S_1, S_2...S_nf)。
 * 4. 【鲁棒性增强】支持多起点寻优 (useMultiStart) 避免局部极小，内置等长降维通道与 L2 正则化约束机制。
 */

#ifndef FITTINGCORE_H
#define FITTINGCORE_H

#include <QObject>
#include <QVector>
#include <QMap>
#include <QFutureWatcher>
#include "modelmanager.h"
#include "fittingsamplingdialog.h"
#include "fittingparameterchart.h"

class FittingCore : public QObject
{
    Q_OBJECT
public:
    explicit FittingCore(QObject *parent = nullptr);

    // ===================== 基础设置接口 =====================
    // 挂载全局模型调度中枢
    void setModelManager(ModelManager* m);

    // 装载现场观测数据（时间 t, 压差 p, 压力导数 d）
    void setObservedData(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);

    // 配置抽样策略（控制计算密度，提升 LM 迭代速度）
    void setSamplingSettings(const QList<SamplingInterval>& intervals, bool enabled);

    // ===================== 拟合控制接口 =====================
    /**
     * @brief 启动自动拟合计算 (异步执行)
     * @param modelType     当前选定的模型类型 (决定底层调用的物理求解器)
     * @param params        待拟合参数集合 (含是否拟合标记、上下限等，包含动态 S_i 阵列)
     * @param weight        目标函数权重 (控制压力差和导数在残差中的比例)
     * @param useLimits     是否在迭代更新步中强制截断超出用户设定边界的参数
     * @param useMultiStart 是否开启广域多起点探索（跳出局部最优）
     */
    void startFit(ModelManager::ModelType modelType, const QList<FitParameter>& params, double weight, bool useLimits = false, bool useMultiStart = false);

    // 中断正在进行的迭代寻优
    void stopFit();

    // ===================== 数据处理与数学算子 =====================
    // 执行数据抽样过滤
    void getLogSampledData(const QVector<double>& srcT, const QVector<double>& srcP, const QVector<double>& srcD,
                           QVector<double>& outT, QVector<double>& outP, QVector<double>& outD);

    /**
     * @brief 计算目标残差向量
     * 包含对观测点计算的偏差，并在多自由度时内置 L2 正则化伪残差，防止参数发散
     */
    QVector<double> calculateResiduals(const QMap<QString, double>& params, ModelManager::ModelType modelType, double weight,
                                       const QVector<double>& t, const QVector<double>& obsP, const QVector<double>& obsD);

    // 计算误差平方和 (SSE)
    double calculateSumSquaredError(const QVector<double>& residuals);

    /**
     * @brief 静态工具：参数进入模型前的预清洗与防崩盘降维
     * @details 拦截并处理 "use_equal_lf" 强制等长约束，以及线源解环境下的 S_i 置零处理
     */
    static QMap<QString, double> preprocessParams(const QMap<QString, double>& rawParams, ModelManager::ModelType type);

    /**
     * @brief 基于观测数据的几何特征 (导数平台、双对数峰谷) 自动推算收敛初值
     */
    static QMap<QString, double> estimateInitialParams(
        const QVector<double>& obsTime,
        const QVector<double>& obsDeltaP,
        const QVector<double>& obsDerivative,
        ModelManager::ModelType modelType);

signals:
    // 每完成一步有效迭代，向 UI 广播一次最新参数和理论曲线
    void sigIterationUpdated(double error, QMap<QString,double> params, QVector<double> t, QVector<double> p, QVector<double> d);

    // 更新进度条
    void sigProgress(int percent);

    // 异步拟合线程结束通知
    void sigFitFinished();

private:
    ModelManager* m_modelManager;
    QVector<double> m_obsTime;
    QVector<double> m_obsDeltaP;
    QVector<double> m_obsDerivative;

    bool m_isCustomSamplingEnabled;
    QList<SamplingInterval> m_customIntervals;

    bool m_stopRequested;
    QFutureWatcher<void> m_watcher;

    // ===================== 内部核心算法实现 =====================
    // 异步任务的包装入口
    void runOptimizationTask(ModelManager::ModelType modelType, QList<FitParameter> fitParams, double weight, bool useLimits, bool useMultiStart);

    // Levenberg-Marquardt (LM) 信赖域算法主体逻辑
    void runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params, double weight, bool useLimits, bool useMultiStart);

    // 中心差分法计算雅可比矩阵 (J)
    QVector<QVector<double>> computeJacobian(const QMap<QString, double>& params, const QVector<double>& baseResiduals,
                                             const QVector<int>& fitIndices, ModelManager::ModelType modelType,
                                             const QList<FitParameter>& currentFitParams, double weight,
                                             const QVector<double>& t, const QVector<double>& obsP, const QVector<double>& obsD);

    // 使用 Cholesky (LDLT) 分解法求解法方程组 (H + lambda*I)*delta = g
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
};

#endif // FITTINGCORE_H
