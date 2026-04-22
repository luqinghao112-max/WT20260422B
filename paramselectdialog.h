/*
 * 文件名: paramselectdialog.h
 * 文件作用: 拟合参数选择与边界配置对话框的头文件。
 * * 功能与架构描述:
 * 1. 承担“自动拟合”前最关键的参数约束配置任务。用户在此决定哪些参数参与反演求解（isFit），
 * 以及这些参数的物理安全边界（下限 min、上限 max）和微调步长（step）。
 * 2. 【核心机制】内置了针对“裂缝条数 (nf)”的动态监听机制。当用户在表格中修改 nf 时，
 * 系统会动态增删对应的阵列参数（如 Lf_i 长度阵列、xf_i 位置阵列、S_i 表皮阵列）。
 * 3. 提供了与底层图版引擎及拟合引擎无缝对接的参数结构体 (FitParameter) 传递接口。
 */

#ifndef PARAMSELECTDIALOG_H
#define PARAMSELECTDIALOG_H

#include <QDialog>
#include <QTableWidget>
#include <QEvent>
#include "fittingparameterchart.h"

namespace Ui {
class ParamSelectDialog;
}

class ParamSelectDialog : public QDialog
{
    Q_OBJECT

public:
    /**
     * @brief 构造函数：初始化参数对话框及表格结构
     * @param params    由上一级传入的初始拟合参数结构体列表
     * @param modelType 当前激活的物理模型枚举，用于生成对应的默认参数和判定模型权限
     * @param fitTime   全局拟合所允许的最大时间截断阈值
     * @param useLimits 标识当前是否强制启用物理上下限约束（用于初始化复选框状态）
     * @param parent    父对象指针
     */
    explicit ParamSelectDialog(const QList<FitParameter>& params,
                               ModelManager::ModelType modelType,
                               double fitTime,
                               bool useLimits,
                               QWidget *parent = nullptr);

    /**
     * @brief 析构函数：负责 UI 树的安全释放
     */
    ~ParamSelectDialog();

    /**
     * @brief 提取接口：供外部 (FittingWidget) 获取用户在表格中修改后的最终参数配置
     * @return 包含最新 value、min、max、isFit 等状态的参数列表
     */
    QList<FitParameter> getUpdatedParams() const;

    /**
     * @brief 提取接口：获取用户设定的最大拟合时间
     */
    double getFittingTime() const;

    /**
     * @brief 提取接口：获取用户是否勾选了“使用上下限控制范围”
     */
    bool getUseLimits() const;

signals:
    /**
     * @brief 信号：当用户点击“智能初值”按钮时，通知外部主引擎根据观测数据特征进行启发式估算
     */
    void estimateInitialParamsRequested();

protected:
    /**
     * @brief 事件过滤器重写：
     * 核心作用是拦截 QDoubleSpinBox (数值调节框) 的鼠标滚轮事件。
     * 防止用户在滚动查看长表格时，鼠标指针不小心划过数值框导致参数被意外修改。
     */
    bool eventFilter(QObject *obj, QEvent *event) override;

    // 授权拟合界面直接访问部分私有方法以实现深度联动
    friend class FittingWidget;
    friend class WT_MultidataFittingWidget;

private:
    Ui::ParamSelectDialog *ui;

public:
    /**
     * @brief 数据采集器：遍历 UI 表格中的所有控件，将文本框、复选框的最新状态写回 m_params 内存副本
     */
    void collectData();

private:
    ModelManager::ModelType m_modelType; // 当前模型的 ID 缓存
    QList<FitParameter> m_params;        // 当前对话框维护的参数列表内存副本

    /**
     * @brief 核心渲染引擎：负责清空旧表格并根据 m_params 重新生成列宽、表头、及各种内嵌控件
     */
    void initTable();

    /**
     * @brief UI 美化函数：根据该行参数是否被勾选为“拟合变量”，动态变换整行的背景色（首尾列除外）
     * @param row   目标行号
     * @param isFit 是否参与拟合
     */
    void updateRowAppearance(int row, bool isFit);

private slots:
    // --- 底部控制栏标准槽函数 ---
    void onConfirm();                 // 用户点击“确定”：保存数据并关闭弹窗
    void onCancel();                  // 用户点击“取消”：丢弃更改并关闭弹窗
    void onResetParams();             // 恢复出厂默认参数，并覆盖当前表格
    void onAutoLimits();              // 根据当前数值自动扩缩上下限（例如 min = val/10, max = val*10）
    void onEstimateInitialParams();   // 触发启发式初值估算请求

    /**
     * @brief 【核心动态引擎】：专门用于监听“裂缝条数 (nf)”修改的异步槽函数。
     * 当 nf 发生改变时，动态增加或削减 Lf_i (半长)、xf_i (位置)、S_i (表皮) 的参数行，
     * 避免多段模型参数面板的冗余或缺失。
     */
    Q_INVOKABLE void onNfChanged();
};

#endif // PARAMSELECTDIALOG_H
