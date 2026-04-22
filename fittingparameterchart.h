/*
 * 文件名: fittingparameterchart.h
 * 文件作用:
 * 1. 声明拟合参数表格的管理器类，负责 UI 表格与底层参数数据(FitParameter)的绑定。
 * 2. 明确定义了 FitParameter 结构体，用于在 UI 和拟合引擎之间传递参数的名称、显示名、当前值、上下限等信息。
 * 3. 声明了动态重构参数阵列的方法 syncFractureParams，用于在裂缝条数改变时动态调整参数列表。
 */

#ifndef FITTINGPARAMETERCHART_H
#define FITTINGPARAMETERCHART_H

#include <QObject>
#include <QTableWidget>
#include <QList>
#include <QMap>
#include <QString>
#include <QTimer>
#include "modelmanager.h" // 仅包含中枢管理器以使用 ModelType 枚举

// =========================================================================
// 全局物理参数结构体定义
// 在此处明确定义，提供给上层 UI 和下层 FittingCore 使用，杜绝未声明报错
// =========================================================================
struct FitParameter {
    QString name;           // 参数内部键名 (如 "kf", "Lf_1", "xf_2", "S_3")
    QString displayName;    // 参数中文显示名称 (如 "裂缝渗透率", "第1段裂缝半长")
    double value;           // 当前设定数值
    double min;             // LM非线性拟合约束下限
    double max;             // LM非线性拟合约束上限
    double step;            // 界面滚轮调节步长
    bool isFit;             // 是否参与拟合运算
    bool isVisible;         // 是否在界面上显示
};

class FittingParameterChart : public QObject
{
    Q_OBJECT
public:
    /**
     * @brief 构造函数，接管 UI 上的 TableWidget
     * @param parentTable 关联的 QTableWidget 指针
     * @param parent 父对象
     */
    explicit FittingParameterChart(QTableWidget *parentTable, QObject *parent = nullptr);

    /**
     * @brief 设置模型管理器中枢，用于获取物理量默认值
     * @param m ModelManager 指针
     */
    void setModelManager(ModelManager *m);

    /**
     * @brief 切换当前试井模型，保留同名参数的值，并调整新特有参数
     * @param newType 新的模型枚举类型
     */
    void switchModel(ModelManager::ModelType newType);

    /**
     * @brief 获取当前表格配置的参数列表 (平铺的标量列表，适配 LM 引擎)
     * @return 包含所有表格内参数配置信息的 FitParameter 列表
     */
    QList<FitParameter> getParameters() const;

    /**
     * @brief 从外部注入并覆盖当前的参数列表，刷新表格
     * @param p 外部传入的参数列表
     */
    void setParameters(const QList<FitParameter> &p);

    /**
     * @brief 将当前模型的参数重置为默认值
     * @param type 目标模型类型
     * @param preserveStates 是否保留当前的“参与拟合”、“可见性”勾选状态
     */
    void resetParams(ModelManager::ModelType type, bool preserveStates = true);

    /**
     * @brief 提取表格中的原始文本(字符串字典)，用于外部映射校验
     * @return 包含键名为参数名，值为参数显示文本的映射表
     */
    QMap<QString, QString> getRawParamTexts() const;

    /**
     * @brief 将 UI 表格用户的修改，同步提取并保存至内部的 m_params 列表
     */
    void updateParamsFromTable();

    /**
     * @brief 静态工具：基于物理意义自动评估并收紧各参数的拟合上下限及滚轮步长
     * @param params 需要被调整上下限的参数列表
     */
    static void adjustLimits(QList<FitParameter>& params);

    /**
     * @brief 触发自动限幅
     */
    void autoAdjustLimits();

    /**
     * @brief 静态工具：生成指定模型的默认参数模板
     * @param type 目标模型类型
     * @param overrideNf (可选) 强制指定的裂缝条数，用于动态展开多缝长/多位置/多表皮参数
     * @return 生成的包含默认值的参数列表
     */
    static QList<FitParameter> generateDefaultParams(ModelManager::ModelType type, int overrideNf = -1);

    /**
     * @brief 静态工具：获取参数的中文名称与单位展示信息
     * @param n 底层参数键名
     * @param cN 输出: 中文显示名称
     * @param s 输出: 缩写 (原参数名)
     * @param uS 输出: 单位字符串的内部表示 (同原参数名)
     * @param u 输出: 实际显示的单位字符串
     */
    static void getParamDisplayInfo(const QString &n, QString &cN, QString &s, QString &uS, QString &u);

signals:
    /**
     * @brief 鼠标滚轮修改数值时触发的防抖信号，通知界面刷新图版
     */
    void parameterChangedByWheel();

protected:
    /**
     * @brief 事件过滤器，用于拦截和处理鼠标滚轮修改数值，防止误操作
     */
    bool eventFilter(QObject *w, QEvent *e) override;

private slots:
    /**
     * @brief 滚轮防抖定时器超时响应
     */
    void onWheelDebounceTimeout();

    /**
     * @brief 表格项被用户编辑后触发的槽函数
     * @param item 被编辑的表格项
     */
    void onTableItemChanged(QTableWidgetItem *item);

private:
    /**
     * @brief 清空并根据当前的 m_params 重新填充表格内容
     */
    void refreshParamTable();

    /**
     * @brief 将单个参数对象转换为一行数据，并插入到表格中
     * @param p 待插入的参数对象
     * @param serialNo 自动递增的序号
     * @param highlight 是否高亮显示该行（针对参与拟合的参数）
     */
    void addRowToTable(const FitParameter& p, int& serialNo, bool highlight);

    /**
     * @brief 核心新增机制：当裂缝条数 nf 发生改变时，动态伸缩 Lf_i, xf_i, S_i 等相关多段阵列参数的行数
     * @param newNf 新的裂缝条数
     */
    void syncFractureParams(int newNf);

private:
    QTableWidget *m_table;               // 关联的 UI 表格组件
    ModelManager *m_modelManager;        // 中枢管理器指针
    QList<FitParameter> m_params;        // 核心参数列表缓存
    QTimer *m_wheelTimer;                // 滚轮防抖定时器
    ModelManager::ModelType m_modelType; // 当前绑定的模型类型
};

#endif // FITTINGPARAMETERCHART_H
