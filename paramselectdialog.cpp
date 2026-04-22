/*
 * 文件名: paramselectdialog.cpp
 * 文件作用: 参数选择配置对话框的具体实现。
 * * 核心重构与优化记录:
 * 1. 【事件过滤器优化】：通过 eventFilter 彻底锁死了 QDoubleSpinBox 的滚轮事件，消除误触隐患。
 * 2. 【多段独立表皮支持】：在 onNfChanged() 方法中，补齐了 S_1 ~ S_nf 的动态解析、注入与销毁逻辑。
 * 3. 【异步安全重绘】：巧妙利用 Qt::QueuedConnection 配合 QMetaObject::invokeMethod 处理 nf 修改事件，
 * 防止在编辑框尚未失去焦点时强行销毁控件导致底层 C++ 野指针崩溃。
 */

#include "paramselectdialog.h"
#include "ui_paramselectdialog.h"
#include "fittingparameterchart.h"
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QDebug>
#include <QMetaObject>
#include "standard_messagebox.h"

/**
 * @brief 内部增强型组件：SmartDoubleSpinBox
 * 继承自 QDoubleSpinBox。原生组件在显示浮点数时总会带着长长的一串 0 (如 50.000000)。
 * 重写 textFromValue 方法，利用 'g' 格式化字符智能抹除多余的 0，使表格清爽整洁。
 */
class SmartDoubleSpinBox : public QDoubleSpinBox {
public:
    explicit SmartDoubleSpinBox(QWidget* parent = nullptr) : QDoubleSpinBox(parent) {}
    QString textFromValue(double value) const override {
        return QString::number(value, 'g', decimals() == 0 ? 10 : decimals());
    }
};

// 构造函数实现
ParamSelectDialog::ParamSelectDialog(const QList<FitParameter> &params, ModelManager::ModelType modelType, double fitTime, bool useLimits, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ParamSelectDialog),
    m_modelType(modelType),
    m_params(params) // 严格按照声明顺序初始化，消除 -Wreorder 编译警告
{
    ui->setupUi(this);
    this->setWindowTitle("拟合参数配置");

    // 初始化全局控制区的数值与状态
    ui->spinTimeMax->setValue(fitTime);
    ui->chkUseLimits->setChecked(useLimits);

    // 信号槽挂载
    connect(ui->btnOk, &QPushButton::clicked, this, &ParamSelectDialog::onConfirm);
    connect(ui->btnCancel, &QPushButton::clicked, this, &ParamSelectDialog::onCancel);
    connect(ui->btnResetDefaults, &QPushButton::clicked, this, &ParamSelectDialog::onResetParams);
    connect(ui->btnAutoLimits, &QPushButton::clicked, this, &ParamSelectDialog::onAutoLimits);

    // 动态嗅探并注入“智能初值”按钮 (防止旧版本 ui 文件缺失该控件)
    QPushButton* estimateButton = this->findChild<QPushButton*>(QStringLiteral("btnEstimateInitial"));
    if (!estimateButton && ui->horizontalLayout_Tools) {
        estimateButton = new QPushButton(QStringLiteral("智能初值"), this);
        estimateButton->setObjectName(QStringLiteral("btnEstimateInitial"));
        ui->horizontalLayout_Tools->insertWidget(4, estimateButton);
    }

    if (estimateButton) {
        connect(estimateButton, &QPushButton::clicked, this, &ParamSelectDialog::onEstimateInitialParams);
        estimateButton->setVisible(true);
        estimateButton->setEnabled(true);
    }

    ui->btnCancel->setAutoDefault(false);

    // 初始化并渲染完整的参数表格
    initTable();
}

ParamSelectDialog::~ParamSelectDialog()
{
    delete ui;
}

// 滚轮事件吞噬器
bool ParamSelectDialog::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::Wheel) {
        // 若当前焦点对象是数值调节框，拦截该事件，不向下传递
        if (qobject_cast<QAbstractSpinBox*>(obj)) {
            return true;
        }
    }
    return QDialog::eventFilter(obj, event);
}

// 恢复默认参数逻辑
void ParamSelectDialog::onResetParams()
{
    if (!Standard_MessageBox::question(this, "确认", "确定要重置为该模型的默认参数吗？")) {
        return;
    }
    // 请求工具类生成全新的默认参数组，并根据量级自动配置上下限，最后触发界面重绘
    m_params = FittingParameterChart::generateDefaultParams(m_modelType);
    FittingParameterChart::adjustLimits(m_params);
    initTable();
}

// 自动拓宽上下限逻辑
void ParamSelectDialog::onAutoLimits()
{
    collectData(); // 先将用户填写的最新 value 同步到内存
    FittingParameterChart::adjustLimits(m_params); // 重新计算 min/max
    initTable();   // 刷新 UI
    Standard_MessageBox::info(this, "提示", "参数上下限及滚轮步长已根据当前值更新。");
}

void ParamSelectDialog::onEstimateInitialParams()
{
    emit estimateInitialParamsRequested();
}

// UI 行高亮染色逻辑
void ParamSelectDialog::updateRowAppearance(int row, bool isFit)
{
    QString yellowHex = "#FFFFCC"; // 柔和的浅黄色提醒
    for(int col = 0; col < ui->tableWidget->columnCount(); ++col) {
        QTableWidgetItem* item = ui->tableWidget->item(row, col);
        if(item) {
            // 首列(拟合勾选)和尾列(显示勾选)不染色，保持控件原本色调
            bool shouldHighlight = isFit && (col != 0 && col != 7);
            if(shouldHighlight) {
                item->setBackground(QColor(yellowHex));
            } else {
                item->setData(Qt::BackgroundRole, QVariant()); // 清除背景色
            }
        }
    }
}

// =========================================================================
// 核心模块：表格渲染器 (initTable)
// =========================================================================
void ParamSelectDialog::initTable()
{
    ui->tableWidget->clear();
    QStringList headers;
    // 制定新的列头架构，将“是否拟合”前置，“是否显示”后置，符合人类视线流
    headers << "拟合变量" << "当前数值" << "单位" << "参数名称" << "下限" << "上限" << "滚轮步长" << "显示";

    ui->tableWidget->setColumnCount(headers.size());
    ui->tableWidget->setHorizontalHeaderLabels(headers);
    ui->tableWidget->setRowCount(m_params.size());
    ui->tableWidget->setAlternatingRowColors(true);

    // CSS 美化表头和单元格网格
    ui->tableWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { background-color: #f0f0f0; border: 1px solid #dcdcdc; padding: 4px; font-weight: bold; }");
    ui->tableWidget->verticalHeader()->setVisible(false);
    ui->tableWidget->verticalHeader()->setDefaultSectionSize(40);
    ui->tableWidget->setStyleSheet("QTableWidget { border: 1px solid #dcdcdc; gridline-color: #e0e0e0; selection-background-color: transparent; } QTableWidget::item { padding: 4px; }");

    QString checkBoxStyle = "QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; } QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; } QCheckBox::indicator:hover { border-color: #0078d7; }";
    QString transparentSpinStyle = "QDoubleSpinBox { background-color: transparent; border: none; margin: 0px; }";
    QString transparentWidgetStyle = "QWidget { background-color: transparent; }";

    // 逐行解析 m_params 并映射至 UI
    for(int i = 0; i < m_params.size(); ++i) {
        FitParameter p = m_params[i];

        bool isEta12 = (p.name == "eta12");
        bool isNf = (p.name == "nf");

        for (int col = 0; col < 8; ++col) {
            if (!ui->tableWidget->item(i, col)) {
                ui->tableWidget->setItem(i, col, new QTableWidgetItem());
            }
        }

        // --- 第 0 列: 拟合使能开关 ---
        QWidget* pWidgetFit = new QWidget();
        pWidgetFit->setStyleSheet(transparentWidgetStyle);
        QHBoxLayout* pLayoutFit = new QHBoxLayout(pWidgetFit);
        QCheckBox* chkFit = new QCheckBox();
        chkFit->setChecked(isEta12 ? false : p.isFit);
        chkFit->setStyleSheet(checkBoxStyle);
        pLayoutFit->addWidget(chkFit);
        pLayoutFit->setAlignment(Qt::AlignCenter);
        pLayoutFit->setContentsMargins(0,0,0,0);

        // 无因次裂缝长度(LfD)及衍生量(eta12)禁止作为独立拟合变量
        if (p.name == "LfD" || isEta12) {
            chkFit->setEnabled(false);
            chkFit->setChecked(false);
        }
        ui->tableWidget->setCellWidget(i, 0, pWidgetFit);

        // --- 第 1 列: 当前数值 (Value) ---
        SmartDoubleSpinBox* spinVal = new SmartDoubleSpinBox();
        spinVal->setStyleSheet(transparentSpinStyle);

        if (isNf) {
            spinVal->setRange(1.0, 9e9);
            spinVal->setDecimals(0); // 裂缝条数必须为整数

            // 【重要：事件绑定】监听 nf 数值的修改动作，采用 QueuedConnection 以策安全
            connect(spinVal, &QDoubleSpinBox::editingFinished, this, [this]() {
                QMetaObject::invokeMethod(this, "onNfChanged", Qt::QueuedConnection);
            });
        } else {
            spinVal->setRange(-9e9, 9e9);
            spinVal->setDecimals(10);
        }
        spinVal->setValue(p.value);

        if(isEta12) spinVal->setEnabled(false);
        else spinVal->installEventFilter(this); // 安装滚轮吞噬器

        ui->tableWidget->setCellWidget(i, 1, spinVal);

        // --- 第 2 列: 物理单位 ---
        QString dummy, dummy2, dummy3, unitStr;
        FittingParameterChart::getParamDisplayInfo(p.name, dummy, dummy2, dummy3, unitStr);
        if (unitStr == "无因次" || unitStr == "小数") unitStr = "-";

        QTableWidgetItem* unitItem = ui->tableWidget->item(i, 2);
        unitItem->setText(unitStr);
        unitItem->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
        unitItem->setFlags(unitItem->flags() & ~Qt::ItemIsEditable);

        // --- 第 3 列: 参数名称与备注 ---
        QString displayNameFull = QString("%1 (%2)").arg(p.displayName).arg(p.name);
        if (isEta12) displayNameFull += " [只读衍生]";

        QTableWidgetItem* nameItem = ui->tableWidget->item(i, 3);
        nameItem->setText(displayNameFull);
        nameItem->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
        nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable);
        nameItem->setData(Qt::UserRole, p.name);

        // --- 第 4 列: 寻优下限 (Min) ---
        SmartDoubleSpinBox* spinMin = new SmartDoubleSpinBox();
        spinMin->setStyleSheet(transparentSpinStyle);
        spinMin->setRange(isNf ? 1.0 : -9e9, 9e9);
        spinMin->setDecimals(isNf ? 0 : 10);
        spinMin->setValue(p.min);
        if(isEta12) spinMin->setEnabled(false); else spinMin->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 4, spinMin);

        // --- 第 5 列: 寻优上限 (Max) ---
        SmartDoubleSpinBox* spinMax = new SmartDoubleSpinBox();
        spinMax->setStyleSheet(transparentSpinStyle);
        spinMax->setRange(isNf ? 1.0 : -9e9, 9e9);
        spinMax->setDecimals(isNf ? 0 : 10);
        spinMax->setValue(p.max);
        if(isEta12) spinMax->setEnabled(false); else spinMax->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 5, spinMax);

        // --- 第 6 列: 手动滚轮步长 (Step) ---
        SmartDoubleSpinBox* spinStep = new SmartDoubleSpinBox();
        spinStep->setStyleSheet(transparentSpinStyle);
        spinStep->setRange(isNf ? 1.0 : 0.0, 10000.0);
        spinStep->setDecimals(isNf ? 0 : 10);
        spinStep->setValue(p.step);
        if(isEta12) spinStep->setEnabled(false); else spinStep->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 6, spinStep);

        // --- 第 7 列: 界面显示可见性开关 ---
        QWidget* pWidgetVis = new QWidget();
        pWidgetVis->setStyleSheet(transparentWidgetStyle);
        QHBoxLayout* pLayoutVis = new QHBoxLayout(pWidgetVis);
        QCheckBox* chkVis = new QCheckBox();
        chkVis->setChecked(p.isVisible);
        chkVis->setStyleSheet(checkBoxStyle);
        pLayoutVis->addWidget(chkVis);
        pLayoutVis->setAlignment(Qt::AlignCenter);
        pLayoutVis->setContentsMargins(0, 0, 0, 0);
        ui->tableWidget->setCellWidget(i, 7, pWidgetVis);

        // --- 核心联动逻辑 ---
        // 约束：如果某参数被勾选为“参与拟合”，那么它在主界面必须“强制显示”，不可隐藏
        connect(chkFit, &QCheckBox::checkStateChanged, this, [this, chkVis, i](Qt::CheckState state){
            bool isFit = (state == Qt::Checked);
            if (isFit) {
                chkVis->setChecked(true);
                chkVis->setEnabled(false); // 锁死显示开关
                chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
            } else {
                chkVis->setEnabled(true);  // 解锁显示开关
                chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; } QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; } QCheckBox::indicator:hover { border-color: #0078d7; }");
            }
            this->updateRowAppearance(i, isFit);
        });

        // 触发表格行着色初始化
        if (p.isFit && !isEta12) {
            chkVis->setChecked(true);
            chkVis->setEnabled(false);
            chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
        }
        updateRowAppearance(i, p.isFit);
    }

    ui->tableWidget->resizeColumnsToContents();
    ui->tableWidget->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch); // 参数名列自适应拉伸
}

// 内存同步函数
void ParamSelectDialog::collectData()
{
    for(int i = 0; i < ui->tableWidget->rowCount(); ++i) {
        if(i >= m_params.size()) break;

        QWidget* wFit = ui->tableWidget->cellWidget(i, 0);
        if (wFit) {
            QCheckBox* cb = wFit->findChild<QCheckBox*>();
            if(cb) m_params[i].isFit = cb->isChecked();
        }

        QDoubleSpinBox* spinVal = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 1));
        if(spinVal) m_params[i].value = spinVal->value();

        QDoubleSpinBox* spinMin = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 4));
        if(spinMin) m_params[i].min = spinMin->value();

        QDoubleSpinBox* spinMax = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 5));
        if(spinMax) m_params[i].max = spinMax->value();

        QDoubleSpinBox* spinStep = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 6));
        if(spinStep) m_params[i].step = spinStep->value();

        QWidget* wVis = ui->tableWidget->cellWidget(i, 7);
        if (wVis) {
            QCheckBox* cb = wVis->findChild<QCheckBox*>();
            if(cb) m_params[i].isVisible = cb->isChecked();
        }
    }
}

// =========================================================================
// 【核心动态裂缝阵列控制器】: onNfChanged()
// 作用：根据最新的裂缝条数 targetNf，自动增删 Lf_i, xf_i 和 S_i 参数。
// =========================================================================
void ParamSelectDialog::onNfChanged()
{
    // 1. 同步内存，防止 UI 尚未存入 m_params 的修改被冲掉
    collectData();

    // 2. 提取当前面板上的目标裂缝条数
    double currentNf = 1.0;
    for (const auto& p : m_params) {
        if (p.name == "nf") {
            currentNf = p.value;
            break;
        }
    }
    int targetNf = std::max(1, (int)std::round(currentNf));

    // 3. 统计当前结构体中已存在的各种裂缝阵列标识，以推断所需拓扑
    int existingLfCount = 0;
    bool hasXf = false;
    bool hasS  = false;
    double defaultLf = 50.0;
    double defaultS  = 1.0;

    for (const auto& p : m_params) {
        if (p.name.startsWith("Lf_")) existingLfCount++;
        if (p.name == "xf_1") hasXf = true;
        if (p.name == "S_1")  hasS  = true;

        // 抓取第1段的参数作为后续新建段的默认模板，提升用户体验
        if (p.name == "Lf_1") defaultLf = p.value;
        if (p.name == "S_1")  defaultS  = p.value;
    }

    // 若当前条数与目标一致，证明未发生实质物理拓扑变动，直接返回
    if (existingLfCount == targetNf) {
        return;
    }

    // ---------------------------------------------------------------------
    // 步骤 A: 【截断清理】从列表底部向上遍历，删除标号大于 targetNf 的所有越界阵列参数
    // ---------------------------------------------------------------------
    for (int i = m_params.size() - 1; i >= 0; --i) {
        QString name = m_params[i].name;
        int idx = 0;

        // 智能解析尾部索引编号 (如 "Lf_12" 提取出 12)
        if (name.startsWith("Lf_")) idx = name.mid(3).toInt();
        else if (name.startsWith("xf_")) idx = name.mid(3).toInt();
        else if (name.startsWith("S_")) idx = name.mid(2).toInt();

        // 斩断多余的参数尾巴
        if (idx > targetNf) {
            m_params.removeAt(i);
        }
    }

    // ---------------------------------------------------------------------
    // 步骤 B: 【动态注入】定位各变量群的插入点，并补齐缺失索引至 targetNf
    // ---------------------------------------------------------------------

    // (B-1) 补齐裂缝半长 Lf_i
    int insertPosLf = m_params.size();
    for (int i = 0; i < m_params.size(); ++i) {
        if (m_params[i].name.startsWith("Lf_")) insertPosLf = i + 1; // 永远插在同类项最后
    }

    for (int i = 1; i <= targetNf; ++i) {
        QString lfName = QString("Lf_%1").arg(i);
        bool found = false;
        for (const auto& p : m_params) { if (p.name == lfName) { found = true; break; } }

        if (!found) {
            FitParameter pLf;
            pLf.name = lfName;
            pLf.displayName = QString("第%1段裂缝半长").arg(i);
            pLf.value = defaultLf;
            pLf.min = 1.0;
            pLf.max = 2000.0;
            pLf.step = 10.0;
            pLf.isFit = true; // 裂缝长度默认开启参与拟合
            pLf.isVisible = true;
            m_params.insert(insertPosLf++, pLf);
        }
    }

    // (B-2) 若该模型支持非均布，则补齐位置阵列 xf_i
    if (hasXf) {
        int insertPosXf = m_params.size();
        for (int i = 0; i < m_params.size(); ++i) {
            if (m_params[i].name.startsWith("xf_")) insertPosXf = i + 1;
        }
        for (int i = 1; i <= targetNf; ++i) {
            QString xfName = QString("xf_%1").arg(i);
            bool found = false;
            for (const auto& p : m_params) { if (p.name == xfName) { found = true; break; } }

            if (!found) {
                FitParameter pXf;
                pXf.name = xfName;
                pXf.displayName = QString("第%1段裂缝位置").arg(i);
                pXf.value = 500.0;
                pXf.min = 0.0;
                pXf.max = 5000.0;
                pXf.step = 50.0;
                pXf.isFit = false; // 物理几何拓扑位置极易引起雅可比奇异，默认不准拟合
                pXf.isVisible = true;
                m_params.insert(insertPosXf++, pXf);
            }
        }
    }

    // (B-3) 若该模型支持多段独立表皮，则补齐表皮阵列 S_i
    if (hasS) {
        int insertPosS = m_params.size();
        for (int i = 0; i < m_params.size(); ++i) {
            if (m_params[i].name.startsWith("S_")) insertPosS = i + 1;
        }
        for (int i = 1; i <= targetNf; ++i) {
            QString sName = QString("S_%1").arg(i);
            bool found = false;
            for (const auto& p : m_params) { if (p.name == sName) { found = true; break; } }

            if (!found) {
                FitParameter pS;
                pS.name = sName;
                pS.displayName = QString("第%1段表皮系数").arg(i);
                pS.value = defaultS;
                pS.min = -5.0; // 支持负表皮 (酸化解堵效果)
                pS.max = 50.0;
                pS.step = 0.5;
                pS.isFit = true; // 各段污染分布默认允许拟合寻找真值
                pS.isVisible = true;
                m_params.insert(insertPosS++, pS);
            }
        }
    }

    // 4. 重建表单模型渲染 UI
    initTable();
}

QList<FitParameter> ParamSelectDialog::getUpdatedParams() const { return m_params; }
double ParamSelectDialog::getFittingTime() const { return ui->spinTimeMax->value(); }
bool ParamSelectDialog::getUseLimits() const { return ui->chkUseLimits->isChecked(); }

// 确认修改，同步到内存并退出
void ParamSelectDialog::onConfirm() {
    collectData();
    accept();
}

// 放弃修改，直接退出
void ParamSelectDialog::onCancel() {
    reject();
}
