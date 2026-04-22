/*
 * 文件名: fittingparameterchart.cpp
 * 作用与功能:
 * 1. 负责参数配置表格的动态渲染与修改监听。
 * 2. 实现了动态标量展开逻辑，摒弃了易出错的字符串数组。
 * 3. 通过 syncFractureParams 方法，当裂缝条数 nf 变化时，自动增减相应的多段参数行（如 Lf_i, xf_i, S_i），并为 LM 非线性反演引擎提供支持。
 * 4. 提供静态方法进行参数名向中文、单位的转换，及初始默认参数的构建与自适应极限配置。
 */

#include "fittingparameterchart.h"
#include <QHeaderView>
#include <QTableWidgetItem>
#include <QDebug>
#include <QBrush>
#include <QColor>
#include <QRegularExpression>
#include <QWheelEvent>
#include <cmath>

// =========================================================================
// 构造函数与事件过滤：初始化表格、定时器及滚轮修改支持
// =========================================================================

FittingParameterChart::FittingParameterChart(QTableWidget *parentTable, QObject *parent)
    : QObject(parent),
    m_table(parentTable),
    m_modelManager(nullptr),
    m_modelType(ModelManager::Model_1)
{
    // 配置滚轮修改的防抖定时器
    m_wheelTimer = new QTimer(this);
    m_wheelTimer->setSingleShot(true);
    m_wheelTimer->setInterval(200);
    connect(m_wheelTimer, &QTimer::timeout, this, &FittingParameterChart::onWheelDebounceTimeout);

    if(m_table) {
        // 初始化表格基础属性
        QStringList headers;
        headers << "序号" << "参数名称" << "数值" << "单位";
        m_table->setColumnCount(headers.size());
        m_table->setHorizontalHeaderLabels(headers);

        m_table->horizontalHeader()->setStyleSheet(
            "QHeaderView::section { background-color: #E0E0E0; color: black; font-weight: bold; border: 1px solid #A0A0A0; }"
            );
        m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
        m_table->horizontalHeader()->setStretchLastSection(true);

        m_table->setColumnWidth(0, 40);
        m_table->setColumnWidth(1, 160);
        m_table->setColumnWidth(2, 80);

        m_table->setSelectionBehavior(QAbstractItemView::SelectRows);
        m_table->setAlternatingRowColors(false);
        m_table->verticalHeader()->setVisible(false);

        // 安装事件过滤器拦截滚轮修改
        m_table->viewport()->installEventFilter(this);
        // 绑定参数数值的手动编辑信号
        connect(m_table, &QTableWidget::itemChanged, this, &FittingParameterChart::onTableItemChanged);
    }
}

bool FittingParameterChart::eventFilter(QObject *w, QEvent *e)
{
    // 拦截表格的鼠标滚轮事件，实现快捷微调参数
    if (w == m_table->viewport() && e->type() == QEvent::Wheel) {
        QWheelEvent *we = static_cast<QWheelEvent*>(e);
        QTableWidgetItem *item = m_table->itemAt(we->position().toPoint());

        if (item && item->column() == 2) { // 仅允许在“数值”列滚动
            QTableWidgetItem *kItem = m_table->item(item->row(), 1);
            if (!kItem) {
                return false;
            }

            QString pName = kItem->data(Qt::UserRole).toString();
            // 无因次缝长等衍生计算量，禁止滚轮修改
            if (pName == "LfD" || pName == "eta12" || pName == "M12") {
                return true;
            }

            // 遍历参数缓存更新相应的值
            for (auto &p : m_params) {
                if (p.name == pName) {
                    QString txt = item->text();

                    // 防止在未平铺开的旧有阵列表示形式上滚动
                    if (txt.contains(',') || txt.contains(QChar(0xFF0C))) {
                        return false;
                    }

                    bool ok;
                    double cur = txt.toDouble(&ok);

                    if (ok) {
                        double nV = cur + (we->angleDelta().y() / 120) * p.step;

                        // 保证裂缝条数为整数且大于等于1
                        if (p.name == "nf") {
                            nV = qMax(1.0, std::round(nV));
                        }

                        // 保证不超过设定的限幅
                        if (p.max > p.min) {
                            nV = qBound(p.min, nV, p.max);
                        }

                        // 设置文本显示
                        if (p.name == "nf") {
                            item->setText(QString::number(nV));
                        } else {
                            item->setText(QString::number(nV, 'g', 6));
                        }

                        p.value = nV;
                        m_wheelTimer->start();
                        return true;
                    }
                }
            }
        }
    }
    return QObject::eventFilter(w, e);
}

void FittingParameterChart::onWheelDebounceTimeout()
{
    // 防抖发送信号，通知图版更新
    emit parameterChangedByWheel();
}

// =========================================================================
// 用户手动编辑及相关参数联动同步
// =========================================================================

void FittingParameterChart::onTableItemChanged(QTableWidgetItem *item)
{
    if (!item || item->column() != 2) return;

    QTableWidgetItem *kItem = m_table->item(item->row(), 1);
    if (!kItem) return;

    QString key = kItem->data(Qt::UserRole).toString();

    // -------------------------------------------------------------
    // 【核心机制】：当裂缝条数(nf)被修改时，启动阵列重构，动态生成缺失的 Lf_i, xf_i, S_i 变量
    // -------------------------------------------------------------
    if (key == "nf") {
        double val = item->text().toDouble();
        int newNf = qMax(1, (int)std::round(val));

        int oldNf = 1;
        for (const auto& p : m_params) {
            if (p.name == "nf") {
                oldNf = (int)p.value;
                break;
            }
        }

        if (newNf != oldNf) {
            syncFractureParams(newNf);
        }

        if (!m_wheelTimer->isActive()) {
            m_wheelTimer->start();
        }
        return;
    }

    // -------------------------------------------------------------
    // 【衍生参数联动】：若修改了裂缝全长(L)或单端半长(Lf)，更新衍生参数 LfD
    // 注意，此功能仅针对传统 0-107 (等长) 模型有效。
    // -------------------------------------------------------------
    if (key == "L" || key == "Lf") {
        double valL = 0.0;
        double valLf = 0.0;
        QTableWidgetItem* iLfD = nullptr;

        for(int i = 0; i < m_table->rowCount(); ++i) {
            QTableWidgetItem* k = m_table->item(i, 1);
            QTableWidgetItem* v = m_table->item(i, 2);

            if(k && v) {
                QString cK = k->data(Qt::UserRole).toString();
                if (cK == "L") {
                    valL = v->text().toDouble();
                } else if (cK == "Lf") {
                    valLf = v->text().toDouble();
                } else if (cK == "LfD") {
                    iLfD = v;
                }
            }
        }

        if (valL > 1e-9 && iLfD) {
            double nLfD = valLf / valL;
            m_table->blockSignals(true);
            iLfD->setText(QString::number(nLfD, 'g', 6));
            m_table->blockSignals(false);

            for(auto& p : m_params) {
                if(p.name == "LfD") {
                    p.value = nLfD;
                    break;
                }
            }
        }
        if (!m_wheelTimer->isActive()) {
            m_wheelTimer->start();
        }
    }
}

// =========================================================================
// 表格数据的重构、初始化及限幅调整
// =========================================================================

void FittingParameterChart::syncFractureParams(int newNf)
{
    int t = (int)m_modelType;

    // 如果是 0-107 均布等长模型，不存在多条参数展开，直接返回
    if (t < 108) {
        return;
    }

    // 1. 缓存旧的参数状态
    QMap<QString, FitParameter> oldParams;
    for (const auto& p : m_params) {
        oldParams[p.name] = p;
    }

    // 2. 利用 overrideNf 获取包含了全新条数阵列的平铺参数模板
    m_params = generateDefaultParams(m_modelType, newNf);

    // 3. 将缓存中的旧数据原封不动还给新列表，确保用户的输入和已配置好的约束不会在阵列拉伸时丢失
    for (auto& p : m_params) {
        if (oldParams.contains(p.name)) {
            p.value = oldParams[p.name].value;
            p.isFit = oldParams[p.name].isFit;
            p.isVisible = oldParams[p.name].isVisible;
            p.min = oldParams[p.name].min;
            p.max = oldParams[p.name].max;
            p.step = oldParams[p.name].step;
        }
    }

    // 4. 对新增（不包含在缓存中）的阵列元素重新计算上下限
    autoAdjustLimits();

    // 5. 重新渲染整个 UI 参数表
    refreshParamTable();
}

void FittingParameterChart::setModelManager(ModelManager *m)
{
    m_modelManager = m;
}

QList<FitParameter> FittingParameterChart::generateDefaultParams(ModelManager::ModelType type, int overrideNf)
{
    QList<FitParameter> params;
    // 从中枢字典获取基础参数映射
    QMap<QString, double> defs = ModelManager::getDefaultParameters(type);

    // 辅助 lambda 函数：添加单个参入到最终输出列表
    auto addParam = [&](QString name, bool isFitDefault, double fallbackValue = 0.0) {
        FitParameter p;
        p.name = name;
        p.value = defs.contains(name) ? defs.value(name) : fallbackValue;
        p.isFit = isFitDefault;
        p.isVisible = true;
        p.min = p.max = p.step = 0;

        QString dummy;
        getParamDisplayInfo(p.name, p.displayName, dummy, dummy, dummy);
        params.append(p);
    };

    // 第一级核心宏观物理常数
    addParam("phi", false, 0.05);
    addParam("h", false, 20.0);
    addParam("rw", false, 0.1);
    addParam("mu", false, 0.5);
    addParam("B", false, 1.2);
    addParam("Ct", false, 5e-4);
    addParam("q", false, 50.0);

    int t = static_cast<int>(type);
    int tBase = t % 108;

    addParam("kf", true, 10.0);
    addParam("k2", true, 10.0);
    addParam("L", true, 1000.0);

    // 确定当前构建逻辑中所基于的裂缝条数
    int nf = (overrideNf > 0) ? overrideNf : (int)defs.value("nf", 9);

    // ============================================
    // 【核心平铺化展开】: 针对模型 > 108 进行多段参数展开
    // 对于每个独立的参数 (Lf_1, xf_1, S_1...) 分配独立的标量名与初值
    // ============================================

    // 1. 独立裂缝长度展开
    if (t >= 108) {
        for (int i = 1; i <= nf; ++i) {
            addParam(QString("Lf_%1").arg(i), true, 50.0);
        }
    } else {
        addParam("Lf", true, 50.0);
    }

    // 2. 独立裂缝位置绝对坐标展开
    if (t >= 216) {
        double L_total = defs.value("L", 1000.0);
        for (int i = 1; i <= nf; ++i) {
            double default_pos = 0.05 * L_total + 0.9 * L_total * (i - 1) / std::max(1, nf - 1);
            addParam(QString("xf_%1").arg(i), true, default_pos);
        }
    }

    // 3. 将总条数参数放入表单，此为拓扑强约束，通常不直接拟合
    addParam("nf", false, nf);
    addParam("rm", true, 50000.0);

    if (defs.contains("re")) {
        addParam("re", true, 20000.0);
    }

    // 4. 多重介质储容与窜流系数展开
    if (tBase >= 72 && tBase <= 107) {
        addParam("omega_f1", true, 0.02);
        addParam("omega_v1", true, 0.01);
        addParam("lambda_m1", true, 4e-4);
        addParam("lambda_v1", true, 1e-4);
        if (tBase - 72 < 24) {
            addParam("omega_f2", true, 0.008);
            addParam("lambda_m2", true, 1e-7);
        }
    } else if (!(tBase >= 7 && tBase <= 12)) {
        addParam("omega1", true, 0.1);
        addParam("omega2", true, 0.001);
        addParam("lambda1", true, 2e-3);
        addParam("lambda2", true, 1e-3);
    }

    // 5. 井筒与表皮：应对独立表皮进行扩展
    int sType = tBase % 4;
    if (sType != 1) { // 剔除无井储的情况
        addParam("cD", true, 0.1);
        if (t >= 216) {
            // 对每段裂缝赋予单独表皮初值
            for (int i = 1; i <= nf; ++i) {
                addParam(QString("S_%1").arg(i), true, 1.0);
            }
        } else {
            addParam("S", true, 1.0);
        }
    }

    if (sType == 2 || sType == 3) {
        addParam("alpha", false, 0.1);
        addParam("C_phi", false, 1e-4);
    }
    addParam("gamaD", false, 0.006);

    // 等长强制降维开关，用于提升多自变量场景收敛性，默认闭合不拟合
    addParam("use_equal_lf", false, 0.0);

    // 6. 追加展示性质辅助变量 (仅针对 0-107 模型)
    if (t < 108) {
        FitParameter lfd;
        lfd.name = "LfD";
        lfd.displayName = "无因次缝长 LfD";
        lfd.value = defs.value("Lf", 50.0) / defs.value("L", 1000.0);
        lfd.isFit = false;
        lfd.isVisible = true;
        lfd.step = 0;
        params.append(lfd);
    }

    return params;
}

void FittingParameterChart::adjustLimits(QList<FitParameter>& params)
{
    for(auto& p : params) {
        // 对只读/展示变量跳过操作
        if(p.name == "LfD" || p.name == "M12" || p.name == "eta12") continue;

        double val = p.value;
        if (std::abs(val) > 1e-15) {
            p.min = val > 0 ? val * 0.1 : val * 10.0;
            p.max = val > 0 ? val * 10.0 : val * 0.1;
        } else {
            p.min = 0.0;
            p.max = 1.0;
        }

        // --- 强制物理限定 ---
        if (p.name == "phi" || p.name.startsWith("omega")) {
            p.max = qMin(p.max, 1.0);
            p.min = qMax(p.min, 0.0001);
        }

        // xf 绝对坐标拦截负数
        if (p.name.startsWith("xf_")) {
            p.min = 0.0;
            p.max = 20000.0;
        }
        else if (p.name == "kf" || p.name == "k2" || p.name == "km" || p.name == "L" || p.name == "Lf" ||
                 p.name.startsWith("Lf_") || p.name == "rm" || p.name == "re" ||
                 p.name.startsWith("lambda") || p.name == "h" || p.name == "rw" ||
                 p.name == "mu" || p.name == "B" || p.name == "Ct" ||
                 p.name == "cD" || p.name == "q" || p.name == "alpha" || p.name == "C_phi")
        {
            p.min = qMax(p.min, qMax(std::abs(val) * 0.01, 1e-6));
        }

        // 独立表皮系数范围
        else if (p.name == "S" || p.name.startsWith("S_")) {
            p.min = -5.0;
            p.max = 50.0;
        }

        if (p.name == "nf") {
            p.min = qMax(std::ceil(p.min), 1.0);
            p.max = std::floor(p.max);
            p.step = 1.0;
        }

        // 步长智能对数级计算
        if (p.max - p.min > 1e-20 && p.name != "nf") {
            double rs = (p.max - p.min) / 20.0;
            double mag = std::pow(10.0, std::floor(std::log10(rs)));
            p.step = qMax(std::round(rs / mag * 10.0) / 10.0, 0.1) * mag;
        } else if (p.name != "nf") {
            p.step = 0.1;
        }
    }
}

void FittingParameterChart::resetParams(ModelManager::ModelType type, bool preserveStates)
{
    m_modelType = type;
    QMap<QString, QPair<bool, bool>> bkp;

    if (preserveStates) {
        for(const auto& p : m_params) {
            bkp[p.name] = {p.isFit, p.isVisible};
        }
    }

    m_params = generateDefaultParams(type);

    if (preserveStates) {
        for(auto& p : m_params) {
            if(bkp.contains(p.name)) {
                p.isFit = (p.name == "LfD") ? false : bkp[p.name].first;
                p.isVisible = bkp[p.name].second;
            }
        }
    }

    autoAdjustLimits();
    refreshParamTable();
}

void FittingParameterChart::autoAdjustLimits()
{
    adjustLimits(m_params);
}

QList<FitParameter> FittingParameterChart::getParameters() const
{
    return m_params;
}

void FittingParameterChart::setParameters(const QList<FitParameter> &p)
{
    m_params = p;
    refreshParamTable();
}

void FittingParameterChart::switchModel(ModelManager::ModelType newType)
{
    m_modelType = newType;
    QMap<QString, double> old;
    for(const auto& p : m_params) {
        old.insert(p.name, p.value);
    }

    resetParams(newType, false);

    for(auto& p : m_params) {
        if(old.contains(p.name)) {
            p.value = old[p.name];
        }
    }

    autoAdjustLimits();

    // 根据切换后继承过来的 L ，智能收紧相关派生参数界限
    double cL = 1000.0;
    for(const auto& p : m_params) {
        if(p.name == "L") {
            cL = p.value;
        }
    }

    for(auto& p : m_params) {
        if(p.name == "rm") {
            p.min = qMax(p.min, cL);
            p.value = qMax(p.value, p.min);
        }
        if(p.name == "LfD") {
            double cLf = 20.0;
            for(const auto& pp : m_params) {
                if(pp.name == "Lf") cLf = pp.value;
            }
            if(cL > 1e-9) {
                p.value = cLf / cL;
            }
        }
    }
    refreshParamTable();
}

void FittingParameterChart::updateParamsFromTable()
{
    if(!m_table) return;

    for(int i = 0; i < m_table->rowCount(); ++i) {
        QTableWidgetItem* iK = m_table->item(i, 1);
        if(!iK) continue;

        QTableWidgetItem* iV = m_table->item(i, 2);
        QString k = iK->data(Qt::UserRole).toString();

        // 统一标量字符串双精度转型
        double v = iV->text().toDouble();

        if (k == "nf") {
            v = qMax(1.0, std::round(v));
        }

        for(auto& p : m_params) {
            if(p.name == k) {
                p.value = v;
                break;
            }
        }
    }
}

QMap<QString, QString> FittingParameterChart::getRawParamTexts() const
{
    QMap<QString, QString> res;
    if(!m_table) return res;

    for(int i = 0; i < m_table->rowCount(); ++i) {
        QTableWidgetItem* k = m_table->item(i, 1);
        QTableWidgetItem* v = m_table->item(i, 2);
        if (k && v) {
            res.insert(k->data(Qt::UserRole).toString(), v->text());
        }
    }
    return res;
}

// =========================================================================
// 渲染功能核心：将内部配置数据按顺序反映到 QTableWidget 界面
// =========================================================================

void FittingParameterChart::refreshParamTable()
{
    if(!m_table) return;

    m_table->blockSignals(true);
    m_table->setRowCount(0);

    int no = 1;
    // 分为两层输出：被设置为"拟合(isFit)"的在顶端高亮显示，作为固定参数的不高亮在底部显示
    for(const auto& p : m_params) {
        if(p.isVisible && p.isFit) {
            addRowToTable(p, no, true);
        }
    }

    for(const auto& p : m_params) {
        if(p.isVisible && !p.isFit) {
            addRowToTable(p, no, false);
        }
    }

    m_table->blockSignals(false);
}

void FittingParameterChart::addRowToTable(const FitParameter& p, int& serialNo, bool highlight)
{
    int r = m_table->rowCount();
    m_table->insertRow(r);

    // LfD, M12 等不参与计算的仅作为参考，给出灰度背景区分
    QColor bg = (p.name == "LfD" || p.name == "M12" || p.name == "eta12") ? QColor(245, 245, 245) : (highlight ? QColor(255, 255, 224) : Qt::white);

    QTableWidgetItem* i0 = new QTableWidgetItem(QString::number(serialNo++));
    i0->setFlags(i0->flags() & ~Qt::ItemIsEditable);
    i0->setTextAlignment(Qt::AlignCenter);
    i0->setBackground(bg);
    m_table->setItem(r, 0, i0);

    QTableWidgetItem* i1 = new QTableWidgetItem(p.displayName);
    i1->setFlags(i1->flags() & ~Qt::ItemIsEditable);
    i1->setData(Qt::UserRole, p.name);
    i1->setBackground(bg);

    if(highlight) {
        QFont f = i1->font();
        f.setBold(true);
        i1->setFont(f);
    }
    m_table->setItem(r, 1, i1);

    QTableWidgetItem* i2 = nullptr;
    if (p.name == "nf") {
        i2 = new QTableWidgetItem(QString::number(qMax(1.0, std::round(p.value))));
    } else {
        i2 = new QTableWidgetItem(QString::number(p.value, 'g', 6));
    }

    i2->setBackground(bg);
    if(highlight) {
        QFont f = i2->font();
        f.setBold(true);
        i2->setFont(f);
    }
    if (p.name == "LfD" || p.name == "M12" || p.name == "eta12") {
        i2->setFlags(i2->flags() & ~Qt::ItemIsEditable);
        i2->setForeground(QBrush(Qt::darkGray));
    }
    m_table->setItem(r, 2, i2);

    QString d, u;
    getParamDisplayInfo(p.name, d, d, d, u);

    QTableWidgetItem* i3 = new QTableWidgetItem((u == "无因次" || u == "小数") ? "-" : u);
    i3->setFlags(i3->flags() & ~Qt::ItemIsEditable);
    i3->setBackground(bg);
    m_table->setItem(r, 3, i3);
}

// =========================================================================
// 核心中文与单位翻译官
// 兼容动态正则映射，为 "S_i" 阵列指派中文显示说明
// =========================================================================
void FittingParameterChart::getParamDisplayInfo(const QString &n, QString &cN, QString &s, QString &uS, QString &u)
{
    if(n == "kf") { cN = "内区渗透率 kf"; u = "mD"; }
    else if(n == "k2" || n == "km") { cN = "外区渗透率 k2"; u = "mD"; }
    else if(n == "L") { cN = "水平井长 L"; u = "m"; }
    else if(n == "Lf") { cN = "统一裂缝半长 Lf"; u = "m"; }
    else if(n == "rm") { cN = "复合内区半径 rm"; u = "m"; }
    else if(n == "omega1") { cN = "内区储容比 ω₁"; u = "无因次"; }
    else if(n == "omega2") { cN = "外区储容比 ω₂"; u = "无因次"; }
    else if(n == "lambda1") { cN = "内区窜流系数 λ₁"; u = "无因次"; }
    else if(n == "lambda2") { cN = "外区窜流系数 λ₂"; u = "无因次"; }
    else if(n == "omega_f1") { cN = "内区裂缝储容比 ωf₁"; u = "无因次"; }
    else if(n == "omega_v1") { cN = "内区溶洞储容比 ωv₁"; u = "无因次"; }
    else if(n == "lambda_m1") { cN = "内区基质窜流系数 λm₁"; u = "无因次"; }
    else if(n == "lambda_v1") { cN = "内区溶洞窜流系数 λv₁"; u = "无因次"; }
    else if(n == "omega_f2") { cN = "外区裂缝储容比 ωf₂"; u = "无因次"; }
    else if(n == "lambda_m2") { cN = "外区基质窜流系数 λm₂"; u = "无因次"; }
    else if(n == "re") { cN = "外区定压/封闭半径 re"; u = "m"; }
    else if(n == "nf") { cN = "裂缝条数 nf"; u = "条"; }
    else if(n == "use_equal_lf") { cN = "防发散强制等长降维开关"; u = "无因次"; }

    // 【核心阵列适配】：对多维度动态阵列给出友善的下标拼接命名
    else if(n.startsWith("Lf_")) {
        int idx = n.mid(3).toInt();
        cN = QString("第%1段裂缝半长 Lf_%2").arg(idx).arg(idx);
        u = "m";
    }
    else if(n.startsWith("xf_")) {
        int idx = n.mid(3).toInt();
        cN = QString("第%1段裂缝绝对位置 xf_%2").arg(idx).arg(idx);
        u = "m";
    }
    else if(n.startsWith("S_")) {
        int idx = n.mid(2).toInt();
        cN = QString("第%1段独立表皮系数 S_%2").arg(idx).arg(idx);
        u = "无因次";
    }

    else if(n == "h") { cN = "有效厚度 h"; u = "m"; }
    else if(n == "rw") { cN = "井筒半径 rw"; u = "m"; }
    else if(n == "phi") { cN = "孔隙度 φ"; u = "小数"; }
    else if(n == "mu") { cN = "流体粘度 μ"; u = "mPa·s"; }
    else if(n == "B") { cN = "原油体积系数 B"; u = "无因次"; }
    else if(n == "Ct") { cN = "系统综合压缩系数 Ct"; u = "MPa⁻¹"; }
    else if(n == "q") { cN = "测试地面恒量产液 q"; u = "m³/d"; }
    else if(n == "cD") { cN = "无因次井储系数 cD"; u = "无因次"; }
    else if(n == "S") { cN = "统一全局表皮系数 S"; u = "无因次"; }
    else if(n == "gamaD") { cN = "无因次压敏渗透模量 γD"; u = "无因次"; }
    else if(n == "LfD") { cN = "无因次裂缝半长 LfD"; u = "无因次"; }
    else if(n == "alpha") { cN = "变井储相态重分时间参数 α"; u = "h"; }
    else if(n == "C_phi") { cN = "变井储相态重分压力参数 Cφ"; u = "MPa"; }
    else { cN = n; u = ""; }

    s = uS = n;
}
