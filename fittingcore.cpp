/*
 * 文件名: fittingcore.cpp
 * 文件作用: 试井拟合核心算法实现
 * * 核心升级记录:
 * 1. 【多表皮追踪】：在雅可比计算与 LM 增量步更新时，加入对 pName.startsWith("S") 的识别，
 * 确保所有表皮变量采用线性步长和线性扰动（因为表皮可为负），而不是传统的对数步长。
 * 2. 【多起点探测升维】：在广域搜索中，针对独立表皮 S_i 专门定制了区间扰动策略，防止初始化坍塌。
 * 3. 【无井筒模型保护】：在 preprocessParams 中增加逻辑，如果是线源解模型，自动将其余多段 S_i 强行置零，屏蔽自由度。
 */

#include "fittingcore.h"
#include "modelparameter.h"
#include <QtConcurrent>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>

FittingCore::FittingCore(QObject *parent)
    : QObject(parent), m_modelManager(nullptr), m_isCustomSamplingEnabled(false), m_stopRequested(false)
{
    // 绑定线程结束事件
    connect(&m_watcher, &QFutureWatcher<void>::finished, this, &FittingCore::sigFitFinished);
}

void FittingCore::setModelManager(ModelManager *m) { m_modelManager = m; }
void FittingCore::setObservedData(const QVector<double> &t, const QVector<double> &p, const QVector<double> &d) {
    m_obsTime = t; m_obsDeltaP = p; m_obsDerivative = d;
}
void FittingCore::setSamplingSettings(const QList<SamplingInterval> &intervals, bool enabled) {
    m_customIntervals = intervals; m_isCustomSamplingEnabled = enabled;
}

void FittingCore::startFit(ModelManager::ModelType modelType, const QList<FitParameter> &params, double weight, bool useLimits, bool useMultiStart) {
    if (m_watcher.isRunning()) return;
    m_stopRequested = false;
    m_watcher.setFuture(QtConcurrent::run([this, modelType, params, weight, useLimits, useMultiStart](){
        runOptimizationTask(modelType, params, weight, useLimits, useMultiStart);
    }));
}

void FittingCore::stopFit() { m_stopRequested = true; }

// ======================= 数据抽样算子 =======================
void FittingCore::getLogSampledData(const QVector<double>& srcT, const QVector<double>& srcP, const QVector<double>& srcD, QVector<double>& outT, QVector<double>& outP, QVector<double>& outD)
{
    outT.clear(); outP.clear(); outD.clear();
    if (srcT.isEmpty()) return;

    struct DataPoint {
        double t, p, d;
        bool operator<(const DataPoint& other) const { return t < other.t; }
        bool operator==(const DataPoint& other) const { return std::abs(t - other.t) < 1e-9; }
    };
    QVector<DataPoint> points;

    if (!m_isCustomSamplingEnabled) {
        int targetCount = 80;
        if (srcT.size() <= targetCount) {
            outT = srcT; outP = srcP; outD = srcD;
            return;
        }
        double tMin = srcT.first() <= 1e-10 ? 1e-4 : srcT.first();
        double tMax = srcT.last();
        double logMin = log10(tMin), logMax = log10(tMax);
        double step = (logMax - logMin) / (targetCount - 1);

        int currentIndex = 0;
        for (int i = 0; i < targetCount; ++i) {
            double targetT = pow(10, logMin + i * step);
            double minDiff = 1e30;
            int bestIdx = currentIndex;
            while (currentIndex < srcT.size()) {
                double diff = std::abs(srcT[currentIndex] - targetT);
                if (diff < minDiff) { minDiff = diff; bestIdx = currentIndex; }
                else break;
                currentIndex++;
            }
            currentIndex = bestIdx;
            points.append({srcT[bestIdx], (bestIdx<srcP.size()?srcP[bestIdx]:0.0), (bestIdx<srcD.size()?srcD[bestIdx]:0.0)});
        }
    } else {
        if (m_customIntervals.isEmpty()) {
            outT = srcT; outP = srcP; outD = srcD;
            return;
        }
        for (const auto& interval : m_customIntervals) {
            double tStart = interval.tStart, tEnd = interval.tEnd;
            int count = interval.count;
            if (count <= 0) continue;

            auto itStart = std::lower_bound(srcT.begin(), srcT.end(), tStart);
            auto itEnd = std::upper_bound(srcT.begin(), srcT.end(), tEnd);
            int idxStart = std::distance(srcT.begin(), itStart);
            int idxEnd = std::distance(srcT.begin(), itEnd);

            if (idxStart >= srcT.size() || idxStart >= idxEnd) continue;

            double subMin = srcT[idxStart], subMax = srcT[idxEnd - 1];
            if (subMin <= 1e-10) subMin = 1e-4;

            double logMin = log10(subMin), logMax = log10(subMax);
            double step = (count > 1) ? (logMax - logMin) / (count - 1) : 0;

            int subCurrentIdx = idxStart;
            for (int i = 0; i < count; ++i) {
                double targetT = (count == 1) ? subMin : pow(10, logMin + i * step);
                double minDiff = 1e30;
                int bestIdx = subCurrentIdx;
                while (subCurrentIdx < idxEnd) {
                    double diff = std::abs(srcT[subCurrentIdx] - targetT);
                    if (diff < minDiff) { minDiff = diff; bestIdx = subCurrentIdx; }
                    else break;
                    subCurrentIdx++;
                }
                subCurrentIdx = bestIdx;
                if (bestIdx < srcT.size()) {
                    points.append({srcT[bestIdx], (bestIdx<srcP.size()?srcP[bestIdx]:0.0), (bestIdx<srcD.size()?srcD[bestIdx]:0.0)});
                }
            }
        }
    }

    std::sort(points.begin(), points.end());
    auto last = std::unique(points.begin(), points.end());
    points.erase(last, points.end());

    for (const auto& p : points) { outT.append(p.t); outP.append(p.p); outD.append(p.d); }
}

// ======================= 预清洗与降维保护 =======================
QMap<QString, double> FittingCore::preprocessParams(const QMap<QString, double>& rawParams, ModelManager::ModelType type)
{
    QMap<QString, double> processed = rawParams;
    ModelParameter* mp = ModelParameter::instance();

    auto getSafeParam = [&](const QString& key, double mpVal, double defaultVal) {
        if (rawParams.contains(key)) return rawParams[key];
        if (std::abs(mpVal) > 1e-15) return mpVal;
        return defaultVal;
    };

    double phi = getSafeParam("phi", mp->getPhi(), 0.05);
    double h   = getSafeParam("h",   mp->getH(),   20.0);
    double Ct  = getSafeParam("Ct",  mp->getCt(),  5e-4);
    double mu  = getSafeParam("mu",  mp->getMu(),  0.5);
    double B   = getSafeParam("B",   mp->getB(),   1.05);
    double q   = getSafeParam("q",   mp->getQ(),   5.0);
    double rw  = getSafeParam("rw",  mp->getRw(),  0.1);

    double nf  = getSafeParam("nf",  mp->getNf(),  9.0);
    int nf_val = std::max(1, (int)std::round(nf));

    processed["phi"] = phi; processed["h"] = h; processed["Ct"] = Ct;
    processed["mu"] = mu;   processed["B"] = B; processed["q"] = q;
    processed["rw"] = rw;   processed["nf"] = nf_val;

    double L = processed.value("L", 0.0);
    if (L < 1e-9) { L = 1000.0; processed["L"] = L; }

    int typeId = (int)type;
    bool isNonEqual = (typeId >= 108);

    if (!isNonEqual) {
        if (processed.contains("Lf")) processed["LfD"] = processed["Lf"] / L;
        else processed["LfD"] = 0.0;
    }

    // 【核心约束1】全缝等长强制绑定
    if (processed.value("use_equal_lf", 0.0) > 0.5) {
        double baseLf = processed.value("Lf_1", 50.0);
        for (int i = 1; i <= nf_val; ++i) {
            QString key = QString("Lf_%1").arg(i);
            if (processed.contains(key)) processed[key] = baseLf;
        }
    }

    if (!processed.contains("M12") && processed.contains("km")) processed["M12"] = processed["km"];

    // 【核心约束2】对无井筒效应模型执行独立表皮置零，封死无意义的求解维度
    int tBase = isNonEqual ? typeId - 108 : typeId;
    int storageType = tBase % 4;
    bool hasStorage = (storageType != 1);

    if (hasStorage) {
        if (!processed.contains("cD")) processed["cD"] = 0.1;
    } else {
        processed["cD"] = 0.0;
        processed["S"] = 0.0;
        for (int i = 1; i <= nf_val; ++i) {
            QString sKey = QString("S_%1").arg(i);
            processed[sKey] = 0.0;
        }
    }

    int gI = tBase % 12;
    if (gI < 4) { if (!processed.contains("re")) processed["re"] = 20000.0; }

    return processed;
}

void FittingCore::runOptimizationTask(ModelManager::ModelType modelType, QList<FitParameter> fitParams, double weight, bool useLimits, bool useMultiStart) {
    runLevenbergMarquardtOptimization(modelType, fitParams, weight, useLimits, useMultiStart);
}

// ======================= LM 核心优化调度 =======================
void FittingCore::runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params, double weight, bool useLimits, bool useMultiStart) {
    if(m_modelManager) m_modelManager->setHighPrecision(false);

    QVector<int> fitIndices;
    for(int i=0; i<params.size(); ++i) {
        if(params[i].isFit && params[i].name != "LfD" && params[i].name != "use_equal_lf") {
            fitIndices.append(i);
        }
    }
    int nParams = fitIndices.size();

    QMap<QString, double> currentParamMap;
    for(const auto& p : params) currentParamMap.insert(p.name, p.value);

    QVector<double> fitT, fitP, fitD;
    getLogSampledData(m_obsTime, m_obsDeltaP, m_obsDerivative, fitT, fitP, fitD);

    QVector<double> residuals = calculateResiduals(currentParamMap, modelType, weight, fitT, fitP, fitD);
    double currentSSE = calculateSumSquaredError(residuals);

    // --- 广域探索（多起点随机扰动防局部最优） ---
    if (useMultiStart && nParams > 0 && !m_stopRequested) {
        QMap<QString, double> bestParamMap = currentParamMap;
        double bestSSE = currentSSE;

        std::srand(static_cast<unsigned>(std::time(nullptr)));
        const int numTrials = 5;

        for (int trial = 0; trial < numTrials && !m_stopRequested; ++trial) {
            QMap<QString, double> trialMap = currentParamMap;
            for (int i = 0; i < nParams; ++i) {
                int pIdx = fitIndices[i];
                QString pName = params[pIdx].name;
                double val = currentParamMap[pName];
                double logPert = (double(std::rand()) / RAND_MAX) * 2.0 - 1.0;

                // 【独立表皮识别】S_i 支持线性变动
                if (pName == "nf") {
                    trialMap[pName] = std::max(1.0, std::round(val + logPert * 5));
                } else if (pName.startsWith("S")) {
                    trialMap[pName] = val + logPert * 5.0; // 表皮线性扰动
                } else if (val > 1e-12) {
                    trialMap[pName] = val * std::pow(10.0, logPert);
                }

                if (useLimits) {
                    trialMap[pName] = qMax(params[pIdx].min, qMin(trialMap[pName], params[pIdx].max));
                }
            }

            QVector<double> trialRes = calculateResiduals(trialMap, modelType, weight, fitT, fitP, fitD);
            double trialSSE = calculateSumSquaredError(trialRes);

            if (trialSSE < bestSSE && !std::isnan(trialSSE)) {
                bestSSE = trialSSE;
                bestParamMap = trialMap;
            }
        }
        currentParamMap = bestParamMap;
        currentSSE = bestSSE;
        residuals = calculateResiduals(currentParamMap, modelType, weight, fitT, fitP, fitD);
    }

    QMap<QString, double> solverParams = preprocessParams(currentParamMap, modelType);
    ModelCurveData curve = m_modelManager->calculateTheoreticalCurve(modelType, solverParams);
    emit sigIterationUpdated(currentSSE/residuals.size(), currentParamMap, std::get<0>(curve), std::get<1>(curve), std::get<2>(curve));

    if(nParams == 0) { emit sigFitFinished(); return; }

    double lambda = 0.1;
    int maxIter = 80;
    double prevSSE = currentSSE;
    int stagnationCount = 0;
    int consecutiveFailures = 0;

    for(int iter = 0; iter < maxIter; ++iter) {
        if(m_stopRequested) break;
        if (!residuals.isEmpty() && (currentSSE / residuals.size()) < 1e-4) break;

        emit sigProgress(iter * 100 / maxIter);

        QVector<QVector<double>> J = computeJacobian(currentParamMap, residuals, fitIndices, modelType, params, weight, fitT, fitP, fitD);
        if(m_stopRequested) break;

        int nRes = residuals.size();
        QVector<QVector<double>> H(nParams, QVector<double>(nParams, 0.0));
        QVector<double> g(nParams, 0.0);

        for(int k=0; k<nRes; ++k) {
            for(int i=0; i<nParams; ++i) {
                g[i] += J[k][i] * residuals[k];
                for(int j=0; j<=i; ++j) H[i][j] += J[k][i] * J[k][j];
            }
        }
        for(int i=0; i<nParams; ++i) {
            for(int j=i+1; j<nParams; ++j) H[i][j] = H[j][i];
        }

        bool stepAccepted = false;
        // 内部搜索循环：若步长引起残差恶化，则增大阻尼并收缩步长，最多试错 20 次
        for(int tryIter=0; tryIter<20; ++tryIter) {
            if(m_stopRequested) break;

            QVector<QVector<double>> H_lm = H;
            for(int i=0; i<nParams; ++i) H_lm[i][i] += lambda * (1.0 + std::abs(H[i][i]));

            QVector<double> negG(nParams);
            for(int i=0; i<nParams; ++i) negG[i] = -g[i];

            QVector<double> delta = solveLinearSystem(H_lm, negG);
            if (delta.isEmpty()) { lambda *= 10.0; continue; }

            QMap<QString, double> trialMap = currentParamMap;

            for(int i=0; i<nParams; ++i) {
                int pIdx = fitIndices[i];
                QString pName = params[pIdx].name;
                double oldVal = currentParamMap[pName];

                // 【重要修复】：多段表皮特征识别，阻截非法的负值取对数运算
                bool isLog = (oldVal > 1e-12 && !pName.startsWith("S") && pName != "nf");

                double newVal;
                if(isLog) newVal = pow(10.0, log10(oldVal) + delta[i]);
                else newVal = oldVal + delta[i];

                if (std::isnan(newVal) || std::isinf(newVal)) newVal = oldVal;
                if (useLimits) newVal = qMax(params[pIdx].min, qMin(newVal, params[pIdx].max));

                trialMap[pName] = newVal;
            }

            QVector<double> newRes = calculateResiduals(trialMap, modelType, weight, fitT, fitP, fitD);
            double newSSE = calculateSumSquaredError(newRes);

            if(newSSE < currentSSE && !std::isnan(newSSE)) {
                currentSSE = newSSE;
                currentParamMap = trialMap;
                residuals = newRes;
                lambda = std::max(1e-7, lambda / 10.0);
                stepAccepted = true;

                double relChange = (prevSSE > 1e-20) ? std::abs(prevSSE - currentSSE) / prevSSE : 0.0;
                if (relChange < 1e-4) stagnationCount++; else stagnationCount = 0;
                prevSSE = currentSSE;
                consecutiveFailures = 0;

                QMap<QString, double> trialSolverParams = preprocessParams(trialMap, modelType);
                ModelCurveData iterCurve = m_modelManager->calculateTheoreticalCurve(modelType, trialSolverParams);
                emit sigIterationUpdated(currentSSE/nRes, currentParamMap, std::get<0>(iterCurve), std::get<1>(iterCurve), std::get<2>(iterCurve));
                break;
            } else {
                lambda *= 10.0;
            }
        }

        if(!stepAccepted) {
            consecutiveFailures++;
            if (consecutiveFailures >= 3) break;
            lambda = std::min(lambda * 100.0, 1e10);
        }
        if (stagnationCount >= 2) break;
    }

    if(m_modelManager) m_modelManager->setHighPrecision(true);
    if (currentParamMap.contains("nf")) currentParamMap["nf"] = std::max(1.0, std::round(currentParamMap["nf"]));

    QMap<QString, double> finalSolverParams = preprocessParams(currentParamMap, modelType);
    ModelCurveData finalCurve = m_modelManager->calculateTheoreticalCurve(modelType, finalSolverParams);
    emit sigIterationUpdated(currentSSE/residuals.size(), currentParamMap, std::get<0>(finalCurve), std::get<1>(finalCurve), std::get<2>(finalCurve));
}

// ======================= 目标函数与雅可比计算 =======================
QVector<double> FittingCore::calculateResiduals(const QMap<QString, double>& params, ModelManager::ModelType modelType, double weight,
                                                const QVector<double>& t, const QVector<double>& obsP, const QVector<double>& obsD) {
    if(!m_modelManager || t.isEmpty()) return QVector<double>();

    QMap<QString, double> solverParams = preprocessParams(params, modelType);
    ModelCurveData res = m_modelManager->calculateTheoreticalCurve(modelType, solverParams, t);
    const QVector<double>& pCal = std::get<1>(res);
    const QVector<double>& dpCal = std::get<2>(res);

    QVector<double> r;
    double wp = weight, wd = 1.0 - weight;

    auto safeResidual = [](double obs, double calc, double w) -> double {
        if (obs <= 1e-10) return 0.0;
        double calcEff = calc;
        if (std::isnan(calcEff) || std::isinf(calcEff) || calcEff <= 0.0) calcEff = 1e-20;
        else if (calcEff < 1e-20) calcEff = 1e-20;
        return (std::log(obs) - std::log(calcEff)) * w;
    };

    int count = qMin((int)obsP.size(), (int)pCal.size());
    for(int i=0; i<count; ++i) r.append(safeResidual(obsP[i], pCal[i], wp));
    int dCount = qMin((int)obsD.size(), (int)dpCal.size());
    dCount = qMin(dCount, count);
    for(int i=0; i<dCount; ++i) r.append(safeResidual(obsD[i], dpCal[i], wd));

    // 【正则化约束体系】：非等长模式下施加基于方差的 L2 伪残差镇压项
    if (params.value("use_equal_lf", 0.0) < 0.5) {
        QVector<double> lf_vals;
        for (auto it = params.constBegin(); it != params.constEnd(); ++it) {
            if (it.key().startsWith("Lf_") && it.key() != "LfD") lf_vals.append(it.value());
        }
        if (lf_vals.size() > 1) {
            double mean_lf = 0.0;
            for (double v : lf_vals) mean_lf += v;
            mean_lf /= lf_vals.size();
            double penalty_weight = 0.005; // 刚性因数，镇定多维振荡
            for (double v : lf_vals) r.append((v - mean_lf) * penalty_weight);
        }
    }
    return r;
}

QVector<QVector<double>> FittingCore::computeJacobian(const QMap<QString, double>& params, const QVector<double>& baseResiduals,
                                                      const QVector<int>& fitIndices, ModelManager::ModelType modelType,
                                                      const QList<FitParameter>& currentFitParams, double weight,
                                                      const QVector<double>& t, const QVector<double>& obsP, const QVector<double>& obsD) {
    int nRes = baseResiduals.size();
    int nParams = fitIndices.size();
    QVector<QVector<double>> J(nRes, QVector<double>(nParams, 0.0));

    for(int j = 0; j < nParams; ++j) {
        if (m_stopRequested) break;

        int idx = fitIndices[j];
        QString pName = currentFitParams[idx].name;
        double val = params.value(pName);

        // 【独立表皮识别】S_i 取消对数差分，使用线性微扰步长
        bool isLog = (val > 1e-12 && !pName.startsWith("S") && pName != "nf");

        double h;
        QMap<QString, double> pPlus = params;
        QMap<QString, double> pMinus = params;

        if (pName == "nf") {
            h = 0.51;
            pPlus[pName] = val + h; pMinus[pName] = val - h;
        } else if(isLog) {
            h = 0.01;
            double valLog = log10(val);
            pPlus[pName] = pow(10.0, valLog + h); pMinus[pName] = pow(10.0, valLog - h);
        } else {
            h = 1e-4;
            pPlus[pName] = val + h; pMinus[pName] = val - h;
        }

        QVector<double> rPlus = this->calculateResiduals(pPlus, modelType, weight, t, obsP, obsD);
        QVector<double> rMinus = this->calculateResiduals(pMinus, modelType, weight, t, obsP, obsD);

        if(rPlus.size() == nRes && rMinus.size() == nRes) {
            for(int i=0; i<nRes; ++i) J[i][j] = (rPlus[i] - rMinus[i]) / (2.0 * h);
        }
    }
    return J;
}

QVector<double> FittingCore::solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b) {
    int n = b.size();
    if (n == 0) return QVector<double>();
    Eigen::MatrixXd matA(n, n);
    Eigen::VectorXd vecB(n);
    for (int i = 0; i < n; ++i) {
        vecB(i) = b[i];
        for (int j = 0; j < n; ++j) matA(i, j) = A[i][j];
    }
    Eigen::VectorXd x = matA.ldlt().solve(vecB);
    QVector<double> res(n);
    for (int i = 0; i < n; ++i) res[i] = x(i);
    return res;
}

double FittingCore::calculateSumSquaredError(const QVector<double>& residuals) {
    double sse = 0.0;
    for(double v : residuals) sse += v*v;
    return sse;
}

// ======================= 智能初值启发器 =======================
QMap<QString, double> FittingCore::estimateInitialParams(
    const QVector<double>& obsTime,
    const QVector<double>& obsDeltaP,
    const QVector<double>& obsDerivative,
    ModelManager::ModelType modelType)
{
    QMap<QString, double> estimated;
    ModelParameter* mp = ModelParameter::instance();

    double q = mp->getQ(), mu = mp->getMu(), B = mp->getB(), h = mp->getH();
    double phi = mp->getPhi(), Ct = mp->getCt(), L = mp->getL();

    int n = qMin(qMin(obsTime.size(), obsDeltaP.size()), obsDerivative.size());
    if (n < 10) return estimated;

    // 提取导数中段平台估计渗透率 (kf)
    int startIdx = n * 3 / 10, endIdx = n * 7 / 10;
    QVector<double> plateauDerivs;
    for (int i = startIdx; i < endIdx; ++i) {
        if (obsDerivative[i] > 1e-10) plateauDerivs.append(obsDerivative[i]);
    }
    double derivPlateau = 0.0;
    if (!plateauDerivs.isEmpty()) {
        std::sort(plateauDerivs.begin(), plateauDerivs.end());
        derivPlateau = plateauDerivs[plateauDerivs.size() / 2];
    }

    double kf_est = 10.0;
    if (derivPlateau > 1e-15 && q > 0 && mu > 0 && B > 0 && h > 0) {
        kf_est = 0.921 * q * mu * B / (h * derivPlateau);
        kf_est = qBound(0.01, kf_est, 100000.0);
    }
    estimated["kf"] = kf_est;

    int storageType = (int)modelType % 4;
    // 早期数据提取井储 (cD)
    if (storageType != 1) {
        int earlyEnd = qMin(n / 5, 20);
        double bestC = 10.0, bestSlopeErr = 1e30;
        for (int i = 1; i < earlyEnd; ++i) {
            if (obsTime[i] > 1e-10 && obsDeltaP[i] > 1e-10 && obsTime[i-1] > 1e-10 && obsDeltaP[i-1] > 1e-10) {
                double slope = (std::log10(obsDeltaP[i]) - std::log10(obsDeltaP[i-1])) / (std::log10(obsTime[i]) - std::log10(obsTime[i-1]));
                double slopeErr = std::abs(slope - 1.0);
                if (slopeErr < bestSlopeErr) {
                    bestSlopeErr = slopeErr;
                    bestC = q * B * obsTime[i] / (24.0 * obsDeltaP[i]);
                }
            }
        }
        double L_half = L / 2.0;
        double denom = phi * h * Ct * L_half * L_half;
        double cD_est = (denom > 1e-20) ? 0.159 * bestC / denom : 0.1;
        estimated["cD"] = qBound(1e-4, cD_est, 10.0);
    }

    // 从最高峰提取表皮系数 (S)
    if (storageType != 1 && derivPlateau > 1e-15) {
        double maxDeriv = 0.0;
        for (int i = 0; i < n; ++i) { if (obsDerivative[i] > maxDeriv) maxDeriv = obsDerivative[i]; }
        double S_est = 1.0;
        if (maxDeriv > derivPlateau * 1.1) S_est = 0.5 * (maxDeriv / derivPlateau - 1.0);
        S_est = qBound(-5.0, S_est, 50.0);

        // 宏观映射：如果是模型 07~09，将宏观估算表皮赋予所有多段裂缝
        estimated["S"] = S_est;
        int nf_est = std::max(1, (int)std::round(mp->getNf()));
        for(int i = 1; i <= nf_est; i++) {
            estimated[QString("S_%1").arg(i)] = S_est;
        }
    }

    // 双孔凹陷特征提取 (omega / lambda)
    double minDeriv = 1e30;
    int minIdx = -1;
    for (int i = startIdx; i < endIdx; ++i) {
        if (obsDerivative[i] > 1e-10 && obsDerivative[i] < minDeriv) { minDeriv = obsDerivative[i]; minIdx = i; }
    }

    int modelId = (int)modelType;
    bool isTriple = (modelId >= 72 && modelId <= 107) || (modelId >= 180 && modelId <= 215) || (modelId >= 288);

    if (derivPlateau > 1e-15 && minDeriv < derivPlateau * 0.9 && minIdx > 0) {
        double omega_est = minDeriv / derivPlateau;
        omega_est = qBound(0.001, omega_est, 0.5);

        if (isTriple) {
            estimated["omega_f1"] = omega_est;
            estimated["omega_v1"] = omega_est * 0.5;
        } else {
            estimated["omega1"] = omega_est;
        }

        if (kf_est > 0 && phi > 0 && mu > 0 && Ct > 0 && L > 0) {
            double tDip = obsTime[minIdx];
            double tD_dip = 3.6 * kf_est * tDip / (phi * mu * Ct * std::pow(L/2.0, 2.0));
            double lambda_est = (tD_dip > 1e-20) ? 1.0 / tD_dip : 1e-3;
            lambda_est = qBound(1e-9, lambda_est, 1.0);
            if (isTriple) {
                estimated["lambda_m1"] = lambda_est;
                estimated["lambda_v1"] = lambda_est * 0.1;
            } else {
                estimated["lambda1"] = lambda_est;
            }
        }
    }

    // 启发式裂缝初始长度推算
    double Lf_est = (L > 0) ? qBound(10.0, L / 20.0, L / 2.0) : 50.0;
    estimated["Lf"] = Lf_est;
    int nf_est = std::max(1, (int)std::round(mp->getNf()));
    for(int i = 1; i <= nf_est; i++) {
        estimated[QString("Lf_%1").arg(i)] = Lf_est;
    }

    return estimated;
}
