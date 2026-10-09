// =============================================================================
//  include/yjk/analysis/ResponseSpectrum.h  ——  设计反应谱与振型分解反应谱法
//
//  GB 50011-2010《建筑抗震设计规范》：
//    · 5.1.4 地震影响系数最大值 αmax 与特征周期 Tg（表 5.1.4-1/2）
//    · 5.1.5 地震影响系数曲线（上升段 / 平台段 / 曲线下降段 / 直线下降段）
//    · 5.2.2 振型分解反应谱法：SRSS（式 5.2.2-2）或 CQC（式 5.2.2-3）
//    · 5.2.2 参与质量不小于总质量 90% 的模态数检查
//
//  【单位】沿用内核 kN-m 体系：
//    α 无量纲；参与质量 m* = γ² 为 t；g = 9.81 m/s²；
//    基底剪力 V = α·m*·g 为 kN。
//
//  【与模态分析的数据流】
//    输入 ModalResult（ModalAnalysis 的产出，含 ω/T/振型/γ/m*），
//    本模块只负责"谱 → 地震作用"，不再碰 K/M 组装；
//    但等效地震力的空间分布需要质量矩阵（MPC 折减后）与方程号空间振型，
//    因此内部经 ModalAnalysis::assembleMass 取 M，从 shapes 抽出主自由度振型。
// =============================================================================
#pragma once

#include <string>
#include <vector>

#include "yjk/analysis/ModalAnalysis.h"
#include "yjk/math/SparseMatrix.h"
#include "yjk/model/Model.h"

namespace yjk {

// -----------------------------------------------------------------------------
//  GB 50011-2010 设计反应谱（地震影响系数曲线，规范图 5.1.5）
// -----------------------------------------------------------------------------
struct DesignSpectrum {
  double alphaMax{0.16};   // 水平地震影响系数最大值（表 5.1.4-1，8 度 0.20g 多遇）
  double tg{0.40};         // 特征周期 s（表 5.1.4-2，二类场地第一组）
  double zeta{0.05};       // 阻尼比（钢筋混凝土取 0.05，钢结构 0.02~0.05）

  // 曲线形状参数（规范 5.1.5）：
  //   γ  —— 衰减指数，0 < ζ < 1
  //   η1 —— 直线下降段斜率调整系数（< 0 时取 0）
  //   η2 —— 阻尼调整系数（< 0.55 时取 0.55）
  double gammaExp() const;   // γ = 0.9 + (0.05 − ζ)/(0.5 + 5ζ)
  double eta1() const;       // η1 = 0.02 + (0.05 − ζ)/(4 + 32ζ)
  double eta2() const;       // η2 = 1 + (0.05 − ζ)/(0.08 + 1.6ζ)

  // 地震影响系数 α(T)。按规范曲线分四段：
  //   0 < T ≤ 0.1        上升段：α = [0.45 + 10(η2 − 0.45)T]·αmax
  //   0.1 < T ≤ Tg       平台段：α = η2·αmax
  //   Tg < T ≤ 5Tg       曲线下降段：α = (Tg/T)^γ · η2·αmax
  //   5Tg < T ≤ 6       直线下降段：α = [η2·0.2^γ − η1(T − 5Tg)]·αmax
  //   T > 6              超曲线终点的长周期按 T = 6 截断（规范曲线定义至 6s）
  double alpha(double T) const;
};

// -----------------------------------------------------------------------------
//  振型分解反应谱法结果
// -----------------------------------------------------------------------------
struct SpectrumResult {
  bool ok{false};
  std::string message;

  int nmodes{0};                     // 参与组合的模态数

  // 各模态：地震影响系数与三个方向的基底剪力（kN）
  std::vector<double> alphaJ;
  std::vector<double> vX, vY, vZ;

  // 组合后的基底剪力（kN）
  double baseX{0.0}, baseY{0.0}, baseZ{0.0};

  // 等效地震力空间分布（【节点序 × 6】索引；方向分量外为 0，kN）——
  // 已按 SRSS/CQC 完成模态组合，可直接作为静力工况施加回结构。
  std::vector<double> forceX, forceY, forceZ;

  // 参与质量比（从 ModalResult 透传）与 ≥90% 检查（规范 5.2.2）
  double massRatioX{0.0}, massRatioY{0.0}, massRatioZ{0.0};
  bool massOkX{false}, massOkY{false}, massOkZ{false};
  bool massOk{false};
};

// -----------------------------------------------------------------------------
//  振型分解反应谱法求解器
// -----------------------------------------------------------------------------
class ResponseSpectrum {
 public:
  explicit ResponseSpectrum(Model& m) : model_(m) {}

  struct Options {
    int method{0};             // 0 = SRSS（式 5.2.2-2），1 = CQC（式 5.2.2-3）
    bool checkMass{true};      // 参与质量 ≥ 90% 检查（不满足仅在结果中标记）
    double verticalScale{0.65};   // 竖向（Z）地震影响系数最大值倍数（规范 5.3.1）
  };

  // 由模态结果计算反应谱地震作用。
  //   mr   —— ModalAnalysis::solve() 的结果（质量归一的振型与 γ/m*）
  //   spec —— 设计反应谱
  // 基底剪力、等效地震力的组合方式（SRSS/CQC）由 opt.method 决定；
  // 竖向地震作用用 spec.alphaMax · opt.verticalScale 计算。
  // 注意：两个重载而非默认参数 —— Options 的默认成员初始化器不能在
  // 封闭类成员函数默认参数中求值（与 SubspaceIteration 相同的原因）。
  SpectrumResult solve(const ModalResult& mr, const DesignSpectrum& spec,
                       const Options& opt);
  SpectrumResult solve(const ModalResult& mr, const DesignSpectrum& spec);

  // ---- 静态工具（供结果输出 / 测试对标）----
  // 模态耦联系数 ρ_jk（规范式 5.2.2-3）：β = ω_k / ω_j
  static double cqcRho(double zeta, double omegaJ, double omegaK);
  // SRSS 组合：√(Σ V_j²)
  static double srss(const std::vector<double>& v);
  // CQC 组合：√(ΣΣ ρ_jk V_j V_k)
  static double cqc(double zeta, const std::vector<double>& omega,
                    const std::vector<double>& v);

 private:
  Model& model_;
};

}  // namespace yjk