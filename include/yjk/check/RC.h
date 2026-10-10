// =============================================================================
//  include/yjk/check/RC.h —— 钢筋混凝土构件正截面承载力验算模块（P4-T1）
//
//  依据：GB 50010-2010《混凝土结构设计规范》第 6 章
//    · §6.2.6  受压区混凝土等效矩形应力图形（α1/β1/εcu）
//    · §6.2.7  界限相对受压区高度 ξb
//    · §6.2.10 正截面受弯承载力（单筋 / 双筋矩形截面）
//    · §6.2.14 受压区高度 x < 2as' 时的受弯承载力
//    · §6.2.15 轴心受压构件正截面承载力（稳定系数 φ）
//    · §6.2.17 矩形截面偏心受压构件正截面承载力（大/小偏心判别）
//    · §6.4.3 / 6.4.4 受扭塑性抵抗矩 Wt 与纯扭承载力
//    · §6.4.8  剪扭构件承载力（混凝土受扭承载力降低系数 βt）
//    · §8.5.1  受弯构件纵向受拉钢筋最小配筋率
//
//  单位约定（与 include/yjk/model/Section.h 一致）：
//    长度 m、力 kN、应力 kPa、弯矩 kN·m、钢筋面积 m²。
//  依赖：仅 C++17 标准库（无 Qt / 外部数值库）。
//
//  设计取向：纯函数输入输出，与后处理层解耦——
//    调用方（后处理 / 测试）把内力 M、N、T 与截面、材料、配筋参数
//    以标量传入，拿到结构化的验算结果；便于用规范算例与手算对拍。
// =============================================================================
#ifndef YJK_CHECK_RC_H
#define YJK_CHECK_RC_H

namespace yjk {
namespace check {

// -----------------------------------------------------------------------------
//  混凝土受压区应力图形参数（GB 50010 §6.2.6 表 6.2.6 / §6.2.1)
// -----------------------------------------------------------------------------
struct StressBlockParams {
  double alpha1 = 1.0;  // 矩形应力图受压区高度换算系数（C50 及以下 = 1.0）
  double beta1 = 0.8;   // 等效矩形应力图高度系数（C50 及以下 = 0.8）
  double epsCu = 0.0033; // 混凝土极限压应变（C50 = 0.0033）
};

// 按混凝土抗压强度设计值 fc (kPa) 取应力图形参数。
// 覆盖 C30~C80（C50 及以下按规范直接取 1.0/0.8/0.0033，
// 更高等级在规范表值间线性内插）。
StressBlockParams stressBlockFromFc(double fc);

// 界限相对受压区高度 ξb（GB 50010 §6.2.7-1）
//   ξb = β1 / (1 + fy / (Es·εcu))
// fy、Es 单位 kPa。
double relativeBoundaryHeight(double fy, double Es, const StressBlockParams& sb);

// -----------------------------------------------------------------------------
//  受弯：矩形截面正截面受弯承载力（GB 50010 §6.2.10 / §6.2.14）
// -----------------------------------------------------------------------------
struct FlexuralInput {
  double b;      // 截面宽度 (m)
  double h;      // 截面高度 (m)
  double fc;     // 混凝土轴心抗压强度设计值 (kPa)
  double ft;     // 混凝土轴心抗拉强度设计值 (kPa)
  double fy;     // 纵向受拉钢筋抗拉强度设计值 (kPa)
  double Es;     // 钢筋弹性模量 (kPa)
  double as;     // 受拉钢筋合力点至受拉边缘距离 (m)
  double as_p;   // 受压钢筋合力点至受压边缘距离 (m)
};

enum class FlexuralStatus {
  Ok,             // 配筋在强度与最小配筋率约束内
  OverReinforced, // 超筋：相对受压区高度超过 ξb，单筋截面不足
  BelowMinReinf   // 少筋：所需 As 小于最小配筋率，按 AsMin 控制
};

struct FlexuralResult {
  double x = 0.0;     // 混凝土受压区高度 (m)
  double xi = 0.0;    // 相对受压区高度 x/h0
  double xiB = 0.0;   // 界限相对受压区高度
  double As = 0.0;    // 受拉钢筋计算面积 (m²)
  double AsMin = 0.0; // 最小配筋面积 (m²)
  FlexuralStatus status = FlexuralStatus::Ok;
};

// 单筋矩形截面：已知弯矩 M (kN·m)，求受拉钢筋 As。
//   As = α1·fc·b·x / fy，其中 x 由 αs = M/(α1·fc·b·h0²) 解二次式。
//   ξ > ξb 时置 OverReinforced（返回按 ξb 对应的配筋，须改双筋或加大截面）。
FlexuralResult flexuralSingle(const FlexuralInput& in, double M);

// 双筋矩形截面：已知受压钢筋 As'（m²）与弯矩 M，求受拉钢筋 As。
//   M1 = M - fy'·As'·(h0 - as')，由 M1 解 x（§6.2.10-3/4）；
//   x < 2as' 时按 §6.2.14 取 As = M / (fy·(h0 - as'))。
//   输入 As' 仍不足以平衡（ξ > ξb）时置 OverReinforced。
FlexuralResult flexuralDouble(const FlexuralInput& in, double M, double As_p);

// 受弯构件纵向受拉钢筋最小配筋率（GB 50010 §8.5.1 表 8.5.1）
//   ρmin = max(0.2%, 0.45·ft/fy)
double minReinfRatio(double ft, double fy);

// -----------------------------------------------------------------------------
//  受压：轴心受压（§6.2.15）与矩形截面偏心受压（§6.2.17）
// -----------------------------------------------------------------------------
struct ColumnInput {
  double b;      // 截面宽度 (m)
  double h;      // 截面高度 (m)
  double fc;     // 混凝土轴心抗压强度设计值 (kPa)
  double ft;     // 混凝土轴心抗拉强度设计值 (kPa)
  double fy;     // 纵向钢筋抗压（拉）强度设计值 (kPa)
  double Es;     // 钢筋弹性模量 (kPa)
  double as;     // 受拉侧钢筋合力点至近边距离 (m)
  double as_p;   // 受压侧钢筋合力点至近边距离 (m)
};

// 轴心受压构件正截面承载力（GB 50010 §6.2.15）
//   Nu = 0.9·φ·(fc·A + fy'·As')
// l0 为构件计算长度 (m)；φ 按 l0/b 查表 6.2.15（其间线性内插）。
double axialCapacity(const ColumnInput& in, double l0, double As_p);

// 矩形截面偏心受压：给定轴力 N (kN) 与弯矩 M (kN·m)，核算正截面承载力。
//   返回该配筋（As、As'）可承受的界限轴力 Nu；Nu >= N 即满足。
//   大偏心 (ξ <= ξb)：由轴向力平衡 x = N/(α1·fc·b)，
//     Nu = α1·fc·b·x + fy'·As' - fy·As；
//   小偏心 (ξ > ξb)：ξ 按规范式迭代求解，Nu 相应公式。
//   e = ei + h/2 - as，ei = e0 + ea，e0 = M/N，
//   ea = max(20mm, h/30)（§6.2.5 附加偏心距）。
struct ColumnResult {
  double Nu = 0.0;     // 正截面受压承载力（与 N 同号比较）(kN)
  double x = 0.0;      // 受压区高度 (m)
  double xi = 0.0;     // 相对受压区高度
  double xiB = 0.0;    // 界限相对受压区高度
  bool largeEccentric = true; // true=大偏心受压, false=小偏心受压
  bool satisfied = false;     // Nu >= N
};
ColumnResult columnCapacity(const ColumnInput& in, double N, double M,
                            double As, double As_p);

// -----------------------------------------------------------------------------
//  受扭：受扭塑性抵抗矩、纯扭承载力（§6.4.3/6.4.4）与剪扭承载力（§6.4.8）
// -----------------------------------------------------------------------------
struct TorsionInput {
  double b;      // 截面宽度 (m)
  double h;      // 截面高度 (m)
  double h0;     // 截面有效高度 (m)（受剪与剪扭用）
  double fc;     // 混凝土轴心抗压强度设计值 (kPa)
  double ft;     // 混凝土轴心抗拉强度设计值 (kPa)
  double fyv;    // 箍筋抗拉强度设计值 (kPa)
  double fy;     // 受扭纵筋抗拉强度设计值 (kPa)
  double c;      // 纵筋形心至截面边缘的混凝土保护层厚度 (m)（求 bcor/hcor）
};

// 矩形截面受扭塑性抵抗矩 Wt = b²·(3h - b)/6（GB 50010 §6.4.3-1）(m³)
double torsionPlasticModulus(double b, double h);

// 纯扭构件受扭承载力（GB 50010 §6.4.4-1）
//   Tu = 0.35·ft·Wt + 1.2·√ζ·fyv·Ast1·Acor/s
//   ζ = fy·Astl·s / (fyv·Ast1·ucor)   （§6.4.4-3，取 0.6<=ζ<=1.7）
// Ast1 单肢抗扭箍筋面积 (m²)、Astl 受扭纵筋总面积 (m²)、s 箍筋间距 (m)。
// 返回纯扭承载力 Tu (kN·m) 与所用 ζ。
struct TorsionResult {
  double Wt = 0.0;   // 受扭塑性抵抗矩 (m³)
  double Tu = 0.0;   // 受扭承载力 (kN·m)
  double zeta = 1.0; // 受扭纵筋与箍筋配筋强度比 ζ
};
TorsionResult pureTorsion(const TorsionInput& in, double Ast1, double Astl,
                          double s);

// 剪扭构件（GB 50010 §6.4.8）：
//   混凝土受扭承载力降低系数 βt = 1.5 / (1 + 0.5·V·Wt/(T·b·h0))，
//   0.5 <= βt <= 1.0（§6.4.8-2）。
//   受扭部分 Tu = 0.35·βt·ft·Wt + 1.2·√ζ·fyv·Ast1·Acor/s（§6.4.8-3）
//   受剪部分 Vu = (1.5 - βt)·0.7·ft·b·h0 + fyv·Asv·h0/s（§6.4.8-4）
// Asv 为受剪箍筋（同一截面各肢）总截面面积 (m²)；其余同 pureTorsion。
struct ShearTorsionResult {
  double Wt = 0.0;
  double zeta = 1.0;  // 受扭纵筋与箍筋配筋强度比 ζ
  double betaT = 1.0; // 混凝土受扭承载力降低系数 βt
  double Tu = 0.0;    // 剪扭共同作用下受扭承载力 (kN·m)
  double Vu = 0.0;    // 剪扭共同作用下受剪承载力 (kN)
};
ShearTorsionResult shearTorsion(const TorsionInput& in, double V, double T,
                                double Ast1, double Astl, double Asv, double s);

}  // namespace check
}  // namespace yjk

#endif  // YJK_CHECK_RC_H