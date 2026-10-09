// =============================================================================
//  include/yjk/analysis/LoadCombination.h  ——  荷载组合（标准/基本/地震）
//
//  【为什么"结果级线性组合"是正确且唯一的做法】
//  线弹性静力分析中，荷载 → 位移/内力是线性映射：
//        u = K⁻¹f，  F_section = k·u − f_eq
//  于是组合荷载 Σαᵢ·fᵢ 的解就是各工况解的同一线性组合（逐分量相加）：
//        u_c = Σαᵢ·uᵢ，  F_c = Σαᵢ·Fᵢ
//  无需重新组装、重新求解 —— 对既有工况结果做线性组合即可。
//  这要求各工况结果来自【同一模型、同一拓扑】（自由度编号一致），
//  本类在 combine() 入口做长度一致性校验，不一致直接返回 !ok。
//
//  【规范依据】
//  · GB 50009-2012《建筑结构荷载规范》3.2 承载能力极限状态设计组合：
//      ─ 可变荷载控制：1.2G + 1.4Q（+ Ψc·其他可变）
//      ─ 永久荷载控制：1.35G + 1.4·Ψc·Q
//      标准组合（正常使用）：G + Q
//  · GB 50011-2010《建筑抗震设计规范》5.4.1 地震作用组合：
//      S = γG·S_GE + γEh·S_Ehk + γEv·S_Evk
//      γG = 1.2（重力荷载代表值），γEh = 1.3（水平地震），γEv = 0.5（竖向）
//      重力荷载代表值 S_GE = 恒载 + ΨE·活载，ΨE = 0.5（住宅/办公）
//
//  【组合结果的语义】
//  组合后的 StaticResult 与原工况结果同构（可作 PostProcessor 输入），
//  统计量（maxDisplacement/maxUxy/maxUz/maxRotation/maxReaction 及极值定位）
//  从组合后的向量重新扫描；residual/seconds 在组合下无定义，置 0/取首工况。
// =============================================================================
#pragma once

#include <initializer_list>
#include <map>
#include <string>
#include <vector>

#include "yjk/analysis/StaticAnalysis.h"

namespace yjk {

// -----------------------------------------------------------------------------
//  规范分项系数与组合值系数（GB 50009-2012 / GB 50011-2010）
// -----------------------------------------------------------------------------
namespace combo_coef {
inline constexpr double kPermanentVar   = 1.2;    // 可变荷载控制时恒载分项系数 γG
inline constexpr double kVariable       = 1.4;    // 可变荷载分项系数 γQ
inline constexpr double kPermanentFixed = 1.35;   // 永久荷载控制时恒载分项系数 γG
inline constexpr double kCombFactor     = 0.7;    // 楼面活载组合值系数 Ψc（5.1.1）
inline constexpr double kGamG           = 1.2;    // 地震组合重力代表值分项 γG
inline constexpr double kGamEh          = 1.3;    // 水平地震作用分项 γEh
inline constexpr double kGamEv          = 0.5;    // 竖向地震作用分项 γEv
inline constexpr double kPsiE           = 0.5;    // 重力代表值中活载 ΨE（5.1.3）
}

// -----------------------------------------------------------------------------
//  输入：一个独立求解出的荷载工况
// -----------------------------------------------------------------------------
struct DesignCase {
  std::string name;                  // 工况名（"D" 恒、"L" 活、"Ex" 水平地震…）
  const StaticResult* result{nullptr};
  double factor{1.0};                // 该工况在本组合中的系数 αᵢ

  DesignCase() = default;
  DesignCase(const std::string& n, const StaticResult& r, double f = 1.0)
      : name(n), result(&r), factor(f) {}
};

// -----------------------------------------------------------------------------
//  组合定义：一条组合 = {名字, 工况名 → 系数}
// -----------------------------------------------------------------------------
struct ComboDef {
  std::string name;
  std::map<std::string, double> factors;
};

// -----------------------------------------------------------------------------
//  荷载组合求解器
// -----------------------------------------------------------------------------
class LoadCombination {
 public:
  // ---- 核心 ----
  // 一条组合：c = Σᵢ factorᵢ · resultᵢ（逐分量线性组合）。
  // 要求各工况结果同长度同拓扑；不满足返回 !ok 并附原因。
  static StaticResult combine(const std::string& name,
                              const std::vector<DesignCase>& cases);

  // 多条组合：对 caseResults 按 defs 逐条执行 combine。
  // 组合定义里出现但输入中没有的工况名会被忽略（系数视为 0）。
  static std::vector<StaticResult> generate(
      const std::map<std::string, StaticResult>& caseResults,
      const std::vector<ComboDef>& defs);

  // ---- 包络 ----
  // 对多条组合结果逐分量取极大/极小，得到 envMax / envMin。
  // 统计极值（maxDisplacement 等）为各组合统计值的包络。
  static void envelope(const std::vector<StaticResult>& combos,
                       StaticResult& envMax, StaticResult& envMin);

  // ---- 规范组合定义便捷构造 ----
  // 标准组合（正常使用）：G + Q
  static ComboDef standard(const std::string& g, const std::string& q);
  // 基本组合（承载能力）可变荷载控制：1.2G + 1.4Q（+1.4·Ψc·额外可变 w）
  static ComboDef basicVariable(const std::string& g, const std::string& q,
                                const std::string& w = "");
  // 基本组合（承载能力）永久荷载控制：1.35G + 1.4·Ψc·Q（+1.4·Ψc·w）
  static ComboDef basicPermanent(const std::string& g, const std::string& q,
                                 const std::string& w = "");
  // 地震组合：1.2·(G + 0.5·Q) + 1.3·Eh（+0.5·Ev 若给出且非空）
  static ComboDef seismic(const std::string& g, const std::string& q,
                          const std::string& eh, const std::string& ev = "");

  // ---- 直接构造 ----
  static ComboDef def(const std::string& name,
                      std::initializer_list<std::pair<const char*, double>> fs);
};

}  // namespace yjk