// =============================================================================
//  include/yjk/post/PostProcessor.h  ——  后处理层
//
//  分析层给出的只有两样东西：位移向量 u 和支座反力。
//  工程师要的却是：内力图、应力云图、配筋控制截面、层间位移角、剪重比……
//  这一层负责把"求解器的输出"翻译成"工程师能读的结果"。
//
//  【这一层最容易出错、后果也最直接的三件事】
//
//  ① 杆端内力必须带固端修正（k·u − f_eq）
//     不做修正时，两端固定梁受均布荷载会得到全零内力。
//     已在 StaticAnalysis 里修好，本层直接用 beamForces[e].ef。
//
//  ② 弯矩沿杆长必须【二次】插值，不能线性
//     有分布荷载时 M(ξ) 是抛物线，线性插值在跨中会差 25%。
//     M(ξ) = M_i + L·[V_i·ξ + (V_j−V_i)·ξ²/2]   （由 M' = V 积分）
//
//  ③ 弯矩与剪力的配对不能搞错
//     V_y → 绕 z 弯 → 配 M_z；  V_z → 绕 y 弯 → 配 M_y
//     配错对内力图看着"差不多"，是最隐蔽的一类错误。
//
//  【应力恢复】
//     σ = N/A ± My·hz/Iy ± Mz·hy/Iz
//  需要截面的边缘纤维距离 hy / hz —— 见 SectionProperties 的说明。
// =============================================================================
#pragma once

#include <array>
#include <string>
#include <vector>

#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/model/Model.h"

namespace yjk {
namespace post {

// -----------------------------------------------------------------------------
//  构件类别
// -----------------------------------------------------------------------------
enum class MemberKind { Column, Beam, Brace, Other };
const char* memberKindName(MemberKind k);

// -----------------------------------------------------------------------------
//  梁沿长某个截面的内力与应力
// -----------------------------------------------------------------------------
struct BeamStation {
  double xi{0.0};                       // 0 = i 端，1 = j 端
  double N{0}, Vy{0}, Vz{0}, T{0}, My{0}, Mz{0};
  double vx{0}, vy{0}, vz{0};           // 剪力【全局】分量 kN（楼层剪力要用矢量累加）
  double sigMax{0}, sigMin{0};          // 截面边缘最大/最小正应力 kPa（拉 +，压 −）
  double tau{0};                        // 剪应力（矩形 1.5·V/A 近似）kPa
  double vonMises{0};
};

struct BeamResult {
  Id id{-1};
  Id ni{-1}, nj{-1};
  MemberKind kind{MemberKind::Other};
  double L{0.0};
  std::array<double, 3> dir{0, 0, 1};   // 杆轴方向余弦（全局）

  BeamElement3D::EndForces ends{};      // 端部截面内力（已含固端修正）
  std::vector<BeamStation> stations;

  // 全杆极值
  double sigMax{0}, sigMin{0}, sigAbsMax{0};
  double vonMisesMax{0};
  double util{0.0};                     // 应力比 = |σ|max / 设计强度（0 表示无强度信息）
  double weight{0.0};                   // 自重 kN

  // 控制截面（工程上真正拿去配筋的三个位置）
  double Mmax{0};                       // 沿杆最大弯矩（绝对值）
  double xiMmax{0};
  double Vmax{0};                       // 沿杆最大剪力（绝对值）
  double Nmax{0};                       // 沿杆最大轴力（绝对值，压力为负）
  double Nmin{0};
};

// -----------------------------------------------------------------------------
//  壳单元结果
// -----------------------------------------------------------------------------
struct ShellResult {
  Id id{-1};
  std::array<Id, 4> nodes{-1, -1, -1, -1};
  double area{0.0}, thickness{0.0};
  ShellElement4::MembraneForces mf{};   // 单元中心（未平滑）
  double weight{0.0};                   // 自重 kN
};

// -----------------------------------------------------------------------------
//  节点场（云图用）
//
//  【为什么必须做节点平滑】
//  壳单元的应力是【分片常数】（单元内一个值），直接上色会看到一格一格的台阶。
//  绕节点做面积加权平均后才是连续云图。
//  不做平滑不是"难看"而已 —— 读最大值时会读到台阶的极值，
//  而这个极值通常出现在应力集中处，是离散误差不是真值。
// -----------------------------------------------------------------------------
struct NodeField {
  std::string name;
  std::string unit;
  std::vector<double> v;                // 按节点序
  double vmin{0.0}, vmax{0.0};
  double vminAbs{0.0}, vmaxAbs{0.0};
};

// -----------------------------------------------------------------------------
//  楼层指标
// -----------------------------------------------------------------------------
struct StoryResult {
  int story{0};
  double z{0.0};            // 该层楼面标高
  double h{0.0};            // 层高（本层标高 − 下层标高）
  int nodeCount{0};

  double maxUxy{0.0};       // 该层最大水平位移 m
  double avgUxy{0.0};       // 该层平均水平位移 m
  double maxUz{0.0};        // 该层最大竖向位移 m

  double drift{0.0};        // 【最大】层间位移角 Δu/h（无量纲）
  double driftAvg{0.0};     // 平均层间位移角
  Id driftNode{-1};
  double driftRatio{0.0};   // 位移比 = 最大层间位移 / 平均层间位移

  double shear{0.0};        // 楼层剪力 kN（该层所有柱的剪力之和）
  double weight{0.0};       // 该层以上总重力荷载 kN
  double shearWeightRatio{0.0};   // 剪重比 V/W
};

// -----------------------------------------------------------------------------
//  全模型包络
// -----------------------------------------------------------------------------
struct Envelope {
  double maxDisp{0.0};        Id maxDispNode{-1};
  double maxUxy{0.0};         Id maxUxyNode{-1};
  double maxUz{0.0};          Id maxUzNode{-1};
  double maxRotation{0.0};    Id maxRotationNode{-1};
  double maxVonMises{0.0};    Id maxVonMisesNode{-1};
  double maxReaction{0.0};    Id maxReactionNode{-1};
  double maxReactionZ{0.0};   Id maxReactionZNode{-1};

  double totalWeight{0.0};    // 总重力荷载 kN
  double baseShear{0.0};      // 基底剪力 kN
  double shearWeightRatio{0.0};
  double maxDrift{0.0};       int maxDriftStory{-1};
  double maxShearRatio{0.0};  int maxShearRatioStory{-1};
  double totalLoadZ{0.0};     // 竖向荷载合计（含反力校核用）
};

// -----------------------------------------------------------------------------
//  后处理器
// -----------------------------------------------------------------------------
class PostProcessor {
 public:
  PostProcessor(const Model& m, const StaticResult& r) : model_(m), res_(r) {}

  // stations：每根杆沿长分几段（含两端），奇数保证跨中落在采样点上
  // columnCos：|杆轴的 z 方向余弦| 大于该值判为柱
  void compute(int stations = 21, double columnCos = 0.7);

  const std::vector<BeamResult>& beamResults() const { return beams_; }
  const std::vector<ShellResult>& shellResults() const { return shells_; }
  const std::vector<StoryResult>& stories() const { return stories_; }
  const std::vector<NodeField>& nodeFields() const { return fields_; }
  const Envelope& envelope() const { return env_; }

  const NodeField* field(const std::string& name) const;
  std::vector<std::string> fieldNames() const;

  // 杆端内力最大的若干根（报告用）
  std::vector<const BeamResult*> topBeamsByStress(int n) const;
  std::vector<const BeamResult*> columns() const;

  // 文本报告
  std::string report() const;

  // 一句话摘要（界面状态栏用）
  std::string summary() const;

 private:
  void computeBeams(int nStation, double columnCos);
  void computeShells();
  void computeNodeFields();
  void computeStories(double columnCos);
  void computeEnvelope();

  const Model& model_;
  const StaticResult& res_;

  std::vector<BeamResult> beams_;
  std::vector<ShellResult> shells_;
  std::vector<StoryResult> stories_;
  std::vector<NodeField> fields_;
  Envelope env_;
  std::vector<double> storyLoad_;    // 各层施加的向下荷载 kN（computeStories 填）
};

}  // namespace post
}  // namespace yjk
