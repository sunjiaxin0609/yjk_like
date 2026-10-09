// =============================================================================
//  include/yjk/analysis/StaticAnalysis.h  ——  线性静力分析
//
//  流程（每一步都必须能被单独检查，否则出错时无法定位）：
//    ① 自由度编号（模型层已完成）
//    ② 组装全局刚度 —— 单元局部 → 全局
//    ③ 约束处理 —— 自由度消去（【不是】罚项）
//    ④ 组装荷载向量 —— 节点荷载 + 单元等效节点荷载
//    ⑤ 求解
//    ⑥ 回填约束自由度（给零或给沉降值）
//    ⑦ 内力恢复
//
//  【为什么用自由度消去而不是罚项】
//  罚项法（在对角线加 1e20）会带来两个问题：
//    · 条件数被污染，稀疏 LDL^T 在病态矩阵上精度损失严重
//    · 反力无法准确恢复（因为罚项本身也"出力"）
//  自由度消去后，约束自由度直接从方程里删掉，反力由
//  "被消去自由度对应的原方程"精确恢复。
//
//  【多点约束（MPC）走的是同一套思路】
//  刚性楼板把整层节点绑成刚体，从属自由度不是"被固定"，
//  而是"等于若干主自由度的线性组合"：
//        u_slave = Σ c_m · u_master
//  组装时把它按系数折到主自由度上（等价于 K_r = Tᵀ K T、F_r = Tᵀ F），
//  求解后【再折回来】得到从属位移。全程不引入任何罚项，
//  条件数与原问题一致，约束力仍可从"被消去那一行的原方程"精确恢复。
// =============================================================================
#pragma once

#include <map>
#include <string>
#include <vector>

#include "yjk/math/Solver.h"
#include "yjk/model/Model.h"

namespace yjk {

// -----------------------------------------------------------------------------
//  分析结果
// -----------------------------------------------------------------------------
struct StaticResult {
  bool ok{false};
  std::string message;

  // 【索引约定】u 与 reaction 都按【节点序 × 6】：节点 n 的第 k 个分量在 [n*6+k]。
  // 不用模型内部的重编号号 —— 理由见 Model.h 中 Node::baseDof 的说明。
  std::vector<double> u;          // 位移（m / rad），约束自由度为 0 或给定的支座沉降
  std::vector<double> reaction;   // 支座反力（kN / kN·m），仅约束自由度非零

  // 单元内力（按单元索引）
  struct BeamForce {
    BeamElement3D::EndForces ef;
    double xi{0.0};               // 取值位置（0 = i 端，1 = j 端）
  };
  std::vector<BeamForce> beamForces;

  struct ShellForce {
    ShellElement4::MembraneForces mf;
    double area{0.0};
    Vec3 normal{0, 0, 1};
  };
  std::vector<ShellForce> shellForces;

  // 统计
  Id ndof{0};          // 自由度数（求解的方程数）
  Id nconstrained{0};
  size_t nnz{0};
  size_t nnzL{0};
  double seconds{0.0};
  double maxDisplacement{0.0};
  double maxUxy{0.0};
  double maxUz{0.0};
  double maxRotation{0.0};
  double maxReaction{0.0};
  double residual{0.0};      // ‖Ku−f‖∞ / ‖f‖∞

  // 极值定位（后处理高亮用）
  Id maxDispNode{-1};
  Id maxDispElem{-1};
};

// -----------------------------------------------------------------------------
//  线性静力求解器
// -----------------------------------------------------------------------------
class StaticAnalysis {
 public:
  // 约束自由度的处理方式
  enum class Support {
    Zero,     // 固定自由度位移 = 0
    Given,    // 用节点上给定的沉降值（supportDisplacement）
  };

  explicit StaticAnalysis(Model& m) : model_(m) {}

  // 支座位移（用于地基沉降、预位移）。索引 = 节点自由度号。
  void setSupportDisplacement(Id globalDof, double value) {
    supportDisp_[globalDof] = value;
    supportMode_ = Support::Given;
  }

  // 求解
  StaticResult solve();

  // 组装（供测试单独检查：刚度矩阵本身是否正确）
  struct Assembly {
    SymSparseMatrix K;
    std::vector<double> f;      // 按【节点序 × 6】索引（含约束自由度上的荷载）
    Id ndof{0};
    size_t triplets{0};

    // ---- 自由度展开表（MPC 的核心）----
    //
    //   expand[gdof] = {(方程号, 系数), ...}
    //      · 自由自由度 → {(自己的方程号, 1)}
    //      · 固定自由度 → 空
    //      · 从属自由度 → 主自由度的展开
    //
    //  【为什么要缓存在这里，而不是各自现算】
    //  组装（刚度）、荷载压缩、位移回填三处都要用同一个展开表。
    //  三处各建一份的话，只要有一处口径不同（比如忘了跳过固定自由度），
    //  得到的方程是自洽的、残差也正常，唯独位移是错的 —— 最难查的一类。
    struct Expand { Id d[4]; double c[4]; int n{0}; };
    std::vector<Expand> expand;
  };
  bool assemble(Assembly& out, std::string* why = nullptr);

  // 自由度消去：返回 (新号 → 原号) 映射
  static std::vector<Id> buildDofMap(const Model& m, Id nTotalDof);

 private:
  // 把单元局部自由度向量散装到全局向量（按 dof 映射）
  static void scatter(const Element& el, const Model& m,
                      const std::vector<double>& local,
                      std::vector<double>& global);

  Model& model_;
  std::map<Id, double> supportDisp_;
  Support supportMode_{Support::Zero};
};

}  // namespace yjk
