// =============================================================================
//  include/yjk/analysis/ModalAnalysis.h  ——  模态分析（无阻尼自由振动）
//
//  求解广义特征值问题
//        K φ = ω² M φ
//  其中 K 为结构刚度矩阵（含 MPC 折减，与静力分析完全一致），
//  M 为质量矩阵（对角集中质量 + MPC 折减后的耦合块）。
//
//  【质量矩阵为什么不是纯对角】
//  集中质量假设：每个单元的质量分配到节点上，只作用于三个平动自由度
//  （梁 1/2 到两端，壳 1/4 到四角；转动惯量为零）。
//  但刚性楼板（MPC）把从属节点的平动自由度折到主自由度上之后，
//  Tᵀ M T 会生成主自由度之间的耦合项：
//        ux_i = UX − Θ·dy_i,  uy_i = UY + Θ·dx_i
//  代入质量项 ½m(üx²+üy²) 后出现：
//        M(UX,UX)=Σm, M(UY,UY)=Σm,
//        M(UX,Θ)=−Σm·dy, M(UY,Θ)=+Σm·dx,       ← 质心偏离参考点时的耦合
//        M(Θ,Θ)=Σm(dx²+dy²)                     ← 楼板绕竖轴的转动惯量
//  这些耦合项是"扭转模态"的物理来源，不能省 ——
//  这也是数学层 SubspaceIteration 支持稀疏质量矩阵的原因。
//
//  【数据流的唯一约定 —— 与静力分析共用一个展开表】
//  K、质量折减、位移回填（从属自由度位移由主自由度算出）三处
//  必须使用【同一个】StaticAnalysis::Assembly::expand 表，
//  否则任何一处的 MPC 口径不同都会导致模态频率错而无人察觉。
//  因此本类直接调用 StaticAnalysis::assemble() 拿它产出的 K 与 expand，
//  不再私自组装刚度。
// =============================================================================
#pragma once

#include <string>
#include <vector>

#include "yjk/math/EigenSolver.h"
#include "yjk/math/SparseMatrix.h"
#include "yjk/model/Model.h"
#include "yjk/analysis/StaticAnalysis.h"

namespace yjk {

// -----------------------------------------------------------------------------
//  模态分析结果
// -----------------------------------------------------------------------------
struct ModalResult {
  bool ok{false};
  std::string message;

  Id ndof{0};            // 自由度数（方程数）
  int nmodes{0};         // 实际求解的模态数（≤ 请求值，受有效质量自由度限制）
  int iter{0};
  double residual{0.0};  // 最差模态残差
  double seconds{0.0};
  size_t nnzK{0};
  size_t nnzM{0};

  // 【索引】以下数组长度 = nmodes，按固有频率升序
  std::vector<double> omega;                    // 圆频率 rad/s
  std::vector<double> freq;                     // 频率 Hz
  std::vector<double> period;                   // 周期 s
  std::vector<std::vector<double>> shapes;      // 振型 [m][节点序 × 6]
  // 质量归一（φᵀMφ = 1）。从属自由度位移已按 MPC 展开回填。

  // 模态参与系数 γ（三个地震方向）与模态参与质量 m* = γ²。
  // 定义：γ_j = φ_jᵀ M r_dir，r_dir 为方向单位向量（X/Y/Z 平动=1）。
  // 物理意义：第 j 阶模态在给定方向地震下"贡献多少有效质量"。
  std::vector<double> gammaX, gammaY, gammaZ;
  std::vector<double> mstarX, mstarY, mstarZ;

  // 方向总质量（rᵀMr）与参与质量累计、占比（≥90% 是规范要求）
  double massX{0.0}, massY{0.0}, massZ{0.0};
  double totalMassX{0.0}, totalMassY{0.0}, totalMassZ{0.0};
  double ratioX{0.0}, ratioY{0.0}, ratioZ{0.0};
};

// -----------------------------------------------------------------------------
//  模态求解器
// -----------------------------------------------------------------------------
class ModalAnalysis {
 public:
  explicit ModalAnalysis(Model& m) : model_(m) {}

  struct Options {
    int nmodes{0};       // 0 = 自动：(层数×3) 与 9 取大，再受有效质量自由度钳制
    double tol{1e-9};    // 子空间迭代残差容差
    int maxIter{300};
    int workDim{0};      // 0 = 自动（min(nmodes+8, 2·nmodes, 有效质量自由度)）
  };

  // 组装质量矩阵（供测试单独检查，与折减后的 K 同一套 expand 表）
  //
  // 单位：质量用 t（吨，= kN·s²/m）。单元自重返回的是 kN（重量），
  // 需 ÷ g=9.81 m/s² 才是质量 —— 这一步最容易漏。
  struct Assembly {
    SymSparseMatrix M;   // ndof × ndof，含 MPC 折减耦合项
    Id ndof{0};
    size_t nnz{0};
    double totalMass{0.0};    // Σ 单元质量 (t)，含被约束掉的部分（仅统计用）
  };
  bool assembleMass(Assembly& out, std::string* why = nullptr);

  // 求解前 nmodes 阶模态（或根据 Options 自动确定阶数）。
  // 注意：两个重载而非默认参数 —— Options 的默认成员初始化器不能在
  // 封闭类成员函数默认参数中求值（与 SubspaceIteration 相同的原因）。
  ModalResult solve();
  ModalResult solve(const Options& opt);

  // ---- 静态工具（供反应谱 / 结果输出复用）----
  // 方向单位向量（方程号索引）：comp=0/1/2 → X/Y/Z 平动自由度置 1。
  // 从属自由度质量已折入主自由度，方向向量只需在【自由】方程号上置 1。
  static std::vector<double> dirVector(const Model& m, int comp);

  // 模态参与系数 γ = φᵀ M r（φ 质量归一）
  static double participation(const SymSparseMatrix& M,
                              const std::vector<double>& phi,
                              const std::vector<double>& r);

 private:
  Model& model_;
};

}  // namespace yjk