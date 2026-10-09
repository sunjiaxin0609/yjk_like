// =============================================================================
//  yjk/math/EigenSolver.h  ——  广义特征值求解器
//
//                    K φ = λ M φ        （λ = ω²，无阻尼自由振动）
//
//  提供两条路线（与线性方程组 Solver.h 的"多路并行"思路一致）：
//   1. DenseGeneralizedEigen  稠密广义雅可比 —— 小规模 / 单元测试 / Ritz 投影
//   2. SubspaceIteration      子空间迭代（Bathe）—— 稀疏 K + 对角集中质量 M，
//                              大规模工程模型的标准路径
//
//  【为什么不引外部特征值库】
//  工程结构分析的模态提取其实只需要"前 p 阶最低特征对"（规范验算用到的
//  也是最低的十几阶或几十阶），这正是子空间迭代的高效区间；
//  而 LANZOS/QR 的优势在"求大量或内部特征值"，结构分析用不到。
//  零外部依赖的前提下，子空间迭代 + 稠密广义雅可比是正确性最可控的组合。
//
//  【为什么质量矩阵用"对角集中质量"而不是一致质量】
//   · 工程惯例（盈建科/PKPM/ETABS 均默认 lumped）：每节点集中平动质量，
//     转动惯量为 0 —— 一致质量需要数值积分，且高阶精度没有工程意义
//   · 对角 M 使"K⁻¹ M 逆迭代"的每次迭代只是一次稀疏前代回代，
//     计算量与二次静力求解相当，代价可控
//   · M 半正定（零质量转动自由度）在算法层面是安全的：M 内积只统计
//     正质量自由度，零质量自由度的振型分量由刚度耦合（静力凝聚）自动确定
//
//  【单位约定】
//  沿用内核 kN-m 体系：K 为 kN/m（kN/rad 等），M 取"吨" t（= kN·s²/m，
//  由单元重量 kN ÷ 9.81 得到），则 λ = ω² 单位 s⁻²，周期 T = 2π/ω (s)。
// =============================================================================
#pragma once

#include <string>
#include <vector>

#include "yjk/math/Solver.h"
#include "yjk/math/SparseMatrix.h"

namespace yjk {

// -----------------------------------------------------------------------------
//  特征对：value = λ = ω²（rad²/s²），vector 为振型（与输入自由度数一致）
// -----------------------------------------------------------------------------
struct EigenPair {
  double value{0.0};
  std::vector<double> vector;
};

// -----------------------------------------------------------------------------
//  1. 稠密广义雅可比：解 K φ = λ M φ（K 对称，M 对称正定）
//
//  算法：Cholesky 转化 + 经典雅可比旋转
//    · M = L Lᵀ（正定 ⇒ 可分解），A = L⁻¹ K L⁻ᵀ 为对称矩阵
//    · 对 A 做雅可比旋转求全部特征对 (μ, z)：K φ = λ M φ ⇔ A z = μ z，
//      其中 φ = L⁻ᵀ z
//  用途：单元测试对标（小系统解析解）、子空间迭代中的 Ritz 投影。
// -----------------------------------------------------------------------------
class DenseGeneralizedEigen {
 public:
  // K、M 为 n×n 稠密矩阵（行主序，只读对称部分即可）。
  // 返回按 λ 升序的全部特征对；M 必须对称正定（半正定用 SubspaceIteration）。
  static bool solve(const std::vector<double>& K, const std::vector<double>& M,
                    int n, std::vector<EigenPair>& out, double tol = 1e-10,
                    int maxSweeps = 100);
};

// -----------------------------------------------------------------------------
//  2. 子空间迭代（Bathe & Wilson）
//
//  思路：用 K⁻¹M 反复作用在"工作子空间"上，把高频成分压掉，
//  再通过 Rayleigh–Ritz 投影求出近似特征对，迭代至收敛。
//
//    X_{k+1} = Q 其中 Q = M-正交化( K⁻¹ (M X_k) )
//    Ritz: K_q = QᵀKQ, M_q = QᵀMQ（小稠密），解稠密广义问题
//    收敛判据：前 nmodes 阶 Ritz 值残差 ||Kφ − λMφ||∞ / ||Kφ||∞ 足够小
//
//  前置条件：K 对称正定（结构已消除刚体自由度），M 对角 ≥ 0。
//  工作子空间维数自动取 q = max(nmodes + 8, 2·nmodes)，封顶到自由度数。
// -----------------------------------------------------------------------------
class SubspaceIteration {
 public:
  struct Options {
    // 显式构造：嵌套 struct 的 default member initializer 不能用于
    // 封闭类成员函数的默认参数（C++ 标准限制），故把默认值放这里。
    Options() = default;
    int maxIter{300};      // 迭代上限（每轮 1 次 K 分解 + q 次回代）
    double tol{1e-9};      // 残差相对容差
    int workDim{0};        // 0 = 自动 q = max(nmodes+8, 2·nmodes)
    unsigned seed{20261008u};  // 初始向量伪随机种子（固定 ⇒ 结果可复现）
  };

  struct Report {
    bool ok{false};
    std::string message;
    int iter{0};           // 实际迭代轮数
    double residual{0.0};  // 最差一阶的最终残差
    int n{0};              // 自由度数
    int q{0};              // 工作子空间维数
    double seconds{0.0};
    Id nnzK{0};
  };

  // 求解前 nmodes 阶最小特征对（λ 升序）。K 稀疏对称正定，M 对角（n 长度，可含 0）。
  // 注意：两个重载而非默认参数 —— Options/Report 是嵌套类，其 default member
  // initializer 不能在封闭类成员函数默认参数中求值（C++ [class.mem] 限制）。
  static bool solve(const SymSparseMatrix& K, const std::vector<double>& M,
                    int nmodes, std::vector<EigenPair>& out);
  static bool solve(const SymSparseMatrix& K, const std::vector<double>& M,
                    int nmodes, std::vector<EigenPair>& out,
                    Options opt, Report* rep);

  // 稀疏对称质量矩阵版本：M 可以是【对角集中质量 + MPC 折叠后的耦合块】。
  // 需要它而不是对角版的原因：刚性楼板（多点约束）把从属节点质量折到主自由度时，
  // 若从属自由度带位置系数（如 ux = UX − θ·dy），TᵀMT 会在主自由度之间产生
  // 非对角耦合项（UX–θ、UY–θ、θ–θ 转动惯量）——这是楼板扭转模态的物理来源，
  // 对角近似会把它丢掉。本实现全部 M 相关运算（内积/正交化/Ritz 投影/残差）
  // 走稀疏矩阵乘，对角 M 只是 nnz = n 的特例。
  static bool solve(const SymSparseMatrix& K, const SymSparseMatrix& M,
                    int nmodes, std::vector<EigenPair>& out);
  static bool solve(const SymSparseMatrix& K, const SymSparseMatrix& M,
                    int nmodes, std::vector<EigenPair>& out,
                    Options opt, Report* rep);

  // Rayleigh 商：直接从向量数值上验算 λ ≈ φᵀKφ / φᵀMφ（诊断辅助，对角 M）
  static double rayleigh(const SymSparseMatrix& K, const std::vector<double>& M,
                         const std::vector<double>& phi);

  // 广义特征值残差范数 ‖Kφ − λMφ‖∞ / ‖Kφ‖∞（对角 M）
  static double residual(const SymSparseMatrix& K, const std::vector<double>& M,
                         const EigenPair& ep);

  // 稀疏 M 版本（结果与对角版一致，供模态/反应谱层与测试使用）
  static double rayleigh(const SymSparseMatrix& K, const SymSparseMatrix& M,
                         const std::vector<double>& phi);
  static double residual(const SymSparseMatrix& K, const SymSparseMatrix& M,
                         const EigenPair& ep);
};

}  // namespace yjk