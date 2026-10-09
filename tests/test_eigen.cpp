// =============================================================================
//  tests/test_eigen.cpp  ——  广义特征值求解器验证
//
//  P0 模态分析的地基：K φ = λ M φ（λ = ω²）。本文件直击求解器本身，
//  用三组独立证据锁定正确性：
//   ① 解析解对标 —— 弹簧-质量链有闭式特征值/振型（固定-自由链）
//   ② 稠密法与稀疏子空间迭代互验 —— 两种算法互相独立
//   ③ 不变量校验 —— φᵀMφ = δ、φᵀKφ = diag(λ)（质量/刚度双正交）
//      以及"零质量自由度"的静力凝聚行为（无质量 dof 上残差必为 0）
//
//  覆盖要点：
//   · 2 自由度系统（特征方程可手解）：λ = 3 ± √5
//   · 固定-自由链 n=3 / n=30：λ_k = 2 − 2cos((2k−1)π/(2n+1))，k=1..n
//   · 零质量 dof：静力凝聚后无质量分量残差严格为 0
//   · 收敛报告 Report（迭代轮数、残差、子空间维数）
// =============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "yjk/math/EigenSolver.h"

using namespace yjk;

static int g_fail = 0;
static int g_pass = 0;

static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}

static void checkNear(double got, double want, double relTol, const char* what) {
  const double err = std::abs(got - want);
  const double rel = std::abs(want) > 1e-30 ? err / std::abs(want) : err;
  if (rel <= relTol) {
    ++g_pass;
    std::printf("  [ OK ] %s  (got %.12g, want %.12g, rel %.2e)\n", what, got, want, rel);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %.12g, want %.12g, rel %.2e)\n", what, got, want, rel);
  }
}

static void checkInt(long long got, long long want, const char* what) {
  if (got == want) { ++g_pass; std::printf("  [ OK ] %s  (got %lld)\n", what, got); }
  else { ++g_fail; std::printf("  [FAIL] %s  (got %lld, want %lld)\n", what, got, want); }
}

// -----------------------------------------------------------------------------
//  构造"固定-自由"弹簧质量链：墙 —— k —— m1 —— k —— m2 —— ... —— k —— mn
//  闭式解（标准教科书）：
//     λ_k = 2k₀/m · (1 − cos((2k−1)π/(2n+1)))，k = 1..n
//  本文件统一取 k₀ = m = 1 ⇒ λ_k = 2 − 2cos((2k−1)π/(2n+1))
// -----------------------------------------------------------------------------
namespace {

// 链的稀疏刚度矩阵 K（对称）与非零质量对角向量 M
// 固定-自由链（左端墙，右端自由，共 n 个质量、n 段弹簧 k=1, m=1）：
//   对角 = [2, 2, ..., 2, 1]（末质量只连一条弹簧），
//   副对角（上三角非对角）= -1（每对相邻质量一条弹簧）。
void chainSystem(int n, SymSparseMatrix& K, std::vector<double>& M) {
  TripletAssembler ta;
  for (int i = 0; i < n; ++i) {
    ta.add(i, i, (i == n - 1) ? 1.0 : 2.0);
    if (i >= 1) ta.add(i, i - 1, -1.0);
  }
  K.buildFrom(ta.triplets(), n);
  M.assign(static_cast<size_t>(n), 1.0);
}

double analyticLambda(int n, int k) {  // k: 1-based mode index
  const double pi = std::acos(-1.0);
  return 2.0 - 2.0 * std::cos((2.0 * k - 1.0) * pi / (2.0 * n + 1.0));
}

// 把稀疏 K 转成稠密行主序（DenseGeneralizedEigen 用），M 对角转稠密
void toDense(const SymSparseMatrix& K, const std::vector<double>& M,
             std::vector<double>& Kd, std::vector<double>& Md) {
  const int n = static_cast<int>(K.size());
  Kd.assign(static_cast<size_t>(n) * n, 0.0);
  Md.assign(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) Kd[static_cast<size_t>(i) * n + j] = K.at(i, j);
    Md[static_cast<size_t>(i) * n + i] = M[static_cast<size_t>(i)];
  }
}

// φᵀXφ 双线性型（X 稠密）
double bilinear(const std::vector<double>& X, const std::vector<double>& phi, int n) {
  double s = 0.0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) s += phi[static_cast<size_t>(i)] * X[static_cast<size_t>(i) * n + j] * phi[static_cast<size_t>(j)];
  return s;
}

}  // namespace

// -----------------------------------------------------------------------------
//  1. 两自由度系统：m1=m2=1, k1=k2=2（串联，两端固定）
//     K = [[4,-2],[-2,2]]，M = I
//     特征方程 (4−λ)(2−λ) − 4 = 0 ⇒ λ = 3 ± √5
// -----------------------------------------------------------------------------
static void testTwoDof() {
  std::printf("\n== 二自由度弹簧系统：λ = 3 ± √5 ==\n");
  const int n = 2;
  std::vector<double> Kd{4.0, -2.0, -2.0, 2.0}, Md{1.0, 0.0, 0.0, 1.0};

  // 稠密法
  std::vector<EigenPair> dense;
  check(DenseGeneralizedEigen::solve(Kd, Md, n, dense), "稠密广义雅可比求解成功");
  if (dense.empty()) { std::printf("  稠密法失败，中止该组\n"); return; }
  checkInt(static_cast<long long>(dense.size()), n, "返回全部特征对");
  checkNear(dense[0].value, 3.0 - std::sqrt(5.0), 1e-9, "λ1");
  checkNear(dense[1].value, 3.0 + std::sqrt(5.0), 1e-9, "λ2");

  // 稀疏子空间迭代
  SymSparseMatrix Ks;
  std::vector<double> Ms = {1.0, 1.0};
  {
    TripletAssembler ta;
    ta.add(0, 0, 4.0); ta.add(1, 1, 2.0); ta.add(1, 0, -2.0);
    Ks.buildFrom(ta.triplets(), n);
  }
  SubspaceIteration::Report rep;
  std::vector<EigenPair> sub;
  if (!SubspaceIteration::solve(Ks, Ms, n, sub, {}, &rep)) {
    std::printf("  FAIL 子空间迭代求解: %s\n", rep.message.c_str());
    return;
  }
  checkInt(rep.iter <= 5 ? 1 : 0, 1, "小系统收敛快（≤5 轮）");
  checkNear(sub[0].value, dense[0].value, 1e-9, "子空间 λ1 == 稠密 λ1");
  checkNear(sub[1].value, dense[1].value, 1e-9, "子空间 λ2 == 稠密 λ2");

  // 不变量：φᵀMφ = 1、φᵀKφ = λ、φᵀKφ/(φᵀMφ) = λ（Rayleigh 商）
  const std::vector<double> Kdv = Kd, Mdv = Md;
  for (int m = 0; m < n; ++m) {
    checkNear(bilinear(Mdv, sub[static_cast<size_t>(m)].vector, n), 1.0, 1e-9,
              "质量归一 φᵀMφ = 1");
    const double kqq = bilinear(Kdv, sub[static_cast<size_t>(m)].vector, n);
    checkNear(kqq / bilinear(Mdv, sub[static_cast<size_t>(m)].vector, n),
              sub[static_cast<size_t>(m)].value, 1e-9, "Rayleigh 商 = λ");
  }
  // φ₁ᵀMφ₂ = 0（M-正交）
  double cross = 0.0;
  for (int i = 0; i < n; ++i)
    cross += sub[0].vector[static_cast<size_t>(i)] * sub[1].vector[static_cast<size_t>(i)];
  checkNear(cross, 0.0, 1e-9, "振型 M-正交 φ₁ᵀMφ₂ = 0");
}

// -----------------------------------------------------------------------------
//  2. 固定-自由链 n=3：λ_k = 2 − 2cos((2k−1)π/7)
// -----------------------------------------------------------------------------
static void testChainThree() {
  std::printf("\n== 三自由度固定-自由链：λ_k = 2 − 2cos((2k−1)π/7) ==\n");
  const int n = 3;
  SymSparseMatrix K;
  std::vector<double> M;
  chainSystem(n, K, M);

  std::vector<double> Kd, Md;
  toDense(K, M, Kd, Md);

  std::vector<EigenPair> dense;
  check(DenseGeneralizedEigen::solve(Kd, Md, n, dense), "稠密法求解成功");
  if (dense.empty()) { std::printf("  稠密法失败，中止该组\n"); return; }
  for (int k = 1; k <= n; ++k)
    checkNear(dense[static_cast<size_t>(k - 1)].value, analyticLambda(n, k), 1e-9,
              "稠密 λ_k 对标解析解");

  SubspaceIteration::Report rep;
  std::vector<EigenPair> sub;
  if (!SubspaceIteration::solve(K, M, n, sub, {}, &rep)) {
    std::printf("  FAIL 子空间迭代求解: %s\n", rep.message.c_str());
    return;
  }
  checkInt(rep.q, n, "子空间维数 = 自由度（小系统全空间）");
  for (int m = 0; m < n; ++m)
    checkNear(sub[static_cast<size_t>(m)].value, dense[static_cast<size_t>(m)].value, 1e-9,
              "子空间 λ 与稠密一致");

  // 残差（子空间方法自带判据）应远小于 1e-7
  for (int m = 0; m < n; ++m)
    check(SubspaceIteration::residual(K, M, sub[static_cast<size_t>(m)]) < 1e-9,
          "残差 ||Kφ−λMφ||∞/||Kφ||∞ < 1e-9");
}

// -----------------------------------------------------------------------------
//  3. 大链 n=30，取前 5 阶：稠密 vs 子空间 vs 解析（子空间需实际迭代）
// -----------------------------------------------------------------------------
static void testChainLarge() {
  std::printf("\n== 三十自由度链，前 5 阶 ==\n");
  const int n = 30, nmodes = 5;
  SymSparseMatrix K;
  std::vector<double> M;
  chainSystem(n, K, M);

  std::vector<double> Kd, Md;
  toDense(K, M, Kd, Md);
  std::vector<EigenPair> dense;
  check(DenseGeneralizedEigen::solve(Kd, Md, n, dense), "稠密法全谱求解成功");
  if (dense.empty()) { std::printf("  稠密法失败，中止该组\n"); return; }
  for (int k = 1; k <= nmodes; ++k)
    checkNear(dense[static_cast<size_t>(k - 1)].value, analyticLambda(n, k), 1e-7,
              "稠密前 5 阶对标解析解");

  SubspaceIteration::Report rep;
  std::vector<EigenPair> sub;
  if (!SubspaceIteration::solve(K, M, nmodes, sub, {}, &rep)) {
    std::printf("  FAIL 子空间迭代（30 自由度）: %s\n", rep.message.c_str());
    return;
  }
  checkInt(rep.q, nmodes + 8, "自动子空间维数 q = nmodes+8 = 13");
  check(rep.residual < 1e-9, "报告残差 < 1e-9");
  checkInt(rep.iter >= 2 ? 1 : 0, 1, "大系统确实发生迭代（≥2 轮）");
  for (int m = 0; m < nmodes; ++m)
    checkNear(sub[static_cast<size_t>(m)].value, dense[static_cast<size_t>(m)].value, 1e-7,
              "子空间前 5 阶与稠密一致");
}

// -----------------------------------------------------------------------------
//  4. 零质量自由度：M = diag(1,1,0)，n=3 固定-自由链且第三个 dof 无质量
//     无质量 dof 上残差必须严格为 0（静力凝聚行为）
// -----------------------------------------------------------------------------
static void testZeroMassDof() {
  std::printf("\n== 零质量自由度（静力凝聚行为）==\n");
  SymSparseMatrix K;
  std::vector<double> M;
  chainSystem(3, K, M);
  M[2] = 0.0;                      // 第三个质量 = 0

  SubspaceIteration::Report rep;
  std::vector<EigenPair> sub;
  if (!SubspaceIteration::solve(K, M, 2, sub, {}, &rep)) {
    std::printf("  FAIL 零质量 dof 迭代: %s\n", rep.message.c_str());
    return;
  }

  // 前两阶 λ 应与"两质量系统 + 刚度凝聚"一致：
  // 无质量 dof 静力凝聚后 K* = [[2,-1],[-1,1]] ⇒ λ = (3±√5)/2
  checkNear(sub[0].value, (3.0 - std::sqrt(5.0)) / 2.0, 1e-7, "凝聚后 λ1 = (3−√5)/2");
  checkNear(sub[1].value, (3.0 + std::sqrt(5.0)) / 2.0, 1e-7, "凝聚后 λ2 = (3+√5)/2");

  // 无质量自由度上的方程残差严格为零：行 2（dof index 2）的 Kφ − λMφ = 0
  for (const auto& ep : sub) {
    std::vector<double> Kx;
    K.multiply(ep.vector, Kx);
    // 第 3 行（i=2）：Kx[2] − λ·0·φ[2] 必须 ≈ 0（刚度行平衡决定该分量）
    checkNear(Kx[2], 0.0, 1e-9, "零质量行残差 = 0（静力凝聚）");
  }
}

// -----------------------------------------------------------------------------
//  5. 频率与周期换算：ω = √λ，T = 2π/ω（单位体系 sanity）
// -----------------------------------------------------------------------------
static void testFreqPeriod() {
  std::printf("\n== 频率/周期换算 ==\n");
  // 单自由度：m = 1 t，k = 1 kN/m ⇒ ω = 1 rad/s，T = 2π s ≈ 6.2832
  SymSparseMatrix K;
  std::vector<double> M{1.0};
  {
    TripletAssembler ta;
    ta.add(0, 0, 1.0);
    K.buildFrom(ta.triplets(), 1);
  }
  std::vector<EigenPair> sub;
  SubspaceIteration::Report rep;
  if (!SubspaceIteration::solve(K, M, 1, sub, {}, &rep)) {
    std::printf("  FAIL 单自由度求解: %s\n", rep.message.c_str());
    return;
  }
  checkNear(std::sqrt(sub[0].value), 1.0, 1e-9, "ω = √λ = 1 rad/s");
  const double pi = std::acos(-1.0);
  checkNear(2.0 * pi / std::sqrt(sub[0].value), 2.0 * pi, 1e-9, "T = 2π s");
}

int main() {
  std::printf("=== test_eigen: 广义特征值求解器 ===\n");
  testTwoDof();
  testChainThree();
  testChainLarge();
  testZeroMassDof();
  testFreqPeriod();

  std::printf("\n通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}