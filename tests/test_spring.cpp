// =============================================================================
//  tests/test_spring.cpp  ——  弹簧单元（SpringElement3D）验证
//
//  验证策略：对标解析解，覆盖——
//   1. 单元级：局部刚度对角形式（k、-k）、k=0 释放标记、全局变换后
//      平动刚体零内力 + 旋转投影刚度；
//   2. 整体级：轴向弹簧支座 u = P/k、内力 -P（拉为正）；
//   3. 串联弹簧：等效刚度 ke = k1·k2/(k1+k2)、中间节点位移 = P/k2、
//      两弹簧内力均 -P；
//   4. 弹簧并联支座：u = P/(k1+k2)、各弹簧内力按刚度分配；
//   5. 45° 倾斜弹簧：全局 X 向等效刚度 = k·cos²45°（Tᵀ Kc T 非轴对齐验证）；
//   6. k=0 释放语义：无荷载自动锁定（dof=-1、不奇异）；带节点荷载被 check()
//      拦截（防荷载被"零刚度自由度"吞掉）。
//
//  弹簧取"拉为正"符号约定：力 +X 推向固定端 → 弹簧受压 → N = -P。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <array>

#include "yjk/element/SpringElement3D.h"
#include "yjk/model/Model.h"
#include "yjk/analysis/StaticAnalysis.h"

using namespace yjk;

static int g_pass = 0, g_fail = 0;
static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}
static void checkNear(double got, double want, double relTol, const char* what) {
  const double rel = std::abs(want) > 1e-30 ? std::abs(got - want) / std::abs(want)
                                            : std::abs(got - want);
  if (rel <= relTol) {
    ++g_pass;
    std::printf("  [ OK ] %s  (got %.8g, want %.8g, rel %.2e)\n", what, got, want, rel);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %.8g, want %.8g, rel %.2e)\n", what, got, want, rel);
  }
}

// ===========================================================================
//  1. 单元级：局部刚度、k=0 释放、平动刚体、旋转投影
// ===========================================================================
static void testLocal() {
  std::printf("\n== 1. 单元级 ==\n");

  // (a) 局部刚度对角结构：K(d,d) = kd，K(d,6+d) = -kd
  {
    std::array<double, 6> ks{};
    ks[0] = 100.0; ks[1] = 200.0; ks[2] = 300.0;
    ks[3] = 400.0; ks[4] = 500.0; ks[5] = 600.0;
    SpringElement3D sp;
    check(sp.init({0, 0, 0}, {1, 0, 0}, ks), "初始化成功");
    check(sp.isValid(), "弹簧单元有效");
    const Mat12 K = sp.localStiffness();
    double worst = 0.0;
    for (int d = 0; d < 6; ++d) {
      const double kd = ks[static_cast<size_t>(d)];
      worst = std::max(worst, std::abs(K(d, d) - kd));
      worst = std::max(worst, std::abs(K(d, 6 + d) + kd));
      worst = std::max(worst, std::abs(K(6 + d, d) + kd));
      worst = std::max(worst, std::abs(K(6 + d, 6 + d) - kd));
    }
    checkNear(worst, 0.0, 1e-15, "局部刚度对角 k 结构（Kdd=kd, Kd,6+d=-kd）");
    double off = 0.0;
    for (int i = 0; i < 12; ++i)
      for (int j = 0; j < 12; ++j) {
        const bool diagPair = (i / 6 == j / 6);
        if (!diagPair && !(i % 6 == j % 6)) off = std::max(off, std::abs(K(i, j)));
      }
    checkNear(off, 0.0, 1e-15, "不同方向间无耦合（副对角线以外为 0）");
  }

  // (b) k=0 释放标记
  {
    std::array<double, 6> ks{};
    ks[0] = 0.0; ks[1] = 5.0; ks[2] = 0.0; ks[3] = 7.0; ks[4] = 0.0; ks[5] = 9.0;
    SpringElement3D sp;
    check(sp.init({0, 0, 0}, {1, 0, 0}, ks), "初始化成功");
    check(sp.isReleasedComp(0) && !sp.isReleasedComp(1) && sp.isReleasedComp(2) &&
          !sp.isReleasedComp(3) && sp.isReleasedComp(4) && !sp.isReleasedComp(5),
          "k=0 方向标记为释放");
  }

  // (c) 平动刚体模式：任意朝向弹簧对整体平动必须零内力
  {
    auto checkRigidTranslate = [&](const Vec3& pi, const Vec3& pj, const Vec3& up,
                                   const char* title) {
      std::array<double, 6> ks{};
      for (double& kd : ks) kd = 1.0;         // 六向全开
      SpringElement3D sp;
      if (!sp.init(pi, pj, ks, up)) { check(false, title); return; }
      if (!sp.isValid()) { check(false, title); return; }
      const Mat12 K = sp.stiffnessGlobal();
      const double kmax = K.maxAbs();
      double worst = 0.0;
      for (int k = 0; k < 3; ++k) {           // 3 个平动
        Vec12d v{};
        for (int n = 0; n < 2; ++n) v[n * 6 + k] = 1.0;
        const Vec12d f = K * v;
        for (double x : f) worst = std::max(worst, std::abs(x));
      }
      checkNear(worst / kmax, 0.0, 1e-14, title);
    };
    checkRigidTranslate({0, 0, 0}, {1, 0, 0}, {0, 0, -1}, "平动刚体零内力（X 向弹簧）");
    checkRigidTranslate({0, 0, 0}, {0, 0, 1}, {0, -1, 0}, "平动刚体零内力（Z 向弹簧）");
    checkRigidTranslate({1, 2, 3}, {4, 6, 5}, {0, 1, 0}, "平动刚体零内力（任意空间弹簧）");
  }

  // (d) 旋转投影：弹簧连 (0,0,0)→(1,1,0)，局部 x' 沿 (1,1)/√2。
  //     单元级验证：给 n0 一个沿全局 X 的单位位移，产生的全局力
  //     = k·(投影)²·êx = k·cos²45°（轴向弹簧的投影刚度）。
  {
    const double k = 400.0;
    const double kPi = std::atan2(0.0, -1.0);
    const double c = std::cos(kPi / 4.0);
    std::array<double, 6> ks{};
    ks[0] = k;
    SpringElement3D sp;
    check(sp.init({0, 0, 0}, {1, 1, 0}, ks, {0, 0, 1}), "45° 弹簧初始化成功");
    check(sp.isValid(), "45° 弹簧有效");
    const Mat12 K = sp.stiffnessGlobal();
    Vec12d v{};
    v[0] = 1.0;                               // n0 沿全局 X 位移 1
    const Vec12d f = K * v;
    // 弹簧只有轴向刚度：n0 位移投影到轴向 = cos45°，轴力 = k·cos45°，
    // 反作用到全局 X = k·cos45°·cos45° = k·cos²45°
    const double fxd = f[0];
    checkNear(fxd, k * c * c, 1e-9, "旋转投影：全局 X 向等效刚度 = k·cos²45°");
    checkNear(f[1], k * c * c, 1e-9, "旋转投影：全局 Y 向等效刚度 = k·cos²45°（对称）");
  }
}

// ===========================================================================
//  2. 轴向弹簧支座：节点 0 --k--> 固定基础，受 +X 力 → u = P/k、N = -P
// ===========================================================================
static void testAxial() {
  std::printf("\n== 2. 轴向弹簧支座 ==\n");
  const double P = 50.0, k = 200.0;
  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id n1 = m.addNode({1, 0, 0}, 1);      // 基础节点（全固定）
  std::array<double, 6> ks{};
  ks[0] = k;
  m.addSpring(n0, n1, ks);
  m.fixNode(n1, true, true, true, true, true, true);
  // 弹簧只约束轴向，其余方向由约束稳住（否则零刚度自由度）
  m.fixNode(n0, false, true, true, true, true, true);
  m.addNodeForce(n0, {P, 0, 0});

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
  checkNear(r.residual, 0.0, 1e-12, "平衡残差 ~ 机器精度");

  const double u = r.u[static_cast<size_t>(n0) * 6 + 0];
  checkNear(u, P / k, 1e-9, "轴向位移 = P/k");
  // 力 +X 推向固定端 → 弹簧受压 → N = -P（拉为正）
  checkNear(r.beamForces[0].ef.N, -P, 1e-9, "弹簧内力 = -P（压缩）");
}

// ===========================================================================
//  3. 串联弹簧：n0--k1--n1--k2--n2（n2 固定，n0 受拉力）
//     ke = k1k2/(k1+k2)、u0 = P/ke、u1 = P/k2、两弹簧内力均 -P
// ===========================================================================
static void testSeries() {
  std::printf("\n== 3. 串联弹簧 ==\n");
  const double P = 30.0, k1 = 150.0, k2 = 300.0;
  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id n1 = m.addNode({1, 0, 0}, 1);
  const Id n2 = m.addNode({2, 0, 0}, 1);
  std::array<double, 6> ksa{}, ksb{};
  ksa[0] = k1; ksb[0] = k2;
  m.addSpring(n0, n1, ksa);
  m.addSpring(n1, n2, ksb);
  m.fixNode(n2, true, true, true, true, true, true);
  m.fixNode(n0, false, true, true, true, true, true);
  m.fixNode(n1, false, true, true, true, true, true);
  m.addNodeForce(n0, {P, 0, 0});

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
  checkNear(r.residual, 0.0, 1e-12, "平衡残差 ~ 机器精度");

  const double ke = k1 * k2 / (k1 + k2);
  checkNear(r.u[static_cast<size_t>(n0) * 6 + 0], P / ke, 1e-9, "u0 = P/ke（串联等效刚度）");
  checkNear(r.u[static_cast<size_t>(n1) * 6 + 0], P / k2, 1e-9, "u1 = P/k2（k2 先承担）");
  checkNear(r.beamForces[0].ef.N, -P, 1e-9, "弹簧 1 内力 = -P（压缩）");
  checkNear(r.beamForces[1].ef.N, -P, 1e-9, "弹簧 2 内力 = -P（压缩）");
}

// ===========================================================================
//  4. 并联弹簧支座：n0 同时连两根弹簧到固定基础，各受 +X 力
//     u = P/(k1+k2)、内力按刚度分配 N1 = -P·k1/(k1+k2)
// ===========================================================================
static void testParallel() {
  std::printf("\n== 4. 并联弹簧支座 ==\n");
  const double P = 60.0, k1 = 100.0, k2 = 300.0;
  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id n1 = m.addNode({0, 1, 0}, 1);      // 基础 1
  const Id n2 = m.addNode({0, 2, 0}, 1);      // 基础 2
  std::array<double, 6> ksa{}, ksb{};
  ksa[0] = k1; ksb[0] = k2;
  m.addSpring(n0, n1, ksa);                   // 竖向弹簧（沿 Y）
  m.addSpring(n0, n2, ksb);
  m.fixNode(n1, true, true, true, true, true, true);
  m.fixNode(n2, true, true, true, true, true, true);
  // n0 仅 uy 自由（弹簧沿 Y 约束该方向）
  m.fixNode(n0, true, false, true, true, true, true);
  m.addNodeForce(n0, {0, P, 0});

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
  checkNear(r.residual, 0.0, 1e-12, "平衡残差 ~ 机器精度");

  const double u = r.u[static_cast<size_t>(n0) * 6 + 1];
  checkNear(u, P / (k1 + k2), 1e-9, "u = P/(k1+k2)（并联等效刚度）");
  checkNear(r.beamForces[0].ef.N, -P * k1 / (k1 + k2), 1e-9, "弹簧 1 内力按刚度分配");
  checkNear(r.beamForces[1].ef.N, -P * k2 / (k1 + k2), 1e-9, "弹簧 2 内力按刚度分配");
}

// ===========================================================================
//  5. 倾斜弹簧（45°）：全局 X 等效刚度 = k·cos²45°、轴力 = -P/cos45°
// ===========================================================================
static void testTilt() {
  std::printf("\n== 5. 45° 倾斜弹簧 ==\n");
  const double P = 20.0, k = 400.0;
  const double kPi = std::atan2(0.0, -1.0);
  const double c = std::cos(kPi / 4.0);
  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id n1 = m.addNode({c, c, 0}, 1);
  std::array<double, 6> ks{};
  ks[0] = k;
  m.addSpring(n0, n1, ks, Vec3{0, 0, 1});
  m.fixNode(n1, true, true, true, true, true, true);
  // 45° 弹簧约束轴向 (1,1)/√2；垂直方向 (面内) 是机构 → 固定 uy，只留 ux
  m.fixNode(n0, false, true, true, true, true, true);
  m.addNodeForce(n0, {P, 0, 0});

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
  checkNear(r.residual, 0.0, 1e-12, "平衡残差 ~ 机器精度");

  const double ux = r.u[static_cast<size_t>(n0) * 6 + 0];
  checkNear(ux, P / (k * c * c), 1e-9, "ux = P/(k·cos²45°)（投影刚度）");
  // n0 沿 +X 位移投影到局部轴为负（压缩）→ N = -P/cos45°
  checkNear(r.beamForces[0].ef.N, -P / c, 1e-9, "轴力 = -P/cos45°（压缩）");
}

// ===========================================================================
//  6. k=0 释放语义：
//     (a) 无荷载方向被 DofNumbering 自动锁定（dof=-1）→ 不奇异；
//     (b) 该方向带节点荷载 → check() 必须拦截（防荷载被吞）
// ===========================================================================
static void testZeroStiff() {
  std::printf("\n== 6. k=0 释放语义 ==\n");
  const double k = 500.0;

  // (a) 无荷载 → 零刚度自由度自动锁定
  {
    Model m;
    const Id n0 = m.addNode({0, 0, 0}, 1);
    const Id n1 = m.addNode({1, 0, 0}, 1);
    std::array<double, 6> ks{};
    ks[0] = k;
    m.addSpring(n0, n1, ks);
    m.fixNode(n1, true, true, true, true, true, true);
    m.assignDofs();
    bool locked = true;
    for (int d = 1; d < 6; ++d)
      if (m.node(n0).dof[d] >= 0) { locked = false; break; }
    check(locked, "k=0 方向无荷载 → 自动锁定（dof=-1）");
    StaticAnalysis sa(m);
    const auto r = sa.solve();
    check(r.ok, "零刚度自动锁定下求解不奇异");
    if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
    checkNear(r.u[static_cast<size_t>(n0) * 6 + 0], 0.0, 1e-12, "无荷载 → ux = 0");
  }

  // (b) 带节点荷载 → check() 拦截
  {
    const double P = 10.0;
    Model m;
    const Id n0 = m.addNode({0, 0, 0}, 1);
    const Id n1 = m.addNode({1, 0, 0}, 1);
    std::array<double, 6> ks{};
    ks[0] = k;
    m.addSpring(n0, n1, ks);
    m.fixNode(n1, true, true, true, true, true, true);
    m.addNodeForce(n0, {0, P, 0});            // 弹簧不提供刚度的 uy 方向
    m.assignDofs();
    auto d = m.check();
    bool found = false;
    for (const auto& e : d.errors)
      if (e.find("无刚度") != std::string::npos) { found = true; break; }
    check(found, "零刚度方向带荷载 → check() 拦截（防荷载被吞）");
  }
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  弹簧单元验证\n");
  std::printf("======================================================\n");

  testLocal();
  testAxial();
  testSeries();
  testParallel();
  testTilt();
  testZeroStiff();

  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}