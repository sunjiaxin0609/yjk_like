// =============================================================================
//  tests/test_floor.cpp  ——  刚性楼板（MPC）+ 板荷载导荷 验证
//
//  这两件事合起来才构成"不建板也能算"：
//    ① 刚性楼板假定  —— 平面内刚度靠它
//    ② 45° 线导荷    —— 楼面荷载的传递路径靠它
//  缺任何一条，模型要么"散"、要么荷载凭空消失，而两者都不报错。
//
//  【本文件里最重要的一条断言】
//    "刚性楼板不得产生虚假刚度"——
//    对称荷载下，带刚性楼板与不带刚性楼板的位移必须一致。
//    这一条抓到过一个真实的严重 bug：组装时"只遍历局部上三角、
//    靠对称性补齐另一半"的技巧在引入 MPC 之后失效，导致每根梁的
//    轴向刚度残留下来（本该精确抵消为 0），结构比不带楼板时
//    【刚 500 倍】—— 而残差、反力平衡、刚体一致性全都正常。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/model/GridMesh.h"

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
//  构造一个"4 柱 + 一圈梁"的单层模型（可选加刚性楼板）
// -----------------------------------------------------------------------------
namespace {

struct Frame {
  Model m;
  Id base[4]{}, top[4]{}, master{-1};
  bool withDiaphragm{false};
};

Frame makeFrame(double b, double h, bool diaphragm) {
  Frame f;
  f.withDiaphragm = diaphragm;
  const auto cs = section::rect(0.5, 0.5);
  const auto bs = section::rect(0.3, 0.6);
  const auto mat = Material::concreteC(30);
  const double xs[4] = {0, b, b, 0}, ys[4] = {0, 0, b, b};
  for (int i = 0; i < 4; ++i) {
    f.base[i] = f.m.addNode({xs[i], ys[i], 0}, 0);
    f.top[i] = f.m.addNode({xs[i], ys[i], h}, 1);
    f.m.addBeam(f.base[i], f.top[i], cs, mat, Vec3{0, -1, 0});
  }
  for (int i = 0; i < 4; ++i)
    f.m.addBeam(f.top[i], f.top[(i + 1) % 4], bs, mat, Vec3{0, 0, -1});
  for (int i = 0; i < 4; ++i) f.m.fixAll(f.base[i]);
  if (diaphragm)
    f.master = f.m.addRigidDiaphragm({f.top[0], f.top[1], f.top[2], f.top[3]}, 1, true);
  return f;
}

double solveTopUx(Model& m, Id probe, const char* what) {
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  if (!r.ok) { std::printf("  [FAIL] %s：求解失败 %s\n", what, r.message.c_str()); ++g_fail; return 0.0; }
  return r.u[static_cast<size_t>(probe) * 6 + 0];
}

}  // namespace

// -----------------------------------------------------------------------------
//  1. 刚性楼板的构造
// -----------------------------------------------------------------------------
static void testConstruction() {
  std::printf("\n== 1. 刚性楼板：构造与自由度编号 ==\n");

  Frame a = makeFrame(6.0, 4.0, false);
  a.m.assignDofs();
  const Id freePlain = a.m.freeDofCount();
  checkInt(freePlain, 24, "无楼板：自由度数 = 8 节点 × 6 − 4 柱脚 × 6 = 24");

  Frame b = makeFrame(6.0, 4.0, true);
  const Node& mn = b.m.node(b.master);
  checkInt(b.master, 8, "主节点建在节点序列末尾");
  checkNear(mn.r.x, 3.0, 1e-14, "主节点在形心 x");
  checkNear(mn.r.y, 3.0, 1e-14, "主节点在形心 y");
  check(mn.fixed[2] && mn.fixed[3] && mn.fixed[4],
        "主节点只保留 (ux, uy, rz)：uz / rx / ry 被固定");
  check(!mn.fixed[0] && !mn.fixed[1] && !mn.fixed[5], "主节点的 (ux, uy, rz) 自由");

  b.m.assignDofs();
  // 9 节点 × 6 = 54；柱脚固定 24；主节点固定 3；从属 4×3 = 12
  checkInt(b.m.freeDofCount(), 15, "有楼板：自由度数 = 54 − 24 − 3 − 12 = 15");

  int nSlave = 0;
  for (Id i = 0; i < b.m.nodeCount(); ++i)
    for (int k = 0; k < 6; ++k)
      if (b.m.node(i).dof[static_cast<size_t>(k)] == kDofSlave) ++nSlave;
  checkInt(nSlave, 12, "从属自由度 12 个（4 个顶层节点 × ux/uy/rz），不占方程号");

  // 从属约束的系数：ux = UX − Θ·dy，uy = UY + Θ·dx
  const Id t0 = b.top[0];                        // 位于 (0,0)，形心 (3,3) ⇒ dx = dy = −3
  const Id mid = b.master;
  const auto& links = b.m.dofLinks();
  auto coeff = [&](Id slaveDof, Id masterDof) -> double {
    for (const auto& L : links) {
      if (L.slave != slaveDof) continue;
      for (const auto& mm : L.masters) if (mm.first == masterDof) return mm.second;
    }
    return 0.0;
  };
  checkNear(coeff(t0 * 6 + 0, mid * 6 + 5), 3.0, 1e-14,
            "ux 的 Θ 系数 = −dy（dy = −3 ⇒ 系数 +3）");
  checkNear(coeff(t0 * 6 + 1, mid * 6 + 5), -3.0, 1e-14,
            "uy 的 Θ 系数 = +dx（dx = −3）");
  checkNear(coeff(t0 * 6 + 0, mid * 6 + 0), 1.0, 1e-14, "ux 对 UX 的系数 = 1");

  // 校核：主节点不该被误报为悬空/零刚度
  b.m.assignDofs();
  const auto diag = b.m.check();
  for (const auto& e : diag.errors) std::printf("       [错误] %s\n", e.c_str());
  check(diag.ok(), "带刚性楼板的模型通过校核（主节点不被误报为悬空节点）");
}

// -----------------------------------------------------------------------------
//  2. 刚性楼板：楼层做平面刚体运动
// -----------------------------------------------------------------------------
static void testRigidBody() {
  std::printf("\n== 2. 刚性楼板：楼层做平面刚体运动 ==\n");

  Frame f = makeFrame(6.0, 4.0, true);
  // 偏心水平力：只在 (6,6) 施加，同时激起平动与扭转
  f.m.addNodeForce(f.top[2], Vec3{100, 0, 0});
  f.m.assignDofs();
  StaticAnalysis sa(f.m);
  const auto r = sa.solve();
  check(r.ok, "偏心水平力求解成功");
  if (!r.ok) return;

  const Vec3 c = f.m.node(f.master).r;
  const double UX = r.u[static_cast<size_t>(f.master) * 6 + 0];
  const double UY = r.u[static_cast<size_t>(f.master) * 6 + 1];
  const double TH = r.u[static_cast<size_t>(f.master) * 6 + 5];
  std::printf("       主节点 UX=%+.6e UY=%+.6e Θ=%+.6e\n", UX, UY, TH);
  check(std::abs(TH) > 1e-12, "偏心荷载激起楼层转动（Θ ≠ 0）");
  check(TH < 0.0, "力作用于 (6,6)、形心 (3,3)：绕 Z 的力矩为负 ⇒ Θ < 0");

  double worst = 0.0, scale = 1e-30;
  for (int i = 0; i < 4; ++i) {
    const double dx = f.m.node(f.top[i]).r.x - c.x;
    const double dy = f.m.node(f.top[i]).r.y - c.y;
    const double ex = UX - TH * dy;
    const double ey = UY + TH * dx;
    const double ax = r.u[static_cast<size_t>(f.top[i]) * 6 + 0];
    const double ay = r.u[static_cast<size_t>(f.top[i]) * 6 + 1];
    worst = std::max(worst, std::max(std::abs(ax - ex), std::abs(ay - ey)));
    scale = std::max(scale, std::max(std::abs(ex), std::abs(ey)));
  }
  std::printf("       刚体运动最大偏差 %.3e（量级 %.3e）\n", worst, scale);
  check(worst / scale < 1e-12, "四个顶层节点的位移严格满足平面刚体运动");

  // 顶层节点的 rz 必须等于楼层转角
  double worstRz = 0.0;
  for (int i = 0; i < 4; ++i)
    worstRz = std::max(worstRz,
                       std::abs(r.u[static_cast<size_t>(f.top[i]) * 6 + 5] - TH) / std::abs(TH));
  check(worstRz < 1e-12, "顶层节点的 rz = 楼层转角 Θ");
}

// -----------------------------------------------------------------------------
//  3. ★ 刚性楼板不得产生虚假刚度（回归测试）
// -----------------------------------------------------------------------------
static void testNoSpuriousStiffness() {
  std::printf("\n== 3. ★ 刚性楼板不产生虚假刚度（对称荷载对照） ==\n");

  // 对称水平力：四个角点各 25 kN，合计 100 kN。
  // 此时顶层四个节点的 ux 本来就相等，所以"带楼板"与"不带楼板"
  // 的位移必须一致 —— 楼板只是把等值约束写出来，不该改变任何刚度。
  auto run = [](bool diaphragm) -> double {
    Frame f = makeFrame(6.0, 4.0, diaphragm);
    for (int i = 0; i < 4; ++i) f.m.addNodeForce(f.top[i], Vec3{25, 0, 0});
    return solveTopUx(f.m, f.top[0], "对称荷载");
  };
  const double d0 = run(false);
  const double d1 = run(true);
  std::printf("       无楼板 ux = %.9e    有楼板 ux = %.9e\n", d0, d1);
  check(d0 > 0.0, "无楼板：位移为正（沿荷载方向）");

  // 【判据】两者必须一致。旧实现（只遍历局部上三角）在这里会差约 500 倍。
  checkNear(d1 / d0, 1.0, 1e-9, "★ 对称荷载下带/不带刚性楼板的位移完全一致");

  // 用位移反推侧移刚度，与"固定-固定柱"的手算量级对账
  Frame f = makeFrame(6.0, 4.0, false);
  const double EI = Material::concreteC(30).E * std::pow(0.5, 4) / 12.0;
  const double kFF = 12.0 * EI / (4.0 * 4.0 * 4.0);      // 单柱 固定-固定
  const double kCant = 3.0 * EI / (4.0 * 4.0 * 4.0);     // 单柱 悬臂
  const double K = 100.0 / d0;
  std::printf("       反推侧移刚度 K = %.5g kN/m（4 柱固定-固定 %.5g / 悬臂 %.5g）\n",
              K, 4 * kFF, 4 * kCant);
  // 有梁参与框架作用 ⇒ 落在"悬臂"与"固定-固定"之间（含 Timoshenko 剪切会略软）
  check(K > 4 * kCant * 0.85 && K < 4 * kFF * 1.05,
        "侧移刚度落在 4 柱悬臂~固定-固定之间（含框架作用与剪切变形）");
  (void)f;

  // 同样的对照放在"带梁加密"的轴网模型上（更接近真实结构）
  auto runGrid = [](bool diaphragm) -> double {
    AxisGrid g;
    g.addAxisX(0.0); g.addAxisX(6.0);
    g.addAxisY(0.0); g.addAxisY(6.0);
    g.addStory(0.0); g.addStory(3.6);
    Model m;
    GridMesh::Spec sp;
    sp.columnSection = section::rect(0.5, 0.5);
    sp.beamSection = section::rect(0.3, 0.6);
    sp.buildSlabs = false;
    GridMesh mesh(g, sp);
    mesh.generate(m);
    for (Id n : mesh.baseNodes()) m.fixAll(n);
    if (diaphragm) mesh.addRigidDiaphragms(m, 1, true);
    for (Id ix = 0; ix < 2; ++ix)
      for (Id jy = 0; jy < 2; ++jy) {
        const Id n = mesh.nodeAtAxis(ix, jy, 1);
        if (n >= 0) m.addNodeForce(n, Vec3{25, 0, 0});
      }
    m.assignDofs();
    StaticAnalysis sa(m);
    const auto r = sa.solve();
    if (!r.ok) return 0.0;
    // 取四个柱顶节点的 ux 平均值
    double s = 0.0;
    for (Id ix = 0; ix < 2; ++ix)
      for (Id jy = 0; jy < 2; ++jy) s += r.u[static_cast<size_t>(mesh.nodeAtAxis(ix, jy, 1)) * 6 + 0];
    return s / 4.0;
  };
  const double g0 = runGrid(false);
  const double g1 = runGrid(true);
  std::printf("       轴网模型：无楼板 ux(平均) = %.9e    有楼板 = %.9e\n", g0, g1);
  checkNear(g1 / g0, 1.0, 1e-9, "★ 轴网模型上同样完全一致");
}

// -----------------------------------------------------------------------------
//  4. 刚性楼板下的整体平衡
// -----------------------------------------------------------------------------
static void testEquilibrium() {
  std::printf("\n== 4. 刚性楼板下的整体平衡 ==\n");

  Frame f = makeFrame(6.0, 4.0, true);
  f.m.addNodeForce(f.top[0], Vec3{30, 10, -50});
  f.m.addNodeForce(f.top[2], Vec3{-15, 20, 0});
  f.m.assignDofs();
  StaticAnalysis sa(f.m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) return;

  // -------------------------------------------------------------------------
  //  【正确做法：用"刚体位移模式"做投影，而不是手写 λ×p 的力矩公式】
  //
  //  对任意刚体模式 v（平动或绕原点转动），都有 vᵀ·K = 0（K 的零空间），
  //  于是   Σ_i [ 平动分量·(F_i + R_i) + 转动分量·(M_i + RM_i) ] = 0
  //
  //  这里 v 的"转动分量"不只是 rz —— 刚体绕 z 转 Θ 时，节点的位移是
  //  (a + ω×p)，而节点的转角分量是 ω。只写 Σ(p×R) 会漏掉
  //  "施加弯矩 × 转角模式分量" 这一项，判据就不完整（本项目踩过）。
  //
  //  用完整的刚体模式做投影，判据对任意荷载（含集中力矩）都成立。
  // -------------------------------------------------------------------------
  double ref = 1.0;
  for (Id i = 0; i < f.m.nodeCount(); ++i)
    for (int k = 0; k < 3; ++k) ref = std::max(ref, std::abs(f.m.node(i).force[k]));

  // 平动：Σ 全部反力 = −Σ 全部外力
  for (int k = 0; k < 3; ++k) {
    double s = 0.0;
    for (Id i = 0; i < f.m.nodeCount(); ++i)
      s += f.m.node(i).force[k] + r.reaction[static_cast<size_t>(i) * 6 + k];
    std::printf("       平动 %d: ΣF + ΣR = %+.3e\n", k, s);
    check(std::abs(s) < 1e-9 * ref, "平动平衡 ΣR = −ΣF");
  }

  // -------------------------------------------------------------------------
  //  转动平衡：必须只用【真正的支座】的反力
  //
  //  【为什么不能把从属自由度的"反力"算进来】
  //  从属自由度上的那三个数是刚性楼板【内部】的约束力 ——
  //  楼板把这层节点拽住的那份力。它由楼板自身承受，不传给地基。
  //  把内部力当成支座反力去做整体力矩平衡，结果当然不闭合。
  //  （实测：用全部反力时绕 Y 轴差 120 kN·m，用真正的支座时差 1e-14。）
  //
  //  判据：以固定支座的反力为外力，与全部外荷载一起做整体力矩平衡。
  // -------------------------------------------------------------------------
  for (int k = 0; k < 3; ++k) {
    Vec3 wv{0, 0, 0};
    wv[k] = 1.0;
    double s = 0.0;
    for (Id i = 0; i < f.m.nodeCount(); ++i) {
      const bool isSupport =
          (f.m.node(i).dof[0] == kDofFixed) || (f.m.node(i).dof[1] == kDofFixed) ||
          (f.m.node(i).dof[2] == kDofFixed);
      if (isSupport) {
        const Vec3 t = cross(wv, f.m.node(i).r);
        const Vec3 R{r.reaction[static_cast<size_t>(i) * 6 + 0],
                     r.reaction[static_cast<size_t>(i) * 6 + 1],
                     r.reaction[static_cast<size_t>(i) * 6 + 2]};
        const Vec3 RM{r.reaction[static_cast<size_t>(i) * 6 + 3],
                      r.reaction[static_cast<size_t>(i) * 6 + 4],
                      r.reaction[static_cast<size_t>(i) * 6 + 5]};
        s += dot(t, R) + dot(wv, RM);
      }
      // 外荷载（全部节点）
      const Vec3 t = cross(wv, f.m.node(i).r);
      s += dot(t, f.m.node(i).force) + dot(wv, f.m.node(i).moment);
    }
    std::printf("       关于 e%d 的力矩: Σ 支座 + Σ 外荷载 = %+.3e\n", k, s);
    check(std::abs(s) < 1e-9 * ref * 20.0, "转动平衡：支座反力矩 + 外荷载力矩 = 0");
  }
  check(r.residual < 1e-10, "解残差 < 1e-10");
}

// -----------------------------------------------------------------------------
//  5. 刚性楼板的逐层开关
// -----------------------------------------------------------------------------
static void testPerStorySwitch() {
  std::printf("\n== 5. 刚性楼板：逐层开关 ==\n");

  AxisGrid g;
  g.addAxisX(0.0); g.addAxisX(6.0);
  g.addAxisY(0.0); g.addAxisY(6.0);
  g.addStory(0.0); g.addStory(3.6); g.addStory(7.2);
  Model m;
  GridMesh::Spec sp;
  sp.buildSlabs = false;
  GridMesh mesh(g, sp);
  mesh.generate(m);

  const Id made = mesh.addRigidDiaphragms(m, 2, true);
  checkInt(made, 1, "fromStory = 2 ⇒ 只做 1 层刚性楼板");
  checkInt(static_cast<long long>(m.diaphragms().size()), 1, "记录 1 个刚性楼板");
  if (!m.diaphragms().empty())
    checkInt(m.diaphragms()[0].story, 2, "作用在第 2 层");

  // coupleRz = false：只约束面内平动
  Model m2;
  GridMesh mesh2(g, sp);
  mesh2.generate(m2);
  mesh2.addRigidDiaphragms(m2, 1, false);
  m2.assignDofs();
  int nSlave = 0;
  for (Id i = 0; i < m2.nodeCount(); ++i)
    for (int k = 0; k < 6; ++k)
      if (m2.node(i).dof[static_cast<size_t>(k)] == kDofSlave) ++nSlave;
  // 第 1、2 层各 4 个节点（无剖分），每节点 ux/uy 两项 → 8×2 = 16
  checkInt(nSlave, 16, "coupleRz = false ⇒ 只约束 ux / uy，不约束 rz");
}

// -----------------------------------------------------------------------------
//  6. 45° 线导荷：守恒与分配形状
// -----------------------------------------------------------------------------
namespace {

struct BayLoads {
  int nSeg{0};
  double total{0.0};
  double maxAlongX{0.0};     // X 向梁上的最大线荷载强度
  double maxAlongY{0.0};
  int segsAlongX{0};
  int segsAlongY{0};
  double plateau{0.0};       // X 向梁上"恒定段"的数量（梯形平台）
};

BayLoads collectBayLoads(const Model& m) {
  BayLoads r;
  for (const auto& e : m.elements()) {
    if (e->type() != ElementType::Beam3D) continue;
    const auto* b = static_cast<const BeamElement*>(e.get());
    const Vec3 ax = b->core().localX();
    const bool alongX = std::abs(ax.x) > 0.7;
    (void)0;
    for (const auto& s : b->loadSegments()) {
      const double q1 = -s.q1.z, q2 = -s.q2.z;    // 向下为正
      if (q1 <= 0.0 && q2 <= 0.0) continue;
      ++r.nSeg;
      r.total += 0.5 * (q1 + q2) * std::abs(s.x2 - s.x1);
      // 【平台段检测与梁的方向无关】曾把它写在"沿 X"分支里，
      // 而 4×6 板格的梯形恰好在 Y 向梁上 —— 于是统计出 0。
      // 统计代码的分支位置也是会写错的。
      if (std::abs(q1 - q2) < 1e-12 && q1 > 0.0) r.plateau += std::abs(s.x2 - s.x1);
      if (alongX) {
        ++r.segsAlongX;
        r.maxAlongX = std::max(r.maxAlongX, std::max(q1, q2));
      } else {
        ++r.segsAlongY;
        r.maxAlongY = std::max(r.maxAlongY, std::max(q1, q2));
      }
    }
  }
  return r;
}

struct BayCase {
  BayLoads loads;
  double reaction{0.0};
  double want{0.0};
  bool solved{false};
};

BayCase runBay(double LX, double LY, double q) {
  BayCase out;
  AxisGrid g;
  g.addAxisX(0.0); g.addAxisX(LY);
  g.addAxisY(0.0); g.addAxisY(LX);
  g.addStory(0.0); g.addStory(3.6);
  Model m;
  GridMesh::Spec sp;
  sp.columnSection = section::rect(0.5, 0.5);
  sp.beamSection = section::rect(0.3, 0.6);
  sp.buildSlabs = false;
  GridMesh mesh(g, sp);
  mesh.generate(m);
  for (Id n : mesh.baseNodes()) m.fixAll(n);
  mesh.addRigidDiaphragms(m, 1, true);

  GridMesh::SlabLoadOptions opt;
  opt.q = q;
  opt.path = GridMesh::SlabLoadOptions::Path::ToBeams;
  mesh.applySlabLoad(m, opt);

  out.loads = collectBayLoads(m);
  out.want = q * LX * LY;
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  out.solved = r.ok;
  if (r.ok)
    for (Id n : mesh.baseNodes()) out.reaction += r.reaction[static_cast<size_t>(n) * 6 + 2];
  return out;
}

}  // namespace

static void testDistribution() {
  std::printf("\n== 6. 45° 线导荷：守恒与分配形状 ==\n");

  // (a) 正方形：四边都是三角形，峰值 = q·L/2
  {
    const double q = 5.0;
    const BayCase c = runBay(6.0, 6.0, q);
    std::printf("       正方形 6×6: %d 段，总荷载 %.6f kN\n", c.loads.nSeg, c.loads.total);
    checkNear(c.loads.total, c.want, 1e-12, "正方形：导到梁上的总荷载 = q·面积");
    checkNear(c.loads.maxAlongX, q * 6.0 / 2.0, 1e-12, "正方形：X 向梁峰值 = q·L/2 = 15");
    checkNear(c.loads.maxAlongY, q * 6.0 / 2.0, 1e-12, "正方形：Y 向梁峰值 = q·L/2 = 15");
    check(c.solved, "正方形：求解成功");
    checkNear(c.reaction, c.want, 1e-10, "正方形：竖向反力合计 = 总荷载");
  }

  // (b) 4×6 双向板：短边三角形、长边梯形，峰值都是 q·L短/2
  {
    const double q = 5.0;
    const BayCase c = runBay(4.0, 6.0, q);
    std::printf("       矩形 4×6: %d 段，总荷载 %.6f kN，X 向 %d 段 / Y 向 %d 段\n",
                c.loads.nSeg, c.loads.total, c.loads.segsAlongX, c.loads.segsAlongY);
    checkNear(c.loads.total, c.want, 1e-12, "矩形 4×6：总荷载 = q·面积（守恒）");
    checkNear(c.loads.maxAlongX, q * 4.0 / 2.0, 1e-12, "短边（4 m）三角形峰值 = q·L短/2 = 10");
    checkNear(c.loads.maxAlongY, q * 4.0 / 2.0, 1e-12, "长边（6 m）梯形峰值 = q·L短/2 = 10");
    checkNear(c.loads.plateau, 4.0, 1e-9,
              "两条长边的梯形平台段合计 = 2 × (L长 − L短) = 4 m");
    check(c.solved, "矩形 4×6：求解成功");
    checkNear(c.reaction, c.want, 1e-10, "矩形 4×6：竖向反力合计 = 总荷载");
  }

  // (c) 单向板：长短边比 ≥ 2 ⇒ 荷载全部由长边方向的两根梁承担
  {
    const double q = 5.0;
    const BayCase c = runBay(4.0, 8.0, q);
    std::printf("       单向板 4×8: %d 段，总荷载 %.6f kN，X 向 %d 段 / Y 向 %d 段\n",
                c.loads.nSeg, c.loads.total, c.loads.segsAlongX, c.loads.segsAlongY);
    checkInt(c.loads.segsAlongX, 0, "单向板：短边方向的梁【不分荷载】");
    check(c.loads.segsAlongY > 0, "单向板：荷载全部落在长边方向的梁上");
    checkNear(c.loads.maxAlongY, q * 4.0 / 2.0, 1e-12, "单向板：均布强度 = q·L短/2 = 10");
    checkNear(c.loads.total, c.want, 1e-12, "单向板：总荷载仍守恒");
    check(c.solved, "单向板：求解成功");
    checkNear(c.reaction, c.want, 1e-10, "单向板：竖向反力合计 = 总荷载");
  }

  // (d) 阈值可调：把 oneWayRatio 放到 2.5，4×6（比值 1.5）仍是双向
  {
    AxisGrid g;
    g.addAxisX(0.0); g.addAxisX(6.0);
    g.addAxisY(0.0); g.addAxisY(4.0);
    g.addStory(0.0); g.addStory(3.6);
    Model m;
    GridMesh::Spec sp;
    sp.buildSlabs = false;
    GridMesh mesh(g, sp);
    mesh.generate(m);
    for (Id n : mesh.baseNodes()) m.fixAll(n);
    GridMesh::SlabLoadOptions opt;
    opt.q = 5.0;
    opt.path = GridMesh::SlabLoadOptions::Path::ToBeams;
    opt.oneWayRatio = 10.0;
    const auto rep = mesh.applySlabLoad(m, opt);
    checkInt(rep.oneWayBays, 0, "阈值为 10 时 4×6 板格按双向板处理");
    check(rep.beamSegs > 0, "阈值为 10 时仍有荷载导到梁上");
  }

  // (e) 缺边梁必须报警，不能静默丢荷载
  {
    AxisGrid g;
    g.addAxisX(0.0); g.addAxisX(6.0);
    g.addAxisY(0.0); g.addAxisY(6.0);
    g.addStory(0.0); g.addStory(3.6);
    g.setBeamAlongX(0, 0, false);              // 抽掉下边的梁
    Model m;
    GridMesh::Spec sp;
    sp.buildSlabs = false;
    GridMesh mesh(g, sp);
    mesh.generate(m);
    for (Id n : mesh.baseNodes()) m.fixAll(n);
    GridMesh::SlabLoadOptions opt;
    opt.q = 5.0;
    opt.path = GridMesh::SlabLoadOptions::Path::ToBeams;
    const auto rep = mesh.applySlabLoad(m, opt);
    check(!rep.warnings.empty(), "缺少边梁时给出警告（不静默丢荷载）");
    if (!rep.warnings.empty())
      std::printf("       [警告] %s\n", rep.warnings[0].c_str());
  }
}

// -----------------------------------------------------------------------------
//  7. 不建板 + 导荷 + 刚性楼板：整条链路
// -----------------------------------------------------------------------------
static void testEndToEnd() {
  std::printf("\n== 7. 不建板 + 导荷 + 刚性楼板：端到端 ==\n");

  AxisGrid g;
  g.addAxisX(0.0); g.addAxisX(6.0); g.addAxisX(12.0);
  g.addAxisY(0.0); g.addAxisY(6.0); g.addAxisY(12.0);
  g.addStory(0.0); g.addStory(3.6); g.addStory(7.2);
  Model m;
  GridMesh::Spec sp;
  sp.columnSection = section::rect(0.5, 0.5);
  sp.beamSection = section::rect(0.3, 0.6);
  sp.buildSlabs = false;                       // 不建板
  GridMesh mesh(g, sp);
  const auto rr = mesh.generate(m);
  check(rr.ok, "网格生成成功");
  checkInt(rr.slabs, 0, "不建板：板单元数 = 0");
  std::printf("       %s\n", rr.message().c_str());

  for (Id n : mesh.baseNodes()) m.fixAll(n);
  checkInt(mesh.addRigidDiaphragms(m, 1, true), 2, "两层刚性楼板");

  GridMesh::SlabLoadOptions opt;
  opt.q = 4.0;
  const auto rep = mesh.applySlabLoad(m, opt);
  std::printf("       %s\n", rep.message().c_str());
  for (const auto& w : rep.warnings) std::printf("       [警告] %s\n", w.c_str());
  check(rep.warnings.empty(), "导荷无警告");
  checkNear(rep.total, 4.0 * 12.0 * 12.0 * 2.0, 1e-12, "两层楼面总荷载 = 4 × 144 × 2 = 1152 kN");

  // 逐节点施加顶层水平力（有刚性楼板时会自动折到主节点）
  int cnt = 0;
  for (Id i = 0; i < m.nodeCount(); ++i) {
    if (m.node(i).diaphragmMaster) continue;
    if (m.node(i).story == 2) { m.addNodeForce(i, Vec3{10, 0, 0}); ++cnt; }
  }
  const double totalFx = 10.0 * cnt;
  std::printf("       顶层 %d 个节点各 10 kN ⇒ 总水平力 %.1f kN\n", cnt, totalFx);

  m.assignDofs();
  const auto diag = m.check();
  for (const auto& e : diag.errors) std::printf("       [错误] %s\n", e.c_str());
  check(diag.ok(), "模型校核通过");

  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) { std::printf("       %s\n", r.message.c_str()); return; }
  std::printf("       自由度 %d，残差 %.2e\n", r.ndof, r.residual);
  check(r.residual < 1e-10, "残差 < 1e-10");

  double sRz = 0, sRx = 0;
  for (Id i = 0; i < m.nodeCount(); ++i) {
    sRz += r.reaction[static_cast<size_t>(i) * 6 + 2];
    sRx += r.reaction[static_cast<size_t>(i) * 6 + 0];
  }
  checkNear(sRz, rep.total, 1e-10, "竖向反力合计 = 楼面总荷载");
  checkNear(std::abs(sRx), totalFx, 1e-10, "水平反力合计 = 顶层总水平力");

  // 【关键】顶层水平力必须原封不动地传到各层柱上 ——
  // 只折了刚度、没折荷载（或反过来）时，位移会偏小而残差仍正常。
  {
    // -----------------------------------------------------------------------
    //  顶层柱的水平剪力【按矢量求和】，不能按 Σ|V| 求和
    //
    //  竖向荷载下框架会"内缩"：柱因梁端弯矩不平衡而出现自平衡的水平剪力
    //  （实测角柱 Vy = ∓7.29 kN、边柱 ∓14.1 kN，Σ 恰好为零）。
    //  按 Σ|V| 累加会把这组自平衡量当成剪力计入 —— 实测把 90 kN 算成 130.5 kN。
    //  本项目的后处理层就是为此改用矢量求和的，测里必须保持一致。
    // -----------------------------------------------------------------------
    double sx = 0.0, sy = 0.0;
    for (Id e = 0; e < m.elementCount(); ++e) {
      const Element& el = *m.elements()[static_cast<size_t>(e)];
      if (el.type() != ElementType::Beam3D) continue;
      const auto* be = static_cast<const BeamElement*>(&el);
      if (std::abs(be->core().localX().z) < 0.7) continue;
      const int sa_ = m.node(el.nodes()[0]).story, sb = m.node(el.nodes()[1]).story;
      if (std::max(sa_, sb) != 2) continue;
      const auto& ef = r.beamForces[static_cast<size_t>(e)].ef;
      // 柱的局部：ex = +Z，ey = −Y，ez = +X ⇒ Vz 沿 +X，Vy 沿 −Y
      sx += ef.Vz;
      sy += -ef.Vy;
    }
    std::printf("       顶层柱水平剪力 ΣVx = %.6f kN，ΣVy = %.3e kN（总水平力 %.1f）\n",
                sx, sy, totalFx);
    checkNear(sx, totalFx, 1e-10, "顶层柱 X 向剪力合计 = 施加的总水平力");
    check(std::abs(sy) < 1e-8 * totalFx,
          "顶层柱 Y 向剪力自平衡（竖向荷载下的框架内缩，Σ 为零）");
  }
}

// -----------------------------------------------------------------------------
//  8. 梁的梯形分布荷载（导荷的直接下游）
// -----------------------------------------------------------------------------
static void testTrapezoidBeamLoad() {
  std::printf("\n== 8. 梁上的三角形 / 梯形分布荷载 ==\n");

  const auto sec = section::rect(0.3, 0.6);
  const auto mat = Material::concreteC(30);

  // (a) 悬臂 + 三角形荷载（固定端最大）：V = q0(L−x)²/2L，M = −q0(L−x)³/6L
  {
    const double L = 6.0, q0 = 12.0;
    const int n = 48;
    Model m;
    std::vector<Id> ns;
    for (int i = 0; i <= n; ++i) ns.push_back(m.addNode({L * i / n, 0, 0}, 0));
    auto qf = [&](double x) { return q0 * (1.0 - x / L); };
    for (int i = 0; i < n; ++i) {
      const double x0 = L * i / n, x1 = L * (i + 1) / n;
      m.addBeam(ns[i], ns[i + 1], sec, mat, Vec3{0, 0, -1})
          ->addLineLoadSegment(0.0, x1 - x0, Vec3{0, 0, -qf(x0)}, Vec3{0, 0, -qf(x1)});
    }
    m.fixAll(ns[0]);
    m.assignDofs();
    StaticAnalysis sa(m);
    const auto r = sa.solve();
    check(r.ok, "悬臂 + 三角形荷载求解成功");
    if (!r.ok) return;

    double worstV = 0, worstM = 0, ref = 0;
    for (int k = 0; k <= 12; ++k) {
      const double x = L * k / 12.0;
      const int e = std::min<int>(static_cast<int>(x / (L / n)), n - 1);
      const auto* el = static_cast<const BeamElement*>(m.elements()[static_cast<size_t>(e)].get());
      std::vector<double> u2(12, 0.0), f2;
      for (int i = 0; i < 12; ++i) u2[i] = r.u[el->nodes()[i / 6] * 6 + i % 6];
      el->equivalentLoads(f2);
      Vec12d uu{}, ff{};
      for (int i = 0; i < 12; ++i) { uu[i] = u2[i]; ff[i] = f2[i]; }
      const double xe = L * e / n;
      const auto st = el->atStation(uu, ff, (x - xe) / (L / n));
      const double Va = q0 * (L - x) * (L - x) / (2 * L);
      const double Ma = -q0 * (L - x) * (L - x) * (L - x) / (6 * L);
      worstV = std::max(worstV, std::abs(st.Vy - Va));
      worstM = std::max(worstM, std::abs(st.Mz - Ma));
      ref = std::max(ref, std::abs(Ma));
    }
    std::printf("       三角形荷载：V 最大偏差 %.2e，M 最大偏差 %.2e（相对 %.2e）\n",
                worstV, worstM, worstM / ref);
    check(worstV / (q0 * L / 2.0) < 1e-6, "三角形荷载：剪力图对标解析解");
    check(worstM / ref < 1e-3, "三角形荷载：弯矩图对标解析解（三次曲线）");
  }

  // (b) 段荷载退化为均布时必须与均布解析式逐位一致
  {
    const double L = 4.0, w = 3.75;
    Model m;
    Id a = m.addNode({0, 0, 0}, 0), b = m.addNode({L, 0, 0}, 0);
    auto* be = m.addBeam(a, b, sec, mat, Vec3{0, 0, -1});
    be->addLineLoadSegment(0.0, L, Vec3{0, 0, -w}, Vec3{0, 0, -w});
    m.fixAll(a);
    m.assignDofs();
    StaticAnalysis sa(m);
    const auto r = sa.solve();
    const auto& ef = r.beamForces[0].ef;
    // 悬臂梁【固定端】的截面剪力 = 全部荷载 = wL。
    // 写成 wL/2 是简支梁的支座反力 —— 梁端剪力与支座反力不是一回事。
    checkNear(ef.Vy, w * L, 1e-9, "均布段荷载：固定端剪力 = wL");
    checkNear(ef.Mz, -w * L * L / 2.0, 1e-9, "均布段荷载：固端弯矩 = −wL²/2");
    // 与"旧 setLineLoad 路径"走的是同一段代码，所以不另设期望值，
    // 只核对量级正确 —— 真正的保护是上面这条解析对标。
  }

  // (c) 内力图在单元交接处必须连续（弯矩跳变 = 图错）
  {
    const double L = 8.0, Ls = 4.0, p = 10.0;
    const double a = Ls / 2, bb = L - Ls / 2;
    const int n = 8;
    Model m;
    std::vector<Id> ns;
    for (int i = 0; i <= n; ++i) ns.push_back(m.addNode({L * i / n, 0, 0}, 0));
    auto qf = [&](double x) {
      if (x <= a) return p * x / a;
      if (x >= bb) return p * (L - x) / a;
      return p;
    };
    for (int i = 0; i < n; ++i) {
      const double x0 = L * i / n, x1 = L * (i + 1) / n;
      m.addBeam(ns[i], ns[i + 1], sec, mat, Vec3{0, 0, -1})
          ->addLineLoadSegment(0.0, x1 - x0, Vec3{0, 0, -qf(x0)}, Vec3{0, 0, -qf(x1)});
    }
    m.fixNode(ns[0], true, true, true, true, false, false);
    m.fixNode(ns[n], true, true, true, false, false, false);
    m.assignDofs();
    StaticAnalysis sa(m);
    const auto r = sa.solve();

    double worstSelf = 0, worstJump = 0;
    for (int e = 0; e < n; ++e) {
      const auto* el = static_cast<const BeamElement*>(m.elements()[static_cast<size_t>(e)].get());
      std::vector<double> u2(12, 0.0), f2;
      for (int i = 0; i < 12; ++i) u2[i] = r.u[el->nodes()[i / 6] * 6 + i % 6];
      el->equivalentLoads(f2);
      Vec12d uu{}, ff{};
      for (int i = 0; i < 12; ++i) { uu[i] = u2[i]; ff[i] = f2[i]; }
      const auto s1 = el->atStation(uu, ff, 1.0);
      worstSelf = std::max(worstSelf,
                           std::abs(s1.Mz - r.beamForces[static_cast<size_t>(e)].ef.Mzj));
      if (e + 1 < n) {
        const auto* e2 = static_cast<const BeamElement*>(m.elements()[static_cast<size_t>(e + 1)].get());
        std::vector<double> u3(12, 0.0), f3;
        for (int i = 0; i < 12; ++i) u3[i] = r.u[e2->nodes()[i / 6] * 6 + i % 6];
        e2->equivalentLoads(f3);
        Vec12d uu3{}, ff3{};
        for (int i = 0; i < 12; ++i) { uu3[i] = u3[i]; ff3[i] = f3[i]; }
        worstJump = std::max(worstJump, std::abs(s1.Mz - e2->atStation(uu3, ff3, 0.0).Mz));
      }
    }
    std::printf("       梯形荷载：单元内自洽 %.2e，节点处跨单元跳变 %.2e\n", worstSelf, worstJump);
    check(worstSelf < 1e-9, "内力图在单元端部与 endForces 自洽");
    check(worstJump < 1e-9, "★ 弯矩图在单元交接处连续（线性插值会跳变）");
  }
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  刚性楼板 + 板荷载导荷 验证\n");
  std::printf("======================================================\n");
  testConstruction();
  testRigidBody();
  testNoSpuriousStiffness();
  testEquilibrium();
  testPerStorySwitch();
  testDistribution();
  testEndToEnd();
  testTrapezoidBeamLoad();

  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}
