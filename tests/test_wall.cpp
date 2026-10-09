// =============================================================================
//  tests/test_wall.cpp  ——  剪力墙功能验证（P2：墙进入轴网/网格/脚本/自重/
//                             楼层剪力/UI/测试 整条链路）
//
//  六节：
//    1. 节点重合 —— 墙节点必须与梁/柱/板共用同一份 nodeAt 结果。
//       （否则墙顶水平力传不到楼盖，模型"看起来正常"但结果全错 ——
//       与"板边节点必须与梁节点重合"是同一类陷阱。）
//    2. 悬臂墙解析解 —— 竖直墙当面内悬臂梁，顶部位移对标
//       Timoshenko 解（弯曲 + 剪切修正），验证壳的"面内抗弯"正确。
//    3. 双墙并联 —— 两片相同墙并联，刚度翻倍 ⇒ 同力位移减半；
//       每片墙底水平反力 = P/2（并联按刚度分配）。
//    4. 墙底剪力对账 —— 楼层剪力（柱+竖直墙底的水平合力）
//       与基底水平反力矢量合闭合（呼应内核 T5 的 400 kN 对账实验）。
//    5. 模态 —— 墙自重进质量矩阵（总质量 = ρ·t·b·H/1000 t）；
//       悬臂墙面内弯曲第一阶频率对标连续梁解析解。
//    6. 端到端 —— 跑 examples/frame_wall.yjk 脚本（4 层电梯井框架）：
//       墙壳数量、第 1 层层剪力 = 基底水平反力、剪重比 > 0。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "yjk/analysis/ModalAnalysis.h"
#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/io/ModelScript.h"
#include "yjk/model/GridMesh.h"
#include "yjk/model/Model.h"
#include "yjk/post/PostProcessor.h"

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
    std::printf("  [ OK ] %s  (got %.10g, want %.10g, rel %.2e)\n", what, got, want, rel);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %.10g, want %.10g, rel %.2e)\n", what, got, want, rel);
  }
}

static void checkInt(long long got, long long want, const char* what) {
  if (got == want) { ++g_pass; std::printf("  [ OK ] %s  (got %lld)\n", what, got); }
  else { ++g_fail; std::printf("  [FAIL] %s  (got %lld, want %lld)\n", what, got, want); }
}

// 把轴网变成"纯墙"模型：关掉所有柱 / 梁 / 板开关（墙开关默认关，只开目标墙）
static void wallsOnly(AxisGrid& g) {
  for (Id ix = 0; ix < g.nx(); ++ix)
    for (Id jy = 0; jy < g.ny(); ++jy)
      g.setColumn(ix, jy, false);
  // 注意签名：setBeamAlongX(ix, jy) = X 轴线 ix、Y 向跨度 jy；
  //          setBeamAlongY(jy, ix) = Y 轴线 jy、X 向跨度 ix
  for (Id ix = 0; ix < g.nx(); ++ix)
    for (Id jy = 0; jy + 1 < g.ny(); ++jy)
      g.setBeamAlongX(ix, jy, false);
  for (Id jy = 0; jy < g.ny(); ++jy)
    for (Id ix = 0; ix + 1 < g.nx(); ++ix)
      g.setBeamAlongY(jy, ix, false);
  for (Id ix = 0; ix + 1 < g.nx(); ++ix)
    for (Id jy = 0; jy + 1 < g.ny(); ++jy)
      g.setSlab(ix, jy, false);
  for (Id s = 0; s < g.nStories(); ++s) g.setSlabStory(s, false);
}

// -----------------------------------------------------------------------------
//  1. 节点重合：墙节点必须与梁 / 柱 / 板共用同一份 axis 节点
// -----------------------------------------------------------------------------
static void testCoincidentNodes() {
  std::printf("\n== 1. 墙节点与梁板节点重合 ==\n");

  // 3×3 轴网、2 层，电梯井两片正交墙全开，柱梁板照常生成
  AxisGrid g;
  g.addAxisX(0.0); g.addAxisX(6.0); g.addAxisX(12.0);
  g.addAxisY(0.0); g.addAxisY(6.0); g.addAxisY(12.0);
  g.addStory(0.0); g.addStory(3.6); g.addStory(7.2);
  g.setDivX(0, 2); g.setDivX(1, 2);      // 每跨 X 向剖 2
  g.setDivY(0, 2); g.setDivY(1, 2);      // 每跨 Y 向剖 2
  // 电梯井：沿 X 轴线 1 的两跨墙（ix=1, jy=0/1）+ 沿 Y 轴线 1 的两跨墙（jy=1, ix=0/1）
  g.setWallAlongX(1, 0, true);
  g.setWallAlongX(1, 1, true);
  g.setWallAlongY(1, 0, true);
  g.setWallAlongY(1, 1, true);

  Model m;
  GridMesh::Spec sp;
  sp.columnSection = section::rect(0.5, 0.5);
  sp.beamSection = section::rect(0.3, 0.6);
  sp.slab = ShellProperties{};           // 默认厚 0.2
  sp.wall.thickness = 0.25;
  GridMesh mesh(g, sp);
  const auto rr = mesh.generate(m);
  check(rr.ok, "网格生成成功");
  if (!rr.ok) { std::printf("       %s\n", rr.error.c_str()); return; }
  std::printf("       %s\n", rr.message().c_str());

  // 每层：X 面 2 跨 × 每跨 2 段 = 4 个壳，Y 面同 = 4，共 8；2 层 → 16
  checkInt(rr.walls, 16, "墙壳单元数量：2 层 ×（X 面 4 + Y 面 4）= 16");

  // ★ 关键：墙开关开/关两模型的实建节点总数必须相等 ——
  //   墙没有"额外添加"任何节点，墙顶/墙底与梁板共用 nodeAt 的同一批节点。
  //   （注意：不能直接跟 mesh.nNodes() 比 —— 那是 5×5×3 的网格槽位数，
  //   其中有 4 个槽位未被任何构件引用，节点按需惰性创建。）
  {
    AxisGrid gNowall;                             // 与 g 同轴网、同剖分，仅不开墙
    gNowall.addAxisX(0.0); gNowall.addAxisX(6.0); gNowall.addAxisX(12.0);
    gNowall.addAxisY(0.0); gNowall.addAxisY(6.0); gNowall.addAxisY(12.0);
    gNowall.addStory(0.0); gNowall.addStory(3.6); gNowall.addStory(7.2);
    gNowall.setDivX(0, 2); gNowall.setDivX(1, 2);
    gNowall.setDivY(0, 2); gNowall.setDivY(1, 2);
    Model mNowall;
    GridMesh meshNowall(gNowall, sp);
    const auto rrNo = meshNowall.generate(mNowall);
    check(rrNo.ok, "关墙对照网格生成成功");
    checkInt(rrNo.walls, 0, "关墙对照：墙壳数 = 0");
    checkInt(static_cast<long long>(mNowall.nodeCount()),
             static_cast<long long>(m.nodeCount()),
             "开墙模型节点数 = 关墙模型节点数（墙零新增节点）");
  }

  // 用 mergeCoincidentNodes 复核：不存在任何坐标重合的重复节点
  const MergeReport mr = mergeCoincidentNodes(m);
  checkInt(mr.merged, 0, "无坐标重合的重复节点（生成器已去重）");

  // 抽查墙单元：每片墙都是竖直壳；底层 8 片底面 2 节点在 story 0，
  // 第 2 层 8 片底面 2 节点在 story 1（story 取节点所在楼面层号）
  int wallElems = 0, bottomOk = 0, lev1Ok = 0;
  for (Id e = 0; e < m.elementCount(); ++e) {
    const Element& el = *m.elements()[static_cast<size_t>(e)];
    if (el.type() != ElementType::Shell4) continue;
    const auto* sh = static_cast<const ShellElement*>(&el);
    if (!sh->isWall()) continue;
    ++wallElems;
    const std::vector<Id>& nds = sh->nodes();
    double zmin = 1e300;
    for (Id nd : nds) zmin = std::min(zmin, m.node(nd).r.z);
    int bottom = 0;
    const Id wantStory = zmin < 1e-9 ? 0 : 1;     // 底层墙 → story 0，上层墙 → story 1
    for (Id nd : nds)
      if (std::abs(m.node(nd).r.z - zmin) < 1e-9 &&
          m.node(nd).story == wantStory) ++bottom;
    if (bottom == 2) {
      if (wantStory == 0) ++bottomOk; else ++lev1Ok;
    }
  }
  checkInt(wallElems, 16, "墙单元总数 16");
  checkInt(bottomOk, 8, "底层 8 片墙底面 2 节点都在 story 0");
  checkInt(lev1Ok, 8, "第 2 层 8 片墙底面 2 节点都在 story 1");

  // 满轴线墙开关默认关 —— 不设即无墙（回归保护：不影响既有框架模型）
  AxisGrid g0;
  g0.addAxisX(0.0); g0.addAxisX(6.0);
  g0.addAxisY(0.0); g0.addAxisY(6.0);
  g0.addStory(0.0); g0.addStory(3.6);
  Model m0;
  GridMesh::Spec sp0;
  GridMesh mesh0(g0, sp0);
  const auto rr0 = mesh0.generate(m0);
  checkInt(rr0.walls, 0, "未设任何墙开关 ⇒ 墙壳数 = 0（默认关）");
}

// -----------------------------------------------------------------------------
//  2. 悬臂墙解析解：竖直墙当面内悬臂梁
//
//  一片"沿 X 的墙"（宽 b、高 H、厚 t），底部全固定，顶部沿 +X 施加水平力 P。
//  墙在 X-Z 平面内做悬臂梁弯曲（面内），跨越来越多层（story）细分。
//  对标 Timoshenko 解：δ = P·H³/(3EI) + P·H/(κGA)
//    I = t·b³/12（绕 Y 轴，即墙平面法向轴）
//    A = b·t，κ = 5/6，G = E/(2(1+ν))
// -----------------------------------------------------------------------------
static void testCantileverWall() {
  std::printf("\n== 2. 悬臂墙解析解（面内弯曲 + 剪切） ==\n");

  const double b = 3.0, H = 3.6, t = 0.25;
  const double E = 3.0e7, nu = 0.2;          // kPa
  const double P = 100.0;                    // kN
  const int nSeg = 16;                       // 沿高度分段

  const double I = t * b * b * b / 12.0;
  const double A = b * t;
  const double G = E / (2.0 * (1.0 + nu));
  const double dFlex = P * H * H * H / (3.0 * E * I);
  const double dShear = P * H / ((5.0 / 6.0) * G * A);
  const double want = dFlex + dShear;
  std::printf("       解析：δ弯曲 = %.6e m，δ剪切 = %.6e m，δ合计 = %.6e m\n",
              dFlex, dShear, want);

  AxisGrid g;
  g.addAxisX(0.0);                             // 悬臂墙所在轴线（y = 0）
  g.addAxisX(100.0);                           // 哑轴线：仅满足轴网 ≥2×2 约束，不设构件
  g.addAxisY(0.0); g.addAxisY(b);              // 宽 b
  for (int k = 0; k <= nSeg; ++k) g.addStory(H * k / nSeg);
  wallsOnly(g);                                // 纯墙模型：关柱梁板
  g.setWallAlongX(0, 0, true);                 // 开目标墙

  Model m;
  GridMesh::Spec sp;
  sp.wall.thickness = t;
  sp.wall.E = E;
  sp.wall.nu = nu;
  GridMesh mesh(g, sp);
  const auto rr = mesh.generate(m);
  check(rr.ok, "悬臂墙网格生成成功");
  if (!rr.ok) { std::printf("       %s\n", rr.error.c_str()); return; }
  checkInt(rr.walls, nSeg, "悬臂墙剖分为 nSeg 个壳单元");

  // 固定底面两个角点节点（z = 0 层）
  const Id b0 = mesh.nodeAtAxis(0, 0, 0);
  const Id b1 = mesh.nodeAtAxis(0, 1, 0);
  const Id t0 = mesh.nodeAtAxis(0, 0, nSeg);
  const Id t1 = mesh.nodeAtAxis(0, 1, nSeg);
  m.fixAll(b0);
  m.fixAll(b1);
  m.addNodeForce(t0, Vec3{P / 2, 0, 0});
  m.addNodeForce(t1, Vec3{P / 2, 0, 0});

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "悬臂墙求解成功");
  if (!r.ok) { std::printf("       %s\n", r.message.c_str()); return; }

  const double got = r.u[static_cast<size_t>(t0) * 6 + 0];
  std::printf("       数值：δ顶 = %.6e m\n", got);

  // 【容差 0.09 的由来】nSeg=16 单列（宽向 1 单元）Q4 膜在面内弯曲下
  // 偏刚约 8.3%，两个来源叠加：
  //   · 单元级：面内弯曲位移场 ux∝z² 双线性膜表达不了，宽向 1 列时
  //     弯曲变形受限（上刚）。数值实验显示 nW=1 单列恒比解析刚 ~8%，
  //     而 nW=4~8 中心点可收到 ~3.5% 以内 —— 属 Q4 膜固有插值误差，
  //     不是壳/板叠加引入。
  //   · 理论级：本墙 H/b = 1.2 属深墙，Timoshenko 细梁解（κ=5/6、
  //     平截面假定）对深墙本身就有几 % 偏差。
  // 故本断言验证"面内抗弯机制正确、量级闭合"（<10%），
  // 精确收敛性由下面独立的 16 vs 32 段自检保证。
  checkNear(got, want, 0.09, "墙顶水平位移对标 Timoshenko 悬臂解（含剪切，±9%）");

  // 网格收敛自检：nSeg 从 16 → 32 位移变化 < 0.5%。
  // 若单元实现有缺陷（沙漏/锁定/装配错位），细分必然大幅漂移；
  // 该断言防止"解析对照蒙混过关、网格却发散"的假通过。
  {
    const int nSeg2 = 32;
    AxisGrid g2;
    g2.addAxisX(0.0); g2.addAxisX(100.0);
    g2.addAxisY(0.0); g2.addAxisY(b);
    for (int k = 0; k <= nSeg2; ++k) g2.addStory(H * k / nSeg2);
    wallsOnly(g2);
    g2.setWallAlongX(0, 0, true);
    Model m2;
    GridMesh::Spec sp2;
    sp2.wall.thickness = t;
    sp2.wall.E = E;
    sp2.wall.nu = nu;
    GridMesh mesh2(g2, sp2);
    mesh2.generate(m2);
    const Id b02 = mesh2.nodeAtAxis(0, 0, 0), b12 = mesh2.nodeAtAxis(0, 1, 0);
    const Id t02 = mesh2.nodeAtAxis(0, 0, nSeg2), t12 = mesh2.nodeAtAxis(0, 1, nSeg2);
    m2.fixAll(b02); m2.fixAll(b12);
    m2.addNodeForce(t02, Vec3{P / 2, 0, 0}); m2.addNodeForce(t12, Vec3{P / 2, 0, 0});
    m2.assignDofs();
    StaticAnalysis sa2(m2);
    const auto r2 = sa2.solve();
    check(r2.ok, "收敛自检：nSeg=32 求解成功");
    if (r2.ok) {
      const double got2 = r2.u[static_cast<size_t>(t02) * 6 + 0];
      std::printf("       nSeg=16: %.6e，nSeg=32: %.6e，变化 %.4f%%\n",
                  got, got2, std::abs(got2 - got) / got * 100.0);
      checkNear(got2, got, 0.005,
                "悬臂墙网格收敛：nSeg 16→32 位移变化 < 0.5%");
    }
  }
}

// -----------------------------------------------------------------------------
//  3. 双墙并联：刚度翻倍、反力平分
// -----------------------------------------------------------------------------
static void testParallelWalls() {
  std::printf("\n== 3. 双墙并联（刚度翻倍 / 反力平分） ==\n");

  const double b = 3.0, H = 3.6, t = 0.25, gap = 2.0;
  const double P = 100.0;                     // 总水平力
  const int nSeg = 8;

  auto buildSingle = [&](int nWall) -> double {
    AxisGrid g;
    g.addAxisX(0.0); g.addAxisX(gap);         // 两片墙的轴线（y = 0 / gap）
    g.addAxisY(0.0); g.addAxisY(b);
    for (int k = 0; k <= nSeg; ++k) g.addStory(H * k / nSeg);
    wallsOnly(g);                              // 纯墙模型：关柱梁板
    g.setWallAlongX(0, 0, true);
    if (nWall == 2) g.setWallAlongX(1, 0, true);
    Model m;
    GridMesh::Spec sp;
    sp.wall.thickness = t;
    GridMesh mesh(g, sp);
    const auto rr = mesh.generate(m);
    if (!rr.ok) return -1.0;
    // 固定所有底面角点节点
    for (int iy = 0; iy < nWall; ++iy) {
      m.fixAll(mesh.nodeAtAxis(iy, 0, 0));
      m.fixAll(mesh.nodeAtAxis(iy, 1, 0));
    }
    // 顶层两个节点各 P/2（两片墙顶各 2 节点 → 每片 P/2）
    for (int iy = 0; iy < nWall; ++iy) {
      m.addNodeForce(mesh.nodeAtAxis(iy, 0, nSeg), Vec3{P / (2 * nWall), 0, 0});
      m.addNodeForce(mesh.nodeAtAxis(iy, 1, nSeg), Vec3{P / (2 * nWall), 0, 0});
    }
    m.assignDofs();
    StaticAnalysis sa(m);
    const auto r = sa.solve();
    if (!r.ok) return -1.0;
    // 顶部位移取墙 0 的顶节点
    return r.u[static_cast<size_t>(mesh.nodeAtAxis(0, 0, nSeg)) * 6 + 0];
  };

  const double d1 = buildSingle(1);
  const double d2 = buildSingle(2);
  std::printf("       单墙 δ = %.6e m，双墙 δ = %.6e m\n", d1, d2);
  check(d1 > 0 && d2 > 0, "两种模型求解成功且位移为正");
  checkNear(d2 / d1, 0.5, 0.02, "双墙并联 ⇒ 总刚度翻倍 ⇒ 位移减半");

  // 反力平分：重新建双墙模型，核对每片墙底水平反力 = P/2
  {
    AxisGrid g;
    g.addAxisX(0.0); g.addAxisX(gap);
    g.addAxisY(0.0); g.addAxisY(b);
    for (int k = 0; k <= nSeg; ++k) g.addStory(H * k / nSeg);
    wallsOnly(g);
    g.setWallAlongX(0, 0, true);
    g.setWallAlongX(1, 0, true);
    Model m;
    GridMesh::Spec sp;
    sp.wall.thickness = t;
    GridMesh mesh(g, sp);
    mesh.generate(m);
    for (int iy = 0; iy < 2; ++iy) {
      m.fixAll(mesh.nodeAtAxis(iy, 0, 0));
      m.fixAll(mesh.nodeAtAxis(iy, 1, 0));
    }
    for (int iy = 0; iy < 2; ++iy) {
      m.addNodeForce(mesh.nodeAtAxis(iy, 0, nSeg), Vec3{P / 4, 0, 0});
      m.addNodeForce(mesh.nodeAtAxis(iy, 1, nSeg), Vec3{P / 4, 0, 0});
    }
    m.assignDofs();
    StaticAnalysis sa(m);
    const auto r = sa.solve();
    check(r.ok, "双墙模型（反力校验）求解成功");
    if (!r.ok) return;
    double rx0 = 0, rx1 = 0;
    for (Id nd : {mesh.nodeAtAxis(0, 0, 0), mesh.nodeAtAxis(0, 1, 0)})
      rx0 += r.reaction[static_cast<size_t>(nd) * 6 + 0];
    for (Id nd : {mesh.nodeAtAxis(1, 0, 0), mesh.nodeAtAxis(1, 1, 0)})
      rx1 += r.reaction[static_cast<size_t>(nd) * 6 + 0];
    std::printf("       墙0 底反力 = %+.6f kN，墙1 底反力 = %+.6f kN（总 %.0f kN）\n",
                rx0, rx1, P);
    checkNear(rx0, -P / 2, 1e-9, "墙0：底水平反力 = −P/2（按刚度平分）");
    checkNear(rx1, -P / 2, 1e-9, "墙1：底水平反力 = −P/2（按刚度平分）");
  }
}

// -----------------------------------------------------------------------------
//  4. 墙底剪力对账：楼层剪力 = 基底水平反力矢量合
//
//  复刻 examples/mini_wall.yjk：2×2 单层、两片正交墙、顶层 +Y 100 kN/节点
//  （共 400 kN）。第 1 层层剪力（柱 + 竖直墙底）必须闭合到 400.00 ——
//  这就是内核 T5 修复"墙底水平内力差一个负号"那条 bug 的实验。
// -----------------------------------------------------------------------------
static void testBaseShear() {
  std::printf("\n== 4. 墙底剪力对账（层剪力 = 基底水平反力） ==\n");

  AxisGrid g;
  g.addAxisX(0.0); g.addAxisX(6.0);
  g.addAxisY(0.0); g.addAxisY(6.0);
  g.addStory(0.0); g.addStory(3.6);
  g.setWallAlongX(0, 0, true);
  g.setWallAlongY(0, 0, true);
  Model m;
  GridMesh::Spec sp;
  sp.columnSection = section::rect(0.5, 0.5);
  sp.beamSection = section::rect(0.3, 0.6);
  sp.buildSlabs = false;
  sp.wall.thickness = 0.25;
  GridMesh mesh(g, sp);
  const auto rr = mesh.generate(m);
  check(rr.ok, "mini 模型网格生成成功");
  if (!rr.ok) { std::printf("       %s\n", rr.error.c_str()); return; }

  for (Id n : mesh.baseNodes()) m.fixAll(n);

  // 顶层每节点 +Y 100 kN（顶层 4 节点 → 400 kN）
  for (Id jx = 0; jx < g.nx(); ++jx)
    for (Id jy = 0; jy < g.ny(); ++jy)
      m.addNodeForce(mesh.nodeAtAxis(jx, jy, 1), Vec3{0, 100, 0});

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "mini 模型求解成功");
  if (!r.ok) { std::printf("       %s\n", r.message.c_str()); return; }

  // 基底水平反力矢量合
  double rx = 0, ry = 0;
  for (Id n : mesh.baseNodes()) {
    rx += r.reaction[static_cast<size_t>(n) * 6 + 0];
    ry += r.reaction[static_cast<size_t>(n) * 6 + 1];
  }
  const double base = std::hypot(rx, ry);
  std::printf("       基底水平反力 = Σ(Rx,Ry) = (%+.4f, %+.4f) ⇒ %.4f kN\n",
              rx, ry, base);
  checkNear(base, 400.0, 1e-9, "基底水平反力矢量合 = 400 kN");

  post::PostProcessor pp(m, r);
  pp.compute();
  check(pp.stories().size() >= 2, "楼层结果含基础层 + 第 1 层");
  if (pp.stories().size() >= 2) {
    // stories()[0] 是基础层（恒 0），第 1 层层剪力是 stories()[1]
    const double v = pp.stories()[1].shear;
    std::printf("       第 1 层层剪力 = %.4f kN\n", v);
    checkNear(v, 400.0, 1e-9, "★ 第 1 层层剪力（柱+墙）= 400 kN = 基底反力");
  }
}

// -----------------------------------------------------------------------------
//  5. 模态：墙自重进质量矩阵，悬臂墙基频对标解析解
// -----------------------------------------------------------------------------
static void testWallModal() {
  std::printf("\n== 5. 墙模态（质量与基频） ==\n");

  const double b = 3.0, H = 3.6, t = 0.25;
  const double E = 3.0e7, nu = 0.2;
  const double rho = 2500.0;                  // kg/m³
  const int nSeg = 12;

  const double I = t * b * b * b / 12.0;
  const double A = b * t;
  const double mTot = rho * A * H / 1000.0;   // t（吨）

  // 悬臂墙底只固定 2 个角点（见下方 fixAll），这两个节点的集中质量
  // 不进入"可动"质量统计（totalMassX = rᵀ·M·r，固定 dof 置零）。
  // 每个壳单元质量均分到 4 节点：mNode = mTot / nSeg / 4。
  const double mFixed2 = 2.0 * (mTot / nSeg / 4.0);   // 2 个固定底角节点质量

  AxisGrid g;
  g.addAxisX(0.0);
  g.addAxisX(100.0);                          // 哑轴线：仅满足轴网 ≥2×2 约束
  g.addAxisY(0.0); g.addAxisY(b);
  for (int k = 0; k <= nSeg; ++k) g.addStory(H * k / nSeg);
  wallsOnly(g);
  g.setWallAlongX(0, 0, true);
  Model m;
  GridMesh::Spec sp;
  sp.wall.thickness = t;
  sp.wall.E = E;
  sp.wall.nu = nu;
  sp.wall.density = rho;
  GridMesh mesh(g, sp);
  const auto rr = mesh.generate(m);
  check(rr.ok, "墙模态网格生成成功");
  if (!rr.ok) { std::printf("       %s\n", rr.error.c_str()); return; }
  checkInt(rr.walls, nSeg, "模态墙剖分为 nSeg 个壳单元");
  std::fprintf(stderr, "    [dbg] modal: nNodes=%lld b0=%lld b1=%lld\n",
               static_cast<long long>(m.nodeCount()),
               static_cast<long long>(mesh.nodeAtAxis(0, 0, 0)),
               static_cast<long long>(mesh.nodeAtAxis(0, 1, 0)));

  m.fixAll(mesh.nodeAtAxis(0, 0, 0));
  m.fixAll(mesh.nodeAtAxis(0, 1, 0));
  m.assignDofs();
  std::fprintf(stderr, "    [dbg] modal: freeDof=%lld\n",
               static_cast<long long>(m.freeDofCount()));

  ModalAnalysis ma(m);
  ModalResult mr = ma.solve(ModalAnalysis::Options{0, 1e-9, 300, 0});
  check(mr.ok, "墙模态求解成功");
  if (!mr.ok) { std::printf("       %s\n", mr.message.c_str()); return; }
  std::printf("       模态数 %d，参与质量 X = %.8f t（解析全质量 %.8f − 固定角点 %.6f）\n",
              mr.nmodes, mr.totalMassX, mTot, mFixed2);
  checkNear(mr.totalMassX, mTot - mFixed2, 1e-9,
            "X 向参与质量 = 墙自重/9.81 − 固定角点集中质量（质量只来自墙）");

  // 悬臂梁（面内弯曲，忽略剪切与转动惯量）：ω₁ = (1.875)²·√(EI/(ρA·H⁴))
  const double EI = E * I;                    // kPa·m⁴ = kN·m²（N·m² 同量纲）
  const double rhoA = rho * A;                // kg/m
  const double wAnal = std::pow(1.875, 2.0) * std::sqrt(EI * 1000.0 / (rhoA * std::pow(H, 4)));
  // 注意单位：EI 以 kN·m²（=1000 N·m²）计 → 乘 1000 转 N·m²

  // 【关键：不能拿 omega[0] 对标面内弯曲】
  // 本墙底只固定两个角点、宽 3m 高 3.6m，最低频模态是面外（板）弯曲
  // （ω≈68 rad/s，面外刚度远小于面内）；面内弯曲基频在
  // "X 向参与质量最大"的模态（nSeg=12 时 ~552 rad/s）。
  int imax = 0;
  for (int im = 1; im < mr.nmodes; ++im)
    if (mr.mstarX[static_cast<size_t>(im)] > mr.mstarX[static_cast<size_t>(imax)])
      imax = im;
  const double wNum = mr.omega[static_cast<size_t>(imax)];
  std::printf("       面内弯曲基频 ω₁ = %.4f rad/s（模态 %d，解析 %.4f）\n",
              wNum, imax + 1, wAnal);
  // 深墙（H/b=1.2）剪切与转动惯量效应显著，欧拉-伯努利解析值偏高；
  // nSeg=12 集中质量亦略偏柔。取方向性判据 [0.6, 1.1] × 解析：
  check(wNum > 0.6 * wAnal && wNum < 1.1 * wAnal,
        "悬臂墙面内弯曲基频落在解析解 ±40% 内（深墙剪切软化 + 集中质量离散）");
}

// -----------------------------------------------------------------------------
//  6. 端到端：examples/frame_wall.yjk
// -----------------------------------------------------------------------------
static void testEndToEnd(const char* scriptPath) {
  std::printf("\n== 6. 端到端：frame_wall.yjk ==\n");

  io::ModelScript sc;
  check(sc.parseFile(scriptPath), "解析 frame_wall.yjk 成功");
  if (!sc.errors().empty() && sc.errors().size() > 0) {
    // parseFile 失败时 errors 非空；这里只做非致命诊断
    if (!sc.errors().empty())
      std::printf("       [解析提示] %s\n", sc.errors()[0].c_str());
  }
  Model m;
  std::string info;
  check(sc.build(m, &info), "按脚本建模型成功");
  if (!info.empty()) std::printf("       %s\n", info.c_str());

  // 墙壳数量：电梯井 = 3 层 ×（沿 X 轴 1 两侧壁 2 跨×2 剖分 + 沿 Y 轴 1 同）？
  // 脚本：wallx 1 0 / 1 1；wally 1 0 / 1 1，divX/divY = 2，层数 = 3（4 标高）：
  //   每层 X 面 4 壳 + Y 面 4 壳 = 8，×3 层 = 24
  int nWall = 0;
  for (Id e = 0; e < m.elementCount(); ++e) {
    const Element& el = *m.elements()[static_cast<size_t>(e)];
    if (el.type() != ElementType::Shell4) continue;
    if (static_cast<const ShellElement&>(el).isWall()) ++nWall;
  }
  checkInt(nWall, 24, "墙壳单元数 = 24（3 层 × 每层 8）");

  // 求解 + 后处理
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "框架+墙求解成功");
  if (!r.ok) { std::printf("       %s\n", r.message.c_str()); return; }

  post::PostProcessor pp(m, r);
  pp.compute();
  const auto& env = pp.envelope();
  const auto& stories = pp.stories();
  check(!stories.empty(), "有楼层结果");
  if (stories.empty()) return;

  // 基底水平反力矢量合（应 ≈ 顶层 Y 向风载 15 kN × 顶层节点数）
  double rx = 0, ry = 0;
  for (Id i = 0; i < m.nodeCount(); ++i) {
    rx += r.reaction[static_cast<size_t>(i) * 6 + 0];
    ry += r.reaction[static_cast<size_t>(i) * 6 + 1];
  }
  const double base = std::hypot(rx, ry);
  // stories()[0] 是基础层（恒 0），取第 1 层 stories()[1]
  const double v1 = stories[1].shear;
  std::printf("       基底水平反力 = %.4f kN，第 1 层层剪力 = %.4f kN\n", base, v1);
  checkNear(v1, base, 1e-6, "第 1 层层剪力 = 基底水平反力矢量合（T5 结论）");

  std::printf("       总重 W = %.2f kN，剪重比 V₁/W = %.4f\n",
              env.totalWeight, env.shearWeightRatio);
  check(env.shearWeightRatio > 0.0, "剪重比 > 0（墙自重已入分母）");
  check(env.shearWeightRatio < 0.5, "剪重比量级合理（< 0.5）");
}

int main(int argc, char** argv) {
  std::printf("======================================================\n");
  std::printf("  yjk_like  剪力墙功能验证\n");
  std::printf("======================================================\n");
  testCoincidentNodes();
  testCantileverWall();
  testParallelWalls();
  testBaseShear();
  testWallModal();
  if (argc > 1) testEndToEnd(argv[1]);
  else std::printf("\n== 6. 端到端：跳过（未给脚本路径） ==\n");

  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}