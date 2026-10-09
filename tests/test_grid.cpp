// =============================================================================
//  tests/test_grid.cpp  ——  轴网建模（GridMesh）验证
//
//  验证原则：
//    ① 拓扑正确：节点/柱/梁/板数量可手算核对
//    ② 【最关键】板边节点与梁节点必须重合 —— 否则荷载传不到梁上且不报错
//    ③ 力学正确：生成的模型能求解，且总反力 = 总荷载
//    ④ 参数错误必须被拒绝，不能静默生成畸形模型
// =============================================================================
#include <cmath>
#include <cstdio>
#include <set>
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
  if (got == want) {
    ++g_pass;
    std::printf("  [ OK ] %s  (got %lld)\n", what, got);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %lld, want %lld)\n", what, got, want);
  }
}

// 统计每个节点被哪些类型的单元使用
static void usageByType(const Model& m, std::vector<char>& byBeam, std::vector<char>& byShell) {
  byBeam.assign(static_cast<size_t>(m.nodeCount()), 0);
  byShell.assign(static_cast<size_t>(m.nodeCount()), 0);
  for (const auto& e : m.elements()) {
    for (Id n : e->nodes()) {
      if (e->type() == ElementType::Beam3D) byBeam[static_cast<size_t>(n)] = 1;
      if (e->type() == ElementType::Shell4) byShell[static_cast<size_t>(n)] = 1;
    }
  }
}

static GridMesh::Spec defaultSpec() {
  GridMesh::Spec s;
  s.columnSection = section::rect(0.5, 0.5);
  s.beamSection = section::rect(0.3, 0.6);
  s.columnMaterial = Material::concreteC(30);
  s.beamMaterial = Material::concreteC(30);
  s.slab.thickness = 0.12;
  s.slab.E = 3.0e7;
  s.slab.nu = 0.2;
  return s;
}

// =============================================================================
// 1. 轴网拓扑：数量可手算核对
// =============================================================================
static void testTopology() {
  std::printf("\n== 1. 轴网拓扑（数量手算核对） ==\n");

  // 3 条 X 向轴线（y=0,6,12）× 3 条 Y 向轴线（x=0,5,10）
  // 2 个楼层：z=0（地面）、z=4（屋面）
  // 每跨板剖分 2×2
  AxisGrid g;
  for (double y : {0.0, 6.0, 12.0}) g.addAxisX(y);
  for (double x : {0.0, 5.0, 10.0}) g.addAxisY(x);
  g.addStory(0.0);
  g.addStory(4.0);
  g.setDivX(0, 2);
  g.setDivX(1, 2);
  g.setDivY(0, 2);
  g.setDivY(1, 2);

  Model m;
  GridMesh mesh(g, defaultSpec());
  const auto r = mesh.generate(m);
  check(r.ok, "网格生成成功");
  if (!r.ok) { std::printf("       %s\n", r.error.c_str()); return; }
  std::printf("       %s\n", r.message().c_str());

  // ---- 手算 ----
  // 细分后 X 向点数 = 1 + 2 + 2 = 5；Y 向点数 = 5
  // 节点：地面层 5×5 = 25？不对 —— 地面层没有板，
  //   只有柱脚(9) + 地面层梁(沿 X 的 3 条线 × (2+2) 段 = 3×4=12 根梁? )
  //
  // 逐项手算：
  //   柱：3×3 = 9 个交点，每个 1 层 → 9 根
  //   梁（每层）：
  //     沿 X：3 条 X 向轴线 × 2 跨 × 2 段 = 12 根
  //     沿 Y：3 条 Y 向轴线 × 2 跨 × 2 段 = 12 根
  //     每层 24 根；地面层与屋面层各 24 → 48 根
  //   板：仅屋面层（底层默认无板）：2×2 板格 × (2×2) 单元 = 16 个
  checkInt(r.columns, 9, "柱数 = 交点数 3×3 × 层数 1 = 9");
  checkInt(r.beams, 48, "梁数 = (沿X 3×2×2 + 沿Y 3×2×2) × 2 层 = 48");
  checkInt(r.slabs, 16, "板数 = 2×2 板格 × 2×2 剖分 × 1 层 = 16");

  // 轴线交点坐标
  const Id n00 = mesh.nodeAtAxis(0, 0, 0);
  check(n00 >= 0, "nodeAtAxis 能查到交点节点");
  if (n00 >= 0) {
    checkNear(m.node(n00).r.x, 0.0, 1e-15, "交点 (0,0) 的 x 坐标");
    checkNear(m.node(n00).r.y, 0.0, 1e-15, "交点 (0,0) 的 y 坐标");
    checkNear(m.node(n00).r.z, 0.0, 1e-15, "交点 (0,0) 底层 z 坐标");
  }
  const Id n11 = mesh.nodeAtAxis(1, 1, 1);
  if (n11 >= 0) {
    checkNear(m.node(n11).r.x, 5.0, 1e-15, "交点 (1,1) 的 x 坐标");
    checkNear(m.node(n11).r.y, 6.0, 1e-15, "交点 (1,1) 的 y 坐标");
    checkNear(m.node(n11).r.z, 4.0, 1e-15, "交点 (1,1) 屋面层 z 坐标");
  }
}

// =============================================================================
// 2. 【最关键】板边节点与梁节点重合
//
//    这是轴网建模唯一"错了不报错"的地方：板剖分后板边中间节点若没有梁，
//    板边在两柱之间就是悬空的，荷载传不到梁上，而模型看起来完全正常。
// =============================================================================
static void testSlabBeamSharedNodes() {
  std::printf("\n== 2. 板边节点与梁节点重合（最关键） ==\n");

  AxisGrid g;
  for (double y : {0.0, 6.0, 12.0}) g.addAxisX(y);
  for (double x : {0.0, 5.0, 10.0}) g.addAxisY(x);
  g.addStory(0.0);
  g.addStory(4.0);
  g.setDivX(0, 3);
  g.setDivX(1, 3);
  g.setDivY(0, 3);
  g.setDivY(1, 3);

  Model m;
  GridMesh mesh(g, defaultSpec());
  const auto r = mesh.generate(m);
  check(r.ok, "网格生成成功");
  if (!r.ok) return;

  std::vector<char> byBeam, byShell;
  usageByType(m, byBeam, byShell);

  // 找出"被壳单元使用、但不被任何梁使用"的节点
  int shellOnly = 0;
  for (Id i = 0; i < m.nodeCount(); ++i)
    if (byShell[i] && !byBeam[i]) ++shellOnly;
  std::printf("       只被板使用、不被梁使用的节点数 = %d\n", shellOnly);

  // 板内部节点（不在任何轴线/细分边界上）确实只连板 —— 这是正常的。
  // 判据应该是【板边界上的节点必须也被梁连接】。
  //
  // 板边界 = 落在 X 向轴线上的（q == yAxisIndex(ix)）或 Y 向轴线上的（p == xAxisIndex(jy)），
  // 或者更严格：直接检查每个壳单元的 4 个节点中，
  // 位于【板格边界】上的那些节点是否被梁连接。
  //
  // 最简洁且严格的做法：沿 X 向轴线的梁必须是"逐段"的 ——
  // 一跨 divX=3 就必须生成 3 根梁，而不是 1 根通长梁。
  // 用梁长验证：跨长 5.0，剖分 3 → 每段 5/3 ≈ 1.667
  // 水平梁长度必须等于"跨长/剖分数"，而不是整跨：
  //   沿 X 的跨长 5.0，剖分 3 → 段长 5/3 ≈ 1.6667
  //   沿 Y 的跨长 6.0，剖分 3 → 段长 6/3 = 2.0
  // 【注意两个方向的跨长不同，别只按一个方向算 —— 这里错过一次】
  int nSegX = 0, nSegY = 0, nFullSpan = 0;
  for (const auto& e : m.elements()) {
    if (e->type() != ElementType::Beam3D) continue;
    const double L = e->length();
    if (std::abs(L - 5.0 / 3.0) < 1e-9) ++nSegX;
    else if (std::abs(L - 2.0) < 1e-9) ++nSegY;
    else if (std::abs(L - 5.0) < 1e-9 || std::abs(L - 6.0) < 1e-9) ++nFullSpan;
  }
  // 沿 X：3 条 X 向轴线 × 2 跨 × 3 段 × 2 层 = 36；沿 Y 同理 36
  checkInt(nSegX, 36, "沿 X 的梁段数 = 36（3 线 × 2 跨 × 3 段 × 2 层）");
  checkInt(nSegY, 36, "沿 Y 的梁段数 = 36");
  checkInt(nFullSpan, 0, "不存在未被打断的通长梁（0 根 = 全部逐段生成）");

  // 再用"悬空板边"直接判据复核：
  // 对每个壳单元，其 4 个节点若有落在网格边界（p 是轴线索引 或 q 是轴线索引）上的，
  // 该节点必须也被梁使用。
  g.buildPoints();
  std::set<Id> axisP, axisQ;
  for (Id jy = 0; jy < g.ny(); ++jy) axisP.insert(g.xAxisIndex(jy));
  for (Id ix = 0; ix < g.nx(); ++ix) axisQ.insert(g.yAxisIndex(ix));

  // 反查节点 → (p,q)
  std::map<Id, std::pair<Id, Id>> pqOf;
  for (Id p = 0; p < g.nXPoints(); ++p)
    for (Id q = 0; q < g.nYPoints(); ++q)
      for (Id s = 0; s < g.nStories(); ++s) {
        const Id n = mesh.nodeAt(p, q, s);
        if (n >= 0) pqOf[n] = {p, q};
      }

  int danglingEdge = 0;
  for (const auto& e : m.elements()) {
    if (e->type() != ElementType::Shell4) continue;
    for (Id n : e->nodes()) {
      auto it = pqOf.find(n);
      if (it == pqOf.end()) continue;
      const Id p = it->second.first, q = it->second.second;
      const bool onBoundary = axisP.count(p) || axisQ.count(q);
      if (onBoundary && !byBeam[static_cast<size_t>(n)]) ++danglingEdge;
    }
  }
  checkInt(danglingEdge, 0, "板边界上不存在只连板、不连梁的节点");
}

// =============================================================================
// 3. 力学验证：生成的模型能算，且总反力 = 总荷载
// =============================================================================
static void testMechanics() {
  std::printf("\n== 3a. 力学验证：荷载沿轴网传到柱脚 ==\n");

  // 3×3 轴线 → 2×2 板格，每格 1 个壳单元（div=1），2 层。
  //
  // 【为什么这里用 div=1】壳单元当前只有【膜】刚度，板弯曲（MITC4）尚未实现，
  // 所以壳单元内部节点的 uz 没有刚度。div=1 时板的 4 个角都落在轴线交点上、
  // 由柱与梁提供竖向刚度，模型可解。div>1 的情形见 3b。
  AxisGrid g;
  for (double y : {0.0, 6.0, 12.0}) g.addAxisX(y);
  for (double x : {0.0, 6.0, 12.0}) g.addAxisY(x);
  g.addStory(0.0);
  g.addStory(3.5);

  Model m;
  GridMesh mesh(g, defaultSpec());
  const auto r = mesh.generate(m);
  check(r.ok, "网格生成成功");
  if (!r.ok) { std::printf("       %s\n", r.error.c_str()); return; }
  std::printf("       %s\n", r.message().c_str());

  const auto bases = mesh.baseNodes();
  checkInt(static_cast<long long>(bases.size()), 9, "柱脚节点数 = 9");
  for (Id n : bases) m.fixAll(n);

  const double q = 5.0;                 // kPa，向下
  int nSlab = 0;
  for (auto& e : m.elements()) {
    if (e->type() != ElementType::Shell4) continue;
    static_cast<ShellElement*>(e.get())->setTransversePressure(-q);
    ++nSlab;
  }
  checkInt(nSlab, 4, "板单元数 = 2×2 = 4");

  m.assignDofs();
  const auto diag = m.check();
  check(diag.ok(), "模型校核通过");
  for (const auto& e : diag.errors) std::printf("       [错误] %s\n", e.c_str());

  StaticAnalysis sa(m);
  const auto res = sa.solve();
  check(res.ok, "静力求解成功");
  if (!res.ok) { std::printf("       %s\n", res.message.c_str()); return; }

  // 总荷载 = q × 平面面积 = 5 × 12×12 = 720 kN
  const double totalLoad = q * 12.0 * 12.0;
  double sumRz = 0.0;
  for (Id n : bases) sumRz += res.reaction[static_cast<size_t>(n) * 6 + 2];
  std::printf("       板面总荷载 = %.2f kN，柱脚竖向反力合计 = %.6f kN\n", totalLoad, sumRz);
  checkNear(sumRz, totalLoad, 1e-9, "柱脚竖向反力合计 = 板面总荷载");

  check(res.residual < 1e-9, "解残差 < 1e-9");
  std::printf("       自由度 %d，nnz=%zu，nnzL=%zu，残差 %.2e\n",
              res.ndof, res.nnz, res.nnzL, res.residual);

  // 中柱受力应大于角柱（板格把荷载往中间聚）
  const Id corner = mesh.nodeAtAxis(0, 0, 0);
  const Id center = mesh.nodeAtAxis(1, 1, 0);
  if (corner >= 0 && center >= 0) {
    const double rc = res.reaction[static_cast<size_t>(corner) * 6 + 2];
    const double rm = res.reaction[static_cast<size_t>(center) * 6 + 2];
    std::printf("       角柱 %.4f kN，中柱 %.4f kN\n", rc, rm);
    check(rm > rc, "中柱反力大于角柱（荷载向中部聚集）");
  }
}

// -----------------------------------------------------------------------------
// 3b. 楼板剖分 —— MITC4 板弯曲到位后的验证
//
//  【这一节原来是"已知缺口"的断言】
//  在板弯曲（MITC4）实现之前，这里断言的是"板剖分后内部节点的 uz 无刚度、
//  模型被判为不可解（检出 4 个自由度）"，用来把缺口固定住、不让它被静默放过。
//  MITC4 到位后该断言【必须替换】—— 留着它等于把"功能缺失"当成正确行为。
//  断言要跟着功能往前走，不能当纪念品。
//
//  现在这里断言的是真实物理：
//    ① 剖分后的板可解（不再有零刚度自由度）
//    ② 竖向反力仍与荷载平衡（新增的板刚度必须自洽）
//    ③ 剖分使板变柔（单格单元对板是严重偏刚的），且加密后趋于收敛
static void testSlabSubdivision() {
  std::printf("\n== 3b. 楼板剖分（MITC4 板弯曲） ==\n");

  const double A = 6.0 * 6.0;
  const double q = 5.0;

  // 返回：是否可解、跨中挠度、柱脚反力合计、壳单元数
  struct R { bool ok{false}; double wc{0}, sumRz{0}; int nShell{0}; };
  auto run = [&](int div) -> R {
    R out;
    AxisGrid g;
    g.addAxisX(0.0); g.addAxisX(6.0);
    g.addAxisY(0.0); g.addAxisY(6.0);
    g.addStory(0.0); g.addStory(3.5);
    if (div > 1) { g.setDivX(0, div); g.setDivY(0, div); }

    Model m;
    GridMesh mesh(g, defaultSpec());
    const auto rr = mesh.generate(m);
    if (!rr.ok) return out;
    for (Id n : mesh.baseNodes()) m.fixAll(n);
    for (auto& e : m.elements())
      if (e->type() == ElementType::Shell4) {
        static_cast<ShellElement*>(e.get())->setTransversePressure(-q);
        ++out.nShell;
      }
    m.assignDofs();
    const auto diag = m.check();
    if (!diag.ok()) {
      for (const auto& e : diag.errors)
        if (e.find("uz 无刚度") != std::string::npos)
          std::printf("       [诊断] %s\n", e.c_str());
      return out;
    }
    StaticAnalysis sa(m);
    const auto res = sa.solve();
    if (!res.ok) return out;
    for (Id n : mesh.baseNodes()) out.sumRz += res.reaction[static_cast<size_t>(n) * 6 + 2];
    // 板格中心节点（细分网格的中间点）
    const Id c = mesh.nodeAt(div / 2, div / 2, 1);
    if (c < 0) return out;
    out.wc = res.u[static_cast<size_t>(c) * 6 + 2];
    out.ok = true;
    return out;
  };

  // ① 剖分后可解（MITC4 到位的直接证据）
  const R r1 = run(1);
  const R r3 = run(3);
  const R r5 = run(5);
  const R r7 = run(7);
  const R r9 = run(9);
  check(r1.ok, "单格楼板可解");
  check(r3.ok, "3×3 剖分楼板可解（内部节点 uz 有刚度 —— MITC4 已实现）");
  check(r5.ok, "5×5 剖分楼板可解");
  std::printf("       壳单元数：1 格 → %d，3×3 → %d，5×5 → %d\n",
              r1.nShell, r3.nShell, r5.nShell);
  if (!r1.ok || !r3.ok || !r5.ok || !r7.ok || !r9.ok) return;

  // ② 竖向平衡仍然成立
  for (const auto& kv : {std::make_pair(1, &r1), std::make_pair(3, &r3), std::make_pair(5, &r5)}) {
    const double want = q * A;
    std::printf("       %d×%d：ΣRz = %.6f（荷载 %.1f，相对差 %.2e）\n",
                kv.first, kv.first, kv.second->sumRz, want,
                std::abs(kv.second->sumRz - want) / want);
    checkNear(kv.second->sumRz, want, 1e-9, "剖分后竖向反力合计 = 板面总荷载");
  }

  // ③ 剖分使板变柔，且加密后收敛
  //
  //  【为什么这里不能用"两次加密相差 < x%"作判据】
  //  这是【四角点支承】的板：角点反力集中在极小的面积上，
  //  薄板理论在角点处是真的应力奇异，挠度收敛天然很慢（近似 O(h·log h)）。
  //  用"区间比值"作判据会误判成"没收敛"。
  //  正确的判据是【增量的衰减趋势】：逐次加密的增量必须单调减小。
  std::printf("       跨中挠度：1 格 %.6e，3×3 %.6e，5×5 %.6e，7×7 %.6e，9×9 %.6e\n",
              r1.wc, r3.wc, r5.wc, r7.wc, r9.wc);
  check(std::abs(r3.wc) > std::abs(r1.wc) * 1.2,
        "3×3 剖分比单格更柔（单格板单元对板面外刚度严重偏刚）");
  check(std::abs(r5.wc) > std::abs(r3.wc) && std::abs(r7.wc) > std::abs(r5.wc) &&
            std::abs(r9.wc) > std::abs(r7.wc),
        "挠度随网格加密单调增大（向真实解趋近）");
  {
    const double i35 = std::abs(r5.wc - r3.wc) / std::abs(r5.wc);
    const double i57 = std::abs(r7.wc - r5.wc) / std::abs(r7.wc);
    const double i79 = std::abs(r9.wc - r7.wc) / std::abs(r9.wc);
    std::printf("       逐次相对增量：3→5 %.1f%%，5→7 %.1f%%，7→9 %.1f%%\n",
                i35 * 100, i57 * 100, i79 * 100);
    check(i79 < i57 && i57 < i35, "逐次增量单调减小（收敛性成立）");
    check(i79 < 0.06, "最后一次加密的相对变化 < 6%（已进入收敛区）");
  }
  check(r5.wc < 0.0 && r9.wc < 0.0, "跨中挠度向下（荷载向下）");

  // ④ 挠度量级核查：与"四角支承方板"的粗略量级比对
  //    四角点支承板的挠度系数约为 (2~4)×10⁻³·q·a⁴/D 量级，
  //    这里只要求落在合理区间（1e-4 ~ 1e-2 · q·a⁴/D），
  //    用于拦住"量级差几百倍"这类错误。
  {
    const double D = 3.0e7 * std::pow(0.12, 3) / (12.0 * (1.0 - 0.2 * 0.2));
    const double ref = q * std::pow(6.0, 4) / D;
    const double ratio = std::abs(r9.wc) / ref;
    std::printf("       |w|/(q·a⁴/D) = %.5f（0.12 m 板，D=%.1f kN·m）\n", ratio, D);
    check(ratio > 1e-4 && ratio < 1e-2, "跨中挠度与 q·a⁴/D 的量级关系合理");
  }
}

// =============================================================================
// 4. 构件开关与参数校验
// =============================================================================
static void testSwitches() {
  std::printf("\n== 4. 构件开关与参数校验 ==\n");

  // 关闭一根梁、一块板、一根柱
  {
    AxisGrid g;
    for (double y : {0.0, 6.0, 12.0}) g.addAxisX(y);
    for (double x : {0.0, 6.0}) g.addAxisY(x);
    g.addStory(0.0);
    g.addStory(4.0);
    // 默认全开时的数量
    Model m0;
    GridMesh mesh0(g, defaultSpec());
    const auto r0 = mesh0.generate(m0);
    check(r0.ok, "默认（全开）生成成功");

    // 关掉 X 向轴线 0 上的跨 0 梁、板格 (0,0)、交点 (0,0) 的柱
    g.setBeamAlongX(0, 0, false);
    g.setSlab(0, 0, false);
    g.setColumn(0, 0, false);
    Model m1;
    GridMesh mesh1(g, defaultSpec());
    const auto r1 = mesh1.generate(m1);
    check(r1.ok, "关闭部分构件后生成成功");
    std::printf("       默认：柱%d 梁%d 板%d ／ 关闭后：柱%d 梁%d 板%d\n",
                r0.columns, r0.beams, r0.slabs, r1.columns, r1.beams, r1.slabs);
    checkInt(r0.columns - r1.columns, 1, "关掉 1 根柱");
    checkInt(r0.slabs - r1.slabs, 1, "关掉 1 块板（1 板格 × 1 剖分）");
    // 沿 X 的梁：跨 0 少 1 段 × 2 层 = 2
    checkInt(r0.beams - r1.beams, 2, "关掉 1 跨梁 × 2 层 = 2 根");
  }

  // 参数错误必须被拒绝
  {
    AxisGrid g;
    g.addAxisX(0.0);
    g.addAxisY(0.0);
    g.addStory(0.0);
    Model m;
    GridMesh mesh(g, defaultSpec());
    const auto r = mesh.generate(m);
    check(!r.ok, "轴线不足 2 条时被拒绝");
    std::printf("       %s\n", r.message().c_str());
  }
  {
    AxisGrid g;
    g.addAxisX(0.0);
    g.addAxisX(6.0);
    g.addAxisY(0.0);
    g.addAxisY(6.0);
    Model m;
    GridMesh mesh(g, defaultSpec());
    const auto r = mesh.generate(m);
    check(!r.ok, "没有楼层标高时被拒绝");
    std::printf("       %s\n", r.message().c_str());
  }
  {
    // 乱序输入
    AxisGrid g;
    g.addAxisX(12.0);
    g.addAxisX(0.0);
    g.addAxisY(0.0);
    g.addAxisY(6.0);
    g.addStory(0.0);
    g.addStory(4.0);
    Model m;
    GridMesh mesh(g, defaultSpec());
    const auto r = mesh.generate(m);
    check(!r.ok, "轴线乱序时被拒绝（要求先 sortAxes 再设开关）");
    std::printf("       %s\n", r.message().c_str());

    check(g.sortAxes(), "sortAxes 报告发生了重排");
    check(g.sorted(), "sortAxes 后已升序");
    Model m2;
    GridMesh mesh2(g, defaultSpec());
    const auto r2 = mesh2.generate(m2);
    check(r2.ok, "排序后生成成功");
  }

  // 剖分数非法
  {
    AxisGrid g;
    g.addAxisX(0.0); g.addAxisX(6.0);
    g.addAxisY(0.0); g.addAxisY(6.0);
    g.addStory(0.0); g.addStory(4.0);
    check(!g.setDivX(0, 0), "剖分数 0 被拒绝");
    check(!g.setDivX(5, 2), "越界的跨度索引被拒绝");
    check(g.setDivX(0, 3), "合法剖分数被接受");
  }
}

// =============================================================================
// 5. mergeCoincidentNodes：手工建模的"看起来连上了"修复
// =============================================================================
static void testMerge() {
  std::printf("\n== 5. mergeCoincidentNodes ==\n");

  const auto mat = Material::concreteC(30);
  const auto sec = section::rect(0.3, 0.6);
  Model m;
  // 三段梁，但中间那段的端点【重复建了两个节点】
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id n1 = m.addNode({2, 0, 0}, 1);
  const Id n1b = m.addNode({2, 0, 0}, 1);        // 与 n1 坐标重合
  const Id n2 = m.addNode({4, 0, 0}, 1);
  m.addBeam(n0, n1, sec, mat);
  m.addBeam(n1b, n2, sec, mat);

  m.assignDofs();
  const auto d0 = m.check();
  bool hasDup = false;
  for (const auto& e : d0.errors)
    if (e.find("坐标重合") != std::string::npos) hasDup = true;
  check(hasDup, "check() 检出坐标重合的节点");

  const auto rep = mergeCoincidentNodes(m, 1e-6);
  std::printf("       合并 %d 个节点，重定向 %d 个单元引用\n", rep.merged, rep.redirected);
  checkInt(rep.merged, 1, "合并掉 1 个重复节点");
  checkInt(rep.redirected, 1, "重定向 1 个单元引用");

  m.assignDofs();
  const auto d1 = m.check();
  bool stillDup = false;
  for (const auto& e : d1.errors)
    if (e.find("坐标重合") != std::string::npos) stillDup = true;
  check(!stillDup, "合并后不再报坐标重合");

  // 合并后结构仍然连通且能算：固定一端，另一端加力
  m.fixAll(n0);
  m.addNodeForce(n2, {0, 0, -10});
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto res = sa.solve();
  check(res.ok, "合并后的模型能求解");
  if (res.ok) {
    const double w = res.u[static_cast<size_t>(n2) * 6 + 2];
    //
    // 【用 Iz 而不是 Iy —— 这里错过一次】
    // 梁沿 X、up=(0,0,-1) ⇒ 局部 y' = −Z（竖向）、z' = +Y。
    // 竖向荷载引起的挠度沿 y'、转动绕 z' ⇒ 用【绕 z' 的惯性矩 Iz】（强轴）。
    // 用 Iy 会差 Iz/Iy = (h/b)² = 4 倍，而且符号方向看不出问题。
    const double exact = -10.0 * 64.0 / (3.0 * mat.E * sec.Iz);   // PL³/3EI
    std::printf("       端部挠度 %.6e，悬臂解析解 %.6e（Iz=%.6g, Iy=%.6g）\n",
                w, exact, sec.Iz, sec.Iy);
    // 两段梁已连通 ⇒ 等效于 4m 悬臂梁；Timoshenko 剪切变形使其略软（约 1.6%）
    checkNear(w / exact, 1.0, 0.03, "合并后 = 4m 悬臂梁（挠度比对解析解）");
  }
}

// =============================================================================
// 6. 多层框架：编号策略与求解规模
// =============================================================================
static void testMultiStory() {
  std::printf("\n== 6. 多层框架 ==\n");

  AxisGrid g;
  for (double y : {0.0, 6.0, 12.0, 18.0}) g.addAxisX(y);      // 4 条
  for (double x : {0.0, 6.0, 12.0, 18.0}) g.addAxisY(x);      // 4 条
  for (double z : {0.0, 3.6, 7.2, 10.8, 14.4}) g.addStory(z); // 5 层（4 层柱）

  Model m;
  GridMesh mesh(g, defaultSpec());
  const auto r = mesh.generate(m);
  check(r.ok, "4×4 轴网 5 层框架生成成功");
  if (!r.ok) { std::printf("       %s\n", r.error.c_str()); return; }
  std::printf("       %s\n", r.message().c_str());
  std::printf("       节点总数 %d\n", m.nodeCount());

  // 手算：柱 4×4×4 层 = 64；板 3×3 板格 × 4 层（底层无板）= 36
  checkInt(r.columns, 64, "柱数 = 16 交点 × 4 层 = 64");
  checkInt(r.slabs, 36, "板数 = 9 板格 × 4 层 = 36");

  for (Id n : mesh.baseNodes()) m.fixAll(n);
  // 屋面加竖向荷载
  for (auto& e : m.elements())
    if (e->type() == ElementType::Shell4)
      static_cast<ShellElement*>(e.get())->setTransversePressure(-3.0);

  m.assignDofs();
  const auto diag = m.check();
  check(diag.ok(), "模型校核通过");
  for (const auto& e : diag.errors) std::printf("       [错误] %s\n", e.c_str());

  StaticAnalysis sa(m);
  const auto res = sa.solve();
  check(res.ok, "静力求解成功");
  if (!res.ok) { std::printf("       %s\n", res.message.c_str()); return; }

  // 板在【除底层外的每个楼面】都有：4 个楼面 × 18×18 m² × 3 kPa
  //
  // 【这里错过一次：只算了"屋面"一层】实际板数 36 = 9 板格 × 4 层，
  // 所以总荷载是 4 倍 —— 求解器给的 3888 是对的，是测试期望写错了。
  const int nSlabLevels = 4;
  const double totalLoad = 3.0 * 18.0 * 18.0 * nSlabLevels;
  double sumRz = 0.0;
  for (Id n : mesh.baseNodes()) sumRz += res.reaction[static_cast<size_t>(n) * 6 + 2];
  std::printf("       %d 个楼面的板，总荷载 = %.2f kN，柱脚竖向反力合计 = %.4f kN\n",
              nSlabLevels, totalLoad, sumRz);
  checkNear(sumRz, totalLoad, 1e-8, "柱脚竖向反力合计 = 各层板面荷载之和");

  std::printf("       自由度 %d，nnz=%zu，nnzL=%zu，残差 %.2e，耗时 %.3f s\n",
              res.ndof, res.nnz, res.nnzL, res.residual, res.seconds);
  check(res.residual < 1e-9, "多层模型残差 < 1e-9");
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  轴网建模（GridMesh）验证\n");
  std::printf("======================================================\n");
  testTopology();
  testSlabBeamSharedNodes();
  testSlabSubdivision();
  testMechanics();
  testSwitches();
  testMerge();
  testMultiStory();
  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}
