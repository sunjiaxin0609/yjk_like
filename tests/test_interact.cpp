// =============================================================================
//  tests/test_interact.cpp  ——  P3 交互建模 T1：内核可变接口 + 序列化往返
//
//  本套件验证"脚本层与交互层共用的同一套内核接口"：
//    · Model::addColumn/addWall/addNodeLoad/removeElement/removeNode/clear
//    · Model::element(id) 按 id 索引查询
//
//  核心验收（round-trip 一致性）：
//    Model 建混合模型（节点 + 梁 + 柱 + 墙 + 壳 + 弹簧 + 约束 + 力/矩/重量）
//     → modelToYjk() 序列化为 .yjk 脚本
//     → ModelScript::parse + build 读回
//     → 节点 / 单元 / 荷载集合逐字段一致（retired 节点重映射后仍一致）。
//
//  这是 T4 的基础：GUI 里手点的构件，将导出成同一种脚本文本；
//  脚本能读回与源模型一致的模型，GUI 与脚本的操作才真正等价。
// =============================================================================
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "yjk/io/ModelScript.h"
#include "yjk/io/YjkWriter.h"
#include "yjk/model/Model.h"

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

// 向量/数组逐项近似比较
template <size_t N>
static void checkArrNear(const std::array<double, N>& got,
                         const std::array<double, N>& want, double tol,
                         const char* what) {
  bool ok = true;
  for (size_t i = 0; i < N; ++i) {
    const double err = std::abs(got[i] - want[i]);
    if (err > tol) { ok = false; break; }
  }
  if (ok) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got {", what);
    for (size_t i = 0; i < N; ++i) std::printf("%s%.10g", i ? ", " : "", got[i]);
    std::printf("}, want {");
    for (size_t i = 0; i < N; ++i) std::printf("%s%.10g", i ? ", " : "", want[i]);
    std::printf("})\n");
  }
}

static void checkVecNear(const Vec3& got, const Vec3& want, double tol, const char* what) {
  const bool ok = std::abs(got.x - want.x) <= tol &&
                  std::abs(got.y - want.y) <= tol &&
                  std::abs(got.z - want.z) <= tol;
  if (ok) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got (%.10g,%.10g,%.10g), want (%.10g,%.10g,%.10g))\n",
                what, got.x, got.y, got.z, want.x, want.y, want.z);
  }
}

// -----------------------------------------------------------------------------
//  1. 可变接口行为：addColumn / addWall / addNodeLoad / remove / clear / element(id)
// -----------------------------------------------------------------------------
static void testMutableApi() {
  std::printf("\n== 1. 内核可变接口行为 ==\n");

  Model m;
  checkInt(m.nodeCount(), 0, "空模型节点数 = 0");
  checkInt(m.elementCount(), 0, "空模型单元数 = 0");

  // 6 个节点：底部矩形 4 角 + 顶部两点
  const Id n0 = m.addNode(Vec3{0, 0, 0}, 0);
  const Id n1 = m.addNode(Vec3{6, 0, 0}, 0);
  const Id n2 = m.addNode(Vec3{6, 6, 0}, 0);
  const Id n3 = m.addNode(Vec3{0, 6, 0}, 0);
  const Id n4 = m.addNode(Vec3{0, 0, 3.6}, 1);
  const Id n5 = m.addNode(Vec3{6, 0, 3.6}, 1);
  checkInt(m.nodeCount(), 6, "添加 6 个节点");
  checkInt(n0, 0, "节点 id 按出现顺序分配（0）");
  checkInt(n5, 5, "节点 id 按出现顺序分配（5）");

  const Material steel = Material::steelQ(355);
  const Material conc = Material::concreteC(30);
  const SectionProperties colSec = section::rect(0.5, 0.5);
  const SectionProperties beamSec = section::rect(0.3, 0.6);

  // addColumn：柱 = 竖向梁，默认 up {0,-1,0}
  BeamElement* col = m.addColumn(n0, n4, colSec, conc);
  check(col != nullptr, "addColumn 成功");
  if (col) {
    check(col->type() == ElementType::Beam3D, "柱是梁单元（Beam3D）");
    checkVecNear(col->upHint(), Vec3{0, -1, 0}, 1e-15, "addColumn 默认 up = {0,-1,0}");
  }

  // addBeam：梁（水平，up 显式给竖直向下 {0,0,-1}）
  BeamElement* beam = m.addBeam(n1, n2, beamSec, steel, Vec3{0, 0, -1});
  check(beam != nullptr, "addBeam 成功");
  if (beam) checkVecNear(beam->upHint(), Vec3{0, 0, -1}, 1e-15, "addBeam up 保持传入值");

  // element(id) 按 id 索引
  check(m.element(0) == col, "element(0) 返回第 0 个单元（柱）");
  check(m.element(1) == beam, "element(1) 返回第 1 个单元（梁）");
  check(m.element(0)->type() == ElementType::Beam3D, "element(0) 类型 Beam3D");

  // addWall：墙 = 壳 + isWall
  ShellProperties wp;
  wp.thickness = 0.25;
  ShellElement* wall = m.addWall({n0, n1, n5, n4}, wp);
  check(wall != nullptr, "addWall 成功");
  if (wall) {
    check(wall->type() == ElementType::Shell4, "墙是壳单元（Shell4）");
    check(wall->isWall(), "addWall 置 isWall = true");
  }

  // addShell：壳（非墙）
  ShellProperties sp;
  sp.thickness = 0.12;
  ShellElement* shell = m.addShell({n1, n2, n3, n0}, sp);
  check(shell != nullptr, "addShell 成功");
  if (shell) check(!shell->isWall(), "addShell 默认 isWall = false");

  // addSpring
  std::array<double, 6> ks{{1e4, 1e4, 1e4, 0.0, 0.0, 1e3}};
  SpringElement* spr = m.addSpring(n3, n2, ks, Vec3{0, 0, 1});
  check(spr != nullptr, "addSpring 成功");
  checkInt(m.elementCount(), 5, "5 个单元（柱/梁/墙/壳/弹簧）");

  // addNodeLoad：力 + 力矩一次施加
  m.addNodeLoad(n4, Vec3{10, 0, -20}, Vec3{0, 5, 0});
  checkVecNear(m.node(n4).force, Vec3{10, 0, -20}, 1e-12, "addNodeLoad 施加力");
  checkVecNear(m.node(n4).moment, Vec3{0, 5, 0}, 1e-12, "addNodeLoad 施加力矩");
  m.addNodeWeight(n4, 30.0);
  checkNear(m.node(n4).weight, 30.0, 1e-12, "addNodeWeight 施加重量");

  // removeElement
  check(m.removeElement(1), "removeElement(1) 成功");
  checkInt(m.elementCount(), 4, "删除后 4 个单元");
  check(m.element(1) == wall, "删除后 element(1) 前移为新墙单元");
  check(m.removeElement(99) == false, "removeElement(99) 越界返回 false");

  // removeNode：被单元引用的节点不能删
  check(m.removeNode(n0) == false, "removeNode(被引用节点) 拒绝");
  // 未被引用的孤立节点可以删（n5 被 wall 引用 → 用一个真正孤立的节点）
  const Id nIso = m.addNode(Vec3{99, 99, 99}, 0);
  check(m.removeNode(nIso), "removeNode(孤立节点) 成功");
  check(m.node(nIso).retired, "被删节点置 retired 标记");
  check(m.removeNode(nIso) == false, "重复删除 retired 节点被拒绝");

  // clear
  m.clear();
  checkInt(m.nodeCount(), 0, "clear 后节点 = 0");
  checkInt(m.elementCount(), 0, "clear 后单元 = 0");
}

// -----------------------------------------------------------------------------
//  2. 序列化往返：混合模型 → .yjk 文本 → 读回，节点/单元/荷载逐项一致
// -----------------------------------------------------------------------------
static void testRoundTrip() {
  std::printf("\n== 2. Model → modelToYjk → ModelScript 往返一致 ==\n");

  // ---- 源模型：手法与 GUI 交互层完全相同（走同一套内核接口）----
  Model src;
  const Id a = src.addNode(Vec3{0, 0, 0}, 0);
  const Id b = src.addNode(Vec3{6, 0, 0}, 0);
  const Id c = src.addNode(Vec3{6, 6, 0}, 0);
  const Id d = src.addNode(Vec3{0, 6, 0}, 0);
  const Id e = src.addNode(Vec3{0, 0, 3.6}, 1);
  const Id f = src.addNode(Vec3{6, 0, 3.6}, 1);
  // 孤立节点 g：将被 removeNode 删除，验证序列化的 retired 重映射
  const Id g = src.addNode(Vec3{50, 50, 50}, 2);

  const Material steel = Material::steelQ(355);
  const Material conc = Material::concreteC(30);

  src.addColumn(a, e, section::rect(0.5, 0.5), conc);                     // 柱
  src.addBeam(b, c, section::rect(0.3, 0.6), steel, Vec3{0, 0, -1});      // 梁
  ShellProperties wp;
  wp.thickness = 0.25;
  wp.E = 3.15e7;
  src.addWall({a, b, f, e}, wp);                                          // 墙
  ShellProperties sp;
  sp.thickness = 0.12;
  src.addShell({b, c, d, a}, sp);                                         // 壳
  src.addSpring(d, c, std::array<double, 6>{{2e4, 2e4, 2e4, 0, 0, 1e3}},
                Vec3{0, 0, 1});                                           // 弹簧

  // 约束：柱脚 a 固定，节点 d 只约束竖向
  src.fixNode(a, true, true, true, true, true, true);
  src.fixNode(d, false, false, true, false, false, false);

  // 荷载：力、力矩、重量
  src.addNodeLoad(e, Vec3{12, -3, -25}, Vec3{0, 4, 1});
  src.addNodeWeight(e, 40.0);
  src.addNodeForce(f, Vec3{0, 0, -10});

  // 删除孤立节点 g：序列化必须跳过它并把单元引用的节点重映射
  check(src.removeNode(g), "源模型删除孤立节点 g");

  // ---- 序列化 ----
  const std::string text = io::modelToYjk(src);
  std::printf("  %zu bytes: %s\n", text.size(),
              text.empty() ? "(空)" : "已写出");
  if (text.empty()) { ++g_fail; return; }

  // ---- 读回 ----
  io::ModelScript script;
  if (!script.parse(text)) {
    ++g_fail;
    std::printf("  [FAIL] 脚本解析失败：\n");
    for (const auto& err : script.errors()) std::printf("       %s\n", err.c_str());
    return;
  }
  Model dst;
  std::string info;
  if (!script.build(dst, &info)) {
    ++g_fail;
    std::printf("  [FAIL] 脚本建模失败：\n");
    for (const auto& err : script.errors()) std::printf("       %s\n", err.c_str());
    return;
  }

  // ---- 节点重映射表：源活跃节点 → 序列化编号 ----
  std::vector<Id> remap(static_cast<size_t>(src.nodeCount()), -1);
  Id next = 0;
  for (Id n = 0; n < src.nodeCount(); ++n)
    if (!src.node(n).retired) remap[static_cast<size_t>(n)] = next++;
  const Id alive = next;                          // 活跃节点数

  checkInt(static_cast<long long>(alive), 6, "源模型活跃节点 = 6");
  checkInt(dst.nodeCount(), 6, "读回节点数 = 6（retired 已跳过）");

  // ---- 逐项比对 ----
  // 触点：重映射编号 + 坐标 + 层号
  for (Id n = 0; n < src.nodeCount(); ++n) {
    const Node& sn = src.node(n);
    if (sn.retired) continue;
    const Id dn = remap[static_cast<size_t>(n)];
    const Node& dnode = dst.node(dn);
    char buf[128];
    std::snprintf(buf, sizeof buf, "节点 %lld 坐标", static_cast<long long>(n));
    checkVecNear(dnode.r, sn.r, 1e-12, buf);
    std::snprintf(buf, sizeof buf, "节点 %lld 层号", static_cast<long long>(n));
    checkInt(dnode.story, sn.story, buf);
    // 约束
    for (int k = 0; k < 6; ++k) {
      std::snprintf(buf, sizeof buf, "节点 %lld 约束分量 %d", static_cast<long long>(n), k);
      check(dnode.fixed[static_cast<size_t>(k)] == sn.fixed[static_cast<size_t>(k)], buf);
    }
    // 荷载
    checkVecNear(dnode.force, sn.force, 1e-12, "节点力一致");
    checkVecNear(dnode.moment, sn.moment, 1e-12, "节点力矩一致");
    checkNear(dnode.weight, sn.weight, 1e-12, "节点重量一致");
  }

  // 单元：类型、节点引用（重映射后）、几何/截面/材料/刚度
  checkInt(dst.elementCount(), src.elementCount(), "单元数量一致");
  for (Id k = 0; k < src.elementCount(); ++k) {
    const Element* se = src.element(k);
    const Element* de = dst.element(k);
    char buf[128];
    std::snprintf(buf, sizeof buf, "单元 %lld 类型", static_cast<long long>(k));
    check(de->type() == se->type(), buf);
    if (de->type() != se->type()) continue;
    std::snprintf(buf, sizeof buf, "单元 %lld 节点引用（重映射）", static_cast<long long>(k));
    bool nsOk = de->nodes().size() == se->nodes().size();
    if (nsOk)
      for (size_t i = 0; i < se->nodes().size(); ++i)
        if (de->nodes()[i] != remap[static_cast<size_t>(se->nodes()[i])]) { nsOk = false; break; }
    check(nsOk, buf);

    switch (se->type()) {
      case ElementType::Beam3D: {
        const auto* sb = static_cast<const BeamElement*>(se);
        const auto* db = static_cast<const BeamElement*>(de);
        checkVecNear(db->upHint(), sb->upHint(), 1e-12, "梁 up 一致");
        const SectionProperties& ss = sb->section();
        const SectionProperties& ds = db->section();
        checkNear(ds.A, ss.A, 1e-12, "梁截面 A 一致");
        checkNear(ds.Iy, ss.Iy, 1e-12, "梁截面 Iy 一致");
        checkNear(ds.Iz, ss.Iz, 1e-12, "梁截面 Iz 一致");
        checkNear(ds.J, ss.J, 1e-12, "梁截面 J 一致");
        checkNear(ds.Asy, ss.Asy, 1e-12, "梁截面 Asy 一致");
        checkNear(ds.Asz, ss.Asz, 1e-12, "梁截面 Asz 一致");
        checkNear(ds.hy, ss.hy, 1e-12, "梁截面 hy 一致");
        checkNear(ds.hz, ss.hz, 1e-12, "梁截面 hz 一致");
        const Material& sm = sb->material();
        const Material& dm = db->material();
        checkNear(dm.E, sm.E, 1e-12, "梁材料 E 一致");
        checkNear(dm.nu, sm.nu, 1e-12, "梁材料 nu 一致");
        checkNear(dm.gamma, sm.gamma, 1e-12, "梁材料 gamma 一致");
        checkNear(dm.fy, sm.fy, 1e-12, "梁材料 fy 一致");
        checkNear(dm.G, sm.G, 1e-12, "梁材料 G 一致");
        break;
      }
      case ElementType::Shell4: {
        const auto* ssh = static_cast<const ShellElement*>(se);
        const auto* dsh = static_cast<const ShellElement*>(de);
        check(dsh->isWall() == ssh->isWall(), "墙标志一致");
        const ShellProperties& ss = ssh->properties();
        const ShellProperties& ds = dsh->properties();
        checkNear(ds.thickness, ss.thickness, 1e-12, "壳厚度一致");
        checkNear(ds.E, ss.E, 1e-12, "壳 E 一致");
        checkNear(ds.nu, ss.nu, 1e-12, "壳 nu 一致");
        checkNear(ds.density, ss.density, 1e-12, "壳 density 一致");
        break;
      }
      case ElementType::Spring3D: {
        const auto* spp = static_cast<const SpringElement*>(se);
        const auto* dpp = static_cast<const SpringElement*>(de);
        checkVecNear(dpp->upHint(), spp->upHint(), 1e-12, "弹簧 up 一致");
        checkArrNear<6>(dpp->stiffness(), spp->stiffness(), 1e-9, "弹簧刚度一致");
        break;
      }
    }
  }

  // 荷载集合：活跃节点每个的力/矩/重全部比对过（上面循环）。
  // 额外确认：序列化文本里出现了三类荷载命令。
  check(text.find("nodeload node") != std::string::npos, "序列化含 nodeload");
  check(text.find("nodemoment node") != std::string::npos, "序列化含 nodemoment");
  check(text.find("nodeweight node") != std::string::npos, "序列化含 nodeweight");
  check(text.find("material ") != std::string::npos, "序列化含 material raw");
  check(text.find("section ") != std::string::npos, "序列化含 section raw");
  check(text.find("wall ") != std::string::npos, "序列化含 wall（墙标志）");
  check(text.find("spring ") != std::string::npos, "序列化含 spring");
}

// -----------------------------------------------------------------------------
//  3. 脚本命令与内核接口等价：GUI 导出的脚本能被脚本引擎原样读回
// -----------------------------------------------------------------------------
static void testScriptEquivalence() {
  std::printf("\n== 3. 脚本命令与内核接口等价 ==\n");

  // 直接写一段"原始构件"脚本（YjkWriter 输出的形态），验证它能被
  // ModelScript 完整读回 —— GUI 另存为 .yjk 产物即这种文本。
  const std::string scriptText =
      "material M0 raw 2.06e8 0.3 78.5 355 490 0 0\n"
      "material M1 raw 3.0e7 0.2 25 0 0 14.3 0\n"
      "section S0 raw 0.25 0.005208333333333333 0.005208333333333333 "
      "0.004069010416666667 0.20833333333333334 0.20833333333333334 "
      "0.0026041666666666665 0.0026041666666666665 0.25 0.25\n"
      "section S1 raw 0.18 0.0054 0.00135 0.0004674063 0.15 0.15 "
      "0.0027 0.000675 0.3 0.15\n"
      "node 0 0 0 0\n"
      "node 6 0 0 0\n"
      "node 6 6 0 0\n"
      "node 0 0 3.6 1\n"
      "column 0 3 S0 M1\n"
      "beam 1 2 S1 M0\n"
      "shell 0 1 3 4 0\n"            // 先用 shell，4 个节点 0 1 3 4（越界→被拦截）
      "fix node 0 1 1 1 1 1 1\n"
      "nodeload node 0 10 0 -20\n"
      "nodemoment node 0 0 5 0\n";

  // 先验证：越界节点号必须被脚本引擎拦截（不能静默建出悬空单元）
  // —— 脚本 1 里 shell 引用节点 4，但只声明了 0..3。
  {
    io::ModelScript bad;
    Model mb;
    const bool p = bad.parse(scriptText);
    const bool b = p && bad.build(mb);
    check(p, "越界脚本 parse 成功（错误在 build 阶段收集）");
    check(!b, "越界节点号被 build 拦截（shell 引用不存在的节点 4）");
    if (b) for (const auto& e : bad.errors()) std::printf("       %s\n", e.c_str());
  }

  // 无论上述 shell 越界是否被拦，下面用合法的墙节点重新验证脚本 → 内核等价
  const std::string script2 =
      "material M0 raw 2.06e8 0.3 78.5 355 490 0 0\n"
      "section S0 raw 0.25 0.005208333333333333 0.005208333333333333 "
      "0.004069010416666667 0.20833333333333334 0.20833333333333334 "
      "0.0026041666666666665 0.0026041666666666665 0.25 0.25\n"
      "node 0 0 0 0\n"
      "node 6 0 0 0\n"
      "node 6 6 0 0\n"
      "node 0 6 0 0\n"
      "node 0 0 3.6 1\n"
      "column 0 4 S0 M0\n"
      "wall 0 1 4 3 0.25\n"           // 墙 0 1 4 3（厚度 0.25）
      "fix node 0 1 1 1 1 1 1\n"
      "nodeload node 0 10 0 -20\n"
      "nodemoment node 0 0 5 0\n"
      "nodeweight node 0 30\n";

  io::ModelScript s2;
  Model m2;
  const bool p2 = s2.parse(script2);
  const bool b2 = p2 && s2.build(m2);
  check(p2, "合法脚本 2 解析成功");
  check(b2, "合法脚本 2 建模成功");
  if (!b2) {
    for (const auto& e : s2.errors()) std::printf("       %s\n", e.c_str());
    ++g_fail;
    return;
  }
  checkInt(m2.nodeCount(), 5, "脚本 2 建出 5 节点");
  checkInt(m2.elementCount(), 2, "脚本 2 建出 2 单元（柱 + 墙）");
  if (m2.elementCount() == 2) {
    check(m2.element(0)->type() == ElementType::Beam3D, "脚本 2 第 0 单元为柱（Beam3D）");
    check(m2.element(1)->type() == ElementType::Shell4, "脚本 2 第 1 单元为墙（Shell4）");
    check(static_cast<ShellElement*>(m2.element(1))->isWall(), "脚本 2 墙标志正确");
    checkVecNear(static_cast<BeamElement*>(m2.element(0))->upHint(),
                 Vec3{0, -1, 0}, 1e-15, "脚本 column 默认 up = {0,-1,0}");
  }
  checkVecNear(m2.node(0).force, Vec3{10, 0, -20}, 1e-12, "脚本 nodeload 生效");
  checkVecNear(m2.node(0).moment, Vec3{0, 5, 0}, 1e-12, "脚本 nodemoment 生效");
  checkNear(m2.node(0).weight, 30.0, 1e-12, "脚本 nodeweight 生效");
  const Node& n0 = m2.node(0);
  check(n0.fixed[0] && n0.fixed[1] && n0.fixed[2] && n0.fixed[3] &&
            n0.fixed[4] && n0.fixed[5],
        "脚本 fix node 全约束生效");
}

// -----------------------------------------------------------------------------
//  4. 清空与重建：交互层"新建工程"语义
// -----------------------------------------------------------------------------
static void testClearRebuild() {
  std::printf("\n== 4. clear 后重建（新建工程语义）==\n");
  Model m;
  m.addNode(Vec3{0, 0, 0});
  m.addNode(Vec3{0, 0, 3.6});
  m.addColumn(0, 1, section::rect(0.5, 0.5), Material::concreteC(30));
  checkInt(m.nodeCount(), 2, "重建前 2 节点");
  checkInt(m.elementCount(), 1, "重建前 1 单元");
  m.clear();
  m.addNode(Vec3{3, 3, 3});
  checkInt(m.nodeCount(), 1, "clear 后可立即重建（旧编号不残留）");
  checkInt(m.elementCount(), 0, "clear 后单元归零");
}

// -----------------------------------------------------------------------------
//  5. T3 属性编辑：setNodeCoord / setBeamProps / setShellProps /
//     setNodeForce/setNodeMoment / Model::clone
//
//  验收对应 GUI 验收"改参后渲染与重算随之更新"：
//    · 改坐标 → 梁长度/质量立即变化（几何重建成功）
//    · 改到零长度 → 被拒绝且坐标回滚（几何退化保护）
//    · 加大截面 → 局部刚度阵变大（能算得更刚）
//    · 改壳厚度 → 厚度读回一致（写回成功）
//    · 节点荷载覆盖式写回
//    · clone 是独立深拷贝（GUI 重算快照不污染主模型）
// -----------------------------------------------------------------------------
static void testPropertyEdit() {
  std::printf("\n== 5. T3 属性编辑 ==");

  Model m;
  const Id n0 = m.addNode(Vec3{0, 0, 0});
  const Id n1 = m.addNode(Vec3{6, 0, 0});
  const Material conc = Material::concreteC(30);
  BeamElement* beam = m.addBeam(n0, n1, section::rect(0.3, 0.6), conc);
  check(beam != nullptr, "addBeam 成功");
  if (!beam) return;
  const double l0 = m.element(0)->length();
  checkNear(l0, 6.0, 1e-12, "梁初始长度 6 m");

  // -- 5a. setNodeCoord：改坐标 → 长度/质量联动 --
  check(m.setNodeCoord(n1, Vec3{3, 0, 0}), "setNodeCoord 移动到 (3,0,0) 成功");
  checkNear(m.element(0)->length(), 3.0, 1e-12, "梁长度随坐标变为 3 m");
  checkNear(m.element(0)->mass(), conc.gamma * 0.3 * 0.6 * 3.0, 1e-9,
            "梁质量随长度更新（=γ·A·L）");
  const double lShort = m.element(0)->length();

  // -- 5b. setNodeCoord：退化回滚 --
  check(!m.setNodeCoord(n1, Vec3{0, 0, 0}), "移到位移处使单元零长度 → 拒绝");
  checkNear(m.element(0)->length(), lShort, 1e-15, "拒绝后长度保持原值（回滚）");
  checkVecNear(m.node(n1).r, Vec3{3, 0, 0}, 1e-15, "拒绝后节点坐标回滚到 (3,0,0)");

  // -- 5c. setBeamProps：加大截面 → 轴向刚度 EA/L 变大 --
  const double kAxial0 = beam->section().A * beam->material().E / m.element(0)->length();
  check(m.setBeamProps(0, section::rect(0.6, 1.0), conc), "setBeamProps 改大截面成功");
  const double kAxial1 = beam->section().A * beam->material().E / m.element(0)->length();
  checkNear(beam->section().A, 0.6, 1e-12, "截面 A 写回 = 0.6 m²");
  check(kAxial1 > kAxial0 * 1.5, "轴向刚度 EA/L 显著变大");

  // -- 5d. setBeamProps：非法材料被拒（E ≤ 0）--
  Material bad = conc;
  bad.E = 0.0;
  check(!m.setBeamProps(0, section::rect(0.3, 0.6), bad),
        "E=0 的材料写回被拒绝");
  checkNear(beam->material().E, conc.E, 1e-12, "拒绝后材料 E 回滚为原值");

  // -- 5e. setShellProps：墙厚翻倍写回 --
  const std::vector<Id> sh = {n0, n1, m.addNode(Vec3{3, 3, 0}), m.addNode(Vec3{0, 3, 0})};
  ShellElement* wall = m.addWall(sh, ShellProperties{});
  check(wall != nullptr, "addWall 成功");
  if (wall) {
    ShellProperties p = wall->properties();
    p.thickness = 0.4;
    check(m.setShellProps(static_cast<Id>(m.elementCount() - 1), p), "setShellProps 厚度写回成功");
    checkNear(wall->properties().thickness, 0.4, 1e-15, "壳厚度读回 = 0.4 m");
    checkNear(wall->mass(), wall->area() * 0.4 * 2500.0 * 9.81 / 1000.0, 1e-9,
              "壳自重随厚度更新（=γ·t·A）");
  }

  // -- 5f. 节点荷载：覆盖式写回 --
  m.setNodeForce(n0, Vec3{10, 0, -20});
  m.setNodeForce(n0, Vec3{5, 0, -5});
  checkVecNear(m.node(n0).force, Vec3{5, 0, -5}, 1e-15, "setNodeForce 覆盖写回（不累加）");
  m.setNodeMoment(n0, Vec3{0, 3, 0});
  m.setNodeMoment(n0, Vec3{0, 0, 7});
  checkVecNear(m.node(n0).moment, Vec3{0, 0, 7}, 1e-15, "setNodeMoment 覆盖写回（不累加）");

  // -- 5g. clone：深拷贝独立，GUI 后台重算快照用 --
  const double lenBefore = m.element(0)->length();
  const Model cp = m.clone();
  checkInt(cp.nodeCount(), m.nodeCount(), "clone 节点数一致");
  checkInt(cp.elementCount(), m.elementCount(), "clone 单元数一致");
  checkVecNear(cp.node(n1).r, m.node(n1).r, 1e-15, "clone 坐标一致");
  checkNear(cp.element(0)->length(), lenBefore, 1e-15, "clone 梁长度一致");

  // 改副本不污染原模型（GUI 线程克隆 → 后台快照求解 → 丢弃）
  check(const_cast<Model&>(cp).setNodeCoord(n1, Vec3{9, 0, 0}),
        "clone 副本可编辑");
  checkNear(cp.element(0)->length(), 9.0, 1e-12, "副本长度变为 9 m");
  checkNear(m.element(0)->length(), lenBefore, 1e-15, "原模型长度不受副本编辑影响");
  checkInt(m.elementCount(), cp.elementCount(), "clone 后原模型完好可继续使用");
}

int main() {
  testMutableApi();
  testRoundTrip();
  testScriptEquivalence();
  testClearRebuild();
  testPropertyEdit();

  std::printf("\n====  test_interact：%d 通过，%d 失败 ====\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}