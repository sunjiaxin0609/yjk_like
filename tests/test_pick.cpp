// =============================================================================
//  tests/test_pick.cpp  ——  射线拾取核心（yjk::interact）单元测试
//
//  覆盖：节点球命中（正对/边缘/容差外）、梁圆柱命中（中点/端点/超容差）、
//  壳四边形命中（两个三角形分块/外部）、多对象遮挡取最近、tMax 截断、退化射线。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "yjk/interact/PickCore.h"

using yjk::Id;
using yjk::Vec3;
using yjk::interact::PickBeam;
using yjk::interact::PickHit;
using yjk::interact::PickKind;
using yjk::interact::PickNode;
using yjk::interact::PickShell;
using yjk::interact::Ray3;
using yjk::interact::pickNearest;
using yjk::interact::pointSegDist;
using yjk::interact::raySegDist;
using yjk::interact::rayTriHit;

static int g_pass = 0;
static int g_fail = 0;

static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}

static void checkNear(double got, double want, double tol, const char* what) {
  if (std::abs(got - want) <= tol) {
    ++g_pass;
    std::printf("  [ OK ] %s  (got %.10g, want %.10g)\n", what, got, want);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %.10g, want %.10g)\n", what, got, want);
  }
}

// 命中类型与 id 的联合断言
static void checkHit(const PickHit& h, PickKind k, Id id, const char* what) {
  if (h.kind == k && h.id == id) {
    ++g_pass;
    std::printf("  [ OK ] %s\n", what);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got kind=%d id=%lld, want kind=%d id=%lld)\n",
                what, static_cast<int>(h.kind), static_cast<long long>(h.id),
                static_cast<int>(k), static_cast<long long>(id));
  }
}

static void testNodePick() {
  std::printf("\n== 节点球命中 ==\n");
  const std::vector<PickNode> nodes{{0, {0, 0, 0}}, {1, {6, 0, 0}}};
  const std::vector<PickBeam> beams;
  const std::vector<PickShell> shells;

  // 正对节点 0：射线 (0,0,10) → 竖直向下穿过节点
  PickHit h = pickNearest(Ray3{Vec3{0, 0, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                          0.8, 0.5);
  checkHit(h, PickKind::Node, 0, "正对节点 0 命中");
  checkNear(h.t, 10.0, 1e-12, "命中距离 = 10");

  // 边缘：水平偏移 0.7 < 0.8 仍命中
  h = pickNearest(Ray3{Vec3{0.7, 0, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.8, 0.5);
  checkHit(h, PickKind::Node, 0, "偏移 0.7（容差 0.8 内）命中");

  // 超出容差：偏移 0.9 > 0.8
  h = pickNearest(Ray3{Vec3{0.9, 0, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.8, 0.5);
  check(h.kind == PickKind::None, "偏移 0.9 超出容差不命中");

  // 射线背向节点（从 z=10 继续往 +z 走）不命中
  h = pickNearest(Ray3{Vec3{0, 0, 10}, Vec3{0, 0, 1}}, nodes, beams, shells,
                  0.8, 0.5);
  check(h.kind == PickKind::None, "背向射线不命中");
}

static void testBeamPick() {
  std::printf("\n== 梁圆柱命中 ==\n");
  // 梁 0：(0,0,0)-(6,0,0)；梁 1：(0,6,0)-(6,6,0)
  const std::vector<PickNode> nodes;
  const std::vector<PickBeam> beams{{0, {0, 0, 0}, {6, 0, 0}}, {1, {0, 6, 0}, {6, 6, 0}}};
  const std::vector<PickShell> shells;

  // 中点正上方
  PickHit h = pickNearest(Ray3{Vec3{3, 0, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                          0.0, 0.5);
  checkHit(h, PickKind::Beam, 0, "梁 0 中点命中");
  checkNear(h.t, 10.0, 1e-12, "命中距离 = 10");

  // 端点上方
  h = pickNearest(Ray3{Vec3{0, 0.3, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.0, 0.5);
  checkHit(h, PickKind::Beam, 0, "梁 0 端点附近命中（s=0 钳制）");

  // 超出容差：横向偏移 0.8 > 0.5
  h = pickNearest(Ray3{Vec3{3, 0.8, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.0, 0.5);
  check(h.kind == PickKind::None, "横向偏移 0.8 超出梁容差不命中");

  // 探测第二根梁（y=6）
  h = pickNearest(Ray3{Vec3{3, 6, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.0, 0.5);
  checkHit(h, PickKind::Beam, 1, "梁 1 中点命中");

  // 斜射：从 (3,-3,5) 指向 +y 方向，应命中线段 [0,6] 上 s=0.5 处
  h = pickNearest(Ray3{Vec3{3, -3, 0}, Vec3{0, 1, 0}}, nodes, beams, shells,
                  0.0, 0.4);
  checkHit(h, PickKind::Beam, 0, "斜向射线命中线段内部");
  checkNear(h.t, 3.0, 1e-12, "斜向命中距离 = 3");
}

static void testShellPick() {
  std::printf("\n== 壳四边形命中 ==\n");
  // 壳 0：z=0 平面上的 6×6 正方形 (0,0,0)(6,0,0)(6,6,0)(0,6,0)
  // 壳 1：z=4 平面上的小四边形
  const std::vector<PickNode> nodes;
  const std::vector<PickBeam> beams;
  PickShell sh0;
  sh0.id = 0;
  sh0.p[0] = {0, 0, 0}; sh0.p[1] = {6, 0, 0};
  sh0.p[2] = {6, 6, 0}; sh0.p[3] = {0, 6, 0};
  PickShell sh1;
  sh1.id = 1;
  sh1.p[0] = {1, 1, 4}; sh1.p[1] = {3, 1, 4};
  sh1.p[2] = {3, 3, 4}; sh1.p[3] = {1, 3, 4};
  const std::vector<PickShell> shells{sh0, sh1};

  // 中心命中 z=0 壳。注意避开壳 1（(1,1)~(3,3) 区域），用左上角区域
  PickHit h = pickNearest(Ray3{Vec3{0.5, 5.5, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                          0.0, 0.0);
  checkHit(h, PickKind::Shell, 0, "壳 0 命中");
  checkNear(h.t, 10.0, 1e-12, "壳命中距离 = 10");

  // 对角线（三角形分块 2 的边）附近：x=4.5,y=4.5 在三角形 (0,0)(6,6)(0,6) 内
  h = pickNearest(Ray3{Vec3{4.5, 4.5, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.0, 0.0);
  checkHit(h, PickKind::Shell, 0, "壳 0 对角分块命中");

  // 四边形外部
  h = pickNearest(Ray3{Vec3{7, 3, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.0, 0.0);
  check(h.kind == PickKind::None, "壳外不命中");

  // 斜射：射线穿过 z=4 壳中心（三角形判交，法向任意）
  h = pickNearest(Ray3{Vec3{2, 2, 0}, Vec3{0, 0, 1}}, nodes, beams, shells,
                  0.0, 0.0);
  checkHit(h, PickKind::Shell, 1, "壳 1 从下方命中");
  checkNear(h.t, 4.0, 1e-12, "壳 1 命中距离 = 4");
}

static void testOcclusion() {
  std::printf("\n== 遮挡与优先级 ==\n");
  // 布局：
  //   壳 3：z=5 平面 6×6 四边形
  //   梁 2：(0,3,3)-(6,3,3)
  //   节点 0：(3,3,4)  —— 在梁正上方的射线路径上
  const std::vector<PickNode> nodes{{0, {3, 3, 4}}};
  const std::vector<PickBeam> beams{{2, {0, 3, 3}, {6, 3, 3}}};
  PickShell sh;
  sh.id = 3;
  sh.p[0] = {0, 0, 5}; sh.p[1] = {6, 0, 5};
  sh.p[2] = {6, 6, 5}; sh.p[3] = {0, 6, 5};
  const std::vector<PickShell> shells{sh};

  // 射线既穿壳又穿梁：壳在 z=5（t=5）先被碰到，梁 z=3（t=7）在后
  PickHit h = pickNearest(Ray3{Vec3{1, 3, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                          0.5, 0.5);
  checkHit(h, PickKind::Shell, 3, "壳挡梁：取最近（壳）");
  checkNear(h.t, 5.0, 1e-12, "壳命中距离 = 5");

  // 去掉壳后，同一射线命中梁（t=7）
  const std::vector<PickShell> noShells;
  h = pickNearest(Ray3{Vec3{1, 3, 10}, Vec3{0, 0, -1}}, nodes, beams, noShells,
                  0.5, 0.5);
  checkHit(h, PickKind::Beam, 2, "无壳时命中梁");
  checkNear(h.t, 7.0, 1e-12, "梁命中距离 = 7");

  // 节点 (3,3,4) 在梁正上方：射线 (3,3,10) 先碰节点（t=6）再碰梁（t=7）
  h = pickNearest(Ray3{Vec3{3, 3, 10}, Vec3{0, 0, -1}}, nodes, beams, noShells,
                  0.5, 0.5);
  checkHit(h, PickKind::Node, 0, "节点挡梁：取最近（节点）");
  checkNear(h.t, 6.0, 1e-12, "节点命中距离 = 6");

  // 同一射线但无节点与梁，只有壳
  h = pickNearest(Ray3{Vec3{2, 2, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.0, 0.0);
  checkHit(h, PickKind::Shell, 3, "仅壳命中");
  checkNear(h.t, 5.0, 1e-12, "壳命中距离 = 5");
}

static void testRanges() {
  std::printf("\n== tMax 截断与退化射线 ==\n");
  const std::vector<PickNode> nodes{{0, {0, 0, 0}}};
  const std::vector<PickBeam> beams;
  const std::vector<PickShell> shells;

  // tMax 截断：实际命中在 t=10，tMax=5 时不命中
  PickHit h = pickNearest(Ray3{Vec3{0, 0, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                          0.8, 0.5, 5.0);
  check(h.kind == PickKind::None, "tMax 截断（5 < 10）不命中");

  h = pickNearest(Ray3{Vec3{0, 0, 10}, Vec3{0, 0, -1}}, nodes, beams, shells,
                  0.8, 0.5, 20.0);
  checkHit(h, PickKind::Node, 0, "tMax 足够时命中");

  // 退化射线（零方向）
  h = pickNearest(Ray3{Vec3{0, 0, 10}, Vec3{0, 0, 0}}, nodes, beams, shells,
                  0.8, 0.5);
  check(h.kind == PickKind::None, "零方向射线不命中");
}

static void testGeomAtoms() {
  std::printf("\n== 几何原子 ==\n");
  // 点到线段距离
  checkNear(pointSegDist({3, 4, 0}, {0, 0, 0}, {6, 0, 0}), 4.0, 1e-12, "点到梁中点距离");
  checkNear(pointSegDist({9, 0, 0}, {0, 0, 0}, {6, 0, 0}), 3.0, 1e-12, "点到端点外距离");
  checkNear(pointSegDist({3, 0, 0}, {0, 0, 0}, {6, 0, 0}), 0.0, 1e-12, "点在线上距离 0");

  // 射线到线段
  double s = -1.0, t = -1.0;
  const bool ok = raySegDist(Ray3{Vec3{3, -2, 0}, Vec3{0, 1, 0}},
                             Vec3{0, 0, 0}, Vec3{6, 0, 0}, &s, &t);
  check(ok, "raySegDist 正常计算");
  checkNear(s, 0.5, 1e-12, "线段最近参数 s = 0.5");
  checkNear(t, 2.0, 1e-12, "射线参数 t = 2");

  // 射线到三角形（z=0 平面）
  double tHit = -1.0;
  const bool hit = rayTriHit(Ray3{Vec3{1, 1, 5}, Vec3{0, 0, -1}},
                             {0, 0, 0}, {3, 0, 0}, {0, 3, 0}, &tHit);
  check(hit, "射线与三角形相交");
  checkNear(tHit, 5.0, 1e-12, "三角形命中 t = 5");
  const bool miss = rayTriHit(Ray3{Vec3{4, 4, 5}, Vec3{0, 0, -1}},
                              {0, 0, 0}, {3, 0, 0}, {0, 3, 0}, nullptr);
  check(!miss, "三角形外不命中");
}

int main() {
  testNodePick();
  testBeamPick();
  testShellPick();
  testOcclusion();
  testRanges();
  testGeomAtoms();

  std::printf("\n====  test_pick：%d 通过，%d 失败 ====\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}