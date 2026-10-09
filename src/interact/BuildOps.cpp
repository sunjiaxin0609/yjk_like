// =============================================================================
//  src/interact/BuildOps.cpp  ——  交互建模操作实现
//
//  所有操作都调用 Model 的 T1 可变接口，本层只负责：
//    · 校验输入（节点存在且未 retired、非零长度、四点不退化）
//    · 把失败原因整理成结构化中文消息（GUI 直接展示，测试逐字断言）
//    · 成功的 ok=true 携带新建单元 id
// =============================================================================
#include "yjk/interact/BuildOps.h"

#include <cmath>

namespace yjk {
namespace interact {

namespace {

// 节点 id 是否可被引用：范围内且未被 mergeCoincidentNodes 删除。
bool nodeUsable(const Model& m, Id n) {
  if (n < 0 || n >= m.nodeCount()) return false;
  return !m.node(n).retired;
}

// 梁/柱共用：校验两点 → 创建。kind 只决定错误消息文案。
BuildOpResult placeBar(Model& m, Id nI, Id nJ, const SectionProperties& sec,
                       const Material& mat, const Vec3& up, const char* what) {
  BuildOpResult r;
  if (!nodeUsable(m, nI) || !nodeUsable(m, nJ)) {
    r.message = std::string(what) + "：节点不存在（已被删除或越界）";
    return r;
  }
  if (nI == nJ) {
    r.message = std::string(what) + "：两个节点不能是同一个（零长度单元）";
    return r;
  }
  BeamElement* e = m.addBeam(nI, nJ, sec, mat, up);
  if (!e) {
    r.message = std::string(what) + "：内核拒绝创建（端点越界或长度为零）";
    return r;
  }
  r.ok = true;
  r.id = static_cast<Id>(m.elementCount() - 1);
  r.message = std::string(what) + "已创建（#" + std::to_string(r.id) + "）";
  return r;
}

// 壳四点是否退化：要求两对角线三角形面积均非零。
// 注意【不做共面检查】—— 脚本 wall/shell 允许翘曲四边形（四个角点
// 不严格共面），ModelScript 与 ShellElement4 都接受，BuildOps 若要求
// 共面会拒绝内核/脚本同意的形状，破坏"GUI 与 .yjk 脚本同操作等价"。
// 只拦真正的退化：三点共线或重合（面积≈0 → 刚度奇异）。
bool shellQuadOk(const Model& m, const std::array<Id, 4>& ns) {
  const Vec3& a = m.coord(ns[0]);
  const Vec3& b = m.coord(ns[1]);
  const Vec3& c = m.coord(ns[2]);
  const Vec3& d = m.coord(ns[3]);
  // 三角形 (a,b,c) 的叉积面积 2 倍
  const Vec3 u = b - a;
  const Vec3 v = c - a;
  const Vec3 n1 = Vec3{u.y * v.z - u.z * v.y,
                       u.z * v.x - u.x * v.z,
                       u.x * v.y - u.y * v.x};
  const double s1 = n1.x * n1.x + n1.y * n1.y + n1.z * n1.z;
  if (s1 < 1e-24) return false;              // 前三线共线或重合
  // 三角形 (a,c,d)
  const Vec3 w = d - a;
  const Vec3 n2 = Vec3{v.y * w.z - v.z * w.y,
                       v.z * w.x - v.x * w.z,
                       v.x * w.y - v.y * w.x};
  const double s2 = n2.x * n2.x + n2.y * n2.y + n2.z * n2.z;
  return !(s2 < 1e-24);
}

}  // namespace

BuildOpResult placeBeam(Model& m, Id nI, Id nJ, const SectionProperties& sec,
                        const Material& mat, const Vec3& up) {
  return placeBar(m, nI, nJ, sec, mat, up, "梁");
}

BuildOpResult placeColumn(Model& m, Id nI, Id nJ, const SectionProperties& sec,
                          const Material& mat) {
  return placeBar(m, nI, nJ, sec, mat, Vec3{0, -1, 0}, "柱");
}

BuildOpResult placeWall(Model& m, const std::array<Id, 4>& ns,
                        const ShellProperties& props) {
  BuildOpResult r;
  for (Id n : ns) {
    if (!nodeUsable(m, n)) {
      r.message = "墙：节点不存在（已被删除或越界）";
      return r;
    }
  }
  // 四点两两不同（重合会在 shellQuadOk 里以"面积为零"拦下，这里提前给明确提示）。
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      if (ns[i] == ns[j]) {
        r.message = "墙：四角节点中有重复（请依次点选四个角点）";
        return r;
      }
  if (!shellQuadOk(m, ns)) {
    r.message = "墙：四个节点共线或不在同一平面（无法成形）";
    return r;
  }
  ShellElement* e = m.addWall(std::vector<Id>(ns.begin(), ns.end()), props);
  if (!e) {
    r.message = "墙：内核拒绝创建";
    return r;
  }
  r.ok = true;
  r.id = static_cast<Id>(m.elementCount() - 1);
  r.message = "墙已创建（#" + std::to_string(r.id) + "）";
  return r;
}

BuildOpResult placeSlabFromNodes(Model& m, const std::vector<Id>& boxed,
                                 const ShellProperties& props) {
  BuildOpResult r;
  if (boxed.size() < 4) {
    r.message = "板：框选范围内节点不足（至少需要 4 个，当前 "
                + std::to_string(boxed.size()) + " 个）";
    return r;
  }
  // 去重 + 校验存在
  std::vector<Id> ids;
  for (Id n : boxed) {
    if (!nodeUsable(m, n)) continue;
    bool dup = false;
    for (Id x : ids)
      if (x == n) { dup = true; break; }
    if (!dup) ids.push_back(n);
  }
  if (ids.size() < 4) {
    r.message = "板：有效节点不足 4 个（其余已删除或重复）";
    return r;
  }

  // 取 (x, y) 包围盒四角的最近节点 —— 框选得到的是一簇节点，
  // "四角最近"给出板的外廓，与脚本里手工点四角建 shell 语义一致。
  double xMin = m.coord(ids[0]).x, xMax = xMin;
  double yMin = m.coord(ids[0]).y, yMax = yMin;
  for (Id n : ids) {
    const Vec3& p = m.coord(n);
    xMin = std::min(xMin, p.x); xMax = std::max(xMax, p.x);
    yMin = std::min(yMin, p.y); yMax = std::max(yMax, p.y);
  }
  const Vec3 corners[4] = {
      {xMin, yMin, 0}, {xMax, yMin, 0}, {xMax, yMax, 0}, {xMin, yMax, 0}};

  std::array<Id, 4> pick;
  for (int k = 0; k < 4; ++k) {
    Id best = -1;
    double bestD = 1e30;
    for (Id n : ids) {
      const Vec3& p = m.coord(n);
      const double d = (p.x - corners[k].x) * (p.x - corners[k].x) +
                       (p.y - corners[k].y) * (p.y - corners[k].y);
      if (d < bestD) { bestD = d; best = n; }
    }
    pick[static_cast<size_t>(k)] = best;
  }
  // 四角节点必须两两不同，且壳不退化
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      if (pick[static_cast<size_t>(i)] == pick[static_cast<size_t>(j)]) {
        r.message = "板：框选节点不构成四角（多个角点最近到同一节点）";
        return r;
      }
  if (!shellQuadOk(m, pick)) {
    r.message = "板：四角节点共线或不在同一平面（无法成形）";
    return r;
  }
  ShellElement* e = m.addShell(std::vector<Id>(pick.begin(), pick.end()), props);
  if (!e) {
    r.message = "板：内核拒绝创建";
    return r;
  }
  r.ok = true;
  r.id = static_cast<Id>(m.elementCount() - 1);
  r.message = "板已创建（#" + std::to_string(r.id) + "）";
  return r;
}

BuildOpResult placeNodeLoad(Model& m, Id n, const Vec3& F, const Vec3& M) {
  BuildOpResult r;
  if (!nodeUsable(m, n)) {
    r.message = "节点荷载：节点不存在（已被删除或越界）";
    return r;
  }
  m.addNodeLoad(n, F, M);
  r.ok = true;
  r.id = n;
  r.message = "节点荷载已施加到节点 #" + std::to_string(n);
  return r;
}

}  // namespace interact
}  // namespace yjk