// =============================================================================
//  src/model/GridMesh.cpp  ——  轴网 → 有限元模型
//  设计说明见 include/yjk/model/GridMesh.h
// =============================================================================
#include "yjk/model/GridMesh.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <utility>

namespace yjk {

// =============================================================================
//  AxisGrid
// =============================================================================

namespace {
// 默认轴线名：X 向用 A/B/C…，Y 向用 1/2/3…
std::string defaultLabel(int k, bool letters) {
  if (!letters) return std::to_string(k + 1);
  std::string s;
  int n = k;
  do {
    s.insert(s.begin(), static_cast<char>('A' + (n % 26)));
    n = n / 26 - 1;
  } while (n >= 0);
  return s;
}
}  // namespace

Id AxisGrid::addAxisX(double y, const std::string& label) {
  const Id id = static_cast<Id>(y_.size());
  y_.push_back(y);
  yLabel_.push_back(label.empty() ? defaultLabel(static_cast<int>(id), true) : label);
  // 开关数组按 (ix*ny+jy) 索引，ny 可能变化，故尺寸在 buildPoints 时重算
  return id;
}

Id AxisGrid::addAxisY(double x, const std::string& label) {
  const Id id = static_cast<Id>(x_.size());
  x_.push_back(x);
  xLabel_.push_back(label.empty() ? defaultLabel(static_cast<int>(id), false) : label);
  return id;
}

Id AxisGrid::addStory(double z, const std::string& label) {
  const Id id = static_cast<Id>(z_.size());
  z_.push_back(z);
  zLabel_.push_back(label.empty() ? ("第" + std::to_string(id + 1) + "层") : label);
  slabStory_.push_back(id == 0 ? 0 : 1);    // 底层（基础顶/地面）默认无板
  return id;
}

// ---- 开关：数组按需扩容，全部默认"开" ----
//
// 【为什么不在 addAxisX/Y 里直接分配】轴线的两个方向是交错添加的
// （先加完 X 向再加 Y 向，或反过来），此时 ny 尚未确定，
// 无法按 (ix*ny+jy) 预分配。统一在第一次访问时按当前尺寸扩容并置 1。
namespace {
void ensure(std::vector<char>& v, size_t need, char init) {
  if (v.size() < need) v.resize(need, init);
}
}  // namespace

bool AxisGrid::setColumn(Id ix, Id jy, bool on) {
  if (ix < 0 || jy < 0 || ix >= nx() || jy >= ny()) return false;
  ensure(col_, static_cast<size_t>(nx()) * ny(), 1);
  col_[static_cast<size_t>(ix) * ny() + jy] = on ? 1 : 0;
  return true;
}
bool AxisGrid::setBeamAlongX(Id ix, Id jy, bool on) {
  if (ix < 0 || jy < 0 || ix >= nx() || jy + 1 >= ny()) return false;
  ensure(beamX_, static_cast<size_t>(nx()) * ny(), 1);
  beamX_[static_cast<size_t>(ix) * ny() + jy] = on ? 1 : 0;
  return true;
}
bool AxisGrid::setBeamAlongY(Id jy, Id ix, bool on) {
  if (jy < 0 || ix < 0 || jy >= ny() || ix + 1 >= nx()) return false;
  ensure(beamY_, static_cast<size_t>(nx()) * ny(), 1);
  beamY_[static_cast<size_t>(ix) * ny() + jy] = on ? 1 : 0;
  return true;
}
bool AxisGrid::setSlab(Id ix, Id jy, bool on) {
  if (ix < 0 || jy < 0 || ix + 1 >= nx() || jy + 1 >= ny()) return false;
  ensure(slab_, static_cast<size_t>(nx()) * ny(), 1);
  slab_[static_cast<size_t>(ix) * ny() + jy] = on ? 1 : 0;
  return true;
}
bool AxisGrid::setSlabStory(Id s, bool on) {
  if (s < 0 || s >= nStories()) return false;
  slabStory_[static_cast<size_t>(s)] = on ? 1 : 0;
  return true;
}
bool AxisGrid::setDivX(Id jy, int n) {
  if (jy < 0 || jy + 1 >= ny() || n < 1) return false;
  divX_.resize(static_cast<size_t>(ny()) - 1, 1);
  divX_[static_cast<size_t>(jy)] = n;
  return true;
}
bool AxisGrid::setDivY(Id ix, int n) {
  if (ix < 0 || ix + 1 >= nx() || n < 1) return false;
  divY_.resize(static_cast<size_t>(nx()) - 1, 1);
  divY_[static_cast<size_t>(ix)] = n;
  return true;
}

bool AxisGrid::hasColumn(Id ix, Id jy) const {
  const size_t k = static_cast<size_t>(ix) * ny() + jy;
  return k < col_.size() ? col_[k] != 0 : true;
}
bool AxisGrid::hasBeamAlongX(Id ix, Id jy) const {
  const size_t k = static_cast<size_t>(ix) * ny() + jy;
  return k < beamX_.size() ? beamX_[k] != 0 : true;
}
bool AxisGrid::hasBeamAlongY(Id jy, Id ix) const {
  const size_t k = static_cast<size_t>(ix) * ny() + jy;
  return k < beamY_.size() ? beamY_[k] != 0 : true;
}
bool AxisGrid::hasSlab(Id ix, Id jy) const {
  const size_t k = static_cast<size_t>(ix) * ny() + jy;
  return k < slab_.size() ? slab_[k] != 0 : true;
}
bool AxisGrid::hasSlabStory(Id s) const {
  return static_cast<size_t>(s) < slabStory_.size() ? slabStory_[static_cast<size_t>(s)] != 0 : true;
}
int AxisGrid::divX(Id jy) const {
  return static_cast<size_t>(jy) < divX_.size() ? divX_[static_cast<size_t>(jy)] : 1;
}
int AxisGrid::divY(Id ix) const {
  return static_cast<size_t>(ix) < divY_.size() ? divY_[static_cast<size_t>(ix)] : 1;
}

bool AxisGrid::sorted() const {
  for (size_t i = 1; i < y_.size(); ++i) if (y_[i] < y_[i - 1]) return false;
  for (size_t j = 1; j < x_.size(); ++j) if (x_[j] < x_[j - 1]) return false;
  for (size_t s = 1; s < z_.size(); ++s) if (z_[s] < z_[s - 1]) return false;
  return true;
}

// 按升序排序坐标，标签跟着一起走。
namespace {
void sortWithLabels(std::vector<double>& v, std::vector<std::string>& lab) {
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  std::stable_sort(idx.begin(), idx.end(),
                   [&](size_t p, size_t q) { return v[p] < v[q]; });
  std::vector<double> v2(v.size());
  std::vector<std::string> l2(lab.size());
  for (size_t i = 0; i < idx.size(); ++i) { v2[i] = v[idx[i]]; l2[i] = lab[idx[i]]; }
  v.swap(v2);
  lab.swap(l2);
}
}  // namespace

// 【重排会失效构件开关 —— 这是 sortAxes 最大的坑】
//
//  col_/beamX_/beamY_/slab_ 按 (ix*ny+jy) 索引，重排 X 向轴线（ix）就把它们打乱了。
//  所以：要么【先排序、再设开关】，要么接受开关被重置。
//  这里采取后者并在返回值里告知（true = 发生过重排），
//  同时把开关清空回默认（全开）—— 宁可回到默认，也不能留下错位的开关。
bool AxisGrid::sortAxes() {
  if (sorted()) return false;
  sortWithLabels(y_, yLabel_);
  sortWithLabels(x_, xLabel_);
  sortWithLabels(z_, zLabel_);
  col_.clear();
  beamX_.clear();
  beamY_.clear();
  slab_.clear();
  return true;
}

// ---- 细分后的一维坐标表 ----
//
//  xs_ = [x0, x0+Δ/div, …, x1, x1+Δ'/div', …, x_{ny-1}]
//  板与梁共用这张表 —— 这就是"板边节点与梁节点必然重合"的实现保证。
void AxisGrid::buildPoints() const {
  xs_.clear();
  xAxisIdx_.assign(static_cast<size_t>(ny()), 0);
  if (ny() == 0) return;
  xs_.push_back(x_[0]);
  for (Id j = 0; j + 1 < ny(); ++j) {
    xAxisIdx_[static_cast<size_t>(j)] = static_cast<Id>(xs_.size()) - 1;
    const int d = divX(j);
    const double a = x_[static_cast<size_t>(j)], b = x_[static_cast<size_t>(j + 1)];
    for (int k = 1; k <= d; ++k) xs_.push_back(a + (b - a) * static_cast<double>(k) / d);
  }
  xAxisIdx_[static_cast<size_t>(ny()) - 1] = static_cast<Id>(xs_.size()) - 1;

  ys_.clear();
  yAxisIdx_.assign(static_cast<size_t>(nx()), 0);
  if (nx() == 0) return;
  ys_.push_back(y_[0]);
  for (Id i = 0; i + 1 < nx(); ++i) {
    yAxisIdx_[static_cast<size_t>(i)] = static_cast<Id>(ys_.size()) - 1;
    const int d = divY(i);
    const double a = y_[static_cast<size_t>(i)], b = y_[static_cast<size_t>(i + 1)];
    for (int k = 1; k <= d; ++k) ys_.push_back(a + (b - a) * static_cast<double>(k) / d);
  }
  yAxisIdx_[static_cast<size_t>(nx()) - 1] = static_cast<Id>(ys_.size()) - 1;
}

// =============================================================================
//  GridMesh
// =============================================================================

std::string GridMesh::Result::message() const {
  if (!ok) return error.empty() ? "网格生成失败" : error;
  return "网格生成成功：节点 " + std::to_string(nodes) + "，柱 " + std::to_string(columns) +
         "，梁 " + std::to_string(beams) + "，板 " + std::to_string(slabs);
}

GridMesh::Result GridMesh::generate(Model& m) {
  Result r;
  if (grid_.nx() < 2 || grid_.ny() < 2) {
    r.error = "轴网至少需要 2 条 X 向轴线和 2 条 Y 向轴线";
    return r;
  }
  if (grid_.nStories() < 1) {
    r.error = "至少需要 1 个楼层标高";
    return r;
  }
  if (!grid_.sorted()) {
    r.error = "轴线/标高未按升序输入，请先调用 AxisGrid::sortAxes() "
              "（重排会重置构件开关，必须在设开关之前排序）";
    return r;
  }

  grid_.buildPoints();

  nXp_ = grid_.nXPoints();
  nYp_ = grid_.nYPoints();
  nSt_ = grid_.nStories();
  const size_t total = static_cast<size_t>(nXp_) * nYp_ * nSt_;
  nodes_.assign(total, -1);

  // 节点：按需创建 + 扁平索引去重。
  //
  // 【为什么用扁平索引而不是坐标哈希去重】
  // 坐标是浮点，用哈希去重需要容差和量化格，边界情况（两个点恰好跨格）
  // 很难处理；而这里每个节点都有唯一的整数网格坐标 (p,q,s)，
  // 用它做键既严格又无容差问题 —— 这正是轴网建模相对手工建模的优势。
  const Id baseNode = m.nodeCount();
  auto nodeAt = [&](Id p, Id q, Id s) -> Id {
    const size_t k = (static_cast<size_t>(p) * nYp_ + q) * nSt_ + s;
    if (nodes_[k] >= 0) return nodes_[k];
    const Id id = m.addNode({grid_.xPt(p), grid_.yPt(q), grid_.storyZ(s)}, static_cast<int>(s));
    nodes_[k] = id;
    return id;
  };

  const Id nx = grid_.nx(), ny = grid_.ny(), ns = grid_.nStories();

  // ---- 柱：交点处，从 z[s] 到 z[s+1] ----
  for (Id ix = 0; ix < nx; ++ix)
    for (Id jy = 0; jy < ny; ++jy) {
      if (!grid_.hasColumn(ix, jy)) continue;
      const Id p = grid_.xAxisIndex(jy), q = grid_.yAxisIndex(ix);
      for (Id s = 0; s + 1 < ns; ++s) {
        const Id a = nodeAt(p, q, s), b = nodeAt(p, q, s + 1);
        if (m.addBeam(a, b, spec_.columnSection, spec_.columnMaterial, spec_.columnUp))
          ++r.columns;
      }
    }

  // ---- 梁：每层楼面，沿轴线逐跨逐段生成 ----
  //
  //  【关键】沿 X 的一跨被剖成 divX(jy) 段 —— 不是一根通长梁。
  //  这样板边的每一个中间节点都落在梁的节点上，荷载才能传下去。
  for (Id s = 0; s < ns; ++s) {
    for (Id ix = 0; ix < nx; ++ix) {
      const Id q = grid_.yAxisIndex(ix);
      for (Id jy = 0; jy + 1 < ny; ++jy) {
        if (!grid_.hasBeamAlongX(ix, jy)) continue;
        const Id p0 = grid_.xAxisIndex(jy), p1 = grid_.xAxisIndex(jy + 1);
        for (Id p = p0; p < p1; ++p) {
          const Id a = nodeAt(p, q, s), b = nodeAt(p + 1, q, s);
          if (m.addBeam(a, b, spec_.beamSection, spec_.beamMaterial, spec_.beamUp))
            ++r.beams;
        }
      }
    }
    for (Id jy = 0; jy < ny; ++jy) {
      const Id p = grid_.xAxisIndex(jy);
      for (Id ix = 0; ix + 1 < nx; ++ix) {
        if (!grid_.hasBeamAlongY(jy, ix)) continue;
        const Id q0 = grid_.yAxisIndex(ix), q1 = grid_.yAxisIndex(ix + 1);
        for (Id q = q0; q < q1; ++q) {
          const Id a = nodeAt(p, q, s), b = nodeAt(p, q + 1, s);
          if (m.addBeam(a, b, spec_.beamSection, spec_.beamMaterial, spec_.beamUp))
            ++r.beams;
        }
      }
    }
  }

  // ---- 板：逐板格逐单元 ----
  //
  //  节点顺序 (p,q) → (p+1,q) → (p+1,q+1) → (p,q+1)：
  //  先沿 +X 再沿 +Y，右手系下法向 = +Z（朝上）。
  //  法向朝上很重要 —— 朝下的话板弯曲内力与面外荷载的符号全反，
  //  而面积和膜内力仍然正确，是最难发现的一类错误。
  //
  //  【不建板时这一步整个跳过】楼面荷载改由 applySlabLoad 按 45° 线导到梁上。
  //  两者是互补的两条路：建板就靠板单元传力，不建板就靠导荷。
  if (spec_.buildSlabs)
  for (Id s = 0; s < ns; ++s) {
    if (!grid_.hasSlabStory(s)) continue;
    for (Id ix = 0; ix + 1 < nx; ++ix)
      for (Id jy = 0; jy + 1 < ny; ++jy) {
        if (!grid_.hasSlab(ix, jy)) continue;
        const Id p0 = grid_.xAxisIndex(jy), p1 = grid_.xAxisIndex(jy + 1);
        const Id q0 = grid_.yAxisIndex(ix), q1 = grid_.yAxisIndex(ix + 1);
        for (Id p = p0; p < p1; ++p)
          for (Id q = q0; q < q1; ++q) {
            std::vector<Id> n4{nodeAt(p, q, s), nodeAt(p + 1, q, s),
                               nodeAt(p + 1, q + 1, s), nodeAt(p, q + 1, s)};
            if (m.addShell(n4, spec_.slab)) ++r.slabs;
          }
      }
  }

  r.nodes = m.nodeCount() - baseNode;
  r.ok = true;
  return r;
}

Id GridMesh::nodeAt(Id p, Id q, Id story) const {
  const size_t k = (static_cast<size_t>(p) * nYp_ + q) * nSt_ + story;
  return k < nodes_.size() ? nodes_[k] : -1;
}

Id GridMesh::nodeAtAxis(Id ix, Id jy, Id story) const {
  grid_.buildPoints();
  return nodeAt(grid_.xAxisIndex(jy), grid_.yAxisIndex(ix), story);
}

std::vector<Id> GridMesh::baseNodes() const {
  grid_.buildPoints();
  std::vector<Id> out;
  for (Id ix = 0; ix < grid_.nx(); ++ix)
    for (Id jy = 0; jy < grid_.ny(); ++jy) {
      if (!grid_.hasColumn(ix, jy)) continue;
      const Id n = nodeAt(grid_.xAxisIndex(jy), grid_.yAxisIndex(ix), 0);
      if (n >= 0) out.push_back(n);
    }
  return out;
}

// =============================================================================
//  刚性楼板
// =============================================================================
Id GridMesh::addRigidDiaphragms(Model& m, Id fromStory, bool coupleRz) const {
  const Id ns = grid_.nStories();
  Id made = 0;
  for (Id s = fromStory; s < ns; ++s) {
    std::vector<Id> ns_;
    for (Id i = 0; i < m.nodeCount(); ++i) {
      const Node& nd = m.node(i);
      if (nd.retired || nd.diaphragmMaster) continue;
      if (nd.story != static_cast<int>(s)) continue;
      ns_.push_back(i);
    }
    if (ns_.size() < 2) continue;
    if (m.addRigidDiaphragm(ns_, static_cast<int>(s), coupleRz) >= 0) ++made;
  }
  return made;
}

// =============================================================================
//  板面荷载导荷（45° 线法）
// =============================================================================
std::string GridMesh::SlabLoadReport::message() const {
  std::ostringstream ss;
  ss << "板格 " << bays << " 个";
  if (oneWayBays > 0) ss << "（其中单向板 " << oneWayBays << " 个）";
  if (shellElems > 0) ss << "，加压到壳单元 " << shellElems << " 块";
  if (beamSegs > 0) ss << "，导到梁上的分布荷载段 " << beamSegs << " 段";
  ss << "，总荷载 " << total << " kN";
  return ss.str();
}

GridMesh::SlabLoadReport GridMesh::applySlabLoad(Model& m,
                                                 const SlabLoadOptions& opt) const {
  SlabLoadReport rep;
  if (opt.q == 0.0) return rep;

  const bool toShells =
      (opt.path == SlabLoadOptions::Path::ToShells) ||
      (opt.path == SlabLoadOptions::Path::Auto && spec_.buildSlabs);

  // ---- 加到壳单元：就是原来的 slabload 行为 ----
  if (toShells) {
    for (auto& e : m.elements()) {
      if (e->type() != ElementType::Shell4) continue;
      auto* sh = static_cast<ShellElement*>(e.get());
      // 脚本里正数表示向下（工程习惯），单元里 -z 为向下
      sh->setTransversePressure(-opt.q);
      rep.total += opt.q * sh->area();
      ++rep.shellElems;
    }
    if (rep.shellElems == 0)
      rep.warnings.push_back("模型里没有板单元，面荷载未生效（需要 grid.slab on）");
    return rep;
  }

  // ---- 导到梁上 ----
  //
  //  节点对 → 梁单元。同一个节点对在同一层只应有一根梁；若有多根
  //  （用户手工加了重复梁），荷载会分摊到每一根上 —— 这是错的，
  //  所以要报出来。
  std::map<std::pair<Id, Id>, std::vector<Id>> beamByPair;
  for (Id e = 0; e < m.elementCount(); ++e) {
    const Element& el = *m.elements()[static_cast<size_t>(e)];
    if (el.type() != ElementType::Beam3D) continue;
    const auto& nd = el.nodes();
    if (nd.size() != 2) continue;
    beamByPair[{std::min(nd[0], nd[1]), std::max(nd[0], nd[1])}].push_back(e);
  }

  grid_.buildPoints();
  const Id nx = grid_.nx(), ny = grid_.ny(), ns = grid_.nStories();

  double onBeams = 0.0;      // 实际加到梁上的总荷载（守恒自检用）

  // 给一条边上的梁段加线性分布荷载。
  //
  //   nodes/tv: 沿边的节点序列与它们的【沿边坐标】（0 起算，m）
  //   I(t):     线荷载强度（kN/m，正 = 向下的量值）
  //   breaks:   荷载折点（沿边坐标）。梁段跨过折点时必须【切开】——
  //             否则一段线性荷载无法表示"三角形顶点"处的折角，
  //             跨中点会被拉平，跨中弯矩直接算错。
  auto addEdgeLoad = [&](const std::vector<Id>& nodes, const std::vector<double>& tv,
                         const std::vector<double>& breaks, auto I) {
    for (size_t i = 0; i + 1 < nodes.size(); ++i) {
      const Id a = nodes[i], b = nodes[i + 1];
      if (a < 0 || b < 0 || a == b) continue;
      const auto it = beamByPair.find({std::min(a, b), std::max(a, b)});
      if (it == beamByPair.end()) continue;
      if (it->second.size() > 1)
        rep.warnings.push_back("节点 " + std::to_string(a) + " 与 " +
                               std::to_string(b) + " 之间有 " +
                               std::to_string(it->second.size()) +
                               " 根梁，荷载被分摊到每一根上");
      auto* beam = dynamic_cast<BeamElement*>(m.elements()[it->second.front()].get());
      if (beam == nullptr || !beam->core().isValid()) continue;
      const double Lb = beam->core().length();
      if (!(Lb > 0.0)) continue;

      // 沿边坐标 → 梁的局部 x（0 在 nodes[i] 端）
      const double t1 = tv[i], t2 = tv[i + 1];
      const double dir = (t2 >= t1) ? 1.0 : -1.0;
      // 在【梁局部坐标】上切分
      std::vector<double> cut{0.0, Lb};
      for (double tb : breaks) {
        const double u = dir * (tb - t1);
        if (u > 1e-9 && u < Lb - 1e-9) cut.push_back(u);
      }
      std::sort(cut.begin(), cut.end());
      for (size_t k = 0; k + 1 < cut.size(); ++k) {
        const double u0 = cut[k], u1 = cut[k + 1];
        if (!(u1 > u0 + 1e-12)) continue;
        const double qa = I(t1 + dir * u0), qb = I(t1 + dir * u1);
        if (qa == 0.0 && qb == 0.0) continue;
        beam->addLineLoadSegment(u0, u1, Vec3{0, 0, -qa}, Vec3{0, 0, -qb});
        ++rep.beamSegs;
        // 【在这里就地累加，不要在事后扫单元】
        // 三角形荷载在两端强度为 0，事后扫描时"两端都非零"的过滤器
        // 会把它们整段丢掉（实测：正方形板格统计出 0 kN）。
        // 在施加点累加则天然覆盖所有段，因为这里的 qa/qb 就是实际值。
        onBeams += 0.5 * (qa + qb) * (u1 - u0);
      }
    }
  };

  for (Id s = 0; s < ns; ++s) {
    if (!grid_.hasSlabStory(s)) continue;
    for (Id ix = 0; ix + 1 < nx; ++ix)
      for (Id jy = 0; jy + 1 < ny; ++jy) {
        if (!grid_.hasSlab(ix, jy)) continue;
        const Id p0 = grid_.xAxisIndex(jy), p1 = grid_.xAxisIndex(jy + 1);
        const Id q0 = grid_.yAxisIndex(ix), q1 = grid_.yAxisIndex(ix + 1);
        const double Lx = grid_.xPt(p1) - grid_.xPt(p0);
        const double Ly = grid_.yPt(q1) - grid_.yPt(q0);
        if (!(Lx > 0.0) || !(Ly > 0.0)) continue;

        // 四边是否都有梁。缺边就无法按 45° 线导荷 ——
        // 静默跳过会让楼面荷载无声无息地消失，所以必须报出来。
        const bool bB = grid_.hasBeamAlongX(ix, jy);
        const bool bT = grid_.hasBeamAlongX(ix + 1, jy);
        const bool bL = grid_.hasBeamAlongY(jy, ix);
        const bool bR = grid_.hasBeamAlongY(jy + 1, ix);
        if (!(bB && bT && bL && bR)) {
          rep.warnings.push_back(
              "第 " + std::to_string(s) + " 层板格 (" + std::to_string(ix) + "," +
              std::to_string(jy) + ") 缺少边梁（下" + (bB ? "有" : "无") + "/上" +
              (bT ? "有" : "无") + "/左" + (bL ? "有" : "无") + "/右" +
              (bR ? "有" : "无") + "），该板格的荷载未导到梁上");
          continue;
        }

        // 沿边的节点序列与沿边坐标
        std::vector<Id> eB, eT, eL, eR;
        std::vector<double> tB, tT, tL, tR;
        for (Id p = p0; p <= p1; ++p) {
          const double t = grid_.xPt(p) - grid_.xPt(p0);
          eB.push_back(nodeAt(p, q0, s)); tB.push_back(t);
          eT.push_back(nodeAt(p, q1, s)); tT.push_back(t);
        }
        for (Id q = q0; q <= q1; ++q) {
          const double t = grid_.yPt(q) - grid_.yPt(q0);
          eL.push_back(nodeAt(p0, q, s)); tL.push_back(t);
          eR.push_back(nodeAt(p1, q, s)); tR.push_back(t);
        }

        const double Ls = std::min(Lx, Ly);
        const double Ll = std::max(Lx, Ly);
        const double peak = opt.q * Ls * 0.5;     // 峰值线荷载 kN/m
        const double half = Ls * 0.5;

        ++rep.bays;
        rep.total += opt.q * Lx * Ly;

        if (Ll / Ls >= opt.oneWayRatio) {
          // ---- 单向板：荷载全部由【长度为 L长 的两条边】承担 ----
          //
          //  板沿短边方向受力，支承它的正是垂直于短向的那两条长边。
          //  短边方向的两条边（长度 L短）不分荷载。
          ++rep.oneWayBays;
          const auto uni = [peak](double) { return peak; };
          if (Lx <= Ly) {
            addEdgeLoad(eL, tL, {}, uni);
            addEdgeLoad(eR, tR, {}, uni);
          } else {
            addEdgeLoad(eB, tB, {}, uni);
            addEdgeLoad(eT, tT, {}, uni);
          }
          continue;
        }

        // ---- 双向板：45° 线（角平分线）法 ----
        //
        //  三角形：I(t) = peak·t/half                (t ≤ half)
        //              = peak·(L−t)/half             (t > half)
        //  梯形：  I(t) = peak·t/half                (t ≤ half)
        //              = peak                        (half < t < L−half)
        //              = peak·(L−t)/half             (t ≥ L−half)
        //
        //  L 取该边自身的长度：短边（L = L短）退化成纯三角形
        //  （half = L/2 ⇒ 平台段长度为 0），长边（L = L长）才是真梯形。
        //  —— 这样写一份函数就够，不需要"这里该用三角形还是梯形"的判断。
        const double midBreak = half;
        const double farBreak = Ll - half;
        const auto trap = [peak, half, Ll](double t) {
          if (t <= half) return peak * t / half;
          if (t >= Ll - half) return peak * (Ll - t) / half;
          return peak;
        };
        const std::vector<double> brk{midBreak, farBreak};
        if (Lx <= Ly) {
          // 长度为 Lx(=L短) 的两条边（下、上）→ 三角形
          const auto triX = [peak, half, Lx](double t) {
            return (t <= half) ? peak * t / half : peak * (Lx - t) / half;
          };
          addEdgeLoad(eB, tB, {half}, triX);
          addEdgeLoad(eT, tT, {half}, triX);
          addEdgeLoad(eL, tL, brk, trap);
          addEdgeLoad(eR, tR, brk, trap);
        } else {
          const auto triY = [peak, half, Ly](double t) {
            return (t <= half) ? peak * t / half : peak * (Ly - t) / half;
          };
          addEdgeLoad(eL, tL, {half}, triY);
          addEdgeLoad(eR, tR, {half}, triY);
          addEdgeLoad(eB, tB, brk, trap);
          addEdgeLoad(eT, tT, brk, trap);
        }
      }
  }

  // ---- 守恒自检 ----
  //
  //  导到梁上的总荷载必须等于 q·Σ板格面积。浮点下不会完全相等，
  //  但相对偏差超过 1e-9 就说明分配公式写错了（曾经踩过：
  //  峰值按各边自身长度取，正方形时总量差一倍）。
  if (rep.total > 0.0) {
    const double dev = std::abs(onBeams - rep.total) / rep.total;
    if (dev > 1e-9) {
      std::ostringstream ss;
      ss << "导荷守恒自检未通过：梁上总荷载 " << onBeams << " kN，板格总荷载 "
         << rep.total << " kN，相对偏差 " << dev;
      rep.warnings.push_back(ss.str());
    }
  }
  return rep;
}

// =============================================================================
//  mergeCoincidentNodes
// =============================================================================
MergeReport mergeCoincidentNodes(Model& m, double tol) {
  MergeReport rep;
  const Id n = m.nodeCount();
  if (n < 2 || !(tol > 0.0)) return rep;

  // 量化到格：tol 为格边长。
  //
  // 【为什么不两两比较】O(n²) 在几千节点时还能接受，几万节点就不行了。
  // 量化哈希是 O(n)。代价是"恰好跨格的两点"可能被漏掉 ——
  // 但重合节点的成因是"同一位置建了两个点"，间距是 0 或 1e-15，
  // 远小于任何合理 tol，跨格风险可忽略。
  auto key = [&](const Vec3& r) -> std::array<long long, 3> {
    return {static_cast<long long>(std::llround(r.x / tol)),
            static_cast<long long>(std::llround(r.y / tol)),
            static_cast<long long>(std::llround(r.z / tol))};
  };
  // 落在格边界附近时，llround 可能把相邻点分到不同格。
  // 保险做法：同时查该点所在格与相邻的 8 个格 —— 用 ±1 偏移构造 27 个候选。
  std::map<std::array<long long, 3>, Id> seen;
  std::vector<Id> keep(static_cast<size_t>(n));
  for (Id i = 0; i < n; ++i) keep[static_cast<size_t>(i)] = i;

  for (Id i = 0; i < n; ++i) {
    const Vec3 r = m.node(i).r;
    const auto k = key(r);
    Id found = -1;
    for (int dx = -1; dx <= 1 && found < 0; ++dx)
      for (int dy = -1; dy <= 1 && found < 0; ++dy)
        for (int dz = -1; dz <= 1 && found < 0; ++dz) {
          auto c = k;
          c[0] += dx; c[1] += dy; c[2] += dz;
          auto it = seen.find(c);
          if (it == seen.end()) continue;
          const Id j = it->second;
          // 必须以真实距离复核，量化只是粗筛
          const Vec3 rj = m.node(j).r;
          if (std::abs(rj.x - r.x) <= tol && std::abs(rj.y - r.y) <= tol &&
              std::abs(rj.z - r.z) <= tol)
            found = j;
        }
    if (found >= 0) {
      keep[static_cast<size_t>(i)] = found;      // i → found（保留编号小的）
      // 荷载转移：不然加在重复节点上的力会被当成"悬空节点荷载"报错
      Node& a = m.node(found);
      Node& b = m.node(i);
      a.force = a.force + b.force;
      a.moment = a.moment + b.moment;
      b.force = Vec3{0, 0, 0};
      b.moment = Vec3{0, 0, 0};
      // 标记退役但【不删除】：结果按节点序 × 6 索引，删除会让编号整体错位。
      // 退役节点随后不被任何单元引用、荷载已清零 ⇒ 由编号器当悬空节点自动锁定。
      b.retired = true;
      ++rep.merged;
    } else {
      seen.emplace(k, i);
    }
  }

  // 重定向所有单元的节点引用
  for (auto& e : m.elements()) {
    for (Id orig : e->nodes()) {
      const Id to = keep[static_cast<size_t>(orig)];
      if (to != orig) { e->remapNode(orig, to); ++rep.redirected; }
    }
  }
  return rep;
}

}  // namespace yjk
