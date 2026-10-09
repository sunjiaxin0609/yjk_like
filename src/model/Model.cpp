// =============================================================================
//  src/model/Model.cpp  ——  模型层实现
// =============================================================================
#include "yjk/model/Model.h"

#include <algorithm>
#include <numeric>
#include <set>
#include <sstream>

namespace yjk {

namespace {
// 该节点是否被梁单元使用
//
// 【为什么必须区分梁与壳】梁单元的 6 个自由度都有物理意义
//（含绕自身轴的扭转 r_x），而膜/板单元里"绕板面法向"的转角
// 不产生任何面内应变 —— 那是零能量自由度。
// 混为一谈会导致：梁-壳混合模型里，梁的扭转自由度被误约束。
bool beamTouches(const std::vector<std::unique_ptr<Element>>& elems, Id nid) {
  for (const auto& e : elems)
    if (e->type() == ElementType::Beam3D)
      for (Id x : e->nodes())
        if (x == nid) return true;
  return false;
}

// 该节点是否被壳单元使用
bool shellTouches(const std::vector<std::unique_ptr<Element>>& elems, Id nid) {
  for (const auto& e : elems)
    if (e->type() == ElementType::Shell4)
      for (Id x : e->nodes())
        if (x == nid) return true;
  return false;
}

// 该节点是否被弹簧单元使用
// 弹簧节点只有 6 个自由度方向的"数字刚度"，与梁一样每个分量都有物理意义
//（只是可能 k=0 被释放）—— 节点属于弹簧单元就不算悬空，不能整体锁定。
bool springTouches(const std::vector<std::unique_ptr<Element>>& elems, Id nid) {
  for (const auto& e : elems)
    if (e->type() == ElementType::Spring3D)
      for (Id x : e->nodes())
        if (x == nid) return true;
  return false;
}
}  // namespace

// -----------------------------------------------------------------------------
//  节点
// -----------------------------------------------------------------------------
Id Model::addNode(const Vec3& r, int story, const std::string& label) {
  Node n;
  n.id = static_cast<Id>(nodes_.size());
  n.r = r;
  n.story = story;
  n.label = label;
  nodes_.push_back(n);
  return n.id;
}

// -----------------------------------------------------------------------------
//  单元
// -----------------------------------------------------------------------------
BeamElement* Model::addBeam(Id nI, Id nJ, const SectionProperties& sec,
                            const Material& mat, const Vec3& up) {
  if (nI < 0 || nJ < 0 || nI >= nodeCount() || nJ >= nodeCount()) return nullptr;
  if (nI == nJ) return nullptr;                       // 零长度单元
  auto e = std::make_unique<BeamElement>(nI, nJ, sec, mat, up);
  BeamElement* raw = e.get();
  raw->build(coord(nI), coord(nJ));
  elems_.push_back(std::move(e));
  return raw;
}

SpringElement* Model::addSpring(Id nI, Id nJ, const std::array<double, 6>& k,
                                const Vec3& up) {
  if (nI < 0 || nJ < 0 || nI >= nodeCount() || nJ >= nodeCount()) return nullptr;
  if (nI == nJ) return nullptr;                       // 零长度单元
  auto e = std::make_unique<SpringElement>(nI, nJ, k, up);
  SpringElement* raw = e.get();
  raw->build(coord(nI), coord(nJ));
  elems_.push_back(std::move(e));
  return raw;
}

ShellElement* Model::addShell(const std::vector<Id>& ns, const ShellProperties& props) {
  if (ns.size() != 4) return nullptr;
  for (Id n : ns)
    if (n < 0 || n >= nodeCount()) return nullptr;
  auto e = std::make_unique<ShellElement>(ns, props);
  ShellElement* raw = e.get();
  std::vector<Vec3> pts;
  pts.reserve(4);
  for (Id n : ns) pts.push_back(coord(n));
  raw->build(pts);
  elems_.push_back(std::move(e));
  return raw;
}

// -----------------------------------------------------------------------------
//  P3 交互建模：可变接口（脚本层与交互层共用）
// -----------------------------------------------------------------------------
BeamElement* Model::addColumn(Id nI, Id nJ, const SectionProperties& sec,
                              const Material& mat, const Vec3& up) {
  // 柱 = 竖向梁。addBeam 已校验节点范围/零长度并初始化几何，
  // 这里只需按柱的语义传 up（默认水平方向，见 Model.h 注释）。
  return addBeam(nI, nJ, sec, mat, up);
}

ShellElement* Model::addWall(const std::vector<Id>& ns, const ShellProperties& props) {
  ShellElement* s = addShell(ns, props);
  if (s) s->setWall(true);               // 墙标志：自重进楼层剪力
  return s;
}

void Model::addNodeLoad(Id n, const Vec3& F, const Vec3& M) {
  if (n < 0 || n >= nodeCount()) return;
  addNodeForce(n, F);
  addNodeMoment(n, M);
}

bool Model::removeElement(Id idx) {
  if (idx < 0 || idx >= elementCount()) return false;
  elems_.erase(elems_.begin() + static_cast<ptrdiff_t>(idx));
  return true;
}

bool Model::removeNode(Id n) {
  if (n < 0 || n >= nodeCount()) return false;
  const Node& nd = nodes_[static_cast<size_t>(n)];
  if (nd.retired) return false;                          // 已删除
  if (nd.diaphragmMaster) return false;                  // 楼板主节点不能删
  for (const auto& e : elems_) {                         // 引用则拒绝
    const auto& ns = e->nodes();
    if (std::find(ns.begin(), ns.end(), n) != ns.end()) return false;
  }
  nodes_[static_cast<size_t>(n)].retired = true;
  return true;
}

void Model::clear() {
  nodes_.clear();
  elems_.clear();
  links_.clear();
  diaphragms_.clear();
  numbering_ = DofNumbering();       // 重置编号策略为默认（StoryDescending）
}

// -----------------------------------------------------------------------------
//  刚性楼板（多点约束）
// -----------------------------------------------------------------------------
//
//  把整层楼的节点绑成"平面内刚体"，只留 3 个自由度：
//      UX（沿 X 平动）、UY（沿 Y 平动）、Θ（绕竖轴转动）
//
//  约束方程（节点 i 相对参考点 (x0,y0) 的偏移为 dx, dy）：
//      ux_i = UX − Θ·dy
//      uy_i = UY + Θ·dx
//      rz_i = Θ                （coupleRz = true 时）
//
//  【这三个式子的来源】平面刚体绕竖轴转 Θ，点 (dx,dy) 的位移增量是
//      Δu = Θ × r ，在二维下为 ( −Θ·dy, +Θ·dx )
//  符号搞反的后果是"扭转方向反了"，而它在对称荷载下完全看不出来 ——
//  必须用一个**偏心荷载**的用例才能验出来。
//
//  【主节点取形心】约束集与参考点位置无关（换参考点只是换几个平移量的数值），
//  所以形心是安全的、也是唯一与用户直觉一致的选择。
//  主节点只保留 (ux, uy, rz)，其余三个自由度直接约束掉 ——
//  楼板的面外刚度本来就不是无穷大，竖向弯曲仍由梁承担。
Id Model::addRigidDiaphragm(const std::vector<Id>& ns, int story, bool coupleRz) {
  Diaphragm D;
  D.story = story;
  D.coupleRz = coupleRz;

  // 形心（只统计有效节点）
  Vec3 c{0, 0, 0};
  int cnt = 0;
  for (Id n : ns) {
    if (n < 0 || n >= nodeCount()) continue;
    if (nodes_[static_cast<size_t>(n)].retired) continue;
    c = c + nodes_[static_cast<size_t>(n)].r;
    ++cnt;
  }
  if (cnt == 0) return kDofFixed;
  c = c * (1.0 / static_cast<double>(cnt));
  D.center = c;

  const Id mid = addNode(c, story, "刚性楼板主节点");
  Node& mn = nodes_[static_cast<size_t>(mid)];
  mn.diaphragmMaster = true;
  mn.fixed[2] = true;      // uz：面外刚度不由刚性楼板假定提供
  mn.fixed[3] = true;      // rx
  mn.fixed[4] = true;      // ry
  D.masterNode = mid;

  const Id mBase = mid * 6;
  for (Id n : ns) {
    if (n < 0 || n >= nodeCount() || n == mid) continue;
    if (nodes_[static_cast<size_t>(n)].retired) continue;
    const Node& nd = nodes_[static_cast<size_t>(n)];
    const double dx = nd.r.x - c.x;
    const double dy = nd.r.y - c.y;
    const Id b = n * 6;
    D.slaves.push_back(n);
    links_.push_back(DofLink{b + 0, {{mBase + 0, 1.0}, {mBase + 5, -dy}}});
    links_.push_back(DofLink{b + 1, {{mBase + 1, 1.0}, {mBase + 5, dx}}});
    if (coupleRz)
      links_.push_back(DofLink{b + 5, {{mBase + 5, 1.0}}});
  }

  diaphragms_.push_back(D);
  return mid;
}

void Model::clearLoads() {
  for (auto& n : nodes_) {
    n.force = {0, 0, 0};
    n.moment = {0, 0, 0};
  }
  for (auto& e : elems_) {
    if (auto* b = dynamic_cast<BeamElement*>(e.get())) {
      b->setSelfWeight(false);
      b->clearLineLoad();
    }
if (auto* s = dynamic_cast<ShellElement*>(e.get())) {
      s->setTransversePressure(0.0);
      s->setMembranePressure(0.0, 0.0);
      s->setSelfWeight(false);
      // 注意：isWall 是构件属性（拓扑），不属于荷载，不清除。
    }
  }
}

// -----------------------------------------------------------------------------
//  自由度编号
//
//  实现见 DofNumbering::assign。三条规则：
//    ① 楼层自上而下（story 大的先编号）
//    ② 同层内按 (y, x) 排序 —— 保证结果确定可复现
//    ③ 节点内先平动后转动
//
//  另外做两件"防呆"：
//    · 标记并自动约束壳单元节点的绕法向转角（零能量自由度）
//    · 检测"完全没有约束"的模型，给出明确警告而不是让求解器报奇异
// -----------------------------------------------------------------------------
Id DofNumbering::assign(std::vector<Node>& nodes,
                        const std::vector<std::unique_ptr<Element>>& elems,
                        const std::vector<DofLink>& links) {
  nfree_ = 0;
  nconstrained_ = 0;
  nauto_ = 0;
  nslave_ = 0;
  order_.clear();
  warn_.clear();

  const Id n = static_cast<Id>(nodes.size());
  if (n == 0) return 0;

  // ---- 第零步：标记从属自由度（多点约束的从端）----
  //
  //  【为什么在编号阶段就要知道】MPC 的从属自由度【不占方程号】：
  //  它没有自己的平衡方程，而是通过系数把贡献折到主自由度上。
  //  若照常给它编号，那一行组装出来是全零（没有任何单元刚度）→ 矩阵奇异，
  //  而求解器只会说"奇异"，看不出是约束系统的问题。
  //
  //  这里只做标记，不解算系数 —— 系数在 StaticAnalysis 组装时使用，
  //  因为那里才需要"全局自由度号 → 方程号"的映射。
  std::vector<char> slave(static_cast<size_t>(n) * 6, 0);
  for (const DofLink& L : links) {
    if (L.slave < 0 || L.slave >= static_cast<Id>(n) * 6) continue;
    if (L.masters.empty()) continue;
    slave[static_cast<size_t>(L.slave)] = 1;
  }

// ---- 第一步：判断每个节点的"绕法向转角"是否为零能量自由度 ----
  //
  //  膜/板单元内，绕板面法向（局部 z'）的转角不产生任何面内应变：
  //  B 矩阵里该转角只出现在钻转刚度项，而多数模型把钻转刚度设为 0。
  //  于是该自由度是"零能量"的 —— 组装出全零行 → 刚度矩阵奇异。
  //
  //  【关键：法向不总是全局 Z】竖直墙的壳法向沿 X 或 Y：
  //    水平板（法向 Z）→ 零能量转角是全局 rz（自由度 5）
  //    沿 X 轴线的墙（法向 Y）→ 零能量转角是全局 ry（自由度 4）
  //    沿 Y 轴线的墙（法向 X）→ 零能量转角是全局 rx（自由度 3）
  //  若像早期实现那样硬编码 rz，竖直墙的 ry/rx 既无梁柱刚度、
  //  又无自动约束 → K 零对角线 → 求解/模态分解失败。
  //
  //  【梁-壳混合模型的关键区别】梁的 6 个自由度都有物理意义
  //  （含绕自身轴的扭转 r_x），绝不能自动约束。
  //  所以判据是"该节点【只】被壳单元使用"，而不是"被壳单元使用"。
  //
  //  【法向冲突的节点不再自动约束】同一节点同时被法向 X 与法向 Y 的
  //  壳使用（如电梯井 L 形墙角）：两片墙互相给对方的面外弯曲转动提供
  //  刚度，rx 与 ry 都不是零能量 —— 只有某个转动方向对【所有】相连壳
  //  都是法向时才是全局零能量。多壳法向不一致时不做自动约束（保持自由，
  //  若确无刚度由 hasStiff 的零刚度检测兜底）。
  //
  //  【为什么不在编号时固定法向方向】法向方向依赖单元的节点顺序
  //  （逆时针才朝上），但【零能量判断只看法向所在的轴】——
  //  法向 ±X/±Y/±Z 归为同一主轴，符号不参与判定，因此顺序无关。
  std::vector<int> shellRotAxis(static_cast<size_t>(n), -1);  // -1=无壳/冲突，0/1/2=X/Y/Z
  for (Id i = 0; i < n; ++i) {
    if (!shellTouches(elems, i) || beamTouches(elems, i)) continue;
    int axis = -1;                 // 已见到的壳法向主轴
    bool conflict = false;
    for (const auto& ep : elems) {
      if (ep->type() != ElementType::Shell4) continue;
      bool touches = false;
      for (Id x : ep->nodes())
        if (x == i) { touches = true; break; }
      if (!touches) continue;
      const Vec3 nrm = ep->normal();
      int a = -1;
      if (std::abs(nrm.x) > 0.999999) a = 0;      // 法向沿 X（Y 向墙）
      else if (std::abs(nrm.y) > 0.999999) a = 1; // 法向沿 Y（X 向墙）
      else if (std::abs(nrm.z) > 0.999999) a = 2; // 法向沿 Z（水平板）
      if (a < 0) { conflict = true; break; }      // 斜壳：无法归入主轴
      if (axis < 0) axis = a;
      else if (axis != a) { conflict = true; break; }
    }
    if (!conflict && axis >= 0) shellRotAxis[static_cast<size_t>(i)] = axis;
    // 冲突或无壳 → 保持 -1：不自动约束，交给 hasStiff 零刚度检测。
  }
  for (Id i = 0; i < n; ++i) {
    if (shellRotAxis[static_cast<size_t>(i)] >= 0)
      nodes[static_cast<size_t>(i)].normalRotationActive = false;
  }

  // ---- 悬空节点（不属于任何单元）的处理 ----
  //
  //  【为什么必须在编号阶段就处理】
  //  悬空节点没有任何单元给它刚度，它的 6 个自由度在组装后是【全零行】。
  //  若照常编号，求解器只会报"刚度矩阵奇异"，用户完全看不出原因 ——
  //  而真正的原因只是"建模时忘了删掉一个多余的节点"。
  //
  //  【处理策略】无荷载的悬空节点 → 自动锁定全部 6 个自由度 + 警告。
  //  这在力学上是严格等价的：该节点无刚度、无荷载，唯一有界的解就是位移为零，
  //  锁定它不改变结构任何其它部分的结果。不静默处理 —— 会给出警告。
  //
  //  【有荷载的悬空节点不能锁定】荷载加在没有单元的节点上无处传递，
  //  锁定等于把荷载悄悄丢掉，比报错更危险。
  //  这种情况保持自由，交给 check() 的"零刚度自由度"检测报成致命错误。
std::vector<char> orphan(static_cast<size_t>(n), 0);
  for (Id i = 0; i < n; ++i) {
    if (beamTouches(elems, i) || shellTouches(elems, i) || springTouches(elems, i)) continue;
    // 刚性楼板主节点不属于任何单元，但它是整个 MPC 的支点，
    // 它的刚度是从属节点折过来的 —— 绝不能当成悬空节点锁掉。
    if (nodes[static_cast<size_t>(i)].diaphragmMaster) continue;
    const Node& nd = nodes[static_cast<size_t>(i)];
    const bool loaded = (nd.force.x != 0.0 || nd.force.y != 0.0 || nd.force.z != 0.0 ||
                         nd.moment.x != 0.0 || nd.moment.y != 0.0 || nd.moment.z != 0.0);
    if (loaded) {
      warn_.push_back("节点 " + std::to_string(i) +
                      " 不属于任何单元，却被施加了荷载 —— 该荷载无处传递，"
                      "请把它加到有单元的节点上");
      continue;
    }
    orphan[static_cast<size_t>(i)] = 1;
    // 不修改 nd.fixed[] —— 那是用户的建模输入，编号器不应反向写回。
    // 锁定在第三步按 orphan[] 完成，效果相同但模型数据保持干净。
    warn_.push_back("节点 " + std::to_string(i) +
                    " 不属于任何单元（悬空节点）—— 已自动锁定其全部自由度，"
                    "建议删除该节点");
  }

// ---- 零刚度（释放）自由度标记 ----
  //
  //  端部释放铰把释放自由度的刚度凝聚为 0；弹簧单元的 k=0 分量同理。
  //  若某 (node, comp) 被至少一个单元引用、但【所有】引用它的单元都
  //  一致释放了该分量 —— 它就是"零刚度自由度"：全局 K 该行全零 → 奇异。
  //
  //  【与悬空节点同构的处理】：
  //  无荷载 → 自动锁定（dof = -1，力学上严格等价：无刚度、无荷载的自由度
  //  唯一有界解就是位移为零）。
  //  有荷载 → 保持自由，交给 check() 的"零刚度自由度检测"报致命错误
  //  （锁定等于把荷载悄悄丢掉，比报错更危险）。
  //
  //  【为什么必须放在编号阶段而不是求解阶段】
  //  释放自由度在单元凝聚后刚度为零，但该自由度可能同时被其它未释放的
  //  单元引用（如梁端释放弯曲、但同一节点还连着柱）—— 那种情况有刚度，
  //  不能锁。只有"没有任何单元提供刚度"时才需要自动锁定。
  std::vector<char> hasStiff(static_cast<size_t>(n) * 6, 0);
  for (const auto& ep : elems) {
    const Element& el = *ep;
    const int nEl = el.localDofs();
    const int dpf = el.dofPerNode();
    const std::vector<Id>& nds = el.nodes();
    for (int a = 0; a < nEl; ++a) {
      const Id nd = nds[static_cast<size_t>(a / dpf)];
      if (nd < 0 || nd >= n) continue;
      const int comp = a % dpf;
      if (!el.isReleased(a / dpf, comp))
        hasStiff[static_cast<size_t>(nd) * 6 + static_cast<size_t>(comp)] = 1;
    }
  }

  // ---- 第二步：确定排序顺序 ----
  order_.resize(n);
  std::iota(order_.begin(), order_.end(), 0);

  if (mode_ == Mode::StoryDescending) {
    std::stable_sort(order_.begin(), order_.end(), [&](Id a, Id b) {
      const Node& A = nodes[static_cast<size_t>(a)];
      const Node& B = nodes[static_cast<size_t>(b)];
      if (A.story != B.story) return A.story > B.story;   // 楼层自上而下
      // 同层内按 (y, x)：先沿 y 后沿 x 展开
      if (std::abs(A.r.y - B.r.y) > 1e-12) return A.r.y < B.r.y;
      if (std::abs(A.r.x - B.r.x) > 1e-12) return A.r.x < B.r.x;
      return a < b;                                        // 稳定、可复现
    });
  }

  // ---- 第三步：按顺序分配编号 ----
  for (Id idx : order_) {
    Node& nd = nodes[static_cast<size_t>(idx)];
    for (int k = 0; k < 6; ++k) {
      if (nd.fixed[k] || orphan[static_cast<size_t>(idx)]) {
        nd.dof[k] = -1;
        if (!nd.fixed[k]) ++nauto_;      // 悬空节点的自动锁定
        ++nconstrained_;
        continue;
      }
// 壳单元节点的绕法向转角：若既没被梁使用、也没显式激活，则自动约束。
      // 判据不是固定 rz —— 而是该节点所有相连壳的法向主轴（第一步已算好）：
      //   水平板 → rz(5)；X 向墙 → rx(3)；Y 向墙 → ry(4)。
      // 零能量自由度 = 绕板面法向的转角，法向沿哪根主轴就约束哪个转动分量。
      // 判据"是否被梁使用"已在第一步计入：被梁使用的节点 shellRotAxis = -1。
      const int rotAxis = shellRotAxis[static_cast<size_t>(idx)];
      const int rotComp = (rotAxis == 0) ? 3 : (rotAxis == 1) ? 4
                                                               : 5;  // X→rx, Y→ry, Z→rz
      if (rotAxis >= 0 && k == rotComp && !nd.normalRotationActive) {
        nd.dof[k] = -1;
        ++nconstrained_;
        ++nauto_;
        continue;
      }
      // 从属自由度：不占方程号（见本函数开头第零步的说明）。
      // 【放在自动约束之后、分配编号之前】—— 顺序有意为之：
      // 用户显式 fixed 的优先级最高，其次是"零能量的自动约束"，
      // 最后才是 MPC。若反过来，一个已经被固定住的自由度会被 MPC
      // 抢走，用户会看到"我明明固定了它，楼层却还在动"。
if (slave[static_cast<size_t>(idx) * 6 + k]) {
        nd.dof[k] = kDofSlave;
        ++nslave_;
        continue;
      }
      // 释放（零刚度）自由度：没有任何单元为它提供刚度（所有连接到的
      // 单元都一致释放了该分量，或该分量在单元里本就是零能量壳法向）
      // → 全局 K 该行全零 → 奇异。无节点荷载时自动锁定（力学上等价：
      // 零刚度 + 零荷载的自由度，唯一有界解就是位移为零；反力恢复时
      // 该行全零 → 反力恒为零，正好符合"释放端内力为零"的物理意义）。
      // 【放在 MPC 从属之后】：被 MPC 引用的自由度即使单元释放了，
      // 也由约束方程提供行为，不能锁定。
      // 【排除刚性楼板主节点】：主节点不被单元引用，但它的刚度由
      // MPC 从属节点折合而来 —— 绝不能当成零刚度锁掉（否则楼层被钉死）。
      // 【排除壳法向已激活的转角】：normalRotationActive 时该转动是用户
      // 显式保留的自由度，由求解器/check 决定其命运，编号阶段不动它。
      const int rotC = (rotAxis == 0) ? 3 : (rotAxis == 1) ? 4 : 5;
      if (hasStiff[static_cast<size_t>(idx) * 6 + static_cast<size_t>(k)] == 0 &&
          !nd.diaphragmMaster &&
          !(rotAxis >= 0 && k == rotC && nd.normalRotationActive)) {
        // 【按分量方向判荷载】—— 不能按"节点整体有无荷载"判断。
        // 例：悬臂梁 j 端释放 ry/rz、荷载加在同节点的 uz 方向 ——
        //   uz 有刚度正常求解；ry/rz 无刚度且【对应方向无荷载】，
        //   应自动锁定（力学等价：零刚度 + 零荷载 → 唯一有界解位移为零）。
        // 只有"该零刚度分量方向本身有荷载"才保持自由，交给 check()
        // 报"无刚度"致命错误（锁定会把荷载悄悄吞掉，比报错更危险）。
        const bool loaded =
            (k < 3 ? (k == 0 ? nd.force.x : k == 1 ? nd.force.y : nd.force.z)
                   : (k == 3 ? nd.moment.x : k == 4 ? nd.moment.y : nd.moment.z)) != 0.0;
if (!loaded) {
          nd.dof[k] = -1;
          ++nconstrained_;
          ++nauto_;
          continue;
        }
      }
      nd.dof[k] = nfree_++;
    }
  }

  // ---- 第四步：防呆检查 ----
  if (nfree_ == 0) {
    warn_.push_back("所有自由度都被约束了 —— 模型没有可动自由度，无法求解");
    return 0;
  }
  if (nconstrained_ == 0) {
    warn_.push_back(
        "模型没有任何约束，刚度矩阵必然奇异。静力分析需要至少 3 个不平行的约束点。");
  }
  return nfree_;
}

bool ShellElement::plateReady() const {
  // 板弯曲需要 MITC4 剪切项；未完成时该单元只能算面内力。
  // 见 ShellElement4::plateStiffness 的说明。
  return false;
}

// -----------------------------------------------------------------------------
//  Model::check
// -----------------------------------------------------------------------------
Model::Diagnostic Model::check() const {
  Diagnostic d;
  const Id nn = nodeCount();
  const Id ne = elementCount();

  if (nn == 0) { d.errors.push_back("模型没有节点"); return d; }
  if (ne == 0) { d.errors.push_back("模型没有单元"); return d; }

  // ---- 节点：重复坐标 ----
  // 重复节点会让两个自由度通过零长度单元耦合，产生"幽灵刚度"或奇异。
  for (Id i = 0; i < nn; ++i) {
    for (Id j = i + 1; j < nn; ++j) {
      const Vec3 a = coord(i), b = coord(j);
      const double d2 = (a.x - b.x) * (a.x - b.x)
                      + (a.y - b.y) * (a.y - b.y)
                      + (a.z - b.z) * (a.z - b.z);
      // 已退役节点（被 mergeCoincidentNodes 合并掉）不参与此项检查 ——
      // 它的坐标本来就与被合并到的节点相同，再报一次"重合"是误导。
      if (node(i).retired || node(j).retired) continue;
      // 刚性楼板主节点也不参与：它【不接任何单元】，重合不产生任何耦合。
      // 而且重合是常态 —— 主节点取楼层形心，而形心往往正好落在某根
      // 柱子或某个轴网交点上（3×3 轴网 + 2×2 剖分时必然重合）。
      if (node(i).diaphragmMaster || node(j).diaphragmMaster) continue;
      if (d2 < 1e-20) {
        std::ostringstream ss;
        ss << "节点 " << i << " 与节点 " << j << " 坐标重合（间距 < 1e-10 m）";
        d.errors.push_back(ss.str());
        if (d.errors.size() > 8) { d.errors.push_back("（其余错误已省略）"); return d; }
      }
    }
  }

  // ---- 单元：几何有效性、节点引用 ----
  std::set<Id> usedNodes;
  for (Id e = 0; e < ne; ++e) {
    const Element& el = *elems_[static_cast<size_t>(e)];
    const std::string tag = "单元 " + std::to_string(e) + "（" + el.typeName() + "）";

    for (Id nid : el.nodes()) {
      if (nid < 0 || nid >= nn) {
        d.errors.push_back(tag + " 引用了不存在的节点 " + std::to_string(nid));
        continue;
      }
      usedNodes.insert(nid);
    }
    if (el.nodes().size() != static_cast<size_t>(el.nodeCount())) {
      d.errors.push_back(tag + " 的节点数不符");
    }

    if (el.type() == ElementType::Beam3D) {
      const auto* b = static_cast<const BeamElement*>(&el);
      if (!b->core().isValid()) {
        d.errors.push_back(tag + " 几何无效：" + b->error());
      } else if (b->length() < 1e-9) {
        d.errors.push_back(tag + " 长度为零");
      }
      if (b->core().usedFallbackUp()) {
        d.warnings.push_back(tag + " 的参考 up 向量与梁轴平行，已自动切换到备用方向 —— "
                                  "局部 y' 方向可能与预期不符，请检查梁方向角");
      }
      const SectionProperties& s = b->section();
      if (!(s.A > 0) || !(s.Iy > 0) || !(s.Iz > 0) || !(s.J > 0)) {
        d.errors.push_back(tag + " 的截面退化（A/Iy/Iz/J 必须全为正）");
      }
} else if (el.type() == ElementType::Shell4) {
      const auto* s = static_cast<const ShellElement*>(&el);
      if (!s->core().isValid()) {
        d.errors.push_back(tag + " 几何无效：" + s->error());
      } else if (s->area() < 1e-12) {
        d.errors.push_back(tag + " 面积为零");
      } else {
        const double n0 = s->core().normal().z;
        // 法向朝下 = 节点顺序是顺时针。
        // 后果不是"面积算错"，而是局部 z' 反向 → 板弯曲内力与面外荷载的符号全反，
        // 而面积、膜内力仍然正确 —— 所以这类错误特别难发现。
        if (n0 < 0.0) {
          d.warnings.push_back(tag + " 节点顺序为顺时针（法向朝下）—— "
                                  "建议改为逆时针，否则局部 z' 朝下，"
                                  "板弯曲内力与面外荷载的符号都会反");
        } else if (std::abs(n0 - 1.0) > 1e-6) {
          std::ostringstream ss;
          ss << tag << " 是倾斜单元（法向 z 分量 = " << n0
             << "）。注意倾斜壳单元的面外刚度不含几何刚度，"
             << "在 P-Δ 或温度作用下的结果不可靠";
          d.warnings.push_back(ss.str());
        }
      }
    } else if (el.type() == ElementType::Spring3D) {
      const auto* s = static_cast<const SpringElement*>(&el);
      if (!s->core().isValid()) {
        d.errors.push_back(tag + " 几何无效：" + s->error());
      } else if (s->length() < 1e-9) {
        d.errors.push_back(tag + " 长度为零");
      }
      // 六个方向刚度全为零 → 单元完全不提供刚度，属于建模错误
      bool allZero = true;
      for (double kd : s->core().stiffness())
        if (!nearlyZero(kd)) { allZero = false; break; }
      if (allZero) {
        d.errors.push_back(tag + " 六个方向刚度全为零 —— 单元不提供任何刚度");
      }
    }
  }

  // ---- 刚性楼板（多点约束）----
  //
  //  【为什么要单独查"主节点没有从属节点"】
  //  一个空的刚性楼板在模型里表现为"多了一个没有任何刚度的节点"，
  //  它的 ux/uy/rz 三个自由度组装出来是全零行 → 矩阵奇异。
  //  而用户看到的信息只会是"奇异"，不知道是自己把楼板节点选空了。
  for (const Diaphragm& D : diaphragms_) {
    if (D.masterNode < 0 || D.masterNode >= nn) continue;
    if (D.slaves.empty()) {
      d.errors.push_back("刚性楼板主节点 " + std::to_string(D.masterNode) +
                         " 没有任何从属节点 —— 该楼层的刚度矩阵会出现全零行");
      continue;
    }
    // 所有从属节点都被固定（或都退役）时，约束集是空的，等于没加楼板
    bool anyFree = false;
    for (Id s : D.slaves) {
      if (s < 0 || s >= nn) continue;
      if (node(s).retired) continue;
      for (const DofLink& L : links_)
        if (L.slave == s * 6 && !L.masters.empty()) { anyFree = true; break; }
      if (anyFree) break;
    }
    if (!anyFree) {
      d.warnings.push_back("刚性楼板主节点 " + std::to_string(D.masterNode) +
                           " 的从属节点全部被约束或全部退役 —— 该刚性楼板实际未生效");
    }
  }

  // ---- 悬空节点 ----
  for (Id i = 0; i < nn; ++i) {
    if (usedNodes.count(i) == 0 && !node(i).retired && !node(i).diaphragmMaster) {
      std::ostringstream ss;
      ss << "节点 " << i << " 不属于任何单元（悬空节点）";
      const auto& nd = node(i);
      if (!nd.label.empty()) ss << "：" << nd.label;
      // 无荷载的悬空节点已在编号阶段自动锁定（见 DofNumbering::assign），
      // 这里只提示用户清理模型；带荷载的悬空节点会被下面的零刚度检测报成错误。
      ss << " —— 无荷载时其自由度已被自动锁定，建议删除该节点";
      d.warnings.push_back(ss.str());
    }
  }

  // ---- 约束充足性（粗判）----
  int nPinned = 0;
  for (Id i = 0; i < nn; ++i) {
    const Node& nd = node(i);
    const int nf = nd.fixed[0] + nd.fixed[1] + nd.fixed[2]
                 + nd.fixed[3] + nd.fixed[4] + nd.fixed[5];
    if (nf >= 3) ++nPinned;
  }
  if (nPinned == 0) {
    d.errors.push_back("没有任何节点被约束 —— 刚度矩阵必然奇异，请先施加支座");
  } else if (nPinned < 3) {
    d.warnings.push_back("被约束（≥3 自由度）的节点少于 3 个，"
                        "可能存在刚体自由度未消除，求解时可能报奇异");
  }

  // ---- 自由度编号 ----
  Id nf = 0;
  for (Id i = 0; i < nn; ++i)
    for (int k = 0; k < 6; ++k)
      if (node(i).dof[k] >= 0) ++nf;
  if (nf == 0 && d.ok()) {
    d.errors.push_back("自由度编号后没有可用自由度（全部被约束？）");
  }

  // ---- 零刚度自由度检测（模型层最有价值的一项检查）----
  //
  //  【为什么必须单独检测】
  //  "自由度存在"不等于"这个自由度有刚度"。
  //  若某个自由度被声明为自由，但没有任何单元给它提供刚度，
  //  组装出的矩阵该行全零 → 刚度矩阵奇异。
  //  求解器只会报"奇异"却不说哪个自由度，
  //  而这类错误的真正原因往往在建模阶段（梁方向角、单元类型不匹配）。
  //
  //  【最常见的成因：梁的扭转自由度接不上】
  //  梁绕自身轴（局部 x'）的扭转刚度是 GJ/L。
  //  两根梁若局部 x' 方向一致，转动自由度就能互相传递；
  //  若一根沿 X、一根沿 Y，它们的扭转轴分别是 X 和 Y ——
  //  共享节点的 r_z 对第一根是"扭转"（有刚度），
  //  对第二根也是"扭转"（有刚度），可以传递；
  //  但 r_y 对第一根是【弯曲】（有刚度）、对第二根也是弯曲，
  //  看似总能传。
  //
  //  真正断链的情况：共享节点处不同单元的局部坐标系使得
  //  某个转动自由度在一个单元里是扭转（GJ/L，很小）
  //  在另一个单元里却是"几乎不受力"的方向。
  //  门式刚架的柱梁连接就是典型：
  //    柱沿 Y、up=(0,−1,0) ⇒ 局部 (x'=+Y, y'=+X, z'=−Z)
  //    梁沿 X、up=(0,0,−1) ⇒ 局部 (x'=+X, y'=−Z, z'=+Y)
  //  共享节点的 r_z：
  //    对柱是绕 z' = 扭转（GJ/L）
  //    对梁是绕 z' = 扭转（GJ/L）
  //  → 数值上能传递，但量级差 3~4 个数量级（梁的 J 比 I 小几个量级），
  //    导致该自由度"数值上几乎是零刚度"，LDL^T 分解后回代精度崩坏。
  //
  //  【检测方法】逐自由度累加其所在行的 K·（单位向量），
  //  即检查该自由度与所有其他自由度的耦合强度是否显著非零。
  {
    // 每个自由度的"刚度贡献度"：Σ_j |K(ij)|
    std::vector<double> stiff(static_cast<size_t>(nf > 0 ? nf : 1), 0.0);
    for (const auto& ep : elems_) {
      const Element& el = *ep;
      LocalMatrix K;
      if (!el.stiffnessLocal(K)) continue;
      const int nEl = el.localDofs();
      const int dpf = el.dofPerNode();
      const std::vector<Id>& nds = el.nodes();
      for (int a = 0; a < nEl; ++a) {
        const Id ga = node(nds[static_cast<size_t>(a / dpf)]).dof[a % dpf];
        if (ga < 0) continue;
        for (int b = 0; b < nEl; ++b) {
          const Id gb = node(nds[static_cast<size_t>(b / dpf)]).dof[b % dpf];
          if (gb < 0) continue;
          const double v = K[static_cast<size_t>(a) * nEl + b];
          stiff[static_cast<size_t>(ga)] += std::abs(v);
        }
      }
    }
    // 全局最大刚度作为量级基准
    double kmax = 0.0;
    for (double v : stiff) kmax = std::max(kmax, v);
    if (kmax > 0.0) {
      // ---- 判据 1：完全无刚度 ----
      const double zeroTol = 1e-14 * kmax;
      // ---- 判据 2：刚度比过大（条件数预警）----
      //
      //  【为什么必须查这一项】
      //  梁单元的"弯曲刚度"（EI/L ≈ 1e7）与"扭转刚度"（GJ/L ≈ 1e5）
      //  相差 2~3 个数量级；再叠加构件长度差异（柱 0.5m vs 梁 1m），
      //  同一根梁上不同自由度的刚度量级可以差 4~5 个数量级。
      //  这在数学上是合法的（刚度矩阵仍然正定），
      //  但会让 LDL^T 的回代精度显著下降 —— 表现是
      //  "求解不报错，残差却到 1e-1 量级"。
      //
      //  门式刚架柱梁连接处就是这种情形：
      //  共享节点的 r_z 对柱是"绕 z' 的扭转"、对梁是"绕 z' 的扭转"，
      //  两者能传递，但扭转刚度 GJ/L 远小于弯曲刚度 EI/L，
      //  于是该自由度在行消去后几乎被消掉，回代时精度崩坏。
      //
      //  经验阈值：刚度比超过 1e8 就值得提醒用户检查梁方向角。
      const double warnRatio = 1e8;
      double kmin = kmax;
      for (Id g = 0; g < nf; ++g) kmin = std::min(kmin, stiff[static_cast<size_t>(g)]);

      int nZero = 0, nWeak = 0;
      for (Id i = 0; i < nn; ++i)
        for (int k = 0; k < 6; ++k) {
          const Id g = node(i).dof[k];
          if (g < 0) continue;
          // 刚性楼板主节点不接任何单元，它的刚度【全部】来自从属节点
          // 通过 MPC 折过来的贡献。用"单元行"去量它必然是零 ——
          // 那不是问题，只是这个量法对它不适用。
          // 真正要防的是"主节点一个从属节点都没有"，已在上面的楼板检查里报了。
          if (node(i).diaphragmMaster) continue;
          const double kg = stiff[static_cast<size_t>(g)];
          if (kg <= zeroTol) {
            std::ostringstream ss;
            ss << "节点 " << i << " 的自由度 "
               << (k == 0 ? "ux" : k == 1 ? "uy" : k == 2 ? "uz"
                       : k == 3 ? "rx" : k == 4 ? "ry" : "rz")
               << " 无刚度（没有任何单元对它提供约束）";
            d.errors.push_back(ss.str());
            ++nZero;
          } else if (kg / kmax < 1.0 / warnRatio) {
            ++nWeak;
          }
        }
      if (nWeak > 0) {
        std::ostringstream ss;
        ss << "有 " << nWeak << " 个自由度的刚度比全模型最大刚度低 "
           << int(std::log10(warnRatio)) << " 个数量级以上。"
           << "这会让求解精度显著下降（残差可能到 1e-1 量级）。"
           << "常见原因：梁的方向角设置使不同构件的局部坐标系不一致，"
           << "导致共享节点的转动自由度难以互相传递。"
           << "建议检查梁的『截面竖向』方向是否统一。";
        d.warnings.push_back(ss.str());
      }
    }
  }

  return d;
}

std::string Model::Diagnostic::summary() const {
  std::ostringstream ss;
  ss << "错误 " << errors.size() << " 项，警告 " << warnings.size() << " 项";
  return ss.str();
}

}  // namespace yjk
