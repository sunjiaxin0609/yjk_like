// =============================================================================
//  include/yjk/model/Model.h  ——  模型层：节点 / 单元 / 荷载 / 约束 / 自由度编号
//
//  这一层回答三个问题：
//    ① 模型里有什么（节点、单元、截面、材料、荷载、约束）
//    ② 自由度怎么编号（直接决定稀疏求解的填充量 —— 见 DofNumbering 说明）
//    ③ 单元荷载怎么等效到节点
//
//  【设计原则：模型层不碰数值算法】
//  它只负责"把工程数据翻译成局部刚度矩阵 + 节点荷载向量"，
//  全局组装、约束消去、求解、后处理都在 analysis 层。
//  这样单元实现可以独立演进，不影响模型组织。
//
//  【自由度映射的唯一约定】
//  单元刚度矩阵永远是"单元局部自由度顺序"（逐节点 dpf 个）。
//  它与全局自由度的关系由 node.dof[] 决定，组装时一次性转换。
//  单元内部永远不知道自己对应的全局号 —— 这样单元可以任意重编号。
// =============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "yjk/element/BeamElement3D.h"
#include "yjk/element/ShellElement4.h"
#include "yjk/element/SpringElement3D.h"
#include "yjk/math/SparseMatrix.h"
#include "yjk/model/Section.h"

namespace yjk {

enum class ElementType { Beam3D, Shell4, Spring3D };

// 单元局部刚度的存储：行主序方阵，n = localDofs()
using LocalMatrix = std::vector<double>;

// -----------------------------------------------------------------------------
//  Node::dof[] 的取值语义
//
//  【为什么要区分 -1 和 -2】
//  两者都"不参与求解"，但反力恢复的方式完全不同：
//    · 固定自由度：反力由"该行的原方程"恢复，位移是已知量（0 或支座沉降）
//    · 从属自由度：位移由主自由度【算出来】，该行方程给出的其实是
//      "约束装置施加的力"（楼板把节点拽住的那份力）
//  把两者混为一谈会得到两种典型错误：从属自由度位移恒为 0（楼层被"钉住"），
//  或者约束力被当成支座反力加进平衡校核（ΣR = −ΣF 不成立）。
// -----------------------------------------------------------------------------
constexpr Id kDofFixed = -1;     // 被约束：用户固定 / 自动锁定 / 悬空节点
constexpr Id kDofSlave = -2;     // 从属于主自由度（多点约束 MPC）

// -----------------------------------------------------------------------------
//  多点约束（MPC）：从属自由度 = Σ 系数 × 主自由度
//
//      u_slave = Σ_m  c_m · u_master(m)
//
//  【为什么用"从属"而不是"罚函数"】
//  罚函数把一个大数加到对角线上去"压住"位移，代价有两个：
//    · 条件数被污染 —— 而本项目数学层已经吃过这个亏：
//      条件数一坏，稀疏 LDL^T 的回代精度就崩（残差从 1e-13 掉到 1e-1）
//    · 反力无法精确恢复 —— 罚项自己也会"出力"，分不清哪部分是约束力
//  从属约束是【精确】的：通过变换矩阵 T 做 K_r = Tᵀ·K·T、F_r = Tᵀ·F，
//  原方程一条不丢，约束力仍能从"被消去那一行的原方程"精确恢复。
//
//  【索引约定】slave / masters 里存的是【全局自由度号】= node*6 + comp，
//  不是方程号。因为 assignDofs() 可以被反复调用（换编号策略），
//  把方程号固化进约束里会让两者失去同步 —— 那是极难查的一类错。
// -----------------------------------------------------------------------------
struct DofLink {
  Id slave{kDofFixed};
  std::vector<std::pair<Id, double>> masters;   // 全局自由度号 → 系数
};

// -----------------------------------------------------------------------------
//  节点
// -----------------------------------------------------------------------------
struct Node {
  Id id{0};
  Vec3 r{0, 0, 0};             // 全局坐标 m

// 节点荷载（由用户逐个添加，此处累加）
  Vec3 force{0, 0, 0};        // kN
  Vec3 moment{0, 0, 0};       // kN·m

  // 节点附加重量（kN，重力荷载代表值折算用，见 gmass 命令）。
  // 与 force 独立：clearLoads() 不清除它；ModalAnalysis 组装质量时
  // 检测到任一节点 weight > 0 就改用"节点重量代表团"路径（自重由
  // 恒载工况的 selfweight 体现，避免重复计重）。
  double weight{0.0};

  // 约束：true = 该自由度被固定
  bool fixed[6]{false, false, false, false, false, false};

  // 建模属性
  int story{0};               // 所属楼层（0 = 地面或未指定）
  std::string label;

  // 已被 mergeCoincidentNodes 合并掉。
  //
  //  【为什么不直接删】结果向量按【节点序 × 6】索引，删掉一个节点会让
  //  后面所有节点号整体错位，历史结果与界面选中状态全部失效。
  //  保留 + 标记是等价且安全的做法：它不再被任何单元引用、荷载已转走、
  //  自由度由编号器自动锁定，等价于"不存在"，但节点号保持稳定。
  bool retired{false};

  // 绕板面法向的转角（rz）是否"激活"。
  // 膜单元内该转动不产生面内应变 → 若不约束则是零能量自由度 → 刚度矩阵奇异。
  // 编号器会据此自动约束（除非用户显式要求激活）。
  bool normalRotationActive{false};

  // 刚性楼板的【主节点】。
  //
  //  【它为什么需要单独的标记】主节点不属于任何单元 —— 它的刚度
  //  完全来自"从属节点通过 MPC 折过来的贡献"。若不标记，
  //  它会被悬空节点逻辑自动锁定、并被零刚度检测报成致命错误，
  //  而实际上它是整个刚性楼板假定的支点。
  bool diaphragmMaster{false};

  // 全局自由度号，-1 = 被约束或未激活。由 DofNumbering 填写。
  Id dof[6]{-1, -1, -1, -1, -1, -1};

  // 【位移向量的索引约定 —— 整个结果体系必须统一】
  //
  //  u / reaction 都按【节点序 × 6】索引：节点 n 的第 k 个分量在 u[n*6+k]。
  //  不用 dof[] 里的重编号号索引，理由是：
  //    · 后处理、结果文件、界面显示全都按"节点 → 分量"访问，用节点序最直观
  //    · dof[] 是模型内部的求解细节，把它泄漏到结果结构里，
  //      一旦编号策略改变（换 numbering 模式）历史结果就无法复用
  //  代价是每次访问要多一层 dof[] 查表 —— 可以忽略。
  Id baseDof() const { return id * 6; }

  Vec3 displacement(const std::vector<double>& u) const {
    Vec3 d{0, 0, 0};
    for (int k = 0; k < 3; ++k) d[k] = u[static_cast<size_t>(id * 6 + k)];
    return d;
  }
  Vec3 rotation(const std::vector<double>& u) const {
    Vec3 d{0, 0, 0};
    for (int k = 0; k < 3; ++k) d[k] = u[static_cast<size_t>(id * 6 + 3 + k)];
    return d;
  }
  // 单个分量（k ∈ [0,6)）
  double u(const std::vector<double>& u, int k) const {
    return u[static_cast<size_t>(id * 6 + k)];
  }
};

// -----------------------------------------------------------------------------
//  单元基类
//
//  抽象出三件事：局部刚度、等效节点荷载、内力恢复。
//  上层组装代码因此不关心具体单元类型 —— 新增单元只需再实现一个子类。
// -----------------------------------------------------------------------------
class Element {
 public:
  virtual ~Element() = default;

// ---- 拓扑 ----
  virtual int nodeCount() const = 0;
  virtual int dofPerNode() const = 0;                 // 梁 6，壳 6
  virtual const std::vector<Id>& nodes() const = 0;
  int localDofs() const { return nodeCount() * dofPerNode(); }

  // 深拷贝（T3 内存重算快照用）。
  // 【为什么必须有】Model 持 unique_ptr<Element>，默认拷贝构造被禁用；
  // 而"GUI 编辑后重算"必须在后台线程拿到一份独立的模型副本
  // （assignDofs 会写 nodes_[].dof[]，与 UI 线程并发读写同一模型是数据竞争）。
  // 各子类成员都是值语义（BeamElement3D 等含普通数组，无堆资源），
  // 用 make_unique<X>(*this) 即可完成逐成员拷贝。
  virtual std::unique_ptr<Element> clone() const = 0;

  // 可变引用 —— 仅供【拓扑修复】使用（见 mergeCoincidentNodes）。
  // 正常建模路径不应改单元的拓扑：改了就等于换了一个单元，
  // 而单元几何（BeamElement3D / ShellElement4）是在 addBeam/addShell 时
  // 按当时的节点坐标初始化的，改节点号不会同步几何 —— 会静默出错。
  virtual std::vector<Id>& nodesMutable() = 0;

  // 把所有指向 from 的引用改成 to。配合 mergeCoincidentNodes 使用。
  void remapNode(Id from, Id to) {
    for (Id& n : nodesMutable())
      if (n == from) n = to;
  }

  // ---- 力学 ----
  // 局部刚度矩阵（行主序，n×n，对称）。返回 false 表示单元无效。
  virtual bool stiffnessLocal(LocalMatrix& K) const = 0;

  // 单元荷载等效到【单元局部自由度】。数组长度 = localDofs()，调用前先清零。
  virtual void equivalentLoads(std::vector<double>& fLocal) const = 0;

  // ---- 元信息 ----
  virtual ElementType type() const = 0;
  virtual std::string typeName() const = 0;
  virtual double length() const { return 0.0; }
  virtual double mass() const { return 0.0; }
  virtual Vec3 normal() const { return {0, 0, 1}; }

  // ---- 端部释放 / 零刚度自由度查询 ----
  //
  // 返回单元的某节点端、某自由度（localNode = 局部节点号，comp = 节点内自由度）
  // 是否被【释放】。释放意味着该单元对这两个全局自由度【不提供任何刚度】。
  //
  // 用途：DofNumbering::assign 在做"零刚度自由度自动锁定"时需要知道
  // 某 (node, comp) 是否被全部引用它的单元一致释放 —— 若是且无荷载，
  // 该自由度就该被自动锁定（否则全局 K 奇异）。
  //
  // 默认实现返回 false（普通单元不释放任何自由度）。
  virtual bool isReleased(int localNode, int comp) const {
    (void)localNode; (void)comp;
    return false;
  }

  // ---- 后处理 ----
  //
  //  杆端内力 = k·u − f_eq。
  //  fEq 是 Element::equivalentLoads() 的输出（全局分量，长度 = localDofs()）；
  //  传空数组即退化为纯 k·u（单元测试对照用）。
  //
  //  【为什么必须显式传 fEq 而不是让单元自己算】
  //  单元自己当然能算，但那样每次恢复内力都要重算一遍等效荷载，
  //  而 StaticAnalysis 在组装阶段已经算过一次了 —— 更重要的是，
  //  显式传参让"有没有做固端修正"在调用点一眼可见，
  //  避免后处理悄悄拿 k·u 当真内力用（这是本项目踩过的坑）。
  virtual BeamElement3D::EndForces endForces(const std::vector<double>& uEl,
                                             const std::vector<double>& fEq) const {
    (void)uEl; (void)fEq;
    return {};                                  // 非梁单元没有杆端内力
  }
  // 便捷重载：不做固端修正（= 纯 k·u）
  BeamElement3D::EndForces endForces(const std::vector<double>& uEl) const {
    return endForces(uEl, std::vector<double>{});
  }
  virtual ShellElement4::MembraneForces shellMembraneForces(const std::vector<double>&) const {
    return {};                                  // 非壳单元没有膜内力
  }
};

// -----------------------------------------------------------------------------
//  梁单元（包装 BeamElement3D）
// -----------------------------------------------------------------------------
class BeamElement : public Element {
 public:
  BeamElement(Id nI, Id nJ, const SectionProperties& sec, const Material& mat,
              const Vec3& up = Vec3{0, 0, -1})
      : nodes_{nI, nJ}, up_(up), sec_(sec), mat_(mat) {}

  bool build(const Vec3& pi, const Vec3& pj) {
    be_.init(pi, pj, sec_, mat_, up_);
    return be_.isValid();
  }
const std::string& error() const { return be_.error(); }
  const BeamElement3D& core() const { return be_; }
  const SectionProperties& section() const { return sec_; }
  const Material& material() const { return mat_; }
  Vec3 upHint() const { return up_; }
  void setUpHint(const Vec3& v) { up_ = v; }

  // T3 属性编辑：改截面与材料。只改成员；调用方（Model::setBeamProps）
  // 负责随后重建几何。be_.init 不触碰 relI_/relJ_（端部释放）与
  // segs_/pt_/selfWeightOn_（荷载），故重 build 安全。
  void setSectionAndMaterial(const SectionProperties& sec, const Material& mat) {
    sec_ = sec;
    mat_ = mat;
  }

  std::unique_ptr<Element> clone() const override {
    // BeamElement 成员全部值语义（BeamElement3D be_ 内部是普通数组），
    // 默认拷贝构造即正确深拷贝。
    return std::make_unique<BeamElement>(*this);
  }

// ---- Element ----
  int nodeCount() const override { return 2; }
  int dofPerNode() const override { return 6; }
  const std::vector<Id>& nodes() const override { return nodes_; }
  std::vector<Id>& nodesMutable() override { return nodes_; }
  ElementType type() const override { return ElementType::Beam3D; }
  std::string typeName() const override { return "梁单元 3D"; }
  double length() const override { return be_.length(); }
  double mass() const override { return be_.mass(); }
  Vec3 normal() const override { return be_.localZ(); }
  bool isReleased(int localNode, int comp) const override {
    return be_.isReleased(localNode, comp);
  }

  // ---- 端部释放（端部铰 / 选择性释放）----
  // localNode: 0 = i 端, 1 = j 端；comp 0..5 = 局部 ux,uy,uz,rx,ry,rz
  void releaseAt(int localNode, int comp) { be_.releaseAt(localNode, comp); }
  void setReleases(const std::array<bool, 6>& iRel, const std::array<bool, 6>& jRel) {
    be_.setReleases(iRel, jRel);
  }
  bool hasReleases() const { return be_.hasReleases(); }
  const std::string& condensationError() const { return be_.condensationError(); }

  bool stiffnessLocal(LocalMatrix& K) const override {
    if (!be_.isValid()) return false;
    const Mat12 ke = be_.stiffness();
    K.resize(144);
    for (int r = 0; r < 12; ++r)
      for (int c = 0; c < 12; ++c) K[static_cast<size_t>(r) * 12 + c] = ke(r, c);
    return true;
  }

void equivalentLoads(std::vector<double>& f) const override {
    f.assign(12, 0.0);
    if (!be_.isValid()) return;
    // 【分布荷载只有一条路径】均布 = "一段全长、两端同强度"的特例。
    // 分开写成两条实现，它们早晚会各自演化 ——
    // 而"均布对、梯形错"这种不一致最难发现（测试用例大多只用均布）。
    for (const auto& s : distributedLoads()) {
      const Vec12d v = be_.equivalentNodalLoadGlobalSegment(s.x1, s.x2, s.q1, s.q2);
      for (int i = 0; i < 12; ++i) f[i] += v[i];
    }
    for (size_t i = 0; i < pt_.size() && i < f.size(); ++i) f[i] += pt_[i];

    // 【端部释放】等效荷载必须与刚度同步凝聚，否则 b 分量的固端弯矩/剪力
    // 会被 scatter 进已锁定的释放自由度 → 反力恢复出现虚假项。
    // 流程：全局 f → 局部 → 凝聚（b 分量清零）→ 回全局。
    if (be_.hasReleases()) {
      Vec12d fgv{};
      for (int i = 0; i < 12 && i < static_cast<int>(f.size()); ++i) fgv[i] = f[i];
      const Vec12d fl = be_.condenseLoad(be_.toLocalVector(fgv));
      const Vec12d fgc = be_.toGlobalForce(fl);
      for (int i = 0; i < 12; ++i) f[i] = fgc[i];
    }
  }

  // 沿梁长的截面内力（含分布荷载的真实积分，见 BeamElement3D::atStation）。
  BeamElement3D::EndForces atStation(const Vec12d& u, const Vec12d& fEq, double xi) const {
    const std::vector<BeamElement3D::DistributedLoad> q = distributedLoads();
    return be_.atStation(u, fEq, xi, &q);
  }

  using Element::endForces;
  BeamElement3D::EndForces endForces(const std::vector<double>& uEl,
                                     const std::vector<double>& fEq) const override {
    Vec12d ue{}, fe{};
    for (int i = 0; i < 12 && i < static_cast<int>(uEl.size()); ++i) ue[i] = uEl[i];
    for (int i = 0; i < 12 && i < static_cast<int>(fEq.size()); ++i) fe[i] = fEq[i];
    return be_.endForces(ue, fe);
  }

  // ---- 荷载设置 ----
  void setSelfWeight(bool on) { selfWeightOn_ = on; }
  bool selfWeight() const { return selfWeightOn_; }

  // 均布线荷载（全局分量）。等价于"一段全长、两端同强度"。
  void setLineLoad(const Vec3& q) {
    segs_.clear();
    segs_.push_back(BeamElement3D::DistributedLoad{0.0, be_.length(), q, q});
  }
  void clearLineLoad() { segs_.clear(); }
  bool hasLineLoad() const { return !segs_.empty(); }
  Vec3 lineLoad() const {
    return segs_.empty() ? Vec3{0, 0, 0} : segs_.front().q1;
  }

  // 线性变化的分段荷载 —— 双向板 45° 导荷产生的三角形/梯形就靠它。
  // x1/x2 是【局部 x】（0 在 i 端），q1/q2 是全局分量。
  void addLineLoadSegment(double x1, double x2, const Vec3& q1, const Vec3& q2) {
    segs_.push_back(BeamElement3D::DistributedLoad{x1, x2, q1, q2});
  }
  const std::vector<BeamElement3D::DistributedLoad>& loadSegments() const { return segs_; }

  // 全部分布荷载（含自重）。atStation 与 equivalentLoads 共用同一份描述，
  // 这是"内力图与固端修正必须一致"的前提。
  std::vector<BeamElement3D::DistributedLoad> distributedLoads() const {
    std::vector<BeamElement3D::DistributedLoad> out = segs_;
    if (selfWeightOn_) {
      const double w = mat_.gamma * sec_.A;
      const Vec3 q{0, 0, -w};
      out.push_back(BeamElement3D::DistributedLoad{0.0, be_.length(), q, q});
    }
    return out;
  }

  // 端部集中力 / 集中力矩（原样传递到指定端，不产生附加弯矩）
  void addPointForce(int end, const Vec3& F) {
    if (end != 0 && end != 1) return;
    appendPt(be_.concentratedForceGlobal(end, F));
  }
  void addPointMoment(int end, const Vec3& M) {
    if (end != 0 && end != 1) return;
    appendPt(be_.concentratedMoment(end, M));
  }
  // 单元内部集中力（按比例分配到两端，无附加弯矩 —— 与"梁段"等效）
  void addMidPointForce(const Vec3& F) {
    Vec12d f{};
    f[0] = f[3] = 0.5 * F.x;
    f[1] = f[4] = 0.5 * F.y;
    f[2] = f[5] = 0.5 * F.z;
    f[6] = f[9] = 0.5 * F.x;
    f[7] = f[10] = 0.5 * F.y;
    f[8] = f[11] = 0.5 * F.z;
    appendPt(f);
  }
bool hasPointLoad() const { return !pt_.empty(); }
  const std::vector<double>& pointLoads() const { return pt_; }
  // T6 全量序列化：等效节点集中力（12 个值）原样写回，与读取互逆。
  void setPointLoads(const std::vector<double>& v) { pt_ = v; }

 private:
  void appendPt(const Vec12d& f) { pt_.insert(pt_.end(), f.begin(), f.end()); }

  std::vector<Id> nodes_{0, 1};
  Vec3 up_{0, 0, -1};
  SectionProperties sec_;
  Material mat_;
  BeamElement3D be_;
  bool selfWeightOn_{false};
  std::vector<BeamElement3D::DistributedLoad> segs_;   // 分布荷载（唯一真源）
  std::vector<double> pt_;      // 集中力/力矩的等效节点值，12 个
};

// -----------------------------------------------------------------------------
//  壳单元（包装 ShellElement4）
// -----------------------------------------------------------------------------
class ShellElement : public Element {
 public:
  ShellElement(const std::vector<Id>& ns, const ShellProperties& props)
      : nodes_(ns), props_(props) {
    if (nodes_.size() != 4) nodes_ = {0, 1, 2, 3};
  }

  bool build(const std::vector<Vec3>& pts) {
    if (pts.size() != 4) { err_ = "壳单元必须有 4 个节点"; return false; }
    Vec3 p[4] = {pts[0], pts[1], pts[2], pts[3]};
    sh_ = ShellElement4(p, props_);
    if (!sh_.isValid()) { err_ = sh_.error(); return false; }
    return true;
  }
const std::string& error() const { return err_.empty() ? sh_.error() : err_; }
  const ShellElement4& core() const { return sh_; }
  const ShellProperties& properties() const { return props_; }

  // T3 属性编辑：改壳属性（厚度 / E / ν / 密度）。只改成员；
  // 调用方（Model::setShellProps）负责随后 rebuild 几何。
  void setProperties(const ShellProperties& props) { props_ = props; }

  std::unique_ptr<Element> clone() const override {
    return std::make_unique<ShellElement>(*this);
  }

  int nodeCount() const override { return 4; }
  int dofPerNode() const override { return 6; }
  const std::vector<Id>& nodes() const override { return nodes_; }
  std::vector<Id>& nodesMutable() override { return nodes_; }
  ElementType type() const override { return ElementType::Shell4; }
  std::string typeName() const override { return "壳单元 4 节点"; }
  double length() const override { return 0.0; }
double mass() const override { return sh_.mass(); }
  Vec3 normal() const override { return sh_.normal(); }

  // ---- 剪力墙标志 ----
  //  GridMesh 生成轴线上的竖直墙壳时置真。用途：
  //    · selfweight 命令只把自重开给墙（水平板不吃，避免向既有
  //      模型注入板重造成回归）；
  //    · 楼层剪力统计把竖直墙下端的水平内力计入（见 PostProcessor）。
  void setWall(bool w) { isWall_ = w; }
  bool isWall() const { return isWall_; }

  // ---- 自重（T4）----
  //  墙自重由壳的 mass() 提供（= γ·t·A，γ 由 density 换算），
  //  以等效节点荷载进入模型，方向永远是全局 -Z —— 不能走
  //  setTransversePressure：那沿壳法向，竖直墙法向是水平的，方向就错了。
  void setSelfWeight(bool on) { selfWeightOn_ = on; }
  bool selfWeight() const { return selfWeightOn_; }

  bool stiffnessLocal(LocalMatrix& K) const override {
    if (!sh_.isValid()) return false;
    const Mat24 ke = sh_.stiffness();
    K.resize(576);
    for (int r = 0; r < 24; ++r)
      for (int c = 0; c < 24; ++c) K[static_cast<size_t>(r) * 24 + c] = ke(r, c);
    return true;
  }

void equivalentLoads(std::vector<double>& f) const override {
    f.assign(24, 0.0);
    if (!sh_.isValid()) return;
    if (pz_ != 0.0) {
      const Vec24d p = sh_.transversePressure(pz_);
      for (int i = 0; i < 24; ++i) f[i] += p[i];
    }
    if (px_ != 0.0 || py_ != 0.0) {
      const Vec24d p = sh_.membranePressure(px_, py_);
      for (int i = 0; i < 24; ++i) f[i] += p[i];
    }
    if (selfWeightOn_) {
      // 自重（kN）均分到 4 个角节点，方向为全局 -Z（重力向下）。
      //  【为什么均分】壳的自重是均匀体积力，四节点等效节点荷载
      //  各取 W/4（对任意平面四边形，常体力在双线性插值下的
      //  等效节点力正是均分）—— Σf = W，剪重比分母对账的前提。
      const double w4 = sh_.mass() / 4.0;
      for (int n = 0; n < 4; ++n) f[static_cast<size_t>(n * 6 + 2)] -= w4;
    }
  }

  ShellElement4::MembraneForces shellMembraneForces(const std::vector<double>& uEl) const override {
    Vec24d ue{};
    for (int i = 0; i < 24 && i < static_cast<int>(uEl.size()); ++i) ue[i] = uEl[i];
    return sh_.membraneForces(ue);
  }

void setTransversePressure(double p) { pz_ = p; }
  void setMembranePressure(double px, double py) { px_ = px; py_ = py; }
  // T6 全量序列化：单元荷载读回（与 set 互逆）
  double transversePressure() const { return pz_; }
  double membranePressureX() const { return px_; }
  double membranePressureY() const { return py_; }
  double area() const { return sh_.area(); }
  bool plateReady() const;

private:
  std::vector<Id> nodes_;
  ShellProperties props_;
  ShellElement4 sh_;
  std::string err_;
  double pz_{0.0};
  double px_{0.0}, py_{0.0};
  bool selfWeightOn_{false};
  bool isWall_{false};
};

// -----------------------------------------------------------------------------
//  弹簧单元（包装 SpringElement3D）
//
//  6 个方向（局部 ux,uy,uz,rx,ry,rz）各自独立的两节点弹簧：
//    k=0 → 该方向释放（配合 DofNumbering 零刚度自动锁定）
//    k→∞ → 近似刚接
// -----------------------------------------------------------------------------
class SpringElement : public Element {
 public:
  SpringElement(Id nI, Id nJ, const std::array<double, 6>& k,
                const Vec3& up = Vec3{0, 0, -1})
      : nodes_{nI, nJ}, k_(k), up_(up) {}

  bool build(const Vec3& pi, const Vec3& pj) {
    sp_.init(pi, pj, k_, up_);
    return sp_.isValid();
  }
  const std::string& error() const { return sp_.error(); }
  const SpringElement3D& core() const { return sp_; }
  const std::array<double, 6>& stiffness() const { return k_; }
  Vec3 upHint() const { return up_; }
  void setUpHint(const Vec3& v) { up_ = v; }

  std::unique_ptr<Element> clone() const override {
    return std::make_unique<SpringElement>(*this);
  }

  // ---- Element ----
  int nodeCount() const override { return 2; }
  int dofPerNode() const override { return 6; }
  const std::vector<Id>& nodes() const override { return nodes_; }
  std::vector<Id>& nodesMutable() override { return nodes_; }
  ElementType type() const override { return ElementType::Spring3D; }
  std::string typeName() const override { return "弹簧单元 3D"; }
  double length() const override { return sp_.length(); }
  double mass() const override { return 0.0; }

  // k_d ≈ 0 → 该方向释放
  bool isReleased(int localNode, int comp) const override {
    (void)localNode;
    return sp_.isReleasedComp(comp);
  }

  bool stiffnessLocal(LocalMatrix& K) const override {
    if (!sp_.isValid()) return false;
    const Mat12 ke = sp_.stiffnessGlobal();
    K.resize(144);
    for (int r = 0; r < 12; ++r)
      for (int c = 0; c < 12; ++c) K[static_cast<size_t>(r) * 12 + c] = ke(r, c);
    return true;
  }

  void equivalentLoads(std::vector<double>& f) const override {
    f.assign(12, 0.0);   // 弹簧不承受分布荷载
  }

  using Element::endForces;
  BeamElement3D::EndForces endForces(const std::vector<double>& uEl,
                                     const std::vector<double>& fEq) const override {
    (void)fEq;
    Vec12d ue{};
    for (int i = 0; i < 12 && i < static_cast<int>(uEl.size()); ++i) ue[i] = uEl[i];
    const std::array<double, 6> f = sp_.endForcesLocal(ue);
    BeamElement3D::EndForces ef;
    ef.N = f[0];   ef.Vy = f[1];   ef.Vz = f[2];
    ef.T = f[3];   ef.My = f[4];   ef.Mz = f[5];
    return ef;
  }

 private:
  std::vector<Id> nodes_;
  std::array<double, 6> k_{};
  Vec3 up_{0, 0, -1};
  SpringElement3D sp_;
};

// -----------------------------------------------------------------------------
//  自由度编号器
//
//  【这是模型层唯一真正困难的部分，直接决定求解性能】
//
//  LDL^T 的消元顺序必须满足「对称消元顺序」：若 A(i,j)≠0 且 i≠j，
//  则 new(i) > new(j)，即置换后矩阵仍然是上三角。
//  数学层已经证明 RCM（带宽压缩排序）不满足该性质 —— 施加上去会让
//  置换后第 0 行对应原序最后一行，那一行的非零项全在原矩阵下三角，
//  而稀疏上三角存储无法按行访问下三角 → 整行漏项 → 结果错误且不报错。
//
//  ⇒ 唯一的办法是【编号顺序本身】就是对的消元顺序。
//
//  【盈建科 / PKPM 的工程做法】本实现照此实现：
//    ① 楼层自上而下：顶层先编号
//    ② 同层内：先平动后转动（让同节点的自由度聚在一起，利于带宽）
//    ③ 同一位置组内按 (y, x) 排序，保证结果确定可复现
//
//  为什么"自上而下"有效：上层节点的消元只会向下影响同层与下层，
//  填充被限制在少数层内；而"自下而上"会让每层消元都影响到所有上层，
//  填充迅速吃掉整个上部结构。
//
//  【不做什么】—— 不做自动重编号。
//  自动重编号（AMD / 嵌套剖分）需要实现图消去算法，收益有限；
//  而"建模时按楼层编号"是用户就能做对的事，代价为零。
//  详见 include/yjk/math/Ordering.h 的说明。
// -----------------------------------------------------------------------------
class DofNumbering {
 public:
  enum class Mode {
    InputOrder,        // 按输入顺序（原样透传，调试用）
    StoryDescending,   // 楼层自上而下  ← 默认
  };

  explicit DofNumbering(Mode m = Mode::StoryDescending) : mode_(m) {}

  // 分配自由度编号。返回自由度数。
  // 被约束的自由度保持 dof = -1；从属自由度标成 -2（见 kDofSlave）；
  // 绕法向的转角按需自动约束。
  //
  // 【从属自由度不占方程号】这是 MPC 的关键：它没有自己的方程，
  // 而是通过系数矩阵把自己的贡献折到主自由度上（见 StaticAnalysis::assemble）。
  // 若给它分配了方程号，那行会因为没有刚度而全零 → 矩阵奇异。
  Id assign(std::vector<Node>& nodes, const std::vector<std::unique_ptr<Element>>& elems,
            const std::vector<DofLink>& links = {});

  Id freeDofCount() const { return nfree_; }
  Id constrainedCount() const { return nconstrained_; }
  Id autoConstrainedCount() const { return nauto_; }
  Id slaveCount() const { return nslave_; }
  Mode mode() const { return mode_; }
  const std::vector<Id>& order() const { return order_; }

  // 【诊断】检测到"约束不足"（刚体自由度未消除）时的建议
  const std::vector<std::string>& warnings() const { return warn_; }

 private:
  Mode mode_;
  Id nfree_{0};
  Id nconstrained_{0};
  Id nauto_{0};
  Id nslave_{0};
  std::vector<Id> order_;
  std::vector<std::string> warn_;
};

// -----------------------------------------------------------------------------
//  模型
// -----------------------------------------------------------------------------
class Model {
 public:
  Model() = default;
  explicit Model(DofNumbering::Mode mode) : numbering_(mode) {}

  // ---- 节点 ----
  Id addNode(const Vec3& r, int story = 0, const std::string& label = "");
  Node& node(Id id) { return nodes_[static_cast<size_t>(id)]; }
  const Node& node(Id id) const { return nodes_[static_cast<size_t>(id)]; }
  std::vector<Node>& nodes() { return nodes_; }
  const std::vector<Node>& nodes() const { return nodes_; }
  Id nodeCount() const { return static_cast<Id>(nodes_.size()); }
  Vec3 coord(Id id) const { return nodes_[static_cast<size_t>(id)].r; }

  // ---- 单元 ----
  // 新建梁单元并按给定节点坐标初始化几何
  BeamElement* addBeam(Id nI, Id nJ, const SectionProperties& sec,
                       const Material& mat, const Vec3& up = Vec3{0, 0, -1});
  ShellElement* addShell(const std::vector<Id>& ns, const ShellProperties& props);
  // 新建弹簧单元：k[0..5] = 局部 ux,uy,uz,rx,ry,rz 六方向刚度
  SpringElement* addSpring(Id nI, Id nJ, const std::array<double, 6>& k,
                           const Vec3& up = Vec3{0, 0, -1});

std::vector<std::unique_ptr<Element>>& elements() { return elems_; }
  const std::vector<std::unique_ptr<Element>>& elements() const { return elems_; }
  Id elementCount() const { return static_cast<Id>(elems_.size()); }

  // 按 id 索引查询单元（id = 该单元在 elements() 中的序号）。
  // 等价于 elements()[idx]，但交互层（拾取/属性面板）需要
  // 一个"按 id 访问"的稳定入口 —— 元素创建顺序 = id 分配顺序。
  Element* element(Id idx) { return elems_[static_cast<size_t>(idx)].get(); }
  const Element* element(Id idx) const {
    return elems_[static_cast<size_t>(idx)].get();
  }

  // ---- P3 交互建模：可变接口（脚本层与交互层共用）----

  // 柱 = 竖向梁。与 addBeam 的唯一区别是默认 up（局部 y 参考方向）取
  // 水平方向 {0, -1, 0}，避免与竖直柱轴共线导致 up 回退（BeamElement
  // 的 up 与梁轴平行时会自动切备用方向，柱轴竖直时 {0,0,-1} 恰好共线）。
  BeamElement* addColumn(Id nI, Id nJ, const SectionProperties& sec,
                         const Material& mat,
                         const Vec3& up = Vec3{0, -1, 0});

  // 墙 = 壳单元 + isWall 标志（自重进楼层剪力、参与墙底反力对账）。
  ShellElement* addWall(const std::vector<Id>& ns, const ShellProperties& props);

  // 节点荷载（力 + 力矩）一次性施加；M 缺省为零。
  void addNodeLoad(Id n, const Vec3& F, const Vec3& M = Vec3{0, 0, 0});

  // ---- P3.5 属性编辑（T3）：GUI 改参写回 ----
  //
  // 【写回语义】坐标/截面/材料都是单元的"几何输入"——改完必须重建单元
  // （局部轴系、刚度、质量都固化在 build 时的几何里）。三个 setter 都
  // 遵循"先备份 → 改 → 重建 → 失败回滚"的事务语义：
  //   成功：返回 true，模型已更新（含所有引用该节点的单元重 build）；
  //   失败：恢复原值并重建回旧几何，返回 false，调用方应恢复界面显示。

  // 改节点坐标：更新坐标并重建所有引用该节点的单元几何。
  // 任一元 build 失败（如两节点重合 → 零长度）则回滚整体坐标。
  bool setNodeCoord(Id n, const Vec3& r);

  // 改梁截面 + 材料（合并成一次写回，避免改两遍触发两次几何重建）。
  bool setBeamProps(Id idx, const SectionProperties& sec, const Material& mat);

  // 改壳属性（墙厚 = thickness）。厚度变大 → 刚度与自重随之更新。
  bool setShellProps(Id idx, const ShellProperties& props);

  // 覆盖式设置节点力/力矩（区别于 addNodeLoad 的累加语义：
  // 属性面板编辑"显示的就是当前值"，再点一次保存仍是同一个值）。
  void setNodeForce(Id n, const Vec3& F);
  void setNodeMoment(Id n, const Vec3& M);

  // 深拷贝整个模型（GUI 编辑后"内存重算"的后台线程快照）。
  Model clone() const;

  // 删除第 idx 个单元。后续单元的序号整体前移 —— 交互层删除后
  // 必须刷新选中态（T5 命令模式用"重建选中态"恢复，不依赖旧 id）。
  // 返回是否删除成功（idx 越界返回 false）。
  bool removeElement(Id idx);

  // 删除节点：仅当该节点【未被任何单元引用】且【不是刚性楼板主节点】
  // 时才有效 —— 内部走 retired 标记（见 Node::retired 注释：
  // 结果向量按节点序×6 索引，不能物理删除）。被引用时返回 false，
  // 调用方应先删除相关单元。
  bool removeNode(Id n);

  // 清空整个模型：节点、单元、刚性楼板、自由度编号全部重置。
  void clear();

// ---- 荷载 ----
  void addNodeForce(Id n, const Vec3& F) { nodes_[n].force = nodes_[n].force + F; }
  void addNodeMoment(Id n, const Vec3& M) { nodes_[n].moment = nodes_[n].moment + M; }
  // 节点附加重量（kN，重力荷载代表值折算后注入；clearLoads 不清除）
  void addNodeWeight(Id n, double w) { nodes_[n].weight += w; }
  // 质量源检测：任一节点有附加重量 → ModalAnalysis 用"节点重量"代表团
  bool hasNodalWeight() const {
    for (const Node& nd : nodes_) if (nd.weight > 0.0) return true;
    return false;
  }
  void clearLoads();

  // ---- 约束 ----
  void fixNode(Id n, bool ux, bool uy, bool uz, bool rx, bool ry, bool rz) {
    Node& nd = nodes_[n];
    nd.fixed[0] = ux; nd.fixed[1] = uy; nd.fixed[2] = uz;
    nd.fixed[3] = rx; nd.fixed[4] = ry; nd.fixed[5] = rz;
  }
  void fixAll(Id n) { fixNode(n, true, true, true, true, true, true); }
  // 常见支座
  void restrainPinned(Id n) { fixNode(n, true, true, true, true, false, true); }
  void restrainRoller(Id n)   { fixNode(n, true, true, true, false, false, false); }

  // ---- 刚性楼板（多点约束）----
  //
  //  【工程意义】真实工程的建模习惯是"只建梁柱 + 加刚性楼板假定"，
  //  而不是真去建楼板单元。刚性假定说的是：楼面在自身平面内刚度无穷大，
  //  于是整层楼在平面内只可能做刚体运动 —— 3 个自由度就够了。
  //  这是 YJK/PKPM 最核心的用户体验之一。
  //
  //  【刚体运动只约束面内 3 个自由度】uz / rx / ry 不动 —— 楼板的面外
  //  刚度本来就不是"无穷大"，竖向弯曲仍由梁承担（这正是"不建板"的代价）。
  struct Diaphragm {
    int story{0};
    Id masterNode{kDofFixed};        // 自动创建的主节点
    Vec3 center{0, 0, 0};            // 约束参考点（形心）
    std::vector<Id> slaves;          // 被绑定的节点
    bool coupleRz{true};             // 是否同时约束绕竖轴的转角
  };

  // 把一组节点绑成刚性楼板，返回新建的主节点号。
  //
  //  【主节点位置只是参考点】约束集与参考点无关（换一个参考点只是换了
  //  几个平移量的数值），所以取形心是安全的、也是唯一与用户直觉一致的选择。
  //  主节点只有 (ux, uy, rz) 三个自由度参与求解，其余三个直接约束掉。
  Id addRigidDiaphragm(const std::vector<Id>& nodes, int story, bool coupleRz = true);

  // T6 全量序列化：把一组节点绑到【已存在的主节点】上，重建刚性楼板。
  //
  //  与 addRigidDiaphragm 的区别：它不新建主节点 —— 主节点在序列化文本里
  //  就是普通 node（带坐标与固定约束），读回后由本接口把该节点升级为主节点
  //  并重建 DofLink。参考点取主节点坐标（= 原形心），links 系数与源一致。
  Id attachRigidDiaphragm(Id master, const std::vector<Id>& slaves,
                          int story, bool coupleRz = true);

  const std::vector<DofLink>& dofLinks() const { return links_; }
  const std::vector<Diaphragm>& diaphragms() const { return diaphragms_; }
  void clearDiaphragms() { links_.clear(); diaphragms_.clear(); }

  // ---- 编号 ----
  Id assignDofs(DofNumbering::Mode m) {
    numbering_ = DofNumbering(m);
    return numbering_.assign(nodes_, elems_, links_);
  }
  Id assignDofs() { return numbering_.assign(nodes_, elems_, links_); }
  Id freeDofCount() const { return numbering_.freeDofCount(); }
  const std::vector<Id>& dofOrder() const { return numbering_.order(); }
  const std::vector<std::string>& dofWarnings() const { return numbering_.warnings(); }

  // ---- 校核 ----
  // 检查几何退化、重复节点、悬空单元、约束不足等。
  // 这些问题在建模阶段最容易犯，而求解器只会报"奇异"却不会说清原因。
  struct Diagnostic {
    std::vector<std::string> errors;     // 必须修
    std::vector<std::string> warnings;   // 建议修
    bool ok() const { return errors.empty(); }
    std::string summary() const;
  };
  Diagnostic check() const;

  double totalMass() const {
    double m = 0.0;
    for (const auto& e : elems_) m += e->mass();
    return m;
  }

 private:
  std::vector<Node> nodes_;
  std::vector<std::unique_ptr<Element>> elems_;
  std::vector<DofLink> links_;
  std::vector<Diaphragm> diaphragms_;
  DofNumbering numbering_{DofNumbering::Mode::StoryDescending};
};

}  // namespace yjk
