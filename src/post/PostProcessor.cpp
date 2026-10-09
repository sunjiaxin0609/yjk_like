// =============================================================================
//  src/post/PostProcessor.cpp
// =============================================================================
#include "yjk/post/PostProcessor.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>

namespace yjk {
namespace post {

namespace {

constexpr double kPi = 3.14159265358979323846;

// 壳单元质量接口给的是 kg（density 单位 kg/m³），转成重力荷载 kN
double shellWeight(double area, double thickness, double density) {
  return area * thickness * density * 9.81 / 1000.0;
}

// σ = N/A ± My·hz/Iy ± Mz·hy/Iz
//
// 【符号约定】轴力拉为正；弯曲项取绝对值叠加，得到截面两条边缘纤维的极值。
// 双向弯矩同时作用时，"最大"一定出现在两个角点之一，
// 所以四个角点里取最大/最小 —— 这里直接叠加绝对值等价于取角点极值。
struct FiberStress {
  double sigMax{0.0}, sigMin{0.0};
};

FiberStress fiberStress(double N, double My, double Mz, const SectionProperties& s) {
  FiberStress r;
  if (!(s.A > 0.0)) return r;
  const double axial = N / s.A;
  const double by = (s.Iy > 0.0) ? std::abs(My) * s.halfWidthZ() / s.Iy : 0.0;
  const double bz = (s.Iz > 0.0) ? std::abs(Mz) * s.halfDepthY() / s.Iz : 0.0;
  r.sigMax = axial + by + bz;
  r.sigMin = axial - by - bz;
  return r;
}

std::string fixed(double v, int w, int p) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%*.*f", w, p, v);
  return std::string(buf);
}

std::string sci(double v, int p = 3) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.*e", p, v);
  return std::string(buf);
}

// 层间位移角写成 1/N 的形式（工程师读这个）
std::string driftAsOneOver(double d) {
  if (!(d > 0.0)) return "  —  ";
  const double n = 1.0 / d;
  char buf[64];
  std::snprintf(buf, sizeof(buf), "1/%5.0f", n);
  return std::string(buf);
}

}  // namespace

const char* memberKindName(MemberKind k) {
  switch (k) {
    case MemberKind::Column: return "柱";
    case MemberKind::Beam:   return "梁";
    case MemberKind::Brace:  return "斜撑";
    default:                 return "杆";
  }
}

// -----------------------------------------------------------------------------
//  梁：沿长分段内力 + 应力
// -----------------------------------------------------------------------------
void PostProcessor::computeBeams(int nStation, double columnCos) {
  beams_.clear();
  if (nStation < 2) nStation = 2;

  const Id ne = model_.elementCount();
  for (Id e = 0; e < ne; ++e) {
    const Element& base = *model_.elements()[static_cast<size_t>(e)];
    if (base.type() != ElementType::Beam3D) continue;
    const BeamElement& el = static_cast<const BeamElement&>(base);
    const BeamElement3D& be = el.core();
    if (be.length() <= 0.0) continue;

    // 单元局部位移与等效节点荷载（atStation 需要固端修正）
    std::vector<double> uEl(12, 0.0), fEq;
    const std::vector<Id>& nds = el.nodes();
    for (int a = 0; a < 12; ++a) {
      const Id nd = nds[static_cast<size_t>(a / 6)];
      const int comp = a % 6;
      const size_t gi = static_cast<size_t>(nd * 6 + comp);
      if (gi < res_.u.size()) uEl[static_cast<size_t>(a)] = res_.u[gi];
    }
    el.equivalentLoads(fEq);
    Vec12d uu{}, ff{};
    for (int i = 0; i < 12; ++i) {
      uu[static_cast<size_t>(i)] = uEl[static_cast<size_t>(i)];
      ff[static_cast<size_t>(i)] = (i < static_cast<int>(fEq.size())) ? fEq[static_cast<size_t>(i)] : 0.0;
    }

    BeamResult br;
    br.id = e;
    br.ni = nds[0];
    br.nj = nds[1];
    br.L = be.length();
    br.dir = {be.localX().x, be.localX().y, be.localX().z};
    br.ends = res_.beamForces[static_cast<size_t>(e)].ef;
    br.weight = el.mass();

    const double cz = std::abs(br.dir[2]);
    if (cz > columnCos)                 br.kind = MemberKind::Column;
    else if (cz < 0.2)                  br.kind = MemberKind::Beam;
    else                                br.kind = MemberKind::Brace;

    const SectionProperties& sec = el.section();
    const double fy = el.material().fy;
    const double fc = el.material().fc;
    const double fDesign = fy > 0.0 ? fy : (fc > 0.0 ? fc : 0.0);

    br.stations.reserve(static_cast<size_t>(nStation));
    double vmMax = 0.0, sigAbsMax = 0.0;
    double mMax = 0.0, xiMMax = 0.0, vMax = 0.0, nMax = -1e300, nMin = 1e300;
    for (int k = 0; k < nStation; ++k) {
      const double xi = static_cast<double>(k) / static_cast<double>(nStation - 1);
      // 【走 wrapper 的 atStation】它会把分布荷载（梯形/三角形）一并传下去，
      // 使内力图按真实荷载积分 —— 直接调 core 的版本只知道"端值"，会把
      // 三角形荷载的二次剪力图画成直线。
      const auto ef = el.atStation(uu, ff, xi);

      BeamStation st;
      st.xi = xi;
      st.N = ef.N; st.Vy = ef.Vy; st.Vz = ef.Vz;
      st.T = ef.T; st.My = ef.My; st.Mz = ef.Mz;

      // 剪力转全局：V_global = Vy·ey + Vz·ez
      //
      //  【为什么必须留全局分量】楼层剪力要按【矢量】累加：
      //  自重会让框架"内缩"，产生一组自平衡的水平反力（Σ 为零）。
      //  若按 Σ|V| 累加，这组自平衡量会被当成剪力计入，
      //  实测把 40 kN 的基底剪力算成 40.58 kN（虚高 1.5%）。
      //  按矢量累加则自动相消 —— 这才是"楼层剪力"的定义。
      {
        const Vec3& ey = be.localY();
        const Vec3& ez = be.localZ();
        st.vx = ef.Vy * ey.x + ef.Vz * ez.x;
        st.vy = ef.Vy * ey.y + ef.Vz * ez.y;
        st.vz = ef.Vy * ey.z + ef.Vz * ez.z;
      }

      const FiberStress fs = fiberStress(ef.N, ef.My, ef.Mz, sec);
      st.sigMax = fs.sigMax;
      st.sigMin = fs.sigMin;
      st.tau = (sec.A > 0.0) ? 1.5 * std::hypot(ef.Vy, ef.Vz) / sec.A : 0.0;
      const double sAbs = std::max(std::abs(st.sigMax), std::abs(st.sigMin));
      st.vonMises = std::sqrt(sAbs * sAbs + 3.0 * st.tau * st.tau);

      br.stations.push_back(st);

      vmMax = std::max(vmMax, st.vonMises);
      sigAbsMax = std::max(sigAbsMax, sAbs);
      br.sigMax = std::max(br.sigMax, st.sigMax);
      br.sigMin = std::min(br.sigMin, st.sigMin);

      const double mm = std::hypot(ef.My, ef.Mz);
      if (mm > mMax) { mMax = mm; xiMMax = xi; }
      vMax = std::max(vMax, std::hypot(ef.Vy, ef.Vz));
      nMax = std::max(nMax, ef.N);
      nMin = std::min(nMin, ef.N);
    }
    br.sigAbsMax = sigAbsMax;
    br.vonMisesMax = vmMax;
    br.util = (fDesign > 0.0) ? sigAbsMax / fDesign : 0.0;
    br.Mmax = mMax;
    br.xiMmax = xiMMax;
    br.Vmax = vMax;
    br.Nmax = nMax;
    br.Nmin = nMin;
    beams_.push_back(std::move(br));
  }
}

void PostProcessor::computeShells() {
  shells_.clear();
  const Id ne = model_.elementCount();
  for (Id e = 0; e < ne; ++e) {
    const Element& base = *model_.elements()[static_cast<size_t>(e)];
    if (base.type() != ElementType::Shell4) continue;
    const ShellElement& sh = static_cast<const ShellElement&>(base);

    ShellResult sr;
    sr.id = e;
    for (int i = 0; i < 4; ++i) sr.nodes[static_cast<size_t>(i)] = sh.nodes()[static_cast<size_t>(i)];
    sr.area = sh.area();
    sr.thickness = sh.properties().thickness;
    sr.mf = res_.shellForces[static_cast<size_t>(e)].mf;
    sr.weight = shellWeight(sr.area, sr.thickness, sh.properties().density);
    shells_.push_back(sr);
  }
}

// -----------------------------------------------------------------------------
//  节点场
// -----------------------------------------------------------------------------
void PostProcessor::computeNodeFields() {
  fields_.clear();
  const Id nn = model_.nodeCount();
  const size_t n = static_cast<size_t>(nn);
  if (nn <= 0) return;

  // 【必须先 reserve】addField 返回的是 vector 元素的引用，
  // 后面再 push_back 会触发重新分配，前面的引用全部失效 ——
  // 表现为"写到已释放内存"，是最难排查的一类崩溃。
  // reserve 之后容量固定，引用才安全。
  fields_.reserve(16);

  auto addField = [&](const std::string& name, const std::string& unit) -> NodeField& {
    fields_.push_back(NodeField{});
    NodeField& f = fields_.back();
    f.name = name;
    f.unit = unit;
    f.v.assign(n, 0.0);
    return f;
  };

  // ---- 位移 ----
  NodeField& fMag = addField("位移幅值", "m");
  NodeField& fux = addField("ux", "m");
  NodeField& fuy = addField("uy", "m");
  NodeField& fuz = addField("uz", "m");
  NodeField& fRot = addField("转角幅值", "rad");
  for (Id i = 0; i < nn; ++i) {
    const Vec3 d = model_.node(i).displacement(res_.u);
    const Vec3 r = model_.node(i).rotation(res_.u);
    const double m = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    fMag.v[static_cast<size_t>(i)] = m;
    fux.v[static_cast<size_t>(i)] = d.x;
    fuy.v[static_cast<size_t>(i)] = d.y;
    fuz.v[static_cast<size_t>(i)] = d.z;
    fRot.v[static_cast<size_t>(i)] = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
  }

  // ---- 应力：绕节点平均 ----
  //
  // 壳给"单元中心常值"，梁给"端部值"。两者量纲一致（kPa）但精度不同，
  // 这里做等权平均 —— 简单平均在网格疏密不均时会偏，
  // 但对工程读数足够，且不会因为权重口径不一致引入新错误。
  std::vector<double> accSx(n, 0.0), accSy(n, 0.0), accTxy(n, 0.0);
  std::vector<double> accVm(n, 0.0), cntVm(n, 0.0);
  std::vector<double> cntS(n, 0.0);

  for (const ShellResult& sr : shells_) {
    for (int i = 0; i < 4; ++i) {
      const size_t nd = static_cast<size_t>(sr.nodes[static_cast<size_t>(i)]);
      if (nd >= n) continue;
      accSx[nd] += sr.mf.sx;
      accSy[nd] += sr.mf.sy;
      accTxy[nd] += sr.mf.txy;
      accVm[nd] += sr.mf.vonMises;
      cntS[nd] += 1.0;
      cntVm[nd] += 1.0;
    }
  }
  for (const BeamResult& br : beams_) {
    if (br.stations.empty()) continue;
    const double vmI = br.stations.front().vonMises;
    const double vmJ = br.stations.back().vonMises;
    const size_t a = static_cast<size_t>(br.ni), b = static_cast<size_t>(br.nj);
    if (a < n) { accVm[a] += vmI; cntVm[a] += 1.0; }
    if (b < n) { accVm[b] += vmJ; cntVm[b] += 1.0; }
  }

  NodeField& fSx = addField("sx", "kPa");
  NodeField& fSy = addField("sy", "kPa");
  NodeField& fTxy = addField("txy", "kPa");
  NodeField& fVm = addField("von Mises 应力", "kPa");
  for (size_t i = 0; i < n; ++i) {
    if (cntS[i] > 0.0) {
      fSx.v[i] = accSx[i] / cntS[i];
      fSy.v[i] = accSy[i] / cntS[i];
      fTxy.v[i] = accTxy[i] / cntS[i];
    }
    if (cntVm[i] > 0.0) fVm.v[i] = accVm[i] / cntVm[i];
  }

  // ---- 支座反力 ----
  NodeField& fRz = addField("竖向反力", "kN");
  for (Id i = 0; i < nn; ++i)
    fRz.v[static_cast<size_t>(i)] = res_.reaction[static_cast<size_t>(i) * 6 + 2];

  for (NodeField& f : fields_) {
    double mn = 1e300, mx = -1e300, mnAbs = 1e300, mxAbs = 0.0;
    for (double v : f.v) {
      mn = std::min(mn, v);
      mx = std::max(mx, v);
      const double a = std::abs(v);
      mnAbs = std::min(mnAbs, a);
      mxAbs = std::max(mxAbs, a);
    }
    if (mn > mx) { mn = mx = 0.0; mnAbs = 0.0; }
    f.vmin = mn; f.vmax = mx; f.vminAbs = mnAbs; f.vmaxAbs = mxAbs;
  }
}

// -----------------------------------------------------------------------------
//  楼层指标
// -----------------------------------------------------------------------------
void PostProcessor::computeStories(double columnCos) {
  stories_.clear();

  const Id nn = model_.nodeCount();
  std::set<int> storySet;
  for (Id i = 0; i < nn; ++i) {
    if (model_.node(i).retired) continue;
    storySet.insert(model_.node(i).story);
  }
  if (storySet.empty()) return;

  std::vector<int> ids(storySet.begin(), storySet.end());
  stories_.resize(ids.size());
  std::vector<std::vector<Id>> members(ids.size());
  for (size_t s = 0; s < ids.size(); ++s) stories_[s].story = ids[s];

  for (Id i = 0; i < nn; ++i) {
    const Node& nd = model_.node(i);
    if (nd.retired) continue;
    for (size_t s = 0; s < ids.size(); ++s)
      if (ids[s] == nd.story) { members[s].push_back(i); break; }
  }

  // 标高：取该层节点 z 的平均（同一楼面的节点 z 应当完全一致）
  for (size_t s = 0; s < ids.size(); ++s) {
    double z = 0.0;
    for (Id n : members[s]) z += model_.node(n).r.z;
    stories_[s].z = members[s].empty() ? 0.0 : z / static_cast<double>(members[s].size());
    stories_[s].nodeCount = static_cast<int>(members[s].size());
    stories_[s].h = (s > 0) ? stories_[s].z - stories_[s - 1].z : 0.0;
  }

  // 水平位移
  for (size_t s = 0; s < stories_.size(); ++s) {
    double mx = 0.0, sum = 0.0, mz = 0.0;
    for (Id n : members[s]) {
      const Vec3 d = model_.node(n).displacement(res_.u);
      const double h = std::hypot(d.x, d.y);
      mx = std::max(mx, h);
      sum += h;
      mz = std::max(mz, std::abs(d.z));
    }
    stories_[s].maxUxy = mx;
    stories_[s].avgUxy = members[s].empty() ? 0.0 : sum / static_cast<double>(members[s].size());
    stories_[s].maxUz = mz;
  }

  // ---- 层间位移角 ----
  //
  //  【必须与"正下方的同一个点"比较】
  //  用本层最大位移减下层最大位移是错的 —— 两个最大值通常不在同一榀，
  //  相减得到的是毫无意义的量。
  //  这里按水平坐标配对（同一根柱的上下节点），找不到配对才退回层平均值。
  for (size_t s = 1; s < stories_.size(); ++s) {
    const double h = stories_[s].h;
    if (!(h > 1e-9)) continue;
    double mx = 0.0, sum = 0.0;
    int cnt = 0;
    Id node = -1;
    for (Id n : members[s]) {
      const Vec3& p = model_.node(n).r;
      // 在下层找同 (x, y) 的节点
      Id below = -1;
      double best = 1e-3;
      for (Id m : members[s - 1]) {
        const Vec3& q = model_.node(m).r;
        const double d = std::hypot(p.x - q.x, p.y - q.y);
        if (d < best) { best = d; below = m; }
      }
      if (below < 0) continue;
      const Vec3 du = model_.node(n).displacement(res_.u);
      const Vec3 dl = model_.node(below).displacement(res_.u);
      const double du1 = std::hypot(du.x - dl.x, du.y - dl.y) / h;
      mx = std::max(mx, du1);
      sum += du1;
      ++cnt;
      if (du1 >= mx) node = n;
    }
    if (cnt == 0) {   // 退回层平均
      const double d = std::abs(stories_[s].avgUxy - stories_[s - 1].avgUxy) / h;
      mx = d; sum = d; cnt = 1;
    }
    stories_[s].drift = mx;
    stories_[s].driftAvg = sum / static_cast<double>(cnt);
    stories_[s].driftNode = node;
    stories_[s].driftRatio = (stories_[s].driftAvg > 1e-12)
                                 ? stories_[s].maxUxy / stories_[s].avgUxy : 0.0;
  }

  // ---- 楼层剪力 ----
  //
  //  定义：第 s 层的层剪力 = 该层所有柱（连接 s−1 与 s）的【水平剪力矢量之和的模】。
  //
  //  【必须用矢量和，不能用 Σ|V|】
  //  自重让框架"内缩"会产生一组自平衡的水平力（一半柱子推 +X、一半推 −X，
  //  合力为零）。Σ|V| 会把它们当成剪力计入，实测 40 kN 的基底剪力被算成
  //  40.58 kN。矢量和自动相消。
  //
  //  校核：最底一层的层剪力 = 所有支座水平反力矢量和的模。
  std::vector<double> sx(stories_.size(), 0.0), sy(stories_.size(), 0.0);
  for (const BeamResult& br : beams_) {
    if (br.kind != MemberKind::Column) continue;
    const int sa = model_.node(br.ni).story;
    const int sb = model_.node(br.nj).story;
    const int hi = std::max(sa, sb);
    if (hi == std::min(sa, sb)) continue;                 // 层间柱才算
    for (size_t s = 0; s < stories_.size(); ++s)
      if (stories_[s].story == hi) {
        // 取下端的水平剪力（下标端 = 楼层较低的一端）
        const BeamStation& st = (sa < sb) ? br.stations.front() : br.stations.back();
        sx[s] += st.vx;
        sy[s] += st.vy;
        break;
      }
  }
  for (size_t s = 0; s < stories_.size(); ++s)
    stories_[s].shear = std::hypot(sx[s], sy[s]);
  (void)columnCos;

  // ---- 各层重力荷载（自上而下累加）----
  //
  //  【为什么必须累加【荷载向量】而不是【构件自重】】
  //  重力荷载代表值 = 自重 + 附加恒载 + 活载。
  //  楼板面荷载、节点荷载同样是重力荷载，只统计构件自重会漏掉一大块
  //  （本例实测：自重 3175 kN，加上 4 kPa 楼面荷载后是 4903 kN，差 35%）。
  //  剪重比 V/W 的分母漏掉面荷载 ⇒ 剪重比被【高估】，
  //  而剪重比是抗震验算的【下限】指标 —— 高估即偏不安全。
  //  所以这里直接扫一遍单元等效荷载 + 节点荷载，而不是用 mass()。
  std::vector<double> wOf(stories_.size(), 0.0);
  auto storyIndex = [&](int id) -> size_t {
    for (size_t s = 0; s < ids.size(); ++s) if (ids[s] == id) return s;
    return ids.size();
  };
  storyLoad_ = wOf;                      // 供 computeEnvelope 复用
  for (const auto& ep : model_.elements()) {
    std::vector<double> fe;
    ep->equivalentLoads(fe);
    const int dpn = ep->dofPerNode();
    const auto& nds = ep->nodes();
    for (size_t a = 0; a < nds.size(); ++a) {
      const double fz = fe[a * static_cast<size_t>(dpn) + 2];
      if (!(fz < 0.0)) continue;                        // 只统计向下（含 NaN 保护）
      const size_t s = storyIndex(model_.node(nds[a]).story);
      if (s < storyLoad_.size()) storyLoad_[s] += -fz;
    }
  }
  for (Id i = 0; i < model_.nodeCount(); ++i) {
    const Node& nd = model_.node(i);
    if (nd.retired || !(nd.force.z < 0.0)) continue;
    const size_t s = storyIndex(nd.story);
    if (s < storyLoad_.size()) storyLoad_[s] += -nd.force.z;
  }
  double above = 0.0;
  for (size_t s = stories_.size(); s-- > 0;) {
    above += storyLoad_[s];
    stories_[s].weight = above;
    stories_[s].shearWeightRatio = (above > 0.0) ? stories_[s].shear / above : 0.0;
  }
}

// -----------------------------------------------------------------------------
void PostProcessor::computeEnvelope() {
  env_ = Envelope{};
  const Id nn = model_.nodeCount();
  for (Id i = 0; i < nn; ++i) {
    const Node& nd = model_.node(i);
    if (nd.retired) continue;
    const Vec3 d = nd.displacement(res_.u);
    const Vec3 r = nd.rotation(res_.u);
    const double dm = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    const double h = std::hypot(d.x, d.y);
    const double rm = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    if (dm > env_.maxDisp) { env_.maxDisp = dm; env_.maxDispNode = i; }
    if (h > env_.maxUxy) { env_.maxUxy = h; env_.maxUxyNode = i; }
    if (std::abs(d.z) > std::abs(env_.maxUz)) { env_.maxUz = d.z; env_.maxUzNode = i; }
    if (rm > env_.maxRotation) { env_.maxRotation = rm; env_.maxRotationNode = i; }

    const size_t b = static_cast<size_t>(i * 6);
    for (int k = 0; k < 6; ++k) {
      const double rv = std::abs(res_.reaction[b + static_cast<size_t>(k)]);
      if (rv > env_.maxReaction) { env_.maxReaction = rv; env_.maxReactionNode = i; }
    }
    const double rz = std::abs(res_.reaction[b + 2]);
    if (rz > env_.maxReactionZ) { env_.maxReactionZ = rz; env_.maxReactionZNode = i; }
    env_.totalLoadZ += res_.reaction[b + 2];
  }

  // 总重力荷载 = 各层施加荷载之和（含面荷载、节点荷载，不只是构件自重）
  for (double w : storyLoad_) env_.totalWeight += w;

  // 基底剪力 = 所有支座水平反力【矢量和】的模（与最底层层剪力互为校核）
  //
  //  同样不能按 Σ|R| 累加：自重引起的自平衡水平反力会混进来。
  double bx = 0.0, by = 0.0;
  for (Id i = 0; i < nn; ++i) {
    const size_t b = static_cast<size_t>(i * 6);
    bx += res_.reaction[b];
    by += res_.reaction[b + 1];
  }
  env_.baseShear = std::hypot(bx, by);
  // 剪重比分母优先取【实际施加的竖向荷载】（= 竖向反力合计），
  // 而不是构件自重合计 —— 两者不等时说明有构件自重没导成荷载
  // （例如楼板没设面荷载），用后者会低估剪重比。
  const double W = (std::abs(env_.totalLoadZ) > 1e-9) ? std::abs(env_.totalLoadZ)
                                                      : env_.totalWeight;
  env_.shearWeightRatio = (W > 0.0) ? env_.baseShear / W : 0.0;

  for (const StoryResult& sr : stories_) {
    if (sr.drift > env_.maxDrift) { env_.maxDrift = sr.drift; env_.maxDriftStory = sr.story; }
    if (sr.shearWeightRatio > env_.maxShearRatio) {
      env_.maxShearRatio = sr.shearWeightRatio;
      env_.maxShearRatioStory = sr.story;
    }
  }

  const NodeField* vm = field("von Mises 应力");
  if (vm) {
    for (Id i = 0; i < nn; ++i) {
      const double v = vm->v[static_cast<size_t>(i)];
      if (v > env_.maxVonMises) { env_.maxVonMises = v; env_.maxVonMisesNode = i; }
    }
  }
}

void PostProcessor::compute(int stations, double columnCos) {
  computeBeams(stations, columnCos);
  computeShells();
  computeNodeFields();
  computeStories(columnCos);
  computeEnvelope();
}

// -----------------------------------------------------------------------------
const NodeField* PostProcessor::field(const std::string& name) const {
  for (const NodeField& f : fields_)
    if (f.name == name) return &f;
  return nullptr;
}

std::vector<std::string> PostProcessor::fieldNames() const {
  std::vector<std::string> out;
  for (const NodeField& f : fields_) out.push_back(f.name);
  return out;
}

std::vector<const BeamResult*> PostProcessor::topBeamsByStress(int n) const {
  std::vector<const BeamResult*> v;
  for (const BeamResult& b : beams_) v.push_back(&b);
  std::stable_sort(v.begin(), v.end(), [](const BeamResult* a, const BeamResult* b) {
    return a->sigAbsMax > b->sigAbsMax;
  });
  if (static_cast<int>(v.size()) > n) v.resize(static_cast<size_t>(n));
  return v;
}

std::vector<const BeamResult*> PostProcessor::columns() const {
  std::vector<const BeamResult*> v;
  for (const BeamResult& b : beams_)
    if (b.kind == MemberKind::Column) v.push_back(&b);
  return v;
}

// -----------------------------------------------------------------------------
std::string PostProcessor::summary() const {
  std::ostringstream os;
  os << "自由度 " << res_.ndof << "，单元 " << beams_.size() + shells_.size()
     << "，最大位移 " << sci(env_.maxDisp, 3) << " m";
  if (env_.maxDrift > 0.0)
    os << "，最大层间位移角 " << driftAsOneOver(env_.maxDrift);
  os << "，残差 " << sci(res_.residual, 2);
  return os.str();
}

std::string PostProcessor::report() const {
  std::ostringstream os;
  os << "======================================================\n";
  os << "  yjk_like  静力分析结果\n";
  os << "======================================================\n";
  os << "求解：自由度 " << res_.ndof << "，约束 " << res_.nconstrained
     << "，nnz " << res_.nnz << "，nnzL " << res_.nnzL << "\n";
  os << "      残差 " << sci(res_.residual, 3) << "，耗时 " << fixed(res_.seconds, 0, 3) << " s\n\n";

  os << "-- 位移包络 --\n";
  os << "  最大位移         " << sci(env_.maxDisp, 4) << " m    节点 " << env_.maxDispNode << "\n";
  os << "  最大水平位移     " << sci(env_.maxUxy, 4) << " m    节点 " << env_.maxUxyNode << "\n";
  os << "  最大竖向位移     " << sci(env_.maxUz, 4) << " m    节点 " << env_.maxUzNode << "\n";
  os << "  最大转角         " << sci(env_.maxRotation, 4) << " rad  节点 " << env_.maxRotationNode << "\n";
  os << "  最大 von Mises   " << fixed(env_.maxVonMises / 1000.0, 10, 2) << " MPa  节点 " << env_.maxVonMisesNode << "\n\n";

  if (!stories_.empty()) {
    os << "-- 楼层指标 --\n";
    os << "  层   标高    层高   最大水平位移   层间位移角   位移比   楼层剪力   剪重比\n";
    for (const StoryResult& s : stories_) {
      os << "  " << fixed(static_cast<double>(s.story), 2, 0)
         << "  " << fixed(s.z, 6, 2)
         << "  " << fixed(s.h, 5, 2)
         << "  " << sci(s.maxUxy, 3) << " m"
         << "  " << driftAsOneOver(s.drift)
         << "  " << fixed(s.driftRatio, 6, 2)
         << "  " << fixed(s.shear, 9, 2) << " kN"
         << "  " << fixed(s.shearWeightRatio, 7, 4) << "\n";
    }
    os << "  最大层间位移角 " << driftAsOneOver(env_.maxDrift)
       << "（第 " << env_.maxDriftStory << " 层）\n";
    os << "  基底剪力 " << fixed(env_.baseShear, 0, 2) << " kN，总重力荷载 "
       << fixed(env_.totalWeight, 0, 2) << " kN，剪重比 "
       << fixed(env_.shearWeightRatio, 0, 4) << "\n\n";
  }

  os << "-- 构件应力 TOP 10 --\n";
  os << "  单元  类别  长度   |σ|max(MPa)   vonMises(MPa)   控制弯矩(kN·m)  最大剪力(kN)  轴力(kN)\n";
  const auto top = topBeamsByStress(10);
  for (const BeamResult* b : top) {
    os << "  " << fixed(static_cast<double>(b->id), 4, 0)
       << "  " << memberKindName(b->kind)
       << "  " << fixed(b->L, 5, 2)
       << "  " << fixed(b->sigAbsMax / 1000.0, 11, 2)
       << "  " << fixed(b->vonMisesMax / 1000.0, 13, 2)
       << "  " << fixed(b->Mmax, 15, 2)
       << "  " << fixed(b->Vmax, 12, 2)
       << "  " << fixed(b->Nmin, 10, 2) << "\n";
  }
  os << "\n";

  // ---------------------------------------------------------------------
  //  【平衡校核 —— 三条独立的账必须对上】
  //
  //    施加的竖向荷载 ΣW   由【荷载装配】统计（自重 + 面荷载 + 节点荷载）
  //    竖向反力合计  ΣRz  由【K·u − f】恢复（完全不同的代码路径）
  //
  //  两条路径独立，能对上就说明：荷载没漏、反力公式没错。
  //  对不上时 90% 是荷载漏项（面荷载没导进去、自重开关没开），
  //  剩下 10% 是反力恢复又重复减了一次荷载。
  // ---------------------------------------------------------------------
  os << "-- 平衡校核 --\n";
  os << "  施加竖向荷载 " << fixed(env_.totalWeight, 0, 4) << " kN"
     << "，竖向反力合计 " << fixed(env_.totalLoadZ, 0, 4) << " kN";
  const double close = (env_.totalWeight > 1e-9)
      ? std::abs(env_.totalLoadZ - env_.totalWeight) / env_.totalWeight : 0.0;
  os << "，闭合差 " << sci(close, 3)
     << (close < 1e-8 ? "  [闭合]" : "  [!! 不闭合：荷载漏项或反力重复减荷载]") << "\n";
  os << "  基底剪力 " << fixed(env_.baseShear, 0, 4) << " kN\n";
  os << "  残差 " << sci(res_.residual, 3) << (res_.residual < 1e-8 ? "  [满足]" : "  [!! 偏大]") << "\n";

  return os.str();
}

}  // namespace post
}  // namespace yjk
