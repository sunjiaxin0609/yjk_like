// =============================================================================
//  src/analysis/ModalAnalysis.cpp  ——  模态分析实现
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>

#include "yjk/analysis/ModalAnalysis.h"

namespace yjk {

namespace {
constexpr double kGravity = 9.81;   // m/s²，把单元自重(kN)换算成质量(t)

// 质量折减：对角集中质量 → Tᵀ M T（与刚度折减共用 expand 表）。
// 【与刚度的区别】刚度项 k(a,b) 是两个自由度各展开一次的内积（双和）；
//  质量项 m 只挂在一个自由度 gd 上，展开后 M_r(p,q) = m·c_p·c_q 是
//  【无序对】的唯一值 —— 遍历 p≤q 即可，绝不能像刚度那样遍历全对
//  （那样 (p,q) 与 (q,p) 上三角化后落进同一槽位累加 ⇒ 质量翻倍）。
void foldMass(const std::vector<double>& agg,
              const std::vector<StaticAnalysis::Assembly::Expand>& expand,
              TripletAssembler& ta) {
  for (size_t gd = 0; gd < agg.size(); ++gd) {
    const double m = agg[gd];
    if (m == 0.0) continue;
    const StaticAnalysis::Assembly::Expand& e = expand[gd];
    if (e.n == 0) continue;                       // 固定自由度：无动能，质量丢弃
    if (e.n == 1) {
      // 自由自由度（d[0] = 自身方程号，c = 1）：对角项
      ta.add(e.d[0], e.d[0], m * e.c[0] * e.c[0]);
      continue;
    }
    // 从属自由度：TᵀMT 折到主自由度（p ≤ q，见函数头注释）
    for (int p = 0; p < e.n; ++p)
      for (int q = p; q < e.n; ++q)
        ta.add(e.d[p], e.d[q], m * e.c[p] * e.c[q]);
  }
}

// 质量聚合向量（全局自由度号索引，节点序 × 6）。
//
// 【质量源的两条路径 —— GB 50011 5.1.3 重力荷载代表值】
//  重力荷载代表值 = 恒载 + ψ·活载（ψ = 组合值系数，楼面 0.5 等）。
//  工程软件通常直接让用户把"折算后的各节点重量"注入模型，而不是
//  把恒载工况的等效节点荷载再反着拆回单元 —— 自重已作为线荷载进入
//  恒载工况的 f，若这里再叠 el.mass() 会重复计重。
//
//  因此：
//    · 任一节点有附加重量（gmass 折算）→ 走"节点重量代表团"路径，
//      总质量 = Σ node.weight / g，按节点分配（自重不再重复计入，
//      因为恒载工况里已含 selfweight）；
//    · 全部节点 weight = 0（未使用 gmass）→ 走原有"单元自重"路径，
//      保持既有模型与测试零回归。
std::vector<double> buildMassAggregate(const Model& m, double& total) {
  const Id nTotal = static_cast<Id>(m.nodeCount() * 6);
  std::vector<double> agg(static_cast<size_t>(nTotal), 0.0);
  total = 0.0;
  if (m.hasNodalWeight()) {
    for (Id n = 0; n < m.nodeCount(); ++n) {
      const double w = m.node(n).weight;        // kN（重量）
      if (w <= 0.0) continue;
      const double mt = w / kGravity;           // t（质量）
      total += mt;
      for (int k = 0; k < 3; ++k)               // 三个平动方向都有惯性
        agg[static_cast<size_t>(n * 6 + k)] += mt;
    }
  } else {
    for (const auto& ep : m.elements()) {
      const Element& el = *ep;
      const double weightKn = el.mass();        // kN（重量）
      if (weightKn <= 0.0) continue;
      const double mTot = weightKn / kGravity;  // t（质量）
      total += mTot;
      const int nn = el.nodeCount();
      const double mNode = mTot / static_cast<double>(nn);
      for (Id nd : el.nodes())
        for (int k = 0; k < 3; ++k)
          agg[static_cast<size_t>(nd * 6 + k)] += mNode;
    }
  }
  return agg;
}
}  // namespace

// -----------------------------------------------------------------------------
//  质量组装
// -----------------------------------------------------------------------------
bool ModalAnalysis::assembleMass(Assembly& out, std::string* why) {
  // 1) 先走一遍静力组装 —— 拿 expand 表（与 K 同一套 MPC 口径）
  StaticAnalysis sa(model_);
  StaticAnalysis::Assembly as;
  if (!sa.assemble(as, why)) return false;

  // 2) 质量聚合 → 节点平动自由度（全局自由度号索引）
  //    质量源：gmass 折算的节点重量（恒载+ψ活载）优先，否则单元自重。
  double total = 0.0;
  std::vector<double> agg = buildMassAggregate(model_, total);

  // 3) MPC 折减 → 自由方程号空间（与 K 同维）
  TripletAssembler ta;
  ta.reserve(model_.elementCount() * 16 + 64);
  foldMass(agg, as.expand, ta);

  out.M.buildFrom(ta.triplets(), as.ndof);
  out.ndof = as.ndof;
  out.nnz = out.M.nnz();
  out.totalMass = total;
  return true;
}

// -----------------------------------------------------------------------------
//  模态求解
// -----------------------------------------------------------------------------
ModalResult ModalAnalysis::solve(const Options& opt) {
  ModalResult res;
  const auto t0 = std::chrono::steady_clock::now();

  // 1) 刚度（含 MPC）+ 展开表 —— 与静力分析共用同一来源
  StaticAnalysis sa(model_);
  StaticAnalysis::Assembly as;
  std::string why;
  if (!sa.assemble(as, &why)) {
    res.message = "刚度组装失败：" + why;
    return res;
  }
  res.ndof = as.ndof;
  res.nnzK = as.K.nnz();
  if (as.ndof == 0) {
    res.message = "没有可求解的自由度";
    return res;
  }
  const std::vector<StaticAnalysis::Assembly::Expand>& expand = as.expand;

  // 2) 质量矩阵
  //    质量源：gmass 折算的节点重量（恒载+ψ活载）优先，否则单元自重。
  //    （与 assembleMass 同一逻辑，保证反应谱取 M 的一致性。）
  double totalMassAgg = 0.0;
  std::vector<double> agg = buildMassAggregate(model_, totalMassAgg);
  TripletAssembler ta;
  ta.reserve(model_.elementCount() * 16 + 64);
  foldMass(agg, expand, ta);
  SymSparseMatrix M;
  M.buildFrom(ta.triplets(), as.ndof);
  res.nnzM = M.nnz();

  // 3) 确定模态数：自动 = max(9, 3×层数)，并受"有效质量自由度"钳制
  int heavy = 0;
  for (Id i = 0; i < as.ndof; ++i)
    if (M.at(i, i) > 1e-30) ++heavy;
  if (heavy == 0) {
    res.message = "质量矩阵全零（模型没有质量：请检查单元材料与截面）";
    return res;
  }
  int maxStory = 0;
  for (const Node& n : model_.nodes())
    maxStory = std::max(maxStory, n.story);
  int nmodes = opt.nmodes > 0 ? opt.nmodes : std::max(9, 3 * maxStory);
  nmodes = std::max(1, std::min(nmodes, heavy));
  res.nmodes = nmodes;

  // 4) 子空间迭代解广义特征问题 Kφ = ω²Mφ
  SubspaceIteration::Options so;
  so.tol = opt.tol;
  so.maxIter = opt.maxIter;
  so.workDim = opt.workDim;
  SubspaceIteration::Report rp;
  std::vector<EigenPair> pairs;
  if (!SubspaceIteration::solve(as.K, M, nmodes, pairs, so, &rp)) {
    res.message = "模态求解失败：" + rp.message;
    return res;
  }
  res.iter = rp.iter;
  res.residual = rp.residual;
  res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  // 5) 结果整理：频率/周期/振型回填
  res.omega.resize(static_cast<size_t>(nmodes));
  res.freq.resize(static_cast<size_t>(nmodes));
  res.period.resize(static_cast<size_t>(nmodes));
  res.shapes.resize(static_cast<size_t>(nmodes));
  const Id nNode = model_.nodeCount();

  for (int m = 0; m < nmodes; ++m) {
    const double lam = pairs[static_cast<size_t>(m)].value;   // ω²
    const double w = std::sqrt(std::max(lam, 0.0));
    res.omega[static_cast<size_t>(m)] = w;
    res.freq[static_cast<size_t>(m)] = w / (2.0 * 3.141592653589793);
    res.period[static_cast<size_t>(m)] = (w > 1e-300) ? (2.0 * 3.141592653589793 / w) : 0.0;

    // 自由方程号 → 节点序 × 6
    std::vector<double>& sh = res.shapes[static_cast<size_t>(m)];
    sh.assign(static_cast<size_t>(nNode * 6), 0.0);
    for (Id n = 0; n < nNode; ++n)
      for (int k = 0; k < 6; ++k) {
        const Id g = model_.node(n).dof[k];
        if (g >= 0) sh[static_cast<size_t>(n * 6 + k)] =
                        pairs[static_cast<size_t>(m)].vector[static_cast<size_t>(g)];
      }
    // 从属自由度按展开表折回（与静力位移回填同一公式）
    for (Id n = 0; n < nNode; ++n)
      for (int k = 0; k < 6; ++k) {
        const Id gd = n * 6 + k;
        if (model_.node(n).dof[k] != kDofSlave) continue;
        const auto& e = expand[static_cast<size_t>(gd)];
        double v = 0.0;
        for (int p = 0; p < e.n; ++p)
          v += e.c[p] * pairs[static_cast<size_t>(m)].vector[static_cast<size_t>(e.d[p])];
        sh[static_cast<size_t>(gd)] = v;
      }
  }

  // 6) 模态参与系数与参与质量（三个方向）
  res.gammaX.assign(static_cast<size_t>(nmodes), 0.0);
  res.gammaY.assign(static_cast<size_t>(nmodes), 0.0);
  res.gammaZ.assign(static_cast<size_t>(nmodes), 0.0);
  res.mstarX.assign(static_cast<size_t>(nmodes), 0.0);
  res.mstarY.assign(static_cast<size_t>(nmodes), 0.0);
  res.mstarZ.assign(static_cast<size_t>(nmodes), 0.0);

  // 方向向量 r（自由方程号）：该方向平动自由度置 1。
  // 从属自由度不单独置 1 —— 它们的质量已折入主自由度，主自由度天然代表整层。
  std::vector<double> rX = dirVector(model_, 0);
  std::vector<double> rY = dirVector(model_, 1);
  std::vector<double> rZ = dirVector(model_, 2);
  std::vector<double> MrX, MrY, MrZ;
  M.multiply(rX, MrX);
  M.multiply(rY, MrY);
  M.multiply(rZ, MrZ);

  res.totalMassX = 0.0;
  res.totalMassY = 0.0;
  res.totalMassZ = 0.0;
  for (Id i = 0; i < as.ndof; ++i) {
    res.totalMassX += rX[static_cast<size_t>(i)] * MrX[static_cast<size_t>(i)];
    res.totalMassY += rY[static_cast<size_t>(i)] * MrY[static_cast<size_t>(i)];
    res.totalMassZ += rZ[static_cast<size_t>(i)] * MrZ[static_cast<size_t>(i)];
  }

  for (int m = 0; m < nmodes; ++m) {
    const std::vector<double>& phi = pairs[static_cast<size_t>(m)].vector;
    res.gammaX[static_cast<size_t>(m)] = participation(M, phi, rX);
    res.gammaY[static_cast<size_t>(m)] = participation(M, phi, rY);
    res.gammaZ[static_cast<size_t>(m)] = participation(M, phi, rZ);
    res.mstarX[static_cast<size_t>(m)] = res.gammaX[static_cast<size_t>(m)] * res.gammaX[static_cast<size_t>(m)];
    res.mstarY[static_cast<size_t>(m)] = res.gammaY[static_cast<size_t>(m)] * res.gammaY[static_cast<size_t>(m)];
    res.mstarZ[static_cast<size_t>(m)] = res.gammaZ[static_cast<size_t>(m)] * res.gammaZ[static_cast<size_t>(m)];
    res.massX += res.mstarX[static_cast<size_t>(m)];
    res.massY += res.mstarY[static_cast<size_t>(m)];
    res.massZ += res.mstarZ[static_cast<size_t>(m)];
  }
  res.ratioX = (res.totalMassX > 1e-300) ? res.massX / res.totalMassX : 0.0;
  res.ratioY = (res.totalMassY > 1e-300) ? res.massY / res.totalMassY : 0.0;
  res.ratioZ = (res.totalMassZ > 1e-300) ? res.massZ / res.totalMassZ : 0.0;

  res.residual = rp.residual;
  res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  res.ok = true;
  res.message = "求解成功";
  return res;
}

// -----------------------------------------------------------------------------
//  静态工具
// -----------------------------------------------------------------------------
std::vector<double> ModalAnalysis::dirVector(const Model& m, int comp) {
  std::vector<double> r(static_cast<size_t>(m.freeDofCount()), 0.0);
  for (Id n = 0; n < m.nodeCount(); ++n) {
    const Id g = m.node(n).dof[comp];   // comp ∈ {0,1,2} = ux/uy/uz
    if (g >= 0) r[static_cast<size_t>(g)] = 1.0;
  }
  return r;
}

double ModalAnalysis::participation(const SymSparseMatrix& M,
                                    const std::vector<double>& phi,
                                    const std::vector<double>& r) {
  const size_t n = phi.size();
  if (n == 0 || r.size() != n) return 0.0;
  std::vector<double> Mr;
  M.multiply(r, Mr);
  double s = 0.0;
  for (size_t i = 0; i < n; ++i) s += phi[i] * Mr[i];
  return s;
}

}  // namespace yjk