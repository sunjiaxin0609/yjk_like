// =============================================================================
//  src/analysis/LoadCombination.cpp  ——  荷载组合实现
// =============================================================================
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "yjk/analysis/LoadCombination.h"

namespace yjk {

namespace {

// 组合后从 u/reaction 重新扫描统计极值（不依赖 Model，纯向量运算）。
void recomputeStats(StaticResult& r) {
  r.maxDisplacement = 0.0;
  r.maxUxy = 0.0;
  r.maxUz = 0.0;
  r.maxRotation = 0.0;
  r.maxReaction = 0.0;
  r.maxDispNode = -1;
  r.maxDispElem = -1;

  const size_t nNode = r.u.size() / 6;
  for (size_t n = 0; n < nNode; ++n) {
    const double dx = r.u[n * 6 + 0];
    const double dy = r.u[n * 6 + 1];
    const double dz = r.u[n * 6 + 2];
    const double rx = r.u[n * 6 + 3];
    const double ry = r.u[n * 6 + 4];
    const double rz = r.u[n * 6 + 5];
    const double dm = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (dm > r.maxDisplacement) {
      r.maxDisplacement = dm;
      r.maxDispNode = static_cast<Id>(n);
    }
    r.maxUxy = std::max(r.maxUxy, std::hypot(dx, dy));
    r.maxUz = std::max(r.maxUz, std::abs(dz));
    r.maxRotation =
        std::max(r.maxRotation, std::sqrt(rx * rx + ry * ry + rz * rz));
  }
  for (double v : r.reaction)
    r.maxReaction = std::max(r.maxReaction, std::abs(v));
}

// 从组合后的面内应力重算主应力与 von Mises（sx/sy/txy 线性组合后
// 主应力是导出量，不能直接相加 —— 必须重算）。
void recomputePrincipal(ShellElement4::MembraneForces& mf) {
  const double sx = mf.sx, sy = mf.sy, txy = mf.txy;
  const double c = 0.5 * (sx + sy);
  const double rad =
      std::sqrt(0.25 * (sx - sy) * (sx - sy) + txy * txy);
  mf.s1 = c + rad;
  mf.s2 = c - rad;
  mf.vonMises = std::sqrt(std::max(
      sx * sx - sx * sy + sy * sy + 3.0 * txy * txy, 0.0));
}

// 以 r0 为模板初始化一个"零组合结果"（拓扑字段取 r0，数值归零）。
void initFrom(const StaticResult& r0, StaticResult& out) {
  out.ok = false;
  out.ndof = r0.ndof;
  out.nconstrained = r0.nconstrained;
  out.nnz = r0.nnz;
  out.nnzL = r0.nnzL;
  out.seconds = 0.0;                 // 组合是纯代数，不计耗时
  out.residual = 0.0;                // 组合结果没有单一残差意义
  out.u.assign(r0.u.size(), 0.0);
  out.reaction.assign(r0.reaction.size(), 0.0);
  out.beamForces.assign(r0.beamForces.size(), StaticResult::BeamForce{});
  out.shellForces.assign(r0.shellForces.size(), StaticResult::ShellForce{});
}

}  // namespace

// -----------------------------------------------------------------------------
//  核心：逐分量线性组合
// -----------------------------------------------------------------------------
StaticResult LoadCombination::combine(const std::string& name,
                                      const std::vector<DesignCase>& cases) {
  StaticResult out;
  out.message = name.empty() ? "组合" : name;

  // 过滤空引用 / 未收敛结果 / 零系数 / 空名
  std::vector<const DesignCase*> refs;
  for (const DesignCase& c : cases) {
    if (!c.result || !c.result->ok || c.factor == 0.0) continue;
    if (c.name.empty()) continue;
    refs.push_back(&c);
  }
  if (refs.empty()) {
    out.message += "：没有可用的工况结果";
    return out;
  }

  // 结构一致性：同模型同拓扑 ⇒ 位移/反力/单元结果长度一致
  const StaticResult& r0 = *refs[0]->result;
  for (const DesignCase* cp : refs) {
    const StaticResult& r = *cp->result;
    if (r.u.size() != r0.u.size() || r.reaction.size() != r0.reaction.size() ||
        r.beamForces.size() != r0.beamForces.size() ||
        r.shellForces.size() != r0.shellForces.size()) {
      out.message += "：工况结果长度不一致（不同模型/拓扑？）";
      return out;
    }
  }

  initFrom(r0, out);

  // 位移（节点序 × 6）
  for (size_t i = 0; i < out.u.size(); ++i) {
    double s = 0.0;
    for (const DesignCase* cp : refs)
      s += cp->factor * cp->result->u[i];
    out.u[i] = s;
  }
  // 支座反力
  for (size_t i = 0; i < out.reaction.size(); ++i) {
    double s = 0.0;
    for (const DesignCase* cp : refs)
      s += cp->factor * cp->result->reaction[i];
    out.reaction[i] = s;
  }

  // 梁端内力：EndForces 为 12 个连续 double（i/j 端各 6 个）
  static_assert(sizeof(BeamElement3D::EndForces) == 12 * sizeof(double),
                "EndForces 必须为 12 个连续标量");
  for (size_t e = 0; e < out.beamForces.size(); ++e) {
    out.beamForces[e].xi = r0.beamForces[e].xi;
    double* pd = &out.beamForces[e].ef.N;
    for (int k = 0; k < 12; ++k) {
      double s = 0.0;
      for (const DesignCase* cp : refs)
        s += cp->factor * (&cp->result->beamForces[e].ef.N)[k];
      pd[k] = s;
    }
  }

  // 壳内力：sx/sy/txy 线性组合，主应力与 von Mises 重算
  for (size_t e = 0; e < out.shellForces.size(); ++e) {
    out.shellForces[e].area = r0.shellForces[e].area;
    out.shellForces[e].normal = r0.shellForces[e].normal;
    ShellElement4::MembraneForces& m = out.shellForces[e].mf;
    m.sx = m.sy = m.txy = 0.0;
    for (const DesignCase* cp : refs) {
      const ShellElement4::MembraneForces& mi = cp->result->shellForces[e].mf;
      m.sx += cp->factor * mi.sx;
      m.sy += cp->factor * mi.sy;
      m.txy += cp->factor * mi.txy;
    }
    recomputePrincipal(m);
  }

  recomputeStats(out);
  out.ok = true;
  out.message = name.empty() ? "组合" : name;   // 保持原消息（首工况名不混入）
  return out;
}

// -----------------------------------------------------------------------------
//  多条组合
// -----------------------------------------------------------------------------
std::vector<StaticResult> LoadCombination::generate(
    const std::map<std::string, StaticResult>& caseResults,
    const std::vector<ComboDef>& defs) {
  std::vector<StaticResult> out;
  for (const ComboDef& d : defs) {
    std::vector<DesignCase> cs;
    for (const auto& kv : d.factors) {
      auto it = caseResults.find(kv.first);
      if (it == caseResults.end() || kv.second == 0.0) continue;
      cs.emplace_back(kv.first, it->second, kv.second);
    }
    StaticResult sr = combine(d.name, cs);
    if (sr.ok) out.push_back(std::move(sr));
  }
  return out;
}

// -----------------------------------------------------------------------------
//  包络：逐分量取极大/极小
// -----------------------------------------------------------------------------
void LoadCombination::envelope(const std::vector<StaticResult>& combos,
                               StaticResult& envMax, StaticResult& envMin) {
  envMax = StaticResult{};
  envMin = StaticResult{};
  if (combos.empty()) return;

  const StaticResult& r0 = combos[0];
  envMax = r0;
  envMin = r0;
  envMax.ok = envMin.ok = true;
  envMax.message = "包络（极大）";
  envMin.message = "包络（极小）";
  envMax.seconds = envMin.seconds = 0.0;
  envMax.residual = envMin.residual = 0.0;

  for (size_t c = 1; c < combos.size(); ++c) {
    const StaticResult& r = combos[c];
    const size_t nu = std::min(envMax.u.size(), r.u.size());
    for (size_t i = 0; i < nu; ++i) {
      envMax.u[i] = std::max(envMax.u[i], r.u[i]);
      envMin.u[i] = std::min(envMin.u[i], r.u[i]);
    }
    const size_t nr = std::min(envMax.reaction.size(), r.reaction.size());
    for (size_t i = 0; i < nr; ++i) {
      envMax.reaction[i] = std::max(envMax.reaction[i], r.reaction[i]);
      envMin.reaction[i] = std::min(envMin.reaction[i], r.reaction[i]);
    }
    const size_t nb = std::min(envMax.beamForces.size(), r.beamForces.size());
    for (size_t e = 0; e < nb; ++e) {
      double* pmax = &envMax.beamForces[e].ef.N;
      double* pmin = &envMin.beamForces[e].ef.N;
      for (int k = 0; k < 12; ++k) {
        const double v = (&r.beamForces[e].ef.N)[k];
        pmax[k] = std::max(pmax[k], v);
        pmin[k] = std::min(pmin[k], v);
      }
    }
    const size_t ns = std::min(envMax.shellForces.size(), r.shellForces.size());
    for (size_t e = 0; e < ns; ++e) {
      const auto& m = r.shellForces[e].mf;
      envMax.shellForces[e].mf.sx = std::max(envMax.shellForces[e].mf.sx, m.sx);
      envMin.shellForces[e].mf.sx = std::min(envMin.shellForces[e].mf.sx, m.sx);
      envMax.shellForces[e].mf.sy = std::max(envMax.shellForces[e].mf.sy, m.sy);
      envMin.shellForces[e].mf.sy = std::min(envMin.shellForces[e].mf.sy, m.sy);
      envMax.shellForces[e].mf.txy = std::max(envMax.shellForces[e].mf.txy, m.txy);
      envMin.shellForces[e].mf.txy = std::min(envMin.shellForces[e].mf.txy, m.txy);
    }
  }

  // 包络后主应力/统计量重算（sx/sy/txy 的 max/min 组合不再对应同一应力态，
  // 主应力与 von Mises 从各分量的包络值重算，保证自洽可显示）。
  for (StaticResult* env : {&envMax, &envMin}) {
    for (auto& sf : env->shellForces) recomputePrincipal(sf.mf);
    recomputeStats(*env);
  }

  // 统计极值 = 各组合统计值的包络（位移幅值等为标量极值，无正负）
  for (const StaticResult& r : combos) {
    envMax.maxDisplacement = std::max(envMax.maxDisplacement, r.maxDisplacement);
    envMax.maxUxy = std::max(envMax.maxUxy, r.maxUxy);
    envMax.maxUz = std::max(envMax.maxUz, r.maxUz);
    envMax.maxRotation = std::max(envMax.maxRotation, r.maxRotation);
    envMax.maxReaction = std::max(envMax.maxReaction, r.maxReaction);
    envMax.maxDispNode = std::max(envMax.maxDispNode, r.maxDispNode);
    envMin.maxDisplacement = std::max(envMin.maxDisplacement, r.maxDisplacement);
    envMin.maxUxy = std::max(envMin.maxUxy, r.maxUxy);
    envMin.maxUz = std::max(envMin.maxUz, r.maxUz);
    envMin.maxRotation = std::max(envMin.maxRotation, r.maxRotation);
    envMin.maxReaction = std::max(envMin.maxReaction, r.maxReaction);
    envMin.maxDispNode = std::max(envMin.maxDispNode, r.maxDispNode);
  }
}

// -----------------------------------------------------------------------------
//  便捷构造
// -----------------------------------------------------------------------------
ComboDef LoadCombination::def(
    const std::string& name,
    std::initializer_list<std::pair<const char*, double>> fs) {
  ComboDef d;
  d.name = name;
  for (const auto& kv : fs) d.factors[kv.first] = kv.second;
  return d;
}

ComboDef LoadCombination::standard(const std::string& g, const std::string& q) {
  return def(g + "+" + q, {{g.c_str(), 1.0}, {q.c_str(), 1.0}});
}

ComboDef LoadCombination::basicVariable(const std::string& g,
                                        const std::string& q,
                                        const std::string& w) {
  ComboDef d = def("1.2" + g + "+1.4" + q,
                   {{g.c_str(), combo_coef::kPermanentVar},
                    {q.c_str(), combo_coef::kVariable}});
  if (!w.empty())
    d.factors[w] = combo_coef::kVariable * combo_coef::kCombFactor;
  return d;
}

ComboDef LoadCombination::basicPermanent(const std::string& g,
                                         const std::string& q,
                                         const std::string& w) {
  ComboDef d = def("1.35" + g + "+1.4x0.7" + q,
                   {{g.c_str(), combo_coef::kPermanentFixed},
                    {q.c_str(), combo_coef::kVariable * combo_coef::kCombFactor}});
  if (!w.empty())
    d.factors[w] = combo_coef::kVariable * combo_coef::kCombFactor;
  return d;
}

ComboDef LoadCombination::seismic(const std::string& g, const std::string& q,
                                  const std::string& eh, const std::string& ev) {
  ComboDef d = def("1.2(G+0.5Q)+1.3" + eh + (ev.empty() ? "" : "+0.5" + ev),
                   {{g.c_str(), combo_coef::kGamG},
                    {q.c_str(), combo_coef::kGamG * combo_coef::kPsiE},
                    {eh.c_str(), combo_coef::kGamEh}});
  if (!ev.empty()) d.factors[ev] = combo_coef::kGamEv;
  return d;
}

}  // namespace yjk