// =============================================================================
//  tests/test_combination.cpp  ——  荷载组合（LoadCombination）单元测试
//
//  模型：单根竖直悬臂柱（长度 L=4m，b×b=0.3×0.3 方形截面）
//       底节点固定、顶节点自由；局部坐标由 setUpHint({1,0,0}) 指定，
//       方形截面保证 Iy=Iz、Asy=Asz，弯曲方向对手算无歧义。
//  工况（均为 x-y 平面内确定性荷载，用手算解析解对标）：
//     D  ：全柱 x 向均布水平线荷载 q=10 kN/m
//     L  ：顶节点 x 向集中力 50 kN
//     W  ：顶节点 x 向集中力 -30 kN（与 D/L 反号，制造包络张力）
//     Ex ：顶节点 y 向集中力 40 kN（水平地震，校验 x/y 平面独立组合）
//
//  Timoshenko 悬臂梁端部挠度解析解（弯曲项 + 剪切项）：
//     均布 q：w = qL⁴/(8EI) + qL²/(2·G·As)
//     集中 P：w = PL³/(3EI) + PL/(G·As)
//  固端弯矩幅值：均布  |M| = qL²/2；集中  |M| = PL。
//
//  覆盖：
//  A. 便捷构造的规范系数核对（standard/basicVariable/basicPermanent/
//     basicVariable+w/seismic/seismic+ev）
//  B. 单条组合：C1=1.2D+1.4L —— 位移解析对照、固端弯矩幅值、
//     支座反力平衡、逐分量线性自洽（combo == Σ factor·case）、
//     未收敛结果被过滤
//  C. 多条组合 generate：含缺失工况容错（忽略、不影响其余组合）
//  D. 包络 envelope：envMax/envMin 逐分量 == 独立 max/min 扫描；
//     包络统计量 >= 各组合统计量
//  E. 拓扑不一致保护：不同长度结果直接拒绝
// =============================================================================
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "yjk/analysis/LoadCombination.h"
#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/model/Model.h"
#include "yjk/model/Section.h"

using namespace yjk;

static int g_fail = 0;
static int g_pass = 0;

static void check(bool ok, const char* what, const std::string& detail = "") {
  if (ok) {
    ++g_pass;
    std::printf("  [ OK ] %s%s%s\n", what, detail.empty() ? "" : "  ",
                detail.c_str());
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s%s%s\n", what, detail.empty() ? "" : "  ",
                detail.c_str());
  }
}

static bool relErr(double got, double want, double tol, double& err) {
  err = std::abs(got - want) / std::max(std::abs(want), 1e-300);
  return err <= tol;
}

// -----------------------------------------------------------------------------
//  模型 / 解析值
// -----------------------------------------------------------------------------
namespace {

const double L = 4.0;              // m
const double B = 0.3;              // 方形截面边长 m
const double E = 3.15e7;           // kPa —— 必须与 concreteC(30) 完全一致！
const double NU = 0.2;
const double I = B * B * B * B / 12.0;          // 6.75e-4 m⁴
const double A = B * B;                          // 0.09 m²
const double As = 5.0 / 6.0 * A;                 // 0.075 m²
const double G = E / (2.0 * (1.0 + NU));         // 1.3125e7 kPa
const double EI = E * I;                         // kN·m²
const double GAs = G * As;                       // kN

const double q = 10.0;              // kN/m（x 向）
const double P = 50.0;              // kN（x 向）
const double PW = -30.0;            // kN（x 向，风）
const double PE = 40.0;             // kN（y 向，水平地震）

// 端部挠度（m，正方向 = 荷载方向）
double dispUniform(double F) { return F * L * L * L * L / (8.0 * EI) + F * L * L / (2.0 * GAs); }
double dispPoint(double F)  { return F * L * L * L / (3.0 * EI) + F * L / GAs; }

const double uD   = dispUniform(q);                       // 0.0158878
const double uL   = dispPoint(P);                         // 0.0528882
const double uW   = dispPoint(PW);                        // -0.0317330
const double uEx  = dispPoint(PE);                        // 0.0423106 (y)

// 固端弯矩幅值（kN·m）
const double mD  = q * L * L / 2.0;                       // 80
const double mL  = P * L;                                 // 200
const double mW  = std::abs(PW) * L;                      // 120
const double mEx = PE * L;                                // 160

}  // namespace

// -----------------------------------------------------------------------------
//  构建悬臂柱模型并施加某一工况
// -----------------------------------------------------------------------------
static StaticResult solveCase(const char* which) {
  Model m;
  const Id base = m.addNode({0, 0, 0}, 0, "base");
  const Id top  = m.addNode({0, 0, L}, 1, "top");
  SectionProperties sec = section::rect(B, B);
  BeamElement* el = m.addBeam(base, top, sec, Material::concreteC(30));
  el->setUpHint({1, 0, 0});          // 局部 y' = 全局 +x，z' = 全局 +y
  m.fixAll(base);
  m.assignDofs();

  const std::string w(which);
  if (w == "D") {
    el->setLineLoad({q, 0, 0});      // 全柱 x 向均布
  } else if (w == "L") {
    m.addNodeForce(top, {P, 0, 0});
  } else if (w == "W") {
    m.addNodeForce(top, {PW, 0, 0});
  } else if (w == "Ex") {
    m.addNodeForce(top, {0, PE, 0});
  }

  StaticAnalysis sa(m);
  return sa.solve();
}

// 独立逐分量扫描包络（测试自己的 max/min，不依赖被测实现）
static void manualEnvelope(const std::vector<StaticResult>& cs,
                           StaticResult& eMax, StaticResult& eMin) {
  const StaticResult& r0 = cs[0];
  eMax = r0;
  eMin = r0;
  for (size_t c = 1; c < cs.size(); ++c) {
    const StaticResult& r = cs[c];
    for (size_t i = 0; i < eMax.u.size(); ++i) {
      eMax.u[i] = std::max(eMax.u[i], r.u[i]);
      eMin.u[i] = std::min(eMin.u[i], r.u[i]);
    }
    for (size_t i = 0; i < eMax.reaction.size(); ++i) {
      eMax.reaction[i] = std::max(eMax.reaction[i], r.reaction[i]);
      eMin.reaction[i] = std::min(eMin.reaction[i], r.reaction[i]);
    }
    for (size_t e = 0; e < eMax.beamForces.size(); ++e) {
      double* pmax = &eMax.beamForces[e].ef.N;
      double* pmin = &eMin.beamForces[e].ef.N;
      for (int k = 0; k < 12; ++k) {
        const double v = (&r.beamForces[e].ef.N)[k];
        pmax[k] = std::max(pmax[k], v);
        pmin[k] = std::min(pmin[k], v);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//  A. 便捷构造系数
// -----------------------------------------------------------------------------
static void test_defs() {
  std::printf("\n== A. 便捷构造的规范系数 ==\n");

  ComboDef c;

  c = LoadCombination::standard("D", "L");
  check(c.name == "D+L" && c.factors.at("D") == 1.0 && c.factors.at("L") == 1.0,
        "standard = G + Q（1.0/1.0）");

  c = LoadCombination::basicVariable("D", "L");
  check(c.factors.at("D") == combo_coef::kPermanentVar &&
            c.factors.at("L") == combo_coef::kVariable &&
            c.factors.size() == 2,
        "basicVariable = 1.2G + 1.4Q");

  c = LoadCombination::basicPermanent("D", "L");
  check(c.factors.at("D") == combo_coef::kPermanentFixed &&
            c.factors.at("L") == combo_coef::kVariable * combo_coef::kCombFactor,
        "basicPermanent = 1.35G + 1.4·0.7Q");

  c = LoadCombination::basicVariable("D", "L", "W");
  check(c.factors.at("D") == 1.2 && c.factors.at("L") == 1.4 &&
            c.factors.at("W") == 1.4 * 0.7 && c.factors.size() == 3,
        "basicVariable + 风：w 系数 = 1.4·Ψc = 0.98");

  c = LoadCombination::seismic("D", "L", "Ex");
  check(c.factors.at("D") == 1.2 && c.factors.at("L") == 1.2 * 0.5 &&
            c.factors.at("Ex") == 1.3 && c.factors.size() == 3,
        "seismic = 1.2G + 0.6Q + 1.3Eh（ΨE=0.5）");

  c = LoadCombination::seismic("D", "L", "Ex", "Ev");
  check(c.factors.at("Ev") == 0.5 && c.factors.size() == 4,
        "seismic + 竖向地震：Ev 系数 = 0.5");
}

// -----------------------------------------------------------------------------
//  B. 单条组合
// -----------------------------------------------------------------------------
static void test_single_combo() {
  std::printf("\n== B. 单条组合 C1 = 1.2D + 1.4L ==\n");

  const StaticResult rD  = solveCase("D");
  const StaticResult rL  = solveCase("L");
  check(rD.ok && rL.ok, "D / L 工况求解成功", rD.message + " / " + rL.message);
  if (!rD.ok || !rL.ok) return;

  // 工况结果自身对标解析解（正 x 力 → 正 x 位移）
  const double gotD = rD.u[1 * 6 + 0];
  const double gotL = rL.u[1 * 6 + 0];
  double err = 0;
  check(relErr(gotD, uD, 1e-6, err), "D 顶节点 x 位移 == qL⁴/8EI + qL²/2GAs",
        "got " + std::to_string(gotD) + " want " + std::to_string(uD) +
            " err " + std::to_string(err));
  check(relErr(gotL, uL, 1e-6, err), "L 顶节点 x 位移 == PL³/3EI + PL/GAs",
        "got " + std::to_string(gotL) + " want " + std::to_string(uL) +
            " err " + std::to_string(err));

  // 组合 C1 = 1.2D + 1.4L
  std::vector<DesignCase> cases;
  cases.emplace_back("D", rD, combo_coef::kPermanentVar);
  cases.emplace_back("L", rL, combo_coef::kVariable);
  const StaticResult c1 = LoadCombination::combine("1.2D+1.4L", cases);
  check(c1.ok, "combine 返回 ok", c1.message);

  // 1) 位移解析对照
  const double wantU = 1.2 * uD + 1.4 * uL;
  check(relErr(c1.u[6], wantU, 1e-6, err),
        "C1 顶 x 位移 == 1.2·D + 1.4·L（解析）",
        "got " + std::to_string(c1.u[6]) + " want " + std::to_string(wantU));

  // 2) 固端弯矩幅值（x 向弯曲 → 局部 z' 轴弯矩 Mz）
  const double wantM = 1.2 * mD + 1.4 * mL;
  // 如果方向确实落在 Mz 上则校验幅值；否则提示字段位置
  double gotMz = std::abs(c1.beamForces[0].ef.Mz);
  double gotMy = std::abs(c1.beamForces[0].ef.My);
  bool mzCheck = relErr(gotMz, wantM, 1e-6, err);
  bool myCheck = relErr(gotMy, wantM, 1e-6, err);
  const std::string mWhat =
      "C1 固端弯矩幅值 == 1.2·qL²/2 + 1.4·PL（Mz=" +
      std::to_string(gotMz) + " My=" + std::to_string(gotMy) + " want " +
      std::to_string(wantM) + "）";
  check(mzCheck || myCheck, mWhat.c_str(),
        mzCheck ? "(落在 Mz)" : "(落在 My)");
  check(mzCheck && std::abs(gotMy) < 1e-9 * wantM ||
            myCheck && std::abs(gotMz) < 1e-9 * wantM,
        "固端弯矩集中在一个弯曲分量上（另一分量 ≈ 0）");

  // 3) 支座反力平衡：C1 Rx = -(1.2·qL + 1.4·P)
  const double wantRx = -(1.2 * q * L + 1.4 * P);
  check(relErr(c1.reaction[0], wantRx, 1e-6, err),
        "C1 支座水平反力 == -(1.2qL + 1.4P)",
        "got " + std::to_string(c1.reaction[0]) + " want " +
            std::to_string(wantRx));

  // 4) 逐分量线性自洽：combo == Σ factor·case（u + reaction + beamForces 全部）
  bool lin = true;
  for (size_t i = 0; i < c1.u.size(); ++i) {
    const double w = 1.2 * rD.u[i] + 1.4 * rL.u[i];
    if (std::abs(c1.u[i] - w) > 1e-9 * std::max(std::abs(w), 1e-9)) lin = false;
  }
  for (size_t i = 0; i < c1.reaction.size(); ++i) {
    const double w = 1.2 * rD.reaction[i] + 1.4 * rL.reaction[i];
    if (std::abs(c1.reaction[i] - w) > 1e-9 * std::max(std::abs(w), 1e-9))
      lin = false;
  }
  const double* p0 = &c1.beamForces[0].ef.N;
  const double* pD = &rD.beamForces[0].ef.N;
  const double* pL = &rL.beamForces[0].ef.N;
  for (int k = 0; k < 12; ++k) {
    const double w = 1.2 * pD[k] + 1.4 * pL[k];
    if (std::abs(p0[k] - w) > 1e-9 * std::max(std::abs(w), 1e-9)) lin = false;
  }
  check(lin, "C1 逐分量 == 1.2·D + 1.4·L（u/reaction/beamForces 全分量）");

  // 5) 统计量重算：maxDisplacement 应为顶节点位移幅值
  check(relErr(c1.maxDisplacement, std::abs(wantU), 1e-6, err),
        "C1 maxDisplacement == 顶节点位移幅值",
        "got " + std::to_string(c1.maxDisplacement));

  // 6) 未收敛结果被过滤
  StaticResult bad = rD;
  bad.ok = false;
  std::vector<DesignCase> cBad;
  cBad.emplace_back("D", bad, 1.2);
  cBad.emplace_back("L", rL, 1.4);
  StaticResult cb = LoadCombination::combine("bad", cBad);
  check(cb.ok && std::abs(cb.u[6] - 1.4 * uL) < 1e-9,
        "未收敛工况被过滤，组合退化为剩余工况");
  std::vector<DesignCase> cNone;
  cNone.emplace_back("X", bad, 1.2);
  StaticResult cn = LoadCombination::combine("none", cNone);
  check(!cn.ok && cn.message.find("没有可用") != std::string::npos,
        "全部无效时返回 !ok 并说明原因", cn.message);
}

// -----------------------------------------------------------------------------
//  C. generate 多条组合 + 缺失工况容错
// -----------------------------------------------------------------------------
static void test_generate() {
  std::printf("\n== C. generate 多条组合 ==\n");

  const StaticResult rD  = solveCase("D");
  const StaticResult rL  = solveCase("L");
  const StaticResult rW  = solveCase("W");
  const StaticResult rEx = solveCase("Ex");

  std::map<std::string, StaticResult> cases;
  cases["D"] = rD;
  cases["L"] = rL;
  cases["W"] = rW;
  cases["Ex"] = rEx;

  std::vector<ComboDef> defs;
  defs.push_back(LoadCombination::standard("D", "L"));
  defs.push_back(LoadCombination::basicVariable("D", "L", "W"));
  defs.push_back(LoadCombination::seismic("D", "L", "Ex"));
  // 组合定义引用不存在的工况 "Wind" → 应被忽略，不影响其余组合
  defs.push_back(LoadCombination::def("with-missing", {{"D", 1.0}, {"Wind", 0.9}}));

  const std::vector<StaticResult> out = LoadCombination::generate(cases, defs);
  check(out.size() == 4, "4 条定义全部生成", "got " + std::to_string(out.size()));

  if (out.size() == 4) {
    double err = 0;
    check(relErr(out[0].u[6], uD + uL, 1e-6, err),
          "standard：顶 x 位移 == D + L",
          std::to_string(out[0].u[6]));
    check(relErr(out[1].u[6], 1.2 * uD + 1.4 * uL + 1.4 * 0.7 * uW, 1e-6, err),
          "basicVariable+W：顶 x 位移 == 1.2D + 1.4L + 0.98W");
    check(std::abs(out[1].u[6] - (1.2 * uD + 1.4 * uL + 0.98 * uW)) < 1e-6,
          "（风工况取负向荷载，包络张力成立）");
    check(relErr(out[2].u[6], 1.2 * uD + 0.6 * uL, 1e-6, err) &&
              relErr(out[2].u[7], 1.3 * uEx, 1e-6, err),
          "seismic：顶 x = 1.2D+0.6L，顶 y = 1.3Ex（x/y 平面独立）",
          "x=" + std::to_string(out[2].u[6]) + " y=" + std::to_string(out[2].u[7]));
    check(relErr(out[3].u[6], uD, 1e-6, err),
          "缺失工况被忽略：with-missing == D", "got " + std::to_string(out[3].u[6]));
  }

  // 空定义 / 空结果集
  check(LoadCombination::generate(cases, {}).empty(), "空定义列表返回空");
  std::map<std::string, StaticResult> emptyMap;
  check(LoadCombination::generate(emptyMap, defs).empty(),
        "空结果集：依赖工况的组合全部被跳过，返回空");
}

// -----------------------------------------------------------------------------
//  D. 包络
// -----------------------------------------------------------------------------
static void test_envelope() {
  std::printf("\n== D. 包络 ==\n");

  const StaticResult rD  = solveCase("D");
  const StaticResult rL  = solveCase("L");
  const StaticResult rW  = solveCase("W");
  const StaticResult rEx = solveCase("Ex");

  std::map<std::string, StaticResult> cases;
  cases["D"] = rD;
  cases["L"] = rL;
  cases["W"] = rW;
  cases["Ex"] = rEx;

  std::vector<ComboDef> defs;
  defs.push_back(LoadCombination::basicVariable("D", "L"));
  defs.push_back(LoadCombination::basicPermanent("D", "L"));
  defs.push_back(LoadCombination::basicVariable("D", "L", "W"));
  defs.push_back(LoadCombination::seismic("D", "L", "Ex"));

  const std::vector<StaticResult> combos = LoadCombination::generate(cases, defs);
  check(combos.size() == 4, "4 条组合生成成功");

  StaticResult envMax, envMin;
  LoadCombination::envelope(combos, envMax, envMin);

  // 1) 逐分量与独立扫描对比
  StaticResult refMax, refMin;
  manualEnvelope(combos, refMax, refMin);
  bool sameU = true, sameR = true, sameB = true;
  for (size_t i = 0; i < envMax.u.size(); ++i)
    if (envMax.u[i] != refMax.u[i] || envMin.u[i] != refMin.u[i]) sameU = false;
  for (size_t i = 0; i < envMax.reaction.size(); ++i)
    if (envMax.reaction[i] != refMax.reaction[i] ||
        envMin.reaction[i] != refMin.reaction[i]) sameR = false;
  for (size_t e = 0; e < envMax.beamForces.size(); ++e) {
    for (int k = 0; k < 12; ++k) {
      if ((&envMax.beamForces[e].ef.N)[k] != (&refMax.beamForces[e].ef.N)[k] ||
          (&envMin.beamForces[e].ef.N)[k] != (&refMin.beamForces[e].ef.N)[k])
        sameB = false;
    }
  }
  check(sameU, "包络 u 逐分量 == 独立 max/min 扫描");
  check(sameR, "包络 reaction 逐分量 == 独立 max/min 扫描");
  check(sameB, "包络 beamForces 逐分量 == 独立 max/min 扫描");

  // 2) 手算对照：顶节点 x 位移包络
  const double c1u = 1.2 * uD + 1.4 * uL;
  const double c2u = 1.35 * uD + 0.98 * uL;
  const double c3u = 1.2 * uD + 1.4 * uL + 0.98 * uW;
  const double c4u = 1.2 * uD + 0.6 * uL;
  const double wantMax = std::max({c1u, c2u, c3u, c4u});
  const double wantMin = std::min({c1u, c2u, c3u, c4u});
  double err = 0;
  check(relErr(envMax.u[6], wantMax, 1e-6, err),
        "包络极大 x 位移 == max(4 条组合)",
        "got " + std::to_string(envMax.u[6]) + " want " +
            std::to_string(wantMax));
  check(relErr(envMin.u[6], wantMin, 1e-6, err),
        "包络极小 x 位移 == min(4 条组合)",
        "got " + std::to_string(envMin.u[6]) + " want " +
            std::to_string(wantMin));
  check(relErr(envMax.u[7], 1.3 * uEx, 1e-6, err),
        "包络极大 y 位移 == 1.3Ex（仅地震组合贡献）");

  // 3) 包络统计量 >= 各组合统计量（单调性）
  bool mono = true;
  for (const StaticResult& c : combos) {
    if (envMax.maxDisplacement + 1e-12 < c.maxDisplacement) mono = false;
    if (envMax.maxUxy + 1e-12 < c.maxUxy) mono = false;
    if (envMax.maxReaction + 1e-12 < c.maxReaction) mono = false;
  }
  check(mono, "包络统计量覆盖每条组合（max ≥ 各组合）");

  // 4) 空输入
  StaticResult e0, e1;
  LoadCombination::envelope({}, e0, e1);
  check(!e0.ok && !e1.ok, "空输入包络返回 !ok");
}

// -----------------------------------------------------------------------------
//  E. 拓扑不一致保护
// -----------------------------------------------------------------------------
static void test_mismatch() {
  std::printf("\n== E. 拓扑不一致保护 ==\n");

  const StaticResult rL = solveCase("L");
  const StaticResult rEx = solveCase("Ex");

  // 不同模型（柱长不同）→ 自由度长度一致但几何不同：
  // 用两个不同长度柱模型的工况结果强行组合 —— 长度一致时线性组合数学上
  // 总是成立，这不是本保护的目标；真正要拦的是"长度已不一致"。
  // 构造一个截断的 StaticResult 模拟不同模型：
  StaticResult shortR = rL;
  shortR.u.resize(shortR.u.size() + 6);   // 节点数不同
  std::vector<DesignCase> cases;
  cases.emplace_back("L", rL, 1.4);
  cases.emplace_back("Ex", shortR, 1.3);
  StaticResult r = LoadCombination::combine("bad-topology", cases);
  check(!r.ok && r.message.find("长度不一致") != std::string::npos,
        "长度不一致 ⟹ 拒绝组合", r.message);

  // 单独一条合法工况仍可组合
  std::vector<DesignCase> single;
  single.emplace_back("Ex", rEx, combo_coef::kGamEh);
  StaticResult r2 = LoadCombination::combine("1.3Ex", single);
  double err = 0;
  check(r2.ok && relErr(r2.u[7], 1.3 * uEx, 1e-6, err),
        "单一工况组合（1.3Ex）仍正确", r2.message);
}

int main() {
  std::printf("=== test_combination: 荷载组合 ===\n");
  test_defs();
  test_single_combo();
  test_generate();
  test_envelope();
  test_mismatch();
  std::printf("\n通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}