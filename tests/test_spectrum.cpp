// =============================================================================
//  tests/test_spectrum.cpp  ——  设计反应谱与振型分解反应谱法单元测试
//
//  覆盖：
//  A. 反应谱曲线四段解析值（上升段/平台段/曲线下降段/直线下降段，
//     与 5Tg 过渡点连续性、T>6s 截断）+ 阻尼调整系数 γ/η1/η2 下限。
//  B. 单自由度解析解：单柱轴向模型，V = α(T)·m*·g，等效地震力顶节点
//     集中力 = 基底剪力（SRSS 单模态下闭合）。
//  C. 两自由度 SRSS：两柱串列轴向，逐模态 V_j = α(T_j)·γ_j²·g 对标，
//     组合值与 SqrtΣV_j² 对标（同时验证 m* = γ² 换算与组合逻辑）。
//  D. CQC：ρ_jk 数值解（β=1→1、远频→0、规范式闭合形式），
//     baseZ = √(V1²+V2²+2ρV1V2) 对标，method=1 与静态 cqc() 一致。
//  E. 参与质量 ≥ 90% 检查：单轴向模型（X/Y 方向无可动质量应豁免）、
//     刚性楼板模型取 1 阶模态（Z 向轴向质量未参与 → 检查应失败）、
//     取全模态通过；checkMass=false 时不阻断。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "yjk/analysis/ResponseSpectrum.h"
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

static Material makeConcrete() {
  Material m;
  m.gamma = 25.0;   // kN/m³
  return m;
}

static constexpr double kG = 9.81;

// -----------------------------------------------------------------------------
//  A. 反应谱曲线解析值
// -----------------------------------------------------------------------------
static void test_spectrum_curve() {
  std::printf("\n== A. 反应谱曲线（ζ = 0.05 默认谱） ==\n");
  const DesignSpectrum sp;   // αmax=0.16, Tg=0.40, ζ=0.05

  double err = 0.0;
  // 形状参数：γ = 0.9（ζ=0.05），η1 = 0.02，η2 = 1.0
  check(relErr(sp.gammaExp(), 0.9, 1e-12, err), "γ(ζ=0.05) = 0.9");
  check(relErr(sp.eta1(), 0.02, 1e-12, err), "η1(ζ=0.05) = 0.02");
  check(relErr(sp.eta2(), 1.0, 1e-12, err), "η2(ζ=0.05) = 1.0");

  // 上升段：α = [0.45 + 10(η2−0.45)T]·αmax，T=0.01
  const double aUp = (0.45 + 10.0 * (1.0 - 0.45) * 0.01) * 0.16;
  check(relErr(sp.alpha(0.01), aUp, 1e-12, err), "上升段 α(0.01)",
        "got " + std::to_string(sp.alpha(0.01)));

  // T=0.1 上升段末 = 平台起点：α = η2·αmax = 0.16
  check(relErr(sp.alpha(0.1), 0.16, 1e-12, err), "平台起点 α(0.1) = 0.16");

  // 平台段：T = Tg = 0.4
  check(relErr(sp.alpha(0.4), 0.16, 1e-12, err), "平台段 α(Tg) = 0.16");

  // 曲线下降段：α = (Tg/T)^γ · η2·αmax，T = 1.0
  const double aDown = std::pow(0.4 / 1.0, 0.9) * 0.16;
  check(relErr(sp.alpha(1.0), aDown, 1e-12, err), "曲线下降段 α(1.0)",
        "got " + std::to_string(sp.alpha(1.0)));

  // 5Tg = 2.0：曲线段末 = 直线段起点，α = 0.2^γ·η2·αmax
  const double a5 = std::pow(0.2, 0.9) * 0.16;
  check(relErr(sp.alpha(2.0), a5, 1e-12, err), "5Tg 过渡点 α(2.0)",
        "got " + std::to_string(sp.alpha(2.0)));

  // 直线下降段：α = [η2·0.2^γ − η1(T − 5Tg)]·αmax，T = 3.0
  const double aLine = (std::pow(0.2, 0.9) - 0.02 * (3.0 - 2.0)) * 0.16;
  check(relErr(sp.alpha(3.0), aLine, 1e-12, err), "直线下降段 α(3.0)",
        "got " + std::to_string(sp.alpha(3.0)));

  // T > 6s：按 T = 6 截断（规范曲线终点）
  const double a6 = (std::pow(0.2, 0.9) - 0.02 * (6.0 - 2.0)) * 0.16;
  check(relErr(sp.alpha(8.0), a6, 1e-12, err), "T=8 按曲线终点 T=6 截断",
        "got " + std::to_string(sp.alpha(8.0)));
  check(std::abs(sp.alpha(8.0) - sp.alpha(6.0)) < 1e-15, "α(8) == α(6)");

  // 阻尼调整：ζ = 0.02（钢结构）→ η2 提高
  DesignSpectrum st;
  st.zeta = 0.02;
  const double eta2s = 1.0 + (0.05 - 0.02) / (0.08 + 1.6 * 0.02);
  check(relErr(st.eta2(), eta2s, 1e-12, err), "η2(ζ=0.02) 提高",
        "got " + std::to_string(st.eta2()));

  // 大阻尼下限：ζ = 0.4 → η2 公式值 < 0.55 应钳到 0.55；
  //                 η1 公式值 < 0 应钳到 0
  DesignSpectrum sbig;
  sbig.zeta = 0.4;
  const double eta1raw = 0.02 + (0.05 - 0.4) / (4.0 + 32.0 * 0.4);   // < 0
  const double eta2raw = 1.0 + (0.05 - 0.4) / (0.08 + 1.6 * 0.4);    // < 0.55
  check(eta1raw < 0.0 && sbig.eta1() == 0.0, "η1(ζ=0.4) 钳到 0",
        "raw " + std::to_string(eta1raw));
  check(eta2raw < 0.55 && std::abs(sbig.eta2() - 0.55) < 1e-12,
        "η2(ζ=0.4) 钳到 0.55");
}

// -----------------------------------------------------------------------------
//  B. 单自由度解析解
// -----------------------------------------------------------------------------
static void test_single_dof() {
  std::printf("\n== B. 单自由度解析解 ==\n");
  Model m;
  const double L = 4.0, A = 0.04, E = 3.0e7;
  SectionProperties sec;
  sec.A = A;
  sec.Iy = sec.Iz = 1.3333e-3;
  sec.J = 0.0;
  const Material mat = makeConcrete();
  const Id n0 = m.addNode({0, 0, 0}, 0, "base");
  const Id n1 = m.addNode({0, 0, L}, 1, "top");
  BeamElement* el = m.addBeam(n0, n1, sec, mat);
  el->setUpHint({1, 0, 0});
  m.fixAll(n0);
  m.fixNode(n1, true, true, false, true, true, true);   // 只留 uz
  m.assignDofs();

  ModalAnalysis ma(m);
  ModalResult mr = ma.solve(ModalAnalysis::Options{1, 1e-10, 300, 0});
  check(mr.ok, mr.ok ? "模态求解成功" : ("失败：" + mr.message).c_str());
  if (!mr.ok) return;

  const double mTop = 0.5 * mat.gamma * sec.A * L / kG;   // t
  const double kAx = E * sec.A / L;
  const double wantW = std::sqrt(kAx / mTop);
  double err = 0.0;
  check(relErr(mr.omega[0], wantW, 5e-3, err), "基频对标解析解");

  // 反应谱：水平方向谱 → Z 向取 verticalScale = 1.0（竖谱与水平同值，便于对标）
  DesignSpectrum spec;
  ResponseSpectrum rs(m);
  ResponseSpectrum::Options opt;
  opt.method = 0;             // SRSS
  opt.verticalScale = 1.0;
  const SpectrumResult r = rs.solve(mr, spec, opt);
  check(r.ok, r.ok ? "反应谱求解成功" : ("失败：" + r.message).c_str());
  if (!r.ok) return;

  // 单自由度：T = 2π/ω，α = α(T)，V = α·m*·g = α·mTop·g
  const double T = 2.0 * 3.141592653589793 / mr.omega[0];
  const double alpha = spec.alpha(T);
  const double wantV = alpha * mTop * kG;
  check(relErr(r.alphaJ[0], alpha, 1e-12, err), "α_j = α(T_j)",
        "got " + std::to_string(r.alphaJ[0]));
  check(relErr(r.baseZ, wantV, 1e-9, err), "基底剪力 V = α·m*·g",
        "got " + std::to_string(r.baseZ) + " want " + std::to_string(wantV));

  // 等效地震力：顶节点 uz 集中力 = 基底剪力（单模态 SRSS 闭合）
  const double fTop = r.forceZ[static_cast<size_t>(n1 * 6 + 2)];
  check(relErr(fTop, wantV, 1e-9, err), "顶节点等效地震力 = 基底剪力",
        "got " + std::to_string(fTop) + " want " + std::to_string(wantV));
  const double fBot = r.forceZ[static_cast<size_t>(n0 * 6 + 2)];
  check(std::abs(fBot) < 1e-12, "底部固定节点等效地震力 = 0");

  // 质量参与：Z 向 ratio = 1；X/Y 向无可动质量应豁免（不误报失败）
  check(r.massRatioZ > 0.999, "Z 向参与质量比 = 1",
        "got " + std::to_string(r.massRatioZ));
  check(r.massOk, "参与质量检查通过（Z=1，X/Y 豁免）",
        r.message.c_str());
}

// -----------------------------------------------------------------------------
//  C. 两自由度 SRSS
// -----------------------------------------------------------------------------
static void test_two_dof_srss() {
  std::printf("\n== C. 两自由度 SRSS ==\n");
  Model m;
  const double L = 2.0, A = 0.04, E = 3.0e7;
  SectionProperties sec;
  sec.A = A;
  sec.Iy = sec.Iz = 1.3333e-3;
  sec.J = 0.0;
  const Material mat = makeConcrete();
  const Id n0 = m.addNode({0, 0, 0}, 0);
  const Id n1 = m.addNode({0, 0, L}, 1);
  const Id n2 = m.addNode({0, 0, 2 * L}, 2);
  BeamElement* e1 = m.addBeam(n0, n1, sec, mat);
  e1->setUpHint({1, 0, 0});
  BeamElement* e2 = m.addBeam(n1, n2, sec, mat);
  e2->setUpHint({1, 0, 0});
  m.fixAll(n0);
  m.fixNode(n1, true, true, false, true, true, true);
  m.fixNode(n2, true, true, false, true, true, true);
  m.assignDofs();

  ModalAnalysis ma(m);
  ModalResult mr = ma.solve(ModalAnalysis::Options{2, 1e-10, 300, 0});
  check(mr.ok, mr.ok ? "模态求解成功" : ("失败：" + mr.message).c_str());
  if (!mr.ok) return;

  DesignSpectrum spec;
  ResponseSpectrum rs(m);
  ResponseSpectrum::Options opt;
  opt.method = 0;
  opt.verticalScale = 1.0;
  const SpectrumResult r = rs.solve(mr, spec, opt);
  check(r.ok, r.ok ? "反应谱求解成功" : ("失败：" + r.message).c_str());
  if (!r.ok) return;

  // 逐模态：V_j = α(T_j)·γ_j²·g，对标实现输出的 vZ[j]
  double err = 0.0;
  double srssWant = 0.0;
  for (int j = 0; j < 2; ++j) {
    const double Tj = 2.0 * 3.141592653589793 / mr.omega[static_cast<size_t>(j)];
    const double aj = spec.alpha(Tj);
    const double gj = mr.gammaZ[static_cast<size_t>(j)];
    const double Vj = aj * gj * gj * kG;
    check(relErr(r.vZ[static_cast<size_t>(j)], Vj, 1e-9, err),
          ("模态 " + std::to_string(j + 1) + " 基底剪力 V = α·γ²·g").c_str(),
          "got " + std::to_string(r.vZ[static_cast<size_t>(j)]) +
              " want " + std::to_string(Vj));
    srssWant += Vj * Vj;
  }
  srssWant = std::sqrt(srssWant);
  check(relErr(r.baseZ, srssWant, 1e-9, err), "SRSS 组合 ΣV_j² 平方根",
        "got " + std::to_string(r.baseZ) + " want " + std::to_string(srssWant));
  // 与静态 srss() 一致
  check(std::abs(r.baseZ - ResponseSpectrum::srss(r.vZ)) <
            1e-9 * r.baseZ,
        "solve 与静态 srss() 一致");
  // 参与质量守恒：Σm* = 方向总质量（先验证手算常量 Mz）
  const double mCol = mat.gamma * sec.A * L / kG;
  const double Mz = mCol + 0.5 * mCol;
  check(std::abs(mr.massZ - Mz) < 1e-6 * Mz, "Σm* = 方向总质量",
        "got " + std::to_string(mr.massZ) + " want " + std::to_string(Mz));
  check(mr.ratioZ > 0.999, "Z 向参与质量比 = 1");
  check(r.massOk, "两模态全取时参与质量检查通过");
}

// -----------------------------------------------------------------------------
//  D. CQC 模态耦联
// -----------------------------------------------------------------------------
static void test_cqc() {
  std::printf("\n== D. CQC 模态耦联 ==\n");
  double err = 0.0;

  // ρ_jk 特殊值：β = ω_k/ω_j = 1 → ρ = 1（重频完全耦联）
  const double rhoSame = ResponseSpectrum::cqcRho(0.05, 2.0, 2.0);
  check(std::abs(rhoSame - 1.0) < 1e-12, "β=1（重频）→ ρ = 1",
        "got " + std::to_string(rhoSame));

  // 远频：β → ∞，ρ → 0
  const double rhoFar = ResponseSpectrum::cqcRho(0.05, 1.0, 100.0);
  check(rhoFar < 1e-3 && rhoFar > 0.0, "远频 β=100 → ρ 接近 0",
        "got " + std::to_string(rhoFar));

  // 规范式闭合形式：β = ω2/ω1 = 3（两自由度链频率比 ≈ 2.62 附近）
  //   手动用通用式重算（避免与实现同源 bug 互证）：
  //   ρ = 8ζ²(1+β)β^1.5 / ((1−β²)² + 4ζ²β(1+β)²)
  const double zeta = 0.05, w1 = 10.0, w2 = 25.0, b = w2 / w1;
  const double z2 = zeta * zeta;
  const double den = (1 - b * b) * (1 - b * b) +
                     4.0 * z2 * b * (1 + b) * (1 + b);
  const double rhoWant = 8.0 * z2 * (1 + b) * std::pow(b, 1.5) / den;
  check(relErr(ResponseSpectrum::cqcRho(zeta, w1, w2), rhoWant, 1e-12, err),
        "ρ_jk 对标规范式（β=2.5）",
        "got " + std::to_string(ResponseSpectrum::cqcRho(zeta, w1, w2)));

  // CQC 组合 = √(V1²+V2²+2ρV1V2)
  const double v1 = 3.0, v2 = 1.2;
  const double want = std::sqrt(v1 * v1 + v2 * v2 +
                                2.0 * rhoWant * v1 * v2);
  // cqc(zeta, omega=[10,25], value=[3,1.2])，ρ12 = cqcRho(ζ,10,25)
  const std::vector<double> om{10.0, 25.0};
  const std::vector<double> vv{3.0, 1.2};
  const double got = ResponseSpectrum::cqc(zeta, om, vv);
  check(relErr(got, want, 1e-12, err), "cqc() = √(ΣΣ ρ_jk V_j V_k)",
        "got " + std::to_string(got) + " want " + std::to_string(want));

  // 单模态 CQC = SRSS 退化
  const std::vector<double> om1{7.0};
  const std::vector<double> vv1{5.0};
  check(std::abs(ResponseSpectrum::cqc(zeta, om1, vv1) - 5.0) < 1e-12,
        "单模态 CQC 退化 = |V|");

  // 与 solve(method=1) 对标：复用两柱链模型
  Model m;
  const double L = 2.0, A = 0.04, E = 3.0e7;
  SectionProperties sec;
  sec.A = A;
  sec.Iy = sec.Iz = 1.3333e-3;
  sec.J = 0.0;
  const Material mat = makeConcrete();
  const Id n0 = m.addNode({0, 0, 0}, 0);
  const Id n1 = m.addNode({0, 0, L}, 1);
  const Id n2 = m.addNode({0, 0, 2 * L}, 2);
  BeamElement* e1 = m.addBeam(n0, n1, sec, mat);
  e1->setUpHint({1, 0, 0});
  BeamElement* e2 = m.addBeam(n1, n2, sec, mat);
  e2->setUpHint({1, 0, 0});
  m.fixAll(n0);
  m.fixNode(n1, true, true, false, true, true, true);
  m.fixNode(n2, true, true, false, true, true, true);
  m.assignDofs();

  ModalAnalysis ma(m);
  ModalResult mr = ma.solve(ModalAnalysis::Options{2, 1e-10, 300, 0});
  if (!mr.ok) return;
  DesignSpectrum spec;
  ResponseSpectrum rs(m);
  ResponseSpectrum::Options opt;
  opt.method = 1;   // CQC
  opt.verticalScale = 1.0;
  const SpectrumResult r = rs.solve(mr, spec, opt);
  if (!r.ok) return;
  // 静态 cqc()：ω 数组已是升序（ω1 < ω2）
  const double wantCqc =
      ResponseSpectrum::cqc(spec.zeta, mr.omega, r.vZ);
  check(std::abs(r.baseZ - wantCqc) < 1e-9 * wantCqc,
        "solve(CQC) 与静态 cqc() 一致",
        "got " + std::to_string(r.baseZ) + " want " + std::to_string(wantCqc));
  // CQC 不小于 SRSS（近频时 CQC ≥ SRSS）
  const double srssV = ResponseSpectrum::srss(r.vZ);
  check(r.baseZ >= srssV - 1e-9 * srssV, "CQC ≥ SRSS（近频耦联）",
        "cqc " + std::to_string(r.baseZ) + " srss " + std::to_string(srssV));
}

// -----------------------------------------------------------------------------
//  E. 参与质量 ≥ 90% 检查
// -----------------------------------------------------------------------------
static void test_mass_check() {
  std::printf("\n== E. 参与质量 ≥ 90% 检查 ==\n");
  double err = 0.0;

  // 刚性楼板单层两柱：X/Y 单自由度（ratio=1），Z 向两支轴向模态。
  // 只取 1 阶（X 平动）→ Z 向参与质量 ≈ 0 → 检查必须失败；
  // 取全模态 → 通过。
  Model m;
  const double L = 3.0, A = 0.04;
  SectionProperties sec;
  sec.A = A;
  sec.Iy = sec.Iz = 1.3333e-3;
  sec.J = 0.0;
  const Material mat = makeConcrete();
  const Id b0 = m.addNode({0, 0, 0}, 0);
  const Id b1 = m.addNode({0, 3, 0}, 0);
  const Id t0 = m.addNode({0, 0, 3}, 1);
  const Id t1 = m.addNode({0, 3, 3}, 1);
  BeamElement* c1 = m.addBeam(b0, t0, sec, mat);
  c1->setUpHint({1, 0, 0});
  BeamElement* c2 = m.addBeam(b1, t1, sec, mat);
  c2->setUpHint({1, 0, 0});
  m.fixAll(b0);
  m.fixAll(b1);
  m.addRigidDiaphragm({t0, t1}, 1, true);
  m.assignDofs();

  ModalAnalysis ma(m);
  ModalResult mr1 = ma.solve(ModalAnalysis::Options{1, 1e-10, 300, 0});
  check(mr1.ok, mr1.ok ? "模态（1 阶）求解成功" : ("失败：" + mr1.message).c_str());
  if (!mr1.ok) return;
  // 1 阶模态只取到 UX/UY/RZ 平动-扭转混合模态（UX 与 RZ 严格简并，
  // 特征向量在简并子空间内旋转，X 向参与质量被分摊）——
  // Z 向轴向模态完全没有参与，因此必然存在参与质量不足的方向。
  check(mr1.totalMassZ > 1e-9, "Z 向存在可动质量（轴向）");
  check(mr1.ratioX < 0.8999 || mr1.ratioZ < 0.8999,
        "1 阶模态参与质量不足（X 分摊或 Z 缺失）",
        "ratioX " + std::to_string(mr1.ratioX) +
            " ratioZ " + std::to_string(mr1.ratioZ));

  DesignSpectrum spec;
  ResponseSpectrum rs(m);
  ResponseSpectrum::Options opt;
  opt.method = 0;
  opt.checkMass = true;
  opt.verticalScale = 1.0;
  const SpectrumResult rb = rs.solve(mr1, spec, opt);
  check(rb.ok, rb.ok ? "反应谱求解（1 阶模态）成功" : ("失败：" + rb.message).c_str());
  check(!rb.massOk, "模态数不足：参与质量检查应失败",
        rb.message.c_str());
  check(rb.massOkZ == false, "Z 向轴向质量未参与 → massOkZ = false");

  // checkMass = false：不阻断结果
  opt.checkMass = false;
  const SpectrumResult rb2 = rs.solve(mr1, spec, opt);
  check(rb2.massOk, "checkMass=false 时检查不阻断");

  // 取全部模态：参与质量齐全 → 通过
  ModalResult mrAll = ma.solve(ModalAnalysis::Options{6, 1e-10, 300, 0});
  check(mrAll.ok && mrAll.ratioZ > 0.999, "全模态 Z 向参与质量比 = 1",
        "got " + std::to_string(mrAll.ratioZ));
  opt.checkMass = true;
  const SpectrumResult rg = rs.solve(mrAll, spec, opt);
  check(rg.ok && rg.massOk, "全模态取齐后参与质量检查通过",
        rg.message.c_str());

  // 单轴向模型（test B 已覆盖）：X/Y 无可动质量豁免 —— 在这里再用
  // 两柱链模型验证"方向无可动质量时 massOkX/Y = true"
  Model m2;
  const double L2 = 2.0;
  SectionProperties sec2;
  sec2.A = A;
  sec2.Iy = sec2.Iz = 1.3333e-3;
  sec2.J = 0.0;
  const Material mat2 = makeConcrete();
  const Id a0 = m2.addNode({0, 0, 0}, 0);
  const Id a1 = m2.addNode({0, 0, L2}, 1);
  const Id a2 = m2.addNode({0, 0, 2 * L2}, 2);
  BeamElement* p1 = m2.addBeam(a0, a1, sec2, mat2);
  p1->setUpHint({1, 0, 0});
  BeamElement* p2 = m2.addBeam(a1, a2, sec2, mat2);
  p2->setUpHint({1, 0, 0});
  m2.fixAll(a0);
  m2.fixNode(a1, true, true, false, true, true, true);
  m2.fixNode(a2, true, true, false, true, true, true);
  m2.assignDofs();
  ModalAnalysis ma2(m2);
  ModalResult mr2 = ma2.solve(ModalAnalysis::Options{2, 1e-10, 300, 0});
  if (!mr2.ok) return;
  ResponseSpectrum rs2(m2);
  ResponseSpectrum::Options opt2;
  opt2.method = 0;
  opt2.checkMass = true;
  opt2.verticalScale = 1.0;
  const SpectrumResult rz = rs2.solve(mr2, spec, opt2);
  check(rz.massOkX && rz.massOkY, "轴向模型 X/Y 无可动质量 → massOkX/Y 豁免",
        rz.message.c_str());
  check(rz.massOk, "轴向模型全向通过参与质量检查");
  (void)err;
}

int main() {
  std::printf("=== test_spectrum: 反应谱与振型分解反应谱法 ===\n");
  test_spectrum_curve();
  test_single_dof();
  test_two_dof_srss();
  test_cqc();
  test_mass_check();
  std::printf("\n通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}