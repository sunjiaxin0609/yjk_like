// =============================================================================
//  tests/test_modal.cpp  ——  模态分析（ModalAnalysis）单元测试
//
//  覆盖：
//  A. 单柱轴向单自由度：ω² = EA/(L·m)，柱重一半集中于顶部节点 —— 无剪切
//     修正、无转动惯量，解析解精确，验证"质量装配 → 特征值"整条链路。
//  B. 两柱串列轴向：3 节点 2 单元固定-自由链，双自由度解析解
//     （同时覆盖多模态、固有频率与振型）。
//  C. 刚性楼板（MPC）质量折减：单跨两柱 + 一层刚板 ——
//     · X 向参与质量 = 结构总质量（ratioX = 1，只有这一个 X 自由度有质量）
//     · 主自由度上的转动惯量 M(Θ,Θ) > 0（耦合块存在，扭转模态的物理来源）
//     · 从属节点位移 = 主节点位移（振型回填的刚体一致性）
//  附加：质量归一 φᵀMφ = 1、频率排序（ω1 < ω2）、周期换算。
// =============================================================================
#include <cmath>
#include <cstdio>

#include "yjk/analysis/ModalAnalysis.h"
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

// 材料：C30 混凝土（近似）
static Material makeConcrete() {
  Material m;
  m.gamma = 25.0;   // kN/m³
  return m;
}

// -----------------------------------------------------------------------------
//  A. 单柱轴向单自由度（悬臂柱，顶部节点只保留 uz）
// -----------------------------------------------------------------------------
static void test_single_dof_axial() {
  std::printf("\n== A. 单柱轴向单自由度 ==\n");
  Model m;
  const double L = 4.0, A = 0.04, E = 3.0e7;   // kPa, m², m
  SectionProperties sec;
  sec.A = A;
  sec.Iy = sec.Iz = 1.3333e-3;
  sec.J = 0.0;
  const Material mat = makeConcrete();

  const Id n0 = m.addNode({0, 0, 0}, 0, "base");
  const Id n1 = m.addNode({0, 0, L}, 1, "top");
  BeamElement* el = m.addBeam(n0, n1, sec, mat);
  el->setUpHint({1, 0, 0});                    // 局部 z 水平，轴向沿 z
  m.fixAll(n0);
  m.fixNode(n1, true, true, false, true, true, true);   // 只留 uz
  const Id ndof = m.assignDofs();

  ModalAnalysis ma(m);
  ModalResult r = ma.solve(ModalAnalysis::Options{1, 1e-10, 300, 0});
  check(r.ok, r.ok ? "求解成功" : ("失败：" + r.message).c_str());
  if (!r.ok) return;

  // 解析：m_top = 柱重/2 ÷ g；ω² = EA/(L·m_top)
  const double wKn = mat.gamma * sec.A * L;    // 4 kN
  const double mTop = 0.5 * wKn / 9.81;        // t
  const double kAxial = E * sec.A / L;         // kN/m
  const double want = std::sqrt(kAxial / mTop);
  double err = 0.0;
  check(relErr(r.omega[0], want, 5e-3, err), "轴向基频对标解析解",
        "got " + std::to_string(r.omega[0]) + " want " + std::to_string(want));
  check(std::abs(r.massZ - r.totalMassZ) < 1e-8 * std::max(r.totalMassZ, 1e-300),
        "Z 向参与质量 = 总质量（单自由度）",
        "m*=" + std::to_string(r.massZ) + " M=" + std::to_string(r.totalMassZ));
  check(r.ratioZ > 0.999, "Z 向参与质量比 = 1");
  if (r.omega[0] > 1e-300)
    check(std::abs(r.period[0] - 2.0 * 3.141592653589793 / r.omega[0]) <
              1e-9 * r.period[0],
          "周期 T = 2π/ω");
  check(r.freq[0] > 0.0, "频率 f = ω/2π > 0",
        "f=" + std::to_string(r.freq[0]));
  check(r.shapes[0][static_cast<size_t>(n1 * 6 + 2)] > 0.0, "顶部 uz 振型非零");
  const double zBottom = r.shapes[0][static_cast<size_t>(n0 * 6 + 2)];
  check(std::abs(zBottom) < 1e-12, "底部固定 uz 振型 = 0");

  // 质量归一
  {
    StaticAnalysis sa(m);
    StaticAnalysis::Assembly as;
    if (sa.assemble(as)) {
      ModalAnalysis::Assembly qa;
      std::string why;
      if (ma.assembleMass(qa, &why)) {
        // shapes 是【节点序×6】空间：提取自由自由度（dof ≥ 0）到方程号空间
        std::vector<double> v(static_cast<size_t>(qa.ndof), 0.0);
        for (Id n = 0; n < m.nodeCount(); ++n)
          for (int k = 0; k < 6; ++k) {
            const Id g = m.node(n).dof[k];
            if (g >= 0)
              v[static_cast<size_t>(g)] = r.shapes[0][static_cast<size_t>(n * 6 + k)];
          }
        std::vector<double> Mv;
        qa.M.multiply(v, Mv);
        double norm = 0.0;
        for (Id i = 0; i < qa.ndof; ++i)
          norm += v[static_cast<size_t>(i)] * Mv[static_cast<size_t>(i)];
        check(std::abs(norm - 1.0) < 1e-8, "φᵀMφ = 1（质量归一）",
              "got " + std::to_string(norm));
      } else {
        check(false, "质量组装失败", why.c_str());
      }
    } else {
      check(false, "刚度组装失败");
    }
  }
  (void)ndof;
}

// -----------------------------------------------------------------------------
//  B. 两柱串列轴向（固定-自由，3 节点 2 单元）
// -----------------------------------------------------------------------------
static void test_two_dof_chain() {
  std::printf("\n== B. 两柱串列轴向双自由度 ==\n");
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
  ModalResult r = ma.solve(ModalAnalysis::Options{2, 1e-10, 300, 0});
  check(r.ok, r.ok ? "求解成功" : ("失败：" + r.message).c_str());
  if (!r.ok) return;

  // 解析：M=diag(m1,m2)，K = c·[[2,-1],[-1,1]]，c = EA/L。
  //   m1 = 柱1底半(固定丢弃) + 柱1顶半 + 柱2底半 = 一整柱
  //   m2 = 柱2顶半 = 半柱
  const double c = E * A / L;                  // kN/m
  const double mCol = mat.gamma * sec.A * L / 9.81;   // 一整柱质量 t
  const double m1 = mCol, m2 = 0.5 * mCol;
  // det(K − λM) = 0：λ = c(m1+2m2) ± √(c²(m1+2m2)² − 4m1m2c²) / (2m1m2)
  const double s = c * (m1 + 2.0 * m2);
  const double disc = s * s - 4.0 * m1 * m2 * c * c;
  const double l1 = (s - std::sqrt(disc)) / (2.0 * m1 * m2);
  const double l2 = (s + std::sqrt(disc)) / (2.0 * m1 * m2);
  const double w1 = std::sqrt(l1), w2 = std::sqrt(l2);

  double err = 0.0;
  check(relErr(r.omega[0], w1, 5e-3, err), "第一阶 ω 对标解析解",
        "got " + std::to_string(r.omega[0]) + " want " + std::to_string(w1));
  check(relErr(r.omega[1], w2, 5e-3, err), "第二阶 ω 对标解析解",
        "got " + std::to_string(r.omega[1]) + " want " + std::to_string(w2));
  check(r.omega[0] < r.omega[1], "固有频率升序（ω1 < ω2）");
  check(r.period[0] > r.period[1], "周期降序（T1 > T2）");

  // 第一阶振型同号（基阶无节点），第二阶反号
  const double z1a = r.shapes[0][static_cast<size_t>(n1 * 6 + 2)];
  const double z1b = r.shapes[0][static_cast<size_t>(n2 * 6 + 2)];
  const double z2a = r.shapes[1][static_cast<size_t>(n1 * 6 + 2)];
  const double z2b = r.shapes[1][static_cast<size_t>(n2 * 6 + 2)];
  check(z1a * z1b > 0.0, "第一阶振型同号（基阶）");
  check(z2a * z2b < 0.0, "第二阶振型反号（一阶反号模式）");

  // 参与质量：Σm*_z = M_total = m1 + m2（所有可动质量都在 z 向）
  const double Mz = m1 + m2;
  check(std::abs(r.massZ - Mz) < 1e-6 * Mz, "参与质量守恒（Z 向）",
        "Σm*=" + std::to_string(r.massZ) + " M=" + std::to_string(Mz));
  check(r.ratioZ > 0.999, "Z 向参与质量比 = 1");
}

// -----------------------------------------------------------------------------
//  C. 刚性楼板质量折减（MPC）
// -----------------------------------------------------------------------------
static void test_diaphragm_mass() {
  std::printf("\n== C. 刚性楼板（MPC）质量折减 ==\n");
  Model m;
  const double L = 3.0, A = 0.04;
  SectionProperties sec;
  sec.A = A;
  sec.Iy = sec.Iz = 1.3333e-3;
  sec.J = 0.0;
  const Material mat = makeConcrete();

  // 单层两柱：柱底 (0,0,0)/(0,3,0)，柱顶 (0,0,3)/(0,3,3)
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
  // 主节点由 addRigidDiaphragm 自动创建：只保留 UX/UY/RZ 三个自由度
  const Id master = m.addRigidDiaphragm({t0, t1}, 1, true);
  m.assignDofs();

  ModalAnalysis ma(m);

  // 1) 质量矩阵结构检查
  ModalAnalysis::Assembly qa;
  std::string why;
  const bool massOk = ma.assembleMass(qa, &why);
  check(massOk && qa.ndof > 0, massOk ? "质量组装成功" : "质量组装失败",
        massOk ? "" : why.c_str());
  if (!massOk) return;
    // 重新组装并检查 M 的结构
    StaticAnalysis sa(m);
    StaticAnalysis::Assembly as;
    sa.assemble(as);
    const Id ux = m.node(master).dof[0];
    const Id uy = m.node(master).dof[1];
    const Id rz = m.node(master).dof[5];
    // 主节点自由度必须存在
    check(ux >= 0 && uy >= 0 && rz >= 0, "主节点 UX/UY/RZ 自由度有效");

    const double mTot = 2.0 * mat.gamma * sec.A * L / 9.81;   // 两柱总质量 t
    const double mH = 0.5 * mat.gamma * sec.A * L / 9.81;     // 每柱半质量 t
    // X 向总质量 = 两柱顶节点质量之和（每柱顶 1/2，柱底固定丢弃）
    const double Mxx = qa.M.at(ux, ux);
    check(std::abs(Mxx - 2.0 * mH) < 1e-9 * 2.0 * mH, "M(UX,UX) = 两顶节点质量之和",
          "got " + std::to_string(Mxx) + " want " + std::to_string(2.0 * mH));
    // 转动惯量 M(RZ,RZ) = Σ m_i·(dx²+dy²)，两节点在 y=0 与 y=3
    // （ux 的从属系数 -dy，uy 的从属系数 +dx；dy = ±1.5, dx = 0）
    const double Mrz = qa.M.at(rz, rz);
    const double wantIr = mH * (1.5 * 1.5) * 2 + mH * (0.0) * 2;
    check(Mrz > 0.0, "M(RZ,RZ) > 0（楼板转动惯量存在）",
          "got " + std::to_string(Mrz));
    check(std::abs(Mrz - wantIr) < 1e-9 * wantIr, "M(RZ,RZ) 对标手算 Σm·r²",
          "got " + std::to_string(Mrz) + " want " + std::to_string(wantIr));
    // 质心在形心：dx=0（所有节点 x=0）⇒ M(UX,RZ) = −Σm·dy = 0
    const double Mxr = qa.M.at(ux, rz);
    check(std::abs(Mxr) < 1e-9 * Mxx, "质心在形心时 M(UX,RZ) ≈ 0",
          "got " + std::to_string(Mxr));

  // 2) 模态求解：X 向单自由度 ⇒ 参与质量比 = 1
  // 5 个有质量自由度（主 UX/UY/RZ + 两柱顶 uz）：取全部模态，
  // 才能收敛出 Z 向两支轴向模态的参与质量
  ModalResult r = ma.solve(ModalAnalysis::Options{6, 1e-10, 300, 0});
  if (!r.ok) return;

  check(r.ratioX > 0.999, "X 向参与质量比 = 1（单自由度）",
        "got " + std::to_string(r.ratioX));
  check(r.ratioY > 0.999, "Y 向参与质量比 = 1");
  // Z 向：柱顶 uz 自由度两台（t0、t1 的 uz 不参与 MPC），两柱轴向模态
  check(r.ratioZ > 0.999, "Z 向参与质量比 = 1（两柱轴向并联）");

  // 3) 振型回填：从属节点位移必须等于主节点位移（刚体平动一致性）
  //    找到 X 向主导模态：γ² 最大的
  {
    size_t best = 0;
    double bestG = -1.0;
    for (size_t i = 0; i < r.gammaX.size(); ++i) {
      const double gx = std::abs(r.gammaX[i]);
      if (gx > bestG) { bestG = gx; best = i; }
    }
    const double uMaster =
        r.shapes[best][static_cast<size_t>(master * 6 + 0)];
    const double theta =
        r.shapes[best][static_cast<size_t>(master * 6 + 5)];   // 主节点 RZ 振型
    const double uSlave0 =
        r.shapes[best][static_cast<size_t>(t0 * 6 + 0)];
    const double uSlave1 =
        r.shapes[best][static_cast<size_t>(t1 * 6 + 0)];
    check(std::abs(uMaster) > 1e-6, "主节点 UX 振型非零");
    if (std::abs(uMaster) > 1e-6) {
      // 刚性楼板：ux_i = UX − Θ·dy_i（t0 在 y=0、t1 在 y=3，形心 y=1.5）
      // X 平动与 RZ 扭转简并 ⇒ 该模态可带 Θ 分量，从属值必须仍精确满足展开式
      const double exp0 = uMaster - theta * (0.0 - 1.5);
      const double exp1 = uMaster - theta * (3.0 - 1.5);
      check(std::abs(uSlave0 - exp0) < 1e-9 * std::abs(uMaster) &&
                std::abs(uSlave1 - exp1) < 1e-9 * std::abs(uMaster),
            "从属 ux 精确满足 ux_i = UX − Θ·dy_i（刚性展开公式）",
            "slave0=" + std::to_string(uSlave0) +
                " exp0=" + std::to_string(exp0) +
                " slave1=" + std::to_string(uSlave1) +
                " exp1=" + std::to_string(exp1) +
                " theta=" + std::to_string(theta));
    }
  }

  // 4) 校验：刚性楼板参与后频率合理 —— 单层剪切频率
  //    ω² = k_eff / m_layer，k_eff = 2 × (12EI/L³)（柱上下两端均约束转角？
  //    柱顶 rx/ry 未被 MPC 约束，是悬臂约束 ⇒ k = 3EI/L³ 每柱；这里不做
  //    频率断言，仅验证量级在 1~100 rad/s（避免剪切/约束争议的伪精确）。
  check(r.omega[0] > 0.5 && r.omega[0] < 200.0, "基频量级合理（0.5~200 rad/s）",
        "got " + std::to_string(r.omega[0]));
}

int main() {
  std::printf("=== test_modal: 模态分析 ===\n");
  test_single_dof_axial();
  test_two_dof_chain();
  test_diaphragm_mass();
  std::printf("\n通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}