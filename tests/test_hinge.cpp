// =============================================================================
//  tests/test_hinge.cpp  ——  端部释放（铰）验证
//
//  验证策略：对标结构力学解析解，覆盖三类工况——
//   1. 单元级：凝聚刚度的对称性、平动刚体零内力、凝聚保持刚度不变性；
//   2. 整体级（简支铰接梁）：铰端弯矩 = 0、跨中弯矩 = PL/4、
//      跨中挠度 = PL³/48EI、端部剪力 = ±P/2、求解残差 ~ 机器精度；
//   3. 倾斜/旋转用例：45° 斜梁铰端【局部】My/Mz = 0，
//      跨中弯矩 = P·L·cos45°/4（验证 Tᵀ Kc T 变换在非轴对齐下成立）。
//
//  端部释放的语义：releaseAt(node, comp) 释放局部坐标系下该端自由度，
//  对应自由度从单元刚度做静力凝聚（Guyan），被释放自由度处内力恒为 0。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "yjk/element/BeamElement3D.h"
#include "yjk/model/Model.h"
#include "yjk/analysis/StaticAnalysis.h"

using namespace yjk;

static int g_pass = 0, g_fail = 0;
static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}
static void checkNear(double got, double want, double relTol, const char* what) {
  const double rel = std::abs(want) > 1e-30 ? std::abs(got - want) / std::abs(want)
                                            : std::abs(got - want);
  if (rel <= relTol) {
    ++g_pass;
    std::printf("  [ OK ] %s  (got %.8g, want %.8g, rel %.2e)\n", what, got, want, rel);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %.8g, want %.8g, rel %.2e)\n", what, got, want, rel);
  }
}

// ===========================================================================
//  1. 单元级：凝聚刚度
// ===========================================================================
static void testCondensation() {
  std::printf("\n== 1. 凝聚刚度 ==\n");
  const Material mat = Material::concreteC(30);
  const SectionProperties sec = section::rect(0.3, 0.5);

  // (a) 对称性：Kc 必须对称（凝聚是 K 的 Schur 补，对称性必须保持）
  {
    BeamElement3D be({0, 0, 0}, {6, 0, 0}, sec, mat, Vec3{0, 0, -1});
    be.releaseAt(0, 5);            // i 端释放 rz
    be.releaseAt(1, 4);            // j 端释放 ry
    check(be.isValid(), "带释放的梁单元有效");
    check(be.condensationOk(), "凝聚成功（Kbb 非奇异）");
    const Mat12 Kc = be.condensedStiffness();
    double asym = 0.0;
    for (int i = 0; i < 12; ++i)
      for (int j = 0; j < 12; ++j) asym = std::max(asym, std::abs(Kc(i, j) - Kc(j, i)));
    checkNear(asym / Kc.maxAbs(), 0.0, 1e-12, "凝聚刚度对称（Kc = Kcᵀ）");
  }

  // (b) 平动刚体模式：凝聚刚度对整体平动必须零内力。
  //     平动模式下两端转角均为 0，释放自由度（转角）上无运动分量，
  //     因此凝聚后刚度对纯平动应严格无内力 —— 这是单元级最直接的完整性判据。
  {
    auto checkRigidTranslate = [&](const Vec3& pi, const Vec3& pj, const Vec3& up,
                                   const char* title) {
      BeamElement3D be(pi, pj, sec, mat, up);
      be.releaseAt(0, 4); be.releaseAt(0, 5);   // i 端铰
      be.releaseAt(1, 4); be.releaseAt(1, 5);   // j 端铰
      if (!be.isValid() || !be.condensationOk()) { check(false, title); return; }
      const Mat12 K = be.stiffness();           // 内部自动取凝聚刚度
      const double kmax = K.maxAbs();
      double worst = 0.0;
      for (int k = 0; k < 3; ++k) {             // 3 个平动刚体模式
        Vec12d v{};
        for (int n = 0; n < 2; ++n) v[n * 6 + k] = 1.0;
        const Vec12d f = K * v;
        for (double x : f) worst = std::max(worst, std::abs(x));
      }
      checkNear(worst / kmax, 0.0, 1e-10, title);
    };
    checkRigidTranslate({0, 0, 0}, {4, 0, 0}, {0, 0, -1}, "铰接梁平动刚体零内力（X 向）");
    checkRigidTranslate({0, 0, 0}, {0, 4, 0}, {0, 0, -1}, "铰接梁平动刚体零内力（Y 向）");
    checkRigidTranslate({1, 2, 3}, {4, 6, 5}, {0, 1, 0}, "铰接梁平动刚体零内力（任意空间杆）");
  }

  // (c) 凝聚保持 a 子空间刚度不变性：
  //     对"没有 b 分量运动"的位移 v（v_b = 0），凝聚刚度与原刚度在 a 子空间
  //     必须给出相同内力：Kc·v|_a = Kaa·v_a（凝聚是精确缩并，非近似）。
  {
    BeamElement3D be({0, 0, 0}, {6, 0, 0}, sec, mat, Vec3{0, 0, -1});
    be.releaseAt(1, 5);                       // 只释放 j 端 rz
    const Mat12 K0 = be.localStiffness();     // 原始 12×12
    const Mat12 Kc = be.condensedStiffness();
    // 取一个不触及释放自由度的位移：j 端绕 y 转角（a 自由度）
    Vec12d v{};
    v[6 + 4] = 1.0;                           // j 端 ry = 1
    const Vec12d fK = K0 * v;                 // 原刚度内力（b 分量 = 释放自由度 j-rz）
    const Vec12d fKc = Kc * v;
    double diffA = 0.0, diffB = 0.0, maxF = 0.0;
    for (int i = 0; i < 12; ++i) {
      maxF = std::max(maxF, std::abs(fK[static_cast<size_t>(i)]));
      if (i == 11) diffB = std::max(diffB, std::abs(fK[static_cast<size_t>(i)] - fKc[static_cast<size_t>(i)]));  // 释放自由度 11 = j·rz
      else diffA = std::max(diffA, std::abs(fK[static_cast<size_t>(i)] - fKc[static_cast<size_t>(i)]));
    }
    checkNear(diffA / std::max(maxF, 1e-30), 0.0, 1e-10,
              "凝聚刚度在保留自由度上与原刚度一致");
  }
}

// =============================================================================
//  2. 简支铰接梁（整体模型）：两端释放 ry/rz + 跨中集中力
//
//  解析解（含剪切变形，Timoshenko）：
//    · 铰端弯矩 = 0
//    · 跨中弯矩 = PL/4
//    · 端部剪力 = ±P/2
//    · 跨中挠度 = PL³/48EI·(1+Φ)，Φ = 12EI/(G·As·L²)
//    （竖向荷载沿局部 y'，绕 z' 弯曲 → I = Iz、As = Asz）
// ===========================================================================
static void testSimplySupportedHinged() {
  std::printf("\n== 2. 简支铰接梁（两端释放 ry/rz）==\n");
  const double L = 6.0, P = 100.0;            // kN, kN·m
  const auto sec = section::rect(0.3, 0.5);   // b=0.3, h=0.5 → Iz = b·h³/12
  const auto mat = Material::concreteC(30);
  const double Iz = sec.Iz;
  checkNear(Iz, 0.3 * 0.5 * 0.5 * 0.5 / 12.0, 1e-12, "截面 Iz = b·h³/12");

  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id nc = m.addNode({L / 2, 0, 0}, 1);
  const Id n1 = m.addNode({L, 0, 0}, 1);
  BeamElement* b1 = m.addBeam(n0, nc, sec, mat, Vec3{0, 0, -1});
  BeamElement* b2 = m.addBeam(nc, n1, sec, mat, Vec3{0, 0, -1});
  // 两端铰：释放绕局部 y/z 的转角（强/弱轴两向都放，端弯矩恒 0）
  b1->releaseAt(0, 4); b1->releaseAt(0, 5);
  b2->releaseAt(1, 4); b2->releaseAt(1, 5);

  // 简支：n0 三平动固定 + rx 固定（防整杆扭转刚体）；n1 滚轴同 rx 固定
  m.fixNode(n0, true, true, true, true, false, true);
  m.fixNode(n1, true, true, true, false, false, true);
  m.addNodeForce(nc, {0, 0, -P});             // 跨中 -Z 集中力

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
  checkNear(r.residual, 0.0, 1e-12, "平衡残差 ~ 机器精度");

  // 跨中挠度（含剪切修正）：w = PL³/48EIz·(1+Φz)，Φz = 12EIz/(G·Asz·L²)
  const double G = mat.shearModulus();
  const double phi = 12.0 * mat.E * Iz / (G * sec.Asz * L * L);
  const double w = -r.u[static_cast<size_t>(nc) * 6 + 2];   // 取幅值
  checkNear(w, P * L * L * L / (48.0 * mat.E * Iz) * (1.0 + phi), 1e-9,
            "跨中挠度 = PL³/48EI·(1+Φ)（含剪切）");

  // 铰端弯矩必须为 0：elem0 i 端（n0 铰）、elem1 j 端（n1 铰）
  const double mI0 = r.beamForces[0].ef.Mz;
  const double mJ1 = r.beamForces[1].ef.Mzj;
  checkNear(mI0, 0.0, 1e-9, "铰端 Mz(n0) = 0");
  checkNear(mJ1, 0.0, 1e-9, "铰端 Mz(n1) = 0");

  // 跨中弯矩 = PL/4（elem1 的 i 端 = 跨中）
  const double mMid = r.beamForces[1].ef.Mz;
  checkNear(mMid, P * L / 4.0, 1e-9, "跨中弯矩 = PL/4");

  // 端部剪力 = ±P/2
  checkNear(r.beamForces[0].ef.Vy, P / 2.0, 1e-9, "梁 1 端部剪力 = +P/2");
  checkNear(r.beamForces[1].ef.Vyj, -P / 2.0, 1e-9, "梁 2 端部剪力 = -P/2");
}

// ===========================================================================
//  3. 悬臂释放：i 端全固支 + j 端释放转角（与无释放悬臂一致）
//     —— 验证"释放一个本来就自由的自由度"不破坏其余刚度；
//        同时验证零刚度自由度按分量方向判荷载自动锁定
//        （j 端 uz 有荷载、ry/rz 无荷载 → 释放转角被锁定，不误报"无刚度"）
//  解析解：w = PL³/3EI·(1+Φ/4)，Φ = 12EI/(G·As·L²)（悬臂端集中力含剪切）
// ===========================================================================
static void testCantileverRelease() {
  std::printf("\n== 3. 悬臂梁（j 端释放 ry/rz）==\n");
  const double L = 4.0, P = 50.0;
  const auto sec = section::rect(0.3, 0.5);
  const auto mat = Material::concreteC(30);
  const double Iz = sec.Iz;

  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id n1 = m.addNode({L, 0, 0}, 1);
  BeamElement* b = m.addBeam(n0, n1, sec, mat, Vec3{0, 0, -1});
  // 释放自由端转角：悬臂自由端本来弯矩为 0，释放不应改变任何结果
  b->releaseAt(1, 4); b->releaseAt(1, 5);
  // i 端全固（含 rx，防扭转刚体）
  m.fixNode(n0, true, true, true, true, true, true);
  m.addNodeForce(n1, {0, 0, -P});             // 自由端 -Z 集中力

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功（同节点 uz 荷载不误锁释放转角）");
  if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
  checkNear(r.residual, 0.0, 1e-12, "平衡残差 ~ 机器精度");

  // 悬臂端挠度（含剪切修正）：w = PL³/3EI·(1+Φ/4)
  const double G = mat.shearModulus();
  const double phi = 12.0 * mat.E * Iz / (G * sec.Asz * L * L);
  const double w = -r.u[static_cast<size_t>(n1) * 6 + 2];
  checkNear(w, P * L * L * L / (3.0 * mat.E * Iz) * (1.0 + phi / 4.0), 1e-9,
            "自由端挠度 = PL³/3EI·(1+Φ/4)（含剪切）");
  const double mM0 = r.beamForces[0].ef.Mz;   // i 端弯矩
  checkNear(mM0, -P * L, 1e-9, "固支端弯矩 = -PL（符号约定：拉上压下的受弯杆件）");
}

// ===========================================================================
//  4. 倾斜（45°，非轴对齐）铰接梁：铰端局部 My/Mz = 0、跨中弯矩解析值
//     —— 验证 Tᵀ Kc T 变换路径在旋转坐标系下成立
// ===========================================================================
static void testTiltedHinged() {
  std::printf("\n== 4. 45° 倾斜铰接梁 ==\n");
  const double L = 6.0, P = 100.0;
  const double kPi = std::atan2(0.0, -1.0);
  const double c = std::cos(kPi / 4.0);
  const auto sec = section::rect(0.3, 0.5);
  const auto mat = Material::concreteC(30);

  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 1);
  const Id nc = m.addNode({L * c / 2, 0, -L * c / 2}, 1);    // 跨中（XZ 平面，45°）
  const Id n1 = m.addNode({L * c, 0, -L * c}, 1);
  BeamElement* b1 = m.addBeam(n0, nc, sec, mat, Vec3{0, 1, 0});
  BeamElement* b2 = m.addBeam(nc, n1, sec, mat, Vec3{0, 1, 0});
  b1->releaseAt(0, 4); b1->releaseAt(0, 5);
  b2->releaseAt(1, 4); b2->releaseAt(1, 5);

  m.fixNode(n0, true, true, true, true, false, true);
  m.fixNode(n1, true, true, true, true, false, true);
  m.addNodeForce(nc, {0, 0, -P});             // 全局 -Z 跨中集中力

  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) { std::printf("      %s\n", r.message.c_str()); return; }
  checkNear(r.residual, 0.0, 1e-12, "平衡残差 ~ 机器精度");

  // 铰端局部 My/Mz 必须 = 0（局部坐标下弯矩 = 0，而非全局）
  const double mI_My = r.beamForces[0].ef.My, mI_Mz = r.beamForces[0].ef.Mz;
  const double mJ_My = r.beamForces[1].ef.Myj, mJ_Mz = r.beamForces[1].ef.Mzj;
  checkNear(mI_My, 0.0, 1e-9, "铰端 n0 局部 My = 0");
  checkNear(mI_Mz, 0.0, 1e-9, "铰端 n0 局部 Mz = 0");
  checkNear(mJ_My, 0.0, 1e-9, "铰端 n1 局部 My = 0");
  checkNear(mJ_Mz, 0.0, 1e-9, "铰端 n1 局部 Mz = 0");

  // 跨中弯矩（elem1 i 端 = 跨中）= P·L·cos45°/4（竖向力沿局部 z 的分量 × L/4）
  const double mMid = r.beamForces[1].ef.My;   // 绕局部 y = 全局 Y 轴
  checkNear(mMid, P * L * c / 4.0, 1e-9, "跨中弯矩 = P·L·cos45°/4");
}

// ===========================================================================
//  5. 旋转用例：绕轴向旋转梁 + 释放验证（三根梁绕 x 轴不同朝向）
//
//  解析解（含剪切修正）：竖向力 F = (0,0,-P) 在局部坐标分解
//    fy = F·ey，fz = F·ez
//    wy = fy·L³/(48E·Iz)·(1+Φz)，Φz = 12EIz/(G·Asz·L²)   —— 绕 z' 弯曲
//    wz = fz·L³/(48E·Iy)·(1+Φy)，Φy = 12EIy/(G·Asy·L²)   —— 绕 y' 弯曲
//    全局 Z 位移解析值 = wy·ey.z + wz·ez.z，与实际节点位移逐点对照。
//  这验证了"弯曲平面随 up 旋转而旋转"（up=-Y 时荷载沿 z'，弱轴弯曲）。
// ===========================================================================
static void testRotatedHinged() {
  std::printf("\n== 5. 绕轴旋转的铰接梁 ==\n");
  const double L = 6.0, P = 100.0;
  const auto sec = section::rect(0.3, 0.5);
  const auto mat = Material::concreteC(30);
  const double Iz = sec.Iz, Iy = sec.Iy;
  const double G = mat.shearModulus();
  const double phiZ = 12.0 * mat.E * Iz / (G * sec.Asz * L * L);
  const double phiY = 12.0 * mat.E * Iy / (G * sec.Asy * L * L);

  // 朝向矩阵：梁沿 +X，up 依次取 -Z / +Z / -Y
  // 竖直荷载 -Z（与 up 不同向 → 弯曲平面随 up 旋转）
  struct Case { Vec3 up; const char* name; };   // Vec3 为聚合类型，可列表初始化
  const Case cases[] = {
      {{0, 0, -1}, "up = -Z（基准）"},
      {{0, 0, 1}, "up = +Z（旋转 180°）"},
      {{0, -1, 0}, "up = -Y（绕轴旋转 90°，弱轴弯曲）"},
  };
  for (const auto& cs : cases) {
    Model m;
    const Id n0 = m.addNode({0, 0, 0}, 1);
    const Id nc = m.addNode({L / 2, 0, 0}, 1);
    const Id n1 = m.addNode({L, 0, 0}, 1);
    BeamElement* b1 = m.addBeam(n0, nc, sec, mat, cs.up);
    BeamElement* b2 = m.addBeam(nc, n1, sec, mat, cs.up);
    b1->releaseAt(0, 4); b1->releaseAt(0, 5);
    b2->releaseAt(1, 4); b2->releaseAt(1, 5);

    m.fixNode(n0, true, true, true, true, false, true);
    m.fixNode(n1, true, true, true, true, false, true);
    m.addNodeForce(nc, {0, 0, -P});

    m.assignDofs();
    StaticAnalysis sa(m);
    const auto r = sa.solve();
    if (!r.ok) { check(false, "求解成功"); continue; }

    // 解析解：F 分解到局部 y'/z'，两方向弯曲独立、线性叠加，再投影回全局 Z
    const Vec3 ey = b1->core().localY(), ez = b1->core().localZ();
    const double fy = (0.0) * ey.x + (0.0) * ey.y + (-P) * ey.z;
    const double fz = (0.0) * ez.x + (0.0) * ez.y + (-P) * ez.z;
    const double wy = fy * L * L * L / (48.0 * mat.E * Iz) * (1.0 + phiZ);
    const double wz = fz * L * L * L / (48.0 * mat.E * Iy) * (1.0 + phiY);
    const double uZana = wy * ey.z + wz * ez.z;    // 全局 Z 位移（向下为负）

    char buf[128];
    std::snprintf(buf, sizeof(buf), "跨中挠度含剪切（%s）", cs.name);
    checkNear(r.u[static_cast<size_t>(nc) * 6 + 2], uZana, 1e-9, buf);
    // 铰端弯矩恒 0
    checkNear(r.beamForces[0].ef.Mz, 0.0, 1e-9, "铰端 Mz = 0");
    checkNear(r.beamForces[1].ef.Mzj, 0.0, 1e-9, "铰端 Mzj = 0");
  }
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  端部释放（铰）验证\n");
  std::printf("======================================================\n");

  testCondensation();
  testSimplySupportedHinged();
  testCantileverRelease();
  testTiltedHinged();
  testRotatedHinged();

  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}