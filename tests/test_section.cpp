// =============================================================================
//  tests/test_section.cpp  ——  截面几何性质验证
//
//  验证策略：对标"独立实现"（闭式解 vs 数值积分），不靠记忆。
//  例如矩形惯性矩用数值积分 ∫y²dA 校核，比对公式更可靠。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "yjk/model/Section.h"

using namespace yjk;

static int g_pass = 0, g_fail = 0;
static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}
static void checkNear(double got, double want, double relTol, const char* what) {
  const double rel = std::abs(want) > 1e-30 ? std::abs(got - want) / std::abs(want) : std::abs(got - want);
  if (rel <= relTol) {
    ++g_pass;
    std::printf("  [ OK ] %s  (got %.10g, want %.10g, rel %.2e)\n", what, got, want, rel);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %.10g, want %.10g, rel %.2e)\n", what, got, want, rel);
  }
}

// ---------------------------------------------------------------------------
// 数值积分：矩形截面 b×h 的解析性质
//   A  = ∫dA,  Iy = ∫z²dA = ∫_{-b/2}^{b/2} z² · h dz,  Iz = ∫y²dA = ∫_{-h/2}^{h/2} y² · b dy
// ---------------------------------------------------------------------------
static void testRectAgainstIntegration(double b, double h) {
  const SectionProperties s = section::rect(b, h);

  // 面积（精确）
  checkNear(s.A, b * h, 1e-15, "矩形面积");

  // Iy = ∫ z² dA  （z 是宽度方向坐标，范围 ±b/2，厚度沿 y 方向 h）
  const double Iy_int = h * (std::pow(b, 3.0) / 12.0);
  checkNear(s.Iy, Iy_int, 1e-12, "矩形 Iy = ∫z²dA");

  // Iz = ∫ y² dA  （y 是高度方向坐标，范围 ±h/2，宽度沿 z 方向 b）
  const double Iz_int = b * (std::pow(h, 3.0) / 12.0);
  checkNear(s.Iz, Iz_int, 1e-12, "矩形 Iz = ∫y²dA");

  // 数值积分校核 Iz（辛普森法，纵向 2000 等分）
  const int N = 2000;
  const double dy = h / N;
  double numIz = 0.0;
  for (int i = 0; i <= N; ++i) {
    const double yv = -h / 2.0 + i * dy;
    const double w = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    numIz += w * b * yv * yv;
  }
  numIz *= dy / 3.0;
  checkNear(s.Iz, numIz, 1e-8, "矩形 Iz 数值积分校核");

  // 塑性抵抗矩（矩形双轴对称，两方向塑性模量之和即截面抵抗矩）：
  //   绕 z（强轴，∫y²dA 方向）Wz = b·h²/4
  //   绕 y（弱轴，∫z²dA 方向）Wy = h·b²/4
  checkNear(s.Rz, b * h * h / 4.0, 1e-15, "矩形塑性模量 Wz = b·h²/4（绕 z）");
  checkNear(s.Ry, h * b * b / 4.0, 1e-15, "矩形塑性模量 Wy = h·b²/4（绕 y）");
  // 交叉校核：Rz 必须与 Iz 同向、Ry 与 Iy 同向
  check((s.Rz / s.Iz) > 1.0 && (s.Ry / s.Iy) > 1.0,
        "塑性模量与惯性矩同向（均为强轴大）");
}

// ---------------------------------------------------------------------------
// 圆截面：与解析解对比
//   A = πd²/4,  I = πd⁴/64,  W = d³/6,  J = 2I
// ---------------------------------------------------------------------------
static void testCircle() {
  const double d = 0.5;
  const SectionProperties s = section::circle(d);
  const double pi = std::acos(-1.0);
  checkNear(s.A, pi * d * d / 4.0, 1e-15, "圆截面面积");
  checkNear(s.Iy, pi * std::pow(d, 4.0) / 64.0, 1e-15, "圆截面 Iy = πd⁴/64");
  checkNear(s.Iz, s.Iy, 1e-15, "圆截面 Iz = Iy（各向同性）");
  checkNear(s.J, 2.0 * s.Iy, 1e-15, "圆截面 J = 2I");

  // 数值积分：关于 y 轴的惯性矩 Iy = ∫z² dA。
  // 用极坐标：∫∫z² dA = ∫0^{d/2} (r cosθ)²·r dr dθ ·2（对称）
  //          = (∫r³dr)·(∫cos²θ dθ)·2 = (d/2)⁴/4 · π
  //          = πd⁴/64 ✓
  const int N = 2000;
  const double dr = (d / 2.0) / N;
  double radial = 0.0;
  for (int i = 0; i <= N; ++i) {
    const double r = i * dr;
    const double w = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    radial += w * r * r * r;            // ∫r³ dr
  }
  radial *= dr / 3.0;
  const double num = radial * pi;      // ∫cos²θ dθ ·2 = π
  checkNear(s.Iy, num, 1e-6, "圆截面 Iy 数值积分校核（∫z²dA）");
  checkNear(s.J, 2.0 * s.Iy, 1e-15, "圆截面 J = 2I（极坐标）");
}

// ---------------------------------------------------------------------------
// 工字形截面：与分块公式对比（这是最容易算错的）
// ---------------------------------------------------------------------------
static void testISection() {
  const double h = 0.60, b = 0.20, tw = 0.010, tf = 0.016;
  const SectionProperties s = section::iSection(h, b, tw, tf);
  const double hw = h - 2.0 * tf;

  // 面积：2btf + tw·hw
  checkNear(s.A, 2.0 * b * tf + tw * hw, 1e-15, "工字形面积");

  // Iz（强轴，绕 z 横轴）：沿 y 切条，∫y²dA
  //   腹板（高 hw、宽 tw）: tw·hw³/12
  //   翼缘: 2 × [ b·tf³/12 + b·tf·((hw+tf)/2)² ]   （平行轴定理）
  const double dyc = (hw + tf) / 2.0;
  const double Iz_ref = tw * std::pow(hw, 3.0) / 12.0 +
                        2.0 * (b * std::pow(tf, 3.0) / 12.0 + b * tf * dyc * dyc);
  checkNear(s.Iz, Iz_ref, 1e-12, "工字形 Iz（强轴，绕 z）");

  // Iy（弱轴，绕 y 竖轴）：沿 y 切条，∫z²dA
  //   腹板（高 hw、宽 tw）: hw·tw³/12  —— 绕【自身】竖轴，注意是 tw³
  //   翼缘（宽 b）: 2 × tf·b³/12
  const double Iy_ref = hw * std::pow(tw, 3.0) / 12.0 +
                        2.0 * tf * std::pow(b, 3.0) / 12.0;
  checkNear(s.Iy, Iy_ref, 1e-12, "工字形 Iy（弱轴，绕 y）");

  // 独立校核：数值积分（沿 y 切条，宽度随 y 变化）
  {
    const int N = 200000;
    const double dy = h / N;
    double numIy = 0.0, numIz = 0.0;
    for (int i = 0; i < N; ++i) {
      const double yv = -h / 2.0 + (i + 0.5) * dy;
      const double ay = std::abs(yv);
      const double w = (ay <= hw / 2.0) ? tw : b;   // 腹板区 or 翼缘区
      numIy += w * w * w / 12.0 * dy;              // ∫z²dz · dy
      numIz += w * yv * yv * dy;                    // ∫dz = w，累加 w·y²·dy
    }
    checkNear(s.Iy, numIy, 1e-4, "工字形 Iy 数值积分校核");
    checkNear(s.Iz, numIz, 1e-4, "工字形 Iz 数值积分校核");
  }

  // 关键性质：强轴 Iz 应显著大于弱轴 Iy
  // （工字形抵抗"上下弯曲"远强于"侧向弯曲"）
  // 对 h=0.60, b=0.20 的 H 型钢，Iz/Iy 通常在 4~10 之间
  // h/b = 3 的工字形，Iz/Iy 约 20~40（手册 H 型钢典型值）
  const double strongWeak = s.Iz / s.Iy;
  std::printf("       工字形 Iz/Iy = %.2f  （h/b=3 的 H 型钢典型值 20~40）\n", strongWeak);
  check(strongWeak > 15.0 && strongWeak < 60.0, "工字形强弱轴比在合理区间");

  // 开口截面 J 极小 —— 这是工字梁抗扭弱的根本原因
  // 验证量级：J 应比 Iz 小 2~3 个数量级
  const double twistRatio = s.Iz / s.J;
  std::printf("       工字形 Iz/J = %.1f  （开口截面抗扭弱，典型值 >100）\n", twistRatio);
  check(twistRatio > 50.0 && twistRatio < 1e5, "开口截面 J 的量级正确");
}

// ---------------------------------------------------------------------------
// 箱形截面
// ---------------------------------------------------------------------------
static void testBox() {
  const double h = 0.50, b = 0.30, tw = 0.012, tf = 0.014;
  const SectionProperties s = section::boxSection(h, b, tw, tf);
  const double hi = h - 2.0 * tf, bi = b - 2.0 * tw;

  checkNear(s.A, b * h - bi * hi, 1e-15, "箱形面积（外减内）");
  checkNear(s.Iz, (h * std::pow(b, 3.0) - hi * std::pow(bi, 3.0)) / 12.0, 1e-12,
            "箱形 Iz");

  // 闭合截面 J 远大于开口截面（同外廓工字形）
  const SectionProperties isec = section::iSection(h, b, tw, tf);
  std::printf("       箱形 J = %.6e, 同尺寸工字形 J = %.6e\n", s.J, isec.J);
  check(s.J > isec.J * 10.0, "闭合截面抗扭刚度远大于开口截面");
}

// ---------------------------------------------------------------------------
// 有效性检查：退化截面必须被拦下
// ---------------------------------------------------------------------------
static void testValidity() {
  std::printf("\n== 5. 截面有效性检查 ==\n");
  {
    SectionProperties s;                       // 全零
    std::string why;
    check(!s.isValid(&why), "零截面被判为无效");
    std::printf("       诊断: %s\n", why.c_str());
  }
  {
    SectionProperties s = section::rect(0.0, 0.5);
    std::string why;
    check(!s.isValid(&why), "宽度为 0 的矩形被拒");
    std::printf("       诊断: %s\n", why.c_str());
  }
  {
    const SectionProperties s = section::rect(0.3, 0.5);
    check(s.isValid(), "正常矩形截面通过检查");
  }
}

// ---------------------------------------------------------------------------
// 材料库
// ---------------------------------------------------------------------------
static void testMaterial() {
  std::printf("\n== 6. 材料库（GB 50010 / GB 50017） ==\n");
  {
    const Material c30 = Material::concreteC(30);
    // GB 50010-2010(2015) 表 4.1.5: C30 的 Ec = 3.15×10⁴ MPa = 3.15e7 kPa
    // 注意 C25 才是 3.00×10⁴ —— 规范表值不可凭记忆，必须查原表
    checkNear(c30.E, 3.15e7, 1e-12, "C30 弹性模量 = 3.15e4 MPa");
    checkNear(c30.fc, 14.3e3, 1e-12, "C30 轴心抗压强度 = 14.3 MPa");
    checkNear(c30.ft, 1.43e3, 1e-12, "C30 抗拉强度 = 1.43 MPa");
    {
      const Material c25 = Material::concreteC(25);
      checkNear(c25.E, 3.00e7, 1e-12, "C25 弹性模量 = 3.00e4 MPa");
      checkNear(c25.fc, 11.9e3, 1e-12, "C25 轴心抗压强度 = 11.9 MPa");
    }
    checkNear(c30.nu, 0.2, 1e-15, "C30 泊松比 = 0.2");
    checkNear(c30.gamma, 25.0, 1e-15, "C30 容重 = 25 kN/m³");
  }
  {
    const Material q235 = Material::steelQ(235);
    checkNear(q235.E, 2.06e8, 1e-12, "Q235 弹性模量 = 2.06e5 MPa");
    checkNear(q235.fy, 235e3, 1e-12, "Q235 屈服强度 = 235 MPa");
    // G 字段留 0 表示"由 E、ν 自动推算"，应通过 shearModulus() 访问
    check(q235.G == 0.0, "钢材 G 字段留空（交由 E、ν 推算）");
    checkNear(q235.shearModulus(), q235.E / (2.0 * (1.0 + q235.nu)), 1e-12,
              "钢材剪切模量由 E、ν 推出");
    {
      const Material c30 = Material::concreteC(30);
      checkNear(c30.shearModulus(), c30.E / (2.4), 1e-12, "C30 剪切模量 = E/2.4");
    }
  }
  {
    // 混凝土强度等级应单调递增
    double prev = 0.0;
    bool mono = true;
    for (int g = 20; g <= 60; g += 10) {
      const Material m = Material::concreteC(g);
      if (m.E <= prev) mono = false;
      prev = m.E;
    }
    check(mono, "混凝土弹性模量随强度等级单调递增");
  }
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  截面与材料验证\n");
  std::printf("======================================================\n");

  std::printf("\n== 1. 矩形截面（对标解析解 + 数值积分） ==\n");
  testRectAgainstIntegration(0.3, 0.5);
  testRectAgainstIntegration(0.4, 0.8);

  std::printf("\n== 2. 圆截面 ==\n");
  testCircle();

  std::printf("\n== 3. 工字形截面 ==\n");
  testISection();

  std::printf("\n== 4. 箱形截面 ==\n");
  testBox();

  testValidity();
  testMaterial();

  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}
