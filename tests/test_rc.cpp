// =============================================================================
//  tests/test_rc.cpp —— 钢筋混凝土正截面承载力验算模块（P4-T1）
//
//  依据：GB 50010-2010《混凝土结构设计规范》第 6 章 / §8.5.1。
//  判据全部为【手算锚点】：期望值与规范公式的独立解析推导逐值对比，
//  不是和实现代码的"同一公式再算一遍"互证。
//
//  覆盖：
//    · 应力图形参数 α1/β1/εcu（C30 常数段、C55+ 内插）与界限 ξb
//    · 单筋受弯：正常 / 少筋（按最小配筋率控制）/ 超筋（返回界限配筋）
//    · 双筋受弯：正常（ξ≤ξb）、As' 已足够（§6.2.14 对受压筋取矩）、
//                x<2as'（§6.2.14）、As' 不足超筋
//    · 轴压：稳定系数 φ 查表与线性内插
//    · 偏心受压：大偏心（对称/非对称配筋）、小偏心（对称/非对称配筋）
//    · 纯扭 / 剪扭：Wt、ζ、βt、Tu、Vu 及 ζ/βt 夹取边界
//    · 非法输入降级（返回不可用状态而非 NaN/崩溃）
//
//  单位约定（与 include/yjk/check/RC.h 一致）：m / kN / kPa。
// =============================================================================
#include <cmath>
#include <cstdio>

#include "yjk/check/RC.h"

using namespace yjk::check;

static int g_fail = 0;
static int g_pass = 0;

static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}
static void checkNear(double got, double want, double relTol, const char* what) {
  const double err = std::abs(got - want);
  const double rel = std::abs(want) > 1e-30 ? err / std::abs(want) : err;
  if (rel <= relTol) { ++g_pass;
    std::printf("  [ OK ] %s  (got %.10g, want %.10g, rel %.2e)\n",
                what, got, want, rel);
  } else { ++g_fail;
    std::printf("  [FAIL] %s  (got %.10g, want %.10g, rel %.2e)\n",
                what, got, want, rel);
  }
}

// -----------------------------------------------------------------------------
//  1. 应力图形参数、界限相对受压区高度、最小配筋率
// -----------------------------------------------------------------------------
static void testStressBlock() {
  std::printf("\n== 1. 应力图形参数 / ξb / 最小配筋率 ==\n");
  // C30：α1=1.0, β1=0.8, εcu=0.0033（规范表 6.2.6 常数段）
  auto sb30 = stressBlockFromFc(14300.0);
  checkNear(sb30.alpha1, 1.0, 1e-12, "C30 α1");
  checkNear(sb30.beta1, 0.8, 1e-12, "C30 β1");
  checkNear(sb30.epsCu, 0.0033, 1e-12, "C30 εcu");
  // C60：α1=0.98, β1=0.78, εcu=0.0032（规范表值，直接命中等级点）
  auto sb60 = stressBlockFromFc(27500.0);
  checkNear(sb60.alpha1, 0.98, 1e-12, "C60 α1");
  checkNear(sb60.beta1, 0.78, 1e-12, "C60 β1");
  checkNear(sb60.epsCu, 0.0032, 1e-12, "C60 εcu");
  // C55~C60 之间线性内插（fc=25.3→26.4，t=0.5）
  auto sb55_ = stressBlockFromFc(26400.0);
  checkNear(sb55_.alpha1, 0.985, 1e-9, "C55.5 α1 内插");
  checkNear(sb55_.beta1, 0.785, 1e-9, "C55.5 β1 内插");
  checkNear(sb55_.epsCu, 0.003225, 1e-9, "C55.5 εcu 内插");

  // ξb（HRB400 fy=360MPa, Es=2.0e5MPa, C30）
  //   ξb = β1/(1+fy/(Es·εcu)) = 0.8/(1+360/(2e5·0.0033)) = 0.517647…
  double xiB = relativeBoundaryHeight(360000.0, 2.0e8, sb30);
  checkNear(xiB, 0.5176470588, 1e-9, "ξb HRB400/C30");

  // 最小配筋率：ρmin = max(0.2%, 0.45·ft/fy)
  //   C30 ft=1.43MPa：0.45·1.43/360 = 0.17875% < 0.2% → 0.2%
  checkNear(minReinfRatio(1430.0, 360000.0), 0.002, 1e-12, "ρmin=C30/HRB400");
  //   ft=1.43, fy=300：0.45·1.43/300 = 0.2145% > 0.2% → 0.2145%
  checkNear(minReinfRatio(1430.0, 300000.0), 0.002145, 1e-12, "ρmin HRB335 控制");
}

// -----------------------------------------------------------------------------
//  2. 单筋矩形截面受弯（GB 50010 §6.2.10）
//     b×h=300×600, C30(fc=14.3MPa/ft=1.43MPa), HRB400, as=40mm → h0=560mm
// -----------------------------------------------------------------------------
static FlexuralInput beam300x600() {
  FlexuralInput in;
  in.b = 0.3; in.h = 0.6;
  in.fc = 14300.0; in.ft = 1430.0;
  in.fy = 360000.0; in.Es = 2.0e8;
  in.as = 0.04; in.as_p = 0.04;
  return in;
}

static void testFlexuralSingle() {
  std::printf("\n== 2. 单筋受弯 ==\n");
  auto in = beam300x600();

  // M=300kN·m：αs=M/(α1·fc·b·h0²)=0.222992, ξ=1-√(1-2αs)=0.255677,
  // x=ξ·h0=0.143179, As=α1·fc·b·x/fy=1706.2mm² > AsMin=360mm²
  auto r1 = flexuralSingle(in, 300.0);
  check(r1.status == FlexuralStatus::Ok, "单筋 M=300 正常");
  checkNear(r1.xi, 0.255676541, 1e-6, "单筋 M=300 ξ");
  checkNear(r1.x, 0.143178863, 1e-6, "单筋 M=300 x");
  checkNear(r1.As, 1706.2148e-6, 1e-6, "单筋 M=300 As");
  checkNear(r1.AsMin, 360.0e-6, 1e-12, "单筋 AsMin=0.2%·b·h");

  // M=40kN·m：As=201.5mm² < AsMin=360mm² → 少筋，按最小配筋率返回
  auto r2 = flexuralSingle(in, 40.0);
  check(r2.status == FlexuralStatus::BelowMinReinf, "单筋 M=40 少筋判别");
  checkNear(r2.As, 360.0e-6, 1e-12, "单筋 M=40 As=AsMin");

  // M=600kN·m：ξ=0.6713 > ξb=0.5176 → 超筋；
  // 返回界限配筋 As=α1·fc·b·ξb·h0/fy=3454.4mm² 作参考
  auto r3 = flexuralSingle(in, 600.0);
  check(r3.status == FlexuralStatus::OverReinforced, "单筋 M=600 超筋判别");
  checkNear(r3.xi, 0.5176470588, 1e-9, "单筋超筋回退 ξ=ξb");
  checkNear(r3.As, 3454.4314e-6, 1e-6, "单筋超筋界限配筋");

  // 非法输入（M<0）降级为不可用，不崩溃
  auto r4 = flexuralSingle(in, -10.0);
  check(r4.status == FlexuralStatus::OverReinforced, "单筋 非法输入降级");
  checkNear(r4.As, 0.0, 1e-12, "单筋 非法输入 As=0");
}

// -----------------------------------------------------------------------------
//  3. 双筋矩形截面受弯（GB 50010 §6.2.10 / §6.2.14）
//     同截面，As'=4Φ20=1256.6mm²，fy'=fy=360MPa
// -----------------------------------------------------------------------------
static void testFlexuralDouble() {
  std::printf("\n== 3. 双筋受弯 ==\n");
  auto in = beam300x600();
  const double As_p = 1256.6371e-6;  // 4Φ20
  // M1 = M - fy'·As'·(h0-as') = 500 - 360e3·1.25664e-3·0.52 = 264.758kN·m
  // ξ=0.221277 < ξb, x=0.123915 > 2as'=0.08
  // As = (α1·fc·b·x + fy'·As')/fy = 2733.3mm²
  auto r1 = flexuralDouble(in, 500.0, As_p);
  check(r1.status == FlexuralStatus::Ok, "双筋 M=500 正常");
  checkNear(r1.x, 0.123915249, 1e-6, "双筋 M=500 x");
  checkNear(r1.As, 2733.2938e-6, 1e-6, "双筋 M=500 As");

  // M=200kN·m：M1<0（As' 本身已足够）→ §6.2.14
  //   As = M/(fy·(h0-as')) = 200/(360e3·0.52) = 1068.4mm²
  auto r2 = flexuralDouble(in, 200.0, As_p);
  check(r2.status == FlexuralStatus::Ok, "双筋 M=200 (M1<0) 正常");
  checkNear(r2.As, 1068.3761e-6, 1e-6, "双筋 M=200 §6.2.14 As");

  // M=300kN·m：x=0.0276 < 2as'=0.08 → §6.2.14 对受压筋合力点取矩
  //   As = M/(fy·(h0-as')) = 300/(360e3·0.52) = 1602.6mm²
  auto r3 = flexuralDouble(in, 300.0, As_p);
  check(r3.status == FlexuralStatus::Ok, "双筋 M=300 (x<2as') 正常");
  checkNear(r3.As, 1602.5641e-6, 1e-6, "双筋 M=300 x<2as' As");

  // M=800kN·m：M1=564.758, ξ=0.5995 > ξb → As' 不足，超筋
  auto r4 = flexuralDouble(in, 800.0, As_p);
  check(r4.status == FlexuralStatus::OverReinforced, "双筋 M=800 超筋判别");
  checkNear(r4.xi, 0.5176470588, 1e-9, "双筋超筋回退 ξ=ξb");

  // As'=0 退化为单筋
  auto r5 = flexuralDouble(in, 300.0, 0.0);
  auto s5 = flexuralSingle(in, 300.0);
  checkNear(r5.As, s5.As, 1e-9, "双筋 As'=0 退化一致");
}

// -----------------------------------------------------------------------------
//  4. 轴心受压（GB 50010 §6.2.15）
//     b×h=400×500, C30, HRB400, as=as'=40mm, As'=4Φ20=1256mm²
// -----------------------------------------------------------------------------
static ColumnInput column400x500() {
  ColumnInput in;
  in.b = 0.4; in.h = 0.5;
  in.fc = 14300.0; in.ft = 1430.0;
  in.fy = 360000.0; in.Es = 2.0e8;
  in.as = 0.04; in.as_p = 0.04;
  return in;
}

static void testAxial() {
  std::printf("\n== 4. 轴心受压 ==\n");
  auto in = column400x500();
  const double As_p = 1256.6371e-6;

  // l0=4.0m → l0/b=10 → φ=0.98（表 6.2.15 直接命中）
  // Nu = 0.9·φ·(fc·A + fy'·As') = 0.9·0.98·(14300·0.2 + 360e3·1.25664e-3)
  double nu10 = axialCapacity(in, 4.0, As_p);
  double want10 = 0.9 * 0.98 * (14300.0 * 0.2 + 360000.0 * As_p);
  checkNear(nu10, want10, 1e-9, "轴压 l0/b=10 φ=0.98");

  // l0=3.6m → l0/b=9 → φ 在 (8→1.00, 10→0.98) 间线性内插 = 0.99
  double nu9 = axialCapacity(in, 3.6, As_p);
  double want9 = 0.9 * 0.99 * (14300.0 * 0.2 + 360000.0 * As_p);
  checkNear(nu9, want9, 1e-9, "轴压 l0/b=9 φ 内插=0.99");

  // As'=0（仅混凝土贡献）
  double nu0 = axialCapacity(in, 4.0, 0.0);
  checkNear(nu0, 0.9 * 0.98 * 14300.0 * 0.2, 1e-9, "轴压 As'=0");
}

// -----------------------------------------------------------------------------
//  5. 偏心受压（GB 50010 §6.2.17）
//     同 column400x500，h0=460mm
// -----------------------------------------------------------------------------
static void testEccentric() {
  std::printf("\n== 5. 偏心受压 ==\n");
  auto in = column400x500();
  const double As = 1256.6371e-6;   // 4Φ20（受拉侧）
  const double As_p = 1256.6371e-6; // 4Φ20（受压侧）

  // ---- 大偏心：N=800, M=250（对称配筋）----
  // e0=M/N=0.3125, ea=max(0.02, h/30)=0.02, ei=0.3325, e=ei+h/2-as=0.5425
  // 0.5A·x² + A(e-h0)·x + [C(e-h0+as') - Cp·e] = 0
  // A=α1·fc·b=5720, C=Cp=fy'·As'=452.389
  // x=0.1881308, ξ=0.408980 < ξb, Nu=A·x + C - Cp = 1076.11
  auto r1 = columnCapacity(in, 800.0, 250.0, As, As_p);
  check(r1.largeEccentric, "大偏心判别（对称）");
  checkNear(r1.x, 0.1881308381, 1e-6, "大偏心 x（对称）");
  checkNear(r1.Nu, 1076.108394, 1e-5, "大偏心 Nu（对称）");
  check(r1.satisfied, "大偏心满足 N<=Nu");
  checkNear(r1.xiB, 0.5176470588, 1e-6, "大偏心 ξb 透出");

  // ---- 小偏心：N=1200, M=150（对称配筋）----
  // e0=0.125, ei=0.145, e=0.355；大偏试探 x=0.3834 → ξ=0.833>ξb → 转小偏
  // denom=(ξb-β1)·h0=-0.129882, 解得 x=0.5010711, ξ=1.089285
  // σs=fy·(ξ-β1)/(ξb-β1)=-368878kPa, Nu=A·x+C-σs·As=3782.01
  auto r2 = columnCapacity(in, 1200.0, 150.0, As, As_p);
  check(!r2.largeEccentric, "小偏心判别（对称）");
  checkNear(r2.x, 0.5010711378, 1e-6, "小偏心 x（对称）");
  checkNear(r2.Nu, 3782.012381, 1e-5, "小偏心 Nu（对称）");
  check(r2.satisfied, "小偏心满足 N<=Nu");

  // ---- 大偏心（非对称：As=4Φ20, As'=2Φ20=628.3mm²）----
  // C=fy'·As'=226.195, Cp=452.389；x=0.2054748, ξ=0.446684
  // Nu=A·x + C - Cp = 949.12
  auto r3 = columnCapacity(in, 800.0, 250.0, As, 628.3185e-6);
  check(r3.largeEccentric, "大偏心判别（非对称）");
  checkNear(r3.x, 0.2054747553, 1e-6, "大偏心 x（非对称）");
  checkNear(r3.Nu, 949.1209046, 1e-5, "大偏心 Nu（非对称）");

  // ---- 小偏心（非对称：As=2Φ20, As'=4Φ20）----
  // x=0.4322931, ξ=0.939768>ξb, σs=-178049kPa, Nu=3037.07
  auto r4 = columnCapacity(in, 1200.0, 150.0, 628.3185e-6, As_p);
  check(!r4.largeEccentric, "小偏心判别（非对称）");
  checkNear(r4.x, 0.4322930566, 1e-6, "小偏心 x（非对称）");
  checkNear(r4.Nu, 3037.074247, 1e-5, "小偏心 Nu（非对称）");

  // 非法输入（M=0, N=0）→ Nu=0 且不满足
  auto r5 = columnCapacity(in, 0.0, 0.0, As, As_p);
  checkNear(r5.Nu, 0.0, 1e-12, "偏心受压 非法输入 Nu=0");
  check(!r5.satisfied, "偏心受压 非法输入不满足");
}

// -----------------------------------------------------------------------------
//  6. 受扭（GB 50010 §6.4）
//     b×h=400×600, C30, HRB400 fy=fyv=360MPa, c=40mm, h0=560mm
// -----------------------------------------------------------------------------
static TorsionInput torsion400x600() {
  TorsionInput in;
  in.b = 0.4; in.h = 0.6; in.h0 = 0.56;
  in.fc = 14300.0; in.ft = 1430.0;
  in.fyv = 360000.0; in.fy = 360000.0;
  in.c = 0.04;
  return in;
}

static void testTorsion() {
  std::printf("\n== 6. 受扭 ==\n");
  auto in = torsion400x600();

  // Wt = b²·(3h-b)/6 = 0.16·1.4/6 = 0.0373333…m³
  checkNear(torsionPlasticModulus(0.4, 0.6), 0.03733333333333333, 1e-12, "Wt");

  // ---- 纯扭：单肢 φ8@100 (Ast1=50.3mm²), 受扭纵筋 4Φ16 (Astl=804mm²) ----
  // bcor=0.32, hcor=0.52, ucor=1.68, Acor=0.1664
  // ζ=fy·Astl·s/(fyv·Ast1·ucor)=804·0.1/(50.3·1.68)=0.951434
  // Tu=0.35·ft·Wt + 1.2·√ζ·fyv·Ast1·Acor/s = 18.685 + 35.269 = 53.954kN·m
  auto pt = pureTorsion(in, 50.3e-6, 804.0e-6, 0.1);
  checkNear(pt.zeta, 0.9514342516, 1e-9, "纯扭 ζ");
  checkNear(pt.Tu, 53.95443880, 1e-8, "纯扭 Tu");

  // ζ 夹取下界：Astl 很小 → ζ_raw=0.118 → 夹到 0.6
  auto ptLo = pureTorsion(in, 50.3e-6, 100.0e-6, 0.1);
  checkNear(ptLo.zeta, 0.6, 1e-12, "纯扭 ζ 下界夹取");
  // ζ 夹取上界：Astl 很大 → ζ_raw=2.367 → 夹到 1.7
  auto ptHi = pureTorsion(in, 50.3e-6, 2000.0e-6, 0.1);
  checkNear(ptHi.zeta, 1.7, 1e-12, "纯扭 ζ 上界夹取");

  // ---- 剪扭：V=400, T=60, 双肢 φ8@100 (Asv=100.6mm²) ----
  // βt=1.5/(1+0.5·V·Wt/(T·b·h0))=1.5/(1+0.55556)=0.964286
  // Tu=0.35·βt·ft·Wt + 箍筋项 = 53.287kN·m
  // Vu=(1.5-βt)·0.7·ft·b·h0 + fyv·Asv·h0/s = 322.93kN
  auto st = shearTorsion(in, 400.0, 60.0, 50.3e-6, 804.0e-6,
                         2 * 50.3e-6, 0.1);
  checkNear(st.betaT, 0.9642857143, 1e-9, "剪扭 βt");
  checkNear(st.Tu, 53.28710547, 1e-8, "剪扭 Tu");
  checkNear(st.Vu, 322.9296, 1e-8, "剪扭 Vu");

  // βt 夹取下界：V 很大 → βt_raw=0.085 → 0.5
  auto stLo = shearTorsion(in, 2000.0, 10.0, 50.3e-6, 804.0e-6,
                           2 * 50.3e-6, 0.1);
  checkNear(stLo.betaT, 0.5, 1e-12, "剪扭 βt 下界夹取");
  // βt 夹取上界：T 很大 → βt_raw=1.4999 → 1.0
  auto stHi = shearTorsion(in, 1.0, 1000.0, 50.3e-6, 804.0e-6,
                           2 * 50.3e-6, 0.1);
  checkNear(stHi.betaT, 1.0, 1e-12, "剪扭 βt 上界夹取");
}

int main() {
  std::printf("test_rc —— 钢筋混凝土正截面承载力验算（P4-T1）\n");
  testStressBlock();
  testFlexuralSingle();
  testFlexuralDouble();
  testAxial();
  testEccentric();
  testTorsion();

  std::printf("\npass=%d fail=%d\n", g_pass, g_fail);
  return g_fail ? 1 : 0;
}