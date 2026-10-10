// =============================================================================
//  src/check/RC.cpp —— 钢筋混凝土构件正截面承载力验算实现（P4-T1）
// =============================================================================
#include "yjk/check/RC.h"

#include <algorithm>
#include <cmath>

namespace yjk {
namespace check {

namespace {

// 混凝土强度等级参数表（GB 50010-2010 表 6.2.6 / 规范配套参数）
// fc 单位为 MPa（转换在调用处处理）。C50 及以下取规范固定值；
// C55~C80 的 α1/β1 取自规范表 6.2.6，εcu 按 §6.2.1 注线性内插。
struct GradePoint {
  double fcMpa;    // 混凝土轴心抗压强度设计值 (MPa)
  double alpha1;
  double beta1;
  double epsCu;
};

constexpr GradePoint kGrades[] = {
    {14.3, 1.00, 0.80, 0.00330},  // C30
    {16.7, 1.00, 0.80, 0.00330},  // C35
    {19.1, 1.00, 0.80, 0.00330},  // C40
    {21.1, 1.00, 0.80, 0.00330},  // C45
    {23.1, 1.00, 0.80, 0.00330},  // C50
    {25.3, 0.99, 0.79, 0.00325},  // C55
    {27.5, 0.98, 0.78, 0.00320},  // C60
    {29.7, 0.97, 0.77, 0.00315},  // C65
    {31.8, 0.96, 0.76, 0.00310},  // C70
    {33.8, 0.95, 0.75, 0.00305},  // C75
    {35.9, 0.94, 0.74, 0.00300},  // C80
};

// 在相邻等级之间线性内插
StressBlockParams interpolate(const GradePoint& lo, const GradePoint& hi,
                              double t) {
  StressBlockParams p;
  p.alpha1 = lo.alpha1 + (hi.alpha1 - lo.alpha1) * t;
  p.beta1 = lo.beta1 + (hi.beta1 - lo.beta1) * t;
  p.epsCu = lo.epsCu + (hi.epsCu - lo.epsCu) * t;
  return p;
}

// 由弯矩求解受压区相对高度 ξ = 1 - sqrt(1 - 2αs)
// 输入 αs = M/(α1·fc·b·h0²)；返回合法解或 -1（判别式非法 / 超筋不适用）
double xiFromAlphaS(double alphaS) {
  const double rad = 1.0 - 2.0 * alphaS;
  if (rad <= 0.0) return -1.0;  // 判别式非法：截面过小或弯矩过大
  return 1.0 - std::sqrt(rad);
}

// 有效高度
inline double effectiveDepth(const FlexuralInput& in) {
  return in.h - in.as;
}

}  // namespace

// -----------------------------------------------------------------------------
//  应力图形参数
// -----------------------------------------------------------------------------
StressBlockParams stressBlockFromFc(double fc) {
  const double fcMpa = fc / 1000.0;
  constexpr int kN = static_cast<int>(sizeof(kGrades) / sizeof(kGrades[0]));
  if (fcMpa <= kGrades[0].fcMpa) return interpolate(kGrades[0], kGrades[0], 0.0);
  if (fcMpa >= kGrades[kN - 1].fcMpa)
    return interpolate(kGrades[kN - 1], kGrades[kN - 1], 0.0);
  for (int i = 0; i + 1 < kN; ++i) {
    if (fcMpa <= kGrades[i + 1].fcMpa) {
      const double t = (fcMpa - kGrades[i].fcMpa) /
                       (kGrades[i + 1].fcMpa - kGrades[i].fcMpa);
      return interpolate(kGrades[i], kGrades[i + 1], t);
    }
  }
  return interpolate(kGrades[kN - 1], kGrades[kN - 1], 0.0);
}

double relativeBoundaryHeight(double fy, double Es,
                              const StressBlockParams& sb) {
  if (Es <= 0.0 || sb.epsCu <= 0.0) return 0.0;
  return sb.beta1 / (1.0 + fy / (Es * sb.epsCu));
}

// -----------------------------------------------------------------------------
//  最小配筋率（GB 50010 §8.5.1 表 8.5.1）
// -----------------------------------------------------------------------------
double minReinfRatio(double ft, double fy) {
  if (fy <= 0.0) return 0.0;
  return std::max(0.002, 0.45 * ft / fy);
}

// -----------------------------------------------------------------------------
//  单筋矩形截面受弯
// -----------------------------------------------------------------------------
FlexuralResult flexuralSingle(const FlexuralInput& in, double M) {
  FlexuralResult r;
  if (in.b <= 0.0 || in.h <= 0.0 || in.fc <= 0.0 || in.fy <= 0.0 ||
      in.as < 0.0 || M < 0.0) {
    r.status = FlexuralStatus::OverReinforced;  // 非法输入：按不可用处理
    return r;
  }
  const StressBlockParams sb = stressBlockFromFc(in.fc);
  r.xiB = relativeBoundaryHeight(in.fy, in.Es, sb);
  const double h0 = effectiveDepth(in);
  if (h0 <= 0.0) {
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }
  const double asMin = minReinfRatio(in.ft, in.fy) * in.b * in.h;
  r.AsMin = asMin;

  const double denom = sb.alpha1 * in.fc * in.b * h0 * h0;
  if (denom <= 0.0) {
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }
  const double alphaS = M / denom;
  const double xi = xiFromAlphaS(alphaS);
  if (xi < 0.0) {  // 2αs >= 1：弯矩超出截面上限（必然超筋）
    // 返回界限配筋供参考
    r.x = r.xiB * h0;
    r.xi = r.xiB;
    r.As = (sb.alpha1 * in.fc * in.b * r.x) / in.fy;
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }
  r.xi = xi;
  r.x = xi * h0;

  // 超筋判断：ξ > ξb 时单筋截面不足，须双筋或加大截面
  if (xi > r.xiB) {
    r.x = r.xiB * h0;
    r.xi = r.xiB;
    r.As = (sb.alpha1 * in.fc * in.b * r.x) / in.fy;
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }

  r.As = (sb.alpha1 * in.fc * in.b * r.x) / in.fy;
  if (r.As < asMin) {
    r.As = asMin;
    r.status = FlexuralStatus::BelowMinReinf;
  } else {
    r.status = FlexuralStatus::Ok;
  }
  return r;
}

// -----------------------------------------------------------------------------
//  双筋矩形截面受弯
// -----------------------------------------------------------------------------
FlexuralResult flexuralDouble(const FlexuralInput& in, double M,
                              double As_p) {
  FlexuralResult r;
  if (in.b <= 0.0 || in.h <= 0.0 || in.fc <= 0.0 || in.fy <= 0.0 ||
      in.as < 0.0 || in.as_p < 0.0 || M < 0.0 || As_p < 0.0) {
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }
  const StressBlockParams sb = stressBlockFromFc(in.fc);
  r.xiB = relativeBoundaryHeight(in.fy, in.Es, sb);
  const double h0 = effectiveDepth(in);
  if (h0 <= 0.0 || in.h - in.as_p <= 0.0) {
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }
  r.AsMin = minReinfRatio(in.ft, in.fy) * in.b * in.h;

  // §6.2.10-4：M1 = M - fy'·As'·(h0 - as')
  const double M1 = M - in.fy * As_p * (h0 - in.as_p);
  if (M1 < 0.0) {
    // 受压钢筋本身已足以抵抗弯矩：按 §6.2.14 最小配筋控制
    r.As = (M / (in.fy * (h0 - in.as_p)));
    if (r.As < r.AsMin) r.As = r.AsMin;
    r.status = r.As > r.AsMin + 1e-12 ? FlexuralStatus::Ok
                                      : FlexuralStatus::BelowMinReinf;
    return r;
  }

  const double denom = sb.alpha1 * in.fc * in.b * h0 * h0;
  if (denom <= 0.0) {
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }
  const double alphaS1 = M1 / denom;
  const double xi = xiFromAlphaS(alphaS1);
  if (xi < 0.0 || xi > r.xiB) {
    // 输入 As' 不足以把 ξ 压到界限内：超筋
    r.x = r.xiB * h0;
    r.xi = r.xiB;
    r.As = (sb.alpha1 * in.fc * in.b * r.x + in.fy * As_p) / in.fy;
    r.status = FlexuralStatus::OverReinforced;
    return r;
  }
  r.xi = xi;
  r.x = xi * h0;

  if (r.x < 2.0 * in.as_p) {
    // §6.2.14：受压区高度过小，按受压钢筋合力点取矩
    r.As = M / (in.fy * (h0 - in.as_p));
  } else {
    // §6.2.10-3：α1·fc·b·x + fy'·As' = fy·As
    r.As = (sb.alpha1 * in.fc * in.b * r.x + in.fy * As_p) / in.fy;
  }
  if (r.As < r.AsMin) {
    r.As = r.AsMin;
    r.status = FlexuralStatus::BelowMinReinf;
  } else {
    r.status = FlexuralStatus::Ok;
  }
  return r;
}

// -----------------------------------------------------------------------------
//  轴心受压（GB 50010 §6.2.15）
//   Nu = 0.9·φ·(fc·A + fy'·As')；φ 按 l0/b 查表 6.2.15 线性内插。
//   注：表值按矩形截面 b 边（短边）控制；圆形取 d。
// -----------------------------------------------------------------------------
namespace {

// 表 6.2.15：l0/b 与稳定系数 φ
constexpr double kPhiL0B[] = {0, 4, 6, 8, 10, 12, 14, 16, 18, 20,
                              22, 24, 26, 28, 30, 34, 38, 40, 42, 46, 48, 50};
constexpr double kPhiVal[] = {1.00, 1.00, 1.00, 1.00, 0.98, 0.95, 0.92, 0.87,
                              0.81, 0.75, 0.70, 0.65, 0.60, 0.56, 0.52, 0.44,
                              0.36, 0.32, 0.29, 0.26, 0.25, 0.24};
constexpr int kPhiN =
    static_cast<int>(sizeof(kPhiVal) / sizeof(kPhiVal[0]));

// 稳定系数 φ（l0/b 线性内插；超出上限取边缘值）
double stabilityFactor(double l0, double b) {
  if (b <= 0.0) return 1.0;
  const double rb = l0 / b;
  if (rb <= kPhiL0B[0]) return kPhiVal[0];
  if (rb >= kPhiL0B[kPhiN - 1]) return kPhiVal[kPhiN - 1];
  for (int i = 0; i + 1 < kPhiN; ++i) {
    if (rb <= kPhiL0B[i + 1]) {
      const double t = (rb - kPhiL0B[i]) / (kPhiL0B[i + 1] - kPhiL0B[i]);
      return kPhiVal[i] + (kPhiVal[i + 1] - kPhiVal[i]) * t;
    }
  }
  return kPhiVal[kPhiN - 1];
}

// 附加偏心距 ea = max(20mm, h/30)（GB 50010 §6.2.5）
inline double extraEccentricity(double h) {
  return std::max(0.02, h / 30.0);
}

}  // namespace

double axialCapacity(const ColumnInput& in, double l0, double As_p) {
  if (in.b <= 0.0 || in.h <= 0.0 || in.fc <= 0.0 || As_p < 0.0) return 0.0;
  const double phi = stabilityFactor(l0, in.b);
  const double A = in.b * in.h;
  return 0.9 * phi * (in.fc * A + in.fy * As_p);
}

// -----------------------------------------------------------------------------
//  矩形截面偏心受压（GB 50010 §6.2.17）
//
//  给定荷载（N, M）与配筋（As, As'），求该配筋在【同一偏心距 e 上】可承受
//  的轴力 Nu；Nu >= N 即满足。
//    e0 = M/N,  ea = max(20mm, h/30),  ei = e0 + ea
//    e  = ei + h/2 - as（轴力至受拉钢筋合力点距离，6.2.17-3）
//  记 A = α1·fc·b, C = fy'·As', Cp = fy·As, h0 = h - as, Hs = h - as'：
//
//  大偏心 (ξ<=ξb 且 x>=2as')：由 (6.2.17-1)(6.2.17-2) 消去 Nu：
//    A/2·x² + A·(e-h0)·x + [C·(e-h0+as') - Cp·e] = 0，取正根；
//    Nu = A·x + C - Cp。
//  x < 2as'（受压筋不屈服）：按规范取受压钢筋合力点矩心
//    e' = ei + h/2 - as'，Nu·e' = Cp·(h0 - as')。
//  小偏心 (ξ>ξb)：σs = fy·(ξ-β1)/(ξb-β1)（6.2.17-7），代入消去 Nu：
//    A/2·x² + [A·(e-h0) - e·Cp/denom]·x
//      + [C·(e - Hs) + e·Cp·β1/denom] = 0，denom = h0·(ξb-β1)，
//    取正根；Nu = A·x + C - σs·As。
// -----------------------------------------------------------------------------
ColumnResult columnCapacity(const ColumnInput& in, double N, double M,
                            double As, double As_p) {
  ColumnResult r;
  if (in.b <= 0.0 || in.h <= 0.0 || in.fc <= 0.0 || in.fy <= 0.0 ||
      N <= 0.0 || M < 0.0 || As < 0.0 || As_p < 0.0) {
    return r;  // 非法输入：Nu=0, satisfied=false
  }
  const StressBlockParams sb = stressBlockFromFc(in.fc);
  r.xiB = relativeBoundaryHeight(in.fy, in.Es, sb);
  const double h0 = in.h - in.as;
  const double Hs = in.h - in.as_p;
  if (h0 <= 0.0 || Hs <= 0.0) return r;

  const double e0 = M / N;
  const double ei = e0 + extraEccentricity(in.h);
  const double e = ei + in.h / 2.0 - in.as;

  const double A = sb.alpha1 * in.fc * in.b;  // α1·fc·b
  const double C = in.fy * As_p;              // fy'·As'
  const double Cp = in.fy * As;               // fy·As

  // ---- 大偏心试探解 ----
  const double ca = 0.5 * A;
  const double cb = A * (e - h0);
  const double cc = C * (e - h0 + in.as_p) - Cp * e;
  double x = 0.0;
  bool largeEcc = false;
  bool smallSolved = false;
  double disc = cb * cb - 4.0 * ca * cc;
  if (disc >= 0.0) {
    x = (-cb + std::sqrt(disc)) / (2.0 * ca);  // 正根
    if (x > 0.0 && x >= 2.0 * in.as_p && x / h0 <= r.xiB) largeEcc = true;
  }

  if (largeEcc) {
    // 大偏心：Nu 由轴力平衡 (6.2.17-1) 给出
    r.x = x;
    r.xi = x / h0;
    r.largeEccentric = true;
    r.Nu = A * x + C - Cp;
  } else if (disc >= 0.0 && x > 0.0 && x < 2.0 * in.as_p) {
    // x < 2as'：受压钢筋不屈服，对受压钢筋合力点取矩（§6.2.17 条文）
    const double eP = ei + in.h / 2.0 - in.as_p;
    if (eP > 0.0) {
      r.x = x;
      r.xi = x / h0;
      r.largeEccentric = true;
      r.Nu = Cp * (h0 - in.as_p) / eP;
    }
  } else {
    // ---- 小偏心：σs = fy·(ξ-β1)/(ξb-β1)，denom = h0·(ξb-β1) ----
    const double denom = (r.xiB - sb.beta1) * h0;
    if (std::abs(denom) > 1e-12 && std::abs(ca) > 1e-12) {
      const double a2 = ca;                                   // A/2
      const double b1 = A * (e - h0) - e * Cp / denom;
      const double c0 = C * (e - Hs) + e * Cp * sb.beta1 / denom;
      const double d2 = b1 * b1 - 4.0 * a2 * c0;
      if (d2 >= 0.0) {
        double xSmall = (-b1 + std::sqrt(d2)) / (2.0 * a2);
        if (xSmall <= 0.0) xSmall = (-b1 - std::sqrt(d2)) / (2.0 * a2);
        if (xSmall > 0.0) {
          x = xSmall;
          const double xi = x / h0;
          const double sigS = in.fy * (xi - sb.beta1) / (r.xiB - sb.beta1);
          r.x = x;
          r.xi = xi;
          r.largeEccentric = false;
          r.Nu = A * x + C - sigS * As;
          smallSolved = r.Nu > 0.0;
        }
      }
    }
    if (!smallSolved) {
      // 迭代解失败（数值退化）：回退到全截面受压控制（保守）
      r.x = in.h;
      r.xi = in.h / h0;
      r.largeEccentric = false;
      r.Nu = A * in.h + C - in.fy * As;  // σs = fy（ξ 很大时下界）
    }
  }
  r.satisfied = r.Nu >= N;
  return r;
}

// -----------------------------------------------------------------------------
//  受扭（GB 50010 §6.4）
// -----------------------------------------------------------------------------
namespace {

// 箍筋内表面所围核心混凝土尺寸（§6.4.3-2）
inline void coreRect(double b, double h, double c, double& bcor,
                     double& hcor, double& ucor, double& acor) {
  bcor = b - 2.0 * c;
  hcor = h - 2.0 * c;
  if (bcor <= 0.0 || hcor <= 0.0) {
    bcor = hcor = 0.0;
    ucor = acor = 0.0;
    return;
  }
  ucor = 2.0 * (bcor + hcor);
  acor = bcor * hcor;
}

// 受扭纵筋与箍筋配筋强度比 ζ（§6.4.4-3），夹在 0.6~1.7（§6.4.4-2）
double torsionZeta(double fy, double Astl, double s, double fyv, double Ast1,
                   double ucor) {
  if (fyv <= 0.0 || Ast1 <= 0.0 || ucor <= 0.0 || s <= 0.0) return 1.0;
  const double z = fy * Astl * s / (fyv * Ast1 * ucor);
  return std::max(0.6, std::min(1.7, z));
}

}  // namespace

double torsionPlasticModulus(double b, double h) {
  if (b <= 0.0 || h <= 0.0) return 0.0;
  return b * b * (3.0 * h - b) / 6.0;  // §6.4.3-1
}

TorsionResult pureTorsion(const TorsionInput& in, double Ast1, double Astl,
                          double s) {
  TorsionResult r;
  r.Wt = torsionPlasticModulus(in.b, in.h);
  if (r.Wt <= 0.0 || in.ft <= 0.0 || in.fyv <= 0.0 || Ast1 <= 0.0 ||
      s <= 0.0) {
    return r;
  }
  double bcor, hcor, ucor, acor;
  coreRect(in.b, in.h, in.c, bcor, hcor, ucor, acor);
  if (acor <= 0.0) return r;
  r.zeta = torsionZeta(in.fy, Astl, s, in.fyv, Ast1, ucor);
  // §6.4.4-1：Tu = 0.35·ft·Wt + 1.2·√ζ·fyv·Ast1·Acor/s
  r.Tu = 0.35 * in.ft * r.Wt +
         1.2 * std::sqrt(r.zeta) * in.fyv * Ast1 * acor / s;
  return r;
}

ShearTorsionResult shearTorsion(const TorsionInput& in, double V, double T,
                                double Ast1, double Astl, double Asv,
                                double s) {
  ShearTorsionResult r;
  r.Wt = torsionPlasticModulus(in.b, in.h);
  if (r.Wt <= 0.0 || in.ft <= 0.0 || in.fyv <= 0.0 || in.h0 <= 0.0 ||
      Ast1 <= 0.0 || s <= 0.0) {
    return r;
  }
  double bcor, hcor, ucor, acor;
  coreRect(in.b, in.h, in.c, bcor, hcor, ucor, acor);
  r.zeta = torsionZeta(in.fy, Astl, s, in.fyv, Ast1, ucor);
  // 混凝土受扭承载力降低系数 βt（§6.4.8-2），夹在 0.5~1.0
  if (T > 0.0 && V > 0.0) {
    const double beta =
        1.5 / (1.0 + 0.5 * V * r.Wt / (T * in.b * in.h0));
    r.betaT = std::max(0.5, std::min(1.0, beta));
  } else {
    r.betaT = 1.0;
  }
  // 受扭部分（§6.4.8-3）：Tu = 0.35·βt·ft·Wt + 1.2·√ζ·fyv·Ast1·Acor/s
  r.Tu = 0.35 * r.betaT * in.ft * r.Wt +
         1.2 * std::sqrt(r.zeta) * in.fyv * Ast1 * acor / s;
  // 受剪部分（§6.4.8-4）：Vu = (1.5-βt)·0.7·ft·b·h0 + fyv·Asv·h0/s
  r.Vu = (1.5 - r.betaT) * 0.7 * in.ft * in.b * in.h0 +
         in.fyv * Asv * in.h0 / s;
  return r;
}

}  // namespace check
}  // namespace yjk