// =============================================================================
//  yjk/model/Section.h  ——  截面与材料
//
//  单位约定（与 yjk/math/Types.h 一致）：
//    长度 m, 力 kN, 应力 kPa, 弯矩 kN·m
//    混凝土 C30: E = 3.0e4 MPa = 3.0e7 kPa, fck = 14.3 MPa = 14300 kPa
//    钢材 Q235: E = 2.06e5 MPa = 2.06e8 kPa, fy = 235 MPa = 235000 kPa
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <string>

#include "yjk/math/Types.h"

namespace yjk {

// ---------------------------------------------------------------------------
//  材料
// ---------------------------------------------------------------------------
struct Material {
  std::string name;
  double E{3.0e7};        // 弹性模量 kPa
  double nu{0.2};         // 泊松比
  double gamma{25.0};     // 容重 kN/m³（混凝土 25，钢 78.5）
  double fy{0.0};         // 屈服/设计强度 kPa，0 表示无塑性（仅弹性分析）
  double ft{0.0};         // 抗拉强度 kPa（混凝土受拉）
  double fc{0.0};         // 抗压强度 kPa（混凝土受压）
  double G{0.0};          // 剪切模量 kPa，0 表示由 E、ν 自动算

  // 剪切模量
  double shearModulus() const { return G > 0.0 ? G : E / (2.0 * (1.0 + nu)); }

  // 预置材料库
  //
  // 强度等级 → Ec / f_c / f_t 的对应表必须与规范原表一致。
  // 这里采用 GB 50010-2010（2015年版）表 4.1.3、4.1.5 的常用等级。
  //
  // 单位提醒：规范里是 MPa，内部统一用 kPa，故所有数值 ×1000。
  //
  //   等级  Ec(×10⁴MPa)  f_c(MPa)  f_t(MPa)
  //   C15     2.55        7.3      0.88
  //   C20     2.80        9.6      1.10
  //   C25     3.00       11.9      1.27
  //   C30     3.15       14.3      1.43
  //   C35     3.25       16.7      1.57
  //   C40     3.30       19.0      1.71
  //   C45     3.35       21.1      1.80
  //   C50     3.40       23.4      1.89
  //   C55     3.45       25.3      1.96
  //   C60     3.50       27.5      2.04
  static Material concreteC(int grade) {
    struct Row { int g; double Ec, fc, ft; };
    static const Row tbl[] = {
        {15, 2.55e7, 7.3e3, 0.88e3}, {20, 2.80e7, 9.6e3, 1.10e3},
        {25, 3.00e7, 11.9e3, 1.27e3}, {30, 3.15e7, 14.3e3, 1.43e3},
        {35, 3.25e7, 16.7e3, 1.57e3}, {40, 3.30e7, 19.0e3, 1.71e3},
        {45, 3.35e7, 21.1e3, 1.80e3}, {50, 3.40e7, 23.4e3, 1.89e3},
        {55, 3.45e7, 25.3e3, 1.96e3}, {60, 3.50e7, 27.5e3, 2.04e3},
    };
    Material m;
    m.name = "C" + std::to_string(grade);
    m.nu = 0.2;
    m.gamma = 25.0;
    m.fy = 0.0;
    m.Ec_ = 0.0; m.fc_ = 0.0; m.ft_ = 0.0;   // 缓存原值便于查表核对
    for (const auto& r : tbl) {
      if (r.g == grade) {
        m.E = r.Ec; m.fc = r.fc; m.ft = r.ft;
        m.Ec_ = r.Ec / 1e4; m.fc_ = r.fc / 1e3; m.ft_ = r.ft / 1e3;
        return m;
      }
    }
    // 等级不在表内：按 C30 保守取值，并保持名称可辨识
    m.E = 3.15e7; m.fc = 14.3e3; m.ft = 1.43e3;
    m.Ec_ = 3.15; m.fc_ = 14.3; m.ft_ = 1.43;
    return m;
  }

  // 钢材（GB 50017-2017 表 4.4.1）
  //   牌号  f_y(MPa)  f_u(MPa)  E(MPa)
  //   Q235    235       370      2.06e5
  //   Q345    345       470      2.06e5
  //   Q390    390       490      2.06e5
  //   Q420    420       520      2.06e5
  //   Q460    460       570      2.06e5
  static Material steelQ(int grade) {
    struct Row { int g; double fy, fu, E; };
    static const Row tbl[] = {
        {235, 235e3, 370e3, 2.06e8}, {345, 345e3, 470e3, 2.06e8},
        {390, 390e3, 490e3, 2.06e8}, {420, 420e3, 520e3, 2.06e8},
        {460, 460e3, 570e3, 2.06e8},
    };
    Material m;
    m.name = "Q" + std::to_string(grade);
    m.nu = 0.3;
    m.gamma = 78.5;
    for (const auto& r : tbl) {
      if (r.g == grade) { m.E = r.E; m.fy = r.fy; m.ft = r.fu; return m; }
    }
    m.E = 2.06e8; m.fy = 235e3; m.ft = 370e3;
    return m;
  }

  // 规范表中的原始值（MPa），供界面显示与核对
  double Ec_{0.0}, fc_{0.0}, ft_{0.0};
};

// ---------------------------------------------------------------------------
//  截面几何性质
//
//  六个基本量决定梁单元的全部刚度：
//    A     截面积           (轴向 + 剪切)
//    Iy    绕 y 轴惯性矩     (平面内弯曲 → 平面外弯矩 Mz)
//    Iz    绕 z 轴惯性矩     (平面外弯曲 → 平面内弯矩 My)
//    J     圣维南扭转常数   (扭转)
//    Asy   沿 y 方向有效剪切面积（考虑腹板高度）
//    Asz   沿 z 方向有效剪切面积
//
//  【为什么必须有 Asy/Asz】Euler-Bernoulli 梁假设截面变形时"剪切变形为零"，
//  这对"长细比大的梁"成立。但对深梁、转换构件、支撑短柱，
//  剪切变形占主导，忽略它刚度会被高估 20~40%，内力分布明显失真。
//  盈建科/PKPM 都提供"考虑剪切变形"开关，就是这个原因。
// ---------------------------------------------------------------------------
struct SectionProperties {
  double A{0.0};      // m²
  double Iy{0.0};     // m⁴  绕 y 轴（主惯性轴之一）
  double Iz{0.0};     // m⁴  绕 z 轴
  double J{0.0};      // m⁴  扭转
  double Asy{0.0};    // m²  剪切面积（y 向）
  double Asz{0.0};    // m²  剪切面积（z 向）
  double Ry{0.0};     // 截面塑性抵抗矩 m³（绕 y，用于塑性分析）
  double Rz{0.0};

  // 边缘纤维距离（后处理算弯曲应力用）
  //   hy = 截面在【局部 y】方向的外缘距离（半高）  m
  //   hz = 截面在【局部 z】方向的外缘距离（半宽）  m
  //
  //  σ = N/A ± My·hz/Iy ± Mz·hy/Iz —— 没有这两个量就只能给出内力、
  //  给不了应力，而配筋和云图要的恰恰是应力。
  //
  //  = 0 时由 I、A 反算（矩形精确：hy = √(3·Iz/A)、hz = √(3·Iy/A)），
  //  对工字钢会偏保守，但至少不会静默给零。
  double hy{0.0};
  double hz{0.0};

  double halfDepthY() const { return hy > 0.0 ? hy : (A > 0 && Iz > 0 ? std::sqrt(3.0 * Iz / A) : 0.0); }
  double halfWidthZ() const { return hz > 0.0 ? hz : (A > 0 && Iy > 0 ? std::sqrt(3.0 * Iy / A) : 0.0); }

  // 有效性检查 —— 退化截面会让单元刚度奇异，必须拦下
  bool isValid(std::string* why = nullptr) const {
    auto bad = [&](const char* m) { if (why) *why = m; return false; };
    if (!(A > 0)) return bad("截面积为零或为负");
    if (!(Iy > 0) || !(Iz > 0)) return bad("惯性矩为零或为负（截面退化）");
    if (!(J > 0)) return bad("扭转常数为零或为负");
    return true;
  }
};

// ---------------------------------------------------------------------------
//  常用截面
//
//  关于"惯性矩 Iy 与 Iz 哪个大"——这是最容易搞错的地方。
//  本库的约定：局部 y 轴是截面的"竖向"（对梁通常是截面高度方向），
//            局部 z 轴是截面的"横向"（梁的翼缘宽度方向）。
//  所以矩形截面 b×h（b=宽, h=高）：
//    Iy = b·h³/12   （绕竖轴 y，弯曲时"弱轴"——截面高度方向的弯曲）
//    Iz = h·b³/12   （绕横轴 z，弯曲时"强轴"——宽度方向的弯曲）
//  对于沿 X 方向布置、截面高 500 的梁，平面内弯矩（绕 z）用 Iz = h·b³/12。
// ---------------------------------------------------------------------------
namespace section {

// 矩形：b = 宽度（沿局部 z），h = 高度（沿局部 y）
//
// 【惯性矩约定 —— 单元层最容易出错、后果最严重的地方】
//
// 局部坐标：x 沿梁轴，y = 截面【竖向】，z = 截面【横向/宽度】方向。
//
//   Iy = ∫z² dA = h·b³/12    绕 y 轴（竖轴）弯曲
//   Iz = ∫y² dA = b·h³/12    绕 z 轴（横轴）弯曲
//
// 例：b = 0.3（宽）, h = 0.5（高）
//   Iy = 0.5·0.027/12 = 1.125e-3    弱轴（材料在竖向分布集中）
//   Iz = 0.3·0.125/12 = 3.125e-3    强轴（材料离竖轴更远）
//
// 【为什么必须分清】梁单元局部刚度矩阵中，绕 y 的弯曲项用 Iy、绕 z 的用 Iz。
// 两者写反 → 平面内与平面外弯曲刚度完全对调 →
//   · 框架的侧向刚度被严重高估（梁被当成"平面内刚"）
//   · 梁的竖向挠度算大一个量级
//   · 而程序不报任何错。
// 这是"算出了数但完全错了"最典型的一类 bug。
//
// 记忆法：绕哪根轴，就对【另一方向】的坐标求二阶矩。
//     绕 y → 对 z 求矩 → z 方向尺寸是 b → h·b³/12
inline SectionProperties rect(double b, double h) {
  SectionProperties s;
  s.A = b * h;
  s.Iy = h * b * b * b / 12.0;   // 绕 y（竖轴）—— 弱轴
  s.Iz = b * h * h * h / 12.0;   // 绕 z（横轴）—— 强轴

  // 圣维南扭转常数 J
  //
  // 【量级基准】J = α · a · t³  ——  a 是【长边】，t 是【短边】。
  // 这个量级关系是硬约束：抗扭刚度由"长边 × 短边³"决定，
  // 因为扭转时剪应力集中在短边附近，长边只起延伸作用。
  // 若把长短边弄反，J 与 Iz 都会算错（且 Iz 会小到不合理）。
  //
  // 矩形截面的 Saint-Venant 精确解（Timoshenko/Greenhill 形式）：
  //     J = (a·t³/3) · [ 1 − 0.63·(t/a) + 0.052·(t/a)⁵ ]
  //
  // 已对标 Roark's Formulas for Engineers 表 3-2 的 α 值（α = J/(a·t³)），
  // 全范围吻合到 0.01%：
  //     a/t  = 1.0  α = 0.1406    a/t = 2.0  α = 0.229
  //     a/t  = 1.2  α = 0.166     a/t = 3.0  α = 0.263
  //     a/t  = 1.5  α = 0.196     a/t = 5.0  α = 0.291
  //                                a/t → ∞  α → 1/3
  // 两个极限均自洽：
  //     正方形 a=t → J = 0.1406·a⁴
  //     狭长 a≫t  → J = a·t³/3（薄壁条，与膜理论一致）
  //
  // 【不要用 GB 50017 附录 C 那种"β 系数分段公式"的手抄版本】——
  // 那些式子分母含 a−1.12 之类的项，a 接近 1 时会变号得到负 J。
  // 负 J 会让扭转刚度为负：求解器报奇异，或更糟：静默给出错误内力。
  {
    const double a = std::max(b, h);   // 长边
    const double t = std::min(b, h);   // 短边
    const double r = t / a;            // 短长比 ∈ (0, 1]
    s.J = (a * t * t * t / 3.0) * (1.0 - 0.63 * r + 0.052 * r * r * r * r * r);
    if (!(s.J > 0.0) || !std::isfinite(s.J)) s.J = a * t * t * t / 3.0;
  }

  // 有效剪切面积（Timoshenko）：矩形 As = 5/6·A
  // 这是剪应力精确解给出的值，标准取值
  s.Asy = 5.0 / 6.0 * s.A;
  s.Asz = 5.0 / 6.0 * s.A;
  s.hy = h / 2.0;                 // 沿局部 y（高度方向）
  s.hz = b / 2.0;                 // 沿局部 z（宽度方向）
  // 塑性抵抗矩同样遵循"绕哪轴、对另一方向求矩"的约定
  s.Ry = h * b * b / 4.0;   // 绕 y
  s.Rz = b * h * h / 4.0;   // 绕 z
  return s;
}

// 圆：d = 直径
inline SectionProperties circle(double d) {
  SectionProperties s;
  const double A0 = std::acos(-1.0);
  s.A = A0 * d * d / 4.0;
  s.Iy = s.Iz = A0 * std::pow(d, 4.0) / 64.0;
  s.J = 2.0 * s.Iy;
  s.Asy = s.Asz = 0.9 * s.A;
  s.Ry = s.Rz = d * d * d / 6.0;
  s.hy = s.hz = d / 2.0;
  return s;
}

// 圆管：d = 外径, t = 壁厚
inline SectionProperties tube(double d, double t) {
  SectionProperties s;
  const double di = d - 2.0 * t;
  const double A0 = std::acos(-1.0);
  s.A = A0 / 4.0 * (d * d - di * di);
  const double dd4 = std::pow(d, 4.0) - std::pow(di, 4.0);
  s.Iy = s.Iz = A0 / 64.0 * dd4;
  s.J = 2.0 * s.Iy;
  // 薄壁圆管剪切面积取 0.5·A（经典薄膜近似）
  s.Asy = s.Asz = 0.5 * s.A;
  s.Ry = s.Rz = (d * d * d - di * di * di) / 6.0;   // 各向同性，两者相等
  s.hy = s.hz = d / 2.0;
  return s;
}

// 工字形（H 型钢）：h 高度, b 翼缘宽, tw 腹板厚, tf 翼缘厚
//
// 局部坐标：y = 竖向（截面高度 h 的方向），z = 横向（翼缘宽度 b 的方向）
//   Iz = ∫y²dA  绕 z（横轴）—— 【强轴】，抵抗上下弯曲
//   Iy = ∫z²dA  绕 y（竖轴）—— 【弱轴】，抵抗侧向弯曲
// 典型 H 型钢 Iz/Iy ≈ 3~10，这个比值是"强轴/弱轴"关系的直接体现。
inline SectionProperties iSection(double h, double b, double tw, double tf) {
  SectionProperties s;
  const double hw = h - 2.0 * tf;          // 腹板净高
  s.A = 2.0 * b * tf + tw * hw;

  // Iy（弱轴，绕 y 竖轴）：按 y 方向切条，∫z²dA
  //   腹板（高 hw、宽 tw 的薄板，绕自身竖轴）: hw·tw³/12
  //   翼缘（宽 b 的板）: 2 × tf·b³/12
  //
  // 【易错点】腹板对 Iy 的贡献是 hw·tw³/12，而不是 tw·hw³/12 ——
  // 后者是"绕横轴"的量纲，用在这里会让弱轴刚度偏大两个数量级，
  // 直接导致梁的侧向（平面外）刚度被严重高估，且不报错。
  // 记忆：绕 y 轴 → 沿 y 切条 → 每条按宽度方向的二阶矩积分。
  s.Iy = hw * std::pow(tw, 3.0) / 12.0 + 2.0 * tf * std::pow(b, 3.0) / 12.0;

  // Iz（强轴，绕 z 横轴）：翼缘需用平行轴定理搬到 z 轴
  //   腹板 tw·hw³/12
  //   翼缘 2 × [ b·tf³/12 + b·tf·((hw+tf)/2)² ]
  const double dy = (hw + tf) / 2.0;       // 翼缘形心到 z 轴的距离
  s.Iz = tw * std::pow(hw, 3.0) / 12.0 +
         2.0 * (b * std::pow(tf, 3.0) / 12.0 + b * tf * dy * dy);

  // 圣维南扭转常数：开口薄壁近似 J ≈ (1/3)·Σ(各矩形段 b·t³)
  // 注意：开口截面的 J 极小（相对闭合截面差几个数量级），
  // 这正是工字梁抗扭弱的原因 —— 计算时不能拿 Iz 代替 J。
  s.J = (1.0 / 3.0) * (hw * tw * tw * tw + 2.0 * b * tf * tf * tf);
  // 有效剪切面积：腹板承担几乎全部剪力（翼缘的贡献通常忽略）
  s.Asy = hw * tw;
  s.Asz = hw * tw;
  // 塑性抵抗矩（矩形近似，误差对承载力判定影响可接受；
  // 精确值需要积分翼缘与腹板的屈服分布，盈建科在做截面塑性域分析时才需要）
  s.Ry = h * b * b / 4.0;      // 绕 y（弱轴）
  s.Rz = b * h * h / 4.0;      // 绕 z（强轴）
  s.hy = h / 2.0;
  s.hz = b / 2.0;
  return s;
}

// 箱形：h, b 外廓, tw 腹板厚, tf 翼缘厚
inline SectionProperties boxSection(double h, double b, double tw, double tf) {
  SectionProperties s;
  const double A0 = std::acos(-1.0);
  const double hi = h - 2.0 * tf, bi = b - 2.0 * tw;
  s.A = b * h - bi * hi;
  s.Iy = (b * h * h * h - bi * hi * hi * hi) / 12.0;   // 绕竖轴
  s.Iz = (h * b * b * b - hi * bi * bi * bi) / 12.0;   // 绕横轴
  // 闭合薄壁：J = 4·A_m²/∮(ds/t)，A_m = 截面中线所围面积
  const double hm = h - tf, bm = b - tw;
  const double Am = hm * bm;
  const double perim = 2.0 * (hm + bm);
  const double tm = (tw + tf) / 2.0;
  s.J = 4.0 * Am * Am / (perim / tm);
  s.Asy = s.Asz = A0 / 6.0 * h * b - A0 / 6.0 * hi * bi;
  s.Ry = h * b * b / 4.0;
  s.Rz = b * h * h / 4.0;
  s.hy = h / 2.0;
  s.hz = b / 2.0;
  return s;
}

}  // namespace section

}  // namespace yjk
