// =============================================================================
//  apps/ui/ColorMap.h  ——  云图色标
//
//  【为什么要单独一个模块】
//  云图的"值 → 颜色"映射看起来只有一行代码，实际上有三个必须做对的点，
//  做错了图还能看，但会读出错误的结论：
//
//  ① 色标两端的值必须由【真实数据】决定，不能固定成 0 ~ max
//     受弯构件的应力有正有负。若把 vmin 钉死为 0，受压侧的负应力全部
//     落在色标最左端，整根构件一个颜色 —— 用户会以为"这根构件没受力"。
//     正确做法：范围取 [min(0, dataMin), max(0, dataMax)]，
//     这样 0 值永远有一个确定的位置，正负都能读出来。
//
//  ② 映射必须钳位。数值微超出范围（浮点误差、极值点）不该产生越界颜色
//     或数组越界。
//
//  ③ 离散色阶 vs 连续渐变：工程上读图习惯是【离散色阶】（YJK 用 10 档左右），
//     因为"哪一档"能直接对应一个数值区间。连续渐变好看但读不出量级。
//     这里两个都提供：实体云图用连续，杆件用离散档。
// =============================================================================
#pragma once

#include <QColor>
#include <QString>
#include <algorithm>
#include <cmath>

namespace ui {

class ColorMap {
 public:
  // 把 t ∈ [0,1] 映射到 rainbow 色标（蓝→青→绿→黄→红）。
  //
  // 分段线性插值，4 个区间 5 个锚点。不用 HSV 直接转，因为
  // HSV 的彩虹在青-绿段变化太快，视觉上会有一段"糊在一起"。
  static QColor rainbow(double t) {
    t = std::clamp(t, 0.0, 1.0);
    static const double pos[5] = {0.0, 0.25, 0.5, 0.75, 1.0};
    static const int rgb[5][3] = {
        {0x2B, 0x5C, 0xD8},   // 蓝
        {0x21, 0xC7, 0xC7},   // 青
        {0x4C, 0xD1, 0x37},   // 绿
        {0xF2, 0xE1, 0x2B},   // 黄
        {0xE0, 0x2E, 0x2E},   // 红
    };
    for (int i = 0; i < 4; ++i) {
      if (t <= pos[i + 1]) {
        const double u = (t - pos[i]) / (pos[i + 1] - pos[i]);
        return QColor(static_cast<int>(rgb[i][0] + u * (rgb[i + 1][0] - rgb[i][0]) + 0.5),
                      static_cast<int>(rgb[i][1] + u * (rgb[i + 1][1] - rgb[i][1]) + 0.5),
                      static_cast<int>(rgb[i][2] + u * (rgb[i + 1][2] - rgb[i][2]) + 0.5));
      }
    }
    return QColor(0xE0, 0x2E, 0x2E);
  }

  // 数值 → 颜色。range 的两个端点必须已经包含 0（见文件头 ①）。
  static QColor of(double v, double lo, double hi) {
    const double span = hi - lo;
    if (!(span > 0.0)) return rainbow(0.5);
    return rainbow((v - lo) / span);
  }

  // 用于图例：把 [lo,hi] 切成 n 档，返回第 i 档（0 = 最低）的颜色。
  static QColor band(double lo, double hi, int i, int n) {
    if (n <= 0) return rainbow(0.5);
    const double t = (i + 0.5) / n;
    (void)lo; (void)hi;
    return rainbow(t);
  }

  // 在已给区间内取一个"好看的"刻度步长：1/2/5 × 10^k。
  // 图例上出现 0.3333 这种刻度是不可接受的 —— 工程师没法心算。
  static double niceStep(double span, int wantBands) {
    if (!(span > 0.0) || wantBands <= 0) return 1.0;
    const double raw = span / wantBands;
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double n = raw / mag;
    double m = 1.0;
    if (n > 5.0)      m = 10.0;
    else if (n > 2.0) m = 5.0;
    else if (n > 1.0) m = 2.0;
    return m * mag;
  }

  // 带单位的数值格式化（云图图例用）。自动切换 m/k、kPa/MPa 这类量级。
  static QString label(double v) {
    const double a = std::abs(v);
    if (a == 0.0) return "0";
    if (a >= 1e5) return QString::number(v / 1e3, 'f', 0) + "k";
    if (a >= 1e3) return QString::number(v / 1e3, 'f', 1) + "k";
    if (a >= 100) return QString::number(v, 'f', 0);
    if (a >= 1)   return QString::number(v, 'f', 2);
    if (a >= 1e-3) return QString::number(v, 'f', 4);
    return QString::number(v, 'e', 1);
  }
};

}  // namespace ui
