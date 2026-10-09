// =============================================================================
//  apps/ui/Theme.h  ——  界面配色与视觉常量
//
//  【为什么把配色集中在一个头文件】
//  结构软件的界面配色不是"好看不好看"的问题，是有工程含义的：
//    · 混凝土构件与钢构件的默认色不同，扫一眼就能区分结构体系
//    · 云图色标必须与工程师看惯的 YJK/PKPM/ANSYS 一致（蓝→红 = 小→大）
//    · 选中色、悬停色要能从一堆构件里一眼挑出来
//  散落在各个 widget 里改一次就要满项目找，所以集中定义。
//
//  【深色视图 / 浅色窗体的理由】
//  三维视口用深色底不是审美偏好：构件用彩色云图上色时，
//  浅底会让黄色/青色失去对比度，深底才能让整个色标都读得出来。
//  这是 CAD/CAE 的通行做法（AutoCAD、ANSYS、YJK 的视口都是深色）。
//  菜单、对话框、表格保持系统浅色，保证与 Windows 其它软件一致。
// =============================================================================
#pragma once

#include <QColor>
#include <QFont>
#include <QString>

namespace ui {

// -----------------------------------------------------------------------------
//  三维视口
// -----------------------------------------------------------------------------
struct Viewport {
  // 背景：自上而下的渐变。纯色底会让空间感消失，看不出构件的前后关系。
  static QColor bgTop()    { return QColor(0x2B, 0x35, 0x43); }
  static QColor bgBottom() { return QColor(0x16, 0x1B, 0x22); }

  // 地面参考网格
  static QColor gridMinor() { return QColor(0x3A, 0x46, 0x55); }
  static QColor gridMajor() { return QColor(0x50, 0x60, 0x74); }

  // 坐标轴三色（X 红 / Y 绿 / Z 蓝，与工业界一致）
  static QColor axisX() { return QColor(0xE0, 0x4A, 0x4A); }
  static QColor axisY() { return QColor(0x4C, 0xAF, 0x50); }
  static QColor axisZ() { return QColor(0x42, 0x8B, 0xE0); }

  static QColor text()      { return QColor(0xD8, 0xDE, 0xE6); }
  static QColor textDim()   { return QColor(0x93, 0x9F, 0xAD); }
};

// -----------------------------------------------------------------------------
//  构件默认色（未上云图时）
//
//  按材料区分：混凝土偏灰暖，钢偏冷灰 —— 这是结构专业的习惯，
//  因为混凝土框架与钢框架的设计流程差别很大，能一眼分清很重要。
// -----------------------------------------------------------------------------
struct Member {
  static QColor concrete()  { return QColor(0xB8, 0xB2, 0xA6); }
  static QColor steel()     { return QColor(0x8E, 0x9A, 0xA8); }
  static QColor slab()      { return QColor(0x6E, 0x7E, 0x92); }
  static QColor wall()      { return QColor(0x7A, 0x86, 0x78); }

  // 交互反馈
  //
  // 【选中色为什么用品红而不是琥珀】
  // 琥珀 #FFC107 与云图色标里的"黄"几乎同色，选中一根黄色构件完全看不出。
  // 品红 #FF2FD0 不在 rainbow 的任何一段上（rainbow 走 蓝→青→绿→黄→红，
  // 没有紫/品红），因此在任何档位上都有对比度。
  static QColor hover()     { return QColor(0x4F, 0xC3, 0xF7); }   // 青：#0 在 rainbow 里但亮度高、且只在悬停瞬间出现
  static QColor selected()  { return QColor(0xFF, 0x2F, 0xD0); }   // 品红：rainbow 之外
  static QColor highlight() { return QColor(0xFF, 0x70, 0x43); }   // 橙红：极值定位
  static QColor node()      { return QColor(0x9A, 0xA8, 0xB8); }
  static QColor support()   { return QColor(0xFF, 0xD5, 0x4F); }
};

// -----------------------------------------------------------------------------
//  云图色标
//
//  与 YJK / PKPM / ANSYS 的 "rainbow" 色标一致：蓝→青→绿→黄→红。
//  【不要换成 viridis 之类的现代配色】—— 工程师看惯了 rainbow，
//  换色会让"红色=危险"这个直觉失效，是可用性倒退。
// -----------------------------------------------------------------------------
struct Contour {
  // 5 个锚点（蓝 - 青 - 绿 - 黄 - 红）
  static QColor c0() { return QColor(0x2B, 0x5C, 0xD8); }
  static QColor c1() { return QColor(0x21, 0xC7, 0xC7); }
  static QColor c2() { return QColor(0x4C, 0xD1, 0x37); }
  static QColor c3() { return QColor(0xF2, 0xE1, 0x2B); }
  static QColor c4() { return QColor(0xE0, 0x2E, 0x2E); }
};

// -----------------------------------------------------------------------------
//  结果状态色（验算指标用）
// -----------------------------------------------------------------------------
struct Status {
  static QColor ok()   { return QColor(0x2E, 0x7D, 0x32); }   // 绿：满足
  static QColor warn() { return QColor(0xEF, 0x6C, 0x00); }   // 橙：接近限值
  static QColor fail() { return QColor(0xC6, 0x28, 0x28); }   // 红：超限
  static QColor info() { return QColor(0x15, 0x65, 0xC0); }
};

// -----------------------------------------------------------------------------
//  字号
// -----------------------------------------------------------------------------
inline int fsSmall()  { return 11; }
inline int fsNormal() { return 12; }
inline int fsTitle()  { return 13; }

// 等宽字体 —— 数值表格必须等宽，否则小数点对不齐，
// 一列数字扫下来看不出量级差异。这不是审美问题。
QFont monoFont(int pt = 9);
QFont uiFont(int pt = 10);

// 统一的数值格式化：结构工程里数量级跨度极大（1e-6 m ~ 1e5 kN），
// 固定 4 位小数会把小量显示成 0.0000，把大量显示成 123456.7890 不好读。
QString fmtNum(double v, int sig = 4);

}  // namespace ui
