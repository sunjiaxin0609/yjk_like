// =============================================================================
//  apps/ui/SceneData.cpp
// =============================================================================
#include "SceneData.h"

#include <algorithm>
#include <cmath>

#include "ColorMap.h"
#include "Theme.h"
#include "yjk/element/BeamElement3D.h"
#include "yjk/element/ShellElement4.h"

namespace ui {

namespace {

inline QVector3D toQ(const Vec3& v) {
  return QVector3D(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z));
}

// 从基向量与两个半尺寸生成一个长方体（8 点 → 6 面）。
//
// 【顶点顺序必须是"外法向右手序"】
// 渲染器靠叉积算面法线做光照，顺序反了法线就朝里，那一面会变成黑色 ——
// 表现为"实体模型有些面是黑的"，很容易被误认为光照参数不对。
void addBox(const QVector3D& c0, const QVector3D& c1,
            const QVector3D& u, const QVector3D& w,
            double hu, double hw,
            const QColor& fill, const QColor& edge,
            int elem, int elemType, double value, bool hasValue,
            std::vector<Face>& out) {
  if (hu <= 0.0 || hw <= 0.0) return;
  const QVector3D du = u * static_cast<float>(hu);
  const QVector3D dw = w * static_cast<float>(hw);

  // 两端各 4 个角点：(-u,-w), (+u,-w), (+u,+w), (-u,+w)
  QVector3D q[8];
  const QVector3D e[4] = {-du - dw, du - dw, du + dw, -du + dw};
  for (int i = 0; i < 4; ++i) {
    q[i] = c0 + e[i];
    q[i + 4] = c1 + e[i];
  }

  auto push = [&](int a, int b, int c, int d) {
    Face f;
    f.p[0] = q[a]; f.p[1] = q[b]; f.p[2] = q[c]; f.p[3] = q[d];
    f.fill = fill;
    f.edge = edge;
    f.elem = elem;
    f.elemType = elemType;
    f.value = value;
    f.hasValue = hasValue;
    out.push_back(f);
  };
  // i 端封头（法向 −axis）、j 端封头（法向 +axis）
  push(3, 2, 1, 0);
  push(4, 5, 6, 7);
  // 四个侧面
  push(0, 1, 5, 4);
  push(1, 2, 6, 5);
  push(2, 3, 7, 6);
  push(3, 0, 4, 7);
}

}  // namespace

const char* memberFieldName(MemberField f) {
  switch (f) {
    case MemberField::Moment:      return "弯矩 M";
    case MemberField::Shear:       return "剪力 V";
    case MemberField::Axial:       return "轴力 N";
    case MemberField::Stress:      return "正应力 σ";
    case MemberField::Utilization: return "应力比";
    default:                       return "";
  }
}

QVector3D SceneBuilder::nodePos(Id n, double scale) const {
  const Node& nd = m_.node(n);
  QVector3D p = toQ(nd.r);
  if (scale != 0.0 && res_.u.size() >= static_cast<size_t>(m_.nodeCount()) * 6) {
    const size_t b = static_cast<size_t>(nd.id * 6);
    p += QVector3D(static_cast<float>(res_.u[b + 0] * scale),
                   static_cast<float>(res_.u[b + 1] * scale),
                   static_cast<float>(res_.u[b + 2] * scale));
  }
  return p;
}

double SceneBuilder::autoDeformScale() const {
  if (res_.u.size() < static_cast<size_t>(m_.nodeCount()) * 6) return 1.0;
  double maxD = 0.0;
  for (Id i = 0; i < m_.nodeCount(); ++i) {
    const size_t b = static_cast<size_t>(i * 6);
    const double dx = res_.u[b], dy = res_.u[b + 1], dz = res_.u[b + 2];
    maxD = std::max(maxD, std::sqrt(dx * dx + dy * dy + dz * dz));
  }
  if (maxD <= 0.0) return 1.0;

  // 模型对角尺度
  QVector3D lo, hi;
  bool first = true;
  for (const Node& nd : m_.nodes()) {
    const QVector3D p = toQ(nd.r);
    if (first) { lo = hi = p; first = false; continue; }
    lo.setX(std::min(lo.x(), p.x())); lo.setY(std::min(lo.y(), p.y())); lo.setZ(std::min(lo.z(), p.z()));
    hi.setX(std::max(hi.x(), p.x())); hi.setY(std::max(hi.y(), p.y())); hi.setZ(std::max(hi.z(), p.z()));
  }
  const QVector3D d = hi - lo;
  const double diag = std::sqrt(static_cast<double>(d.x() * d.x() + d.y() * d.y() + d.z() * d.z()));
  if (diag <= 0.0) return 1.0;

  // 目标：最大位移占对角线的 4%。再大就会"看不出真实形状"，
  // 再小则肉眼分辨不出，这是 CAE 软件通行的比例。
  return 0.04 * diag / maxD;
}

double SceneBuilder::memberValue(const post::BeamResult& br, MemberField f) const {
  switch (f) {
    case MemberField::Moment: return br.Mmax;
    case MemberField::Shear:  return br.Vmax;
    case MemberField::Axial:  return std::max(std::abs(br.Nmax), std::abs(br.Nmin));
    case MemberField::Stress: return br.sigAbsMax;
    case MemberField::Utilization: return br.util;
    default: return 0.0;
  }
}

// -----------------------------------------------------------------------------
//  主构建
// -----------------------------------------------------------------------------
Scene SceneBuilder::build(const RenderOptions& opt) const {
  Scene s;

  const double scale = opt.deformed ? (opt.autoScale ? autoDeformScale() : opt.deformScale) : 0.0;
  s.appliedScale = scale;

  // 最大位移（状态栏）
  if (res_.u.size() >= static_cast<size_t>(m_.nodeCount()) * 6) {
    for (Id i = 0; i < m_.nodeCount(); ++i) {
      const size_t b = static_cast<size_t>(i * 6);
      s.maxDisp = std::max(s.maxDisp, std::sqrt(res_.u[b] * res_.u[b] +
                                                res_.u[b + 1] * res_.u[b + 1] +
                                                res_.u[b + 2] * res_.u[b + 2]));
    }
  }

  // ---- 云图取值 ----
  const std::vector<post::NodeField> kNoFields;
  const std::vector<post::BeamResult> kNoBeams;
  const std::vector<post::NodeField>& nf = pp_ ? pp_->nodeFields() : kNoFields;
  const std::vector<post::BeamResult>& brs = pp_ ? pp_->beamResults() : kNoBeams;
  const bool useNode = opt.contourOn && opt.nodeField >= 0 &&
                       opt.nodeField < static_cast<int>(nf.size());
  const bool useMember = opt.contourOn && !useNode && opt.memberField != MemberField::None;

  std::vector<double> shellVal;    // 逐壳单元
  std::vector<double> beamVal;     // 逐梁单元
  double vLo = 0.0, vHi = 0.0;

  if (useNode) {
    s.valueLabel = nf[static_cast<size_t>(opt.nodeField)].name;
    s.valueUnit = nf[static_cast<size_t>(opt.nodeField)].unit;

    // 【关键】范围必须包含 0。位移分量、应力都有正负，
    // 把范围钉死在 [dataMin, dataMax] 会让 0 落在色标中间某个不确定位置，
    // 用户没法判断"这个颜色到底是拉还是压"。包含 0 后，
    // 蓝→红的中点永远对应 0，正负一目了然。
    double a = 0.0, b = 0.0;
    bool first = true;
    for (double v : nf[static_cast<size_t>(opt.nodeField)].v) {
      if (first) { a = b = v; first = false; continue; }
      a = std::min(a, v); b = std::max(b, v);
    }
    vLo = std::min(a, 0.0);
    vHi = std::max(b, 0.0);
    if (vHi - vLo <= 0.0) { vLo -= 1.0; vHi += 1.0; }

    shellVal.assign(m_.elementCount(), 0.0);
    beamVal.assign(m_.elementCount(), 0.0);
    // 壳面：4 节点平均
    for (Id e = 0; e < m_.elementCount(); ++e) {
      const Element& el = *m_.elements()[static_cast<size_t>(e)];
      if (el.type() != ElementType::Shell4) continue;
      double sum = 0.0;
      for (Id n : el.nodes()) sum += nf[static_cast<size_t>(opt.nodeField)].v[static_cast<size_t>(n)];
      shellVal[static_cast<size_t>(e)] = sum / 4.0;
    }
    // 梁：两端平均（杆件很长时这是简化，构件级云图更精确）
    for (const post::BeamResult& br : brs) {
      const double a2 = nf[static_cast<size_t>(opt.nodeField)].v[static_cast<size_t>(br.ni)];
      const double b2 = nf[static_cast<size_t>(opt.nodeField)].v[static_cast<size_t>(br.nj)];
      beamVal[static_cast<size_t>(br.id)] = 0.5 * (a2 + b2);
    }
  } else if (useMember) {
    s.valueLabel = memberFieldName(opt.memberField);
    s.valueUnit = (opt.memberField == MemberField::Utilization) ? "" : "kN";
    if (opt.memberField == MemberField::Moment || opt.memberField == MemberField::Shear) s.valueUnit = "kN";
    if (opt.memberField == MemberField::Moment) s.valueUnit = "kN·m";
    if (opt.memberField == MemberField::Stress) s.valueUnit = "kPa";

    beamVal.assign(m_.elementCount(), 0.0);
    for (const post::BeamResult& br : brs)
      beamVal[static_cast<size_t>(br.id)] = memberValue(br, opt.memberField);

    // 构件级云图都是非负量（绝对值/应力比）⇒ 从 0 起。
    // 【从 0 起是必须的】若从 min 起，受力最小的构件会显示成红色，
    // 与"红色=危险"的直觉完全相反。
    vLo = 0.0;
    vHi = 0.0;
    for (double v : beamVal) vHi = std::max(vHi, v);
    if (vHi <= 0.0) vHi = 1.0;
    shellVal.assign(m_.elementCount(), 0.0);
  }

  s.vLo = vLo;
  s.vHi = vHi;

  auto colorOf = [&](double v, bool has) -> QColor {
    if (!has) return Member::concrete();
    return ColorMap::of(v, vLo, vHi);
  };

  // ---- 构件 ----
  const int hoverType = opt.hoverType;
  const int hoverElem = opt.hoverElem;
  const int extremeElem = opt.extremeElem;

  for (Id e = 0; e < m_.elementCount(); ++e) {
    const Element& el = *m_.elements()[static_cast<size_t>(e)];
    const bool hovered = (hoverElem == static_cast<int>(e)) && hoverType != 0;
    const bool isSel = std::find(opt.selected.begin(), opt.selected.end(), static_cast<int>(e)) != opt.selected.end();
    const bool isExtreme = (extremeElem == static_cast<int>(e)) && opt.extremeType != 0;

    // 交互色的优先级：悬停 > 选中 > 极值 > 云图/材质
    //
    // 【顺序不能颠倒】悬停是瞬时的，必须最醒目；极值标记是常驻的，
    // 让它压过悬停会导致"鼠标指到哪都看不出高亮"。
    auto finalColor = [&](QColor base) {
      if (hovered) return Member::hover();
      if (isSel) return Member::selected();
      if (isExtreme) return Member::highlight();
      return base;
    };

    if (el.type() == ElementType::Beam3D) {
      const auto& be = static_cast<const BeamElement&>(el);
      const QVector3D pa = nodePos(el.nodes()[0], scale);
      const QVector3D pb = nodePos(el.nodes()[1], scale);

      const bool hasVal = useNode || useMember;
      const double v = hasVal ? beamVal[static_cast<size_t>(e)] : 0.0;
      const QColor base = hasVal ? colorOf(v, true) : Member::concrete();
      const QColor col = finalColor(base);

      // 拾取用中心线（无论是否实体显示都要有）
      s.pickBeams.push_back({static_cast<int>(e), pa, pb});

      // 叠加高亮轮廓
      if (hovered || isSel || isExtreme) {
        const QColor hl = hovered ? Member::hover() : (isSel ? Member::selected() : Member::highlight());
        s.highlight.push_back({pa, pb, hl, hovered ? 5.0 : 4.0, static_cast<int>(e), 1});
      }

      if (opt.solid && be.core().isValid()) {
        const auto& sec = be.section();
        const double hy = sec.halfDepthY();
        const double hz = sec.halfWidthZ();
        // 截面尺寸为 0 ⇒ 退化成线，而不是画一个零厚度盒子（会 z-fighting）
        if (hy > 1e-6 && hz > 1e-6) {
          const auto& c = be.core();
          addBox(pa, pb, toQ(c.localY()), toQ(c.localZ()), hy, hz, col,
                 QColor(0, 0, 0, 0), static_cast<int>(e), 1, v, hasVal, s.faces);
        } else {
          s.lines.push_back({pa, pb, col, 2.0, static_cast<int>(e), 1});
        }
      } else {
        s.lines.push_back({pa, pb, col, hovered ? 3.5 : 2.0, static_cast<int>(e), 1});
      }

      // 端点必须计入包围盒，否则线框模式下节点在模型外时包围盒偏小
      s.addBounds(pa);
      s.addBounds(pb);
    }
  }

  for (Id e = 0; e < m_.elementCount(); ++e) {
    const Element& el = *m_.elements()[static_cast<size_t>(e)];
    if (el.type() != ElementType::Shell4) continue;
    const bool hovered = (hoverElem == static_cast<int>(e)) && hoverType == 2;
    const bool isSel = std::find(opt.selected.begin(), opt.selected.end(), static_cast<int>(e)) != opt.selected.end();
    const bool isExtreme = (extremeElem == static_cast<int>(e)) && opt.extremeType == 2;

// 壳单元取值：节点场用四节点平均；构件级云图（弯矩/剪力/轴力）只对梁有意义
    // —— 此时板/墙用中性色，不能跟着涂成最低档的蓝：那会被读成"板的弯矩为零"，
    // 而实际上板根本不在这个云图的统计范围内。
    const bool shellHasVal = useNode;
    const double v = shellHasVal ? shellVal[static_cast<size_t>(e)] : 0.0;
    // 未上云图时：竖墙用墙色（与楼板区分——墙是抗侧力构件，扫一眼视图
    // 就能看出核心筒/电梯井在哪），水平板用板色。
    const bool isWallElem = static_cast<const ShellElement&>(el).isWall();
    QColor base = shellHasVal ? ColorMap::of(v, vLo, vHi)
                              : (isWallElem ? Member::wall() : Member::slab());
    if (hovered) base = Member::hover();
    else if (isSel) base = Member::selected();
    else if (isExtreme) base = Member::highlight();

    Face f;
    for (int i = 0; i < 4; ++i) {
      f.p[i] = nodePos(el.nodes()[static_cast<size_t>(i)], scale);
      s.addBounds(f.p[i]);
    }
    f.fill = base;
    f.elem = static_cast<int>(e);
    f.elemType = 2;
    f.value = v;
    f.hasValue = shellHasVal;
    // 板面描边。不描边时整片楼板是一块纯色，看不出被剖分成几块；
    // 剖分信息对判断网格质量很重要，所以淡描边始终打开。
    f.edge = QColor(0, 0, 0, 90);

    if (hovered || isSel || isExtreme) {
      const QColor hl = hovered ? Member::hover() : (isSel ? Member::selected() : Member::highlight());
      for (int i = 0; i < 4; ++i)
        s.highlight.push_back({f.p[i], f.p[(i + 1) % 4], hl, 3.5, static_cast<int>(e), 2});
    }

    if (opt.solid) {
      s.faces.push_back(f);
    } else {
      // 线框模式：板只画四条边，不填色。
      // 【必须真的不填色】否则"线框"只是把梁换成线、板还是实心，
      // 内部的梁全被楼板挡住 —— 用户切到线框恰恰就是想看被挡住的构件。
      for (int i = 0; i < 4; ++i)
        s.lines.push_back({f.p[i], f.p[(i + 1) % 4], base, 1.2, static_cast<int>(e), 2});
    }
  }

  if (opt.showGrid) addGrid(opt, s);
  if (opt.showSupports) addSupports(opt, scale, s);
  if (opt.showNodes) addNodes(opt, scale, s);

  if (!s.boundsValid) { s.lo = QVector3D(0, 0, 0); s.hi = QVector3D(1, 1, 1); }
  return s;
}

// -----------------------------------------------------------------------------
//  节点
// -----------------------------------------------------------------------------
void SceneBuilder::addNodes(const RenderOptions& opt, double scale, Scene& s) const {
  const double sz = std::max(3.0, (s.hi - s.lo).length() * 0.5);
  for (const Node& nd : m_.nodes()) {
    if (nd.retired) continue;
    const bool hovered = (opt.hoverType == 3 && opt.hoverElem == static_cast<int>(nd.id));
    Point p;
    p.p = nodePos(nd.id, scale);
    p.color = hovered ? Member::hover() : Member::node();
    p.size = hovered ? 7.0 : 4.5;
    p.node = static_cast<int>(nd.id);
    s.points.push_back(p);
    (void)sz;
  }
}

// -----------------------------------------------------------------------------
//  支座与约束符号
// -----------------------------------------------------------------------------
void SceneBuilder::addConstraintGlyph(const Node& nd, const QVector3D& p, double sz, Scene& s) const {
  // 工程制图的约束符号：
  //   三向固定 → 一个立方体墩（简化为三向十字线 + 底部横线）
  //   竖向约束（滚动）→ 底下一个小三角
  // 这里不用图标贴图，直接用线段画，好处是随视角天然透视正确。
  const bool ux = nd.fixed[0], uy = nd.fixed[1], uz = nd.fixed[2];
  if (!ux && !uy && !uz) return;
  const QColor c = Member::support();

  // 底部水平十字（表示支承面）
  s.lines.push_back({p + QVector3D(-sz, 0, 0), p + QVector3D(sz, 0, 0), c, 2.0, -1, 0});
  s.lines.push_back({p + QVector3D(0, -sz, 0), p + QVector3D(0, sz, 0), c, 2.0, -1, 0});

  // 竖杆
  s.lines.push_back({p, p - QVector3D(0, 0, sz * 1.2f), c, 2.0, -1, 0});

  // 约束方向指示：哪个方向固定就画一根短粗线
  const double t = sz * 1.6;
  if (uz) s.lines.push_back({p - QVector3D(0, 0, sz * 1.2f), p - QVector3D(0, 0, sz * 1.2f) + QVector3D(t, 0, 0), c, 1.5, -1, 0});
  if (uz) s.lines.push_back({p - QVector3D(0, 0, sz * 1.2f), p - QVector3D(0, 0, sz * 1.2f) - QVector3D(t, 0, 0), c, 1.5, -1, 0});
  if (ux) s.lines.push_back({p, p + QVector3D(0, sz * 1.4f, 0), Member::support(), 1.5, -1, 0});
  if (uy) s.lines.push_back({p, p + QVector3D(sz * 1.4f, 0, 0), Member::support(), 1.5, -1, 0});
}

void SceneBuilder::addSupports(const RenderOptions& opt, double scale, Scene& s) const {
  const double sz = std::max(0.05, (s.hi - s.lo).length() * 0.012);
  for (const Node& nd : m_.nodes()) {
    if (nd.retired) continue;
    if (!nd.fixed[0] && !nd.fixed[1] && !nd.fixed[2]) continue;
    addConstraintGlyph(nd, nodePos(nd.id, scale), sz, s);
  }
  (void)opt;
}

// -----------------------------------------------------------------------------
//  地面参考网格
//
//  【网格必须画在结构底面标高，而不是 z = 0】
//  多层地下室或带高差的模型里，z=0 可能悬在半空。
//  画在底面标高，网格才真正起"地面"的参照作用。
// -----------------------------------------------------------------------------
void SceneBuilder::addGrid(const RenderOptions& opt, Scene& s) const {
  if (!s.boundsValid) return;
  const QVector3D c = (s.lo + s.hi) * 0.5f;
  const double dx = static_cast<double>(s.hi.x() - s.lo.x());
  const double dy = static_cast<double>(s.hi.y() - s.lo.y());
  const double span = std::max({dx, dy, 1.0}) * 1.6;
  const double step = ColorMap::niceStep(span, 10);
  const double z = s.lo.z();
  const double cx = static_cast<double>(c.x()), cy = static_cast<double>(c.y());
  const int n = static_cast<int>(std::ceil(span / step / 2.0)) + 1;

  // 每 5 格加粗一条，形成"主格/次格"。没有主次之分时，
  // 一片均匀的网格线读不出尺度 —— 数格子数不出米数。
  const double majorEvery = step * 5.0;
  for (int i = -n; i <= n; ++i) {
    const double x = cx + i * step;
    const double y = cy + i * step;
    const bool majX = (std::abs(std::remainder(x - cx, majorEvery)) < step * 0.25);
    const bool majY = (std::abs(std::remainder(y - cy, majorEvery)) < step * 0.25);
    s.lines.push_back({QVector3D(static_cast<float>(x), static_cast<float>(cy - n * step), static_cast<float>(z)),
                       QVector3D(static_cast<float>(x), static_cast<float>(cy + n * step), static_cast<float>(z)),
                       Viewport::gridMinor(), majX ? 1.4 : 0.8, -1, 0});
    s.lines.push_back({QVector3D(static_cast<float>(cx - n * step), static_cast<float>(y), static_cast<float>(z)),
                       QVector3D(static_cast<float>(cx + n * step), static_cast<float>(y), static_cast<float>(z)),
                       Viewport::gridMinor(), majY ? 1.4 : 0.8, -1, 0});
  }
  (void)opt;
}

}  // namespace ui
