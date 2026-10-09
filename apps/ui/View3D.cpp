// =============================================================================
//  apps/ui/View3D.cpp  ——  三维视口实现
//
//  渲染流程（每帧）：
//    ① 背景渐变
//    ② 构建绘制项：投影所有面/线/点到屏幕，同时算视空间深度
//    ③ 按深度排序（远 → 近），画家算法逐项绘制
//    ④ 叠加层：色标、坐标轴三角、信息文字
//
//  【为什么是整个场景排序，而不是分开画面和线】
//  分开画会让地面参考网格"透过"实体模型显示出来，空间关系立刻错乱。
//  这是软件渲染里最经典的一个错误，而且看起来像"模型是透明的"，
//  很容易被误判成透明度参数问题。
// =============================================================================
#include "View3D.h"

#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QPixmap>
#include <QPolygonF>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>

#include "ColorMap.h"
#include "Theme.h"
#include "yjk/interact/PickCore.h"

namespace ui {

namespace {

inline QVector3D toQ(const Vec3& v) {
  return QVector3D(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z));
}

// 固定的双光源方向（世界坐标）。
// 【为什么不用"跟随相机的头灯"】头灯会让所有朝向相机的面亮度一致，
// 体块之间的转折完全消失，实体模型看起来像一片平面色块。
// 固定光源才能保持"光从斜上方来"的空间感。
const QVector3D kLight1 = QVector3D(0.42f, -0.55f, 0.72f).normalized();
const QVector3D kLight2 = QVector3D(-0.55f, 0.45f, 0.30f).normalized();

QColor shadeColor(const QColor& base, const QVector3D& n) {
  const double d1 = std::max(0.0f, QVector3D::dotProduct(n, kLight1));
  const double d2 = std::max(0.0f, QVector3D::dotProduct(n, kLight2));
  // 环境光取 0.42：再低背光面就全黑了，看不到截面轮廓。
  double k = 0.42 + 0.52 * d1 + 0.20 * d2;
  k = std::clamp(k, 0.30, 1.30);
  return QColor(std::min(255, static_cast<int>(base.red() * k + 0.5)),
                std::min(255, static_cast<int>(base.green() * k + 0.5)),
                std::min(255, static_cast<int>(base.blue() * k + 0.5)),
                base.alpha());
}

// 像素容差(px) → 世界容差(m)。与拾取共用：视口中心向上 1 px 的
// 世界位移即"每像素对应多少米"，对透视/正交投影都正确。
double worldPerPixAt(const Camera& cam, const QRect& vp) {
  const QVector3D c0 = cam.unproject(0, 0, vp);
  const QVector3D c1 = cam.unproject(0, 1, vp);
  return std::max(static_cast<double>((c1 - c0).length()), 1e-6);
}

struct DrawItem {
  double depth{0.0};          // 视空间 z（越负越远）
  int kind{0};                // 0 = 面，1 = 线，2 = 点
  int index{0};
  QPolygonF screen;
  QPointF a, b;
  QColor color;
  double width{1.0};
  double size{4.0};
};

}  // namespace

// =============================================================================
//  构造 / 数据源
// =============================================================================
View3D::View3D(QWidget* parent) : QWidget(parent) {
  setFocusPolicy(Qt::StrongFocus);
  setMouseTracking(true);
  setMinimumSize(400, 300);
  setAttribute(Qt::WA_OpaquePaintEvent);

  animTimer_.setInterval(40);   // 25 fps 对变形动画足够，再高只是浪费 CPU
  connect(&animTimer_, &QTimer::timeout, this, &View3D::onAnimTick);
}

void View3D::setModel(const Model* m) {
  model_ = m;
  result_ = nullptr;
  post_ = nullptr;
  sel_.clear();
  selType_ = 0;
  opt_.hoverElem = -1;
  opt_.hoverType = 0;
  opt_.selected.clear();
  opt_.extremeElem = -1;
  opt_.extremeType = 0;
  opt_.contourOn = false;
  opt_.nodeField = -1;
  opt_.memberField = MemberField::None;
  opt_.deformed = false;
  dirty_ = true;
  update();
}

void View3D::setResult(const StaticResult* r) {
  result_ = r;
  dirty_ = true;
  update();
}

void View3D::setPost(const post::PostProcessor* pp) {
  post_ = pp;
  dirty_ = true;
  update();
}

void View3D::clearPost() {
  post_ = nullptr;
  result_ = nullptr;
  opt_.contourOn = false;
  opt_.nodeField = -1;
  opt_.memberField = MemberField::None;
  opt_.deformed = false;
  opt_.extremeElem = -1;
  opt_.extremeType = 0;
  dirty_ = true;
  update();
}

void View3D::setOptions(const RenderOptions& o) {
  if (o.deformed && !opt_.deformed) animPhase_ = M_PI;   // 刚打开时立刻满变形
  if (!o.deformed) animPhase_ = 0.0;
  opt_ = o;
  dirty_ = true;
  update();
}

void View3D::refresh() {
  dirty_ = true;
  update();
}

void View3D::setHighlightExtreme(int elem, int type) {
  opt_.extremeElem = elem;
  opt_.extremeType = type;
  dirty_ = true;
  update();
}

void View3D::clearSelection() {
  sel_.clear();
  selType_ = 0;
  opt_.selected.clear();
  selNode_ = -1;
  opt_.selectedNodes.clear();
  dirty_ = true;
  update();
  emit selectionChanged();
}

void View3D::selectElement(int elem, int type) {
  sel_.clear();
  opt_.selected.clear();
  selType_ = 0;
  selNode_ = -1;
  opt_.selectedNodes.clear();
  if (elem >= 0 && type > 0) {
    sel_.push_back(elem);
    opt_.selected.push_back(elem);
    selType_ = type;
  }
  dirty_ = true;
  update();
  emit selectionChanged();
}

void View3D::setAnimating(bool on) {
  if (animating_ == on) return;
  animating_ = on;
  if (on) {
    opt_.deformed = true;
    animPhase_ = M_PI;
    animTimer_.start();
  } else {
    animTimer_.stop();
    animPhase_ = M_PI;
  }
  dirty_ = true;
  update();
  emit viewChanged();
}

void View3D::onAnimTick() {
  animPhase_ += 0.07;
  if (animPhase_ > 2.0 * M_PI) animPhase_ -= 2.0 * M_PI;
  dirty_ = true;
  update();
}

void View3D::setPreset(Camera::Preset p) {
  cam_.setPreset(p);
  proj_.valid = false;
  dirty_ = true;
  update();
  emit viewChanged();
}

void View3D::setOrtho(bool on) {
  cam_.setOrtho(on);
  proj_.valid = false;
  dirty_ = true;
  update();
  emit viewChanged();
}

bool View3D::ortho() const { return cam_.ortho(); }

void View3D::frameAll() {
  rebuild();
  // 只调视距，不动朝向 —— 否则"正立面 + Home"会跳回轴测
  if (scene_.boundsValid) cam_.fitBounds(scene_.lo, scene_.hi);
  proj_.valid = false;
  update();
  emit viewChanged();
}

void View3D::resetView() {
  rebuild();
  if (scene_.boundsValid) cam_.frameBounds(scene_.lo, scene_.hi);
  proj_.valid = false;
  update();
  emit viewChanged();
}

QRect View3D::viewport() const {
  // 右侧 100 px 给色标，底部 34 px 给坐标轴三角。
  // 叠加信息不能压到构件上 —— 压住的那个位置往往正好是要看的节点。
  return QRect(2, 2, std::max(1, width() - 104), std::max(1, height() - 36));
}

// =============================================================================
//  场景重建
// =============================================================================
void View3D::rebuild() {
  if (!dirty_) return;
  dirty_ = false;

  if (!model_ || model_->nodeCount() == 0) {
    scene_ = Scene();
    return;
  }

  const StaticResult& r = result_ ? *result_ : empty_;
  SceneBuilder probe(*model_, r, post_);

  RenderOptions o = opt_;

  // 变形比例：先把"自动系数"算出来变成显式系数，再乘动画相位。
  // 顺序不能反 —— 若先乘相位再自动定标，每一帧的自动系数都不同，
  // 变形会变成抖动的噪声而不是平滑摆动。
  if (o.deformed) {
    double s = o.autoScale ? probe.autoDeformScale() : o.deformScale;
    if (!(s > 0.0)) s = 1.0;
    o.autoScale = false;
    o.deformScale = s * animFactor();
  }
  appliedScale_ = o.deformed ? o.deformScale : 0.0;

  scene_ = probe.build(o);
}

double View3D::animFactor() const {
  if (!animating_) return 1.0;
  // 余弦把相位 [0,2π] 映到 [0,1] 且两端导数为零 —— 极值处不会"急停急起"
  return 0.06 + 0.94 * 0.5 * (1.0 - std::cos(animPhase_));
}

void View3D::ensureProj() {
  const QRect vp = viewport();
  if (proj_.valid && proj_.rect == vp) return;
  proj_.rect = vp;
  proj_.vp = cam_.viewProjection(vp);
  proj_.valid = true;
}

// =============================================================================
//  绘制
// =============================================================================
void View3D::paintEvent(QPaintEvent*) {
  rebuild();
  ensureProj();

  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing, true);

  const QRect vp = viewport();
  drawBackground(p, vp);

  // 【不要用 QPainter::setWindow/setViewport 做投影】
  // 投影已经在 drawScene 里手工算到屏幕像素了。再叠一层窗口变换
  // 等于把变换做了两遍 —— 表现为"图形跑到画布外面"或"缩得只剩一点"。
  p.save();
  p.setClipRect(vp);
  drawScene(p, vp);
  p.restore();

  drawAxisTriad(p, vp);
  drawColorBar(p, vp);
  drawHud(p, vp);
}

void View3D::drawBackground(QPainter& p, const QRect& vp) {
  QLinearGradient g(0, 0, 0, height());
  g.setColorAt(0.0, Viewport::bgTop());
  g.setColorAt(1.0, Viewport::bgBottom());
  p.fillRect(rect(), g);
  (void)vp;
}

void View3D::drawScene(QPainter& p, const QRect& vp) {
  if (scene_.faces.empty() && scene_.lines.empty() && scene_.points.empty()) {
    p.setPen(Viewport::textDim());
    p.setFont(uiFont(11));
    p.drawText(vp, Qt::AlignCenter,
               model_ ? QStringLiteral("模型为空 —— 打开一个 .yjk 建模脚本")
                      : QStringLiteral("未加载模型\n\n文件 → 打开建模脚本(*.yjk)"));
    return;
  }

  const QMatrix4x4& vpMat = proj_.vp;
  const QMatrix4x4 view = cam_.view();
  const QVector3D eye = cam_.eye();

  // 手动投影到【屏幕像素】。不用 QPainter 的窗口变换，
  // 因为深度排序、拾取、HUD 都要用同一套屏幕坐标，手工算一次最省心。
  auto toScreen = [&](const QVector4D& c) {
    const double ndcX = c.x() / c.w(), ndcY = c.y() / c.w();
    return QPointF(vp.x() + (ndcX * 0.5 + 0.5) * vp.width(),
                   vp.y() + (1.0 - (ndcY * 0.5 + 0.5)) * vp.height());
  };

  std::vector<DrawItem> items;
  items.reserve(scene_.faces.size() + scene_.lines.size() + scene_.points.size());

  // ---- 面 ----
  for (size_t i = 0; i < scene_.faces.size(); ++i) {
    const Face& f = scene_.faces[i];
    if (f.n < 3) continue;

    QPolygonF poly;
    bool behind = false;
    double depth = 0.0;
    for (int k = 0; k < f.n; ++k) {
      const QVector4D c = vpMat * QVector4D(f.p[k], 1.0f);
      if (c.w() <= 1e-5f) { behind = true; break; }
      poly.append(toScreen(c));
      const QVector4D v = view * QVector4D(f.p[k], 1.0f);
      depth += v.z();
    }
    if (behind) continue;
    depth /= f.n;

    // 背面剔除：只对【闭合体】（梁的长方体）做。
    // 壳单元是单面构件，剔除会让从下方看楼板时整块板消失 —— 那是个 bug 不是优化。
    QVector3D nrm(0, 0, 0);
    if (f.elemType == 1) {
      const QVector3D e1 = f.p[1] - f.p[0];
      const QVector3D e2 = f.p[2] - f.p[0];
      nrm = QVector3D::crossProduct(e1, e2);
      if (nrm.lengthSquared() < 1e-20f) continue;
      nrm.normalize();
      const QVector3D ctr = (f.p[0] + f.p[1] + f.p[2] + f.p[3]) * 0.25f;
      if (QVector3D::dotProduct(nrm, eye - ctr) <= 0.0f) continue;   // 背向相机
    }

    DrawItem it;
    it.kind = 0;
    it.index = static_cast<int>(i);
    it.depth = depth;
    it.screen = poly;
    it.color = (f.elemType == 1 && nrm.lengthSquared() > 0.5f) ? shadeColor(f.fill, nrm) : f.fill;
    items.push_back(std::move(it));
  }

  // ---- 线 ----
  for (size_t i = 0; i < scene_.lines.size(); ++i) {
    const Line& l = scene_.lines[i];
    const QVector4D ca = vpMat * QVector4D(l.a, 1.0f);
    const QVector4D cb = vpMat * QVector4D(l.b, 1.0f);
    if (ca.w() <= 1e-5f || cb.w() <= 1e-5f) continue;
    const QVector4D va = view * QVector4D(l.a, 1.0f);
    const QVector4D vb = view * QVector4D(l.b, 1.0f);
    DrawItem it;
    it.kind = 1;
    it.index = static_cast<int>(i);
    it.depth = 0.5 * (static_cast<double>(va.z()) + static_cast<double>(vb.z()));
    it.a = toScreen(ca);
    it.b = toScreen(cb);
    it.color = l.color;
    it.width = l.width;
    items.push_back(std::move(it));
  }

  // ---- 点 ----
  for (size_t i = 0; i < scene_.points.size(); ++i) {
    const Point& pt = scene_.points[i];
    const QVector4D c = vpMat * QVector4D(pt.p, 1.0f);
    if (c.w() <= 1e-5f) continue;
    const QVector4D v = view * QVector4D(pt.p, 1.0f);
    DrawItem it;
    it.kind = 2;
    it.index = static_cast<int>(i);
    it.depth = v.z();
    it.a = toScreen(c);
    it.color = pt.color;
    it.size = pt.size;
    items.push_back(std::move(it));
  }

  // ---- 深度排序：远（z 更负）先画 ----
  std::stable_sort(items.begin(), items.end(),
                   [](const DrawItem& a, const DrawItem& b) { return a.depth < b.depth; });

  // 构件少时开抗锯齿（边缘干净），多了就关（否则填充率是瓶颈）
  const bool aa = items.size() < 4000;
  p.setRenderHint(QPainter::Antialiasing, aa);

  QPen pen;
  pen.setCosmetic(true);
  for (const DrawItem& it : items) {
    if (it.kind == 0) {
      // ---- 面 ----
      const Face& f = scene_.faces[static_cast<size_t>(it.index)];
      p.setBrush(it.color);
      // 淡描边：把相邻面的接缝勾出来。不描边时同色相邻面会糊成一片，
      // 看不出这是一根梁还是两根拼在一起。
      if (f.edge.alpha() > 0) {
        pen.setColor(f.edge);
        pen.setWidthF(0.7);
        p.setPen(pen);
      } else {
        p.setPen(Qt::NoPen);
      }
      p.drawPolygon(it.screen);
    } else if (it.kind == 1) {
      // ---- 线 ----
      pen.setColor(it.color);
      pen.setWidthF(it.width);
      pen.setCapStyle(Qt::RoundCap);
      p.setPen(pen);
      p.setBrush(Qt::NoBrush);
      p.drawLine(it.a, it.b);
    } else {
      // ---- 点 ----
      p.setPen(Qt::NoPen);
      p.setBrush(it.color);
      const double h = it.size * 0.5;
      p.drawRect(QRectF(it.a.x() - h, it.a.y() - h, it.size, it.size));
    }
  }

  // ---- 叠加高亮层：选中 / 悬停 / 极值构件的轮廓 ----
  //
  // 刻意画在【所有几何之后】，不参与深度排序。
  // 理由是选择高亮的语义是"告诉我选中的东西在哪"，而不是"如实遮挡"：
  // 被前面构件挡住的选中构件如果不可见，用户会以为没选中。
  // 先描一圈深色再描亮色，保证在深色背景和浅色云图上都有对比度。
  if (!scene_.highlight.empty()) {
    for (int pass = 0; pass < 2; ++pass) {
      for (const Line& l : scene_.highlight) {
        const QVector4D ca = vpMat * QVector4D(l.a, 1.0f);
        const QVector4D cb = vpMat * QVector4D(l.b, 1.0f);
        if (ca.w() <= 1e-5f || cb.w() <= 1e-5f) continue;
        if (pass == 0) {
          pen.setColor(QColor(0, 0, 0, 190));
          pen.setWidthF(l.width + 3.5);
        } else {
          pen.setColor(l.color);
          pen.setWidthF(l.width);
        }
        pen.setCapStyle(Qt::RoundCap);
        p.setPen(pen);
        p.setBrush(Qt::NoBrush);
        p.drawLine(toScreen(ca), toScreen(cb));
      }
    }
  }

  p.setRenderHint(QPainter::Antialiasing, true);
}

// -----------------------------------------------------------------------------
//  坐标轴三角（左下角）
//
//  【必须画这个】三维视图没有它，用户无法判断当前朝哪个方向看 ——
//  轴测视角下"左上"和"右上"在视觉上是对称的，方向感全靠这个三角。
// -----------------------------------------------------------------------------
void View3D::drawAxisTriad(QPainter& p, const QRect& vp) {
  const int ox = 42;
  const int oy = vp.height() - 18;
  const int len = 26;

  // 只取相机的旋转部分（不含平移），这样三角只表示朝向不表示位置
  const QVector3D right = cam_.rightDir();
  const QVector3D up = cam_.upDir();
  auto proj = [&](const QVector3D& w) {
    // 世界方向 → 屏幕偏移：世界 +X 在屏幕上的投影 = (dot(X,right), -dot(X,up))
    return QPointF(ox + len * QVector3D::dotProduct(w, right),
                   oy - len * QVector3D::dotProduct(w, up));
  };

  const QPointF o(ox, oy);
  struct Ax { QVector3D d; QColor c; const char* n; };
  const Ax axes[3] = {
      {QVector3D(1, 0, 0), Viewport::axisX(), "X"},
      {QVector3D(0, 1, 0), Viewport::axisY(), "Y"},
      {QVector3D(0, 0, 1), Viewport::axisZ(), "Z"},
  };

  p.setRenderHint(QPainter::Antialiasing, true);
  for (const Ax& a : axes) {
    const QPointF e = proj(a.d);
    p.setPen(QPen(a.c, 1.8));
    p.drawLine(o, e);
    p.setPen(a.c);
    p.setFont(uiFont(8));
    p.drawText(e + QPointF(3, 3), a.n);
  }
  p.setPen(QPen(Viewport::textDim(), 1.0));
  p.setBrush(Viewport::textDim());
  p.drawEllipse(o, 2.0, 2.0);
}

// -----------------------------------------------------------------------------
//  色标
// -----------------------------------------------------------------------------
void View3D::drawColorBar(QPainter& p, const QRect& vp) {
  (void)vp;        // 色标固定在视口最右侧，与视口矩形无关
  if (!scene_.hasContour()) return;

  const int barW = 18;
  const int x = width() - barW - 42;
  const int top = 46;
  const int bot = height() - 78;
  if (bot <= top + 40) return;

  p.setRenderHint(QPainter::Antialiasing, true);

  // 连续渐变色条
  QLinearGradient g(0, top, 0, bot);
  for (int i = 0; i <= 32; ++i) {
    const double t = 1.0 - i / 32.0;      // 顶部 = 最大值
    g.setColorAt(i / 32.0, ColorMap::rainbow(t));
  }
  p.setPen(QPen(QColor(255, 255, 255, 60), 1));
  p.setBrush(g);
  p.drawRect(x, top, barW, bot - top);

  // 刻度：用 1/2/5 的整数步长，工程师能心算的量级
  const double span = scene_.vHi - scene_.vLo;
  const double step = ColorMap::niceStep(span, 6);
  p.setFont(uiFont(8));
  const QFontMetricsF fm(p.font());
  for (double v = std::ceil(scene_.vLo / step) * step; v <= scene_.vHi + step * 1e-6; v += step) {
    const double t = (v - scene_.vLo) / span;
    const int y = static_cast<int>(bot - t * (bot - top));
    p.setPen(QColor(255, 255, 255, 110));
    p.drawLine(x + barW, y, x + barW + 4, y);
    const QString lab = ColorMap::label(v);
    p.setPen(Viewport::text());
    p.drawText(QPointF(x + barW + 7, y + fm.height() * 0.35), lab);
  }

  // 标题
  p.setFont(uiFont(9));
  p.setPen(Viewport::text());
  const QString title = QString::fromStdString(scene_.valueLabel);
  p.drawText(QRectF(x - 46, top - 22, barW + 92, 18), Qt::AlignCenter, title);
  if (!scene_.valueUnit.empty()) {
    p.setFont(uiFont(8));
    p.setPen(Viewport::textDim());
    p.drawText(QRectF(x - 46, top - 10, barW + 92, 14), Qt::AlignCenter,
               QString::fromStdString(scene_.valueUnit));
  }
}

// -----------------------------------------------------------------------------
//  信息层
// -----------------------------------------------------------------------------
void View3D::drawHud(QPainter& p, const QRect& vp) {
  p.setFont(uiFont(9));

  // 左上：模型统计
  QStringList lines;
  if (model_) {
    lines << QStringLiteral("节点 %1 · 单元 %2").arg(model_->nodeCount()).arg(model_->elementCount());
  }
  if (scene_.maxDisp > 0.0) {
    lines << QStringLiteral("最大位移 %1 m").arg(fmtNum(scene_.maxDisp, 3));
  }
  if (opt_.deformed && appliedScale_ > 0.0) {
    lines << QStringLiteral("变形放大 ×%1").arg(fmtNum(appliedScale_, 3));
  }

  if (!lines.isEmpty()) {
    const QFontMetricsF fm(p.font());
    double w = 0;
    for (const QString& s : lines) w = std::max(w, fm.horizontalAdvance(s));
    const QRectF box(10, 10, w + 20, fm.height() * lines.size() + 12);
    p.setPen(Qt::NoPen);
    p.setBrush(QColor(0, 0, 0, 90));
    p.drawRoundedRect(box, 4, 4);
    p.setPen(Viewport::text());
    double y = box.top() + 6 + fm.ascent();
    for (const QString& s : lines) {
      p.drawText(QPointF(box.left() + 10, y), s);
      y += fm.height();
    }
  }

  // 底部：鼠标提示
  if (!hint_.isEmpty()) {
    const QFontMetricsF fm(p.font());
    const QRectF box(vp.left() + 90, vp.bottom() - 22,
                     fm.horizontalAdvance(hint_) + 18, fm.height() + 8);
    p.setPen(Qt::NoPen);
    p.setBrush(QColor(0, 0, 0, 110));
    p.drawRoundedRect(box, 4, 4);
    p.setPen(Viewport::text());
    p.drawText(box.adjusted(9, 2, -9, -2), Qt::AlignVCenter | Qt::AlignLeft, hint_);
  }

  // 空模型的引导文字
  if ((!model_ || model_->elementCount() == 0) && scene_.faces.empty()) {
    p.setPen(Viewport::textDim());
    p.setFont(uiFont(11));
    p.drawText(vp, Qt::AlignCenter,
               QStringLiteral("未加载模型\n\n文件 → 打开建模脚本（*.yjk），或直接点「开始计算」用示例模型试跑"));
  }
}

// =============================================================================
//  拾取
//
//  【为什么换成射线拾取】旧实现把每条梁投影回屏幕，再用"点到线段距离"判定，
//  对梁/壳尚可，但节点（点）在屏幕上投影后依然是点，永远点不中；
//  而且投影判定无法正确区分前后遮挡，放大后密集区域的命中很随意。
//  现在用世界空间射线 + yjk::interact::pickNearest 纯函数（在测试里独立验证
//  过）：View3D 只负责"像素 → 射线"的换算，命中逻辑与脚本层共用一套。
//
//  【容差必须随缩放变化】固定 px 容差在缩得很小的时候会把整片结构
//  都算作"命中"，放得很大时又点不中；按视距 + fov + 视口高换算成
//  世界单位容差，效果等价于"屏幕上固定 6 px 的容差带"。
// =============================================================================
int View3D::pickAt(const QPoint& pos, int* type) const {
  *type = 0;
  if (!model_) return -1;

  const QRect vp = viewport();
  if (!vp.contains(pos)) return -1;

  // 像素 → 世界空间射线：相机眼睛是原点，指向光标反投影点。
  const QVector3D eye = cam_.eye();
  const QVector3D dir = cam_.unproject(pos.x(), pos.y(), vp) - eye;
  const double tMax = cam_.radius() * 12.0 + 1.0;  // 与投影 far 平面一致

  // 屏幕容差(px) → 世界容差(m)：由 unproject 差分直接得到，
  // 对透视/正交投影都自动正确，不需要知道 fov 的具体值。
  const double worldPerPix = worldPerPixAt(cam_, vp);
  const double nodeR = 6.0 * worldPerPix;
  const double beamR = 6.0 * worldPerPix;

  // 组装可拾取图元（节点/梁中心线/壳四边形），都带变形。
  std::vector<yjk::interact::PickNode> nodes;
  nodes.reserve(scene_.pickNodes.size());
  for (const Scene::PickNode& pn : scene_.pickNodes) {
    nodes.push_back({static_cast<yjk::Id>(pn.id),
                     {pn.p.x(), pn.p.y(), pn.p.z()}});
  }
  std::vector<yjk::interact::PickBeam> beams;
  beams.reserve(scene_.pickBeams.size());
  for (const Scene::PickBeam& pb : scene_.pickBeams) {
    beams.push_back({static_cast<yjk::Id>(pb.elem),
                     {pb.a.x(), pb.a.y(), pb.a.z()},
                     {pb.b.x(), pb.b.y(), pb.b.z()}});
  }
  std::vector<yjk::interact::PickShell> shells;
  for (const Face& f : scene_.faces) {
    if (f.elemType != 2) continue;
    yjk::interact::PickShell sh;
    sh.id = static_cast<yjk::Id>(f.elem);
    for (int i = 0; i < 4; ++i) sh.p[i] = {f.p[i].x(), f.p[i].y(), f.p[i].z()};
    shells.push_back(sh);
  }

  const yjk::interact::PickHit hit = yjk::interact::pickNearest(
      yjk::interact::Ray3({eye.x(), eye.y(), eye.z()},
                          {dir.x(), dir.y(), dir.z()}),
      nodes, beams, shells, nodeR, beamR, tMax);
  if (hit.kind == yjk::interact::PickKind::None) return -1;
  *type = (hit.kind == yjk::interact::PickKind::Node) ? 3
        : (hit.kind == yjk::interact::PickKind::Beam) ? 1 : 2;
  return static_cast<int>(hit.id);
}

// =============================================================================
//  交互
// =============================================================================
void View3D::mousePressEvent(QMouseEvent* e) {
  lastPos_ = e->pos();
  setFocus();

  // 中键 = 平移，这是 CAD/CAE 的通用约定。
  // 也支持 Shift + 左键，因为很多笔记本没有中键。
  if (e->button() == Qt::MiddleButton ||
      (e->button() == Qt::LeftButton && (e->modifiers() & Qt::ShiftModifier))) {
    panning_ = true;
    setCursor(Qt::ClosedHandCursor);
  } else if (e->button() == Qt::LeftButton) {
    rotating_ = true;
    setCursor(Qt::SizeAllCursor);
  }
}

void View3D::mouseMoveEvent(QMouseEvent* e) {
  const QPoint d = e->pos() - lastPos_;

  if (rotating_) {
    cam_.orbit(d.x(), d.y());
    proj_.valid = false;
    lastPos_ = e->pos();
    update();
    emit viewChanged();
    return;
  }
  if (panning_) {
    cam_.pan(d.x(), d.y(), viewport());
    proj_.valid = false;
    lastPos_ = e->pos();
    update();
    emit viewChanged();
    return;
  }

  // 悬停高亮
  int t = 0;
  const int el = pickAt(e->pos(), &t);
  if (el != opt_.hoverElem || t != opt_.hoverType) {
    opt_.hoverElem = el;
    opt_.hoverType = t;
    dirty_ = true;
    update();
    emit hoverChanged(el, t);
  }

  // 提示文字
  QString h;
if (el >= 0 && t == 1) {
    h = QStringLiteral("梁 #%1").arg(el);
    if (post_ && el < static_cast<int>(post_->beamResults().size())) {
      const auto& br = post_->beamResults()[static_cast<size_t>(el)];
      h += QStringLiteral("  L=%1m  M=%2kN·m  V=%3kN  N=%4kN")
               .arg(fmtNum(br.L, 3), fmtNum(br.Mmax, 3), fmtNum(br.Vmax, 3),
                    fmtNum(0.5 * (br.Nmax + br.Nmin), 3));
    }
  } else if (el >= 0 && t == 2) {
    h = QStringLiteral("壳 #%1").arg(el);
    if (post_ && el < static_cast<int>(post_->shellResults().size())) {
      const auto& sr = post_->shellResults()[static_cast<size_t>(el)];
      h += QStringLiteral("  面积=%1m²").arg(fmtNum(sr.area, 3));
    }
  } else if (el >= 0 && t == 3) {
    // 节点：从拾取图元里找坐标（含变形），显示 id 与三向坐标。
    for (const Scene::PickNode& pn : scene_.pickNodes) {
      if (pn.id == el) {
        h = QStringLiteral("节点 #%1  (%2, %3, %4) m")
                .arg(el)
                .arg(fmtNum(pn.p.x(), 3), fmtNum(pn.p.y(), 3), fmtNum(pn.p.z(), 3));
        break;
      }
    }
    if (h.isEmpty())
      h = QStringLiteral("节点 #%1").arg(el);
  }
  if (h != hint_) { hint_ = h; update(); }
}

void View3D::mouseReleaseEvent(QMouseEvent* e) {
  const QPoint d = e->pos() - lastPos_;
  const bool wasRotating = rotating_;
  const bool wasPanning = panning_;
  rotating_ = false;
  panning_ = false;
  unsetCursor();

  // 点击（几乎没拖动）才算拾取。
  // 【必须判断这个】否则每次旋转结束都会顺手选中一个构件，
  // 用户会觉得"选中的东西总是自己变"。
  if (e->button() == Qt::LeftButton && !wasPanning && wasRotating &&
      d.manhattanLength() < 4) {
    int t = 0;
    const int el = pickAt(e->pos(), &t);
if (el >= 0) {
      if (t == 3) {
        // 节点拾取：与单元选择互斥 —— 选中节点时清空单元选择。
        sel_.clear();
        selType_ = 0;
        opt_.selected.clear();
        selNode_ = el;
        opt_.selectedNodes.assign({el});
        dirty_ = true;
        update();
        emit nodePicked(el);
        emit selectionChanged();
        return;
      }
      sel_.clear();
      opt_.selected.clear();
      selType_ = t;
      sel_.push_back(el);
      opt_.selected.push_back(el);
      selNode_ = -1;
      opt_.selectedNodes.clear();
      dirty_ = true;
      update();
      emit elementPicked(el, t);
      emit selectionChanged();
      return;
    }
    // 点空白处 = 取消选择
    if (!sel_.empty() || selNode_ >= 0) {
      sel_.clear();
      selType_ = 0;
      opt_.selected.clear();
      selNode_ = -1;
      opt_.selectedNodes.clear();
      dirty_ = true;
      update();
      emit selectionChanged();
    }
  }
}

void View3D::wheelEvent(QWheelEvent* e) {
  const double steps = e->angleDelta().y() / 120.0;
  if (steps == 0.0) return;

  // 以光标为锚点缩放：光标指的位置在屏幕上保持不动。
  // 【这个细节决定了"能不能放大看一个节点"】不锚定的话，
  // 每次滚轮都要重新把想看的位置拖回中间，密集框架下极其烦人。
  const QPointF pos = e->position();
  if (viewport().contains(pos.toPoint())) {
    const QVector3D anchor = cam_.unproject(pos.x(), pos.y(), viewport());
    cam_.zoom(steps, &anchor);
  } else {
    cam_.zoom(steps, nullptr);
  }
  proj_.valid = false;
  update();
  emit viewChanged();
}

void View3D::keyPressEvent(QKeyEvent* e) {
  switch (e->key()) {
    case Qt::Key_Home:  frameAll(); return;
    case Qt::Key_F:     frameAll(); return;
    case Qt::Key_1:     setPreset(Camera::Preset::Iso);    return;
    case Qt::Key_2:     setPreset(Camera::Preset::Front);  return;
    case Qt::Key_3:     setPreset(Camera::Preset::Left);   return;
    case Qt::Key_4:     setPreset(Camera::Preset::Top);    return;
    case Qt::Key_O:     setOrtho(!ortho());                return;
    case Qt::Key_Space: setAnimating(!animating());        return;
    case Qt::Key_Escape: clearSelection();                 return;
    default: break;
  }
  QWidget::keyPressEvent(e);
}

void View3D::leaveEvent(QEvent*) {
  if (opt_.hoverElem >= 0) {
    opt_.hoverElem = -1;
    opt_.hoverType = 0;
    hint_.clear();
    dirty_ = true;
    update();
    emit hoverChanged(-1, 0);
  }
}

void View3D::resizeEvent(QResizeEvent* e) {
  proj_.valid = false;
  QWidget::resizeEvent(e);
}

QImage View3D::grabImage() {
  rebuild();
  ensureProj();
  QImage img(size() * devicePixelRatioF(), QImage::Format_ARGB32_Premultiplied);
  img.setDevicePixelRatio(devicePixelRatioF());
  img.fill(Qt::transparent);
  QPainter p(&img);
  p.setRenderHint(QPainter::Antialiasing, true);
  const QRect vp = viewport();
  drawBackground(p, vp);
  p.save();
  p.setClipRect(vp);
  drawScene(p, vp);
  p.restore();
  drawAxisTriad(p, vp);
  drawColorBar(p, vp);
  drawHud(p, vp);
  return img;
}

}  // namespace ui
