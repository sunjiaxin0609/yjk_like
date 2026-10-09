// =============================================================================
//  apps/ui/Camera.h  ——  三维视图相机（轨道式）
//
//  【为什么用"轨道相机"而不是"自由飞行相机"】
//  结构软件的视角操作有明确的工程语义：
//    · 旋转 = 绕着结构转，不是绕着相机转 —— 结构始终在视野中心
//    · 平移 = 平移视点，但保持视角不变
//    · 缩放 = 改变视距，透视感随之变化
//  自由飞行相机的缺点是"转着转着就找不到结构了"，工程师必须常驻一个
//  "全览"操作把它拉回来。轨道相机天然不会丢目标。
//
//  【视图状态可保存/恢复】
//  结构设计是一个反复看同一个立面、同一个角度的过程。相机参数必须能
//  一键回到预设视角（前/后/左/右/俯/轴测），否则每次都要重新找角度。
// =============================================================================
#pragma once

#include <QMatrix4x4>
#include <QRect>
#include <QVector3D>

#include <algorithm>
#include <cmath>

namespace ui {

class Camera {
 public:
  // ---- 预设视角 ----
  enum class Preset { Iso, Front, Back, Left, Right, Top, Bottom };

  Camera() { reset(); }

  void reset() {
    target_ = QVector3D(0, 0, 0);
    radius_ = 10.0;
    yaw_ = -55.0;       // 绕 Z 轴（工程惯例：Z 向上，所以方位角绕 Z）
    pitch_ = 24.0;      // 与水平面夹角
    ortho_ = false;
    fovY_ = 32.0;
  }

  // 结构跨度（用于自动定视距）——模型换一个大小完全不同的工程时，
  // 相机距离不跟着变就会"要么一个点，要么跑到模型内部"。
  //
  // 【视距必须由视场角反算，不能直接取对角线的倍数】
  // 在距离 r 处，纵向可见高度 = 2·r·tan(fovY/2)。
  // 要让对角线 diag 装得下，需 r = diag / (2·tan(fovY/2))。
  // 早期版本直接取 r = 0.67·diag（当 fovY=32° 时可见高度只有 0.39·diag），
  // 结果模型涨出屏幕 2.5 倍，看上去像"相机钻进了楼板里"。
  // 只调视距与视点中心，【不动朝向】。
  //
  // 【这一点非常关键】"缩放到全览"（Home / Zoom Extents）在 CAD 里
  // 是【保持当前视角】的。如果它顺手把朝向也归零，用户刚切到正立面、
  // 一按 Home 就跳回轴测，预设视角功能等于废掉。
  void fitBounds(const QVector3D& lo, const QVector3D& hi, double fill = 1.55) {
    target_ = (lo + hi) * 0.5f;
    const QVector3D d = hi - lo;
    const double diag = std::sqrt(static_cast<double>(d.x() * d.x() + d.y() * d.y() + d.z() * d.z()));
    const double halfFov = fovY_ * M_PI / 360.0;
    const double t = std::max(std::tan(halfFov), 1e-3);
    radius_ = std::max(diag * fill / (2.0 * t), 1e-3);
  }

  // 定视距 + 归位到默认轴测朝向。用于【首次加载模型】和【重置视角】。
  void frameBounds(const QVector3D& lo, const QVector3D& hi, double fill = 1.55) {
    fitBounds(lo, hi, fill);
    // 默认稍微俯视，这样能同时看到平面布置与竖向高度
    yaw_ = -55.0;
    pitch_ = 24.0;
  }

  void setPreset(Preset p) {
    switch (p) {
      // 轴测：同时能看到三个方向，是"总览"视角。
      // 方位角取 -55°、仰角 24° 是工程图里最常用的轴测方向（类似西南等轴测）。
      case Preset::Iso:    yaw_ = -55.0; pitch_ = 24.0;  break;
      // 正立面：从 -Y 方向看（工程立面图的习惯视向）
      case Preset::Front:  yaw_ = -90.0; pitch_ = 0.0;   break;
      case Preset::Back:   yaw_ =  90.0; pitch_ = 0.0;   break;
      case Preset::Left:   yaw_ = 180.0; pitch_ = 0.0;   break;
      case Preset::Right:  yaw_ =   0.0; pitch_ = 0.0;   break;
      // 俯视（平面图）：pitch 限到 89.5° 而不是 90°，
      // 因为正 90° 时"上方向"退化，视图会突然翻面。
      case Preset::Top:    yaw_ = -90.0; pitch_ = 89.5;  break;
      case Preset::Bottom: yaw_ = -90.0; pitch_ = -89.5; break;
    }
  }

  // ---- 交互 ----
  // 轨道旋转。限位在 ±89.5°，理由同上：越过极点会让视图翻转。
  void orbit(double dxPix, double dyPix) {
    yaw_ -= dxPix * 0.35;
    pitch_ = std::clamp(pitch_ + dyPix * 0.35, -89.5, 89.5);
    while (yaw_ > 180.0) yaw_ -= 360.0;
    while (yaw_ < -180.0) yaw_ += 360.0;
  }

  // 平移：把屏幕像素位移换算成世界位移。
  // 关键是【按视距缩放】—— 拉得远时拖动同样像素要移动更多距离，
  // 否则远处的小位移会让结构飞出屏幕，近处又完全拖不动。
  void pan(double dxPix, double dyPix, const QRect& vp) {
    if (vp.height() <= 0) return;
    const double worldPerPix = 2.0 * radius_ * std::tan(fovY_ * M_PI / 360.0) /
                               static_cast<double>(vp.height());
    const QVector3D right = rightDir();
    const QVector3D up = upDir();
    target_ += (-dxPix * worldPerPix) * right + (dyPix * worldPerPix) * up;
  }

  // 缩放：乘法不是加法。加法缩放在接近结构时会一步穿过，
  // 乘法缩放（每格 ×1.1）在整个尺度范围上手感一致。
  void zoom(double steps, const QVector3D* anchor = nullptr) {
    const double f = std::pow(1.12, -steps);
    if (!anchor) { radius_ = std::max(radius_ * f, 1e-4); return; }
    // 以光标处为锚点缩放：光标所指的点在屏幕上保持不动。
    // 这是三维软件的标准手感，能"盯着一个节点放大看"。
    const QVector3D a = *anchor;
    target_ = a + (target_ - a) * static_cast<float>(f);
    radius_ = std::max(radius_ * f, 1e-4);
  }

  void setOrtho(bool on) { ortho_ = on; }
  bool ortho() const { return ortho_; }

  // ---- 派生量 ----
  QVector3D eye() const {
    const double ry = yaw_ * M_PI / 180.0, rp = pitch_ * M_PI / 180.0;
    const double cp = std::cos(rp);
    return target_ + static_cast<float>(radius_) *
                         QVector3D(static_cast<float>(cp * std::cos(ry)),
                                   static_cast<float>(cp * std::sin(ry)),
                                   static_cast<float>(std::sin(rp)));
  }
  QVector3D forward() const { return (target_ - eye()).normalized(); }
  QVector3D rightDir() const {
    return QVector3D::crossProduct(forward(), QVector3D(0, 0, 1)).normalized();
  }
  // 工程惯例：Z 轴向上。若视线正好与 Z 平行（俯视图），
  // cross 会退化为零向量，此时改用 Y 作为参考轴。
  QVector3D upDir() const {
    QVector3D r = QVector3D::crossProduct(forward(), QVector3D(0, 0, 1));
    if (r.lengthSquared() < 1e-12f) r = QVector3D::crossProduct(forward(), QVector3D(0, 1, 0));
    r.normalize();
    return QVector3D::crossProduct(r, forward()).normalized();
  }

  double radius() const { return radius_; }
  QVector3D target() const { return target_; }

  QMatrix4x4 view() const {
    QMatrix4x4 m;
    m.lookAt(eye(), target_, upDir());
    return m;
  }

  // 投影矩阵。近远平面按视距动态取 —— 固定 0.1/1000 在
  // 小模型（米级）上会把深度精度全浪费在远处，导致面片闪烁(z-fighting)。
  QMatrix4x4 projection(const QRect& vp) const {
    QMatrix4x4 m;
    const double aspect = vp.height() > 0
                              ? static_cast<double>(vp.width()) / vp.height()
                              : 1.0;
    const double near_ = std::max(radius_ * 0.005, 1e-3);
    const double far_ = radius_ * 12.0 + 1.0;
    if (ortho_) {
      const double h = radius_ * std::tan(fovY_ * M_PI / 360.0);
      m.ortho(-h * aspect, h * aspect, -h, h, -far_, far_);
    } else {
      m.perspective(fovY_, aspect, near_, far_);
    }
    return m;
  }

  QMatrix4x4 viewProjection(const QRect& vp) const { return projection(vp) * view(); }

  // 世界 → 屏幕（含视口映射）。返回 false 表示点在相机后方。
  bool project(const QVector3D& w, const QRect& vp, QMatrix4x4& vpMat,
               QPointF* out, double* depthOut = nullptr) const {
    const QVector4D c = vpMat * QVector4D(w, 1.0f);
    if (std::abs(static_cast<double>(c.w())) < 1e-12) return false;
    const double ndcX = c.x() / c.w();
    const double ndcY = c.y() / c.w();
    const double ndcZ = c.z() / c.w();
    out->setX(vp.x() + (ndcX * 0.5 + 0.5) * vp.width());
    out->setY(vp.y() + (1.0 - (ndcY * 0.5 + 0.5)) * vp.height());
    if (depthOut) *depthOut = ndcZ;
    return true;
  }

  // 把屏幕坐标反投影成一个世界点，用于"以光标为锚点缩放"。
  // 取近平面上的点即可 —— 只用来定方向，不需要精确的交点。
  QVector3D unproject(double sx, double sy, const QRect& vp) const {
    const double ndcX = vp.width() > 0 ? 2.0 * (sx - vp.x()) / vp.width() - 1.0 : 0.0;
    const double ndcY = vp.height() > 0 ? 1.0 - 2.0 * (sy - vp.y()) / vp.height() : 0.0;
    const double aspect = vp.height() > 0
                              ? static_cast<double>(vp.width()) / vp.height()
                              : 1.0;
    const double h = std::tan(fovY_ * M_PI / 360.0);
    return eye() + static_cast<float>(radius_) * forward() +
           static_cast<float>(ndcX * h * aspect * radius_) * rightDir() +
           static_cast<float>(ndcY * h * radius_) * upDir();
  }

 private:
  QVector3D target_{0, 0, 0};
  double radius_{10.0};
  double yaw_{-55.0};
  double pitch_{24.0};
  bool ortho_{false};
  double fovY_{32.0};
};

}  // namespace ui
