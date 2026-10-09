// =============================================================================
//  apps/ui/View3D.h  ——  三维视口
//
//  【为什么用 QPainter 软件渲染，而不是 QOpenGLWidget】
//  这是本项目唯一一个"选了看起来更笨的方案"的地方，理由要说清楚：
//
//  ① 可用性优先。OpenGL 在远程桌面、虚拟机、老旧集显上会退到 1.1 软实现，
//     而现代 Qt 的 QOpenGLWidget 要求 3.3 core profile —— 直接白屏。
//     结构工程师恰恰经常在远程桌面上跑设计软件。白屏是致命的，
//     性能差一点只是体验问题。
//  ② 本模型的规模用不着 GPU。一个 3×3 柱网 4 层的框架，实体显示也就
//     一千多个四边形；到一万个构件时 QPainter 仍在 60fps 量级。
//     真要上十万构件，正确的做法是 LOD + 分块裁剪，不是简单换个后端。
//  ③ 软件渲染让深度排序、拾取、叠加文字共用同一个坐标系，
//     不会出现"GPU 画的几何和 CPU 画的标注对不上"这类经典麻烦。
//
//  实现的是标准的小型软件三维管线：
//      世界 → 视图 → 投影 → 视口 → 画家算法排序 → 双光源着色 → QPainter
// =============================================================================
#pragma once

#include <QColor>
#include <QImage>
#include <QMatrix4x4>
#include <QPoint>
#include <QRect>
#include <QString>
#include <QTimer>
#include <QWidget>

#include <memory>
#include <vector>

#include "Camera.h"
#include "SceneData.h"
#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/model/Model.h"
#include "yjk/post/PostProcessor.h"

namespace ui {

class View3D : public QWidget {
  Q_OBJECT
 public:
  explicit View3D(QWidget* parent = nullptr);

  // 数据源。三者必须同时有效才算"有结果"；
  // 只有 model 时显示线框 —— 用户确认"模型建对了没"靠的就是这一步。
  void setModel(const Model* m);
  void setResult(const StaticResult* r);
  void setPost(const post::PostProcessor* pp);
  void clearPost();

  const RenderOptions& options() const { return opt_; }
  void setOptions(const RenderOptions& o);

  void refresh();
  void frameAll();                 // 缩放到全览（保持当前视角）
  void resetView();                // 全览 + 回到默认轴测朝向
  void setPreset(Camera::Preset p);
  void setOrtho(bool on);
  bool ortho() const;

  // 变形动画：在 0 与放大系数之间往复。
  // 【为什么要动起来】静力变形图看不出"哪里在动"。动起来之后，
  // 薄弱层、扭转效应一眼就能看出来 —— 这是抗震设计最常用的判断手段。
  void setAnimating(bool on);
  bool animating() const { return animating_; }

  void setHighlightExtreme(int elem, int type);
  void clearSelection();
  void selectElement(int elem, int type);
  const std::vector<int>& selection() const { return sel_; }
  int selectionType() const { return selType_; }

  QImage grabImage();      // 导出截图

 signals:
  void elementPicked(int elem, int type);   // type: 1 = 梁，2 = 壳
  void nodePicked(int node);
  void selectionChanged();
  void viewChanged();
  void hoverChanged(int elem, int type);

 protected:
  void paintEvent(QPaintEvent*) override;
  void resizeEvent(QResizeEvent*) override;
  void mousePressEvent(QMouseEvent*) override;
  void mouseMoveEvent(QMouseEvent*) override;
  void mouseReleaseEvent(QMouseEvent*) override;
  void wheelEvent(QWheelEvent*) override;
  void keyPressEvent(QKeyEvent*) override;
  void leaveEvent(QEvent*) override;

 private slots:
  void onAnimTick();

 private:
  QRect viewport() const;
  void rebuild();
  void ensureProj();
  double animFactor() const;

  void drawBackground(QPainter& p, const QRect& vp);
  void drawScene(QPainter& p, const QRect& vp);
  void drawAxisTriad(QPainter& p, const QRect& vp);
  void drawColorBar(QPainter& p, const QRect& vp);
  void drawHud(QPainter& p, const QRect& vp);

  int pickAt(const QPoint& pos, int* type) const;

  const Model* model_{nullptr};
  const StaticResult* result_{nullptr};
  const post::PostProcessor* post_{nullptr};

  StaticResult empty_;         // result_ 为空时的替身，避免到处判空

  Camera cam_;
  RenderOptions opt_;
  Scene scene_;
  bool dirty_{true};
  double appliedScale_{0.0};

  struct ProjCache {
    QMatrix4x4 vp;
    QRect rect;
    bool valid{false};
  };
  ProjCache proj_;

  QPoint lastPos_;
  bool rotating_{false};
  bool panning_{false};

std::vector<int> sel_;
  int selType_{0};
  int selNode_{-1};   // 选中的节点 id（-1 = 未选中节点）

  bool animating_{false};
  double animPhase_{M_PI};
  QTimer animTimer_;

  QString hint_;
};

}  // namespace ui
