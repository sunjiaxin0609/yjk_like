// =============================================================================
//  apps/ui/ModelTreePanel.h  ——  左侧模型树
//
//  【为什么模型树是结构软件里用得最多的面板】
//  三维视口里"点中一根梁"在密集框架里很难做准（前后遮挡、像素太细），
//  而模型树是【确定性】的：要看第 3 层第 5 根柱，就从树上点它。
//  所以树不是装饰，它是三维视口的互补入口，两者必须双向联动：
//  树上点 → 视口高亮；视口点 → 树上展开并选中。
// =============================================================================
#pragma once

#include <QTreeWidget>
#include <QWidget>

#include <vector>

#include "SceneData.h"

namespace ui {

class ModelTreePanel : public QWidget {
  Q_OBJECT
 public:
  explicit ModelTreePanel(QWidget* parent = nullptr);

  void setModel(const Model* m);
  void setPost(const post::PostProcessor* pp);

  // 从视口反选：展开到该构件并在树上高亮。
  // blockSignal 防止"树上选中 → 通知视口 → 视口再通知树"的死循环。
  void selectElement(int elem, int type);

 signals:
  void elementActivated(int elem, int type);
  void nodeActivated(int node);

 private:
  void rebuild();

  QTreeWidget* tree_{nullptr};
  const Model* model_{nullptr};
  const post::PostProcessor* post_{nullptr};
  bool block_{false};
};

}  // namespace ui
