// =============================================================================
//  apps/ui/PropertyPanel.h  ——  属性面板
//
//  回答一个问题："我刚点中的这个构件到底是什么、受力多大。"
//
//  【为什么属性要用"分组表"而不是一行行 label】
//  一根梁的属性有 20 多项（几何、截面、材料、内力、应力）。
//  平铺成一列时，找"这根的混凝土等级是多少"要扫过十几行数字；
//  分组之后（几何 / 截面 / 内力 / 验算）每一组只有 4~6 行，眼睛能跳着看。
// =============================================================================
#pragma once

#include <QWidget>

#include "SceneData.h"

class QTreeWidget;
class QTreeWidgetItem;
class QLabel;

namespace ui {

class PropertyPanel : public QWidget {
  Q_OBJECT
 public:
  explicit PropertyPanel(QWidget* parent = nullptr);

  void setModel(const Model* m);
  void setResult(const StaticResult* r);
  void setPost(const post::PostProcessor* pp);

  void showElement(int elem, int type);
  void showNode(int node);
  void clear();

 private:
  void rebuild();

  QLabel* title_{nullptr};
  QTreeWidget* tree_{nullptr};

  const Model* model_{nullptr};
  const StaticResult* result_{nullptr};
  const post::PostProcessor* post_{nullptr};

  int curElem_{-1};
  int curType_{0};
  int curNode_{-1};
};

}  // namespace ui
