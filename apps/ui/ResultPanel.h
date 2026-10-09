// =============================================================================
//  apps/ui/ResultPanel.h  ——  结果面板
//
//  内容按"工程师实际翻结果的顺序"组织：
//    ① 楼层指标 —— 抗震审图第一眼看的东西（层间位移角、剪重比）
//    ② 整体指标 —— 总重、基底剪力、最大位移
//    ③ 构件控制 —— 哪几根构件最危险
//    ④ 文本报告 —— 可以整段复制进计算书
//  每个表都【可排序 + 可导出为 CSV】，因为真正的使用场景是
//  "按层间位移角排序看薄弱层"、"把构件内力表发给同事核对"。
// =============================================================================
#pragma once

#include <QWidget>

#include <vector>

#include "SceneData.h"

class QTabWidget;
class QTableWidget;
class QPlainTextEdit;
class QLabel;

namespace ui {

class ResultPanel : public QWidget {
  Q_OBJECT
 public:
  explicit ResultPanel(QWidget* parent = nullptr);

  void setModel(const Model* m);
  void setResult(const StaticResult* r);
  void setPost(const post::PostProcessor* pp);

  void refresh();

 signals:
  // 在表里双击某一行 → 视口高亮该构件
  void elementActivated(int elem, int type);
  void storyActivated(int story);

 private:
  void buildStoryTable();
  void buildEnvelope();
  void buildMemberTable();
  void buildReport();

  QTabWidget* tabs_{nullptr};
  QTableWidget* story_{nullptr};
  QTableWidget* envelope_{nullptr};
  QTableWidget* members_{nullptr};
  QPlainTextEdit* report_{nullptr};
  QLabel* summary_{nullptr};

  const Model* model_{nullptr};
  const StaticResult* result_{nullptr};
  const post::PostProcessor* post_{nullptr};
};

}  // namespace ui
