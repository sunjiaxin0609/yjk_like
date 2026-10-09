// =============================================================================
//  apps/ui/PropertyPanel.h  ——  属性面板
//
//  回答一个问题："我刚点中的这个构件到底是什么、受力多大。"
//
//  【为什么属性要用"分组表"而不是一行行 label】
//  一根梁的属性有 20 多项（几何、截面、材料、内力、应力）。
//  平铺成一列时，找"这根的混凝土等级是多少"要扫过十几行数字；
//  分组之后（几何 / 截面 / 内力 / 验算）每一组只有 4~6 行，眼睛能跳着看。
//
//  【T3 可编辑化】坐标 / 截面 / 材料 / 墙厚 / 荷载值可以双击修改。
//  可编辑行通过 data(1, kFieldRole) 存 PropField 键；itemChanged 时
//  解析文本 → 写回 Model（内核 setter 自带事务回滚）→ 发 modelEdited()。
//  编辑期间 MainWindow 会刷新树/视口并按内存模型重算。
// =============================================================================
#pragma once

#include <QtGlobal>

#include <QWidget>

#include "SceneData.h"

class QTreeWidget;
class QTreeWidgetItem;
class QLabel;

namespace ui {

// 可编辑字段的键（存在 QTreeWidgetItem::data(1, kFieldRole)）。
// 【为什么用枚举而不是字符串】itemChanged 每敲一键都要分派一次，
// 字符串比较放在热路径里既慢又易写错拼写；用枚举 switch 分发代价为零。
enum class PropField : int {
  None = 0,
  // 节点坐标
  NodeX, NodeY, NodeZ,
  // 节点荷载（覆盖式写回，见 Model::setNodeForce/setNodeMoment）
  NodeFx, NodeFy, NodeFz,
  NodeMx, NodeMy, NodeMz,
  // 梁：截面（矩形等效宽/高）+ 材料
  BeamB, BeamH,
  BeamE, BeamNu, BeamGamma,
  // 壳：厚度（墙厚）+ E/ν/密度
  ShellT, ShellE, ShellNu, ShellDensity,
};
constexpr int kFieldRole = Qt::UserRole + 100;

class PropertyPanel : public QWidget {
  Q_OBJECT
 public:
  explicit PropertyPanel(QWidget* parent = nullptr);

  void setModel(Model* m);          // 非 const：编辑要写回模型
  void setResult(const StaticResult* r);
  void setPost(const post::PostProcessor* pp);

  void showElement(int elem, int type);
  void showNode(int node);
  void clear();

  // busy 期间禁用编辑（后台重算时 UI 线程不能写模型）
  void setEditEnabled(bool on);
  bool editEnabled() const { return editEnabled_; }

  // 当前选中对象（MainWindow 树刷新后恢复选中）
  int currentElem() const { return curElem_; }
  int currentType() const { return curType_; }
  int currentNode() const { return curNode_; }

 signals:
  // 一次成功的属性编辑（已写回 Model）。MainWindow 据此刷新树/视图并重算。
  void modelEdited();
  // 编辑被拒绝（格式错误 / 数值越界 / 写回失败）：面板已恢复原值，消息用于日志。
  void editRejected(const QString& why);

 private:
  void rebuild();
  bool applyEdit(PropField f, double v);
  void onItemChanged(QTreeWidgetItem* item, int column);

  QLabel* title_{nullptr};
  QTreeWidget* tree_{nullptr};

  Model* model_{nullptr};
  const StaticResult* result_{nullptr};
  const post::PostProcessor* post_{nullptr};

  int curElem_{-1};
  int curType_{0};
  int curNode_{-1};
  bool editEnabled_{true};
  bool rebuilding_{false};    // rebuild 期间 setText 也触发 itemChanged，用它防环
};

}  // namespace ui