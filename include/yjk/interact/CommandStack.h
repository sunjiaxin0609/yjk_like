// =============================================================================
//  yjk/interact/CommandStack.h  ——  撤销/重做命令栈（纯 C++，无 Qt 依赖）
//
//  【为什么用整模型快照而不是差分命令】
//    · Model::removeElement 是"按序号删除、后续前移"的语义，删一个单元后
//      所有后继 id 整体错位；属性编辑是"重建几何"的事务操作。要为每种
//      操作写差分 undo 既脆又容易漏（漏一条就把模型改成不可逆的坏状态）。
//    · 本项目的模型量级是几百节点/单元，clone() 深拷贝是 µs 级；
//      而 MainWindow 持有 std::unique_ptr<Model> —— 撤销/重做只需换指针，
//      O(1) 且绝不产生部分修改状态。
//    · 撤销栈 LRU 深度有限（默认 50 步），超出丢弃最老记录。
//
//  【语义】
//    push(before, after)         记录"一次已完成模型变更"（调用方负责先改
//                                再把前后快照交进来；本类不改模型）。
//    undo(cur)                   把 cur 指向的当前模型换回"操作前"状态，
//                                当前模型（= 操作后）转入重做栈。
//    redo(cur)                   对称：恢复"操作后"状态。
//    clear()                     新文件/清空模型时历史全部作废。
//
//  与 BuildOps 的关系：BuildOps 回答"拾取后建什么"，本类回答"建了如何
//  反悔"。两者都做成无 Qt 纯函数，GUI 与脚本层共用、可脱离界面单测。
// =============================================================================
#pragma once

#include <cstddef>
#include <deque>
#include <memory>

#include "yjk/model/Model.h"

namespace yjk {
namespace interact {

class CommandStack {
 public:
  explicit CommandStack(std::size_t maxDepth = 50);

  // 记录一次模型变更。before/after 各做一次深拷贝存入栈；
  // 栈满则丢弃最老记录（FIFO，保证栈深 ≤ maxDepth）。
  void push(const Model& before, const Model& after);

  bool canUndo() const { return !undo_.empty(); }
  bool canRedo() const { return !redo_.empty(); }
  std::size_t undoDepth() const { return undo_.size(); }
  std::size_t redoDepth() const { return redo_.size(); }
  std::size_t maxDepth() const { return maxDepth_; }

  // 撤销：cur 被替换为栈顶操作的"前"快照；原 cur 转入重做栈。
  // 无可撤销记录返回 false，cur 不变。
  bool undo(std::unique_ptr<Model>& cur);

  // 重做：cur 被替换为最近撤销操作的"后"快照。无可重做记录返回 false。
  bool redo(std::unique_ptr<Model>& cur);

  // 清空历史（打开新文件 / 新建工程时调用）。
  void clear();

 private:
  struct Snapshot {
    std::unique_ptr<Model> before;  // 操作前
    std::unique_ptr<Model> after;   // 操作后
  };
  std::deque<Snapshot> undo_;
  std::deque<Snapshot> redo_;
  std::size_t maxDepth_;
};

}  // namespace interact
}  // namespace yjk