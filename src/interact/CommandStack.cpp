// =============================================================================
//  src/interact/CommandStack.cpp  ——  撤销/重做命令栈实现
// =============================================================================
#include "yjk/interact/CommandStack.h"

namespace yjk {
namespace interact {

CommandStack::CommandStack(std::size_t maxDepth)
    : maxDepth_(maxDepth > 0 ? maxDepth : 1) {}

void CommandStack::push(const Model& before, const Model& after) {
  Snapshot s;
  s.before = std::make_unique<Model>(before.clone());
  s.after = std::make_unique<Model>(after.clone());
  undo_.push_back(std::move(s));
  // 栈深有限：超出丢弃最老记录（FIFO）。新操作使重做历史失效。
  if (undo_.size() > maxDepth_) undo_.pop_front();
  redo_.clear();
}

bool CommandStack::undo(std::unique_ptr<Model>& cur) {
  if (undo_.empty()) return false;
  Snapshot s = std::move(undo_.back());
  undo_.pop_back();
  // 当前模型即"该操作后状态"：连同快照中的 after 一起转入重做栈，
  // 保证 redo 后 undo 历史仍然正确（redo 可用后再次 undo 能回到 before）。
  Snapshot r;
  r.before = std::move(cur);
  r.after = std::move(s.after);
  redo_.push_back(std::move(r));
  cur = std::move(s.before);
  return true;
}

bool CommandStack::redo(std::unique_ptr<Model>& cur) {
  if (redo_.empty()) return false;
  Snapshot s = std::move(redo_.back());
  redo_.pop_back();
  // 对称：当前模型是"操作前状态"，连同快照中的 after 压回撤销栈。
  Snapshot u;
  u.before = std::move(cur);
  u.after = std::move(s.after);
  undo_.push_back(std::move(u));
  cur = std::move(s.before);
  return true;
}

void CommandStack::clear() {
  undo_.clear();
  redo_.clear();
}

}  // namespace interact
}  // namespace yjk