// =============================================================================
//  apps/ui/ModelTreePanel.cpp
// =============================================================================
#include "ModelTreePanel.h"

#include <QFont>
#include <QHeaderView>
#include <QLabel>
#include <QVBoxLayout>

#include <map>

#include "Theme.h"
#include "yjk/element/BeamElement3D.h"

namespace ui {

namespace {

// 挂在树项上的数据。用 role 存"点击后视口要做什么"。
constexpr int kRoleKind = Qt::UserRole;       // 0 = 无，1 = 单元，2 = 节点
constexpr int kRoleId = Qt::UserRole + 1;
constexpr int kRoleType = Qt::UserRole + 2;   // 单元：1 梁 2 壳

QTreeWidgetItem* addChild(QTreeWidgetItem* parent, const QStringList& cols,
                          int kind, int id, int type) {
  auto* it = new QTreeWidgetItem(parent, cols);
  it->setData(0, kRoleKind, kind);
  it->setData(0, kRoleId, id);
  it->setData(0, kRoleType, type);
  return it;
}

QString beamKindName(const post::BeamResult* br, const Vec3& dir) {
  if (br) {
    switch (br->kind) {
      case post::MemberKind::Column: return QStringLiteral("柱");
      case post::MemberKind::Beam:   return QStringLiteral("梁");
      case post::MemberKind::Brace:  return QStringLiteral("撑");
      default: break;
    }
  }
  // 没有后处理时按几何判断：杆轴接近竖直就是柱
  return std::abs(dir.z) > 0.7 ? QStringLiteral("柱") : QStringLiteral("梁");
}

}  // namespace

ModelTreePanel::ModelTreePanel(QWidget* parent) : QWidget(parent) {
  auto* lay = new QVBoxLayout(this);
  lay->setContentsMargins(0, 0, 0, 0);
  lay->setSpacing(0);

  tree_ = new QTreeWidget(this);
  tree_->setHeaderLabels({QStringLiteral("对象"), QStringLiteral("信息")});
  tree_->header()->setStretchLastSection(true);
  tree_->header()->resizeSection(0, 190);
  tree_->setAlternatingRowColors(true);
  tree_->setUniformRowHeights(true);
  tree_->setFont(uiFont(9));
  lay->addWidget(tree_);

  connect(tree_, &QTreeWidget::itemClicked, this, [this](QTreeWidgetItem* it, int) {
    if (block_ || !it) return;
    const int kind = it->data(0, kRoleKind).toInt();
    const int id = it->data(0, kRoleId).toInt();
    const int type = it->data(0, kRoleType).toInt();
    if (kind == 2) emit nodeActivated(id);
    else if (kind == 1) emit elementActivated(id, type);
  });
}

void ModelTreePanel::setModel(const Model* m) {
  model_ = m;
  rebuild();
}

void ModelTreePanel::setPost(const post::PostProcessor* pp) {
  post_ = pp;
  rebuild();
}

void ModelTreePanel::rebuild() {
  block_ = true;
  tree_->clear();
  if (!model_) { block_ = false; return; }

  const Id ne = model_->elementCount();
  const Id nn = model_->nodeCount();

  // ---- 概览 ----
  auto* rootModel = new QTreeWidgetItem(tree_, {QStringLiteral("模型"), QString()});
  rootModel->setExpanded(true);
  addChild(rootModel, {QStringLiteral("节点"), QString::number(nn)}, 0, -1, 0);
  addChild(rootModel, {QStringLiteral("单元"), QString::number(ne)}, 0, -1, 0);

  // 楼层是从节点属性里统计出来的（模型层没有单独的"楼层表"）
  std::map<int, int> storyNodes, storyElems;
  for (const Node& nd : model_->nodes()) storyNodes[nd.story]++;
  for (Id e = 0; e < ne; ++e) {
    const auto& el = *model_->elements()[static_cast<size_t>(e)];
    int top = -1000000;
    for (Id n : el.nodes()) top = std::max(top, model_->node(n).story);
    storyElems[top]++;
  }
  const int nStory = static_cast<int>(storyNodes.size());
  addChild(rootModel, {QStringLiteral("楼层数"), QString::number(nStory)}, 0, -1, 0);

  // ---- 材料与截面 ----
  //
  // 【为什么要在这里列出，而不是只放在属性面板里】
  // 截面用错是建模阶段最贵的一类错误 —— 构件全建对了、荷载也对，
  // 但用了 200×400 而不是 300×600，结果差一倍且没有任何报错。
  // 在树上把每个截面各用了几根杆摆出来，"只用了 1 根的那个截面"
  // 一眼就能发现是误建。
  auto* rootMat = new QTreeWidgetItem(tree_, {QStringLiteral("材料与截面"), QString()});
  std::map<std::string, int> matUse;
  std::map<std::string, int> secUse;
  for (Id e = 0; e < ne; ++e) {
    const auto& el = *model_->elements()[static_cast<size_t>(e)];
    if (el.type() != ElementType::Beam3D) continue;
    const auto& be = static_cast<const BeamElement&>(el);
    const auto& sec = be.section();
    const auto& mat = be.material();
    // 用 A、Iy、Iz 三元组当截面指纹（SectionProperties 没有名字字段）
    char buf[160];
    std::snprintf(buf, sizeof(buf), "A=%.4g m²  Iy=%.4g  Iz=%.4g", sec.A, sec.Iy, sec.Iz);
    secUse[buf]++;
    matUse[mat.name.empty() ? std::string("(未命名)") : mat.name]++;
  }
  for (const auto& kv : matUse) {
    addChild(rootMat, {QStringLiteral("材料 "), QString::fromStdString(kv.first) +
                       QStringLiteral(" · %1 根杆").arg(kv.second)}, 0, -1, 0);
  }
  for (const auto& kv : secUse) {
    addChild(rootMat, {QStringLiteral("截面 "), QString::fromStdString(kv.first) +
                       QStringLiteral(" · %1 根").arg(kv.second)}, 0, -1, 0);
  }

// ---- 构件（按 类别 → 楼层 分组）----
  //
  // 【为什么按类别分组而不是按编号平铺】
  // 一个 4 层 3×3 框架就有 87 个单元，平铺出来是一堵墙，找不到东西。
  // 按"柱/梁/板/墙"分组后，"所有柱"是一个可整体检查的集合，
  // 这是结构工程师核对模型的实际方式。
  //
  // 【为什么墙和板分开两组】
  // 墙壳与楼板壳在模型树里必须能一眼区分：isWall() 的壳是抗侧力竖构件
  // （电梯井、核心筒），承担楼层剪力的主体；水平板壳只扛竖向荷载。
  // 混在一组里，选中一片电梯井墙却显示在"楼板"分组下，会误导画图顺序。
  auto* rootCol = new QTreeWidgetItem(tree_, {QStringLiteral("柱"), QString()});
  auto* rootBeam = new QTreeWidgetItem(tree_, {QStringLiteral("梁"), QString()});
  auto* rootSlab = new QTreeWidgetItem(tree_, {QStringLiteral("楼板"), QString()});
  auto* rootWall = new QTreeWidgetItem(tree_, {QStringLiteral("剪力墙"), QString()});

  std::map<int, QTreeWidgetItem*> colByStory, beamByStory, slabByStory, wallByStory;
  int nCol = 0, nBeam = 0, nSlab = 0, nWall = 0;

  auto groupOf = [&](QTreeWidgetItem* root, std::map<int, QTreeWidgetItem*>& byStory, int story) {
    auto it = byStory.find(story);
    if (it != byStory.end()) return it->second;
    auto* g = new QTreeWidgetItem(root, {QStringLiteral("%1 层").arg(story), QString()});
    g->setExpanded(story >= nStory - 1);       // 只展开顶层，避免一打开就是几百行
    byStory[story] = g;
    return g;
  };

  for (Id e = 0; e < ne; ++e) {
    const auto& el = *model_->elements()[static_cast<size_t>(e)];
    if (el.type() == ElementType::Beam3D) {
      const auto& be = static_cast<const BeamElement&>(el);
      const auto& dir = be.core().localX();
      const post::BeamResult* br = nullptr;
      if (post_ && e < static_cast<Id>(post_->beamResults().size()))
        br = &post_->beamResults()[static_cast<size_t>(e)];
      const QString kind = beamKindName(br, dir);
      const int sa = model_->node(el.nodes()[0]).story;
      const int sb = model_->node(el.nodes()[1]).story;
      const int top = std::max(sa, sb);
      const QString info = QStringLiteral("#%1  L=%2m").arg(e).arg(fmtNum(be.length(), 3));
      if (kind == QStringLiteral("柱")) {
        addChild(groupOf(rootCol, colByStory, top), {info, kind}, 1, static_cast<int>(e), 1);
        ++nCol;
      } else {
        addChild(groupOf(rootBeam, beamByStory, top), {info, kind}, 1, static_cast<int>(e), 1);
        ++nBeam;
      }
} else {
      // 壳单元：isWall() 的竖墙进"剪力墙"分组，否则是水平楼板。
      // 二者同为 Shell4 但承载角色完全不同（抗侧 vs 竖向传荷），
      // 树上的分组与名称必须反映这一点。
      const auto& sh = static_cast<const ShellElement&>(el);
      const int top = model_->node(el.nodes()[0]).story;
      const QString info = QStringLiteral("#%1  %2m²").arg(e).arg(fmtNum(sh.area(), 3));
      if (sh.isWall()) {
        addChild(groupOf(rootWall, wallByStory, top), {info, QStringLiteral("墙")},
                 1, static_cast<int>(e), 2);
        ++nWall;
      } else {
        addChild(groupOf(rootSlab, slabByStory, top), {info, QStringLiteral("板")},
                 1, static_cast<int>(e), 2);
        ++nSlab;
      }
    }
  }

  rootCol->setText(1, QString::number(nCol));
  rootBeam->setText(1, QString::number(nBeam));
  rootSlab->setText(1, QString::number(nSlab));
  rootWall->setText(1, QString::number(nWall));
  // 数量少的类别默认展开，数量多的收起 —— 打开面板第一眼要能看到全貌
  rootCol->setExpanded(nCol <= 30);
  rootBeam->setExpanded(nBeam <= 30);
  rootSlab->setExpanded(false);
  rootWall->setExpanded(nWall <= 30);

  block_ = false;
}

void ModelTreePanel::selectElement(int elem, int type) {
  block_ = true;
  QTreeWidgetItemIterator it(tree_);
  for (; *it; ++it) {
    QTreeWidgetItem* item = *it;
    if (item->data(0, kRoleKind).toInt() == 1 &&
        item->data(0, kRoleId).toInt() == elem &&
        item->data(0, kRoleType).toInt() == type) {
      // 展开所有祖先，否则选中项藏在收起的节点里，看起来"没反应"
      for (QTreeWidgetItem* p = item->parent(); p; p = p->parent()) p->setExpanded(true);
      tree_->setCurrentItem(item);
      tree_->scrollToItem(item);
      break;
    }
  }
  block_ = false;
}

}  // namespace ui
