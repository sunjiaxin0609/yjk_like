// =============================================================================
//  apps/ui/MainWindow.cpp
// =============================================================================
#include "MainWindow.h"

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QImage>
#include <QInputDialog>
#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QStatusBar>
#include <QThread>
#include <QToolBar>
#include <QVBoxLayout>

#include <QElapsedTimer>
#include <QTime>

#include <QMimeData>

#include <algorithm>
#include <cmath>
#include <utility>

#include "ModelTreePanel.h"
#include "PropertyPanel.h"
#include "ResultPanel.h"
#include "Theme.h"
#include "View3D.h"
#include "yjk/interact/BuildOps.h"
#include "yjk/post/ResultExport.h"

namespace ui {

namespace sec = yjk::section;
using yjk::Material;
using yjk::SectionProperties;
using yjk::ShellProperties;

namespace {

// -----------------------------------------------------------------------------
//  后台计算的结果包
//
//  【为什么要打包整体搬运】PostProcessor 内部持有 Model 与 StaticResult 的
//  【引用】。如果先搬 Model、再搬 Result，中间那一刻 post 的引用指向的是
//  已经失效的旧对象。三者必须一起构造、一起移动，且都放在堆上
//  （unique_ptr 的移动不会移动对象本身，引用始终有效）。
// -----------------------------------------------------------------------------
struct Bundle {
  std::unique_ptr<Model> model;
  std::unique_ptr<StaticResult> result;
  std::unique_ptr<post::PostProcessor> post;

  std::vector<std::string> errors, warnings;
  std::string info;          // 建模摘要
  std::string diagErrors;    // 模型校核（多条合并成文本）
  std::string diagWarnings;
  std::string solveMessage;
  int nStories{0};

  bool parsed{false};
  bool buildOk{false};
  bool checkOk{false};
  bool solveOk{false};

  double buildMs{0.0}, solveMs{0.0};
};

}  // namespace

// =============================================================================
MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
  setWindowTitle(QStringLiteral("yjk_like  ——  结构分析"));
  resize(1480, 900);

view_ = new View3D(this);
  setCentralWidget(view_);

  cmds_ = std::make_unique<yjk::interact::CommandStack>(50);   // T5 撤销/重做栈

  buildActions();
  buildMenus();
  buildToolBar();
  buildDocks();
  buildStatusBar();

  // ---- 视口 ↔ 面板 的双向联动 ----
  //
  // 【必须双向且防环】视口点选 → 树上高亮；树上点选 → 视口高亮。
  // 两条链路都会触发"选中变化"，不加防护会无限递归。
  // 这里靠各自的 block 标志断开：视口发信号时树只更新显示不回调，
  // 反之亦然（ModelTreePanel 内部已用 block_ 处理）。
  connect(view_, &View3D::elementPicked, this, [this](int e, int t) {
    prop_->showElement(e, t);
    tree_->selectElement(e, t);
  });
  connect(view_, &View3D::nodePicked, this, [this](int n) { prop_->showNode(n); });
connect(view_, &View3D::selectionChanged, this, [this] {
    if (view_->selection().empty()) prop_->clear();
    // T5：有选中对象才允许 Delete
    if (actDelete_)
      actDelete_->setEnabled(!view_->selection().empty() || view_->selectedNode() >= 0);
  });
  connect(view_, &View3D::viewChanged, this, [this] { syncOptionsToView(); });

  connect(tree_, &ModelTreePanel::elementActivated, this, [this](int e, int t) {
    view_->selectElement(e, t);
    prop_->showElement(e, t);
  });
  connect(tree_, &ModelTreePanel::nodeActivated, this, [this](int n) {
    prop_->showNode(n);
    view_->selectElement(-1, 0);
  });

connect(results_, &ResultPanel::elementActivated, this, [this](int e, int t) {
    view_->selectElement(e, t);
    prop_->showElement(e, t);
    tree_->selectElement(e, t);
    setWindowTitle(QStringLiteral("yjk_like  ——  结构分析  [定位到单元 #%1]").arg(e));
  });

  // ---- T3 属性编辑链路 ----
  // 属性面板双击改参 -> 写回 Model（面板内部已做事务回滚）
  // -> modelEdited：结果失效 + 刷新视图/树 + 内存重算
  // -> editRejected：面板已恢复原值，这里只负责记日志
connect(prop_, &PropertyPanel::modelEdited, this, &MainWindow::onModelEdited);
  connect(prop_, &PropertyPanel::editRejected, this,
          [this](const QString& why) { log(why, true); });

  // ---- T4 创建工具：视口事件 → 状态机 ----
  connect(view_, &View3D::buildNodeHit, this, &MainWindow::onBuildNodeHit);
  connect(view_, &View3D::buildBoxNodes, this, &MainWindow::onBuildBoxNodes);
  connect(view_, &View3D::buildEscaped, this, &MainWindow::onBuildEscaped);
  connect(view_, &View3D::buildToolChanged, this, &MainWindow::onBuildToolChanged);

  log(QStringLiteral("就绪。打开一个 .yjk 建模脚本，或直接点「开始计算」用内置示例试跑。"));
  log(QStringLiteral("创建工具：B 梁 · C 柱 · W 墙 · P 板 · L 节点荷载（Esc 退出）"));
  setBusy(false);
}

MainWindow::~MainWindow() = default;

// =============================================================================
//  动作 / 菜单 / 工具栏
// =============================================================================
void MainWindow::buildActions() {
  actOpen_ = new QAction(QStringLiteral("打开建模脚本..."), this);
  actOpen_->setShortcut(QKeySequence::Open);
  connect(actOpen_, &QAction::triggered, this, &MainWindow::onOpen);

  actReload_ = new QAction(QStringLiteral("重新加载"), this);
  actReload_->setShortcut(QKeySequence::Refresh);
  connect(actReload_, &QAction::triggered, this, &MainWindow::onReload);

  actCheck_ = new QAction(QStringLiteral("模型校核"), this);
  connect(actCheck_, &QAction::triggered, this, &MainWindow::onCheck);

  actRun_ = new QAction(QStringLiteral("开始计算"), this);
  actRun_->setShortcut(Qt::Key_F5);
  connect(actRun_, &QAction::triggered, this, &MainWindow::onRun);

  actSolid_ = new QAction(QStringLiteral("实体显示"), this);
  actSolid_->setCheckable(true);
  actSolid_->setChecked(true);

  actDeform_ = new QAction(QStringLiteral("变形显示"), this);
  actDeform_->setCheckable(true);

  actAnimate_ = new QAction(QStringLiteral("变形动画"), this);
  actAnimate_->setCheckable(true);

  actOrtho_ = new QAction(QStringLiteral("正交投影"), this);
  actOrtho_->setCheckable(true);

  actNodes_ = new QAction(QStringLiteral("显示节点"), this);
  actNodes_->setCheckable(true);
  actSupports_ = new QAction(QStringLiteral("显示支座"), this);
  actSupports_->setCheckable(true);
  actSupports_->setChecked(true);
  actGrid_ = new QAction(QStringLiteral("显示参考网格"), this);
  actGrid_->setCheckable(true);
  actGrid_->setChecked(true);

actContour_ = new QAction(QStringLiteral("显示云图"), this);
  actContour_->setCheckable(true);

  // ---- T4 创建工具（B/C/W/P/L，WidgetShortcut → 输入框内不误触）----
  auto mkTool = [this](const QString& text, const char* key, const QString& tip,
                       View3D::BuildTool t) {
    QAction* a = new QAction(text, view_);   // 归 view_ 所有：快捷键只在视口有焦点时生效
    a->setCheckable(true);
    a->setShortcut(QKeySequence(QString::fromUtf8(key)));
    a->setShortcutContext(Qt::WidgetShortcut);
    a->setToolTip(tip);
    connect(a, &QAction::triggered, this,
            [this, t] { onToolTriggered(t); });
    return a;
  };
  actToolBeam_ = mkTool(QStringLiteral("画梁 (B)"), "B",
                        QStringLiteral("两点定梁：依次点选起点、终点节点"),
                        View3D::BuildTool::Beam);
  actToolColumn_ = mkTool(QStringLiteral("画柱 (C)"), "C",
                          QStringLiteral("两点定柱：依次点选底、顶节点"),
                          View3D::BuildTool::Column);
  actToolWall_ = mkTool(QStringLiteral("画墙 (W)"), "W",
                        QStringLiteral("四角定墙：底左→底右→顶右→顶左"),
                        View3D::BuildTool::Wall);
  actToolSlab_ = mkTool(QStringLiteral("画板 (P)"), "P",
                        QStringLiteral("框选定板：拖拽框选板所在区域的角点"),
                        View3D::BuildTool::Slab);
  actToolLoad_ = mkTool(QStringLiteral("节点荷载 (L)"), "L",
                        QStringLiteral("点击节点施加集中力/力矩"),
                        View3D::BuildTool::Load);
  toolGroup_ = new QActionGroup(this);
  toolGroup_->setExclusive(true);
  for (QAction* a : {actToolBeam_, actToolColumn_, actToolWall_, actToolSlab_,
                     actToolLoad_})
    toolGroup_->addAction(a);

  // ---- T5 撤销 / 重做 / 删除 ----
  // 撤销/重做是全局编辑动作（WindowShortcut）；删除归 view_ 且 WidgetShortcut：
  // 只在视口有焦点时生效，避免属性面板输入框里按 Delete 误删整个选中对象。
  actUndo_ = new QAction(QStringLiteral("撤销"), this);
  actUndo_->setShortcut(QKeySequence::Undo);
  actUndo_->setEnabled(false);
  connect(actUndo_, &QAction::triggered, this, &MainWindow::onUndo);

  actRedo_ = new QAction(QStringLiteral("重做"), this);
  actRedo_->setShortcut(QKeySequence::Redo);
  actRedo_->setEnabled(false);
  connect(actRedo_, &QAction::triggered, this, &MainWindow::onRedo);

  actDelete_ = new QAction(QStringLiteral("删除选中"), view_);
  actDelete_->setShortcut(QKeySequence::Delete);
  actDelete_->setShortcutContext(Qt::WidgetShortcut);
  actDelete_->setEnabled(false);
  actDelete_->setToolTip(QStringLiteral("删除选中的构件/节点（可撤销）"));
  connect(actDelete_, &QAction::triggered, this, &MainWindow::onDeleteSelected);
}

void MainWindow::buildMenus() {
  // ---- 文件 ----
  auto* file = menuBar()->addMenu(QStringLiteral("文件(&F)"));
  file->addAction(actOpen_);
  file->addAction(actReload_);
  file->addSeparator();
  file->addAction(QStringLiteral("导出结果 JSON..."), this, &MainWindow::onExportJson);
  file->addAction(QStringLiteral("导出可视化 HTML..."), this, &MainWindow::onExportHtml);
  file->addAction(QStringLiteral("导出文本报告..."), this, &MainWindow::onExportReport);
  file->addAction(QStringLiteral("导出视图截图..."), this, &MainWindow::onExportShot);
file->addSeparator();
  file->addAction(QStringLiteral("退出"), this, &QWidget::close);

  // ---- 编辑（T5 撤销/重做/删除）----
  auto* edit = menuBar()->addMenu(QStringLiteral("编辑(&E)"));
  edit->addAction(actUndo_);
  edit->addAction(actRedo_);
  edit->addSeparator();
  edit->addAction(actDelete_);

  // ---- 视图 ----
  auto* view = menuBar()->addMenu(QStringLiteral("视图(&V)"));
  view->addAction(actSolid_);
  view->addAction(actDeform_);
  view->addAction(actAnimate_);
  view->addSeparator();
  view->addAction(actContour_);
  view->addSeparator();
  view->addAction(actNodes_);
  view->addAction(actSupports_);
  view->addAction(actGrid_);
  view->addSeparator();
  view->addAction(actOrtho_);
  view->addSeparator();

  auto* presets = view->addMenu(QStringLiteral("预设视角"));
  struct P { const char* n; Camera::Preset p; QKeySequence k; };
  const P ps[] = {
      {"轴测图", Camera::Preset::Iso, QKeySequence(Qt::Key_1)},
      {"正立面", Camera::Preset::Front, QKeySequence(Qt::Key_2)},
      {"侧立面", Camera::Preset::Left, QKeySequence(Qt::Key_3)},
      {"俯视图（平面）", Camera::Preset::Top, QKeySequence(Qt::Key_4)},
      {"背立面", Camera::Preset::Back, QKeySequence()},
      {"右立面", Camera::Preset::Right, QKeySequence()},
      {"仰视图", Camera::Preset::Bottom, QKeySequence()},
  };
  for (const P& p : ps) {
    auto* a = presets->addAction(QString::fromUtf8(p.n));
    if (!p.k.isEmpty()) a->setShortcut(p.k);
    const Camera::Preset pr = p.p;
    connect(a, &QAction::triggered, this, [this, pr] {
      view_->setPreset(pr);
      if (cmbPreset_) cmbPreset_->setCurrentIndex(0);
    });
  }
  presets->addSeparator();
  auto* aFit = presets->addAction(QStringLiteral("全览 (Home)"));
  aFit->setShortcut(Qt::Key_Home);
  connect(aFit, &QAction::triggered, this, [this] { view_->frameAll(); });

// ---- 建模（T4 创建工具）----
  auto* build = menuBar()->addMenu(QStringLiteral("建模(&M)"));
  build->addAction(actToolBeam_);
  build->addAction(actToolColumn_);
  build->addAction(actToolWall_);
  build->addAction(actToolSlab_);
  build->addAction(actToolLoad_);
  build->addSeparator();
  build->addAction(QStringLiteral("退出创建工具 (Esc)"), this, [this] {
    view_->setBuildTool(View3D::BuildTool::None);
  });

  // 视口右键菜单 = 同一组创建工具（常用操作少点几下）
  view_->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(view_, &QWidget::customContextMenuRequested, this,
          [this](const QPoint& pos) {
            QMenu menu(this);
            menu.addAction(actToolBeam_);
            menu.addAction(actToolColumn_);
            menu.addAction(actToolWall_);
            menu.addAction(actToolSlab_);
            menu.addAction(actToolLoad_);
            menu.addSeparator();
            menu.addAction(QStringLiteral("退出创建工具 (Esc)"), this, [this] {
              view_->setBuildTool(View3D::BuildTool::None);
            });
            menu.exec(view_->mapToGlobal(pos));
          });

  // ---- 计算 ----
  auto* calc = menuBar()->addMenu(QStringLiteral("计算(&C)"));
  calc->addAction(actCheck_);
  calc->addAction(actRun_);

  // ---- 帮助 ----
  auto* help = menuBar()->addMenu(QStringLiteral("帮助(&H)"));
  help->addAction(QStringLiteral("快捷键与操作说明"), this, [this] {
    QMessageBox::information(this, QStringLiteral("操作说明"),
        QStringLiteral(
            "【三维视图操作】\n"
            "  左键拖动          旋转（轨道）\n"
            "  中键 / Shift+左键 平移\n"
            "  滚轮              以光标为锚点缩放\n"
            "  左键单击          选中构件 / 节点\n"
            "  单击空白          取消选择\n\n"
            "【快捷键】\n"
            "  1 / 2 / 3 / 4     轴测 / 正立面 / 侧立面 / 俯视\n"
            "  Home 或 F         缩放到全览\n"
            "  O                 正交 / 透视切换\n"
            "  空格              变形动画 开 / 关\n"
            "  Esc               取消选择\n"
            "  F5                开始计算\n\n"
            "【怎么看结果】\n"
            "  云图范围【必然包含 0】：色标中点对应 0 值，\n"
            "  蓝色一侧为负（压），红色一侧为正（拉）。"));
  });
  help->addAction(QStringLiteral("关于"), this, &MainWindow::onAbout);
}

void MainWindow::buildToolBar() {
  auto* tb = addToolBar(QStringLiteral("主工具栏"));
  tb->setMovable(false);
  tb->setFont(uiFont(9));
  tb->setIconSize(QSize(18, 18));

  tb->addAction(actOpen_);
  tb->addAction(actReload_);
  tb->addSeparator();
  tb->addAction(actCheck_);
  tb->addAction(actRun_);
  tb->addSeparator();
tb->addAction(actSolid_);
  tb->addAction(actDeform_);
  tb->addAction(actAnimate_);
  tb->addSeparator();

  // T4 创建工具
  tb->addWidget(new QLabel(QStringLiteral(" 创建 "), tb));
  tb->addAction(actToolBeam_);
  tb->addAction(actToolColumn_);
  tb->addAction(actToolWall_);
  tb->addAction(actToolSlab_);
  tb->addAction(actToolLoad_);
  tb->addSeparator();

  // 云图变量
  tb->addWidget(new QLabel(QStringLiteral(" 云图 "), tb));
  actContour_->setChecked(false);
  tb->addAction(actContour_);
  cmbField_ = new QComboBox(tb);
  cmbField_->setFont(uiFont(9));
  cmbField_->setMinimumWidth(190);
  cmbField_->setToolTip(QStringLiteral("节点级云图（壳面按节点值插值）"));
  tb->addWidget(cmbField_);

  cmbMemberField_ = new QComboBox(tb);
  cmbMemberField_->setFont(uiFont(9));
  cmbMemberField_->setMinimumWidth(120);
  cmbMemberField_->setToolTip(QStringLiteral("构件级云图（整根杆一个颜色）"));
  cmbMemberField_->addItem(QStringLiteral("（无）"), static_cast<int>(MemberField::None));
  cmbMemberField_->addItem(QStringLiteral("弯矩 M"), static_cast<int>(MemberField::Moment));
  cmbMemberField_->addItem(QStringLiteral("剪力 V"), static_cast<int>(MemberField::Shear));
  cmbMemberField_->addItem(QStringLiteral("轴力 N"), static_cast<int>(MemberField::Axial));
  cmbMemberField_->addItem(QStringLiteral("正应力 σ"), static_cast<int>(MemberField::Stress));
  cmbMemberField_->addItem(QStringLiteral("应力比"), static_cast<int>(MemberField::Utilization));
  tb->addWidget(cmbMemberField_);

  tb->addSeparator();
  tb->addWidget(new QLabel(QStringLiteral(" 视角 "), tb));
  cmbPreset_ = new QComboBox(tb);
  cmbPreset_->setFont(uiFont(9));
  cmbPreset_->addItems({QStringLiteral("轴测"), QStringLiteral("正立面"),
                        QStringLiteral("侧立面"), QStringLiteral("俯视")});
  cmbPreset_->setCurrentIndex(-1);
  connect(cmbPreset_, QOverload<int>::of(&QComboBox::activated), this, [this](int i) {
    const Camera::Preset ps[4] = {Camera::Preset::Iso, Camera::Preset::Front,
                                  Camera::Preset::Left, Camera::Preset::Top};
    if (i >= 0 && i < 4) view_->setPreset(ps[i]);
  });
  tb->addWidget(cmbPreset_);
  tb->addAction(actOrtho_);

  // 视图选项变化统一走一个槽，避免每个动作都写一遍
  const QList<QAction*> opts = {actSolid_, actDeform_, actAnimate_, actOrtho_,
                                actNodes_, actSupports_, actGrid_, actContour_};
  for (QAction* a : opts)
    connect(a, &QAction::toggled, this, &MainWindow::onViewOptionsChanged);
  connect(cmbField_, QOverload<int>::of(&QComboBox::activated), this,
          &MainWindow::onViewOptionsChanged);
  connect(cmbMemberField_, QOverload<int>::of(&QComboBox::activated), this,
          &MainWindow::onViewOptionsChanged);
}

void MainWindow::buildDocks() {
  // ---- 左：模型树 ----
  auto* dTree = new QDockWidget(QStringLiteral("模型导航"), this);
  dTree->setObjectName(QStringLiteral("dockTree"));
  tree_ = new ModelTreePanel(dTree);
  dTree->setWidget(tree_);
  addDockWidget(Qt::LeftDockWidgetArea, dTree);

  // ---- 右：属性 + 结果 + 消息（用 dock 竖排）----
  auto* dProp = new QDockWidget(QStringLiteral("属性"), this);
  dProp->setObjectName(QStringLiteral("dockProp"));
  prop_ = new PropertyPanel(dProp);
  dProp->setWidget(prop_);
  addDockWidget(Qt::RightDockWidgetArea, dProp);

  auto* dRes = new QDockWidget(QStringLiteral("计算结果"), this);
  dRes->setObjectName(QStringLiteral("dockResult"));
  results_ = new ResultPanel(dRes);
  dRes->setWidget(results_);
  addDockWidget(Qt::RightDockWidgetArea, dRes);

  auto* dLog = new QDockWidget(QStringLiteral("消息"), this);
  dLog->setObjectName(QStringLiteral("dockLog"));
  log_ = new QPlainTextEdit(dLog);
  log_->setReadOnly(true);
  log_->setFont(monoFont(9));
  log_->setMaximumBlockCount(2000);   // 长会话不无限吃内存
  dLog->setWidget(log_);
  addDockWidget(Qt::BottomDockWidgetArea, dLog);

  // 横向：结果面板给得宽一点 —— 楼层指标表有 9 列，窄了全是省略号。
  resizeDocks({dTree, dProp, dRes}, {240, 320, 500}, Qt::Horizontal);
  // 纵向：属性面板要高。一根梁的属性有 20 多项，矮了连"材料"那一组
  // 都露不出来，用户会以为属性面板是空的。
  resizeDocks({dProp, dRes}, {560, 380}, Qt::Vertical);
  resizeDocks({dLog}, {140}, Qt::Vertical);
}

void MainWindow::buildStatusBar() {
  lblModel_ = new QLabel(QStringLiteral("未加载模型"), this);
  lblResult_ = new QLabel(QStringLiteral("—"), this);
  progress_ = new QProgressBar(this);
  progress_->setMaximumWidth(160);
  progress_->setVisible(false);
  progress_->setRange(0, 0);          // 不确定进度：求解耗时无法预估

  lblModel_->setFont(uiFont(9));
  lblResult_->setFont(monoFont(9));

  statusBar()->addWidget(lblModel_, 1);
  statusBar()->addPermanentWidget(lblResult_);
  statusBar()->addPermanentWidget(progress_);
}

// =============================================================================
//  视图选项
// =============================================================================
void MainWindow::onViewOptionsChanged() {
  RenderOptions o = view_->options();
  o.solid = actSolid_->isChecked();
  o.showNodes = actNodes_->isChecked();
  o.showSupports = actSupports_->isChecked();
  o.showGrid = actGrid_->isChecked();
  o.deformed = actDeform_->isChecked();
  o.contourOn = actContour_->isChecked();
  // 【不能在 index = -1 时读 currentData()】
  // 此时返回的是无效 QVariant，toInt() 静默给 0 —— 与"真的选中了第 0 个场"
  // 完全无法区分。后果是"撤销选择节点场"时会连带把构件云图也关掉。
  o.nodeField = (cmbField_->currentIndex() >= 0 && cmbField_->currentData().isValid())
                    ? cmbField_->currentData().toInt()
                    : -1;
  o.memberField = static_cast<MemberField>(
      cmbMemberField_->currentIndex() >= 0 ? cmbMemberField_->currentData().toInt() : 0);

  // 节点场与构件场是互斥的两种上色方式。
  // 【必须显式互斥】两个同时勾上时，代码里"节点场优先"，
  // 用户切了构件场却看不到变化，会以为功能坏了。
  if (o.nodeField >= 0 && o.memberField != MemberField::None) {
    o.memberField = MemberField::None;
    cmbMemberField_->setCurrentIndex(0);
  }
  // 没有后处理结果时不该允许开云图
  if (o.contourOn && !post_) {
    o.contourOn = false;
    actContour_->setChecked(false);
    log(QStringLiteral("还没有计算结果，云图不可用。"), true);
  }

  view_->setOptions(o);
  if (actAnimate_->isChecked() != view_->animating())
    view_->setAnimating(actAnimate_->isChecked());
  if (actOrtho_->isChecked() != view_->ortho()) view_->setOrtho(actOrtho_->isChecked());
}

void MainWindow::syncOptionsToView() {
  // 视口从内部改了状态（例如按空格切了动画），把工具栏同步过来。
  // 【不加防护会死循环】必须用 blockSignals。
  const bool an = view_->animating();
  const bool ortho = view_->ortho();
  if (actAnimate_->isChecked() != an) {
    actAnimate_->blockSignals(true);
    actAnimate_->setChecked(an);
    actAnimate_->blockSignals(false);
  }
  if (actOrtho_->isChecked() != ortho) {
    actOrtho_->blockSignals(true);
    actOrtho_->setChecked(ortho);
    actOrtho_->blockSignals(false);
  }
}

// =============================================================================
//  数据发布
// =============================================================================
void MainWindow::publishModel() {
  view_->setModel(model_.get());
  tree_->setModel(model_.get());
  prop_->setModel(model_.get());
  results_->setModel(model_.get());

  if (model_) {
    // 楼层数从节点属性统计（模型层没有单独的楼层表）
    int mx = 0;
    for (const Node& nd : model_->nodes()) mx = std::max(mx, nd.story);
    lblModel_->setText(QStringLiteral("节点 %1 · 单元 %2 · 楼层 %3  ·  %4")
                           .arg(model_->nodeCount())
                           .arg(model_->elementCount())
                           .arg(mx + 1)
                           .arg(scriptPath_.isEmpty() ? QStringLiteral("(未保存脚本)")
                                                      : QFileInfo(scriptPath_).fileName()));
  } else {
    lblModel_->setText(QStringLiteral("未加载模型"));
  }
}

void MainWindow::publishResults(bool autostyle) {
  view_->setResult(result_.get());
  view_->setPost(post_.get());
  prop_->setResult(result_.get());
  prop_->setPost(post_.get());
  results_->setResult(result_.get());
  results_->setPost(post_.get());

  refreshFieldCombo();

  // 自动切到"变形 + 位移云图"：
  // 算完之后用户的第一件事就是看变形，让他们再点四下菜单是浪费。
  // 【autostyle=false 用于属性编辑触发的内存重算】编辑过程中用户视角
  // 正盯着某一根梁，此时自动切云图 + frameAll 会把视口拽走，
  // 打断"改一处 → 看一处"的工作流。
  if (autostyle && result_ && result_->ok) {
    actDeform_->setChecked(true);
    actContour_->setChecked(true);
    if (cmbField_->count() > 0) cmbField_->setCurrentIndex(0);
    actSolid_->setChecked(true);
    onViewOptionsChanged();
    view_->frameAll();
  }

  if (result_ && result_->ok) {
    lblResult_->setText(QStringLiteral("残差 %1 · 自由度 %2 · 求解 %3 s")
                            .arg(fmtNum(result_->residual, 2))
                            .arg(result_->ndof)
                            .arg(fmtNum(result_->seconds, 3)));
  } else {
    lblResult_->setText(QStringLiteral("—"));
  }
}

void MainWindow::clearResults() {
  result_.reset();
  post_.reset();
  view_->clearPost();
  prop_->setResult(nullptr);
  prop_->setPost(nullptr);
  results_->setResult(nullptr);
  results_->setPost(nullptr);
  tree_->setPost(nullptr);
  lblResult_->setText(QStringLiteral("—"));
  actContour_->setChecked(false);
  actDeform_->setChecked(false);
  actAnimate_->setChecked(false);
}

void MainWindow::refreshFieldCombo() {
  const QString keep = cmbField_->currentText();
  cmbField_->blockSignals(true);
  cmbField_->clear();
  if (post_) {
    for (size_t i = 0; i < post_->nodeFields().size(); ++i) {
      const auto& f = post_->nodeFields()[i];
      QString label = QString::fromStdString(f.name);
      if (!f.unit.empty()) label += QStringLiteral(" (%1)").arg(QString::fromStdString(f.unit));
      cmbField_->addItem(label, static_cast<int>(i));
    }
    // 默认选"位移幅值"：算完第一眼看的是变形
    const int idx = cmbField_->findText(keep, Qt::MatchContains);
    cmbField_->setCurrentIndex(idx >= 0 ? idx : 0);
  } else {
    cmbField_->addItem(QStringLiteral("（无结果）"), -1);
  }
  cmbField_->blockSignals(false);
}

void MainWindow::log(const QString& s, bool error) {
  (void)error;      // 保留参数：将来要按级别上色（QPlainTextEdit 支持富文本追加）
  if (!log_) return;
  // 时间戳用 HH:mm:ss：长会话里"哪条消息是刚才那条"只能靠时间分辨
  const QString stamp = QTime::currentTime().toString(QStringLiteral("HH:mm:ss"));
  log_->appendPlainText(QStringLiteral("[%1] %2").arg(stamp, s));
}

// =============================================================================
//  打开 / 重新加载
// =============================================================================
bool MainWindow::openScript(const QString& path, bool runImmediately) {
  io::ModelScript ms;
  if (!ms.parseFile(path.toStdString())) {
    QString msg = QStringLiteral("脚本解析失败：\n");
    for (const auto& e : ms.errors()) msg += QString::fromStdString(e) + "\n";
    QMessageBox::warning(this, QStringLiteral("打开失败"), msg);
    log(QStringLiteral("解析失败 %1").arg(path), true);
    for (const auto& e : ms.errors()) log(QString("    ") + QString::fromStdString(e), true);
    return false;
  }

  auto m = std::make_unique<Model>();
  std::string info;
  if (!ms.build(*m, &info)) {
    QString msg = QStringLiteral("建模失败：\n");
    for (const auto& e : ms.errors()) msg += QString::fromStdString(e) + "\n";
    QMessageBox::warning(this, QStringLiteral("建模失败"), msg);
    log(QStringLiteral("建模失败 %1").arg(path), true);
    for (const auto& e : ms.errors()) log(QString("    ") + QString::fromStdString(e), true);
    return false;
  }

clearResults();
  model_ = std::move(m);
  cmds_->clear();                    // 新模型：历史命令栈失效，全部清空
  updateCmdState();
  scriptPath_ = path;
  scriptErrors_ = ms.errors();
  scriptWarnings_ = ms.warnings();
  buildInfo_ = info;
  scriptLoaded_ = true;

  publishModel();
  view_->frameAll();

  log(QStringLiteral("已加载 %1").arg(QFileInfo(path).fileName()));
  log(QString("    ") + QString::fromStdString(info));
  for (const auto& w : scriptWarnings_) log(QStringLiteral("    [警告] ") + QString::fromStdString(w), true);

  setWindowTitle(QStringLiteral("yjk_like  ——  %1").arg(QFileInfo(path).fileName()));

  if (runImmediately) onRun();
  return true;
}

void MainWindow::onOpen() {
  const QString path = QFileDialog::getOpenFileName(
      this, QStringLiteral("打开建模脚本"), QString(),
      QStringLiteral("yjk 建模脚本 (*.yjk);;文本文件 (*.txt);;所有文件 (*)"));
  if (path.isEmpty()) return;
  openScript(path, false);
}

void MainWindow::onReload() {
  if (scriptPath_.isEmpty()) {
    QMessageBox::information(this, QStringLiteral("重新加载"),
                             QStringLiteral("当前模型不是从脚本加载的。"));
    return;
  }
  openScript(scriptPath_, false);
}

// =============================================================================
//  计算
// =============================================================================
// =============================================================================
//  T4：交互建模（工具栏/菜单/右键/快捷键 B C W P L）
//
//  状态机：工具激活 → 视口单击节点（buildNodeHit）→ 攒够点即commit。
//    梁 / 柱   2 点
//    墙        4 点（底左→底右→顶右→顶左，与脚本 wall <a> <b> <c> <d> 同构）
//    板        Slab 工具下拖拽框选，视口发 buildBoxNodes
//    节点荷载  1 点，直接弹对话框输入力/矩
//  提交走 yjk::interact::BuildOps（无 Qt 纯函数，与脚本层共用内核接口），
//  成功后 refreshAfterBuild()：结果失效 + 刷新视图/树 + 若有结果则内存重算。
// =============================================================================
void MainWindow::onToolTriggered(View3D::BuildTool t) {
  // 再次点同一个工具 = 退出（CAD 惯例：Esc / 再点一次都能退出）
  if (buildTool_ == t) {
    view_->setBuildTool(View3D::BuildTool::None);
    return;
  }
  buildTool_ = t;
  pendingNodes_.clear();
  view_->setBuildTool(t);
  QString tip;
  switch (t) {
    case View3D::BuildTool::Beam:
      tip = QStringLiteral("梁工具：依次点击起点、终点节点；Esc 退出");
      break;
    case View3D::BuildTool::Column:
      tip = QStringLiteral("柱工具：依次点击底、顶节点；Esc 退出");
      break;
    case View3D::BuildTool::Wall:
      tip = QStringLiteral("墙工具：依次点击底左→底右→顶右→顶左；Esc 退出");
      break;
    case View3D::BuildTool::Slab:
      tip = QStringLiteral("板工具：拖拽框选板四周节点；Esc 退出");
      break;
    case View3D::BuildTool::Load:
      tip = QStringLiteral("荷载工具：点击节点施加荷载；Esc 退出");
      break;
    default:
      break;
  }
  if (!tip.isEmpty()) statusBar()->showMessage(tip, 4000);
}

void MainWindow::onBuildToolChanged(View3D::BuildTool t) {
  buildTool_ = t;
  pendingNodes_.clear();
  // 同步工具栏/菜单勾选态（Esc 退出后归位）
  for (QAction* a :
       {actToolBeam_, actToolColumn_, actToolWall_, actToolSlab_, actToolLoad_}) {
    const bool on = (a == actToolBeam_ && t == View3D::BuildTool::Beam) ||
                    (a == actToolColumn_ && t == View3D::BuildTool::Column) ||
                    (a == actToolWall_ && t == View3D::BuildTool::Wall) ||
                    (a == actToolSlab_ && t == View3D::BuildTool::Slab) ||
                    (a == actToolLoad_ && t == View3D::BuildTool::Load);
    a->setChecked(on);
  }
}

void MainWindow::gatherBuildNode(int node) {
  if (buildTool_ == View3D::BuildTool::None) return;
  if (node < 0) {   // 点空白：不丢弃已选点？不 —— 建工具里点空就是点空，保持现状
    return;
  }
  // 去重（同一节点连点两次没意义，视为一次点击）
  for (int x : pendingNodes_)
    if (x == node) return;
  pendingNodes_.push_back(node);

  const size_t need =
      (buildTool_ == View3D::BuildTool::Wall) ? 4
      : (buildTool_ == View3D::BuildTool::Beam ||
         buildTool_ == View3D::BuildTool::Column)
          ? 2
          : 1;   // Load / Slab 立即处理
  log(QStringLiteral("已选节点 #%1（%2/%3）")
          .arg(node)
          .arg(pendingNodes_.size())
          .arg(need));

  if (pendingNodes_.size() >= need) {
    commitBuild(std::vector<Id>(pendingNodes_.begin(), pendingNodes_.end()),
                buildTool_);
    pendingNodes_.clear();
  }
}

void MainWindow::commitBuild(const std::vector<Id>& ids, View3D::BuildTool t) {
  if (!model_) return;
  // T5：操作前快照 —— 命令栈记录"这一次建模"的撤销节点
  const Model before = model_->clone();
  // 默认属性：C30 混凝土 + 矩形截面（用户之后可在属性面板细调 ——
  // T3 已打通"改参即写回模型 + 自动重算"链路）。
  const Material mat = Material::concreteC(30);
  const SectionProperties sec = sec::rect(0.3, 0.6);
  const ShellProperties shell;   // 默认厚度等内核取默认值

  yjk::interact::BuildOpResult r;
  switch (t) {
    case View3D::BuildTool::Beam:
      r = yjk::interact::placeBeam(*model_, ids[0], ids[1], sec, mat);
      break;
    case View3D::BuildTool::Column:
      r = yjk::interact::placeColumn(*model_, ids[0], ids[1], sec, mat);
      break;
    case View3D::BuildTool::Wall:
      r = yjk::interact::placeWall(*model_, {ids[0], ids[1], ids[2], ids[3]}, shell);
      break;
    case View3D::BuildTool::Slab:
      r = yjk::interact::placeSlabFromNodes(*model_, ids, shell);
      break;
    case View3D::BuildTool::Load: {
      // 节点荷载：弹框输入力/矩（6 分量），取消则不施加
      const Id n = ids[0];
      bool ok = false;
      const double fx = QInputDialog::getDouble(
          this, QStringLiteral("节点荷载"), QStringLiteral("集中力 Fx (kN) [0]"),
          0.0, -1e6, 1e6, 2, &ok);
      if (!ok) return;
      const double fy = QInputDialog::getDouble(
          this, QStringLiteral("节点荷载"), QStringLiteral("集中力 Fy (kN) [0]"),
          0.0, -1e6, 1e6, 2, &ok);
      if (!ok) return;
      const double fz = QInputDialog::getDouble(
          this, QStringLiteral("节点荷载"), QStringLiteral("集中力 Fz (kN) [0，向下为负]"),
          0.0, -1e6, 1e6, 2, &ok);
      if (!ok) return;
      const double mx = QInputDialog::getDouble(
          this, QStringLiteral("节点荷载"), QStringLiteral("力矩 Mx (kN·m) [0]"),
          0.0, -1e6, 1e6, 2, &ok);
      if (!ok) return;
      const double my = QInputDialog::getDouble(
          this, QStringLiteral("节点荷载"), QStringLiteral("力矩 My (kN·m) [0]"),
          0.0, -1e6, 1e6, 2, &ok);
      if (!ok) return;
      const double mz = QInputDialog::getDouble(
          this, QStringLiteral("节点荷载"), QStringLiteral("力矩 Mz (kN·m) [0]"),
          0.0, -1e6, 1e6, 2, &ok);
      if (!ok) return;
      r = yjk::interact::placeNodeLoad(*model_, n, Vec3{fx, fy, fz},
                                       Vec3{mx, my, mz});
      break;
    }
    default:
      return;
  }

  if (r.ok) {
    log(QString::fromStdString(r.message));
    cmds_->push(before, *model_);                  // T5：可撤销的模型变更
    refreshAfterBuild();
    updateCmdState();
  } else {
    log(QStringLiteral("创建失败：%1").arg(QString::fromStdString(r.message)), true);
    statusBar()->showMessage(QStringLiteral("创建失败：") +
                                 QString::fromStdString(r.message),
                             5000);
  }
}

void MainWindow::refreshAfterBuild() {
  if (!model_) return;
  // 复用 T3 编辑链路：结果失效 + 刷新视图/树 + 自动重算（若有结果）
  onModelEdited();
}

void MainWindow::onBuildNodeHit(int node) { gatherBuildNode(node); }

void MainWindow::onBuildBoxNodes(const std::vector<int>& nodes) {
  if (buildTool_ != View3D::BuildTool::Slab) return;
  log(QStringLiteral("框选 %1 个节点").arg(nodes.size()));
  commitBuild(std::vector<Id>(nodes.begin(), nodes.end()), buildTool_);
}

void MainWindow::onBuildEscaped() {
  pendingNodes_.clear();
  log(QStringLiteral("已退出创建工具"));
}

// =============================================================================
//  T5 删除 + 撤销/重做
//
//  【快照式命令栈】所有模型变更（T4 创建、这里的删除）都以
//  push(操作前, 操作后) 记录；撤销/重做只换 model_ 指针，不重放任何
//  差分计算 —— 与"结果失效铁律③"天然一致：换完模型立即 refreshAfterBuild()。
// =============================================================================
void MainWindow::updateCmdState() {
  if (!actUndo_ || !actRedo_) return;
  actUndo_->setEnabled(cmds_ && cmds_->canUndo());
  actRedo_->setEnabled(cmds_ && cmds_->canRedo());
}

void MainWindow::onDeleteSelected() {
  if (!model_ || !view_) return;
  const auto& sel = view_->selection();
  const int selNode = view_->selectedNode();
  if (sel.empty() && selNode < 0) return;

  const Model before = model_->clone();
  int removed = 0;

  if (!sel.empty()) {
    // removeElement 是"按序号删除、后续前移"：从大到小删，前移不误伤。
    std::vector<int> ids(sel.begin(), sel.end());
    std::sort(ids.rbegin(), ids.rend());
    for (int id : ids) {
      if (model_->removeElement(id)) {
        ++removed;
      } else {
        log(QStringLiteral("删除失败：单元 #%1（越界或不可删除）").arg(id), true);
      }
    }
  }
  if (selNode >= 0) {
    if (model_->removeNode(selNode)) {
      ++removed;
    } else {
      log(QStringLiteral("删除失败：节点 #%1 仍被构件引用（先删除引用它的构件）")
              .arg(selNode),
          true);
    }
  }

  if (removed > 0) {
    cmds_->push(before, *model_);          // 记录一次可撤销的模型变更
    view_->clearSelection();
    log(QStringLiteral("已删除 %1 个对象（Ctrl+Z 可撤销）").arg(removed));
    refreshAfterBuild();
    updateCmdState();
  } else {
    statusBar()->showMessage(QStringLiteral("没有可删除的对象"), 3000);
  }
}

void MainWindow::onUndo() {
  if (!model_ || !cmds_ || !cmds_->canUndo()) return;
  cmds_->undo(model_);
  view_->clearSelection();
  log(QStringLiteral("已撤销（剩余可撤销 %1 步）").arg(cmds_->undoDepth()));
  refreshAfterBuild();                     // 结果失效 + 刷新视图/树 + 重算
  updateCmdState();
}

void MainWindow::onRedo() {
  if (!model_ || !cmds_ || !cmds_->canRedo()) return;
  cmds_->redo(model_);
  view_->clearSelection();
  log(QStringLiteral("已重做（剩余可重做 %1 步）").arg(cmds_->redoDepth()));
  refreshAfterBuild();
  updateCmdState();
}

void MainWindow::onCheck() {
  if (!model_) {
    QMessageBox::information(this, QStringLiteral("模型校核"),
                             QStringLiteral("还没有模型。先打开一个 .yjk 脚本。"));
    return;
  }
  model_->assignDofs();
  const Model::Diagnostic d = model_->check();
  QString msg = QString::fromStdString(d.summary());
  if (d.ok() && d.warnings.empty()) {
    QMessageBox::information(this, QStringLiteral("模型校核"),
                             QStringLiteral("校核通过，未发现问题。"));
  } else {
    QMessageBox::information(this, QStringLiteral("模型校核"), msg);
  }
  log(QStringLiteral("模型校核：") + QString::fromStdString(d.summary()));
  for (const auto& e : d.errors) log(QStringLiteral("    [错误] ") + QString::fromStdString(e), true);
  for (const auto& w : d.warnings) log(QStringLiteral("    [警告] ") + QString::fromStdString(w));
}

void MainWindow::onRun() {
  if (busy_) return;

  // 没有模型时用内置示例，让用户"打开就能看到东西"。
  // 【这不是偷懒】第一印象决定用户会不会继续用：空白窗口 + 要求先找文件，
  // 大多数人到这一步就走了。
  if (!model_ || scriptPath_.isEmpty()) {
    const QString sample = QFileDialog::getOpenFileName(
        this, QStringLiteral("选择建模脚本（取消则用内置示例 frame3x4.yjk）"),
        QApplication::applicationDirPath() + QStringLiteral("/../examples"),
        QStringLiteral("yjk 建模脚本 (*.yjk)"));
    if (!sample.isEmpty()) {
      if (!openScript(sample, true)) return;
      return;
    }
    // 内置示例
    const QString builtin = QCoreApplication::applicationDirPath() +
                            QStringLiteral("/examples/frame3x4.yjk");
    if (QFileInfo::exists(builtin)) {
      if (!openScript(builtin, true)) return;
      return;
    }
    QMessageBox::information(this, QStringLiteral("开始计算"),
                             QStringLiteral("还没有模型。\n\n"
                                            "请从「文件 → 打开建模脚本」载入一个 .yjk 文件，"
                                            "例如项目里的 examples/frame3x4.yjk。"));
    return;
  }

  const QString path = scriptPath_;
  auto b = std::make_shared<Bundle>();
  setBusy(true, QStringLiteral("计算中…"));
  log(QStringLiteral("开始计算 %1").arg(QFileInfo(path).fileName()));

  // ---- 后台线程 ----
  QThread* th = QThread::create([b, path]() {
    QElapsedTimer t;
    t.start();

    io::ModelScript ms;
    if (!ms.parseFile(path.toStdString())) {
      b->errors = ms.errors();
      return;
    }
    b->parsed = true;

    b->model = std::make_unique<Model>();
    if (!ms.build(*b->model, &b->info)) {
      b->errors = ms.errors();
      return;
    }
    b->warnings = ms.warnings();
    b->buildOk = true;
    b->buildMs = t.elapsed();
    t.restart();

    b->model->assignDofs();
    const Model::Diagnostic d = b->model->check();
    b->checkOk = d.ok();
    for (const auto& e : d.errors) b->diagErrors += e + "\n";
    for (const auto& w : d.warnings) b->diagWarnings += w + "\n";

    // 校核不过就不求解 —— 与其让求解器报"刚度矩阵奇异"，
    // 不如直接把"哪个节点悬空/哪个单元退化"告诉用户。
    if (!b->checkOk) return;

    StaticAnalysis sa(*b->model);
    auto sr = std::make_unique<StaticResult>(sa.solve());
    b->solveMessage = sr->message;
    b->solveOk = sr->ok;
    b->solveMs = t.elapsed();
    b->result = std::move(sr);
    if (!b->solveOk) return;

    b->post = std::make_unique<post::PostProcessor>(*b->model, *b->result);
    b->post->compute(21, 0.7);
  });

  connect(th, &QThread::finished, this, [this, b, th]() {
    th->deleteLater();
    setBusy(false);

    if (!b->parsed || !b->buildOk) {
      log(QStringLiteral("建模失败"), true);
      QString msg = QStringLiteral("建模失败：\n");
      for (const auto& e : b->errors) {
        log(QStringLiteral("    ") + QString::fromStdString(e), true);
        msg += QString::fromStdString(e) + "\n";
      }
      QMessageBox::warning(this, QStringLiteral("建模失败"), msg);
      return;
    }

    clearResults();
    model_ = std::move(b->model);
    buildInfo_ = b->info;
    publishModel();

    for (const auto& w : b->warnings) log(QStringLiteral("    [警告] ") + QString::fromStdString(w), true);

    if (!b->checkOk) {
      log(QStringLiteral("模型校核未通过，已中止求解"), true);
      for (const auto& line : QString::fromStdString(b->diagErrors).split('\n'))
        if (!line.isEmpty()) log(QStringLiteral("    [错误] ") + line, true);
      QMessageBox::warning(this, QStringLiteral("模型校核未通过"),
                           QStringLiteral("发现问题，已中止求解（求解器在这种情况下只会报"
                                          "\"刚度矩阵奇异\"，无法定位）：\n\n") +
                               QString::fromStdString(b->diagErrors));
      return;
    }
    for (const auto& line : QString::fromStdString(b->diagWarnings).split('\n'))
      if (!line.isEmpty()) log(QStringLiteral("    [警告] ") + line);

    if (!b->solveOk) {
      log(QStringLiteral("求解失败：") + QString::fromStdString(b->solveMessage), true);
      QMessageBox::warning(this, QStringLiteral("求解失败"),
                           QString::fromStdString(b->solveMessage));
      return;
    }

    result_ = std::move(b->result);
    post_ = std::move(b->post);
    tree_->setPost(post_.get());
    publishResults();

    log(QStringLiteral("求解完成：自由度 %1，残差 %2，耗时 %3 s（建模 %4 ms + 求解 %5 ms）")
            .arg(result_->ndof)
            .arg(fmtNum(result_->residual, 2))
            .arg(fmtNum(result_->seconds, 3))
            .arg(fmtNum(b->buildMs, 3))
            .arg(fmtNum(b->solveMs, 3)));
    log(QStringLiteral("    ") + QString::fromStdString(post_->summary()));
    statusBar()->showMessage(QStringLiteral("计算完成"), 4000);
  });

th->start();
}

// -----------------------------------------------------------------------------
//  T3：属性面板编辑后
// -----------------------------------------------------------------------------
void MainWindow::onModelEdited() {
  if (!model_) return;

  // 铁律③：模型变了，结果必须失效 —— 用户绝不能看着旧云图做新决策
  const bool hadResult = (result_ != nullptr);
  clearResults();

  // 刷新渲染与模型树（血缘关系：编辑的是同一个 model_ 对象，
  // 不能走 setModel —— 它会清空视口的结果/选中态，见 View3D 注释）
  view_->refresh();
  tree_->setModel(model_.get());

  // 恢复选中：刷新树后选中态丢失，把刚才编辑的对象重新高亮
  // 【不能用 view_->selectElement(-1, 0)】那会触发"空选中 →
  // 面板 clear"的联动链，把属性面板的当前对象清掉，下一次编辑就断了。
  if (prop_->currentElem() >= 0) {
    view_->selectElement(prop_->currentElem(), prop_->currentType());
    tree_->selectElement(prop_->currentElem(), prop_->currentType());
  }
  // 节点场景：视口选中态在 refresh() 后依然保留，无需干预；
  // 树没有节点级高亮接口，保持原样即可。

  // 上次重算还没结束（例如编辑器销毁时又提交了一次修改）：
  // 记一笔，等当前重算结束用最新模型补跑，别把这次编辑吞掉。
  if (busy_) {
    pendingRecompute_ = true;
    return;
  }

  // 之前有结果 → 用内存模型重算（onRun 会从磁盘重解析，把编辑丢掉）
  if (hadResult) recomputeFromMemory();
}

void MainWindow::recomputeFromMemory() {
  if (!model_ || busy_) return;

  // UI 线程做深拷贝快照：后台线程 assignDofs/solve 会写 nodes_[].dof[]，
  // 与 UI 线程并发读写同一个 Model 是数据竞争（Model.h 有专门注释）。
  // clone 保持节点/单元顺序，故求解结果的 u 索引与主线程 model_ 对齐。
  auto b = std::make_shared<Bundle>();
  b->model = std::make_unique<Model>(model_->clone());
  setBusy(true, QStringLiteral("属性已修改，重算中…"));
  prop_->setEditEnabled(false);

  QThread* th = QThread::create([b]() {
    b->model->assignDofs();
    const Model::Diagnostic d = b->model->check();
    b->checkOk = d.ok();
    for (const auto& e : d.errors) b->diagErrors += e + "\n";
    for (const auto& w : d.warnings) b->diagWarnings += w + "\n";
    if (!b->checkOk) return;

    StaticAnalysis sa(*b->model);
    auto sr = std::make_unique<StaticResult>(sa.solve());
    b->solveMessage = sr->message;
    b->solveOk = sr->ok;
    b->result = std::move(sr);
  });

  connect(th, &QThread::finished, this, [this, b, th]() {
    th->deleteLater();
    setBusy(false);
    prop_->setEditEnabled(true);

    if (!b->checkOk || !b->solveOk) {
      log(QStringLiteral("重算未通过（结果保持失效，模型仍生效）："), true);
      if (!b->diagErrors.empty()) log(QString::fromStdString(b->diagErrors).trimmed(), true);
      if (!b->solveMessage.empty())
        log(QString::fromStdString(b->solveMessage), true);
      QMessageBox::warning(this, QStringLiteral("重算未通过"),
                           b->checkOk ? QString::fromStdString(b->solveMessage)
                                      : QStringLiteral("模型校核未通过，结果已失效：\n\n") +
                                            QString::fromStdString(b->diagErrors));
    } else {
      // 注意：PostProcessor 持有 Model 的【引用】，
      // 不能引用后台线程里的快照（快照随 b 析构后悬垂）。
      // 主线程 model_ 此刻与快照内容一致（clone + 顺序保持），
      // 用它重建 post 即可安全存活。
      result_ = std::move(b->result);
      post_ = std::make_unique<post::PostProcessor>(*model_, *result_);
      post_->compute(21, 0.7);
      tree_->setPost(post_.get());
      publishResults(false);      // 不自动切视角/云图，别打断编辑工作流

      log(QStringLiteral("重算完成：自由度 %1，残差 %2（属性编辑后自动重算）")
              .arg(result_->ndof)
              .arg(fmtNum(result_->residual, 2)));
      statusBar()->showMessage(QStringLiteral("已重算"), 3000);
    }

    // 重算进行中又被编辑过 → 用最新模型补跑一次（见 onModelEdited）
    if (pendingRecompute_) {
      pendingRecompute_ = false;
      recomputeFromMemory();
    }
  });

  th->start();
}

void MainWindow::setBusy(bool busy, const QString& what) {
  busy_ = busy;
  progress_->setVisible(busy);
  actRun_->setEnabled(!busy);
  actCheck_->setEnabled(!busy);
  actOpen_->setEnabled(!busy);
  actReload_->setEnabled(!busy);
  if (busy) {
    statusBar()->showMessage(what);
    QApplication::setOverrideCursor(Qt::BusyCursor);
  } else {
    QApplication::restoreOverrideCursor();
  }
}

// =============================================================================
//  导出
// =============================================================================
void MainWindow::onExportJson() {
  if (!post_) {
    QMessageBox::information(this, QStringLiteral("导出"), QStringLiteral("还没有计算结果。"));
    return;
  }
  const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("导出 JSON"),
                                                    QStringLiteral("result.json"),
                                                    QStringLiteral("JSON (*.json)"));
  if (path.isEmpty()) return;
  post::ExportOptions opt;
  opt.stationsPerMember = 21;
  const QString stem = path.left(path.lastIndexOf('.'));
  const auto r = post::exportAll(stem.toStdString(), *model_, *result_, *post_, opt);
  if (r.ok) {
    log(QStringLiteral("已导出：%1").arg(QString::fromStdString(r.jsonPath)));
    statusBar()->showMessage(QStringLiteral("导出完成"), 3000);
  } else {
    QMessageBox::warning(this, QStringLiteral("导出失败"), QString::fromStdString(r.error));
  }
}

void MainWindow::onExportHtml() {
  if (!post_) {
    QMessageBox::information(this, QStringLiteral("导出"), QStringLiteral("还没有计算结果。"));
    return;
  }
  const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("导出可视化 HTML"),
                                                    QStringLiteral("viewer.html"),
                                                    QStringLiteral("HTML (*.html)"));
  if (path.isEmpty()) return;
  post::ExportOptions opt;
  const QString stem = path.left(path.lastIndexOf('.'));
  const auto r = post::exportAll(stem.toStdString(), *model_, *result_, *post_, opt);
  if (r.ok) {
    log(QStringLiteral("已导出：%1").arg(QString::fromStdString(r.htmlPath)));
    statusBar()->showMessage(QStringLiteral("导出完成"), 3000);
  } else {
    QMessageBox::warning(this, QStringLiteral("导出失败"), QString::fromStdString(r.error));
  }
}

void MainWindow::onExportReport() {
  if (!post_) {
    QMessageBox::information(this, QStringLiteral("导出"), QStringLiteral("还没有计算结果。"));
    return;
  }
  const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("导出文本报告"),
                                                    QStringLiteral("report.txt"),
                                                    QStringLiteral("文本 (*.txt)"));
  if (path.isEmpty()) return;
  QFile f(path);
  if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QMessageBox::warning(this, QStringLiteral("导出失败"), f.errorString());
    return;
  }
  QTextStream ts(&f);
  ts.setEncoding(QStringConverter::Utf8);
  ts << QString::fromStdString(post_->report());
  log(QStringLiteral("已导出报告：%1").arg(path));
}

void MainWindow::onExportShot() {
  const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("导出视图截图"),
                                                    QStringLiteral("view.png"),
                                                    QStringLiteral("PNG (*.png)"));
  if (path.isEmpty()) return;
  const QImage img = view_->grabImage();
  if (img.isNull() || !img.save(path)) {
    QMessageBox::warning(this, QStringLiteral("导出失败"), QStringLiteral("无法写入图片"));
    return;
  }
  log(QStringLiteral("已导出截图：%1").arg(path));
}

QImage MainWindow::renderToImage() { return view_->grabImage(); }

// -----------------------------------------------------------------------------
//  自检截图
//
//  把"界面能不能正确渲染"变成一条可自动化的断言。
//  手工点界面既慢又不可复现，而渲染出图是可以进 CI 的。
// -----------------------------------------------------------------------------
bool MainWindow::applyShotSettings(const QString& preset, const QString& contour, int selectElem) {
  if (!model_) return false;

  // 选中一个构件：用于验证"视口 ↔ 属性面板 ↔ 模型树"三者的联动
  if (selectElem >= 0 && selectElem < model_->elementCount()) {
    const int t = (model_->elements()[static_cast<size_t>(selectElem)]->type() == ElementType::Beam3D)
                      ? 1 : 2;
    view_->selectElement(selectElem, t);
    prop_->showElement(selectElem, t);
    tree_->selectElement(selectElem, t);
  }
  if (!model_) {
    log(QStringLiteral("截图失败：没有模型"), true);
    return false;
  }

  // 视角
  if (preset == QStringLiteral("front"))      view_->setPreset(Camera::Preset::Front);
  else if (preset == QStringLiteral("left"))  view_->setPreset(Camera::Preset::Left);
  else if (preset == QStringLiteral("top"))   view_->setPreset(Camera::Preset::Top);
  else if (preset == QStringLiteral("right")) view_->setPreset(Camera::Preset::Right);
  else                                        view_->setPreset(Camera::Preset::Iso);

  // 云图变量
  if (post_) {
    actContour_->blockSignals(true);
    cmbField_->blockSignals(true);
    cmbMemberField_->blockSignals(true);

    actContour_->setChecked(true);
    cmbField_->setCurrentIndex(-1);
    cmbMemberField_->setCurrentIndex(0);

    MemberField mf = MemberField::None;
    int nodeField = -1;
    QString want = contour;
    if (want == QStringLiteral("moment"))      mf = MemberField::Moment;
    else if (want == QStringLiteral("shear"))  mf = MemberField::Shear;
    else if (want == QStringLiteral("axial"))  mf = MemberField::Axial;
    else if (want == QStringLiteral("stress")) mf = MemberField::Stress;
    else if (want == QStringLiteral("util"))   mf = MemberField::Utilization;
    else {
      // 节点场：按名字里的关键字找（ux / uy / uz / sx / sy / vm / disp）
      const QStringList keys = {QStringLiteral("位移幅值"), QStringLiteral("ux"),
                                QStringLiteral("uy"), QStringLiteral("uz"),
                                QStringLiteral("sx"), QStringLiteral("sy"),
                                QStringLiteral("von Mises"), QStringLiteral("竖向反力")};
      QString key = keys.value(0);
      if (want == QStringLiteral("ux")) key = QStringLiteral("ux");
      else if (want == QStringLiteral("uy")) key = QStringLiteral("uy");
      else if (want == QStringLiteral("uz")) key = QStringLiteral("uz");
      else if (want == QStringLiteral("sx")) key = QStringLiteral("sx");
      else if (want == QStringLiteral("sy")) key = QStringLiteral("sy");
      else if (want == QStringLiteral("vm")) key = QStringLiteral("von Mises");
      for (int i = 0; i < cmbField_->count(); ++i)
        if (cmbField_->itemText(i).contains(key)) { nodeField = i; break; }
      if (nodeField >= 0) cmbField_->setCurrentIndex(nodeField);
    }
    if (mf != MemberField::None) {
      for (int i = 0; i < cmbMemberField_->count(); ++i)
        if (cmbMemberField_->itemData(i).toInt() == static_cast<int>(mf)) {
          cmbMemberField_->setCurrentIndex(i);
          cmbField_->setCurrentIndex(-1);
          break;
        }
    }

    actContour_->blockSignals(false);
    cmbField_->blockSignals(false);
    cmbMemberField_->blockSignals(false);
    onViewOptionsChanged();
  }

  view_->frameAll();
  return true;
}

bool MainWindow::shootTo(const QString& pngPath, const QString& preset, const QString& contour,
                         int selectElem, bool fullWindow) {
  if (!applyShotSettings(preset, contour, selectElem)) {
    log(QStringLiteral("截图失败：没有模型"), true);
    return false;
  }

  // fullWindow：连面板一起抓（自检用）
  // 否则只抓三维视口（出图用：报告里要的是干净的构件图，不带工具栏）
  const QImage img = fullWindow ? grab().toImage() : view_->grabImage();
  if (img.isNull() || !img.save(pngPath)) {
    log(QStringLiteral("截图写入失败：%1").arg(pngPath), true);
    return false;
  }
  log(QStringLiteral("截图完成：%1 (%2×%3)")
          .arg(pngPath).arg(img.width()).arg(img.height()));
  return true;
}

void MainWindow::onAbout() {
  QMessageBox::about(this, QStringLiteral("关于 yjk_like"),
      QStringLiteral(
          "<b>yjk_like</b> —— 结构有限元分析与设计原型<br><br>"
          "内核：C++17，零外部依赖<br>"
          "界面：Qt %1<br><br>"
          "<b>分层结构</b><br>"
          "① 数学层：稀疏 LDL^T 直接法（消元树符号分析）<br>"
          "② 单元层：3D 梁（Timoshenko）、4 节点壳（膜）<br>"
          "③ 模型层：截面库 · 节点/单元/荷载/约束 · 自由度编号<br>"
          "④ 分析层：线性静力<br>"
          "⑤ 后处理：内力/应力恢复 · 楼层指标 · 云图<br>"
          "⑥ 界面层：本程序（软件三维渲染，无需 GPU）<br><br>"
          "<b>已验证</b>：全部判据对照解析解，见项目 NOTES.md")
          .arg(QString::fromLatin1(qVersion())));
}

}  // namespace ui
