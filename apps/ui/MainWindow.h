// =============================================================================
//  apps/ui/MainWindow.h  ——  主窗口
//
//  【界面层的三条铁律】
//
//  ① 界面层不计算，只搬运。
//     所有力学结果都来自 yjk_core（StaticAnalysis / PostProcessor）。
//     界面里出现的任何数字都必须能在命令行 `yjk run` 的输出里找到同一个值
//     —— 否则一旦对不上，无法判断是内核错了还是界面显示错了。
//
//  ② 求解必须离开 UI 线程。
//     大模型求解是秒级到分钟级的。放在 UI 线程里，窗口会白屏假死，
//     Windows 会给标题栏加"无响应"，用户的第一反应是强杀进程。
//
//  ③ 模型改了就必须让结果失效。
//     "改完模型但忘记重算"是结构软件最危险的场景：用户看着一张
//     旧结果的云图做设计决策。所以任何会改模型的操作都要先 clearResults()。
// =============================================================================
#pragma once

#include <QMainWindow>
#include <QString>

#include <memory>
#include <string>
#include <vector>

#include "SceneData.h"
#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/io/ModelScript.h"
#include "yjk/model/Model.h"
#include "yjk/post/PostProcessor.h"

class QAction;
class QCheckBox;
class QComboBox;
class QLabel;
class QPlainTextEdit;
class QProgressBar;

namespace ui {

// 内核命名空间的别名。界面代码里到处写 yjk:: 只会淹没真正要读的信息。
namespace io = yjk::io;
using yjk::StaticAnalysis;

class View3D;
class ModelTreePanel;
class PropertyPanel;
class ResultPanel;

class MainWindow : public QMainWindow {
  Q_OBJECT
 public:
  explicit MainWindow(QWidget* parent = nullptr);
  ~MainWindow() override;

  // 打开建模脚本（用于命令行参数 / 拖放）
  bool openScript(const QString& path, bool runImmediately = false);

  // 无界面自检用：把当前视图渲染成 PNG
  QImage renderToImage();

  // 自检截图：设定视角与云图变量后渲染并存盘。
  // 走显式方法而不是 QMetaObject::invokeMethod(名字字符串) ——
  // 名字拼错只会在运行时静默失败，而显式方法编译期就能查出来。
  bool shootTo(const QString& pngPath, const QString& preset, const QString& contour,
               int selectElem = -1, bool fullWindow = false);

 private:
  // 应用自检截图前的视图设置（视角 / 云图 / 选中）。
  // 【必须与"截哪一块"分开】否则一旦选"截整窗"，
  // 设置这一步就被跳过，命令行参数全部失效 —— 这是个真实踩过的坑。
  bool applyShotSettings(const QString& preset, const QString& contour, int selectElem);

 public:

 private slots:
  void onOpen();
  void onReload();
  void onRun();
  void onCheck();
  void onExportJson();
  void onExportHtml();
  void onExportReport();
  void onExportShot();
  void onAbout();
  void onViewOptionsChanged();

 private:
  void buildActions();
  void buildMenus();
  void buildToolBar();
  void buildDocks();
  void buildStatusBar();

  // 把当前 model_ 交给各面板 / 视口
  void publishModel();
  // autostyle=true 时自动切变形+云图+frameAll；属性编辑重算时传 false，
  // 不打断编辑工作流（视角/云图状态保持）
  void publishResults(bool autostyle = true);
  void clearResults();
  void setBusy(bool busy, const QString& what = QString());

  void log(const QString& s, bool error = false);
  void refreshFieldCombo();
  void syncOptionsToView();

  // ---- T3 属性编辑 -> 内存重算 ----
  // 属性面板改参后：结果失效 -> 刷新视图/树 -> 若之前有结果则后台重算。
  // 【为什么不能复用 onRun】onRun 每次都从磁盘 .yjk 重新解析，
  // 会把 GUI 里刚做的编辑整个丢掉。
  void onModelEdited();
  void recomputeFromMemory();
  // 后台重算进行中又收到模型编辑时置位：当前重算结束后再补跑一次，
  // 保证「重算结果 = 最后一次编辑后的模型」而不是吞掉中间修改。
  bool pendingRecompute_{false};

  // ---- 数据（界面持有，生命周期覆盖所有面板）----
  std::unique_ptr<Model> model_;
  std::unique_ptr<StaticResult> result_;
  std::unique_ptr<post::PostProcessor> post_;

  QString scriptPath_;
  QString outDir_;
  std::string buildInfo_;
  std::vector<std::string> scriptErrors_, scriptWarnings_;
  bool scriptLoaded_{false};

  // ---- 子控件 ----
  View3D* view_{nullptr};
  ModelTreePanel* tree_{nullptr};
  PropertyPanel* prop_{nullptr};
  ResultPanel* results_{nullptr};
  QPlainTextEdit* log_{nullptr};

  // ---- 工具栏控件 ----
  QAction* actOpen_{nullptr};
  QAction* actReload_{nullptr};
  QAction* actRun_{nullptr};
  QAction* actCheck_{nullptr};
  QAction* actSolid_{nullptr};
  QAction* actDeform_{nullptr};
  QAction* actAnimate_{nullptr};
  QAction* actOrtho_{nullptr};
  QAction* actNodes_{nullptr};
  QAction* actSupports_{nullptr};
  QAction* actGrid_{nullptr};
  QAction* actContour_{nullptr};
  QComboBox* cmbField_{nullptr};
  QComboBox* cmbMemberField_{nullptr};
  QComboBox* cmbPreset_{nullptr};

  QLabel* lblModel_{nullptr};
  QLabel* lblResult_{nullptr};
  QProgressBar* progress_{nullptr};

  bool busy_{false};
};

}  // namespace ui