// =============================================================================
//  apps/ui/main.cpp  ——  Qt 界面入口
//
//  两种运行方式：
//    ① 常规：双击运行，从菜单打开 .yjk 脚本
//    ② 自检：yjk_ui --script a.yjk --shot out.png [--preset iso] [--contour disp]
//       无界面渲染一张图后退出。用途是【让界面本身可被回归测试】——
//       手工点界面既慢又不可复现，而"能不能渲染出图"是可以自动化断言的。
// =============================================================================
#include <QApplication>
#include <QCommandLineParser>
#include <QFileInfo>
#include <QImage>
#include <QSurfaceFormat>
#include <QTimer>
#include <cstdio>

#include "MainWindow.h"
#include "Theme.h"

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  app.setApplicationName(QStringLiteral("yjk_like"));
  app.setApplicationDisplayName(QStringLiteral("yjk_like 结构分析"));
  app.setOrganizationName(QStringLiteral("yjk_like"));

  // 全局字体：结构软件的表格数字很多，默认字体在中文 Windows 下偏大，
  // 一个 9 列的结果表会立刻放不下。统一调到 9pt。
  app.setFont(ui::uiFont(9));

  QCommandLineParser parser;
  parser.setApplicationDescription(QStringLiteral("yjk_like 结构分析 Qt 界面"));
  parser.addHelpOption();
  QCommandLineOption optScript({"s", "script"}, QStringLiteral("启动时打开并计算 .yjk 脚本"), "file");
  QCommandLineOption optShot({"", "shot"}, QStringLiteral("渲染一张 PNG 后退出（自检用）"), "png");
  QCommandLineOption optPreset({"", "preset"},
                               QStringLiteral("截图视角：iso|front|left|top"), "name",
                               QStringLiteral("iso"));
  QCommandLineOption optContour({"", "contour"},
                                QStringLiteral("截图云图变量：disp|ux|uy|uz|sx|sy|vm|moment|shear|axial"),
                                "name", QStringLiteral("disp"));
  QCommandLineOption optSize({"", "size"}, QStringLiteral("截图尺寸 WxH"), "WxH",
                             QStringLiteral("1440x900"));
  QCommandLineOption optFull({"", "window"},
                             QStringLiteral("截【整个窗口】（含面板），默认只截三维视口"));
  QCommandLineOption optSelect({"", "select"}, QStringLiteral("预先选中某构件（自检联动用）"),
                               "id", QStringLiteral("-1"));
  parser.addOption(optScript);
  parser.addOption(optShot);
  parser.addOption(optPreset);
  parser.addOption(optContour);
  parser.addOption(optSize);
  parser.addOption(optFull);
  parser.addOption(optSelect);
  parser.process(app);

  ui::MainWindow w;

  if (parser.isSet(optScript)) {
    if (!w.openScript(parser.value(optScript), true)) {
      std::printf("打开脚本失败\n");
      return 2;
    }
  }

  if (parser.isSet(optShot)) {
    // ---- 自检模式：布局完成后渲染一张图，然后退出 ----
    //
    // 【为什么必须 show() 再截图】View3D 的视口尺寸来自实际布局，
    // 不 show 的话控件尺寸还是默认值，截出来的比例和真实窗口对不上。
    // 延迟 350 ms 是为了等布局与首帧绘制稳定（Dock 的尺寸分配是异步的）。
    const QStringList wh = parser.value(optSize).split('x');
    int W = 1440, H = 900;
    if (wh.size() == 2) { W = wh[0].toInt(); H = wh[1].toInt(); }
    w.resize(W, H);
    w.show();

    const QString shot = parser.value(optShot);
    const QString preset = parser.value(optPreset).toLower();
    const QString contour = parser.value(optContour).toLower();
    const bool full = parser.isSet(optFull);
    const int sel = parser.value(optSelect).toInt();
    QTimer::singleShot(350, &app, [&w, shot, preset, contour, full, sel]() {
      // 两种截图用途不同：
      //   · 视口图 = 出图（报告里插的干净图，只有构件和云图）
      //   · 窗口图 = 自检（要连面板一起看，才能发现布局/表格的问题）
      // 【两者共用同一个 applyShotSettings】，所以 --preset/--contour/--select
      // 在两种模式下都生效。
      const bool ok = w.shootTo(shot, preset, contour, sel, full);
      std::printf(ok ? "截图完成：%s\n" : "截图失败：%s\n",
                  shot.toLocal8Bit().constData());
      QCoreApplication::exit(ok ? 0 : 3);
    });
    return app.exec();
  }

  w.show();
  return app.exec();
}
