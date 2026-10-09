// =============================================================================
//  apps/ui/SceneData.h  ——  渲染场景构建
//
//  职责：把（模型 + 求解结果 + 后处理结果）翻译成一堆【世界坐标的三角形面
//  与线段】。渲染器只认这些，不认梁单元、壳单元、荷载 —— 这样渲染代码
//  不需要知道有限元的任何事，单元的演进不会波及界面。
//
//  【为什么梁要画成"盒子"而不是一条线】
//  线框看起来快，但会丢掉两个工程师非常依赖的信息：
//    · 构件截面方向的相对关系（强轴/弱轴朝向），梁柱对不齐一眼就能看出
//    · 遮挡关系。线框下前后构件叠在一起，密集框架根本读不出空间关系
//  所以默认走实体（每根梁 6 个四边形面），需要看内部时才切回线框。
//
//  【变形与云图必须共用同一套节点坐标】
//  这是最容易出错的地方：如果变形只在骨架上做、云图的插值却用原始坐标，
//  两者会错位，看上去是"云图没跟着变形"。所以 Scene 里所有几何
//  （构件、节点、支座、网格线）都在构建时就已经把变形加进去了。
// =============================================================================
#pragma once

#include <QColor>
#include <QVector3D>

#include <array>
#include <string>
#include <vector>

#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/model/Model.h"
#include "yjk/post/PostProcessor.h"

namespace ui {

// 内核类型在 yjk 命名空间。界面代码里到处写 yjk:: 只会淹没真正的信息，
// 这里一次性引入（只引入类型名，不引入任何符号）。
using yjk::BeamElement;
using yjk::Element;
using yjk::ElementType;
using yjk::Id;
using yjk::Model;
using yjk::Node;
using yjk::ShellElement;
using yjk::StaticResult;
using yjk::Vec3;

// 命名空间别名，不是 using 声明 —— 后处理层里有 BeamResult / ShellResult
// 这类名字，用 using 会把它们全部拉进 ui，与界面自己的类型撞名。
namespace post = yjk::post;

// -----------------------------------------------------------------------------
//  云图变量
// -----------------------------------------------------------------------------
enum class MemberField {
  None = 0,
  Moment,        // 弯矩 M（沿杆绝对值最大）
  Shear,         // 剪力 V
  Axial,         // 轴力 N
  Stress,        // 正应力 σ（截面边缘极值绝对值）
  Utilization,   // 应力比
};

const char* memberFieldName(MemberField f);

// -----------------------------------------------------------------------------
//  渲染选项
// -----------------------------------------------------------------------------
struct RenderOptions {
  bool solid{true};            // 实体显示；false = 线框
  bool showNodes{false};
  bool showSupports{true};
  bool showGrid{true};
  bool showConstraints{true};  // 约束符号

  // 云图：构件级 与 节点级 二选一
  MemberField memberField{MemberField::None};
  int nodeField{-1};           // >= 0 时按节点场上色（优先于 memberField）
  bool contourOn{false};

  // 变形
  bool deformed{false};
  double deformScale{0.0};     // 0 由 autoScale 自动定
  bool autoScale{true};

  // 高亮
  int hoverElem{-1};
  int hoverType{0};            // 0 无 1 梁 2 壳
  std::vector<int> selected;   // 选中的单元索引
  int selectedType{0};
  int extremeElem{-1};         // 包络极值所在构件（用醒目色）
  int extremeType{0};
};

// -----------------------------------------------------------------------------
//  渲染图元
// -----------------------------------------------------------------------------
struct Face {
  QVector3D p[4];
  int n{4};
  QColor fill;
  QColor edge{0, 0, 0, 0};      // 透明 = 不描边
  int elem{-1};                 // 所属单元索引（-1 = 非构件）
  int elemType{0};              // 1 = 梁，2 = 壳
  double value{0.0};            // 云图值（图例与拾取提示用）
  bool hasValue{false};
};

struct Line {
  QVector3D a, b;
  QColor color;
  double width{1.0};
  int elem{-1};
  int elemType{0};
};

// 节点标记。单独一类而不是塞进 Face，
// 因为渲染器对点的处理方式完全不同（画小方块、不参与光照、不参与深度排序的背面剔除）。
struct Point {
  QVector3D p;
  QColor color;
  double size{4.0};
  int node{-1};
};

struct Scene {
  std::vector<Face> faces;
  std::vector<Line> lines;
  std::vector<Point> points;

  // 叠加高亮层：选中 / 悬停 / 极值构件的轮廓。
  //
  // 【为什么必须单独一层，而不是靠改填充色】
  // 云图把整个 rainbow 色标都用满了，任何"醒目的填充色"都会撞上云图里的
  // 某一档 —— 选中色选黄就撞黄色档、选红就撞红色档。实测把选中色设成
  // 琥珀 #FFC107 时，选中一根黄色构件完全看不出来。
  // 正确的做法是【轮廓】：亮色描边画在场景最上层，与填充色无关。
  // 代价是"被挡住的选中构件也会透出来"，但这恰恰是选择高亮想要的行为
  // —— 用户需要知道"我选中的东西在哪"，哪怕它在别的构件后面。
  std::vector<Line> highlight;

  QVector3D lo{0, 0, 0}, hi{0, 0, 0};   // 包围盒
  bool boundsValid{false};

  // 云图范围（已按"包含 0"修正，见 ColorMap.h）
  double vLo{0.0}, vHi{1.0};
  std::string valueLabel;               // 图例标题
  std::string valueUnit;

  // 变形信息（状态栏显示）
  double appliedScale{1.0};
  double maxDisp{0.0};

  // 拾取用的构件中心线（世界坐标，含变形）。
  // 壳单元直接按投影后的四边形做点包含判断，不需要额外存。
  struct PickBeam {
    int elem{-1};
    QVector3D a, b;
  };
  std::vector<PickBeam> pickBeams;

  bool hasContour() const { return !valueLabel.empty(); }

  void addBounds(const QVector3D& p) {
    if (!boundsValid) { lo = hi = p; boundsValid = true; return; }
    lo.setX(std::min(lo.x(), p.x())); lo.setY(std::min(lo.y(), p.y())); lo.setZ(std::min(lo.z(), p.z()));
    hi.setX(std::max(hi.x(), p.x())); hi.setY(std::max(hi.y(), p.y())); hi.setZ(std::max(hi.z(), p.z()));
  }
};

// -----------------------------------------------------------------------------
//  构建器
// -----------------------------------------------------------------------------
class SceneBuilder {
 public:
  // res / pp 都可以为空指针：模型建好但还没算的时候，
  // 界面仍然要能把构件线框画出来（这是用户确认"模型建对了没"的关键一步）。
  SceneBuilder(const Model& m, const StaticResult& r, const post::PostProcessor* pp)
      : m_(m), res_(r), pp_(pp) {}

  Scene build(const RenderOptions& opt) const;

  // 自动变形放大系数：让最大位移占模型对角尺度的 4%。
  //
  // 【为什么要自动】位移量级（1e-3 m）与模型尺度（1e1 m）差 4 个数量级，
  // 放大系数用 1 时结构性变形完全看不见，用户会以为"没算"；
  // 固定成 100 又会在小模型上糊成一团。按尺度自适应才是对的。
  double autoDeformScale() const;

  // 变形后的节点坐标（Scene 与状态栏共用同一份换算，避免两处不一致）
  QVector3D nodePos(Id n, double scale) const;

 private:
  void addNodes(const RenderOptions& opt, double scale, Scene& s) const;
  void addSupports(const RenderOptions& opt, double scale, Scene& s) const;
  void addGrid(const RenderOptions& opt, Scene& s) const;
  void addConstraintGlyph(const Node& nd, const QVector3D& p, double sz, Scene& s) const;

  // 构件级云图取值
  double memberValue(const post::BeamResult& br, MemberField f) const;

  const Model& m_;
  const StaticResult& res_;
  const post::PostProcessor* pp_;      // 可空
};

}  // namespace ui
