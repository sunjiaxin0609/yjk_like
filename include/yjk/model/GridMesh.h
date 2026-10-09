// =============================================================================
//  yjk/model/GridMesh.h  ——  轴网 → 有限元模型 的自动生成器
//
//  【为什么必须有这一层】
//  手工建模型时，一块楼板要列 4 个节点、一条梁要列 2 个节点。
//  一栋 6 层 3×4 跨的楼就有几千个节点 —— 手工列不现实，也必须保证
//  "板边节点与梁节点是同一个节点"。这是盈建科/PKPM 用轴网建模的根本原因。
//
//  【最关键的一条工程约束：板边节点必须与梁节点重合】
//
//  若把一跨板剖成 3×3 个壳单元，而沿该跨的梁仍是"柱到柱一根通长梁"：
//     · 板边的中间节点【只被壳单元使用】，不连任何梁；
//     · 这些节点竖向自由、只有板给它刚度 ⇒ 板边在两柱之间实际上是"悬空"的；
//     · 荷载通过板传到边上后【传不到梁上】，计算结果完全错；
//     · 求解器不会报任何错 —— 模型看起来完全正常。
//
//  ⇒ 所以必须【先确定每条轴线的一维细分表】，梁与板共用同一套细分点。
//    本实现的做法是：细分挂在【跨度】上（divX[j] / divY[i]），
//    梁沿该跨生成 divX[j] 段，板的边也落在同样的点上 —— 天然重合。
//
//  【坐标与编号约定】
//    · X 向轴线：沿 X 延伸，位置由 y 确定（编号习惯 A / B / C）
//    · Y 向轴线：沿 Y 延伸，位置由 x 确定（编号习惯 1 / 2 / 3）
//    · 交点 (ix, jy) 位于 (x[jy], y[ix])
//    · 板格 (ix, jy)：X 向跨度 ix（y[ix]→y[ix+1]）× Y 向跨度 jy（x[jy]→x[jy+1]）
//    · 楼层 s 的楼面标高 z[s]；柱从 z[s] 连到 z[s+1]
//    · 壳单元节点按 (p,q)→(p+1,q)→(p+1,q+1)→(p,q+1) 排列，
//      在 +X→+Y 的手系下法向朝 +Z（朝上），符合 ShellElement4 的约定
//    · 墙（剪力墙）：位于【轴线上】的竖直壳，跨楼层 z[s]→z[s+1] 一段。
//      沿 X 的墙在 X 向轴线 ix 上（y 固定），跨 jy；沿 Y 的墙在 Y 向轴线 jy
//      上（x 固定），跨 ix。节点顺序与法向约定见 generate() 内的注释。
//
//  【不做什么】不生成荷载、不生成支座。这两件事依赖具体工程判断，
//  由用户在生成后按 nodeAtAxis() 查到的节点号自行施加。
// =============================================================================
#pragma once

#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <vector>

#include "yjk/model/Model.h"

namespace yjk {

// -----------------------------------------------------------------------------
//  轴网：纯几何 + 拓扑，不含任何有限元概念
// -----------------------------------------------------------------------------
class AxisGrid {
 public:
  // ---- 轴线（按升序添加；generate 时会自动排序并给出警告）----
  Id addAxisX(double y, const std::string& label = "");   // 沿 X 延伸，y 定位
  Id addAxisY(double x, const std::string& label = "");   // 沿 Y 延伸，x 定位
  Id addStory(double z, const std::string& label = "");   // 楼面标高（含基础顶/地面）

  Id nx() const { return static_cast<Id>(y_.size()); }        // X 向轴线数
  Id ny() const { return static_cast<Id>(x_.size()); }        // Y 向轴线数
  Id nStories() const { return static_cast<Id>(z_.size()); }

  double axisX(Id i) const { return y_[static_cast<size_t>(i)]; }
  double axisY(Id j) const { return x_[static_cast<size_t>(j)]; }
  double storyZ(Id s) const { return z_[static_cast<size_t>(s)]; }
  const std::string& axisXLabel(Id i) const { return yLabel_[static_cast<size_t>(i)]; }
  const std::string& axisYLabel(Id j) const { return xLabel_[static_cast<size_t>(j)]; }

  // ---- 构件开关 ----
  //  除墙外所有开关默认"开"（楼板默认不含底层，见 setSlabStory 的说明）。
  //
  //  【为什么墙默认"关"】柱/梁/板是"框架楼盖必然有"的构件；墙则不是
  //  每个轴网都有的 —— 若默认开，既有框架模型会凭空多出满轴线的墙，
  //  刚度、质量、自重全变。所以墙的开关语义与其他构件相反：**不设即无墙**。
  //
  //  索引越界时返回 false 且不改变状态 —— 建模阶段的参数错误不该抛异常。
  bool setColumn(Id ix, Id jy, bool on);          // 交点 (ix,jy) 是否有柱
  bool setBeamAlongX(Id ix, Id jy, bool on);      // 沿 X 的梁：X 向轴线 ix，跨 jy
  bool setBeamAlongY(Id jy, Id ix, bool on);      // 沿 Y 的梁：Y 向轴线 jy，跨 ix
  bool setSlab(Id ix, Id jy, bool on);            // 板格 (ix,jy) 是否有板
  bool setSlabStory(Id s, bool on);               // 第 s 层楼面是否有板
  bool setWallAlongX(Id ix, Id jy, bool on);      // 沿 X 的墙：X 向轴线 ix，跨 jy
  bool setWallAlongY(Id jy, Id ix, bool on);      // 沿 Y 的墙：Y 向轴线 jy，跨 ix
  bool setDivX(Id jy, int n);                     // Y 向跨度 jy 沿 X 的剖分数（≥1）
  bool setDivY(Id ix, int n);                     // X 向跨度 ix 沿 Y 的剖分数（≥1）

  bool hasColumn(Id ix, Id jy) const;
  bool hasBeamAlongX(Id ix, Id jy) const;
  bool hasBeamAlongY(Id jy, Id ix) const;
  bool hasSlab(Id ix, Id jy) const;
  bool hasSlabStory(Id s) const;
  bool hasWallAlongX(Id ix, Id jy) const;
  bool hasWallAlongY(Id jy, Id ix) const;
  int divX(Id jy) const;
  int divY(Id ix) const;

  // ---- 派生：细分后的一维坐标表（梁与板共用，保证节点重合）----
  //
  //  xPts: x[j] 之间按 divX[j] 等分后展开成一维表
  //  xAxisIndex(j) = 轴线 j 在该表中的下标
  //
  //  声明为 const：表是纯派生量，缓存在 mutable 成员里，
  //  这样 GridMesh 持有 const AxisGrid& 也能随时调用。
  void buildPoints() const;
  Id nXPoints() const { return static_cast<Id>(xs_.size()); }
  Id nYPoints() const { return static_cast<Id>(ys_.size()); }
  double xPt(Id p) const { return xs_[static_cast<size_t>(p)]; }
  double yPt(Id q) const { return ys_[static_cast<size_t>(q)]; }
  Id xAxisIndex(Id j) const { return xAxisIdx_[static_cast<size_t>(j)]; }
  Id yAxisIndex(Id i) const { return yAxisIdx_[static_cast<size_t>(i)]; }

  // 排序（升序）。返回是否发生了重排 —— 重排说明输入顺序有误。
  bool sortAxes();
  bool sorted() const;

 private:
  std::vector<double> y_, x_, z_;                  // X向轴线(y) / Y向轴线(x) / 楼面标高
  std::vector<std::string> yLabel_, xLabel_, zLabel_;
  std::vector<char> col_, beamX_, beamY_, slab_;   // 开关，索引 = ix*ny+jy
  std::vector<char> wallX_, wallY_;                // 墙开关，同一索引（默认关）
  std::vector<char> slabStory_;
  std::vector<int> divX_, divY_;
  // 细分后的一维坐标表（buildPoints 的缓存，逻辑上是 const）
  mutable std::vector<double> xs_, ys_;
  mutable std::vector<Id> xAxisIdx_, yAxisIdx_;
};

// -----------------------------------------------------------------------------
//  网格生成器：轴网 + 构件规格 → 模型
// -----------------------------------------------------------------------------
class GridMesh {
 public:
  // 构件规格。柱与梁用梁单元，板与墙用壳单元。
  struct Spec {
    SectionProperties columnSection, beamSection;
    Material columnMaterial, beamMaterial;
    ShellProperties slab;
    ShellProperties wall;           // 剪力墙规格（默认厚 0.2 混凝土壳）

    // 局部 up 方向（决定 Iy/Iz 哪个是平面内 —— 见 BeamElement3D 的说明）
    //
    //  柱：杆轴竖直，up 必须取水平方向，否则投影为零、单元退化。
    //  梁：up 取 (0,0,-1)，即截面的"竖向"朝下，符合结构习惯。
    //      水平梁不论沿 X 还是沿 Y，(0,0,-1) 都合法。
    Vec3 columnUp{0, -1, 0};
    Vec3 beamUp{0, 0, -1};

    // 是否生成楼板【壳单元】。
    //
    //  · true（默认）—— 建板。楼面剖分成壳单元，荷载直接加压在板上，
    //    楼板的面外刚度也参与工作。
    //  · false —— **不建板**。这是真实工程的常规做法：楼板单元数量巨大、
    //    且楼板本身不参与抗侧（平面内刚度无穷大是很好的近似）。
    //    此时"楼面荷载怎么传到梁上"必须靠 applySlabLoad 按 45° 线导荷，
    //    "平面内刚度"靠刚性楼板假定（Model::addRigidDiaphragm）提供。
    //
    //  【这两件事必须成对】只关掉板而不导荷，楼面荷载会凭空消失 ——
    //  而模型仍然可解、残差正常，只有"总荷载对不上"这一条能暴露它。
    bool buildSlabs{true};

    // 是否生成墙【壳单元】。轴网上没有墙开关时（AxisGrid 默认无墙），
    // 此开关不产生任何单元。
    bool buildWalls{true};
  };

  struct Result {
    bool ok{false};
    std::string error;
    Id nodes{0}, columns{0}, beams{0}, slabs{0}, walls{0};
    std::string message() const;
  };

  GridMesh(const AxisGrid& grid, const Spec& spec) : grid_(grid), spec_(spec) {}

  // 把节点与单元写进模型。
  // 【前置】调用前模型可以已有内容 —— 生成器只在需要时新增节点，
  // 且用坐标去重，不会与已有节点重复。
  Result generate(Model& m);

  // ---- 生成后查询（施加支座/荷载必须靠它）----
  Id nodeAt(Id p, Id q, Id story) const;                 // 细分网格坐标 → 节点号
  Id nodeAtAxis(Id ix, Id jy, Id story) const;           // 轴线交点 → 节点号
  std::vector<Id> baseNodes() const;                     // 底层柱脚节点
  Id nNodes() const { return static_cast<Id>(nodes_.size()); }

  // ---- 刚性楼板 ----
  //
  // 逐层把该层节点绑成刚性楼板。fromStory 之前的楼层跳过
  // （嵌固层/基础顶通常整体固定，不需要刚性假定）。
  // 返回新建的主节点数。
  Id addRigidDiaphragms(Model& m, Id fromStory = 1, bool coupleRz = true) const;

  // ---- 板面荷载导荷 ----
  //
  // 【为什么需要它】一旦决定"不建板"，楼面荷载就失去了传递路径 ——
  // 荷载得有路走到梁上。本函数按【45° 线（角平分线）法】把每个板格上的
  // 均布面荷载分到四边的梁上：
  //
  //      短边方向的梁 → 三角形荷载
  //      长边方向的梁 → 梯形荷载
  //      峰值均为 q·L短/2
  //
  // 【峰值为什么取 q·L短/2 而不是 q·L该边/2】
  // 45° 线把板格分成四块，每块的"从属宽度"在边中部是 L短/2
  // （按角平分线几何算出来就是这个值）。取这个值才【守恒】：
  //      2·(½·L短·qL短/2) + 2·[qL短/2·(L长−L短) + ½·L短·qL短/2] = q·L短·L长
  // 若按各边自身长度取峰值，荷载会不守恒（正方形时差一倍）。
  //
  // 【长短边比 ≥ oneWayRatio 时按单向板】荷载全部由长边方向的两根梁承担
  // （板沿短向受力，支承它的正是那两条长边）。GB 50010 规定 ≥3，
  // 行业习惯取 2 —— 取 2 是偏安全的（长边梁分到的更多）。
  struct SlabLoadOptions {
    enum class Path {
      Auto,      // 有板单元 → 压在板上；不建板 → 导到梁上
      ToShells,  // 强制加压在板单元上
      ToBeams,   // 强制按 45° 线导到梁上
    };
    Path path{Path::Auto};
    double q{0.0};              // 面荷载 kPa，正 = 向下（与 slabload 命令一致）
    double oneWayRatio{2.0};    // 长短边比达到该值时按单向板
  };

  struct SlabLoadReport {
    int bays{0};            // 处理的板格数
    int oneWayBays{0};      // 其中按单向板处理的
    int shellElems{0};      // 加压力到壳单元的数量
    int beamSegs{0};        // 生成的梁上分布荷载段数
    double total{0.0};      // 施加的总荷载 kN（可与竖向反力合计对账）
    std::vector<std::string> warnings;
    std::string message() const;
  };
  SlabLoadReport applySlabLoad(Model& m, const SlabLoadOptions& opt) const;

 private:
  const AxisGrid& grid_;
  Spec spec_;
  std::vector<Id> nodes_;      // 扁平索引 (p * nYPoints + q) * nStories + s
  Id nXp_{0}, nYp_{0}, nSt_{0};
};

// -----------------------------------------------------------------------------
//  合并坐标重合的节点
//
//  【为什么需要】用户手工建模时最容易在两个构件"看起来连上了"的地方
//  各建一个节点。Model::check() 会把它报成"坐标重合"，但只报不改。
//
//  【合并策略】保留编号小的那个，把荷载累加过去，
//  并把所有单元里的引用重定向。被合并掉的节点随后不再被任何单元使用、
//  且荷载已清零 ⇒ 由 DofNumbering 的悬空节点逻辑自动锁定，不产生自由度。
//
//  【为什么不真的删除节点】结果向量按【节点序 × 6】索引，
//  删除会让历史节点号全部错位。保留 + 自动锁定是等价且安全的做法。
// -----------------------------------------------------------------------------
struct MergeReport {
  int merged{0};               // 被合并掉的节点数
  int redirected{0};           // 被重定向的单元引用数
};
MergeReport mergeCoincidentNodes(Model& m, double tol = 1e-6);

}  // namespace yjk
