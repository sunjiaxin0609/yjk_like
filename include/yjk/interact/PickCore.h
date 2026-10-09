// =============================================================================
//  yjk/interact/PickCore.h  ——  三维射线拾取核心（纯函数，无 Qt 依赖）
//
//  【为什么要有这么一层】
//  交互层（View3D）负责把鼠标位置换算成世界空间射线，渲染层（SceneData）
//  负责把内核模型翻译成可拾取图元（节点点/梁中心线/壳四边形），而"射线
//  到底命中了什么"是一个纯粹的几何问题 —— 输入一组图元 + 一条射线，
//  输出最近的命中。把它做成无 Qt、无渲染状态的纯函数后：
//    · 可以在 build-p1（w64devkit，无 Qt）上直接单测，不依赖 GUI 环境；
//    · GUI 与脚本层共用同一套命中逻辑，不会出现"界面点得中、脚本算不到"；
//    · 后续命令对象（加/删构件）与拾取解耦，各自可测。
//
//  命中语义（图形学标准约定）：
//    · 节点：以 p 为心、nodeR 为半径的球。命中判据 = 射线到球心距离 ≤ nodeR；
//    · 梁：以中心线线段为轴、beamR 为半径的圆柱（含两端）。命中判据 =
//      射线到线段的最短距离 ≤ beamR；
//    · 壳：四边形（按对角线拆成两个三角形），射线与三角形求交。
//  重叠时按射线参数 t 取最近者（t = 命中点到 o 的距离，d 为单位向量）。
//
//  单位约定跟随内核：长度 m。
// =============================================================================
#pragma once

#include <vector>

#include "yjk/math/Types.h"

namespace yjk {
namespace interact {

// 世界空间射线。d 应为单位向量（pickNearest 内部会再次归一化兜底）。
struct Ray3 {
  Vec3 o;
  Vec3 d;
  constexpr Ray3() = default;
  constexpr Ray3(const Vec3& origin, const Vec3& dir) : o(origin), d(dir) {}
};

// ---- 可拾取图元 ------------------------------------------------------------
struct PickNode {
  Id id{-1};
  Vec3 p{0, 0, 0};
};

struct PickBeam {
  Id id{-1};
  Vec3 a{0, 0, 0};
  Vec3 b{0, 0, 0};
};

struct PickShell {
  Id id{-1};
  Vec3 p[4];
};

enum class PickKind { None = 0, Node, Beam, Shell };

struct PickHit {
  PickKind kind{PickKind::None};
  Id id{-1};
  double t{0.0};      // 命中点到射线原点的距离（越大越远）
  double dist{0.0};   // 命中判据用的最小距离（节点/梁的容差判定，壳=0）
};

// ---- 几何原子（也可以直接单测）---------------------------------------------
// 点到线段的最短距离
double pointSegDist(const Vec3& p, const Vec3& a, const Vec3& b);

// 射线到线段的最短距离。命中时把线段的最近参数 s∈[0,1] 与射线参数 t 写回。
// 返回 false 表示射线与线段所在的平行直线平行（几何上退化）。
bool raySegDist(const Ray3& ray, const Vec3& a, const Vec3& b,
                double* sOut = nullptr, double* tOut = nullptr);

// 射线与三角形求交（Möller–Trumbore）。命中返回 true 并写回 t（>0）。
bool rayTriHit(const Ray3& ray, const Vec3& v0, const Vec3& v1, const Vec3& v2,
               double* tOut = nullptr);

// ---- 主拾取 -----------------------------------------------------------------
// 在一组图元中找射线最近命中。tMax 为有效射程（防极远处无关图元干扰；
// 传 0/负值 = 不限幅）。
// 优先级：仅按 t 最近；t 相同时 节点 > 梁 > 壳（数值上几乎不可能相等，
// 此序仅作防御）。
PickHit pickNearest(const Ray3& ray,
                    const std::vector<PickNode>& nodes,
                    const std::vector<PickBeam>& beams,
                    const std::vector<PickShell>& shells,
                    double nodeR, double beamR, double tMax = 0.0);

}  // namespace interact
}  // namespace yjk