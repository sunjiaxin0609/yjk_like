// =============================================================================
//  yjk/interact/BuildOps.h  ——  交互建模操作（纯函数，无 Qt 依赖）
//
//  【与 PickCore 的关系】
//  T2 的 PickCore 回答"鼠标点在哪里"（拾取），本模块回答"拾取后干什么"
//  （建模）。同样做成无 Qt、无渲染状态的纯函数后：
//    · 可以在 build-p1（w64devkit，无 Qt）上直接单测，不依赖 GUI 环境；
//    · GUI 与脚本层共用同一套建模逻辑 —— GUI 里手点的构件，
//      导出的 .yjk 与脚本写的构件逐项等价（T4 验收点）；
//    · 后续命令对象（撤销/重做，T5）直接在 BuildOps 之上包装，
//      回退时拿到的是同一套结构化错误消息。
//
//  【所有操作都走 Model 的 T1 可变接口】
//  placeBeam/placeColumn/placeWall/placeSlabFromNodes/placeNodeLoad
//  分别只调用 Model::addBeam/addColumn/addWall/addShell/addNodeLoad ——
//  不做任何几何上的额外秘密。校验（节点存在、非零长度、四点不退化）
//  在这里做而不是散落 GUI，是为了让错误消息结构化、可单测。
//
//  【单位约定跟随内核】长度 m、力 kN、力矩 kN·m。
// =============================================================================
#pragma once

#include <array>
#include <string>
#include <vector>

#include "yjk/math/Types.h"
#include "yjk/model/Model.h"
#include "yjk/model/Section.h"

namespace yjk {
namespace interact {

// 一次建模操作的结果。
// ok=false 时 message 是面向用户的结构化中文错误消息，GUI 直接展示。
struct BuildOpResult {
  bool ok{false};
  Id id{-1};                 // 新建单元（或荷载作用的节点）id
  std::string message;
};

// 两点定梁：nI → nJ。upHint 缺省 = 梁平面 {0,0,-1}（与脚本 beam 一致）。
BuildOpResult placeBeam(Model& m, Id nI, Id nJ, const SectionProperties& sec,
                        const Material& mat, const Vec3& up = Vec3{0, 0, -1});

// 两点定柱：nI → nJ（竖向），走 addColumn（默认 up {0,-1,0}）。
BuildOpResult placeColumn(Model& m, Id nI, Id nJ, const SectionProperties& sec,
                          const Material& mat);

// 沿轴线定墙：底边两点 + 顶边两点，共 4 节点（逆时针：底左、底右、顶右、顶左）。
// 与脚本 wall <a> <b> <c> <d> 完全同构；校验四点不退化（不共线、面积非零）。
BuildOpResult placeWall(Model& m, const std::array<Id, 4>& ns,
                        const ShellProperties& props);

// 框选区域定板：boxed = 框选收集到的节点 id（可能有冗余内部节点）。
// 取这些节点 (x, y) 包围盒的四角最近节点建四边形板，与脚本 shell 等价。
// 框内节点数 < 4 或四角退化时返回结构化错误。
BuildOpResult placeSlabFromNodes(Model& m, const std::vector<Id>& boxed,
                                 const ShellProperties& props);

// 点击定节点荷载：单点力 + 力矩（缺省零），走 addNodeLoad。
BuildOpResult placeNodeLoad(Model& m, Id n, const Vec3& F,
                            const Vec3& M = Vec3{0, 0, 0});

}  // namespace interact
}  // namespace yjk