// =============================================================================
//  include/yjk/io/YjkWriter.h  ——  Model → .yjk 脚本文本序列化
//
//  【为什么需要它】
//  · T1 验收：Model 建模型 → 导出 → 再读入，节点/单元/荷载集合完全一致。
//  · T6 另存为 .yjk：GUI 里改完模型后要把全量数据落成脚本，与 ModelScript
//    解析互逆 —— 本序列化器就是那条"出"的链路（ModelScript::parse+build
//    是"进"的链路）。
//
//  【输出格式 = 原始构件命令 + 单元级荷载 + 刚性楼板】
//  序列化产物使用"原始构件命令"（node / beam / column / shell / wall /
//  spring / fix / nodeload / nodemoment / nodeweight / material raw /
//  section raw）以及单元级荷载/释放/楼板命令（release / beamsw /
//  beamseg / beampoint / shellsw / shellp / diaphragm.bind），不用轴网
//  命令 —— 轴网是生成式描述（几乘几跨几层），不能无损表达任意增删改后的
//  模型；原始构件命令是"一份数据"，逐条列出每个节点/单元/荷载，读回后
//  与源模型逐字段一致（T6 全量序列化：自重开关、线荷载段、等效节点荷载、
//  端部释放、板面压/膜压、刚性楼板 DofLink 全部覆盖）。
//
//  编号语义：
//  · 节点：按 addNode 顺序分配 id = 序号。序列化时跳过 retired 节点并
//    **重映射编号**（活跃节点按 0..n-1 输出，单元节点引用同步转换），
//    这样即使模型里有被合并/删除过的节点也能正确读回。
//  · 单元：保持 elements() 顺序输出；柱（addColumn）与梁（addBeam）都是
//    BeamElement，序列化统一输出为 beam 命令并带显式 up —— 读回后单元
//    类型/几何/截面/材料/up 与源模型逐字段一致（柱没有独立类型标志）。
//  · 墙：ShellElement::isWall() 为 true 时输出 wall，否则输出 shell。
//  · 单元级荷载/释放序号 = elements() 中对应类型的出现计数（只数该类型，
//    与 ModelScript `release beam <序号>` 的 seen 计数口径一致）。
//  · 刚性楼板：主节点作为普通 node 输出（fixed[2..4] 以 fix 命令表达），
//    每个 Diaphragm 输出一条 diaphragm.bind —— 读回后由
//    Model::attachRigidDiaphragm 重建，DofLink 系数与源模型逐位一致。
//
//  【数字精度】所有 double 以 17 位有效数字输出，保证按位读回一致。
// =============================================================================
#pragma once

#include <string>

#include "yjk/model/Model.h"

namespace yjk {
namespace io {

// 把模型序列化为 .yjk 脚本文本（原始构件命令）。
// 返回的文本可直接交给 ModelScript::parse + build 读回。
std::string modelToYjk(const Model& m);

}  // namespace io
}  // namespace yjk