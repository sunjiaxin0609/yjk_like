// =============================================================================
//  include/yjk/io/YjkWriter.h  ——  Model → .yjk 脚本文本序列化
//
//  【为什么需要它】
//  · T1 验收：Model 建模型 → 导出 → 再读入，节点/单元/荷载集合完全一致。
//  · T6 另存为 .yjk：GUI 里改完模型后要把全量数据落成脚本，与 ModelScript
//    解析互逆 —— 本序列化器就是那条"出"的链路（ModelScript::parse+build
//    是"进"的链路）。
//
//  【输出格式 = 原始构件命令】
//  序列化产物只使用"原始构件命令"（node / beam / column / shell / wall /
//  spring / fix / nodeload / nodemoment / nodeweight / material raw /
//  section raw），不使用轴网命令 —— 轴网是生成式描述（几乘几跨几层），
//  不能无损表达任意增删改后的模型；原始构件命令是"一份数据"，逐条列出
//  每个节点/单元/荷载，读回后与源模型逐字段一致。
//
//  编号语义：
//  · 节点：按 addNode 顺序分配 id = 序号。序列化时跳过 retired 节点并
//    **重映射编号**（活跃节点按 0..n-1 输出，单元节点引用同步转换），
//    这样即使模型里有被合并/删除过的节点也能正确读回。
//  · 单元：保持 elements() 顺序输出；柱（addColumn）与梁（addBeam）都是
//    BeamElement，序列化统一输出为 beam 命令并带显式 up —— 读回后单元
//    类型/几何/截面/材料/up 与源模型逐字段一致（柱没有独立类型标志）。
//  · 墙：ShellElement::isWall() 为 true 时输出 wall，否则输出 shell。
//
//  【未序列化的数据】（T6 全量序列化时补齐）
//  · 刚性楼板（Diaphragm / DofLink）：约束属于分析装配的一部分，本轮
//    T1 验收范围是"节点/单元/荷载集合"，round-trip 测试模型不建楼板。
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