// =============================================================================
//  include/yjk/io/ModelScript.h  ——  .yjk 文本建模脚本
//
//  【为什么要有文本建模脚本，而不是只做 C++ API】
//    · 工程师描述一个工程只需要几十行，不需要写 main() 和编译
//    · 脚本是【数据】不是【代码】→ 可以版本管理、diff、批量生成参数化方案
//    · 命令行 `yjk run a.yjk` 就能从模型直接走到结果文件
//
//  【格式设计取舍】
//  选了"扁平命令 + 点分命名空间"，没有做缩进块。
//  缩进块（YAML 式）可读性更好，但要处理缩进层级、列表项、续行，
//  解析器的复杂度是这里的两倍以上，而出错信息反而更难定位。
//  一行一命令 ⇒ 报错可以精确到行号，这是建模时最需要的。
//
//  语法：
//    # 注释（行内任意位置起）
//    material   <名>  concrete <牌号>          C30→concrete 30
//    material   <名>  steel    <牌号>
//    section    <名>  rect   <b> <h>
//    section    <名>  circle <d>
//    section    <名>  tube   <d> <t>
//    section    <名>  i      <h> <b> <tw> <tf>
//    section    <名>  box    <h> <b> <tw> <tf>
//
//    grid.axisX   <y...>                    X 向轴线的 y 坐标
//    grid.axisY   <x...>                    Y 向轴线的 x 坐标
//    grid.story   <z...>                    楼面标高（第一个为基底）
//    grid.divX    <n...>                    每个 Y 向跨度沿 X 的剖分数
//    grid.divY    <n...>                    每个 X 向跨度沿 Y 的剖分数
//    grid.column  <截面名> <材料名>
//    grid.beam    <截面名> <材料名>
//    grid.slabthick <t>
//    grid.slaboff  <层号...>                这些层不布板（默认底层不布）
//    grid.nocolumn <ix> <jy>                关掉交点柱
//    grid.nobeamx  <ix> <jy>                关掉沿 X 的梁
//    grid.nobeamy  <jy> <ix>                关掉沿 Y 的梁
//    grid.noslab   <ix> <jy>                关掉板格
//
//    support base  fixed|pinned|roller      柱脚
//    support node  <id> fixed|pinned|roller
//    support story <层>  fixed|pinned|roller
//
//    selfweight on|off                      所有梁柱自重
//    slabload   <q>                         板面荷载 kPa（正数 = 向下）
//    beamload   <qz>                        梁上竖向均布 kN/m（正数 = 向下）
//    nodeload   node  <id>    <fx> <fy> <fz>
//    nodeload   story <层>    <fx> <fy> <fz>
//    nodeload   top           <fx> <fy> <fz>
//
//    case <工况名>                   切换荷载工况（缺省工况 "D"）
//    modal [nmodes]                  请求模态分析（缺省/auto/0 = 自动取值
//                                    max(9, 层数×3)，GB 50011 每方向≥9）
//    spectrum <alphaMax> <Tg> [srss|cqc] [verticalK]
//                                    请求反应谱分析；第 4 参 verticalK 为
//                                    竖向地震影响系数放大系数（默认 0.65，
//                                    GB 50011 5.3.1 取水平地震作用的 65%）
//    gmass <恒载工况> <活载工况> [ψ=0.5]
//                                    GB 50011 5.1.3 重力荷载代表值转质量：
//                                    恒载 + ψ·活载 折算为节点重量（需 case 命令）
//    combo <kind> <名> <case...>     荷载组合定义
//      kind ∈ standard|basicvar|basicperm|seismic
//      standard( D L ) = D + L            （正常使用标准组合）
//      basicvar( D L [W] ) = 1.2D+1.4L(+1.4·0.7·W)   （基本组合·可变控制）
//      basicperm( D L [W] ) = 1.35D+1.4·0.7·L(+1.4·0.7·W)（基本组合·永久控制）
//      seismic( D L Ex [Ev] ) = 1.2(D+0.5L)+1.3Ex(+0.5Ev) （地震组合）
//
//    out <输出文件名前缀>
//
//  【工况与动态分析】
//  脚本中若出现 case/modal/spectrum/combo 中任意命令，build() 只建立
//  几何与约束（不含任何荷载）；CLI 随后对每个工况调用 applyCase() 施加
//  该工况的荷载并分别求解，再执行模态/反应谱/组合/包络与导出。
//  不写这些命令时行为与旧版完全一致：所有荷载同属缺省工况 "D"。
// =============================================================================
#pragma once

#include <map>
#include <functional>
#include <string>
#include <vector>

#include "yjk/model/GridMesh.h"
#include "yjk/model/Model.h"

namespace yjk {
namespace io {

class ModelScript {
 public:
  // 解析。失败时 errors() 给出带行号的原因。
  bool parse(const std::string& text);
  bool parseFile(const std::string& path);

  // 按脚本建模型（几何 + 约束 + 支座；见文件头"工况与动态分析"）。
  // info 可选，接收一行人类可读的建模摘要。
  // build 会往 errors_/warnings_ 里追加诊断，所以不是 const
  bool build(Model& m, std::string* info = nullptr);

  // ---- 荷载工况 ----
  // 通常为主线内没有 case 命令 → 返回空（CLI 视为唯一的缺省工况 "D"）。
  // 该工况下的全部荷载已由 build() 直接施加到 m 上（旧行为）。
  bool hasCases() const { return !caseOrder_.empty(); }
  const std::vector<std::string>& caseNames() const { return caseOrder_; }

  // 对已 build 好的模型施加指定工况的荷载（供多工况分析逐 case 调用）。
  // 需要 build() 先生成网格（脚本含轴网时自动）；失败把原因推进 errors_。
  bool applyCase(Model& m, const std::string& name, std::string* info = nullptr);

  // ---- 动态分析请求（modal/spectrum/combo/gmass）----
  bool hasDynamic() const { return hasModalCmd_ || hasSpectrum_ || hasGmass_; }
  // 是否显式请求模态分析（modal 命令出现即 true）
  bool hasModal() const { return hasModalCmd_; }
  int modalNmodes() const { return modalNmodes_; }    // 0 = 自动（ModalAnalysis 侧决定）
  bool hasSpectrum() const { return hasSpectrum_; }
  double spectrumAlphaMax() const { return alphaMax_; }
  double spectrumTg() const { return tg_; }
  int spectrumMethod() const { return method_; }       // 0 = SRSS, 1 = CQC
  double spectrumVerticalScale() const { return verticalScale_; } // 竖向 k (GB 50011 5.3.1)

  // ---- gmass：GB 50011 5.1.3 重力荷载代表值转质量 ----
  // 恒载 + ψ·活载 折算为节点重量（恒载工况自重经 selfweight 体现，
  // 转质量后 ModalAnalysis 改用"节点重量代表团"，不再叠加单元自重，
  // 避免重复计重 —— 见 ModalAnalysis.cpp buildMassAggregate）。
  bool hasGmass() const { return hasGmass_; }
  const std::string& gmassD() const { return gmassD_; }   // 恒载工况名
  const std::string& gmassL() const { return gmassL_; }   // 活载工况名
  double gmassPsi() const { return gmassPsi_; }           // 重力组合值系数（默认 0.5）

  struct ComboRequest {
    std::string kind;                 // standard/basicvar/basicperm/seismic
    std::string name;
    std::vector<std::string> cases;   // 引用的工况名（seismic 可含 Ex/Ey/Ev）
  };
  const std::vector<ComboRequest>& combos() const { return combos_; }

  const std::vector<std::string>& errors() const { return errors_; }
  const std::vector<std::string>& warnings() const { return warnings_; }
  const std::string& outStem() const { return outStem_; }

  std::string describe() const;      // 回显解析到的命令（调试用）

 private:
  struct Command {
    std::string name;
    std::vector<std::string> args;
    int line{0};
  };

  std::vector<Command> cmds_;         // 原始命令流（describe/debug 用）
  std::map<std::string, std::vector<Command>> caseLoads_;  // 工况 → 荷载命令
  std::vector<std::string> caseOrder_;                     // 工况出现顺序
  std::string currentCase_{"D"};
  bool sawCaseCmd_{false};            // 是否出现过 case 命令

  // 透传给 CLI 的网格/约束状态：CLI 只有 Model，没有 grid/spec，
  // 因此由 build() 把生成的网格暂存，applyCase() 复用同一份网格与板导荷设置。
  // GridMesh 持有 const AxisGrid& → 必须值持有 AxisGrid（否则 build 返回后悬垂）。
  AxisGrid grid_;                                // build 填充的轴网（值持有）
  std::unique_ptr<GridMesh> mesh_;               // build 生成，applyCase 复用
  GridMesh::SlabLoadOptions slabOpt_;            // slabload.to/ratio 全局设置
  GridMesh::Spec meshSpec_;                      // 网格构件规格（随 grid_ 缓存）
  bool hasGrid_{false};                          // 脚本是否定义了轴网
  bool meshBuilt_{false};                        // mesh_ 已生成（applyCase 前置）

  int modalNmodes_{0};
  bool hasModalCmd_{false};          // 脚本是否出现 modal 命令（auto 也算）
  bool hasSpectrum_{false};
  double alphaMax_{0.16};
  double tg_{0.40};
  int method_{0};
  double verticalScale_{0.65};          // 竖向地震影响系数放大系数

  // gmass：重力荷载代表值（恒载 + ψ·活载）转质量的请求参数
  bool hasGmass_{false};
  std::string gmassD_, gmassL_;
  double gmassPsi_{0.5};
  std::vector<ComboRequest> combos_;

  std::vector<std::string> errors_, warnings_;
  std::string outStem_{"yjk_out"};

  // ---- 内部工具（build 与 applyCase 共用）----
  // 用成员 grid_/meshSpec_ 生成网格并缓存到 mesh_（applyCase 复用同一网格）。
  bool ensureMeshCached(Model& m, std::string* info);
  // 执行单条荷载命令（selfweight/slabload/beamload/nodeload），
  // 错误进 errors_，警告进 warnings_。返回 false = 无法处理。
  bool applyLoadCmd(Model& m, const Command& c, std::string* info);
};

}  // namespace io
}  // namespace yjk
