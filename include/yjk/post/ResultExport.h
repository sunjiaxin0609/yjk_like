// =============================================================================
//  include/yjk/post/ResultExport.h  ——  结果导出（JSON）与可视化界面生成
//
//  【为什么要自己写 JSON 而不引第三方库】
//  内核的定位是"零外部依赖"。为一个只写不读的格式引入 nlohmann/json
//  会带来几万行头文件依赖，而这里需要的只是"把几个数组按固定 schema 写出去"。
//  手写 200 行足够，而且 schema 完全可控。
//
//  【界面层为什么生成 HTML 而不是 Qt】
//    · 内核可以在无 GUI 的环境（CI / 服务器）里跑，生成的结果文件
//      用浏览器打开就能看 —— 不需要编译 Qt
//    · 依赖为零：Canvas 2D 手绘三维投影，不引 three.js（离线可用）
//    · 结果文件自带数据，可以邮件/微信直接发，不丢上下文
// =============================================================================
#pragma once

#include <string>
#include <vector>

#include "yjk/analysis/ModalAnalysis.h"
#include "yjk/analysis/ResponseSpectrum.h"
#include "yjk/post/PostProcessor.h"

namespace yjk {
namespace post {

struct ExportOptions {
  int stationsPerMember{21};
  bool includeStations{true};     // 是否导出每根杆沿长的分段内力（内力图要用）
  bool includeNodeFields{true};
  int precision{6};
};

// 生成结果 JSON（UTF-8，无 BOM）
std::string toJson(const Model& m, const StaticResult& r, const PostProcessor& pp,
                   const ExportOptions& opt = {}, const std::string& title = "yjk_like 分析结果");

bool writeJson(const std::string& path, const Model& m, const StaticResult& r,
               const PostProcessor& pp, const ExportOptions& opt = {},
               const std::string& title = "yjk_like 分析结果", std::string* err = nullptr);

// 生成自带数据的 HTML 后处理界面（离线可用）
//   json 会被整体内嵌到 <script> 里，替换模板中的 __YJK_DATA__ 占位符。
std::string htmlViewer(const std::string& json, const std::string& title = "yjk_like 后处理");

bool writeHtmlViewer(const std::string& path, const std::string& json,
                     const std::string& title = "yjk_like 后处理", std::string* err = nullptr);

// 纯文本报告（写文件用）
bool writeReport(const std::string& path, const PostProcessor& pp, std::string* err = nullptr);

// 便捷入口：一次写出 json / html / txt
struct ExportFiles {
  std::string jsonPath, htmlPath, txtPath;
  bool ok{false};
  std::string error;
};
ExportFiles exportAll(const std::string& stem, const Model& m, const StaticResult& r,
                      const PostProcessor& pp, const ExportOptions& opt = {});

// -----------------------------------------------------------------------------
//  动态分析导出（模态 / 反应谱 / 多工况与荷载组合）
//
//  静态单工况（旧接口）只输出"一个工况"的结果；动态流程要同时呈现：
//    · 模态：固有频率 / 周期 / 参与质量
//    · 反应谱：基底剪力 / 参与质量 ≥ 90% 检查
//    · 各工况与各组合：层剪力（stories）、组合后内力（杆端/跨中极值）
//  DynamicExport 只持有指针，不拷贝结果 —— 调用方（CLI）保证其生命周期。
//  每个工况/组合的内部后处理（PostProcessor）在导出函数内部完成。
// -----------------------------------------------------------------------------
struct DynamicExport {
  struct NamedResult {
    std::string name;
    const StaticResult* r{nullptr};
  };

  const ModalResult* modal{nullptr};        // 可选
  const SpectrumResult* spectrum{nullptr};   // 可选
  // 重力代表值（恒载 + ψ·活载，GB 50011 5.1.3），按层存放：
  //   first  = 层号，second = 该层及以上合计重力代表值 kN。
  // 供各工况/组合的剪重比作分母（V / W）。为空时无剪重比基准。
  std::vector<std::pair<int, double>> gravityByStory;
  std::vector<NamedResult> cases;            // 静力工况（D/L/W/Ex/Ey…，按 caseNames 顺序）
  std::vector<NamedResult> combos;           // 组合结果
  const StaticResult* envMax{nullptr};       // 可选：包络极值工况
  const StaticResult* envMin{nullptr};
};

// 动态 JSON（UTF-8，无 BOM）。包含：
//   modal / spectrum / cases 与 combos 的层剪力表、combos 的杆件内力极值表、
//   envelope 包络段。位移按 节点序×6（与静态 toJson 一致）。
std::string toJsonDynamic(const Model& m, const DynamicExport& dyn,
                          const ExportOptions& opt = {},
                          const std::string& title = "yjk_like 动态分析结果");

// 纯文本报告：模态表 → 反应谱摘要 → 工况/组合层剪力 → 组合内力极值 → 包络
std::string dynamicReport(const Model& m, const DynamicExport& dyn,
                          const ExportOptions& opt = {});

struct DynamicExportFiles {
  std::string jsonPath, txtPath;
  bool ok{false};
  std::string error;
};
// 一次写出 <stem>.dynamic.json 与 <stem>.dynamic.txt（不覆盖静态三件套）
DynamicExportFiles exportDynamic(const std::string& stem, const Model& m,
                                 const DynamicExport& dyn,
                                 const ExportOptions& opt = {});

}  // namespace post
}  // namespace yjk
