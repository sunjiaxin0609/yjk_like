// =============================================================================
//  src/post/ResultExport.cpp  ——  JSON 导出 / HTML 界面生成 / 文本报告落盘
//
//  【数据格式的两条硬约束】
//
//  ① 位移按【节点序 × 6】：节点 n 的第 k 个分量在 u[n*6+k]。
//     与 StaticResult、Model::Node 完全一致，别在导出层重编号。
//
//  ② JSON 里不能出现 NaN / Inf。
//     求解发散或退化单元会产生 NaN，若不拦下，JSON.parse 会直接抛异常，
//     界面表现为"白屏"而不是"这里算错了" —— 最难排查的一类故障。
//     所以所有浮点写出前统一过 finite() 过滤。
//
//  【为什么节点/分段用数组而不是对象】
//  一个实际工程模型有几万个节点、几十万条记录。
//  对象格式 {"x":1.0,"y":2.0,...} 的键名会占掉 60% 以上的体积，
//  而这里的 schema 是固定且已知的，键名是纯粹的冗余。
//  改成定长数组后体积约 1/3，解析也更快。
//  代价是"看 JSON 不知道哪一列是什么" —— 用顶层的 _schema 字段补偿。
// =============================================================================
#include "yjk/post/ResultExport.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

#include "yjk/post/ViewerTemplate.h"

namespace yjk {
namespace post {
namespace {

// -----------------------------------------------------------------------------
//  数值写出
// -----------------------------------------------------------------------------
std::string num(double v, int prec) {
  // JSON 不承认 NaN / Infinity / -Infinity
  if (!std::isfinite(v)) return "0";
  std::ostringstream ss;
  ss << std::setprecision(prec) << v;
  return ss.str();
}

std::string esc(const std::string& s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (char c : s) {
    switch (c) {
      case '"':  out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n";  break;
      case '\r': out += "\\r";  break;
      case '\t': out += "\\t";  break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned>(c) & 0xff);
          out += buf;
        } else {
          out += c;
        }
    }
  }
  return out;
}

std::string str(const std::string& s) { return "\"" + esc(s) + "\""; }

struct Arr {
  std::string s{"["};
  bool first{true};
  void d(double v, int p) {
    if (!first) s += ",";
    first = false;
    s += num(v, p);
  }
  void i(long long v) {
    if (!first) s += ",";
    first = false;
    s += std::to_string(v);
  }
  std::string end() const { return s + "]"; }
};

bool writeFile(const std::string& path, const std::string& content, std::string* err) {
  std::ofstream of(path, std::ios::binary | std::ios::trunc);
  if (!of) {
    if (err) *err = "无法打开文件写入：" + path;
    return false;
  }
  of.write(content.data(), static_cast<std::streamsize>(content.size()));
  of.flush();
  if (!of) {
    if (err) *err = "写入失败（磁盘满？）：" + path;
    return false;
  }
  return true;
}

}  // namespace

// -----------------------------------------------------------------------------
//  JSON
//
//  schema 说明（定长数组的含义）—— 界面与后续读者都以此为准
//    node    : [x, y, z, story, fixedMask, ux, uy, uz, rx, ry, rz, Rx, Ry, Rz]
//    station : [xi, N, Vy, Vz, T, My, Mz, sigMax, sigMin, vonMises]
//    beam    : [id, ni, nj, kind, L, sigMax, sigMin, vmMax, util, Mmax, Vmax,
//               Nmax, Nmin, weight, stations[]]
//    shell   : [id, n0, n1, n2, n3, area, t, sx, sy, txy, s1, s2, vm, weight]
// -----------------------------------------------------------------------------
std::string toJson(const Model& m, const StaticResult& r, const PostProcessor& pp,
                   const ExportOptions& opt, const std::string& title) {
  const int p = opt.precision;
  std::string s;
  s.reserve(1 << 20);

  s += "{\n";
  s += "  \"title\": " + str(title) + ",\n";
  s += "  \"_schema\": {\n";
  s += "    \"node\": [\"x\",\"y\",\"z\",\"story\",\"fixedMask\",\"ux\",\"uy\",\"uz\","
       "\"rx\",\"ry\",\"rz\",\"Rx\",\"Ry\",\"Rz\"],\n";
  s += "    \"station\": [\"xi\",\"N\",\"Vy\",\"Vz\",\"T\",\"My\",\"Mz\","
       "\"sigMax\",\"sigMin\",\"vonMises\"],\n";
  s += "    \"beam\": [\"id\",\"ni\",\"nj\",\"kind\",\"L\",\"sigMax\",\"sigMin\","
       "\"vmMax\",\"util\",\"Mmax\",\"Vmax\",\"Nmax\",\"Nmin\",\"weight\",\"stations\"],\n";
  s += "    \"shell\": [\"id\",\"n0\",\"n1\",\"n2\",\"n3\",\"area\",\"t\","
       "\"sx\",\"sy\",\"txy\",\"s1\",\"s2\",\"vm\",\"weight\"],\n";
  s += "    \"note\": \"位移与反力均按 节点序×6；fixedMask 位0..5 = ux..rz\"\n";
  s += "  },\n";

  // ---- meta ----
  s += "  \"meta\": {\n";
  s += "    \"ok\": " + std::string(r.ok ? "true" : "false") + ",\n";
  s += "    \"message\": " + str(r.message) + ",\n";
  s += "    \"ndof\": " + std::to_string(r.ndof) + ",\n";
  s += "    \"nconstrained\": " + std::to_string(r.nconstrained) + ",\n";
  s += "    \"nnz\": " + std::to_string(r.nnz) + ",\n";
  s += "    \"nnzL\": " + std::to_string(r.nnzL) + ",\n";
  s += "    \"seconds\": " + num(r.seconds, 4) + ",\n";
  s += "    \"residual\": " + num(r.residual, 6) + ",\n";
  s += "    \"nodeCount\": " + std::to_string(m.nodeCount()) + ",\n";
  s += "    \"elementCount\": " + std::to_string(m.elementCount()) + "\n";
  s += "  },\n";

  // ---- nodes ----
  s += "  \"nodes\": [\n";
  const Id nn = m.nodeCount();
  for (Id i = 0; i < nn; ++i) {
    const Node& nd = m.node(i);
    int mask = 0;
    for (int k = 0; k < 6; ++k) if (nd.fixed[k]) mask |= (1 << k);
    Arr a;
    a.d(nd.r.x, p); a.d(nd.r.y, p); a.d(nd.r.z, p);
    a.i(nd.story); a.i(mask);
    for (int k = 0; k < 6; ++k) a.d(r.u[static_cast<size_t>(i) * 6 + k], p);
    for (int k = 0; k < 3; ++k) a.d(r.reaction[static_cast<size_t>(i) * 6 + k], p);
    s += "    " + a.end();
    if (i + 1 < nn) s += ",";
    s += "\n";
  }
  s += "  ],\n";

  // ---- beams ----
  s += "  \"beams\": [\n";
  const auto& beams = pp.beamResults();
  for (size_t bi = 0; bi < beams.size(); ++bi) {
    const BeamResult& b = beams[bi];
    Arr a;
    a.i(b.id); a.i(b.ni); a.i(b.nj);
    a.i(static_cast<int>(b.kind));
    a.d(b.L, p);
    a.d(b.sigMax, p); a.d(b.sigMin, p); a.d(b.vonMisesMax, p); a.d(b.util, p);
    a.d(b.Mmax, p); a.d(b.Vmax, p); a.d(b.Nmax, p); a.d(b.Nmin, p);
    a.d(b.weight, p);
    // 分段：嵌套数组
    {
      std::string inner{"["};
      if (opt.includeStations) {
        const size_t take = static_cast<size_t>(opt.stationsPerMember);
        const size_t ns = b.stations.size();
        for (size_t k = 0; k < ns; ++k) {
          // 采样：优先取端点与跨中
          if (ns > take && k % ((ns + take - 1) / take) != 0 && k + 1 != ns) continue;
          const BeamStation& st = b.stations[k];
          Arr t;
          t.d(st.xi, 4);
          t.d(st.N, p); t.d(st.Vy, p); t.d(st.Vz, p); t.d(st.T, p);
          t.d(st.My, p); t.d(st.Mz, p);
          t.d(st.sigMax, p); t.d(st.sigMin, p); t.d(st.vonMises, p);
          if (inner.size() > 1) inner += ",";
          inner += t.end();
        }
      }
      inner += "]";
      if (!a.first) a.s += ",";
      a.first = false;
      a.s += inner;
    }
    s += "    " + a.end();
    if (bi + 1 < beams.size()) s += ",";
    s += "\n";
  }
  s += "  ],\n";

  // ---- shells ----
  s += "  \"shells\": [\n";
  const auto& shells = pp.shellResults();
  for (size_t si = 0; si < shells.size(); ++si) {
    const ShellResult& q = shells[si];
    Arr a;
    a.i(q.id);
    for (int k = 0; k < 4; ++k) a.i(q.nodes[static_cast<size_t>(k)]);
    a.d(q.area, p); a.d(q.thickness, p);
    a.d(q.mf.sx, p); a.d(q.mf.sy, p); a.d(q.mf.txy, p);
    a.d(q.mf.s1, p); a.d(q.mf.s2, p); a.d(q.mf.vonMises, p);
    a.d(q.weight, p);
    s += "    " + a.end();
    if (si + 1 < shells.size()) s += ",";
    s += "\n";
  }
  s += "  ],\n";

  // ---- node fields（云图）----
  s += "  \"fields\": [\n";
  if (opt.includeNodeFields) {
    const auto& fs = pp.nodeFields();
    for (size_t fi = 0; fi < fs.size(); ++fi) {
      const NodeField& f = fs[fi];
      s += "    {\"name\":" + str(f.name) + ",\"unit\":" + str(f.unit) +
           ",\"min\":" + num(f.vmin, p) + ",\"max\":" + num(f.vmax, p) + ",\"v\":[";
      for (size_t i = 0; i < f.v.size(); ++i) {
        if (i) s += ",";
        s += num(f.v[i], p);
      }
      s += "]}";
      if (fi + 1 < fs.size()) s += ",";
      s += "\n";
    }
  }
  s += "  ],\n";

  // ---- stories ----
  s += "  \"stories\": [\n";
  const auto& st = pp.stories();
  for (size_t i = 0; i < st.size(); ++i) {
    const StoryResult& q = st[i];
    s += "    {\"story\":" + std::to_string(q.story) +
         ",\"z\":" + num(q.z, p) + ",\"h\":" + num(q.h, p) +
         ",\"nodeCount\":" + std::to_string(q.nodeCount) +
         ",\"maxUxy\":" + num(q.maxUxy, p) + ",\"avgUxy\":" + num(q.avgUxy, p) +
         ",\"maxUz\":" + num(q.maxUz, p) +
         ",\"drift\":" + num(q.drift, p) + ",\"driftAvg\":" + num(q.driftAvg, p) +
         ",\"driftNode\":" + std::to_string(q.driftNode) +
         ",\"driftRatio\":" + num(q.driftRatio, p) +
         ",\"shear\":" + num(q.shear, p) + ",\"weight\":" + num(q.weight, p) +
         ",\"shearWeightRatio\":" + num(q.shearWeightRatio, p) + "}";
    if (i + 1 < st.size()) s += ",";
    s += "\n";
  }
  s += "  ],\n";

  // ---- envelope ----
  const Envelope& e = pp.envelope();
  s += "  \"envelope\": {\n";
  s += "    \"maxDisp\":" + num(e.maxDisp, p) + ",\"maxDispNode\":" + std::to_string(e.maxDispNode) + ",\n";
  s += "    \"maxUxy\":" + num(e.maxUxy, p) + ",\"maxUxyNode\":" + std::to_string(e.maxUxyNode) + ",\n";
  s += "    \"maxUz\":" + num(e.maxUz, p) + ",\"maxUzNode\":" + std::to_string(e.maxUzNode) + ",\n";
  s += "    \"maxVonMises\":" + num(e.maxVonMises, p) + ",\"maxVonMisesNode\":" +
       std::to_string(e.maxVonMisesNode) + ",\n";
  s += "    \"maxReaction\":" + num(e.maxReaction, p) + ",\"maxReactionNode\":" +
       std::to_string(e.maxReactionNode) + ",\n";
  s += "    \"maxReactionZ\":" + num(e.maxReactionZ, p) + ",\"maxReactionZNode\":" +
       std::to_string(e.maxReactionZNode) + ",\n";
  s += "    \"totalWeight\":" + num(e.totalWeight, p) +
       ",\"baseShear\":" + num(e.baseShear, p) +
       ",\"shearWeightRatio\":" + num(e.shearWeightRatio, p) + ",\n";
  s += "    \"maxDrift\":" + num(e.maxDrift, p) + ",\"maxDriftStory\":" +
       std::to_string(e.maxDriftStory) + "\n";
  s += "  },\n";

  s += "  \"memberKindNames\": [\"柱\",\"梁\",\"支撑\",\"其它\"],\n";
  s += "  \"summary\": " + str(pp.summary()) + "\n";
  s += "}\n";
  return s;
}

bool writeJson(const std::string& path, const Model& m, const StaticResult& r,
               const PostProcessor& pp, const ExportOptions& opt,
               const std::string& title, std::string* err) {
  return writeFile(path, toJson(m, r, pp, opt, title), err);
}

// -----------------------------------------------------------------------------
//  HTML 界面
//
//  【为什么用替换占位符而不是把 JSON 直接拼进模板字符串】
//  JSON 里可能含 "</script>"（节点 label 里完全可能出现），
//  直接拼接会提前闭合脚本标签 → 界面白屏。
//  这里把模板里的 __YJK_DATA__ 替换成 JSON，并对 "</" 做转义，
//  使其不可能闭合标签。
// -----------------------------------------------------------------------------
namespace {
std::string jsonSafeForScript(const std::string& j) {
  std::string out;
  out.reserve(j.size() + 64);
  for (size_t i = 0; i < j.size(); ++i) {
    if (j[i] == '<' && i + 1 < j.size() && (j[i + 1] == '/' )) {
      out += "<\\/";
      ++i;
    } else {
      out += j[i];
    }
  }
  return out;
}
}  // namespace

std::string htmlViewer(const std::string& json, const std::string& title) {
  std::string html = kViewerTemplate;
  const std::string key = "__YJK_DATA__";
  const size_t pos = html.find(key);
  if (pos == std::string::npos) return html;         // 模板被改坏了：原样返回
  html.replace(pos, key.size(), jsonSafeForScript(json));

  const std::string tkey = "__YJK_TITLE__";
  size_t p2 = html.find(tkey);
  while (p2 != std::string::npos) {
    html.replace(p2, tkey.size(), esc(title));
    p2 = html.find(tkey, p2 + title.size());
  }
  return html;
}

bool writeHtmlViewer(const std::string& path, const std::string& json,
                     const std::string& title, std::string* err) {
  return writeFile(path, htmlViewer(json, title), err);
}

bool writeReport(const std::string& path, const PostProcessor& pp, std::string* err) {
  return writeFile(path, pp.report(), err);
}

ExportFiles exportAll(const std::string& stem, const Model& m, const StaticResult& r,
                      const PostProcessor& pp, const ExportOptions& opt) {
  ExportFiles out;
  const std::string title = "yjk_like 分析结果  [" + stem + "]";
  out.jsonPath = stem + ".json";
  out.htmlPath = stem + ".html";
  out.txtPath = stem + ".txt";

  const std::string json = toJson(m, r, pp, opt, title);
  std::string err;
  if (!writeFile(out.jsonPath, json, &err)) { out.error = err; return out; }
  if (!writeFile(out.htmlPath, htmlViewer(json, title), &err)) { out.error = err; return out; }
  if (!writeFile(out.txtPath, pp.report(), &err)) { out.error = err; return out; }
  out.ok = true;
  return out;
}

// =============================================================================
//  动态分析导出（模态 / 反应谱 / 多工况 / 荷载组合）
// =============================================================================
namespace {

// 楼层结果 → JSON 对象（与静态 toJson 的 stories 段同构）
std::string storiesJson(const PostProcessor& pp, int p) {
  std::string s;
  const auto& st = pp.stories();
  for (size_t i = 0; i < st.size(); ++i) {
    const StoryResult& q = st[i];
    if (i) s += ",";
    s += "{\"story\":" + std::to_string(q.story) +
         ",\"z\":" + num(q.z, p) + ",\"h\":" + num(q.h, p) +
         ",\"nodeCount\":" + std::to_string(q.nodeCount) +
         ",\"maxUxy\":" + num(q.maxUxy, p) + ",\"avgUxy\":" + num(q.avgUxy, p) +
         ",\"maxUz\":" + num(q.maxUz, p) +
         ",\"drift\":" + num(q.drift, p) + ",\"driftAvg\":" + num(q.driftAvg, p) +
         ",\"driftNode\":" + std::to_string(q.driftNode) +
         ",\"driftRatio\":" + num(q.driftRatio, p) +
         ",\"shear\":" + num(q.shear, p) + ",\"weight\":" + num(q.weight, p) +
         ",\"shearWeightRatio\":" + num(q.shearWeightRatio, p) + "}";
  }
  return s;
}

// 杆件内力极值表（组合后内力）：[id, kind, L, Mmax, Vmax, Nmax, Nmin, sigMax]
std::string beamExtremesJson(const PostProcessor& pp, int p) {
  std::string s;
  const auto& beams = pp.beamResults();
  for (size_t i = 0; i < beams.size(); ++i) {
    const BeamResult& b = beams[i];
    if (i) s += ",";
    Arr a;
    a.i(b.id); a.i(static_cast<int>(b.kind)); a.d(b.L, p);
    a.d(b.Mmax, p); a.d(b.Vmax, p); a.d(b.Nmax, p); a.d(b.Nmin, p);
    a.d(b.sigMax, p);
    s += a.end();
  }
  return s;
}

// 单个静力结果的摘要（多工况/组合共用）
std::string resultSummaryJson(const StaticResult& r, int p) {
  std::string s;
  s += "{\"ok\":" + std::string(r.ok ? "true" : "false") +
       ",\"ndof\":" + std::to_string(r.ndof) +
       ",\"maxDisplacement\":" + num(r.maxDisplacement, p) +
       ",\"maxUxy\":" + num(r.maxUxy, p) +
       ",\"maxUz\":" + num(r.maxUz, p) +
       ",\"maxRotation\":" + num(r.maxRotation, p) +
       ",\"maxReaction\":" + num(r.maxReaction, p) +
       ",\"residual\":" + num(r.residual, p) +
       ",\"seconds\":" + num(r.seconds, p) + "}";
  return s;
}

}  // namespace

std::string toJsonDynamic(const Model& m, const DynamicExport& dyn,
                          const ExportOptions& opt, const std::string& title) {
  const int p = opt.precision;
  std::string s;
  s.reserve(1 << 20);

  s += "{\n";
  s += "  \"title\": " + str(title) + ",\n";
  s += "  \"_schema\": \"yjk.dynamic\",\n";
  s += "  \"_schemaNote\": \"cases/combos 各自携带 stories（层剪力）与摘要；"
       "combos 另含 beamExtremes=[id,kind,L,Mmax,Vmax,Nmax,Nmin,sigMax]；"
       "gravityByStory=[story,W] 为剪重比 V/W 的重力代表值分母\",\n";

  // ---- gravityByStory：重力代表值（剪重比分母）----
  s += "  \"gravityByStory\": [";
  for (size_t i = 0; i < dyn.gravityByStory.size(); ++i) {
    if (i) s += ",";
    s += "[" + std::to_string(dyn.gravityByStory[i].first) + "," +
         num(dyn.gravityByStory[i].second, p) + "]";
  }
  s += "],\n";

  // ---- modal ----
  if (dyn.modal) {
    const ModalResult& q = *dyn.modal;
    s += "  \"modal\": {\n";
    s += "    \"ok\": " + std::string(q.ok ? "true" : "false") + ",\n";
    s += "    \"message\": " + str(q.message) + ",\n";
    s += "    \"nmodes\": " + std::to_string(q.nmodes) + ",\n";
    s += "    \"ndof\": " + std::to_string(q.ndof) + ",\n";
    s += "    \"iter\": " + std::to_string(q.iter) + ",\n";
    s += "    \"residual\": " + num(q.residual, 6) + ",\n";
    s += "    \"seconds\": " + num(q.seconds, 4) + ",\n";
    const char* tab = "    ";
    auto vecJson = [&](const std::vector<double>& v, const char* key) {
      s += tab; s += key; s += ": [";
      for (size_t i = 0; i < v.size(); ++i) { if (i) s += ","; s += num(v[i], p); }
      s += "],\n";
    };
    vecJson(q.freq, "\"freq\"");
    vecJson(q.period, "\"period\"");
    vecJson(q.omega, "\"omega\"");
    vecJson(q.gammaX, "\"gammaX\"");
    vecJson(q.gammaY, "\"gammaY\"");
    vecJson(q.gammaZ, "\"gammaZ\"");
    vecJson(q.mstarX, "\"mstarX\"");
    vecJson(q.mstarY, "\"mstarY\"");
    vecJson(q.mstarZ, "\"mstarZ\"");
    s += "    \"massX\": " + num(q.massX, p) + ",\n";
    s += "    \"massY\": " + num(q.massY, p) + ",\n";
    s += "    \"massZ\": " + num(q.massZ, p) + ",\n";
    s += "    \"totalMassX\": " + num(q.totalMassX, p) + ",\n";
    s += "    \"totalMassY\": " + num(q.totalMassY, p) + ",\n";
    s += "    \"totalMassZ\": " + num(q.totalMassZ, p) + ",\n";
    s += "    \"ratioX\": " + num(q.ratioX, p) + ",\n";
    s += "    \"ratioY\": " + num(q.ratioY, p) + ",\n";
    s += "    \"ratioZ\": " + num(q.ratioZ, p) + "\n";
    s += "  },\n";
  } else {
    s += "  \"modal\": null,\n";
  }

  // ---- spectrum ----
  if (dyn.spectrum) {
    const SpectrumResult& q = *dyn.spectrum;
    s += "  \"spectrum\": {\n";
    s += "    \"ok\": " + std::string(q.ok ? "true" : "false") + ",\n";
    s += "    \"message\": " + str(q.message) + ",\n";
    s += "    \"nmodes\": " + std::to_string(q.nmodes) + ",\n";
    s += "    \"baseX\": " + num(q.baseX, p) + ",\n";
    s += "    \"baseY\": " + num(q.baseY, p) + ",\n";
    s += "    \"baseZ\": " + num(q.baseZ, p) + ",\n";
    s += "    \"massRatioX\": " + num(q.massRatioX, p) + ",\n";
    s += "    \"massRatioY\": " + num(q.massRatioY, p) + ",\n";
    s += "    \"massRatioZ\": " + num(q.massRatioZ, p) + ",\n";
    s += "    \"massOkX\": " + std::string(q.massOkX ? "true" : "false") + ",\n";
    s += "    \"massOkY\": " + std::string(q.massOkY ? "true" : "false") + ",\n";
    s += "    \"massOkZ\": " + std::string(q.massOkZ ? "true" : "false") + ",\n";
    s += "    \"massOk\": " + std::string(q.massOk ? "true" : "false") + "\n";
    s += "  },\n";
  } else {
    s += "  \"spectrum\": null,\n";
  }

  // ---- cases（各工况：摘要 + 层剪力）----
  s += "  \"cases\": [\n";
  for (size_t i = 0; i < dyn.cases.size(); ++i) {
    const auto& nc = dyn.cases[i];
    PostProcessor pp(m, *nc.r);
    pp.compute();
    s += "    {\"name\":" + str(nc.name) + ",";
    s += "\"summary\":" + resultSummaryJson(*nc.r, p) + ",";
    s += "\"stories\":[" + storiesJson(pp, p) + "]}";
    if (i + 1 < dyn.cases.size()) s += ",";
    s += "\n";
  }
  s += "  ],\n";

  // ---- combos（组合：摘要 + 层剪力 + 杆内力极值）----
  s += "  \"combos\": [\n";
  for (size_t i = 0; i < dyn.combos.size(); ++i) {
    const auto& nc = dyn.combos[i];
    PostProcessor pp(m, *nc.r);
    pp.compute();
    s += "    {\"name\":" + str(nc.name) + ",";
    s += "\"summary\":" + resultSummaryJson(*nc.r, p) + ",";
    s += "\"stories\":[" + storiesJson(pp, p) + "],";
    s += "\"beamExtremes\":[" + beamExtremesJson(pp, p) + "]}";
    if (i + 1 < dyn.combos.size()) s += ",";
    s += "\n";
  }
  s += "  ],\n";

  // ---- envelope（组合包络）----
  const auto envStoryJson = [&](const StaticResult* r) -> std::string {
    PostProcessor pp(m, *r);
    pp.compute();
    return storiesJson(pp, p);
  };
  if (dyn.envMax && dyn.envMin) {
    s += "  \"envelope\": {\n";
    s += "    \"max\": {\"summary\":" + resultSummaryJson(*dyn.envMax, p) +
         ",\"stories\":[" + envStoryJson(dyn.envMax) + "]},\n";
    s += "    \"min\": {\"summary\":" + resultSummaryJson(*dyn.envMin, p) +
         ",\"stories\":[" + envStoryJson(dyn.envMin) + "]}\n";
    s += "  },\n";
  } else {
    s += "  \"envelope\": null,\n";
  }

  s += "  \"nodeCount\": " + std::to_string(m.nodeCount()) + ",\n";
  s += "  \"memberKindNames\": [\"柱\",\"梁\",\"支撑\",\"其它\"],\n";
  s += "  \"summary\": \"动态分析：模态/反应谱/多工况组合\"\n";
  s += "}\n";
  return s;
}

std::string dynamicReport(const Model& m, const DynamicExport& dyn,
                          const ExportOptions& opt) {
  const int p = opt.precision;
  std::ostringstream os;

  // 重力代表值（GB 50011 5.1.3 恒载 + ψ·活载）：各层以上合计 W 为剪重比
  // V/W 的分母。无基准时回退用各工况自身荷载统计的比值。
  std::map<int, double> gravWeight;   // story → 该层及以上重力代表值 kN
  for (const auto& kv : dyn.gravityByStory) gravWeight[kv.first] = kv.second;
  const auto ratio = [&gravWeight](const StoryResult& q) -> double {
    const auto it = gravWeight.find(q.story);
    if (it != gravWeight.end() && it->second > 1e-9) return q.shear / it->second;
    return q.shearWeightRatio;
  };

  const auto rep = [&](const std::vector<DynamicExport::NamedResult>& v,
                       const char* what) {
    for (const auto& nc : v) {
      if (!nc.r || !nc.r->ok) { os << what << " " << nc.name << ": 求解失败\n"; continue; }
      PostProcessor pp(m, *nc.r);
      pp.compute();
      os << what << " " << nc.name << ": 最大位移 " << nc.r->maxDisplacement << " m  "
         << "Uxy " << nc.r->maxUxy << "  反力 " << nc.r->maxReaction << " kN\n";
      for (const StoryResult& q : pp.stories())
        os << "    层 " << q.story << ": 剪力 " << q.shear << " kN  位移角 "
           << q.drift << "  剪重比 " << ratio(q)
           << "  重力W " << (gravWeight.count(q.story)
                                 ? gravWeight.at(q.story)
                                 : q.weight)
           << " kN  maxUz " << q.maxUz << " m\n";
    }
  };

  os << "================================================================\n";
  os << " yjk_like 动态分析结果（模态 / 反应谱 / 多工况 / 荷载组合）\n";
  os << "================================================================\n";

  if (!gravWeight.empty()) {
    os << "\n[重力代表值]（GB 50011 5.1.3 恒载 + ψ·活载，剪重比 V/W 分母）\n";
    double totalW = 0.0;
    for (const auto& kv : dyn.gravityByStory)
      if (kv.second > totalW) totalW = kv.second;
    for (const auto& kv : dyn.gravityByStory)
      os << "  层 " << kv.first << ": 及以上 W = " << kv.second << " kN\n";
    os << "  整楼重力代表值 W = " << totalW << " kN\n";
  }

  if (dyn.modal) {
    const ModalResult& q = *dyn.modal;
    os << "\n[模态分析]  " << (q.ok ? "成功" : "失败：" + q.message)
       << "  ndof=" << q.ndof << "  迭代 " << q.iter << "  残差 " << q.residual << "\n";
    if (q.ok) {
      os << "  阶   频率Hz      周期s       ωrad/s     参与质量比 X / Y / Z\n";
      for (int j = 0; j < q.nmodes; ++j) {
        os << "  " << (j + 1) << "   " << q.freq[j] << "   " << q.period[j]
           << "   " << q.omega[j] << "   "
           << q.gammaX[j] << " / " << q.gammaY[j] << " / " << q.gammaZ[j] << "\n";
      }
      os << "  总质量 X=" << q.totalMassX << " Y=" << q.totalMassY
         << " Z=" << q.totalMassZ << " t\n";
      os << "  参与质量累计  X " << q.ratioX * 100.0 << "%   Y " << q.ratioY * 100.0
         << "%   Z " << q.ratioZ * 100.0 << "%\n";
    }
  }

  if (dyn.spectrum) {
    const SpectrumResult& q = *dyn.spectrum;
    os << "\n[反应谱]  " << (q.ok ? "成功" : "失败：" + q.message) << "\n";
    // 谱参数不在 SpectrumResult 里 —— 由 CLI 侧在 message 里带出，这里打印结果重点
    os << "  基底剪力  X=" << q.baseX << "  Y=" << q.baseY
       << "  Z=" << q.baseZ << " kN\n";
    os << "  参与质量比  X " << q.massRatioX * 100.0 << "%"
       << (q.massOkX ? " ✓" : " ✗") << "   Y " << q.massRatioY * 100.0 << "%"
       << (q.massOkY ? " ✓" : " ✗") << "   Z " << q.massRatioZ * 100.0 << "%"
       << (q.massOkZ ? " ✓" : " ✗") << "\n";
    os << "  参与质量 ≥90% 检查: " << (q.massOk ? "通过" : "未满足（请增加模态阶数）") << "\n";
  }

  os << "\n[静力工况]（层剪力 / 层间位移角 / 竖向位移）\n";
  rep(dyn.cases, "工况");

  os << "\n[荷载组合]（层剪力与组合后控制内力）\n";
  for (const auto& nc : dyn.combos) {
    if (!nc.r || !nc.r->ok) { os << "组合 " << nc.name << ": 求解失败\n"; continue; }
    PostProcessor pp(m, *nc.r);
    pp.compute();
    os << "\n  组合 " << nc.name << ": 最大位移 " << nc.r->maxDisplacement
       << " m  Uxy " << nc.r->maxUxy << "  反力 " << nc.r->maxReaction << " kN\n";
    for (const StoryResult& q : pp.stories())
      os << "    层 " << q.story << ": 剪力 " << q.shear << " kN  位移角 "
         << q.drift << "  剪重比 " << ratio(q)
         << "  重力W " << (gravWeight.count(q.story) ? gravWeight.at(q.story)
                                                     : q.weight)
         << " kN\n";
    os << "    组合后内力极值（Mmax / Vmax / Nmax，绝对值）:\n";
    for (const BeamResult& b : pp.beamResults())
      os << "      杆 " << b.id << "(" << memberKindName(b.kind) << ") "
         << "Mmax=" << b.Mmax << "  Vmax=" << b.Vmax << "  Nmax=" << b.Nmax
         << "  Nmin=" << b.Nmin << "  σmax=" << b.sigMax << " kPa\n";
  }

  if (dyn.envMax && dyn.envMin) {
    PostProcessor ppx(m, *dyn.envMax);
    ppx.compute();
    os << "\n[包络]  maxDisp=" << dyn.envMax->maxDisplacement
       << "  maxUxy=" << dyn.envMax->maxUxy
       << "  maxUz=" << dyn.envMax->maxUz
       << "  maxReaction=" << dyn.envMax->maxReaction << "\n";
  }
  os << "\n================================================================\n";
  return os.str();
}

DynamicExportFiles exportDynamic(const std::string& stem, const Model& m,
                                 const DynamicExport& dyn, const ExportOptions& opt) {
  DynamicExportFiles out;
  const std::string title = "yjk_like 动态分析结果  [" + stem + "]";
  const std::string json = toJsonDynamic(m, dyn, opt, title);
  std::string err;
  out.jsonPath = stem + ".dynamic.json";
  out.txtPath = stem + ".dynamic.txt";
  if (!writeFile(out.jsonPath, json, &err)) { out.error = err; return out; }
  if (!writeFile(out.txtPath, dynamicReport(m, dyn, opt), &err)) {
    out.error = err;
    return out;
  }
  out.ok = true;
  return out;
}

}  // namespace post
}  // namespace yjk
