// =============================================================================
//  apps/yjk.cpp  ——  命令行入口
//
//    yjk run   <模型.yjk> [-o 输出前缀]      建模 → 校核 → 求解 → 后处理 → 导出
//    yjk check <模型.yjk>                    只做建模与校核（不求解）
//    yjk help
//
//  【为什么命令行而不是先做 GUI】
//  求解器正确性的验证靠"批量跑算例、对解析解"，这个过程必须能脚本化。
//  GUI 在这种工作流里是负担。命令行到位后，GUI 只是一个壳。
//
//  【脚本里的动态分析命令（case/modal/spectrum/combo）】
//  .yjk 脚本一旦出现这些命令，build() 只建立几何与约束（不施加荷载），
//  本程序随后按下列流水线执行：
//    ① 每个工况：clearLoads() → applyCase() 施加工况荷载 → 静力求解
//    ② modal：模态分析（频率/周期/参与质量）
//    ③ spectrum：反应谱 → 等效地震力 → 作为 Ex/Ey(/Ez) 独立工况再静力求解
//    ④ combo：对既有工况结果做线性组合 + 包络
//    ⑤ 导出：<前缀>{,.json,.html,.txt}（主工况三件套，同旧版）
//            + <前缀>.dynamic.json / .dynamic.txt（模态/谱/层剪力/组合内力）
//
//  【退出码】
//    0 成功   1 脚本/建模错误   2 校核未通过   3 求解失败   4 用法错误
//  明确的退出码让批量脚本能区分"模型写错了"和"算不出来"。
// =============================================================================
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "yjk/analysis/LoadCombination.h"
#include "yjk/analysis/ModalAnalysis.h"
#include "yjk/analysis/ResponseSpectrum.h"
#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/io/ModelScript.h"
#include "yjk/model/Model.h"
#include "yjk/post/PostProcessor.h"
#include "yjk/post/ResultExport.h"

namespace {

void usage() {
  std::printf(
      "yjk_like —— 结构有限元分析内核\n"
      "\n"
      "用法:\n"
      "  yjk run   <模型.yjk> [-o 输出前缀]    建模 → 校核 → 求解 → 后处理 → 导出\n"
      "  yjk check <模型.yjk>                  只建模与校核，不求解\n"
      "  yjk help\n"
      "\n"
      "脚本可含动态分析命令:\n"
      "  case <工况名>                    切换后续荷载所属工况（缺省 D）\n"
      "  modal <nmodes>                   请求模态分析\n"
      "  spectrum <alphaMax> <Tg> [srss|cqc]  请求反应谱分析\n"
      "  gmass <恒载工况> <活载工况> [ψ=0.5]  重力荷载代表值(恒+ψ活)转节点质量\n"
      "  combo <kind> <名> <case...>      荷载组合（standard/basicvar/basicperm/seismic）\n"
      "\n"
      "输出（默认前缀 = 脚本里的 out 命令，缺省 yjk_out）:\n"
      "  <前缀>.json           主工况结果数据\n"
      "  <前缀>.html           后处理界面（浏览器直接打开，离线可用）\n"
      "  <前缀>.txt            主工况文本报告\n"
      "  <前缀>.dynamic.json   动态分析数据（模态/反应谱/层剪力/组合内力）\n"
      "  <前缀>.dynamic.txt    动态分析文本报告\n"
      "\n"
      "退出码: 0 成功 / 1 脚本错误 / 2 校核未通过 / 3 求解失败 / 4 用法错误\n");
}

// 脚本 combo 命令 → LoadCombination 组合定义（规范分项系数见 LoadCombination.h）
yjk::ComboDef buildComboDef(const yjk::io::ModelScript::ComboRequest& req) {
  using namespace yjk::combo_coef;
  yjk::ComboDef d;
  d.name = req.name;
  const auto& c = req.cases;
  if (req.kind == "standard") {
    d.factors[c[0]] = 1.0;
    d.factors[c[1]] = 1.0;
  } else if (req.kind == "basicvar") {
    d.factors[c[0]] = kPermanentVar;                       // 1.2G
    d.factors[c[1]] = kVariable;                           // 1.4Q
    if (c.size() > 2) d.factors[c[2]] = kVariable * kCombFactor;  // 1.4·0.7·W
  } else if (req.kind == "basicperm") {
    d.factors[c[0]] = kPermanentFixed;                     // 1.35G
    d.factors[c[1]] = kVariable * kCombFactor;             // 1.4·0.7·Q
    if (c.size() > 2) d.factors[c[2]] = kVariable * kCombFactor;  // 1.4·0.7·W
  } else {  // seismic
    d.factors[c[0]] = kGamG;                               // 1.2G
    d.factors[c[1]] = kGamG * kPsiE;                       // 1.2·0.5·Q
    if (c.size() > 2) d.factors[c[2]] = kGamEh;            // 1.3Ex
    if (c.size() > 3) d.factors[c[3]] = kGamEv;            // 0.5Ev
  }
  return d;
}

int cmdRun(const std::string& path, const std::string& stemOverride) {
  yjk::io::ModelScript sc;
  if (!sc.parseFile(path)) {
    for (const auto& e : sc.errors()) std::printf("[脚本错误] %s\n", e.c_str());
    return 1;
  }

  yjk::Model m;
  std::string info;
  if (!sc.build(m, &info)) {
    for (const auto& e : sc.errors()) std::printf("[建模错误] %s\n", e.c_str());
    return 1;
  }
  for (const auto& w : sc.warnings()) std::printf("[警告] %s\n", w.c_str());
  if (!info.empty()) std::printf("%s\n", info.c_str());

  m.assignDofs();
  const yjk::Model::Diagnostic d = m.check();
  for (const auto& w : d.warnings) std::printf("[警告] %s\n", w.c_str());
  if (!d.ok()) {
    for (const auto& e : d.errors) std::printf("[校核] %s\n", e.c_str());
    return 2;
  }

  const std::string stem = stemOverride.empty() ? sc.outStem() : stemOverride;
  yjk::post::ExportOptions opt;

  // 没有 case/modal/spectrum/combo 命令 → 旧版单工况路径（行为不变）
  if (!sc.hasCases() && !sc.hasDynamic()) {
    yjk::StaticAnalysis sa(m);
    const yjk::StaticResult r = sa.solve();
    if (!r.ok) {
      std::printf("[求解失败] %s\n", r.message.c_str());
      return 3;
    }
    yjk::post::PostProcessor pp(m, r);
    pp.compute();
    std::printf("%s", pp.report().c_str());
    const auto f = yjk::post::exportAll(stem, m, r, pp, opt);
    if (!f.ok) {
      std::printf("[导出失败] %s\n", f.error.c_str());
      return 3;
    }
    std::printf("\n已写出:\n  %s\n  %s\n  %s\n",
                f.jsonPath.c_str(), f.htmlPath.c_str(), f.txtPath.c_str());
    return 0;
  }

  // ================== 动态分析流程 ==================
  std::string modalDesc = sc.hasModal()
      ? (sc.modalNmodes() > 0 ? std::to_string(sc.modalNmodes()) + " 阶"
                              : std::string("自动(max(9,3×层))"))
      : std::string("无");
  std::printf("[动态分析] 工况 %zu 个，模态 %s，反应谱 %s，组合 %zu 条\n",
              sc.caseNames().size(), modalDesc.c_str(),
              sc.hasSpectrum() ? "有" : "无", sc.combos().size());

  // ① 静力工况：无 case 命令时缺省工况 "D"（荷载已由 build 施加）
  std::vector<std::string> order =
      sc.hasCases() ? sc.caseNames() : std::vector<std::string>{"D"};
  std::map<std::string, yjk::StaticResult> caseResults;
  yjk::StaticAnalysis sa(m);
  std::string caseInfo;
  // gmass 反算出的各层重力代表值（story → 该层及以上合计 W，剪重比分母）
  std::vector<std::pair<int, double>> gravityByStory;
  for (const auto& name : order) {
    caseInfo.clear();
    if (sc.hasCases()) {
      m.clearLoads();
      if (!sc.applyCase(m, name, &caseInfo)) {
        for (const auto& e : sc.errors()) std::printf("[工况错误] %s\n", e.c_str());
        return 1;
      }
      if (!caseInfo.empty()) std::printf("    %s\n", caseInfo.c_str());
    }
    yjk::StaticResult r = sa.solve();
    if (!r.ok) {
      std::printf("[求解失败] 工况 %s: %s\n", name.c_str(), r.message.c_str());
      return 3;
    }
    caseResults[name] = r;
    std::printf("  工况 %-3s 最大位移 %-10.6g m  Uxy %-9.5g  反力 %-9.4g kN\n",
                name.c_str(), r.maxDisplacement, r.maxUxy, r.maxReaction);
  }

  // ② gmass：重力荷载代表值转质量（GB 50011 5.1.3）
  //    节点重量 = 恒载工况竖向节点力 + ψ·活载工况竖向节点力。
  //    恒载工况的 selfweight 把自重作为线荷载进入节点等效荷载，
  //    故折算出的节点重量已含自重 —— ModalAnalysis 检测到节点重量
  //    后改用"节点重量代表团"路径，不再叠加单元自重（避免重复计重）。
  yjk::ModalResult mr;
  if (sc.hasGmass()) {
    std::vector<double> w(static_cast<size_t>(m.nodeCount()), 0.0);
    auto collect = [&](const std::string& cname, double fac) {
      m.clearLoads();
      caseInfo.clear();
      if (!sc.applyCase(m, cname, &caseInfo)) {
        for (const auto& e : sc.errors()) std::printf("[工况错误] %s\n", e.c_str());
        return false;
      }
      if (!caseInfo.empty()) std::printf("    %s\n", caseInfo.c_str());
      yjk::StaticAnalysis::Assembly as;
      std::string why;
      if (!sa.assemble(as, &why)) {
        std::printf("[gmass失败] 工况 %s 组装失败：%s\n", cname.c_str(), why.c_str());
        return false;
      }
      // f 按【节点序 × 6】索引，压缩 [2] 是竖向分量（z 轴向上，向下为负）
      for (yjk::Id i = 0; i < m.nodeCount(); ++i) {
        const double fz = as.f[static_cast<size_t>(i) * 6 + 2];
        if (fz < 0.0) w[static_cast<size_t>(i)] += fac * (-fz);
      }
      return true;
    };
    bool gmassOk = collect(sc.gmassD(), 1.0) && collect(sc.gmassL(), sc.gmassPsi());
    m.clearLoads();
    if (gmassOk) {
      double wSum = 0.0;
      int cnt = 0;
      for (yjk::Id i = 0; i < m.nodeCount(); ++i)
        if (w[static_cast<size_t>(i)] > 0.0) {
          m.addNodeWeight(i, w[static_cast<size_t>(i)]);
          wSum += w[static_cast<size_t>(i)];
          ++cnt;
        }
      std::printf("  gmass: 恒载(%s) + %.2f·活载(%s) → 节点重量 %d 个，合计 %.3f kN\n",
                  sc.gmassD().c_str(), sc.gmassPsi(), sc.gmassL().c_str(), cnt, wSum);

      // 重力代表值工况 G：把节点重量作为竖向节点荷载静力求解，
      // 供各工况/组合的剪重比作分母（剪重比 V/W，W = 该层以上重力代表值）。
      // G 不进入 order/caseResults（不参与组合与包络），仅取各层重力 W。
      m.clearLoads();
      for (yjk::Id i = 0; i < m.nodeCount(); ++i)
        if (w[static_cast<size_t>(i)] > 0.0)
          m.addNodeForce(i, yjk::Vec3{0.0, 0.0, -w[static_cast<size_t>(i)]});
      {
        const yjk::StaticResult g = sa.solve();
        if (g.ok) {
          // 必须在 clearLoads 之前从模型取各层重力（PostProcessor 读模型荷载）
          yjk::post::PostProcessor pg(m, g);
          pg.compute();
          for (const auto& q : pg.stories())
            gravityByStory.push_back({q.story, q.weight});
          std::printf("  重力代表值 W 合计 %.3f kN（剪重比基准）\n", wSum);
        } else {
          std::printf("[警告] 重力代表值工况 G 求解失败：%s\n", g.message.c_str());
        }
      }
      m.clearLoads();
    }
  }

  // ② 模态分析
  if (sc.hasModal() || sc.hasGmass()) {
    yjk::ModalAnalysis ma(m);
    yjk::ModalAnalysis::Options mo;
    mo.nmodes = sc.modalNmodes();       // 0 = 自动（max(9, 层数×3)）
    mr = ma.solve(mo);
    if (mr.ok) {
      std::printf("  模态: %d 阶完成（ndof %d，迭代 %d，残差 %g）\n",
                  mr.nmodes, mr.ndof, mr.iter, mr.residual);
      std::printf("    参与质量 X %.1f%%  Y %.1f%%  Z %.1f%%（总质量 %.3f t）\n",
                  mr.ratioX * 100.0, mr.ratioY * 100.0, mr.ratioZ * 100.0,
                  mr.massX);
    } else {
      std::printf("[模态失败] %s\n", mr.message.c_str());
    }
  }

  // ③ 反应谱 → 等效地震力工况（Ex/Ey/Ez，独立静力求解）
  yjk::SpectrumResult sr;
  if (sc.hasSpectrum() && mr.ok) {
    yjk::DesignSpectrum spec;
    spec.alphaMax = sc.spectrumAlphaMax();
    spec.tg = sc.spectrumTg();
    yjk::ResponseSpectrum::Options so;
    so.method = sc.spectrumMethod();
    so.verticalScale = sc.spectrumVerticalScale();
    yjk::ResponseSpectrum rs(m);
    sr = rs.solve(mr, spec, so);
    if (sr.ok) {
      std::printf("  反应谱: 基底剪力 X=%-9.4g Y=%-9.4g Z=%-9.4g kN\n",
                  sr.baseX, sr.baseY, sr.baseZ);
      std::printf("    参与质量 X %.1f%%%s  Y %.1f%%%s（≥90%% 检查 %s）\n",
                  sr.massRatioX * 100.0, sr.massOkX ? "✓" : "✗",
                  sr.massRatioY * 100.0, sr.massOkY ? "✓" : "✗",
                  sr.massOk ? "通过" : "未通过");
    } else {
      std::printf("[反应谱失败] %s\n", sr.message.c_str());
    }
    // 等效地震力 → 节点荷载 → 静力求解，作为独立工况（组合引用 Ex/Ey）
    const auto applySpectrum = [&](const std::vector<double>& f,
                                   const char* name) {
      m.clearLoads();
      for (yjk::Id i = 0; i < m.nodeCount(); ++i) {
        yjk::Vec3 v{0, 0, 0};
        bool any = false;
        for (int k = 0; k < 3; ++k) {
          if (m.node(i).dof[k] < 0) continue;   // 固定/从属自由度质量已折入主自由度
          const double x = f[static_cast<size_t>(i) * 6 + k];
          if (x != 0.0) { v[k] = x; any = true; }
        }
        if (any) m.addNodeForce(i, v);
      }
      const yjk::StaticResult r = sa.solve();
      if (!r.ok) {
        std::printf("[求解失败] 工况 %s: %s\n", name, r.message.c_str());
        return;
      }
      order.push_back(name);
      caseResults[name] = r;
      std::printf("  工况 %-3s 最大位移 %-10.6g m\n", name, r.maxDisplacement);
    };
    applySpectrum(sr.forceX, "Ex");
    applySpectrum(sr.forceY, "Ey");
    applySpectrum(sr.forceZ, "Ez");
  }

  // ④ 荷载组合 + 包络
  std::vector<yjk::ComboDef> defs;
  for (const auto& cr : sc.combos()) defs.push_back(buildComboDef(cr));
  std::vector<yjk::StaticResult> combos;
  if (!defs.empty()) {
    combos = yjk::LoadCombination::generate(caseResults, defs);
    if (combos.size() != defs.size())
      std::printf("[组合警告] 生成 %zu 条 / 定义 %zu 条（引用缺失工况者被忽略）\n",
                  combos.size(), defs.size());
    for (size_t i = 0; i < combos.size(); ++i) {
      std::printf("  组合 %-12s 最大位移 %-10.6g m  Uxy %-9.5g\n",
                  defs[i].name.c_str(), combos[i].maxDisplacement,
                  combos[i].maxUxy);
    }
  }
  yjk::StaticResult envMax, envMin;
  if (!combos.empty()) yjk::LoadCombination::envelope(combos, envMax, envMin);

  // ⑤ 导出：主工况三件套 + 动态 json/txt
  const yjk::StaticResult& mainR = caseResults.at(order.front());
  yjk::post::PostProcessor pp(m, mainR);
  pp.compute();
  std::printf("%s", pp.report().c_str());
  const auto f = yjk::post::exportAll(stem, m, mainR, pp, opt);
  if (!f.ok) {
    std::printf("[导出失败] %s\n", f.error.c_str());
    return 3;
  }

  yjk::post::DynamicExport dyn;
  if (mr.ok) dyn.modal = &mr;
  if (sr.ok) dyn.spectrum = &sr;
  dyn.gravityByStory = gravityByStory;
  for (const auto& name : order)
    dyn.cases.push_back({name, &caseResults.at(name)});
  for (size_t i = 0; i < combos.size(); ++i)
    dyn.combos.push_back({defs[i].name, &combos[i]});
  if (!combos.empty()) {
    dyn.envMax = &envMax;
    dyn.envMin = &envMin;
  }
  const auto fd = yjk::post::exportDynamic(stem, m, dyn, opt);
  if (!fd.ok) {
    std::printf("[导出失败] %s\n", fd.error.c_str());
    return 3;
  }

  std::printf("\n已写出:\n  %s\n  %s\n  %s\n  %s\n  %s\n",
              f.jsonPath.c_str(), f.htmlPath.c_str(), f.txtPath.c_str(),
              fd.jsonPath.c_str(), fd.txtPath.c_str());
  return 0;
}

int cmdCheck(const std::string& path) {
  yjk::io::ModelScript sc;
  if (!sc.parseFile(path)) {
    for (const auto& e : sc.errors()) std::printf("[脚本错误] %s\n", e.c_str());
    return 1;
  }
  yjk::Model m;
  std::string info;
  if (!sc.build(m, &info)) {
    for (const auto& e : sc.errors()) std::printf("[建模错误] %s\n", e.c_str());
    return 1;
  }
  for (const auto& w : sc.warnings()) std::printf("[警告] %s\n", w.c_str());
  std::printf("%s\n", info.c_str());
  m.assignDofs();
  const yjk::Model::Diagnostic d = m.check();
  for (const auto& w : d.warnings) std::printf("[警告] %s\n", w.c_str());
  if (!d.ok()) {
    for (const auto& e : d.errors) std::printf("[校核] %s\n", e.c_str());
    return 2;
  }
  std::printf("校核通过：节点 %d，单元 %d，自由度 %d\n",
              m.nodeCount(), m.elementCount(), m.freeDofCount());
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) { usage(); return 4; }
  const std::string cmd = argv[1];
  if (cmd == "help" || cmd == "-h" || cmd == "--help") { usage(); return 0; }
  if (cmd != "run" && cmd != "check") { usage(); return 4; }

  std::string path;
  std::string stem;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-o" && i + 1 < argc) { stem = argv[++i]; continue; }
    if (!a.empty() && a[0] == '-') { std::printf("未知选项: %s\n", a.c_str()); return 4; }
    if (path.empty()) path = a;
    else { std::printf("多余的参数: %s\n", a.c_str()); return 4; }
  }
  if (path.empty()) { std::printf("缺少模型文件路径\n"); return 4; }

  return (cmd == "run") ? cmdRun(path, stem) : cmdCheck(path);
}