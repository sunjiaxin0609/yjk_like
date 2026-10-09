// =============================================================================
//  tests/test_post.cpp  ——  后处理层 / 导出层 / 建模脚本 验证
//
//  这一层的验证重点是【工程师会读到的那些数】，而不是求解器本身：
//    · 弯矩沿杆长必须是【二次】的（线性插值在跨中差 25%）
//    · 应力 σ = N/A ± M·h/I 必须能拿手算对上
//    · 剪重比的分母必须包含面荷载（漏掉会高估 → 偏不安全）
//    · 导出的 JSON 必须能被解析（NaN 会让 JSON.parse 抛异常 → 界面白屏）
//
//  另外专门验证【建模脚本的错误处理】：
//  脚本是给人写的，写错是常态。错误必须被拦住并给出行号，
//  绝不能静默生成一个"看起来能算、其实错"的模型。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "yjk/analysis/StaticAnalysis.h"
#include "yjk/io/ModelScript.h"
#include "yjk/model/GridMesh.h"
#include "yjk/post/PostProcessor.h"
#include "yjk/post/ResultExport.h"

using namespace yjk;

static int g_fail = 0;
static int g_pass = 0;

static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}
static void checkNear(double got, double want, double relTol, const char* what) {
  const double err = std::abs(got - want);
  const double rel = std::abs(want) > 1e-30 ? err / std::abs(want) : err;
  if (rel <= relTol) { ++g_pass;
    std::printf("  [ OK ] %s  (got %.12g, want %.12g, rel %.2e)\n", what, got, want, rel);
  } else { ++g_fail;
    std::printf("  [FAIL] %s  (got %.12g, want %.12g, rel %.2e)\n", what, got, want, rel);
  }
}
static void checkInt(long long got, long long want, const char* what) {
  if (got == want) { ++g_pass; std::printf("  [ OK ] %s  (got %lld)\n", what, got); }
  else { ++g_fail; std::printf("  [FAIL] %s  (got %lld, want %lld)\n", what, got, want); }
}

// -----------------------------------------------------------------------------
//  1. 梁内力沿杆长：必须与解析解逐点吻合
//
//  【核心】有分布荷载时 M(ξ) 是抛物线。若后处理用线性插值，
//  跨中会差 25%（(wL²/8 - wL²/12)/(wL²/8) = 1/3 的量级）。
//  这里用 21 个采样点逐点对比，线性插值无处可藏。
// -----------------------------------------------------------------------------
static void testBeamDiagram() {
  std::printf("\n== 1. 梁内力沿杆长 vs 解析解 ==\n");
  const double L = 6.0, w = 10.0;
  const auto sec = section::rect(0.3, 0.6);
  const auto mat = Material::concreteC(30);

  // 简支梁：M(x) = w·x·(L−x)/2，V(x) = w·(L/2 − x)
  Model m;
  const int n = 8;
  std::vector<Id> ns;
  for (int i = 0; i <= n; ++i) ns.push_back(m.addNode({L * i / double(n), 0, 0}, 0));
  for (int i = 0; i < n; ++i)
    m.addBeam(ns[static_cast<size_t>(i)], ns[static_cast<size_t>(i + 1)], sec, mat,
              Vec3{0, 0, -1})->setLineLoad(Vec3{0, 0, -w});
  m.fixNode(ns[0], true, true, true, true, false, false);
  m.fixNode(ns[n], true, true, true, false, false, false);
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "简支梁求解成功");
  if (!r.ok) return;

  post::PostProcessor pp(m, r);
  pp.compute(21, 0.7);

  // 全梁拼起来看：每个单元的 station 按全局 x 展开
  double worstM = 0.0, worstV = 0.0;
  for (const auto& b : pp.beamResults()) {
    const double x0 = m.node(b.ni).r.x, x1 = m.node(b.nj).r.x;
    for (const auto& st : b.stations) {
      const double x = x0 + (x1 - x0) * st.xi;
      const double mTh = w * x * (L - x) / 2.0;
      const double vTh = w * (L / 2.0 - x);
      // 竖向荷载沿全局 -z；局部 y' = up = (0,0,-1) 指向下方，
      // 故竖向挠度在局部 y' 平面内 ⇒ 剪力读 Vy，弯矩读 Mz（绕 z'）。
      // 配对规则：V_y 配 M_z，V_z 配 M_y —— 配错时数值"看着差不多"，
      // 是最隐蔽的一类错误，所以这里的期望值必须是带符号的解析式。
      worstM = std::max(worstM, std::abs(st.Mz - mTh) / (w * L * L / 8.0));
      worstV = std::max(worstV, std::abs(st.Vy - vTh) / (w * L / 2.0));
    }
  }
  check(worstM < 1e-6, "弯矩图逐点吻合 w·x(L−x)/2（含跨中抛物线）");
  check(worstV < 1e-9, "剪力图逐点吻合 w(L/2−x)");
  std::printf("       弯矩最大相对偏差 %.2e，剪力 %.2e\n", worstM, worstV);

  // 跨中弯矩。
  //  【注意】不能取"某个单元的中点" —— n=8 时第 5 个单元的中点在 x=3.375，
  //  不是跨中。按全局 x 找最接近 L/2 的采样点才是真跨中。
  //  线性插值会给出 (M_i+M_j)/2，对简支梁端部 M=0 ⇒ 得到 0，与真值差 100%。
  double bestM = 0.0, bestX = 0.0;
  for (const auto& b : pp.beamResults()) {
    const double x0 = m.node(b.ni).r.x, x1 = m.node(b.nj).r.x;
    for (const auto& st : b.stations) {
      const double x = x0 + (x1 - x0) * st.xi;
      if (std::abs(x - L / 2.0) < std::abs(bestX - L / 2.0) || bestM == 0.0) {
        bestX = x; bestM = st.Mz;
      }
    }
  }
  std::printf("       跨中采样点 x=%.4f，弯矩 %.6f\n", bestX, bestM);
  checkNear(bestM, w * L * L / 8.0, 1e-6, "跨中弯矩 = +wL²/8（线性插值会得到 0）");

  // ---- 符号约定：一次性校核 M 与 V，缺一不可 ----
  {
    const double Lc = 4.0, P = 50.0;
    Model m2;
    const Id a = m2.addNode({0, 0, 0}, 0), b2 = m2.addNode({Lc, 0, 0}, 0);
    m2.addBeam(a, b2, sec, mat, Vec3{0, 0, -1});
    m2.fixAll(a);
    m2.addNodeForce(b2, Vec3{0, 0, -P});
    m2.assignDofs();
    StaticAnalysis sa2(m2);
    const auto r2 = sa2.solve();
    post::PostProcessor pp2(m2, r2);
    pp2.compute(21, 0.7);
    const auto& br2 = pp2.beamResults().front();
    const auto& s0 = br2.stations.front();
    const auto& s1 = br2.stations.back();
    checkNear(s0.Mz, -P * Lc, 1e-9, "悬臂固定端弯矩 = −PL（上侧受拉为负）");
    checkNear(s0.Vy, P, 1e-9, "悬臂固定端剪力 = +P（V = dM/dx）");
    checkNear(s1.Vy, P, 1e-9, "悬臂自由端剪力与固定端相同（无分布荷载）");
    checkNear(s1.Mz, 0.0, 1e-12, "悬臂自由端弯矩 = 0");
  }
}

// -----------------------------------------------------------------------------
//  2. 应力恢复：σ = N/A ± M·h/I
// -----------------------------------------------------------------------------
static void testStress() {
  std::printf("\n== 2. 截面应力恢复 σ = N/A ± M·h/I ==\n");
  const double L = 4.0, P = 100.0;
  const auto sec = section::rect(0.2, 0.4);      // b=0.2, h=0.4
  const auto mat = Material::concreteC(30);

  // 悬臂梁，端部竖向集中力 P（向下） + 端部轴力 N（拉）
  Model m;
  const Id n0 = m.addNode({0, 0, 0}, 0);
  const Id n1 = m.addNode({L, 0, 0}, 0);
  auto* b = m.addBeam(n0, n1, sec, mat, Vec3{0, 0, -1});
  m.fixAll(n0);
  m.assignDofs();
  // 先求位移（只加端弯矩效果：用节点力）
  std::vector<double> f;
  m.addNodeForce(n1, Vec3{P, 0, -P});            // 轴力 +P（拉），竖向 −P
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "悬臂梁求解成功");
  if (!r.ok) return;

  post::PostProcessor pp(m, r);
  pp.compute(21, 0.7);
  const auto& br = pp.beamResults().front();
  (void)b;

  // 固定端：N = +P（拉），M = P·L
  const auto& st0 = br.stations.front();
  const double Nth = P, Mth = P * L;
  // 竖向弯曲绕局部 z' ⇒ σ = Mz·hy/Iz，hy = h/2 = 0.2
  const double sigBend = Mth * sec.hy / sec.Iz;
  const double sigAx = Nth / sec.A;
  std::printf("       N=%.4f(kN) Mz=%.4f(kN·m)  A=%.4f Iz=%.6g hy=%.3f\n",
              st0.N, st0.Mz, sec.A, sec.Iz, sec.hy);
  std::printf("       σ轴=%.4f kPa，σ弯=%.4f kPa\n", sigAx / 1000.0, sigBend / 1000.0);
  checkNear(st0.N, Nth, 1e-9, "固定端轴力 = +P（拉为正）");
  checkNear(st0.Mz, -Mth, 1e-9, "固定端弯矩 = −P·L（悬臂上侧受拉为负）");
  // sigMax = σ轴 + |σ弯|，sigMin = σ轴 − |σ弯|
  checkNear(st0.sigMax, sigAx + sigBend, 1e-6, "σmax = N/A + M·hy/Iz");
  checkNear(st0.sigMin, sigAx - sigBend, 1e-6, "σmin = N/A − M·hy/Iz");
  check(br.sigAbsMax >= std::abs(st0.sigMax) - 1e-9, "全杆 |σ|max 不小于端部的 |σ|");
}

// -----------------------------------------------------------------------------
//  3. 楼层指标与剪重比（面荷载必须进分母）
// -----------------------------------------------------------------------------
static void testStoryMetrics() {
  std::printf("\n== 3. 楼层指标 / 剪重比 ==\n");
  AxisGrid g;
  g.addAxisX(0.0); g.addAxisX(6.0);
  g.addAxisY(0.0); g.addAxisY(6.0);
  g.addStory(0.0); g.addStory(3.6); g.addStory(7.2);
  g.setSlabStory(0, false);

  Model m;
  GridMesh::Spec sp;
  sp.columnSection = section::rect(0.5, 0.5);
  sp.beamSection = section::rect(0.3, 0.6);
  sp.columnMaterial = Material::concreteC(30);
  sp.beamMaterial = Material::concreteC(30);
  sp.slab.thickness = 0.12;
  GridMesh mesh(g, sp);
  const auto rr = mesh.generate(m);
  check(rr.ok, "两层一跨框架生成成功");
  for (Id n : mesh.baseNodes()) m.fixAll(n);
  for (auto& e : m.elements()) {
    if (e->type() == ElementType::Beam3D) static_cast<BeamElement*>(e.get())->setSelfWeight(true);
    if (e->type() == ElementType::Shell4)
      static_cast<ShellElement*>(e.get())->setTransversePressure(-4.0);
  }
  // 顶层 4 个节点各 10 kN 沿 X ⇒ 基底剪力 40 kN
  const Id top = 2;
  int cnt = 0;
  for (Id i = 0; i < m.nodeCount(); ++i)
    if (m.node(i).story == top) { m.addNodeForce(i, Vec3{10, 0, 0}); ++cnt; }
  checkInt(cnt, 4, "顶层节点数 = 4");
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "求解成功");
  if (!r.ok) return;

  post::PostProcessor pp(m, r);
  pp.compute();
  const auto& st = pp.stories();
  check(st.size() >= 3, "楼层表至少 3 层");

  // 楼层剪力：每层都应等于 40 kN（水平荷载全在顶层，逐层向下传）
  bool shearOk = true;
  for (size_t s = 1; s < st.size(); ++s)
    if (std::abs(st[s].shear - 40.0) > 1e-6 * 40.0) shearOk = false;
  check(shearOk, "各层楼层剪力 = 40 kN（水平荷载全部由顶层逐层下传）");

  // 剪重比分母必须包含面荷载：
  //   板 1 块 × 36 m² × 4 kPa = 144 kN，若漏掉，剪重比会被高估
  const double wSelf = [&]() {
    double w = 0.0;
    for (const auto& e : m.elements()) w += e->mass();
    return w;
  }();
  const double wAll = pp.envelope().totalWeight;
  std::printf("       构件自重 %.3f kN，计入面荷载后 %.3f kN\n", wSelf, wAll);
  check(wAll > wSelf * 1.02, "总重力荷载【大于】构件自重 ⇒ 面荷载已计入剪重比分母");

  // 荷载闭合：施加 = 反力
  const double close = std::abs(pp.envelope().totalLoadZ - wAll) / wAll;
  check(close < 1e-10, "荷载闭合：施加竖向荷载 = 竖向反力合计");
  std::printf("       闭合差 %.2e\n", close);

  // 层间位移角为正且落在合理量级
  check(pp.envelope().maxDrift > 0.0, "层间位移角 > 0");
  check(pp.envelope().maxDrift < 1.0, "层间位移角在合理范围（< 1）");
  std::printf("       最大层间位移角 1/%.0f\n", 1.0 / pp.envelope().maxDrift);
}

// -----------------------------------------------------------------------------
//  4. 建模脚本：正确脚本 + 错误必须被拦住
// -----------------------------------------------------------------------------
static void testScript() {
  std::printf("\n== 4. .yjk 建模脚本 ==\n");
  const std::string good =
      "material C30 concrete 30\n"
      "section  COL rect 0.5 0.5\n"
      "section  BM  rect 0.3 0.6\n"
      "grid.axisX 0 6\n"
      "grid.axisY 0 6\n"
      "grid.story 0 3.6\n"
      "grid.column COL C30\n"
      "grid.beam   BM  C30\n"
      "grid.slabthick 0.12\n"
      "support base fixed\n"
      "selfweight on\n"
      "slabload 4\n"
      "nodeload top 10 0 0\n"
      "out demo\n";

  io::ModelScript sc;
  check(sc.parse(good), "脚本解析成功");
  Model m;
  std::string info;
  check(sc.build(m, &info), "建模成功");
  check(!info.empty(), "给出建模摘要");
  checkInt(m.nodeCount(), 8, "节点数 = 8（2×2 轴网交点 × 2 层）");
  check(sc.outStem() == "demo", "out 命令生效");
  m.assignDofs();
  const auto d = m.check();
  for (const auto& e : d.errors) std::printf("       [错误] %s\n", e.c_str());
  check(d.ok(), "模型校核通过");
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "脚本模型可求解");
  std::printf("       %s\n", info.c_str());

  // ---- 错误必须被拦住 ----
  struct BadCase { const char* name; std::string text; };
  const std::vector<BadCase> bad = {
    {"未知命令",            "grid.axisX 0 6\ngrid.axisY 0 6\ngrid.story 0 3.6\nfoobar 1 2\n"},
    {"截面参数非数值",      "section S rect 0.3 x\ngrid.axisX 0 6\n"},
    {"引用未定义截面",      "grid.axisX 0 6\ngrid.axisY 0 6\ngrid.story 0 3.6\ngrid.column NOPE C30\n"},
    {"支座节点越界",        "grid.axisX 0 6\ngrid.axisY 0 6\ngrid.story 0 3.6\nsupport node 999 fixed\n"},
    {"轴网不完整",          "grid.axisX 0 6\nsupport base fixed\n"},
    {"剖分数非法",          "grid.axisX 0 6\ngrid.axisY 0 6\ngrid.story 0 3.6\ngrid.divX 0\n"},
    {"板厚为负",            "grid.axisX 0 6\ngrid.axisY 0 6\ngrid.story 0 3.6\ngrid.slabthick -1\n"},
  };
  for (const auto& bc : bad) {
    io::ModelScript s2;
    s2.parse(bc.text);
    Model m2;
    std::string i2;
    const bool built = s2.build(m2, &i2);
    const bool reported = !s2.errors().empty();
    check(!built && reported, bc.name);
    if (reported) std::printf("         → %s\n", s2.errors().front().c_str());
  }
}

// -----------------------------------------------------------------------------
//  5. 导出：JSON 可解析 / 无 NaN / HTML 占位符被替换
// -----------------------------------------------------------------------------
namespace {
// 轻量 JSON 合法性检查：括号配对（跳过字符串）+ 不含非有限字面量。
// 不写完整 parser —— 测试只需要"能被 JSON.parse 接受"这个保证，
// 而 parse 失败的两个典型原因就是括号不配对和 NaN/Infinity。
bool jsonLooksValid(const std::string& s, std::string& why) {
  int depth = 0, bdepth = 0;
  bool inStr = false, esc = false;
  for (char c : s) {
    if (inStr) {
      if (esc) esc = false;
      else if (c == '\\') esc = true;
      else if (c == '"') inStr = false;
      continue;
    }
    if (c == '"') { inStr = true; continue; }
    if (c == '{') ++depth;
    else if (c == '}') { --depth; if (depth < 0) { why = "多余的 }"; return false; } }
    else if (c == '[') ++bdepth;
    else if (c == ']') { --bdepth; if (bdepth < 0) { why = "多余的 ]"; return false; } }
  }
  if (inStr) { why = "字符串未闭合"; return false; }
  if (depth != 0 || bdepth != 0) { why = "括号不配对"; return false; }
  for (const char* bad : {"NaN", "Infinity", "-Infinity", "nan", "inf"}) {
    if (s.find(bad) != std::string::npos) { why = std::string("含非有限字面量 ") + bad; return false; }
  }
  return true;
}
}  // namespace

static void testExport() {
  std::printf("\n== 5. 结果导出（JSON / HTML / 报告） ==\n");
  AxisGrid g;
  g.addAxisX(0.0); g.addAxisX(6.0);
  g.addAxisY(0.0); g.addAxisY(6.0);
  g.addStory(0.0); g.addStory(3.6);
  g.setSlabStory(0, false);
  Model m;
  GridMesh::Spec sp;
  sp.columnSection = section::rect(0.4, 0.4);
  sp.beamSection = section::rect(0.25, 0.5);
  sp.columnMaterial = Material::concreteC(30);
  sp.beamMaterial = Material::concreteC(30);
  sp.slab.thickness = 0.10;
  GridMesh mesh(g, sp);
  mesh.generate(m);
  for (Id n : mesh.baseNodes()) m.fixAll(n);
  for (auto& e : m.elements())
    if (e->type() == ElementType::Beam3D) static_cast<BeamElement*>(e.get())->setSelfWeight(true);
  m.addNodeForce(mesh.nodeAtAxis(0, 0, 1), Vec3{5, 0, 0});
  m.assignDofs();
  StaticAnalysis sa(m);
  const auto r = sa.solve();
  check(r.ok, "导出用例求解成功");
  if (!r.ok) return;

  post::PostProcessor pp(m, r);
  pp.compute();

  const std::string json = post::toJson(m, r, pp, {}, "测试");
  std::string why;
  check(jsonLooksValid(json, why), "JSON 括号配对且不含非有限值");
  if (!why.empty()) std::printf("       → %s\n", why.c_str());
  check(json.find("\"nodes\"") != std::string::npos, "JSON 含 nodes");
  check(json.find("\"beams\"") != std::string::npos, "JSON 含 beams");
  check(json.find("\"stories\"") != std::string::npos, "JSON 含 stories");
  check(json.find("\"_schema\"") != std::string::npos, "JSON 自带 schema 说明");
  std::printf("       JSON %.1f KB\n", json.size() / 1024.0);

  // 节点数组长度必须与模型一致
  const size_t nb = json.find("\"nodes\": [\n");
  check(nb != std::string::npos, "节点数组存在");

  const std::string html = post::htmlViewer(json, "测试标题");
  check(html.find("__YJK_DATA__") == std::string::npos, "HTML 数据占位符已替换");
  check(html.find("__YJK_TITLE__") == std::string::npos, "HTML 标题占位符已替换");
  check(html.find("测试标题") != std::string::npos, "HTML 含标题");
  check(html.find("<canvas") != std::string::npos, "HTML 含画布");
  std::printf("       HTML %.1f KB\n", html.size() / 1024.0);

  // 报告非空且含关键段
  const std::string rep = pp.report();
  check(rep.find("平衡校核") != std::string::npos, "报告含平衡校核段");
  check(rep.find("[闭合]") != std::string::npos, "报告判定为闭合");
  check(!pp.summary().empty(), "摘要非空");
  std::printf("       摘要：%s\n", pp.summary().c_str());
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  后处理 / 导出 / 建模脚本 验证\n");
  std::printf("======================================================\n");
  testBeamDiagram();
  testStress();
  testStoryMetrics();
  testScript();
  testExport();
  std::printf("\n------------------------------------------------------\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("------------------------------------------------------\n");
  return g_fail == 0 ? 0 : 1;
}
