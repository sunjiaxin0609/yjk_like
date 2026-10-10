// =============================================================================
//  src/io/ModelScript.cpp  ——  .yjk 建模脚本解析与建模
//
//  【错误处理原则：收集而不是抛出】
//  建模脚本是给人写的，一次写对的概率不高。
//  遇到错误就抛异常/中断，用户只能"改一个、跑一次"；
//  收集全部错误一次报出来，可以一次改完。
//  所以所有解析与建模函数都返回 bool + 往 errors_ 里追加，不抛异常。
//
//  【数值解析必须校验】
//  strtod 对 "abc" 返回 0 且不报错 —— 不校验的话，
//  用户写成 "rect 0.3 x" 会得到 b=0.3、h=0 的退化截面，
//  一路算下去得到"结果很奇怪但没报错"的模型。
//  这类错误最难发现，必须在这里拦住。
// =============================================================================
#include "yjk/io/ModelScript.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

#include "yjk/model/GridMesh.h"

namespace yjk {
namespace io {
namespace {

std::vector<std::string> splitWs(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream is(s);
  std::string t;
  while (is >> t) out.push_back(t);
  return out;
}

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

// 严格数值解析：必须整串消费完，且结果有限
bool toNum(const std::string& s, double& v) {
  if (s.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const double x = std::strtod(s.c_str(), &end);
  if (end != s.c_str() + s.size()) return false;
  if (errno == ERANGE || !std::isfinite(x)) return false;
  v = x;
  return true;
}

bool toInt(const std::string& s, int& v) {
  double d = 0;
  if (!toNum(s, d)) return false;
  if (std::abs(d - std::floor(d)) > 1e-12) return false;
  v = static_cast<int>(d);
  return true;
}

}  // namespace

// -----------------------------------------------------------------------------
bool ModelScript::parse(const std::string& text) {
  cmds_.clear();
  errors_.clear();
  warnings_.clear();
  caseLoads_.clear();
  caseOrder_.clear();
  currentCase_ = "D";
  sawCaseCmd_ = false;

  int lineNo = 0;
  std::istringstream in(text);
  std::string raw;
  while (std::getline(in, raw)) {
    ++lineNo;
    if (!raw.empty() && raw.back() == '\r') raw.pop_back();

    // 去注释
    const size_t h = raw.find('#');
    if (h != std::string::npos) raw = raw.substr(0, h);

    const std::vector<std::string> tk = splitWs(raw);
    if (tk.empty()) continue;

    Command c;
    c.name = lower(tk[0]);
    c.args.assign(tk.begin() + 1, tk.end());
    c.line = lineNo;

    // ---- 工况与动态分析请求：解析期收集，不参与几何建模 ----
    if (c.name == "case") {
      if (c.args.empty()) errors_.push_back("第 " + std::to_string(lineNo) +
                                             " 行：case 需要一个工况名");
      else {
        currentCase_ = c.args[0];
        sawCaseCmd_ = true;
        if (std::find(caseOrder_.begin(), caseOrder_.end(), currentCase_) ==
            caseOrder_.end())
          caseOrder_.push_back(currentCase_);
      }
      continue;
    }
    if (c.name == "modal") {
      // modal            —— 自动（max(9, 层数×3)，GB 50011 每方向 ≥9）
      // modal auto | 0   —— 同上
      // modal <n>        —— 显式取前 n 阶（n ≥ 1）
      hasModalCmd_ = true;
      if (c.args.empty()) {
        modalNmodes_ = 0;
      } else if (lower(c.args[0]) == "auto" || c.args[0] == "0") {
        modalNmodes_ = 0;
      } else {
        int n = 0;
        if (!toInt(c.args[0], n) || n < 1)
          errors_.push_back("第 " + std::to_string(lineNo) +
                            " 行：modal 参数应为正整数、0 或 auto（0/auto = 自动取值）");
        else modalNmodes_ = n;
      }
      continue;
    }
    if (c.name == "spectrum") {
      if (c.args.size() < 2 || !toNum(c.args[0], alphaMax_) ||
          !toNum(c.args[1], tg_) || alphaMax_ <= 0.0 || tg_ <= 0.0) {
        errors_.push_back("第 " + std::to_string(lineNo) +
                          " 行：spectrum 需要 <alphaMax> <Tg> [srss|cqc]"
                          "（正数）");
        hasSpectrum_ = false;
      } else {
        hasSpectrum_ = true;
        method_ = 0;
        if (c.args.size() >= 3 && lower(c.args[2]) == "cqc") method_ = 1;
        else if (c.args.size() >= 3 && lower(c.args[2]) != "srss")
          warnings_.push_back("第 " + std::to_string(lineNo) + " 行：谱组合方式 '" +
                              c.args[2] + "' 未知（srss/cqc），按 SRSS 处理");
        // 第 4 参：竖向地震影响系数放大系数 k（竖向谱取水平谱 × k，
        // GB 50011 5.3.1 竖向地震作用取水平地震作用 0.65），默认 0.65
        if (c.args.size() >= 4) {
          double k = 0.0;
          if (!toNum(c.args[3], k) || k <= 0.0)
            errors_.push_back("第 " + std::to_string(lineNo) +
                              " 行：spectrum 第 4 参（竖向系数 k）应为正数，"
                              "得到 '" + c.args[3] + "'");
          else verticalScale_ = k;
        }
      }
      continue;
    }
    if (c.name == "gmass") {
      // gmass <恒载工况> <活载工况> [ψ=0.5]
      // GB 50011 5.1.3 重力荷载代表值转质量：恒载 + ψ·活载。
      if (c.args.size() < 2) {
        errors_.push_back("第 " + std::to_string(lineNo) +
                          " 行：gmass 需要 <恒载工况> <活载工况> [ψ=0.5]");
        continue;
      }
      double psi = 0.5;
      if (c.args.size() >= 3 && !toNum(c.args[2], psi)) {
        errors_.push_back("第 " + std::to_string(lineNo) + " 行：组合值系数 '" +
                          c.args[2] + "' 不是合法数值");
        continue;
      }
      if (!(psi > 0.0) || psi > 1.0) {
        errors_.push_back("第 " + std::to_string(lineNo) +
                          " 行：组合值系数必须在 (0, 1] 区间（规范取 0.5 等）");
        continue;
      }
      gmassD_ = c.args[0];
      gmassL_ = c.args[1];
      gmassPsi_ = psi;
      hasGmass_ = true;
      continue;
    }
    if (c.name == "combo") {
      // combo <standard|basicvar|basicperm|seismic> <组合名> <工况名...>
      if (c.args.size() < 3) {
        errors_.push_back("第 " + std::to_string(lineNo) +
                          " 行：combo 需要 <kind> <组合名> <工况名...>（至少一个工况）");
        continue;
      }
      ComboRequest r;
      r.kind = lower(c.args[0]);
      r.name = c.args[1];
      r.cases.assign(c.args.begin() + 2, c.args.end());
      if (r.kind != "standard" && r.kind != "basicvar" &&
          r.kind != "basicperm" && r.kind != "seismic") {
        errors_.push_back("第 " + std::to_string(lineNo) + " 行：组合类型 '" +
                          r.kind + "' 未知（standard/basicvar/basicperm/seismic）");
        continue;
      }
      combos_.push_back(std::move(r));
      continue;
    }

    if (c.name == "out") {
      if (c.args.empty()) errors_.push_back("第 " + std::to_string(lineNo) +
                                            " 行：out 需要一个文件名前缀");
      else outStem_ = c.args[0];
    }

    // 多工况模式下，荷载命令按当前工况分组存储。
    // 【为什么放在 parse 阶段】case 命令只在此处可见：它被 consume 不进
    // cmds_，build() 遍历 cmds_ 时 currentCase_ 保持 parse 结束时的值，
    // 若在 build 阶段分组会把所有荷载误归到最后声明的工况。
    if (sawCaseCmd_ && (c.name == "selfweight" || c.name == "slabload" ||
                        c.name == "beamload" || c.name == "nodeload" ||
                        c.name == "nodemoment")) {
      caseLoads_[currentCase_].push_back(c);
      continue;
    }

    cmds_.push_back(std::move(c));
  }

  // gmass 必须与 case 命令配套：恒/活荷载分属不同工况才能折算。
  // （没有 case 命令时所有荷载同属缺省工况 "D"，无法区分恒载与活载。）
  if (hasGmass_ && !sawCaseCmd_) {
    errors_.push_back("gmass 需要 case 命令：恒载与活载必须分属不同工况（如 case D / case L）");
  } else if (hasGmass_) {
    auto defined = [&](const std::string& n) {
      return std::find(caseOrder_.begin(), caseOrder_.end(), n) != caseOrder_.end();
    };
    if (!defined(gmassD_))
      errors_.push_back("gmass 引用的恒载工况 '" + gmassD_ + "' 未定义（先用 case 声明）");
    if (!defined(gmassL_))
      errors_.push_back("gmass 引用的活载工况 '" + gmassL_ + "' 未定义（先用 case 声明）");
  }
  return errors_.empty();
}

bool ModelScript::parseFile(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) {
    errors_.push_back("无法打开脚本文件：" + path);
    return false;
  }
  std::ostringstream ss;
  ss << f.rdbuf();
  return parse(ss.str());
}

std::string ModelScript::describe() const {
  std::ostringstream os;
  for (const Command& c : cmds_) {
    os << "  " << c.line << ": " << c.name;
    for (const auto& a : c.args) os << " " << a;
    os << "\n";
  }
  return os.str();
}

// -----------------------------------------------------------------------------
//  建模
//
//  【为什么用 map 存命名截面/材料，而不是立刻建成对象】
//  脚本里 section 定义可以出现在 grid.column 之后（人写脚本不会care顺序）。
//  若边解析边建，后定义的截面就引用不到。
//  所以先把命令全部收下，build 时按顺序执行；截面/材料表在执行中逐步填充，
//  引用发生在"用到它的那一刻"，天然支持任意书写顺序。
// -----------------------------------------------------------------------------
bool ModelScript::build(Model& m, std::string* info) {
  errors_.clear();
  warnings_.clear();

  std::map<std::string, Material> mats;
  std::map<std::string, SectionProperties> secs;

  // 默认值：不写 material/section 也能跑
  mats["C30"] = Material::concreteC(30);
  secs["B"] = section::rect(0.3, 0.6);
  secs["C"] = section::rect(0.5, 0.5);

  // 重置成员缓存状态（脚本可能被重复 build —— 例如 CLI 多工况流程）
  grid_ = AxisGrid();
  meshSpec_ = GridMesh::Spec();
  meshSpec_.columnSection = secs["C"];
  meshSpec_.beamSection = secs["B"];
  meshSpec_.columnMaterial = mats["C30"];
  meshSpec_.beamMaterial = mats["C30"];
  hasGrid_ = false;
  mesh_.reset();
  meshBuilt_ = false;
  slabOpt_ = GridMesh::SlabLoadOptions();
  slabOpt_.path = GridMesh::SlabLoadOptions::Path::Auto;

  // 面荷载的去处与单向板判定（slabload.to / slabload.ratio 可改）
  // 网格必须在使用节点查询前生成（applyLoadCmd / support / diaphragm 共用）
  auto ensureMesh = [&]() -> bool { return ensureMeshCached(m, info); };

  auto findMat = [&](const std::string& n, Material& out, int line) -> bool {
    auto it = mats.find(n);
    if (it == mats.end()) {
      errors_.push_back("第 " + std::to_string(line) + " 行：未定义的材料 '" + n + "'");
      return false;
    }
    out = it->second;
    return true;
  };
  auto findSec = [&](const std::string& n, SectionProperties& out, int line) -> bool {
    auto it = secs.find(n);
    if (it == secs.end()) {
      errors_.push_back("第 " + std::to_string(line) + " 行：未定义的截面 '" + n + "'");
      return false;
    }
    out = it->second;
    return true;
  };
  auto need = [&](const Command& c, size_t n) -> bool {
    if (c.args.size() < n) {
      errors_.push_back("第 " + std::to_string(c.line) + " 行：" + c.name +
                        " 需要至少 " + std::to_string(n) + " 个参数，实际 " +
                        std::to_string(c.args.size()) + " 个");
      return false;
    }
    return true;
  };

  for (const Command& c : cmds_) {
    const auto& a = c.args;

    // ---------------- 材料 ----------------
    if (c.name == "material") {
      if (!need(c, 3)) continue;
      const std::string kind = lower(a[1]);
      if (kind == "raw") {
        // material <名> raw <E> <nu> <gamma> <fy> <ft> <fc> <G>
        // 原始材料：任意字段直接给；G ≤ 0 表示由 E、ν 自动算（与
        // Material 默认一致）。交互层写入的任意材料都走这条命令。
        if (!need(c, 9)) continue;
        double E = 0, nu = 0, gamma = 0, fy = 0, ft = 0, fc = 0, G = 0;
        bool ok = true;
        const char* names[7] = {"E", "ν", "γ", "fy", "ft", "fc", "G"};
        double* dst[7] = {&E, &nu, &gamma, &fy, &ft, &fc, &G};
        for (int i = 0; i < 7; ++i) {
          if (!toNum(a[static_cast<size_t>(i + 2)], *dst[i])) {
            errors_.push_back("第 " + std::to_string(c.line) + " 行：材料参数 " +
                              names[i] + " '" + a[static_cast<size_t>(i + 2)] +
                              "' 不是合法数值");
            ok = false;
          }
        }
        if (!ok) continue;
        if (!(E > 0.0) || !(nu >= 0.0 && nu < 0.5) || !(gamma >= 0.0)) {
          errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：材料参数不合法（需 E>0，0≤ν<0.5，γ≥0）");
          continue;
        }
        Material rawMat;
        rawMat.E = E; rawMat.nu = nu; rawMat.gamma = gamma;
        rawMat.fy = fy; rawMat.ft = ft; rawMat.fc = fc; rawMat.G = G;
        mats[a[0]] = rawMat;
        continue;
      }
      int grade = 0;
      if (!toInt(a[2], grade)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：牌号 '" + a[2] + "' 不是整数");
        continue;
      }
      if (kind == "concrete" || kind == "c") mats[a[0]] = Material::concreteC(grade);
      else if (kind == "steel" || kind == "q" || kind == "s") mats[a[0]] = Material::steelQ(grade);
      else errors_.push_back("第 " + std::to_string(c.line) + " 行：材料类型 '" + a[1] +
                             "' 未知（应为 concrete / steel / raw）");
      continue;
    }

    // ---------------- 截面 ----------------
    if (c.name == "section") {
      if (!need(c, 2)) continue;
      const std::string kind = lower(a[1]);
      std::vector<double> p;
      bool ok = true;
      for (size_t i = 2; i < a.size(); ++i) {
        double v = 0;
        if (!toNum(a[i], v)) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：截面参数 '" + a[i] +
                            "' 不是合法数值");
          ok = false;
          break;
        }
        p.push_back(v);
      }
      if (!ok) continue;
      if (kind == "rect" && p.size() == 2) secs[a[0]] = section::rect(p[0], p[1]);
      else if (kind == "circle" && p.size() == 1) secs[a[0]] = section::circle(p[0]);
      else if (kind == "tube" && p.size() == 2) secs[a[0]] = section::tube(p[0], p[1]);
      else if (kind == "i" && p.size() == 4) secs[a[0]] = section::iSection(p[0], p[1], p[2], p[3]);
      else if (kind == "box" && p.size() == 4) secs[a[0]] = section::boxSection(p[0], p[1], p[2], p[3]);
      else if (kind == "raw" && p.size() == 10) {
        // section <名> raw <A> <Iy> <Iz> <J> <Asy> <Asz> <Ry> <Rz> <hy> <hz>
        // 原始截面：任意几何直接给（交互层读取-修改-写回任意截面都用这条）。
        SectionProperties s;
        s.A = p[0]; s.Iy = p[1]; s.Iz = p[2]; s.J = p[3];
        s.Asy = p[4]; s.Asz = p[5]; s.Ry = p[6]; s.Rz = p[7];
        s.hy = p[8]; s.hz = p[9];
        std::string why;
        if (!s.isValid(&why)) errors_.push_back("第 " + std::to_string(c.line) +
                                                " 行：截面 '" + a[0] + "' 无效：" + why);
        else secs[a[0]] = s;
      }
      else errors_.push_back("第 " + std::to_string(c.line) + " 行：截面 '" + a[0] + "' 的型式 '" +
                             a[1] + "' 与参数个数（" + std::to_string(p.size()) +
                             "）不匹配（rect 2 / circle 1 / tube 2 / i 4 / box 4 / raw 10）");
      continue;
    }

    // ---------------- 原始构件（P3：交互层与脚本层共用同一套内核接口） ----------------
    // 不经轴网，直接从显式节点/单元建模型。node 必须先行声明（按出现顺序
    // 分配 id=0,1,2,…），beam/column/shell/wall/fix/nodeweight 引用这些 id。
    // YjkWriter 序列化的脚本（另存为 .yjk）就只包含这一组命令，不写轴网。
    if (c.name == "node") {
      if (!need(c, 3)) continue;
      double x = 0, y = 0, z = 0;
      int story = 0;
      if (!toNum(a[0], x) || !toNum(a[1], y) || !toNum(a[2], z)) {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：node 坐标不是三个合法数值");
        continue;
      }
      if (a.size() >= 4 && !toInt(a[3], story)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：层号 '" + a[3] +
                          "' 不是整数");
        continue;
      }
      m.addNode(Vec3{x, y, z}, story);
      continue;
    }
    if (c.name == "beam" || c.name == "column") {
      // beam <i> <j> <截面> <材料> [upx upy upz]
      // column 与 beam 同构：柱 = 竖向梁，只是无 up 时默认 {0,-1,0}。
      if (!need(c, 4)) continue;
      Id i = 0, j = 0;
      int ni = 0, nj = 0;
      const bool okI = toInt(a[0], ni) && ni >= 0;
      const bool okJ = toInt(a[1], nj) && nj >= 0;
      if (!okI || !okJ) {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：单元端点号必须是 ≥0 的整数");
        continue;
      }
      if (ni >= m.nodeCount() || nj >= m.nodeCount()) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：单元端点号 (" +
                          std::to_string(ni) + "," + std::to_string(nj) +
                          ") 越界（当前节点 0.." + std::to_string(m.nodeCount() - 1) + "）");
        continue;
      }
      i = static_cast<Id>(ni); j = static_cast<Id>(nj);
      SectionProperties sec;
      Material mat;
      if (!findSec(a[2], sec, c.line) || !findMat(a[3], mat, c.line)) continue;
      Vec3 up = (c.name == "column") ? Vec3{0, -1, 0} : Vec3{0, 0, -1};
      if (a.size() >= 7) {
        double ux = 0, uy = 0, uz = 0;
        if (!toNum(a[4], ux) || !toNum(a[5], uy) || !toNum(a[6], uz)) {
          errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：up 向量不是三个合法数值");
          continue;
        }
        up = Vec3{ux, uy, uz};
      }
      BeamElement* b = (c.name == "column")
                           ? m.addColumn(i, j, sec, mat, up)
                           : m.addBeam(i, j, sec, mat, up);
      if (!b) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：" + c.name +
                          " 创建失败（端点越界或长度为零）");
        continue;
      }
      continue;
    }
    if (c.name == "shell" || c.name == "wall") {
      // shell <a> <b> <c> <d> <厚> [E nu density]   墙 = shell + isWall
      if (!need(c, 5)) continue;
      std::vector<Id> ns;
      bool nodeOk = true;
      for (size_t k = 0; k < 4; ++k) {
        int v = 0;
        if (!toInt(a[k], v) || v < 0 || v >= m.nodeCount()) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号 '" + a[k] +
                            "' 越界（当前节点 0.." + std::to_string(m.nodeCount() - 1) + "）");
          nodeOk = false;
          break;
        }
        ns.push_back(static_cast<Id>(v));
      }
      if (!nodeOk) continue;
      double t = 0;
      if (!toNum(a[4], t) || !(t > 0.0)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：厚度必须为正数");
        continue;
      }
      ShellProperties props;
      props.thickness = t;
      if (a.size() >= 6) {
        double E = 0;
        if (!toNum(a[5], E) || !(E > 0.0)) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：E 必须为正数");
          continue;
        }
        props.E = E;
      }
      if (a.size() >= 7) {
        double nu = 0;
        if (!toNum(a[6], nu) || !(nu >= 0.0 && nu < 0.5)) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：ν 必须在 [0,0.5)");
          continue;
        }
        props.nu = nu;
      }
      if (a.size() >= 8) {
        double d = 0;
        if (!toNum(a[7], d) || !(d > 0.0)) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：density 必须为正数");
          continue;
        }
        props.density = d;
      }
      if (c.name == "wall") m.addWall(ns, props);
      else m.addShell(ns, props);
      continue;
    }
    if (c.name == "fix") {
      // fix node <id> <ux> <uy> <uz> <rx> <ry> <rz>（0/1，任意组合约束）
      if (!need(c, 8)) continue;
      if (lower(a[0]) != "node") {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：fix 的目标只能是 node（当前只支持显式节点）");
        continue;
      }
      int id = 0;
      if (!toInt(a[1], id) || id < 0 || id >= m.nodeCount()) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号 " + a[1] +
                          " 越界（有效范围 0.." + std::to_string(m.nodeCount() - 1) + "）");
        continue;
      }
      std::array<bool, 6> fx{};
      bool ok = true;
      for (int k = 0; k < 6; ++k) {
        const std::string s = lower(a[static_cast<size_t>(k + 2)]);
        if (s == "1" || s == "on" || s == "true") fx[static_cast<size_t>(k)] = true;
        else if (s == "0" || s == "off" || s == "false") fx[static_cast<size_t>(k)] = false;
        else {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：约束分量 '" +
                            a[static_cast<size_t>(k + 2)] + "' 应为 0/1");
          ok = false;
        }
      }
      if (!ok) continue;
      m.fixNode(static_cast<Id>(id), fx[0], fx[1], fx[2], fx[3], fx[4], fx[5]);
      continue;
    }
    if (c.name == "nodeweight") {
      // nodeweight node <id> <w> —— 节点附加重量（建模属性，不分工况）
      if (!need(c, 3)) continue;
      if (lower(a[0]) != "node") {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：nodeweight 的目标只能是 node");
        continue;
      }
      int id = 0;
      if (!toInt(a[1], id) || id < 0 || id >= m.nodeCount()) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号 " + a[1] +
                          " 越界（有效范围 0.." + std::to_string(m.nodeCount() - 1) + "）");
        continue;
      }
      double w = 0;
      if (!toNum(a[2], w) || !(w > 0.0)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：重量必须为正数");
        continue;
      }
      m.addNodeWeight(static_cast<Id>(id), w);
      continue;
    }

    // ---------------- 单元级荷载 / 释放（T6 全量序列化） ----------------
    // 序号 = 该类单元在 elements() 中的计数（与 `release beam <序号>` 的
    // seen 计数口径一致：只数对应类型，跳过其他单元）。必须出现在
    // 单元命令之后（序列化文本天然满足：单元块 → 荷载块 → 约束块）。
    if (c.name == "beamsw") {
      // beamsw <序号> <0|1> —— 该梁的自重开关
      if (!need(c, 2)) continue;
      int id = 0;
      if (!toInt(a[0], id) || id < 0) { errors_.push_back("第 " + std::to_string(c.line) +
                                        " 行：梁序号 '" + a[0] + "' 不是非负整数"); continue; }
      const std::string v = lower(a[1]);
      if (v != "0" && v != "1" && v != "on" && v != "off" && v != "true" && v != "false") {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：beamsw 参数应为 0/1"); continue;
      }
      const bool on = (v == "1" || v == "on" || v == "true");
      int seen = 0;
      bool hit = false;
      for (auto& e : m.elements()) {
        if (e->type() != ElementType::Beam3D) continue;
        if (seen++ != id) continue;
        static_cast<BeamElement*>(e.get())->setSelfWeight(on);
        hit = true;
        break;
      }
      if (!hit) errors_.push_back("第 " + std::to_string(c.line) +
                                  " 行：梁序号 " + a[0] + " 越界（梁单元共 " +
                                  std::to_string(seen) + " 根）");
      continue;
    }
    if (c.name == "beamseg") {
      // beamseg <序号> <x1> <x2> <q1x q1y q1z q2x q2y q2z>
      if (!need(c, 9)) continue;
      int id = 0;
      if (!toInt(a[0], id) || id < 0) { errors_.push_back("第 " + std::to_string(c.line) +
                                        " 行：梁序号 '" + a[0] + "' 不是非负整数"); continue; }
      double v[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      bool ok = true;
      for (int i = 0; i < 8; ++i) {
        if (!toNum(a[static_cast<size_t>(i + 1)], v[static_cast<size_t>(i)])) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：beamseg 参数 '" +
                            a[static_cast<size_t>(i + 1)] + "' 不是合法数值");
          ok = false;
        }
      }
      if (!ok) continue;
      int seen = 0;
      bool hit = false;
      for (auto& e : m.elements()) {
        if (e->type() != ElementType::Beam3D) continue;
        if (seen++ != id) continue;
        static_cast<BeamElement*>(e.get())->addLineLoadSegment(
            v[0], v[1], Vec3{v[2], v[3], v[4]}, Vec3{v[5], v[6], v[7]});
        hit = true;
        break;
      }
      if (!hit) errors_.push_back("第 " + std::to_string(c.line) +
                                  " 行：梁序号 " + a[0] + " 越界（梁单元共 " +
                                  std::to_string(seen) + " 根）");
      continue;
    }
    if (c.name == "beampoint") {
      // beampoint <序号> <12 个等效节点荷载值>
      if (!need(c, 13)) continue;
      int id = 0;
      if (!toInt(a[0], id) || id < 0) { errors_.push_back("第 " + std::to_string(c.line) +
                                        " 行：梁序号 '" + a[0] + "' 不是非负整数"); continue; }
      std::vector<double> pt;
      pt.reserve(12);
      bool ok = true;
      for (int i = 0; i < 12; ++i) {
        double x = 0;
        if (!toNum(a[static_cast<size_t>(i + 1)], x)) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：beampoint 分量 '" +
                            a[static_cast<size_t>(i + 1)] + "' 不是合法数值");
          ok = false;
        }
        pt.push_back(x);
      }
      if (!ok) continue;
      int seen = 0;
      bool hit = false;
      for (auto& e : m.elements()) {
        if (e->type() != ElementType::Beam3D) continue;
        if (seen++ != id) continue;
        static_cast<BeamElement*>(e.get())->setPointLoads(pt);
        hit = true;
        break;
      }
      if (!hit) errors_.push_back("第 " + std::to_string(c.line) +
                                  " 行：梁序号 " + a[0] + " 越界（梁单元共 " +
                                  std::to_string(seen) + " 根）");
      continue;
    }
    if (c.name == "shellsw") {
      // shellsw <序号> <0|1> —— 该壳（含墙）的自重开关
      if (!need(c, 2)) continue;
      int id = 0;
      if (!toInt(a[0], id) || id < 0) { errors_.push_back("第 " + std::to_string(c.line) +
                                        " 行：壳序号 '" + a[0] + "' 不是非负整数"); continue; }
      const std::string v = lower(a[1]);
      if (v != "0" && v != "1" && v != "on" && v != "off" && v != "true" && v != "false") {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：shellsw 参数应为 0/1"); continue;
      }
      const bool on = (v == "1" || v == "on" || v == "true");
      int seen = 0;
      bool hit = false;
      for (auto& e : m.elements()) {
        if (e->type() != ElementType::Shell4) continue;
        if (seen++ != id) continue;
        static_cast<ShellElement*>(e.get())->setSelfWeight(on);
        hit = true;
        break;
      }
      if (!hit) errors_.push_back("第 " + std::to_string(c.line) +
                                  " 行：壳序号 " + a[0] + " 越界（壳单元共 " +
                                  std::to_string(seen) + " 个）");
      continue;
    }
    if (c.name == "shellp") {
      // shellp <序号> <pz> [px py] —— 板面压（kPa，单元内存储值）与膜压
      // 注意：pz 直接写入单元（不做工程正负号转换）—— 序列化输出的是
      // 单元内存储的原值，与 YjkWriter 输出互逆。
      if (!need(c, 2)) continue;
      int id = 0;
      if (!toInt(a[0], id) || id < 0) { errors_.push_back("第 " + std::to_string(c.line) +
                                        " 行：壳序号 '" + a[0] + "' 不是非负整数"); continue; }
      double pz = 0, px = 0, py = 0;
      if (!toNum(a[1], pz)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：面压 '" + a[1] +
                          "' 不是合法数值"); continue;
      }
      if (a.size() >= 3 && !toNum(a[2], px)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：膜压 px '" + a[2] +
                          "' 不是合法数值"); continue;
      }
      if (a.size() >= 4 && !toNum(a[3], py)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：膜压 py '" + a[3] +
                          "' 不是合法数值"); continue;
      }
      int seen = 0;
      bool hit = false;
      for (auto& e : m.elements()) {
        if (e->type() != ElementType::Shell4) continue;
        if (seen++ != id) continue;
        ShellElement* s = static_cast<ShellElement*>(e.get());
        s->setTransversePressure(pz);
        if (a.size() >= 3) s->setMembranePressure(px, py);
        hit = true;
        break;
      }
      if (!hit) errors_.push_back("第 " + std::to_string(c.line) +
                                  " 行：壳序号 " + a[0] + " 越界（壳单元共 " +
                                  std::to_string(seen) + " 个）");
      continue;
    }

    // ---------------- 轴网 ----------------
    if (c.name == "grid.axisx") {
      for (const auto& s : a) {
        double v = 0;
        if (!toNum(s, v)) { errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：'" + s + "' 不是合法坐标"); continue; }
        grid_.addAxisX(v);
      }
      hasGrid_ = true;
      continue;
    }
    if (c.name == "grid.axisy") {
      for (const auto& s : a) {
        double v = 0;
        if (!toNum(s, v)) { errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：'" + s + "' 不是合法坐标"); continue; }
        grid_.addAxisY(v);
      }
      hasGrid_ = true;
      continue;
    }
    if (c.name == "grid.story") {
      for (const auto& s : a) {
        double v = 0;
        if (!toNum(s, v)) { errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：'" + s + "' 不是合法标高"); continue; }
        grid_.addStory(v);
      }
      hasGrid_ = true;
      continue;
    }
    if (c.name == "grid.divx" || c.name == "grid.divy") {
      for (size_t i = 0; i < a.size(); ++i) {
        int n = 1;
        if (!toInt(a[i], n) || n < 1) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：剖分数 '" + a[i] +
                            "' 必须是 ≥1 的整数");
          continue;
        }
        const bool okk = (c.name == "grid.divx") ? grid_.setDivX(static_cast<Id>(i), n)
                                                 : grid_.setDivY(static_cast<Id>(i), n);
        if (!okk) warnings_.push_back("第 " + std::to_string(c.line) + " 行：第 " +
                                      std::to_string(i) + " 个剖分数越界（跨度数不足），已忽略");
      }
      continue;
    }
    if (c.name == "grid.column") {
      if (!need(c, 2)) continue;
      findSec(a[0], meshSpec_.columnSection, c.line);
      findMat(a[1], meshSpec_.columnMaterial, c.line);
      continue;
    }
    if (c.name == "grid.beam") {
      if (!need(c, 2)) continue;
      findSec(a[0], meshSpec_.beamSection, c.line);
      findMat(a[1], meshSpec_.beamMaterial, c.line);
      continue;
    }
    if (c.name == "grid.slabthick") {
      if (!need(c, 1)) continue;
      double t = 0.12;
      if (!toNum(a[0], t) || !(t > 0.0)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：板厚必须为正数");
        continue;
      }
      meshSpec_.slab.thickness = t;
      continue;
    }
    if (c.name == "grid.slab") {
      // grid.slab on|off —— 是否生成楼板【壳单元】。
      //
      //  【与 grid.slaboff 的区别】两者看起来都像"关掉板"，语义完全不同：
      //    grid.slab off    —— 这一层【不建板单元】：荷载改走导荷到梁，
      //                        面内刚度靠刚性楼板假定。板还是有的。
      //    grid.slaboff 1 2 —— 这几层【没有楼面】：既不建板也不导荷。
      //  把它俩混起来会得到"楼面荷载凭空消失"的模型 —— 而模型照样能解。
      if (!need(c, 1)) continue;
      const std::string v = lower(a[0]);
      if (v == "on" || v == "1" || v == "true") meshSpec_.buildSlabs = true;
      else if (v == "off" || v == "0" || v == "false") meshSpec_.buildSlabs = false;
      else { errors_.push_back("第 " + std::to_string(c.line) +
                               " 行：grid.slab 的参数应为 on / off"); continue; }
      if (mesh_) {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：grid.slab 必须写在任何需要生成网格的命令之前"
                          "（它决定网格怎么生成）");
      }
      continue;
    }
    if (c.name == "grid.slaboff") {
      for (const auto& s : a) {
        int k = 0;
        if (!toInt(s, k)) { errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：层号 '" + s + "' 不是整数"); continue; }
        grid_.setSlabStory(static_cast<Id>(k), false);
      }
      continue;
    }
    if (c.name == "grid.wall") {
      // grid.wall on|off —— 是否生成墙【壳单元】。
      //
      //  【与 AxisGrid 墙开关的配合】AxisGrid 默认无墙（不设即无墙），
      //  所以 grid.wall on 本身不会凭空生成墙；必须配合 grid.wallx /
      //  grid.wally 先在轴网上开墙段。grid.wall off 用于整模型关墙。
      if (!need(c, 1)) continue;
      const std::string v = lower(a[0]);
      if (v == "on" || v == "1" || v == "true") meshSpec_.buildWalls = true;
      else if (v == "off" || v == "0" || v == "false") meshSpec_.buildWalls = false;
      else { errors_.push_back("第 " + std::to_string(c.line) +
                               " 行：grid.wall 的参数应为 on / off"); continue; }
      if (mesh_) {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：grid.wall 必须写在任何需要生成网格的命令之前"
                          "（它决定网格怎么生成）");
      }
      continue;
    }
    if (c.name == "grid.wallthick") {
      if (!need(c, 1)) continue;
      double t = 0.2;
      if (!toNum(a[0], t) || !(t > 0.0)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：墙厚必须为正数");
        continue;
      }
      meshSpec_.wall.thickness = t;
      continue;
    }
    if (c.name == "grid.wallx" || c.name == "grid.wally" ||
        c.name == "grid.nowallx" || c.name == "grid.nowally") {
      // 墙段开关（索引一律 0 基，越界只警告不报错，与 nobeamx 一致）：
      //   grid.wallx  <ix> <jy>   在 X 向轴线 ix、Y 向跨度 jy 加墙
      //   grid.wally  <jy> <ix>   在 Y 向轴线 jy、X 向跨度 ix 加墙
      //   grid.nowallx <ix> <jy>  移除该墙段
      //   grid.nowally <jy> <ix>  移除该墙段
      if (!need(c, 2)) continue;
      int i0 = 0, i1 = 0;
      if (!toInt(a[0], i0) || !toInt(a[1], i1)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：索引必须是整数");
        continue;
      }
      bool okk = false;
      if (c.name == "grid.wallx") okk = grid_.setWallAlongX(i0, i1, true);
      else if (c.name == "grid.wally") okk = grid_.setWallAlongY(i0, i1, true);
      else if (c.name == "grid.nowallx") okk = grid_.setWallAlongX(i0, i1, false);
      else okk = grid_.setWallAlongY(i0, i1, false);
      if (!okk) warnings_.push_back("第 " + std::to_string(c.line) + " 行：索引 (" +
                                    std::to_string(i0) + "," + std::to_string(i1) +
                                    ") 越界，已忽略");
      continue;
    }
    if (c.name == "diaphragm") {
      // diaphragm on [起始层]   —— 加刚性楼板假定（默认从第 1 层起）
      // diaphragm off           —— 不加（默认）
      //
      //  【为什么默认从第 1 层起】第 0 层是嵌固层/基础顶，通常整体固定，
      //  给它加刚性假定是多余的；强行加上还会把柱脚节点绑成一个刚体，
      //  掩盖"某个柱脚漏了约束"这类建模错误。
      if (!need(c, 1)) continue;
      const std::string v = lower(a[0]);
      if (v == "off" || v == "0" || v == "false") continue;
      if (v != "on" && v != "1" && v != "true") {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：diaphragm 的参数应为 on / off（on 后可跟起始层号）");
        continue;
      }
      if (!ensureMesh()) continue;
      int from = 1;
      if (a.size() >= 2 && !toInt(a[1], from)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：起始层号 '" + a[1] +
                          "' 不是整数");
        continue;
      }
      const Id made = mesh_->addRigidDiaphragms(m, static_cast<Id>(std::max(0, from)), true);
      if (info) *info += "；刚性楼板 " + std::to_string(made) + " 层";
      continue;
    }
    if (c.name == "diaphragm.bind") {
      // diaphragm.bind <story> <coupleRz 0/1> <masterId> <slaveIds...>
      // 把 slave 节点绑到【已存在】的主节点上（T6 全量序列化读回路径）。
      // 与 `diaphragm on` 的区别：不新建主节点、不依赖网格 ——
      // 主节点在序列化文本里就是普通 node（带固定约束），由本命令
      // 把它升级为主节点并重建 DofLink（系数与源模型逐位一致）。
      if (!need(c, 4)) continue;
      int story = 0;
      if (!toInt(a[0], story)) { errors_.push_back("第 " + std::to_string(c.line) +
                                 " 行：层号 '" + a[0] + "' 不是整数"); continue; }
      const std::string cr = lower(a[1]);
      if (cr != "0" && cr != "1" && cr != "on" && cr != "off" &&
          cr != "true" && cr != "false") {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：coupleRz 应为 0/1"); continue;
      }
      const bool coupleRz = (cr == "1" || cr == "on" || cr == "true");
      int mid = 0;
      if (!toInt(a[2], mid) || mid < 0 || mid >= m.nodeCount()) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：主节点号 '" + a[2] +
                          "' 越界（当前节点 0.." + std::to_string(m.nodeCount() - 1) + "）");
        continue;
      }
      std::vector<Id> slaves;
      bool ok = true;
      for (size_t i = 3; i < a.size(); ++i) {
        int v = 0;
        if (!toInt(a[i], v) || v < 0 || v >= m.nodeCount()) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：从属节点号 '" + a[i] +
                            "' 越界（当前节点 0.." + std::to_string(m.nodeCount() - 1) + "）");
          ok = false;
        } else {
          slaves.push_back(static_cast<Id>(v));
        }
      }
      if (!ok) continue;
      const Id made = m.attachRigidDiaphragm(static_cast<Id>(mid), slaves, story, coupleRz);
      if (made == kDofFixed) errors_.push_back("第 " + std::to_string(c.line) +
                                  " 行：diaphragm.bind 主节点 " + a[2] + " 无效");
      else if (info) *info += "；刚性楼板绑定 " + std::to_string(slaves.size()) + " 节点";
      continue;
    }
    if (c.name == "grid.nocolumn" || c.name == "grid.nobeamx" ||
        c.name == "grid.nobeamy" || c.name == "grid.noslab") {
      if (!need(c, 2)) continue;
      int i0 = 0, i1 = 0;
      if (!toInt(a[0], i0) || !toInt(a[1], i1)) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：索引必须是整数");
        continue;
      }
      bool okk = false;
      if (c.name == "grid.nocolumn") okk = grid_.setColumn(i0, i1, false);
      else if (c.name == "grid.nobeamx") okk = grid_.setBeamAlongX(i0, i1, false);
      else if (c.name == "grid.nobeamy") okk = grid_.setBeamAlongY(i0, i1, false);
      else okk = grid_.setSlab(i0, i1, false);
      if (!okk) warnings_.push_back("第 " + std::to_string(c.line) + " 行：索引 (" +
                                    std::to_string(i0) + "," + std::to_string(i1) +
                                    ") 越界，已忽略");
      continue;
    }

    // ---------------- 支座 ----------------
    if (c.name == "support") {
      if (!need(c, 2)) continue;
      //
      //  参数顺序统一为  support <目标> <目标号> <类型>
      //    support base  fixed          （柱脚没有"号"）
      //    support node  5   pinned
      //    support story 2   fixed
      //  类型永远是最后一个参数 —— 目标有无编号都保持同一读法，
      //  避免"有时 a[1] 是类型、有时 a[1] 是编号"这种最易写错的接口。
      const std::string kind = lower(a.back());
      auto apply = [&](Id n) {
        if (kind == "fixed") m.fixAll(n);
        else if (kind == "pinned") m.restrainPinned(n);
        else if (kind == "roller") m.restrainRoller(n);
        else errors_.push_back("第 " + std::to_string(c.line) + " 行：支座类型 '" + a[1] +
                               "' 未知（fixed / pinned / roller）");
      };
      const std::string tgt = lower(a[0]);
      if (tgt == "base") {
        if (!ensureMesh()) continue;
        const auto bases = mesh_->baseNodes();
        if (bases.empty()) warnings_.push_back("第 " + std::to_string(c.line) +
                                               " 行：没有找到柱脚节点（是否关掉了全部柱？）");
        for (Id n : bases) apply(n);
      } else if (tgt == "node") {
        if (!need(c, 3)) continue;
        int id = 0;
        if (!toInt(a[1], id) || id < 0 || id >= m.nodeCount()) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号 " + a[1] +
                            " 越界（有效范围 0.." + std::to_string(m.nodeCount() - 1) + "）");
          continue;
        }
        apply(static_cast<Id>(id));
      } else if (tgt == "story") {
        if (!need(c, 3)) continue;
        if (!ensureMesh()) continue;
        int s = 0;
        if (!toInt(a[1], s)) { errors_.push_back("第 " + std::to_string(c.line) +
                               " 行：层号不是整数"); continue; }
        int cnt = 0;
        for (Id i = 0; i < m.nodeCount(); ++i) {
          // 刚性楼板主节点是计算过程的产物，不该被用户施加支座
          if (m.node(i).diaphragmMaster) continue;
          if (m.node(i).story == s) { apply(i); ++cnt; }
        }
        if (cnt == 0) warnings_.push_back("第 " + std::to_string(c.line) + " 行：第 " +
                                          std::to_string(s) + " 层没有节点");
      } else {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：支座目标 '" + a[0] +
                          "' 未知（base / node / story）");
      }
      continue;
    }

    // ---------------- 弹簧单元 ----------------
    //  spring <节点i> <节点j> <kx> <ky> <kz> <krx> <kry> <krz> [upx upy upz]
    //  两节点之间 6 个方向的独立弹簧（局部 ux,uy,uz,rx,ry,rz）。
    //    k=0      → 该方向释放（配合 DofNumbering 零刚度自动锁定）
    //    k 很大   → 近似刚接（比"端部释放"更强的连接自由度）
    //  刚度单位：平动 kN/m，转动 kN·m/rad。
    //  可选的 up 向量定局部 y' 方向（缺省 {0,0,-1}，与梁单元相同）。
    //  注意：弹簧必须引用已存在的节点（addNode 或网格生成后）。
    if (c.name == "spring") {
      if (!need(c, 8)) continue;
      // 节点号引用的是"已经存在的节点"（网格生成或库 API 添加）。
      // 脚本里节点由轴网生成 —— 若尚未生成网格则先生成，否则节点号必然越界。
      if (hasGrid_ && !mesh_ && !ensureMesh()) continue;
      auto nodeIndex = [&](const std::string& s, Id& out) -> bool {
        int v = 0;
        if (!toInt(s, v) || v < 0 || v >= m.nodeCount()) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号 " + s +
                            " 越界（有效范围 0.." + std::to_string(m.nodeCount() - 1) + "）");
          return false;
        }
        out = static_cast<Id>(v);
        return true;
      };
      Id ni = 0, nj = 0;
      if (!nodeIndex(a[0], ni) || !nodeIndex(a[1], nj)) continue;
      std::array<double, 6> ks{};
      bool kOk = true;
      for (int i = 0; i < 6; ++i) {
        if (!toNum(a[static_cast<size_t>(i + 2)], ks[static_cast<size_t>(i)]) ||
            ks[static_cast<size_t>(i)] < 0.0) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：弹簧刚度 '" +
                            a[static_cast<size_t>(i + 2)] + "' 必须是 ≥0 的数");
          kOk = false;
        }
      }
      if (!kOk) continue;
      Vec3 up{0, 0, -1};
      if (a.size() >= 11) {
        double ux = 0, uy = 0, uz = 0;
        if (!toNum(a[8], ux) || !toNum(a[9], uy) || !toNum(a[10], uz)) {
          errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：up 向量不是三个合法数值");
          continue;
        }
        up = Vec3{ux, uy, uz};
      }
      const Vec3 pi = m.coord(ni), pj = m.coord(nj);
      m.addSpring(ni, nj, ks, up);
      if (pi.x == pj.x && pi.y == pj.y && pi.z == pj.z) {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：弹簧两端节点坐标相同（长度为零）");
      }
      continue;
    }

    // ---------------- 端部释放（铰） ----------------
    //  release <目标> <参数...> <分量...>
    //    release beam  0  i  rx ry      （第 0 根梁的 i 端释放 rx、ry）
    //    release node  5  rx ry rz      （节点 5 上连接的所有梁，都在该端
    //                                     释放 rx、ry、rz）
    //  分量名：ux uy uz rx ry rz（单元局部坐标）。
    //  release 必须在单元建好之后使用（"建好"= addBeam / mesh 生成之后）。
    if (c.name == "release") {
      if (!need(c, 2)) continue;
      // 分量名 → 节点内自由度号；未知返回 -1
      auto compIndex = [](const std::string& s) -> int {
        if (s == "ux") return 0;
        if (s == "uy") return 1;
        if (s == "uz") return 2;
        if (s == "rx") return 3;
        if (s == "ry") return 4;
        if (s == "rz") return 5;
        return -1;
      };
      const std::string tgt = lower(a[0]);
      // 分量从第几个参数开始：beam = 3（跳过 序号 和 端名），node = 2
      const size_t compStart = (tgt == "beam") ? 3 : 2;
      std::vector<int> comps;
      bool compErr = false;
      for (size_t t = compStart; t < a.size(); ++t) {
        const int cc = compIndex(lower(a[t]));
        if (cc < 0) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：未知分量 '" + a[t] +
                            "'（应为 ux uy uz rx ry rz）");
          compErr = true;
        } else {
          comps.push_back(cc);
        }
      }
      if (compErr || comps.empty()) {
        if (!compErr) errors_.push_back("第 " + std::to_string(c.line) +
                                        " 行：release 至少需要一个分量（ux uy uz rx ry rz）");
        continue;
      }
      if (tgt == "beam") {
        // release beam <序号> <i|j> <分量...>
        if (!need(c, 4)) continue;
        int id = 0;
        if (!toInt(a[1], id) || id < 0) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：梁序号 " + a[1] +
                            " 不是非负整数");
          continue;
        }
        const std::string end = lower(a[2]);
        int localNode = -1;
        if (end == "i" || end == "0") localNode = 0;
        else if (end == "j" || end == "1") localNode = 1;
        else {
          errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：梁端应为 i 或 j（也可写 0 / 1）");
          continue;
        }
        int seen = 0;
        bool applied = false;
        for (auto& e : m.elements()) {
          if (e->type() != ElementType::Beam3D) continue;
          if (seen++ != id) continue;
          BeamElement* be = static_cast<BeamElement*>(e.get());
          for (int cc : comps) be->releaseAt(localNode, cc);
          applied = true;
          break;
        }
        if (!applied) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：梁序号 " + a[1] +
                            " 越界（模型里梁单元（含轴网生成的）共 " +
                            std::to_string(seen) + " 根）");
        }
      } else if (tgt == "node") {
        // release node <节点号> <分量...> —— 该节点上连接的每根梁的那一端
        if (!need(c, 3)) continue;
        int nid = 0;
        if (!toInt(a[1], nid) || nid < 0 || nid >= m.nodeCount()) {
          errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号 " + a[1] +
                            " 越界（有效范围 0.." + std::to_string(m.nodeCount() - 1) + "）");
          continue;
        }
        int hit = 0;
        for (auto& e : m.elements()) {
          if (e->type() != ElementType::Beam3D) continue;
          const auto& ns = e->nodes();
          int localNode = -1;
          if (ns[0] == nid) localNode = 0;
          else if (ns[1] == nid) localNode = 1;
          if (localNode < 0) continue;
          BeamElement* be = static_cast<BeamElement*>(e.get());
          for (int cc : comps) be->releaseAt(localNode, cc);
          ++hit;
        }
        if (hit == 0) {
          warnings_.push_back("第 " + std::to_string(c.line) + " 行：节点 " + a[1] +
                              " 上没有任何梁单元，release 未生效");
        }
      } else {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：release 目标 '" + a[0] +
                          "' 未知（beam / node）");
      }
      continue;
    }

    // ---------------- 荷载 ----------------
    if (c.name == "selfweight" || c.name == "slabload" ||
        c.name == "beamload" || c.name == "nodeload" ||
        c.name == "nodemoment") {
      if (sawCaseCmd_) {
        // 多工况模式：荷载不在此处施加，按工况分组，等 applyCase() 重放
        caseLoads_[currentCase_].push_back(c);
        continue;
      }
      if (!applyLoadCmd(m, c, info) && errors_.empty())
        errors_.push_back("第 " + std::to_string(c.line) + " 行：荷载命令 '" +
                          c.name + "' 处理失败");
      continue;
    }
    if (c.name == "slabload.to") {
      // slabload.to auto|shell|beam —— 面荷载的去处
      if (!need(c, 1)) continue;
      const std::string v = lower(a[0]);
      if (v == "auto") slabOpt_.path = GridMesh::SlabLoadOptions::Path::Auto;
      else if (v == "shell" || v == "shells") slabOpt_.path = GridMesh::SlabLoadOptions::Path::ToShells;
      else if (v == "beam" || v == "beams") slabOpt_.path = GridMesh::SlabLoadOptions::Path::ToBeams;
      else { errors_.push_back("第 " + std::to_string(c.line) +
                               " 行：slabload.to 的参数应为 auto / shell / beam"); continue; }
      continue;
    }
    if (c.name == "slabload.ratio") {
      if (!need(c, 1)) continue;
      double r = 2.0;
      if (!toNum(a[0], r) || !(r >= 1.0)) {
        errors_.push_back("第 " + std::to_string(c.line) +
                          " 行：长短边比阈值必须是不小于 1 的数");
        continue;
      }
      slabOpt_.oneWayRatio = r;
      continue;
    }
    if (c.name == "out") continue;      // 已在 parse 里处理

    errors_.push_back("第 " + std::to_string(c.line) + " 行：未知命令 '" + c.name + "'");
  }

  // 多工况模式：预生成网格缓存，applyCase() 施加板荷载/顶节点荷载时复用
  // （生成网格会改动 m，所以必须在建模阶段完成，不能在 applyCase 里补）。
  if (sawCaseCmd_ && hasGrid_ && !mesh_) ensureMesh();

  return errors_.empty();
}

// -----------------------------------------------------------------------------
//  网格缓存（build 与 applyCase 共用）
// -----------------------------------------------------------------------------
bool ModelScript::ensureMeshCached(Model& m, std::string* info) {
  if (mesh_) return true;
  if (!hasGrid_) {
    errors_.push_back("支座/荷载命令出现在轴网定义之前（缺少 grid.axisX/axisY/story）");
    return false;
  }
  if (grid_.nx() < 2 || grid_.ny() < 2 || grid_.nStories() < 2) {
    errors_.push_back("轴网不完整：至少需要 2 条 X 向轴线、2 条 Y 向轴线、2 个标高");
    return false;
  }
  mesh_ = std::make_unique<GridMesh>(grid_, meshSpec_);
  const GridMesh::Result r = mesh_->generate(m);
  if (!r.ok) {
    errors_.push_back("网格生成失败：" + r.error);
    return false;
  }
  meshBuilt_ = true;
  if (info) *info = "网格：" + r.message();
  return true;
}

// -----------------------------------------------------------------------------
//  单条荷载命令的执行（build 的直接施加 与 applyCase 的重放 共用）
// -----------------------------------------------------------------------------
bool ModelScript::applyLoadCmd(Model& m, const Command& c, std::string* info) {
  const auto& a = c.args;

  if (c.name == "selfweight") {
    const bool on = !(a.size() >= 1 && lower(a[0]) == "off");
    for (auto& e : m.elements()) {
      if (e->type() != ElementType::Beam3D) continue;
      static_cast<BeamElement*>(e.get())->setSelfWeight(on);
    }
    // 剪力墙自重（T4）：只对墙壳开关，水平板壳不吃自重 ——
    // 板的重量由 slabload 面荷载体现，若把板也包进来，既有
    // 轴网模型会凭空多出一大块重量（板自重），破坏 18 个测试。
    for (auto& e : m.elements()) {
      if (e->type() != ElementType::Shell4) continue;
      ShellElement* s = static_cast<ShellElement*>(e.get());
      if (s->isWall()) s->setSelfWeight(on);
    }
    return true;
  }

  if (c.name == "slabload") {
    if (a.empty()) { errors_.push_back("第 " + std::to_string(c.line) +
                     " 行：slabload 需要一个面荷载数值"); return true; }
    double q = 0;
    if (!toNum(a[0], q)) { errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：面荷载不是合法数值"); return true; }
    // 有轴网 → 网格已由 ensureMeshCached 生成（build 尾部或此处）；
    // 否则（手工模型）只能加压在已有的壳单元上。
    if (!mesh_ && hasGrid_) ensureMeshCached(m, info);
    if (mesh_) {
      slabOpt_.q = q;
      const auto rep = mesh_->applySlabLoad(m, slabOpt_);
      if (info) *info += "；" + rep.message();
      for (const auto& w : rep.warnings)
        warnings_.push_back("第 " + std::to_string(c.line) + " 行：" + w);
      return true;
    }
    int n = 0;
    for (auto& e : m.elements()) {
      if (e->type() != ElementType::Shell4) continue;
      // 脚本里正数表示向下（工程习惯），单元里 -z 为向下
      static_cast<ShellElement*>(e.get())->setTransversePressure(-q);
      ++n;
    }
    if (n == 0) warnings_.push_back("第 " + std::to_string(c.line) +
                                    " 行：模型里没有板单元，slabload 未生效");
    return true;
  }

  if (c.name == "beamload") {
    if (a.empty()) { errors_.push_back("第 " + std::to_string(c.line) +
                     " 行：beamload 需要一个线荷载数值"); return true; }
    double q = 0;
    if (!toNum(a[0], q)) { errors_.push_back("第 " + std::to_string(c.line) +
                            " 行：线荷载不是合法数值"); return true; }
    int n = 0;
    for (auto& e : m.elements()) {
      if (e->type() != ElementType::Beam3D) continue;
      auto* b = static_cast<BeamElement*>(e.get());
      // 只加在水平构件上：柱上加竖向均布没有意义且会污染轴力
      if (std::abs(b->core().localX().z) > 0.7) continue;
      b->setLineLoad(Vec3{0, 0, -q});
      ++n;
    }
    if (n == 0) warnings_.push_back("第 " + std::to_string(c.line) +
                                    " 行：模型里没有水平梁，beamload 未生效");
    return true;
  }

  if (c.name == "nodeload") {
    if (a.size() < 3) { errors_.push_back("第 " + std::to_string(c.line) +
                       " 行：nodeload 需要 3 个力分量"); return true; }
    Vec3 f{0, 0, 0};
    // nodeload <目标> <目标号> <fx> <fy> <fz>；top 没有目标号
    const size_t base = (lower(a[0]) == "top") ? 1u : 2u;
    if (a.size() < base + 3) {
      errors_.push_back("第 " + std::to_string(c.line) + " 行：nodeload 需要 3 个力分量");
      return true;
    }
    bool ok = true;
    for (int k = 0; k < 3; ++k)
      if (!toNum(a[base + static_cast<size_t>(k)], f[k])) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：力分量 '" +
                          a[base + static_cast<size_t>(k)] + "' 不是合法数值");
        ok = false;
      }
    if (!ok) return true;
    const std::string tgt = lower(a[0]);
    if (tgt == "node") {
      int id = 0;
      if (!toInt(a[1], id) || id < 0 || id >= m.nodeCount()) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号越界");
        return true;
      }
      m.addNodeForce(static_cast<Id>(id), f);
    } else if (tgt == "story") {
      int s = 0;
      if (!toInt(a[1], s)) { errors_.push_back("第 " + std::to_string(c.line) +
                             " 行：层号不是整数"); return true; }
      int cnt = 0;
      for (Id i = 0; i < m.nodeCount(); ++i) {
        if (m.node(i).diaphragmMaster) continue;   // 见下面 top 分支的说明
        if (m.node(i).story == s) { m.addNodeForce(i, f); ++cnt; }
      }
      if (cnt == 0) warnings_.push_back("第 " + std::to_string(c.line) + " 行：第 " +
                                        std::to_string(s) + " 层没有节点");
    } else if (tgt == "top") {
      // 网格在 build() 阶段已预生成（applyCase 调用点在建模之后）
      int top = -1;
      for (Id i = 0; i < m.nodeCount(); ++i) top = std::max(top, m.node(i).story);
      int cnt = 0;
      for (Id i = 0; i < m.nodeCount(); ++i) {
        // 【必须跳过刚性楼板主节点】它按楼层分类，但它是计算过程的产物
        // ——对它施加"整层水平力"会让荷载翻倍：主节点的力与从属节点的力
        // 最终都落在同一个方程上。而这个错误在残差、反力平衡里都看不出来。
        if (m.node(i).diaphragmMaster) continue;
        if (m.node(i).story == top) { m.addNodeForce(i, f); ++cnt; }
      }
      if (cnt == 0) warnings_.push_back("第 " + std::to_string(c.line) +
                                        " 行：顶层没有节点");
    } else {
      errors_.push_back("第 " + std::to_string(c.line) + " 行：nodeload 目标 '" + a[0] +
                        "' 未知（node / story / top）");
    }
    return true;
  }

  if (c.name == "nodemoment") {
    // nodemoment <目标> <目标号> <mx> <my> <mz>；top 没有目标号。
    // 与 nodeload 同构（P3：节点力矩荷载，施加到节点力矩）。
    if (a.size() < 3) { errors_.push_back("第 " + std::to_string(c.line) +
                       " 行：nodemoment 需要 3 个力矩分量"); return true; }
    Vec3 g{0, 0, 0};
    const size_t base = (lower(a[0]) == "top") ? 1u : 2u;
    if (a.size() < base + 3) {
      errors_.push_back("第 " + std::to_string(c.line) + " 行：nodemoment 需要 3 个力矩分量");
      return true;
    }
    bool ok = true;
    for (int k = 0; k < 3; ++k)
      if (!toNum(a[base + static_cast<size_t>(k)], g[k])) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：力矩分量 '" +
                          a[base + static_cast<size_t>(k)] + "' 不是合法数值");
        ok = false;
      }
    if (!ok) return true;
    const std::string tgt = lower(a[0]);
    if (tgt == "node") {
      int id = 0;
      if (!toInt(a[1], id) || id < 0 || id >= m.nodeCount()) {
        errors_.push_back("第 " + std::to_string(c.line) + " 行：节点号越界");
        return true;
      }
      m.addNodeMoment(static_cast<Id>(id), g);
    } else if (tgt == "story") {
      int s = 0;
      if (!toInt(a[1], s)) { errors_.push_back("第 " + std::to_string(c.line) +
                             " 行：层号不是整数"); return true; }
      int cnt = 0;
      for (Id i = 0; i < m.nodeCount(); ++i) {
        if (m.node(i).diaphragmMaster) continue;
        if (m.node(i).story == s) { m.addNodeMoment(i, g); ++cnt; }
      }
      if (cnt == 0) warnings_.push_back("第 " + std::to_string(c.line) + " 行：第 " +
                                        std::to_string(s) + " 层没有节点");
    } else if (tgt == "top") {
      int top = -1;
      for (Id i = 0; i < m.nodeCount(); ++i) top = std::max(top, m.node(i).story);
      int cnt = 0;
      for (Id i = 0; i < m.nodeCount(); ++i) {
        if (m.node(i).diaphragmMaster) continue;
        if (m.node(i).story == top) { m.addNodeMoment(i, g); ++cnt; }
      }
      if (cnt == 0) warnings_.push_back("第 " + std::to_string(c.line) +
                                        " 行：顶层没有节点");
    } else {
      errors_.push_back("第 " + std::to_string(c.line) + " 行：nodemoment 目标 '" + a[0] +
                        "' 未知（node / story / top）");
    }
    return true;
  }

  return false;
}

// -----------------------------------------------------------------------------
//  工况荷载重放（多工况分析流程：对每个工况调用一次）
//
//  【前置】build() 已经建好几何/约束/支座并生成网格缓存（mesh_）。
//  【语义】只施加工况名下的荷载命令，不改动几何与约束；
//  同一工况重复调用会把荷载叠加（CLI 保证先 clearLoads() 再调用）。
// -----------------------------------------------------------------------------
bool ModelScript::applyCase(Model& m, const std::string& name, std::string* info) {
  const auto it = caseLoads_.find(name);
  if (it == caseLoads_.end()) {
    errors_.push_back("applyCase：工况 '" + name + "' 未定义");
    return false;
  }
  for (const Command& c : it->second) {
    if (!applyLoadCmd(m, c, info) && errors_.empty())
      errors_.push_back("第 " + std::to_string(c.line) + " 行：荷载命令 '" +
                        c.name + "' 处理失败");
  }
  return errors_.empty();
}

}  // namespace io
}  // namespace yjk
