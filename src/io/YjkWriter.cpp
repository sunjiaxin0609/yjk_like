// =============================================================================
//  src/io/YjkWriter.cpp  ——  Model → .yjk 脚本文本序列化
// =============================================================================
#include "yjk/io/YjkWriter.h"

#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

namespace yjk {
namespace io {

namespace {

// 材料按力学字段去重。name 是显示名（concreteC 生成 "C30" 等），
// 序列化时用 "M<序号>" 重命名，读回后力学属性一致即可。
struct MatKey {
  double E, nu, gamma, fy, ft, fc, G;
  bool operator<(const MatKey& o) const {
    if (E != o.E) return E < o.E;
    if (nu != o.nu) return nu < o.nu;
    if (gamma != o.gamma) return gamma < o.gamma;
    if (fy != o.fy) return fy < o.fy;
    if (ft != o.ft) return ft < o.ft;
    if (fc != o.fc) return fc < o.fc;
    return G < o.G;
  }
};
MatKey matKey(const Material& m) {
  return {m.E, m.nu, m.gamma, m.fy, m.ft, m.fc, m.G};
}

// 截面按全部几何性质去重。
struct SecKey {
  double A, Iy, Iz, J, Asy, Asz, Ry, Rz, hy, hz;
  bool operator<(const SecKey& o) const {
    if (A != o.A) return A < o.A;
    if (Iy != o.Iy) return Iy < o.Iy;
    if (Iz != o.Iz) return Iz < o.Iz;
    if (J != o.J) return J < o.J;
    if (Asy != o.Asy) return Asy < o.Asy;
    if (Asz != o.Asz) return Asz < o.Asz;
    if (Ry != o.Ry) return Ry < o.Ry;
    if (Rz != o.Rz) return Rz < o.Rz;
    if (hy != o.hy) return hy < o.hy;
    return hz < o.hz;
  }
};
SecKey secKey(const SectionProperties& s) {
  return {s.A, s.Iy, s.Iz, s.J, s.Asy, s.Asz, s.Ry, s.Rz, s.hy, s.hz};
}

const char* boolStr(bool b) { return b ? "1" : "0"; }

}  // namespace

std::string modelToYjk(const Model& m) {
  std::ostringstream os;
  os << std::setprecision(17);

  // ---- 材料 / 截面去重命名表 ----
  // 全部在单元引用之前输出 —— ModelScript build() 边执行边填充
  // 命名表，`beam S0 M0` 执行到那一刻表里必须已有定义。
  std::map<MatKey, std::string> matNames;
  std::map<SecKey, std::string> secNames;
  auto matName = [&](const Material& mt) {
    const MatKey k = matKey(mt);
    auto it = matNames.find(k);
    if (it != matNames.end()) return it->second;
    const std::string name = "M" + std::to_string(matNames.size());
    matNames.emplace(k, name);
    return name;
  };
  auto secName = [&](const SectionProperties& sp) {
    const SecKey k = secKey(sp);
    auto it = secNames.find(k);
    if (it != secNames.end()) return it->second;
    const std::string name = "S" + std::to_string(secNames.size());
    secNames.emplace(k, name);
    return name;
  };

  // 第一次遍历：收集所有被单元引用的材料/截面，保证命名表先行输出。
  for (const auto& e : m.elements()) {
    switch (e->type()) {
      case ElementType::Beam3D: {
        const auto* b = static_cast<const BeamElement*>(e.get());
        matName(b->material());
        secName(b->section());
        break;
      }
      case ElementType::Shell4:
        // 壳的"材料"内嵌在 ShellProperties 里（E/nu/density），
        // 序列化为 shell 命令的内联参数，不引用材料表。
        break;
      case ElementType::Spring3D:
        break;                          // 弹簧无材料/截面
    }
  }

  // ---- 输出材料命令 ----
  for (const auto& kv : matNames) {
    Material m_;
    m_.E = kv.first.E;
    m_.nu = kv.first.nu;
    m_.gamma = kv.first.gamma;
    m_.fy = kv.first.fy;
    m_.ft = kv.first.ft;
    m_.fc = kv.first.fc;
    m_.G = kv.first.G;
    os << "material " << kv.second << " raw "
       << m_.E << " " << m_.nu << " " << m_.gamma << " "
       << m_.fy << " " << m_.ft << " " << m_.fc << " "
       << m_.G << "\n";
  }

  // ---- 输出截面命令 ----
  for (const auto& kv : secNames) {
    os << "section " << kv.second << " raw "
       << kv.first.A << " " << kv.first.Iy << " " << kv.first.Iz << " "
       << kv.first.J << " " << kv.first.Asy << " " << kv.first.Asz << " "
       << kv.first.Ry << " " << kv.first.Rz << " " << kv.first.hy << " "
       << kv.first.hz << "\n";
  }

  // ---- 节点重映射 ----
  // 跳过 retired 节点；活跃节点按 0..n-1 重新编号（addNode 顺序分配 id）。
  std::vector<Id> remap(static_cast<size_t>(m.nodeCount()), kDofFixed);
  Id next = 0;
  for (Id n = 0; n < m.nodeCount(); ++n) {
    if (m.node(n).retired) continue;
    remap[static_cast<size_t>(n)] = next++;
  }

  for (Id n = 0; n < m.nodeCount(); ++n) {
    if (m.node(n).retired) continue;
    const Node& nd = m.node(n);
    os << "node " << nd.r.x << " " << nd.r.y << " " << nd.r.z
       << " " << nd.story << "\n";
  }

  // ---- 单元 ----
  for (const auto& e : m.elements()) {
    std::vector<std::string> ns;
    for (Id nid : e->nodes())
      ns.push_back(std::to_string(remap[static_cast<size_t>(nid)]));
    switch (e->type()) {
      case ElementType::Beam3D: {
        const auto* b = static_cast<const BeamElement*>(e.get());
        const std::string mn = matName(b->material());
        const std::string sn = secName(b->section());
        const Vec3 up = b->upHint();
        // 统一输出为 beam 命令并带显式 up —— 柱与梁在单元层无区别，
        // 读回后 up 保持一致（Column 的语义由 up 方向体现）。
        os << "beam " << ns[0] << " " << ns[1] << " " << sn << " " << mn
           << " " << up.x << " " << up.y << " " << up.z << "\n";
        break;
      }
      case ElementType::Shell4: {
        const auto* s = static_cast<const ShellElement*>(e.get());
        const ShellProperties& p = s->properties();
        const std::string kind = s->isWall() ? "wall" : "shell";
        os << kind << " " << ns[0] << " " << ns[1] << " " << ns[2] << " "
           << ns[3] << " " << p.thickness << " " << p.E << " " << p.nu
           << " " << p.density << "\n";
        break;
      }
      case ElementType::Spring3D: {
        const auto* sp = static_cast<const SpringElement*>(e.get());
        const std::array<double, 6>& k = sp->stiffness();
        const Vec3 up = sp->upHint();
        os << "spring " << ns[0] << " " << ns[1];
        for (int i = 0; i < 6; ++i) os << " " << k[static_cast<size_t>(i)];
        os << " " << up.x << " " << up.y << " " << up.z << "\n";
        break;
      }
    }
  }

  // ---- 约束 ----
  for (Id n = 0; n < m.nodeCount(); ++n) {
    if (m.node(n).retired) continue;
    const Node& nd = m.node(n);
    if (!nd.fixed[0] && !nd.fixed[1] && !nd.fixed[2] && !nd.fixed[3] &&
        !nd.fixed[4] && !nd.fixed[5])
      continue;
    os << "fix node " << remap[static_cast<size_t>(n)];
    for (int k = 0; k < 6; ++k) os << " " << boolStr(nd.fixed[k]);
    os << "\n";
  }

  // ---- 节点荷载 ----
  for (Id n = 0; n < m.nodeCount(); ++n) {
    if (m.node(n).retired) continue;
    const Node& nd = m.node(n);
    const Id rn = remap[static_cast<size_t>(n)];
    if (nd.force.x != 0.0 || nd.force.y != 0.0 || nd.force.z != 0.0)
      os << "nodeload node " << rn << " " << nd.force.x << " " << nd.force.y
         << " " << nd.force.z << "\n";
    if (nd.moment.x != 0.0 || nd.moment.y != 0.0 || nd.moment.z != 0.0)
      os << "nodemoment node " << rn << " " << nd.moment.x << " " << nd.moment.y
         << " " << nd.moment.z << "\n";
    if (nd.weight != 0.0)
      os << "nodeweight node " << rn << " " << nd.weight << "\n";
  }

  return os.str();
}

}  // namespace io
}  // namespace yjk