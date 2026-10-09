// =============================================================================
//  apps/ui/PropertyPanel.cpp
// =============================================================================
#include "PropertyPanel.h"

#include <QFont>
#include <QHeaderView>
#include <QLabel>
#include <QTreeWidget>
#include <QVBoxLayout>

#include <cmath>

#include "Theme.h"
#include "yjk/element/BeamElement3D.h"
#include "yjk/element/ShellElement4.h"

namespace ui {

// PropertyPanel 用到的内核类型里，BeamElement/ShellElement/Element/Model/
// StaticResult/Vec3 已由 SceneData.h using 引入；截面/材料/壳属性/工厂
// 不在那个名单里，这里补齐（只引类型与工厂命名空间，不引符号）。
using yjk::Material;
using yjk::SectionProperties;
using yjk::ShellProperties;
namespace section = yjk::section;

namespace {

// 数值行：左右两列，单位写在标签里。
// 【单位必须写在标签里】结构工程里同一个量有多个单位在流通
// （kN / kN·m / kPa / MPa / m / mm），不写单位的话
// "M=45.2" 到底是 kN·m 还是 kN·mm 只能靠猜，而差了 1000 倍。
QTreeWidgetItem* kv(QTreeWidgetItem* parent, const QString& k, const QString& v,
                    const QColor& vColor = QColor()) {
  auto* it = new QTreeWidgetItem(parent, {k, v});
  it->setFont(0, uiFont(9));
  it->setFont(1, monoFont(9));
  if (vColor.isValid()) it->setForeground(1, vColor);
  return it;
}

// 可编辑数值行：在 kv 基础上挂字段键 + 编辑 flag。
// 【为什么颜色用 info】结构软件的属性面板里，"能改"比"当前值"更重要
// —— 用户扫一眼就知道哪些是输入、哪些是只读结果。
QTreeWidgetItem* kvEditable(QTreeWidgetItem* parent, const QString& k, const QString& v,
                            PropField f) {
  auto* it = kv(parent, k, v, Status::info());
  it->setData(1, kFieldRole, static_cast<int>(f));
  it->setFlags(it->flags() | Qt::ItemIsEditable);
  return it;
}

QTreeWidgetItem* group(QTreeWidget* tree, const QString& name) {
  auto* g = new QTreeWidgetItem(tree, {name, QString()});
  QFont f = uiFont(9);
  f.setBold(true);
  g->setFont(0, f);
  g->setExpanded(true);
  return g;
}

}  // namespace

PropertyPanel::PropertyPanel(QWidget* parent) : QWidget(parent) {
  auto* lay = new QVBoxLayout(this);
  lay->setContentsMargins(4, 4, 4, 4);
  lay->setSpacing(3);

  title_ = new QLabel(QStringLiteral("未选中对象"), this);
  QFont tf = uiFont(9);
  tf.setBold(true);
  title_->setFont(tf);
  title_->setStyleSheet(QStringLiteral("color:%1;").arg(Status::info().name()));
  lay->addWidget(title_);

  tree_ = new QTreeWidget(this);
  tree_->setColumnCount(2);
  tree_->setHeaderLabels({QStringLiteral("项"), QStringLiteral("值")});
  tree_->header()->setStretchLastSection(true);
  tree_->header()->resizeSection(0, 108);
  tree_->setRootIsDecorated(true);
  tree_->setAlternatingRowColors(true);
  tree_->setUniformRowHeights(true);
  tree_->setFont(uiFont(9));
  // 双击进入编辑；回车提交，Esc 恢复原值由 Qt 内部处理
  tree_->setEditTriggers(QAbstractItemView::DoubleClicked |
                         QAbstractItemView::EditKeyPressed);
  lay->addWidget(tree_, 1);

  // 编辑提交 → 解析 → 写回 → 刷新。rebuilding_ 防"rebuild 里 setText
  // 又触发 itemChanged"的死循环（QTreeWidget::itemChanged 在数据变化时
  // 就会发，不只是用户编辑完）。
  connect(tree_, &QTreeWidget::itemChanged, this, &PropertyPanel::onItemChanged);
}

void PropertyPanel::setModel(Model* m) { model_ = m; rebuild(); }
void PropertyPanel::setResult(const StaticResult* r) { result_ = r; rebuild(); }
void PropertyPanel::setPost(const post::PostProcessor* pp) { post_ = pp; rebuild(); }

void PropertyPanel::setEditEnabled(bool on) { editEnabled_ = on; rebuild(); }

void PropertyPanel::clear() {
  curElem_ = -1;
  curType_ = 0;
  curNode_ = -1;
  rebuild();
}

void PropertyPanel::showElement(int elem, int type) {
  curElem_ = elem;
  curType_ = type;
  curNode_ = -1;
  rebuild();
}

void PropertyPanel::showNode(int node) {
  curNode_ = node;
  curElem_ = -1;
  curType_ = 0;
  rebuild();
}

// -----------------------------------------------------------------------------
//  编辑处理
// -----------------------------------------------------------------------------
void PropertyPanel::onItemChanged(QTreeWidgetItem* item, int column) {
  if (rebuilding_ || !model_ || !item) return;
  if (column != 1) return;                            // 只编辑"值"列

  const QVariant fv = item->data(1, kFieldRole);
  if (!fv.isValid()) return;                          // 只读行
  const PropField f = static_cast<PropField>(fv.toInt());
  if (f == PropField::None) return;

  bool okNum = false;
  const double v = item->text(1).toDouble(&okNum);
  if (!okNum || !std::isfinite(v)) {
    emit editRejected(QStringLiteral("数值无效：“%1”（已恢复原值）").arg(item->text(1)));
    rebuild();                                        // 恢复显示原值
    return;
  }

  if (!applyEdit(f, v)) {
    rebuild();                                        // 写回失败（含内核回滚），恢复原值
    return;
  }
  rebuild();                                          // 成功：刷新联动值（如截面 A/Iy/J 随 b/h 变）
  emit modelEdited();
}

bool PropertyPanel::applyEdit(PropField f, double v) {
  // ---- 节点 ----
  if (curNode_ >= 0 && curNode_ < model_->nodeCount()) {
    switch (f) {
      case PropField::NodeX: case PropField::NodeY: case PropField::NodeZ: {
        Vec3 r = model_->node(curNode_).r;
        if (f == PropField::NodeX) r.x = v;
        else if (f == PropField::NodeY) r.y = v;
        else r.z = v;
        if (!model_->setNodeCoord(curNode_, r)) {
          emit editRejected(QStringLiteral("坐标修改使单元退化（零长度/共线），已回滚"));
          return false;
        }
        return true;
      }
      case PropField::NodeFx: case PropField::NodeFy: case PropField::NodeFz: {
        Vec3 F = model_->node(curNode_).force;
        if (f == PropField::NodeFx) F.x = v;
        else if (f == PropField::NodeFy) F.y = v;
        else F.z = v;
        model_->setNodeForce(curNode_, F);
        return true;
      }
      case PropField::NodeMx: case PropField::NodeMy: case PropField::NodeMz: {
        Vec3 M = model_->node(curNode_).moment;
        if (f == PropField::NodeMx) M.x = v;
        else if (f == PropField::NodeMy) M.y = v;
        else M.z = v;
        model_->setNodeMoment(curNode_, M);
        return true;
      }
      default: break;
    }
  }

  // ---- 单元 ----
  if (curElem_ < 0 || curElem_ >= model_->elementCount()) return false;
  Element* el = model_->element(curElem_);
  if (!el) return false;

  if (el->type() == ElementType::Beam3D) {
    auto* be = static_cast<BeamElement*>(el);
    const SectionProperties sec0 = be->section();
    const Material mat0 = be->material();

    // 截面：矩形等效 b/h 编辑。
    // 【为什么用 b/h 而不是 A/Iy/Iz】用户写"300×500"直觉是宽×高；
    // 直接改 A/Iy/Iz 需要同时维护四个耦合的量，任何一想到处错。
    // 对矩形截面 b=2·hz、h=2·hy 精确互逆；非矩形截面（工字钢等）
    // 显示的是等效外廓，编辑后按矩形重建（1-1 编辑语义，可接受）。
    switch (f) {
      case PropField::BeamB: case PropField::BeamH: {
        double b = 2.0 * sec0.halfWidthZ();
        double h = 2.0 * sec0.halfDepthY();
        if (f == PropField::BeamB) b = v; else h = v;
        if (!(b > 0) || !(h > 0) || !std::isfinite(b) || !std::isfinite(h)) {
          emit editRejected(QStringLiteral("截面尺寸必须为正数"));
          return false;
        }
        if (!model_->setBeamProps(curElem_, section::rect(b, h), mat0)) {
          emit editRejected(QStringLiteral("截面写回失败"));
          return false;
        }
        return true;
      }
      case PropField::BeamE: case PropField::BeamNu: case PropField::BeamGamma: {
        Material m = mat0;
        if (f == PropField::BeamE) m.E = v;
        else if (f == PropField::BeamNu) m.nu = v;
        else m.gamma = v;
        if (!(m.E > 0) || !std::isfinite(m.E)) { emit editRejected(QStringLiteral("弹性模量必须为正")); return false; }
        if (!(m.nu > -1.0 && m.nu < 0.5)) { emit editRejected(QStringLiteral("泊松比应为 (-1, 0.5)")); return false; }
        if (!(m.gamma > 0) || !std::isfinite(m.gamma)) { emit editRejected(QStringLiteral("容重必须为正")); return false; }
        if (!model_->setBeamProps(curElem_, sec0, m)) {
          emit editRejected(QStringLiteral("材料写回失败"));
          return false;
        }
        return true;
      }
      default: break;
    }
  } else if (el->type() == ElementType::Shell4) {
    auto* sh = static_cast<ShellElement*>(el);
    ShellProperties p = sh->properties();
    switch (f) {
      case PropField::ShellT: p.thickness = v; break;
      case PropField::ShellE: p.E = v; break;
      case PropField::ShellNu: p.nu = v; break;
      case PropField::ShellDensity: p.density = v; break;
      default: break;
    }
    std::string why;
    if (!p.isValid(&why)) {
      emit editRejected(QString::fromStdString(why));
      return false;
    }
    if (!model_->setShellProps(curElem_, p)) {
      emit editRejected(QStringLiteral("壳属性写回失败"));
      return false;
    }
    return true;
  }

  return false;
}

// -----------------------------------------------------------------------------
//  重建
// -----------------------------------------------------------------------------
void PropertyPanel::rebuild() {
  rebuilding_ = true;                        // 防 itemChanged 死循环
  tree_->clear();

  const bool editable = editEnabled_ && model_ != nullptr;
  auto editableOr = [&](QTreeWidgetItem* g, const QString& k, const QString& v, PropField f) {
    return editable ? kvEditable(g, k, v, f) : kv(g, k, v);
  };

  if (!model_) {
    title_->setText(QStringLiteral("未加载模型"));
    rebuilding_ = false;
    return;
  }

  // ---- 节点 ----
  if (curNode_ >= 0 && curNode_ < model_->nodeCount()) {
    const Node& nd = model_->node(curNode_);
    title_->setText(QStringLiteral("节点 #%1").arg(nd.id));

    auto* g0 = group(tree_, QStringLiteral("坐标"));
    editableOr(g0, QStringLiteral("X (m)"), fmtNum(nd.r.x, 5), PropField::NodeX);
    editableOr(g0, QStringLiteral("Y (m)"), fmtNum(nd.r.y, 5), PropField::NodeY);
    editableOr(g0, QStringLiteral("Z (m)"), fmtNum(nd.r.z, 5), PropField::NodeZ);
    kv(g0, QStringLiteral("楼层"), QString::number(nd.story));

    auto* g1 = group(tree_, QStringLiteral("约束"));
    static const char* kn[6] = {"ux", "uy", "uz", "rx", "ry", "rz"};
    QStringList fixed;
    for (int k = 0; k < 6; ++k)
      if (nd.fixed[k]) fixed << kn[k];
    kv(g1, QStringLiteral("固定自由度"),
       fixed.isEmpty() ? QStringLiteral("自由") : fixed.join(","));
    if (!nd.retired) kv(g1, QStringLiteral("状态"), QStringLiteral("有效"));
    else kv(g1, QStringLiteral("状态"), QStringLiteral("已合并（退役）"), Status::warn());

    auto* g2 = group(tree_, QStringLiteral("荷载"));
    editableOr(g2, QStringLiteral("Fx (kN)"), fmtNum(nd.force.x, 4), PropField::NodeFx);
    editableOr(g2, QStringLiteral("Fy (kN)"), fmtNum(nd.force.y, 4), PropField::NodeFy);
    editableOr(g2, QStringLiteral("Fz (kN)"), fmtNum(nd.force.z, 4), PropField::NodeFz);
    editableOr(g2, QStringLiteral("Mx (kN·m)"), fmtNum(nd.moment.x, 4), PropField::NodeMx);
    editableOr(g2, QStringLiteral("My (kN·m)"), fmtNum(nd.moment.y, 4), PropField::NodeMy);
    editableOr(g2, QStringLiteral("Mz (kN·m)"), fmtNum(nd.moment.z, 4), PropField::NodeMz);
    if (editable)
      kv(g2, QStringLiteral("提示"), QStringLiteral("双击数值可修改"), Status::warn());

    if (result_ && result_->ok && result_->u.size() >= static_cast<size_t>(model_->nodeCount()) * 6) {
      const size_t b = static_cast<size_t>(curNode_ * 6);
      auto* g3 = group(tree_, QStringLiteral("位移"));
      // 位移用 mm 显示：结构工程的位移量级是毫米级，
      // 用 m 显示会得到一串 0.00xxxx，读起来要数零。
      kv(g3, QStringLiteral("ux (mm)"), fmtNum(result_->u[b + 0] * 1e3, 4));
      kv(g3, QStringLiteral("uy (mm)"), fmtNum(result_->u[b + 1] * 1e3, 4));
      kv(g3, QStringLiteral("uz (mm)"), fmtNum(result_->u[b + 2] * 1e3, 4));
      kv(g3, QStringLiteral("水平合位移 (mm)"),
         fmtNum(std::sqrt(result_->u[b] * result_->u[b] + result_->u[b + 1] * result_->u[b + 1]) * 1e3, 4));

      auto* g4 = group(tree_, QStringLiteral("支座反力"));
      bool any = false;
      static const char* rn[6] = {"Rx (kN)", "Ry (kN)", "Rz (kN)", "Mx (kN·m)", "My (kN·m)", "Mz (kN·m)"};
      for (int k = 0; k < 6; ++k) {
        const double v = result_->reaction[b + k];
        if (std::abs(v) > 1e-12) any = true;
        kv(g4, QString::fromUtf8(rn[k]), fmtNum(v, 4));
      }
      if (!any) kv(g4, QStringLiteral("(无)"), QStringLiteral("该节点不是支座"));
    }
    rebuilding_ = false;
    return;
  }

  // ---- 单元 ----
  if (curElem_ < 0 || curElem_ >= model_->elementCount()) {
    title_->setText(QStringLiteral("未选中对象"));
    auto* g = group(tree_, QStringLiteral("提示"));
    kv(g, QStringLiteral("操作"), QStringLiteral("在三维视图或模型树中点选构件"));
    rebuilding_ = false;
    return;
  }

  const Element& el = *model_->elements()[static_cast<size_t>(curElem_)];
  const bool isBeam = (el.type() == ElementType::Beam3D);
  title_->setText(QStringLiteral("%1  #%2")
                      .arg(isBeam ? QStringLiteral("梁")
                                  : QStringLiteral("壳"))
                      .arg(curElem_));

  auto* g0 = group(tree_, QStringLiteral("拓扑"));
  QString ns;
  for (size_t i = 0; i < el.nodes().size(); ++i) {
    if (i) ns += " → ";
    ns += QString::number(el.nodes()[i]);
  }
  kv(g0, QStringLiteral("节点"), ns);
  if (isBeam) kv(g0, QStringLiteral("长度 (m)"), fmtNum(el.length(), 4));

  if (isBeam) {
    const auto& be = static_cast<const BeamElement&>(el);
    const auto& sec = be.section();
    const auto& mat = be.material();
    // 矩形等效宽/高（对矩形截面精确，非矩形显示外廓）
    const double b = 2.0 * sec.halfWidthZ();
    const double h = 2.0 * sec.halfDepthY();

    auto* g1 = group(tree_, QStringLiteral("截面"));
    editableOr(g1, QStringLiteral("宽 b (m)"), fmtNum(b, 5), PropField::BeamB);
    editableOr(g1, QStringLiteral("高 h (m)"), fmtNum(h, 5), PropField::BeamH);
    kv(g1, QStringLiteral("A (m²)"), fmtNum(sec.A, 4));
    kv(g1, QStringLiteral("Iy (m⁴)"), fmtNum(sec.Iy, 4));
    kv(g1, QStringLiteral("Iz (m⁴)"), fmtNum(sec.Iz, 4));
    kv(g1, QStringLiteral("J (m⁴)"), fmtNum(sec.J, 4));

    auto* g2 = group(tree_, QStringLiteral("材料"));
    kv(g2, QStringLiteral("名称"), QString::fromStdString(mat.name));
    editableOr(g2, QStringLiteral("E (kPa)"), fmtNum(mat.E, 4), PropField::BeamE);
    editableOr(g2, QStringLiteral("ν"), fmtNum(mat.nu, 3), PropField::BeamNu);
    editableOr(g2, QStringLiteral("容重 (kN/m³)"), fmtNum(mat.gamma, 3), PropField::BeamGamma);
    if (mat.fc > 0) kv(g2, QStringLiteral("f_c (kPa)"), fmtNum(mat.fc, 3));
    if (mat.fy > 0) kv(g2, QStringLiteral("f_y (kPa)"), fmtNum(mat.fy, 3));

    auto* g3 = group(tree_, QStringLiteral("荷载"));
    kv(g3, QStringLiteral("自重"), be.selfWeight() ? QStringLiteral("开") : QStringLiteral("关"));
    if (be.hasLineLoad()) {
      const Vec3 q = be.lineLoad();
      kv(g3, QStringLiteral("均布线荷载"),
         QStringLiteral("(%1, %2, %3) kN/m").arg(fmtNum(q.x, 3), fmtNum(q.y, 3), fmtNum(q.z, 3)));
    }
    if (editable)
      kv(g3, QStringLiteral("提示"), QStringLiteral("双击数值可修改"), Status::warn());

    // 结果
    if (post_ && curElem_ < static_cast<int>(post_->beamResults().size())) {
      const post::BeamResult& br = post_->beamResults()[static_cast<size_t>(curElem_)];

      auto* g4 = group(tree_, QStringLiteral("控制截面内力"));
      kv(g4, QStringLiteral("M_max (kN·m)"), fmtNum(br.Mmax, 4),
         Status::info());
      kv(g4, QStringLiteral("位置 ξ"), fmtNum(br.xiMmax, 3));
      kv(g4, QStringLiteral("V_max (kN)"), fmtNum(br.Vmax, 4));
      kv(g4, QStringLiteral("N 范围 (kN)"),
         fmtNum(br.Nmin, 4) + " ~ " + fmtNum(br.Nmax, 4));

      auto* g5 = group(tree_, QStringLiteral("应力"));
      kv(g5, QStringLiteral("σ_max (kPa)"), fmtNum(br.sigMax, 4));
      kv(g5, QStringLiteral("σ_min (kPa)"), fmtNum(br.sigMin, 4));
      kv(g5, QStringLiteral("σ_abs (kPa)"), fmtNum(br.sigAbsMax, 4));
      kv(g5, QStringLiteral("von Mises (kPa)"), fmtNum(br.vonMisesMax, 4));

      auto* g6 = group(tree_, QStringLiteral("验算"));
      if (br.util > 0.0) {
        const QColor c = br.util > 1.0 ? Status::fail()
                                        : (br.util > 0.85 ? Status::warn() : Status::ok());
        kv(g6, QStringLiteral("应力比"), fmtNum(br.util, 3), c);
        kv(g6, QStringLiteral("结论"),
           br.util > 1.0 ? QStringLiteral("超限") : QStringLiteral("满足"), c);
      } else {
        kv(g6, QStringLiteral("应力比"), QStringLiteral("无强度信息"));
      }
      kv(g6, QStringLiteral("自重 (kN)"), fmtNum(br.weight, 4));
    }
  } else {
    const auto& sh = static_cast<const ShellElement&>(el);
    const auto& pr = sh.properties();
    auto* g1 = group(tree_, QStringLiteral("几何"));
    kv(g1, QStringLiteral("面积 (m²)"), fmtNum(sh.area(), 4));
    const Vec3 n = sh.normal();
    kv(g1, QStringLiteral("法向"),
       QStringLiteral("(%1, %2, %3)").arg(fmtNum(n.x, 2), fmtNum(n.y, 2), fmtNum(n.z, 2)));
    auto* g2 = group(tree_, QStringLiteral("属性"));
    editableOr(g2, QStringLiteral("厚度 (m)"), fmtNum(pr.thickness, 4), PropField::ShellT);
    editableOr(g2, QStringLiteral("E (kPa)"), fmtNum(pr.E, 4), PropField::ShellE);
    editableOr(g2, QStringLiteral("ν"), fmtNum(pr.nu, 3), PropField::ShellNu);
    editableOr(g2, QStringLiteral("密度 (kg/m³)"), fmtNum(pr.density, 4), PropField::ShellDensity);
    kv(g2, QStringLiteral("自重 (kN)"), fmtNum(sh.mass(), 4));
    if (editable)
      kv(g2, QStringLiteral("提示"), QStringLiteral("双击数值可修改"), Status::warn());
    if (!sh.plateReady()) {
      auto* g3 = group(tree_, QStringLiteral("提示"));
      kv(g3, QStringLiteral("板弯曲"), QStringLiteral("未实现（仅膜刚度）"), Status::warn());
    }
  }

  rebuilding_ = false;
}

}  // namespace ui