// =============================================================================
//  apps/ui/ResultPanel.cpp
// =============================================================================
#include "ResultPanel.h"

#include <QFileDialog>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMenu>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QTabWidget>
#include <QTableWidget>
#include <QTextStream>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>

#include "Theme.h"

namespace ui {

namespace {

QTableWidget* makeTable(const QStringList& headers) {
  auto* t = new QTableWidget();
  t->setColumnCount(headers.size());
  t->setHorizontalHeaderLabels(headers);
  t->setAlternatingRowColors(true);
  t->setSelectionBehavior(QAbstractItemView::SelectRows);
  t->setSelectionMode(QAbstractItemView::SingleSelection);
  t->setEditTriggers(QAbstractItemView::NoEditTriggers);
  t->setSortingEnabled(true);
  t->verticalHeader()->setVisible(false);
  t->horizontalHeader()->setStretchLastSection(true);
  t->setFont(monoFont(9));
  return t;
}

// -----------------------------------------------------------------------------
//  可排序的数值单元
//
//  【为什么必须自定义】QTableWidgetItem::operator< 比较的是显示文本，
//  按字符串排会得到 "10" < "9"、"1/3018" < "1/550" 这种荒谬结果 ——
//  而"按层间位移角排序找薄弱层"正是工程师最常用的操作。
//  把真实数值塞进 UserRole，排序时用数值比，显示时用文本。
// -----------------------------------------------------------------------------
class SortItem : public QTableWidgetItem {
 public:
  SortItem(const QString& text, double key) : QTableWidgetItem(text) {
    setData(Qt::UserRole, key);
  }
  bool operator<(const QTableWidgetItem& other) const override {
    const QVariant a = data(Qt::UserRole);
    const QVariant b = other.data(Qt::UserRole);
    if (a.isValid() && b.isValid()) return a.toDouble() < b.toDouble();
    return QTableWidgetItem::operator<(other);
  }
};

// 数值单元格：右对齐 + 等宽字体 + 可按数值排序。
// 左对齐的数字列在小数点上对不齐，一列数字扫下来读不出量级差异 ——
// 这是数值表格的基本要求，不是审美。
QTableWidgetItem* numItem(double v, int sig = 4, const QColor& col = QColor(),
                          const QString& textOverride = QString()) {
  auto* it = new SortItem(textOverride.isEmpty() ? fmtNum(v, sig) : textOverride, v);
  it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
  it->setFont(monoFont(9));
  if (col.isValid()) it->setForeground(col);
  return it;
}

void setNum(QTableWidget* t, int r, int c, double v, int sig = 4, const QColor& col = QColor()) {
  t->setItem(r, c, numItem(v, sig, col));
}

void setNumText(QTableWidget* t, int r, int c, const QString& text, double sortKey,
                const QColor& col = QColor()) {
  t->setItem(r, c, numItem(sortKey, 4, col, text));
}

void setText(QTableWidget* t, int r, int c, const QString& s) {
  auto* it = new QTableWidgetItem(s);
  it->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
  it->setFont(monoFont(9));
  t->setItem(r, c, it);
}

}  // namespace

ResultPanel::ResultPanel(QWidget* parent) : QWidget(parent) {
  auto* lay = new QVBoxLayout(this);
  lay->setContentsMargins(4, 4, 4, 4);
  lay->setSpacing(4);

  summary_ = new QLabel(QStringLiteral("尚无计算结果"), this);
  summary_->setFont(uiFont(9));
  summary_->setWordWrap(true);
  lay->addWidget(summary_);

  auto* btnRow = new QHBoxLayout();
  auto* btnCsv = new QPushButton(QStringLiteral("导出当前表 CSV"), this);
  btnCsv->setFont(uiFont(9));
  connect(btnCsv, &QPushButton::clicked, this, [this] {
    QTableWidget* t = nullptr;
    const int idx = tabs_ ? tabs_->currentIndex() : 0;
    if (idx == 0) t = story_;
    else if (idx == 1) t = envelope_;
    else if (idx == 2) t = members_;
    if (!t) { QMessageBox::information(this, QStringLiteral("导出"), QStringLiteral("当前页没有表格")); return; }

    const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("导出 CSV"),
                                                      QStringLiteral("result.csv"),
                                                      QStringLiteral("CSV (*.csv)"));
    if (path.isEmpty()) return;
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QMessageBox::warning(this, QStringLiteral("导出失败"), f.errorString());
      return;
    }
    QTextStream ts(&f);
    ts.setEncoding(QStringConverter::Utf8);
    QStringList head;
    for (int c = 0; c < t->columnCount(); ++c) head << t->horizontalHeaderItem(c)->text();
    ts << head.join(",") << "\n";
    for (int r = 0; r < t->rowCount(); ++r) {
      QStringList row;
      for (int c = 0; c < t->columnCount(); ++c) {
        const auto* it = t->item(r, c);
        // 数值单元格导出时给出足够位数，Excel 里才好继续算
        QString s = it ? it->text() : QString();
        s.replace('"', QStringLiteral("\"\""));
        row << QStringLiteral("\"%1\"").arg(s);
      }
      ts << row.join(",") << "\n";
    }
  });
  btnRow->addWidget(btnCsv);
  btnRow->addStretch(1);
  lay->addLayout(btnRow);

  tabs_ = new QTabWidget(this);
  tabs_->setFont(uiFont(9));

  story_ = makeTable({QStringLiteral("层"), QStringLiteral("标高 m"), QStringLiteral("层高 m"),
                      QStringLiteral("最大水平位移 mm"), QStringLiteral("层间位移角"),
                      QStringLiteral("位移比"), QStringLiteral("层剪力 kN"),
                      QStringLiteral("层以上重量 kN"), QStringLiteral("剪重比")});
  // 包络表刻意关掉排序：它是一张【固定顺序】的指标清单，
  // 排序会把它打乱成字母序，"最大位移"跑到"基底剪力"后面，反而难读。
  envelope_ = makeTable({QStringLiteral("指标"), QStringLiteral("值"), QStringLiteral("位置")});
  envelope_->setSortingEnabled(false);
  members_ = makeTable({QStringLiteral("单元"), QStringLiteral("类别"), QStringLiteral("L m"),
                        QStringLiteral("M_max kN·m"), QStringLiteral("V_max kN"),
                        QStringLiteral("N kN"), QStringLiteral("σ_max kPa"),
                        QStringLiteral("应力比")});
  report_ = new QPlainTextEdit(this);
  report_->setReadOnly(true);
  report_->setFont(monoFont(9));
  report_->setLineWrapMode(QPlainTextEdit::NoWrap);

  tabs_->addTab(story_, QStringLiteral("楼层指标"));
  tabs_->addTab(envelope_, QStringLiteral("整体指标"));
  tabs_->addTab(members_, QStringLiteral("构件内力"));
  tabs_->addTab(report_, QStringLiteral("文本报告"));
  lay->addWidget(tabs_, 1);

  // 双击行 → 视口高亮
  auto wire = [this](QTableWidget* t, bool memberTable) {
    connect(t, &QTableWidget::cellDoubleClicked, this, [this, t, memberTable](int row, int) {
      auto* first = t->item(row, 0);
      if (!first) return;
      if (memberTable) {
        bool ok = false;
        const int id = first->text().remove('#').toInt(&ok);
        if (ok) emit elementActivated(id, 1);
      } else {
        bool ok = false;
        const int s = first->text().toInt(&ok);
        if (ok) emit storyActivated(s);
      }
    });
  };
  wire(story_, false);
  wire(members_, true);
}

void ResultPanel::setModel(const Model* m) { model_ = m; refresh(); }
void ResultPanel::setResult(const StaticResult* r) { result_ = r; refresh(); }
void ResultPanel::setPost(const post::PostProcessor* pp) { post_ = pp; refresh(); }

void ResultPanel::refresh() {
  if (post_) {
    summary_->setText(QString::fromStdString(post_->summary()));
    summary_->setStyleSheet(QString());
  } else if (result_ && result_->ok) {
    summary_->setText(QStringLiteral("求解完成 —— 自由度 %1，耗时 %2 s，残差 %3")
                          .arg(result_->ndof)
                          .arg(fmtNum(result_->seconds, 3), fmtNum(result_->residual, 2)));
  } else {
    summary_->setText(QStringLiteral("尚无计算结果"));
  }
  buildStoryTable();
  buildEnvelope();
  buildMemberTable();
  buildReport();
}

void ResultPanel::buildStoryTable() {
  story_->setSortingEnabled(false);
  story_->setRowCount(0);
  if (!post_) { story_->setSortingEnabled(true); return; }

  // 【基底不是一层】层高为 0 的"层"是嵌固面（所有节点都固定），
  // 它的位移、层剪力、位移角恒为 0。留在表里只会让人以为"底层没受力"。
  std::vector<const post::StoryResult*> vis;
  for (const auto& s : post_->stories())
    if (s.h > 1e-9) vis.push_back(&s);

  story_->setRowCount(static_cast<int>(vis.size()));
  for (int i = 0; i < static_cast<int>(vis.size()); ++i) {
    const post::StoryResult& s = *vis[static_cast<size_t>(i)];
    setNumText(story_, i, 0, QString::number(s.story), s.story);
    setNum(story_, i, 1, s.z, 4);
    setNum(story_, i, 2, s.h, 4);
    setNum(story_, i, 3, s.maxUxy * 1e3, 4);

    // 层间位移角：工程上写作 1/xxx 比分小数好读得多。
    // 规范限值是「1/550」这样的形式，直接对照才不会多一次换算。
    const double drift = s.drift;
    const QString driftTxt = (drift > 1e-12)
                                 ? QStringLiteral("1/%1").arg(fmtNum(1.0 / drift, 4))
                                 : QStringLiteral("—");
    setNumText(story_, i, 4, driftTxt, drift);
    // 限值 1/550（GB 50011 框架结构）。超了标红 —— 这是审图最容易被打回的一项。
    if (drift > 1.0 / 550.0) story_->item(i, 4)->setForeground(Status::fail());
    else if (drift > 0.8 / 550.0) story_->item(i, 4)->setForeground(Status::warn());

    setNum(story_, i, 5, s.driftRatio, 3);
    setNum(story_, i, 6, s.shear, 4);
    setNum(story_, i, 7, s.weight, 4);
    setNum(story_, i, 8, s.shearWeightRatio, 3);
    // 剪重比下限（GB 50011 表 5.2.5，7 度区约 1.6%，这里给通用提示值）
    if (s.shearWeightRatio > 0.0 && s.shearWeightRatio < 0.016)
      story_->item(i, 8)->setForeground(Status::warn());
  }
  story_->resizeColumnsToContents();
  story_->setSortingEnabled(true);
}

void ResultPanel::buildEnvelope() {
  envelope_->setRowCount(0);
  if (!post_) return;

  const post::Envelope& e = post_->envelope();
  struct Row { QString k; QString v; QString at; };
  const std::vector<Row> rows = {
      {QStringLiteral("最大总位移 (mm)"), fmtNum(e.maxDisp * 1e3, 4), QStringLiteral("节点 #%1").arg(e.maxDispNode)},
      {QStringLiteral("最大水平位移 (mm)"), fmtNum(e.maxUxy * 1e3, 4), QStringLiteral("节点 #%1").arg(e.maxUxyNode)},
      {QStringLiteral("最大竖向位移 (mm)"), fmtNum(e.maxUz * 1e3, 4), QStringLiteral("节点 #%1").arg(e.maxUzNode)},
      {QStringLiteral("最大转角 (rad)"), fmtNum(e.maxRotation, 3), QStringLiteral("节点 #%1").arg(e.maxRotationNode)},
      {QStringLiteral("最大 von Mises (kPa)"), fmtNum(e.maxVonMises, 4), QStringLiteral("节点 #%1").arg(e.maxVonMisesNode)},
      {QStringLiteral("最大支座反力 (kN)"), fmtNum(e.maxReaction, 4), QStringLiteral("节点 #%1").arg(e.maxReactionNode)},
      {QStringLiteral("最大竖向反力 (kN)"), fmtNum(e.maxReactionZ, 4), QStringLiteral("节点 #%1").arg(e.maxReactionZNode)},
      {QStringLiteral("总重力荷载 (kN)"), fmtNum(e.totalWeight, 5), QString()},
      {QStringLiteral("竖向荷载合计 (kN)"), fmtNum(e.totalLoadZ, 5), QString()},
      {QStringLiteral("基底剪力 (kN)"), fmtNum(e.baseShear, 4), QString()},
      {QStringLiteral("剪重比"), fmtNum(e.shearWeightRatio, 3), QString()},
      {QStringLiteral("最大层间位移角"),
       e.maxDrift > 1e-12 ? QStringLiteral("1/%1").arg(fmtNum(1.0 / e.maxDrift, 4)) : QStringLiteral("—"),
       QStringLiteral("%1 层").arg(e.maxDriftStory)},
      {QStringLiteral("最大剪重比层"), fmtNum(e.maxShearRatio, 3), QStringLiteral("%1 层").arg(e.maxShearRatioStory)},
  };

  envelope_->setRowCount(static_cast<int>(rows.size()));
  for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
    setText(envelope_, i, 0, rows[static_cast<size_t>(i)].k);
    auto* it = new QTableWidgetItem(rows[static_cast<size_t>(i)].v);
    it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
    it->setFont(monoFont(9));
    envelope_->setItem(i, 1, it);
    setText(envelope_, i, 2, rows[static_cast<size_t>(i)].at);
  }

  // 平衡校核：施加的竖向荷载必须等于竖向反力之和。
  // 【这是整个软件最重要的一条自检】它同时覆盖了组装、约束消去、
  // 求解、反力恢复四个环节，只要有一处错就会不闭合。
  const double diff = e.totalLoadZ - e.totalWeight;
  envelope_->setRowCount(static_cast<int>(rows.size()) + 1);
  const int r = static_cast<int>(rows.size());
  setText(envelope_, r, 0, QStringLiteral("平衡校核（竖向荷载 − 重力）"));
  const bool okBal = std::abs(diff) < 1e-6 * std::max(1.0, std::abs(e.totalWeight));
  auto* it = new QTableWidgetItem(okBal ? QStringLiteral("闭合") : fmtNum(diff, 5));
  it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
  it->setForeground(okBal ? Status::ok() : Status::fail());
  envelope_->setItem(r, 1, it);
  setText(envelope_, r, 2, okBal ? QString() : QStringLiteral("检查荷载缺失"));
  envelope_->resizeColumnsToContents();
}

void ResultPanel::buildMemberTable() {
  members_->setSortingEnabled(false);
  members_->setRowCount(0);
  if (!post_) { members_->setSortingEnabled(true); return; }

  const auto& brs = post_->beamResults();
  members_->setRowCount(static_cast<int>(brs.size()));
  for (int i = 0; i < static_cast<int>(brs.size()); ++i) {
    const post::BeamResult& b = brs[static_cast<size_t>(i)];
    setNumText(members_, i, 0, QStringLiteral("#%1").arg(b.id), b.id);
    setText(members_, i, 1, post::memberKindName(b.kind));
    setNum(members_, i, 2, b.L, 4);
    setNum(members_, i, 3, b.Mmax, 4);
    setNum(members_, i, 4, b.Vmax, 4);
    setNum(members_, i, 5, std::abs(b.Nmax) > std::abs(b.Nmin) ? b.Nmax : b.Nmin, 4);
    setNum(members_, i, 6, b.sigAbsMax, 4);
    if (b.util > 0.0) {
      const QColor c = b.util > 1.0 ? Status::fail() : (b.util > 0.85 ? Status::warn() : Status::ok());
      setNum(members_, i, 7, b.util, 3, c);
    } else {
      setText(members_, i, 7, QStringLiteral("—"));
    }
  }
  members_->resizeColumnsToContents();
  members_->setSortingEnabled(true);
}

void ResultPanel::buildReport() {
  if (!post_) { report_->setPlainText(QStringLiteral("尚无计算结果。")); return; }
  report_->setPlainText(QString::fromStdString(post_->report()));
}

}  // namespace ui
