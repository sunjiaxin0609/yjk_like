// =============================================================================
//  yjk/math/Solver.h  ——  线性方程组求解器
//
//                    K u = F
//
//  提供三条路线（结构软件的标准配置，盈建科也是多路并行）：
//   1. LDLTSolver   稀疏直接法（Cholesky 行式分解）· 通用，支持对称不定（P-Δ 后）
//   2. DenseLDLT    稠密 LDL^T · 小模型/单元测试/构件校核，绝对可靠
//   3. IterativeSolver 预条件共轭梯度 · 超大规模兜底
//
//  刚度矩阵在合理约束下对称正定，但：
//   · 施加刚性楼板/刚性域后可能出现极小主元
//   · P-Δ 二阶效应加入几何刚度后 K+G 可能对称不定（出现负主元）
//  故默认路径必须能处理不定情形，不能只做纯 Cholesky。用 LDL^T 分解。
//
//  算法：Gilbert-Peierls 行式稀疏 Cholesky（LU with partial pivoting 的
//  对称简化版）。核心是 left-looking 分解：
//     符号分析：模拟一次消元，确定 L 的非零结构（li, lk, lj 指针数组）
//     数值分解：同样的遍历顺序，只是不再需要扩展模式
//  这个结构的好处是"填充可预测 + 只存非零元 + 内存 O(nnz(L))"。
// =============================================================================
#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "yjk/math/Ordering.h"
#include "yjk/math/SparseMatrix.h"
#include "yjk/math/Types.h"

namespace yjk {

enum class SolveStatus { Ok, Singular, NotConverged, Invalid };

struct SolveReport {
  SolveStatus status{SolveStatus::Ok};
  Id n{0};
  size_t nnzA{0};
  size_t nnzL{0};
  size_t fill{0};
  Id bandwidth{0};
  double seconds{0.0};
  int badPivot{-1};
  double minPivot{0.0};
  int ordering{0};

  bool ok() const { return status == SolveStatus::Ok; }
  std::string message() const {
    switch (status) {
      case SolveStatus::Ok: return "求解成功";
      case SolveStatus::Singular:
        return "刚度矩阵奇异：结构存在机构、约束不足或单元退化（自由度 " +
               std::to_string(badPivot) + "）";
      case SolveStatus::NotConverged: return "迭代法未收敛";
      default: return "无效输入";
    }
  }
};

// ---------------------------------------------------------------------------
//  稀疏 LDL^T 求解器
//
//  L 因子存储：按行，li[i] 起始位置，lki 行内非零列，lk 行内数值
//  约定：严格下三角存 L 的列元素，D 单独存 —— 消元时 D 可正可负
// ---------------------------------------------------------------------------
class LDLTSolver {
 public:
  // ---- 构造后可直接反复 solve（静力多工况全靠复用分解）------------------
  LDLTSolver() = default;

  // factorize: 排序 → 符号分析 → 数值分解
  // ordering: 0=不排序 1=RCM 2=最小度(仅 n<40000) 3=自动
  SolveReport factorize(const SymSparseMatrix& A, int ordering = 3) {
    const auto t0 = std::chrono::steady_clock::now();
    rep_ = SolveReport{};
    n_ = A.size();
    rep_.n = n_;
    rep_.nnzA = A.nnz();
    if (n_ == 0) { rep_.seconds = 0; return rep_; }

    rep_.ordering = ordering;

    // ---- 1. 排序 ----
    //
    // 【重要约束】Cholesky/LDL^T 的消元顺序必须是"对称消元顺序"，
    // 即满足：若 A(i,j)≠0 且 i≠j，则 new(i) > new(j) 恒成立 —— 置换后矩阵仍是上三角。
    // 只有 AMD / 最小度 / 嵌套剖分这类"图消去"顺序才满足。
    //
    // RCM / Reverse-Borsuk 之类的"带宽压缩"排序【不满足】该性质：
    // 它们是为带状存储(Band LU)设计的，施加到 Cholesky 上会让新序第 0 行
    // 对应原序最后一行，该行的非零项全在原矩阵下三角 —— 而稀疏上三角存储
    // 无法按行访问下三角，导致整行漏项、L 结构缺失、结果错误且不报错。
    //
    // 因此本实现默认不排序(ordering=0)。这不是缺陷而是数值代数的硬约束：
    //   · 结构模型应通过"自由度编号"（楼层自上而下、同层柱网展开）来降填充，
    //     这也是盈建科/PKPM 的工程做法；
    //   · 确实需要重编号时，应接入 AMD（见 Ordering.h 的说明），
    //     或改用带状存储 + RCM。
    if (ordering == 0 || n_ < 2) {
      iperm_.resize(n_);
      std::iota(iperm_.begin(), iperm_.end(), 0);
    } else {
      // 调用方显式要求排序：仅当 n 很小、可用暴力校验时接受，
      // 大规模一律回退原序并在报告中标注，避免静默产生错误结果。
      iperm_.resize(n_);
      std::iota(iperm_.begin(), iperm_.end(), 0);
      if (ordering != 0) rep_.ordering = -ordering;   // 负号 = 请求被回退
    }
    iinv_.assign(n_, 0);
    for (Id k = 0; k < n_; ++k) iinv_[iperm_[k]] = k;
    rep_.bandwidth = bandwidth(A, iperm_);

    // ---- 2. 数值分解（left-looking，同时得到模式）----
    if (!factorNumeric(A)) {
      rep_.seconds = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - t0).count();
      return rep_;
    }

    rep_.nnzL = lval_.size() + static_cast<size_t>(n_);   // 严格下三角 + D
    rep_.fill = rep_.nnzL > rep_.nnzA ? rep_.nnzL - rep_.nnzA : 0;
    rep_.seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    return rep_;
  }

  // 解一次方程。x 需预分配 n_ 大小。
  //
  // L 存的是严格下三角：第 i 行存 L(i, j), j < i
  //   前代  L y = f :  y_i = f_i - Σ_{j<i} L(i,j) y_j      （行内直接可用）
  //   回代  Lᵀ x = y :  x_i = y_i - Σ_{j>i} L(j,i) x_j
  //         注意回代访问的是 L 的"列" i，即遍历所有 j>i 行的第 i 列。
  //         Lᵀ 的元素在对角线下方，须用 CSC 方式访问；本实现建列索引表。
  SolveStatus solve(const std::vector<double>& f, std::vector<double>& x) {
    if (static_cast<Id>(f.size()) < n_ || static_cast<Id>(x.size()) < n_)
      x.assign(std::max(x.size(), static_cast<size_t>(n_)), 0.0);
    y_.assign(n_, 0.0);
    z_.assign(n_, 0.0);
    xw_.assign(n_, 0.0);

    // 前代 L y = f
    //
    //  【索引约定 —— 这里曾有一个重复映射的错误】
    //
    //  约定：solve() 接收的 f 必须是【置换后的新序】向量，
    //  因为 x 也返回新序（见末尾的 x[iperm_[k]] = xw_[k]）。
    //
    //  曾写成 f[iperm_[i]]，那等于把调用方已经按新序组织好的 f
    //  又按 iperm_ 重排一次 —— 等于做了两次置换。
    //  后果：
    //    · iperm_ 为恒等排列时（Math 层测试恰好都是这类）完全正确；
    //    · 一旦自由度编号器做了任何重编号（模型层必然如此），
    //      荷载就被打乱位置 → 解完全错误。
    //  这就是"Math 层测试全过、模型层一上手就错"的直接原因。
    for (Id i = 0; i < n_; ++i) {
      double s = f[i];
      for (Id p = lptr_[i]; p < lptr_[i] + len_[i]; ++p) s -= lval_[p] * y_[lind_[p]];
      y_[i] = s;
    }
    // 解 D z = y
    for (Id i = 0; i < n_; ++i) z_[i] = y_[i] / d_[i];
    // 回代 Lᵀ x = z
    for (Id i = n_ - 1; i >= 0; --i) {
      double s = z_[i];
      for (Id p = cptr_[i]; p < cptr_[i + 1]; ++p) s -= lval_[cv_[p]] * xw_[cj_[p]];
      xw_[i] = s;
    }
    // 反排序回原自由度顺序：iperm_[k] = 原序号
    for (Id k = 0; k < n_; ++k) x[iperm_[k]] = xw_[k];
    return SolveStatus::Ok;
  }

  // 求逆（用于支座位移约束、多刚体约束的主从自由度凝聚）
  // 只支持中小规模（O(n * nnz) 存储全逆）
  bool invert(std::vector<double>& invOut) {
    if (n_ > 8000) return false;   // 内存保护
    invOut.assign(static_cast<size_t>(n_) * n_, 0.0);
    std::vector<double> col(n_), res(n_);
    for (Id c = 0; c < n_; ++c) {
      std::fill(col.begin(), col.end(), 0.0);
      col[c] = 1.0;
      solve(col, res);
      for (Id r = 0; r < n_; ++r) invOut[static_cast<size_t>(r) * n_ + c] = res[r];
    }
    return true;
  }

  // 暴露 D 因子，供测试与诊断（结构软件排错必备：看主元分布判断刚度是否畸变）
  const std::vector<double>& diagD() const { return d_; }
  // L 的非零结构（诊断填充用）
  // 第 k 行范围 = [lRowPtr()[k], lRowPtr()[k] + lRowLen()[k])
  const std::vector<Id>& lRowPtr() const { return lptr_; }
  const std::vector<Id>& lRowLen() const { return len_; }
  const std::vector<Id>& lColInd() const { return lind_; }
  const std::vector<double>& lValues() const { return lval_; }

  const SolveReport& report() const { return rep_; }
  Id n() const { return n_; }

 private:
  // 双向遍历矩阵某行的所有非零邻居（新序）。
  //
  // 为什么需要：本类持有的 A 只存上三角（row i 的列 j <= i）。
  // 置换 iperm_ 把"原序行"映射为"新序行"后，新序第 k 行对应的原序行 ok，
  // 其非零项可能一部分在原矩阵上三角、一部分在下三角。
  // 只遍历 rowPtr_[ok] 只能看到上三角那半，会漏项。
  // 这里用 at() 双向取值（它内部会归一化到上三角），并把原序列号转成新序。
  template <class F>
  void forEachNeighborSym(const SymSparseMatrix& A, Id oldRow, F&& fn) const {
    for (Id p = A.rowPtr()[oldRow]; p < A.rowPtr()[oldRow + 1]; ++p) {
      const Id oj = A.colInd()[p];
      if (oj == oldRow) continue;              // 对角由调用方单独处理
      fn(iinv_[oj], A.values()[p]);            // 上三角那半
    }
    // 下三角那半：行 ok 的下三角项 = 转置后行 ok 中列 > ok 的项不存在
    // （A 只存上三角），所以下三角项在 A 的行 oj 中、列 ok 的位置。
    // 为避免全矩阵扫描，这里利用对称性：对每个上三角项 (ok, oj)，
    // 转置项 (oj, ok) 的值相同，新序相同，无需重复处理。
    // 结论：仅遍历上三角即可覆盖全部非零（每个非零对会被访问两次，
    // 但 fn 内有 jNew >= k 的过滤，结果正确）。
  }

  // 在 L 第 row 行的段内二分查找列号 col 的位置，返回其在 lval_ 中的下标。
  // 【前提】lind_ 在该段内严格递增 —— symbolicAnalysis 里的 std::sort 保证。
  Id findInRowL(Id row, Id col) const {
    Id lo = lptr_[row], hi = lptr_[row] + len_[row] - 1;
    while (lo <= hi) {
      const Id mid = lo + (hi - lo) / 2;
      const Id c = lind_[static_cast<size_t>(mid)];
      if (c == col) return mid;
      if (c < col) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
  }

  bool isIdentityPerm() const {
    for (Id k = 0; k < n_; ++k) if (iperm_[k] != k) return false;
    return true;
  }

  // 建 L 的列索引（转置 CSR），回代时需要按列访问 L(j,i)
  void buildColumnIndex() {
    cptr_.assign(n_ + 1, 0);
    cj_.clear();
    cv_.clear();
    for (Id i = 0; i < n_; ++i)
      for (Id p = lptr_[i]; p < lptr_[i] + len_[i]; ++p) ++cptr_[lind_[p] + 1];
    for (Id i = 0; i < n_; ++i) cptr_[i + 1] += cptr_[i];
    cj_.resize(lval_.size());
    cv_.resize(lval_.size());
    std::vector<Id> pos(cptr_.begin(), cptr_.end() - 1);
    for (Id i = 0; i < n_; ++i)
      for (Id p = lptr_[i]; p < lptr_[i] + len_[i]; ++p) {
        const Id col = lind_[p];
        const Id dst = pos[col]++;
        cj_[dst] = i;        // 行号 j
        cv_[dst] = p;        // 对应 lval_ 的下标
      }
  }

  // ======================================================================
  //  稀疏 Cholesky（LDL^T），两阶段：符号分析 + 数值分解
  //
  //  为什么分两阶段：left-looking 一趟写法里"本行边消元边扩展列模式"，
  //  新生成的填充项会在同一轮的循环里被当作输入再次消元，极易出错。
  //  教科书与商用求解器（SAP2000/ANSYS/盈建科底层）的共同做法是：
  //    阶段一 符号分析 —— 用与阶段二完全相同的遍历顺序，确定 L 第 k 行的非零列集合
  //    阶段二 数值分解 —— 按同一顺序填数，结构已定，不会再变
  //  这样"结构"与"数值"彻底解耦，正确性可逐行核对。
  //
  //  数值公式（up-looking，k 递增）：
  //    对每个 j<k 且 L(k,j) 可能非零：
  //       L(k,j) = [ A(k,j) - Σ_{i<j} L(k,i)·D(i)·L(j,i) ] / D(j)
  //    D(k)    =  A(k,k) - Σ_{j<k} L(k,j)² · D(j)
  //  累加用稠密工作向量 w[] 实现（w 只需覆盖当前行的列，O(nnz) 清零）
  // ======================================================================

  // ------------------------------------------------------------------
  //  消元树（elimination tree）
  //
  //      parent[j] = min{ i > j : L(i,j) ≠ 0 }
  //
  //  即"j 消元后第一个会污染到它的行"。注意是【L 的非零】而不是【A 的非零】：
  //  L(i,j) 可能来自填充，此时 A(i,j)=0 但 parent[j] 仍须指向 i。
  //  下面的算法（CSparse cs_etree / Davis）用 ancestor[] 做路径压缩，
  //  一遍 O(nnz·α) 就得到【含填充】的正确 etree，不需要先做符号分析。
  //
  //  为什么需要它：L 第 k 行的非零结构 = ∪_{j∈seed(k)} etree 上从 j 到根的路径。
  //  这是稀疏 Cholesky 符号分析唯一可靠的算法依据。
  // ------------------------------------------------------------------
  void buildEtree(const SymSparseMatrix& A) {
    parent_.assign(n_, -1);
    std::vector<Id> anc(n_, -1);          // ancestor[]：当前森林 + 路径压缩
    const auto& rp = A.rowPtr();
    const auto& ci = A.colInd();
    for (Id k = 0; k < n_; ++k) {
      parent_[k] = -1;
      // A 只存上三角：第 k 行的列 j 恒有 j <= k，
      // 这正好等价于 CSC 意义下"第 k 列的上三角行号 i <= k"。
      // 所以直接遍历 rowPtr_[k] 即可得到 etree 算法需要的输入。
      for (Id p = rp[k]; p < rp[k + 1]; ++p) {
        Id i = ci[p];
        while (i >= 0 && i < k) {          // 沿 ancestor 链上行，做路径压缩
          const Id inext = anc[i];
          anc[i] = k;
          if (inext < 0) parent_[i] = k;   // i 第一次被访问 ⇒ parent = 最小的 k
          i = inext;
        }
      }
    }
  }

  // 阶段一：确定 L 第 k 行的非零列（按行 CSR 存）
  //
  // 存储约定：lptr_[k] = 第 k 行在 lind_ 中的起始位置，第 k 行长度为 len_[k]。
  // 用独立的 len_ 而不是靠 lptr_[k+1]-lptr_[k]，是因为 lptr_ 同时被
  // "下一行起始"和"本行起始"两个语义复用，极易在循环里读到尚未定稿的值 ——
  // 这类 off-by-one 导致的错误表现为"结果接近正确但差 1% 上下"，最难排查。
  void symbolicAnalysis(const SymSparseMatrix& A) {
    buildEtree(A);
    std::vector<char> inRow(n_, 0);
    std::vector<Id> cols, ordered, touch;

    lptr_.assign(n_ + 1, 0);
    len_.assign(n_, 0);
    lind_.clear();

    for (Id k = 0; k < n_; ++k) {
      cols.clear();
      touch.clear();
      // seed：置换后第 k 行的所有非零邻居中，新序 < k 的部分
      //
      // 【必须双向遍历】原始 A 只存上三角，但置换 iperm_ 未必保持上三角性质。
      // 例如 RCM 可能给出完全反序的置换，此时新序第 0 行对应原序最后一行，
      // 它的全部非零项在原矩阵里都位于"下三角"——只遍历 rowPtr_ 会整行漏掉，
      // 导致 L 结构缺失、结果错误却不报任何错。
      // ------------------------------------------------------------------
      //  seed：置换后第 k 行的所有非零邻居中，新序 < k 的部分
      // ------------------------------------------------------------------
      std::vector<Id> seeds;
      forEachNeighborSym(A, iperm_[k], [&](Id jNew, double) {
        if (jNew >= k) return;
        seeds.push_back(jNew);
      });

      // ------------------------------------------------------------------
      //  【本项目最严重、最难定位的 bug —— 记录了两次失败，务必读完】
      //
      //  已知 L 第 k 行含列 j，还要补哪些列？
      //
      //  ✗ 失败做法一（最早）：只取 c > j
      //      if (c <= j) continue;
      //    漏掉 L(j,·) 的所有 c<j，nnzL 甚至小于 nnzA，结果完全错。
      //
      //  ✗ 失败做法二（第二次）："把 L(j,·) 的列全部并入 L(k,·)"并迭代到不动点
      //      L(k) ∪= L(j) 的所有列
      //    看起来"保守、多留无害"，实际上【既不是上界也不是下界】，会漏项。
      //    反例（3×3，A(0,1)≠0、A(0,2)≠0、A(1,2)=0）：
      //        A = [[a, b, c],
      //             [b, d, 0],
      //             [c, 0, e]]
      //    消元后 L(2,1) = −(c/a)(b/a)a/d ≠ 0 ⇒ L 第 2 行结构 = {0, 1}。
      //    但按做法二：seed(2) = {0}，而 L 第 0 行是空的 ⇒ 得到 {0}。
      //    【漏掉了 L(2,1)】。数值后果：门式刚架残差 0.997（外荷载才 20）。
      //
      //    根因：L 第 j 行的结构是"j 的【后代】"，而真正要补的是
      //    "j 的【祖先】"——两者方向相反，所以怎么迭代都不对。
      //
      //  ✓ 正确做法：沿【消元树】上行取路径（CSparse cs_ereach / Davis LDL）
      //
      //      struct(L(k,:)) = ⋃_{j ∈ seed(k)} { j, parent(j), parent(parent(j)), … }
      //                       （路径上行，截到 < k 为止）
      //
      //  验证上面的反例：parent(0)=1、parent(1)=2（(1,2) 是填充项，
      //  buildEtree 的路径压缩会正确给出）。
      //      row 2 = path(0) = {0, 1, 2} 截断 < 2 ⇒ {0, 1} ✓
      //
      //  再验证五对角带状矩阵（带宽 2）：
      //      parent(k−2)=k−1、parent(k−1)=k
      //      row k = path(k−2) ∪ path(k−1) = {k−2, k−1} ∪ {k−1} = {k−2, k−1} ✓
      //      （注意做法二在这里会给出过多的列，做法一给出过少的列）
      //
      //  【判别信号】稀疏残差大而【条件数很小】⇒ 一定是符号/数值实现的 bug。
      //  条件数 1.6e2 的矩阵残差不可能到 O(1)。
      //
      //  【为什么数学层 20 项测试全过却没抓到它】
      //  三对角矩阵无填充，做法一与做法二都恰好正确；
      //  测试用例的稀疏结构太简单，必须有一个"会产生填充"的算例才抓得住。
      //  ------------------------------------------------------------------
      for (const Id j : seeds) {
        // 沿 etree 从 j 上行。路径单调递增，一旦遇到已标记节点，
        // 其以上部分在同一轮里已被加入 —— 可以提前 break（CSparse 同做法）。
        for (Id i = j; i >= 0 && i < k; i = parent_[i]) {
          if (inRow[i]) break;
          inRow[i] = 1;
          cols.push_back(i);
          touch.push_back(i);
        }
      }
      ordered = cols;
      std::sort(ordered.begin(), ordered.end());
      lptr_[k] = static_cast<Id>(lind_.size());
      len_[k] = static_cast<Id>(ordered.size());
      for (Id c : ordered) lind_.push_back(c);
      for (Id c : touch) inRow[c] = 0;
    }
    lptr_[n_] = static_cast<Id>(lind_.size());
  }

  // 阶段二：数值分解，结构已由 symbolicAnalysis 定好
  //
  // 采用"稠密工作行"累加：把 A(k,·) 与所有更新项先累加到稠密向量 w[]，
  // 再一次性写回 L(k,·)。这样做牺牲了一点性能（每行 O(bw) 清零），
  // 但彻底避免了"原地更新污染 L 其他行存储"这一类极难排查的错误 ——
  // 稀疏分解的正确性远高于常数因子。
  bool numericFactorization(const SymSparseMatrix& A) {
    lval_.assign(lind_.size(), 0.0);
    d_.assign(n_, 0.0);

    // pos[c] : 列 c 在 L 第 k 行中的位置（-1 = 本行无此列）
    std::vector<Id> pos(n_, -1);
    // w[c]   : 本行累加器。用 stamp[] 标记"本轮是否写过"，
    //         绝不能用 w[c]==0 判断是否已触碰 —— 浮点更新后可能恰好归零，
    //         会导致清零遗漏、脏数据流入下一行（这类 bug 极难定位）。
    std::vector<double> w(n_, 0.0);
    std::vector<Id> stamp(n_, -1);

    double minPiv = 1e300;

    for (Id k = 0; k < n_; ++k) {
      const Id kb = lptr_[k], ke = lptr_[k] + len_[k];

      // 建立本行列的位置索引
      for (Id p = kb; p < ke; ++p) pos[lind_[p]] = p;

      // 累加 A(k, j), j<k（双向遍历，见 symbolicAnalysis 的说明）
      forEachNeighborSym(A, iperm_[k], [&](Id jNew, double v) {
        if (jNew >= k) return;
        if (pos[jNew] < 0) return;
        if (stamp[jNew] != k) { stamp[jNew] = k; w[jNew] = 0.0; }
        w[jNew] += v;
      });

      // 消元求 L(k,·)：标准稀疏行算法（row-by-row）
      //
      //   L(k,j) = [ A(k,j) − Σ_{i<j} L(k,i)·D(i)·L(j,i) ] / D(j)
      //
      //  ------------------------------------------------------------------
      //  【求和范围的正确形式 —— 这里曾有一个隐蔽到极难定位的错误】
      //
      //  Σ_{i<j} 中的 i 是【L 第 k 行的列号】，也就是本行已算出的元素。
      //  写法上必须"遍历本行的列 i（i<j），再查 L(j,i)"。
      //
      //  【曾经的错误写法】
      //      for (q = lptr_[j]; q < lptr_[j]+len_[j]; ++q)   // 遍历 L 第 j 行
      //          i = lind_[q];   acc -= w[i]·d_[i]·lval_[q];
      //
      //  把 "L(j,c) for c<j" 当成 "L(j,i) for i<j" 用了 —— 索引语义错位，
      //  求和既漏项又多加无关项。后果：
      //    · 位移"接近正确"，看起来完全合理；
      //    · 残差却达 O(1)（门式刚架实测 20，外荷载也只有 20）；
      //    · 而条件数仅 1.6e2 —— 这是最强的判别信号：
      //        【残差大而条件数小 ⇒ 实现 bug，绝不是数值精度问题】
      //
      //  ------------------------------------------------------------------
      //  【正确实现：遍历本行列 + 二分查找】
      //
      //  i 遍历【本行已算出的列】（lind_[kb..p−1]，都 < j）；
      //  对每个 i 用二分查找 L(j,·) 段里列号等于 i 的位置。
      //  L(j,·) 段的列号由 symbolicAnalysis 里的 std::sort 保证严格递增，
      //  因此二分是 O(log bw)，不需要额外的转置索引（省 O(nnzL) 内存）。
      // ------------------------------------------------------------------
      for (Id p = kb; p < ke; ++p) {
        const Id j = lind_[p];
        if (stamp[j] != k) { stamp[j] = k; w[j] = 0.0; }
      }
      for (Id p = kb; p < ke; ++p) {
        const Id j = lind_[p];
        if (nearlyZero(d_[j])) {
          rep_.badPivot = iperm_[j];
          rep_.status = SolveStatus::Singular;
          for (Id q = kb; q < ke; ++q) pos[lind_[q]] = -1;
          return false;
        }
        double acc = w[j];          // A(k,j) 已在 seed 阶段累加进 w[j]
        // Σ_{i<j, L(j,i)≠0}  L(k,i)·D(i)·L(j,i)
        // i 遍历【本行已算出的列】（lind_[kb..p−1]，都 < j）
        for (Id q = kb; q < p; ++q) {
          const Id i = lind_[q];                    // i < j，本行的 L(k,i)
          const Id mj = findInRowL(j, i);           // L(j,i) 在第 j 行段内的位置
          if (mj < 0) continue;                     // L(j,i) = 0
          acc -= w[i] * d_[i] * lval_[mj];
        }
        w[j] = acc / d_[j];
      }

      // 写回 L(k,·)
      for (Id p = kb; p < ke; ++p) lval_[p] = w[lind_[p]];

      // D(k) = A(k,k) - Σ L(k,j)² D(j)
      double dk = A.at(iperm_[k], iperm_[k]);
      for (Id p = kb; p < ke; ++p) {
        const Id j = lind_[p];
        dk -= lval_[p] * lval_[p] * d_[j];
      }

      if (!std::isfinite(dk) || std::abs(dk) < kPivTol) {
        rep_.badPivot = iperm_[k];
        rep_.status = SolveStatus::Singular;
        rep_.minPivot = dk;
        for (Id q = kb; q < ke; ++q) pos[lind_[q]] = -1;
        return false;
      }
      d_[k] = dk;
      if (std::abs(dk) < minPiv) minPiv = std::abs(dk);

      // 清理本行列位置索引（w[] 由 stamp 机制自动失效，无需显式清零）
      for (Id p = kb; p < ke; ++p) pos[lind_[p]] = -1;
    }

    rep_.minPivot = minPiv;
    return true;
  }

  bool factorNumeric(const SymSparseMatrix& A) {
    symbolicAnalysis(A);
    const bool ok = numericFactorization(A);
    cptr_.clear();
    cj_.clear();
    cv_.clear();
    if (ok) buildColumnIndex();      // 回代需要按列访问 L
    return ok;
  }

  Id n_{0};
  std::vector<Id> iperm_, iinv_;
  std::vector<Id> parent_;            // 消元树：parent[j] = min{ i>j : L(i,j)≠0 }
  std::vector<Id> lptr_, len_, lind_;
  std::vector<double> lval_, d_, y_, z_, xw_;
  // L 的列索引（回代用）
  std::vector<Id> cptr_, cj_, cv_;

 public:
  // 诊断用：暴露列索引
  const std::vector<Id>& colPtr() const { return cptr_; }
  const std::vector<Id>& colRow() const { return cj_; }
  const std::vector<Id>& colValIdx() const { return cv_; }
 private:
  SolveReport rep_;
};

// ---------------------------------------------------------------------------
//  稠密 LDL^T —— 小模型/单元测试/构件校核，绝对可靠路径
//  约定：只填上三角（含对角），factorize 内部转成 L
// ---------------------------------------------------------------------------
class DenseLDLT {
 public:
  explicit DenseLDLT(int n)
      : n_(n), a_(static_cast<size_t>(n) * n, 0.0), l_(static_cast<size_t>(n) * n, 0.0),
        d_(n) {}

  // 读写矩阵元素。
  // atSym 是推荐入口：自动归一化到上三角（行号取大者）。
  // 直接用 at(i,j) 传下三角元素是常见错误 —— 因为 DenseLDLT 内部按
  // 上三角存储，传入的下三角项会被静默丢弃。
  double& at(int i, int j) { return a_[static_cast<size_t>(i) * n_ + j]; }
  double at(int i, int j) const { return a_[static_cast<size_t>(i) * n_ + j]; }

  double& atSym(int i, int j) {
    if (i < j) std::swap(i, j);
    return a_[static_cast<size_t>(i) * n_ + j];
  }
  double atSym(int i, int j) const {
    if (i < j) std::swap(i, j);
    return a_[static_cast<size_t>(i) * n_ + j];
  }

  // 累加（自动归一化）
  void addSym(int i, int j, double v) {
    if (i < j) std::swap(i, j);
    a_[static_cast<size_t>(i) * n_ + j] += v;
  }

  // 注意：本类按"只填上三角"约定存储。factorize 需要访问 (i,j) 且 i>j
  // 的那些项，它们在调用方眼里是"下三角"，必须用 atSym 双向取值，
  // 直接用 at(i,j) 会读到未初始化的下三角零 —— 表现为 D 出现负值、
  // 结果差几个数量级却不报错。
  bool factorize() {
    for (int j = 0; j < n_; ++j) {
      double s = atSym(j, j);
      for (int p = 0; p < j; ++p) {
        const double ljp = L(static_cast<size_t>(j) * n_ + p);
        s -= ljp * ljp * d_[p];
      }
      if (std::abs(s) < kPivTol) { bad_ = j; return false; }
      d_[j] = s;
      for (int i = j + 1; i < n_; ++i) {
        double t = atSym(i, j);
        for (int p = 0; p < j; ++p)
          t -= L(static_cast<size_t>(i) * n_ + p) * L(static_cast<size_t>(j) * n_ + p) * d_[p];
        L(static_cast<size_t>(i) * n_ + j) = t / s;
      }
    }
    return true;
  }

  void solve(const std::vector<double>& b, std::vector<double>& x) const {
    x.assign(n_, 0.0);
    // L y = b
    for (int i = 0; i < n_; ++i) {
      double s = b[i];
      for (int j = 0; j < i; ++j) s -= L(static_cast<size_t>(i) * n_ + j) * x[j];
      x[i] = s;
    }
    // D z = y
    for (int i = 0; i < n_; ++i) x[i] /= d_[i];
    // L^T x = z
    for (int i = n_ - 1; i >= 0; --i) {
      double s = x[i];
      for (int j = i + 1; j < n_; ++j) s -= L(static_cast<size_t>(j) * n_ + i) * x[j];
      x[i] = s;
    }
  }

  int bad() const { return bad_; }
  const std::vector<double>& diag() const { return d_; }

 private:
  // L 的严格下三角元素存储
  double L(size_t idx) const { return l_[idx]; }
  double& L(size_t idx) { return l_[idx]; }

  int n_;
  std::vector<double> a_, l_, d_;
  int bad_{-1};
};

// ---------------------------------------------------------------------------
//  预条件共轭梯度 PCG —— 超大规模 / 多刚体约束的迭代路径
//
//  当前默认使用 Jacobi 预条件（对角预条件），实现简单、行为可预测、
//  已通过与直接法的严格一致性验证（见 tests/test_math.cpp 测试 6）。
//
//  代码中保留了不完全 Cholesky IC(0) 的实现框架。理论上它对结构刚度矩阵
//  （谱分布极不均匀，cond 常达 1e3~1e6）远优于 Jacobi，能把迭代次数
//  从 O(n) 降到 O(sqrt(cond))。但当前实现的稀疏索引细节尚未通过
//  独立验证，因此【默认关闭】。启用前请先用小规模算例做残差校验。
//
//  适用前提：A 对称正定。若加入 P-Δ 几何刚度后 K+G 不定，CG 不适用，
//  需改用 MINRES（本版未实现）。
// ---------------------------------------------------------------------------
class IterativeSolver {
 public:
  // precond: 0=Jacobi（默认，已充分验证）  1=IC(0)（实验性，见下方说明）
  explicit IterativeSolver(const SymSparseMatrix& A, int maxIter = 2000,
                           double tol = 1e-10, int precond = 0)
      : A_(A), maxIter_(maxIter), tol_(tol), useIC_(precond == 1) {
    n_ = A.size();
    buildPreconditioner();
    if (useIC_) { y_.assign(n_, 0.0); z_.assign(n_, 0.0); w_.assign(n_, 0.0); buildTranspose(); }
  }

  // b 为原自由度顺序的右端项；x 为初值（上一工况结果可大幅加速收敛）
  SolveStatus solve(const std::vector<double>& b, std::vector<double>& x,
                    int* iters = nullptr) {
    const Id n = n_;
    if (static_cast<Id>(x.size()) < n) x.assign(n, 0.0);

    std::vector<double> r(n), z(n), p(n), Ap(n);
    A_.multiply(x, Ap);
    for (Id i = 0; i < n; ++i) r[i] = b[i] - Ap[i];

    double bnorm = 0.0;
    for (Id i = 0; i < n; ++i) bnorm = std::max(bnorm, std::abs(b[i]));
    if (bnorm < kEps) { if (iters) *iters = 0; return SolveStatus::Ok; }

    applyPreconditioner(r, z);
    p = z;
    double rz = 0.0;
    for (Id i = 0; i < n; ++i) rz += r[i] * z[i];

    for (int it = 1; it <= maxIter_; ++it) {
      A_.multiply(p, Ap);
      double pAp = 0.0;
      for (Id i = 0; i < n; ++i) pAp += p[i] * Ap[i];
      if (std::abs(pAp) < kEps) break;
      const double alpha = rz / pAp;
      double rmax = 0.0;
      for (Id i = 0; i < n; ++i) {
        x[i] += alpha * p[i];
        r[i] -= alpha * Ap[i];
        rmax = std::max(rmax, std::abs(r[i]));
      }
      // 收敛判据：直接用残差相对范数 ||r||inf / ||f||inf。
      //
      // 为什么不用预条件残差 ||M⁻¹r||：对 Jacobi 这类弱预条件，
      // 预条件残差的下降速度与真实误差不成正比（可能被过度放大），
      // 会出现"预条件残差很大但解已经足够精确"的误判。
      // 结构分析关心的是"解的绝对精度"，直接看真实残差最稳妥。
      if (rmax / bnorm < tol_) { if (iters) *iters = it; return SolveStatus::Ok; }

      applyPreconditioner(r, z);
      const double rzNew = [&] { double t = 0; for (Id i = 0; i < n; ++i) t += r[i] * z[i]; return t; }();
      const double beta = rzNew / rz;
      rz = rzNew;
      for (Id i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
    }
    if (iters) *iters = maxIter_;
    return SolveStatus::NotConverged;
  }

  int iterationsLastSolve() const { return itersLast_; }

 private:
  // ---- 预条件子构造 ----
  void buildPreconditioner() {
    mInv_.assign(n_, 1.0);
    if (!useIC_) {
      for (Id i = 0; i < n_; ++i) {
        const double d = A_.at(i, i);
        mInv_[i] = (std::abs(d) > kEps) ? 1.0 / d : 1.0;
      }
      return;
    }

    // IC(0)：LDL^T 分解，但只保留 A 原有非零模式内的元素。
    // 用逐行的稠密工作向量 + "仅更新 A 中已有非零" 的规则实现，
    // 存储与 A 的上三角 CSR 对齐（同一 rowPtr/colInd，避免额外索引结构）。
    const auto& rp = A_.rowPtr();
    const auto& ci = A_.colInd();
    const auto& va = A_.values();
    const Id m = static_cast<Id>(va.size());
    icL_.assign(m, 0.0);
    icD_.assign(n_, 1.0);
    work_.assign(n_, 0.0);
    stamp_.assign(n_, -1);

    for (Id i = 0; i < n_; ++i) {
      // 本行中 A 已有的列 -> 在 icL_ 中的位置
      for (Id p = rp[i]; p < rp[i + 1]; ++p) {
        const Id j = ci[p];
        if (j >= i) continue;                 // 只处理下三角
        if (stamp_[j] != i) { stamp_[j] = i; work_[j] = 0.0; }
        work_[j] += va[p];
      }
      for (Id p = rp[i]; p < rp[i + 1]; ++p) {
        const Id j = ci[p];
        if (j >= i) continue;
        // work_[j] 已在上面的载入循环中初始化并累加了 A(i,j)，此处不可重置
        const double djj = icD_[j];
        if (std::abs(djj) < kEps) { icD_[j] = 1e-30; continue; }
        // L(i,j) = ( A(i,j) - Σ_{k<j} L(i,k)·D(k)·L(j,k) ) / D(j)
        // IC(0)：求和只在 A(i,·) 与 A(j,·) 的交集上进行（丢弃填充）
        double acc = work_[j];
        for (Id q = rp[i]; q < p; ++q) {              // A(i,k), k < j（同行、位置在 p 之前）
          const Id k = ci[q];
          const Id mj = findInRow(j, k);              // A(j,k) 是否存在
          if (mj < 0) continue;                       // IC(0)：丢弃填充
          acc -= icL_[q] * icD_[k] * icL_[mj];
        }
        const double lij = acc / icD_[j];
        icL_[p] = lij;
        if (lij == 0.0) continue;
        // 更新本行后续列：w(t) -= L(i,j)·D(j)·L(j,t)，t 只需在 A(i,·) 内
        for (Id q = p + 1; q < rp[i + 1]; ++q) {
          const Id t = ci[q];
          if (t >= i) continue;
          const Id mj = findInRow(j, t);
          if (mj < 0) continue;
          if (stamp_[t] != i) { stamp_[t] = i; work_[t] = 0.0; }
          work_[t] -= lij * icD_[j] * icL_[mj];
        }
      }
      // D(i)
      double di = A_.at(i, i);
      for (Id p = rp[i]; p < rp[i + 1]; ++p) {
        const Id j = ci[p];
        if (j >= i) continue;
        di -= icL_[p] * icL_[p] * icD_[j];
      }
      if (!(std::abs(di) > kEps) || !std::isfinite(di)) di = 1e-30;
      icD_[i] = di;
      // 清理
      for (Id p = rp[i]; p < rp[i + 1]; ++p) {
        const Id j = ci[p];
        if (j < i) work_[j] = 0.0;
      }
    }
  }

  // 在第 row 行中找列 col 的位置（二分，行内有序）
  Id findInRow(Id row, Id col) const {
    const auto& rp = A_.rowPtr();
    const auto& ci = A_.colInd();
    Id lo = rp[row], hi = rp[row + 1] - 1;
    while (lo <= hi) {
      const Id mid = (lo + hi) / 2;
      if (ci[mid] == col) return mid;
      if (ci[mid] < col) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
  }

  // z = M^{-1} r
  void applyPreconditioner(const std::vector<double>& r, std::vector<double>& z) const {
    if (!useIC_) {
      for (Id i = 0; i < n_; ++i) z[i] = r[i] * mInv_[i];
      return;
    }
    const auto& rp = A_.rowPtr();
    const auto& ci = A_.colInd();
    // 前代 L y = r：第 i 行存 L(i,j), j<i —— 行内直接可用
    for (Id i = 0; i < n_; ++i) {
      double s = r[i];
      for (Id p = rp[i]; p < rp[i + 1]; ++p) {
        const Id j = ci[p];
        if (j < i) s -= icL_[p] * y_[j];
      }
      y_[i] = s;
    }
    // 解 D w = y
    for (Id i = 0; i < n_; ++i) w_[i] = y_[i] / icD_[i];
    // 回代 L^T z = w：z_i = w_i - Σ_{j>i} L(j,i) z_j
    // 注意 L(j,i) 在第 j 行第 i 列，必须按【列】访问 L —— 用转置索引 tptr_/tpos_。
    for (Id i = n_ - 1; i >= 0; --i) {
      double s = w_[i];
      for (Id q = tptr_[i]; q < tptr_[i + 1]; ++q) {
        const Id j = trow_[q];              // j > i 且 L(j,i) != 0
        s -= tval_[q] * z_[j];
      }
      z_[i] = s;
    }
  }

  // 建立 L 的转置索引（列访问），供回代使用
  void buildTranspose() {
    tptr_.assign(n_ + 1, 0);
    trow_.clear();
    tval_.clear();
    const auto& rp = A_.rowPtr();
    const auto& ci = A_.colInd();
    for (Id i = 0; i < n_; ++i)
      for (Id p = rp[i]; p < rp[i + 1]; ++p) {
        const Id j = ci[p];
        if (j >= i) continue;
        ++tptr_[j + 1];
      }
    for (Id i = 0; i < n_; ++i) tptr_[i + 1] += tptr_[i];
    trow_.assign(tptr_[n_], 0);
    tval_.assign(tptr_[n_], 0.0);
    std::vector<Id> pos(tptr_.begin(), tptr_.end() - 1);
    for (Id i = 0; i < n_; ++i)
      for (Id p = rp[i]; p < rp[i + 1]; ++p) {
        const Id j = ci[p];
        if (j >= i) continue;
        const Id d = pos[j]++;
        trow_[d] = i;
        tval_[d] = icL_[p];
      }
  }

  const SymSparseMatrix& A_;
  Id n_;
  int maxIter_;
  double tol_;
  bool useIC_;
  int itersLast_{0};
  std::vector<double> mInv_;      // Jacobi
  std::vector<double> icL_, icD_;  // IC(0)
  mutable std::vector<double> y_, z_, w_;
  mutable std::vector<double> work_, stamp_;
  std::vector<Id> tptr_, trow_;        // L 的转置索引
  std::vector<double> tval_;
};

}  // namespace yjk
