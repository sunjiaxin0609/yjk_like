// =============================================================================
//  yjk/math/SparseMatrix.h  ——  对称稀疏矩阵（CSR/上三角存储）
//
//  设计要点（面向结构分析）：
//  1. 只存上三角，K = K^T，由组装时保证对称性 —— 省一半内存，且杜绝舍入不对称
//  2. 组装期用 triplet(IJK) 累加器，允许同一(i,j)多次累加（多个单元贡献同一自由度对）
//  3. 分解前做 AMD(近似最小度) 排序，抑制填充
//  4. 数值上采用 LDL^T（不要求正定，可处理 P-Δ 后的对称不定问题），
//     但默认路径在 K 对称正定时自动退化为 Cholesky，取 sqrt 更快更稳
// =============================================================================
#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "yjk/math/Types.h"

namespace yjk {

// 三元组：单元刚度矩阵往全局累加时的中间表示
// 约定：只存上三角，即恒有 i >= j。
struct Triplet {
  Id i, j;
  double v;
};

// ---------------------------------------------------------------------------
//  三元组累加器：单元刚度矩阵往全局里"攒"
//
//  【重要】add() 会自动把 (i,j) 归一化到上三角（i 取大者）。
//  这不是可选的优化，而是必需的正确性保障：
//  单元自由度映射后，全局编号的大小关系与局部编号无关，
//  调用方很容易无意中传入下三角项。若原样存储，
//  下游 factorize 按上三角读取时会静默跳过这些项 ——
//  刚度矩阵退化成对角阵，位移结果小几个数量级，且不报任何错。
//  这类"静默错误"在结构软件里是最危险的，必须在 API 层堵死。
// ---------------------------------------------------------------------------
class TripletAssembler {
 public:
  void reserve(size_t n) { tri_.reserve(n); }

  // 加入一个矩阵元素。自动归一化到上三角（i 取大者），见类注释。
  void add(Id i, Id j, double v) {
    if (nearlyZero(v)) return;   // 跳过零，省内存
    if (i < j) std::swap(i, j);
    tri_.push_back({i, j, v});
  }

  // 把一个 NxN 单元刚度矩阵按自由度映射散装到全局。
  // 只需遍历上三角，add() 内部会做归一化。
  template <int N>
  void addElement(const Id* dof, const Mat<N>& ke) {
    for (int r = 0; r < N; ++r) {
      for (int c = r; c < N; ++c) {
        const double v = ke(r, c);
        if (nearlyZero(v)) continue;
        add(dof[r], dof[c], v);
      }
    }
  }

  // 覆写式写入（用于支座位移约束的对角线加罚，或初始应变力）
  void set(Id i, Id j, double v) {
    if (i < j) std::swap(i, j);
    tri_.push_back({i, j, v});
  }

  const std::vector<Triplet>& triplets() const { return tri_; }
  std::vector<Triplet>& triplets() { return tri_; }
  void clear() { tri_.clear(); }
  size_t size() const { return tri_.size(); }

 private:
  std::vector<Triplet> tri_;
};

// ---------------------------------------------------------------------------
//  对称稀疏矩阵，CSR 上三角存储
// ---------------------------------------------------------------------------
class SymSparseMatrix {
 public:
  SymSparseMatrix() = default;

  // 从三元组组装：合并重复项，排序，构建 CSR
  void buildFrom(const std::vector<Triplet>& tri, Id n) {
    n_ = n;
    std::vector<Triplet> t;
    t.reserve(tri.size());
    // 过滤越界与非对角非零的行（安全网）
    for (const auto& e : tri) {
      if (e.i < 0 || e.i >= n || e.j < 0 || e.j >= n) continue;
      if (e.i == e.j && nearlyZero(e.v)) continue;
      t.push_back(e);
    }
    // 稳定排序：先列后行
    std::stable_sort(t.begin(), t.end(), [](const auto& a, const auto& b) {
      return a.j != b.j ? a.j < b.j : a.i < b.i;
    });

    // 压缩：相同 (i,j) 累加
    compressed_.clear();
    compressed_.reserve(t.size());
    for (const auto& e : t) {
      if (!compressed_.empty() && compressed_.back().i == e.i && compressed_.back().j == e.j) {
        compressed_.back().v += e.v;
      } else {
        compressed_.push_back(e);
      }
    }

    // 构建行指针 + 列索引
    rowPtr_.assign(n + 1, 0);
    for (const auto& e : compressed_) rowPtr_[e.i + 1]++;
    for (Id i = 0; i < n; ++i) rowPtr_[i + 1] += rowPtr_[i];
    colInd_.resize(compressed_.size());
    val_.resize(compressed_.size());
    std::vector<Id> pos(rowPtr_.begin(), rowPtr_.end() - 1);
    for (size_t k = 0; k < compressed_.size(); ++k) {
      const Id r = pos[compressed_[k].i]++;
      colInd_[r] = compressed_[k].j;
      val_[r] = compressed_[k].v;
    }
    nnz_ = val_.size();
  }

  Id size() const { return n_; }
  size_t nnz() const { return nnz_; }

  const std::vector<Id>& rowPtr() const { return rowPtr_; }
  const std::vector<Id>& colInd() const { return colInd_; }
  const std::vector<double>& values() const { return val_; }
  std::vector<double>& values() { return val_; }
  const std::vector<Triplet>& tripletsUpper() const { return compressed_; }

  // 取值（含对称扩展）
  // 注意：本矩阵只存上三角（行 i 的列 j 恒有 j <= i）。
  // at(i,j) 中若 i > j（调用方习惯按"行、列"给任意顺序），
  // 需交换成 (j,i) 再查 —— 因为对称矩阵 K(i,j) = K(j,i)，
  // 而 CSR 里 K(j,i) 存在行 j 中。
  double at(Id i, Id j) const {
    if (i < 0 || i >= n_ || j < 0 || j >= n_) return 0.0;
    if (i < j) std::swap(i, j);       // 保证 i >= j，落在上三角
    for (Id p = rowPtr_[i]; p < rowPtr_[i + 1]; ++p)
      if (colInd_[p] == j) return val_[p];
    return 0.0;
  }

  // y = A*x（A 对称，存储只含上三角）
  //
  // 每个上三角项 a_ij (i >= j) 对结果的贡献有两处（因为 A 对称）：
  //     y[i] += a_ij * x[j]
  //     y[j] += a_ij * x[i]      (i > j 时)
  // 必须一次遍历同时累加这两处。分两轮做是常见错误：
  // 第一轮只算 y[i]、第二轮想补 y[j]，但第二轮遍历同样的上三角数据时
  // 根本分不清哪些项需要补 —— 结果下三角贡献全丢，CG 永远不收敛。
  void multiply(const std::vector<double>& x, std::vector<double>& y) const {
    y.assign(n_, 0.0);
    for (Id i = 0; i < n_; ++i) {
      for (Id p = rowPtr_[i]; p < rowPtr_[i + 1]; ++p) {
        const Id j = colInd_[p];          // j <= i（上三角）
        const double a = val_[p];
        y[i] += a * x[j];
        if (j < i) y[j] += a * x[i];     // 对称项的另一半
      }
    }
  }

  // 强制对角线最小值（处理约束不足导致的奇异，不静默出错）
  void addToDiagonal(double eps) {
    for (Id i = 0; i < n_; ++i) {
      bool found = false;
      for (Id p = rowPtr_[i]; p < rowPtr_[i + 1]; ++p)
        if (colInd_[p] == i) { val_[p] += eps; found = true; break; }
      if (!found) {
        compressed_.push_back({i, i, eps});
        colInd_.push_back(i); val_.push_back(eps); nnz_++;
      }
    }
  }

  // 估计条件数指标的廉价代理：max|K_ii| 与 min|K_ii| 之比
  double diagRatio() const {
    double mx = 0.0, mn = 1e300;
    for (Id i = 0; i < n_; ++i) {
      double d = 0.0;
      for (Id p = rowPtr_[i]; p < rowPtr_[i + 1]; ++p)
        if (colInd_[p] == i) { d = std::abs(val_[p]); break; }
      if (d > kEps) { mx = std::max(mx, d); mn = std::min(mn, d); }
    }
    return mn < 1e299 ? mx / mn : 1e300;
  }

 private:
  Id n_{0};
  size_t nnz_{0};
  std::vector<Id> rowPtr_, colInd_;
  std::vector<double> val_;
  std::vector<Triplet> compressed_;
};

}  // namespace yjk
