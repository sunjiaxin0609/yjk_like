// =============================================================================
//  src/math/EigenSolver.cpp
// =============================================================================
#include "yjk/math/EigenSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>

namespace yjk {

namespace {

// ---------------------------------------------------------------------------
//  稠密小工具（n 很小，O(n³) 足够；只为可读性与确定性，不为性能）
// ---------------------------------------------------------------------------
inline double matAt(const std::vector<double>& A, int n, int i, int j) {
  return A[static_cast<size_t>(i) * n + j];
}
inline double& matAt(std::vector<double>& A, int n, int i, int j) {
  return A[static_cast<size_t>(i) * n + j];
}

// 对 M 做 Cholesky：M = L Lᵀ，L 严格下三角 + 对角（返回下三角按列主序）
// 失败（非正定）返回 false。
bool cholesky(const std::vector<double>& M, int n, std::vector<double>& L) {
  L.assign(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = matAt(M, n, i, j);
      for (int k = 0; k < j; ++k) s -= matAt(L, n, i, k) * matAt(L, n, j, k);
      if (i == j) {
        if (s <= 1e-300) return false;   // M 必须正定
        matAt(L, n, i, i) = std::sqrt(s);
      } else {
        matAt(L, n, i, j) = s / matAt(L, n, j, j);
      }
    }
  }
  return true;
}

// 求解 L X = B（L 为 Cholesky 下三角因子，B 按列主序，X 输出按列主序，B 可复用为 X）
void solveL(const std::vector<double>& L, int n, std::vector<double>& B) {
  for (int c = 0; c < n; ++c)
    for (int r = 0; r < n; ++r) {
      double acc = matAt(B, n, r, c);
      for (int k = 0; k < r; ++k)
        acc -= matAt(L, n, r, k) * matAt(B, n, k, c);
      matAt(B, n, r, c) = acc / matAt(L, n, r, r);
    }
}

// 求解 Lᵀ φ = z（上三角回代），φ 覆盖 z
void solveLt(const std::vector<double>& L, int n, std::vector<double>& z) {
  for (int i = n - 1; i >= 0; --i) {
    double acc = z[static_cast<size_t>(i)];
    for (int j = i + 1; j < n; ++j)
      acc -= matAt(L, n, j, i) * z[static_cast<size_t>(j)];
    z[static_cast<size_t>(i)] = acc / matAt(L, n, i, i);
  }
}

// ---------------------------------------------------------------------------
//  经典雅可比旋转：对称矩阵 A（n×n 稠密）对角化，返回特征对
//  A = Z diag(μ) Zᵀ，Z 累积旋转（列 = 特征向量）
// ---------------------------------------------------------------------------
struct JacobiResult {
  bool ok{false};
  int sweeps{0};
  double off{0.0};
};

JacobiResult jacobi(std::vector<double>& A, int n, std::vector<double>& Z) {
  JacobiResult res;
  Z.assign(static_cast<size_t>(n) * n, 0.0);
  if (n == 0) { res.ok = true; return res; }
  for (int i = 0; i < n; ++i) matAt(Z, n, i, i) = 1.0;

  for (int sweep = 0; sweep < 100; ++sweep) {
    // 一轮扫描：所有上三角对 (p,q)
    double off2 = 0.0;
    int nrot = 0;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) {
        const double apq = matAt(A, n, p, q);
        if (std::abs(apq) < 1e-300) continue;
        // 旋转角：tan(2θ) = 2apq / (aqq − app)
        const double app = matAt(A, n, p, p), aqq = matAt(A, n, q, q);
        const double tau = (aqq - app) / (2.0 * apq);
        const double t = (tau >= 0.0) ? 1.0 / (tau + std::sqrt(1.0 + tau * tau))
                                      : 1.0 / (tau - std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t), s = c * t;
        // 旋转 A：只更新 p,q 行/列
        for (int k = 0; k < n; ++k) {
          if (k == p || k == q) continue;
          const double akp = matAt(A, n, k, p), akq = matAt(A, n, k, q);
          matAt(A, n, k, p) = c * akp - s * akq;
          matAt(A, n, p, k) = matAt(A, n, k, p);
          matAt(A, n, k, q) = s * akp + c * akq;
          matAt(A, n, q, k) = matAt(A, n, k, q);
        }
        const double appNew = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        const double aqqNew = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        matAt(A, n, p, p) = appNew;
        matAt(A, n, q, q) = aqqNew;
        matAt(A, n, p, q) = 0.0;
        matAt(A, n, q, p) = 0.0;
        // 累积 Z
        for (int k = 0; k < n; ++k) {
          const double zkp = matAt(Z, n, k, p), zkq = matAt(Z, n, k, q);
          matAt(Z, n, k, p) = c * zkp - s * zkq;
          matAt(Z, n, k, q) = s * zkp + c * zkq;
        }
        ++nrot;
      }
    // 收敛判据：非对角平方和相对值
    off2 = 0.0;
    double diag2 = 0.0;
    for (int p = 0; p < n; ++p)
      for (int q = 0; q < n; ++q) {
        const double v = matAt(A, n, p, q);
        if (p == q) diag2 += v * v; else off2 += v * v;
      }
    res.off = std::sqrt(off2);
    if (res.off <= 1e-12 * std::sqrt(diag2) || (nrot == 0)) {
      res.ok = true;
      res.sweeps = sweep + 1;
      return res;
    }
  }
  res.ok = false;   // 100 轮未收敛（实际不会发生）
  return res;
}

// 确定性的伪随机（固定种子 ⇒ 子空间迭代结果可复现）
struct Rng {
  unsigned s;
  explicit Rng(unsigned seed) : s(seed) {}
  double next() {          // [0,1)
    s = s * 1664525u + 1013904223u;
    return (s >> 8) * (1.0 / 16777216.0);
  }
};

}  // namespace

// -----------------------------------------------------------------------------
//  DenseGeneralizedEigen
// -----------------------------------------------------------------------------
bool DenseGeneralizedEigen::solve(const std::vector<double>& K,
                                  const std::vector<double>& M, int n,
                                  std::vector<EigenPair>& out, double tol,
                                  int maxSweeps) {
  out.clear();
  if (n <= 0) return true;

  // 1. M = L Lᵀ
  std::vector<double> L;
  if (!cholesky(M, n, L)) return false;

  // 2. A = L⁻¹ K L⁻ᵀ
  //    B = L⁻¹K（按列解 L B = K）
  std::vector<double> B = K;
  solveL(L, n, B);
  //    C = L⁻¹ Bᵀ ⇒ A = Cᵀ = B L⁻ᵀ
  std::vector<double> BT(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) matAt(BT, n, i, j) = matAt(B, n, j, i);
  solveL(L, n, BT);
  std::vector<double> A(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) matAt(A, n, i, j) = matAt(BT, n, j, i);

  // 3. 雅可比对角化 A = Z diag(μ) Zᵀ
  std::vector<double> Zt;
  auto jac = jacobi(A, n, Zt);
  if (!jac.ok) return false;

  // 4. 特征对 (μ, z)：φ = L⁻ᵀ z，并按 φᵀMφ = 1 归一
  out.resize(static_cast<size_t>(n));
  for (int c = 0; c < n; ++c) {
    EigenPair& ep = out[static_cast<size_t>(c)];
    ep.value = matAt(A, n, c, c);
    ep.vector.assign(static_cast<size_t>(n), 0.0);
    for (int r = 0; r < n; ++r)
      ep.vector[static_cast<size_t>(r)] = matAt(Zt, n, r, c);
    solveLt(L, n, ep.vector);
    // 归一化 φᵀMφ = 1
    double mnorm = 0.0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        mnorm += ep.vector[static_cast<size_t>(i)] * matAt(M, n, i, j) * ep.vector[static_cast<size_t>(j)];
    if (mnorm > 0.0 && std::isfinite(mnorm)) {
      const double s = 1.0 / std::sqrt(mnorm);
      for (double& v : ep.vector) v *= s;
    }
  }
  std::sort(out.begin(), out.end(),
            [](const EigenPair& a, const EigenPair& b) { return a.value < b.value; });
  (void)tol; (void)maxSweeps;
  return true;
}

// -----------------------------------------------------------------------------
//  SubspaceIteration
// -----------------------------------------------------------------------------
bool SubspaceIteration::solve(const SymSparseMatrix& K, const std::vector<double>& M,
                              int nmodes, std::vector<EigenPair>& out) {
  return solve(K, M, nmodes, out, Options(), nullptr);
}

// 对角质量向量版本：构造对角稀疏矩阵后统一走稀疏版本（单一实现，避免两处漂移）
bool SubspaceIteration::solve(const SymSparseMatrix& K, const std::vector<double>& M,
                              int nmodes, std::vector<EigenPair>& out,
                              Options opt, Report* rep) {
  std::vector<Triplet> tri;
  tri.reserve(M.size());
  for (size_t i = 0; i < M.size(); ++i)
    if (M[i] != 0.0)
      tri.push_back(Triplet{static_cast<Id>(i), static_cast<Id>(i), M[i]});
  SymSparseMatrix Ms;
  Ms.buildFrom(tri, static_cast<Id>(M.size()));
  return solve(K, Ms, nmodes, out, opt, rep);
}

bool SubspaceIteration::solve(const SymSparseMatrix& K, const SymSparseMatrix& M,
                              int nmodes, std::vector<EigenPair>& out) {
  return solve(K, M, nmodes, out, Options(), nullptr);
}

bool SubspaceIteration::solve(const SymSparseMatrix& K, const SymSparseMatrix& M,
                              int nmodes, std::vector<EigenPair>& out,
                              Options opt, Report* rep) {
  out.clear();
  Report r;
  const auto t0 = std::chrono::steady_clock::now();
  r.n = static_cast<int>(K.size());

  if (r.n == 0 || nmodes <= 0 || nmodes > r.n) {
    r.message = "自由度或模态数不合法";
    if (rep) *rep = r;
    return false;
  }
  if (static_cast<int>(M.size()) < r.n) {
    r.message = "质量矩阵维度不足";
    if (rep) *rep = r;
    return false;
  }

  // 工作子空间维数 q
  int q = opt.workDim > 0 ? opt.workDim : std::max(nmodes + 8, 2 * nmodes);
  q = std::min(q, r.n);
  r.q = q;
  r.nnzK = K.nnz();

  // 有效质量自由度：M 对角 > 0 的自由度 —— 决定 M 内积的有效子空间。
  // 对半正定 M，对角为 0 ⇒ 整行/列恒为 0（|M_ij| ≤ √(M_ii·M_jj)），
  // 因此只用对角就能确定有效支撑，无需扫描耦合项。
  std::vector<int> heavy;
  heavy.reserve(r.n);
  for (int i = 0; i < r.n; ++i)
    if (M.at(i, i) > 1e-30) heavy.push_back(i);
  if (heavy.empty()) {
    r.message = "质量矩阵全零（没有可激发模态的平动质量）";
    if (rep) *rep = r;
    return false;
  }

  // 工作子空间维数不超过有效质量自由度：无质量 dof 不贡献谱维数
  // （逆迭代 Y=K⁻¹MQ 会自动把无质量 dof 的静力凝聚位移算进 Q，
  //   因此 Ritz 投影只取带质量子空间的 q 维即可，无质量分量随 Q 携带）。
  if (q > static_cast<int>(heavy.size())) {
    q = static_cast<int>(heavy.size());
  }
  r.q = q;
  if (nmodes > q) {
    r.message = "请求模态数超过有效质量自由度";
    if (rep) *rep = r;
    return false;
  }

  // 初始工作向量 X（n×q）：第 0 列叠全 1 激发最低阶，其余取"质量最大的
  // 自由度"单位向量 + 确定性随机分量（保证线性独立且不与目标子空间正交）。
  std::vector<double> X(static_cast<size_t>(r.n) * q, 0.0);
  {
    std::sort(heavy.begin(), heavy.end(),
              [&](int a, int b) { return M.at(a, a) > M.at(b, b); });
    for (int j = 0; j < q; ++j) {
      const int h = heavy[static_cast<size_t>(j) % heavy.size()];
      X[static_cast<size_t>(j) * r.n + h] = M.at(h, h);   // 单位向量（质量加权）
      Rng rng(opt.seed + static_cast<unsigned>(j) * 101u);
      for (int i = 0; i < r.n; ++i)
        X[static_cast<size_t>(j) * r.n + i] += 1e-3 * M.at(i, i) * rng.next();
    }
    if (r.n > 0) for (int i = 0; i < r.n; ++i) X[static_cast<size_t>(i)] += M.at(i, i);  // 第 0 列叠 M 对角
  }

  // K 的 LDL^T 分解只做一次，之后每轮 q 次回代
  LDLTSolver solver;
  const auto repF = solver.factorize(K, 0);
  if (!repF.ok()) {
    r.message = "刚度矩阵分解失败（模态分析）";
    if (rep) *rep = r;
    return false;
  }

  // 工作向量：Y(n×q) 与 Q(n×q)。第 j 列 = vec[j*n + i]
  std::vector<double> Q(static_cast<size_t>(r.n) * q, 0.0);
  auto orthonormalize = [&](const std::vector<double>& V, std::vector<double>& Vout) {
    // MGS（带两次重正交，避免丢正交性导致的迭代停滞），M-内积走稀疏矩阵乘
    std::vector<double> Mv(r.n, 0.0);
    Vout = V;
    for (int j = 0; j < q; ++j) {
      for (int repass = 0; repass < 2; ++repass) {
        for (int k = 0; k < j; ++k) {
          double* vj = &Vout[static_cast<size_t>(j) * r.n];
          const double* qk = &Vout[static_cast<size_t>(k) * r.n];
          M.multiply(std::vector<double>(vj, vj + r.n), Mv);
          double alpha = 0.0;
          for (int h : heavy) alpha += Mv[static_cast<size_t>(h)] * qk[h];
          for (int i = 0; i < r.n; ++i) vj[i] -= alpha * qk[i];
        }
      }
      double* vj = &Vout[static_cast<size_t>(j) * r.n];
      std::vector<double> Mvj(r.n, 0.0);
      M.multiply(std::vector<double>(vj, vj + r.n), Mvj);
      double nrm = 0.0;
      for (int h : heavy) nrm += vj[h] * Mvj[h];
      if (nrm > 1e-300) {
        nrm = 1.0 / std::sqrt(nrm);
        for (int i = 0; i < r.n; ++i) vj[i] *= nrm;
        continue;
      }
      // 数值退化：用伪随机方向补齐（保证有 M-范数）
      Rng rng(opt.seed + static_cast<unsigned>(j) * 1009u + 7u);
      for (int i = 0; i < r.n; ++i) vj[i] = rng.next() - 0.5;
      for (int repass = 0; repass < 2; ++repass)
        for (int k = 0; k < j; ++k) {
          const double* qk = &Vout[static_cast<size_t>(k) * r.n];
          M.multiply(std::vector<double>(vj, vj + r.n), Mv);
          double alpha = 0.0;
          for (int h : heavy) alpha += Mv[static_cast<size_t>(h)] * qk[h];
          for (int i = 0; i < r.n; ++i) vj[i] -= alpha * qk[i];
        }
      M.multiply(std::vector<double>(vj, vj + r.n), Mvj);
      nrm = 0.0;
      for (int h : heavy) nrm += vj[h] * Mvj[h];
      if (nrm <= 1e-300) {
        // 已无关紧：给一个 {1,0,...} 归一
        vj[0] = 1.0;
        for (int i = 1; i < r.n; ++i) vj[i] = 0.0;
        continue;
      }
      nrm = 1.0 / std::sqrt(nrm);
      for (int i = 0; i < r.n; ++i) vj[i] *= nrm;
    }
  };

  std::vector<double> Y(static_cast<size_t>(r.n) * q, 0.0);
  std::vector<double> b(r.n, 0.0), Yc(r.n, 0.0);
  std::vector<double> Kq(static_cast<size_t>(q) * q, 0.0);
  std::vector<double> Mq(static_cast<size_t>(q) * q, 0.0);
  std::vector<EigenPair> ritz;

  orthonormalize(X, Q);
  double worstRes = 1e300;

  for (int iter = 1; iter <= opt.maxIter; ++iter) {
    // Y = K⁻¹ M Q：每列一次前代+回代（M 稀疏乘）
    for (int j = 0; j < q; ++j) {
      const double* qj = &Q[static_cast<size_t>(j) * r.n];
      M.multiply(std::vector<double>(qj, qj + r.n), b);
      if (solver.solve(b, Yc) != SolveStatus::Ok) {
        r.message = "逆迭代回代失败（奇异？）";
        if (rep) *rep = r;
        return false;
      }
      for (int i = 0; i < r.n; ++i) Y[static_cast<size_t>(j) * r.n + i] = Yc[static_cast<size_t>(i)];
    }
    orthonormalize(Y, Q);   // M-正交归一化 ⇒ 新工作子空间

    // Ritz 投影：Kq = QᵀKQ、Mq = QᵀMQ（小稠密）
    std::vector<double> qcol(r.n, 0.0), Kqj(r.n, 0.0), Mqj(r.n, 0.0);
    for (int j = 0; j < q; ++j) {
      const double* qj = &Q[static_cast<size_t>(j) * r.n];
      std::copy(qj, qj + r.n, qcol.begin());
      K.multiply(qcol, Kqj);
      M.multiply(qcol, Mqj);
      double* KqjRow = &matAt(Kq, q, j, 0);
      double* MqjRow = &matAt(Mq, q, j, 0);
      for (int k = 0; k < q; ++k) {
        const double* qk = &Q[static_cast<size_t>(k) * r.n];
        // 刚度内积 Kq(j,k) = (K·q_j)ᵀ·q_k（K 对称 ⇒ 与 (K·q_k)ᵀ·q_j 相等）
        double kjk = 0.0, mjk = 0.0;
        for (int i = 0; i < r.n; ++i) {
          kjk += Kqj[static_cast<size_t>(i)] * qk[i];
          mjk += Mqj[static_cast<size_t>(i)] * qk[i];
        }
        KqjRow[k] = kjk;
        MqjRow[k] = mjk;
      }
    }

    // 小广义问题 Kq z = μ Mq z（q ≤ 几百，稠密直接解）
    if (!DenseGeneralizedEigen::solve(Kq, Mq, q, ritz, 1e-12, 200)) {
      r.message = "Ritz 投影求解失败";
      if (rep) *rep = r;
      return false;
    }

    // 新迭代向量 X = Q · Z（Z = Ritz 向量矩阵），本质保持 M-正交。
    // 注意：先保留本轮投影基 Q 用于收敛判据，判据通过后再更新 Q。
    std::vector<double> Xnew(static_cast<size_t>(r.n) * q, 0.0);
    for (int j = 0; j < q; ++j)
      for (int k = 0; k < q; ++k) {
        const double zjk = ritz[static_cast<size_t>(k)].vector[static_cast<size_t>(j)];
        const double* qk = &Q[static_cast<size_t>(k) * r.n];
        double* xj = &Xnew[static_cast<size_t>(j) * r.n];
        for (int i = 0; i < r.n; ++i) xj[i] += zjk * qk[i];
      }

    // 收敛判据：前 nmodes 阶残差。
    // 注意：ritz[m].vector 是 Q 基下的坐标 z_m（q 维），必须还原为物理
    // 振型 φ_m = Q·z_m（n 维）后才能求 Kφ − λMφ 的残差。
    worstRes = 0.0;
    r.iter = iter;
    r.residual = 0.0;
    std::vector<double> phi(r.n, 0.0);
    bool converged = true;
    for (int m = 0; m < nmodes; ++m) {
      const EigenPair& rp = ritz[static_cast<size_t>(m)];
      std::fill(phi.begin(), phi.end(), 0.0);
      for (int k = 0; k < q; ++k) {
        const double zkm = rp.vector[static_cast<size_t>(k)];
        const double* qk = &Q[static_cast<size_t>(k) * r.n];
        for (int i = 0; i < r.n; ++i) phi[static_cast<size_t>(i)] += zkm * qk[i];
      }
      const double res = residual(K, M, {rp.value, phi});
      worstRes = std::max(worstRes, res);
      if (!(res <= opt.tol)) converged = false;
    }
    r.residual = worstRes;
    if (converged) {
      r.ok = true;
      r.message = "收敛";
      r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      // 输出物理振型 φ_m = Q·z_m（Q 仍是本轮投影基，未被覆盖）
      out.clear();
      out.reserve(static_cast<size_t>(nmodes));
      for (int m = 0; m < nmodes; ++m) {
        const EigenPair& rp = ritz[static_cast<size_t>(m)];
        EigenPair ep;
        ep.value = rp.value;
        ep.vector.assign(static_cast<size_t>(r.n), 0.0);
        for (int k = 0; k < q; ++k) {
          const double zkm = rp.vector[static_cast<size_t>(k)];
          const double* qk = &Q[static_cast<size_t>(k) * r.n];
          for (int i = 0; i < r.n; ++i)
            ep.vector[static_cast<size_t>(i)] += zkm * qk[i];
        }
        out.push_back(std::move(ep));
      }
      if (rep) *rep = r;
      return true;
    }

    // 未收敛：更新工作基，进入下一轮
    orthonormalize(Xnew, Q);
  }

  r.message = "达到最大迭代次数未收敛（残差 " + std::to_string(worstRes) + "）";
  if (rep) *rep = r;
  return false;
}

double SubspaceIteration::rayleigh(const SymSparseMatrix& K, const std::vector<double>& M,
                                   const std::vector<double>& phi) {
  std::vector<double> Kx;
  K.multiply(phi, Kx);
  double num = 0.0, den = 0.0;
  for (size_t i = 0; i < phi.size(); ++i) {
    num += phi[i] * Kx[i];
    den += M[i] * phi[i] * phi[i];
  }
  return (den > 1e-300) ? num / den : 0.0;
}

double SubspaceIteration::residual(const SymSparseMatrix& K, const std::vector<double>& M,
                                   const EigenPair& ep) {
  std::vector<double> Ax;
  K.multiply(ep.vector, Ax);
  double num = 0.0, den = 0.0;
  const size_t n = ep.vector.size();
  for (size_t i = 0; i < n; ++i) {
    num = std::max(num, std::abs(Ax[i] - ep.value * M[i] * ep.vector[i]));
    den = std::max(den, std::abs(Ax[i]));
  }
  return (den > 1e-300) ? num / den : num;
}

double SubspaceIteration::rayleigh(const SymSparseMatrix& K, const SymSparseMatrix& M,
                                   const std::vector<double>& phi) {
  std::vector<double> Kx, Mx;
  K.multiply(phi, Kx);
  M.multiply(phi, Mx);
  double num = 0.0, den = 0.0;
  for (size_t i = 0; i < phi.size(); ++i) {
    num += phi[i] * Kx[i];
    den += phi[i] * Mx[i];
  }
  return (den > 1e-300) ? num / den : 0.0;
}

double SubspaceIteration::residual(const SymSparseMatrix& K, const SymSparseMatrix& M,
                                   const EigenPair& ep) {
  std::vector<double> Ax, Mx;
  K.multiply(ep.vector, Ax);
  M.multiply(ep.vector, Mx);
  double num = 0.0, den = 0.0;
  const size_t n = ep.vector.size();
  for (size_t i = 0; i < n; ++i) {
    num = std::max(num, std::abs(Ax[i] - ep.value * Mx[i]));
    den = std::max(den, std::abs(Ax[i]));
  }
  return (den > 1e-300) ? num / den : num;
}

}  // namespace yjk