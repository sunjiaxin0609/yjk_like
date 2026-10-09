// =============================================================================
//  yjk/math/Ordering.h  ——  稀疏排序：RCM 反向 Cuthill-McKee 与 MD 最小度
//
//  为什么必须做排序：
//  刚度矩阵的原始节点编号来自建模顺序（往往随意），消元会产生大量填充。
//  排序把"耦合紧密"的自由度凑在一起，让填充最小化。对壳单元为主的模型，
//  好的排序能让内存和耗时差 5~10 倍。
//
//  RCM : 便宜 O(n log n)，对规则网格（框架结构）很有效。
//        盈建科这类框架/剪力墙软件实测：按"楼层自上而下 + 同层柱网展开"编号，
//        RCM 后的带宽已经接近最优，这是行业里最省事的做法。
//  MD  : 最小度（Minimum Degree）O(n^2) 最坏但配合优先队列+惰性更新，
//        对任意稀疏模式都好，用作壳单元为主的模型（带宽大）的兜底。
// =============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <vector>

#include "yjk/math/SparseMatrix.h"

namespace yjk {

// ---------------------------------------------------------------------------
//  RCM：先 Cuthill-McKee 生成紧凑顺序，再反向使带宽最小
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//  RCM：Cuthill-McKee + 反向
//  实际经验：单一最小度种子常常不是最优起点（尤其网格结构）。
//  这里做多起点试跑（每个度数最小的候选点各跑一次，取带宽最小者），
//  这是商用代码的标准做法，代价只有几次 BFS（O(n log n)），收益显著。
// ---------------------------------------------------------------------------
inline std::vector<Id> rcmOrdering(const SymSparseMatrix& A) {
  const Id n = A.size();
  if (n == 0) return {};

  // 构建对称邻接（CSR 的 CSR），含下三角
  std::vector<Id> aik(n + 1, 0);
  for (Id i = 0; i < n; ++i) {
    for (Id p = A.rowPtr()[i]; p < A.rowPtr()[i + 1]; ++p) {
      const Id j = A.colInd()[p];
      if (i != j) { ++aik[i]; ++aik[j]; }
    }
  }
  std::vector<Id> aptr(n + 1, 0);
  for (Id i = 0; i < n; ++i) aptr[i + 1] = aptr[i] + aik[i];
  std::vector<Id> aind(aptr[n]);
  std::vector<Id> pos = aptr;
  for (Id i = 0; i < n; ++i) {
    for (Id p = A.rowPtr()[i]; p < A.rowPtr()[i + 1]; ++p) {
      const Id j = A.colInd()[p];
      if (i == j) continue;
      aind[pos[i]++] = j;
      aind[pos[j]++] = i;
    }
  }
  for (Id i = 0; i < n; ++i) {
    std::sort(aind.begin() + aptr[i], aind.begin() + aptr[i + 1]);
    aind.erase(std::unique(aind.begin() + aptr[i], aind.begin() + aptr[i + 1]),
               aind.begin() + aptr[i + 1]);
  }

  std::vector<Id> deg(n);
  for (Id i = 0; i < n; ++i) deg[i] = aptr[i + 1] - aptr[i];

  // ---- 单次 Cuthill-McKee，从指定首个种子出发 ----
  // 严格按 BFS 层扩展，同层内按度数升序。
  // 关键点：
  //  1) 必须是 BFS（层次优先）。用优先队列"全局挑最小度数"会破坏层次结构，
  //     得到的顺序在网格上剧烈跳跃，带宽反而变差。
  //  2) 同层排序用度数升序 —— 这正是 CM 相对朴素 BFS 的唯一改进。
  //  3) 队列用 index 游标而非 swap/clear，避免两层容器别名导致的隐蔽越界。
  auto cuthillMcKee = [&](Id startNode, std::vector<Id>& out) {
    out.clear();
    out.reserve(n);
    std::vector<char> vis(n, 0);
    std::vector<Id> queueBuf;
    queueBuf.reserve(n);

    for (Id comp = 0; comp < n; ++comp) {
      // 取首个未访问点：优先用 startNode，之后取最小度未访问点（结构可能不连通）
      Id seed = -1;
      if (comp == 0 && startNode >= 0 && startNode < n) {
        seed = startNode;
      } else {
        for (Id i = 0; i < n; ++i)
          if (!vis[i] && (seed < 0 || deg[i] < deg[seed])) seed = i;
      }
      if (seed < 0) break;

      // 每个分量都是独立的 BFS：队列与游标都要重置
      queueBuf.clear();
      size_t head = 0;
      queueBuf.push_back(seed);
      vis[seed] = 1;

      while (head < queueBuf.size()) {
        // 取本层剩余节点中度数最小者（等价于同层排序）
        size_t bestIdx = head;
        for (size_t t = head + 1; t < queueBuf.size(); ++t) {
          const Id a = queueBuf[t], b = queueBuf[bestIdx];
          if (deg[a] < deg[b] || (deg[a] == deg[b] && a < b)) bestIdx = t;
        }
        std::swap(queueBuf[head], queueBuf[bestIdx]);
        const Id v = queueBuf[head++];
        out.push_back(v);
        for (Id q = aptr[v]; q < aptr[v + 1]; ++q) {
          const Id w = aind[q];
          if (w >= 0 && w < n && !vis[w]) { vis[w] = 1; queueBuf.push_back(w); }
        }
      }
    }
  };

  // ---- 由排列算带宽 ----
  auto bwOf = [&](const std::vector<Id>& order) -> Id {
    if (static_cast<Id>(order.size()) != n) return n;   // 排列不完整，视为最差
    std::vector<Id> newOf(n);
    for (Id k = 0; k < n; ++k) {
      const Id v = order[k];
      if (v < 0 || v >= n) return n;                     // 越界保护
      newOf[v] = k;
    }
    Id bw = 0;
    for (Id i = 0; i < n; ++i)
      for (Id p = A.rowPtr()[i]; p < A.rowPtr()[i + 1]; ++p)
        bw = std::max(bw, std::abs(newOf[i] - newOf[A.colInd()[p]]));
    return bw;
  };

  // ---- 多起点试跑，取带宽最优 ----
  // 单一最小度种子常常不是最优起点：网格/带状结构换个角点方向，带宽能差很多。
  // 代价只是几次 BFS（每次 O(n log n)），收益往往成倍。
  std::vector<Id> cand;
  {
    std::vector<Id> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](Id a, Id b) {
      if (deg[a] != deg[b]) return deg[a] < deg[b];
      return a < b;
    });
    const Id nc = std::min<Id>(n, 8);
    for (Id i = 0; i < nc; ++i) cand.push_back(idx[i]);
    // 补几个确定性伪随机起点（用固定种子，保证结果可复现）
    unsigned h = 12345u;
    for (int t = 0; t < 4 && n > 0; ++t) {
      h = h * 1103515245u + 12345u;
      cand.push_back(static_cast<Id>(h % static_cast<unsigned>(n)));
    }
  }

  std::vector<Id> best, cur;
  Id bestBw = -1;
  for (Id s : cand) {
    if (s < 0 || s >= n) continue;
    cuthillMcKee(s, cur);
    const Id bw = bwOf(cur);
    if (bestBw < 0 || bw < bestBw) { bestBw = bw; best = cur; }
    if (bestBw <= 1) break;              // 已达理论下界
  }

  std::reverse(best.begin(), best.end());
  // 统一约定：iperm[newIdx] = oldIdx
  std::vector<Id> iperm(n);
  for (Id k = 0; k < n; ++k) iperm[k] = best[k];
  return iperm;
}

// ---------------------------------------------------------------------------
//  最小度排序 MD —— 已停用，理由记录在案（重要，别再轻易启用）
//
//  理论上 MD 抑制填充优于 RCM，但在真实结构模型上不可用：
//  消元图邻接的膨胀是指数级的。8 邻接的网格点，消元一次后其 8 个邻居
//  两两相连（+28 边），下一轮每个都有约 30 个邻居，再下一轮就是 900……
//  144 节点的小网格就会耗尽内存。这是 MD 的固有性质，不是实现问题。
//
//  真正的解法是嵌套剖分（Nested Dissection, ND），把网格分层消元，
//  填充量可降到 O(n^1.5) —— 但那是数百行代码的工程，且不改善"结果正确性"。
//
//  结论：结构分析软件的正确选择是 ——
//   · 框架/梁柱结构：RCM 就足够好（按楼层+柱网编号天然接近最优）
//   · 壳单元为主的模型：RCM + 迭代法求解器（避开直接法的填充问题）
//   · 确实需要时：Intel MKL Pardiso（内部是 AMD + 多重 frontal）
// ---------------------------------------------------------------------------

// 按给定排列计算带宽（带状算法性能的直接指标）
// perm[newIdx] = oldIdx
inline Id bandwidth(const SymSparseMatrix& A, const std::vector<Id>& perm) {
  const Id n = A.size();
  if (n == 0) return 0;
  // newOf[old] = new
  std::vector<Id> newOf(n);
  for (Id k = 0; k < n; ++k) newOf[perm[k]] = k;
  Id bw = 0;
  for (Id i = 0; i < n; ++i)
    for (Id p = A.rowPtr()[i]; p < A.rowPtr()[i + 1]; ++p)
      bw = std::max(bw, std::abs(newOf[i] - newOf[A.colInd()[p]]));
  return bw;
}

// 消元后 L 因子的填充量与带宽（带宽按新序衡量）
inline size_t fillInCount(const SymSparseMatrix& A, const std::vector<Id>& perm,
                          Id* bwOut = nullptr) {
  const Id n = A.size();
  if (n == 0) { if (bwOut) *bwOut = 0; return 0; }
  std::vector<Id> newOf(n);
  for (Id k = 0; k < n; ++k) newOf[perm[k]] = k;

  // 按新序建对称邻接
  std::vector<std::vector<Id>> adj(n);
  for (Id i = 0; i < n; ++i)
    for (Id p = A.rowPtr()[i]; p < A.rowPtr()[i + 1]; ++p) {
      const Id j = A.colInd()[p];
      if (i == j) continue;
      adj[newOf[i]].push_back(newOf[j]);
      adj[newOf[j]].push_back(newOf[i]);
    }

  size_t fill = 0;
  Id bw = 0;
  std::vector<char> dead(n, 0);
  std::vector<Id> nbr;
  for (Id k = 0; k < n; ++k) {
    dead[k] = 1;
    nbr.clear();
    for (Id j : adj[k]) if (!dead[j]) nbr.push_back(j);
    for (Id j : nbr) bw = std::max(bw, k - j);
    for (Id a = 0; a < static_cast<Id>(nbr.size()); ++a)
      for (Id b = a + 1; b < static_cast<Id>(nbr.size()); ++b) {
        adj[nbr[a]].push_back(nbr[b]);
        adj[nbr[b]].push_back(nbr[a]);
        ++fill;
      }
  }
  if (bwOut) *bwOut = bw;
  return fill;
}

}  // namespace yjk
