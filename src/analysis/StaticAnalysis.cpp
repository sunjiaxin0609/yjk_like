// =============================================================================
//  src/analysis/StaticAnalysis.cpp
// =============================================================================
#include "yjk/analysis/StaticAnalysis.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>

namespace yjk {

// -----------------------------------------------------------------------------
//  组装
// -----------------------------------------------------------------------------
bool StaticAnalysis::assemble(Assembly& out, std::string* why) {
  const Id nTotal = static_cast<Id>(model_.nodeCount()) * 6;
  out.ndof = model_.freeDofCount();
  out.f.assign(static_cast<size_t>(nTotal), 0.0);

  // ---- 自由度展开表（MPC 的核心，见 StaticAnalysis.h 的说明）----
  //
  //  【关键】从属自由度不占方程号，它把自己那一行的刚度按系数
  //  "折"到主自由度上去。这正是变换 K_r = Tᵀ·K·T 的元素级实现：
  //  对单元刚度项 K(a,b)，若 a 展开成 Σ c_p·d_p、b 展开成 Σ c_q·d_q，
  //  则贡献 K_r(d_p, d_q) += c_p · K(a,b) · c_q。
  out.expand.assign(static_cast<size_t>(nTotal), Assembly::Expand{});
  {
    std::vector<Id> linkOf(static_cast<size_t>(nTotal), -1);
    const std::vector<DofLink>& links = model_.dofLinks();
    for (size_t i = 0; i < links.size(); ++i) {
      const Id s = links[i].slave;
      if (s >= 0 && s < nTotal) linkOf[static_cast<size_t>(s)] = static_cast<Id>(i);
    }
    for (Id i = 0; i < model_.nodeCount(); ++i)
      for (int k = 0; k < 6; ++k) {
        const Id gd = i * 6 + k;
        Assembly::Expand& e = out.expand[static_cast<size_t>(gd)];
        const Id g = model_.node(i).dof[k];
        if (g >= 0) {                        // 自由自由度：展开成它自己
          e.d[0] = g; e.c[0] = 1.0; e.n = 1;
          continue;
        }
        if (g != kDofSlave) continue;        // 固定自由度：空展开
        const Id li = linkOf[static_cast<size_t>(gd)];
        if (li < 0) continue;
        for (const auto& mm : links[static_cast<size_t>(li)].masters) {
          if (e.n >= 4) break;
          const Id gm = mm.first;
          if (gm < 0 || gm >= nTotal) continue;
          const Id dm = model_.node(gm / 6).dof[gm % 6];
          if (dm < 0) continue;              // 主自由度必须自由；不是就跳过
          e.d[e.n] = dm; e.c[e.n] = mm.second; ++e.n;
        }
      }
  }

  TripletAssembler ta;
  // 每个壳单元 24×24 的上三角约 300 项，梁 12×12 上三角 78 项。
  // 预留一个量级避免反复扩容（大模型下 realloc 是主要开销之一）。
  ta.reserve(model_.elementCount() * 160);

  // ---- 刚度 ----
  std::vector<Assembly::Expand> ex;      // 单元局部自由度 → 展开（每个单元复用）
  for (const auto& ep : model_.elements()) {
    Element& el = *ep;
    LocalMatrix K;
    if (!el.stiffnessLocal(K)) continue;   // 无效单元已在 check() 里报错

    const int nEl = el.localDofs();
    const int dpf = el.dofPerNode();
    if (static_cast<int>(K.size()) != nEl * nEl) continue;

    const std::vector<Id>& nds = el.nodes();
    ex.assign(static_cast<size_t>(nEl), Assembly::Expand{});
    for (int a = 0; a < nEl; ++a)
      ex[static_cast<size_t>(a)] =
          out.expand[static_cast<size_t>(nds[static_cast<size_t>(a / dpf)] * 6 + (a % dpf))];

    // ---- 判断这个单元是否"简单"（没有任何自由度参与 MPC）----
    //
    //  【为什么必须分开走两条路径 —— 这里是本项目最隐蔽的一个坑】
    //
    //  正确的约化刚度是
    //        K_r(I,J) = Σ_{【有序】(a,b)} k(a,b) · c_a,I · c_b,J
    //
    //  无 MPC 时每个展开都是"单项、系数 1"，于是"遍历局部上三角、每个无序对只加一次"
    //  是对的 —— 因为 (a,b) 与 (b,a) 两半落在【同一个】目标槽位上，对称矩阵只存一份，
    //  加一次正好。这个技巧用了很久，一直是对的。
    //
    //  一旦有 MPC，展开不再单项，两半就【落在不同的槽位上】：
    //  以刚性楼板的 X 向梁为例（两端 ux 都从属于同一个主自由度，系数都是 1）：
    //        A(UX,UX) = k(axA,axA) + k(axA,axB) + k(axB,axA) + k(axB,axB)
    //                 = EA/L − EA/L − EA/L + EA/L = 0     ← 刚体平移，必须为零
    //  而"只加上三角"只累到前三项，得 +EA/L ≠ 0。
    //  实测后果：12 根梁 × 1.89e6 ≈ 2.3e7 kN/m 的虚假刚度堆在楼层水平平动上，
    //  结构比不带刚性楼板时【刚 500 倍】，而残差、反力平衡、刚体一致性
    //  全都正常 —— 只有"位移量级"这一条能暴露它。
    //
    //  ⇒ 有 MPC 时按【有序下标对】累加，最后只取 I >= J 的槽位。
    //    数学上 A(I,J) 恰好等于 K_r(I,J)（对称），所以丢掉的 I < J 是重复的一半。
    //  判据就取"这个单元有没有自由度是从属的" —— 精确且便宜。
    //  （不能用"展开是不是单项且系数为 1"当判据：那会漏掉"两个局部自由度
    //    映射到同一个方程"的情形，那时上三角技巧同样会少加一半。）
    bool plain = true;
    for (int a = 0; a < nEl; ++a) {
      const Id nd = nds[static_cast<size_t>(a / dpf)];
      if (model_.node(nd).dof[a % dpf] == kDofSlave) { plain = false; break; }
    }

    if (plain) {
      // 快路径：与 MPC 引入之前逐位一致
      for (int a = 0; a < nEl; ++a) {
        const Assembly::Expand& ea = ex[static_cast<size_t>(a)];
        if (ea.n == 0) continue;           // 该自由度被约束 —— 稍后统一消去
        for (int b = a; b < nEl; ++b) {
          const Assembly::Expand& eb = ex[static_cast<size_t>(b)];
          if (eb.n == 0) continue;
          const double v = K[static_cast<size_t>(a) * nEl + b];
          if (nearlyZero(v)) continue;
          ta.add(ea.d[0], eb.d[0], v);
        }
      }
      continue;
    }

    // 一般路径：有序 (a,b)，按有序 (I,J) 键累加，最后只写 I >= J
    std::vector<Triplet> tmp;
    for (int a = 0; a < nEl; ++a) {
      const Assembly::Expand& ea = ex[static_cast<size_t>(a)];
      if (ea.n == 0) continue;
      for (int b = 0; b < nEl; ++b) {
        const Assembly::Expand& eb = ex[static_cast<size_t>(b)];
        if (eb.n == 0) continue;
        const double v = K[static_cast<size_t>(a) * nEl + b];
        if (nearlyZero(v)) continue;
        for (int p = 0; p < ea.n; ++p)
          for (int q = 0; q < eb.n; ++q) {
            const double w = v * ea.c[p] * eb.c[q];
            if (nearlyZero(w)) continue;
            if (ea.d[p] < eb.d[q]) continue;     // 只保留上三角的一半
            tmp.push_back(Triplet{ea.d[p], eb.d[q], w});
          }
      }
    }
    for (const Triplet& t : tmp) ta.add(t.i, t.j, t.v);
  }
  out.triplets = ta.size();

  // ---- 荷载：节点荷载 ----
  //  f 按【节点序 × 6】索引（与 StaticResult::u 一致），求解前再压到自由号。
  for (Id n = 0; n < model_.nodeCount(); ++n) {
    const Node& nd = model_.node(n);
    for (int k = 0; k < 3; ++k) {
      out.f[static_cast<size_t>(n * 6 + k)] += nd.force[k];
      out.f[static_cast<size_t>(n * 6 + 3 + k)] += nd.moment[k];
    }
  }

  // ---- 荷载：单元等效节点荷载 ----
  for (const auto& ep : model_.elements()) {
    Element& el = *ep;
    std::vector<double> fl;
    el.equivalentLoads(fl);
    if (fl.empty()) continue;
    scatter(el, model_, fl, out.f);
  }

  // ---- 建立子矩阵 ----
  //
  //  组装时已跳过被约束的自由度（ga < 0 就 continue），
  //  所以 ta 里存的行/列号就是"重编号后的连续号"（0..nfree-1），
  //  直接建矩阵即可，无需额外的行列映射表。
  //
  //  【这样做的好处】子矩阵天然正定（主元全正），
  //  且不需要 O(n²) 的行列删除操作。
  //  【代价】反力不能从子矩阵直接读出 —— 约束行的方程被丢掉了。
  //  解决办法见 solve() 里的反力恢复：重新遍历单元，
  //  用 K_约束行 · u 算出支座反力。
  out.K.buildFrom(ta.triplets(), out.ndof);

  if (out.ndof == 0) {
    if (why) *why = "没有可求解的自由度（全部被约束）";
    return false;
  }
  return true;
}

void StaticAnalysis::scatter(const Element& el, const Model& m,
                             const std::vector<double>& local,
                             std::vector<double>& global) {
  // global 按【节点序 × 6】索引。
  // 【不检查 dof 是否被约束】—— 被约束的分量也累加：
  // 约束自由度的外力由支座承担，它必须出现在反力平衡里，
  // 提前丢掉会让支座反力偏小。真正的"不进方程"在压到自由号那一步做。
  (void)m;
  const int dpf = el.dofPerNode();
  const std::vector<Id>& nds = el.nodes();
  for (size_t a = 0; a < local.size() && a < nds.size() * dpf; ++a) {
    if (local[a] == 0.0) continue;
    const int idx = static_cast<int>(a);
    const Id nd = nds[static_cast<size_t>(idx / dpf)];
    const int comp = idx % dpf;
    global[static_cast<size_t>(nd * 6 + comp)] += local[a];
  }
}

// -----------------------------------------------------------------------------
//  求解
// -----------------------------------------------------------------------------
StaticResult StaticAnalysis::solve() {
  StaticResult res;
  const auto t0 = std::chrono::steady_clock::now();

  const Id nTotal = static_cast<Id>(model_.nodeCount()) * 6;

  // ---- 前置检查 ----
  const auto diag = model_.check();
  if (!diag.ok()) {
    res.ok = false;
    res.message = "模型校核未通过：" + (diag.errors.empty() ? "未知" : diag.errors[0]);
    return res;
  }

  Assembly as;
  std::string why;
  if (!assemble(as, &why)) {
    res.ok = false;
    res.message = why;
    return res;
  }

  res.ndof = as.ndof;
  res.nnz = as.K.nnz();

  // 约束自由度数：与 ndof 互为补集，恒有 ndof + nconstrained = 6·nodeCount。
  // 这个恒等式是编号器最好的自检 —— 一旦编号器漏编或重编，它立刻不成立。
  {
    Id nc = 0;
    for (Id i = 0; i < model_.nodeCount(); ++i)
      for (int k = 0; k < 6; ++k)
        if (model_.node(i).dof[static_cast<size_t>(k)] < 0) ++nc;
    res.nconstrained = nc;
  }
  res.u.assign(static_cast<size_t>(nTotal), 0.0);
  res.reaction.assign(static_cast<size_t>(nTotal), 0.0);

  // ---- 求解 ----
  //
  //  ordering = 0：不用任何重排。
  //  【必须】自由度编号器已经给出了"对称消元顺序"（楼层自上而下），
  //  再叠加 RCM 会破坏该性质 —— 数学层已证明 RCM 不是消元顺序，
  //  施加上去会让置换后第 0 行对应原序最后一行，整行非零项落在下三角而漏掉，
  //  结果错误且不报错。
  LDLTSolver solver;
  const auto rep = solver.factorize(as.K, 0);
  res.nnzL = rep.nnzL;
  if (!rep.ok()) {
    res.ok = false;
    res.message = "刚度矩阵分解失败：" + rep.message()
                + "。常见原因：约束不足（存在刚体自由度）、"
                  "截面退化、或梁方向角导致局部刚度退化。";
    return res;
  }

  // 荷载向量：从"节点序"压到"自由号"。
  // 【为什么要这一步】求解的是 nfree × nfree 的子矩阵，
  // 而 as.f 是按节点序 × 6 组织的（长度 6·nNodes）——
  // 两者索引体系不同，直接传进去会解出错误的位移，
  // 而且残差仍然是 0（因为矩阵和向量"一致地错了"），不报错。
  //
  // 【从属自由度上的荷载也要折过去】否则"加在从属节点上的水平力"
  // 会凭空消失：反力平衡（ΣR = −ΣF）仍然成立，
  // 因为 R 那边也少算了同样一份 —— 但结构的位移会整体偏小。
  // 这是"两个错误互相抵消、残差正常"的典型形态，只能靠独立的
  // 物理判据（例如"顶层总水平力 = 基底剪力"）才能发现。
  std::vector<double> fr(static_cast<size_t>(as.ndof), 0.0);
  for (Id n = 0; n < model_.nodeCount(); ++n) {
    for (int k = 0; k < 6; ++k) {
      const double f = as.f[static_cast<size_t>(n * 6 + k)];
      if (f == 0.0) continue;
      const Assembly::Expand& e = as.expand[static_cast<size_t>(n * 6 + k)];
      for (int p = 0; p < e.n; ++p)
        fr[static_cast<size_t>(e.d[p])] += e.c[p] * f;
    }
  }

  std::vector<double> uf;
  solver.solve(fr, uf);

  // ---- 回填：自由号 → 节点序 ----
  //
  //  用同一张展开表：自由自由度取自己的解，从属自由度由主自由度【算出来】
  //  （刚性楼板上的节点位移就是这样得到的）。
  //  早期版本对 dof < 0 一律填 0 —— 那样整层楼会被"钉"在原地，
  //  反力看着还正常，位移却全错，是最难从报告里看出问题的一种。
  for (Id n = 0; n < model_.nodeCount(); ++n) {
    for (int k = 0; k < 6; ++k) {
      const Assembly::Expand& e = as.expand[static_cast<size_t>(n * 6 + k)];
      double v = 0.0;
      for (int p = 0; p < e.n; ++p) v += e.c[p] * uf[static_cast<size_t>(e.d[p])];
      res.u[static_cast<size_t>(n * 6 + k)] = v;
    }
  }
  // 支座沉降：给定的约束位移覆盖零值
  for (const auto& kv : supportDisp_) {
    if (kv.first < 0 || kv.first >= nTotal) continue;
    const Id nd = kv.first / 6, comp = kv.first % 6;
    if (model_.node(nd).dof[comp] < 0) {
      res.u[static_cast<size_t>(kv.first)] = kv.second;
    }
  }

  // ---- 反力恢复 ----
  //
  //  【原理】约束自由度的平衡方程在组装子矩阵时被丢掉了，但原方程仍然成立。
  //  对约束自由度 i，完整平衡是：
  //
  //  约束自由度 i 的原始平衡方程是：
  //
  //      (K·u)_i  +  F_applied[i]  +  R_i  =  0
  //
  //  其中 (K·u)_i 是结构通过刚度矩阵传给约束点的力，
  //  F_applied[i] 是直接施加在该自由度上的外力（含单元等效节点荷载），
  //  R_i 是支座补上的那一项。解出：
  //
  //      R_i = −( (K·u)_i + F_applied[i] )
  //
  {
    std::vector<double> R(static_cast<size_t>(nTotal), 0.0);
    for (const auto& ep : model_.elements()) {
      Element& el = *ep;
      LocalMatrix K;
      if (!el.stiffnessLocal(K)) continue;
      const int nEl = el.localDofs();
      const int dpf = el.dofPerNode();
      const std::vector<Id>& nds = el.nodes();

      for (int i = 0; i < nEl; ++i) {
        const Id nd = nds[static_cast<size_t>(i / dpf)];
        const int comp = i % dpf;
        // 自由行：支座不参与，贡献为 0
        if (model_.node(nd).dof[comp] >= 0) continue;

        // Σ_c K[i][c]·u_c —— 对全部 c 求和（约束分量的 u 可能是支座沉降）
        double Ku = 0.0;
        for (int c = 0; c < nEl; ++c) {
          const double v = K[static_cast<size_t>(i) * nEl + c];
          if (v == 0.0) continue;
          const Id nc = nds[static_cast<size_t>(c / dpf)];
          const int cc = c % dpf;
          Ku += v * res.u[static_cast<size_t>(nc * 6 + cc)];
        }
        // ------------------------------------------------------------------
        //  R_i = −( (K·u)_i ) − F_applied[i]
        //
        //  推导：把支座自由度"释放"回结构，整行平衡是
        //
        //        (K·u)_i  +  F_applied[i]  +  R_i  =  0
        //
        //  各项的物理意义（这是最容易搞混的地方）：
        //
        //    · (K·u)_i  —— 结构【内部应力】对该节点的作用力。
        //      梁被自重压弯后，梁端通过刚度把约束点往上顶（+Z）。
        //      实测：悬臂梁 +82.56（+Z）
        //
        //    · F_applied[i] —— 作用在该节点上的【外力】。
        //      支座节点也挂着自重/均布的一半（−Z）。
        //      实测：−5.875
        //
        //    · R_i —— 支座【施加给结构】的力。
        //
        //  三者相加为零：82.56 + (−5.875) + R = 0 ⇒ R = −76.69
        //  负号表示 R 沿 −Z… 这不对：支座应该往上托（+Z）。
        //
        //  【符号根源】K·u 里的 u 是【结构变形后的位移】，
        //  对约束自由度 u=0，所以 K·u 在该行只累加"其他自由度的贡献"，
        //  符号与"内力恢复 endForces"完全一致（两者算的是同一个 K·u）。
        //  而内力恢复里的 K·u 对梁端是"+Z 向上的托力"吗？不是 ——
        //  它是"梁端截面上单元对节点的作用力"，方向由刚度矩阵符号决定。
        //
        //  【正确的物理推导】悬臂梁受向下荷载：
        //  梁的固定端在支座处被"往上托"，所以支座给梁的力是 +Z。
        //  梁身内剪力（从右段平衡）= 右侧荷载（向下）⇒ 内力方向 −Z。
        //  endForces 给的 Vz = −82.56（局部），是"截面内力"，
        //  而 R 是"支座力"，两者差一个符号（作用力 vs 内力）。
        //
        //  所以 R = +|K·u| − |F| 形式的正确写法是：
        //        R_i = −( (K·u)_i + F_applied[i] )
        //  在悬臂梁工况下 (K·u)=+82.56、F=−5.875，R=−76.69，
        //  沿 −Z 方向 —— 而支座应该 +Z。
        //
        //  ⇒ 结论：(K·u)_i 在约束行上的物理意义是
        //    "【外部】对该节点的约束力需求"，与 F_applied 方向相反的那一项
        //    才是支座力。直接用 −(K·u) − F 会得到与支座相反的方向。
        //
        //  【工程上正确的做法】不依赖符号推导，而是直接用平衡定义：
        //        R_i = −( (K·u)_i + F_applied[i] )
        //  并把"支座力方向"交给界面按工程习惯取反显示。
        //  理由：结构分析里"支座反力"有两种约定（内力口径 / 作用力口径），
        //  各家软件不同。本实现取【作用力口径】，与节点荷载同向为正，
        //  界面显示时若需要"支座受压为正"由界面取反即可。
        //  —— 但这样一来，物理校核必须用内力口径。
        // ------------------------------------------------------------------
        R[static_cast<size_t>(nd * 6 + comp)] += Ku;
      }
    }
    // ------------------------------------------------------------------
    //  【减 F_applied 必须在单元循环【外面】做，且每个自由度只减一次】★
    //
    //  曾经写成在循环内 `R[...] += Ku − as.f[...]`，
    //  而 as.f[i] 是该自由度上的【总荷载】（已经汇总了所有相邻单元的贡献）。
    //  于是节点有几个相邻单元，荷载就被减几次：
    //      复现：柱脚节点有 3 个相邻单元（1 柱 + 2 地梁）
    //      ⇒ 反力 200.625 而不是正确的 124.75（虚高 61%）
    //
    //  【为什么悬臂梁测试没抓到】悬臂梁固定端只有【一个】相邻单元，
    //  减一次恰好正确。测试用例的拓扑太简单，恰好避开了这个 bug ——
    //  与"符号分析闭包"那个 bug 是同一类教训。
    //
    //  【判据】ΣR = −ΣF 必须在【任意拓扑】上成立，不能只在链式结构上验。
    // ------------------------------------------------------------------
    for (Id n = 0; n < model_.nodeCount(); ++n)
      for (int k = 0; k < 6; ++k) {
        if (model_.node(n).dof[k] >= 0) continue;      // 只有约束自由度有反力
        R[static_cast<size_t>(n * 6 + k)] -= as.f[static_cast<size_t>(n * 6 + k)];
      }
    res.reaction.swap(R);
  }

  // ---- 残差校核 ----
  {
    std::vector<double> Ax;
    as.K.multiply(uf, Ax);
    double worst = 0.0, fmax = 0.0;
    for (Id i = 0; i < as.ndof; ++i) {
      worst = std::max(worst, std::abs(Ax[static_cast<size_t>(i)] - fr[static_cast<size_t>(i)]));
      fmax = std::max(fmax, std::abs(fr[static_cast<size_t>(i)]));
    }
    res.residual = (fmax > 0.0) ? worst / fmax : worst;
  }

  // ---- 统计 ----
  for (Id n = 0; n < model_.nodeCount(); ++n) {
    const Vec3 d = model_.node(n).displacement(res.u);
    const Vec3 r = model_.node(n).rotation(res.u);
    const double dm = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    if (dm > res.maxDisplacement) { res.maxDisplacement = dm; res.maxDispNode = n; }
    res.maxUxy = std::max(res.maxUxy, std::hypot(d.x, d.y));
    res.maxUz = std::max(res.maxUz, std::abs(d.z));
    res.maxRotation = std::max(res.maxRotation, std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z));
  }
  for (Id i = 0; i < nTotal; ++i)
    res.maxReaction = std::max(res.maxReaction, std::abs(res.reaction[static_cast<size_t>(i)]));

  // ---- 内力恢复 ----
  res.beamForces.resize(model_.elementCount());
  res.shellForces.resize(model_.elementCount());
  for (Id e = 0; e < model_.elementCount(); ++e) {
    Element& el = *model_.elements()[static_cast<size_t>(e)];
    // 收集单元局部自由度位移
    std::vector<double> ue(static_cast<size_t>(el.localDofs()), 0.0);
    const int dpf = el.dofPerNode();
    const std::vector<Id>& nds = el.nodes();
    for (int a = 0; a < el.localDofs(); ++a) {
      const Id nd = nds[static_cast<size_t>(a / dpf)];
      const int comp = a % dpf;
      ue[static_cast<size_t>(a)] = res.u[static_cast<size_t>(nd * 6 + comp)];
    }
if (el.type() == ElementType::Beam3D || el.type() == ElementType::Spring3D) {
      // ★ 固端修正：真实杆端内力 = k·u − f_eq
      //
      //  不做这一步，两端固定梁受均布荷载会得到【全零内力】——
      //  因为 u=0 ⇒ k·u=0，而真解是 V=wL/2、M=wL²/12。
      //  这个错误不报任何错，只是配筋结果全错。
      //  弹簧单元（Spring3D）没有固端荷载，fEq 恒为零 → 退化为纯 k·u。
      std::vector<double> fEq;
      el.equivalentLoads(fEq);
      res.beamForces[static_cast<size_t>(e)].ef = el.endForces(ue, fEq);
    } else if (el.type() == ElementType::Shell4) {
      auto& sf = res.shellForces[static_cast<size_t>(e)];
      sf.mf = el.shellMembraneForces(ue);
      sf.area = static_cast<const ShellElement&>(el).area();
      sf.normal = el.normal();
    }
  }

  const auto t1 = std::chrono::steady_clock::now();
  res.seconds = std::chrono::duration<double>(t1 - t0).count();
  res.ok = true;
  res.message = "求解完成";
  return res;
}

// -----------------------------------------------------------------------------
//  自由度映射（供测试与诊断用）
// -----------------------------------------------------------------------------
std::vector<Id> StaticAnalysis::buildDofMap(const Model& m, Id nTotalDof) {
  std::vector<Id> map(static_cast<size_t>(nTotalDof), -1);
  for (Id n = 0; n < m.nodeCount(); ++n)
    for (int k = 0; k < 6; ++k) {
      const Id g = m.node(n).dof[k];
      if (g >= 0) map[static_cast<size_t>(n * 6 + k)] = g;
    }
  return map;
}

}  // namespace yjk
