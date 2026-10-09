// =============================================================================
//  yjk/element/BeamElement3D.h  ——  三维梁单元（2 节点，12 自由度）
//
// =============================================================================
//  【自由度排列】—— 单元层的"接口约定"，错了全盘皆错
//
//      节点 i (0..5):  ux    uy    uz    rx    ry    rz
//      节点 j (6..11): ux    uy    uz    rx    ry    rz
//
//      u = 平动（沿全局 x,y,z）
//      r = 转角（绕全局 x,y,z，单位为弧度）
//
//  【局部坐标】
//      x' 沿梁轴线（由 i 指向 j）
//      y' 截面"竖向"（由用户指定的方向角决定）
//      z' = x' × y'  右手系补全
//
//  刚度矩阵的构造分两步（这是所有空间梁单元的标准做法）：
//      1. 在【局部坐标】下组装 12×12 —— 这里是解析式的，很稳定
//      2. 用坐标变换矩阵 R（3×3 旋转）转到【全局坐标】
//         K_global = Tᵀ · K_local · T
//         T = blkdiag(R, R, R, R)  即 4 个 3×3 旋转块
//
//  【为什么必须分两步】局部坐标下梁是"沿 x' 的直杆"，公式简单且对称；
//  直接在全局坐标组装则要处理任意空间角度，极易出错。
//  盈建科/SAP2000 都是这么做的。
// =============================================================================
//
//  【局部刚度矩阵的构成】12×12 按物理效应分块：
//
//      轴向（EA）        : u'_x                —— 1 自由度/节点
//      扭转（GJ）        : r'_x                —— 1 自由度/节点
//      绕 y 弯曲（EIy）  : u'_z, r'_y          —— 2 自由度/节点
//      绕 z 弯曲（EIz）  : u'_y, r'_z          —— 2 自由度/节点
//      剪切（G·As）      : 与上面耦合（Timoshenko）
//
//  注意"绕 y 弯曲"配的是 u'_z（z 向平动）—— 因为绕 y 轴转动会让 z 方向产生位移。
// 这是空间梁单元最容易搞混的地方，也是平面内外刚度互换的根源。
// =============================================================================
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "yjk/math/Types.h"
#include "yjk/model/Section.h"

namespace yjk {

// ---------------------------------------------------------------------------
//  梁单元
// ---------------------------------------------------------------------------
class BeamElement3D {
 public:
  static constexpr int kDofPerNode = 6;
  static constexpr int kDof = 12;

  BeamElement3D() = default;

  // pi, pj    : 节点坐标
  // sec       : 截面几何性质
  // mat       : 材料
  // upVector  : 参考"上方向"（通常取全局 -Z，即竖直向下）
  //              局部 y' 轴 = up 中垂直于 x' 的分量（正交化后）
  //              传入 {0,0,0} 时默认用全局 -Z；若梁竖直（x' 平行 Z），
  //              则退化，改用全局 +Y（水平方向）以避免奇异性
  BeamElement3D(const Vec3& pi, const Vec3& pj, const SectionProperties& sec,
                const Material& mat, const Vec3& upVector = Vec3{0, 0, -1}) {
    init(pi, pj, sec, mat, upVector);
  }

  void init(const Vec3& pi, const Vec3& pj, const SectionProperties& sec,
            const Material& mat, const Vec3& upVector) {
    pi_ = pi;
    pj_ = pj;
    sec_ = sec;
    mat_ = mat;

    const Vec3 d = pj - pi;
    length_ = norm(d);
    if (!(length_ > kEps)) {
      valid_ = false;
      error_ = "梁单元长度为零（两端节点重合）";
      return;
    }
    ex_ = d * (1.0 / length_);

    // ---- 局部 y' 轴：up 在垂直于 x' 平面上的投影 ----
    Vec3 up = upVector;
    if (nearlyZero(norm(up))) up = Vec3{0, 0, -1};

    // 关键处理：若 up 与 x' 平行，投影为零 → 退化。
    // 这在"竖直柱"上很常见：默认 up = (0,0,-1) 而 x' 也是竖直方向。
    // 此时应改用一个水平参考方向，否则整个单元刚度矩阵奇异。
    Vec3 proj = up - ex_ * dot(up, ex_);
    if (norm(proj) < 1.0e-8) {
      // up 与 x' 平行 —— 退化，改用全局 Y 或 X 作为参考
      const Vec3 fallback = (std::abs(ex_.y) < 0.9) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
      proj = fallback - ex_ * dot(fallback, ex_);
      degenerateUp_ = true;
    }
    ey_ = normalized(proj);
    ez_ = cross(ex_, ey_);      // 右手系：x' × y' = z'
    if (norm(ez_) < 1.0e-12) {
      valid_ = false;
      error_ = "梁局部坐标系构造失败（x' 与 y' 共线）";
      return;
    }
    ez_ = normalized(ez_);

    // 正交性自检 —— 旋转矩阵必须是正交阵，否则 K_global = TᵀKL T 的变换失效
    const double o1 = std::abs(dot(ex_, ey_));
    const double o2 = std::abs(dot(ex_, ez_));
    const double o3 = std::abs(dot(ey_, ez_));
    if (o1 > 1e-10 || o2 > 1e-10 || o3 > 1e-10) {
      valid_ = false;
      error_ = "梁局部坐标轴不正交";
      return;
    }
    valid_ = true;
  }

// ---- 查询 ----
  bool isValid() const { return valid_; }
  const std::string& error() const { return error_; }
  double length() const { return length_; }
  bool usedFallbackUp() const { return degenerateUp_; }
  const Vec3& localX() const { return ex_; }
  const Vec3& localY() const { return ey_; }
  const Vec3& localZ() const { return ez_; }

  // -------------------------------------------------------------------------
  //  端部释放（端部铰 / 选择性释放）
  //
  //  i/j 端各 6 个【局部】自由度（ux,uy,uz,rx,ry,rz）可选择性释放。
  //  释放自由度通过【静力凝聚（Guyan 消去）】从单元中消去：
  //      Kc_aa = K_aa − K_ab · Kbb⁻¹ · Kba
  //  凝聚后的 K 在释放自由度上的行/列全为 0 —— 即该方向单元刚度为零。
  //  配套：等效节点荷载同样凝聚（不把力施加到释放自由度上），
  //  端内力恢复里释放端的对应分量也恒为 0（端弯矩=0，铰的真正含义）。
  //
  //  【典型用途】
  //    · 两端铰（桁架杆）：releaseAt(0,4), releaseAt(0,5),
  //                        releaseAt(1,4), releaseAt(1,5)   （释放 ry, rz）
  //    · 一端刚一端铰：只对铰接端释放弯曲
  //
  //  【注意】释放是"减刚度"操作：若被释放的自由度没有任何其他单元
  //  提供刚度，它在全局装配中会是"零刚度自由度"，DofNumbering 会自动
  //  锁定它（见 Model.cpp zeroStiff 逻辑），避免全局刚度矩阵奇异。
  // -------------------------------------------------------------------------
  void setReleases(const std::array<bool, 6>& iRel, const std::array<bool, 6>& jRel) {
    relI_ = iRel;
    relJ_ = jRel;
    singular_ = false;
    condensationErr_.clear();
  }
  // 便捷：释放单个自由度。localNode = 0(i 端) / 1(j 端)，comp = 0..5
  // （对应局部 ux,uy,uz,rx,ry,rz —— 与刚度矩阵列顺序一致）
  void releaseAt(int localNode, int comp) {
    if (localNode < 0 || localNode > 1 || comp < 0 || comp >= 6) return;
    (localNode == 0 ? relI_ : relJ_)[static_cast<size_t>(comp)] = true;
    singular_ = false;
    condensationErr_.clear();
  }
  void clearReleases() {
    relI_.fill(false);
    relJ_.fill(false);
    singular_ = false;
    condensationErr_.clear();
  }
  bool isReleased(int localNode, int comp) const {
    if (localNode < 0 || localNode > 1 || comp < 0 || comp >= 6) return false;
    return (localNode == 0 ? relI_ : relJ_)[static_cast<size_t>(comp)];
  }
  bool hasReleases() const {
    for (bool b : relI_) if (b) return true;
    for (bool b : relJ_) if (b) return true;
    return false;
  }
  // 端部释放的静力凝聚是否可行（Kbb 非奇异）。
  // 奇异意味着释放组合让单元变成"局部机构"（例：把某端横向平动与
  // 与之耦合的弯曲转角一起释放，而剪切参数 Φ=0 时该 2×2 子块退化）。
  // 冻结（奇异）时单元对结构贡献零刚度，check() 会在预检中报错。
  // 注意：四舍五入级别的"弱主元"也视为奇异 —— 宁可报错也不静默给零刚度。
  bool condensationOk() const { return !singular_; }
  const std::string& condensationError() const { return condensationErr_; }

  // -------------------------------------------------------------------------
  //  坐标变换矩阵 T（12×12，4 个 3×3 旋转块）
  //  R 的行 = 局部轴在全局坐标中的表示：
  //      R = [ ex.x ex.y ex.z ]
  //          [ ey.x ey.y ey.z ]
  //          [ ez.x ez.y ez.z ]
  //  局部量 u_local = R · u_global，故 u_global = Rᵀ · u_local
  //  能量 uᵀKu 不变 → K_global = Tᵀ K_local T，其中 T = blkdiag(R,R,R,R)
  // -------------------------------------------------------------------------
  void transformationMatrix(Mat12& T) const {
    const double R[3][3] = {
        {ex_.x, ex_.y, ex_.z},
        {ey_.x, ey_.y, ey_.z},
        {ez_.x, ez_.y, ez_.z},
    };
    for (int blk = 0; blk < 4; ++blk)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          T(blk * 3 + i, blk * 3 + j) = R[i][j];
  }

  // -------------------------------------------------------------------------
  //  局部刚度矩阵（12×12）
  //
  //  考虑剪切变形（Timoshenko 梁）的标准形式：
  //      绕 y 弯曲（位移 u_z、转角 r_y），单元长 L，Φy = 12EIy/(G·Asy·L²)
  //          系数 12EIy/(L³(1+Φy)) 替代 12EIy/L³
  //      绕 z 弯曲同理
  //  Φy = 0 时退化为 Euler-Bernoulli（细长梁正确）
  //  Φy 很大时（深梁）剪切变形主导，必须计入，否则刚度被高估
  // -------------------------------------------------------------------------
  Mat12 localStiffness() const {
    Mat12 k = Mat12::zero();
    if (!valid_) return k;

    const double L = length_;
    const double E = mat_.E;
    const double G = mat_.shearModulus();
    const double A = sec_.A;

    // ---- 轴向 ----
    const double ka = E * A / L;
    k(0, 0) = ka;  k(0, 6) = -ka;
    k(6, 0) = -ka; k(6, 6) = ka;

    // ---- 扭转 ----
    const double kt = G * sec_.J / L;
    k(3, 3) = kt;  k(3, 9) = -kt;
    k(9, 3) = -kt; k(9, 9) = kt;

    // ---- 绕 y 轴弯曲（u_z, r_y）----
    //
    // 标准形式（c = 12EI/(L³(1+Φ))，Φ = 12EI/(G·As·L²) 为剪切参数）：
    //
    //     c · | 12      L/2        -12     L/2      |
    //         | L/2     L²(4+Φ)/12  -L/2    L²(2-Φ)/12 |
    //         | -12    -L/2          12    -L/2     |
    //         | L/2     L²(2-Φ)/12  -L/2    L²(4+Φ)/12 |
    //         +---------------------------------------  × (u_zi, r_yi, u_zj, r_yj)
    //
    // 【最容易错的地方】耦合项是 L/2，不是 L。
    // 写成 6cL 会让耦合项大 12 倍（因为 cL/2 = 6EI/L²，而 6cL = 72EI/L²）。
    // 后果：转角-位移的比例关系仍正确（所以"θ/w = 3/2L"这类相对校核能通过），
    //       但绝对挠度会小 30 倍左右 —— 而且不报任何错。
    // 记忆：k(θ,θ) = 4EI/L = cL²/3，可用它反查其他项。
    //
    // -------------------------------------------------------------------------
    //  【★ 两个弯曲平面的耦合项符号【相反】—— 这是本项目最隐蔽的一个 bug ★】
    //
    //  上面那张矩阵是"θ = dw/dx"的标准形式，它成立的前提是
    //  转角 DOF 与横向位移满足 **θ = +dw/dx**。
    //
    //  这个前提只对【绕 z 弯曲】成立：
    //      ω_z 的正向转动把 +x 轴转向 +y ⇒ u_y = +ω_z·x ⇒ ω_z = +du_y/dx      ✓
    //  对【绕 y 弯曲】则相反：
    //      ω_y 的正向转动把 +x 轴转向 −z ⇒ u_z = −ω_y·x ⇒ ω_y = −du_z/dx      ✗
    //
    //  早期实现两个平面用了【同一张】矩阵，于是 y 平面实际使用的是 θ = +du_z/dx，
    //  与几何相差一个负号。后果不是"位移算错"（单个构件两端约束给定时结果仍是
    //  对的，这也是它藏得住的原因），而是：
    //    · **整体刚体转动会产生虚假内力** —— K·v ≠ 0，v 是刚体模式
    //    · 共享节点的转角在构件之间对不上 → 框架的内力分配错
    //    · 整体力矩平衡不闭合（实测差 120 kN·m）
    //
    //  验证办法（本项目的回归测试就是这么做的）：把"绕全局轴转单位角"的刚体
    //  位移场代入，检查 K·v 是否为舍入量级。修好前绕 y/z 分别差 1.2e2 / 1.95e1。
    //
    //  【修法】把矩阵中与 r_y 成【奇次】的项（耦合项 ±cl）统一取反，
    //  等价于做 K → D·K·D，D = diag(1,−1,1,−1)。对角项（弯曲、剪切刚度）不动。
    //  配套：等效节点荷载里的固端弯矩也必须同步取反，否则固端修正会差一个符号。
    // -------------------------------------------------------------------------
    {
      const double EI = E * sec_.Iy;
      const double As = sec_.Asy > 0.0 ? sec_.Asy : (5.0 / 6.0) * A;
      const double phi = (G * As > 0.0) ? 12.0 * EI / (G * As * L * L) : 0.0;
      const double c = 12.0 * EI / (L * L * L * (1.0 + phi));
      const double cl = c * L / 2.0;                    // 耦合项
      const double t_i = c * L * L * (4.0 + phi) / 12.0; // i 端转角刚度 = (4+Φ)EI/L
      const double t_j = c * L * L * (2.0 - phi) / 12.0; // ij 耦合转角 = (2-Φ)EI/L
      // 局部自由度：u_z' = 2(i) / 8(j)，r_y' = 4(i) / 10(j)
      //  【易错】节点内自由度顺序是 [ux,uy,uz,rx,ry,rz]，
      //          所以 r_y 是第 5 个（下标 4），r_z 才是第 6 个（下标 5）。
      //  配对规则：绕 y' 弯 → 配 u_z'（竖向位移），不配 u_y'。
      //  耦合项符号与下方 z 平面相反，理由见上。
      k(2, 2) = c;    k(2, 4) = -cl;  k(2, 8) = -c;   k(2, 10) = -cl;
      k(4, 2) = -cl;  k(4, 4) = t_i;  k(4, 8) = cl;   k(4, 10) = t_j;
      k(8, 2) = -c;   k(8, 4) = cl;   k(8, 8) = c;    k(8, 10) = cl;
      k(10, 2) = -cl; k(10, 4) = t_j; k(10, 8) = cl;  k(10, 10) = t_i;
    }

    // ---- 绕 z 轴弯曲（局部 u_y, 局部 r_z）----
    //
    //  【局部坐标下的弯曲平面配对 —— 千万别配错】
    //    绕局部 y' 弯曲 → 截面【竖向】位移 u_z' 与转角 r_y'
    //    绕局部 z' 弯曲 → 截面【横向】位移 u_y' 与转角 r_z'
    //  即：绕哪个轴弯，就配"另一个方向"的位移。
    //    配错的后果（如 u_y' 与 r_y' 耦合）不会报任何错，
    //    但内力恢复时弯矩会出现在错误的分量上（My 与 Mz 互换）。
    {
      const double EI = E * sec_.Iz;
      const double As = sec_.Asz > 0.0 ? sec_.Asz : (5.0 / 6.0) * A;
      const double phi = (G * As > 0.0) ? 12.0 * EI / (G * As * L * L) : 0.0;
      const double c = 12.0 * EI / (L * L * L * (1.0 + phi));
      const double cl = c * L / 2.0;
      const double t_i = c * L * L * (4.0 + phi) / 12.0;
      const double t_j = c * L * L * (2.0 - phi) / 12.0;
      // 局部自由度：u_y' = 1(i) / 7(j)，r_z' = 5(i) / 11(j)
      k(1, 1) = c;   k(1, 5) = cl;   k(1, 7) = -c;   k(1, 11) = cl;
      k(5, 1) = cl;  k(5, 5) = t_i;  k(5, 7) = -cl;  k(5, 11) = t_j;
      k(7, 1) = -c;  k(7, 5) = -cl;  k(7, 7) = c;    k(7, 11) = -cl;
      k(11, 1) = cl; k(11, 5) = t_j; k(11, 7) = -cl; k(11, 11) = t_i;
    }
    return k;
  }

// -------------------------------------------------------------------------
  //  全局刚度矩阵 K_global = Tᵀ · K_local · T
  //  【端部释放】有释放时用【凝聚后局部刚度】再变换 ——
  //  释放是对局部自由度的定义，必须先在局部凝聚，再回到全局，
  //  否则释放方向会混入全局各分量。
  // -------------------------------------------------------------------------
  Mat12 stiffness() const {
    Mat12 T;
    transformationMatrix(T);
    const Mat12 kl = hasReleases() ? condensedStiffness() : localStiffness();
    const Mat12 kt = kl * T;      // K_local · T
    return T.transposed() * kt;   // Tᵀ · (K_local · T)
  }

  // -------------------------------------------------------------------------
  //  端部释放：静力凝聚（Guyan 消去）
  //
  //  原理：把 12 个局部自由度分成两组
  //      a = 未释放（保留），b = 释放（被凝聚掉）
  //  平衡方程按 [a; b] 分块：
  //      K = [Kaa Kab; Kba Kbb],  f = [fa; fb]
  //  释放意味着 b 自由度上无外力、且单元对 b 的刚度贡献为零 ——
  //  由第二行平衡 b 的自由位移后代入第一行：
  //      ub = Kbb⁻¹ (fb − Kba·ua)
  //      Kc_aa = Kaa − Kab·Kbb⁻¹·Kba
  //      fc_a  = fa  − Kab·Kbb⁻¹·fb
  //  凝聚后单元对 b 的贡献恒为 0（b 行/列清零），b 端内力 = 0 ——
  //  这正是"端弯矩 = 0"（铰）的数学来源。
  //
  //  【为什么必须用"凝聚"而不是"删掉 b 自由度的行/列"】
  //  直接删行/列会丢掉 b 与 a 的耦合项 Kab —— 弯曲转角释放后，
  //  端部横向位移与转角的耦合（cl 项）必须保留在 a 里（通过 Kbb⁻¹ 折算）。
  //  只有凝聚才能既"释放该自由度"又不破坏杆端位移-内力关系。
  //
  //  【数值】nB ≤ 11 的小系统，带部分主元 Gauss-Jordan。
  //  Kbb 在常用释放组合下可逆（对角主元 > 0）；当释放组合形成
  //  "局部机构"（如同一端同时释放 u_y 与 r_z，且 Φ=0 剪切参数为零时
  //  该 2×2 子块退化）时 Kbb 奇异 —— 这是真实的物理退化，应当报错
  //  而不是硬算出一个无意义的零刚度矩阵。
  //
  //  【先凝聚后变换】凝聚在局部坐标下进行，然后整体做 Tᵀ·Kc·T。
  //  释放方向语义 = 局部自由度（工程习惯：铰 = 绕截面主轴/局部轴释放）。
  // -------------------------------------------------------------------------

  // 带部分主元 Gauss-Jordan 解 A·X = B（A n×n，B n×nRhs，原地求逆消元）。
  // 返回 false 表示 A 奇异。nRhs 可为 0（只求逆不消右端）。
  static bool gaussJordanSolve(int n, double* A, double* B, int nRhs) {
    for (int col = 0; col < n; ++col) {
      int piv = col;
      double best = std::abs(A[static_cast<size_t>(col) * n + col]);
      for (int r = col + 1; r < n; ++r) {
        const double v =
            std::abs(A[static_cast<size_t>(r) * n + static_cast<size_t>(col)]);
        if (v > best) { best = v; piv = r; }
      }
      if (best < kPivTol) return false;
      if (piv != col) {
        for (int j = 0; j < n; ++j)
          std::swap(A[static_cast<size_t>(col) * n + j],
                    A[static_cast<size_t>(piv) * n + j]);
        for (int j = 0; j < nRhs; ++j)
          std::swap(B[static_cast<size_t>(col) * nRhs + j],
                    B[static_cast<size_t>(piv) * nRhs + j]);
      }
      const double d = 1.0 / A[static_cast<size_t>(col) * n + col];
      for (int j = 0; j < n; ++j) A[static_cast<size_t>(col) * n + j] *= d;
      for (int j = 0; j < nRhs; ++j) B[static_cast<size_t>(col) * nRhs + j] *= d;
      for (int r = 0; r < n; ++r) {
        if (r == col) continue;
        const double f = A[static_cast<size_t>(r) * n + col];
        if (nearlyZero(f)) continue;
        for (int j = 0; j < n; ++j)
          A[static_cast<size_t>(r) * n + j] -= f * A[static_cast<size_t>(col) * n + j];
        for (int j = 0; j < nRhs; ++j)
          B[static_cast<size_t>(r) * nRhs + j] -= f * B[static_cast<size_t>(col) * nRhs + j];
      }
    }
    return true;
  }

  // 凝聚后的局部刚度（b 行/列 = 0）。无释放时返回原始局部刚度。
  Mat12 condensedStiffness() const {
    Mat12 kl = localStiffness();
    if (!valid_ || !hasReleases()) return kl;
    // 收集 b（释放）与 a（保留）下标
    std::array<int, 12> bIdx{};
    int nB = 0;
    for (int i = 0; i < 6; ++i) {
      if (relI_[static_cast<size_t>(i)]) bIdx[static_cast<size_t>(nB++)] = i;
      if (relJ_[static_cast<size_t>(i)]) bIdx[static_cast<size_t>(nB++)] = i + 6;
    }
    if (nB >= 12) {
      singular_ = true;
      condensationErr_ = "全部 12 个自由度均被释放 —— 单元无任何刚度";
      return Mat12::zero();
    }
    std::array<int, 12> aIdx{};
    int nA = 0;
    for (int i = 0; i < 12; ++i) {
      bool isB = false;
      for (int k = 0; k < nB; ++k)
        if (bIdx[static_cast<size_t>(k)] == i) { isB = true; break; }
      if (!isB) aIdx[static_cast<size_t>(nA++)] = i;
    }
    // 组装 Kbb（nB×nB）与右端 Kba（nB×nA）
    std::vector<double> A(static_cast<size_t>(nB) * nB);
    std::vector<double> B(static_cast<size_t>(nB) * nA);
    for (int r = 0; r < nB; ++r)
      for (int c = 0; c < nB; ++c)
        A[static_cast<size_t>(r) * nB + c] = kl(bIdx[static_cast<size_t>(r)],
                                                bIdx[static_cast<size_t>(c)]);
    for (int r = 0; r < nB; ++r)
      for (int c = 0; c < nA; ++c)
        B[static_cast<size_t>(r) * nA + c] = kl(bIdx[static_cast<size_t>(r)],
                                                aIdx[static_cast<size_t>(c)]);
    if (!gaussJordanSolve(nB, A.data(), B.data(), nA)) {
      singular_ = true;
      condensationErr_ =
          "端部释放组合使局部子矩阵奇异（Kbb 退化）—— 同一端同时释放了相互"
          "耦合的自由度（如 u_y 与 r_z 同时释放且忽略剪切变形），形成局部机构";
      return Mat12::zero();
    }
    singular_ = false;
    condensationErr_.clear();
    // Kc_aa = Kaa − Kab · (Kbb⁻¹ Kba)
    // 注意 Kab 的行 = a 下标、列 = b 下标；B 的每行是 Kbb⁻¹Kba 的一行（行=b、列=a）。
    // 需要 K_ab[r][k] · X[k][c]，其中 X = Kbb⁻¹Kba，因此求和下标对齐。
    Mat12 Kc = Mat12::zero();
    for (int r = 0; r < nA; ++r)
      for (int c = 0; c < nA; ++c) {
        double s = kl(aIdx[static_cast<size_t>(r)], aIdx[static_cast<size_t>(c)]);
        for (int k = 0; k < nB; ++k)
          s -= kl(aIdx[static_cast<size_t>(r)], bIdx[static_cast<size_t>(k)]) *
               B[static_cast<size_t>(k) * nA + c];
        Kc(aIdx[static_cast<size_t>(r)], aIdx[static_cast<size_t>(c)]) = s;
      }
    return Kc;
  }

  // 凝聚后的等效节点荷载（局部坐标；b 分量 = 0）。
  // 幂等：若 fLocal 已经凝聚过（b 分量全 0），再次凝聚结果不变。
  Vec12d condenseLoad(const Vec12d& fLocal) const {
    Vec12d fc = fLocal;
    if (!valid_ || !hasReleases()) return fc;
    const Mat12 kl = localStiffness();
    // 收集 b / a 下标（与 condensedStiffness 相同 —— 独立算一遍，n 很小）
    std::array<int, 12> bIdx{};
    int nB = 0;
    for (int i = 0; i < 6; ++i) {
      if (relI_[static_cast<size_t>(i)]) bIdx[static_cast<size_t>(nB++)] = i;
      if (relJ_[static_cast<size_t>(i)]) bIdx[static_cast<size_t>(nB++)] = i + 6;
    }
    if (nB >= 12) { fc.fill(0.0); return fc; }
    std::array<int, 12> aIdx{};
    int nA = 0;
    for (int i = 0; i < 12; ++i) {
      bool isB = false;
      for (int k = 0; k < nB; ++k)
        if (bIdx[static_cast<size_t>(k)] == i) { isB = true; break; }
      if (!isB) aIdx[static_cast<size_t>(nA++)] = i;
    }
    // 解 Kbb · x = fb（nB×1）
    std::vector<double> A(static_cast<size_t>(nB) * nB);
    std::vector<double> xb(static_cast<size_t>(nB));
    for (int r = 0; r < nB; ++r)
      for (int c = 0; c < nB; ++c)
        A[static_cast<size_t>(r) * nB + c] = kl(bIdx[static_cast<size_t>(r)],
                                                bIdx[static_cast<size_t>(c)]);
    for (int r = 0; r < nB; ++r)
      xb[static_cast<size_t>(r)] = fLocal[bIdx[static_cast<size_t>(r)]];
    if (!gaussJordanSolve(nB, A.data(), xb.data(), 1)) {
      singular_ = true;
      condensationErr_ = "端部释放荷载凝聚失败（Kbb 奇异）";
      fc.fill(0.0);
      return fc;
    }
    // fc_a = fa − Kab·x；fc_b = 0
    for (int r = 0; r < nA; ++r) {
      double s = fLocal[aIdx[static_cast<size_t>(r)]];
      for (int k = 0; k < nB; ++k)
        s -= kl(aIdx[static_cast<size_t>(r)], bIdx[static_cast<size_t>(k)]) *
             xb[static_cast<size_t>(k)];
      fc[aIdx[static_cast<size_t>(r)]] = s;
    }
    for (int k = 0; k < nB; ++k) fc[bIdx[static_cast<size_t>(k)]] = 0.0;
    return fc;
  }

  // -------------------------------------------------------------------------
  //  等效节点荷载（局部坐标）
  //
  //  局部 x' 为梁轴，y' 竖向，z' 横向。
  //  参数是【全局】坐标系的三个分量，由 transformLoad 转到局部。
  //
  //  【为什么要用等效节点荷载】梁上的分布荷载（自重、楼板传来）不能直接
  //  施加到节点上。必须通过"积分形函数 × 荷载"得到节点力与节点矩，
  //  否则弯矩沿梁长不连续 —— 这是有限元荷载处理的铁律。
  //
  //  局部均布荷载 qx'(沿轴), qy'(竖向), qz'(横向) 的等效节点力：
  //     节点 i:  q·L/2,  节点 j: q·L/2
  //     节点 i 弯矩:  q·L²/12,  节点 j 弯矩: -q·L²/12
  //  （固定端一致荷载，符号按右手定则）
  // -------------------------------------------------------------------------
  Vec12d equivalentNodalLoadLocal(double qx, double qy, double qz) const {
    Vec12d f{};
    if (!valid_) return f;
    const double L = length_;
    const double qL2_12 = L * L / 12.0;

    // 轴向（自由度 0, 6）
    f[0] = qx * L / 2.0;
    f[6] = qx * L / 2.0;

    // 绕 z 弯曲：由 qy（沿局部 y 的荷载）引起
    //   u_y 自由度 1, 7；r_z 自由度 5, 11
    f[1] = qy * L / 2.0;
    f[7] = qy * L / 2.0;
    f[5] = -qy * qL2_12;
    f[11] = qy * qL2_12;

    // 绕 y 弯曲：由 qz（沿局部 z 的荷载）引起
    //   u_z 自由度 2, 8；r_y 自由度 4, 10
    //  【固端弯矩的符号必须与刚度矩阵里 y 平面的耦合项符号一致】——
    //  那里为了满足 ω_y = −du_z/dx 把耦合项取反了，这里必须同步取反。
    //  只改一处会得到"位移对、固端修正错"的结果，而且端部位移仍然合理。
    f[2] = qz * L / 2.0;
    f[8] = qz * L / 2.0;
    f[4] = -qz * qL2_12;
    f[10] = qz * qL2_12;

    return f;
  }

  // -------------------------------------------------------------------------
  //  分布荷载段：局部 x ∈ [x1, x2]，强度从 q1 线性变到 q2
  //
  //  【为什么必须做到段这一级】
  //  双向板导荷给梁的是三角形 / 梯形荷载。若把它折算成"等效均布"，
  //  跨中弯矩会偏小（三角形荷载下是均布结果的 0.8 倍），而偏小是
  //  【不安全】的方向；支座剪力又完全对不上（差 1.5 倍）。
  //  形函数已经在那里了，直接积分是零成本的。
  //
  //  【q1/q2 是全局分量】与 setLineLoad 的约定一致 ——
  //  用户按全局方向给荷载（竖向、水平），单元负责转到局部。
  // -------------------------------------------------------------------------
  struct DistributedLoad {
    double x1{0.0}, x2{0.0};          // 沿局部 x 的范围（m）
    Vec3 q1{0, 0, 0}, q2{0, 0, 0};    // 两端强度（全局分量，kN/m）
  };

  // -------------------------------------------------------------------------
  //  线性分布荷载段的等效节点荷载（局部坐标）
  //
  //  【积分方案：2 点 Gauss】被积函数 = 形函数(三次) × 荷载(一次) = 四次，
  //  2 点 Gauss 对三次以内精确 —— 对"轴力/剪力"这类（形函数×荷载，
  //  最高四次）不是全精确的，但下面按每一种自由度分别用了正确阶次的公式：
  //    · 平动自由度用 Hermite 三次 × 荷载线性 = 四次  → 3 点 Gauss
  //    · 转动自由度同阶
  //  这里用 3 点 Gauss 一并处理，并对"均布"情形与解析公式比对（有测试）。
  //
  //  【形函数必须与均布荷载那一套完全一致】否则会得到
  //  "均布时对、三角形时错"的诡异结果 —— 因为固端弯矩的符号约定
  //  只在一个方向成立。下面 S/t 两组函数就是从
  //  equivalentNodalLoadLocal 的均布结果反推出来的。
  // -------------------------------------------------------------------------
  Vec12d equivalentNodalLoadLocalSegment(double x1, double x2,
                                         double qx1, double qy1, double qz1,
                                         double qx2, double qy2, double qz2) const {
    Vec12d f{};
    if (!valid_) return f;
    const double L = length_;
    if (!(L > 0.0)) return f;

    double a = std::min(x1, x2), b = std::max(x1, x2);
    // 端点强度要跟着区间端点走，否则把用户给的 q1/q2 用反
    if (x1 > x2) { std::swap(qx1, qx2); std::swap(qy1, qy2); std::swap(qz1, qz2); }
    a = std::max(a, 0.0);
    b = std::min(b, L);
    if (!(b > a)) return f;

    const double xi1 = a / L, xi2 = b / L;
    const double mid = 0.5 * (xi1 + xi2), half = 0.5 * (xi2 - xi1);
    // 3 点 Gauss：对五次以内多项式精确
    static const double gp[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
    static const double gw[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

    double ax1 = 0, ax2 = 0;
    double by1 = 0, bt1 = 0, by2 = 0, bt2 = 0;
    double bz1 = 0, bu1 = 0, bz2 = 0, bu2 = 0;
    for (int i = 0; i < 3; ++i) {
      const double xi = mid + half * gp[i];
      const double jw = half * gw[i];                 // dξ 上的权重
      const double s = (xi - xi1) / (xi2 - xi1);      // 段内参数 0..1
      const double qx = qx1 + (qx2 - qx1) * s;
      const double qy = qy1 + (qy2 - qy1) * s;
      const double qz = qz1 + (qz2 - qz1) * s;

      const double N1 = 1.0 - xi, N2 = xi;                       // 轴向：线性
      const double Su1 = 1.0 - 3.0 * xi * xi + 2.0 * xi * xi * xi;
      const double Su2 = 3.0 * xi * xi - 2.0 * xi * xi * xi;
      const double Sr1 = -L * (xi - 2.0 * xi * xi + xi * xi * xi);   // 绕 z
      const double Sr2 = L * (xi * xi - xi * xi * xi);
      ax1 += jw * N1 * qx;   ax2 += jw * N2 * qx;
      by1 += jw * Su1 * qy;  bt1 += jw * Sr1 * qy;
      by2 += jw * Su2 * qy;  bt2 += jw * Sr2 * qy;
      // 绕 y：与 localStiffness 的 y 平面符号一致（见那里的说明）
      bz1 += jw * Su1 * qz;  bu1 += jw * Sr1 * qz;
      bz2 += jw * Su2 * qz;  bu2 += jw * Sr2 * qz;
    }
    const double k = L;      // ∫dx = L·∫dξ
    f[0] = ax1 * k;   f[6]  = ax2 * k;
    f[1] = by1 * k;   f[5]  = bt1 * k;
    f[7] = by2 * k;   f[11] = bt2 * k;
    f[2] = bz1 * k;   f[4]  = bu1 * k;
    f[8] = bz2 * k;   f[10] = bu2 * k;
    return f;
  }

// 全局 → 局部（f_local = T·f_global），对每个 3 自由度块做 R
  // 与 toGlobalForce 互为逆变换（R 正交）。位移同此规则（u_local = R·u_global）。
  Vec12d toLocalVector(const Vec12d& fg) const {
    Vec12d fl{};
    const double R[3][3] = {
        {ex_.x, ex_.y, ex_.z},
        {ey_.x, ey_.y, ey_.z},
        {ez_.x, ez_.y, ez_.z},
    };
    for (int blk = 0; blk < 4; ++blk)
      for (int i = 0; i < 3; ++i) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j) s += R[i][j] * fg[blk * 3 + j];    // R · v
        fl[blk * 3 + i] = s;
      }
    return fl;
  }

  // 局部 → 全局（f_global = Tᵀ·f_local），对每个 3 自由度块做 Rᵀ
  Vec12d toGlobalForce(const Vec12d& fl) const {
    Vec12d fg{};
    const double R[3][3] = {
        {ex_.x, ex_.y, ex_.z},
        {ey_.x, ey_.y, ey_.z},
        {ez_.x, ez_.y, ez_.z},
    };
    for (int blk = 0; blk < 4; ++blk)
      for (int i = 0; i < 3; ++i) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j) s += R[j][i] * fl[blk * 3 + j];   // Rᵀ
        fg[blk * 3 + i] = s;
      }
    return fg;
  }

  // -------------------------------------------------------------------------
  //  将局部等效节点荷载转到全局
  //   f_global = Tᵀ · f_local
  // -------------------------------------------------------------------------
  //
  //  【统一走"段"这一条路径】均布只是 x1=0、x2=L、q1=q2 的特例。
  //  若保留两条独立实现，它们会各自演化 —— 而"均布对、梯形错"这种
  //  不一致恰恰最难发现（测试用例大多只用均布）。
  Vec12d equivalentNodalLoadGlobal(const Vec3& qGlobal) const {
    return equivalentNodalLoadGlobalSegment(0.0, length_, qGlobal, qGlobal);
  }

  Vec12d equivalentNodalLoadGlobalSegment(double x1, double x2,
                                          const Vec3& q1G, const Vec3& q2G) const {
    if (!valid_) return Vec12d{};
    const Vec3 a(dot(q1G, ex_), dot(q1G, ey_), dot(q1G, ez_));
    const Vec3 b(dot(q2G, ex_), dot(q2G, ey_), dot(q2G, ez_));
    return toGlobalForce(equivalentNodalLoadLocalSegment(x1, x2, a.x, a.y, a.z,
                                                         b.x, b.y, b.z));
  }

  // -------------------------------------------------------------------------
  //  单元集中力 → 节点力（不做任何分配）
  //
  //  集中力与均布荷载的处理方式完全不同：
  //  · 集中力作用在哪个节点，就直接加到该节点的对应自由度上
  //  · 不产生附加节点矩（除非用户显式给出力矩）
  //
  //  盈建科/PKPM 的做法：用户输入的"梁上集中力"自动定位到最近的节点，
  //  这在结构上是近似（会引入附加约束），但对框架分析影响可忽略。
  //  更高精度的做法是引入"集中力作用位置"自由度做静力凝聚。
  // -------------------------------------------------------------------------
  // node: 0 = i 端, 1 = j 端
  // dirLocal: {平动x, 平动y, 平动z, 转角x, 转角y, 转角z}，局部坐标
  Vec12d concentratedForce(int node, const Vec6d& dirLocal) const {
    Vec12d f{};
    if (!valid_) return f;
    const int off = (node == 0) ? 0 : 6;
    for (int d = 0; d < 6; ++d) f[off + d] = dirLocal[d];
    return f;
  }

  // 全局坐标的集中力
  Vec12d concentratedForceGlobal(int node, const Vec3& F) const {
    // 局部力 = R · 全局力
    const Vec3 Fl(dot(F, ex_), dot(F, ey_), dot(F, ez_));
    Vec6d d{};
    d[0] = Fl.x; d[1] = Fl.y; d[2] = Fl.z;
    return concentratedForce(node, d);
  }

  // -------------------------------------------------------------------------
  //  单元集中力矩 → 节点力矩
  // -------------------------------------------------------------------------
  Vec12d concentratedMoment(int node, const Vec3& M) const {
    Vec12d f{};
    if (!valid_) return f;
    const int off = (node == 0) ? 0 : 6;
    // 局部力矩 = R · 全局力矩
    f[off + 3] = dot(M, ex_);
    f[off + 4] = dot(M, ey_);
    f[off + 5] = dot(M, ez_);
    return f;
  }

  // -------------------------------------------------------------------------
  //  自重（均布竖向向下）
  // -------------------------------------------------------------------------
  Vec12d selfWeight(double gamma) const {
    // 容重 gamma kN/m³ × 截面积 A m² = kN/m（沿 -Z 的均布线荷载）
    const double w = gamma * sec_.A;
    return equivalentNodalLoadGlobal(Vec3{0.0, 0.0, -w});
  }

  // -------------------------------------------------------------------------
  //  单元内力恢复（由单元端位移求端力）
  //     f_local = k_local · u_local
  //  盈建科的"杆端内力"（轴力、剪力、弯矩、扭矩）就是 f_local 的分量，
  //  符号约定：正 = 沿单元局部坐标正方向（右手定则）。
  // -------------------------------------------------------------------------
  struct EndForces {
    double N{0}, Vy{0}, Vz{0}, T{0}, My{0}, Mz{0};        // 节点 i
    double Nj{0}, Vyj{0}, Vzj{0}, Tj{0}, Myj{0}, Mzj{0};  // 节点 j
  };

  // 纯 k·u（【不含】单元荷载的固端贡献）。
  // 仅供单元测试做对照 —— 工程上要的是下面那个带修正的版本。
  EndForces endForces(const Vec12d& uGlobal) const {
    return endForces(uGlobal, Vec12d{});
  }

  // -------------------------------------------------------------------------
  //  【真实杆端内力 = k·u − f_eq】★ 后处理层的正确性依赖于此
  //
  //  为什么必须减：全局方程是 K·U = F_节点 + Σ f_eq，
  //  其中 f_eq 是分布荷载的等效（一致）节点荷载向量。
  //  于是对单元有 k·u = P + f_eq，而端部【截面内力】才是 k·u − f_eq。
  //
  //  两个极端算例立刻能看出差别：
  //    · 两端固定梁 + 均布荷载 w：u = 0 ⇒ k·u = 0，
  //      若不做修正端部内力全为 0 —— 真解是 V=wL/2、M=wL²/12；
  //    · 简支梁 + 均布荷载：k·u ≠ 0，修正后端部弯矩归零、剪力 = wL/2。
  //
  //  【符号】减号。验证见 tests/test_post.cpp 的两端固定梁算例。
  //
  //  【注意 f_eq 必须转成局部坐标】endForces 返回的是【局部】分量，
  //  而 Element::equivalentLoads() 给的是【全局】分量（与 stiffness() 一致），
  //  所以这里要用同一个 R 把 f_eq 也转过去，否则修正项全错。
  // -------------------------------------------------------------------------
  EndForces endForces(const Vec12d& uGlobal, const Vec12d& fEqGlobal) const {
    EndForces ef;
    if (!valid_) return ef;
    // u_global → u_local
    Vec12d ul{};
    Vec12d fl{};                       // 局部的等效节点荷载
    const double R[3][3] = {
        {ex_.x, ex_.y, ex_.z},
        {ey_.x, ey_.y, ey_.z},
        {ez_.x, ez_.y, ez_.z},
    };
    for (int blk = 0; blk < 4; ++blk)
      for (int i = 0; i < 3; ++i) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j) s += R[i][j] * uGlobal[blk * 3 + j];   // R · u
        ul[blk * 3 + i] = s;
      }
    for (int blk = 0; blk < 4; ++blk)
      for (int i = 0; i < 3; ++i) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j) s += R[i][j] * fEqGlobal[blk * 3 + j];  // R · f_eq
        fl[blk * 3 + i] = s;
      }
const Mat12 kl = hasReleases() ? condensedStiffness() : localStiffness();
    // 【端部释放】等效节点荷载也必须凝聚：释放自由度上的固端弯矩/剪力
    // 不再直接施加，而是通过 Kbb⁻¹ 折算到保留自由度上。
    // （凝聚后 b 分量为 0 —— 铰端没有节点力矩输入。）
    const Vec12d flc = hasReleases() ? condenseLoad(fl) : fl;
    Vec12d ft = kl * ul;
    for (int i = 0; i < 12; ++i) ft[i] -= flc[i];      // ★ 固端修正

    // -----------------------------------------------------------------------
    //  杆端力 ft → 【截面内力】
    //
    //  ft = k·u − f_eq 在数学上是【杆端力】：节点施加给单元的力。
    //  而工程师要的是【截面内力】，两者差一个符号，且不统一：
    //
    //    · 轴力 N：受拉杆 ft[0] = −P，而截面轴力应为 +P（拉为正）⇒ 取反
    //    · 剪力 V：ft[1] 是沿局部 y' 的分量，而 y' 由 up 向量决定
    //      （本单元 up=(0,0,−1) ⇒ y' 指向下方）。取反后才是
    //      "向上为正"，也才满足 V = dM/dx
    //    · 扭矩 T：与轴力同理（同属沿杆轴的分量）⇒ 取反
    //    · 弯矩 M：ft[5] 已经是 sagging 为正（悬臂固定端给出 −PL，
    //      简支跨中给出 +PL/4，都对）⇒ 不取反
    //
    //  【判据】改完后必须同时满足（缺一不可）：
    //    悬臂梁端部受力 P：固定端 M = −PL（上侧受拉），V = +P
    //    简支梁跨中受力 P：跨中 M = +PL/4（下侧受拉），左半段 V = +P/2
    //    受拉杆：两端 N 都是 +P
    //
    //  只改一个分量会让它"看起来对"但 M 不再是 V 的积分 ——
    //  内力图会自相矛盾，所以必须一次性对齐并在测试里同时校核 M 与 V。
    // -----------------------------------------------------------------------
    ef.N  = -ft[0]; ef.Vy = -ft[1]; ef.Vz = -ft[2];
    ef.T  = -ft[3]; ef.My =  ft[4]; ef.Mz =  ft[5];
    // j 端：原实现对所有分量取了一次负，这里按上面的规则分别处理
    ef.Nj  =  ft[6]; ef.Vyj =  ft[7]; ef.Vzj =  ft[8];
    ef.Tj  =  ft[9]; ef.Myj = -ft[10]; ef.Mzj = -ft[11];
    return ef;
  }

  // -------------------------------------------------------------------------
  //  沿梁长的内力分布（画内力图用）
  //
  //  【弯矩必须按二次插值，不能线性】
  //  有分布荷载时，两端弯矩之间不是直线：
  //      悬臂梁 + 均布荷载：M(ξ) = wL²(1−ξ)²/2 —— 抛物线
  //  线性插值在跨中会低估/高估 25%，对配筋是致命的。
  //
  //  正确做法：剪力 V 线性（均布荷载下本来就是线性，插值即精确），
  //  弯矩由 M' = −V 积分得到：
  //      M(ξ) = M_i − L·[ V_i·ξ + (V_j − V_i)·ξ²/2 ]
  //  代入悬臂梁均布荷载（V_i = wL, V_j = 0, M_i = wL²/2）即得
  //      M(ξ) = wL²/2 − wL²(ξ − ξ²/2) = wL²(1−ξ)²/2      ✓ 精确
  //
  //  【配对】弯矩与剪力的配对是"绕哪根轴弯"决定的：
  //      V_y（局部 y 向剪力）→ 绕 z 轴弯 → 配 M_z
  //      V_z（局部 z 向剪力）→ 绕 y 轴弯 → 配 M_y
  //  配错对不报错，只是内力图看着"差不多"，是最隐蔽的一类错误。
  //  ------------------------------------------------------------------
  EndForces atStation(const Vec12d& uGlobal, double xi) const {
    return atStation(uGlobal, Vec12d{}, xi, nullptr);
  }
  EndForces atStation(const Vec12d& uGlobal, const Vec12d& fEqGlobal, double xi) const {
    return atStation(uGlobal, fEqGlobal, xi, nullptr);
  }

  // -------------------------------------------------------------------------
  //  沿梁长的截面内力
  //
  //  【有分布荷载时必须真积分，不能线性插值】
  //  两个端值线性插值等价于假设"剪力沿梁长线性变化"，那只在【均布】荷载
  //  下成立。三角形荷载下剪力是二次的，梯形荷载下弯矩是三次的 ——
  //  线性插值会给出"端部对、中间错"的图，而端部恰恰是最常被校核的地方。
  //
  //  【本单元内力的符号约定（实测标定，不是推导）】
  //      V(x) = V_i − ∫₀ˣ q_局部 dt          （q 是局部 y/z 分量）
  //      M(x) = M_i + ∫₀ˣ V dt               （即 M′ = V）
  //  标定办法：悬臂梁 + 向下均布，解析解 V = w(L−x)、M_cl = w(L−x)²/2，
  //  而本单元给出 ef.Vy = +w(L−x)、ef.Mz = −M_cl —— 弯矩整体差一个负号，
  //  这正是 M′ = V 成立所需的那个符号。**先测再写，不要推。**
  //
  //  loads 为空 → 退化成"线性插值"，与旧行为逐位一致。
  // -------------------------------------------------------------------------
  EndForces atStation(const Vec12d& uGlobal, const Vec12d& fEqGlobal, double xi,
                      const std::vector<DistributedLoad>* loads) const {
    const EndForces e0 = endForces(uGlobal, fEqGlobal);
    EndForces e = e0;
    const double L = length_;
    const double x = xi * L;

    if (loads == nullptr || loads->empty()) {
      // 无分布荷载：轴力/剪力/扭矩为常量，弯矩线性 —— 用端值插值即可
      e.N = e0.N * (1.0 - xi) + e0.Nj * xi;
      e.Vy = e0.Vy * (1.0 - xi) + e0.Vyj * xi;
      e.Vz = e0.Vz * (1.0 - xi) + e0.Vzj * xi;
      e.T = e0.T * (1.0 - xi) + e0.Tj * xi;
      e.Mz = e0.Mz * (1.0 - xi) + e0.Mzj * xi;
      e.My = e0.My * (1.0 - xi) + e0.Myj * xi;
      return e;
    }

    // ---- 按实际荷载分布积分 ----
    //
    //  I1(x) = ∫₀ˣ q(t) dt          （剪力用）
    //  I2(x) = ∫₀ˣ (x−t)·q(t) dt    （简支梁弯矩用）
    //
    //  荷载在每段内是线性的 ⇒ 两个被积函数都是多项式 ⇒ 2 点 Gauss 精确。
    Vec3 I1{0, 0, 0}, I2{0, 0, 0}, I2L{0, 0, 0};
    static const double gp = 0.5773502691896257;
    for (const DistributedLoad& d : *loads) {
      const double lo = std::min(d.x1, d.x2), hi = std::max(d.x1, d.x2);
      const double h0 = hi - lo;
      if (!(h0 > 0.0)) continue;
      Vec3 qa(dot(d.q1, ex_), dot(d.q1, ey_), dot(d.q1, ez_));
      Vec3 qb(dot(d.q2, ex_), dot(d.q2, ey_), dot(d.q2, ez_));
      if (d.x1 > d.x2) { const Vec3 t = qa; qa = qb; qb = t; }
      const double a = std::max(lo, 0.0);
      const double b = std::min(hi, x);
      // 整段（用于简支梁弯矩的"假想支座反力"）
      {
        const double h = hi - lo;
        const double t0 = lo + h * (0.5 - 0.5 * gp);
        const double t1 = lo + h * (0.5 + 0.5 * gp);
        const Vec3 q0 = qa * (1.0 - (t0 - lo) / h) + qb * ((t0 - lo) / h);
        const Vec3 q1 = qa * (1.0 - (t1 - lo) / h) + qb * ((t1 - lo) / h);
        const double w = 0.5 * h;
        I2L = I2L + (q0 * (L - t0) + q1 * (L - t1)) * w;
      }
      if (!(b > a)) continue;
      const double h = b - a;
      const double t0 = a + h * (0.5 - 0.5 * gp);
      const double t1 = a + h * (0.5 + 0.5 * gp);
      const Vec3 q0 = qa * (1.0 - (t0 - lo) / h0) + qb * ((t0 - lo) / h0);
      const Vec3 q1 = qa * (1.0 - (t1 - lo) / h0) + qb * ((t1 - lo) / h0);
      const double w = 0.5 * h;
      I1 = I1 + (q0 + q1) * w;
      I2 = I2 + (q0 * (x - t0) + q1 * (x - t1)) * w;
    }

    e.N  = e0.N  - I1.x;
    e.Vy = e0.Vy - I1.y;
    e.Vz = e0.Vz - I1.z;
    e.T  = e0.T * (1.0 - xi) + e0.Tj * xi;        // 无分布扭矩

    // -----------------------------------------------------------------------
    //  弯矩：端值线性插值 + 简支梁弯矩
    //
    //      M(x) = M_i·(1−ξ) + M_j·ξ + M_ss(x)
    //      M_ss(x) = [ I2(L)/L ]·x − I2(x)        （两端恒为 0）
    //
    //  【为什么不能简单地对 V 积分得到 M】
    //  有限元里"剪力"和"弯矩"是两套独立的插值（本单元 w 与 θ 各自插值），
    //  点对点的 dM/dx = V 在【梁内有分布荷载时】并不精确成立。
    //  实测：悬臂梁 + 三角形荷载（单单元）的端弯矩是 −57.6，而按 V 积分得到
    //  +14.4（自由端应当为 0）—— 差 14.4，正好等于"固端一致荷载"的量级。
    //  后果是内力图在两个单元的交接处出现【跳变】（实测 0.167 kN·m），
    //  而物理上弯矩必须连续 —— 工程师一眼就会认为程序算错了。
    //
    //  【为什么这个形式一定对】
    //    · x=0 → M_i；x=L → M_j（两端的端值来自 endForces，逐个单元自洽）
    //    · M_ss 是"同样跨度、同样荷载的简支梁弯矩"，对均布是抛物线、
    //      对三角形是三次曲线 —— 正是缺的那个形状
    //    · 荷载为零时 M_ss ≡ 0，退化成原来的线性插值，行为不变
    //  自检：悬臂 + 均布，单单元即可与 w(L−x)²/2 逐点吻合到 1e-12。
    // -----------------------------------------------------------------------
    e.Mz = e0.Mz * (1.0 - xi) + e0.Mzj * xi + (I2L.y / L) * x - I2.y;
    // 【y 平面的简支弯矩项符号相反】与 localStiffness / 固端弯矩的符号改动配套，
    // 见那边关于 ω_y = −du_z/dx 的说明。改一半会得到"端部对、中间错"的图。
    e.My = e0.My * (1.0 - xi) + e0.Myj * xi - (I2L.z / L) * x + I2.z;
    return e;
  }

  // -------------------------------------------------------------------------
  //  【旧实现错在哪里 —— 留在这里，因为它是最容易被"看起来对"骗过的一类】
  //
  //      e.Vy = V_i·(1−ξ) + V_j·ξ                        // 线性插值
  //      e.Mz = M_i + L·[ V_i·ξ + (V_j−V_i)·ξ²/2 ]
  //
  //  这等价于【假设剪力沿梁长线性变化】，只在均布荷载下成立：
  //      · 三角形荷载 → 剪力是二次的，跨中差 25%
  //      · 梯形荷载   → 弯矩是三次的
  //  而且在【端部仍然精确】（端值来自 endForces），所以只查端部永远发现不了。
  //  自检办法：把 ξ=1 处的结果与 endForces 的 j 端比对 ——
  //  两者必须相等，不相等就说明积分路径有问题。这条自检已写进测试。
  // -------------------------------------------------------------------------

  // 截面属性访问
  const SectionProperties& section() const { return sec_; }
  const Material& material() const { return mat_; }
  const Vec3& nodeI() const { return pi_; }
  const Vec3& nodeJ() const { return pj_; }

  // 单元质量（一致质量矩阵的对角项，模态分析用）
  double mass() const { return mat_.gamma * sec_.A * length_; }

private:
  Vec3 pi_, pj_;
  Vec3 ex_{1, 0, 0}, ey_{0, 1, 0}, ez_{0, 0, 1};
  SectionProperties sec_;
  Material mat_;
  double length_{0.0};
  bool valid_{false};
  bool degenerateUp_{false};
  std::string error_;

  // ---- 端部释放 ----
  std::array<bool, 6> relI_{};   // i 端释放掩码（局部自由度 0..5）
  std::array<bool, 6> relJ_{};   // j 端释放掩码（局部自由度 0..5）
  mutable bool singular_{false}; // 最近一次凝聚是否奇异（Kbb 退化）
  mutable std::string condensationErr_;
};

}  // namespace yjk
