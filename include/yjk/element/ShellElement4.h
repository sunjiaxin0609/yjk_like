// =============================================================================
//  yjk/element/ShellElement4.h  ——  四节点壳单元（膜 + 板，MITC4）
//
// =============================================================================
//  【自由度排列】—— 与梁单元一致，这是全软件的接口约定
//
//      节点 i (0..5):   ux    uy    uz    rx    ry    rz
//      节点 j (6..11):  ux    uy    uz    rx    ry    rz
//
//      壳单元把 6 个自由度分成两组独立工作：
//        膜（平面内）：ux, uy          —— 2 个/节点
//        板（平面外）：uz, rx, ry      —— 3 个/节点
//        rz：绕板面法向的转角
//            【注意】它既不属于膜（膜里扭转被抑制），
//            也不属于板（板里绕法向转动无弯曲）。
//            盈建科的处理是把它作为"钻转刚度"的自由度单独赋予一个小刚度。
//            本实现采用标准做法：给一个可调的小抗扭刚度，
//            用户可设为 0（则该自由度被约束）。
//
//  【坐标约定】
//      壳单元本身是平面的，局部坐标 x-y 在板面内，z 为板面法向。
//      局部系原点取在单元形心。
//      节点坐标逆时针排列（保证雅可比为正）。
//
//  【刚度组成】
//      K_shell = Tᵀ · (K_membrane + K_plate) · T
//
//      K_membrane : 平面应力/平面应变 —— 楼板的面内受力
//      K_plate    : Mindlin-Reissner 板弯曲（含剪切）—— 楼板的面外变形
//
//  【平面应力 vs 平面应变】
//      楼板、剪力墙通常按【平面应力】（薄壁构件，σ_z ≈ 0）
//      但当单元位于楼板中部、上下均有楼板时，可按【平面应变】（约束更强）
//      盈建科默认平面应力，并在属性里提供切换。
//      公式差异只在 D 矩阵：
//        平面应力：D = E/(1-ν²) · [1  ν  0; ν  1  0; 0 0 (1-ν)/2]
//        平面应变：D = E/((1+ν)(1-2ν)) · [1-ν ν 0; ν 1-ν 0; 0 0 (1-2ν)/2]
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <string>

#include "yjk/math/Types.h"
#include "yjk/model/Section.h"

namespace yjk {

// ---------------------------------------------------------------------------
//  壳单元的截面参数
// ---------------------------------------------------------------------------
struct ShellProperties {
  double thickness{0.2};        // 厚度 m
  double E{3.0e7};              // 弹性模量 kPa
  double nu{0.2};               // 泊松比
  bool planeStrain{false};      // false=平面应力（薄板）, true=平面应变
  double drillingStiffness{0.0};// 绕法向的抗扭刚度系数（0 = 约束该自由度）
                                     // >0 时 K_rz = coeff · E·t·a²·a²/a
  double membraneThinningFactor{1.0};  // 膜单元的泊松比折减（考虑混凝土开裂后软化）
  double density{2500.0};        // 密度 kg/m³（模态分析用；混凝土 2500，钢 7850）
                                     // 与 Material.gamma（容重 kN/m³）的区别：
                                     // gamma = density × 9.81 / 1000，这里不混用，
                                     // 避免"单位是 kN 还是 N"的隐式换算出错。

  // 平面内弹性矩阵 D（3×3，应力-应变：εx, εy, γxy）
  Mat3 membraneD() const {
    Mat3 D = Mat3::zero();
    if (planeStrain) {
      const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
      D(0, 0) = c * (1.0 - nu);
      D(1, 1) = c * (1.0 - nu);
      D(0, 1) = D(1, 0) = c * nu;
      D(2, 2) = c * (1.0 - 2.0 * nu) / 2.0;
    } else {
      const double c = E / (1.0 - nu * nu);
      D(0, 0) = c;
      D(1, 1) = c;
      D(0, 1) = D(1, 0) = c * nu * membraneThinningFactor;
      D(2, 2) = c * (1.0 - nu * membraneThinningFactor) / 2.0;
    }
    return D;
  }

  // 剪切模量（板弯曲的剪切部分用）
  double G() const { return E / (2.0 * (1.0 + nu)); }

  // 弯曲刚度 D_b = E·t³/[12(1-ν²)]  （每单位宽度）
  double bendingRigidity() const {
    return E * thickness * thickness * thickness / (12.0 * (1.0 - nu * nu));
  }
  // 膜弯曲刚度 D_m = E·t/(1-ν²)
  double membraneRigidity() const { return E * thickness / (1.0 - nu * nu); }

  // 剪切刚度（Mindlin 板）：κ·G·t，其中 κ 是剪切修正系数
  // 取 5/6 是标准值（与 Timoshenko 梁一致）
  double shearRigidity() const { return 5.0 / 6.0 * G() * thickness; }

  bool isValid(std::string* why = nullptr) const {
    if (!(thickness > 0)) { if (why) *why = "壳单元厚度为零或为负"; return false; }
    if (!(E > 0)) { if (why) *why = "弹性模量为零或为负"; return false; }
    if (!(nu > -1.0 && nu < 0.5)) { if (why) *why = "泊松比超出物理范围 (-1, 0.5)"; return false; }
    return true;
  }
};

// ---------------------------------------------------------------------------
//  四节点壳单元
// ---------------------------------------------------------------------------
class ShellElement4 {
 public:
  static constexpr int kDof = 24;

  // p[0..3] : 节点坐标，逆时针
  // 默认构造：单元处于"未初始化"状态，isValid() 返回 false。
  // 存在的意义是让上层容器（std::vector、可选成员）能持有未构建的实例，
  // 随后用赋值运算符替换。ShellElement4 是纯值语义，无堆资源，安全。
  ShellElement4() = default;

  ShellElement4(const Vec3 p[4], const ShellProperties& props) {
    init(p, props);
  }

  void init(const Vec3 p[4], const ShellProperties& props) {
    props_ = props;
    for (int i = 0; i < 4; ++i) p_[i] = p[i];

    if (!props_.isValid(&error_)) { valid_ = false; return; }

    // ---- 建立局部坐标 ----
    // 原点：形心
    for (int i = 0; i < 4; ++i) centroid_ = centroid_ + p_[i];
    centroid_ = centroid_ * 0.25;

    // x 轴：节点 0 → 节点 1
    Vec3 ex = p_[1] - p_[0];
    const double l01 = norm(ex);
    if (!(l01 > kEps)) { valid_ = false; error_ = "壳单元边 01 长度为零"; return; }
    ex = ex * (1.0 / l01);

    // 【这里曾经漏掉 ex_ = ex; —— 是一个潜伏了很久的 bug】
    //
    //  只写了局部变量 ex，没有写回成员 ex_。于是 ex_ 一直是默认值 (1,0,0)。
    //  后果：
    //    · 局部 y、z 轴是正确的（ey 由 ez × ex 用局部变量算，ez 由法向算），
    //      只有 ex_ 错 —— 局部坐标系呈"半个正确"的畸形状态
    //    · local_[i].x = d·ex_ 退化成全局 x 分量，
    //      而 local_[i].y 又是正确的局部 y ⇒ 两个坐标来自不同的基
    //    · 轴对齐的矩形单元上 ex 恰好就是 (1,0,0)，一切正常 ⇒ 测试全过
    //    · 斜交轴网、旋转后的板、任意方向的墙 —— 刚度矩阵与内力全错，
    //      且不报错、量级也正常（刚体校核在旋转后才会失败）
    //  是"测试用例没离开过舒适区"的典型：所有既有用例都建在全局坐标轴上。
    ex_ = ex;

    // z 轴：单元法向（由右手定则，保证局部系与全局右手系一致）
    Vec3 nrm = cross(p_[1] - p_[0], p_[2] - p_[0]);
    const double area2 = norm(nrm);
    if (!(area2 > kEps)) { valid_ = false; error_ = "壳单元前三点共线（面积为零）"; return; }
    ez_ = nrm * (1.0 / area2);

    // y 轴：ez × ex —— 用【局部变量 ex】而不是 ex_，
    // 这样即使将来 ex_ 的赋值又被误删，ey 仍然是正交的（错误不会雪崩）
    Vec3 ey = cross(ez_, ex);
    if (!(norm(ey) > 1e-10)) { valid_ = false; error_ = "壳单元局部坐标构造失败"; return; }
    ey_ = normalized(ey);

    // ---- 局部坐标（形心为原点）----
    for (int i = 0; i < 4; ++i) {
      const Vec3 d = p_[i] - centroid_;
      local_[i] = Vec3{dot(d, ex_), dot(d, ey_), dot(d, ez_)};
    }

    // ---- 检查逆时针（雅可比应为正）----
    // 用叉积判断：节点顺序若为顺时针，需要提示用户
    const Vec3 e1 = local_[1] - local_[0];
    const Vec3 e2 = local_[3] - local_[0];
    const double cr = e1.x * e2.y - e1.y * e2.x;
    if (cr < 0) {
      // 逆时针判定失败。这不是致命错误（局部坐标已按ez=ex×ey定好），
      // 但会导致雅可比符号不一致，需要警告。
      clockwise_ = true;
    }

    // ---- 面积 ----
    //
    // 【关键】必须用【全局坐标】的鞋带公式，不能用局部坐标。
    // 局部坐标是单元自身的平面（z' = 0），用局部坐标算出的面积是
    // "投影面积"，对倾斜单元会偏小 —— 例如把 2×3 的板绕 x 轴倾斜，
    // 投影面积变成 2×3·cos(θ)，而真实面积不变。
    // 后果：倾斜楼板的膜应力与弯曲应力都按错误面积缩放，偏小可达 30%，
    // 且程序不报错。盈建科/PKPM 都用全局坐标的面积。
    //
    // 抗畸变优于三角剖分（Newell 法对凹四边形更稳健）。
    // Newell 公式：S = 1/2 · Σ (V_i × V_{i+1})
    // 【符号约定】V = P − 原点。对平面多边形，|S| 就是真实面积。
    // 常见错误：写成 1/2 Σ V_i × V_{i−1}（顺序反了）会让四个三角形
    // 首尾相消、面积偏小 —— 倾斜单元上尤其明显（水平单元看不出来，
    // 因为此时退化到熟悉的鞋带公式）。
    Vec3 sum{0, 0, 0};
    for (int i = 0; i < 4; ++i) {
      const Vec3 a = p_[i] - centroid_;
      const Vec3 b = p_[(i + 1) % 4] - centroid_;
      sum = sum + cross(a, b);
    }
    area_ = 0.5 * norm(sum);
    if (!(area_ > kEps)) { valid_ = false; error_ = "壳单元面积为零"; return; }

    // ---- 一阶导数（常应变四边形 CST 型，用于膜部分）----
    // 对线性形函数，导数为常数，只依赖 4 个点的坐标
    computeShapeGradients();

    valid_ = true;
  }

  bool isValid() const { return valid_; }
  const std::string& error() const { return error_; }
  double area() const { return area_; }
  double thickness() const { return props_.thickness; }
  const ShellProperties& properties() const { return props_; }
  const Vec3& normal() const { return ez_; }

  // -------------------------------------------------------------------------
  //  4 节点双线性形函数及其导数
  //
  //  自然坐标 (ξ, η) ∈ [-1,1]²：
  //    N1 = ¼(1-ξ)(1-η),  N2 = ¼(1+ξ)(1-η)
  //    N3 = ¼(1+ξ)(1+η),  N4 = ¼(1-ξ)(1+η)
  //
  //  对线性（常应变）单元，形函数导数为常数：
  //    ∂N1/∂x = -y/(4A),  ∂N1/∂y = x/(4A)     —— A 为单元面积
  //  这是 CST 的经典结果：导数只与对边坐标有关，与本节点无关。
  // -------------------------------------------------------------------------
  static void shapeFunctions(double xi, double eta, double N[4]) {
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
  }

  // 形函数对 x, y 的导数（常数）
  const double* dNdx() const { return dNdx_; }
  const double* dNdy() const { return dNdy_; }

  // -------------------------------------------------------------------------
  //  膜刚度（平面内）8×8
  //
  //  B 矩阵（3×8）：
  //    ε = B·u,  ε = [εx, εy, γxy]ᵀ
  //    B_i = | ∂Ni/∂x   0        ∂Ni/∂y |
  //          | 0         ∂Ni/∂x  ∂Ni/∂y |
  //
  //  K_m = t · ∫ Bᵀ·D·B dA   （常应变 ⇒ 积分 = A·(单位)）
  // -------------------------------------------------------------------------
  //  膜刚度（平面内）8×8
  //
  //    应变：ε = B·u,  ε = [εx, εy, γxy]ᵀ
  //    B_i = | ∂Ni/∂x   0       ∂Ni/∂y |      （节点 i 的 3×2 块）
  //          | 0         ∂Ni/∂x ∂Ni/∂y |
  //          | ∂Ni/∂y   ∂Ni/∂x 0       |
  //
  //    K_m = t·A·Bᵀ·D·B     （常应变 ⇒ ∫dA = A）
  //
  //  ---------------------------------------------------------------------
  //  【易错点 1】B 的第三行（γxy）是 {∂Ni/∂y, ∂Ni/∂x}。
  //
  //  因为 γxy = ∂u_x/∂y + ∂u_y/∂x，两项用的是【不同自由度】配【不同方向】的导数：
  //      ∂u_x/∂y → u_x 配 ∂Ni/∂y
  //      ∂u_y/∂x → u_y 配 ∂Ni/∂x
  //
  //  误写成 {∂Ni/∂y, ∂Ni/∂y}（两个都配 ∂/∂y）是最常见的错误，
  //  后果是刚体旋转校核通不过（残余与转角 θ 成正比，约 0.4%），
  //  但平移校核仍通过 —— 所以只做"平移校核"会漏掉这个 bug。
  //
  //  ---------------------------------------------------------------------
  //  【易错点 2】必须是【双重节点循环】 i × j。
  //  只算对角块（i = j）会导致刚体平移产生内力。
  Mat8 membraneStiffness() const {
    Mat8 K = Mat8::zero();
    if (!valid_) return K;

    const Mat3 D = props_.membraneD();
    const double t = props_.thickness;
    const double fac = t * area_;

    // 预存 B_i 的 3×2 块
    double Bi[4][3][2];
    for (int i = 0; i < 4; ++i) {
      const double dx = dNdx_[i], dy = dNdy_[i];
      Bi[i][0][0] = dx; Bi[i][0][1] = 0.0;   // εx = ∂u_x/∂x
      Bi[i][1][0] = 0.0; Bi[i][1][1] = dy;   // εy = ∂u_y/∂y  ← 注意是 dy！
      Bi[i][2][0] = dy; Bi[i][2][1] = dx;   // γxy = ∂u_x/∂y + ∂u_y/∂x
    }

    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        for (int a = 0; a < 2; ++a) {
          for (int b = 0; b < 2; ++b) {
            double v = 0.0;
            for (int r = 0; r < 3; ++r)
              for (int s = 0; s < 3; ++s) {
                const double d = D(r, s);
                if (nearlyZero(d)) continue;
                v += Bi[i][r][a] * d * Bi[j][s][b];
              }
            K(2 * i + a, 2 * j + b) += fac * v;
          }
        }
      }
    }
    return K;
  }

  // -------------------------------------------------------------------------
  //  板刚度（平面外）12×12 —— MITC4
  //
  //  局部板自由度：节点 i → (w, θx, θy) = (3i, 3i+1, 3i+2)
  //
  //  ============================  符号约定（必须先定死）  ============================
  //
  //  取法线在两个竖直平面内的转角作为独立场：
  //      βx ≡  θy      （法线在 x–z 面内转动，正号指向 +x）
  //      βy ≡ −θx      （法线在 y–z 面内转动，正号指向 +y）
  //
  //  于是
  //      曲率：κx = ∂βx/∂x =  ∂θy/∂x
  //            κy = ∂βy/∂y = −∂θx/∂y
  //            κxy= ∂βx/∂y + ∂βy/∂x = ∂θy/∂y − ∂θx/∂x
  //      横向剪切：γxz = ∂w/∂x + βx = ∂w/∂x + θy
  //                γyz = ∂w/∂y + βy = ∂w/∂y − θx
  //
  //  【为什么必须先定死符号】κy 的定义里那个负号，与 Mx/My 的公式是绑在一起的。
  //  单独看刚度矩阵（只用到 κᵀ D κ 这个二次型）根本看不出符号定义，
  //  但内力恢复（Mx = D₀(κx + νκy)）会因此整体错位 ——
  //  表现为"位移算对了、弯矩图是错的"，且形似正确极难发现。
  //  本单元用一条独立判据把它锁死（见下节"验证判据"）。
  //
  //  ============================  剪切锁定与 MITC4  ============================
  //
  //  【问题】纯弯曲时中性层的横向剪切应变必须为零。
  //  4 节点双线性单元无法精确表示 w = c·x² 这样的抛物挠曲，
  //  位移插值给出的 γxz 是一个与坐标成正比的"寄生剪切"。
  //  而厚度越薄，剪切刚度 κGt 相对弯曲刚度 Et³/12 越占优
  //  （比值 10G/(E t²)，t=0.1 m、E=3e7 时约 42 倍），
  //  于是单元被这个虚假的剪切"锁住"，挠度比理论值小 1~2 个数量级。
  //
  //  【Dvorkin & Bathe 1980 的解法 —— 假定自然应变（ANS）】
  //  不去直接插值 γ 的两个分量，而是：
  //    ① 先在每条边的中点取【共变剪切应变】
  //         γ_rz = ∂w/∂r + (∂x/∂r)·θy − (∂y/∂r)·θx
  //         γ_sz = ∂w/∂s + (∂x/∂s)·θy − (∂y/∂s)·θx
  //       （共变分量 = 笛卡尔分量在自然坐标基上的投影：γ_rz = g_r·(γxz, γyz)）
  //    ② 用【线性插值】把 4 个采样值铺满单元：
  //         γ_rz(r,s) = ½(1−s)·γ_rz(0,−1) + ½(1+s)·γ_rz(0,+1)
  //         γ_sz(r,s) = ½(1−r)·γ_sz(−1,0) + ½(1+r)·γ_sz(+1,0)
  //    ③ 用高斯点处的 Jacobian 逆把共变分量换回笛卡尔分量
  //
  //  【为什么绑定点取边中点就对了 —— 这可以手算验证】
  //  取矩形单元 x = (a/2)r，纯弯曲模式 w = c·x²、θy = −2cx、θx = 0
  //  （真实剪切为零）。位移插值给出分片双线性的 w，于是
  //      γ_rz = (a²c/2)·r        ← 寄生剪切，正是锁定的来源
  //  在边中点 (0, ±1) 处：γ_rz(0,±1) = (a²c/2)·0 = 0，
  //  线性插值后处处为零 —— 寄生剪切被完全消除。
  //  而常剪切模式 w = γ₀·x 下：γ_rz(0,±1) = γ₀a/2，插值后仍为 γ₀a/2，
  //  换算回笛卡尔分量恰好得到 γ₀ —— 常剪切也被精确表示。
  //  ⇒ 这两个模式都能精确表示，正是"不锁定"的充要条件。
  //
  //  【已实测失败的做法 —— 记录在此，避免重复踩】
  //    ✗ B 矩阵的 γ 行直接用 ∂N_i/∂x 与 N_i（体积插值，即普通低阶单元）
  //      → 剪切项大 100 倍，挠度小约 1500 倍；t 降到 0.02 时挠度甚至变负
  //    ✗ 只把 θ 的插值降阶、保留 ∂w/∂x 为常数 → 同样无改善
  //    ✗ 绑定点取"相邻边中点的中点" → 只对平行四边形成立，畸变单元不成立
  //    ✗ 绑定点改成 (r=±1, s=0) 采样 γ_rz（即按 r 方向插值）
  //      → 纯弯曲下退化为 γ_rz = (a²c/2)·r，与位移插值一模一样，等于没做
  //      （采样点必须落在"另一方向"的两端，插值方向不能弄反）
  //
  //  ============================  验证判据  ============================
  //  ① 四边简支方板均布荷载，w_max = 0.00406·q·a⁴/D；
  //     t 从 0.2 降到 0.01 时 w/w_理论 必须稳定在 0.9~1.05（不随 t 崩塌）
  //  ② 厚板（t/a = 1/10）因剪切变形【比薄板理论更软】—— 这是物理，不是误差
  //  ③ 8×8 网格跨中挠度与薄板理论偏差 < 15%
  //  ④ 悬臂板端部挠度对标 w = q·a⁴/(8D)（单位宽板带，可比 Timoshenko 解）
  //  见 tests/test_shell.cpp 的 testPlateBending。
  //  =========================================================================
  Mat12 plateStiffness() const {
    Mat12 K = Mat12::zero();
    if (!valid_) return K;

    const double t = props_.thickness;
    const double nu = props_.nu;
    // 弯曲刚度系数 D₀ = E t³/[12(1−ν²)]，弯曲弹性矩阵 = D₀·[[1,ν,0],[ν,1,0],[0,0,(1−ν)/2]]
    const double D0 = props_.E * t * t * t / (12.0 * (1.0 - nu * nu));
    const double Ds = props_.shearRigidity();       // 5/6·G·t
    if (!(D0 > 0) || !(Ds > 0)) return K;

    double Db[3][3] = {{1.0, nu, 0.0},
                       {nu, 1.0, 0.0},
                       {0.0, 0.0, (1.0 - nu) * 0.5}};

    // 2×2 高斯积分。
    // 【为什么不能只用 1 点】弯曲应变场是双线性的，1 点积分会让 κ 只有
    // 3 个独立采样 → 12 个板自由度里有 4 个以上得不到约束 → 出现零能模态。
    static const double g2[2] = {-0.57735026918962584, 0.57735026918962584};

    for (int ga = 0; ga < 2; ++ga)
      for (int gb = 0; gb < 2; ++gb) {
        const double r = g2[ga], s = g2[gb];

        double N[4], dNdr[4], dNds[4];
        shapeNatural(r, s, N, dNdr, dNds);

        // Jacobian：J = ∂(x,y)/∂(r,s)
        double xr = 0, yr = 0, xs = 0, ys = 0;
        for (int i = 0; i < 4; ++i) {
          xr += dNdr[i] * local_[i].x;  yr += dNdr[i] * local_[i].y;
          xs += dNds[i] * local_[i].x;  ys += dNds[i] * local_[i].y;
        }
        const double detJ = xr * ys - yr * xs;
        if (std::abs(detJ) < 1e-30) continue;
        const double inv = 1.0 / detJ;

        // 笛卡尔导数（标准链式法则）
        double dNdx[4], dNdy[4];
        for (int i = 0; i < 4; ++i) {
          dNdx[i] = (ys * dNdr[i] - yr * dNds[i]) * inv;
          dNdy[i] = (-xs * dNdr[i] + xr * dNds[i]) * inv;
        }

        // ---------------- 弯曲项（3×12 B 矩阵） ----------------
        double Bb[3][12] = {};
        for (int j = 0; j < 4; ++j) {
          Bb[0][3 * j + 2] =  dNdx[j];        // κx  =  ∂θy/∂x
          Bb[1][3 * j + 1] = -dNdy[j];        // κy  = −∂θx/∂y
          Bb[2][3 * j + 2] =  dNdy[j];        // κxy =  ∂θy/∂y − ∂θx/∂x
          Bb[2][3 * j + 1] = -dNdx[j];
        }
        for (int i = 0; i < 12; ++i)
          for (int j = 0; j < 12; ++j) {
            double v = 0.0;
            for (int p = 0; p < 3; ++p) {
              if (Bb[p][i] == 0.0) continue;
              for (int q = 0; q < 3; ++q) {
                const double d = Db[p][q];
                if (d == 0.0 || Bb[q][j] == 0.0) continue;
                v += Bb[p][i] * d * Bb[q][j];
              }
            }
            if (v != 0.0) K(i, j) += D0 * v * detJ;
          }

        // ---------------- 剪切项（MITC4 假定应变场） ----------------
        // 四个绑定点：各边中点在自然坐标下的位置
        double rz[2][12], sz[2][12];    // [0] = 负端，[1] = 正端
        covariantShearRow(0.0, -1.0, rz[0], sz[0]);   // 下边中点
        covariantShearRow(0.0,  1.0, rz[1], sz[1]);   // 上边中点
        double rzL[12], szL[12], rzR[12], szR[12];
        covariantShearRow(-1.0, 0.0, rzL, szL);       // 左边中点
        covariantShearRow( 1.0, 0.0, rzR, szR);       // 右边中点

        // 【插值方向不能弄反】
        //   γ_rz 的采样点在 s 的两端 ⇒ 沿 s 线性插值
        //   γ_sz 的采样点在 r 的两端 ⇒ 沿 r 线性插值
        const double wA = 0.5 * (1.0 - s), wC = 0.5 * (1.0 + s);
        const double wD = 0.5 * (1.0 - r), wB = 0.5 * (1.0 + r);

        double Bxc[12] = {}, Byc[12] = {};
        for (int i = 0; i < 12; ++i) {
          const double grz = wA * rz[0][i] + wC * rz[1][i];
          const double gsz = wD * szL[i] + wB * szR[i];
          // 共变 → 笛卡尔：(γxz, γyz) = J⁻¹·(γ_rz, γ_sz)
          Bxc[i] = ( ys * grz - yr * gsz) * inv;
          Byc[i] = (-xs * grz + xr * gsz) * inv;
        }
        for (int i = 0; i < 12; ++i)
          for (int j = 0; j < 12; ++j)
            K(i, j) += Ds * (Bxc[i] * Bxc[j] + Byc[i] * Byc[j]) * detJ;
      }
    return K;
  }

  // -------------------------------------------------------------------------
  //  总刚度 K(24×24)
  //    DOF：节点 i → (ux, uy, uz, rx, ry, rz) = 6i .. 6i+5
  //
  //  【两步走：先在局部坐标组装，再做变换】
  //
  //  膜与板的刚度都是在【单元局部坐标】里推出来的（局部 x-y 在板面内），
  //  而组装到整体结构时用的是【全局自由度】。两者之间差一个旋转：
  //
  //      K_global = Tᵀ · K_local · T,      u_local = T · u_global
  //
  //  T 是每节点 6×6 的对角块 diag(R, R)，R 的三行是局部基向量
  //  （R 把全局分量投到局部：u = ex·d, v = ey·d, w = ez·d）。
  //  转动分量用同一个 R —— 转角是赝矢量，在【正常】旋转（det R = +1）下
  //  与位移向量同规律变换。
  //
  //  【这是一个真实存在过的缺陷】
  //  早期实现把全局 ux,uy 直接当局部的 u,v 用（即隐含 T = I）。
  //  后果只在"局部轴与全局轴不平行"时暴露：
  //    · 轴对齐的矩形网格 —— T 是符号置换，结果恰好正确（所以测试全过）
  //    · 斜交轴网、倾斜屋面、任意方向的墙 —— 刚度矩阵是错的，
  //      且错得不显眼：位移量级正常、平衡也成立，只是数值不对
  //  内力恢复函数一直是对的（它显式做了投影），
  //  于是出现"位移错、应力按正确公式从错误的位移里算"这种更难查的局面。
  // -------------------------------------------------------------------------
  Mat24 stiffness() const {
    Mat24 K = Mat24::zero();
    if (!valid_) return K;
    const Mat8 Km = membraneStiffness();
    const Mat12 Kb = plateStiffness();

    // 膜：局部 (2i+a) → 局部 24 自由度 (6i+a)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        for (int a = 0; a < 2; ++a)
          for (int b = 0; b < 2; ++b) {
            const double v = Km(2 * i + a, 2 * j + b);
            if (v != 0.0) K(6 * i + a, 6 * j + b) += v;
          }

    // 板：局部 (3i+a) → 局部 24 自由度 (6i+2+a)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) {
            const double v = Kb(3 * i + a, 3 * j + b);
            if (v != 0.0) K(6 * i + 2 + a, 6 * j + 2 + b) += v;
          }

    // 钻转刚度：绕板面法向（rz）
    // 【为什么需要】膜单元内"绕法向"的转动被完全抑制（面内应变与它无关），
    // 若不给刚度，该自由度就是一个"零能量"自由度 → 刚度矩阵奇异。
    // 盈建科/PKPM 的处理：给一个很小的正刚度，工程上等效于"约束"。
    if (props_.drillingStiffness > 0.0) {
      const double kr = props_.drillingStiffness * props_.E
                      * props_.thickness * area_;
      for (int i = 0; i < 4; ++i) {
        K(6 * i + 5, 6 * i + 5) += kr;
        for (int j = 0; j < 4; ++j)
          if (i != j) K(6 * i + 5, 6 * j + 5) -= kr / 3.0;
      }
    }

    return toGlobalDofs(K);
  }

  // -------------------------------------------------------------------------
  //  荷载：面外均布荷载（沿单元法向，kN/m²，正为沿 +ez）
  //
  //  【返回的是【全局】分量】—— 这一点必须说清楚。
  //  荷载本身是沿局部法向施加的，但上层组装（StaticAnalysis::scatter）
  //  直接把它当作全局自由度上的力来散装。若不在这里做变换，
  //  倾斜板的面荷载方向就错了（轴对齐单元恰好看不出来）。
  //  与刚度矩阵一样：局部推出来的东西，出这一层之前必须转回全局。
  // -------------------------------------------------------------------------
  Vec24d transversePressure(double pz) const {
    Vec24d f{};
    if (!valid_) return f;
    for (int i = 0; i < 4; ++i) f[6 * i + 2] = pz * area_ / 4.0;
    return toGlobalVector(f);
  }

  //  面内均布压力（局部 x, y 方向，kN/m²，正为压）
  Vec24d membranePressure(double px, double py) const {
    Vec24d f{};
    if (!valid_) return f;
    for (int i = 0; i < 4; ++i) {
      f[6 * i + 0] = -px * area_ / 4.0;
      f[6 * i + 1] = -py * area_ / 4.0;
    }
    return toGlobalVector(f);
  }

  // -------------------------------------------------------------------------
  //  面内内力恢复（单元中心应力）
  //
  //  ε = B·u  （先转到局部坐标），σ = D·ε
  //  返回平面内应力与主应力、von Mises
  // -------------------------------------------------------------------------
  struct MembraneForces {
    double sx{0}, sy{0}, txy{0};     // 单元局部坐标系下的应力
    double s1{0}, s2{0};             // 主应力
    double vonMises{0};
  };

  MembraneForces membraneForces(const Vec24d& uGlobal) const {
    MembraneForces r;
    if (!valid_) return r;

    // 转到局部（正交基 ⇒ 直接投影）
    double ul[8];
    for (int i = 0; i < 4; ++i) {
      const double ux = uGlobal[6 * i + 0];
      const double uy = uGlobal[6 * i + 1];
      const double uz = uGlobal[6 * i + 2];
      ul[2 * i + 0] = ux * ex_.x + uy * ex_.y + uz * ex_.z;
      ul[2 * i + 1] = ux * ey_.x + uy * ey_.y + uz * ey_.z;
    }
    double eps[3] = {0, 0, 0};
    for (int i = 0; i < 4; ++i) {
      // εx = ∂u_x/∂x,  εy = ∂u_y/∂y,  γxy = ∂u_x/∂y + ∂u_y/∂x
      eps[0] += dNdx_[i] * ul[2 * i + 0];
      eps[1] += dNdy_[i] * ul[2 * i + 1];        // ← dNdy 不是 dNdx
      eps[2] += dNdy_[i] * ul[2 * i + 0] + dNdx_[i] * ul[2 * i + 1];
    }
    const Mat3 D = props_.membraneD();
    r.sx = D(0, 0) * eps[0] + D(0, 1) * eps[1];
    r.sy = D(0, 1) * eps[0] + D(1, 1) * eps[1];
    r.txy = D(2, 2) * eps[2];

    const double c = 0.5 * (r.sx + r.sy);
    const double rad = std::sqrt(0.25 * (r.sx - r.sy) * (r.sx - r.sy) + r.txy * r.txy);
    r.s1 = c + rad;
    r.s2 = c - rad;
    r.vonMises = std::sqrt(r.sx * r.sx - r.sx * r.sy + r.sy * r.sy
                          + 3.0 * r.txy * r.txy);
    return r;
  }

  // -------------------------------------------------------------------------
  //  板弯曲内力：Mx, My, Mxy（kN·m/m）与 Qx, Qy（kN/m）
  // -------------------------------------------------------------------------
  struct PlateForces {
    double Mx{0}, My{0}, Mxy{0};
    double Qx{0}, Qy{0};
  };

  //  符号约定见 plateStiffness 的文件头（βx = θy，βy = −θx）。
  //  本构：
  //      Mx  = D₀(κx + νκy)
  //      My  = D₀(κy + νκx)
  //      Mxy = D₀(1−ν)/2 · κxy
  //      Qx  = κGt·γxz,  Qy = κGt·γyz
  //
  //  【这里曾经有两个错误，都是"形似正确"的类型】
  //  ① Mx 与 My 的公式写反了：原实现 Mx = −D₀·∂θx/∂y = D₀·κy，
  //     而 κy 属于 My。位移能算对（刚度矩阵只用到二次型，与符号约定无关），
  //     但弯矩报告出来是错的 —— 且因为是"两个方向互换"，
  //     在正方形板块上 Mx = My，根本看不出来。
  //  ② Mxy 被直接写死为 0（注释说"简化"）。
  //     缺了扭矩项，Mx/My 在斜交板上会明显偏小，
  //     更关键的是"主弯矩方向"完全错了 —— 而配筋是按主弯矩方向布置的。
  //
  //  【为什么剪力要从 MITC4 假定场取，而不是直接用 ∂w/∂x + θy】
  //  刚度矩阵用的是假定应变场。如果内力恢复换成位移插值场，
  //  两者不一致，会出现 Q ≠ dM/dx（平衡校核失败），
  //  而且薄板下 Q 会带着寄生剪切 —— 恰恰是 MITC4 要消除的那个量。
  // -------------------------------------------------------------------------
  PlateForces plateForces(const Vec24d& uGlobal) const {
    PlateForces r;
    if (!valid_) return r;

    // ---- 全局 → 局部 ----
    double w[4], tx[4], ty[4];
    for (int i = 0; i < 4; ++i) {
      const double ux = uGlobal[6 * i + 0];
      const double uy = uGlobal[6 * i + 1];
      const double uz = uGlobal[6 * i + 2];
      const double rx = uGlobal[6 * i + 3];
      const double ry = uGlobal[6 * i + 4];
      const double rz = uGlobal[6 * i + 5];
      w[i]  = ux * ez_.x + uy * ez_.y + uz * ez_.z;
      tx[i] = rx * ex_.x + ry * ex_.y + rz * ex_.z;   // 局部 θx
      ty[i] = rx * ey_.x + ry * ey_.y + rz * ey_.z;   // 局部 θy
    }

    // ---- 形心 (r = s = 0) 处的形函数与笛卡尔导数 ----
    double N[4], dr[4], ds[4];
    shapeNatural(0.0, 0.0, N, dr, ds);
    double xr = 0, yr = 0, xs = 0, ys = 0;
    for (int i = 0; i < 4; ++i) {
      xr += dr[i] * local_[i].x;  yr += dr[i] * local_[i].y;
      xs += ds[i] * local_[i].x;  ys += ds[i] * local_[i].y;
    }
    const double detJ = xr * ys - yr * xs;
    const double inv = (std::abs(detJ) > 1e-30) ? 1.0 / detJ : 0.0;
    double dNdx[4], dNdy[4];
    for (int i = 0; i < 4; ++i) {
      dNdx[i] = (ys * dr[i] - yr * ds[i]) * inv;
      dNdy[i] = (-xs * dr[i] + xr * ds[i]) * inv;
    }

    // ---- 曲率 ----
    double kx = 0, ky = 0, kxy = 0;
    for (int i = 0; i < 4; ++i) {
      kx  +=  dNdx[i] * ty[i];                              //  ∂θy/∂x
      ky  -=  dNdy[i] * tx[i];                              // −∂θx/∂y
      kxy +=  dNdy[i] * ty[i] - dNdx[i] * tx[i];            //  ∂θy/∂y − ∂θx/∂x
    }
    const double t = props_.thickness;
    const double nu = props_.nu;
    const double D0 = props_.E * t * t * t / (12.0 * (1.0 - nu * nu));
    r.Mx  = D0 * (kx + nu * ky);
    r.My  = D0 * (ky + nu * kx);
    r.Mxy = D0 * (1.0 - nu) * 0.5 * kxy;

    // ---- 横向剪力：与刚度矩阵同一套 MITC4 假定应变场 ----
    const double v[12] = {w[0], tx[0], ty[0], w[1], tx[1], ty[1],
                          w[2], tx[2], ty[2], w[3], tx[3], ty[3]};
    double rzA[12], szA[12], rzC[12], szC[12], rzD[12], szD[12], rzB[12], szB[12];
    covariantShearRow(0.0, -1.0, rzA, szA);
    covariantShearRow(0.0,  1.0, rzC, szC);
    covariantShearRow(-1.0, 0.0, rzD, szD);
    covariantShearRow( 1.0, 0.0, rzB, szB);

    // 形心处四个插值权的值都是 ½（½(1−0)）
    double grz = 0.0, gsz = 0.0;
    for (int k = 0; k < 12; ++k) {
      grz += 0.5 * (rzA[k] + rzC[k]) * v[k];
      gsz += 0.5 * (szD[k] + szB[k]) * v[k];
    }
    const double Ds = props_.shearRigidity();
    r.Qx = Ds * ( ys * grz - yr * gsz) * inv;
    r.Qy = Ds * (-xs * grz + xr * gsz) * inv;
    return r;
  }

  // ---- 查询 ----
  const Vec3& localX() const { return ex_; }
  const Vec3& localY() const { return ey_; }
  const Vec3& localZ() const { return ez_; }
  const Vec3& centroid() const { return centroid_; }
  const Vec3* nodes() const { return p_; }

  // 自重 kN（与其它单元的 mass() 同单位）
  //
  //  【这里曾漏了密度】原实现返回 area×thickness —— 那是【体积 m³】，
  //  既不是质量也不是重量。后果是总重量统计差 1000 倍以上，
  //  剪重比（V/W）被严重低估，而剪重比是抗震的【下限】指标。
  //  density 单位是 kg/m³ ⇒ 乘 9.81/1000 换成 kN/m³（2500 → 24.5，≈ 混凝土容重 25）。
  double mass() const { return area_ * props_.thickness * props_.density * 9.81 / 1000.0; }

 private:
  // 局部 → 全局（向量形式）：f_g = Tᵀ · f_l
  //
  //  与刚度矩阵同一个 T：T 是每节点 6×6 的 diag(R, R)。
  //  荷载是"余向量"（力的分量），所以用 Tᵀ 而不是 T ——
  //  这正是虚功原理的要求：f_lᵀ·u_l = f_lᵀ·(T·u_g) = (Tᵀ·f_l)ᵀ·u_g。
  Vec24d toGlobalVector(const Vec24d& fl) const {
    const double R[3][3] = {{ex_.x, ex_.y, ex_.z},
                            {ey_.x, ey_.y, ey_.z},
                            {ez_.x, ez_.y, ez_.z}};
    Vec24d fg{};
    for (int n = 0; n < 4; ++n) {
      for (int c = 0; c < 3; ++c) {
        // 平动：f_g[p] = Σ_c R[c][p] · f_l[c]
        fg[6 * n + c] = R[0][c] * fl[6 * n + 0] + R[1][c] * fl[6 * n + 1] +
                        R[2][c] * fl[6 * n + 2];
        fg[6 * n + 3 + c] = R[0][c] * fl[6 * n + 3] + R[1][c] * fl[6 * n + 4] +
                            R[2][c] * fl[6 * n + 5];
      }
    }
    return fg;
  }

  // 局部 → 全局：K_g = Tᵀ · K_l · T
  //
  //  T 是每节点 6×6 的 diag(R, R)，R 的三行是局部基向量（ex, ey, ez）。
  //  这里不显式构造 24×24 的 T（那是 4608 个 double，白占栈），
  //  而是按"每节点 6 自由度块"直接做两次缩并：
  //      K_g(gi,gj) = Σ_c Σ_d Tn(p,c) · K_l(li,lj) · Tn(q,d)
  Mat24 toGlobalDofs(const Mat24& Kl) const {
    const double R[3][3] = {{ex_.x, ex_.y, ex_.z},
                            {ey_.x, ey_.y, ey_.z},
                            {ez_.x, ez_.y, ez_.z}};
    // Tn[局部][全局] = ∂l_局部 / ∂g_全局
    //
    //  l = Tn · g  ——【行是局部、列是全局】。这个方向很容易搞反，
    //  而搞反的后果是算成 T·K·Tᵀ 而不是 Tᵀ·K·T：
    //    · 轴对齐单元上 T 是符号置换且 T = Tᵀ ⇒ 看不出问题
    //    · 一旦有真实旋转，刚度矩阵就是一个"转错方向"的矩阵，
    //      表现为刚体转动产生巨大虚假内力（残差 O(1)，不是舍入误差）
    double Tn[6][6] = {};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        Tn[a][b] = R[a][b];              // 局部平动 ← 全局平动
        Tn[3 + a][3 + b] = R[a][b];      // 局部转动 ← 全局转动
      }

    // K_g(i,j) = Σ_{a,b} T(a,i) · K_l(a,b) · T(b,j)
    //   a,b 为【局部】分量下标；i,j 为【全局】分量下标
    Mat24 Kg = Mat24::zero();
    for (int ni = 0; ni < 4; ++ni)
      for (int nj = 0; nj < 4; ++nj)
        for (int iG = 0; iG < 6; ++iG)          // 全局分量
          for (int jG = 0; jG < 6; ++jG) {
            double v = 0.0;
            for (int aL = 0; aL < 6; ++aL) {    // 局部分量
              const double t1 = Tn[aL][iG];
              if (t1 == 0.0) continue;
              for (int bL = 0; bL < 6; ++bL) {
                const double t2 = Tn[bL][jG];
                if (t2 == 0.0) continue;
                const double kk = Kl(6 * ni + aL, 6 * nj + bL);
                if (kk == 0.0) continue;
                v += t1 * kk * t2;
              }
            }
            if (v != 0.0) Kg(6 * ni + iG, 6 * nj + jG) = v;
          }
    return Kg;
  }

  // 双线性形函数及其【自然坐标】导数，一次算齐（板部分要在任意 (r,s) 处求值）
  static void shapeNatural(double r, double s, double N[4], double dNdr[4], double dNds[4]) {
    N[0] = 0.25 * (1.0 - r) * (1.0 - s);
    N[1] = 0.25 * (1.0 + r) * (1.0 - s);
    N[2] = 0.25 * (1.0 + r) * (1.0 + s);
    N[3] = 0.25 * (1.0 - r) * (1.0 + s);
    dNdr[0] = -0.25 * (1.0 - s); dNdr[1] = 0.25 * (1.0 - s);
    dNdr[2] =  0.25 * (1.0 + s); dNdr[3] = -0.25 * (1.0 + s);
    dNds[0] = -0.25 * (1.0 - r); dNds[1] = -0.25 * (1.0 + r);
    dNds[2] =  0.25 * (1.0 + r); dNds[3] =  0.25 * (1.0 - r);
  }

  // -------------------------------------------------------------------------
  //  MITC4 绑定点上的【共变剪切应变】系数行
  //
  //  共变分量 = 笛卡尔剪切向量在自然坐标基上的投影：
  //      γ_rz = g_r·(γxz, γyz) = (∂x/∂r)·γxz + (∂y/∂r)·γyz
  //  代入 γxz = ∂w/∂x + θy、γyz = ∂w/∂y − θx，并注意
  //      ∂w/∂r = (∂x/∂r)(∂w/∂x) + (∂y/∂r)(∂w/∂y)
  //  可得
  //      γ_rz = Σ_i (∂Ni/∂r)·w_i + (∂x/∂r)·θy − (∂y/∂r)·θx
  //  其中 θx、θy 用形函数插值到该点。
  //
  //  【这里没有"½ 的对称化"】网上流传的 MITC4 写法常把这一项写成
  //      γxz(A) = ½[ (θy2+θy4) + (w2−w1)/L1 + (w3−w4)/L3 ]
  //  那种形式是把"边上的 ∂w/∂r"直接写成两端差商、并取两条边的平均。
  //  本实现用形函数在该点的解析导数，与上式在矩形上等价，
  //  但在畸变四边形上不会引入"哪条边配哪个节点"的歧义。
  //
  //  rowRz / rowSz 各 12 项，对应板自由度 (w, θx, θy) × 4 节点。
  // -------------------------------------------------------------------------
  void covariantShearRow(double rp, double sp, double rowRz[12], double rowSz[12]) const {
    double N[4], dr[4], ds[4];
    shapeNatural(rp, sp, N, dr, ds);

    double xr = 0, yr = 0, xs = 0, ys = 0;
    for (int i = 0; i < 4; ++i) {
      xr += dr[i] * local_[i].x;  yr += dr[i] * local_[i].y;
      xs += ds[i] * local_[i].x;  ys += ds[i] * local_[i].y;
    }
    for (int i = 0; i < 12; ++i) { rowRz[i] = 0.0; rowSz[i] = 0.0; }
    for (int j = 0; j < 4; ++j) {
      rowRz[3 * j + 0] = dr[j];             // ∂w/∂r
      rowRz[3 * j + 1] = -yr * N[j];        // −(∂y/∂r)·θx
      rowRz[3 * j + 2] =  xr * N[j];        // +(∂x/∂r)·θy
      rowSz[3 * j + 0] = ds[j];             // ∂w/∂s
      rowSz[3 * j + 1] = -ys * N[j];
      rowSz[3 * j + 2] =  xs * N[j];
    }
  }

  // 4 节点双线性单元的形函数与导数
  //
  //   N1 = (1−ξ)(1−η)/4    N2 = (1+ξ)(1−η)/4
  //   N3 = (1+ξ)(1+η)/4    N4 = (1−ξ)(1+η)/4
  //
  //  【必须用【单元中心】处的导数值 —— 这里曾经差了整 2 倍】
  //
  //  对双线性单元，∂N/∂ξ 在单元内是【线性变化】的：
  //      ∂N1/∂ξ = −(1−η)/4
  //  它在"角点"处的值是 ∓0.5（η = ∓1），在"中心"处的值是 ∓0.25（η = 0）。
  //  常应变（CST 型）的 B 矩阵要的是【中心值】∓0.25。
  //
  //  早期实现取了角点值 ∓0.5，于是：
  //    · dN/dx 大 2 倍 ⇒ 膜应变大 2 倍 ⇒ 膜刚度【大 4 倍】
  //    · 应力恢复也大 2 倍
  //    · 刚体位移校核仍然通过（Σ∂Ni/∂x = 0 与整体缩放无关）
  //    · 而当时的单元测试把"dNdx 幅值 = 1/a"当成了期望值写进断言 ——
  //      错误的值被测试"确认"了，于是这个 bug 一直被掩盖
  //
  //  独立复核办法（不依赖任何推导，只看形函数本身）：对 N_i 在单元中心做
  //  数值微分。矩形 a×b 的正确答案是
  //      ∂N_i/∂x = ∓1/(2a)，  ∂N_i/∂y = ∓1/(2b)
  //  直观理解：沿 x 方向 N_i 从 1 线性降到 0，跨越的是【整个宽度 a】，
  //  所以斜率是 1/a —— 中心差分给出的正是这个值。
  //  角点值 ∓0.5 相当于"用半个单元的宽度去除整个落差"，本身就不成立。
  //
  //  Jacobian 用的 Σ(∂N_i/∂ξ)·x_i 本来就是中心处的求和，
  //  所以两处一致，不存在"两个系数要分开处理"的问题。
  //
  //  校核：矩形单元 a×b 的 dNdx = ∓1/(2a)，dNdy = ∓1/(2b)
  void computeShapeGradients() {
    // 形函数对自然坐标的导数，【取单元中心 ξ = η = 0 处的值】
    static const double dNdXiCtr[4]  = {-0.25, 0.25, 0.25, -0.25};
    static const double dNdEtaCtr[4] = {-0.25, -0.25, 0.25, 0.25};
    // 四个角点的自然坐标（Jacobian 在其上用同一条公式求和，等价于中心值）
    static const double xi_i[4]  = {-1, 1, 1, -1};
    static const double eta_i[4] = {-1, -1, 1, 1};

    // Jacobian：∂x/∂ξ = (1/4)[−x1 + x2 + x3 − x4]，其余同理
    double xxi = 0, xeta = 0, yxi = 0, yeta = 0;
    for (int i = 0; i < 4; ++i) {
      xxi  += xi_i[i]  * local_[i].x / 4.0;
      xeta += eta_i[i] * local_[i].x / 4.0;
      yxi  += xi_i[i]  * local_[i].y / 4.0;
      yeta += eta_i[i] * local_[i].y / 4.0;
    }
    const double detJ = xxi * yeta - xeta * yxi;
    const double inv = (std::abs(detJ) > 1e-30) ? 1.0 / detJ : 0.0;

    // ∂ξ/∂x = y_η/detJ,  ∂η/∂x = −y_ξ/detJ
    // ∂ξ/∂y = −x_η/detJ, ∂η/∂y = x_ξ/detJ
    for (int i = 0; i < 4; ++i) {
      dNdx_[i] = dNdXiCtr[i]  * yeta * inv - dNdEtaCtr[i] * yxi * inv;
      dNdy_[i] = -dNdXiCtr[i] * xeta * inv + dNdEtaCtr[i] * xxi * inv;
    }
  }

  Vec3 p_[4];
  Vec3 local_[4];
  Vec3 centroid_;
  Vec3 ex_{1, 0, 0}, ey_{0, 1, 0}, ez_{0, 0, 1};
  double dNdx_[4] = {0, 0, 0, 0};
  double dNdy_[4] = {0, 0, 0, 0};
  double area_{0.0};
  ShellProperties props_;
  bool valid_{false};
  bool clockwise_{false};
  std::string error_;
};

}  // namespace yjk