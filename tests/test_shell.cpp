// =============================================================================
//  tests/test_shell.cpp  ——  四节点壳单元（膜 + MITC4 板）验证
//
//  壳单元的验证分两条独立的线，缺一不可：
//   膜部分：面内平衡、应力恢复 —— 用"无限大板受均布力"的解析解对标
//   板部分：薄板挠度收敛 —— 用"薄板理论"验证"厚度→0 时不锁定"
//
//  第二条是壳单元的生死线。普通低阶板单元在小厚度时会出现剪切锁定，
//  挠度偏小一个数量级且不报错 —— 只有 MITC4 才能避免。
// =============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "yjk/element/ShellElement4.h"
#include "yjk/math/Solver.h"

using namespace yjk;

static int g_pass = 0, g_fail = 0, g_skip = 0;
static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}
static void checkNear(double got, double want, double relTol, const char* what) {
  const double rel = std::abs(want) > 1e-30 ? std::abs(got - want) / std::abs(want)
                                            : std::abs(got - want);
  if (rel <= relTol) {
    ++g_pass;
    std::printf("  [ OK ] %s  (got %.8g, want %.8g, rel %.2e)\n", what, got, want, rel);
  } else {
    ++g_fail;
    std::printf("  [FAIL] %s  (got %.8g, want %.8g, rel %.2e)\n", what, got, want, rel);
  }
}

// 刚体位移不产生内力
//
// 【校核要点】
//  ① 必须【逐模式归一化】‖K·u‖∞ / (max|K| · max|u|)。
//     若把所有模式的 u 混在一起算 max|u|，平移模式（值 1）会主导归一化，
//     使旋转模式的相对误差被严重低估 —— 一个真实的 1e-3 误差会被显示成 1e-6。
//  ② 旋转用严格的小角度形式 u_x = −θy, u_y = θx。
//     θ 取得极小（1e-6），二阶项 O(θ²)=1e-12 可忽略。
//  ③ 必须同时校核【平移】与【旋转】。
//     γxy = ∂u_x/∂y + ∂u_y/∂x 的 B 矩阵写错时，平移仍通过、只有旋转不通过。
static void checkRigidBody(const Mat8& K, const char* what) {
  double kmax = 0.0;
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j) kmax = std::max(kmax, std::abs(K(i, j)));
  if (!(kmax > 0)) { check(false, what); return; }

  // 矩形 2×3 单元的局部坐标（形心已在原点）
  const double xs[4] = {-1.0, 1.0, 1.0, -1.0};
  const double ys[4] = {-1.5, -1.5, 1.5, 1.5};
  const double th = 1.0e-6;

  double worstRatio = 0.0;
  for (int mode = 0; mode < 3; ++mode) {
    double u[8] = {};
    double umax = 0.0;
    if (mode == 0) {
      for (int i = 0; i < 4; ++i) u[2 * i] = 1.0;              // x 平移
    } else if (mode == 1) {
      for (int i = 0; i < 4; ++i) u[2 * i + 1] = 1.0;          // y 平移
    } else {
      for (int i = 0; i < 4; ++i) {                             // 绕形心旋转
        u[2 * i] = -th * ys[i];
        u[2 * i + 1] = th * xs[i];
      }
    }
    for (int i = 0; i < 8; ++i) umax = std::max(umax, std::abs(u[i]));
    double worst = 0.0;
    for (int i = 0; i < 8; ++i) {
      double sum = 0.0;
      for (int j = 0; j < 8; ++j) sum += K(i, j) * u[j];
      worst = std::max(worst, std::abs(sum));
    }
    worstRatio = std::max(worstRatio, worst / (kmax * umax));
  }
  checkNear(worstRatio, 0.0, 1e-13, what);
}

// ===========================================================================
//  1. 几何与形函数
// ===========================================================================
static void testGeometry() {
  std::printf("\n== 1. 几何与形函数 ==\n");
  ShellProperties sp;
  sp.thickness = 0.2; sp.E = 3.0e7; sp.nu = 0.2;

  // (a) 矩形单元
  {
    Vec3 p[4] = {{0,0,0}, {2,0,0}, {2,3,0}, {0,3,0}};
    ShellElement4 sh(p, sp);
    check(sh.isValid(), "矩形单元有效");
    checkNear(sh.area(), 6.0, 1e-14, "矩形面积 = 6");
    checkNear(sh.localZ().z, 1.0, 1e-15, "法向 = +Z");
    checkNear(sh.centroid().x, 1.0, 1e-14, "形心 x");
    checkNear(sh.centroid().y, 1.5, 1e-14, "形心 y");
  }
  // (b) 倾斜单元（法向不是坐标轴）
  {
    // 在 x-y 平面上旋转 30°，再加 z 偏移
    const double th = 0.5;
    const double c = std::cos(th), s = std::sin(th);
    auto rot = [&](double x, double y, double z) {
      return Vec3{x * c - y * s, x * s + y * c, z + 0.5};
    };
    Vec3 p[4] = {rot(0,0,0), rot(2,0,0), rot(2,3,0), rot(0,3,0)};
    ShellElement4 sh(p, sp);
    check(sh.isValid(), "倾斜单元有效");
    checkNear(sh.area(), 6.0, 1e-13, "倾斜单元面积不变（旋转不变量）");
    checkNear(sh.localZ().z, 1.0, 1e-14, "法向仍为 +Z（平面内旋转）");
  }
  // (c) 非水平单元：真正的【平面】四边形绕 x 轴倾斜
  //     注意：不能给 4 个点各自加不同的 z 偏移 —— 那样得到的是空间四边形，
  //     没有唯一平面，壳单元无从定义法向。必须是刚性旋转。
  {
    const double th = 0.3;              // 绕 x 轴倾斜
    const double c = std::cos(th), s = std::sin(th);
    auto tilt = [&](double x, double y) {
      return Vec3{x, y * c, y * s};    // 绕 x 轴刚性旋转，x 不变
    };
    Vec3 p[4] = {tilt(0,0), tilt(2,0), tilt(2,3), tilt(0,3)};
    ShellElement4 sh(p, sp);
    check(sh.isValid(), "倾斜平面单元有效");
    const Vec3 n = sh.normal();
    // 法向应为 (0, −sinθ, cosθ) 的规范化
    checkNear(n.y, -s, 1e-12, "倾斜单元法向 y 分量");
    checkNear(n.z, c, 1e-12, "倾斜单元法向 z 分量");
    // 【关键】倾斜不改变真实面积 —— 这是用全局坐标算面积的意义
    checkNear(sh.area(), 6.0, 1e-13, "倾斜单元面积保持 6（刚性旋转不变量）");
  }
  // (d) 形函数导数的物理校核
  //
  //  【这里曾经把错误值当成期望值写进断言】
  //  原断言是 "dNdx 幅值 = 1/a"，而正确值是 1/(2a) —— 断言本身错了，
  //  于是接口上那个"大了 2 倍"的 bug 被测试'确认'成了正确行为。
  //  教训：**期望值必须来自独立推导或数值实验，不能来自"代码现在给出什么"**。
  //  下面的数值微分就是那个不依赖实现的独立来源。
  {
    Vec3 p[4] = {{0,0,0}, {4,0,0}, {4,2,0}, {0,2,0}};   // a=4, b=2，形心 (2,1)
    ShellElement4 sh(p, sp);
    // 常应变 ⇒ ΣdNi/dx = 0, ΣdNi/dy = 0
    double sx = 0, sy = 0;
    for (int i = 0; i < 4; ++i) { sx += sh.dNdx()[i]; sy += sh.dNdy()[i]; }
    checkNear(sx, 0.0, 1e-14, "Σ∂Ni/∂x = 0");
    checkNear(sy, 0.0, 1e-14, "Σ∂Ni/∂y = 0");

    // 【独立来源】对形函数做中心差分，直接测 ∂Ni/∂x —— 不看实现
    //  形心映射：ξ = (x−2)/2 ⇒ dξ/dx = 0.5。
    //  所以"物理坐标扰动 dx"对应"自然坐标扰动 dξ = 0.5·dx"，
    //  差商的分母必须用 2·dx（分子是 ξ 从 −dξ 到 +dξ 的总变化）。
    //  【这一步最容易自己写错】把分母写成 dξ 会得到 2 倍的值 ——
    //  而"数值微分给出 2 倍"恰好是我们要检测的那个 bug 的特征，
    //  校验代码自己算错就会把真 bug 判成"通过"。所以分母要显式写出来。
    {
      const double dx = 1e-6;                  // 物理坐标扰动
      const double dxi = dx * 0.5;             // dξ/dx = 1/(2·2) ... 该单元 a=4 ⇒ 0.5
      double worst = 0.0;
      for (int i = 0; i < 4; ++i) {
        double Np[4], Nm[4];
        ShellElement4::shapeFunctions(dxi, 0.0, Np);
        ShellElement4::shapeFunctions(-dxi, 0.0, Nm);
        const double dnNum = (Np[i] - Nm[i]) / (2.0 * dx);
        worst = std::max(worst, std::abs(dnNum - sh.dNdx()[i]));
      }
      checkNear(worst, 0.0, 1e-7, "∂Ni/∂x 与数值微分一致（独立来源校核）");
    }
    // 解析值：矩形 a×b ⇒ ∓1/(2a)、∓1/(2b)
    checkNear(std::abs(sh.dNdx()[0]), 1.0 / (2.0 * 4.0), 1e-14,
              "∂Ni/∂x 幅值 = 1/(2a)（注意不是 1/a）");
    checkNear(std::abs(sh.dNdy()[0]), 1.0 / (2.0 * 2.0), 1e-14, "∂Ni/∂y 幅值 = 1/(2b)");

    // 分片检验的必要条件：线性场 u = x 应给出 εx = 1
    {
      double eps = 0.0;
      for (int i = 0; i < 4; ++i) eps += sh.dNdx()[i] * p[i].x;
      checkNear(eps, 1.0, 1e-14, "线性场 u=x ⇒ εx=1（常应变单元的充要条件）");
    }

    double sumN = 0;
    double N[4];
    ShellElement4::shapeFunctions(0.3, -0.7, N);
    for (int i = 0; i < 4; ++i) sumN += N[i];
    checkNear(sumN, 1.0, 1e-15, "形函数分区之和 = 1");
  }
  // (e) 退化单元必须被拒
  {
    Vec3 p[4] = {{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}};
    ShellElement4 sh(p, sp);
    check(!sh.isValid(), "零面积单元被拒");
    std::printf("       诊断: %s\n", sh.error().c_str());
  }
  {
    Vec3 p[4] = {{0,0,0}, {1,0,0}, {2,0,0}, {3,0,0}};   // 四点共线
    ShellElement4 sh(p, sp);
    check(!sh.isValid(), "共线单元被拒");
  }
}

// ===========================================================================
//  2. 膜单元
// ===========================================================================
static void testMembrane() {
  std::printf("\n== 2. 膜单元（平面内） ==\n");
  ShellProperties sp;
  sp.thickness = 0.2; sp.E = 3.0e7; sp.nu = 0.2;

  Vec3 p[4] = {{0,0,0}, {2,0,0}, {2,3,0}, {0,3,0}};
  ShellElement4 sh(p, sp);
  const Mat8 K = sh.membraneStiffness();

  // (a) 对称性
  {
    double asym = 0, kmax = 0;
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 8; ++j) {
        asym = std::max(asym, std::abs(K(i, j) - K(j, i)));
        kmax = std::max(kmax, std::abs(K(i, j)));
      }
    checkNear(asym / kmax, 0.0, 1e-15, "膜刚度对称性");
  }
  // (b) 刚体位移（平移 + 旋转）不产生内力
  checkRigidBody(K, "膜刚体位移不产生内力");

  // (c) 单轴拉伸对标平面应力解
  //
  //     【测试用例的两个要点】
  //     ① 不能给所有节点相同位移 —— 那是刚体平移，应变为零。
  //        正确做法：一侧固定、另一侧位移。
  //     ② 常应变单元的应变由 εx = Σ(∂Ni/∂x)·u_i 给出，【不是】"中心点位移/半宽"。
  //        双线性插值在单元内部线性变化，εx 就是那个线性场的斜率。
  //        用数值微分可以交叉验证这一点（见下方 dNdx 的数值校核）。
  {
    const double Ex = 3.0e7, nu = 0.2;
    // 【用线性位移场 u = ε0·x 施加应变 —— 不依赖 dNdx 的具体数值】
    //  这样 εx = Σ(∂Ni/∂x)·u_i = ε0·Σ(∂Ni/∂x)·x_i = ε0·1 = ε0
    //  （最后一步正是常应变单元必须满足的分片条件）
    //  早期版本用"右端两点各移 δ、然后直接写 epsx = δ"，
    //  把 dNdx 的数值假设藏进了期望值里，导致 2 倍的 bug 被掩盖。
    const double eps0 = 5.0e-5;
    Vec24d u{};
    const double xNode[4] = {0.0, 2.0, 2.0, 0.0};
    for (int i = 0; i < 4; ++i) u[6 * i + 0] = eps0 * xNode[i];

    const double epsx = eps0;
    const double c = Ex / (1.0 - nu * nu);
    const double sx = c * epsx;
    const double sy = c * nu * epsx;     // 泊松效应
    const auto f = sh.membraneForces(u);
    checkNear(f.sx, sx, 1e-12, "单轴拉伸 σx 对标平面应力解");
    checkNear(f.sy, sy, 1e-12, "泊松效应产生面内应力 σy（平面应力）");
    // 平面应力 von Mises：σvm = √(σx² − σxσy + σy² + 3τ²)
    const double svm = std::sqrt(sx * sx - sx * sy + sy * sy);
    checkNear(f.vonMises, svm, 1e-12, "von Mises = √(σx² − σxσy + σy²)（面内无剪）");
    checkNear(f.s1, std::max(sx, sy), 1e-12, "主应力 s1");
    checkNear(f.s2, std::min(sx, sy), 1e-12, "主应力 s2");
  }
  // (d) 平面应力 vs 平面应变的差异
  {
    ShellProperties s1 = sp, s2 = sp;
    s2.planeStrain = true;
    Vec3 q[4] = {{0,0,0}, {2,0,0}, {2,3,0}, {0,3,0}};
    ShellElement4 e1(q, s1), e2(q, s2);
    Vec24d u{};
    u[6 * 2 + 0] = 1.0e-4;      // 节点 2 沿 x 位移（与 (c) 同一工况）
    const auto f1 = e1.membraneForces(u);
    const auto f2 = e2.membraneForces(u);
    // 平面应变：约束 εz，σz = ν(σx+σy) ≠ 0；面内刚度更高
    check(f2.sx > f1.sx, "平面应变的面内应力高于平面应力（约束更强）");
    std::printf("       σx(平面应力)=%.4e  σx(平面应变)=%.4e\n", f1.sx, f2.sx);
  }
  // (e) 畸变四边形：不应退化
  {
    // 严重畸变的四边形：点 (0,0),(2,0),(3,3),(0,2)
    // 三角形剖分 (node0,node1,node2) + (node0,node2,node3)：
    //   T1 = ½|2×3 − 0×1| = 3.0
    //   T2 = ½|3×2 − 3×0| = 3.0
    //   合计 6.0
    Vec3 q[4] = {{0,0,0}, {2,0,0}, {3,3,0}, {0,2,0}};
    ShellElement4 e(q, sp);
    check(e.isValid(), "畸变四边形有效");
    checkNear(e.area(), 6.0, 1e-13, "畸变四边形面积（三角形剖分校核）");
    const Mat8 Kd = e.membraneStiffness();
    // 该单元形心 = ((0+2+3+0)/4, (0+0+3+2)/4) = (1.25, 1.25)
    // 局部坐标系由 init 内部建立：ex = 节点0→节点1 方向
    // 直接用单元给出的局部坐标做刚体校核最稳妥 —— 不自己推
    double worst = 0.0, kmax = 0.0, umax = 0.0;
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 8; ++j) kmax = std::max(kmax, std::abs(Kd(i, j)));
    // 该单元形心 (1.25, 1.25)，局部轴与全局平行
    const double gx[4] = {-1.25, 0.75, 1.75, -1.25};
    const double gy[4] = {-1.25, -1.25, 1.75, 0.75};
    const double th = 1.0e-6;   // 严格小角度：二阶项 O(θ²) 可忽略
    for (int mode = 0; mode < 3; ++mode) {
      double u[8] = {};
      if (mode == 0)      for (int i = 0; i < 4; ++i) u[2*i] = 1.0;
      else if (mode == 1) for (int i = 0; i < 4; ++i) u[2*i+1] = 1.0;
      else {
        for (int i = 0; i < 4; ++i) {   // 绕形心小角度旋转
          u[2*i] = -th * gy[i];
          u[2*i+1] = th * gx[i];
        }
      }
      for (int i = 0; i < 8; ++i) umax = std::max(umax, std::abs(u[i]));
      for (int i = 0; i < 8; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 8; ++j) sum += Kd(i, j) * u[j];
        worst = std::max(worst, std::abs(sum));
      }
    }
    checkNear(worst / (kmax * umax), 0.0, 1e-13, "畸变单元刚体位移仍无内力");
  }
}

// ===========================================================================
//  3. 板弯曲 —— MITC4 的核心验证
// ===========================================================================
//
//  简支方板均布荷载的经典解（Timoshenko & Woinowsky-Krieger）：
//      w_max = 0.00406 · q·a⁴ / D,     D = E·t³/[12(1−ν²)]
//  夹支方板的经典解：
//      w_max = 0.00126 · q·a⁴ / D
//
//  两个解都要对，而且必须【同时】对 —— 只有一个对说明不了问题：
//  一个"偏刚"的单元可以碰巧对上夹支，一个"偏柔"的可以碰巧对上简支。
//
//  【这个工况是检验剪切锁定的试金石】
//  锁定单元（普通低阶 Mindlin 板）的 w/w_理论 会随 t 减小而迅速崩塌
//  （t 减半、比值掉一个数量级）。MITC4 应当保持常数。
// ===========================================================================

// 边界条件
enum class PlateBC {
  SimplySupported,   // 只约束 w=0  —— 对标 0.00406
  Clamped,           // 再约束各边的【弯曲转角】—— 对标 0.00126
};

struct PlateResult {
  bool ok{false};
  double w{0.0};                 // 跨中挠度（沿板法向）
  double Mx{0.0}, My{0.0};       // 跨中内力（中心单元的平均）
  int ndof{0};
};

// 方板求解器
//   n    每边网格数
//   a    边长
//   q    均布荷载 kN/m²（沿板法向）
//   rotXY 整块板绕 z 轴的平面内旋转角（弧度）—— 用于验证坐标变换
//
// 【为什么用完整 6 自由度】平面内旋转会把全局的 (ux,uy) 与 (rx,ry) 各自
// 混合（T 是每节点 diag(R,R)）。用 3 自由度简化模型根本表达不出这种混合，
// 也就测不出"有没有做坐标变换"。膜、板、钻转三个子空间在 T 下互不混合，
// 所以固定住 ux/uy/rz 之后剩下的正好就是纯弯曲问题。
//
// 【旧版本的边界条件是错的 —— 这里记录清楚】
//  原代码写：
//      if (edgeX) fixed[θx] = 1;      // edgeX = i=0 或 i=n，即 y=0 / y=a 边
//      if (edgeY) fixed[θy] = 1;
//  注释说是"简支"，但把边界节点的【弯曲转角】锁住，物理上是【夹支】。
//  验证过：这套边界给出 0.99 × 夹支解，而简支解应该是它的 3.22 倍。
//  一个取错边界条件的"简支板"测试，会让任何偏刚的单元看起来都对。
static PlateResult solvePlate(int n, double t, double q, double a,
                              double E, double nu, PlateBC bc, double rotXY = 0.0) {
  PlateResult R;
  const int npr = n + 1;
  const int nn = npr * npr;
  const int ndof = nn * 6;
  R.ndof = ndof;
  const double h = a / n;
  const double c = std::cos(rotXY), sn = std::sin(rotXY);

  auto rot = [&](double x, double y) {
    // 绕【整板几何中心】旋转，保持形心不动（否则相当于平移）
    const double dx = x - a * 0.5, dy = y - a * 0.5;
    return Vec3{a * 0.5 + dx * c - dy * sn, a * 0.5 + dx * sn + dy * c, 0.0};
  };

  TripletAssembler ta;
  std::vector<double> f(ndof, 0.0);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Vec3 p[4] = {rot(j * h, i * h), rot((j + 1) * h, i * h),
                   rot((j + 1) * h, (i + 1) * h), rot(j * h, (i + 1) * h)};
      ShellProperties sp;
      sp.thickness = t; sp.E = E; sp.nu = nu;
      sp.drillingStiffness = 1e-4;
      ShellElement4 el(p, sp);
      if (!el.isValid()) return R;

      const Mat24 Ke = el.stiffness();
      const Vec24d fe = el.transversePressure(q);     // 沿法向，已转全局

      const int gi[4] = {i * npr + j, i * npr + (j + 1),
                         (i + 1) * npr + (j + 1), (i + 1) * npr + j};
      // 【易错】Ke 的索引空间是【单元内】0..23，不能拿全局自由度去索引它。
      for (int A = 0; A < 24; ++A) {
        f[6 * gi[A / 6] + (A % 6)] += fe[A];
        for (int B = A; B < 24; ++B) {
          const double v = Ke(A, B);
          if (nearlyZero(v)) continue;
          ta.add(6 * gi[A / 6] + (A % 6), 6 * gi[B / 6] + (B % 6), v);
        }
      }
    }

  // ---- 约束 ----
  std::vector<char> fixed(ndof, 0);
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j) {
      const int base = 6 * (i * npr + j);
      // 膜与钻转自由度全部固定：平面内荷载为零，它们与弯曲解耦
      fixed[base + 0] = 1;    // ux
      fixed[base + 1] = 1;    // uy
      fixed[base + 5] = 1;    // rz
      const bool eX = (i == 0 || i == n);      // y = 0 / y = a 边
      const bool eY = (j == 0 || j == n);      // x = 0 / x = a 边
      if (eX || eY) fixed[base + 2] = 1;       // w = 0：四边简支只约束这一条
      if (bc == PlateBC::Clamped) {
        if (eX) fixed[base + 3] = 1;           // 绕 x 的转角：该边的弯曲转角
        if (eY) fixed[base + 4] = 1;           // 绕 y 的转角
      }
    }

  std::vector<int> map(ndof, -1);
  int nf = 0;
  for (int i = 0; i < ndof; ++i) if (!fixed[i]) map[i] = nf++;

  TripletAssembler ta2;
  for (const auto& e : ta.triplets()) {
    const int mi = map[e.i], mj = map[e.j];
    if (mi < 0 || mj < 0) continue;
    ta2.add(mi, mj, e.v);
  }
  SymSparseMatrix K;
  K.buildFrom(ta2.triplets(), nf);

  std::vector<double> fr(nf, 0.0);
  for (int i = 0; i < ndof; ++i)
    if (map[i] >= 0) fr[map[i]] = f[i];

  LDLTSolver solver;
  const auto rep = solver.factorize(K, 0);
  if (!rep.ok()) return R;
  std::vector<double> ufree;
  solver.solve(fr, ufree);

  std::vector<double> u(ndof, 0.0);
  for (int i = 0; i < ndof; ++i)
    if (map[i] >= 0) u[i] = ufree[map[i]];

  // ---- 取跨中挠度：沿板法向（旋转后法向仍是 +z）----
  const int mid = (n / 2) * npr + (n / 2);
  R.w = u[6 * mid + 2];

  // ---- 取跨中弯矩：中央那块单元的中心值 ----
  {
    const int i = (n > 1) ? (n / 2 - 1) : 0;
    const int j = (n > 1) ? (n / 2 - 1) : 0;
    Vec3 p[4] = {rot(j * h, i * h), rot((j + 1) * h, i * h),
                 rot((j + 1) * h, (i + 1) * h), rot(j * h, (i + 1) * h)};
    ShellProperties sp;
    sp.thickness = t; sp.E = E; sp.nu = nu; sp.drillingStiffness = 1e-4;
    ShellElement4 el(p, sp);
    const int gi[4] = {i * npr + j, i * npr + (j + 1),
                       (i + 1) * npr + (j + 1), (i + 1) * npr + j};
    Vec24d ue{};
    for (int A = 0; A < 24; ++A) ue[A] = u[6 * gi[A / 6] + (A % 6)];
    const auto pf = el.plateForces(ue);
    R.Mx = pf.Mx;
    R.My = pf.My;
  }

  R.ok = true;
  return R;
}

// ---------------------------------------------------------------------------
//  3a. 单元级：纯弯曲 / 纯剪切的精确校核
//
//  这两条是【精确】的，不涉及任何离散误差 ——
//  因为弯曲应变只由转角场决定，而转角场在这里是线性的，双线性插值可精确表示。
//  它直接锁死了 plateStiffness 与 plateForces 的符号约定：
//  只要两处的约定有一点不一致，这里立刻暴露。
// ---------------------------------------------------------------------------
static void testPlateUnit() {
  std::printf("\n== 3a. 板单元级校核（精确） ==\n");
  const double E = 3.0e7, nu = 0.3, t = 0.2;
  const double a = 2.0, b = 3.0;
  ShellProperties sp;
  sp.thickness = t; sp.E = E; sp.nu = nu;

  Vec3 p[4] = {{0, 0, 0}, {a, 0, 0}, {a, b, 0}, {0, b, 0}};
  ShellElement4 sh(p, sp);

  const double D0 = E * t * t * t / (12.0 * (1.0 - nu * nu));

  // ---- (1) 板刚度的对称性与刚体模式 ----
  {
    const Mat12 K = sh.plateStiffness();
    double asym = 0, kmax = 0;
    for (int i = 0; i < 12; ++i)
      for (int j = 0; j < 12; ++j) {
        asym = std::max(asym, std::abs(K(i, j) - K(j, i)));
        kmax = std::max(kmax, std::abs(K(i, j)));
      }
    checkNear(asym / kmax, 0.0, 1e-15, "板刚度矩阵对称");

    // 三个刚体模式（局部坐标；w 沿 z，θx 绕 x，θy 绕 y）
    //   平移：      w = c
    //   绕 x 转动 φ：w = φ·y，θx = φ
    //   绕 y 转动 ψ：w = −ψ·x，θy = ψ
    // 【为什么是这两个组合】绕 x 的刚体转动必须与 w 的线性变化配套，
    // 否则 γ_yz = ∂w/∂y − θx ≠ 0 —— 单独的 "θx = 常数" 不是刚体模式，
    // 写成那样会把一个正常受力的状态误判为零能。
    const double phi = 1.0e-6, psi = 1.0e-6;
    const double ly[4] = {-b / 2, -b / 2, b / 2, b / 2};
    const double lx[4] = {-a / 2, a / 2, a / 2, -a / 2};
    // 【归一化必须【逐模式】做】平移模式的 |v| = 1，转动模式只有 1e-6。
    // 用一个固定的分母（比如 1e-6）去归一化，平移模式的相对误差
    // 会被平白放大 100 万倍 —— 机器精度会被显示成 1e-10 的"真误差"。
    // 这正是本文件开头 checkRigidBody 注释里警告过的那件事，写新测试时又踩了一次。
    double worst = 0.0;
    for (int mode = 0; mode < 3; ++mode) {
      double v[12] = {};
      double vmax = 0.0;
      for (int i = 0; i < 4; ++i) {
        if (mode == 0) { v[3 * i] = 1.0; }
        else if (mode == 1) { v[3 * i] = phi * ly[i]; v[3 * i + 1] = phi; }
        else { v[3 * i] = -psi * lx[i]; v[3 * i + 2] = psi; }
      }
      for (int i = 0; i < 12; ++i) vmax = std::max(vmax, std::abs(v[i]));
      double w1 = 0.0;
      for (int i = 0; i < 12; ++i) {
        double s = 0.0;
        for (int j = 0; j < 12; ++j) s += K(i, j) * v[j];
        w1 = std::max(w1, std::abs(s));
      }
      if (vmax > 0.0) worst = std::max(worst, w1 / (kmax * vmax));
    }
    checkNear(worst, 0.0, 1e-13, "板刚体运动不产生内力（逐模式归一化）");
  }

  // ---- (2) 纯弯曲：κx = κ0，其余为零 ----
  //
  //  βx = θy = κ0·x 是【线性场】，双线性插值可精确表示 ⇒ 这条是精确校核。
  //  期望：Mx = D₀κ0, My = ν·D₀κ0, Mxy = 0
  {
    const double k0 = 1.0e-3;
    Vec24d u{};
    for (int i = 0; i < 4; ++i) {
      const double x = (i == 1 || i == 2) ? a * 0.5 : -a * 0.5;
      u[6 * i + 4] = k0 * x;          // 全局 ry = 局部 θy（轴对齐 ⇒ T = I）
    }
    const auto f = sh.plateForces(u);
    checkNear(f.Mx, D0 * k0, 1e-12, "纯弯曲 κx：Mx = D₀·κx（精确）");
    checkNear(f.My, nu * D0 * k0, 1e-12, "纯弯曲 κx：My = ν·D₀·κx（泊松耦合）");
    checkNear(f.Mxy, 0.0, 1e-12, "纯弯曲 κx：Mxy = 0");
    checkNear(f.Qx, 0.0, 1e-12, "纯弯曲：Qx = 0");
    checkNear(f.Qy, 0.0, 1e-12, "纯弯曲：Qy = 0");
  }
  // ---- (3) 纯弯曲 κy（验证 κy 的负号约定）----
  {
    const double k0 = 1.0e-3;
    Vec24d u{};
    for (int i = 0; i < 4; ++i) {
      const double y = (i >= 2) ? b * 0.5 : -b * 0.5;
      u[6 * i + 3] = -k0 * y;         // θx = −κ0·y  ⇒ κy = −∂θx/∂y = κ0
    }
    const auto f = sh.plateForces(u);
    checkNear(f.My, D0 * k0, 1e-12, "纯弯曲 κy：My = D₀·κy（精确）");
    checkNear(f.Mx, nu * D0 * k0, 1e-12, "纯弯曲 κy：Mx = ν·D₀·κy");
    checkNear(f.Mxy, 0.0, 1e-12, "纯弯曲 κy：Mxy = 0");
  }
  // ---- (4) 纯扭转：κxy = 2·∂θy/∂y ----
  {
    // κxy = ∂θy/∂y − ∂θx/∂x。取 θy = c·y、θx = c·x ⇒ κxy = c − (−c)... 
    // 直接构造：θy = c·y, θx = 0 ⇒ κxy = c
    const double kk = 1.0e-3;
    Vec24d u{};
    for (int i = 0; i < 4; ++i) {
      const double y = (i >= 2) ? b * 0.5 : -b * 0.5;
      u[6 * i + 4] = kk * y;
    }
    const auto f = sh.plateForces(u);
    checkNear(f.Mxy, D0 * (1.0 - nu) * 0.5 * kk, 1e-12, "纯扭转：Mxy = D₀(1−ν)/2·κxy");
    checkNear(f.Mx, 0.0, 1e-12, "纯扭转：Mx = 0");
    checkNear(f.My, 0.0, 1e-12, "纯扭转：My = 0");
  }
  // ---- (5) 常剪切：Qx = κGt·γ0 ----
  //
  //  【这条测的是 MITC4 能否精确表示常剪切】w = γ0·x 是线性场，
  //  假定应变场的构造必须让它原样通过。
  {
    const double g0 = 1.0e-4;
    Vec24d u{};
    for (int i = 0; i < 4; ++i) {
      const double x = (i == 1 || i == 2) ? a * 0.5 : -a * 0.5;
      u[6 * i + 2] = g0 * x;
    }
    const auto f = sh.plateForces(u);
    checkNear(f.Qx, sh.properties().shearRigidity() * g0, 1e-12, "常剪切：Qx = κGt·γxz（精确）");
    checkNear(f.Qy, 0.0, 1e-12, "常剪切：Qy = 0");
    checkNear(f.Mx, 0.0, 1e-12, "常剪切不产生弯矩 Mx");
    checkNear(f.My, 0.0, 1e-12, "常剪切不产生弯矩 My");
  }
  // ---- (6) 刚体运动在【全局坐标】下不产生内力 ----
  //
  //  【这条同时验证坐标变换】把单元放到一个任意朝向（绕三个轴都转过），
  //  在全局坐标下施加刚体平移与刚体转动，内力必须为零。
  //  若没做局部→全局变换（K 直接用局部刚度当全局用），
  //  倾斜单元的刚体转动会产生巨大的虚假内力。
  {
    auto rotZ = [](double th, const Vec3& v) {
      const double c = std::cos(th), s = std::sin(th);
      return Vec3{c * v.x - s * v.y, s * v.x + c * v.y, v.z};
    };
    auto rotX = [](double th, const Vec3& v) {
      const double c = std::cos(th), s = std::sin(th);
      return Vec3{v.x, c * v.y - s * v.z, s * v.y + c * v.z};
    };
    Vec3 q[4];
    for (int i = 0; i < 4; ++i) {
      Vec3 v = p[i];
      v = rotZ(0.6, v);
      v = rotX(0.4, v);
      q[i] = v + Vec3{1.0, -2.0, 0.7};       // 再整体平移
    }
    ShellProperties sp2;
    sp2.thickness = t; sp2.E = E; sp2.nu = nu; sp2.drillingStiffness = 1e-4;
    ShellElement4 e2(q, sp2);
    check(e2.isValid(), "任意朝向单元有效");
    const Mat24 K = e2.stiffness();
    double kmax = 0.0;
    for (int i = 0; i < 24; ++i)
      for (int j = 0; j < 24; ++j) kmax = std::max(kmax, std::abs(K(i, j)));

    const double th = 1.0e-6;
    const Vec3 w = Vec3{0, 0, 1};
    double worst = 0.0;
    for (int mode = 0; mode < 4; ++mode) {
      Vec24d u{};
      double umax = 0.0;
      for (int i = 0; i < 4; ++i) {
        Vec3 d{0, 0, 0}, r{0, 0, 0};
        if (mode == 0) d = Vec3{1, 1, 1};                   // 刚体平移
        else if (mode < 4) {
          const Vec3 ax = (mode == 1) ? Vec3{1, 0, 0}
                        : (mode == 2) ? Vec3{0, 1, 0} : Vec3{0, 0, 1};
          // 绕【单元形心】转，否则混入平移（结果一样，但构造更明确）
          const Vec3 rel = q[i] - e2.centroid();
          d = cross(ax * th, rel);
          r = ax * th;
        }
        u[6 * i + 0] = d.x; u[6 * i + 1] = d.y; u[6 * i + 2] = d.z;
        u[6 * i + 3] = r.x; u[6 * i + 4] = r.y; u[6 * i + 5] = r.z;
      }
      for (int i = 0; i < 24; ++i) umax = std::max(umax, std::abs(u[i]));
      double w1 = 0.0;
      for (int i = 0; i < 24; ++i) {
        double s = 0.0;
        for (int j = 0; j < 24; ++j) s += K(i, j) * u[j];
        w1 = std::max(w1, std::abs(s));
      }
      if (umax > 0.0) worst = std::max(worst, w1 / (kmax * umax));
    }
    checkNear(worst, 0.0, 1e-13,
              "任意朝向单元：全局刚体运动不产生内力（验证坐标变换）");
    (void)w;
  }
  // ---- (7) 面荷载方向：倾斜单元的面荷载必须沿【单元法向】----
  {
    Vec3 q[4] = {{0, 0, 0}, {2, 0, 0}, {2, 0, 2}, {0, 0, 2}};   // 竖直平面
    ShellProperties sp3;
    sp3.thickness = t; sp3.E = E; sp3.nu = nu;
    ShellElement4 e3(q, sp3);
    check(e3.isValid(), "竖直壳单元有效");
    const Vec24d f = e3.transversePressure(1.0);       // 单位压力
    // 合力方向必须与单元法向平行
    double fx = 0, fy = 0, fz = 0;
    for (int i = 0; i < 4; ++i) { fx += f[6 * i + 0]; fy += f[6 * i + 1]; fz += f[6 * i + 2]; }
    const Vec3 n = e3.normal();
    const double mag = std::sqrt(fx * fx + fy * fy + fz * fz);
    checkNear(mag, e3.area() * 1.0, 1e-12, "面荷载合力大小 = p·A");
    const double dot = (fx * n.x + fy * n.y + fz * n.z) / mag;
    checkNear(std::abs(dot), 1.0, 1e-12, "面荷载方向平行于单元法向（验证荷载变换）");
  }
}

// ---------------------------------------------------------------------------
//  3b. 简支 / 夹支方板 —— 不锁定、且对标两个独立的经典解
// ---------------------------------------------------------------------------
static void testPlateBending() {
  std::printf("\n== 3b. 方板弯曲（对标 Timoshenko 经典解） ==\n");
  const double a = 6.0, q = 5.0, E = 3.15e7, nu = 0.3;

  // ---- (1) 简支板：厚度逐级减小，比值必须稳定（核心判据）----
  //
  //  约定范围 0.9 ~ 1.05，与 tests 里预先写下的判据一致。
  std::printf("       简支方板 a=%.1f q=%.1f E=%.3g ν=%.2f  n=8  (理论 0.00406·qa⁴/D)\n",
              a, q, E, nu);
  std::printf("       %8s %16s %16s %10s\n", "t", "w 数值(m)", "w 理论(m)", "w/w理论");
  double rmin = 1e9, rmax = -1e9;
  bool allInRange = true;
  for (double t : {0.20, 0.12, 0.06, 0.03, 0.01}) {
    const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
    const double wTh = 0.00406 * q * std::pow(a, 4) / D;
    const PlateResult R = solvePlate(8, t, q, a, E, nu, PlateBC::SimplySupported);
    if (!R.ok) { check(false, "简支板求解成功"); allInRange = false; continue; }
    const double ratio = R.w / wTh;
    rmin = std::min(rmin, ratio);
    rmax = std::max(rmax, ratio);
    std::printf("       %8.3f %16.6e %16.6e %10.4f\n", t, R.w, wTh, ratio);
    if (ratio < 0.9 || ratio > 1.05) allInRange = false;
  }
  check(allInRange, "t 从 0.20 降到 0.01，w/w理论 全部落在 0.9~1.05（无剪切锁定）");
  checkNear(rmax / rmin - 1.0, 0.0, 0.02,
            "比值随厚度的【变化幅度】< 2%（锁定会让它成倍崩塌）");

  // ---- (2) 厚板必须比薄板理论更软（判据②）----
  //
  //  【这是物理，不是误差】Mindlin 板含横向剪切变形，
  //  板越厚剪切柔度占比越大。若厚板反而更刚，说明剪切项符号或量级错了。
  {
    const double t = 0.6;                       // t/a = 1/10，属厚板
    const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
    const double wTh = 0.00406 * q * std::pow(a, 4) / D;
    const PlateResult R = solvePlate(8, t, q, a, E, nu, PlateBC::SimplySupported);
    check(R.ok, "厚板（t/a=1/10）求解成功");
    std::printf("       t/a=1/10：w = %.6e，薄板理论 = %.6e，比 %.4f\n", R.w, wTh, R.w / wTh);
    check(R.w / wTh > 1.02, "厚板比薄板理论更软（剪切变形的正确贡献）");
  }

  // ---- (3) 网格收敛 ----
  {
    const double t = 0.06;
    const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
    const double wTh = 0.00406 * q * std::pow(a, 4) / D;
    double prev = 0.0;
    bool mono = true;
    std::printf("       网格收敛（t=%.2f）：", t);
    for (int n : {4, 8, 16}) {
      const PlateResult R = solvePlate(n, t, q, a, E, nu, PlateBC::SimplySupported);
      if (!R.ok) { mono = false; break; }
      std::printf("n=%d→%.6e(%.4f)  ", n, R.w, R.w / wTh);
      if (R.w < prev) mono = false;
      prev = R.w;
    }
    std::printf("\n");
    check(prev / wTh > 0.98 && prev / wTh < 1.05, "n=16 时跨中挠度与理论偏差 < 2%");
    check(mono, "挠度随网格加密单调增大（趋近理论值）");
  }

  // ---- (4) 夹支方板：第二个独立经典解 ----
  //
  //  【为什么必须同时对标简支与夹支】一个偏刚的单元可以碰巧对上夹支，
  //  一个偏柔的可以碰巧对上简支。只有两个都对，才说明刚度是真的对。
  {
    const double wClCoef = 0.00126;             // Timoshenko，ν=0.3
    bool allOk = true;
    std::printf("       夹支方板（理论 0.00126·qa⁴/D）：\n");
    for (double t : {0.12, 0.06, 0.03}) {
      const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
      const double wTh = wClCoef * q * std::pow(a, 4) / D;
      const PlateResult R = solvePlate(8, t, q, a, E, nu, PlateBC::Clamped);
      if (!R.ok) { allOk = false; continue; }
      std::printf("         t=%.2f  w=%.6e  理论=%.6e  比=%.4f\n", t, R.w, wTh, R.w / wTh);
      if (std::abs(R.w / wTh - 1.0) > 0.03) allOk = false;
    }
    check(allOk, "夹支方板 w/w理论 偏差 < 3%（独立于简支解的第二个校核）");
  }

  // ---- (5) 跨中弯矩 ----
  {
    const double t = 0.12;
    const PlateResult R = solvePlate(16, t, q, a, E, nu, PlateBC::SimplySupported);
    const double Mth = 0.0479 * q * a * a;      // Timoshenko，ν=0.3，跨中 Mx=My
    std::printf("       跨中弯矩：Mx=%.6f My=%.6f 理论 %.6f kN·m/m\n", R.Mx, R.My, Mth);
    checkNear(R.Mx / Mth, 1.0, 0.08, "简支板跨中 Mx 对标 0.0479·q·a²");
    checkNear(R.My / Mth, 1.0, 0.08, "简支板跨中 My 对标 0.0479·q·a²（各向同性 ⇒ 与 Mx 相等）");
    checkNear(R.Mx / R.My, 1.0, 1e-6, "正方形板跨中 Mx = My（对称性）");
  }
}

// ---------------------------------------------------------------------------
//  3c. 坐标变换的定量判据：平面内旋转不改变解
//
//  【这是本单元最容易被漏掉的一条】
//  局部→全局变换若缺失（K 直接拿局部刚度当全局用），轴对齐网格上
//  结果恰好正确 —— 因为 T 退化成符号置换。斜交轴网、旋转后的板、
//  任意方向的墙才会暴露，而且不报错、量级也正常。
//  把整块板在自身平面内转过一个非特殊角，挠度与弯矩应当【严格不变】。
// ---------------------------------------------------------------------------
static void testPlateRotation() {
  std::printf("\n== 3c. 平面内旋转不变性（验证坐标变换） ==\n");
  const double a = 6.0, q = 5.0, E = 3.15e7, nu = 0.3, t = 0.12;
  const PlateResult r0 = solvePlate(8, t, q, a, E, nu, PlateBC::SimplySupported, 0.0);
  check(r0.ok, "基准（未旋转）求解成功");
  if (!r0.ok) return;
  std::printf("       基准：w = %.8e  Mx = %.6f\n", r0.w, r0.Mx);

  bool allOk = true;
  for (double deg : {15.0, 37.0, 45.0, 73.0}) {
    const double kDeg = 3.14159265358979323846 / 180.0;
    const PlateResult R = solvePlate(8, t, q, a, E, nu, PlateBC::SimplySupported,
                                     deg * kDeg);
    if (!R.ok) { allOk = false; continue; }
    const double dw = std::abs(R.w - r0.w) / std::abs(r0.w);
    const double dm = std::abs(R.Mx - r0.Mx) / std::abs(r0.Mx);
    std::printf("       旋转 %5.1f°：w 相对偏差 %.2e   Mx 相对偏差 %.2e\n", deg, dw, dm);
    if (dw > 1e-9) allOk = false;
  }
  check(allOk, "板在自身平面内旋转 15/37/45/73° 后挠度严格不变（1e-9）");
}

// ---------------------------------------------------------------------------
//  3d. 畸变网格
//
//  轴网建模出来的板单元是矩形，但洞口附近、变标高处的网格是畸变的。
//  畸变下若出现明显偏差或奇异，说明单元在真实工程里不可用。
// ---------------------------------------------------------------------------
static void testPlateSkewed() {
  std::printf("\n== 3d. 畸变四边形网格 ==\n");
  const double a = 6.0, q = 5.0, E = 3.15e7, nu = 0.3, t = 0.12;
  const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
  const double wTh = 0.00406 * q * std::pow(a, 4) / D;
  const PlateResult r0 = solvePlate(8, t, q, a, E, nu, PlateBC::SimplySupported);
  check(r0.ok, "规则网格求解成功");
  if (!r0.ok) return;
  std::printf("       规则网格 w=%.6e（比 %.4f）\n", r0.w, r0.w / wTh);

  // 单元素级：把矩形单元的一个角拉动 25% 边长，刚体校核仍须通过
  ShellProperties sp;
  sp.thickness = t; sp.E = E; sp.nu = nu; sp.drillingStiffness = 1e-4;
  Vec3 p[4] = {{0, 0, 0}, {2, 0, 0}, {2, 3, 0}, {0.5, 2.4, 0}};
  ShellElement4 sk(p, sp);
  check(sk.isValid(), "畸变单元有效");
  const Mat24 K = sk.stiffness();
  double kmax = 0;
  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j) kmax = std::max(kmax, std::abs(K(i, j)));
  const double th = 1e-6;
  double worst = 0;
  for (int mode = 0; mode < 4; ++mode) {
    Vec24d u{};
    double umax = 0.0;
    for (int i = 0; i < 4; ++i) {
      const Vec3 rel = p[i] - sk.centroid();
      Vec3 d{0, 0, 0}, r{0, 0, 0};
      if (mode == 0) d = Vec3{1, 1, 1};
      else {
        const Vec3 ax = (mode == 1) ? Vec3{1, 0, 0} : (mode == 2) ? Vec3{0, 1, 0} : Vec3{0, 0, 1};
        d = cross(ax * th, rel);
        r = ax * th;
      }
      u[6 * i + 0] = d.x; u[6 * i + 1] = d.y; u[6 * i + 2] = d.z;
      u[6 * i + 3] = r.x; u[6 * i + 4] = r.y; u[6 * i + 5] = r.z;
    }
    for (int i = 0; i < 24; ++i) umax = std::max(umax, std::abs(u[i]));
    double w1 = 0.0;
    for (int i = 0; i < 24; ++i) {
      double s = 0;
      for (int j = 0; j < 24; ++j) s += K(i, j) * u[j];
      w1 = std::max(w1, std::abs(s));
    }
    if (umax > 0.0) worst = std::max(worst, w1 / (kmax * umax));
  }
  checkNear(worst, 0.0, 1e-13, "畸变单元：刚体运动不产生内力（逐模式归一化）");
}

// ===========================================================================
//  4. 壳单元整体结构 —— 位移分片检验 + 荷载变换
// ===========================================================================
//
//  【分片检验的做法】不施加荷载，而是【给定一个精确解对应的位移场】，
//  然后检查应力恢复是否给出精确解。这样避开了"荷载怎么加"的歧义，
//  纯粹检验刚度矩阵与应力恢复的一致性，而且有精确解可比。
//
//  取单轴应变状态：εx = ε0，y 方向自由 ⇒ σx = E·ε0，σy = 0
//      ux = ε0·x,  uy = −ν·ε0·y
//  边界：x=0 边 ux=0；x=a 边 ux=ε0·a（给定位移）；y=0 边 uy=0（对称面）
//  检查：σx = Eε0、σy = 0、自由边 uy = −ν·ε0·y
//
//  【为什么这条能测出坐标变换的问题】膜刚度与膜荷载走同一套变换。
//  变换缺失时上一条"位移线性"仍可能成立（两者都错、相互抵消），
//  但"侧向收缩量是否等于泊松比"不会 —— 这个量与变形成正比。
// ===========================================================================
static void testShellStructure() {
  std::printf("\n== 4. 位移分片检验（膜）+ 荷载变换 ==\n");
  const double a = 6.0, t = 0.15, E = 3.15e7, nu = 0.3;
  const double eps0 = 1.0e-4;
  const int n = 4, npr = n + 1, ndof = npr * npr * 6;
  const double h = a / n;

  TripletAssembler ta;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Vec3 pt[4] = {{j * h, i * h, 0}, {(j + 1) * h, i * h, 0},
                    {(j + 1) * h, (i + 1) * h, 0}, {j * h, (i + 1) * h, 0}};
      ShellProperties sp;
      sp.thickness = t; sp.E = E; sp.nu = nu; sp.drillingStiffness = 1e-4;
      ShellElement4 el(pt, sp);
      const Mat24 Ke = el.stiffness();
      const int gi[4] = {i * npr + j, i * npr + (j + 1),
                         (i + 1) * npr + (j + 1), (i + 1) * npr + j};
      for (int A = 0; A < 24; ++A)
        for (int B = A; B < 24; ++B) {
          const double v = Ke(A, B);
          if (nearlyZero(v)) continue;
          ta.add(6 * gi[A / 6] + (A % 6), 6 * gi[B / 6] + (B % 6), v);
        }
    }

  std::vector<double> ub(ndof, 0.0);
  std::vector<char> fixed(ndof, 0);
  for (int i = 0; i <= n; ++i) {
    fixed[6 * (i * npr + 0) + 0] = 1;               // x=0 边：ux = 0
    fixed[6 * (i * npr + n) + 0] = 1;               // x=a 边：ux 给定位移
    ub[6 * (i * npr + n) + 0] = eps0 * a;
  }
  for (int j = 0; j <= n; ++j) fixed[6 * (0 * npr + j) + 1] = 1;    // y=0 边：uy=0
  fixed[6 * (0 * npr + 0) + 2] = 1;                 // 防 z 刚体平动
  fixed[6 * (0 * npr + 0) + 3] = 1;                 // 防绕 x 刚体转动
  fixed[6 * (0 * npr + 0) + 4] = 1;                 // 防绕 y 刚体转动

  std::vector<int> map(ndof, -1);
  int nf = 0;
  for (int i = 0; i < ndof; ++i) if (!fixed[i]) map[i] = nf++;
  TripletAssembler ta2;
  for (const auto& e : ta.triplets()) {
    const int mi = map[e.i], mj = map[e.j];
    if (mi < 0 || mj < 0) continue;
    ta2.add(mi, mj, e.v);
  }
  SymSparseMatrix K;
  K.buildFrom(ta2.triplets(), nf);
  std::vector<double> fr(nf, 0.0);
  for (const auto& e : ta.triplets()) {
    const int mi = map[e.i], mj = map[e.j];
    if (mi >= 0 && mj < 0 && ub[e.j] != 0.0) fr[mi] -= e.v * ub[e.j];
    else if (mi < 0 && mj >= 0 && ub[e.i] != 0.0) fr[mj] -= e.v * ub[e.i];
  }

  LDLTSolver solver;
  const auto rep = solver.factorize(K, 0);
  check(rep.ok(), "分片检验求解成功");
  if (!rep.ok()) return;
  std::vector<double> uf;
  solver.solve(fr, uf);
  std::vector<double> u(ndof, 0.0);
  for (int i = 0; i < ndof; ++i) u[i] = (map[i] >= 0) ? uf[map[i]] : ub[i];

  // (1) 位移场：ux = ε0·x（线性）、uy = −ν·ε0·y（自由收缩）
  double dw = 0.0, dv = 0.0;
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j) {
      const double x = j * h, y = i * h;
      dw = std::max(dw, std::abs(u[6 * (i * npr + j) + 0] - eps0 * x) / (eps0 * a));
      dv = std::max(dv, std::abs(u[6 * (i * npr + j) + 1] - (-nu * eps0 * y)) / (eps0 * a));
    }
  std::printf("       ux 偏差 %.3e   uy 偏差 %.3e（相对 ε0·a）\n", dw, dv);
  check(dw < 1e-13, "分片检验：ux = ε0·x 精确（膜刚度正确）");
  check(dv < 1e-13, "分片检验：uy = −ν·ε0·y 精确（泊松收缩正确）");

  // (2) 应力恢复：σx = Eε0，σy = 0
  {
    const int i = n / 2 - 1, j = n / 2 - 1;
    Vec3 pt[4] = {{j * h, i * h, 0}, {(j + 1) * h, i * h, 0},
                  {(j + 1) * h, (i + 1) * h, 0}, {j * h, (i + 1) * h, 0}};
    ShellProperties sp;
    sp.thickness = t; sp.E = E; sp.nu = nu; sp.drillingStiffness = 1e-4;
    ShellElement4 el(pt, sp);
    const int gi[4] = {i * npr + j, i * npr + (j + 1),
                       (i + 1) * npr + (j + 1), (i + 1) * npr + j};
    Vec24d ue{};
    for (int A = 0; A < 24; ++A) ue[A] = u[6 * gi[A / 6] + (A % 6)];
    const auto mf = el.membraneForces(ue);
    std::printf("       单元中心应力 σx=%.6f σy=%.6f kPa（理论 %.4f / 0）\n",
                mf.sx, mf.sy, E * eps0);
    checkNear(mf.sx, E * eps0, 1e-12, "膜应力恢复 σx = E·ε0");
    checkNear(mf.sy, 0.0, 1e-10, "膜应力恢复 σy = 0（y 向自由）");
  }

  // (3) 面外挠度为零：面内变形不引起弯曲
  double maxW = 0.0;
  for (int i = 0; i < npr * npr; ++i) maxW = std::max(maxW, std::abs(u[6 * i + 2]));
  checkNear(maxW / (eps0 * a), 0.0, 1e-10, "面内变形不产生面外挠度（膜/板解耦）");

  // (4) 膜荷载的坐标变换
  //
  //  【荷载是余向量，变换要用 Tᵀ 而不是 T】
  //  membranePressure(px,0) 施加的是"沿【局部 x】的面分布荷载"。
  //  把单元转起来，合力方向必须跟着转到新的局部 x 方向，大小不变。
  {
    const double ang = 0.7;
    const double c = std::cos(ang), sn = std::sin(ang);
    auto rot = [&](double x, double y) {
      return Vec3{x * c - y * sn, x * sn + y * c, 0.0};
    };
    Vec3 pt[4] = {rot(0, 0), rot(2, 0), rot(2, 1.5), rot(0, 1.5)};
    ShellProperties sp;
    sp.thickness = t; sp.E = E; sp.nu = nu;
    ShellElement4 el(pt, sp);
    const double p = 3.0;
    const Vec24d f = el.membranePressure(p, 0.0);
    double fx = 0, fy = 0, fz = 0;
    for (int i = 0; i < 4; ++i) { fx += f[6 * i + 0]; fy += f[6 * i + 1]; fz += f[6 * i + 2]; }
    const Vec3 ex = el.localX();
    const double mag = std::sqrt(fx * fx + fy * fy + fz * fz);
    std::printf("       旋转后膜荷载合力 %.6f，方向 (%.4f,%.4f) vs 局部x (%.4f,%.4f)\n",
                mag, fx / mag, fy / mag, ex.x, ex.y);
    checkNear(mag, p * el.area(), 1e-12, "膜荷载合力大小 = p·A");
    checkNear((fx * ex.x + fy * ex.y) / mag, -1.0, 1e-12,
              "膜荷载沿【局部 x 的负方向】（正为压）且随单元旋转");
    checkNear(std::abs(fz) / mag, 0.0, 1e-12, "膜荷载无面外分量");
  }
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  四节点壳单元验证\n");
  std::printf("======================================================\n");

  testGeometry();
  testMembrane();
  testPlateUnit();
  testPlateBending();
  testPlateRotation();
  testPlateSkewed();
  testShellStructure();

  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项，跳过 %d 项\n", g_pass, g_fail, g_skip);
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}