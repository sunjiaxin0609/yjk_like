// =============================================================================
//  yjk/math/Types.h  ——  全局基础类型
//  单位约定（和所有有限元教科书一致，内部统一用 kN-m-N，输出时再换算）：
//    长度  m        力    kN       应力  kPa ( = kN/m^2 )
//    线荷载 kN/m    面荷载 kN/m^2   集中力 kN      弯矩 kN·m
//  截面模量等几何量以 m 为单位输入，内部 E = 3.0e7 kPa (C30)，fck=20.1 MPa=20100 kPa
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <limits>

namespace yjk {

using Id = int;

// ---- 精度与容差 -----------------------------------------------------------
// 结构分析里"接近奇异"的判断要比数值分析严格得多：
// 刚度矩阵的条件数随构件刚度比急剧恶化（板/梁刚度比可达 1e4~1e5），
// 所以奇异阈值取相对量级，并且要求 DOF 有明确约束。
constexpr double kEps = 1.0e-14;
constexpr double kPivTol = 1.0e-30;   // 主元绝对下限，低于此认为奇异

inline bool nearlyZero(double v, double tol = 1.0e-12) { return std::abs(v) <= tol; }

// ---- 三维向量 -------------------------------------------------------------
struct Vec3 {
  double x{0.0}, y{0.0}, z{0.0};

  constexpr Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
  Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

  // 按分量下标访问（0=x, 1=y, 2=z）
  //
  // 【为什么需要】自由度编号、荷载组装、后处理全是"按分量循环"的写法：
  //     for (int k = 0; k < 3; ++k) f[dof[k]] += force[k];
  // 若没有 []，每个这样的循环都要写三遍 if，或者用 lambda 取分量 ——
  // 而后者反而更容易出错（梁单元层就踩过"哪个分量对应哪个自由度"的坑）。
  // 统一用 [] 后，"第 k 个平动分量"在代码里永远是 force[k]，
  // 读代码时能直接对应自由度编号的 k，不会串位。
  double& operator[](int i) { return (i == 0) ? x : ((i == 1) ? y : z); }
  double operator[](int i) const { return (i == 0) ? x : ((i == 1) ? y : z); }
  static constexpr int size() { return 3; }
};

inline Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
inline Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
inline Vec3 operator*(Vec3 a, double s) { return a *= s; }
inline Vec3 operator*(double s, Vec3 a) { return a *= s; }
inline Vec3 operator-(const Vec3& a) { return Vec3{-a.x, -a.y, -a.z}; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double norm2(const Vec3& a) { return dot(a, a); }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

inline Vec3 normalized(const Vec3& a) {
  const double n = norm(a);
  return nearlyZero(n) ? Vec3{} : a * (1.0 / n);
}

// ---- 二维向量（平面内单元/2D框架用）--------------------------------------
struct Vec2 {
  double x{0.0}, y{0.0};
  constexpr Vec2() = default;
  constexpr Vec2(double x_, double y_) : x(x_), y(y_) {}
  Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }
  Vec2& operator-=(const Vec2& o) { x -= o.x; y -= o.y; return *this; }
};
inline Vec2 operator+(Vec2 a, const Vec2& b) { return a += b; }
inline Vec2 operator-(Vec2 a, const Vec2& b) { return a -= b; }
inline double dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
inline double cross(const Vec2& a, const Vec2& b) { return a.x * b.y - a.y * b.x; }
inline double norm(const Vec2& a) { return std::sqrt(a.x * a.x + a.y * a.y); }

// ---- 4x4 / 6x6 / 12x12 小矩阵 --------------------------------------------
// 单元层反复用到 12x12 以下的稠密小矩阵运算，自己写比引入外部 BLAS 快得多
// （单元刚度矩阵是热点，内存局部性极好，编译器能向量化）。

template <int N>
struct Mat {
  std::array<double, N * N> a{};   // 行主序 row-major

  double& operator()(int r, int c) { return a[r * N + c]; }
  double operator()(int r, int c) const { return a[r * N + c]; }

  static Mat zero() { return Mat{}; }
  static Mat identity() {
    Mat m;
    for (int i = 0; i < N; ++i) m(i, i) = 1.0;
    return m;
  }
  // 对称/反对称矩阵构造
  static Mat fromSymmetricUpper(const double* up) {
    Mat m;
    int p = 0;
    for (int r = 0; r < N; ++r)
      for (int c = r; c < N; ++c, ++p) { m(r, c) = up[p]; m(c, r) = up[p]; }
    return m;
  }

  Mat& operator+=(const Mat& o) {
    for (int i = 0; i < N * N; ++i) a[i] += o.a[i];
    return *this;
  }
  Mat& operator-=(const Mat& o) {
    for (int i = 0; i < N * N; ++i) a[i] -= o.a[i];
    return *this;
  }
  Mat& operator*=(double s) {
    for (int i = 0; i < N * N; ++i) a[i] *= s;
    return *this;
  }
  // C = this * o
  // 注：不提供 `friend operator*(const Mat&, const Mat&)` ——
  // 它与成员版本在 `a * b` 时会产生二义性（ISO C++ 明确规定为 ill-formed，
  // 尽管 GCC 会选更优的那个）。成员版本已足够。
  Mat operator*(const Mat& o) const {
    Mat c;
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        double s = 0.0;
        for (int k = 0; k < N; ++k) s += (*this)(i, k) * o(k, j);
        c(i, j) = s;
      }
    }
    return c;
  }
  // 注：不提供 `friend operator*(const Mat&, const Mat&)` ——
  // 它与成员版本在 `a * b` 时会产生二义性。
  Mat transposed() const {
    Mat t;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) t(i, j) = (*this)(j, i);
    return t;
  }

  // 逆矩阵（高斯-约当，带列主元）。N<=12 时最稳。
  // 对刚度矩阵而言失败通常意味着单元退化或约束不足，调用方应检查 ok。
  Mat inverse(bool* ok = nullptr) const {
    Mat m = *this;
    Mat inv = Mat::identity();
    for (int col = 0; col < N; ++col) {
      // 选主元
      int piv = col;
      double best = std::abs(m(col, col));
      for (int r = col + 1; r < N; ++r) {
        const double v = std::abs(m(r, col));
        if (v > best) { best = v; piv = r; }
      }
      if (best < kPivTol) { if (ok) *ok = false; return Mat{}; }
      if (piv != col) {
        for (int j = 0; j < N; ++j) { std::swap(m(col, j), m(piv, j)); std::swap(inv(col, j), inv(piv, j)); }
      }
      const double d = 1.0 / m(col, col);
      for (int j = 0; j < N; ++j) { m(col, j) *= d; inv(col, j) *= d; }
      for (int r = 0; r < N; ++r) {
        if (r == col) continue;
        const double f = m(r, col);
        if (nearlyZero(f)) continue;
        for (int j = 0; j < N; ++j) { m(r, j) -= f * m(col, j); inv(r, j) -= f * inv(col, j); }
      }
    }
    if (ok) *ok = true;
    return inv;
  }

  // y = this * v
  std::array<double, N> operator*(const std::array<double, N>& v) const {
    std::array<double, N> y{};
    for (int i = 0; i < N; ++i) {
      double s = 0.0;
      for (int j = 0; j < N; ++j) s += (*this)(i, j) * v[j];
      y[i] = s;
    }
    return y;
  }

  double maxAbs() const {
    double m = 0.0;
    for (double v : a) m = std::max(m, std::abs(v));
    return m;
  }
};

// 常用别名
using Mat2 = Mat<2>;
using Mat3 = Mat<3>;
using Mat4 = Mat<4>;
using Mat6 = Mat<6>;
using Mat8 = Mat<8>;     // 膜单元刚度（4 节点 × 2 面内自由度）
using Mat12 = Mat<12>;   // 板弯曲刚度（4 节点 × 3 面外自由度）
using Mat16 = Mat<16>;   // MITC4 剪切 B 矩阵（8 绑定点 × 2 分量）
using Mat24 = Mat<24>;   // 壳单元总刚度（4 节点 × 6 自由度）
using Vec3d = std::array<double, 3>;
using Vec6d = std::array<double, 6>;
using Vec12d = std::array<double, 12>;
using Vec24d = std::array<double, 24>;

}  // namespace yjk
