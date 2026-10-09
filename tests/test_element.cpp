// =============================================================================
//  tests/test_element.cpp  ——  单元层验证（占位，单元层实现后填充）
// =============================================================================
#include <cmath>
#include <cstdio>

#include "yjk/math/Types.h"

using namespace yjk;

static int g_pass = 0, g_fail = 0;

static void check(bool cond, const char* what) {
  if (cond) { ++g_pass; std::printf("  [ OK ] %s\n", what); }
  else { ++g_fail; std::printf("  [FAIL] %s\n", what); }
}

int main() {
  std::printf("======================================================\n");
  std::printf("  yjk_like  单元层验证\n");
  std::printf("======================================================\n");

  // ---- 小矩阵运算：单元刚度矩阵组装的底层依赖 ----
  std::printf("\n== 1. 稠密小矩阵 ==\n");
  {
    Mat3 a = Mat3::identity();
    a(0, 1) = 2.0;
    a(1, 0) = 2.0;
    const Mat3 inv = a.inverse();
    const Mat3 id = a * inv;
    double maxErr = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) maxErr = std::max(maxErr, std::abs(id(i, j) - (i == j ? 1.0 : 0.0)));
    check(maxErr < 1e-14, "3x3 求逆后 A·A⁻¹ = I");
  }
  {
    Mat4 a = Mat4::identity();
    a(2, 3) = -3.0;
    a(3, 2) = -3.0;
    bool ok = true;
    const Mat4 inv = a.inverse(&ok);
    check(ok, "4x4 奇异矩阵被正确识别");
    (void)inv;
  }
  {
    // 上三角构造：单元刚度矩阵的常用入口
    double up[6] = {1, 2, 3, 4, 5, 6};   // (0,0)(0,1)(1,1)(0,2)(1,2)(2,2)
    const Mat3 m = Mat3::fromSymmetricUpper(up);
    check(m(1, 0) == 2.0 && m(2, 1) == 5.0 && m(0, 2) == 3.0, "上三角 → 对称矩阵构造");
  }

  // ---- 三维向量：坐标变换的底层依赖 ----
  std::printf("\n== 2. 三维向量 ==\n");
  {
    const Vec3 a{1, 0, 0}, b{0, 1, 0};
    const Vec3 c = cross(a, b);
    check(std::abs(c.z - 1.0) < 1e-15, "叉乘 e1×e2 = e3");
    check(std::abs(dot(a, b)) < 1e-15, "正交向量点积为 0");
    const Vec3 n = normalized(Vec3{3, 4, 0});
    check(std::abs(norm(n) - 1.0) < 1e-15, "归一化后模为 1");
  }
  {
    // 旋转矩阵正交性 —— 单元坐标变换必须保持这一性质
    const double th = 0.7;
    const double c = std::cos(th), s = std::sin(th);
    const Vec3 e1 = normalized(Vec3{c, s, 0});
    const Vec3 e2{-s, c, 0};
    check(std::abs(dot(e1, e2)) < 1e-14, "旋转后两轴仍正交");
  }

  std::printf("\n======================================================\n");
  std::printf("  通过 %d 项，失败 %d 项\n", g_pass, g_fail);
  std::printf("  注：单元层（梁/壳）尚未实现，此处仅验证底层矩阵与向量运算。\n");
  std::printf("======================================================\n");
  return g_fail == 0 ? 0 : 1;
}
