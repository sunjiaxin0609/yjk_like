// -----------------------------------------------------------------------------
//  SpringElement3D —— 三维弹簧单元（两节点，6 方向独立刚度）
//
//  用途：
//    · 弹簧支座 / 弹性节点连接（k=0 释放该方向、大 k 近似刚接）
//    · 与 BeamElement3D 的端部释放（铰）配合，实现"部分刚接、部分弹性"
//
//  力学模型：
//    局部坐标系 x' 沿 i→j 方向，y'/z' 由 up 参考向量构造（与梁单元同约定）。
//    每个局部自由度 d∈{0..5}（ux,uy,uz,rx,ry,rz）是一根独立的两节点弹簧：
//        K_local(d,d)   =  k_d
//        K_local(d,6+d) = -k_d       （j 端）
//        K_local(6+d,d) = -k_d
//        K_local(6+d,6+d)= k_d
//    即 力 = k_d·(u_j − u_i)。6 个方向互不耦合。
//    k_d = 0 → 该方向【释放】（单元在该方向不提供任何刚度），
//    DofNumbering 的零刚度锁定会自动接管（若该自由度无其他单元提供刚度）。
//
//  全局刚度：K_global = Tᵀ · K_local · T，T = blkdiag(R,R,R,R)，
//    R 为局部→全局旋转矩阵（u_local = R·u_global，与 BeamElement3D 一致）。
//
//  与梁单元的区别：弹簧没有长度衰减（k 不除以 L），也没有弯曲-平移耦合；
//  它就是在两个节点的对应自由度之间直接拉一根"数字弹簧"。
// -----------------------------------------------------------------------------
#pragma once
#include <array>
#include <cmath>
#include <string>
#include "yjk/math/Types.h"

namespace yjk {

class SpringElement3D {
 public:
  SpringElement3D() = default;

  // k[0..5]：局部 ux,uy,uz,rx,ry,rz 六个方向的刚度（力学单位）：
  //   平动刚度 kN/m（或 kips/in 等一致单位制），转动刚度 kN·m/rad
  bool init(const Vec3& pi, const Vec3& pj, const std::array<double, 6>& k,
            const Vec3& up = Vec3{0, 0, -1}) {
    pi_ = pi;
    pj_ = pj;
    k_ = k;
    valid_ = false;
    error_.clear();

    const Vec3 dx = pj - pi;
    length_ = norm(dx);
    if (length_ < 1.0e-12) {
      error_ = "弹簧单元两端节点重合（长度为零）";
      return false;
    }
    ex_ = dx * (1.0 / length_);

    // 与 BeamElement3D 相同的局部坐标系构造：up 投影到垂直 x' 的平面作为 y'，
    // z' = x' × y'。up 与 x' 平行时退化 → 回退到与全局 Y/X 垂直的方向。
    Vec3 proj = up - ex_ * dot(up, ex_);
    if (norm(proj) < 1.0e-12) {
      const Vec3 fallback = (std::abs(ex_.y) < 0.9) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
      proj = fallback - ex_ * dot(fallback, ex_);
    }
    ey_ = normalized(proj);
    ez_ = cross(ex_, ey_);
    if (norm(ez_) < 1.0e-12) {
      error_ = "弹簧单元局部坐标系构造失败（x' 与 y' 共线）";
      return false;
    }
    ez_ = normalized(ez_);

    for (double kd : k_)
      if (kd < 0.0) {
        error_ = "弹簧刚度不能为负";
        return false;
      }
    valid_ = true;
    return true;
  }

  bool isValid() const { return valid_; }
  const std::string& error() const { return error_; }
  double length() const { return length_; }
  const std::array<double, 6>& stiffness() const { return k_; }
  const Vec3& localX() const { return ex_; }
  const Vec3& localY() const { return ey_; }
  const Vec3& localZ() const { return ez_; }

  // k_d ≈ 0 → 该方向释放（单元不提供刚度）
  bool isReleasedComp(int comp) const { return nearlyZero(k_[static_cast<size_t>(comp)]); }

  // -------------------------------------------------------------------------
  //  局部刚度（12×12）：6 个方向各一根两节点弹簧，互不耦合
  // -------------------------------------------------------------------------
  Mat12 localStiffness() const {
    Mat12 K = Mat12::zero();
    if (!valid_) return K;
    for (int d = 0; d < 6; ++d) {
      const double kd = k_[static_cast<size_t>(d)];
      K(d, d) = kd;
      K(d, 6 + d) = -kd;
      K(6 + d, d) = -kd;
      K(6 + d, 6 + d) = kd;
    }
    return K;
  }

  // -------------------------------------------------------------------------
  //  全局刚度：K_global = Tᵀ K_local T，T = blkdiag(R,R,R,R)
  // -------------------------------------------------------------------------
  Mat12 stiffnessGlobal() const {
    Mat12 T = Mat12::zero();
    const double R[3][3] = {
        {ex_.x, ex_.y, ex_.z},
        {ey_.x, ey_.y, ey_.z},
        {ez_.x, ez_.y, ez_.z},
    };
    for (int blk = 0; blk < 4; ++blk)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) T(blk * 3 + i, blk * 3 + j) = R[i][j];

    const Mat12 kl = localStiffness();
    const Mat12 kt = kl * T;          // K_local · T
    Mat12 Kg = Mat12::zero();
    for (int r = 0; r < 12; ++r)
      for (int c = 0; c < 12; ++c) {
        double s = 0.0;
        for (int t = 0; t < 12; ++t) s += T(t, r) * kt(t, c);   // Tᵀ · (K_local·T)
        Kg(r, c) = s;
      }
    return Kg;
  }

  // -------------------------------------------------------------------------
  //  全局 → 局部（u_local = R·u_global，逐 3 自由度块）
  // -------------------------------------------------------------------------
  Vec12d toLocalVector(const Vec12d& ug) const {
    Vec12d ul{};
    const double R[3][3] = {
        {ex_.x, ex_.y, ex_.z},
        {ey_.x, ey_.y, ey_.z},
        {ez_.x, ez_.y, ez_.z},
    };
    for (int blk = 0; blk < 4; ++blk)
      for (int i = 0; i < 3; ++i) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j) s += R[i][j] * ug[blk * 3 + j];
        ul[blk * 3 + i] = s;
      }
    return ul;
  }

  // -------------------------------------------------------------------------
  //  弹簧内力恢复（i 端局部力）
  //   F_d = k_d·(u_j − u_i)，d = 0..5；j 端内力 = i 端取反。
  //   f_i 的符号约定与梁单元 EndForces 对齐：对 i 端分析体，正 = 拉力。
  // -------------------------------------------------------------------------
  //  uElGlobal：全局位移（12 长度），返回局部 i 端 6 内力
  std::array<double, 6> endForcesLocal(const Vec12d& uElGlobal) const {
    std::array<double, 6> f{};
    if (!valid_) return f;
    const Vec12d ul = toLocalVector(uElGlobal);
    for (int d = 0; d < 6; ++d) {
      const double kd = k_[static_cast<size_t>(d)];
      f[static_cast<size_t>(d)] = kd * (ul[6 + d] - ul[d]);   // k·(u_j − u_i)
    }
    return f;
  }

 private:
  Vec3 pi_{}, pj_{};
  Vec3 ex_{}, ey_{}, ez_{};
  std::array<double, 6> k_{};
  double length_ = 0.0;
  bool valid_ = false;
  std::string error_;
};

}  // namespace yjk