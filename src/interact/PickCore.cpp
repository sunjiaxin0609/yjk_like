// =============================================================================
//  yjk/interact/PickCore.cpp  ——  射线拾取几何实现
// =============================================================================
#include "yjk/interact/PickCore.h"

#include <algorithm>
#include <cmath>

namespace yjk {
namespace interact {

namespace {

// 射线与三角形求交（Möller–Trumbore 的单三角形版本，供四边形拆半使用）。
// 命中判据：t > 0，且重心坐标在 [0,1]。
bool triHit(const Ray3& ray, const Vec3& v0, const Vec3& v1, const Vec3& v2,
            double* tOut) {
  const Vec3 e1 = v1 - v0;
  const Vec3 e2 = v2 - v0;
  const Vec3 p = cross(ray.d, e2);
  const double det = dot(e1, p);
  // 平行（或退化三角形）时无交点；同时要求 det 符号稳定，不搞"双重命中"
  if (std::abs(det) < 1e-14) return false;
  const double inv = 1.0 / det;
  const Vec3 s = ray.o - v0;
  const double u = dot(s, p) * inv;
  if (u < 0.0 || u > 1.0) return false;
  const Vec3 q = cross(s, e1);
  const double v = dot(ray.d, q) * inv;
  if (v < 0.0 || u + v > 1.0) return false;
  const double t = dot(e2, q) * inv;
  if (t <= 0.0) return false;
  if (tOut) *tOut = t;
  return true;
}

}  // namespace

double pointSegDist(const Vec3& p, const Vec3& a, const Vec3& b) {
  const Vec3 ab = b - a;
  const double len2 = norm2(ab);
  if (len2 < 1e-30) return norm(p - a);
  const double s = dot(p - a, ab) / len2;
  const double t = std::clamp(s, 0.0, 1.0);
  return norm(p - (a + t * ab));
}

bool raySegDist(const Ray3& ray, const Vec3& a, const Vec3& b,
                double* sOut, double* tOut) {
  const Vec3 u = b - a;         // 线段方向
  const Vec3 v = ray.d;         // 射线方向（单位）
  const Vec3 w0 = a - ray.o;

  const double aa = dot(u, u);
  const double bb = dot(u, v);
  const double cc = dot(v, v);
  const double dd = dot(u, w0);
  const double ee = dot(v, w0);
  const double denom = aa * cc - bb * bb;

  double s = 0.0;
  if (std::abs(denom) > 1e-14) {
    s = (bb * ee - cc * dd) / denom;
  }
  s = std::clamp(s, 0.0, 1.0);
  // 固定 s 后在射线上求最近点参数（简化公式）
  const double t = dot(w0 + s * u, v) / cc;
  if (sOut) *sOut = s;
  if (tOut) *tOut = t;
  return true;
}

bool rayTriHit(const Ray3& ray, const Vec3& v0, const Vec3& v1, const Vec3& v2,
               double* tOut) {
  return triHit(ray, v0, v1, v2, tOut);
}

PickHit pickNearest(const Ray3& rayIn,
                    const std::vector<PickNode>& nodes,
                    const std::vector<PickBeam>& beams,
                    const std::vector<PickShell>& shells,
                    double nodeR, double beamR, double tMax) {
  Ray3 ray = rayIn;
  const double dl = norm(ray.d);
  if (dl < 1e-30) return PickHit{};                 // 退化射线：无命中
  if (std::abs(dl - 1.0) > 1e-9) ray.d = ray.d * (1.0 / dl);  // 归一化兜底

  const bool clampT = tMax > 0.0;
  PickHit best;
  best.t = clampT ? tMax : 1e300;

  // ---- 节点：射线到球心的最短距离 ----
  // 射线最近点参数 tN = dot(p-o, d)，距离 = |p - (o + tN·d)|
  if (nodeR > 0.0) {
    const double r2 = nodeR * nodeR;
    for (const PickNode& n : nodes) {
      const double t = dot(n.p - ray.o, ray.d);
      if (t < 0.0 || t >= best.t) continue;
      const Vec3 foot = ray.o + t * ray.d;
      const double d2 = norm2(n.p - foot);
      if (d2 <= r2) {
        best.kind = PickKind::Node;
        best.id = n.id;
        best.t = t;
        best.dist = std::sqrt(d2);
        break;    // 同一节点只需命中一次
      }
    }
  }

  // ---- 梁：射线到线段的最近距离（容差圆柱）----
  if (beamR > 0.0) {
    const double r2 = beamR * beamR;
    for (const PickBeam& bm : beams) {
      double s = 0.0, t = 0.0;
      raySegDist(ray, bm.a, bm.b, &s, &t);
      if (t < 0.0 || t >= best.t) continue;
      const Vec3 segPt = bm.a + s * (bm.b - bm.a);
      const Vec3 rayPt = ray.o + t * ray.d;
      const double d2 = norm2(segPt - rayPt);
      if (d2 <= r2) {
        best.kind = PickKind::Beam;
        best.id = bm.id;
        best.t = t;
        best.dist = std::sqrt(d2);
      }
    }
  }

  // ---- 壳：四边形拆两个三角形，最近交点 ----
  for (const PickShell& sh : shells) {
    double t = 0.0;
    if (triHit(ray, sh.p[0], sh.p[1], sh.p[2], &t) && t > 0.0 && t < best.t) {
      best.kind = PickKind::Shell;
      best.id = sh.id;
      best.t = t;
      best.dist = 0.0;
    }
    if (triHit(ray, sh.p[0], sh.p[2], sh.p[3], &t) && t > 0.0 && t < best.t) {
      best.kind = PickKind::Shell;
      best.id = sh.id;
      best.t = t;
      best.dist = 0.0;
    }
  }

  return best;
}

}  // namespace interact
}  // namespace yjk