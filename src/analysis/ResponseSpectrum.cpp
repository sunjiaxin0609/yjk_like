// =============================================================================
//  src/analysis/ResponseSpectrum.cpp  ——  设计反应谱与振型分解反应谱法实现
// =============================================================================
#include <algorithm>
#include <cmath>

#include "yjk/analysis/ResponseSpectrum.h"

namespace yjk {

namespace {
constexpr double kG = 9.81;   // m/s² —— 参与质量(t)·g 才是地震力(kN)
}  // namespace

// -----------------------------------------------------------------------------
//  GB 50011-2010 地震影响系数曲线（规范 5.1.5）
// -----------------------------------------------------------------------------
double DesignSpectrum::gammaExp() const {
  // γ = 0.9 + (0.05 − ζ) / (0.5 + 5ζ)，ζ 为阻尼比
  return 0.9 + (0.05 - zeta) / (0.5 + 5.0 * zeta);
}

double DesignSpectrum::eta1() const {
  // η1 = 0.02 + (0.05 − ζ) / (4 + 32ζ)；小于 0 时取 0
  const double e = 0.02 + (0.05 - zeta) / (4.0 + 32.0 * zeta);
  return (e < 0.0) ? 0.0 : e;
}

double DesignSpectrum::eta2() const {
  // η2 = 1 + (0.05 − ζ) / (0.08 + 1.6ζ)；小于 0.55 时取 0.55
  const double e = 1.0 + (0.05 - zeta) / (0.08 + 1.6 * zeta);
  return (e < 0.55) ? 0.55 : e;
}

double DesignSpectrum::alpha(double T) const {
  if (!(T > 0.0)) return alphaMax;   // 零/负周期防御：按平台段最大值
  const double e2 = eta2();
  const double g  = gammaExp();
  const double e1 = eta1();
  const double t  = std::min(T, 6.0);   // 曲线终点：规范定义至 6s

  if (t <= 0.1) {
    // 上升段：α = [0.45 + 10(η2 − 0.45)T]·αmax
    return (0.45 + 10.0 * (e2 - 0.45) * t) * alphaMax;
  }
  if (t <= tg) {
    // 平台段
    return e2 * alphaMax;
  }
  if (t <= 5.0 * tg) {
    // 曲线下降段：α = (Tg/T)^γ · η2·αmax
    return std::pow(tg / t, g) * e2 * alphaMax;
  }
  // 直线下降段：α = [η2·0.2^γ − η1(T − 5Tg)]·αmax
  const double a = e2 * std::pow(0.2, g) - e1 * (t - 5.0 * tg);
  if (a < 0.0) return 0.0;
  return a * alphaMax;
}

// -----------------------------------------------------------------------------
//  振型分解反应谱法
// -----------------------------------------------------------------------------
double ResponseSpectrum::cqcRho(double zeta, double omegaJ, double omegaK) {
  // ρ_jk = 8ζ²(1+β)β^1.5 / ( (1−β²)² + 4ζ²β(1+β)² )，β = ω_k / ω_j
  // （规范式 5.2.2-3）。β = 1（重频）时 ρ = 1。
  const double b = omegaK / omegaJ;
  if (b <= 0.0 || omegaJ <= 0.0) return 0.0;
  const double z2 = zeta * zeta;
  const double num = 8.0 * z2 * (1.0 + b) * std::pow(b, 1.5);
  const double den = (1.0 - b * b) * (1.0 - b * b) +
                     4.0 * z2 * b * (1.0 + b) * (1.0 + b);
  if (den == 0.0) return 1.0;
  return num / den;
}

double ResponseSpectrum::srss(const std::vector<double>& v) {
  double s = 0.0;
  for (double x : v) s += x * x;
  return std::sqrt(s);
}

double ResponseSpectrum::cqc(double zeta, const std::vector<double>& omega,
                             const std::vector<double>& v) {
  const size_t n = std::min(omega.size(), v.size());
  double s = 0.0;
  for (size_t j = 0; j < n; ++j) {
    for (size_t k = 0; k < n; ++k) {
      s += cqcRho(zeta, omega[j], omega[k]) * v[j] * v[k];
    }
  }
  return std::sqrt(std::max(s, 0.0));
}

SpectrumResult ResponseSpectrum::solve(const ModalResult& mr,
                                       const DesignSpectrum& spec,
                                       const Options& opt) {
  SpectrumResult res;
  if (!mr.ok) {
    res.message = "模态结果不可用（ModalResult.ok = false）";
    return res;
  }
  const int n = mr.nmodes;
  if (n <= 0 || mr.omega.empty()) {
    res.message = "没有可用的模态";
    return res;
  }
  res.nmodes = n;

  // 1) 每个模态的地震影响系数（按周期 T = 2π/ω）
  res.alphaJ.resize(static_cast<size_t>(n));
  for (int j = 0; j < n; ++j) {
    const double t = (mr.omega[static_cast<size_t>(j)] > 1e-300)
                         ? (2.0 * 3.141592653589793 /
                            mr.omega[static_cast<size_t>(j)])
                         : 0.0;
    res.alphaJ[static_cast<size_t>(j)] = spec.alpha(t);
  }

  // 2) 竖向谱：Z 向影响系数取水平谱 × verticalScale（规范 5.3.1）
  const double alphaV = spec.alphaMax * opt.verticalScale;
  const double scaleV = (spec.alphaMax > 1e-300)
                            ? (alphaV / spec.alphaMax)
                            : 0.0;

  // 3) 各模态基底剪力：V_j = α_j · m*_j · g（m* = γ²，kN）
  //    地震系数本身要么全水平、要么 Z 向乘竖向比，故此处按方向缩放。
  res.vX.resize(static_cast<size_t>(n));
  res.vY.resize(static_cast<size_t>(n));
  res.vZ.resize(static_cast<size_t>(n));
  for (int j = 0; j < n; ++j) {
    const double a = res.alphaJ[static_cast<size_t>(j)];
    const double gx = mr.gammaX.empty()
                          ? 0.0
                          : mr.gammaX[static_cast<size_t>(j)];
    const double gy = mr.gammaY.empty()
                          ? 0.0
                          : mr.gammaY[static_cast<size_t>(j)];
    const double gz = mr.gammaZ.empty()
                          ? 0.0
                          : mr.gammaZ[static_cast<size_t>(j)];
    res.vX[static_cast<size_t>(j)] = a * gx * gx * kG;
    res.vY[static_cast<size_t>(j)] = a * gy * gy * kG;
    res.vZ[static_cast<size_t>(j)] = a * scaleV * gz * gz * kG;
  }

  // 4) 组合基底剪力：SRSS 或 CQC
  if (opt.method == 0) {
    res.baseX = srss(res.vX);
    res.baseY = srss(res.vY);
    res.baseZ = srss(res.vZ);
  } else {
    res.baseX = cqc(spec.zeta, mr.omega, res.vX);
    res.baseY = cqc(spec.zeta, mr.omega, res.vY);
    res.baseZ = cqc(spec.zeta, mr.omega, res.vZ);
  }

  // 5) 等效地震力空间分布（节点序 × 6）
  //    由静力分析同款质量组装得到 M（MPC 折减后），
  //    从 shapes 抽出主自由度振型 φ（自由方程号空间），
  //    f_j = α_j·γ_j·M·φ_j 即为第 j 阶模态地震作用；
  //    再按 SRSS/CQC 组合出最终等效静力荷载。
  const Id nNode = model_.nodeCount();
  res.forceX.assign(static_cast<size_t>(nNode * 6), 0.0);
  res.forceY.assign(static_cast<size_t>(nNode * 6), 0.0);
  res.forceZ.assign(static_cast<size_t>(nNode * 6), 0.0);

  // 抽 φ：节点序 → 方程号（自由自由度 g ≥ 0；从属自由度由 MPC 折入主自由度）
  const size_t ndof = static_cast<size_t>(mr.ndof);
  if (ndof > 0) {
    // 取质量矩阵（与模态分析同源：集中质量 + MPC 折减）
    ModalAnalysis ma(model_);
    ModalAnalysis::Assembly qa;
    std::string why;
    if (!ma.assembleMass(qa, &why)) {
      res.message = "反应谱计算：质量组装失败：" + why;
      return res;
    }
    if (static_cast<size_t>(qa.ndof) != ndof) {
      res.message = "反应谱计算：模态与质量矩阵自由度不一致";
      return res;
    }

    // 逐模态组合：先算 f_j = α·γ·Mφ，按方向放入节点序向量
    std::vector<std::vector<double>> fx(static_cast<size_t>(n));
    std::vector<std::vector<double>> fy(static_cast<size_t>(n));
    std::vector<std::vector<double>> fz(static_cast<size_t>(n));
    for (int j = 0; j < n; ++j) {
      fx[static_cast<size_t>(j)].assign(static_cast<size_t>(nNode * 6), 0.0);
      fy[static_cast<size_t>(j)].assign(static_cast<size_t>(nNode * 6), 0.0);
      fz[static_cast<size_t>(j)].assign(static_cast<size_t>(nNode * 6), 0.0);

      // φ（方程号空间）
      std::vector<double> phi(ndof, 0.0);
      for (Id nd = 0; nd < nNode; ++nd)
        for (int k = 0; k < 6; ++k) {
          const Id g = model_.node(nd).dof[k];
          if (g >= 0)
            phi[static_cast<size_t>(g)] =
                mr.shapes[static_cast<size_t>(j)]
                         [static_cast<size_t>(nd * 6 + k)];
        }

      std::vector<double> Mphi;
      qa.M.multiply(phi, Mphi);

      const double aj = res.alphaJ[static_cast<size_t>(j)];
      const double gx = mr.gammaX.empty()
                            ? 0.0
                            : mr.gammaX[static_cast<size_t>(j)];
      const double gy = mr.gammaY.empty()
                            ? 0.0
                            : mr.gammaY[static_cast<size_t>(j)];
      const double gz = mr.gammaZ.empty()
                            ? 0.0
                            : mr.gammaZ[static_cast<size_t>(j)];

      const double cx = aj * gx * kG;
      const double cy = aj * gy * kG;
      const double cz = aj * scaleV * gz * kG;
      if (cx == 0.0 && cy == 0.0 && cz == 0.0) continue;

      for (Id nd = 0; nd < nNode; ++nd)
        for (int k = 0; k < 3; ++k) {   // 地震作用只产生平动力
          const Id g = model_.node(nd).dof[k];
          if (g < 0) continue;          // 固定/从属自由度上的质量已折入主自由度
          const double mphi = Mphi[static_cast<size_t>(g)];
          const size_t idx = static_cast<size_t>(nd * 6 + k);
          fx[static_cast<size_t>(j)][idx] = cx * mphi;
          fy[static_cast<size_t>(j)][idx] = cy * mphi;
          fz[static_cast<size_t>(j)][idx] = cz * mphi;
        }
    }

    // 组合（SRSS/CQC 均按分量逐点组合）
    auto combine = [&](const std::vector<std::vector<double>>& fj,
                       std::vector<double>& out) {
      const size_t len = static_cast<size_t>(nNode * 6);
      for (size_t i = 0; i < len; ++i) {
        if (opt.method == 0) {
          double s = 0.0;
          for (int j = 0; j < n; ++j) {
            const double x = fj[static_cast<size_t>(j)][i];
            s += x * x;
          }
          out[i] = std::sqrt(s);
        } else {
          // CQC：√(ΣΣ ρ_jk f_j f_k)
          double s = 0.0;
          for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
              s += cqcRho(spec.zeta, mr.omega[static_cast<size_t>(j)],
                          mr.omega[static_cast<size_t>(k)]) *
                   fj[static_cast<size_t>(j)][i] *
                   fj[static_cast<size_t>(k)][i];
          out[i] = std::sqrt(std::max(s, 0.0));
        }
      }
    };
    combine(fx, res.forceX);
    combine(fy, res.forceY);
    combine(fz, res.forceZ);
  }

  // 6) 参与质量比（规范 5.2.2：≥ 90%）。
  //    某方向根本没有可动质量（如轴向模型的 X/Y 平动全被约束）时，
  //    该方向不存在地震作用需求，视为满足（豁免检查）。
  //    checkMass=false 为整体开关：关闭时全部方向视为通过（仅记录比值）。
  res.massRatioX = mr.ratioX;
  res.massRatioY = mr.ratioY;
  res.massRatioZ = mr.ratioZ;
  if (opt.checkMass) {
    const bool noX = mr.totalMassX < 1e-9;
    const bool noY = mr.totalMassY < 1e-9;
    const bool noZ = mr.totalMassZ < 1e-9;
    res.massOkX = noX || res.massRatioX > 0.8999;
    res.massOkY = noY || res.massRatioY > 0.8999;
    res.massOkZ = noZ || res.massRatioZ > 0.8999;
  } else {
    res.massOkX = res.massOkY = res.massOkZ = true;
  }
  res.massOk = res.massOkX && res.massOkY && res.massOkZ;
  if (opt.checkMass && !res.massOk)
    res.message = "求解成功，但参与质量 < 90%（需增加模态数）";

  res.ok = true;
  if (res.message.empty()) res.message = "求解成功";
  return res;
}

SpectrumResult ResponseSpectrum::solve(const ModalResult& mr,
                                       const DesignSpectrum& spec) {
  return solve(mr, spec, Options());
}

}  // namespace yjk