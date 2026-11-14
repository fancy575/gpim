#ifndef GPIM_COPULA_H
#define GPIM_COPULA_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace gpim_joint {

// All copulas below are survival copulas.  Arguments are cumulative hazards,
// so the first and second uniform arguments are exp(-x) and exp(-y).
// Long-double arithmetic protects the subtraction in Frank and Plackett.
struct CopulaLogs {
  long double c, du, dv, duv;
};

inline int copula_code(const std::string& name) {
  if (name == "independence") return 0;
  if (name == "gumbel") return 1;
  if (name == "clayton") return 2;
  if (name == "frank") return 3;
  if (name == "plackett") return 4;
  return -1;
}

inline long double log_add(long double a, long double b) {
  if (a == -std::numeric_limits<long double>::infinity()) return b;
  if (b == -std::numeric_limits<long double>::infinity()) return a;
  const long double m = std::max(a, b);
  return m + std::log1p(std::exp(std::min(a, b) - m));
}

inline CopulaLogs copula_logs(long double x, long double y,
                              double coordinate, int family) {
  const long double inf = std::numeric_limits<long double>::infinity();
  const long double nan = std::numeric_limits<long double>::quiet_NaN();
  if (!(x >= 0 && y >= 0) || !std::isfinite(x) || !std::isfinite(y))
    return {nan, nan, nan, nan};
  if (family == 0) return {-x - y, -y, -x, 0};
  const long double k = std::exp(static_cast<long double>(coordinate));
  if (!std::isfinite(k) || !(k > 0)) return {nan, nan, nan, nan};

  if (family == 1) {
    const long double alpha = 1 + k;
    if (alpha == 1) return {-x - y, -y, -x, 0};
    if (x == 0 && y == 0) return {0, 0, 0, inf};
    if (x == 0) return {-y, -inf, 0, -inf};
    if (y == 0) return {-x, 0, -inf, -inf};
    const long double lx = std::log(x), ly = std::log(y);
    const long double lm = std::max(lx, ly);
    const long double lt = lm +
      std::log(std::exp(alpha * (lx - lm)) +
               std::exp(alpha * (ly - lm))) / alpha;
    const long double t = std::exp(lt);
    return {-t,
            -t + k * (lx - lt) + x,
            -t + k * (ly - lt) + y,
            -t + x + y + k * (lx + ly - 2 * lt) +
              std::log1p(k / t)};
  }

  if (family == 2) {
    const long double kx = k * x, ky = k * y;
    long double la;
    if (std::max(kx, ky) < 0.5L) {
      la = std::log1p(std::expm1(kx) + std::expm1(ky));
    } else {
      const long double m = std::max(kx, ky);
      la = m + std::log(std::exp(kx - m) + std::exp(ky - m) -
                       std::exp(-m));
    }
    return {-la / k,
            (k + 1) * x - (1 / k + 1) * la,
            (k + 1) * y - (1 / k + 1) * la,
            std::log1p(k) + (k + 1) * (x + y) - (1 / k + 2) * la};
  }

  const long double u = std::exp(-x), v = std::exp(-y);
  if (family == 3) {
    const long double ku = k * u, kv = k * v;
    const long double a = -std::expm1(-ku);
    const long double b = -std::expm1(-kv);
    const long double d = -std::expm1(-k);
    if (!(a > 0 && b > 0 && d > 0)) return {nan, nan, nan, nan};
    const long double la = std::log(a), lb = std::log(b), ld = std::log(d);
    // d-a*b = exp(-k*u)*(1-exp(-k*v))
    //         + exp(-k*v)*(1-exp(-k*(1-v))).
    const long double second = -std::expm1(-k * (1 - v));
    const long double lden = log_add(-ku + lb,
      second > 0 ? -kv + std::log(second) : -inf);
    const long double lp = la + lb - ld;
    const long double linner = lp < -0.5L ?
      std::log1p(-std::exp(lp)) : lden - ld;
    const long double lc = linner < 0 ? std::log(-linner) - std::log(k) :
      (x == 0 && y == 0 ? 0 : nan);
    return {lc, -ku + lb - lden, -kv + la - lden,
            std::log(k) + ld - ku - kv - 2 * lden};
  }

  if (family == 4) {
    const long double a = 1 + (k - 1) * (u + v);
    const long double disc = k >= 1 ?
      std::pow(1 + (k - 1) * (u - v), 2) + 4 * (k - 1) * v * (1 - u) :
      a * a - 4 * k * (k - 1) * u * v;
    if (!(disc > 0)) return {nan, nan, nan, nan};
    const long double sq = std::sqrt(disc);
    const long double denom = a + sq;
    if (!(denom > 0)) return {nan, nan, nan, nan};
    const long double lc = std::log(2 * k) - x - y - std::log(denom);
    auto derivative = [&](long double w) {
      if (w == 1) return 1.0L;
      if (w == 0) return 0.0L;
      const long double z = a - 2 * k * w;
      if (z < 0) return (1 - z / sq) / 2;
      return 2 * k * w * (1 - w) / (sq * (sq + z));
    };
    const long double cu = derivative(v), cv = derivative(u);
    const long double numerator = 1 + (k - 1) * (u + v - 2 * u * v);
    return {lc, cu > 0 ? std::log(cu) : -inf,
            cv > 0 ? std::log(cv) : -inf,
            std::log(k) + std::log(numerator) - 1.5L * std::log(disc)};
  }
  return {nan, nan, nan, nan};
}

inline long double pattern_value(long double x, long double y, double a,
                                 int family, int death, int nonfatal,
                                 bool conditional = false) {
  const CopulaLogs z = copula_logs(x, y, a, family);
  if (conditional) return nonfatal ? z.dv : z.c;
  if (death && nonfatal) return z.duv - x - y;
  if (death) return z.du - x;
  if (nonfatal) return z.dv - y;
  return z.c;
}

struct LocalScore {
  double value, log_hd, log_ht, dependence;
};

inline LocalScore local_score(double hd, double ht, double a,
                              int family, int death, int nonfatal,
                              bool conditional = false) {
  const long double x = hd, y = ht;
  const long double h = 1.0e-4L;
  auto value = [&](long double xx, long double yy, double aa) {
    return pattern_value(xx, yy, aa, family, death, nonfatal, conditional);
  };
  const long double f = value(x, y, a);
  if (!std::isfinite(f))
    return {static_cast<double>(f), NAN, NAN, NAN};
  if (family == 0) {
    // log C=-x-y and log Cv=-x.  Event density factors cancel in q.
    return {static_cast<double>(f), -hd,
            conditional && nonfatal ? 0.0 : -ht, 0.0};
  }
  auto derivative = [&](int axis) {
    if ((axis == 0 && x == 0) || (axis == 1 && y == 0)) return 0.0;
    auto at = [&](int j) {
      if (axis == 0) return value(x * std::exp(j * h), y, a);
      if (axis == 1) return value(x, y * std::exp(j * h), a);
      return value(x, y, a + static_cast<double>(j * h));
    };
    return static_cast<double>((at(-2) - 8 * at(-1) + 8 * at(1) - at(2)) /
                               (12 * h));
  };
  return {static_cast<double>(f), derivative(0), derivative(1), derivative(2)};
}

} // namespace gpim_joint

#endif
