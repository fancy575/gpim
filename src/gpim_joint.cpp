// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include "gpim_copula.h"

namespace {

struct MarginValue {
  double hazard = 0, log_density_factor = 0;
  double b0_H = 0, b1_H = 0, b0_h = 0, b1_h = 0;
  double beta_H = 0, beta_h = 0;
  int last = 0, event = -1;
  bool valid = true;
};

struct Margin {
  std::string type;
  Rcpp::NumericMatrix X;
  Rcpp::NumericVector grid;
  std::vector<double> par, jumps, cumulative;
  int offset, baseline, p, length;
  bool valid = true;

  Margin(const std::string& type_, Rcpp::NumericMatrix X_,
         Rcpp::NumericVector grid_, Rcpp::NumericVector all, int offset_) :
    type(type_), X(X_), grid(grid_), offset(offset_), p(X_.ncol()) {
    if (type == "cox") baseline = grid.size();
    else if (type == "weibull" || type == "lognormal") baseline = 2;
    else if (type == "exponential") baseline = 1;
    else Rcpp::stop("Unknown joint-model margin.");
    length = baseline + p;
    if (offset + length > all.size())
      Rcpp::stop("The joint-model parameter vector is too short.");
    par.assign(all.begin() + offset, all.begin() + offset + length);
    if (type == "cox") {
      jumps.resize(baseline);
      cumulative.assign(baseline + 1, 0);
      for (int k = 0; k < baseline; ++k) {
        if (!R_finite(grid[k]) || grid[k] < 0 || (k && grid[k] <= grid[k - 1]))
          Rcpp::stop("Cox baseline event grids must be finite, nonnegative and strictly increasing.");
        jumps[k] = std::exp(par[k]);
        cumulative[k + 1] = cumulative[k] + jumps[k];
        if (!(jumps[k] > 0) || !R_finite(cumulative[k + 1])) valid = false;
      }
    }
  }

  MarginValue evaluate(int i, double t, bool event) const {
    MarginValue ans;
    double eta = 0;
    for (int j = 0; j < p; ++j) eta += X(i, j) * par[baseline + j];
    if (type == "cox") {
      ans.last = std::upper_bound(grid.begin(), grid.end(), t) - grid.begin();
      const double risk = std::exp(eta);
      ans.hazard = risk * cumulative[ans.last];
      ans.beta_H = ans.beta_h = 1;
      if (event) {
        if (!ans.last || grid[ans.last - 1] != t) {
          ans.valid = false;
          return ans;
        }
        ans.event = ans.last - 1;
        ans.log_density_factor = par[ans.event] + eta;
      }
      ans.valid = valid && R_finite(risk) && R_finite(ans.hazard) &&
        R_finite(ans.log_density_factor);
      return ans;
    }
    if (t == 0) {
      ans.valid = !event;
      return ans;
    }
    const double lt = std::log(t);
    if (type == "weibull" || type == "exponential") {
      const double shape = type == "weibull" ? std::exp(par[1]) : 1;
      const double log_H = par[0] + shape * lt + eta;
      ans.hazard = std::exp(log_H);
      ans.log_density_factor = par[0] + std::log(shape) + (shape - 1) * lt + eta;
      ans.b0_H = ans.b0_h = 1;
      ans.b1_H = shape * lt;
      ans.b1_h = 1 + shape * lt;
      ans.beta_H = ans.beta_h = 1;
    } else {
      const double sd = std::exp(par[1]);
      const double z = (lt - par[0] - eta) / sd;
      const double ls = R::pnorm(z, 0.0, 1.0, false, true);
      const double lp = R::dnorm(z, 0.0, 1.0, true);
      const double mills = std::exp(lp - ls);
      ans.hazard = -ls;
      const double log_H = ans.hazard > 0 ? std::log(ans.hazard) :
        R::pnorm(z, 0.0, 1.0, true, true);
      const double ratio = std::exp(lp - ls - log_H);
      ans.log_density_factor = lp - par[1] - lt - ls;
      ans.b0_H = -ratio / sd;
      ans.b1_H = -z * ratio;
      ans.b0_h = (z - mills) / sd;
      ans.b1_h = z * z - 1 - z * mills;
      ans.beta_H = ans.b0_H;
      ans.beta_h = ans.b0_h;
    }
    ans.valid = R_finite(ans.hazard) && R_finite(ans.log_density_factor) &&
      R_finite(ans.b0_H) && R_finite(ans.b1_H) &&
      R_finite(ans.b0_h) && R_finite(ans.b1_h);
    return ans;
  }

  void add_total(Rcpp::NumericVector& gradient, std::vector<double>& buckets,
                 const MarginValue& v, int i, double g, double event) const {
    if (type == "cox") {
      if (v.last && cumulative[v.last] > 0)
        buckets[v.last - 1] += g / cumulative[v.last];
      if (event && v.event >= 0) gradient[offset + v.event] += event;
    } else {
      gradient[offset] += g * v.b0_H + event * v.b0_h;
      if (baseline == 2)
        gradient[offset + 1] += g * v.b1_H + event * v.b1_h;
    }
    const double gb = g * v.beta_H + event * v.beta_h;
    for (int j = 0; j < p; ++j)
      gradient[offset + baseline + j] += gb * X(i, j);
  }

  void finish_total(Rcpp::NumericVector& gradient,
                    const std::vector<double>& buckets) const {
    if (type != "cox") return;
    double remaining = 0;
    for (int k = baseline - 1; k >= 0; --k) {
      remaining += buckets[k];
      gradient[offset + k] += jumps[k] * remaining;
    }
  }

  void add_row(Rcpp::NumericMatrix& score, const MarginValue& v, int i,
               double g, double event, double multiplier = 1) const {
    if (type == "cox") {
      if (v.last && cumulative[v.last] > 0) {
        const double factor = multiplier * g / cumulative[v.last];
        for (int k = 0; k < v.last; ++k)
          score(i, offset + k) += factor * jumps[k];
      }
      if (event && v.event >= 0) score(i, offset + v.event) += multiplier * event;
    } else {
      score(i, offset) += multiplier * (g * v.b0_H + event * v.b0_h);
      if (baseline == 2)
        score(i, offset + 1) += multiplier * (g * v.b1_H + event * v.b1_h);
    }
    const double gb = multiplier * (g * v.beta_H + event * v.beta_h);
    for (int j = 0; j < p; ++j)
      score(i, offset + baseline + j) += gb * X(i, j);
  }
};

struct JointData {
  Rcpp::NumericVector u, v;
  Rcpp::LogicalVector death, nonfatal;
  Rcpp::NumericMatrix XD, XT;
  int n;
  JointData(Rcpp::List data) : u(data["u"]), v(data["v"]),
    death(Rcpp::as<Rcpp::LogicalVector>(data["death"])),
    nonfatal(Rcpp::as<Rcpp::LogicalVector>(data["nonfatal"])),
    XD(data["XD"]), XT(data["XT"]), n(u.size()) {
    if (!n || v.size() != n || death.size() != n || nonfatal.size() != n ||
        XD.nrow() != n || XT.nrow() != n)
      Rcpp::stop("Incompatible dimensions in joint-model records.");
    for (int i = 0; i < n; ++i) {
      if (!R_finite(u[i]) || u[i] < 0 || !R_finite(v[i]) || v[i] < 0 ||
          v[i] > u[i] || death[i] == NA_LOGICAL || nonfatal[i] == NA_LOGICAL)
        Rcpp::stop("Invalid joint-model history. Nonfatal events must be recorded by the terminal observation time.");
      if (death[i] && nonfatal[i] && !(v[i] < u[i]))
        Rcpp::stop("A nonfatal event must precede death in the continuous joint model.");
    }
    for (double x : XD) if (!R_finite(x)) Rcpp::stop("Joint-model covariates must be finite.");
    for (double x : XT) if (!R_finite(x)) Rcpp::stop("Joint-model covariates must be finite.");
  }
};

Rcpp::List invalid_result(int n, int p, bool individual) {
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["loglik"] = R_NegInf,
    Rcpp::_["gradient"] = Rcpp::NumericVector(p, NA_REAL),
    Rcpp::_["loglik_i"] = Rcpp::NumericVector(n, NA_REAL));
  if (individual) {
    Rcpp::NumericMatrix scores(n, p);
    std::fill(scores.begin(), scores.end(), NA_REAL);
    out["scores"] = scores;
  }
  return out;
}

Rcpp::NumericVector get_grid(Rcpp::List spec, const char* name) {
  return spec.containsElementNamed(name) ? Rcpp::as<Rcpp::NumericVector>(spec[name]) :
    Rcpp::NumericVector();
}

bool score_finite(const gpim_joint::LocalScore& g) {
  return R_finite(g.value) && R_finite(g.log_hd) &&
    R_finite(g.log_ht) && R_finite(g.dependence);
}

} // namespace

// The Cox option uses the baseline-increment likelihood.  Every observed
// event contributes its log baseline jump, not a derivative of a step function.
// The gradient includes the entire cumulative-hazard contribution of each
// individual, including those without an observed event.  In the counting-
// process representation these are the event and compensator score terms.
// [[Rcpp::export]]
Rcpp::List gpim_joint_cpp(Rcpp::NumericVector par, Rcpp::List data,
                          Rcpp::List spec, bool individual = false) {
  JointData d(data);
  const int family = gpim_joint::copula_code(Rcpp::as<std::string>(spec["copula"]));
  if (family < 0) Rcpp::stop("Unknown joint-model copula.");
  for (double a : par) if (!R_finite(a)) return invalid_result(d.n, par.size(), individual);
  Margin md(Rcpp::as<std::string>(spec["marginD"]), d.XD, get_grid(spec, "gridD"), par, 0);
  Margin mt(Rcpp::as<std::string>(spec["marginT"]), d.XT, get_grid(spec, "gridT"), par, md.length);
  const int dependence_index = md.length + mt.length;
  const int npar = dependence_index + (family != 0);
  if (par.size() != npar) Rcpp::stop("Incorrect joint-model parameter-vector length.");
  if (!md.valid || !mt.valid) return invalid_result(d.n, npar, individual);
  const double a = family ? par[dependence_index] : 0;
  Rcpp::NumericVector gradient(npar), loglik_i(d.n);
  Rcpp::NumericMatrix scores(individual ? d.n : 0, individual ? npar : 0);
  std::vector<double> bucketD(md.baseline, 0), bucketT(mt.baseline, 0);
  double loglik = 0;
  for (int i = 0; i < d.n; ++i) {
    if (i % 128 == 0) Rcpp::checkUserInterrupt();
    const bool dd = d.death[i], dt = d.nonfatal[i];
    const double t = dt ? d.v[i] : d.u[i];
    const MarginValue vd = md.evaluate(i, d.u[i], dd);
    const MarginValue vt = mt.evaluate(i, t, dt);
    if (!vd.valid || !vt.valid) return invalid_result(d.n, npar, individual);
    const gpim_joint::LocalScore g = gpim_joint::local_score(
      vd.hazard, vt.hazard, a, family, dd, dt);
    if (!score_finite(g)) return invalid_result(d.n, npar, individual);
    const double li = g.value + (dd ? vd.log_density_factor : 0) +
      (dt ? vt.log_density_factor : 0);
    if (!R_finite(li)) return invalid_result(d.n, npar, individual);
    loglik_i[i] = li;
    loglik += li;
    md.add_total(gradient, bucketD, vd, i, g.log_hd, dd);
    mt.add_total(gradient, bucketT, vt, i, g.log_ht, dt);
    if (family) gradient[dependence_index] += g.dependence;
    if (individual) {
      md.add_row(scores, vd, i, g.log_hd, dd);
      mt.add_row(scores, vt, i, g.log_ht, dt);
      if (family) scores(i, dependence_index) = g.dependence;
    }
  }
  md.finish_total(gradient, bucketD);
  mt.finish_total(gradient, bucketT);
  for (double g : gradient) if (!R_finite(g)) return invalid_result(d.n, npar, individual);
  Rcpp::List out = Rcpp::List::create(Rcpp::_["loglik"] = loglik,
    Rcpp::_["gradient"] = gradient, Rcpp::_["loglik_i"] = loglik_i);
  if (individual) out["scores"] = scores;
  return out;
}

// Conditional probabilities use only the record restricted at L.  The model
// parameters may have been estimated from the longer fitting record.
// [[Rcpp::export]]
Rcpp::List gpim_joint_q_cpp(Rcpp::NumericVector par, Rcpp::List data,
                            Rcpp::List spec, double L) {
  if (!R_finite(L) || L <= 0) Rcpp::stop("L must be positive and finite.");
  JointData d(data);
  const int family = gpim_joint::copula_code(Rcpp::as<std::string>(spec["copula"]));
  if (family < 0) Rcpp::stop("Unknown joint-model copula.");
  for (double a : par) if (!R_finite(a)) Rcpp::stop("Joint-model parameters must be finite.");
  Margin md(Rcpp::as<std::string>(spec["marginD"]), d.XD, get_grid(spec, "gridD"), par, 0);
  Margin mt(Rcpp::as<std::string>(spec["marginT"]), d.XT, get_grid(spec, "gridT"), par, md.length);
  const int dependence_index = md.length + mt.length;
  const int npar = dependence_index + (family != 0);
  if (par.size() != npar) Rcpp::stop("Incorrect joint-model parameter-vector length.");
  if (!md.valid || !mt.valid) Rcpp::stop("Invalid fitted joint-model baseline.");
  const double a = family ? par[dependence_index] : 0;
  Rcpp::NumericVector q(d.n), logq(d.n);
  Rcpp::NumericMatrix dq(d.n, npar);
  for (int i = 0; i < d.n; ++i) {
    if (i % 128 == 0) Rcpp::checkUserInterrupt();
    // A death recorded exactly at L has min(D,L)=L, matching the pair rule.
    if (d.u[i] >= L) {
      q[i] = 1;
      logq[i] = 0;
      continue;
    }
    if (d.death[i]) {
      q[i] = 0;
      logq[i] = R_NegInf;
      continue;
    }
    const bool dt = d.nonfatal[i];
    const double t = dt ? d.v[i] : d.u[i];
    const MarginValue vl = md.evaluate(i, L, false);
    const MarginValue vu = md.evaluate(i, d.u[i], false);
    const MarginValue vt = mt.evaluate(i, t, false);
    if (!vl.valid || !vu.valid || !vt.valid)
      Rcpp::stop("Nonfinite cumulative hazard in a conditional tie probability.");
    const gpim_joint::LocalScore gl = gpim_joint::local_score(
      vl.hazard, vt.hazard, a, family, 0, dt, true);
    const gpim_joint::LocalScore gu = gpim_joint::local_score(
      vu.hazard, vt.hazard, a, family, 0, dt, true);
    if (!score_finite(gl) || !score_finite(gu))
      Rcpp::stop("The fitted joint model gives an undefined conditional tie probability.");
    const double lq = gl.value - gu.value;
    if (!R_finite(lq) || lq > 1.0e-8)
      Rcpp::stop("Invalid conditional tie probability from the fitted joint model.");
    logq[i] = std::min(0.0, lq);
    q[i] = std::exp(logq[i]);
    md.add_row(dq, vl, i, gl.log_hd, 0, q[i]);
    md.add_row(dq, vu, i, gu.log_hd, 0, -q[i]);
    mt.add_row(dq, vt, i, gl.log_ht - gu.log_ht, 0, q[i]);
    if (family) dq(i, dependence_index) = q[i] * (gl.dependence - gu.dependence);
  }
  return Rcpp::List::create(Rcpp::_["q"] = q, Rcpp::_["logq"] = logq,
                            Rcpp::_["dq"] = dq);
}
