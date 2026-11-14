#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

namespace {
inline void pim_link_values(double eta, int link, double& mu, double& comp,
                            double& prime, double& k, double& kp) {
  if (link == 0) {
    mu = (1.0 + eta) / 2.0; comp = (1.0 - eta) / 2.0;
    prime = 0.5; k = 1.0; kp = 0.0;
  } else if (link == 1) {
    const double e = std::exp(-std::abs(eta));
    mu = eta >= 0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
    comp = eta >= 0 ? e / (1.0 + e) : 1.0 / (1.0 + e);
    prime = mu * comp; k = 1.0; kp = 0.0;
  } else {
    mu = R::pnorm(eta, 0.0, 1.0, true, false);
    comp = R::pnorm(eta, 0.0, 1.0, false, false);
    prime = R::dnorm(eta, 0.0, 1.0, false);
    k = std::exp(R::dnorm(eta, 0.0, 1.0, true) -
                 R::pnorm(eta, 0.0, 1.0, true, true) -
                 R::pnorm(eta, 0.0, 1.0, false, true));
    kp = k * (-eta - k * (comp - mu));
  }
}
}

// Evaluate the pair equation and its derivative. The final pass accumulates
// subject projections and nuisance derivatives without a pair-score matrix.
// [[Rcpp::export]]
Rcpp::List gpim_equation_cpp(Rcpp::NumericMatrix Z, Rcpp::NumericVector beta,
                             Rcpp::NumericVector y, Rcpp::NumericVector weight,
                             int link, bool details, Rcpp::List state) {
  const int m = Z.nrow(), p = Z.ncol();
  if (!m || beta.size() != p || y.size() != m || weight.size() != m ||
      link < 0 || link > 2) Rcpp::stop("Invalid pair-equation dimensions.");
  const int n = Rcpp::as<int>(state["n"]);
  Rcpp::IntegerVector ii, jj, li, lj;
  Rcpp::LogicalVector active;
  Rcpp::NumericVector risk, h0, q, factor;
  Rcpp::NumericMatrix W, ac;
  int ne = 0, nc = 0;
  if (details) {
    ii = state["i"]; jj = state["j"];
    li = state["left_i"]; lj = state["left_j"];
    active = state["derivative"]; risk = state["risk"];
    h0 = state["h0"];
    W = Rcpp::as<Rcpp::NumericMatrix>(state["W"]);
    ac = Rcpp::as<Rcpp::NumericMatrix>(state["a_cum"]);
    q = state["q"]; factor = state["factor"];
    ne = h0.size(); nc = W.ncol();
    if (ii.size() != m || jj.size() != m || li.size() != m ||
        lj.size() != m || active.size() != m || factor.size() != m ||
        risk.size() != n || q.size() != n || W.nrow() != n ||
        ac.nrow() != ne || ac.ncol() != nc)
      Rcpp::stop("Invalid subject or nuisance dimensions.");
  }
  Rcpp::NumericVector score(p), eta(details ? m : 0);
  Rcpp::NumericMatrix jacobian(p, p), rows(details ? n : 0, p),
    buckets(ne, p), gamma(p, nc), loading(details ? n : 0, p);
  double min_mu = R_PosInf, max_mu = R_NegInf, max_eta = 0.0;
  std::vector<double> nuisance_gamma(nc, 0.0);
  for (int r = 0; r < m; ++r) {
    if ((r & 16383) == 0) Rcpp::checkUserInterrupt();
    double lin = 0.0;
    for (int k = 0; k < p; ++k) lin += Z(r, k) * beta[k];
    double mu, comp, prime, mult, multp;
    pim_link_values(lin, link, mu, comp, prime, mult, multp);
    const double resid = y[r] * comp - (1.0 - y[r]) * mu;
    const double sc = weight[r] * mult * resid;
    const double jc = weight[r] * (multp * resid - mult * prime);
    int a = 0, b = 0, ia = 0, ib = 0;
    if (details) {
      eta[r] = lin; min_mu = std::min(min_mu, mu);
      max_mu = std::max(max_mu, mu); max_eta = std::max(max_eta, std::abs(lin));
      a = ii[r] - 1; b = jj[r] - 1; ia = li[r]; ib = lj[r];
      if (a < 0 || b < 0 || a >= n || b >= n || ia < 0 || ib < 0 ||
          ia > ne || ib > ne) Rcpp::stop("Invalid pair subject or time index.");
      if (active[r]) for (int c = 0; c < nc; ++c) {
        nuisance_gamma[c] =
          risk[a] * (ia ? h0[ia-1] * W(a,c) - ac(ia-1,c) : 0.0) +
          risk[b] * (ib ? h0[ib-1] * W(b,c) - ac(ib-1,c) : 0.0);
      }
    }
    for (int k = 0; k < p; ++k) {
      const double value = Z(r,k) * sc;
      score[k] += value / m;
      for (int l = 0; l <= k; ++l)
        jacobian(k,l) += Z(r,k) * Z(r,l) * jc / m;
      if (!details) continue;
      rows(a,k) += value; rows(b,k) += value;
      if (active[r]) {
        if (ia) buckets(ia-1,k) += value * risk[a] / m;
        if (ib) buckets(ib-1,k) += value * risk[b] / m;
        for (int c = 0; c < nc; ++c)
          gamma(k,c) += value * nuisance_gamma[c] / m;
      }
      if (factor[r] != 0.0) {
        const double v = Z(r,k) * mult * resid * factor[r] / m;
        loading(a,k) += v * q[b];
        loading(b,k) += v * q[a];
      }
    }
  }
  for (int k = 0; k < p; ++k)
    for (int l = 0; l < k; ++l) jacobian(l,k) = jacobian(k,l);
  return Rcpp::List::create(
    Rcpp::_["mean_score"] = score, Rcpp::_["mean_jacobian"] = jacobian,
    Rcpp::_["subject_sum"] = rows, Rcpp::_["censor_buckets"] = buckets,
    Rcpp::_["censor_gamma"] = gamma, Rcpp::_["joint_loading"] = loading,
    Rcpp::_["eta"] = eta,
    Rcpp::_["fitted_range"] = Rcpp::NumericVector::create(min_mu,max_mu),
    Rcpp::_["max_abs_eta"] = max_eta);
}

// Additional observed nonfatal comparisons for one nonfatal component.
// Recorded fatal ties and genuine complete-history ties keep the RCW term.
// [[Rcpp::export]]
Rcpp::List gpim_ctw_pairs_cpp(Rcpp::List pairs, Rcpp::NumericVector event_time,
                              double L) {
  Rcpp::NumericVector u = pairs["u"], s = pairs["s"], h = pairs["h"];
  Rcpp::IntegerVector ii = pairs["i"], jj = pairs["j"];
  Rcpp::LogicalVector died = pairs["died"], resolved = pairs["resolved"];
  const int n = u.size(), m = ii.size();
  if (event_time.size() != n) Rcpp::stop("Invalid nonfatal event vector.");
  Rcpp::NumericVector outcome(m,0.5), time(m,0.0);
  Rcpp::LogicalVector copied(m,false), lower(m,false);
  for (int r = 0; r < m; ++r) {
    if ((r & 32767) == 0) Rcpp::checkUserInterrupt();
    const int a = ii[r]-1, b = jj[r]-1;
    if (resolved[r] && (s[r] < L || h[r] == 0.5)) {
      copied[r] = true; outcome[r] = h[r]; time[r] = s[r]; continue;
    }
    const double va = event_time[a] <= u[a] && event_time[a] < L ? event_time[a] : R_PosInf;
    const double vb = event_time[b] <= u[b] && event_time[b] < L ? event_time[b] : R_PosInf;
    const double first = std::min(va,vb);
    if (va != vb && first < L && u[a] >= first && u[b] >= first &&
        !(died[a] && u[a] <= first) && !(died[b] && u[b] <= first)) {
      lower[r] = true; time[r] = first;
      outcome[r] = va > vb ? 1.0 : (va < vb ? 0.0 : 0.5);
    }
  }
  return Rcpp::List::create(Rcpp::_["y"] = outcome,
    Rcpp::_["time"] = time, Rcpp::_["copied"] = copied,
    Rcpp::_["lower"] = lower);
}
