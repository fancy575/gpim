// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>

// Construct comparisons from observed histories in their specified priority order.
// Returned subject indices follow R's one-based convention.
// [[Rcpp::export]]
Rcpp::List gpim_pairs_cpp(Rcpp::NumericVector observed_u,
                         Rcpp::LogicalVector death,
                         Rcpp::NumericMatrix event_times,
                         Rcpp::IntegerVector priorities,
                         double L) {
  const int n = observed_u.size();
  const int q = event_times.ncol();
  if (n < 2 || death.size() != n || event_times.nrow() != n)
    Rcpp::stop("Incompatible subject dimensions in hierarchical comparisons.");
  if (!R_finite(L) || L <= 0)
    Rcpp::stop("The restriction time must be positive and finite.");
  if (priorities.size() != q + 1 || priorities[0] != 1)
    Rcpp::stop("Endpoint priorities must begin with fatal code 1.");
  for (int k = 0; k <= q; ++k)
    if (priorities[k] == NA_INTEGER || priorities[k] <= 0)
      Rcpp::stop("Endpoint priorities must be positive integers.");

  const double pair_count = static_cast<double>(n) * (n - 1) / 2;
  if (pair_count > std::numeric_limits<int>::max())
    Rcpp::stop("Too many subject pairs for the model matrix.");
  const int m = static_cast<int>(pair_count);

  Rcpp::NumericVector u(n);
  Rcpp::LogicalVector died(n), censor(n), complete(n);
  for (int a = 0; a < n; ++a) {
    if (!R_finite(observed_u[a]) || observed_u[a] < 0 ||
        death[a] == NA_LOGICAL)
      Rcpp::stop("Terminal histories contain invalid times or death indicators.");
    u[a] = std::min(observed_u[a], L);
    died[a] = death[a] && observed_u[a] <= L;
    censor[a] = !death[a] && observed_u[a] < L;
    complete[a] = died[a] || u[a] == L;
    for (int k = 0; k < q; ++k)
      if (ISNAN(event_times(a, k)) || event_times(a, k) < 0)
        Rcpp::stop("Nonfatal event times must be nonnegative or positive infinity.");
  }

  // Preserve the names supplied by the subject-level data preparation.
  SEXP u_names = observed_u.attr("names");
  SEXP death_names = death.attr("names");
  SEXP flag_names = Rf_isNull(death_names) ? u_names : death_names;
  if (!Rf_isNull(u_names)) u.attr("names") = u_names;
  if (!Rf_isNull(flag_names)) {
    died.attr("names") = flag_names;
    censor.attr("names") = flag_names;
    complete.attr("names") = flag_names;
  }

  Rcpp::IntegerVector i(m), j(m), priority(m, NA_INTEGER);
  Rcpp::NumericVector s(m), h(m, 0.5);
  Rcpp::LogicalVector resolved(m), ch(m);
  Rcpp::CharacterVector pair_u_names, pair_flag_names;
  if (!Rf_isNull(u_names)) pair_u_names = Rcpp::CharacterVector(m);
  if (!Rf_isNull(flag_names)) pair_flag_names = Rcpp::CharacterVector(m);
  Rcpp::CharacterVector named_u = Rf_isNull(u_names) ?
    Rcpp::CharacterVector() : Rcpp::as<Rcpp::CharacterVector>(u_names);
  Rcpp::CharacterVector named_flags = Rf_isNull(flag_names) ?
    Rcpp::CharacterVector() : Rcpp::as<Rcpp::CharacterVector>(flag_names);

  int r = 0;
  for (int a = 0; a < n - 1; ++a) {
    Rcpp::checkUserInterrupt();
    for (int b = a + 1; b < n; ++b, ++r) {
      i[r] = a + 1;
      j[r] = b + 1;
      s[r] = std::min(u[a], u[b]);
      const bool di = died[a] && u[a] == s[r];
      const bool dj = died[b] && u[b] == s[r];
      resolved[r] = s[r] == L || di || dj;
      ch[r] = complete[a] && complete[b];
      if (!Rf_isNull(u_names)) pair_u_names[r] = named_u[a];
      if (!Rf_isNull(flag_names)) pair_flag_names[r] = named_flags[a];
      if (ch[r] && !resolved[r])
        Rcpp::stop("Internal error in complete-history comparison.");

      if (!resolved[r]) continue;
      if (s[r] < L && di && !dj) {
        h[r] = 0;
        priority[r] = 1;
        continue;
      }
      if (s[r] < L && dj && !di) {
        h[r] = 1;
        priority[r] = 1;
        continue;
      }

      priority[r] = 0;
      for (int k = 0; k < q; ++k) {
        const double ya = std::min(event_times(a, k), u[a]);
        const double yb = std::min(event_times(b, k), u[b]);
        if (ya != yb) {
          h[r] = ya > yb ? 1.0 : 0.0;
          priority[r] = priorities[k + 1];
          break;
        }
      }
    }
  }

  if (!Rf_isNull(u_names)) {
    s.attr("names") = pair_u_names;
    resolved.attr("names") = pair_u_names;
  } else if (!Rf_isNull(flag_names)) {
    resolved.attr("names") = pair_flag_names;
  }
  if (!Rf_isNull(flag_names)) ch.attr("names") = pair_flag_names;

  return Rcpp::List::create(
    Rcpp::_["i"] = i, Rcpp::_["j"] = j, Rcpp::_["u"] = u,
    Rcpp::_["died"] = died, Rcpp::_["censor"] = censor,
    Rcpp::_["complete"] = complete, Rcpp::_["s"] = s,
    Rcpp::_["resolved"] = resolved, Rcpp::_["ch"] = ch,
    Rcpp::_["h"] = h, Rcpp::_["priority"] = priority);
}

