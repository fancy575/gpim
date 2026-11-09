#' @name gpim-package
#' @aliases gpim-package NULL
#' @title Probabilistic Index Regression for Prioritized Survival Outcomes
#' @description Formula-based CH, HT, and DL regression with multiple restriction times, long-format hierarchical outcomes, censoring-adjusted sandwich inference, convergence controls, and coefficient trajectories.
#' @references Fang, X., Cao, Z., Mao, L., and Li, F. (2026). Probabilistic index models for hierarchical survival outcomes. Manuscript under review.
#' @seealso \code{\link{gpim}}, \code{\link{HierSurv}}, \code{\link{hfaction}}, \code{\link{pim_diagnostics}}
#' @useDynLib gpim, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef vcov confint nobs predict
#' @keywords internal
"_PACKAGE"
