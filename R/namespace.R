#' @name gpim-package
#' @aliases gpim-package NULL
#' @title Probabilistic Index Regression for Prioritized Survival Outcomes
#' @description Formula-based CHW, RCW, and CTW regression with multiple restriction times, observed hierarchical event histories, selectable copula and margin models, subject influence covariance estimation, and coefficient trajectories.
#' @references Fang, X., Cao, Z., Mao, L., and Li, F. (2026). Probabilistic index models for hierarchical survival outcomes. Manuscript under review.
#' @seealso \code{\link{gpim}}, \code{\link{HierSurv}}, \code{\link{hfaction}}, \code{\link{pim_joint_model}}, \code{\link{pim_diagnostics}}
#' @useDynLib gpim, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef vcov confint nobs predict logLik influence
#' @keywords internal
"_PACKAGE"
