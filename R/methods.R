# S3 presentation methods. Inference remains on the fitted coefficient scale.

.pim_check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'conf.level' (or 'level') must be one number strictly between 0 and 1.",
         call. = FALSE)
  }
  level
}

.pim_same_time <- function(x, y) {
  abs(x - y) <= sqrt(.Machine$double.eps) * pmax(1, abs(x), abs(y))
}

.pim_result_table <- function(object, conf.level = object$conf.level,
                              method = NULL, link = NULL, L = NULL) {
  if (!inherits(object, "prior_pim")) {
    stop("'object' must be a fitted 'prior_pim' object.", call. = FALSE)
  }
  conf.level <- .pim_check_level(conf.level)
  tab <- object$table
  if (!is.data.frame(tab) || !all(c("L", "method", "link", "term",
      "estimate", "std.error", "converged", "boundary") %in% names(tab))) {
    stop("The fitted object does not contain a complete coefficient table.",
         call. = FALSE)
  }
  if (!is.null(method)) {
    if (!is.character(method) || anyNA(method)) {
      stop("'method' must contain method names.", call. = FALSE)
    }
    method <- unique(toupper(method))
    if (!identical(method, "ALL")) {
      bad <- setdiff(method, unique(as.character(tab$method)))
      if (length(bad)) stop("Unavailable method: ", paste(bad, collapse = ", "),
                            call. = FALSE)
      tab <- tab[tab$method %in% method, , drop = FALSE]
    }
  }
  if (!is.null(link)) {
    if (!is.character(link) || anyNA(link)) {
      stop("'link' must contain link names.", call. = FALSE)
    }
    link <- unique(tolower(link))
    if (!identical(link, "all")) {
      bad <- setdiff(link, unique(as.character(tab$link)))
      if (length(bad)) stop("Unavailable link: ", paste(bad, collapse = ", "),
                            call. = FALSE)
      tab <- tab[tab$link %in% link, , drop = FALSE]
    }
  }
  if (!is.null(L)) {
    if (!is.numeric(L) || !length(L) || any(!is.finite(L))) {
      stop("'L' must be a nonempty finite numeric vector.", call. = FALSE)
    }
    available <- unique(tab$L)
    found <- vapply(L, function(tt) any(.pim_same_time(available, tt)),
                    logical(1))
    if (any(!found)) {
      stop("No fitted model at restriction time(s): ",
           paste(format(L[!found], digits = 8), collapse = ", "),
           ". Use L values stored in the fitted object.", call. = FALSE)
    }
    keep <- vapply(tab$L, function(tt) any(.pim_same_time(tt, L)),
                   logical(1))
    tab <- tab[keep, , drop = FALSE]
  }
  if (!nrow(tab)) stop("No fitted models match the requested selection.",
                       call. = FALSE)
  good <- tab$converged %in% TRUE & !(tab$boundary %in% TRUE) &
    is.finite(tab$estimate) & is.finite(tab$std.error) & tab$std.error > 0
  tab$statistic <- tab$p.value <- tab$conf.low <- tab$conf.high <- NA_real_
  critical <- stats::qnorm((1 + conf.level) / 2)
  tab$statistic[good] <- tab$estimate[good] / tab$std.error[good]
  tab$p.value[good] <- 2 * stats::pnorm(-abs(tab$statistic[good]))
  tab$conf.low[good] <- tab$estimate[good] - critical * tab$std.error[good]
  tab$conf.high[good] <- tab$estimate[good] + critical * tab$std.error[good]
  rownames(tab) <- NULL
  tab
}

.pim_select_parameters <- function(tab, parameters) {
  available <- unique(as.character(tab$term))
  if (is.null(parameters)) return(tab)
  if (is.numeric(parameters)) {
    if (!length(parameters) || anyNA(parameters) ||
        any(parameters != as.integer(parameters)) ||
        any(parameters < 1 | parameters > length(available))) {
      stop("Parameter indices are outside the available coefficient range.",
           call. = FALSE)
    }
    parameters <- available[parameters]
  }
  if (!is.character(parameters) || !length(parameters) || anyNA(parameters)) {
    stop("'parameters' (or 'parm') must contain coefficient names or indices.",
         call. = FALSE)
  }
  bad <- setdiff(parameters, available)
  if (length(bad)) stop("Unknown coefficient: ", paste(bad, collapse = ", "),
                        ". Available coefficients: ",
                        paste(available, collapse = ", "), call. = FALSE)
  tab[tab$term %in% parameters, , drop = FALSE]
}

.pim_fit_label <- function(fit) {
  paste0(fit$method, "/", fit$link, "/L=", format(fit$L, digits = 10,
                                                 trim = TRUE))
}

.pim_selected_fits <- function(object, method = NULL, link = NULL,
                               L = NULL) {
  tab <- .pim_result_table(object, method = method, link = link, L = L)
  keys <- unique(tab[c("L", "method", "link")])
  keep <- vapply(object$fits, function(ff) {
    any(keys$method == ff$method & keys$link == ff$link &
          .pim_same_time(keys$L, ff$L))
  }, logical(1))
  fits <- object$fits[keep]
  if (!length(fits)) stop("Selected model components are missing from the object.",
                          call. = FALSE)
  names(fits) <- make.unique(vapply(fits, .pim_fit_label, character(1)))
  fits
}

#' @rdname prior_pim_methods
#' @order 9
#' @export
trajectory_data <- function(object, parameters = NULL, method = NULL,
                            link = NULL, L = NULL,
                            conf.level = object$conf.level) {
  tab <- .pim_result_table(object, conf.level, method, link, L)
  .pim_select_parameters(tab, parameters)
}

#' @name pim_diagnostics
#' @aliases pim_diagnostics
#' @title Inspect Convergence, Censoring Weights, and Pair Counts
#' @description Return one row for each fitted method, link, and restriction time. The table
#' collects convergence results, the number of available comparisons, and
#' diagnostics for the censoring weights and covariance estimate.
#' @param object A fitted object returned by \code{\link{gpim}}.
#' @param method Optional vector of fitted methods, such as
#'     \code{c("CH", "DL")}. \code{NULL} selects all methods.
#' @param link Optional vector of fitted links. \code{NULL} selects all links.
#' @param L Optional numeric vector selecting restriction times already
#'     fitted. No interpolation or refitting is performed.
#' @details The \code{status} column describes the numerical result and takes the values
#' \code{"converged"}, \code{"nonconverged"}, \code{"boundary"}, or
#' \code{"failed"}. Numerical convergence does not by itself establish that
#' regular Wald inference is available. The \code{inference_available} column
#' also requires an interior estimate and finite, positive standard errors.
#' The \code{requires_attention} column is true when regular inference is
#' unavailable or the covariance matrix is not positive semidefinite.
#'
#' The \code{score_norm} is the Euclidean norm of the mean estimating equation.
#' Its convergence criterion is given by \code{score_tolerance}.
#' \code{step_norm} records the norm of the last accepted coefficient step.
#' The three columns beginning with \code{censoring_} report the iterations,
#' score norm, and score tolerance for the fitted censoring model.
#'
#' All four pair-count columns count unordered pairs of individuals.
#' \code{n_pairs} counts all possible pairs. \code{n_resolved} counts comparisons
#' whose outcome can be determined under the hierarchical rule.
#' \code{n_complete_history} counts pairs with both restricted individual
#' histories observed. \code{n_retained} counts the pairs satisfying the selected
#' method's observation requirement. For HT, the mean response term still uses
#' all pairs, including those without an observed win or loss.
#'
#' \code{min_probability} and \code{max_weight} summarize the observation
#' probabilities and weights for retained pairs. \code{n_truncated} counts
#' retained pairs whose probabilities were at or below a positive
#' \code{weight_floor}. A positive floor is an explicit modification of the
#' estimating equation when truncation occurs. The default floor is \code{0.0001},
#' which limits each inverse probability of censoring weight to \code{10000}.
#' Set \code{weight_floor = 0} in \code{pim_control} to disable truncation.
#'
#' Large weights or few retained pairs can help explain unstable estimates.
#' These diagnostics are descriptive and do not supply a universal threshold
#' for accepting a fit. The \code{message} column records the fitting result.
#' Unavailable diagnostics are returned as \code{NA}.
#' @return A data frame with columns \code{L}, \code{method}, \code{link},
#' \code{status}, \code{converged}, \code{boundary}, \code{variance_ok},
#' \code{inference_available}, \code{requires_attention}, \code{iterations},
#' \code{score_norm}, \code{score_tolerance}, \code{step_norm}, \code{n_pairs},
#' \code{n_resolved}, \code{n_complete_history}, \code{n_retained},
#' \code{min_probability}, \code{max_weight}, \code{n_truncated},
#' \code{covariance_positive_semidefinite}, \code{censoring_iterations},
#' \code{censoring_score_norm}, \code{censoring_score_tolerance}, and
#' \code{message}.
#' @seealso \code{\link{gpim}}, \code{\link{pim_control}},
#'   \code{\link{summary.prior_pim}}
#' @export
pim_diagnostics <- function(object, method = NULL, link = NULL, L = NULL) {
  fits <- .pim_selected_fits(object, method, link, L)
  number <- function(x, default = NA_real_) {
    if (length(x) == 1L && is.numeric(x)) as.numeric(x) else default
  }
  rows <- lapply(fits, function(ff) {
    dd <- ff$diagnostics
    converged <- isTRUE(ff$converged)
    boundary <- isTRUE(ff$boundary)
    variance_ok <- isTRUE(ff$variance_ok)
    status <- ff$status
    if (is.null(status) || length(status) != 1L) {
      status <- if (boundary) "boundary" else if (converged) "converged" else
        if (all(is.finite(ff$coefficients))) "nonconverged" else "failed"
    }
    inference_available <- converged && !boundary && variance_ok &&
      length(ff$se) > 0L && all(is.finite(ff$se) & ff$se > 0)
    covariance_psd <- dd$covariance_positive_semidefinite
    if (length(covariance_psd) != 1L || !is.logical(covariance_psd))
      covariance_psd <- NA
    score_tolerance <- number(ff$score_tolerance, number(object$control$tol))
    cox <- ff$censoring
    data.frame(L = ff$L, method = ff$method, link = ff$link,
      status = status, converged = converged, boundary = boundary,
      variance_ok = variance_ok, inference_available = inference_available,
      requires_attention = !inference_available || identical(covariance_psd, FALSE),
      iterations = as.integer(number(ff$iterations)),
      score_norm = number(ff$score_norm), score_tolerance = score_tolerance,
      step_norm = number(ff$step_norm),
      n_pairs = number(dd$n_pairs), n_resolved = number(dd$n_resolved),
      n_complete_history = number(dd$n_complete_history),
      n_retained = number(dd$n_retained),
      min_probability = number(dd$minimum_observation_probability),
      max_weight = number(dd$maximum_weight), n_truncated = number(dd$n_truncated),
      covariance_positive_semidefinite = covariance_psd,
      censoring_iterations = as.integer(number(cox$iterations)),
      censoring_score_norm = number(cox$score_norm),
      censoring_score_tolerance = number(cox$score_tolerance),
      message = if (length(ff$message)) paste(ff$message, collapse = " ") else "",
      stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  ans
}

#' @name prior_pim_methods
#' @aliases summary.prior_pim
#' @aliases print.prior_pim
#' @aliases print.summary.prior_pim
#' @aliases coef.prior_pim
#' @aliases vcov.prior_pim
#' @aliases confint.prior_pim
#' @aliases nobs.prior_pim
#' @aliases plot.prior_pim
#' @aliases trajectory_data
#' @title Summaries and Coefficient Trajectories for Prioritized PIM Fits
#' @description Summarize fitted coefficients, extract model components, and plot coefficient
#' trajectories over the restriction time. Wald inference includes the dependence
#' between comparisons sharing an individual and estimation of the censoring
#' distribution, as computed by the fitted model.
#' @param object,x A fitted object returned by \code{gpim}.
#'     For \code{print.summary.prior_pim}, an object returned by \code{summary}.
#' @param conf.level,level Confidence level strictly between zero and one.
#' @param method Optional vector selecting fitted methods, such as
#'     \code{c("CH", "DL")}. \code{NULL} selects all methods.
#' @param link Optional vector selecting fitted links. \code{NULL} selects all
#'     links.
#' @param L Optional numeric vector selecting restriction times that were
#'     fitted. No interpolation or refitting is performed.
#' @param digits Number of significant digits used when printing.
#' @param parm,parameters Coefficient names or indices. The default selects
#'     every coefficient.
#' @param conf.int Whether to plot pointwise confidence intervals.
#' @param ncol Number of columns in the panel layout. The default uses at most
#'     two columns.
#' @param ... For \code{plot}, arguments passed to
#'     \code{\link[graphics]{plot.default}} when initializing each panel, such as
#'     \code{xlim}, \code{ylim}, or \code{cex.main}. For the other methods,
#'     additional arguments are currently unused.
#' @details \code{summary} reports the estimate, standard error, Wald statistic, two-sided
#' p-value, and pointwise confidence interval for each coefficient, method, link,
#' and restriction time. The printed summary also reports how many fits
#' converged and how many have regular Wald inference. It lists at most six
#' fits requiring attention. Use \code{\link{pim_diagnostics}} to inspect the
#' full diagnostics table. Inference is omitted for fits that did not converge or
#' that lie on the identity model's admissibility boundary. Nonpositive or
#' unavailable standard errors do not produce a Wald interval.
#'
#' Coefficients remain on their fitted scales. Logit coefficients describe log
#' win odds, probit coefficients are on the standard normal quantile scale,
#' and identity coefficients describe net benefit. No exponentiation is
#' performed by these methods.
#'
#' The covariance uses the sparse correlation of pairs sharing an individual.
#' Let \eqn{\mathcal{P}} contain one oriented comparison for each distinct
#' pair, let \eqn{m=n(n-1)/2}, and let \eqn{\zeta_{ab}} denote a pair score
#' after the correction for estimating the censoring distribution. Define
#' \eqn{S_i} as the sum of \eqn{\zeta_{ab}} over pairs containing individual
#' \eqn{i}. The central matrix in the sandwich estimator is
#' \deqn{\widehat B=m^{-2}\left\{\sum_{i=1}^n S_iS_i^T
#' -\sum_{(a,b)\in\mathcal{P}}\zeta_{ab}\zeta_{ab}^T\right\}.}
#' The subtraction counts the contribution of each individual pair once.
#' The covariance is \eqn{\widehat A^{-1}\widehat B\widehat A^{-T}},
#' where \eqn{\widehat A} is the negative derivative of the mean pair score.
#'
#' \code{plot} draws one panel for every selected coefficient and link. Methods
#' are distinguished by color, line type, and point symbol. The horizontal axis
#' uses the numerical restriction times. Lines and shaded pointwise intervals
#' join adjacent fitted times and are a display convention, not an additional
#' smooth model. Lines break at unsuccessful fits. A boundary estimate can be
#' shown with an open point but has no regular Wald interval. A single fitted
#' time is displayed with points and interval bars.
#'
#' The displayed intervals are pointwise, not simultaneous confidence bands
#' over the restriction-time trajectory. The same individuals contribute to
#' the different restriction times. These methods do not treat the fitted
#' trajectories as independent data.
#' @return \code{summary} returns an object of class \code{summary.prior_pim}, with a
#' tidy coefficient table in \code{coefficients} (also available as \code{table})
#' and the table returned by \code{\link{pim_diagnostics}} in
#' \code{diagnostics}. Print methods return their input
#' invisibly.
#'
#' \code{coef} returns a named coefficient vector when one model is selected,
#' or a coefficient-by-model matrix when several models are selected.
#' \code{vcov} returns a covariance matrix for one model, or a named list of
#' covariance matrices for several models. \code{confint} returns a matrix for
#' one selected model, or a data frame with the method, link, restriction time,
#' coefficient name, and interval limits for several models. \code{nobs} returns
#' the number of independent individuals, rather than the number of event rows
#' or pairwise comparisons.
#'
#' \code{trajectory_data} returns the selected tidy coefficient table with
#' intervals at the requested confidence level. \code{plot} returns this table
#' invisibly after drawing the plots.
#' @seealso \code{\link{gpim}}, \code{\link{pim_diagnostics}}
#' @examples
#' # Synthetic HF-ACTION-like data, with time measured in years.
#' data("hfaction", package = "gpim")
#'
#' # Fit all three methods at the default restriction times.
#' fit <- gpim(
#'   HierSurv(time, status) ~
#'     treatment + I(age / 10) + sex + I(lvef / 5),
#'   data = hfaction,
#'   id = "id",
#'   method = c("CH", "HT", "DL"),
#'   link = "logit",
#'   censor_formula = ~ treatment + I(age / 10) + sex
#' )
#'
#' # Summarize DL estimates at every fitted time.
#' summary(fit, method = "DL")
#'
#' # Select one fitted time and request 90 percent confidence intervals.
#' middle_time <- fit$L[2]
#' summary(fit, method = "DL", L = middle_time, conf.level = 0.90)
#'
#' # Extract coefficients, covariance, intervals, and the subject count.
#' coef(fit, method = "DL", L = middle_time)
#' vcov(fit, method = "DL", L = middle_time)
#' confint(fit, method = "DL", L = middle_time)
#' nobs(fit)
#'
#' # Plot every DL coefficient over time with pointwise intervals.
#' plot(fit, method = "DL")
#'
#' # Compare CH, HT, and DL for the treatment coefficient.
#' plot(fit, parameters = "treatmentExercise", conf.int = FALSE)
#'
#' # Obtain the estimates and intervals used in the DL plots.
#' trajectory_data(fit, method = "DL")
#' @rdname prior_pim_methods
#' @order 1
#' @export
summary.prior_pim <- function(object, conf.level = object$conf.level,
                              method = NULL, link = NULL, L = NULL, ...) {
  tab <- .pim_result_table(object, conf.level, method, link, L)
  diagnostics <- pim_diagnostics(object, method, link, L)
  structure(list(call = object$call, formula = object$formula,
    censor_formula = object$censor_formula, n = object$n,
    n_rows = object$n_rows,
    priorities = if (!is.null(object$priorities)) object$priorities else
      object$data$priorities,
    L = unique(tab$L), method = unique(tab$method),
    link = unique(tab$link), conf.level = conf.level,
    coefficients = tab, table = tab,
    diagnostics = diagnostics), class = "summary.prior_pim")
}

.pim_print_summary <- function(x, digits) {
  cat("Prioritized probabilistic index regression\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call, width.cutoff = 80L), collapse = "\n"), "\n")
  }
  cat("\nIndividuals:", x$n, "  Event records:", x$n_rows, "\n")
  if (length(x$priorities)) {
    cat("Priority order:", paste(x$priorities, collapse = " > "), "\n")
  }
  cat("Restriction times (L):", paste(format(x$L, digits = 7, trim = TRUE),
                                   collapse = ", "), "\n")
  cat("Methods:", paste(x$method, collapse = ", "),
      "  Links:", paste(x$link, collapse = ", "), "\n")
  dg <- x$diagnostics
  nfits <- nrow(dg)
  cat("Convergence:", sum(dg$converged), "of", nfits,
      "fits. Regular Wald inference:", sum(dg$inference_available),
      "of", nfits, "fits.\n\n")
  tt <- x$coefficients
  fnum <- function(y) formatC(y, digits = digits, format = "g")
  shown <- data.frame(L = fnum(tt$L), Method = tt$method, Link = tt$link,
    Term = tt$term, Estimate = fnum(tt$estimate), SE = fnum(tt$std.error),
    Lower = fnum(tt$conf.low), Upper = fnum(tt$conf.high),
    z = fnum(tt$statistic), P = format.pval(tt$p.value, digits = digits,
                                         eps = .Machine$double.eps),
    check.names = FALSE, stringsAsFactors = FALSE)
  print(shown, row.names = FALSE, right = TRUE)
  cat("\n", format(100 * x$conf.level, trim = TRUE),
      "% pointwise Wald confidence intervals.\n", sep = "")
  scales <- c(logit = "log win odds", probit = "probit", identity = "net benefit")
  shown_links <- x$link[x$link %in% names(scales)]
  cat("Scales:", paste(paste(shown_links, scales[shown_links], sep = " = "),
                       collapse = "; "), ".\n", sep = " ")
  bad <- dg[dg$requires_attention, , drop = FALSE]
  if (nrow(bad)) {
    cat("\nFits requiring attention:\n")
    nshow <- min(6L, nrow(bad))
    print(bad[seq_len(nshow), c("L", "method", "link", "status", "message")],
          row.names = FALSE, right = FALSE)
    if (nrow(bad) > nshow)
      cat(nrow(bad) - nshow, "further fit(s) require attention.\n")
    cat("Use pim_diagnostics(fit) for the full diagnostics table.\n")
  }
  invisible(x)
}

#' @rdname prior_pim_methods
#' @order 2
#' @export
print.prior_pim <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  .pim_print_summary(summary.prior_pim(x), digits = digits)
  invisible(x)
}

#' @rdname prior_pim_methods
#' @order 3
#' @export
print.summary.prior_pim <- function(x,
    digits = max(3L, getOption("digits") - 3L), ...) {
  .pim_print_summary(x, digits = digits)
  invisible(x)
}

#' @rdname prior_pim_methods
#' @order 4
#' @export
coef.prior_pim <- function(object, method = NULL, link = NULL, L = NULL, ...) {
  fits <- .pim_selected_fits(object, method, link, L)
  if (length(fits) == 1L) return(fits[[1L]]$coefficients)
  terms <- object$term_names
  if (is.null(terms)) terms <- unique(unlist(lapply(fits, function(ff)
    names(ff$coefficients)), use.names = FALSE))
  result <- matrix(NA_real_, nrow = length(terms), ncol = length(fits),
                   dimnames = list(terms, names(fits)))
  for (j in seq_along(fits)) {
    cc <- fits[[j]]$coefficients
    if (is.null(names(cc))) names(cc) <- terms
    result[match(names(cc), terms), j] <- cc
  }
  result
}

#' @rdname prior_pim_methods
#' @order 5
#' @export
vcov.prior_pim <- function(object, method = NULL, link = NULL, L = NULL, ...) {
  fits <- .pim_selected_fits(object, method, link, L)
  result <- lapply(fits, function(ff) ff$vcov)
  if (length(result) == 1L) result[[1L]] else result
}

#' @rdname prior_pim_methods
#' @order 6
#' @export
confint.prior_pim <- function(object, parm, level = 0.95, method = NULL,
                              link = NULL, L = NULL, ...) {
  tab <- .pim_result_table(object, conf.level = level, method = method,
                           link = link, L = L)
  if (!missing(parm)) tab <- .pim_select_parameters(tab, parm)
  keys <- unique(tab[c("L", "method", "link")])
  if (nrow(keys) == 1L) {
    result <- as.matrix(tab[c("conf.low", "conf.high")])
    rownames(result) <- tab$term
    colnames(result) <- paste0(format(100 * c((1 - level) / 2,
                                 (1 + level) / 2), trim = TRUE), "%")
    return(result)
  }
  tab[c("L", "method", "link", "term", "conf.low", "conf.high")]
}

#' @rdname prior_pim_methods
#' @order 7
#' @export
nobs.prior_pim <- function(object, ...) object$n

.pim_true_runs <- function(good) {
  if (!length(good)) return(list())
  runs <- rle(good)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  lapply(which(runs$values), function(j) seq.int(starts[j], ends[j]))
}

#' @rdname prior_pim_methods
#' @order 8
#' @export
plot.prior_pim <- function(x, parameters = NULL, method = NULL, link = NULL,
                           conf.int = TRUE, conf.level = x$conf.level,
                           ncol = NULL, L = NULL, ...) {
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("'conf.int' must be TRUE or FALSE.", call. = FALSE)
  }
  tab <- trajectory_data(x, parameters = parameters, method = method,
                         link = link, L = L, conf.level = conf.level)
  terms <- unique(as.character(tab$term))
  links <- unique(as.character(tab$link))
  panels <- expand.grid(term = terms, link = links, stringsAsFactors = FALSE)
  npanels <- nrow(panels)
  if (is.null(ncol)) ncol <- min(2L, ceiling(sqrt(npanels)))
  if (!is.numeric(ncol) || length(ncol) != 1L || !is.finite(ncol) ||
      ncol < 1 || ncol != as.integer(ncol)) {
    stop("'ncol' must be a positive integer.", call. = FALSE)
  }
  ncol <- min(as.integer(ncol), npanels)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(ceiling(npanels / ncol), ncol),
    mar = c(4.2, 4.5, 2.5, 1.0), oma = c(0, 0, 0, 0), las = 1)
  palette <- c(CH = "#6683A6", HT = "#BC9455", DL = "#4E9B8F")
  line_type <- c(CH = 2, HT = 3, DL = 1)
  point_type <- c(CH = 16, HT = 17, DL = 15)
  axis_label <- c(logit = "Log win odds", probit = "Probit coefficient",
                  identity = "Net benefit")
  dots <- list(...)
  for (panel in seq_len(npanels)) {
    dd <- tab[tab$term == panels$term[panel] &
                tab$link == panels$link[panel], , drop = FALSE]
    methods <- unique(as.character(dd$method))
    valid <- is.finite(dd$estimate) &
      (dd$converged %in% TRUE | dd$boundary %in% TRUE)
    limits <- c(0, dd$estimate[valid])
    if (conf.int) limits <- c(limits, dd$conf.low[valid], dd$conf.high[valid])
    limits <- limits[is.finite(limits)]
    ylim <- if (length(limits)) range(limits) else c(-1, 1)
    span <- diff(ylim)
    if (span == 0) span <- max(0.1, abs(ylim[1]) * 0.1)
    ylim <- ylim + c(-1, 1) * span * 0.08
    xlim <- range(dd$L)
    xspan <- diff(xlim)
    if (xspan == 0) xspan <- max(0.01, abs(xlim[1]) * 0.1)
    xlim <- xlim + c(-1, 1) * xspan * 0.05
    label <- unname(axis_label[panels$link[panel]])
    if (is.na(label)) label <- panels$link[panel]
    args <- utils::modifyList(list(x = NA_real_, y = NA_real_, type = "n",
      xlim = xlim, ylim = ylim, xlab = "Restriction time (L)", ylab = label,
      main = panels$term[panel], bty = "l"), dots)
    do.call(graphics::plot.default, args)
    graphics::abline(h = 0, col = "#B9BFC5", lty = 3, lwd = 0.8)
    for (mm in methods) {
      md <- dd[dd$method == mm, , drop = FALSE]
      md <- md[order(md$L), , drop = FALSE]
      cc <- unname(palette[mm])
      if (is.na(cc)) cc <- "#606870"
      lt <- unname(line_type[mm]); if (is.na(lt)) lt <- 1
      pp <- unname(point_type[mm]); if (is.na(pp)) pp <- 16
      good <- is.finite(md$estimate) &
        (md$converged %in% TRUE | md$boundary %in% TRUE)
      if (conf.int) {
        ci_good <- good & is.finite(md$conf.low) & is.finite(md$conf.high)
        for (idx in .pim_true_runs(ci_good)) {
          if (length(idx) > 1L) {
            graphics::polygon(c(md$L[idx], rev(md$L[idx])),
              c(md$conf.low[idx], rev(md$conf.high[idx])),
              col = grDevices::adjustcolor(cc, alpha.f = 0.15), border = NA)
          } else {
            graphics::segments(md$L[idx], md$conf.low[idx],
              md$L[idx], md$conf.high[idx], col = cc, lwd = 1.2)
            cap <- diff(xlim) * 0.012
            graphics::segments(md$L[idx] - cap,
              c(md$conf.low[idx], md$conf.high[idx]),
              md$L[idx] + cap,
              c(md$conf.low[idx], md$conf.high[idx]), col = cc, lwd = 1.2)
          }
        }
      }
      yy <- md$estimate
      yy[!good] <- NA_real_
      graphics::lines(md$L, yy, col = cc, lty = lt, lwd = 2)
      pch <- rep(pp, nrow(md))
      pch[md$boundary %in% TRUE] <- 1
      graphics::points(md$L, yy, col = cc, pch = pch, cex = 0.8)
    }
    if (!any(valid)) {
      graphics::text(mean(xlim), mean(ylim), "No valid coefficient estimates",
                     cex = 0.85, col = "#606870")
    }
    graphics::legend("topleft", legend = methods, col = palette[methods],
      lty = line_type[methods], pch = point_type[methods], lwd = 2,
      bty = "n", cex = 0.8, inset = 0.015)
  }
  invisible(tab)
}
