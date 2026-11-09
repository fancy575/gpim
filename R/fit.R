# Public fitting interface for censoring-adjusted pairwise estimating equations.

.pim_control <- function(maxit = 100L, tol = 1e-8, weight_floor = 0.0001,
                        keep_pairs = FALSE, step_tol = 1e-10,
                        max_halving = 20L, max_step = 2,
                        cox_maxit = 100L, cox_tol = 1e-9,
                        cox_score_tol = 1e-6, start = NULL,
                        trace = FALSE, progress = FALSE,
                        on_failure = c("warn", "stop", "quiet")) {
    positive <- function(x, name) {
        if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0)
            stop(name, " must be a positive finite number.", call. = FALSE)
    }
    integer_value <- function(x, name, lower) {
        if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
            x < lower || x > .Machine$integer.max || x != floor(x))
            stop(name, " must be an integer at least ", lower, ".", call. = FALSE)
    }
    logical_value <- function(x, name) {
        if (length(x) != 1L || !is.logical(x) || is.na(x))
            stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    integer_value(maxit, "maxit", 1)
    integer_value(max_halving, "max_halving", 0)
    integer_value(cox_maxit, "cox_maxit", 1)
    positive(tol, "tol")
    positive(step_tol, "step_tol")
    positive(max_step, "max_step")
    positive(cox_tol, "cox_tol")
    positive(cox_score_tol, "cox_score_tol")
    if (cox_tol <= .Machine$double.eps^0.75 || cox_tol >= 1)
        stop("cox_tol must exceed .Machine$double.eps^0.75 and be less than one.",
             call. = FALSE)
    if (length(weight_floor) != 1L || !is.numeric(weight_floor) ||
        !is.finite(weight_floor) || weight_floor < 0 || weight_floor >= 1)
        stop("weight_floor must be in [0, 1).", call. = FALSE)
    logical_value(keep_pairs, "keep_pairs")
    logical_value(trace, "trace")
    logical_value(progress, "progress")
    if (!is.null(start) && (!is.numeric(start) || !length(start) ||
                           !is.null(dim(start)) || any(!is.finite(start))))
        stop("start must be NULL or a finite numeric coefficient vector.", call. = FALSE)
    if (!is.null(names(start)) &&
        (anyNA(names(start)) || any(!nzchar(names(start))) || anyDuplicated(names(start))))
        stop("start coefficient names must be nonempty and unique.", call. = FALSE)
    on_failure <- match.arg(on_failure)
    list(maxit = as.integer(maxit), tol = tol, weight_floor = weight_floor,
         keep_pairs = keep_pairs, step_tol = step_tol,
         max_halving = as.integer(max_halving), max_step = max_step,
         cox_maxit = as.integer(cox_maxit), cox_tol = cox_tol,
         cox_score_tol = cox_score_tol, start = start, trace = trace,
         progress = progress, on_failure = on_failure)
}

#' @rdname gpim
#' @order 3
#' @export
pim_control <- function(...) {
    args <- list(...)
    if (length(args) && (is.null(names(args)) || anyNA(names(args)) ||
                        any(!nzchar(names(args))) || anyDuplicated(names(args))))
        stop("Control arguments must have unique, nonempty names.", call. = FALSE)
    unknown <- setdiff(names(args), names(formals(.pim_control)))
    if (length(unknown))
        stop("Unknown control fields: ", paste(unknown, collapse = ", "),
             call. = FALSE)
    do.call(.pim_control, args)
}

.pim_failed_fit <- function(terms, message, control) {
    p <- length(terms)
    list(coefficients = stats::setNames(rep(NA_real_, p), terms),
         vcov = matrix(NA_real_, p, p, dimnames = list(terms, terms)),
         se = stats::setNames(rep(NA_real_, p), terms),
         converged = FALSE, boundary = FALSE, variance_ok = FALSE, message = message,
         status = "failed", iterations = NA_integer_, score_norm = NA_real_,
         score_tolerance = control$tol, step_norm = NA_real_,
         step_tolerance = control$step_tol, start_used = NULL,
         iteration_history = data.frame(iteration = integer(), score_norm = numeric(),
                                        step_norm = numeric(), halvings = integer()),
         diagnostics = list())
}

.pim_fit_table <- function(fits, level) {
    zcrit <- stats::qnorm((1 + level)/2)
    out <- lapply(fits, function(f) {
        b <- f$coefficients
        s <- f$se
        if (is.null(s)) s <- sqrt(diag(f$vcov))
        valid <- isTRUE(f$converged) && !isTRUE(f$boundary) && isTRUE(f$variance_ok)
        if (!valid) s[] <- NA_real_
        z <- b/s
        data.frame(L = f$L, method = f$method, link = f$link,
                   term = names(b), estimate = as.numeric(b),
                   std.error = as.numeric(s), statistic = as.numeric(z),
                   p.value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
                   conf.low = as.numeric(b - zcrit * s),
                   conf.high = as.numeric(b + zcrit * s),
                   converged = isTRUE(f$converged),
                   boundary = isTRUE(f$boundary), row.names = NULL,
                   stringsAsFactors = FALSE)
    })
    ans <- do.call(rbind, out)
    rownames(ans) <- NULL
    ans
}

#' @name gpim
#' @aliases gpim
#' @aliases prior_pim
#' @aliases pim_control
#' @title Fit Probabilistic Index Regression for Prioritized Survival Outcomes
#' @description Fits censoring-adjusted probabilistic index regression at one or more restriction times. The input contains observed long-format event histories and baseline covariates.
#' @param formula An R formula, for example \code{HierSurv(time, status) ~ x1 + b}. A one-sided formula such as \code{~ x1 + b} uses the columns named by \code{time} and \code{status}. Factors, interactions, and standard formula transformations are supported. A free intercept is removed because the predictors are differences between individuals.
#' @param data A data frame with one row per observed event and one terminal row per individual. Covariates used in either formula must be constant within an individual.
#' @param id Name of the subject identifier column.
#' @param time Name of the column containing observed event or terminal times, used with a one-sided formula. With a two-sided formula, specify this column in \code{HierSurv(time, status)}.
#' @param status Name of the event-code column, used with a one-sided formula. Code 0 denotes terminal censoring, code 1 denotes the terminal fatal event, and codes 2 and above denote nonfatal endpoint types. With a two-sided formula, specify this column in \code{HierSurv(time, status)}.
#' @param L Positive numeric restriction times. The default is the 25th, 50th, and 75th sample quantiles of the individual observed death-or-censoring times. Duplicates are removed and times are sorted. Supplied times may not exceed the largest observed terminal time.
#' @param probs Quantile probabilities used only when \code{L = NULL}.
#' @param method One or more of \code{"CH"}, \code{"HT"}, and \code{"DL"}, or \code{"all"}. The default is \code{"DL"}.
#' @param link One or more of \code{"identity"}, \code{"logit"}, and \code{"probit"}. The centered identity link models conditional net benefit.
#' @param censor_formula One-sided formula for a Cox model of censoring. The default uses the same right-hand side as \code{formula}. Use \code{~ 1} for an unconditional Breslow estimate, or specify different censoring predictors, for example \code{~ w}.
#' @param priorities Optional vector of positive event codes ordered from highest to lowest clinical priority. It must begin with 1 for the fatal event and include every positive event code present in \code{data}. The default is increasing numerical order.
#' @param recurrent \code{"first"} uses the first event within each nonfatal type, with a warning if repeated records occur. \code{"error"} rejects repeated occurrences. This package does not implement a recurrence-count comparison.
#' @param na.action \code{"na.fail"}, the default, rejects missing input. \code{"na.omit"} removes an entire individual's history when a required value is missing.
#' @param conf.level Confidence level for the stored coefficient table.
#' @param control A list from \code{pim_control}.
#' @param ... Named numerical and reporting options for \code{pim_control}. See Details for the options and their defaults.
#' @details \code{gpim} is the main fitting function. \code{prior_pim} is an identical alias retained for existing analysis scripts.
#' \code{\link{HierSurv}} constructs the response from the observed time and event-code columns for use on the left-hand side of the formula.
#'
#' The outcome of comparing individuals \eqn{i} and \eqn{j} is one if \eqn{i} wins, one half for a genuine tie, and zero if \eqn{i} loses. A later event time is favorable. Restricted fatal times are compared first. The next endpoint is compared only when a tie is established on every higher-priority component. If censoring prevents this determination, the pairwise comparison remains indeterminate.
#'
#' CH requires both individual histories through death or the restriction time. HT and DL use comparisons resolved at the earlier restricted fatal time. HT weights the centered outcome and leaves the modeled-mean term unweighted. DL weights the entire regression residual. Their untruncated population equations have the same target under the censoring assumptions.
#'
#' For a difference \eqn{z} between outcome model-matrix rows, write \eqn{\eta=z^T\beta}. The means are \eqn{(1+\eta)/2}, \eqn{\operatorname{expit}(\eta)}, and \eqn{\Phi(\eta)} for identity, logit, and probit, respectively. Each satisfies \eqn{h(-\eta)=1-h(\eta)}. Identity coefficients describe conditional net benefit. Estimates at an identity admissibility boundary do not receive regular Wald inference.
#'
#'
#' Censoring is fitted separately at each restriction time with a Cox model and Breslow ties, using follow-up truncated at that time. Censoring exactly at the restriction time does not prevent observation through the restriction time. Survival weights use the left limit of the cumulative censoring hazard. The sandwich variance includes overlapping-pair dependence and the fitted baseline-hazard and censoring-coefficient martingale contributions. It is not an independent-pairs GLM variance.
#'
#' Individuals must be independent. Censoring must be independent of the event history given the specified censoring covariates, with adequate observation probabilities. Covariates are baseline variables. The interval at each selected time is pointwise. Inference treats the evaluation times as fixed and does not add uncertainty for estimating a default sample quantile. Prespecified times are preferable for confirmatory inference.
#'
#' For the iterative links, convergence requires the score norm to meet \code{tol}. Reaching the iteration limit, a failed line search, or a small step with an inadequate score is recorded separately. Inspect \code{pim_diagnostics(fit)} for convergence, available comparisons, weights, and variance flags. Individual fits retain their message and iteration history. Scale continuous covariates to useful units before imposing extremely small tolerances.
#'
#' Numerical options are supplied through \code{control = pim_control(...)}.
#' All arguments to \code{pim_control} must be named. The available options are
#' \describe{
#' \item{\code{maxit}}{Maximum number of score iterations for nonlinear links. The identity equation uses a closed-form solution and, if required, an internal constrained boundary calculation. Default \code{100L}.}
#' \item{\code{tol}}{Required Euclidean norm of the mean pair score for convergence. Default \code{1e-8}.}
#' \item{\code{weight_floor}}{Truncation on the estimated pair observation probability used in the inverse probability of censoring weight. The weight is \eqn{1/\max(\widehat\pi_{ij}, \epsilon)}, where \eqn{\epsilon} is \code{weight_floor}. The default \code{0.0001} limits the weight to \code{10000}. Set \code{weight_floor = 0} to disable truncation. A positive floor changes the estimating equation when the estimated probability is below the floor. The censoring adjustment in the variance calculation differentiates the weight only when the probability exceeds the floor.}
#' \item{\code{keep_pairs}}{Whether to retain pair-level quantities in the fitted objects. This can substantially increase memory use. Default \code{FALSE}.}
#' \item{\code{step_tol}}{Positive tolerance for detecting a small relative accepted step. A small step without meeting \code{tol} is treated as stagnation, not successful convergence. Default \code{1e-10}.}
#' \item{\code{max_halving}}{Nonnegative integer limiting the number of step halvings in a line search. Zero permits only the initial proposed step. Default \code{20L}.}
#' \item{\code{max_step}}{Positive bound on the Euclidean norm of a proposed Newton step. Default \code{2}.}
#' \item{\code{cox_maxit}}{Positive integer iteration limit for the censoring Cox fit. Default \code{100L}.}
#' \item{\code{cox_tol}}{The \code{eps} argument passed to \code{survival::coxph.control}. It governs the change in partial log likelihood and must exceed \code{.Machine$double.eps^0.75}, the default Cholesky tolerance, and be less than one. This is separate from the PIM score tolerance. Default \code{1e-9}.}
#' \item{\code{cox_score_tol}}{Positive bound on the fitted censoring Cox score norm divided by the number of censoring events. A fit that fails this additional score check is not used to construct weights. Default \code{1e-6}.}
#' \item{\code{start}}{Optional finite coefficient vector for every requested nonlinear fit. An unnamed vector must follow the model-matrix columns. A named vector must contain exactly those names and is reordered to match them. No intercept is included. Identity fits ignore starting values. Default \code{NULL}.}
#' \item{\code{trace}}{Logical. Print the individual PIM iterations. Default \code{FALSE}.}
#' \item{\code{progress}}{Logical. Print progress through the method, link, and time combinations. Default \code{FALSE}.}
#' \item{\code{on_failure}}{\code{"warn"} retains problematic fits with a warning. \code{"quiet"} retains them without a fitting warning. \code{"stop"} stops when a requested fit needs attention, including nonconvergence, a boundary estimate, or invalid or suspect variance. Input errors always stop. Default \code{"warn"}.}
#' }
#' For example, \code{control = pim_control(trace = TRUE)} prints iteration
#' details, and \code{control = pim_control(progress = TRUE)} reports progress
#' through the requested fits.
#' @return An object of class \code{prior_pim}. Its main fields are \code{table}, a tidy coefficient and inference table; \code{fits}, the individual fitted models with covariance matrices and diagnostic information; \code{data}, the validated subject and event data; \code{L}; \code{priorities}; the formulas and call; and \code{messages}. Restriction times are stored in \code{L}, both in each coefficient or diagnostic table and in the fitted object. Numerical failures remain in the object with unavailable inference and an explanatory message. \code{pim_control} returns a control list.
#' @references Fang, X., Cao, Z., Mao, L., and Li, F. (2026). Probabilistic index models for hierarchical survival outcomes. Manuscript under review.
#' @seealso \code{\link{HierSurv}}, \code{\link{hfaction}}, \code{\link{pim_diagnostics}}, \code{\link{summary.prior_pim}}, \code{\link{plot.prior_pim}}
#' @examples
#' # Synthetic data, with time measured in years.
#' # Status 0 is censoring, 1 is death, and 2 is hospitalization.
#' data("hfaction", package = "gpim")
#'
#' # Fit DL with a logit link at the three default restriction times.
#' # Age is measured per 10 years and LVEF per 5 percentage points.
#' fit <- gpim(
#'   HierSurv(time, status) ~
#'     treatment + I(age / 10) + sex + I(lvef / 5),
#'   data = hfaction,
#'   id = "id",
#'   method = "DL",
#'   link = "logit",
#'   censor_formula = ~ treatment + I(age / 10) + sex
#' )
#'
#' # Default times are the 25th, 50th, and 75th percentiles of the
#' # observed death-or-censoring time, counting each individual once.
#' # To use specified times, add L = c(1, 2, 3) to the call above.
#' fit$L
#'
#' # Estimates, standard errors, confidence intervals, and p-values.
#' summary(fit)
#'
#' # One trajectory panel per coefficient, with pointwise intervals.
#' plot(fit)
#' @rdname gpim
#' @order 2
#' @export
prior_pim <- function(formula, data, id = "id", time = "time",
                      status = "status", L = NULL,
                      probs = c(0.25, 0.5, 0.75), method = "DL",
                      link = "logit", censor_formula = NULL,
                      priorities = NULL, recurrent = "first",
                      na.action = "na.fail", conf.level = 0.95,
                      control = pim_control()) {
    cl <- match.call()
    if (!inherits(formula, "formula"))
        stop("formula must be an R formula.", call. = FALSE)
    if (length(conf.level) != 1L || !is.numeric(conf.level) ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
        stop("conf.level must be strictly between zero and one.", call. = FALSE)
    if (!is.list(control)) stop("control must be a list from pim_control().")
    if (length(control) && (is.null(names(control)) || any(!nzchar(names(control))) ||
                            anyNA(names(control)) || anyDuplicated(names(control))))
        stop("control must be a list with unique, nonempty field names.", call. = FALSE)
    unknown <- setdiff(names(control), names(pim_control()))
    if (length(unknown)) stop("Unknown control fields: ", paste(unknown, collapse = ", "))
    control <- do.call(pim_control, utils::modifyList(pim_control(), control))
    if (!is.character(method) || !length(method) || anyNA(method))
        stop("method must be CH, HT, DL, or a vector of these.")
    method <- unique(toupper(method))
    if (identical(method, "ALL")) method <- c("CH", "HT", "DL")
    if (any(!method %in% c("CH", "HT", "DL")))
        stop("Unknown method. Choose CH, HT, DL, or all.", call. = FALSE)
    if (!is.character(link) || !length(link) || anyNA(link))
        stop("link must contain one or more supported link names.", call. = FALSE)
    link <- unique(tolower(link))
    if (any(!link %in% c("identity", "logit", "probit")))
        stop("Unknown link. Choose identity, logit, or probit.",
             call. = FALSE)
    dat <- .pim_prepare_data(formula, data, id = id, time = time,
                             status = status, priorities = priorities,
                             censor_formula = censor_formula,
                             recurrent = recurrent, na.action = na.action)
    if (!is.null(control$start)) {
        if (length(control$start) != ncol(dat$X))
            stop("start must contain one value per model coefficient (",
                 ncol(dat$X), ").", call. = FALSE)
        if (!is.null(names(control$start))) {
            if (!setequal(names(control$start), colnames(dat$X)))
                stop("Named start values must match the model coefficient names: ",
                     paste(colnames(dat$X), collapse = ", "), ".", call. = FALSE)
            control$start <- control$start[colnames(dat$X)]
        }
        names(control$start) <- colnames(dat$X)
    }
    automatic_L <- is.null(L)
    requested_L <- L
    if (automatic_L) {
        if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
            any(probs <= 0 | probs >= 1))
            stop("probs must contain numbers strictly between zero and one.")
        L <- as.numeric(stats::quantile(dat$u, probs = probs, type = 7,
                                            names = FALSE))
        if (any(L <= 0)) {
            warning("Nonpositive default quantiles were omitted. Supply L explicitly if needed.",
                    call. = FALSE)
            L <- L[L > 0]
        }
    }
    if (!is.numeric(L) || !length(L) || any(!is.finite(L)) || any(L <= 0))
        stop("L must be a nonempty vector of positive finite restriction times.",
             call. = FALSE)
    if (any(L > max(dat$u)))
        stop("A restriction time exceeds the largest observed death-or-censoring time.",
             call. = FALSE)
    L <- sort(unique(L))
    fits <- list()
    censoring <- list()
    messages <- character()
    total_fits <- length(L) * length(link) * length(method)
    fit_number <- 0L
    for (L_i in L) {
        if (control$progress)
            message("Preparing comparisons and censoring weights at time ",
                    format(L_i, digits = 7), ".")
        cache <- tryCatch(.pim_prepare_horizon(dat, L_i, control),
                          error = identity)
        if (!inherits(cache, "error")) {
            censoring[[format(L_i, digits = 16)]] <- cache$censoring
            if (is.null(censoring[[format(L_i, digits = 16)]]))
                censoring[[format(L_i, digits = 16)]] <- cache$cox
        }
        for (lnk in link) for (m in method) {
            fit_number <- fit_number + 1L
            if (control$progress || control$trace)
                message("Fit ", fit_number, "/", total_fits, ": ", m, ", ",
                        lnk, ", L = ", format(L_i, digits = 7), ".")
            f <- if (inherits(cache, "error")) {
                .pim_failed_fit(colnames(dat$X), conditionMessage(cache), control)
            } else {
                tryCatch(.pim_fit_cached(cache, m, lnk, control), error = function(e)
                    .pim_failed_fit(colnames(dat$X), conditionMessage(e), control))
            }
            f$L <- L_i
            f$method <- m
            f$link <- lnk
            if (is.null(names(f$coefficients))) names(f$coefficients) <- colnames(dat$X)
            if (is.null(f$se)) f$se <- sqrt(diag(f$vcov))
            if (is.null(f$boundary)) f$boundary <- FALSE
            key <- paste(m, lnk, format(L_i, digits = 16), sep = ":")
            fits[[key]] <- f
            if (!isTRUE(f$converged) || isTRUE(f$boundary) || !isTRUE(f$variance_ok) ||
                any(!is.finite(f$se)) ||
                identical(f$diagnostics$covariance_positive_semidefinite, FALSE)) {
                messages <- c(messages, paste0(key, " -- ", f$message))
                if (control$on_failure == "stop")
                    stop("Fit ", key, " requires attention. ", f$message, call. = FALSE)
            }
            if (control$progress)
                message("  ", f$status, ". ", f$message)
        }
    }
    ans <- list(call = cl, formula = formula, censor_formula = dat$censor_formula,
                n = nrow(dat$X), n_rows = dat$n_rows, L = L, probs = probs,
                automatic_L = automatic_L, requested_L = requested_L,
                method = method, link = link, term_names = colnames(dat$X),
                priorities = dat$priorities, fits = fits, conf.level = conf.level,
                table = .pim_fit_table(fits, conf.level), control = control,
                data = dat, censoring = censoring, diagnostics = lapply(fits, `[[`, "diagnostics"),
                messages = messages, version = "0.0.5")
    class(ans) <- "prior_pim"
    if (length(messages) && control$on_failure == "warn")
        warning(length(messages), " fit(s) require attention. ",
                "Inspect pim_diagnostics(fit) or fit$messages for details.", call. = FALSE)
    ans
}

# Main fitting name. Keep prior_pim available for existing analyses.
#' @rdname gpim
#' @order 1
#' @export gpim
gpim <- prior_pim

