# Joint event models for conditional tie weighting.

#' Specify a Joint Model for Conditional Tie Weighting
#'
#' @description Specifies the joint model used to estimate conditional fatal-tie
#' probabilities for CTW. The built-in models describe one fatal and one
#' nonfatal event, with no nonfatal event after death.
#' @param copula One of \code{"gumbel"}, \code{"clayton"}, \code{"frank"},
#'   \code{"plackett"}, or \code{"independence"}.
#' @param margins A character vector specifying the death and nonfatal margin
#'   functions, in that order. Choices are \code{"weibull"},
#'   \code{"lognormal"}, \code{"exponential"}, and \code{"cox"}.
#'   A single choice is used for both margins. A named vector must use
#'   \code{death} and \code{nonfatal}.
#' @param death_formula,nonfatal_formula Optional one-sided formulas for the
#'   covariates in the two margin functions. Each defaults to the PIM formula's
#'   right-hand side. Use \code{~1} for a model without covariates.
#' @param ... Named fitting controls. See Details.
#' @details
#' Gumbel, Clayton, and Frank describe positive dependence. Plackett allows
#' positive or negative dependence, and independence has no dependence
#' parameter. Weibull, exponential, and Cox margins use proportional hazards.
#' The lognormal margin is an accelerated failure time model. The second margin
#' is a function used to specify the joint distribution on the region where the
#' nonfatal event precedes death. It is not the marginal survival function of
#' the actual nonfatal event, which may never occur before death.
#'
#' All parameters are fitted jointly using the likelihood contributions for
#' both observed events, an observed nonfatal event followed by censoring,
#' death without a preceding nonfatal event, and censoring without either
#' event. A Cox margin has a separate baseline hazard increment at each
#' distinct observed event time for that margin. Both sets of increments are
#' included in the joint fit and in its information matrix. A nonfatal Cox
#' model that treats death as independent censoring is not used as the final
#' joint model.
#'
#' Covariates are centered for the joint fit. Baseline parameters refer to
#' those centered covariates. Regression coefficients retain their original
#' units. Parametric margins require strictly positive event times. Zero
#' censoring times contribute the corresponding boundary probabilities.
#'
#' The fitting controls in \code{...} are \code{maxit = 800L},
#' \code{reltol = 1e-9}, \code{score_tol = 1e-5},
#' \code{hessian_step = 1e-3}, \code{starts = 3L}, \code{start = NULL},
#' and \code{trace = FALSE}. \code{starts} controls the number of deterministic
#' initial dependence values. \code{start}, when supplied, is a complete
#' parameter vector in the order reported by \code{coef(pim_joint_fit(fit))}.
#' \code{hessian_step} is the relative step for differentiating the native
#' score to obtain the full information matrix. No generalized inverse or
#' ridge term is used in the variance calculation.
#'
#' The fitted model is shared across the requested restriction times and links.
#' Conditional tie probabilities use only each subject's history through the
#' relevant restriction time. Inference is conditional on the selected model
#' family and formulas. The constructor does not select a model using the
#' outcomes and then ignore that selection.
#' @return An object of class \code{pim_joint_model}, supplied as
#'   \code{joint_model} in \code{gpim}.
#' @references
#' Chen, Y.-H. (2012). Maximum likelihood analysis of semicompeting risks data
#' with semiparametric regression models. \emph{Lifetime Data Analysis},
#' \bold{18}, 36--57. \doi{10.1007/s10985-011-9202-4}.
#'
#' Arachchige, S. J., Chen, X., and Zhou, Q. M. (2025). Two-stage pseudo maximum
#' likelihood estimation of semiparametric copula-based regression models for
#' semi-competing risks data. \emph{Lifetime Data Analysis}, \bold{31}, 52--75.
#' \doi{10.1007/s10985-024-09640-z}.
#' @seealso \code{\link{gpim}}, \code{\link{pim_joint_fit}}
#' @examples
#' pim_joint_model(copula = "gumbel", margins = "weibull")
#' pim_joint_model(copula = "clayton", margins = c("cox", "cox"),
#'                 death_formula = ~ treatment + age,
#'                 nonfatal_formula = ~ treatment + age)
#' @export
pim_joint_model <- function(copula = "gumbel",
                            margins = c("weibull", "weibull"),
                            death_formula = NULL, nonfatal_formula = NULL, ...) {
    choices <- c("gumbel", "clayton", "frank", "plackett", "independence")
    if (!is.character(copula) || length(copula) != 1L || is.na(copula))
        stop("'copula' must name one supported copula.", call. = FALSE)
    copula <- match.arg(tolower(copula), choices)
    margin_choices <- c("weibull", "lognormal", "exponential", "cox")
    if (!is.character(margins) || !length(margins) %in% c(1L, 2L) || anyNA(margins))
        stop("'margins' must contain one or two margin names.", call. = FALSE)
    if (length(margins) == 1L) margins <- rep(margins, 2L)
    if (!is.null(names(margins))) {
        if (!setequal(names(margins), c("death", "nonfatal")) ||
            anyDuplicated(names(margins)))
            stop("Named 'margins' must use 'death' and 'nonfatal'.", call. = FALSE)
        margins <- margins[c("death", "nonfatal")]
    }
    margins <- tolower(unname(margins))
    if (any(!margins %in% margin_choices))
        stop("Supported margins are weibull, lognormal, exponential, and cox.",
             call. = FALSE)
    names(margins) <- c("death", "nonfatal")
    for (nm in c("death_formula", "nonfatal_formula")) {
        f <- get(nm, inherits = FALSE)
        if (!is.null(f) && (!inherits(f, "formula") || length(f) != 2L))
            stop("'", nm, "' must be NULL or a one-sided formula.", call. = FALSE)
    }
    defaults <- list(maxit = 800L, reltol = 1e-9, score_tol = 1e-5,
                     hessian_step = 1e-3, starts = 3L, start = NULL,
                     trace = FALSE)
    args <- list(...)
    if (length(args) && (is.null(names(args)) || anyNA(names(args)) ||
        any(!nzchar(names(args))) || anyDuplicated(names(args))))
        stop("Joint-model controls must have unique, nonempty names.", call. = FALSE)
    unknown <- setdiff(names(args), names(defaults))
    if (length(unknown))
        stop("Unknown joint-model control(s) ", paste(unknown, collapse = ", "),
             ".", call. = FALSE)
    ctl <- utils::modifyList(defaults, args, keep.null = TRUE)
    for (nm in c("maxit", "starts")) {
        x <- ctl[[nm]]
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
            x != floor(x) || x > .Machine$integer.max)
            stop("Joint-model ", nm, " must be a positive integer.", call. = FALSE)
        ctl[[nm]] <- as.integer(x)
    }
    for (nm in c("reltol", "score_tol", "hessian_step")) {
        x <- ctl[[nm]]
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
            stop("Joint-model ", nm, " must lie strictly between zero and one.",
                 call. = FALSE)
    }
    if (!is.logical(ctl$trace) || length(ctl$trace) != 1L || is.na(ctl$trace))
        stop("Joint-model trace must be TRUE or FALSE.", call. = FALSE)
    if (!is.null(ctl$start) && (!is.numeric(ctl$start) || !length(ctl$start) ||
        !is.null(dim(ctl$start)) || any(!is.finite(ctl$start))))
        stop("Joint-model start must be NULL or a finite numeric vector.", call. = FALSE)
    structure(list(copula = copula, margins = margins,
                   death_formula = death_formula,
                   nonfatal_formula = nonfatal_formula, control = ctl),
              class = "pim_joint_model")
}

#' @export
print.pim_joint_model <- function(x, ...) {
    cat("Joint model for conditional tie weighting\n")
    cat("Copula ", x$copula, "\n", sep = "")
    cat("Death margin ", x$margins[[1L]], "\n", sep = "")
    cat("Nonfatal margin ", x$margins[[2L]], "\n", sep = "")
    invisible(x)
}

.pim_joint_design <- function(dat, formula, label) {
    if (is.null(formula)) {
        x <- dat$X
        terms <- dat$terms
        contrasts <- dat$contrasts
        xlevels <- dat$xlevels
    } else {
        tt <- .pim_terms(formula, dat$subjects,
                         c(dat$id_column, dat$response_variables), label)
        mf <- .pim_make_frame(tt, dat$subjects, label)
        if (any(!stats::complete.cases(mf)))
            stop(label, " covariates contain missing values after sample selection.",
                 call. = FALSE)
        mm <- .pim_model_matrix(mf, label, allow_empty = TRUE)
        x <- mm$matrix
        terms <- mm$terms
        contrasts <- mm$contrasts
        xlevels <- mm$xlevels
    }
    x <- as.matrix(x)
    center <- if (ncol(x)) colMeans(x) else numeric()
    if (ncol(x)) x <- sweep(x, 2L, center, "-")
    storage.mode(x) <- "double"
    list(matrix = x, center = center, terms = terms, contrasts = contrasts,
         xlevels = xlevels)
}

.pim_joint_margin_names <- function(prefix, margin, grid, xnames) {
    baseline <- switch(margin,
        weibull = c("log_rate", "log_shape"),
        exponential = "log_rate", lognormal = c("location", "log_sd"),
        cox = paste0("log_jump", seq_along(grid)))
    regression <- if (length(xnames)) paste0("beta.", xnames) else character()
    paste0(prefix, ".", c(baseline, regression))
}

.pim_joint_margin_start <- function(margin, time, status, X, grid, control) {
    p <- ncol(X)
    if (!any(status))
        stop("A joint margin has no observed events and cannot be estimated.",
             call. = FALSE)
    if (margin != "cox" && any(time[status] <= 0))
        stop("Parametric joint margins require strictly positive event times. ",
             "Use an appropriate event-time specification or a Cox margin.",
             call. = FALSE)
    if (margin == "cox") {
        initial <- tryCatch(.pim_cox(time, status, X, control),
                            error = function(e) NULL)
        beta <- if (is.null(initial)) rep(0, p) else unname(initial$coefficients)
        risk <- if (p) exp(drop(X %*% beta)) else rep(1, length(time))
        number <- tabulate(match(time[status], grid), nbins = length(grid))
        jumps <- vapply(seq_along(grid), function(k)
            number[k] / sum(risk[time >= grid[k]]), numeric(1L))
        return(c(log(jumps), beta))
    }
    # These independent-margin fits provide initial values only. The final
    # estimate always maximizes the semi-competing joint likelihood.
    positive <- time > 0
    yy <- survival::Surv(time[positive], status[positive])
    xx <- X[positive, , drop = FALSE]
    formula <- if (p) yy ~ xx else yy ~ 1
    initial <- tryCatch(suppressWarnings(survival::survreg(
        formula, dist = margin,
        control = survival::survreg.control(maxiter = 100L))),
        error = function(e) NULL)
    if (!is.null(initial) && all(is.finite(stats::coef(initial))) &&
        is.finite(initial$scale) && initial$scale > 0) {
        cf <- unname(stats::coef(initial))
        if (margin == "lognormal")
            return(c(cf[1L], log(initial$scale), cf[-1L]))
        if (margin == "exponential") return(c(-cf[1L], -cf[-1L]))
        shape <- 1 / initial$scale
        return(c(-cf[1L] * shape, log(shape), -cf[-1L] * shape))
    }
    if (margin == "lognormal") {
        z <- log(time[positive])
        spread <- stats::sd(z)
        if (!is.finite(spread) || spread < 0.1) spread <- 1
        return(c(stats::median(z), log(spread), rep(0, p)))
    }
    rate <- sum(status) / sum(time)
    if (margin == "exponential") return(c(log(rate), rep(0, p)))
    c(log(rate), 0, rep(0, p))
}

.pim_joint_information <- function(par, data, spec, step) {
    n <- length(data$u)
    d <- length(par)
    ans <- matrix(0, d, d)
    for (j in seq_len(d)) {
        h <- step * (1 + abs(par[j]))
        plus <- minus <- plus2 <- minus2 <- par
        plus[j] <- plus[j] + h
        minus[j] <- minus[j] - h
        plus2[j] <- plus2[j] + 2 * h
        minus2[j] <- minus2[j] - 2 * h
        gp <- gpim_joint_cpp(plus, data, spec, FALSE)$gradient
        gm <- gpim_joint_cpp(minus, data, spec, FALSE)$gradient
        gp2 <- gpim_joint_cpp(plus2, data, spec, FALSE)$gradient
        gm2 <- gpim_joint_cpp(minus2, data, spec, FALSE)$gradient
        ans[, j] <- -(gm2 - 8 * gm + 8 * gp - gp2) / (12 * h * n)
    }
    asymmetry <- max(abs(ans - t(ans))) / max(1, max(abs(ans)))
    ans <- (ans + t(ans)) / 2
    dimnames(ans) <- list(names(par), names(par))
    attr(ans, "relative_asymmetry") <- asymmetry
    ans
}

.pim_fit_joint <- function(dat, spec, control) {
    if (!inherits(spec, "pim_joint_model"))
        stop("'joint_model' must be an object from pim_joint_model().", call. = FALSE)
    if (ncol(dat$event_times) != 1L)
        stop("Built-in copula models require one fatal and one nonfatal component.",
             call. = FALSE)
    n <- length(dat$u)
    death <- as.logical(dat$death)
    nonfatal <- is.finite(dat$event_times[, 1L])
    v <- ifelse(nonfatal, dat$event_times[, 1L], dat$u)
    if (any(nonfatal & death & v >= dat$u))
        stop("The joint model requires an observed nonfatal event to precede death.",
             call. = FALSE)
    if (any(nonfatal & v > dat$u))
        stop("A nonfatal event cannot follow the end of observation.", call. = FALSE)
    if (!any(death) || !any(nonfatal))
        stop("The built-in joint model requires observed fatal and nonfatal events.",
             call. = FALSE)
    desD <- .pim_joint_design(dat, spec$death_formula, "Death margin")
    desT <- .pim_joint_design(dat, spec$nonfatal_formula, "Nonfatal margin")
    gridD <- if (spec$margins[1L] == "cox") sort(unique(dat$u[death])) else numeric()
    gridT <- if (spec$margins[2L] == "cox") sort(unique(v[nonfatal])) else numeric()
    native_data <- list(u = as.numeric(dat$u), death = as.integer(death),
                        nonfatal = as.integer(nonfatal), v = as.numeric(v),
                        XD = desD$matrix, XT = desT$matrix)
    native_spec <- list(copula = spec$copula,
                        marginD = unname(spec$margins[1L]),
                        marginT = unname(spec$margins[2L]),
                        gridD = gridD, gridT = gridT)
    dstart <- .pim_joint_margin_start(spec$margins[1L], dat$u, death,
                                     desD$matrix, gridD, control)
    tstart <- .pim_joint_margin_start(spec$margins[2L], v, nonfatal,
                                     desT$matrix, gridT, control)
    par <- c(dstart, tstart)
    par_names <- c(.pim_joint_margin_names("death", spec$margins[1L],
                                          gridD, colnames(desD$matrix)),
                   .pim_joint_margin_names("nonfatal", spec$margins[2L],
                                          gridT, colnames(desT$matrix)))
    dependent <- spec$copula != "independence"
    dependence_starts <- switch(spec$copula, gumbel = log(c(1, 0.25, 3)),
        clayton = log(c(2, 0.5, 6)), frank = log(c(5, 2, 10)),
        plackett = log(c(3, 1, 10)), independence = numeric())
    if (dependent) {
        par <- c(par, dependence_starts[1L])
        par_names <- c(par_names, "copula.log_parameter")
    }
    names(par) <- par_names
    ctl <- spec$control
    if (!is.null(ctl$start)) {
        if (length(ctl$start) != length(par))
            stop("Joint-model start has the wrong number of parameters.", call. = FALSE)
        if (!is.null(names(ctl$start))) {
            if (anyDuplicated(names(ctl$start)) || !setequal(names(ctl$start), par_names))
                stop("Named joint-model starting values must match all parameter names.",
                     call. = FALSE)
            par[] <- ctl$start[par_names]
        } else par[] <- ctl$start
    }
    cache <- new.env(parent = emptyenv())
    cache$par <- NULL
    evaluate <- function(x) {
        if (identical(as.numeric(x), cache$par)) return(cache$value)
        z <- tryCatch(gpim_joint_cpp(x, native_data, native_spec, FALSE),
                      error = function(e) NULL)
        if (is.null(z) || !is.finite(z$loglik) ||
            length(z$gradient) != length(x) || any(!is.finite(z$gradient)))
            z <- list(loglik = -Inf, gradient = rep(NA_real_, length(x)))
        cache$par <- as.numeric(x)
        cache$value <- z
        z
    }
    objective <- function(x) -evaluate(x)$loglik / n
    gradient <- function(x) -evaluate(x)$gradient / n
    candidates <- vector("list", if (dependent) ctl$starts else 1L)
    for (k in seq_along(candidates)) {
        start <- par
        if (dependent && k > 1L)
            start[length(start)] <- dependence_starts[
                1L + (k - 1L) %% length(dependence_starts)]
        candidates[[k]] <- tryCatch(stats::nlminb(start, objective, gradient,
            control = list(iter.max = ctl$maxit, eval.max = 2L * ctl$maxit,
                           rel.tol = ctl$reltol, x.tol = 1e-9,
                           trace = if (isTRUE(ctl$trace) || isTRUE(control$trace)) 1L else 0L)),
            error = function(e) list(objective = Inf, convergence = 1L,
                                     message = conditionMessage(e)))
    }
    values <- vapply(candidates, function(x) x$objective, numeric(1L))
    if (!any(is.finite(values)))
        stop("The joint likelihood could not be optimized from the supplied initial values.",
             call. = FALSE)
    best <- candidates[[which.min(values)]]
    par <- stats::setNames(as.numeric(best$par), par_names)
    evaluation <- gpim_joint_cpp(par, native_data, native_spec, TRUE)
    information <- .pim_joint_information(par, native_data, native_spec,
                                         ctl$hessian_step)
    if (any(!is.finite(information)))
        stop("The full joint-model information matrix is not finite.", call. = FALSE)
    factor <- tryCatch(chol(information), error = function(e) NULL)
    if (is.null(factor))
        stop("The full joint-model information matrix is not positive definite. ",
             "Check the selected margin and dependence models.", call. = FALSE)
    scale <- pmax(sqrt(diag(information)), 1e-8)
    score_norm <- max(abs(evaluation$gradient / n) / scale)
    # Refine a nearly converged optimizer solution using the same full score.
    # This is Newton iteration, not a modification of the information matrix.
    refinements <- 0L
    while (score_norm > ctl$score_tol && refinements < 4L) {
        direction <- backsolve(factor, forwardsolve(t(factor),
                              evaluation$gradient / n))
        accepted <- FALSE
        for (halving in 0:12) {
            trial <- par + direction / 2^halving
            z <- evaluate(trial)
            if (is.finite(z$loglik) && z$loglik >= evaluation$loglik - 1e-10) {
                par <- stats::setNames(trial, par_names)
                accepted <- TRUE
                break
            }
        }
        if (!accepted) break
        refinements <- refinements + 1L
        evaluation <- gpim_joint_cpp(par, native_data, native_spec, TRUE)
        information <- .pim_joint_information(par, native_data, native_spec,
                                             ctl$hessian_step)
        factor <- tryCatch(chol(information), error = function(e) NULL)
        if (is.null(factor)) break
        scale <- pmax(sqrt(diag(information)), 1e-8)
        score_norm <- max(abs(evaluation$gradient / n) / scale)
    }
    if (is.null(factor) || !is.finite(score_norm) || score_norm > ctl$score_tol)
        stop("The joint likelihood score did not converge sufficiently for inference. ",
             "Increase maxit or review the joint-model specification.", call. = FALSE)
    if (dependent && spec$copula %in% c("gumbel", "clayton", "frank") &&
        exp(par[length(par)]) < 1e-6)
        stop("The estimated dependence is at the independence boundary. ",
             "Use copula = 'independence' for an independence model.", call. = FALSE)
    scores <- as.matrix(evaluation$scores)
    dimnames(scores) <- list(as.character(dat$id), par_names)
    dependence <- if (!dependent) NULL else
        exp(par[length(par)]) + as.numeric(spec$copula == "gumbel")
    structure(list(coefficients = par, par = par, parameter_names = par_names,
                   scores = scores, information = information,
                   loglik = evaluation$loglik, logLik = evaluation$loglik,
                   loglik_i = evaluation$loglik_i, gradient = evaluation$gradient,
                   convergence = best$convergence, converged = TRUE,
                   message = best$message, score_norm = score_norm,
                   score_tolerance = ctl$score_tol,
                   iterations = best$iterations, refinements = refinements,
                   n = n, id = as.character(dat$id), dependence = dependence,
                   model = spec, native_data = native_data, native_spec = native_spec,
                   designD = desD, designT = desT,
                   centerD = desD$center, centerT = desT$center,
                   baseline_convention = "right-continuous cumulative hazards",
                   variance = "full joint likelihood subject influence"),
              class = "pim_joint_fit")
}

.pim_joint_predict <- function(fit, dat, L) {
    if (!inherits(fit, "pim_joint_fit") || !isTRUE(fit$converged))
        stop("Conditional tie probabilities require a fitted joint model.", call. = FALSE)
    if (!identical(as.character(dat$id), fit$id))
        stop("The joint fit and PIM must use the same individuals in the same order.",
             call. = FALSE)
    if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L <= 0)
        stop("L must be a positive finite number.", call. = FALSE)
    z <- gpim_joint_q_cpp(fit$par, fit$native_data, fit$native_spec, L)
    probability <- matrix(as.numeric(z$q), ncol = 1L,
                          dimnames = list(fit$id, "nonfatal"))
    jacobian <- as.matrix(z$dq)
    dimnames(jacobian) <- list(fit$id, fit$parameter_names)
    if (any(!is.finite(probability)) || any(probability < 0 | probability > 1) ||
        any(!is.finite(jacobian)))
        stop("The fitted joint model returned invalid conditional probabilities or derivatives.",
             call. = FALSE)
    list(probability = probability, jacobian = list(jacobian),
         log_probability = matrix(as.numeric(z$logq), ncol = 1L))
}

#' Extract the Fitted CTW Joint Model
#'
#' @description Extracts the joint model fitted once for the CTW analysis.
#' @param object A fitted \code{prior_pim} object containing CTW, or an object
#'   returned by \code{pim_joint_fit} for the extraction methods.
#' @param x A fitted joint-model object for printing.
#' @param ... Further arguments passed to or from methods.
#' @details The returned object contains the fitted parameter vector, the
#' subject scores, and the negative mean Hessian for all parameters. With Cox
#' margins, these include every estimated baseline hazard increment.
#' \code{vcov} returns their centered subject-influence covariance. The PIM
#' covariance also includes the pair score, the fitted censoring model, and
#' cross-covariances between all contributions.
#'
#' \code{coef} returns the native parameterization, with log rates, log shapes,
#' log standard deviations, log Cox increments, and the transformed dependence
#' parameter identified by their names. For Gumbel the dependence parameter is
#' \eqn{1+\exp(a)}. For Clayton, Frank, and Plackett it is \eqn{\exp(a)}.
#' The fitted dependence parameter is also stored in \code{dependence}.
#' \code{logLik} counts all fitted parameters, including Cox baseline increments,
#' in its degrees-of-freedom attribute. Cox increments and parametric densities
#' use different likelihood measures, so their log likelihoods or AIC values
#' should not be compared directly when choosing between those margin forms.
#' @return \code{pim_joint_fit} returns an object of class \code{pim_joint_fit},
#'   or \code{NULL} when CTW reduces to RCW without fitting a joint model.
#'   \code{coef}, \code{vcov}, and \code{logLik} return coefficients, their
#'   covariance matrix, and the fitted log likelihood, respectively.
#' @seealso \code{\link{pim_joint_model}}, \code{\link{gpim}}
#' @export
pim_joint_fit <- function(object) {
    if (inherits(object, "pim_joint_fit")) return(object)
    if (!inherits(object, "prior_pim"))
        stop("'object' must be a fitted gpim object.", call. = FALSE)
    if (is.null(object$joint))
        stop("This fit does not contain a CTW joint model.", call. = FALSE)
    if (inherits(object$joint, "error"))
        stop("The CTW joint model could not be fitted. ",
             conditionMessage(object$joint), call. = FALSE)
    if (isTRUE(object$joint$not_required)) return(NULL)
    object$joint
}

#' @rdname pim_joint_fit
#' @export
coef.pim_joint_fit <- function(object, ...) object$coefficients

#' @rdname pim_joint_fit
#' @export
vcov.pim_joint_fit <- function(object, ...) {
    influence <- t(solve(object$information, t(object$scores)))
    influence <- sweep(influence, 2L, colMeans(influence), "-")
    ans <- crossprod(influence) / object$n^2
    dimnames(ans) <- list(object$parameter_names, object$parameter_names)
    ans
}

#' @rdname pim_joint_fit
#' @export
logLik.pim_joint_fit <- function(object, ...) {
    structure(object$loglik, df = length(object$coefficients), nobs = object$n,
              class = "logLik")
}

#' @rdname pim_joint_fit
#' @export
print.pim_joint_fit <- function(x, ...) {
    cat("Joint model for conditional tie weighting\n")
    cat("Copula ", x$model$copula, "\n", sep = "")
    cat("Margins ", paste(x$model$margins, collapse = ", "), "\n", sep = "")
    if (!is.null(x$dependence))
        cat("Dependence parameter ", format(x$dependence, digits = 4L), "\n", sep = "")
    cat("Log likelihood ", format(x$loglik, digits = 7L), "\n", sep = "")
    keep <- !grepl("\\.log_jump[0-9]+$|^copula\\.", names(x$coefficients))
    print(x$coefficients[keep], ...)
    invisible(x)
}
