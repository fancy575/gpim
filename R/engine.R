# Numerical implementation of the CHW, RCW, and CTW estimating equations.
# All routines use observed event and censoring histories.

.pim_rowsum <- function(x, group, n) {
  x <- as.matrix(x)
  ans <- matrix(0, nrow = n, ncol = ncol(x))
  if (!length(group) || !ncol(x)) return(ans)
  z <- rowsum(x, group = group, reorder = FALSE)
  ans[as.integer(rownames(z)), ] <- z
  ans
}

.pim_cumsum <- function(x, reverse = FALSE) {
  x <- as.matrix(x)
  if (!nrow(x) || !ncol(x)) return(x)
  if (reverse) x <- x[nrow(x):1L, , drop = FALSE]
  ans <- vapply(seq_len(ncol(x)), function(j) cumsum(x[, j]),
                numeric(nrow(x)))
  ans <- matrix(ans, nrow = nrow(x), ncol = ncol(x))
  if (reverse) ans <- ans[nrow(ans):1L, , drop = FALSE]
  ans
}

.pim_at <- function(x, index) {
  if (is.matrix(x)) {
    rbind(rep(0, ncol(x)), x)[index + 1L, , drop = FALSE]
  } else {
    c(0, x)[index + 1L]
  }
}

.pim_left_index <- function(time, grid) {
  findInterval(time, grid, left.open = TRUE)
}

.pim_cox <- function(u, censor, W, control) {
  n <- length(u)
  W <- as.matrix(W)
  q <- ncol(W)
  center <- if (q) colMeans(W) else numeric()
  Wc <- if (q) sweep(W, 2L, center, "-") else W
  gamma <- stats::setNames(rep(0, q), colnames(W))
  event_time <- sort(unique(u[censor]))
  ne <- length(event_time)
  fit <- NULL
  fit_warnings <- character()
  score_norm <- 0
  if (q && ne) {
    fit <- withCallingHandlers(
      survival::coxph(survival::Surv(u, censor) ~ Wc,
                      ties = "breslow", singular.ok = FALSE,
                      control = survival::coxph.control(
                        iter.max = control$cox_maxit,
                        eps = control$cox_tol, timefix = FALSE),
                      x = TRUE, y = TRUE),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    gamma[] <- stats::coef(fit)
    if (any(!is.finite(gamma)) ||
        any(grepl("infinite|converge", fit_warnings, ignore.case = TRUE))) {
      stop("The censoring Cox model did not have a finite converged solution. ",
           paste(unique(fit_warnings), collapse = " "), call. = FALSE)
    }
  }
  exp_eta <- if (q) exp(drop(Wc %*% gamma)) else rep(1, n)
  if (any(!is.finite(exp_eta))) {
    stop("The fitted censoring risk scores are not finite.", call. = FALSE)
  }
  risk0 <- d_lambda <- numeric(ne)
  xbar <- matrix(0, ne, q)
  information <- matrix(0, q, q)
  event_n <- tabulate(match(u[censor], event_time), nbins = ne)
  for (k in seq_len(ne)) {
    at_risk <- u >= event_time[k]
    er <- exp_eta[at_risk]
    risk_sum <- sum(er)
    risk0[k] <- risk_sum / n
    d_lambda[k] <- event_n[k] / risk_sum
    if (q) {
      ww <- Wc[at_risk, , drop = FALSE]
      xb <- colSums(ww * er) / risk_sum
      xbar[k, ] <- xb
      information <- information + event_n[k] *
        (crossprod(ww, ww * er) / risk_sum - tcrossprod(xb))
    }
  }
  h0 <- cumsum(d_lambda)
  a_cum <- .pim_cumsum(xbar * d_lambda)
  if_gamma <- matrix(0, n, q, dimnames = list(NULL, colnames(W)))
  if (q && ne) {
    inv_info <- tryCatch(solve(information / n), error = function(e) NULL)
    if (is.null(inv_info) || any(!is.finite(inv_info))) {
      stop("The censoring Cox information matrix is singular.", call. = FALSE)
    }
    upper <- findInterval(u, event_time)
    psi <- -exp_eta * (Wc * .pim_at(h0, upper) - .pim_at(a_cum, upper))
    ei <- match(u[censor], event_time)
    psi[censor, ] <- psi[censor, , drop = FALSE] +
      Wc[censor, , drop = FALSE] - xbar[ei, , drop = FALSE]
    if_gamma <- psi %*% inv_info
    colnames(if_gamma) <- colnames(W)
    cox_score <- colSums(psi)
    score_norm <- sqrt(sum(cox_score^2)) / max(1, sum(censor))
    if (!is.finite(score_norm) || score_norm > control$cox_score_tol) {
      stop("The censoring Cox score did not converge sufficiently.", call. = FALSE)
    }
  }
  list(fit = fit, coefficients = gamma, center = center, W = Wc,
       followup = u, status = censor, exp_eta = exp_eta,
       event_time = event_time, d_lambda = d_lambda, h0 = h0,
       risk0 = risk0, xbar = xbar, a_cum = a_cum,
       if_gamma = if_gamma, information = information,
       n_events = sum(censor), no_censoring = !ne,
       iterations = if (is.null(fit)) 0L else fit$iter,
       score_norm = score_norm, score_tolerance = control$cox_score_tol,
       likelihood_tolerance = control$cox_tol,
       warnings = unique(fit_warnings), ties = "breslow",
       survival_convention = "exp(-Lambda_C(t-) * exp(W gamma))")
}

.pim_probability <- function(log_inverse, floor) {
  probability <- exp(-log_inverse)
  if (floor > 0) {
    weight <- exp(pmin(log_inverse, -log(floor)))
    derivative_active <- log_inverse < -log(floor)
  } else {
    weight <- exp(log_inverse)
    derivative_active <- rep(TRUE, length(weight))
  }
  if (any(!is.finite(weight))) {
    stop("An inverse observation probability is not finite. ",
         "Use an earlier restriction time or an explicitly chosen weight_floor.",
         call. = FALSE)
  }
  list(probability = probability, weight = weight,
       derivative_active = derivative_active)
}

.pim_prepare_horizon <- function(dat, L, control) {
  n <- nrow(dat$X)
  if (n < 3L) stop("At least three independent subjects are required.", call. = FALSE)
  pairs <- gpim_pairs_cpp(dat$u, dat$death, dat$event_times, dat$priorities, L)
  ii <- pairs$i
  jj <- pairs$j
  h <- pairs$h
  u <- pairs$u
  died <- pairs$died
  censor <- pairs$censor
  complete <- pairs$complete
  s <- pairs$s
  resolved <- pairs$resolved
  ch <- pairs$ch
  priority <- pairs$priority
  cox <- .pim_cox(u, censor, dat$W, control)
  m <- length(ii)
  rho_ch <- rho_dl <- numeric(m)
  pi_ch <- pi_dl <- rep(NA_real_, m)
  deriv_ch <- deriv_dl <- rep(FALSE, m)
  left_u <- .pim_left_index(u, cox$event_time)
  left_s <- integer(m)
  left_s[resolved] <- .pim_left_index(s[resolved], cox$event_time)
  if (any(ch)) {
    li <- .pim_at(cox$h0, left_u[ii[ch]]) * cox$exp_eta[ii[ch]] +
      .pim_at(cox$h0, left_u[jj[ch]]) * cox$exp_eta[jj[ch]]
    z <- .pim_probability(li, control$weight_floor)
    rho_ch[ch] <- z$weight
    pi_ch[ch] <- z$probability
    deriv_ch[ch] <- z$derivative_active
  }
  if (any(resolved)) {
    li <- .pim_at(cox$h0, left_s[resolved]) *
      (cox$exp_eta[ii[resolved]] + cox$exp_eta[jj[resolved]])
    z <- .pim_probability(li, control$weight_floor)
    rho_dl[resolved] <- z$weight
    pi_dl[resolved] <- z$probability
    deriv_dl[resolved] <- z$derivative_active
  }
  Z <- dat$X[ii, , drop = FALSE] - dat$X[jj, , drop = FALSE]
  list(dat = dat, n = n, m = m, L = L, i = ii, j = jj, Z = Z,
       u = u, died = died, censor = censor, complete = complete,
       s = s, resolved = resolved, ch = ch, h = h, priority = priority,
       cox = cox, rho_ch = rho_ch, rho_dl = rho_dl,
       pi_ch = pi_ch, pi_dl = pi_dl,
       deriv_ch = deriv_ch, deriv_dl = deriv_dl,
       left_u = left_u, left_s = left_s, native_pairs = pairs)
}

.pim_add_ctw <- function(cache, joint, control) {
  if (inherits(joint, "error")) {
    cache$ctw_error <- conditionMessage(joint)
    return(cache)
  }
  if (isTRUE(joint$not_required)) {
    cache$ctw <- list(y = cache$h, weight = cache$rho_dl,
      diagnostic_weight = cache$rho_dl, retained = cache$resolved,
      active = cache$resolved, derivative = cache$deriv_dl,
      left_i = cache$left_s, left_j = cache$left_s,
      probability = cache$pi_dl, q = rep(0,cache$n),
      factor = rep(0,cache$m), jacobian = matrix(0,cache$n,0),
      lower = rep(FALSE,cache$m))
    cache$joint <- joint
    return(cache)
  }
  pred <- .pim_joint_predict(joint, cache$dat, cache$L)
  q <- as.numeric(pred$probability[,1L])
  pair <- gpim_ctw_pairs_cpp(cache$native_pairs,
                            cache$dat$event_times[,1L], cache$L)
  li <- integer(cache$m)
  used <- pair$copied | pair$lower
  li[used] <- .pim_left_index(pair$time[used],cache$cox$event_time)
  base <- weight <- numeric(cache$m)
  prob <- rep(NA_real_,cache$m)
  deriv <- rep(FALSE,cache$m)
  copied <- pair$copied
  weight[copied] <- cache$rho_dl[copied]
  prob[copied] <- cache$pi_dl[copied]
  deriv[copied] <- cache$deriv_dl[copied]
  lower <- pair$lower
  if (any(lower)) {
    h <- .pim_at(cache$cox$h0,li[lower]) *
      (cache$cox$exp_eta[cache$i[lower]] + cache$cox$exp_eta[cache$j[lower]])
    z <- .pim_probability(h,control$weight_floor)
    base[lower] <- z$weight
    weight[lower] <- z$weight * q[cache$i[lower]] * q[cache$j[lower]]
    prob[lower] <- z$probability
    deriv[lower] <- z$derivative_active
  }
  cache$ctw <- list(y = pair$y, weight = weight,
    diagnostic_weight = weight, retained = weight > 0, active = weight > 0,
    derivative = deriv, left_i = li, left_j = li, probability = prob,
    q = q, factor = base, jacobian = pred$jacobian[[1L]], lower = lower)
  cache$joint <- joint
  cache
}

.pim_use <- function(cache, method) {
  if (method == "CTW") {
    if (!is.null(cache$ctw_error)) stop(cache$ctw_error,call.=FALSE)
    if (is.null(cache$ctw)) stop("CTW joint-model information is missing.",call.=FALSE)
    return(cache$ctw)
  }
  if (method == "CHW") {
    return(list(y = cache$h, weight = cache$rho_ch,
      diagnostic_weight = cache$rho_ch, retained = cache$ch,
      active = cache$ch, derivative = cache$deriv_ch,
      left_i = cache$left_u[cache$i], left_j = cache$left_u[cache$j],
      probability = cache$pi_ch))
  }
  if (method != "RCW") stop("Unknown weighting method.",call.=FALSE)
  list(y = cache$h, weight = cache$rho_dl,
       diagnostic_weight = cache$rho_dl, retained = cache$resolved,
       active = cache$resolved, derivative = cache$deriv_dl,
       left_i = cache$left_s, left_j = cache$left_s,
       probability = cache$pi_dl)
}

.pim_equation <- function(beta, cache, use, link, details = FALSE) {
  state <- list(n = cache$n)
  if (details) state <- list(n = cache$n, i = cache$i, j = cache$j,
    left_i = use$left_i, left_j = use$left_j, derivative = use$derivative,
    risk = cache$cox$exp_eta, h0 = cache$cox$h0,
    W = cache$cox$W, a_cum = cache$cox$a_cum,
    q = if (is.null(use$q)) rep(0,cache$n) else use$q,
    factor = if (is.null(use$factor)) rep(0,cache$m) else use$factor)
  gpim_equation_cpp(cache$Z,as.numeric(beta),use$y,use$weight,
    match(link,c("identity","logit","probit"))-1L,details,state)
}

.pim_solve <- function(cache, use, link, control) {
  p <- ncol(cache$Z)
  empty_history <- data.frame(iteration = integer(), score_norm = numeric(),
                              step_norm = numeric(), halvings = integer())
  bad <- list(beta = stats::setNames(rep(NA_real_, p), colnames(cache$Z)),
              converged = FALSE, boundary = FALSE, iterations = 0L,
              step_norm = NA_real_, iteration_history = empty_history,
              start_used = NULL,
              message = "The weighted covariate difference design is singular.")
  if (!any(use$active)) {
    bad$message <- "No comparisons satisfy this method's observation requirement."
    return(bad)
  }
  M <- crossprod(cache$Z, cache$Z * use$weight) / cache$m
  rhs <- drop(crossprod(cache$Z, use$weight * (2 * use$y - 1))) / cache$m
  linear <- tryCatch(drop(solve(M, rhs)), error = function(e) NULL)
  if (is.null(linear) || any(!is.finite(linear))) return(bad)
  if (link == "identity") {
    eta <- drop(cache$Z %*% linear)
    bound <- max(abs(eta))
    beta <- linear
    boundary <- bound >= 1 - 1e-10
    projection_message <- NULL
    if (boundary && p == 1L) {
      beta <- sign(linear) * (1 - 1e-12) / max(abs(cache$Z))
    } else if (boundary) {
      X <- sweep(cache$dat$X, 2L, colMeans(cache$dat$X), "-")
      ui <- rbind(cbind(X, 1), cbind(-X, -1))
      ci <- c(rep(0, cache$n), rep(-1, cache$n))
      objective <- function(theta) {
        b <- theta[seq_len(p)]
        drop(crossprod(b, M %*% b) / 2 - crossprod(rhs, b))
      }
      gradient <- function(theta) {
        c(drop(M %*% theta[seq_len(p)]) - rhs, 0)
      }
      constrained <- tryCatch(stats::constrOptim(
        c(rep(0, p), 0.5), objective, gradient, ui = ui, ci = ci,
        mu = 1e-8, method = "BFGS", outer.iterations = 100,
        # This only supplies an admissible descriptive boundary projection.
        # It does not turn an inadmissible root into a regular converged fit.
        control = list(maxit = 500L, reltol = 1e-11)),
        error = function(e) NULL)
      if (is.null(constrained)) {
        bad$boundary <- TRUE
        bad$message <- "The identity root is inadmissible and constrained fitting failed."
        return(bad)
      }
      beta <- constrained$par[seq_len(p)]
      if (constrained$convergence != 0L)
        projection_message <- paste("The numerical boundary projection did not converge;",
                                    "its coefficient values are provisional.")
    }
    score_norm <- sqrt(sum(.pim_equation(beta, cache, use, link)$mean_score^2))
    # The identity solution is algebraic. The same score check is still applied
    # to an interior root, so an unattainably strict tolerance is not ignored.
    converged <- !boundary && is.finite(score_norm) && score_norm <= control$tol
    result_message <- if (boundary)
      "The identity root reaches or exceeds the admissible boundary. Regular Wald inference is unavailable."
    else if (converged) "Closed-form identity estimating-equation solution."
    else "The closed-form identity solution did not meet the requested score tolerance."
    if (!is.null(projection_message))
      result_message <- paste(result_message, projection_message)
    if (!is.null(control$start))
      result_message <- paste(result_message, "start is ignored for the identity link.")
    if (control$trace)
      message("  Identity solution: score norm = ", format(score_norm, digits = 5),
              if (boundary) "; boundary solution." else ".")
    return(list(beta = stats::setNames(beta, colnames(cache$Z)),
                unconstrained = stats::setNames(linear, colnames(cache$Z)),
                converged = converged, boundary = boundary, iterations = 1L,
                step_norm = NA_real_, start_used = NULL,
                iteration_history = data.frame(iteration = 1L, score_norm = score_norm,
                                               step_norm = NA_real_, halvings = NA_integer_),
                message = result_message))
  }
  beta <- if (!is.null(control$start)) control$start else
    linear / (2 * .pim_link(0, link)$prime)
  start_used <- stats::setNames(beta, colnames(cache$Z))
  current <- .pim_equation(beta, cache, use, link)
  sn <- sqrt(sum(current$mean_score^2))
  history <- data.frame(iteration = 0L, score_norm = sn,
                        step_norm = NA_real_, halvings = NA_integer_)
  if (control$trace)
    message("  Iteration 0: score norm = ", format(sn, digits = 5), ".")
  converged <- FALSE
  message <- "Maximum number of estimating-equation iterations reached."
  iterations <- 0L
  accepted_step_norm <- NA_real_
  for (iteration in seq_len(control$maxit)) {
    if (!is.finite(sn)) {
      message <- "The estimating equation is non-finite."
      break
    }
    # Success is defined by || U(beta) ||_2 <= tol, where U is the sum of
    # pair scores divided by the total number of unordered subject pairs.
    if (sn <= control$tol) {
      converged <- TRUE
      message <- "Estimating-equation score tolerance reached."
      break
    }
    step <- tryCatch(drop(solve(current$mean_jacobian, current$mean_score)),
                     error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      message <- "The estimating-equation Jacobian is singular."
      break
    }
    step_norm <- sqrt(sum(step^2))
    if (step_norm > control$max_step) step <- control$max_step * step / step_norm
    accepted <- FALSE
    for (ls in 0:control$max_halving) {
      candidate_beta <- beta - step * 0.5^ls
      candidate <- .pim_equation(candidate_beta, cache, use, link)
      candidate_norm <- sqrt(sum(candidate$mean_score^2))
      if (is.finite(candidate_norm) && candidate_norm < sn) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      message <- "Estimating-equation line search failed."
      break
    }
    accepted_step_norm <- sqrt(sum((candidate_beta - beta)^2))
    beta <- candidate_beta
    current <- candidate
    sn <- candidate_norm
    iterations <- iteration
    history <- rbind(history,
                     data.frame(iteration = iteration, score_norm = sn,
                                step_norm = accepted_step_norm, halvings = ls))
    if (control$trace)
      base::message("  Iteration ", iteration, ": score norm = ",
                    format(sn, digits = 5), ", step norm = ",
                    format(accepted_step_norm, digits = 5), ", halvings = ", ls, ".")
    if (sn <= control$tol) {
      converged <- TRUE
      message <- "Estimating-equation score tolerance reached."
      break
    }
    # A small accepted step means stagnation, not convergence. This relative
    # threshold does not replace the score condition, even near a flat root.
    if (accepted_step_norm <= control$step_tol * (1 + sqrt(sum(beta^2)))) {
      message <- "The parameter step is small but the score tolerance was not reached."
      break
    }
  }
  list(beta = stats::setNames(beta, colnames(cache$Z)), converged = converged,
       boundary = FALSE, iterations = iterations, step_norm = accepted_step_norm,
       iteration_history = history, start_used = start_used, message = message)
}

.pim_nuisance <- function(cache, eq, use, method) {
  n <- cache$n
  p <- ncol(cache$Z)
  cox <- cache$cox
  ne <- length(cox$event_time)
  nc <- ncol(cox$W)
  baseline <- coefficient <- joint <- matrix(0,n,p)
  if (ne) {
    qstart <- .pim_cumsum(eq$censor_buckets,reverse=TRUE)
    cumulative <- .pim_cumsum(qstart * (cox$d_lambda/cox$risk0))
    upper <- findInterval(cox$followup,cox$event_time)
    baseline <- -.pim_at(cumulative,upper) * cox$exp_eta
    ev <- cox$status
    if (any(ev)) {
      at <- match(cox$followup[ev],cox$event_time)
      baseline[ev,] <- baseline[ev,,drop=FALSE] +
        qstart[at,,drop=FALSE]/cox$risk0[at]
    }
    if (nc) coefficient <- cox$if_gamma %*% t(eq$censor_gamma)
  }
  derivative <- matrix(0,p,0)
  if (method == "CTW" && ncol(use$jacobian)) {
    derivative <- crossprod(eq$joint_loading,use$jacobian)
    # The information includes copula, regression and all baseline parameters.
    # Scores and censoring martingales remain in the same subject row.
    if (any(derivative != 0))
      joint <- cache$joint$scores %*%
        solve(cache$joint$information,t(derivative))
  }
  list(total=baseline+coefficient+joint, baseline=baseline,
       coefficient=coefficient,joint=joint,joint_derivative=derivative)
}

.pim_pair_table <- function(cache) {
  ans <- data.frame(i = cache$i, j = cache$j,
             id_i = cache$dat$id[cache$i], id_j = cache$dat$id[cache$j],
             response = ifelse(cache$resolved, cache$h, NA_real_),
             resolved = cache$resolved, complete_history = cache$ch,
             determining_priority = cache$priority,
             ascertainment_time = ifelse(cache$resolved, cache$s, NA_real_),
             probability_CHW = cache$pi_ch, probability_RCW = cache$pi_dl,
             weight_CHW = cache$rho_ch, weight_RCW = cache$rho_dl,
             stringsAsFactors = FALSE)
  if (!is.null(cache$ctw)) {
    ans$response_CTW <- ifelse(cache$ctw$retained, cache$ctw$y, NA_real_)
    ans$probability_CTW <- cache$ctw$probability
    ans$weight_CTW <- cache$ctw$weight
    ans$conditional_nonfatal_comparison <- cache$ctw$lower
  }
  ans
}

.pim_fit_cached <- function(cache, method, link, control) {
  use <- .pim_use(cache, method)
  solution <- .pim_solve(cache, use, link, control)
  p <- ncol(cache$Z)
  cn <- colnames(cache$Z)
  empty <- matrix(NA_real_, p, p, dimnames = list(cn, cn))
  out <- list(coefficients = solution$beta, vcov = empty,
              se = stats::setNames(rep(NA_real_, p), cn),
              converged = solution$converged, boundary = solution$boundary,
              status = if (solution$boundary) "boundary" else
                if (solution$converged) "converged" else
                if (any(!is.finite(solution$beta))) "failed" else "nonconverged",
              message = solution$message, iterations = solution$iterations,
              score_norm = NA_real_, score_tolerance = control$tol,
              step_norm = solution$step_norm, step_tolerance = control$step_tol,
              iteration_history = solution$iteration_history,
              start_used = solution$start_used, bread = empty, meat = empty,
              influence = matrix(NA_real_, cache$n, p,
                                 dimnames = list(cache$dat$id, cn)),
              censoring = cache$cox, variance_ok = FALSE,
              L = cache$L, method = method, link = link)
  pi <- use$probability
  dw <- use$diagnostic_weight[use$retained]
  observed_pi <- pi[use$retained]
  out$diagnostics <- list(
    n_subjects = cache$n, n_pairs = cache$m,
    n_resolved = sum(cache$resolved), n_complete_history = sum(cache$ch),
    n_retained = sum(use$retained), n_active = sum(use$active),
    n_genuine_ties = sum(cache$resolved & cache$h == 0.5),
    n_censoring_events = cache$cox$n_events,
    n_completed_histories = sum(cache$complete),
    minimum_observation_probability = if (length(observed_pi)) min(observed_pi) else NA_real_,
    maximum_weight = if (length(dw)) max(dw) else NA_real_,
    effective_pairs = if (length(dw) && sum(dw^2) > 0) sum(dw)^2 / sum(dw^2) else 0,
    weight_floor = control$weight_floor,
    n_truncated = sum(use$retained & !use$derivative),
    score_tolerance = control$tol, step_tolerance = control$step_tol,
    start_ignored = link == "identity" && !is.null(control$start),
    censoring_adjustment = "Cox/Breslow martingale contribution",
    variance_method = "Centered subject influence covariance",
    identity_scale = if (link == "identity") "net benefit" else NULL)
  if (isTRUE(control$keep_pairs)) out$pairs <- .pim_pair_table(cache)
  if (any(!is.finite(solution$beta))) return(out)
  eq <- .pim_equation(solution$beta, cache, use, link, details = TRUE)
  out$score_norm <- sqrt(sum(eq$mean_score^2))
  out$bread <- -eq$mean_jacobian
  dimnames(out$bread) <- list(cn, cn)
  out$diagnostics$fitted_probability_range <- eq$fitted_range
  out$diagnostics$max_abs_linear_predictor <- eq$max_abs_eta
  if (link == "identity") out$unconstrained_coefficients <- solution$unconstrained
  if (!out$converged || out$boundary) return(out)
  inv_a <- tryCatch(solve(out$bread), error = function(e) NULL)
  if (is.null(inv_a) || any(!is.finite(inv_a))) {
    out$message <- paste(out$message, "The sandwich bread is singular.")
    return(out)
  }
  nuisance <- .pim_nuisance(cache,eq,use,method)
  pair <- sweep(2 * eq$subject_sum/(cache$n-1L),2L,
                2 * eq$mean_score,"-")
  score_components <- list(pair=pair,
    censoring_baseline=nuisance$baseline,
    censoring_regression=nuisance$coefficient,joint=nuisance$joint)
  out$influence_components <- lapply(score_components,function(x) {
    v <- x %*% t(inv_a)
    v <- sweep(v,2L,colMeans(v),"-")
    dimnames(v) <- list(cache$dat$id,cn)
    v
  })
  out$influence <- Reduce(`+`,out$influence_components)
  score_influence <- pair + nuisance$total
  score_influence <- sweep(score_influence,2L,colMeans(score_influence),"-")
  out$meat <- crossprod(score_influence)/cache$n^2
  out$vcov <- crossprod(out$influence)/cache$n^2
  dimnames(out$meat) <- dimnames(out$vcov) <- list(cn,cn)
  out$nuisance <- nuisance
  out$diagnostics$bread_condition <- kappa(out$bread,exact=TRUE)
  out$diagnostics$nuisance_rms <- sqrt(mean(rowSums(nuisance$total^2)))
  out$diagnostics$baseline_martingale_rms <- sqrt(mean(rowSums(nuisance$baseline^2)))
  out$diagnostics$coefficient_martingale_rms <- sqrt(mean(rowSums(nuisance$coefficient^2)))
  out$diagnostics$joint_model_rms <- sqrt(mean(rowSums(nuisance$joint^2)))
  out$diagnostics$adjusted_score_mean <- colMeans(score_influence)
  out$variance_ok <- all(is.finite(out$vcov)) && all(diag(out$vcov)>=0)
  if (out$variance_ok) {
    out$se <- sqrt(diag(out$vcov))
    eig <- eigen(out$vcov,symmetric=TRUE,only.values=TRUE)$values
    out$diagnostics$covariance_eigenvalues <- eig
    out$diagnostics$covariance_positive_semidefinite <-
      min(eig)>=-1e-10*max(1,max(abs(eig)))
  } else out$message <- paste(out$message,"The subject influence covariance is not finite.")

  out
}

