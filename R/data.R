# Long-format endpoint data and subject-level formula preparation.

#' @name HierSurv
#' @aliases HierSurv
#' @title Construct a Long-Format Prioritized Survival Response
#' @description Constructs the response for prioritized survival outcomes in a \code{gpim} formula. The response contains the observed time and event code for each row of the long-format data.
#' @param time Numeric vector of nonnegative observed event times or terminal follow-up times.
#' @param status Integer vector of event codes, with one entry for each value of \code{time}. Code 0 denotes terminal censoring, code 1 denotes the terminal fatal event, and codes 2 and above identify nonfatal endpoint types.
#' @details Use \code{HierSurv(time, status)} on the left-hand side of the model formula. The fitting function assembles and validates each individual's event history using the subject identifier supplied through \code{id}.
#'
#' Every individual must have exactly one terminal record with status 0 or 1. Nonfatal events have separate rows and cannot occur after the terminal time. The status code identifies the endpoint type, rather than a binary event indicator. If a nonfatal event is not recorded during follow-up, omit its event row. The order of endpoint comparisons is specified by the \code{priorities} argument to \code{gpim}.
#' @return A two-column response object containing the supplied time and status.
#' @seealso \code{\link{gpim}}
#' @export
HierSurv <- function(time, status) {
    if (!is.numeric(time) || !is.numeric(status) ||
        is.matrix(time) || is.matrix(status)) {
        stop("'time' and 'status' must be numeric vectors.", call. = FALSE)
    }
    if (length(time) != length(status)) {
        stop("'time' and 'status' must have the same length.", call. = FALSE)
    }
    if (any(!is.na(time) & (!is.finite(time) | time < 0))) {
        stop("Observed times must be finite and nonnegative.", call. = FALSE)
    }
    if (any(!is.na(status) &
            (!is.finite(status) | status < 0 | status != floor(status) |
             status > .Machine$integer.max))) {
        stop("Event codes must be nonnegative integers.", call. = FALSE)
    }
    ans <- cbind(time = as.numeric(time), status = as.numeric(status))
    class(ans) <- c("HierSurv", "matrix", "array")
    ans
}

.pim_na_action <- function(na.action) {
    if (is.function(na.action)) {
        if (identical(na.action, stats::na.fail)) return("na.fail")
        if (identical(na.action, stats::na.omit)) return("na.omit")
    }
    if (is.character(na.action) && length(na.action) == 1L &&
        !is.na(na.action) && na.action %in% c("na.fail", "na.omit")) {
        return(na.action)
    }
    stop("'na.action' must be 'na.fail' or 'na.omit' (or the corresponding function).",
         call. = FALSE)
}

.pim_rhs_formula <- function(formula) {
    if (!inherits(formula, "formula") || !length(formula) %in% c(2L, 3L)) {
        stop("Supply an R formula, such as HierSurv(time, status) ~ treatment + age.",
             call. = FALSE)
    }
    if (length(formula) == 2L) return(formula)
    ans <- stats::as.formula(call("~", formula[[3L]]), env = environment(formula))
    ans
}

.pim_response <- function(formula, data, time, status) {
    if (length(formula) == 2L) {
        for (nm in c(time, status)) {
            if (!is.character(nm) || length(nm) != 1L || is.na(nm) ||
                !nm %in% names(data)) {
                stop("For a one-sided formula, 'time' and 'status' must name columns in 'data'.",
                     call. = FALSE)
            }
        }
        return(list(response = HierSurv(data[[time]], data[[status]]),
                    variables = c(time, status)))
    }
    lhs <- formula[[2L]]
    is_hier <- is.call(lhs) && (
        identical(lhs[[1L]], as.name("HierSurv")) ||
        (is.call(lhs[[1L]]) && length(lhs[[1L]]) == 3L &&
         identical(lhs[[1L]][[1L]], as.name("::")) &&
         identical(lhs[[1L]][[3L]], as.name("HierSurv"))))
    if (!is_hier) {
        stop("The formula response must be HierSurv(time, status). Alternatively use a one-sided formula.",
             call. = FALSE)
    }
    matched <- tryCatch(match.call(definition = HierSurv, call = lhs),
                        error = function(e) {
                            stop("Use HierSurv(time, status) on the left-hand side: ",
                                 conditionMessage(e), call. = FALSE)
                        })
    if (is.null(matched$time) || is.null(matched$status)) {
        stop("Both 'time' and 'status' are required in HierSurv().", call. = FALSE)
    }
    tt <- eval(matched$time, envir = data, enclos = environment(formula))
    ss <- eval(matched$status, envir = data, enclos = environment(formula))
    response <- HierSurv(tt, ss)
    if (nrow(response) != nrow(data)) {
        stop("The HierSurv response must have one entry for every long-format data row.",
             call. = FALSE)
    }
    list(response = response, variables = all.vars(lhs))
}

.pim_terms <- function(rhs, data, exclude, label) {
    # Expand a dot over candidate baseline columns, excluding the response and ID.
    candidate <- data[, setdiff(names(data), exclude), drop = FALSE]
    ans <- tryCatch(stats::terms(rhs, data = candidate), error = function(e) {
        stop(label, " formula could not be evaluated: ", conditionMessage(e), call. = FALSE)
    })
    if (length(attr(ans, "offset"))) {
        stop("Offsets are not supported in the ", tolower(label),
             " formula. Use explicitly modeled baseline predictors.", call. = FALSE)
    }
    ans
}

.pim_missing_rows <- function(x) {
    if (is.data.frame(x) || is.matrix(x)) return(!stats::complete.cases(x))
    if (is.list(x) && !is.factor(x)) {
        stop("Formula covariates must be atomic vectors or matrix columns, not list columns.",
             call. = FALSE)
    }
    is.na(x)
}

.pim_constant_within <- function(x, groups, variable) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) {
        for (ii in groups) {
            first <- x[ii[1L], , drop = FALSE]
            if (length(ii) > 1L && !all(vapply(ii[-1L], function(j) {
                identical(unname(x[j, , drop = FALSE]), unname(first))
            }, logical(1L)))) {
                stop("Baseline covariate '", variable,
                     "' varies within an individual. Time-varying covariates are not supported.",
                     call. = FALSE)
            }
        }
    } else {
        if (is.list(x)) {
            stop("Baseline covariate '", variable, "' must not be a list column.",
                 call. = FALSE)
        }
        bad <- vapply(groups, function(ii) length(unique(x[ii])) > 1L, logical(1L))
        if (any(bad)) {
            stop("Baseline covariate '", variable,
                 "' varies within an individual. Time-varying covariates are not supported.",
                 call. = FALSE)
        }
    }
    invisible(NULL)
}

.pim_make_frame <- function(tt, subjects, label) {
    tryCatch(stats::model.frame(tt, data = subjects, na.action = stats::na.pass,
                               drop.unused.levels = TRUE),
             error = function(e) {
                 stop(label, " formula could not be evaluated on one row per subject: ",
                      conditionMessage(e), call. = FALSE)
             })
}

.pim_model_matrix <- function(mf, label, allow_empty = FALSE) {
    mm <- tryCatch(stats::model.matrix(attr(mf, "terms"), data = mf),
                   error = function(e) {
                       stop(label, " model matrix could not be constructed: ",
                            conditionMessage(e), call. = FALSE)
                   })
    contrasts <- attr(mm, "contrasts")
    keep <- colnames(mm) != "(Intercept)"
    mm <- mm[, keep, drop = FALSE]
    if (!ncol(mm) && !allow_empty) {
        stop("The outcome formula must contain at least one nonconstant covariate. A free PIM intercept is not identifiable.",
             call. = FALSE)
    }
    if (any(!is.finite(mm))) {
        stop(label, " model matrix contains non-finite values. Check transformations and covariates.",
             call. = FALSE)
    }
    if (ncol(mm)) {
        centered <- sweep(mm, 2L, colMeans(mm), FUN = "-")
        z <- qr(centered, tol = 1e-9)
        if (z$rank < ncol(mm)) {
            offending <- colnames(mm)[z$pivot[seq.int(z$rank + 1L, ncol(mm))]]
            stop(label, " covariate differences are rank deficient. Nonestimable column(s): ",
                 paste(offending, collapse = ", "),
                 ". Remove redundant or constant predictors. Use the usual formula intercept for factor contrasts; the PIM intercept is removed automatically.",
                 call. = FALSE)
        }
    }
    xlevels <- lapply(mf[vapply(mf, function(z) is.factor(z) || is.character(z),
                               logical(1L))], function(z) {
        if (is.factor(z)) levels(z) else levels(factor(z))
    })
    list(matrix = mm, contrasts = contrasts, xlevels = xlevels,
         terms = attr(mf, "terms"))
}

.pim_prepare_data <- function(formula, data, id = "id", time = "time",
                              status = "status", priorities = NULL,
                              censor_formula = NULL, recurrent = "first",
                              na.action = "na.fail") {
    if (!is.data.frame(data) || !nrow(data)) {
        stop("'data' must be a nonempty long-format data frame.", call. = FALSE)
    }
    if (anyDuplicated(names(data))) {
        stop("Column names in 'data' must be unique.", call. = FALSE)
    }
    if (!is.character(id) || length(id) != 1L || is.na(id) || !id %in% names(data)) {
        stop("'id' must be the name of the individual identifier column.", call. = FALSE)
    }
    ids <- data[[id]]
    if (!(is.numeric(ids) || is.character(ids) || is.factor(ids)) ||
        is.matrix(ids) || anyNA(ids) || any(as.character(ids) == "")) {
        stop("Every row must have a nonmissing numeric, character, or factor individual ID.",
             call. = FALSE)
    }
    if (is.numeric(ids) && any(!is.finite(ids))) {
        stop("Numeric individual IDs must be finite.", call. = FALSE)
    }
    if (is.factor(ids)) ids <- as.character(ids)
    na_method <- .pim_na_action(na.action)
    recurrent <- match.arg(recurrent, c("first", "error"))
    rhs <- .pim_rhs_formula(formula)
    rr <- .pim_response(formula, data, time, status)
    response <- rr$response
    if (is.null(censor_formula)) {
        censor_formula <- rhs
    } else if (!inherits(censor_formula, "formula") || length(censor_formula) != 2L) {
        stop("'censor_formula' must be one-sided, such as ~ age or ~ 1.", call. = FALSE)
    }
    excluded <- unique(c(id, rr$variables))
    outcome_terms <- .pim_terms(rhs, data, excluded, "Outcome")
    censor_terms <- .pim_terms(censor_formula, data, excluded, "Censoring")
    used <- unique(c(all.vars(outcome_terms), all.vars(censor_terms)))
    # Variables in the formula must be data columns. This avoids row-length
    # ambiguities when a formula environment contains a long-format vector.
    absent <- setdiff(used, names(data))
    if (length(absent)) {
        stop("Formula covariate(s) must be columns in 'data': ",
             paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    missing_row <- !stats::complete.cases(response)
    for (nm in used) missing_row <- missing_row | .pim_missing_rows(data[[nm]])
    omitted <- unique(ids[missing_row])
    if (length(omitted) && na_method == "na.fail") {
        stop("Missing time, event code, or modeled covariate. Use na.action = 'na.omit' to remove entire individuals.",
             call. = FALSE)
    }
    original_n <- length(unique(ids))
    keep <- !ids %in% omitted
    data <- data[keep, , drop = FALSE]
    response <- response[keep, , drop = FALSE]
    ids <- ids[keep]
    repeat {
        subject_ids <- unique(ids)
        subject_ids <- subject_ids[order(subject_ids)]
        if (length(subject_ids) < 3L) {
            stop("At least three individuals with complete required information are needed.",
                 call. = FALSE)
        }
        subject_index <- match(ids, subject_ids)
        groups <- split(seq_along(ids), factor(subject_index,
                                              levels = seq_along(subject_ids)))
        for (nm in used) .pim_constant_within(data[[nm]], groups, nm)
        first <- match(subject_ids, ids)
        subjects <- data[first, , drop = FALSE]
        rownames(subjects) <- as.character(subject_ids)
        outcome_frame <- .pim_make_frame(outcome_terms, subjects, "Outcome")
        censor_frame <- .pim_make_frame(censor_terms, subjects, "Censoring")
        complete <- stats::complete.cases(outcome_frame) &
                    stats::complete.cases(censor_frame)
        if (all(complete)) break
        if (na_method == "na.fail") {
            stop("A formula transformation produced missing values. Use na.action = 'na.omit' to remove entire individuals.",
                 call. = FALSE)
        }
        removed <- subject_ids[!complete]
        omitted <- unique(c(omitted, removed))
        keep <- !ids %in% removed
        data <- data[keep, , drop = FALSE]
        response <- response[keep, , drop = FALSE]
        ids <- ids[keep]
    }
    times <- response[, "time"]
    codes <- response[, "status"]
    if (is.null(priorities)) {
        priorities <- sort(unique(c(1, codes[codes > 0])))
    } else {
        if (!is.numeric(priorities) || !length(priorities) ||
            anyNA(priorities) || any(!is.finite(priorities)) ||
            any(priorities < 1 | priorities != floor(priorities) |
                priorities > .Machine$integer.max) ||
            priorities[1L] != 1 || anyDuplicated(priorities)) {
            stop("'priorities' must be distinct positive integer codes, with fatal code 1 first.",
                 call. = FALSE)
        }
        unlisted <- setdiff(unique(codes[codes > 0]), priorities)
        if (length(unlisted)) {
            stop("Observed event code(s) missing from 'priorities': ",
                 paste(unlisted, collapse = ", "),
                 ". Include every observed positive event code.", call. = FALSE)
        }
    }
    priorities <- as.integer(priorities)
    terminal <- codes %in% c(0, 1)
    nt <- vapply(groups, function(ii) sum(terminal[ii]), integer(1L))
    if (any(nt != 1L)) {
        offending <- as.character(subject_ids[nt != 1L])
        stop("Each individual must have exactly one terminal record with status 0 (censoring) or 1 (fatal). Check ID(s): ",
             paste(utils::head(offending, 8L), collapse = ", "), ".", call. = FALSE)
    }
    terminal_rows <- vapply(groups, function(ii) ii[terminal[ii]], integer(1L))
    u <- times[terminal_rows]
    death <- codes[terminal_rows] == 1L
    names(u) <- names(death) <- as.character(subject_ids)
    if (any(times > u[subject_index])) {
        bad <- unique(as.character(ids[times > u[subject_index]]))
        stop("An event occurs after that individual's terminal censoring or death record. Check ID(s): ",
             paste(utils::head(bad, 8L), collapse = ", "), ".", call. = FALSE)
    }
    # Preserve a useful one-row record with terminal follow-up columns. Modeled
    # baseline values equal those in the frames already constructed above.
    subjects <- data[terminal_rows, , drop = FALSE]
    rownames(subjects) <- as.character(subject_ids)
    lower <- priorities[-1L]
    event_times <- matrix(Inf, nrow = length(subject_ids), ncol = length(lower),
                          dimnames = list(as.character(subject_ids), as.character(lower)))
    rec_ids <- subject_ids[FALSE]
    rec_codes <- rec_events <- integer()
    for (k in seq_along(lower)) {
        for (i in seq_along(subject_ids)) {
            rows <- groups[[i]][codes[groups[[i]]] == lower[k]]
            if (length(rows)) event_times[i, k] <- min(times[rows])
            if (length(rows) > 1L) {
                rec_ids <- c(rec_ids, subject_ids[i])
                rec_codes <- c(rec_codes, lower[k])
                rec_events <- c(rec_events, length(rows))
            }
        }
    }
    recurrent_counts <- data.frame(id = rec_ids, priority = rec_codes,
                                    records = rec_events,
                                    discarded = rec_events - 1L,
                                    row.names = NULL)
    if (nrow(recurrent_counts)) {
        msg <- paste0("Repeated nonfatal records were found for ",
                      nrow(recurrent_counts), " individual/priority combination(s). ")
        if (recurrent == "error") {
            stop(msg, "The model uses time to the first event of each type. Set recurrent = 'first' to use the first occurrence.",
                 call. = FALSE)
        }
        warning(msg, "Only the first occurrence of each type is used; this is not a recurrent-event comparison rule.",
                call. = FALSE)
    }
    outcome <- .pim_model_matrix(outcome_frame, "Outcome")
    censor <- .pim_model_matrix(censor_frame, "Censoring", allow_empty = TRUE)
    rownames(outcome$matrix) <- rownames(censor$matrix) <- as.character(subject_ids)
    ordered <- order(subject_index, times, ifelse(codes == 0, Inf, codes))
    data <- data[ordered, , drop = FALSE]
    rownames(data) <- NULL
    if (length(omitted)) {
        warning("Removed ", length(omitted),
                " entire individual(s) because required data were missing.",
                call. = FALSE)
    }
    list(X = outcome$matrix, W = censor$matrix, u = u, death = death,
         event_times = event_times, priorities = priorities,
         id = stats::setNames(subject_ids, as.character(subject_ids)),
         subjects = subjects, data = data, formula = formula,
         censor_formula = censor_formula, terms = outcome$terms,
         censor_terms = censor$terms, contrasts = outcome$contrasts,
         censor_contrasts = censor$contrasts, xlevels = outcome$xlevels,
         censor_xlevels = censor$xlevels, n = length(subject_ids),
         n_rows = nrow(data), n_original = original_n,
         omitted_ids = omitted, recurrent_counts = recurrent_counts,
         id_column = id, response_variables = rr$variables,
         na.action = na_method, recurrent = recurrent)
}

#' @name hfaction
#' @aliases hfaction
#' @docType data
#' @title Synthetic Heart Failure Data
#' @description Synthetic death and hospitalization records in the clinical setting of Heart Failure: A Controlled Trial Investigating Outcomes of Exercise Training (HF-ACTION). These data contain no original trial participant records and do not reproduce trial estimates.
#' @format A data frame with 407 rows for 300 subjects. Each subject has one terminal death or censoring record and a separate hospitalization record when observed.
#' \describe{
#' \item{id}{Subject identifier.}
#' \item{time}{Observed time in years.}
#' \item{status}{0 for censoring or administrative end of follow-up, 1 for death, and 2 for first hospitalization.}
#' \item{treatment}{Factor with levels \code{"Usual care"} and \code{"Exercise"}.}
#' \item{age}{Age in years.}
#' \item{sex}{Factor with levels \code{"Female"} and \code{"Male"}.}
#' \item{lvef}{Left ventricular ejection fraction, in percentage points.}
#' \item{bmi}{Body mass index.}
#' \item{nyha}{New York Heart Association functional class, with levels \code{"II"} and \code{"III or IV"}.}
#' }
#' @details Baseline covariates are repeated consistently within a subject. Follow-up continues after hospitalization. An absent hospitalization row means that hospitalization was not observed before the terminal record. The data are identified as synthetic by an attribute.
#' @seealso \code{\link{gpim}}, \code{\link{hierarchy}}
#' @keywords datasets
"hfaction"

#' @name hierarchy
#' @aliases hierarchy
#' @docType data
#' @title Synthetic Data with Three Endpoint Priorities
#' @description Synthetic observed event histories with a fatal event and two types of nonfatal events. These data contain no clinical participant records.
#' @format A data frame for 180 subjects, with event records and one terminal record per subject.
#' \describe{
#' \item{id}{Subject identifier.}
#' \item{time}{Observed event or terminal follow-up time in a common time unit.}
#' \item{status}{0 for censoring, 1 for the fatal event, 2 for the next priority, and 3 for the lowest priority.}
#' \item{x1}{Continuous baseline covariate.}
#' \item{b}{Binary baseline covariate.}
#' \item{w}{Centered binary covariate for the censoring model.}
#' }
#' @details The default comparison order is 1, 2, 3. Each nonfatal type refers to its first observed event. A later event time is favorable. An endpoint is considered only after every higher priority has been established as tied.
#' @seealso \code{\link{gpim}}, \code{\link{hfaction}}
#' @keywords datasets
"hierarchy"
