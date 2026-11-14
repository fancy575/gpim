#' Extract Subject Influence Functions
#'
#' @description Extracts the centered coefficient influence functions used
#' in the covariance calculation for CHW, RCW, or CTW.
#' @param model A fitted object returned by \code{\link{gpim}}.
#' @param method,link,L Optional selections of fitted methods, links, and
#' restriction times. \code{NULL} selects all fitted values.
#' @param components If \code{TRUE}, return the separate pairwise, censoring
#' baseline, censoring regression, and joint-model contributions.
#' @param ... Further arguments passed to or from methods.
#' @details Each matrix has one row per individual and one column per
#' regression coefficient. The contributions are on the coefficient scale,
#' include the inverse estimating-equation derivative, and are centered
#' columnwise. For a selected influence matrix \code{a}, the fitted covariance
#' is \code{crossprod(a) / nobs(model)^2}. The component matrices sum to
#' \code{a}. Their covariances must not be added separately, because that
#' would omit cross-covariances.
#' @return A matrix when one fit is selected, or a named list of matrices
#' for several fits. With \code{components = TRUE}, each selected fit returns
#' a list of its four component matrices.
#' @seealso \code{\link{gpim}}, \code{\link{pim_joint_fit}}
#' @export
influence.prior_pim <- function(model, method = NULL, link = NULL, L = NULL,
                                components = FALSE, ...) {
    if (!is.logical(components) || length(components) != 1L || is.na(components))
        stop("components must be TRUE or FALSE.", call. = FALSE)
    fits <- .pim_selected_fits(model, method, link, L)
    field <- if (components) "influence_components" else "influence"
    ans <- lapply(fits, function(f) {
        if (is.null(f[[field]]))
            stop("Influence functions are not available for the selected fit. ",
                 "Inspect pim_diagnostics(model).", call. = FALSE)
        f[[field]]
    })
    if (length(ans) == 1L) ans[[1L]] else ans
}
