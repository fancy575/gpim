# Inverse links and their derivatives for the pair estimating equations.
# The identity link retains the net-benefit scale used by this package.

.pim_link <- function(eta, link) {
  if (!is.character(link) || length(link) != 1L || is.na(link)) {
    stop("'link' must be a single link name.", call. = FALSE)
  }
  if (link == "identity") {
    return(list(mu = (1 + eta) / 2, complement = (1 - eta) / 2,
                prime = rep(0.5, length(eta)), multiplier = rep(1, length(eta)),
                multiplier_prime = rep(0, length(eta))))
  }
  if (link == "logit") {
    mu <- stats::plogis(eta)
    complement <- stats::plogis(-eta)
    return(list(mu = mu, complement = complement, prime = mu * complement,
                multiplier = rep(1, length(eta)),
                multiplier_prime = rep(0, length(eta))))
  }
  if (link == "probit") {
    mu <- stats::pnorm(eta)
    complement <- stats::pnorm(-eta)
    phi <- stats::dnorm(eta)
    cc <- exp(stats::dnorm(eta, log = TRUE) -
                stats::pnorm(eta, log.p = TRUE) -
                stats::pnorm(-eta, log.p = TRUE))
    return(list(mu = mu, complement = complement, prime = phi, multiplier = cc,
                multiplier_prime = cc * (-eta - cc * (complement - mu))))
  }
  stop("Unknown link. Choose identity, logit, or probit.",
       call. = FALSE)
}
