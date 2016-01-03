#' Logit function
#' @param x input vector
#' @export
logit <- function(x) {
  return(log(x - (1 / x)))
}

#' Inverse logit
#' @param x input vector
#' @export
logitInv <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

#' liklihood function for binomial
#' @param beta vector of coefficients
#' @param X matrix of covariates
#' @param y dichotomous vector of outcomes (0, 1)
#' @export
LL <- function(beta, X, y) {
  probs <- logitInv(X %*% beta)
  return(-sum(y * log(probs) + (1 - y) * log(1 - probs)))
}

#' definition of beta hat. Not used in iterative process
#' @param X matrix of coefficients
#' @param W diagonal... [EXPLAIN]
#' @param z vector of... [EXPLAIN]
#' @export
betaHat <- function(X, W, z) {
  stp1 <- MASS::ginv(t(X) %*% W %*% X)
  stp2 <- stp1 %*% t(X) %*% W %*% z
  return(stp2)
}

#' Variance of coefficients
#' @param X matrix of coefficients
#' @param W diagonal matrix of... [EXPLAIN]
#' @export
varBeta <- function(X, W) {
  return(MASS::ginv(t(X) %*% W %*% X))
}

#' z vector
#' @param beta vector of coefficients
#' @param X matrix of covariates
#' @param y vecto of outcomes
#' @export
zval <- function(beta, X, y) {
  etaHat <- X %*% beta
  muHat <- logitInv(etaHat)
  z <- etaHat * ((y - muHat) / (muHat * (1 - muHat)))
  return(z)
}

#' W matrix
#' @param X matrix of covariates
#' @param beta vector of coefficients
#' @export
Wmat <- function(beta, X) {
  etaHat <- X %*% beta
  muHat <- logitInv(etaHat)
  w <- matrix(data = 0,
              nrow = length(muHat),
              ncol = length(muHat))
  diag(w) <- muHat * (1 - muHat)
  return(w)
}

#' Initialise z
#' @param y vector of outcomes
#' @export
zInit <- function(y) {
  z <- (y + 0.5) / (1 - y + 0.5)
  return(z)
}

#' Initial guess of beta vector
#' @param X matrix of covariates
#' @param y vector of outcomes
betaInit <- function(X, y) {
  return(betaHat(X,
                 diag(x = 1,
                      nrow = nrow(X),
                      ncol = nrow(X)),
                 zInit(y)))
}
