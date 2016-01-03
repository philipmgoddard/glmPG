#' S4 class definition for a linear model
#'
#' @slot model character vector description of model
#' @slot intercept is an intercept included in the model? Takes values TRUE or FALSE
#' @slot outcome numeric vector of the outcome ('y'). Should be binary
#' @slot covariates data frame of the covariates. At present must be numeric
#' @slot modelCoef matrix representation of the fitted model coefficients ('betahats')
#' @export
setClass(
  Class = "GLMP",
  slots = list(
    model = "character",
    intercept = "logical",
    outcome = "numeric",
    covariates = "data.frame",
    modelCoef = "matrix")
)

#' Initialiser for GLMP objects
#' @param .Object object of class GLMP
#' @param model character string representing model
#' @param intercept logical. is intercept included?
#' @param outcome numerical vector of the outcome
#' @param covariates data frame holding numerical covariates
#' @param modelCoef matrix of the model coefficients
#' @export
setMethod(
  f = "initialize",
  signature = "GLMP",
  definition = function(.Object, model, intercept, outcome, covariates,
                        modelCoef) {

    colnames(modelCoef) <- c("betaHat", "stdErr", "z", "Pr(>|z|)")
    ifelse(intercept,
           rownames(modelCoef) <- c("Intercept", names(covariates)),
           rownames(modelCoef) <- names(covariates))

    .Object@model <- model
    .Object@intercept <- intercept
    .Object@outcome <- outcome
    .Object@covariates <- covariates
    .Object@modelCoef <- modelCoef
    return(.Object)
  }
)

#' Constructor for GLMP
#' @param outcome character vector name of column containing outcome ("y")
#' @param covariates character vector of column names containing covariates ('x')
#' @param int logical. is an intercept included?
#' @param df data frame
#' @import MASS
#' @import methods
#' @import stats
#' @export
glmp <- function(outcome, covariates, int = TRUE, df) {

  y <- df[, names(df) == outcome]
  X <- df[, names(df) %in% covariates, drop = FALSE]
  if (int) X <- cbind(1, X)
  nVar <- ncol(X)
  X <- data.matrix(X)

  # initial values and solve for beta coefficients
  # iteratively
  init <- betaInit(X, y)
  beta <- optim(init,
                 fn = LL,
                 X = X,
                 y = y,
                 method = 'Nelder-Mead',
                 hessian = TRUE)$par

  stdErr <- sqrt(diag(varBeta(X, Wmat(X, beta))))
  z <- beta / stdErr
  pVals <- pnorm(abs(z), lower.tail = FALSE) * 2

  return(new(Class = "GLMP",
             model = (paste0(paste0(outcome, " ~ "), paste(covariates, collapse = " + "))),
             intercept = int,
             outcome = y,
             covariates = df[, names(df) %in% covariates, drop = FALSE],
             modelCoef = cbind(beta, stdErr, z, pVals)))
}
