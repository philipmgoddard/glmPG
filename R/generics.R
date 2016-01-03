#' @include class-definition.R
NULL

#' Show generic
#' @param object object of class GLMP
#' @export
setMethod(
  f = "show",
  signature = "GLMP",
  function(object) {
    cat("\nCall: \n")
    cat(object@model)
    cat("\n\nCoefficients: \n \n")
    print(object@modelCoef)
  }
)
