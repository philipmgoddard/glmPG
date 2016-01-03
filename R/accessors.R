#' @include class-definition.R
NULL

#' Set the generic for accessor (getter) for coefficients
#' @param object object of class GLMP
#' @export
setGeneric("getCoef",
           function(object){
             standardGeneric("getCoef")
           })

#' @describeIn getCoef
#' @export
setMethod("getCoef",
          signature = "GLMP",
          function(object){
            out <- object@modelCoef[, 1]
            return(out)
          })
