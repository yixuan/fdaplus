#' Basis Functions for Functional Data
#'
#' The "\code{basis+}" class is the base class for all basis function objects,
#' including \code{\link[=bspline+-class]{bspline+}} and
#' \code{\link[=fourier+-class]{fourier+}}.
#'
#' @slot range A length-two numeric vector to define the interval on which
#'             basis functions can be evaluated.
#' @slot nbasis The total number of basis functions, including the ones that need
#'              to be dropped.
#' @slot dropind Indices of basis functions that need to be dropped.
#' @slot ncoef The actual number of basis functions to represent functional data,
#'             satisfying the condition \code{ncoef == nbasis - length(dropind)}.
#'
#' @seealso \code{\link[=bspline+-class]{bspline+}},
#' \code{\link[=fourier+-class]{fourier+}}.
#'
#' @export
setClass("basis+", slots = c(range = "numeric",
                             nbasis = "numeric",
                             dropind = "numeric",
                             ncoef = "numeric"),
         validity = function(object) {
             if(length(object@range) != 2)
                 return("range must be a vector of length 2")
             if(object@nbasis <= 0)
                 return("nbasis must be positive")
             if(object@ncoef != object@nbasis - length(object@dropind))
                 return("ncoef must be equal to nbasis-length(dropind)")
             if(object@ncoef <= 0)
                 return("number of dropped indices must be smaller than nbasis")
             return(TRUE)
         }
)

#' B-spline Basis Functions
#'
#' The "\code{bspline+}" class represents B-spline basis functions.
#'
#' @section Extends:
#' Class \code{\link[=basis+-class]{basis+}}
#'
#' @slot range A length-two numeric vector to define the interval on which
#'             basis functions can be evaluated. Default is \code{c(0, 1)}.
#' @slot nbasis The total number of basis functions, including the ones that need
#'              to be dropped. Default is 4.
#' @slot dropind Indices of basis functions that need to be dropped. Default is
#'               \code{numeric(0)} (no basis will be dropped).
#' @slot ncoef The actual number of basis functions to represent functional data,
#'             satisfying the condition \code{ncoef == nbasis - length(dropind)}.
#' @slot degree Degree of the B-splines. Default is 3, standing for cubic splines.
#' @slot knots A numeric vector giving the inner knots. Must satisfy
#'             \code{nbasis == degree + length(knots) + 1}.
#'
#' @export
setClass("bspline+", contains = "basis+",
         slots = c(degree = "numeric", knots = "numeric"),
         prototype = list(range = c(0, 1),
                          nbasis = 4,
                          dropind = numeric(0),
                          ncoef = 4,
                          degree = 3,
                          knots = numeric(0)),
         validity = function(object) {
             if(object@nbasis != object@degree + length(object@knots) + 1)
                 return("nbasis must be equal to degree+length(knots)+1")
             return(TRUE)
         }
)

#' Fourier Basis Functions
#'
#' The "\code{fourier+}" class represents Fourier basis functions.
#'
#' @section Extends:
#' Class \code{\link[=basis+-class]{basis+}}
#'
#' @slot range A length-two numeric vector to define the interval on which
#'             basis functions can be evaluated. Default is \code{c(0, 1)}.
#' @slot nbasis The total number of basis functions, including the ones that need
#'              to be dropped. Default is 3.
#' @slot dropind Indices of basis functions that need to be dropped. Default is
#'               \code{numeric(0)} (no basis will be dropped).
#' @slot ncoef The actual number of basis functions to represent functional data,
#'             satisfying the condition \code{ncoef == nbasis - length(dropind)}.
#' @slot period The length of the cycle that Fourier basis functions repeat
#'              themselves. Default is 1.
#'
#' @export
setClass("fourier+", contains = "basis+",
         slots = c(period = "numeric"),
         prototype = list(range = c(0, 1),
                          nbasis = 3,
                          dropind = numeric(0),
                          ncoef = 3,
                          period = 1),
         validity = function(object) {
             if(object@period <= 0)
                 return("period must be positive")
             return(TRUE)
         }
)



#' @describeIn wrap Converting "basisfd" objects
wrap.basisfd = function(obj, ...)
{
    if(obj$type == "bspline")
    {
        res = new("bspline+",
                  range = obj$rangeval,
                  nbasis = obj$nbasis,
                  dropind = as.numeric(obj$dropind),
                  ncoef = obj$nbasis - length(obj$dropind),
                  degree = obj$nbasis - length(obj$params) - 1,
                  knots = as.numeric(obj$params))
    } else if(obj$type == "fourier") {
        res = new("fourier+",
                  range = obj$rangeval,
                  nbasis = obj$nbasis,
                  dropind = as.numeric(obj$dropind),
                  ncoef = obj$nbasis - length(obj$dropind),
                  period = obj$params)
    } else stop("type not supported yet")
    res
}

#' @rdname subsetter-methods
setMethod("[",
          signature(x = "basis+", i = "numeric", j = "missing", drop = "ANY"),
          function(x, i) {
              i = as.integer(i)
              if(any(i > x@ncoef))
                  stop("subscript out of bounds")
              orgind = setdiff(1:x@nbasis, x@dropind)
              newind = orgind[i]
              newdropind = setdiff(1:x@nbasis, newind)
              initialize(x, dropind = newdropind,
                         ncoef = x@nbasis - length(newdropind))
          }
)

## A generic implementation of plot() for basis+ class
## Calls feval()
setMethod("plot", signature(x = "basis+", y = "missing"),
          function(x, y, ...) {
              x0 = seq(x@range[1], x@range[2], length.out = 101)
              args = list(...)
              if(!"type" %in% names(args))
                  args = c(args, type = "l")
              if(!"xlab" %in% names(args))
                  args = c(args, xlab = "t")
              if(!"ylab" %in% names(args))
                  args = c(args, ylab = "Basis function")
              y = feval(x, x0)
              args = c(list(x = x0, y = t(y)), args)
              do.call(graphics::matplot, args)
          }
)

## A generic implementation of inner product for basis+ class
## Calls feval()
setMethod("%*%", signature(x = "basis+", y = "basis+"),
          function(x, y) {
              if(!isTRUE(all.equal(x@range, y@range)))
                  stop("range of x and y must be the same")

              res = .Call("basis_inprod", x, y)
              dim(res) = c(x@ncoef, y@ncoef)
              res
          }
)
