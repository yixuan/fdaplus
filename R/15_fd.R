#' Univariate Functional Data Object
#'
#' A univariate functional data object, of the class \code{fd+}, can be
#' thought of as a vector whose elements are univariate functions
#' that can be expressed as linear combination of a system of basis functions.
#' Functions in one object are assumed to share the same basis, so an
#' \code{fd+} object is composed of one \code{\link[=basis+-class]{basis+}}
#' object and the corresponding coefficients.
#'
#' @slot basis Basis object of class \code{\link[=basis+-class]{basis+}}.
#' @slot coefs A matrix of dimension \code{n} by \code{p} where \code{n}
#'             is the number of functions in this object, and \code{p}
#'             is the number of basis functions in \code{basis}.
#'             The \code{i}-th row of \code{coefs} contains the coefficients
#'             for the \code{i}-th function in the object.
#'
#' @export
setClass("fd+", slots = c(coefs = "matrix", basis = "basis+"),
         validity = function(object) {
             if(!nrow(object@coefs))
                 return("fd+ object must contain at least one function")
             if(ncol(object@coefs) != object@basis@ncoef)
                 return("ncol(coefs) must be equal to the number of basis functions")
             return(TRUE)
         }
)



#' Creating A Univariate Functional Data Object
#'
#' This function constructs an \code{\link[=fd+-class]{fd+}} object
#' that represents a set of univariate functional data.
#'
#' @param coefs A matrix whose each row gives the coefficients of basis functions
#'              for each curve in this object. If \code{coefs} is a vector, it
#'              is treated as a matrix of one row.
#' @param basis Basis of this functional data object. Should be of class
#'              \code{\link[=basis+-class]{basis+}} or one of its subclasses
#'              (\code{\link[=bspline+-class]{bspline+}} and
#'              \code{\link[=fourier+-class]{fourier+}}).
#'
#' @return An \code{\link[=fd+-class]{fd+}} object with the given basis and coefficients.
#' @author Yixuan Qiu <\url{http://statr.me/}>
#' @export
fd_new = function(coefs, basis)
{
    if(is.matrix(coefs))
    {
        return(new("fd+", coefs = coefs, basis = basis))
    } else if(is.numeric(coefs)) {
        coefs = as.numeric(coefs)
        return(new("fd+", coefs = t(coefs), basis = basis))
    } else {
        stop("coefs must a matrix or a vector")
    }
}

#' @describeIn wrap Converting "fd" objects
wrap.fd = function(obj, ...)
{
    new("fd+", coefs = t(obj$coefs), basis = wrap(obj$basis))
}

#' @rdname subsetter-methods
setMethod("[",
          signature(x = "fd+", i = "numeric", j = "missing", drop = "ANY"),
          function(x, i) {
              i = as.integer(i)
              newcoefs = x@coefs[i, , drop = FALSE]
              initialize(x, coefs = newcoefs)
          }
)

#' Combining Multiple Univariate Functional Data Objects
#'
#' This function combines a series of \code{\link[=fd+-class]{fd+}} objects
#' that share the same basis into a single one.
#'
#' @usage \method{c}{`fd+`}(x, ...)
#' @param x,\dots \code{\link[=fd+-class]{fd+}} objects to be combined. They should
#'                have the same basis functions.
#'
#' @return An \code{\link[=fd+-class]{fd+}} object that contains all the functions
#'         in the list.
#' @author Yixuan Qiu <\url{http://statr.me/}>
#' @export
`c.fd+` = function(x, ...)
{
    objs = list(x, ...)
    if(length(objs) < 2)  return(x)
    for(i in 2:length(objs))
    {
        if(!is(objs[[i]], "fd+"))
            stop("all elements must be fd+ objects")
        if(!identical(objs[[i]]@basis, objs[[i - 1]]@basis))
            stop("basis of each fd+ object must be the same")
    }
    newcoefs = do.call(rbind, lapply(objs, function(x) x@coefs))
    initialize(x, coefs = newcoefs)
}

#' @rdname feval-methods
#'
#' @section Method (fd+, numeric):
#' \tabular{lcl}{
#'   \code{f}  \tab - \tab  A \code{\link[=fd+-class]{fd+}} object. \cr
#'   \code{x}  \tab - \tab  A numeric vector.
#' }
#'
#' \code{feval(f, x)} returns a matrix \code{R} of
#' \code{n} rows and \code{length(x)} columns where \code{n} is the number of
#' function curves in \code{f}. \code{R[i, j]}
#' equals the value of the \code{i}-th function evaluated on \code{x[j]}.
setMethod("feval", signature(f = "fd+", x = "numeric"),
          function(f, x, ...) {
              f@coefs %*% feval(f@basis, x)
          }
)

#' @rdname plot-methods
setMethod("plot", signature(x = "fd+", y = "missing"),
          function(x, y, ...) {
              x0 = seq(x@basis@range[1], x@basis@range[2],
                       length.out = 101)
              args = list(...)
              if(!"type" %in% names(args))
                  args = c(args, type = "l")
              if(!"xlab" %in% names(args))
                  args = c(args, xlab = "t")
              if(!"ylab" %in% names(args))
                  args = c(args, ylab = "Functional data")
              y = feval(x, x0)
              args = c(list(x = x0, y = t(y)), args)
              do.call(graphics::matplot, args)
          }
)

#' @rdname inprod-methods
#'
#' @section Method (univariate function vs univariate function):
#' \code{x} and \code{y} are two \code{\link[=fd+-class]{fd+}} objects
#' which we do not require to have the same types of basis functions.
#'
#' Assume that \code{x} contains \code{m} functions and \code{y} contains
#' \code{n} functions, and then \code{x \%*\% y} returns an \code{m} by \code{n}
#' matrix \code{P}, whose element \code{P[i, j]} is the inner product between
#' the \code{i}-th function of \code{x} and the \code{j}-th function of \code{y}.
setMethod("%*%", signature(x = "fd+", y = "fd+"),
          function(x, y) {
              x@coefs %*% (x@basis %*% y@basis) %*% t(y@coefs)
          }
)

#' @rdname inprod-methods
#'
#' @section Method (univariate function vs basis):
#' Assueme that \code{x} is an \code{\link[=fd+-class]{fd+}} object and \code{y}
#' a \code{\link[=basis+-class]{basis+}} object (or vice versa). We do not require
#' \code{x} to have the same basis as \code{y}.
#'
#' If \code{x} contains \code{m} functions and \code{y} contains
#' \code{n} functions, then \code{x \%*\% y} returns an \code{m} by \code{n}
#' matrix \code{P}, whose element \code{P[i, j]} is the inner product between
#' the \code{i}-th function of \code{x} and the \code{j}-th function of \code{y}.
setMethod("%*%", signature(x = "fd+", y = "basis+"),
          function(x, y) {
              x@coefs %*% (x@basis %*% y)
          }
)

#' @rdname inprod-methods
setMethod("%*%", signature(x = "basis+", y = "fd+"),
          function(x, y) {
              (x %*% y@basis) %*% t(y@coefs)
          }
)

## Arithmetic between fd+ and a scalar
#' @rdname arithmetic-methods
#'
#' @section Method (fd+, numeric scalar):
#' An \code{\link[=fd+-class]{fd+}} object can be multiplied or divided by
#' a numeric scalar, which has the effect that all the functions in this object
#' will be scaled by this constant.
setMethod("*", signature(e1 = "fd+", e2 = "numeric"),
          function(e1, e2) {
              initialize(e1, coefs = e2[1] * e1@coefs)
          }
)
#' @rdname arithmetic-methods
setMethod("*", signature(e1 = "numeric", e2 = "fd+"),
          function(e1, e2) {
              initialize(e2, coefs = e1[1] * e2@coefs)
          }
)
#' @rdname arithmetic-methods
setMethod("/", signature(e1 = "fd+", e2 = "numeric"),
          function(e1, e2) {
              initialize(e1, coefs = e1@coefs / e2[1])
          }
)

## Arithmetic between fd+ objects
#' @rdname arithmetic-methods
#'
#' @section Method (fd+, fd+):
#' An \code{\link[=fd+-class]{fd+}} object can be added to or subtracted from
#' another \code{\link[=fd+-class]{fd+}} object, if they have the same basis,
#' and:
#'
#' \itemize{
#'   \item One of them contains only one function curve, so that this single
#'         function will be added to or subtracted from all the functions
#'         in the other object.
#'   \item Or, they have the same number of functions, so that addition and
#'         subtraction will be done element-wisely.
#' }
setMethod("+", signature(e1 = "fd+", e2 = "fd+"),
          function(e1, e2) {
              if(!identical(e1@basis, e2@basis))
                  stop("need to have the same basis functions");
              if(nrow(e1@coefs) != nrow(e2@coefs))
                  stop("need to contain the same number of functions")
              initialize(e1, coefs = e1@coefs + e2@coefs)
          }
)
#' @rdname arithmetic-methods
setMethod("-", signature(e1 = "fd+", e2 = "fd+"),
          function(e1, e2) {
              if(!identical(e1@basis, e2@basis))
                  stop("need to have the same basis functions");
              if(nrow(e1@coefs) != nrow(e2@coefs))
                  stop("need to contain the same number of functions")
              initialize(e1, coefs = e1@coefs - e2@coefs)
          }
)

## Calculate mean function
setMethod("mean", signature(x = "fd+"),
          function(x, ...) {
              initialize(x, coefs = t(colMeans(x@coefs)))
          }
)

## Calculate sd and var function
## Seems difficult to implement

## Calculate bivariate covariance function
setMethod("cov", signature(x = "fd+", y = "missing",
                           use = "ANY", method = "ANY"),
          function(x, y, use, method) {
              if(nrow(x@coefs) <= ncol(x@coefs))
                  stop("number of functions should be greater than the number of basis functions")
              coefs = var(x@coefs)
              new("bifd+", coefs = coefs, sbasis = x@basis, tbasis = x@basis)
          }
)
