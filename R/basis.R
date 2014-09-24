setClass("basis+", slots = c(range = "numeric",
                             nbasis = "numeric",
                             dropind = "numeric",
                             ncoef = "numeric")
)

setClass("bspline+", contains = "basis+",
         slots = c(degree = "numeric", knots = "numeric"),
         prototype = list(range = c(0, 1),
                          nbasis = 4,
                          dropind = numeric(0),
                          ncoef = 4,
                          degree = 3,
                          knots = numeric(0)),
         validity = function(object) {
             if(length(object@range) != 2)
                 return("range must be a vector of length 2")
             if(object@nbasis != object@degree + length(object@knots) + 1)
                 return("nbasis must be equal to degree+length(knots)+1")
             if(object@ncoef != object@nbasis - length(object@dropind))
                 return("ncoef must be equal to nbasis-length(dropind)")
             return(TRUE)
         }
)

## Wrap objects in fda to the new class
wrap = function(obj, ...)
    UseMethod("wrap")


## Evaluate functional data objects
## *** For basis+ class, it returns an p by T matrix
##     p is the number of basis functions, T is the length of x
## *** For fd+ class, it returns an n by T matrix
##     n is the number of functions, T is the length of x
if(!isGeneric("feval"))
    setGeneric("feval", function(f, x, ...) standardGeneric("feval"))

## Plot functional data objects
if(!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

## Inner product of functional data objects
## *** For basis+ class, it returns an p by p matrix
##     p is the number of basis functions
##     P(i, j) = integrate B_i(x) * B_j(x) dx
if(!isGeneric("%*%"))
    setGeneric("%*%")

## Penalty matrix of basis function
## P(i, j) = integrate B^(k)_i(x) * B^(k)_j(x) dx
## k is the order of derivative, given by the "penalty" argument
if(!isGeneric("penmat"))
    setGeneric("penmat", function(basis, penalty, ...)
        standardGeneric("penmat"))


## A generic implementation of plot() for basis+ class
## Will call feval()
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
              if(length(x@dropind))
                  y = y[, -x@dropind]
              args = c(list(x = x0, y = t(y)), args)
              do.call(graphics::matplot, args)
          }
)
