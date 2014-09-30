## coefs is a n by p matrix
## n is the number of functions, p is the number of basis functions
setClass("fd+", slots = c(coefs = "matrix", basis = "basis+"),
         validity = function(object) {
             if(!nrow(object@coefs))
                 return("fd+ object must contain at least one function")
             if(ncol(object@coefs) != object@basis@ncoef)
                 return("ncol(coefs) must be equal to the number of basis functions")
             return(TRUE)
         }
)



wrap.fd = function(obj, ...)
{
    new("fd+", coefs = t(obj$coefs), basis = wrap(obj$basis))
}

## A generic implementation of "[" for fd+ class
setMethod("[",
          signature(x = "fd+", i = "numeric", j = "missing", drop = "ANY"),
          function(x, i, j, drop) {
              i = as.integer(i)
              newcoefs = x@coefs[i, , drop = FALSE]
              initialize(x, coefs = newcoefs)
          }
)

## A generic implementation of feval() for fd+ class
## Will call feval() on the basis
setMethod("feval", signature(f = "fd+", x = "numeric"),
          function(f, x, ...) {
              f@coefs %*% feval(f@basis, x)
          }
)

## A generic implementation of plot() for fd+ class
## Will call feval() on fd+
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

## A generic implementation of "%*%" for fd+ class
## Will call "%*%" on the basis
setMethod("%*%", signature(x = "fd+", y = "fd+"),
          function(x, y) {
              x@coefs %*% (x@basis %*% y@basis) %*% t(y@coefs)
          }
)

## Arithmetic between fd+ and a scalar
setMethod("*", signature(e1 = "fd+", e2 = "numeric"),
          function(e1, e2) {
              initialize(e1, coefs = e2 * e1@coefs)
          }
)
setMethod("*", signature(e1 = "numeric", e2 = "fd+"),
          function(e1, e2) {
              initialize(e2, coefs = e1 * e2@coefs)
          }
)

## Calculate mean function
setMethod("mean", signature(x = "fd+"),
          function(x, ...) {
              initialize(x, coefs = t(colMeans(x@coefs)))
          }
)

## Calculate sd and var function
setMethod("sd", signature(x = "fd+", na.rm = "ANY"),
          function(x, na.rm) {
              initialize(x, coefs = t(apply(x@coefs, 2, sd)))
          }
)
setMethod("var", signature(x = "fd+", y = "missing", na.rm = "ANY", use = "ANY"),
          function(x, y, na.rm, use) {
              initialize(x, coefs = t(apply(x@coefs, 2, var)))
          }
)
