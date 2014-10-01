## coefs is an n by p matrix
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

## A generic implementation of "c" for fd+ class
`c.fd+` = function(x, ..., recursive = FALSE)
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
setMethod("/", signature(e1 = "fd+", e2 = "numeric"),
          function(e1, e2) {
              initialize(e1, coefs = e1@coefs / e2)
          }
)

## Arithmetic between fd+ objects
setMethod("+", signature(e1 = "fd+", e2 = "fd+"),
          function(e1, e2) {
              if(!identical(e1@basis, e2@basis))
                  stop("need to have the same basis functions");
              if(nrow(e1@coefs) != nrow(e2@coefs))
                  stop("need to contain the same number of functions")
              initialize(e1, coefs = e1@coefs + e2@coefs)
          }
)
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
setMethod("sd", signature(x = "fd+", na.rm = "ANY"),
          function(x, na.rm = FALSE) {
              initialize(x, coefs = t(apply(x@coefs, 2, sd)))
          }
)
setMethod("var", signature(x = "fd+", y = "missing", na.rm = "ANY", use = "ANY"),
          function(x, y = NULL, na.rm = FALSE, use) {
              initialize(x, coefs = t(apply(x@coefs, 2, var)))
          }
)

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
