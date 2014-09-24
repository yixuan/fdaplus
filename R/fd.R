setClass("fd+", slots = c(coefs = "matrix", basis = "basis+"),
         validity = function(object) {
             if(ncol(object@coefs) != object@basis@ncoef)
                 return("ncol(coefs) must be equal to the number of basis functions")
             return(TRUE)
         }
)

## A generic implementation of feval() for fd+ class
## Will call feval() on the basis
setMethod("feval", signature(f = "fd+", x = "numeric"),
          function(f, x, ...) {
              feval(f@basis, x) %*% t(f@coefs)
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
              args = c(list(x = x0, y = y), args)
              do.call(graphics::matplot, args)
          }
)