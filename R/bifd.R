## coefs is an p1 by p2 matrix
## p1 == sbasis@ncoef
## p2 == tbasis@ncoef
setClass("bifd+", slots = c(coefs = "matrix",
                            sbasis = "basis+",
                            tbasis = "basis+"),
         validity = function(object) {
             if(nrow(object@coefs) != object@sbasis@ncoef)
                 return("nrow(coefs) must be equal to sbasis@ncoef")
             if(ncol(object@coefs) != object@tbasis@ncoef)
                 return("ncol(coefs) must be equal to tbasis@ncoef")
             return(TRUE)
         }
)



bifd_new = function(coefs, sbasis, tbasis = sbasis)
{
    new("bifd+", coefs = coefs, sbasis = sbasis, tbasis = tbasis)
}

wrap.bifd = function(obj, ...)
{
    new("bifd+", coefs = obj$coefs, sbasis = wrap(obj$sbasis),
        tbasis = wrap(obj$tbasis))
}

## Arithmetic between bifd+ and a scalar
setMethod("*", signature(e1 = "bifd+", e2 = "numeric"),
          function(e1, e2) {
              initialize(e1, coefs = e2 * e1@coefs)
          }
)
setMethod("*", signature(e1 = "numeric", e2 = "bifd+"),
          function(e1, e2) {
              initialize(e2, coefs = e1 * e2@coefs)
          }
)
setMethod("/", signature(e1 = "bifd+", e2 = "numeric"),
          function(e1, e2) {
              initialize(e1, coefs = e1@coefs / e2)
          }
)

## Arithmetic between bifd+ objects
setMethod("+", signature(e1 = "bifd+", e2 = "bifd+"),
          function(e1, e2) {
              if(!identical(e1@sbasis, e2@sbasis) |
                 !identical(e1@tbasis, e2@tbasis))
                  stop("need to have the same basis functions");
              initialize(e1, coefs = e1@coefs + e2@coefs)
          }
)
setMethod("-", signature(e1 = "bifd+", e2 = "bifd+"),
          function(e1, e2) {
              if(!identical(e1@sbasis, e2@sbasis) |
                 !identical(e1@tbasis, e2@tbasis))
                  stop("need to have the same basis functions");
              initialize(e1, coefs = e1@coefs - e2@coefs)
          }
)

## A generic implementation of feval() for bifd+ class
## Will call feval() on the basis
## Return a matrix K(i, j) with
##          K(i, j) = f(x[i], y[j])
setMethod("feval", signature(f = "bifd+", x = "numeric"),
          function(f, x, y, ...) {
              y = as.numeric(y)
              crossprod(feval(f@sbasis, x), f@coefs) %*% feval(f@tbasis, y)
          }
)

## A generic implementation of plot() for bifd+ class
## Will call feval() on bifd+
setMethod("plot", signature(x = "bifd+", y = "missing"),
          function(x, y, ...) {
              x0 = seq(x@sbasis@range[1], x@sbasis@range[2],
                       length.out = 101)
              y0 = seq(x@tbasis@range[1], x@tbasis@range[2],
                       length.out = 101)
              z = feval(x, x0, y0)
              persp(x0, y0, z, ...)
          }
)
setMethod("plot3d", signature(x = "bifd+"),
          function(x, ...) {
              x0 = seq(x@sbasis@range[1], x@sbasis@range[2],
                       length.out = 101)
              y0 = seq(x@tbasis@range[1], x@tbasis@range[2],
                       length.out = 101)
              z = feval(x, x0, y0)
              rgl::persp3d(x0, y0, z, ...)
          }
)