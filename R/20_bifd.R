#' Bivariate Function Using Basis Expansion
#'
#' This class defines bivariate functions that can be expressed by two sets
#' of basis functions and the associated coefficient matrix. It takes the following
#' form:
#'
#' \deqn{K(s,t)=\sum_{i,j} a_{ij}f_i(s)g_j(t)}{K(s, t) = sum_{i,j} a_{ij} * f_i(s) * g_j(t)}
#'
#' Here \eqn{K(s, t)} is the bivariate function, \eqn{f_i} and \eqn{g_j} are
#' two basis systems, and \eqn{A = (a_{ij})} is the coefficient matrix.
#'
#' @slot sbasis,tbasis Basis objects of class \code{\link[=basis+-class]{basis+}},
#'                     not necessarily of the same type (for example, one can be
#'                     \code{\link[=bspline+-class]{bspline+}} and the other be
#'                     \code{\link[=fourier+-class]{fourier+}}).
#' @slot coefs A matrix of dimension \code{m} by \code{n} where \code{m}
#'             is the number of functions in \code{sbasis}, and \code{n}
#'             is the number of functions in \code{tbasis}.
#'
#' @export
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



#' Creating A Bivariate Function Using Basis Expansion
#'
#' This function constructs a \code{\link[=bifd+-class]{bifd+}} object
#' that represents a bivariate function.
#'
#' @slot coefs The coefficient matrix.
#' @slot sbasis,tbasis Basis objects of class \code{\link[=basis+-class]{basis+}},
#'                     not necessarily of the same type (for example, one can be
#'                     \code{\link[=bspline+-class]{bspline+}} and the other be
#'                     \code{\link[=fourier+-class]{fourier+}}).
#'
#' @return An \code{\link[=bifd+-class]{bifd+}} object with the given
#'         bases and coefficients.
#' @author Yixuan Qiu <\url{http://statr.me/}>
#' @export
bifd_new = function(coefs, sbasis, tbasis = sbasis)
{
    new("bifd+", coefs = coefs, sbasis = sbasis, tbasis = tbasis)
}

#' @describeIn wrap Converting "bifd" objects
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

#' @rdname feval-methods
#'
#' @section Method (bifd+, numeric, numeric):
#' \tabular{lcl}{
#'   \code{f}     \tab - \tab  A \code{\link[=bifd+-class]{bifd+}} object. \cr
#'   \code{x, y}  \tab - \tab  Numeric vectors.
#' }
#'
#' \code{feval(f, x, y)} returns a matrix \code{R} of
#' \code{length(x)} rows and \code{length(y)} columns, with \code{R[i, j]}
#' equal to the value of \code{f(x[i], x[j])}.
setMethod("feval", signature(f = "bifd+", x = "numeric"),
          function(f, x, y, ...) {
              y = as.numeric(y)
              crossprod(feval(f@sbasis, x), f@coefs) %*% feval(f@tbasis, y)
          }
)

#' @rdname plot-methods
setMethod("plot", signature(x = "bifd+", y = "missing"),
          function(x, y, ..., engine = c("graphics", "rgl")) {
              x0 = seq(x@sbasis@range[1], x@sbasis@range[2],
                       length.out = 101)
              y0 = seq(x@tbasis@range[1], x@tbasis@range[2],
                       length.out = 101)
              z = feval(x, x0, y0)

              if(engine[1] == "rgl")
                  rgl::persp3d(x0, y0, z, ...)
              else
                  persp(x0, y0, z, ...)
          }
)

## Integral transform (product) of two bivariate functions
setMethod("%*%", signature(x = "bifd+", y = "bifd+"),
          function(x, y) {
              newcoef = x@coefs %*% (x@tbasis %*% y@sbasis) %*% y@coefs
              new("bifd+", coefs = newcoef, sbasis = x@sbasis,
                  tbasis = y@tbasis)
          }
)

## Integral transform (product) on fd+
setMethod("%*%", signature(x = "bifd+", y = "fd+"),
          function(x, y) {
              newcoef = x@coefs %*% (x@tbasis %*% y@basis) %*% t(y@coefs)
              new("fd+", coefs = t(newcoef), basis = x@sbasis)
          }
)
setMethod("%*%", signature(x = "fd+", y = "bifd+"),
          function(x, y) {
              newcoef = x@coefs %*% (x@basis %*% y@sbasis) %*% y@coefs
              new("fd+", coefs = newcoef, basis = y@tbasis)
          }
)

## Square root of a symmetric matrix
sqrtm = function(x)
{
    if(!isSymmetric(x))
        stop("x must be a symmetric matrix")
    .Call("sqrtm", x, PACKAGE = "fdaplus")
}
## Inverse of the square root of a symmetric matrix
sqrtInvm = function(x)
{
    if(!isSymmetric(x))
        stop("x must be a symmetric matrix")
    .Call("sqrtInvm", x, PACKAGE = "fdaplus")
}
## Both of above
sqrtBothm = function(x)
{
    if(!isSymmetric(x))
        stop("x must be a symmetric matrix")
    .Call("sqrtBothm", x, PACKAGE = "fdaplus")
}
## Power of (symmetric) bivariate function
power_bifd = function(x, k)
{
    if(!isSymmetric(x@coefs) |
       !identical(x@sbasis, x@tbasis))
        stop("need a symmetric bivariate function")

    xmat = x@coefs
    w = penmat(x@sbasis, 0)
    wsqrtBoth = sqrtBothm(w)
    wsqrt = wsqrtBoth$sqrt
    wsqrtInv = wsqrtBoth$sqrtInv

    mdecomp = wsqrt %*% xmat %*% wsqrt
    e = eigen(mdecomp)

    newcoef = wsqrtInv %*% e$vectors %*% diag(e$values^k) %*% t(e$vectors) %*% wsqrtInv
    initialize(x, coefs = newcoef)
}
setMethod("^", signature(e1 = "bifd+", e2 = "numeric"),
          function(e1, e2) power_bifd(e1, e2)
)
