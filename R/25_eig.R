## A thin wrapper of the eigen decomposition of bivariate function.
## I do this since I want to use "funs" as the field name
## rather than "harmonics", which is used by the fda package.
## However, it is helpful to be compatible with fda.
## By wrapping the result with an S4 object, I can map result$harmonics
## to result$functions, without introducing an extra field.
setClassUnion("fd+_or_NULL", c("fd+", "NULL"))
setClass("eig+", slots = c(values = "numeric",
                           funs = "fd+_or_NULL"))

## Eigen decomposition
eig_bifd = function(x, k, ...)
{
    if(!isSymmetric(x@coefs) |
           !identical(x@sbasis, x@tbasis))
        stop("need a symmetric bivariate function")
    
    cct = x@coefs
    wboth = sqrtBothm(penmat(x@sbasis, 0))
    wsqrt = wboth$sqrt
    wsqrtInv = wboth$sqrtInv
    mat = wsqrt %*% cct %*% wsqrt
    
    ## values only
    if(k < 1)
    {
        e = eigen(mat, only.values = TRUE)
        funs = NULL
        return(new("eig+", values = e$values, funs = NULL))
    }
    
    k = min(k, nrow(x@coefs))
    e = eigen(mat)
    b = crossprod(e$vectors[, 1:k], wsqrtInv)
    new("eig+",
        values = e$values,
        funs = fd_new(b, x@sbasis))
}

setMethod("eig", signature(x = "bifd+", k = "numeric"), eig_bifd)
setMethod("eig", signature(x = "bifd+", k = "missing"),
          function(x, k, ...) eig_bifd(x, nrow(x@coefs)))

setMethod("$", signature(x = "eig+"),
          function(x, name) {
              if(name == "harmonics")
                  return(x@funs)
              else
                  return(slot(x, name))
          }
)
