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
