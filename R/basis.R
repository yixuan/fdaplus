setClass("basis+", slots = c(range = "numeric",
                             nbasis = "numeric",
                             dropind = "numeric")
)

setClass("bspline+", contains = "basis+",
         slots = c(degree = "numeric", knots = "numeric"),
         prototype = list(range = c(0, 1),
                          nbasis = 4,
                          dropind = numeric(0),
                          degree = 3,
                          knots = numeric(0)),
         validity = function(object) {
             if(length(object@range) != 2)
                 return("range must be a vector of length 2")
             if(object@nbasis != object@degree + length(object@knots) + 1)
                 return("nbasis must be equal to degree+length(knots)+1")
             return(TRUE)
         }
)

## Wrap objects in fda to the new class
wrap = function(obj, ...)
    UseMethod("wrap")


## Evaluate functional data objects
if(!isGeneric("feval"))
    setGeneric("feval", function(f, x, ...) standardGeneric("feval"))

if(!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
