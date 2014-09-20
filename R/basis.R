setClass("basis+", slots = c(rangeval = "numeric",
                             nbasis = "numeric",
                             dropind = "numeric")
)

setClass("bspline+", contains = "basis+",
         slots = c(knots = "numeric"),
         prototype = list(rangeval = c(0, 1),
                          nbasis = 4,
                          dropind = numeric(0),
                          knots = numeric(0)),
         validity = function(object) {
             TRUE
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
