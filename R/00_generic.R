## Wrap objects in fda package to the new class
wrap = function(obj, ...)
    UseMethod("wrap")



## Retrieve a subset of functional data objects
## "[" is already defined in base package as a generic


## Retrieve a field
## "$" is already defined in base package as a generic


## Concatenate a sequence of functional data objects
## "c" is already a generic S3 method


#' Evaluating Various Types of Functional Data Objects
#'
#' @name feval
#' @usage feval(f, x, ...)
#'
#' ## For bifd+
#' # feval(f, x, y, ...)
#' @param f An object of class \code{\link[=basis+-class]{basis+}},
#'        \code{\link[=fd+-class]{fd+}}, or \code{\link[=bifd+-class]{bifd+}}.
#' @param x Values that \code{f} is evaluated on.
#' @param y Only used when \code{f} is of class \code{\link[=bifd+-class]{bifd+}}.
#' @param \dots Additional parameters, currently unused.
#' @return See section \strong{Methods} for the actual implementations of
#'         this generic function.
NULL
## Evaluate functional data objects
## *** For basis+ class, it returns an p by T matrix
##     p is the number of basis functions, T is the length of x
## *** For fd+ class, it returns an n by T matrix
##     n is the number of functions, T is the length of x
## *** For bifd+ class, it returns an n1 by n2 matrix
##     n1 is the length of x, n2 is the length of y
if(!isGeneric("feval"))
    setGeneric("feval", function(f, x, ...) standardGeneric("feval"))



## Plot functional data objects
if(!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
if(!isGeneric("plot3d"))
    setGeneric("plot3d", function(x, ...) standardGeneric("plot3d"))


## Inner product of functional data objects
## *** For basis+ class, it returns an p by p matrix
##     p is the number of basis functions
##     P(i, j) = integrate B_i(x) * B_j(x) dx
if(!isGeneric("%*%"))
    setGeneric("%*%")


## Arithmetic operations
## "+", "-", "*", "/" and "^" are already defined in base package as generics


## Penalty matrix of basis function
## P(i, j) = integrate B^(k)_i(x) * B^(k)_j(x) dx
## k is the order of derivative, given by the "penalty" argument
if(!isGeneric("penmat"))
    setGeneric("penmat", function(basis, penalty, ...)
        standardGeneric("penmat"))


## Mean function
if(!isGeneric("mean"))
    setGeneric("mean")


## Standard deviation and variance function
if(!isGeneric("sd"))
    setGeneric("sd")
if(!isGeneric("var"))
    setGeneric("var")

## Covariance function
if(!isGeneric("cov"))
    setGeneric("cov")

## Eigen decomposition of a (symmetric) bivariate function
if(!isGeneric("eig"))
    setGeneric("eig", function(x, k, ...) standardGeneric("eig"))
