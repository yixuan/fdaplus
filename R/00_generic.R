## Wrap objects in fda package to the new class
wrap = function(obj, ...)
    UseMethod("wrap")



## Retrieve a subset of functional data objects
## "[" is already defined in base package as a generic


## Evaluate functional data objects
## *** For basis+ class, it returns an p by T matrix
##     p is the number of basis functions, T is the length of x
## *** For fd+ class, it returns an n by T matrix
##     n is the number of functions, T is the length of x
if(!isGeneric("feval"))
    setGeneric("feval", function(f, x, ...) standardGeneric("feval"))


## Plot functional data objects
if(!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))


## Inner product of functional data objects
## *** For basis+ class, it returns an p by p matrix
##     p is the number of basis functions
##     P(i, j) = integrate B_i(x) * B_j(x) dx
if(!isGeneric("%*%"))
    setGeneric("%*%")


## Arithmetic operations
## "*" is already defined in base package as a generic


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