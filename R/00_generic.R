#' Converting Objects in the 'fda' Package into their Counterparts in 'fdaplus'
#'
#' This function can convert objects in the \pkg{fda} package
#' into S4 objects that have been redesigned in the \pkg{fdaplus} package.
#'
#' @param obj Object to be converted, whose type should be "\code{basisfd}",
#'            "\code{fd}" or "\code{bifd}".
#' @param \dots Currently unused.
#' @return Below gives the conversion result:
#' \itemize{
#'   \item \code{basisfd} objects will be converted into
#'         \code{\link[=basis+-class]{basis+}} objects
#'   \item \code{fd} objects will be converted into
#'         \code{\link[=fd+-class]{fd+}} objects
#'   \item \code{bifd} objects will be converted into
#'         \code{\link[=bifd+-class]{bifd+}} objects
#' }
#' @author Yixuan Qiu <\url{http://statr.me/}>
#' @export
wrap = function(obj, ...)
    UseMethod("wrap")



#' Subsetting A Functional Data Object
#'
#' This generic function provides methods to subset a functional data
#' object (basis function or univariate function).
#'
#' @name [-methods
#' @rdname subsetter-methods
#'
#' @param x An object of class \code{\link[=basis+-class]{basis+}}
#'          or \code{\link[=fd+-class]{fd+}}.
#' @param i An integer vector giving the indices of functions that are going
#'          to be extracted.
#' @param \dots Additional parameters, currently unused.
#' @return An object of the same class as \code{x} containing the extracted functions.
#' @export
#' @author Yixuan Qiu <\url{http://statr.me/}>
NULL

## "[" is already defined in base package as a generic


## Retrieve a field
## "$" is already defined in base package as a generic


## Concatenate a sequence of functional data objects
## "c" is already a generic S3 method


#' Evaluating A Functional Data Object
#'
#' This generic function provides methods to evaluate different types of
#' functional data (basis function, univariate function, bivariate function)
#' at specified argument values.
#'
#' @name feval-methods
#' @aliases feval
#' @usage feval(f, x, ...)
#'
#' ## For bifd+
#' # feval(f, x, y, ...)
#' @param f An object of class \code{\link[=basis+-class]{basis+}},
#'        \code{\link[=fd+-class]{fd+}}, or \code{\link[=bifd+-class]{bifd+}}.
#' @param x Values that \code{f} is evaluated on.
#' @param y Only used when \code{f} is of class \code{\link[=bifd+-class]{bifd+}},
#'          in which case \code{f} is evaluated on two arguments, \code{x} and
#'          \code{y}.
#' @param \dots Additional parameters, currently unused.
#' @return See \strong{Method (...)} sections for the actual implementations of
#'         this generic function.
#' @export
#' @author Yixuan Qiu <\url{http://statr.me/}>
NULL

if(!isGeneric("feval"))
    setGeneric("feval", function(f, x, ...) standardGeneric("feval"))



#' Visualizing A Functional Data Object
#'
#' This generic function provides methods to visualize different types of
#' functional data (basis function, univariate function, bivariate function).
#'
#' @name plot-methods
#'
#' @param x An object of class \code{\link[=basis+-class]{basis+}},
#'        \code{\link[=fd+-class]{fd+}}, or \code{\link[=bifd+-class]{bifd+}}.
#' @param y Unused and should be missing. For compatibility with the signature
#'        of the generic function.
#' @param \dots Additional graphical parameters passed to \code{matplot()},
#'        \code{persp()} or \code{rgl::persp3d()}.
#'
#' @export
#' @author Yixuan Qiu <\url{http://statr.me/}>
NULL

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
