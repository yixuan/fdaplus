#' Creating B-spline Basis Functions
#'
#' This function constructs a \code{\link[=bspline+-class]{bspline+}} object
#' that represents a series of B-spline basis functions.
#'
#' \code{nbasis}, \code{order} and \code{breaks} are related to each other,
#' satisfying the condition \code{nbasis == order + length(breaks)}. Hence when
#' any one of \code{nbasis} and \code{breaks} is properly specified,
#' the other one can be computed accordingly. For details,
#'
#' \itemize{
#'   \item If both \code{nbasis} and \code{breaks} are \code{NULL} (the default),
#'         \code{nbasis} will be set to \code{order} and there is no inner knot.
#'   \item If \code{nbasis} is \code{NULL} and \code{breaks} is given, then
#'         \code{nbasis} will be calculated using the equation above.
#'   \item If \code{nbasis} is given and \code{breaks} unspecified, a series
#'         of equally spaced inner knots will be calculated.
#'   \item If both \code{nbasis} and \code{breaks} are provided, then the program
#'         will check their validity.
#' }
#'
#' @param range A length-two numeric vector to define the interval on which
#'              basis functions can be evaluated. Default is \code{c(0, 1)}.
#' @param nbasis The total number of basis functions, including the ones that need
#'               to be dropped. See section \strong{Details}.
#' @param order The order of the B-spline functions, which is the degree of
#'              splines plus one. Default is 4, standing for cubic splines.
#' @param breaks A vector of break points (inner knots) to define the B-spline basis.
#'               This does not include the end points, which are defined by
#'               the parameter \code{range}. See section \strong{Details}.
#' @param dropind Indices of basis functions that need to be dropped. Default is
#'                \code{NULL}, meaning no basis will be dropped.
#'
#' @return A \code{\link[=bspline+-class]{bspline+}} object representing the
#'         basis functions.
#' @author Yixuan Qiu <\url{http://statr.me/}>
#' @export
basis_bs = function(range = c(0, 1), nbasis = NULL, order = 4,
                    breaks = NULL, dropind = NULL)
{
    range = as.numeric(range)
    order = as.integer(order)
    dropind = unique(as.integer(dropind))

    if(is.null(nbasis) & is.null(breaks))
    {
        nbasis = order
        knots = numeric(0)
    } else if(is.null(nbasis) & (!is.null(breaks))) {
        breaks = as.numeric(breaks)
        knots = breaks[!breaks %in% range]
        nbasis = order + length(knots)
    } else if((!is.null(nbasis)) & is.null(breaks)) {
        nbasis = as.integer(nbasis)
        if(nbasis < order)
            stop("nbasis must be greater than or equal to order")
        knots = seq(range[1], range[2], length.out = nbasis - order + 2)
        knots = knots[-c(1, length(knots))]
    } else {
        nbasis = as.integer(nbasis)
        breaks = as.numeric(breaks)
        knots = breaks[!breaks %in% range]
        if(nbasis != order + length(knots))
            stop("nbasis must be equal to order plus the number of inner knots")
    }

    new("bspline+",
        range = range,
        nbasis = nbasis,
        dropind = dropind,
        ncoef = nbasis - length(dropind),
        degree = order - 1,
        knots = knots)
}



#' @rdname feval-methods
#'
#' @section Method (bspline+, numeric):
#' \tabular{lcl}{
#'   \code{f}  \tab - \tab  A \code{\link[=bspline+-class]{bspline+}} object. \cr
#'   \code{x}  \tab - \tab  A numeric vector.
#' }
#'
#' \code{feval(f, x)} returns a matrix \code{R} of
#' \code{f@@ncoef} rows and \code{length(x)} columns, with \code{R[i, j]}
#' representing the value of the \code{i}-th basis function evaluated on \code{x[j]}.
setMethod("feval", signature(f = "bspline+", x = "numeric"),
          function(f, x, ...) {
              ord = f@degree + 1
              mat = splineDesign(c(rep(f@range[1], ord),
                                   f@knots,
                                   rep(f@range[2], ord)),
                                 x, ord, ...)
              if(!length(f@dropind))
                  return(t(mat))
              else
                  return(t(mat[, -f@dropind]))
          }
)

#' @rdname plot-methods
setMethod("plot", signature(x = "bspline+", y = "missing"),
          function(x, y, ...) {
              args = list(...)
              if(!"ylab" %in% names(args))
                  args = c(args, ylab = "Bspline basis")
              args = c(list(x = x), args)
              do.call(callNextMethod, args)
          }
)

#' @rdname inner_product
setMethod("%*%", signature(x = "bspline+", y = "bspline+"),
          function(x, y) {
              if(!isTRUE(all.equal(x@range, y@range)))
                  stop("range of x and y must be the same")

              ordx = as.integer(x@degree + 1)
              ordy = as.integer(y@degree + 1)
              allknotsx = as.numeric(c(rep(x@range[1], ordx), x@knots,
                                       rep(x@range[2], ordx)))
              allknotsy = as.numeric(c(rep(y@range[1], ordy), y@knots,
                                       rep(y@range[2], ordy)))
              res = .Call("bspline_inprod", x, y, allknotsx, allknotsy,
                          PACKAGE = "fdaplus")
              dim(res) = c(x@nbasis, y@nbasis)
              if(!length(x@dropind) & !length(y@dropind))
                  return(res)
              else
                  return(res[setdiff(1:nrow(res), x@dropind),
                             setdiff(1:ncol(res), y@dropind)])
          }
)

#' @rdname penmat-methods
setMethod("penmat", signature(basis = "bspline+", penalty = "numeric"),
          function(basis, penalty, ...) {
              ord = as.integer(basis@degree + 1)
              penalty = as.integer(penalty)
              if(penalty < 0)
                  stop("'penalty' must be >= 0")
              if(penalty >= ord)
                  return(matrix(0, basis@ncoef, basis@ncoef))
              allknots = as.numeric(c(rep(basis@range[1], ord),
                                      basis@knots,
                                      rep(basis@range[2], ord)))
              res = .Call("bspline_penmat", basis, allknots, penalty,
                          PACKAGE = "fdaplus")
              if(!length(basis@dropind))
                  return(res)
              else
                  return(res[-basis@dropind, -basis@dropind])
          }
)
