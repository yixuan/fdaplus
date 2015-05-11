#' Creating Fourier Basis Functions
#'
#' This function constructs a \code{\link[=fourier+-class]{fourier+}} object
#' that represents a series of Fourier basis functions.
#'
#' @param range A length-two numeric vector to define the interval on which
#'              basis functions can be evaluated. Default is \code{c(0, 1)}.
#' @param nbasis The total number of basis functions, including the ones that need
#'               to be dropped. \code{nbasis} should be an odd number. If an
#'               even number is specified, it will be increased by 1 automatically.
#'               Default is 3.
#' @param period The length of the cycle that Fourier basis functions repeat
#'               themselves. Default is the width of \code{range}.
#' @param dropind Indices of basis functions that need to be dropped. Default is
#'                \code{NULL}, meaning no basis will be dropped.
#'
#' @return A \code{\link[=fourier+-class]{fourier+}} object representing the
#'         basis functions.
#' @author Yixuan Qiu <\url{http://statr.me/}>
#' @export
basis_fourier = function(range = c(0, 1), nbasis = 3, period = diff(range),
                         dropind = NULL)
{
    range = as.numeric(range)
    nbasis = as.integer(nbasis)
    period = as.numeric(period)
    dropind = unique(as.integer(dropind))

    if(nbasis %% 2 == 0)
    {
        warning("nbasis should be an odd number; will be increased by 1")
        nbasis = nbasis + 1
    }

    new("fourier+",
        range = range,
        nbasis = nbasis,
        dropind = dropind,
        ncoef = nbasis - length(dropind),
        period = period)
}



# Assuming starting point is 0, scale = sqrt(2 / T), omega = 2 * pi / T
# B1(x) = scale / sqrt(2)
# B2(x) = scale * sin(omega * x)
# B3(x) = scale * cos(omega * x)
# B4(x) = scale * sin(omega * 2 * x)
# B5(x) = scale * cos(omega * 2 * x)
# ...
# B_{2k}(x) = scale * sin(omega * k * x)
# B_{2k+1}(x) = scale * cos(omega * k * x)

#' @rdname feval-methods
#'
#' @section Method (fourier+, numeric):
#' \tabular{lcl}{
#'   \code{f}  \tab - \tab  A \code{\link[=fourier+-class]{fourier+}} object. \cr
#'   \code{x}  \tab - \tab  A numeric vector.
#' }
#'
#' \code{feval(f, x)} returns a matrix \code{R} of
#' \code{f@@ncoef} rows and \code{length(x)} columns, with \code{R[i, j]}
#' representing the value of the \code{i}-th basis function evaluated on \code{x[j]}.
setMethod("feval", signature(f = "fourier+", x = "numeric"),
          function(f, x, ...) {
              ind = as.integer(setdiff(seq(f@nbasis), f@dropind))
              x = x - f@range[1]
              .Call("fourier_feval", x, ind, f@period)
          }
)

#' @rdname plot-methods
setMethod("plot", signature(x = "fourier+", y = "missing"),
          function(x, y, ...) {
              args = list(...)
              if(!"ylab" %in% names(args))
                  args = c(args, ylab = "Fourier basis")
              args = c(list(x = x), args)
              do.call(callNextMethod, args)
          }
)

setMethod("%*%", signature(x = "fourier+", y = "fourier+"),
          function(x, y) {
              if(!isTRUE(all.equal(x@range, y@range)))
                  stop("range of x and y must be the same")

              indx = as.integer(setdiff(seq(x@nbasis), x@dropind))
              indy = as.integer(setdiff(seq(y@nbasis), y@dropind))
              .Call("fourier_inprod", x@range - x@range[1],
                    indx, indy, x@period, y@period)
          }
)

setMethod("penmat", signature(basis = "fourier+", penalty = "numeric"),
          function(basis, penalty, ...) {
              penalty = as.integer(penalty)
              if(penalty < 0)
                  stop("'penalty' must be >= 0")
              if(penalty == 0)
                  return(basis %*% basis)

              ind = as.integer(setdiff(seq(basis@nbasis), basis@dropind))
              .Call("fourier_penmat", basis@range - basis@range[1],
                    ind, basis@period, penalty)
          }
)
