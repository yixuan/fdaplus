wrap.basisfd = function(obj, ...)
{
    if(obj$type == "bspline")
        res = new("bspline+",
                  range = obj$rangeval,
                  nbasis = obj$nbasis,
                  dropind = as.numeric(obj$dropind),
                  degree = obj$nbasis - length(obj$params) - 1,
                  knots = as.numeric(obj$params))
}


setMethod("feval", signature(f = "bspline+", x = "numeric"),
          function(f, x, ...) {
              ord = f@degree + 1
              mat = splineDesign(c(rep(f@range[1], ord),
                                   f@knots,
                                   rep(f@range[2], ord)),
                                 x, ord, ...)
              return(mat)
          }
)

setMethod("plot", signature(x = "bspline+", y = "missing"),
          function(x, y, ...) {
              x0 = seq(x@range[1], x@range[2], length.out = 101)
              args = list(...)
              if(!"type" %in% names(args))
                  args = c(args, type = "l")
              if(!"xlab" %in% names(args))
                  args = c(args, xlab = "t")
              if(!"ylab" %in% names(args))
                  args = c(args, ylab = "B-spline basis")
              y = feval(x, x0)
              if(length(x@dropind))
                  y = y[, -x@dropind]
              args = c(list(x = x0, y = y), args)
              do.call(graphics::matplot, args)
          }
)

setMethod("%*%", signature(x = "bspline+", y = "bspline+"),
          function(x, y) {
              if(!isTRUE(all.equal(x@range, y@range)))
                  stop("range of x and y must be the same")
              res = .Call("bspline_inprod", x, y, x@range)
              dim(res) = c(x@nbasis, y@nbasis)
              if(!length(x@dropind) & !length(y@dropind))
                  return(res)
              else
                  return(res[setdiff(1:nrow(res), x@dropind),
                             setdiff(1:ncol(res), y@dropind)])
          }
)
