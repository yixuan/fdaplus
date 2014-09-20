wrap.basisfd = function(obj, ...)
{
    if(obj$type == "bspline")
        res = new("bspline+",
                  rangeval = obj$rangeval,
                  nbasis = obj$nbasis,
                  dropind = as.numeric(obj$dropind),
                  knots = as.numeric(obj$params))
}


setMethod("feval", signature(f = "bspline+", x = "numeric"),
          function(f, x, ...) {
              mat = bs(x, knots = f@knots,
                       degree = f@nbasis - length(f@knots) - 1,
                       intercept = TRUE,
                       Boundary.knots = f@rangeval)
              ## Preserve dim and dimnames
              attributes(mat) = attributes(mat)[c("dim", "dimnames")]
              return(mat)
          }
)

setMethod("plot", signature(x = "bspline+", y = "missing"),
          function(x, y, ...) {
              x0 = seq(x@rangeval[1], x@rangeval[2], length.out = 101)
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
              if(!isTRUE(all.equal(x@rangeval, y@rangeval)))
                  stop("range of x and y must be the same")
              res = .Call("bspline_inprod", x, y, x@rangeval)
              dim(res) = c(x@nbasis, y@nbasis)
              return(res)
          }
)
