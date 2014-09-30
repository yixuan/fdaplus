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
        knots = seq(range[0], range[1], length.out = nbasis - order + 2)
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

setMethod("plot", signature(x = "bspline+", y = "missing"),
          function(x, y, ...) {
              args = list(...)
              if(!"ylab" %in% names(args))
                  args = c(args, ylab = "Bspline basis")
              args = c(list(x = x), args)
              do.call(callNextMethod, args)
          }
)

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
              res = .Call("bspline_inprod", x, y, allknotsx, allknotsy)
              dim(res) = c(x@nbasis, y@nbasis)
              if(!length(x@dropind) & !length(y@dropind))
                  return(res)
              else
                  return(res[setdiff(1:nrow(res), x@dropind),
                             setdiff(1:ncol(res), y@dropind)])
          }
)

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
              res = .Call("bspline_penmat", basis, allknots, penalty)
              if(!length(basis@dropind))
                  return(res)
              else
                  return(res[-basis@dropind, -basis@dropind])
          }
)