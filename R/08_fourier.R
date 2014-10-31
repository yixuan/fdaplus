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
setMethod("feval", signature(f = "fourier+", x = "numeric"),
          function(f, x, ...) {
              scale = sqrt(2 / f@period)
              omega = 2 * pi / f@period
              ind = as.integer(setdiff(seq(f@ncoef), f@dropind))
              x = x - f@range[1]
              .Call("fourier_feval", x, ind, omega, scale)
          }
)
