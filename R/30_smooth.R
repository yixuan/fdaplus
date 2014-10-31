# Page 97 of FDA book
smooth_data = function(x, y, basis, penalty = 2L, lambda = 0)
{
    ## X: T x 1
    ## Y: T x n
    x = as.numeric(x)
    y = as.matrix(y)
    if(length(x) != nrow(y))
        stop("length(x) must be equal to ncol(y)")
    
    if(lambda < 0)  lambda = 0
    
    ## Phi: T x B
    ## beta = (Phi' * Phi + lambda * R)^(-1) * Phi' * y
    Phi = t(feval(basis, x))
    if(lambda == 0)
    {
        coefs = .lm.fit(Phi, y)$coefficients
    } else {
        R = penmat(basis, penalty)
        PhiTPhi = crossprod(Phi)
        PhiTy = crossprod(Phi, y)
        coefs = solve(PhiTPhi + lambda * R, PhiTy)
    }
    
    fd = fd_new(t(coefs), basis)
    return(list(fd = fd))
}