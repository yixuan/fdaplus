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
    ## Hat = Phi * (Phi' * Phi + lambda * R)^(-1) * Phi'
    PhiT = feval(basis, x)
    Phi = t(PhiT)
    nT = nrow(Phi)
    nB = ncol(Phi)
    
    if(lambda == 0)
    {
        res = .lm.fit(Phi, y)
        coefs = res$coefficients
        lambda_df = nB
        resid = as.matrix(res$residuals)
    } else {
        R = penmat(basis, penalty)
        PhiTPhi = crossprod(Phi)
        PhiRPhi = solve(PhiTPhi + lambda * R, PhiT)
        coefs = PhiRPhi %*% y
        lambda_df = sum(PhiT * PhiRPhi)
        resid = y - Phi %*% (PhiRPhi %*% y)
    }
    
    fd = fd_new(t(coefs), basis)
    sse = colSums(resid^2)
    gcv = nT * sse / (nT - lambda_df)^2
    
    return(list(fd = fd, df = lambda_df, gcv = gcv, SSE = sse))
}
