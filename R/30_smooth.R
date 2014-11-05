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
    nT = nrow(Phi)
    nB = ncol(Phi)
    n = ncol(y)
    
    ## If lambda != 0, we use the augmented least square method
    ## Page 89 of FDA book
    ##
    ## R = L' * L, Choleski decomposition of R
    ## However, chol() does not allow semi-positive definite R,
    ## so we first do eigen decomposition R = G * D * G',
    ## and then set L = D^0.5 * G'
    ##
    ## Phi = rbind(Phi, sqrt(lambda) * L)
    ## y = rbind(y, matrix(0, ncol(Phi), ncol(y)))
    if(lambda != 0)
    {
        R = penmat(basis, penalty)
        e = eigen(R)
        L = tcrossprod(diag(x = sqrt(abs(e$values)), nB, nB), e$vectors)
        Phi = rbind(Phi, sqrt(lambda) * L)
        y = rbind(y, matrix(0, nB, n))
    }
    
    res = .lm.fit(Phi, y)
    fd = fd_new(t(res$coefficients), basis)
    QR = res[c("qr", "rank", "qraux", "pivot")]
    class(QR) = "qr"
    Q = qr.Q(QR)
    lambda_df = sum(Q[1:nT, ]^2)
    ## Page 97 of FDA book
    resid = as.matrix(res$residuals)[1:nT, ]
    sse = colSums(resid^2)
    gcv = nT * sse / (nT - lambda_df)^2
    
    return(list(fd = fd, df = lambda_df, gcv = gcv, SSE = sse))
}

