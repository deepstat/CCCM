library(MASS)

est_CCCM <- function(
        Y,
        alpha=0.05
){
    n <- dim(Y)[1]
    d <- dim(Y)[2]
    p <- dim(Y)[3]

    M  <- matrix(0, d, p)
    M2 <- matrix(0, d, d)
    M3 <- matrix(0, d*p, d*d)
    M4 <- matrix(0, d*d, d*d)
    F_ <- matrix(0, d*p, d*p)
    
    for(iter in 1:n){
        M  <- M + Y[iter,,]
        M2 <- M2 + Y[iter,,] %*% t(Y[iter,,])
        M3 <- M3 + c(Y[iter,,]) %*% t(c(Y[iter,,] %*% t(Y[iter,,])))
        M4 <- M4 + c(Y[iter,,] %*% t(Y[iter,,])) %*%
            t(c(Y[iter,,] %*% t(Y[iter,,])))
        F_ <- F_ + c(Y[iter,,]) %*% t(c(Y[iter,,]))
    }
    M  <- M/n
    M2 <- M2/n
    M3 <- M3/n
    M4 <- M4/n
    F_ <- F_/n
    
    Z <- M%*%t(M)
    
    denom_1 <- ((d-2)*sum(diag(M2)))+sum(M2)+(2*(sum(diag(Z)) - sum(Z)))
    
    d_FTCM_d_M  <- 2*(M-(matrix(1, d, d)%*%M))/denom_1
    d_FTCM_d_M2 <- (((d-2)*diag(rep(1,d)))+matrix(1, d, d))/denom_1
    d_FTCM_d_M2 <- d_FTCM_d_M2 -
        ((d*diag(rep(1,d))) - matrix(1, d, d))/((d*sum(diag(M2))) - sum(M2))
    d_FTCM_d_M2 <- d_FTCM_d_M2/2
    
    Var_vecY       <- F_ - (c(M) %*% t(c(M)))
    Cov_vecY_vecYY <- M3 - (c(M) %*% t(c(M2)))
    Var_vecYY      <- M4 - (c(M2) %*% t(c(M2)))
    
    eta <- t(c(d_FTCM_d_M)) %*% Var_vecY %*% c(d_FTCM_d_M)
    eta <- eta +
        (2 * (t(c(d_FTCM_d_M)) %*% Cov_vecY_vecYY %*% c(d_FTCM_d_M2)))
    eta <- eta +
        (t(c(d_FTCM_d_M2)) %*% Var_vecYY %*% c(d_FTCM_d_M2))
    
    CCCM <- (sum(M2) - sum(diag(M2)) + sum(diag(Z)) - sum(Z))/
        (((d-1)*sum(diag(M2))) + sum(diag(Z)) - sum(Z))
    
    CI <- c(-qnorm(1 - alpha/2) * sqrt(eta/n),
            qnorm(1 - alpha/2) * sqrt(eta/n))
    CI <- CI + atanh(CCCM)
    CI <- tanh(CI)
    
    return(list(
        CCCM       =CCCM,
        lower_bound=CI[1],
        upper_bound=CI[2]
    ))
}
