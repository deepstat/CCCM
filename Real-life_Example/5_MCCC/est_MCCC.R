library(MASS)

est_MCCC <- function(
        Y,
        alpha=0.05
){
    n <- dim(Y)[1]
    d <- dim(Y)[2]
    p <- dim(Y)[3]
    
    if(d != 2){
        stop("The number of observers should be two.")
    }

    V_D   <- matrix(0, p, p)
    V_I   <- matrix(0, p, p)
    Sigma <- matrix(0, 2*p*p, 2*p*p)
    Phi_D <- matrix(0, n, p*p)
    Phi_I <- matrix(0, n, p*p)
    
    for(i in 1:n){
        V_D <- V_D + ((Y[i,1,]-Y[i,2,]) %*% t(Y[i,1,]-Y[i,2,]))
        for(j in 1:n){
            if(i != j){
                V_I <- V_I + ((Y[i,1,]-Y[j,2,]) %*% t(Y[i,1,]-Y[j,2,]))
            }
            if(i < j){
                Phi_D[i,] <- Phi_D[i,] +
                    c((Y[i,1,]-Y[i,2,]) %*% t(Y[i,1,]-Y[i,2,])) +
                    c((Y[j,1,]-Y[j,2,]) %*% t(Y[j,1,]-Y[j,2,]))
                Phi_I[i,] <- Phi_I[i,] +
                    c((Y[i,1,]-Y[j,2,]) %*% t(Y[i,1,]-Y[j,2,])) +
                    c((Y[j,1,]-Y[i,2,]) %*% t(Y[j,1,]-Y[i,2,]))
            }
        }
        Phi_D[i,] <- Phi_D[i,]/(2*(n-i))
        Phi_I[i,] <- Phi_I[i,]/(2*(n-i))
    }
    V_D <- V_D/n
    V_I <- V_I/(n*(n-1))
    
    eigen_V_I <- eigen(V_I)
    V_I_U <- eigen_V_I$vectors
    V_I_D <- diag(eigen_V_I$values)
    
    half_inv_V_I <- V_I_U %*% diag(diag(V_I_D)^-.5) %*% t(V_I_U)
    H <- half_inv_V_I %*% V_D %*% half_inv_V_I
    g <- (sum(H^2))^-.5 * c(H)
    r <- 1 - sqrt(sum(H^2)/p)
    dz_dr <- 1/(1-r^2)
    u <- c(c(V_D), c(V_I))
    
    for(i in 1:(n-1)){
        Phi_i <- c(Phi_D[i,], Phi_I[i,])
        Sigma <- Sigma + ((Phi_i - u) %*% t(Phi_i - u))
    }
    Sigma <- 4 * Sigma / (n-1)
    
    Gamma_1 <- kronecker(half_inv_V_I, half_inv_V_I)
    Gamma_2 <- kronecker(half_inv_V_I, diag(rep(1,p))) +
        kronecker(diag(rep(1,p)), half_inv_V_I)
    Gamma_2 <- solve(Gamma_2) %*% (
        kronecker(V_D %*% half_inv_V_I, diag(rep(1,p))) +
            kronecker(diag(rep(1,p)), V_D %*% half_inv_V_I))
    Gamma_2 <- -kronecker(solve(V_I), solve(V_I)) %*% Gamma_2
    Gamma <- cbind(Gamma_1, Gamma_2)
    Sigma_star <- Gamma %*% Sigma %*% t(Gamma)
    var_r <- c(t(g) %*% Sigma_star %*% g)/p
    CI <- c(-qnorm(1-alpha/2), qnorm(1-alpha/2))
    CI <- CI*dz_dr*sqrt(var_r/n)
    CI <- CI + atanh(r)
    CI <- tanh(CI)
    
    return(list(
        MCCC       =r,
        lower_bound=CI[1],
        upper_bound=CI[2]
    ))
}
