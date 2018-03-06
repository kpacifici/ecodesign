
# Estimate detection probability
p_hat <- function(n,Y,a=1,b=1){  
  ok  <- n>0 & Y>0
  n   <- n[ok]
  Y   <- Y[ok]

  p   <- seq(0.001,0.999,.001)
  d   <- dbeta(p,a,b,log=TRUE)
  if(sum(ok)>0){for(j in 1:length(Y)){
     d <- d + dbinom(Y[j],n[j],p,log=TRUE) - 
              log(1-dbinom(0,n[j],p))
  }}
  p <- p[which.max(d)]
return(p)}


# Link function and its derivatives
G   <- function(theta){pnorm(theta)}
Gp  <- function(theta){dnorm(theta)}
Gpp <- function(theta){-theta*dnorm(theta)}

# Compute the Gaussian approximation to the log-likelihood
Laplace <- function(Y,n,det,init){

    M        <- rep(0,length(n))
    V        <- rep(0,length(n))

    these    <- which((n>0) & (Y>0))
    if(length(these)>0){
     ini      <- init[these]
     g        <- G(ini)
     gp       <- Gp(ini)
     gpp      <- Gpp(ini)
     V[these] <- (gp^2-gpp*g)/(g^2)
     M[these] <- V[these]*ini + gp/g
    }

    these    <- which((n>0) & (Y==0))
    if(length(these)>0){
     ini      <- init[these]
     g        <- G(ini)
     gp       <- Gp(ini)
     gpp      <- Gpp(ini)
     q        <- (1-det)^(n[these])-1    
     part     <- q*g+1
     V[these] <- ((q*gp)^2 - q*gpp*part)/(part^2)
     M[these] <- V[these]*ini + q*gp/part
    }

return(cbind(M,V))}

# Compute some matrices used in the posterior approximation
matrices <- function(X,SigInv,gamma,LamInv){
   library(emulator)
   SX       <- SigInv%*%X
   inner    <- solve(t(X)%*%SX+LamInv)
   outer    <- SX%*%inner%*%LamInv%*%gamma
   SXiXS    <- quad.tform(inner,SX)
list(SigInv=SigInv,SXiXS=SXiXS,outer=outer)}


# Compute the posterior of theta
post_theta <- function(n,Sigma,lap,v){
    obs <- which(n>0)
    pre <- which(n==0)
    mn  <- 0*n
    sd  <- 0*n

    H   <- solve(Sigma[obs,obs])
    Q   <- H+diag(lap[obs,2])
    V   <- solve(Q)
    M   <- V%*%lap[obs,1]
    
    mn[obs] <- M
    sd[obs] <- sqrt(diag(V))
 
    S12     <- Sigma[pre,obs]
    SH      <- S12%*%H    
    SQ      <- S12%*%(H%*%V%*%H)    
    mn[pre] <- SH%*%M
    sd[pre] <- sqrt(rowSums(SQ*S12) + v[pre] - rowSums(SH*S12))

return(cbind(mn,sd))}

# Approximate the posterior of z-bar given the posterior of theta
theta2zbar <- function(Y,n,det,mnsd,tau=seq(0.005,0.995,0.005)){
  zbar   <- 1+0*Y
  theta0 <- qnorm(tau)
  for(j in which(Y==0)){
    theta   <- mnsd[j,1] + mnsd[j,2]*theta0
    q       <- (1-det)^n[j]
    g       <- pnorm(theta)
    gt      <- q*g/(q*g+1-g)
    zbar[j] <- mean(gt)
  }
return(zbar)} 

# Get the intial esimates to base the Laplace approximation around
get.init <- function(Y,n,det,fudge=1){
  fudge*ifelse(Y>0,1,-1/((1-det)^n+1))
}

# Compute V(D)
criteria <- function(Y,n,SIGMA,PRIVAR,occ){
 mse <- 0
 for(r in 1:ncol(Y)){
   det_hat <- p_hat(n,Y[,r])
   init    <- get.init(Y[,r],n,det_hat)
   lap     <- Laplace(Y[,r],n,det_hat,init)
   app     <- post_theta(n,SIGMA,lap=lap,PRIVAR)
   zbar    <- theta2zbar(Y[,r],n,det_hat,app)
   mse     <- mse + mean((zbar-occ[,r])^2)/nrow(Y)
 }
return(mse)}

