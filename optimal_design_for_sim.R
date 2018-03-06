rm(list=ls())
library(geoR)

setwd("S:\\Documents\\EcoDesign\\Code_for_ecodesign\\")
source("approx.R")

design <- 1

####################################################
######         OPTIMIZATION SETTINGS           #####
####################################################:
 
 m      <- 20 
 M      <- m^2         # number of total sites
 nocc   <- 5           # number of sampling occasions at each sample site
 p      <- 2           # number of covariates
 reps   <- 1000        # number of data sets for MC integration
reps <- 50 # This is faster as a demo
 rho    <- 3           # parameters  
 nu     <- 0.5 
 sigma  <- 1
 ns     <- 36
 Xtype  <- design

 beta.mn <- c(0,1)
 beta.sd <- c(0.5,0.5)
 det.a	 <- 2
 det.b	 <- 2

####################################################
######          SET-UP SPATIAL DATA            #####
####################################################:
 
 s     <- expand.grid(1:m,1:m)
 r     <- sqrt(rowSums((s-m/2)^2))
 if(Xtype==1){
   X <- exp(-10*(r/m)^2)
   X <- cos(s[,1]) + cos(s[,2])
   X <- ifelse(X>0.01,X/2,0.01)
 }
 if(Xtype==2){
   X <- exp(-10*((r-7)/5)^2)
 }
 X <- qnorm(0.98*X+0.01)

####################################################
######       GENERATE DATA FOR MC APPROX       #####
####################################################:

 d     <- as.matrix(dist(s))
 S     <- matern(d,rho,nu)*sigma^2
 P     <- t(chol(S)) 
 Y     <- matrix(0,M,reps)
 occ   <- matrix(0,M,reps)
 for(r in 1:reps){
   set.seed(r*0820)
   beta    <- rnorm(2,beta.mn,beta.sd)
   det     <- rbeta(1,det.a,det.b)
   theta   <- beta[1]+X*beta[2] + P%*%rnorm(M)
   occ[,r] <- rbinom(M,1,pnorm(theta))
   Y[,r]   <- rbinom(M,nocc,det*occ[,r])
 }


####################################################
######           SYMMETRY/ADJACENCY            #####
####################################################:

 x<-cbind(1,X)

 Q     <- solve(S)
 mat   <- matrices(x,Q,gamma=rep(0,p),LamInv=0.1*diag(p))
 s     <- expand.grid(1:m,1:m)-m/2-.5
 uni   <- unique(abs(s))
 nuni  <- nrow(uni)
 sym   <- list()
 first <- NULL
 for(j in 1:nuni){
    sym[[j]] <- which(abs(s[,1])==uni[j,1] & abs(s[,2])==uni[j,2])
    first    <- c(first,max(sym[[j]]))
 }

 duni   <- as.matrix(dist(uni))
 adj1   <- adj2 <- NULL
 for(i in 1:nuni){for(j in 1:nuni){if(i<j & duni[i,j]==1){
    adj1<-c(adj1,i)
    adj2<-c(adj2,j)
 }}}
 na     <- length(adj1)

####################################################
######           EXCHANGE ALGORITHM            #####
####################################################:

 par(mfrow=c(1,1))
 nt <- 20
 nr <- 10
 n  <- matrix(0,M,nr)
 VZ <- rep(0,nr)


 SIGMA         <- S + 10*x%*%t(x)
 PRIVAR        <- diag(SIGMA) 

 for(r in 1:nr){

   print(paste("Random start",r))

   set.seed(r*0820)
   init  <- sample(1:nuni,round(ns/4),replace=FALSE)
   curn  <- rep(0,M)

   for(j in init){curn[sym[[j]]]<-nocc}

    YY           <- Y
    YY[curn==0,] <- 0
    curVZ        <- criteria(YY,curn,SIGMA,PRIVAR,occ)
    done         <- FALSE

    for(t in 1:nt){if(!done){
     print(paste("Loop number",t))
     prevn <- curn
     for(i in sample(1:na)){
       id1 <- sym[[adj1[i]]]
       id2 <- sym[[adj2[i]]]
       if(mean(curn[id1] == curn[id2])<1){  
         cann             <- curn
         cann[c(id1,id2)] <- cann[c(id2,id1)] 
         YY               <- Y
         YY[cann==0,]     <- 0
         canVZ        <- criteria(YY,cann,SIGMA,PRIVAR,occ)
         if(canVZ<curVZ){
           curVZ <- canVZ
           curn  <- cann
           plot(s,main=round(t+i/na,3),cex=2)
           points(s[curn>0,],pch=19,col=2,cex=2)
           points(s[id1,],pch=2)
           points(s[id2,],pch=2)
         }    
       }
     } 
     done <- mean(prevn==curn)==1
   }} 
   n[,r] <- curn
   VZ[r] <- curVZ
}
