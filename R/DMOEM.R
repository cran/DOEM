#' The DMOEM is an overrelaxation algorithm in distributed manner, which is used to solve the parameter estimation of Poisson mixture model.
#'
#' @param y is a data matrix
#' @param M is the number of subsets
#' @param K is the number of Poisson distribution
#' @param seed is the recommended way to specify seeds
#' @param alpha0 is the initial value of the mixing weight under the EM algorithm
#' @param lambda0 is the initial value of the mean under the EM algorithm
#' @param MOEMalpha0 is the initial value of the mixing weight under the monotonically overrelaxed EM algorithm
#' @param MOEMlambda0 is the initial value of the mean under the monotonically overrelaxed EM algorithm
#' @param omega is the overrelaxation factor
#' @param T is the number of iterations
#' @param epsilon is the threshold value
#'
#' @return DMOEMtime,DMOEMalpha,DMOEMlambda
#' @export
#'

#' @examples
#' library(stats)
#' set.seed(637351)
#' K=5 
#' alpha1=c(rep(1/K,K)) 
#' lambda1=c(1,2,3,4,5) 
#' n=300 
#' U=sample(c(1:n),n,replace=FALSE)
#' y= c(rep(0,n)) 
#' for(i in 1:n){
#' if(U[i]<=0.2*n){
#' y[i] = rpois(1,lambda1[1])} 
#' else if(U[i]>0.2*n & U[i]<=0.4*n){
#' y[i] = rpois(1,lambda1[2])} 
#' else if(U[i]>0.4*n & U[i]<=0.6*n){
#' y[i] = rpois(1,lambda1[3])} 
#' else if(U[i]>0.6*n & U[i]<=0.8*n){
#' y[i] = rpois(1,lambda1[4])}
#' else if(U[i]>0.8*n ){
#' y[i] = rpois(1,lambda1[5])} 
#' }
#' M=5
#' seed=637351
#' set.seed(123) 
#' e=sample(c(1:n),K)
#' alpha0= MOEMalpha0=e/sum(e)
#' lambda0= MOEMlambda0=c(1.5,2.5,3.5,4.5,5.5)
#' omega=0.8
#' T=10
#' epsilon=0.005
#' DMOEM(y,M,K,seed,alpha0,lambda0,MOEMalpha0,MOEMlambda0,omega,T,epsilon)

DMOEM=function(y,M,K,seed,alpha0,lambda0,MOEMalpha0,MOEMlambda0,omega,T,epsilon){
n=length(y)
nm=n/M
L1=L2=matrix(rep(0,M*K),nrow=M)
set.seed(seed)
mr=matrix(sample(c(1:n),n,replace=FALSE),nrow = M,ncol=nm,byrow=TRUE)
time1=system.time(for (m in 1:M) {
y1=y[mr[m,]] 
alpha=alpha0
lambda=lambda0
MOEMalpha=MOEMalpha0 
MOEMlambda=MOEMlambda0
for (step in 1:T){
w=matrix(rep(0, K*nm), nrow = nm)
wy=matrix(rep(0, K*nm), nrow = nm)
oldalpha=alpha 
oldlambda=lambda
oldMOEMalpha=MOEMalpha 
oldMOEMlambda=MOEMlambda
for (k in 1:K){
for (i in 1:nm){
w[i,k]=(oldalpha[k]*exp(-oldlambda[k])* 
oldlambda[k]^y1[i])/sum(oldalpha*exp(-oldlambda)* 
oldlambda^ y1[i])
wy[i,k]=w[i,k]*y1[i]
}
}
for (k in 1:K){
alpha[k]=sum(w[,k])/nm
lambda[k]=sum(wy[,k])/sum(w[,k])
}
r= c(rep(0,K))
for (k in 1:K){
r[k]=alpha[k]/oldMOEMalpha[k]
}
omegat= omega*(min(r))/(1+omega-omega*(min(r))) 
for (k in 1:K){
MOEMalpha[k]=(1+omegat)*alpha[k]-omegat*oldMOEMalpha[k]
MOEMlambda[k]=(1+omega)*lambda[k]-omega*oldMOEMlambda[k]  
}
    if(max(abs(MOEMalpha-oldMOEMalpha))<epsilon &
       max(abs(MOEMlambda-oldMOEMlambda))<epsilon) break 
    cat(
   "step",step,"\n",
   "MOEMalpha",MOEMalpha,"\n",
   "MOEMlambda",MOEMlambda,"\n"  
)
}
L1[m,]=MOEMalpha 
L2[m,]=MOEMlambda
}
)
time=time1/M
alphamao=colSums(L1)/M 
lambdamao=colSums(L2)/M
return(list(DMOEMtime=time,DMOEMalpha=alphamao,DMOEMlambda=lambdamao))
} 