#' The DE-OEM algorithm replaces E-step with stochastic step in distributed manner, which is used to solve the parameter estimation of Poisson mixture model.
#'
#' @param y is a vector
#' @param M is the number of subsets
#' @param K is the number of Poisson distribution
#' @param seed is the recommended way to specify seeds
#' @param alpha0 is the initial value of the mixing weight
#' @param lambda0 is the initial value of the mean
#' @param a represents the power of the reciprocal of the step size
#' @param b indicates that the M-step is not implemented for the first b data points
#'
#' @return DE_OEMtime,DE_OEMalpha,DE_OEMlambda
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
#' alpha0=e/sum(e)
#' lambda0=c(1.5,2.5,3.5,4.5,5.5)
#' a=1
#' b=5
#' DE_OEM(y,M,K,seed,alpha0,lambda0,a,b)

DE_OEM=function(y,M,K,seed,alpha0,lambda0,a,b){
n=length(y)
nm=n/M
L1=L2=matrix(rep(0,M*K),nrow=M)
set.seed(seed)
mr=matrix(sample(c(1:n),n,replace=FALSE),nrow = M,ncol=nm,byrow=TRUE)
time1=system.time(for (m in 1:M) {
y1=y[mr[m,]]
alpha=alpha0
lambda=lambda0
S1=S2= c(rep(0,K))
for (k in 1:K){
S1[k]= (alpha[k]*exp(-lambda[k])* 
lambda[k]^ y1[1])/sum((alpha*exp(-lambda)* lambda^ y1[1]))
S2[k]=(alpha[k]*exp(-lambda[k])* lambda[k]^ 
y1[1])* y1[1]/sum((alpha*exp(-lambda)* lambda^ y1[1]))
}
for (t in 0:(nm-1)) {
oldS1=S1
oldS2=S2
oldalpha=alpha 
oldlambda=lambda
gamma=1/(t+1)^a
w=c(rep(0,K)) 
for (k in 1:K){
w[k]=(oldalpha[k]*exp(-oldlambda[k])* 
oldlambda[k]^ y1[t+1])/sum((oldalpha
*exp(-oldlambda)* oldlambda^ y1[t+1]))
S1[k]= oldS1[k]+ gamma*(w[k]- oldS1[k])
S2[k]= oldS2[k]+ gamma*(w[k]*y1[t+1]- oldS2[k])
if(t<b){
alpha[k]=alpha[k]
lambda[k]=lambda[k]
}
else {
alpha[k]= S1[k] 
lambda[k]= S2[k]/S1[k]
}
}
cat("alpha",alpha,"\n")
cat("lambda",lambda,"\n")
}
L1[m,]=alpha 
L2[m,]=lambda
}
)
time=time1/M
alphamao=colSums(L1)/M 
lambdamao=colSums(L2)/M 
return(list(DE_OEMtime=time,DE_OEMalpha=alphamao,DE_OEMlambda=lambdamao))
}