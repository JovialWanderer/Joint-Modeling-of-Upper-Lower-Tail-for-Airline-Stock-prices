library(mvtnorm)
set.seed(42)
n<-100
p<-3
X<-matrix(rnorm(n*p),ncol=n)
beta_true<-c(2,3,5)
square_sigma_true<-4
X_t<-t(X)
Y<- X_t %*% beta_true+rnorm(n,mean=0,sd=sqrt(square_sigma_true))
identity_mat<-diag(p)
#---------------------------Gibbs Sampling----------------------------------#
#---------------------------Set initial value------------------------------#
initial_beta<-rep(1,p)
initial_sigma<-2
#---------------------------Updating Beta-----------------------------------#
beta.update<-function(Y,X,X_t,square_sig,identity_mat){
  inv_sigma<-((X %*% X_t)/square_sig)+identity_mat*(100)
  sigma_mat<-solve(inv_sigma)
  meu_vec<-(sigma_mat %*% (X %*% Y))
  beta_out<-t(meu_vec)+rmvnorm(1,mean=rep(0,p),sigma = sigma_mat)
  beta_out
}
#----------------------------Updating Sigma---------------------------------#
square_sigma.update<-function(n,a,b,Y,beta,X_t){
  a_new<-n/2+a
  Y_new<-Y-(X_t %*% t(beta))
  sum_data<-(t(Y_new)) %*% Y_new
  b_new<-b+(1/2*sum_data)
  square_sigma_inv<-rgamma(1,shape=a_new,rate=b_new)
  square_sigma<-1/square_sigma_inv
  square_sigma
}
#--------------------------------Gibbs MCMC-----------------------------------#
gibbs_mcmc<-function(Y,beta.init,sigmasq.init,X,X_t,identity_mat,a,b,iters,n,p){
  #chain initiation
  beta<-beta.init
  sigmasq<-sigmasq.init
  
  #define chains
  beta.chain<-matrix(0,nrow=iters,ncol=p)
  sigmasq.chain<-numeric(iters)
  
  #start gibbs mcmc
  for(i in 1:iters){
    beta<-beta.update(Y,X,X_t,sigmasq,identity_mat)
    sigmasq<-square_sigma.update(n,a,b,Y,beta,X_t)
    beta.chain[i,]<-beta
    sigmasq.chain[i]<-sigmasq
  }
  
  #return chains
  out<-list(beta.chain=beta.chain,sigmasq.chain=sigmasq.chain)
  return(out)
}
mcmc.out<-gibbs_mcmc(Y=Y,
                     beta.init = initial_beta,
                     sigmasq.init = initial_sigma,
                     X=X,
                     X_t = X_t,
                     identity_mat = identity_mat,
                     a=0.01,
                     b=0.01,
                     iters = 30000,
                     n=n,
                     p=p
                       )
print(mcmc.out$beta.chain)
print(mcmc.out$sigmasq.chain)
plot(mcmc.out$beta.chain[5:30000,1],xlab="Iteration",ylab="beta_1",type="l")

plot(mcmc.out$beta.chain[5:30000,2],xlab="Iteration",ylab="beta_2",type="l")

plot(mcmc.out$beta.chain[5:30000,3],xlab="Iteration",ylab="beta_3",type="l")

plot(mcmc.out$sigmasq.chain[5:30000],xlab="Iteration",ylab="sigmasq",type="l")
