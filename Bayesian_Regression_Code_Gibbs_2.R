set.seed(42)
n<-100
p<-3
X<-matrix(rnorm(n*p),ncol=n)
beta_true<-c(2,3,5)
square_sigma_true<-4
identity_mat<-diag(p)
X_t<-t(X)
meu_t<- X_t %*% beta_true+rnorm(n,mean=0,sd=sqrt(square_sigma_true))
Y<-rep(0,n)
for(i in 1:n){
  Y[i]<-rnorm(1,mean=meu_t[i],sd=1)
}
meu.update<-function(Y,X_t,sigmasq,beta,n){
  #print(1)
  meu_mean<-(sigmasq*Y+ X_t %*% (beta))/(sigmasq+1)
  identity_n<-diag(n)
  meu_sd<-(sigmasq/(sigmasq+1))*identity_n
  meu_out<-t(t(meu_mean)+rmvnorm(1,mean=rep(0,n),sigma=meu_sd))
  meu_out
}
#meu.update(Y,X_t,initial_sigma,initial_beta,n)
beta.update<-function(meu,X,X_t,square_sig,identity_mat){
  #print(2)
  inv_sigma<-((X %*% X_t)/square_sig)+identity_mat
  sigma_mat<-solve(inv_sigma)
  meu_vec<-(sigma_mat %*% (X %*% (meu)))/(square_sig)
  beta_out<-t(t(meu_vec)+rmvnorm(1,mean=rep(0,p),sigma = sigma_mat))
  beta_out
}
#beta.update(initial_meu,X,X_t,initial_sigma,identity_mat)
square_sigma.update<-function(n,a,b,meu,beta,X_t){
  #print(3)
  a_new<-n/2+a
  meu_new<-(meu)-(X_t %*% (beta))
  sum_data<-(t(meu_new)) %*% meu_new
  b_new<-b+(1/2*sum_data)
  square_sigma_inv<-rgamma(1,shape=a_new,rate=b_new)
  square_sigma<-1/square_sigma_inv
  square_sigma
}
#square_sigma.update(n,0.01,0.01,initial_meu,initial_beta,X_t)
gibbs_mcmc_meu<-function(Y,meu.init,beta.init,sigmasq.init,X,X_t,identity_mat,a,b,iters,n,p){
  #chain initiation
  meu<-meu.init
  beta<-beta.init
  sigmasq<-sigmasq.init
  
  #define chains
  meu.chain<-matrix(0,nrow = iters,ncol = n)
  beta.chain<-matrix(0,nrow=iters,ncol=p)
  sigmasq.chain<-numeric(iters)
  
  #start gibbs mcmc
  for(i in 1:iters){
    meu<-meu.update(Y,X_t,sigmasq,beta,n)
    beta<-beta.update(meu,X,X_t,sigmasq,identity_mat)
    sigmasq<-square_sigma.update(n,a,b,meu,beta,X_t)
    #print(meu)
    #print(beta)
    #print(sigmasq)
    meu.chain[i,]<-meu
    beta.chain[i,]<-beta
    sigmasq.chain[i]<-sigmasq
  }
  
  #return chains
  out<-list(meu.chain=meu.chain,beta.chain=beta.chain,sigmasq.chain=sigmasq.chain)
  return(out)
}
initial_beta<-rep(1,p)
initial_sigma<-2
initial_meu<-X_t %*% initial_beta+rnorm(n,mean=0,sd=sqrt(initial_sigma))
meu_mcmc.out<-gibbs_mcmc_meu(Y=Y,
                             meu.init = initial_meu,
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
plot(meu_mcmc.out$meu.chain[,84],xlab="Iteration",ylab=84,type="l")
plot(meu_mcmc.out$beta.chain[,1],xlab="Iteration",ylab="beta_1",type="l")
