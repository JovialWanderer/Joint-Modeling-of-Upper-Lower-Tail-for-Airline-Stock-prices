#Final Basyesian model
meu_vec<-c()
install.packages("Matrix")
library(Matrix)
blocks<-list()
block_fin<-list()
for(i in 1:length(final.list)){
  meu_vec<-c(meu_vec,final.list[[i]]$mles.final)
  blocks[[i]]<-c(matrix(final.list[[i]]$mles.cov.final,nrow=3,ncol = 3))
}
for(i in 1:length(final.list)){
  mat<-blocks[[i]]
  block_fin[[i]]<-matrix(mat,nrow=3,ncol=3)
}
n<-length(final.list)
p<-2
X_mat<-matrix(rnorm(n*p),ncol=n)
A_mat<-kronecker(X_mat,diag(3))
permutation_mat<-matrix(rep(0,9*n*n),nrow=3*n,ncol=3*n)
for(i in 1:n){
  for(j in 1:n){
    permutation_mat[(j-1)*n+i,(i-1)*3+j]<-1
  }
}
install.packages("MASS")
install.packages("MCMCpack")
library(MCMCpack)
library(MASS)
C_phi_matrix<-function(phi,n){
  c_mat<-matrix(0,nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      c_mat[i,j]<-phi^(abs(i-j))
    }
  }
  return(c_mat)
}
#calculating inverse of c_phi matrix
C_phi_inv<-function(phi,n){
  c_mat_inv<-matrix(0,nrow=n,ncol=n)
  for (i in 1:n){
    if(i==1 || i==n){
      c_mat_inv[i,i]=1
    }
    else{
      c_mat_inv[i,i]<-1+phi^2
    }
    if(i<n){
      c_mat_inv[i,i+1]<-(-phi)
      c_mat_inv[i+1,i]<-(-phi)
    }
  }
  return(c_mat_inv/(1-phi^2))
}
block_diagonal_matrix<-as.matrix(bdiag(block_fin))
meu_update<-function(meu_data,permutation_mat,sigma_meu,sigma_inv,A_mat,beta_curr){
x_beta_mean<-(A_mat %*% beta_curr)
per_meu_std<-permutation_mat %*% sigma_meu %*% t(permutation_mat)
meu_std<-solve(solve(sigma_inv)+solve(per_meu_std))
per_x_beta_mean<-permutation_mat %*% x_beta_mean
meu_mean<- (meu_std %*% ((sigma_inv %*% meu_data)+
             (solve(per_meu_std) %*% per_x_beta_mean)))
meu_updated<-mvrnorm(n=1,mu=meu_mean,Sigma=meu_std)
return(meu_updated)
}
beta_update<-function(A_mat,sigma_meu,permutation_mat,meu){
  meu_per<-permutation_mat %*% meu
  beta_std<-solve(t(A_mat) %*% solve(sigma_meu) %*% A_mat+ (1/10000)*diag(6))
  beta_mean<-beta_std %*% (A_mat %*% solve(sigma_meu) %*% meu_per)
  beta_updated<-mvrnorm(n=1,mu=beta_mean,Sigma=beta_std)
  return(beta_updated)
}
sigma_update<-function(a,b,A_mat,beta,sigma_phi,n){
  nu<-b+3*n
  M<-A_mat %*% beta
  S<-t(meu-M) %*% solve(sigma_phi) %*% (meu-M)
  V_mat<-a*diag(3)+S
  sigma_updated<-riwish(nu,V_mat)
  return(sigma_updated)
}
posterior_trial<-function(Y,meu,X,beta,phi,square_sigma){
  return(likelihood_func(Y,meu,X,beta,phi,square_sigma)+#
           log(prior_beta(beta))+#log of prior of beta
           log(prior_square_sigma(square_sigma))+#log of prior for square of sigma
           log(prior_phi(phi))#log of prior for phi
  )
}
mcmc_algo<-function(meu_data,permutation_mat,sigma.init,
                    sigma_inv,A_mat,beta.init,a,b,phi.init,meu.init){
  meu<-meu.init
  beta<-beta.init
  sigma<-sigma.init
  phi<-phi.init
  #define chains
  meu.chain<-matrix(0,nrow = iters,ncol = 3*n)
  beta.chain<-matrix(0,nrow=iters,ncol=6)
  sigma.chain<-array(0,dim = c(iters,3,3))
  phi.chain<-numeric(iters)
  
  #start gibbs mcmc
  for(i in 1:iters){
    meu<-meu_updated(Y,X_t,sigmasq,beta,n)
    beta<-beta_updated(meu,X,X_t,sigmasq,identity_mat)
    sigma<-sigma_updated(n,a,b,meu,beta,X_t)
    proposed_phi<-phi+rnorm(1,mean=0,sd=0.01)
    
    #log of the proposed posterior is defined as proposed_posterior
    proposed_posterior<-posterior_trial(Y,meu,X,beta,proposed_phi,sigma)
    #log of the current posterior is defined as curr_posterior
    curr_posterior<-posterior_trial(Y,meu,X,beta,phi,sigma)
    #print(c(proposed_posterior,curr_posterior,i))
    #to calculate (proposed posterior)/(current posterior) we subtract 
    #curr_posterior from proposed_posterior and exponentiate it
    alpha<-exp(proposed_posterior-curr_posterior)
    #generate uniform r.v. to check for acceptance 
    a<-runif(1,min=0,max = 1)
    #if true then accept
    if(a < alpha){
      phi<-proposed_phi
      accept<-accept+1
    }
    #print(meu)
    #print(beta)
    #print(sigmasq)
    meu.chain[i,]<-meu
    beta.chain[i,]<-beta
    sigma.chain[i,]<-sigma
    phi.chain[i]<-phi
  }
  
  #return chains
  out<-list(meu.chain=meu.chain,beta.chain=beta.chain,sigmasq.chain=sigmasq.chain)
  return(out)
}