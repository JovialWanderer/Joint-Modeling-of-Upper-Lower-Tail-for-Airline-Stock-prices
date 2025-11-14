set.seed(42)
likelihood_trial_meu<-function(Y,meu,X,beta,square_sigma){
  n<-length(Y)
  exp_sum_y<-sum((Y-meu)^2)
  exp_sum_meu<-sum((meu-X%*% beta)^2)
  return((-(exp_sum_y)/2-exp_sum_meu/(2*square_sigma))-((n/2)*log(square_sigma)))
}
#print(likelihood_trial_meu(Y_meu,initial_meu_meu,X_meu,initial_beta_meu,initial_sigma_meu))
prior_beta_trial_meu<-function(beta){
  p<-length(beta)
  return(dmvnorm(beta,mean=rep(0,p),sigma=10000*diag(p)))
}
prior_square_sigma_trial_meu<-function(sig){
  return(dinvgamma(sig,shape=0.01,scale=0.01))
}
posterior_trial_meu<-function(Y,meu,X,beta,square_sigma){
  return(likelihood_trial_meu(Y,meu,X,beta,square_sigma)+
           log(prior_beta_trial_meu(beta))+
           log(prior_square_sigma_trial_meu(square_sigma)))
}
metropolis_hasting_algo_two_layer<-function(Y,ini_meu,X,ini_beta,ini_sigma,iter){
  p<-length(ini_beta)
  n<-length(ini_meu)
  beta_iter<-matrix(0,nrow=iter,ncol=p)
  sigma_iter<-numeric(iter)
  meu_iter<-matrix(0,nrow=iter,ncol=n)
  accept<-0
  
  beta_curr<-ini_beta
  sigma_curr<-ini_sigma
  meu_curr<-ini_meu
  for(i in 1:iter){
    proposed_beta<-beta_curr+rnorm(p,mean = 0,sd = 0.1)
    proposed_sigma<-sigma_curr+rnorm(1,mean = 0,sd = 0.1)
    proposed_meu<-meu_curr+rnorm(n,mean=0,sd=0.1)
    a1<-posterior_trial_meu(Y,proposed_meu,X,proposed_beta,proposed_sigma)
    a2<-posterior_trial_meu(Y,meu_curr,X,beta_curr,sigma_curr)
    #print(a1-a2)
    alpha<-exp(a1-a2)
    #print(alpha)
    a<-runif(1,min=0,max = 1)
    if(a < alpha){
      beta_curr<-proposed_beta
      sigma_curr<-proposed_sigma
      meu_curr<-proposed_meu
      accept<-accept+1
    }
    beta_iter[i,]<-beta_curr
    sigma_iter[i]<-sigma_curr
    meu_iter[i,]<-meu_curr
  }
  
  return(list(beta=beta_iter,squared_sigma=sigma_iter,meu=meu_iter,acceptance_rate=accept/iter))
}

n<-100
p<-3
X_meu<-matrix(rnorm(n*p),ncol=p)
beta_true_meu<-c(2,3,5)
square_sigma_true_meu<-4
meu_true<- X_meu %*% beta_true_meu+rnorm(n,mean=0,sd=sqrt(square_sigma_true_meu))
Y_meu<-rnorm(n,mean=meu_true,sd=1)
initial_beta_meu<-rep(1,p)
initial_sigma_meu<-2
initial_meu_meu<-rep(0,n)
iter<-100000
library(mvtnorm)
ans_meu<-metropolis_hasting_algo_two_layer(Y_meu,initial_meu_meu,X_meu,initial_beta_meu,initial_sigma_meu,iter)
print(ans_meu$meu)
print(ans_meu$beta)
print(ans_meu$squared_sigma)
print(ans_meu$acceptance_rate)
plot(ans_meu$squared_sigma,xlab="Iteration",ylab="squared sigma",type="l")
#try out with new dataset and do the plots