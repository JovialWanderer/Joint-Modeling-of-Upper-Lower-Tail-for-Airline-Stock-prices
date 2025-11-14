set.seed(42)
library(mvtnorm)
#model 3
#calculation of C_phi matrix
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
#calculating likelihood for prior of Y|meu and prior of meu 
likelihood_func<-function(Y,meu,X,beta,phi,square_sigma){
  n<-length(Y)
  c_mat_inv<-C_phi_inv(phi,n)
  c_mat<-C_phi_matrix(phi,n)
  Y_diff<-Y-meu
  meu_diff<-meu-X%*%beta
  log_y_sum<-t(Y_diff)%*%Y_diff
  log_meu_sum<-t(meu_diff)%*%(((1/square_sigma)*c_mat_inv) %*% (meu_diff))
  modified_c_mat_det<-square_sigma*c_mat
  det_c_mat<-det(modified_c_mat_det)#abs(square_sigma*((modified_c_mat_det)^(n-1)))
  log_likelihood<-(-1/2)*(log_y_sum+log_meu_sum+log(abs(det_c_mat)))
  return(log_likelihood)
}
#print(likelihood_func(Y_true,meu_true_val,X,beta_true,phi_true,square_sigma_true))
#calculating prior for beta
prior_beta<-function(beta){
  p<-length(beta)
  return(dmvnorm(beta,mean=rep(0,p),sigma=10000*diag(p)))
}
#prior for square sigma parameter
prior_square_sigma<-function(sig){
  return(dinvgamma(sig,shape=0.01,scale=0.01))
}
#prior for the phi parameter
prior_phi<-function(phi){
  return(dunif(phi,min = 0,max = 1))
}
#finding the posterior for the data
posterior_trial<-function(Y,meu,X,beta,phi,square_sigma){
  return(likelihood_func(Y,meu,X,beta,phi,square_sigma)+#
           log(prior_beta(beta))+#log of prior of beta
           log(prior_square_sigma(square_sigma))+#log of prior for square of sigma
           log(prior_phi(phi))#log of prior for phi
         )
}
#metropolis hasting algorithm for the parameters
metropolis_hasting_algo_three<-function(Y,ini_meu,X,ini_beta,ini_sigma,ini_phi,iter){
  #length of beta and meu vector
  p<-length(ini_beta)
  n<-length(ini_meu)
  #beta_iter,sigma_iter,... to store the values of beta,sigma,... at every iteration
  beta_iter<-matrix(0,nrow=iter,ncol=p)
  sigma_iter<-numeric(iter)
  meu_iter<-matrix(0,nrow=iter,ncol=n)
  phi_iter<-numeric(iter)
  accept<-0#to find the accpetance.
  #initializing beta,sigma,...
  beta_curr<-ini_beta
  sigma_curr<-ini_sigma
  meu_curr<-ini_meu
  phi_curr<-ini_phi
  for(i in 1:iter){
    #generate proposals
    
    proposed_beta<-beta_curr+rnorm(p,mean = 0,sd =0.1)
    proposed_sigma<-sigma_curr+rnorm(1,mean = 0,sd =0.1)
    proposed_meu<-meu_curr+rnorm(n,mean=0,sd=0.1)
    proposed_phi<-phi_curr+rnorm(1,mean=0,sd=0.01)
    
    #log of the proposed posterior is defined as proposed_posterior
    proposed_posterior<-posterior_trial(Y,proposed_meu,X,proposed_beta,proposed_phi,proposed_sigma)
    #log of the current posterior is defined as curr_posterior
    curr_posterior<-posterior_trial(Y,meu_curr,X,beta_curr,phi_curr,sigma_curr)
    #print(c(proposed_posterior,curr_posterior,i))
    #to calculate (proposed posterior)/(current posterior) we subtract 
    #curr_posterior from proposed_posterior and exponentiate it
    alpha<-exp(proposed_posterior-curr_posterior)
    #generate uniform r.v. to check for acceptance 
    a<-runif(1,min=0,max = 1)
    #if true then accept
    if(a < alpha){
      beta_curr<-proposed_beta
      sigma_curr<-proposed_sigma
      meu_curr<-proposed_meu
      phi_curr<-proposed_phi
      accept<-accept+1
    }
    beta_iter[i,]<-beta_curr
    sigma_iter[i]<-sigma_curr
    meu_iter[i,]<-meu_curr
    phi_iter[i]<-phi_curr
  }
  #return the values for each iterations
  return(list(beta=beta_iter,squared_sigma=sigma_iter,meu=meu_iter,phi=phi_iter,acceptance_rate=accept/iter))
}
#initializing dataset
n<-100
p<-3
initial_beta<-rep(1,p)
initial_square_sigma<-2
initial_meu<-rnorm(n,mean = rep(0,n),sd=1)
initial_phi<-0.3
#true values of the parameters
iter<-100000
phi_true<-0.5
beta_true<-c(2,3,5)
square_sigma_true<-4
X<-matrix(rnorm(n*p),ncol=p)
meu_true<- X%*% beta_true
sigma_mat_true<-square_sigma_true*C_phi_matrix(phi_true,n)
meu_true_val<-t(rmvnorm(1,mean=meu_true,sigma=sigma_mat_true))
Y_true<-t(rmvnorm(1,mean=meu_true_val,sigma=diag(n)))
#algorithmic step
model3_ans<-metropolis_hasting_algo_three(Y_true,initial_meu,X,initial_beta,initial_square_sigma,initial_phi,iter)
#Plots for the dataset
plot(model3_ans$phi,xlab="Iteration100000",ylab="phi",type="l")
plot(model3_ans$beta[,3],xlab = "Iteration100000",ylab="beta3",type="l")
abline(h=mean(model3_ans$beta[,3]),col="red")
print(model3_ans$acceptance_rate)
print(mean(model3_ans$phi))
