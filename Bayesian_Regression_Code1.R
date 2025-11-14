library(MCMCpack)
likelihood_trial<-function(Y,X,beta,square_sigma){
  n<-length(Y)
  exp_sum<-sum((Y-X%*% beta)^2)
  return((-exp_sum/(2*square_sigma))-((n/2)*log(square_sigma)))
}
prior_beta_trial<-function(beta){
  p<-length(beta)
  return(dmvnorm(beta,mean=rep(0,p),sigma=100*diag(p)))
}
prior_square_sigma_trial<-function(sig){
  return(dinvgamma(sig,shape=0.01,scale=0.01))
}
posterior_trial<-function(Y,X,beta,square_sigma){
  return(likelihood_trial(Y,X,beta,square_sigma)+
           log(prior_beta_trial(beta))+
           log(prior_square_sigma_trial(square_sigma)))
}
metropolis_hasting_algo<-function(Y,X,ini_beta,ini_sigma,iter){
  p<-length(ini_beta)
  beta_iter<-matrix(0,nrow=iter,ncol=p)
  sigma_iter<-numeric(iter)
  accept<-0
  
  beta_curr<-ini_beta
  sigma_curr<-ini_sigma
  for(i in 1:iter){
    proposed_beta<-beta_curr+rnorm(p,mean = 0,sd = 0.1)
    proposed_sigma<-sigma_curr+rnorm(1,mean = 0,sd = 0.1)

    alpha<-(exp(posterior_trial(Y,X,proposed_beta,proposed_sigma))/
             exp(posterior_trial(Y,X,beta_curr,sigma_curr)))
    print(alpha)
    a<-runif(1,min=0,max = 1)
    if(a < alpha){
      beta_curr<-proposed_beta
      sigma_curr<-proposed_sigma
      accept<-accept+1
    }
    beta_iter[i,]<-beta_curr
    sigma_iter[i]<-sigma_curr
  }
  
  return(list(beta=beta_iter,squared_sigma=sigma_iter,acceptance_rate=accept/iter))
}

set.seed(42)
n<-100
p<-3
X<-matrix(rnorm(n*p),ncol=p)
beta_true<-c(2,3,5)
square_sigma_true<-4
Y<- X %*% beta_true+rnorm(n,mean=0,sd=sqrt(square_sigma_true))

initial_beta<-rep(1,p)
initial_sigma<-2

iter<-1000
library(mvtnorm)
ans<-metropolis_hasting_algo(Y,X,initial_beta,initial_sigma,iter)

print(ans$beta)
print(ans$squared_sigma)
print(ans$acceptance_rate)
library(ggplot2)
posterior_df <- data.frame(cbind(rep(1:iter), ans$beta))
posterior_df
colnames(posterior_df) <- c("Iteration", paste0("Beta", 1:p))

# Melt the data for ggplot2
install.package(reshape2)
library(reshape2)
posterior_melted <- melt(posterior_df, id.vars = "Iteration")

# Create density plots
density_plots <- ggplot(posterior_melted, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Posterior Density Plots", x = "Parameter Value", y = "Density")

# Plot
print(density_plots)