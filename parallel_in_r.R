
pasymlp<- function(x, del.l, del.u){
  coeff<-del.l*as.numeric(x <= 0)+del.u*as.numeric(x>0)
  out<-coeff*exp(-abs(x) / coeff) / (del.l + del.u)
  out <- ifelse(x <= 0, out, 1 - out)
  return(out)
}#this is the  CDF of AL

dasymlp <- function(x, del.l, del.u){
  coeff <- del.l * as.numeric(x <= 0) + del.u * as.numeric(x > 0)
  out <- exp(-abs(x) / coeff) / (del.l + del.u)
  return(out)
}#this is the PDF of AL

qasymlp <- function(x, del.l, del.u){
  cutoff <- del.l / (del.l + del.u)
  part1 <- (del.l+del.u)/(del.l*as.numeric(x<= cutoff)+del.u*as.numeric(x> cutoff))
  part2 <- ifelse(x<= cutoff, log(x), log(1 -x))
  part3 <- (del.l*as.numeric(x<= cutoff)-del.u*as.numeric(x>cutoff))
  out <- part3*(log(part1) + part2)
  out
}#this is the quantile dunno what is used for

rasymlp <- function(n, del.l, del.u){
  qasymlp(runif(n), del.l, del.u)
}#this is the simulator

#---------------
# Check the functions corresponding to AL(del.l, del.u)
neg.log <- function(params){
  del.l <- 1 / (1 + exp(-params[1]))
  del.u <- 1 / (1 + exp(-params[2]))
  out <- -sum(log(dasymlp(R, del.l, del.u)))
  out
}# mle calculation
# CDF and quantile function of AL(del.l, del.u) + AL(1 - del.l, 1 - del.u)

k1 <- function(del.l, del.u){
  del.l^3 / ((del.l + del.u) * (2 * del.l - 1) * (1 + del.l - del.u))}

k2 <- function(del.l, del.u){
  (del.l - 1)^3 / ((2 * del.l - 1) * (del.l - del.u - 1) * (2 - del.l - del.u))}

k3 <- function(del.l, del.u){
  del.u^3 / ((del.l + del.u) * (2 * del.u - 1) * (del.l - del.u - 1))}

k4 <- function(del.l, del.u){
  (del.u - 1)^3 / ((2 * del.u - 1) * (1 + del.l - del.u) * (2 - del.l - del.u))}

k5 <- function(del.l, del.u){2 / ((1 + 2 * del.u) * (2 * del.u - 3))}

k6 <- function(del.l, del.u){
  (-12 * del.u + 12 * del.u^2 - 5) / ((1 + 2 * del.u)^2 * (2 * del.u - 3)^2)}

k7 <- function(del.l, del.u){4 * del.u^3 / ((1 + 2 * del.u)^2 * (2 * del.u - 1))}

k8 <- function(del.l, del.u){4 * (del.u - 1)^3 / ((2 * del.u - 3)^2 * (2 * del.u - 1))}

k9 <- function(del.l, del.u){4 * del.l^3 / ((1 + 2 * del.l)^2 * (2 * del.l - 1))}

k10 <- function(del.l, del.u){4 * (del.l - 1)^3 / ((2 * del.l - 3)^2 * (2 * del.l - 1))}

k11 <- function(del.l, del.u){2 / ((1 + 2 * del.l) * (2 * del.l - 3))}

k12 <- function(del.l, del.u){
  (-12 * del.l + 12 * del.l^2 - 5) / ((1 + 2 * del.l)^2 * (2 * del.l - 3)^2)}

pasympl.sum <- function(x, del.l, del.u){
  if((del.l == 0.5) & (del.u == 0.5)){
    val <- ifelse(x <= 0, 0.5*(1 - x)*exp(2*x),1 - 0.5*(1 + x)*exp(-2 * x))
  }
  else if((del.l == 0.5) & (del.u != 0.5)){
    val<- ifelse(x<=0, (k5(del.l,del.u)*x - k6(del.l,del.u))*exp(2*x),
                 1 - k7(del.l,del.u)*exp(-x/del.u)-k8(del.l, del.u)*exp(x/(del.u - 1)))
  }
  else if((del.l != 0.5) & (del.u == 0.5)){
    val<-ifelse(x <= 0, k9(del.l, del.u)*exp(x/del.l)+k10(del.l,del.u)*exp(x/(1-del.l)),
                1+(k11(del.l,del.u)*x + k12(del.l,del.u))*exp(-2*x))
  }
  else{
    val<-ifelse(x<= 0,k1(del.l, del.u)*exp(x/del.l) - k2(del.l,del.u)*exp(x/(1-del.l)),
                1+k3(del.l,del.u)*exp(-x/del.u)-k4(del.l,del.u)*exp(x/(del.u - 1)))
  }
  return(val)
}

dasympl.sum <- function(x, del.l, del.u){
  if((del.l == 0.5) & (del.u == 0.5)){
    val <- (1 + abs(x)) * exp(-2 * abs(x)) - 0.5 * exp(-2 * abs(x))
  }
  else if((del.l == 0.5) & (del.u != 0.5)){
    val <- ifelse(x <= 0, exp(2 * x) * (k5(del.l, del.u) * (1 + 2 * x) - 2 * k6(del.l, del.u)),
                  k7(del.l, del.u) / del.u * exp(-x / del.u) + 
                    k8(del.l, del.u) / (1 - del.u) * exp(-x / (1 - del.u)))
  }
  else if((del.l != 0.5) & (del.u == 0.5)){
    val <- ifelse(x <= 0, k9(del.l, del.u) / del.l * exp(x / del.l) + 
                    k10(del.l, del.u) / (1 - del.l) * exp(x / (1 - del.l)),
                  exp(-2*x)*(k11(del.l,del.u)*(1-2*x)-2*k12(del.l,del.u)))
  }
  else{
    val <- ifelse(x <= 0, k1(del.l,del.u)/del.l*exp(x/del.l) - 
                    k2(del.l,del.u)/(1-del.l)*exp(x/(1 - del.l)),
                  -k3(del.l,del.u)/del.u*exp(-x/del.u) + 
                    k4(del.l,del.u)/(1-del.u)*exp(-x/(1-del.u)))
  }
  return(val)
}
qasympl.sum_0.5_0.5 <- function(p){
  solver <- function(x){
    abs(log(ifelse(p<0.5,p,1-p))-log(0.5)-log(1-x)-2*x)
  }
  solver.logit <- function(xx){
    solver(log(xx) - log(1 - xx))
  }
  out<-optimize(solver.logit,c(0, 0.5),tol= 1e-10)$minimum
  out <- log(out)-log(1-out)
  out <- ifelse(p <0.5,out,-out)
  out
}

qasympl.sum_del.l_del.u <- function(p, del.l, del.u){
  cdf0 <-as.numeric(pasympl.sum(0, del.l, del.u))
  #print(length(p))
  if(p <= cdf0){
    solver <- function(x){
      if(del.l > 0.5){
        val <- abs(log(p) - x / del.l - 
                     log(k1(del.l, del.u) - k2(del.l, del.u) * 
                           exp(x * (2 * del.l - 1) / (del.l * (1 - del.l)))))
      }
      else{
        val <- abs(log(p) - x / (1 - del.l) - 
                     log(k1(del.l, del.u) * exp(x * (1 - 2 * del.l) / (del.l * (1 - del.l))) - 
                           k2(del.l, del.u)))}
      val
    }
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(0, cdf0), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)
  }else{
    solver <- function(x){
      if(del.u > 0.5){
        val <- abs(log(1 - p) + x / del.u - 
                     log(-k3(del.l, del.u) + k4(del.l, del.u) *
                           exp(x * (1 - 2 * del.u) / (del.u * (1 - del.u)))))
      }else{
        val <- abs(log(1 - p) + x / (1 - del.u) - 
                     log(-k3(del.l, del.u) * exp(x * (2 * del.u - 1) / (del.u * (1 - del.u))) + 
                           k4(del.l, del.u)))}
      val}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(cdf0, 1), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)}
  out
}

qasympl.sum_0.5_del.u <- function(p, del.u){
  del.l <- 0.5
  cdf0 <- as.numeric(pasympl.sum(0, del.l, del.u))
  if(p <= cdf0){
    solver <- function(x){
      abs(log(p) - 2 * x - log((k5(del.l, del.u) * x) - k6(del.l, del.u)))}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(0, cdf0), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)
  }else{
    solver <- function(x){
      if(del.u > 0.5){
        val <- abs(log(1 - p) + x / del.u - 
                     log(k7(del.l, del.u) + k8(del.l, del.u) * 
                           exp(x * (1 - 2 * del.u) / (del.u * (1 - del.u)))))
      }else{
        val <- abs(log(1 - p) + x / (1 - del.u) - 
                     log(k7(del.l, del.u) * exp(x * (2 * del.u - 1) / (del.u * (1 - del.u))) + 
                           k8(del.l, del.u)))}
      val}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(cdf0, 1), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)}
  out
}

qasympl.sum_del.l_0.5 <- function(p, del.l){
  del.u <- 0.5
  cdf0 <- as.numeric(pasympl.sum(0, del.l, del.u))
  if(p > cdf0){
    solver <- function(x){
      abs(log(1 - p) + 2 * x - log(-k11(del.l, del.u) * x - k12(del.l, del.u)))}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(cdf0, 1), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)
  }else{
    solver <- function(x){
      if(del.l > 0.5){
        val <- abs(log(p) - x / del.l - log(k9(del.l, del.u) + k10(del.l, del.u) * 
                                              exp(x * (2 * del.l - 1) / (del.l * (1 - del.l)))))
      }else{
        val <- abs(log(p) - x / (1 - del.l) - 
                     log(k9(del.l, del.u) * exp(x * (1 - 2 * del.l) / (del.l * (1 - del.l))) + 
                           k10(del.l, del.u)))}
      val}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(0, cdf0), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)}
  out
}

qasympl.sum <- function(p, del.l, del.u){
  if((del.l == 0.5) & (del.u == 0.5)){
    return(qasympl.sum_0.5_0.5(p))
  }else if((del.l == 0.5) & (del.u != 0.5)){
    return(qasympl.sum_0.5_del.u(p, del.u))
  }else if((del.l != 0.5) & (del.u == 0.5)){
    return(qasympl.sum_del.l_0.5(p, del.l))
  }else{
    return(qasympl.sum_del.l_del.u(p, del.l, del.u))
  }
}

neg.log <- function(params){
  del.l <- 1 / (1 + exp(-params[1]))
  del.u <- 1 / (1 + exp(-params[2]))
  
  out <- -sum(log(dasympl.sum(X, del.l, del.u)))
  out}

# bivariate PDF and CDF of X = (X1, X2)

qnorm.pasymlp <- function(x, del.l, del.u){
  coeff <- del.l * as.numeric(x <= 0) + del.u * as.numeric(x > 0)
  out <- coeff * exp(-abs(x) / coeff) / (del.l + del.u)
  out <- ifelse(x <= 0, qnorm(out), -qnorm(out))
  return(out)
}

library(mvtnorm)

pasympl.gauss <- function(x1, x2, del.l, del.u, rho){
  
  den.integrand <- function(x){
    sapply(x, function(xx){
      r <- log(xx) - log(1 - xx)
      q1 <- qnorm.pasymlp(x1 - r, 1 - del.l, 1 - del.u)
      q2 <- qnorm.pasymlp(x2 - r, 1 - del.l, 1 - del.u)
      
      q <- c(q1, q2)
      out <- pmvnorm(lower = rep(-Inf, 2), upper = q, mean = c(0, 0), 
                     corr = matrix(c(1, rho, rho, 1), nrow = 2))[1]
      out <- out * dasymlp(r, del.l, del.u)
      out <- out / (xx * (1 - xx))
      out})}
  
  out <- integrate(den.integrand, lower = 0, upper = 1)$value
  out}

dasympl.gauss <- function(x1, x2, del.l, del.u, rho){
  
  den.integrand <- function(x){
    r <- log(x) - log(1 - x)
    q1 <- qnorm.pasymlp(x1 - r, 1 - del.l, 1 - del.u)
    q2 <- qnorm.pasymlp(x2 - r, 1 - del.l, 1 - del.u)
    
    q <- cbind(q1, q2)
    out <- dmvnorm(q, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2), log = T)
    out <- out + log(dasymlp(x1 - r, 1 - del.l, 1 - del.u)) +
      log(dasymlp(x2 - r, 1 - del.l, 1 - del.u)) + log(dasymlp(r, del.l, del.u))
    out <- out - dnorm(qnorm.pasymlp(x1 - r, 1 - del.l, 1 - del.u), log = T) -
      dnorm(qnorm.pasymlp(x2 - r, 1 - del.l, 1 - del.u), log = T)
    out <- out - log(x) - log(1 - x)
    exp(out)}
  out <- integrate(den.integrand, lower = 0, upper = 1)$value
  
  out}
# maximum likelihood estimation for joint

transform <- list(
  logit = function(x, lower = 0, upper = 1){
    x <- (x - lower) / (upper - lower)
    return(log(x / (1 - x)))
  },
  inv.logit = function(x, lower = 0, upper = 1){
    p <- exp(x) / (1 + exp(x))
    p <- p * (upper - lower) + lower
    return(p)
  }
)
#--------------------
# copula part
c.u1.u2 <- function(u1, u2, del.l, del.u, rho){
  x.1 <- qasympl.sum(u1, del.l, del.u)
  x.2 <- qasympl.sum(u2, del.l, del.u)
  out <- dasympl.gauss(x.1, x.2, del.l, del.u, rho) / 
    (dasympl.sum(x.1, del.l, del.u) * dasympl.sum(x.2, del.l, del.u))
  out
}
#3. negative likelihood
neg.log <- function(params){
  del.l <- transform$inv.logit(params[1], 0.01, 0.99)
  del.u <- transform$inv.logit(params[2], 0.01, 0.99)
  rho <- transform$inv.logit(params[3], -0.99, 0.99)
  out <- sum(sapply(1:nrow(airline_dum), function(i){
    -log(c.u1.u2(airline_dum[1][i],airline_dum[2][i], del.l, del.u, rho))}))
  out
}
n <- dim(airline_matrix)[1]
print(n)
param_vec<-matrix(nrow=568,ncol=3)
neg.log <- function(params){
  del.l <- transform$inv.logit(params[1], 0.01, 0.99)
  del.u <- transform$inv.logit(params[2], 0.01, 0.99)
  rho <- transform$inv.logit(params[3], -0.99, 0.99)
  
  out <- sum(sapply(1:nrow(airline_dum), function(i){
    -log(c.u1.u2(airline_dum[i, 1], airline_dum[i, 2], del.l, del.u, rho))}))
  out}
for(i in 51:(n/2-51)){
  ############################################################
  init <- c(transform$logit(0.2, 0.01, 0.99), 
            transform$logit(0.2, 0.01, 0.99), 
            transform$logit(0.2, -0.99, 0.99))
  airline_dum<-rbind.data.frame(airline_matrix[(i-50):(i+50),],
                                airline_matrix[(620+i-51):(620+i+49),])
  mles.i <- optim(init, neg.log,hessian = T)$par
  mles.i <- c(transform$inv.logit(mles.i[1], 0.01, 0.99),
              transform$inv.logit(mles.i[2], 0.01, 0.99),
              transform$inv.logit(mles.i[3], -0.99, 0.99))
  print(mles.i)
}
print(param_vec)
