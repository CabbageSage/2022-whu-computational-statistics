set.seed(114514)
library(mixtools)
n <- 100
lambda_1 <- 1
lambda_2 <- 10
p <- 0.2

# ---------------
x <- rexpmix(n, c(p, 1-p), c(lambda_1, lambda_2))

v <- rbinom(n, 1, p)
x_1 <- v * rexp(n, lambda_1) + (1-v) * rexp(n, lambda_2)

x_2 <- p * rexp(n, lambda_1) + (1-p) * rexp(n, lambda_2)
# ----------------------
biexpLL <- function(theta, y) {
  # define parameters
  w <- 1/(1+exp(-theta[1]))
  lambda_1 <- exp(theta[2])
  lambda_2 <- exp(theta[3])
  # likelihood function with dexp
  l <- w * dexp((y), rate = 1/lambda_1) + (1 - w) * dexp((y), rate = 1/lambda_2)
  return(- sum(log(l + 1e-9)))
}
o <- optim(par=c(0.5,0.1,0.2),
           fn=biexpLL,
           y=x)

1-1/(1+exp(-o$par[1])) #p=0.3359249
exp(o$par[2]) #lambda1=0.09561656
exp(o$par[3]) # lambda2=0.8987072

#-------------
est <- function(data, ind){
  x <- data[ind]
  o <- optim(par=c(0.5,0.1,0.2),
             fn=biexpLL,
             y=x)
  return(c(1-1/(1+exp(-o$par[1])), exp(o$par[2]), exp(o$par[3])))
}

x_boot <- boot::boot(x, est, 100)
x_boot
#0.09097494
#0.01529895
#0.27952524
# ----------------------------------


boot::boot.ci(x_boot, type=c("norm","basic"))
# ( 0.1529,  0.4992 )   ( 0.1318,  0.5215 )  

# ----------------------
calCI <- function(n,alpha){
  x <- rexpmix(n, c(p, 1-p), c(lambda_1, lambda_2))
  return((n-1) * var(x) / qchisq(alpha, df = n-1))
}
UCL<-replicate(1000,expr=calCI(n=200, alpha=.05))
mean(0.1529 < UCL & UCL < 0.4992) #0.771
mean(0.1318 < UCL & UCL < 0.5215) #0.808

