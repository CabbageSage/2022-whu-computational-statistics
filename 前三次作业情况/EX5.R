####EX5.13

m <- 10000

g <- function(x) {
  ((x^2)*(exp(-x^2/2))/(sqrt(2*pi)))* (x > 1)
}

u <- runif(m)     #using f1
x1<-qnorm(u*(1-pnorm(1))+pnorm(1))
fx1<-exp(-x1^2/2)/((sqrt(2*pi))*(1-pnorm(1)))
fg1 <- g(x1)/fx1
theta.hat1 <- mean(fg1)
se1 <- sd(fg1)/sqrt(m)
theta.hat1
se1

x2<-sqrt(1-2*log(1-u))  ###f2
fx2<-(x2*exp(-x2^2/2))/exp(-1/2)
fg2<-g(x2)/fx2
theta.hat2<-mean(fg2)
se2<-sd(fg2)/sqrt(m)
theta.hat2
se2


###Ex5.15 comparing for stratified importance sampling and importance sampling
set.seed(100)
M<-1000
k <- 5     #number of strata
m <- M / k  #replicates per stratum
N <- 50     #number of times to repeat the estimation
 
                                                                                                 
g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

fgs<-matrix(0,k,m)  ####the j th element is g_j(x_i)/f_j(x_i), i=1,...,m
estimates <- matrix(0, N, 4)
for (i in 1:N) {
  u <- runif(M)     #f3, inverse transform method
  x <- - log(1 - u * (1 - exp(-1)))
  fg <- g(x) / (exp(-x) / (1 - exp(-1)))
  estimates[i, 1] <- mean(fg)    ##hattheta_MC
  estimates[i, 2] <- var(fg)/M
  for (j in 1:k){
    u2 <- runif(m)
    xj <- - log(exp(-(j-1)/5) - u2 * (exp(-(j-1)/5) - exp(-j/5)))
    fgs[j,] <- g(xj) / (exp(-xj) / (exp(-(j-1)/5) - exp(-j/5)))
  }
  estimates[i, 3] <- sum(apply(fgs,1,mean))  ####\hattheta_S
  estimates[i, 4] <- sum(apply(fgs,1,var)/m)
}

results1<-apply(estimates, 2, mean)
results2<-apply(estimates[,c(1,3)], 2, var)  ###SSE

cat("hatthetaMC=",results1[1],"\n")
cat("hatthetaMC.varese=",results1[2],"\n")
cat("hatthetaMC.varsse=",results2[1],"\n")
cat("hatthetaS=",results1[3],"\n")
cat("hatthetaS.varese=",results1[4],"\n")
cat("hatthetaS.varsse=",results2[2],"\n")

(results1[2]-results1[4])/results1[2] ## variance reduction rate for varese
(results2[1]-results2[2])/results2[1] ## variance reduction rate for varsse
