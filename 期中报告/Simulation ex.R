##从已知的总体分布中实现抽样(即得到y1...yn)####
gendata=function(n,mu,sigma){
  data1=rnorm(n,mu,sigma)
  return(data1)
}

##进行统计推断，估计总体均值#####
est=function(data){
  n=length(data)
  T1=mean(data)
  sT1=sd(data)/sqrt(n)  ##统计量的标准差的估计
  return(list(est=T1,sd=sT1))
}


##评价统计量####################
main=function(n,mu,sigma,Nsim){
##Nsim:number of simulation 模拟重复的次数
  estimate=rep(0,Nsim)    ##存放估计量T1(有s次)
  se.est=rep(0,Nsim)      ##存放标准差的估计
  cp.est=rep(0,Nsim)      ##考察coverage probability的估计值(渐近正态性)
  for(i in 1:Nsim){
    data1=gendata(n,mu,sigma)   ##生成Y1...Yn
    result=est(data1)
    estimate[i]=result$est
    se.est[i]=result$sd
  }
  
  cp.est[mu>=estimate-1.96*se.est&mu<=estimate+1.96*se.est]=1
  Bias=mean(estimate)-mu
  SSE=sd(estimate)  ##抽样的标准差(sampling standard estimator)
  ESE=mean(se.est)  ##统计量理论上的标准差的估计量，通常与SSE比较接近
  CP=mean(cp.est)
  return(list(Bias=Bias,SSE=SSE,ESE=ESE,CP=CP))  ##将CP值和1-α进行比较(在0.85以上较好)
  
}

set.seed(250)   ##建议调试阶段不要设置种子
main(100,0,1,1000)

