## -----------------------------------------------------------------------------
dt<-USPersonalExpenditure

## -----------------------------------------------------------------------------
n<-paste(row.names(dt),seq=",")
cat("数据集USPersonalExpenditure收集了美国在1940,1945,1950,1955,1960年",n,"五个方面的消费。")

## -----------------------------------------------------------------------------
y<-as.numeric(colnames(dt))
plot(y,dt[1,],main="USPersonalExpenditure",xlab="年份",ylab="消费(十亿)",ylim = c(0,100),type = "l")
for (i in 2:5) {
  lines(y,dt[i,],col=i)
}
legend("topleft",row.names(dt),col=1:5,lty=1)

## -----------------------------------------------------------------------------
dt

## -----------------------------------------------------------------------------
rrayleigh<-function(n=1000,sigma=1){
  return(sqrt(-2*sigma^2*log(runif(n))))
}

## -----------------------------------------------------------------------------
par(mfrow=c(2,3))
hist(rrayleigh(sigma = 1/4))
hist(rrayleigh(sigma = 1/2))
hist(rrayleigh(sigma = 1))
hist(rrayleigh(sigma = 2))
hist(rrayleigh(sigma = 4))
hist(rrayleigh(sigma = 10))

## -----------------------------------------------------------------------------
rbimodal<-function(n=1000,p){
  return(rnorm(n,mean = sample(c(0,3),n,replace = T,prob = c(p,1-p))))
}
hist(rbimodal(p=0.75))

## -----------------------------------------------------------------------------
hist(rbimodal(p = 0.9))
hist(rbimodal(p = 0.75))
hist(rbimodal(p = 0.6))
hist(rbimodal(p = 0.55))
hist(rbimodal(p = 0.5))
hist(rbimodal(p = 0.45))
hist(rbimodal(p = 0.4))
hist(rbimodal(p = 0.25))
hist(rbimodal(p = 0.1))

## -----------------------------------------------------------------------------
rCompoundPoisson<-function(t,n=1,lambda=1,a=1,b=1){
  Nt<-rpois(n,lambda*t)
  X<-rep(0,n)
  for (i in 1:n) {
    X[i]<-sum(rgamma(Nt[i],a,b))
  }
  return(X)
}

## -----------------------------------------------------------------------------
#设定参数
para<-matrix(0,ncol=3,nrow=16)
colnames(para)=c("lambda","a","b")
para[,1]<-c(rep(0.5,4),rep(1,4),rep(2,4),rep(6,4))
para[,2]<-rep(c(0.5,1,1,3),4)
para[,3]<-rep(c(3,3,1,1),4)
#理论值
t<-10
theo<-matrix(0,ncol=2,nrow=16)
colnames(theo)=c("theo-mean","theo-var")
EY<-para[,2]/para[,3]
EY2<-EY^2+para[,2]/(para[,3]^2)
theo[,1]<-para[,1]*t*EY
theo[,2]<-para[,1]*t*EY2
#模拟值
est<-matrix(0,ncol=2,nrow=16)
colnames(est)=c("est-mean","est-var")
for (k in 1:16) {
  x<-rCompoundPoisson(t=10,n=1000,para[k,1],para[k,2],para[k,3])
  est[k,1]<-mean(x)
  est[k,2]<-var(x)
}
#输出结果
res<-cbind(para,est,theo)
#print(res)
knitr::kable(res)

## -----------------------------------------------------------------------------
x<-seq(0,1,0.1)
m<-1000
u<-runif(m)
cdf<-numeric(length(x))
for (i in 1:length(x)) {
  g<-x[i]^3*u^2*(1-x[i]*u)^2
  cdf[i]<-mean(g)/beta(3,3)
}
Phi<-pbeta(x,3,3)
print(round(rbind(x,cdf,Phi),3))

## -----------------------------------------------------------------------------
rRayleigh<-function(n=1000,sigma=1,anti= TRUE){
  u<-runif(n)
  if (!anti) v<-runif(n) else
    v <- 1-u
  x1<-sqrt(-2*sigma^2*log(u))
  x2<-sqrt(-2*sigma^2*log(v))
  x<-(x1+x2)/2
  return(x)
}

## -----------------------------------------------------------------------------
K<-100
# 对偶变量的方差
sd_an<-0
for (k in 1:K) {
  sd_an<-sd_an+var(rRayleigh(anti= TRUE))
}
sd_an<-sd_an/K
print(sd_an)
# 独立变量的方差
sd_in<-0
for (k in 1:K) {
  sd_in<-sd_in+var(rRayleigh(anti= FALSE))
}
sd_in<-sd_in/K
print(sd_in)
# 对比方差减小比例
print((sd_in-sd_an)/sd_in)

## ----fig.width=10-------------------------------------------------------------
    x <- seq(0, 1, .01)
    w <- 2
    h <- 1/x^4/sqrt(2*pi)*exp(-1/2/x^2) * (x > 0) * (x < 1)
    f1 <- rep(1,length(x))
    f2 <- (1 / pi) / (1 + x^2)
    gs <- c(expression(h(x)==1/(x^4*sqrt(2*pi))*exp(-1/(2*x^2)) *1[(x>0)]* 1[(x < 1)]),
            expression(f[1](x)==1),
            expression(f[2](x)==(1+x^2)^{-1}/pi))
    #for color change lty to col
    par(mfrow=c(1,2))
    #figure (a)
    plot(x, h, type = "l", ylab = "",
         ylim = c(0,2), lwd = w,col=1,main='(A)')
    lines(x, f1, lty = 2, lwd = w,col=2)
    lines(x, f2, lty = 3, lwd = w,col=3)
    legend("topright", legend = gs,
           lty = 1:3, lwd = w, inset = 0.02,col=1:3)

    #figure (b)
    plot(x, h/f1, type = "l", ylab = "",
        ylim = c(0,4.5), lwd = w, lty = 2,col=2,main='(B)')
    lines(x, h/f2, lty = 3, lwd = w,col=3)
    legend("topright", legend = gs[-1],
           lty = 2:3, lwd = w, inset = 0.02,col=2:3)

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- se <- numeric(2)
h <- function(x) {
  1/x^4/sqrt(2*pi)*exp(-1/2/x^2) * (x > 0) * (x < 1)
}

## -----------------------------------------------------------------------------
# f1
x <- runif(m) 
fh <- h(x)
theta.hat[1] <- mean(fh)
se[1] <- sd(fh)
# f2
x <- rexp(m, 1)
fh <- h(x) / exp(-x)
theta.hat[2] <- mean(fh)
se[2] <- sd(fh)
# 结果
rbind(theta.hat,se)

## -----------------------------------------------------------------------------
n<-20
alpha<-0.05
CL <- replicate(1000, expr = {
x <- rchisq(n,2)
c(mean(x)+ qt(alpha/2, df=n-1)*sqrt(var(x))/ sqrt(n),mean(x)- qt(alpha/2, df=n-1)*sqrt(var(x))/ sqrt(n))
} )
#compute the mean to get the confidence level
length(intersect(which(CL[1,]<2),which(CL[2,]>2)))/1000

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values

## -----------------------------------------------------------------------------
for (j in 1:m) {
  x <- rchisq(n, 1)
  ttest <- t.test(x, alternative = "greater", mu = mu0)
  p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
for (j in 1:m) {
  x <- runif(n,0,2)
  ttest <- t.test(x, alternative = "greater", mu = mu0)
  p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
for (j in 1:m) {
  x <- rexp(n, 1)
  ttest <- t.test(x, alternative = "greater", mu = mu0)
  p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
nn <- c(10,20,30,50,100,500)  # 样本容量
alpha <- 0.05                 # 显著性水平
d <- 2                        # 随机变量的维数
b0 <- qchisq(1-alpha,df=d*(d+1)*(d+2)/6)*6/nn  # 每种样本容量临界值向量

# 计算多元样本偏度统计量
mul.sk <- function(x){
  n <- nrow(x) # 样本个数
  xbar <- colMeans(x) 
  sigma.hat <- (n-1)/n*cov(x) # MLE估计
  
  b <- 0
  for(i in 1:nrow(x)){
    for(j in 1:nrow(x)){
      b <- b+((x[i,]-xbar)%*%solve(sigma.hat)%*%(x[j,]-xbar))^3
    }
  }
  return(b/(n^2))
}

# 计算第一类错误的经验估计
library(mvtnorm)
set.seed(200)
p.reject <- vector(mode = "numeric",length = length(nn)) # 保存模拟结果

m <- 100

for(i in 1:length(nn)){
  mul.sktests <- vector(mode = "numeric",length = m)
  for(j in 1:m){
    data <- rmvnorm(nn[i],mean = rep(0,d))
    mul.sktests[j] <- as.integer(mul.sk(data)>b0[i])
  }
  p.reject[i] <- mean(mul.sktests)
}
p.reject

## -----------------------------------------------------------------------------
summ <- rbind(nn,p.reject)
rownames(summ) <- c("n","estimate")
knitr::kable(summ)

## -----------------------------------------------------------------------------
alpha <- 0.1
n <- 30      # 样本大小
m <- 200   # 重复次数
epsilon <- c(seq(0,0.15,0.01),seq(0.15,1,0.05))
N <- length(epsilon)
power <- vector(mode = "numeric",length = N)
b0 <- qchisq(1-alpha,df=d*(d+1)*(d+2)/6)*6/n  #临界值

# 对这列epsilon分别求power
for(j in 1:N){
  e <- epsilon[j]
  mul.sktests <- numeric(m)
  for(i in 1:m){
    # 生成混合分布
    u <- sample(c(1,0),size = n,replace = T,prob = c(1-e,e))
    data1 <- rmvnorm(n,sigma = diag(1,d))
    data2 <- rmvnorm(n,sigma = diag(100,d))
    data <- u*data1+(1-u)*data2
    mul.sktests[i] <- as.integer(mul.sk(data)>b0)
  }
  power[j] <- mean(mul.sktests)
}

# 绘制功效函数
plot(epsilon,power,type="b",xlab=bquote(epsilon),ylim=c(0,1))
abline(h=0.1,lty=3,col="lightblue")
se <- sqrt(power*(1-power)/m)  # 绘制标准误差
lines(epsilon,power-se,lty=3)
lines(epsilon,power+se,lty=3)

## -----------------------------------------------------------------------------
library(bootstrap)
data("scor",package = "bootstrap")
n<-nrow(scor)
lam.hat<-eigen(cov(scor))$value
theta.hat<-lam.hat[1]/sum(lam.hat)
theta.hat

## -----------------------------------------------------------------------------
B<-50
theta.b<-numeric(B)
for (b in 1:B) {
  i<-sample(1:n,size = n,replace = TRUE)
  S<-cov(scor[i,])
  lam<-eigen(S)$value
  theta.b[i]<-lam[1]/sum(lam)
}
bias<-mean(theta.b-theta.hat)
bias
sd<-sd(theta.b)
sd

## -----------------------------------------------------------------------------
theta.j<-numeric(n)
for (i in 1:n) {
  S<-cov(scor[-i,])
  lam<-eigen(S)$value
  theta.j[i]<-lam[1]/sum(lam)
}
bias<-(n-1)*(mean(theta.j)-theta.hat)
bias
se<-sqrt((n-1)*mean((theta.j-mean(theta.j))^2))
se

## -----------------------------------------------------------------------------
library(boot)
theta.boot<-function(dat,ind){
  S<-cov(dat[ind,])
  lam<-eigen(S)$value
  lam[1]/sum(lam)
}
boot.obj<-boot(scor,statistic = theta.boot,R=2000)
print(boot.ci(boot.obj,type = "perc"))

## -----------------------------------------------------------------------------
boot.BCa<-function(x,th0,th,stat,conf=.95){
  x<-as.matrix(x)
  n<-nrow(x)
  N<-1:n
  alpha<-(1+c(-conf,conf))/2
  zalpha<-qnorm(alpha)
  z0<-qnorm(sum(th<th0)/length(th))
  th.jack <- numeric(n)
  for (i in 1:n) {
    J <- N[1:(n-1)]
    th.jack[i] <- stat(x[-i, ], J)
  }
  L <- mean(th.jack) - th.jack
  a <- sum(L^3)/(6 * sum(L^2)^1.5)
  adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
  limits <- quantile(th, adj.alpha, type=6)
  return(list("est"=th0, "BCa"=limits))
}
boot.BCa(scor,th0<-theta.hat,th=theta.b,stat=theta.boot)

## -----------------------------------------------------------------------------
library(boot)
n<-10000
x<-rnorm(n)
theta.boot<-function(dat,ind){
  mean(dat[ind])
}
boot.obj<-boot(x,statistic = theta.boot,R=2000)
ci<-boot.ci(boot.obj,type = c("norm","basic","perc"))
print(ci)

## -----------------------------------------------------------------------------
intervals<-matrix(0,nrow = 3,ncol = 2)
intervals[1,]<-ci$normal[-1]
intervals[2,]<-ci$basic[-3:-1]
intervals[3,]<-ci$percent[-3:-1]

## -----------------------------------------------------------------------------
M<-1000
mont.mean<-numeric(M)
for (i in 1:M) {
  xm<-rnorm(n)
  mont.mean[i]<-mean(xm)
}
cover<-rep(0,3)
left<-rep(0,3)
right<-rep(0,3)

## -----------------------------------------------------------------------------
for (j in 1:3) {
  left.id<-which(mont.mean<intervals[j,1])
  right.id<-which(mont.mean>intervals[j,2])
  left[j]<-length(left.id)
  right[j]<-length(right.id)
  cover[j]<-M-left[j]-right[j]
}
cover/M # empirical coverage rates for the sample mean
left/M #proportion of times that the confidence intervals miss on the left
right/M #porportion of times that the confidence intervals miss on the right

## -----------------------------------------------------------------------------
library(boot)
n<-10000
x<-rchisq(n,5)
theta.boot<-function(dat,ind){
  mean(dat[ind])
}
boot.obj<-boot(x,statistic = theta.boot,R=2000)
ci<-boot.ci(boot.obj,type = c("norm","basic","perc"))
print(ci)

## -----------------------------------------------------------------------------
intervals<-matrix(0,nrow = 3,ncol = 2)
intervals[1,]<-ci$normal[-1]
intervals[2,]<-ci$basic[-3:-1]
intervals[3,]<-ci$percent[-3:-1]

## -----------------------------------------------------------------------------
M<-1000
mont.mean<-numeric(M)
for (i in 1:M) {
  xm<-rchisq(n,5)
  mont.mean[i]<-mean(xm)
}
cover<-rep(0,3)
left<-rep(0,3)
right<-rep(0,3)

## -----------------------------------------------------------------------------
for (j in 1:3) {
  left.id<-which(mont.mean<intervals[j,1])
  right.id<-which(mont.mean>intervals[j,2])
  left[j]<-length(left.id)
  right[j]<-length(right.id)
  cover[j]<-M-left[j]-right[j]
}
cover/M # empirical coverage rates for the sample mean
left/M #proportion of times that the confidence intervals miss on the left
right/M #porportion of times that the confidence intervals miss on the right

## -----------------------------------------------------------------------------
n<-20
x<-sample(1:50,n)
y<-sample(1:50,n)

## -----------------------------------------------------------------------------
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- n*2
set.seed(0)
reps <- numeric(R) #storage for replicates
t0 <- cor(x,y,method = "spearman")
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- cor(x1, y1,method = "spearman")
}
p <- mean(abs(c(t0, reps)) >= abs(t0))

round(c(p,cor.test(x,y,method = "spearman")$p.value),3)


## -----------------------------------------------------------------------------
n<-20
x<-sample(1:50,n)
y<-x*2

## -----------------------------------------------------------------------------
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- n*2
set.seed(0)
reps <- numeric(R) #storage for replicates
t0 <- cor(x,y,method = "spearman")
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- cor(x1, y1,method = "spearman")
}
p <- mean(abs(c(t0, reps)) >= abs(t0))

round(c(p,cor.test(x,y,method = "spearman")$p.value),3)

## ----warning=F----------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]
  n2 <- sizes[2]
  n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0)
  z <- z[ix, ]
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1] 
  block2 <- NN$nn.idx[(n1+1):n,-1] 
  i1 <- sum(block1 < n1 +.5)
  i2 <- sum(block2 > n1+.5) 
  (i1+i2) / (k*n)
}
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
  sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}


## -----------------------------------------------------------------------------
m <- 1000
k<-3
R<-99
set.seed(12345)
n1 <- n2 <- 100
n <- n1+n2
N = c(n1,n2)

## -----------------------------------------------------------------------------
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- rnorm(n1,0,1)
  y <- rnorm(n2,0,1.4)
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12345)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- rnorm(n1,0,1)
  y <- rnorm(n2,1/3,1)
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12345)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- rt(n1,1)
  y <- rt(n2,2)
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12345)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- rnorm(n1,-1,1)+rnorm(n1,1,1)
  y <- rnorm(n2,-1/2,1)+rnorm(n2,1,1)
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12345)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1 <- 10
n2 <- 20
n <- n1+n2
N = c(n1,n2)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- rnorm(n1,0,1)
  y <- rnorm(n2,0,3)
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12345)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
f <- function(x) {
  return(1/pi/(1+x^2))
}
m <- 10000
f1<-function(m){
  x <- numeric(m)
  x[1] <- rnorm(1)
  k <- 0
  u <- runif(m)
  for (i in 2:m) {
    xt <- x[i-1]
    y <- rnorm(1)
    num <- f(y) * dnorm(xt)
    den <- f(xt) * dnorm(y)
    if (u[i] <= num/den){
      x[i] <- y
      } else {
        x[i] <- xt
        k <- k+1     #y is rejected
      }
  }
  return(x)
}
set.seed(0)
x<-f1(m)

## -----------------------------------------------------------------------------
index <- 5000:5500
y1 <- x[index]
plot(index, y1, type="l", main="", ylab="x")

## -----------------------------------------------------------------------------
b <- 1001      #discard the burn-in sample
y <- x[b:m]
a <- ppoints(10)

QC <- qcauchy(a)   #quantiles of Cauchy
Q <- quantile(y, a)
par(mfrow=c(1,2))
    qqplot(QC, Q, main="",
        xlab="Cauchy Quantiles", ylab="Sample Quantiles")
    abline(c(0,0),c(1,1),col='blue',lwd=2)
    hist(y, breaks="scott", main="", xlab="", freq=FALSE,xlim = c(0,5))
    lines(QC, f(QC))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
#generate the chains
m<-15000
k<-4
X <- matrix(0, nrow=k, ncol=m)
for (i in 1:k) {
  set.seed(i*12345)
  X[i, ] <- f1(m)
}
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i,] <- psi[i,] / (1:ncol(psi))
for (i in 1:k){
   if(i==1){
      plot(1:(m-b+1),psi[i,b:m],ylim=c(-0.2,0.3), type="l", xlab='Index', ylab=bquote(phi))
    }else{
      lines(psi[i, b:m], col=i)
    }
}
par(mfrow=c(1,1)) #restore default

## -----------------------------------------------------------------------------
rhat <- rep(0,(m-b+1))
for (j in 1:(m-b+1))  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat, type="l", xlab="", ylab="R",ylim = c(1,3))
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
#initialize constants and parameters
N <- 5000 #length of chain
burn <- 1000 #burn-in length
f2<-function(m,a=2,b=2,n=10){
  X <- matrix(0, m, 2) #the chain, a bivariate sample
  ###### generate the chain #####
  X[1, ] <- c(1, .5) #initialize
  for (i in 2:m) {
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, y)
    x <- X[i, 1]
    X[i, 2] <- rbeta(1,x+a,n-x+b)
    }
  b <- burn + 1
  return(X[b:m, ])
}
x<-f2(N)
cat('Means: ',round(colMeans(x),2))
cat('Standard errors: ',round(apply(x,2,sd),2))
cat('Correlation coefficients: ', round(cor(x[,1],x[,2]),2))
plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
lines(x[,2],col=2,lwd=2)
legend('bottomright',c(expression(x),expression(y)),col=1:2,lwd=2)
plot(x,xlab='x',ylab='y')
hist(x[,1],xlab='Random numbers[x]')
hist(x[,2],xlab='Random numbers[y]')

## -----------------------------------------------------------------------------
m<-20000
k<-4
X <- matrix(0, nrow=k, ncol=m-b+1)
Y <- matrix(0, nrow=k, ncol=m-b+1)
for (i in 1:k) {
  set.seed(i*123)
  ff<-f2(m)
  X[i, ] <- t(ff[,1])
  Y[i, ] <- t(ff[,2])
}
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i,] <- psi[i,] / (1:ncol(psi))
for (i in 1:k){
   if(i==1){
      plot(1:(m-b+1),psi[i,], type="l", xlab='Index', ylab=bquote(phi))
    }else{
      lines(psi[i, ], col=i)
    }
}
par(mfrow=c(1,1)) #restore default

## -----------------------------------------------------------------------------
rhat1 <- rep(0,m-b+1)
for (j in 1:(m-b+1))  rhat1[j] <- Gelman.Rubin(psi[,1:j])

## -----------------------------------------------------------------------------
X<-Y
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i,] <- psi[i,] / (1:ncol(psi))
for (i in 1:k){
   if(i==1){
      plot(1:(m-b+1),psi[i,], type="l", xlab='Index', ylab=bquote(phi))
    }else{
      lines(psi[i, ], col=i)
    }
}
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
rhat2 <- rep(0,m-b+1)
for (j in 1:(m-b+1))  rhat2[j] <- Gelman.Rubin(psi[,1:j])

## -----------------------------------------------------------------------------
plot(rhat1, type="l", xlab="", ylab="R",ylim = c(1,3))
lines(rhat2,col=2)
abline(h=1.2, lty=2)
legend('topright',c(expression(x),expression(y)),col=1:2,lwd=2)

## -----------------------------------------------------------------------------
fk<-function(a,k){
  d<-length(a)
  rk<-(-1)^k/factorial(k)/(2^k)*(sum(a^2))^(k+1)/(2*k+1)/(2*k+2)*gamma((d+1)/2)*gamma(k+3/2)/gamma(k+d/2+1)
  return(rk)
}
a<-rep(1,5)
fk(a,1)
fk(a,2)
fk(a,3)
fk(a,4)
fk(a,5)

## -----------------------------------------------------------------------------
f<-function(a){
  r<-0
  for (i in 1:20) {
    if(abs(fk(a,i))>10^(-5))r<-r+fk(a,i) else break
  }
  return(round(r,5))
}
f(rep(1,5))
f(c(1,2,3))
f(rnorm(10))

## -----------------------------------------------------------------------------
a<-c(1,2)
f(a)

## -----------------------------------------------------------------------------
g0<-function(u,k=4){
  (1+u^2/(k-1))^(-k/2)
}

gleft<-function(a,k=4){
  2*gamma(k/2)/sqrt(pi*(k-1))/gamma((k-1)/2)*integrate(g0,lower = 0,upper = sqrt(a^2*(k-1)/(k-a^2)))$value
}


## ----warning=F----------------------------------------------------------------
k0<-c(4:25,100)
res<-vector(length = 0)
for (k in k0) {
  g<-function(a){
   gleft(a,k)-gleft(a,k+1)
  }
  solution<-uniroot(g,lower = 0.1,upper = 2)
  res[as.character(k)]<-round(unlist(solution),5)

}
res


## ----warning=F----------------------------------------------------------------
S<-function(a,k){
  1-pt(sqrt(a^2*k/(k+1-a^2)),k)
}
k0<-c(4:25,100,500,1000)
res0<-vector(length = 0)
for (k in k0) {
  SS<-function(a){
  S(a,k)-S(a,k-1)
  }
  solution<-uniroot(SS,lower = 0.1,upper = 2)
  res0[as.character(k)]<-round(unlist(solution),5)

}
rbind(res0,c(res,0,0))

## -----------------------------------------------------------------------------
N<-20
l<-rep(2,N)
x0<-c(0.54, 0.48, 0.33, 0.43,  0.91,  0.21, 0.85)
n1<-7
n2<-3
n<-10
eme<-function(la){
  return(1+1/la)
}
emm<-function(xm){
  return(n/(sum(x0)+sum(xm)))
}
for (k in 1:N) {
  xm<-eme(rep(l[k],n2))
  l[k+1]<-emm(xm)
  if(abs(l[k+1]-l[k])<10^(-5))break
}
l[1:(k+1)]#收敛过程
l[k+1]#EM方法结果

## -----------------------------------------------------------------------------
n1/(sum(x0)+n2)

## -----------------------------------------------------------------------------
trims<-c(0,0.1,0.2,0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
mod1<-lapply(formulas,lm,data= mtcars)
lapply(mod1, rsq)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

## -----------------------------------------------------------------------------
mod2<-lapply(bootstraps,lm,formula = mpg ~ disp)
lapply(mod2, rsq)

## -----------------------------------------------------------------------------
dt<-matrix(sample(1:100),nrow = 20)
df<-data.frame(dt)#一个20*5的numeric data frame
df
vapply(df,sd,FUN.VALUE=0)

## -----------------------------------------------------------------------------
dff<-cbind(df,rep("a",20),df+4,rep("b",20))
colnames(dff)<-as.character(1:12) #一个20*12的mixed data frame
dff
id<-vapply(dff, is.numeric, FUN.VALUE = 1)#判断
vapply(dff[id],sd,FUN.VALUE=0)

## ----warning=FALSE------------------------------------------------------------
library(parallel)
mcsapply<-function(X, FUN, ...,mc.cores = 1L, simplify=TRUE, USE.NAMES = TRUE){
  makeCluster(getOption("cl.cores",mc.cores))->cl
  mcl<-parLapply(cl,X,FUN,mc.cores = mc.cores)
  if(simplify==TRUE) mcl<-unlist(mcl)
  if(USE.NAMES==FALSE) row.names(mcl)<-NULL
  return(mcl)
}
#参数取值对比
mcsapply(df,mean)
mcsapply(df,mean,simplify=FALSE)
mcsapply(df,mean,USE.NAMES = FALSE)
#时间对比
system.time(mcsapply(df,mean))
system.time(mcsapply(df,mean,mc.cores = 2))
system.time(mcsapply(df,mean,mc.cores = 4))
system.time(mcsapply(df,mean,mc.cores = 8))

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('
NumericMatrix gibbsc(int N,int burn,int a,int b,int n){
  NumericMatrix X(N,2);
  NumericMatrix XX(N-burn,2);
  float x,y;
  X(0,0)=1;
  X(0,1)=0.5;
  for(int i=1;i<N;i++){
    y=X(i-1,1);
    X(i,0)=rbinom(1,n,y)[0];
    x=X(i,0);
    X(i,1)=rbeta(1,x+a,n-x+b)[0];
  }
  for(int k=0;k<N-burn;k++){
    XX(k,0)=X(k+burn,0);
    XX(k,1)=X(k+burn,1);
  }
  return XX;
}
')
res1<-gibbsc(5000,1000,2,2,10)

## -----------------------------------------------------------------------------
gibbsr<-function(N,burn,a=2,b=2,n=10){
  X <- matrix(0, N, 2) #the chain, a bivariate sample
  ###### generate the chain #####
  X[1, ] <- c(1, .5) #initialize
  for (i in 2:N) {
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, y)
    x <- X[i, 1]
    X[i, 2] <- rbeta(1,x+a,n-x+b)
    }
  b <- burn + 1
  return(X[b:N, ])
}
res2<-gibbsr(5000,1000,2,2,10)

## -----------------------------------------------------------------------------
a <- ppoints(50)
Q1x <- quantile(res1[,1], a)   #quantiles of Cauchy
Q2x <- quantile(res2[,1], a)
Q1y <- quantile(res1[,2], a)   #quantiles of Cauchy
Q2y <- quantile(res2[,2], a)
qqplot(Q1x, Q2x, main="",
        xlab="Rcpp Quantiles of x", ylab="R Quantiles of x")
    abline(c(0,0),c(1,1),col='blue',lwd=2)
qqplot(Q1y, Q2y, main="",
        xlab="Rcpp Quantiles of y", ylab="R Quantiles of y")
    abline(c(0,0),c(1,1),col='blue',lwd=2)

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(Rcpp=gibbsc(5000,1000,2,2,10),R=gibbsr(5000,1000,2,2,10))
summary(ts)[,c(1,3,5,6)]

