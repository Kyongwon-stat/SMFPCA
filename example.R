
n=100# sample size
ntx=10 # number of  time points to sample from each predictor
p=5 # number of dimension
N=50
K=7
M=5

sigma=0.1 # for rnorm when generating sample

# observed time points for balanced
nnt <- ntx * 5
tx <- round(seq(1, nnt, len=ntx))
tx <- matrix(tx, nrow=n, ncol=ntx, byrow=T)/ nnt


X=array(0,c(n,ntx,p))
X.ns=array(0,c(n,ntx,p))
alpha=array(0,c(n,N,p))
for (j in 1:p){
  alpha[,,j]=generate.function(n,nt=ntx,p=N,fun="brownian")[[5]]}
for (k in 1:N){
  for (i in 1:n){
    for(ll in 1:p){
      X[i,,ll]=X[i,,ll]+(alpha[i,k,ll])*fpcawiener(tx[ll,],1)
    }
  }
}

X[,,1]=X[,,1]+rnorm(n*ntx, sd=.5)
X[,,4]=X[,,4]+rnorm(n*ntx, sd=.5)
X[,,2]=X[,,2]+rnorm(n*ntx, sd=.5)
X[,,5]=X[,,1]+X[,,2]+X[,,4]+rnorm(n*ntx, sd=10)
X[,,3]=X[,,3]+rnorm(n*ntx, sd=.5)




# stack X by using rbind
funcx=numeric()
for (k in 1:p){
  funcx=rbind(funcx,(X[,,k]))
}

dimen <- dim(funcx)[[2]]
size <- dim(funcx)[[1]]

# produce splines functions defined on interval [0,1]
databasis=create.fourier.basis(rangeval=c(0,1),nbasis=5)


# In functional data
x.fd=Data2fd(tx[1,],t(funcx),databasis)

xf <- t(x.fd$coefs)

# make data matrix

x <- numeric()
for(k in 1:p){
  x=cbind(x, xf[((k-1)*n+1):(k*n),])
}


############################################################
# result from fca method
############################################################

index <- rep(c(1:p),c(rep(5,p)))
Z <- groupsparsePCA(x,m=5,c(rep(0.6,5)),index)$Z # block different mu

res.plot1 <- fd(Z[c(1:5),],databasis)
res.plot2 <- fd(Z[c(6:10),],databasis)
res.plot3 <- fd(Z[c(11:15),],databasis)
res.plot4 <- fd(Z[c(16:20),],databasis)
res.plot5 <- fd(Z[c(21:25),],databasis)
par(mfrow=c(3,5))
plot(res.plot1)
title("i=1")
plot(res.plot2)
title("i=2")
plot(res.plot3)
title("i=3")
plot(res.plot4)
title("i=4")
plot(res.plot5)
title("i=5")
