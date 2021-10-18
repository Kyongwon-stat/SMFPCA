#################################################
# Functions for FGGM, ploting ROC curves etc
#################################################
library('fda')
library(MASS)
#library(splines)
library(Matrix)
library('BB')
library(MFPCA)

###############################################################
#          SUBROUTINE: Moore-Penrose type power               #
#           keep only the first m eigenvalues                 #
###############################################################
mppower1 = function(matrix,power,m){
eig = eigen(matrix)
eval = eig$values
evec = eig$vectors
tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
t(evec[,1:m])
return(tmp)
}
###############################################################
#          Subroutine: estimate gamma                        #
###############################################################
estgam = function(x){
n = nrow(x)
sig = 0
for(i in 1:(n-1)){
for(j in (i+1):n){
d = x[i,]-x[j,]
sig = sig + sqrt(t(d)%*%d)}}
sig = sig/(n*(n-1)/2)
gamm = 1/(2*(sig^2))
return(c(gamm))}

###################################################
#  Moore-Penrose type power           #
#  Taking power ignoring 0 eigenvalues;           #
#    ignoring criterion=ignore                    #
###################################################
mppower = function(matrix,power,ignore){
eig = eigen(matrix)
eval = eig$values
evec = eig$vectors
m = length(eval[abs(eval)>ignore])
tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
t(evec[,1:m])
return(tmp)
}
#######################################################
# Power of a matrix
######################################
matpower = function(a,alpha){
a = (a + t(a))/2
tmp = eigen(a)
return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
t(tmp$vectors))}
#########################################
# First level functional space H
#Function that generates the data from a RKHS
# using Gaussian radial basic kernel
# and covariance function kernel
#####################################################
# function that generates the eigenfunctions
# for Brownian covariance kernel function
 ###################################
fpcawiener <- function(t, j){
 sqrt(2) * sin((j-1/2)*pi*t)}
 ###################################
# function that generates an iid sample from X(t)
# X is one-dimensional
# using Gaussian radial basis function and Brownian
# covarance kernel functon
# p is only for brownian.
 ###################################
generate.function <- function(n, nt, N=10, gamma=7, p=100, sd=1, fun="brownian"){
if(fun=="gaussian"){
t1 <- 0.1#;t2<-0.6;t3<-0.9
t11 <- 0.9#;t22<-0.5;t33<-0.8
tt <- seq(0,1, len=nt)
x <- matrix(0,n,ntx)   # observed x
ttt <- NULL # it is for the function x
a <- NULL
for(i in 1:n){
xi <- rep(0, nt)
a[[i]] <- rnorm(N)
ttt[[i]] <- runif(N)
for(k in 1:N){
xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)}
x[i,] <- xi}
result <- list(x=x, tt=tt)}
if(fun=="brownian"){
nnt <- nt * 5
tt <-(1:nnt)/nnt
BM <- matrix(0, n, nnt)
a <- matrix(rnorm(n*p,sd=sd),n,p)
for(j in 1: p) {
a[,j] <- a[,j]/((j-1/2)*pi)
for(i in 1:n){
BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)}
}
#unbalanced
x.ub <- matrix(0, n, nt)
t.ub <- matrix(0, n, nt)
for(i in 1:n){
t.ub[i, ] <- sort(sample(1:nnt, nt))
x.ub[i, ] <- BM[i, t.ub[i,]]
}
t.ub <- t.ub / nnt
#balanced
t.b <- round(seq(1, nnt, len=nt))
x.b <- BM[,t.b]
t.b <- matrix(t.b, nrow=n, ncol=nt, byrow=T)/ nnt
result <- list(x.ub=x.ub, x.b=x.b, t.ub=t.ub, t.b=t.b, a=a)}
return(result)}

###############################################################
# function:  Gram matrix for first-level Hilbert space
#            for a given kernel and time set
# tte: is the time set for evaluation
# tto: the time set for observations
# if you just want to compute gram matrix then set tte=tto
#  kern is either brown or gauss
###############################################################
gramt=function(tte,tto,kern){
  ltte=length(tte);ltto=length(tto)
  if (kern=="gauss"){
    a1=matrix(tte^2,ltte,ltto);a2=tte%*%t(tto);a3=t(matrix(tto^2,ltto,ltte))
    a=a1-2*a2+a3
    b1=matrix(tto^2,ltto,ltto);b2=tto%*%t(tto);b3=t(matrix(tto^2,ltto,ltto))
    b=b1-2*b2+b3
    sigma=sum(sqrt(b))/(2*choose(ltto,2));gamma=1/(2*sigma^2)
    ktmat=exp(-gamma*(a))}
  if(kern=="brown"){
    arr=array(0,c(ltte,ltto,2))
    arr[,,1]=matrix(tte,ltte,ltto);arr[,,2]=t(matrix(tto,ltto,ltte))
    ktmat=apply(arr,c(1,2),min)}
  return(t(ktmat))
}


###############################################################
# function: estimates one function
###############################################################
evalx=function(f,tte,tto,ridge,kern){
  kt=gramt(tto,tto,kern)
  scale=eigen(kt)$values[1]
  ktinv=matpower(kt+scale*ridge*diag(nrow(kt)),-1)
  kt1=t(gramt(tte,tto,kern))
  out=kt1%*%ktinv%*%f
  return(c(out))}
###############################################################
# function: estimates a sample of functions
###############################################################
evalxmat=function(ff,tte,ttt,ridge,kern){
  n=dim(ff)[1];ffcoo=numeric()
  for(i in 1:n) ffcoo=rbind(ffcoo,evalx(ff[i,],tte,ttt[i,],ridge,kern))
  return(ffcoo)
}


############################################################
###### Miscellaneous function for fsgm ####################
############################################################
CalGam=function(A){
  if(is.vector(A)){
    n = length(A)}
  else n = dim(A)[1]
  tmp=rowSums(as.matrix(A*A))%*%t(rep(1,n))
  K=as.numeric(tmp+t(tmp)-2*A%*%t(A))
  K=K*(K>=0)
  #tou=as.real(sum(sqrt(K))/(n*(n-1)))
  tou=sum(sqrt(K))/(n*(n-1))
  gam=1/(2*tou^2)
  return(gam)
}

KGaussian=function(gamma,A,B){
  if(is.vector(A)){
    n = length(A)}
  else n = dim(A)[1]
  if(is.vector(B)){
    m = length(B)}
  else m = dim(B)[1]
  tmp_1=rowSums(as.matrix(A*A))%*%matrix(1,1,m)
  tmp_2=rowSums(as.matrix(B*B))%*%matrix(1,1,n)
  K=tmp_1+t(tmp_2)-2*A%*%t(B)
  K=exp(-K*gamma)
  return(K)
}

KBrown=function(tte,tto){
  ltte=length(tte);ltto=length(tto)
  arr=array(0,c(ltte,ltto,2))
  arr[,,1]=matrix(tte,ltte,ltto);arr[,,2]=t(matrix(tto,ltto,ltte))
  ktmat=apply(arr,c(1,2),min)
}



standvec = function(x){
  return((x - mean(x))/sd(x))}


###############################################################
#                      3/8/2017   fapo
# use k_T to generate basis functions; use L2 as the inner
# product; use simpson's rule to compute L2 inner product
# contain both the balanced case and unbalanced case
# for either the balanced or the unbalanced case, use ttt from
# the output of generatex-p30-n300.R; in the balanced case
# ttt contains identical rows; in the unbalanced case, ttt contains
# randomly generated rows.
###############################################################
###############################################################
# function: generate simpson's weights (ns must be odd)
###############################################################
simpson=function(ns){
  wei=rep(0,ns);wei[1]=1;wei[ns]=1
  for(i in 1:((ns-1)/2)) wei[2*i]=4
  for(i in 1:((ns-1)/2-1)) wei[2*i+1]=2
  h=1/(ns-1); wei=wei*(h/3)
  return(wei)
}
###############################################################
# function: operator norm
###############################################################
onorm=function(a) return(eigen(round((a+t(a))/2,8))$values[1])
##############################################################
#      function: trace of a matrix
##############################################################
tr=function(a) return(sum(diag(a)))
###############################################################
#          function: Moore-Penrose type power
#           keep only the first m eigenvalues
###############################################################
mppower1 = function(matrix,power,m){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
    t(evec[,1:m])
  return(tmp)
}
###############################################################
#       function: Moore-Penrose type power
#  ignoring eigenvalues  less than a threshold
###############################################################
mppower = function(matrix,power,ignore){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
    t(evec[,1:m])
  return(tmp)
}
###############################################################
#          function: power of a matrix
###############################################################
matpower = function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))}
###############################################################
# function:  Gram matrix for first-level Hilbert space
#            for a given kernel and time set
# tte: is the time set for evaluation
# tto: the time set for observations
# if you just want to compute gram matrix then set tte=tto
#  kern is either brown or gauss
###############################################################
gramt=function(tte,tto,kern){
  ltte=length(tte);ltto=length(tto)
  if (kern=="gauss"){
    a1=matrix(tte^2,ltte,ltto);a2=tte%*%t(tto);a3=t(matrix(tto^2,ltto,ltte))
    a=a1-2*a2+a3
    b1=matrix(tto^2,ltto,ltto);b2=tto%*%t(tto);b3=t(matrix(tto^2,ltto,ltto))
    b=b1-2*b2+b3
    sigma=sum(sqrt(b))/(2*choose(ltto,2));gamma=1/(2*sigma^2)
    ktmat=exp(-gamma*(a))}
  if(kern=="brown"){
    arr=array(0,c(ltte,ltto,2))
    arr[,,1]=matrix(tte,ltte,ltto);arr[,,2]=t(matrix(tto,ltto,ltte))
    ktmat=apply(arr,c(1,2),min)}
  return(t(ktmat))
}
##############################################################
#   function:  Generalized Cross Validation for expsilon_T
#  et is a grid of values
# xxx are the observed random functions
##############################################################
gcvt=function(xxx,ttt,et,kern){
  n=dim(xxx)[1];nt=dim(xxx)[2];p=dim(xxx)[3]
  nuset=numeric();deset=numeric()
  for(i in 1:n){
    kt=gramt(ttt[i,],ttt[i,],kern)
    scale=onorm(kt)
    ktinv=matpower(kt+et*scale*diag(nt),-1)
    for(j in 1:p) {
      nuset=c(nuset,sum((xxx[i,,j]-kt%*%ktinv%*%xxx[i,,j])^2))
      deset=c(deset,(1-tr(kt%*%ktinv)/nt)^2)}}
  out=sum(nuset/deset)
  return(out)
}
###############################################################
# function: estimates one function
###############################################################
evalx=function(f,tte,tto,ridge,kern){
  kt=gramt(tto,tto,kern)
  scale=eigen(kt)$values[1]
  ktinv=matpower(kt+scale*ridge*diag(nrow(kt)),-1)
  kt1=t(gramt(tte,tto,kern))
  out=kt1%*%ktinv%*%f
  return(c(out))}
###############################################################
# function: estimates a sample of functions
###############################################################
evalxmat=function(ff,tte,ttt,ridge,kern){
  n=dim(ff)[1];ffcoo=numeric()
  for(i in 1:n) ffcoo=rbind(ffcoo,evalx(ff[i,],tte,ttt[i,],ridge,kern))
  return(ffcoo)
}
###############################################################
# function: gram matrix for the second RHKS for a sample of x
# x are the estimated functions at the evaluation points
# kt is the diagonal of simpson weights to compute the L2 inner product
###############################################################
gramx=function(x,kt){
  n=dim(x)[1]
  k2=x%*%kt%*%t(x);k1=t(matrix(diag(k2),n,n));k3=t(k1);k=k1-2*k2+k3
  sigma=sum(sqrt(k))/(2*choose(n,2));gamma=1/(2*sigma^2)
  return(exp(-gamma*(k1-2*k2+k3)))
}
###############################################################
# function: Q = I - J/n
###############################################################
qmat = function(n) return(diag(n)-rep(1,n)%*%t(rep(1,n))/n)
###############################################################
#  function:  K to G=QKQ
###############################################################
cgram=function(K){n=dim(K)[1];Q=qmat(n);return(Q%*%K%*%Q)}
###############################################################
# function: Compute matrix A_i and Lambda_i
###############################################################
aili=function(Gi,ridge){
  n=dim(Gi)[1];scale=eigen(Gi)$values[1];Q=qmat(n)
  mat=Gi+ridge*scale*Q
  Ai=(n^(-1/2))*mppower(Gi,1/2,10^(-7))%*%mppower(mat,-1/2,10^(-7))
  Li=Q-Ai%*%Ai
  return(list(Ai=Ai,Li=Li))
}


###############################################################
# This function is from sparsePCA(Chavent) package. We edit 
# it with our own purpose. 
###############################################################
groupsparsePCA <- function(A, m, lambda, index=1:ncol(A), block=1, mu=1/1:m,post=FALSE,center=TRUE,scale=TRUE,init=NULL)
{
  n <- nrow(A)
  pp <- ncol(A) #number of scalar variables
  if (m > min(n-1,pp))
    stop("m must be smaller than rank(A)",call. = FALSE)
  if (length(lambda)!=m)
    stop("lambda must be a vector of size m",call. = FALSE)
  if ((max(lambda) >1) || (max(lambda) < 0))
    stop("Values in lambda must be in [0,1]",call. = FALSE)
  if (length(index) !=ncol(A))
    stop("the length of index is not correct",call. = FALSE)

  norm2 <- function(x) sqrt(sum(x^2))
  polar <- function(x)
  {
    obj <- svd(x)
    obj$u %*% t(obj$v)
  }

  p <- length(unique(index)) #number of groups variables

  iter_max <- 1000   # maximum number of admissible iterations
  epsilon <- 0.0001  # accuracy of the stopping criterion

  Z <- matrix(0,ncol(A),m)
  rownames(Z) <- colnames(A)
  A <- as.matrix(A)
  index <- as.factor(index)

  if (center==TRUE)
    A <- scale(A,center=TRUE, scale=FALSE) #center data
  if (scale==TRUE)
    A <- scale(A,center=FALSE,scale=TRUE)*sqrt(n/(n-1)) #scale data

  if ((m==1) || (m>1 && block==0)) #single-unit algorithm (deflation is used if m>1)
  {
    B <- A
    X <- matrix(NA,n,m)
    for (comp in 1:m)            #loop on the components
    {
      ai <- lapply(split(data.frame(t(B)),index),t) # list of group variables (matrices)
      i_max <- which.max(lapply(ai,function(x){norm(x,"2")})) #norm of a matrix is the largest singular value
      gamma_max <- norm(ai[[i_max]],type="2")
      gamma <- lambda[comp]*gamma_max
      x <-svd(B)$u[,1]     #initialization point
      f <- matrix(0,iter_max,1)
      iter <- 1
      S <- function(v)
      {
        alpha=sqrt(sum(v^2))
        if (alpha > 0) v/alpha*((alpha>=gamma)*(alpha-gamma))
      }
      repeat
      {
        Ax <- lapply(ai,function(a){crossprod(a,x)})
        tresh <-  lapply(Ax,S)
        f[iter] <- sum(unlist(lapply(tresh,function(v){sum(v^2)})))  #cost function
        if (f[iter]==0) #the sparsity parameter is too high: all entries of the loading vector are zero
          break
        grad <-  B%*%unlist(tresh)
        x <- grad/norm2(grad)
        if (((iter>=2) && ((f[iter]-f[iter-1])/f[iter-1] < epsilon)) || (iter >= iter_max)) #stopping criterion
        {
          if (iter >= iter_max) print("Maximum number of iterations reached")
          break
        }
        iter <- iter+1
      }
      Ax <- lapply(ai,function(a){crossprod(a,x)})
      #pattern <- abs(Ax)-gamma >0
      tresh <-  lapply(Ax,S)
      z <- unlist(tresh)
      if (max(abs(z))>0)
        z <- z/norm2(z)
      #if (post==TRUE) z <- pattern_filling(A,z)
      y <- B%*%z
      B <- B-y%*%t(z)
      Z[,comp]<-z
      X[,comp] <- y/norm2(y)
    }
  }
  if (m>1 && block==1) # block algorithm
  {
    e <- svd(A)
    d <- e$d[1:m]
    u <- e$u[,1:m]
    ai <- lapply(split(data.frame(t(A)),index),t) # list of group variables
    i_max <- which.max(lapply(ai,function(x){norm(x,"2")}))
    gamma_max <- norm(ai[[i_max]],type="2")
    gamma <- lambda*gamma_max*d/d[1]
    #initialization
    if (!is.null(init))
    {
      x <- A %*% init %*% diag(mu ^ 2)
      x <- polar(x)
    } else
      x <- u

    f <- matrix(0,iter_max,1)
    iter <- 1
    S <- function(v)
    {
      alpha <- apply(v,2,norm2)
      sweep(v,2,alpha,"/")%*%diag((alpha >= gamma)*(alpha-gamma))
    }
    repeat
    {
      Ax <- lapply(ai,function(a){crossprod(a,x)})
      tresh <-  lapply(Ax,S)
      #f[iter] <- sum(unlist(lapply(tresh,function(v){sum(v^2)})))  #cost function
      f[iter] <- sum(unlist(lapply(tresh,function(v){sum(apply(v,2,norm2)^2*mu^2)})))  #cost function
      if (f[iter]==0) #the sparsity parameter is too high: all entries of the loading vector are zero
        break
      for(i in seq(length(tresh))) tresh2 <- if(i == 1) tresh[[i]] else rbind(tresh2,tresh[[i]])
      grad <-  2*A%*%tresh2 %*% diag(mu^2)
      x <- polar(grad)
      if (((iter>=2) && ((f[iter]-f[iter-1])/f[iter-1] < epsilon)) || (iter >= iter_max)) #stopping criterion
      {
        if (iter >= iter_max) print("Maximum number of iterations reached")
        break
      }
      iter <- iter+1
    }
    Ax <- lapply(ai,function(a){crossprod(a,x)})
    tresh <-  lapply(Ax,S)
    for(i in seq(length(tresh))) tresh2 <- if(i == 1) tresh[[i]] else rbind(tresh2,tresh[[i]])
    z <- tresh2
    for (j in 1:m)
    {
      if (max(abs(z[,j]))>0)
        Z[,j] <- z[,j]/norm2(z[,j])
    }
    #if (post==TRUE) Z <- pattern_filling(A,Z,mu)
    X=x
  }
  Y <- A %*% Z
  res <- list(Z=Z,Y=Y,gamma=gamma,B=A,X=X)
  class(res) <- "sparsePCA"
  return(res)
}


###############################################################################
#      Example
###############################################################################

#n=100# sample size
#ntx=10 # number of  time points to sample from each predictor
#p=5 # number of dimension
#N=50
#K=7
#M=5

#sigma=0.1 # for rnorm when generating sample

# observed time points for balanced
#nnt <- ntx * 5
#tx <- round(seq(1, nnt, len=ntx))
#tx <- matrix(tx, nrow=n, ncol=ntx, byrow=T)/ nnt


#X=array(0,c(n,ntx,p))
#X.ns=array(0,c(n,ntx,p))
#alpha=array(0,c(n,N,p))
#for (j in 1:p){
#  alpha[,,j]=generate.function(n,nt=ntx,p=N,fun="brownian")[[5]]}
#for (k in 1:N){
#  for (i in 1:n){
#    for(ll in 1:p){
#      X[i,,ll]=X[i,,ll]+(alpha[i,k,ll])*fpcawiener(tx[ll,],1)
#    }
#  }
#}

#X[,,1]=X[,,1]+rnorm(n*ntx, sd=.5)
#X[,,4]=X[,,4]+rnorm(n*ntx, sd=.5)
#X[,,2]=X[,,2]+rnorm(n*ntx, sd=.5)
#X[,,5]=X[,,1]+X[,,2]+X[,,4]+rnorm(n*ntx, sd=10)
#X[,,3]=X[,,3]+rnorm(n*ntx, sd=.5)




# stack X by using rbind
#funcx=numeric()
#for (k in 1:p){
#  funcx=rbind(funcx,(X[,,k]))
#}

#dimen <- dim(funcx)[[2]]
#size <- dim(funcx)[[1]]

# produce splines functions defined on interval [0,1]
#databasis=create.fourier.basis(rangeval=c(0,1),nbasis=5)


# In functional data
#x.fd=Data2fd(tx[1,],t(funcx),databasis)

#xf <- t(x.fd$coefs)

# make data matrix

#x <- numeric()
#for(k in 1:p){
#  x=cbind(x, xf[((k-1)*n+1):(k*n),])
#}


############################################################
# result from SMFPCA method
############################################################

#index <- rep(c(1:p),c(rep(5,p)))
#Z <- groupsparsePCA(x,m=5,c(rep(0.6,5)),index)$Z 

#res.plot1 <- fd(Z[c(1:5),],databasis)
#res.plot2 <- fd(Z[c(6:10),],databasis)
#res.plot3 <- fd(Z[c(11:15),],databasis)
#res.plot4 <- fd(Z[c(16:20),],databasis)
#res.plot5 <- fd(Z[c(21:25),],databasis)
#par(mfrow=c(3,5))
#plot(res.plot1)
#title("i=1")
#plot(res.plot2)
#title("i=2")
#plot(res.plot3)
#title("i=3")
#plot(res.plot4)
#title("i=4")
#plot(res.plot5)
#title("i=5")








