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
# To implement the group sparse maximum variance approach in 
# Chavent and Chavent (2021), the following r package 
# is required.
###############################################################

devtools::install_github("chavent/sparsePCA")






