
d1 <- function(x, y, d, h, t0) {
  2*(x-t0)^2*exp(-(y-x)^2/2/h^2) /d^2/3*(2*pi)^(-0.5)/h
} 

FUN1 <- function(y, d, h, t0, t1) {
  integrate(d1, lower = t0, upper = t1, y = y, d = d, h=h, t0 = t0)$value
}

d2 <- function(x, y, d, h, t1){
  (1 - 4*(x-t1-d/2)^2/3/d^2) * exp(-(y-x)^2/2/h^2)*(2*pi)^(-0.5)/h
}

FUN2 <- function(y, d, h, t1, t2) {
  integrate(d2, lower = t1, upper = t2, y = y, d = d, h=h, t1 = t1)$value
}

d3 <- function(x, y, d, h, t3){
  2 * (x-t3)^2 * exp(-(y-x)^2/2/h^2)/3/d^2*(2*pi)^(-0.5)/h
}

FUN3 <- function(y, d, h, t2, t3) {
  integrate(d3, lower = t2, upper = t3, y = y, d = d, h=h, t3 = t3)$value
}

makegamma <- function(y_values, t0, d, h) {
  t1 <- t0 + d
  t2 <- t1 + d
  t3 <- t2 + d 
  gamma <- sapply(y_values, function(y) {
    FUN1(y, d, h, t0, t1) + FUN2(y, d, h, t1, t2) + FUN3(y, d, h, t2, t3)
  })
  return(gamma)
}


makec=function(y,tk,d,h){
  nc=length(y)
  nm=length(tk)-3
  cvec=1:nm
  for(j in 1:nm){
    cvec[j]=sum(makegamma(y,tk[j],d,h))/nc
  }
  cvec
}


deconmodetest <- function(y, h, v=0.01, lam=NULL, eps1=2, eps2=1, od=1/9, B=100){
  nraw=length(y)
  s1=min(y)
  s2=max(y)
  tk = seq(s1,s2,length.out = round(12*nraw^od))
  d = tk[2]-tk[1]
  m = length(tk)-3
  yp = seq(s1-3*h,s2+3*h,length.out = 2000); np = length(yp)
  yraw = y
  y = y[y>=s1 & y<=s2]
  n = length(y)
  if(is.null(lam)){
    K = kurtosis(y)
    if(K<2){
      lam = 10^2*n^(-od)
    } else if(K>2 & K<5){
      lam = 10^(4-K)*n^(-od)
    } else if(K>5 & K<9){
      lam = 10^(3/2-K/2)*n^(-od)
    } else{
      lam = 10^(-3)*n^(-od)
    }
    lam = v*lam
  }
  eps=ifelse(abs(skewness(y))>0.7,eps2,eps1)
  delta  = matrix(0,nrow=np,ncol=m)
  for(i in 1:m){
    delta[yp>tk[i]&yp<=tk[i+1],i] = 2*(yp[yp>tk[i]&yp<=tk[i+1]]-tk[i])^2/d^2/3
    delta[yp>=tk[i+1]&yp<=tk[i+2],i] = 1-4*(yp[yp>=tk[i+1]&yp<=tk[i+2]]-tk[i+1]-d/2)^2/d^2/3
    delta[yp>tk[i+2]&yp<=tk[i+3],i] = 2*(yp[yp>tk[i+2]&yp<=tk[i+3]]-tk[i+3])^2/d^2/3
  }
  gamma = matrix(0,nrow=np,ncol=m)
  for(j in 1:m){
    gamma[,j] = makegamma(yp,tk[j],d, h)
  }
  hmat = matrix(0,nrow=m,ncol=m)
  for(i in 1:m){
    for(j in 1:i)
      hmat[i,j] = hmat[j,i]=sum(gamma[,i]*gamma[,j])*(yp[2]-yp[1])
  }
  avec = 1:m*0+4*d/3
  wmat = matrix(0,nrow=m,ncol=m-1)
  for(i in 1:(m-1)){
    wmat[i,i] = -1;wmat[i+1,i]=1
  }
  b0 = 1:m*0+3/m/d/4
  cvec = makec(y,tk,d,h)
  D <- matrix(0,m+1,m)
  for(i in 1:m){
    D[i,i]=3
    D[i+1,i]=-3
  }
  for (i in 1:(m-1)) D[i,i+1]=-1
  for(i in 1:(m-1)) D[i+2,i]=1
  D = 4*D/3*sqrt(d)
  ##  get unimodal
  ans1=deconumfit(y,tk,hmat,cvec,b0,wmat,D,lam,eps=eps,od)
  ans2=deconbmfit(y,tk,hmat,cvec,b0,wmat,D,lam)
  t1=as.numeric((ans1$crit-ans2$crit))
  #t2=as.numeric((ans1$crit-ans2$crit)/abs(ans1$crit))
  fhat1=round(delta%*%ans1$bhat,10)
  fhat2=round(delta%*%ans2$bhat,10)
  # ghat1=round(gamma%*%ans1$bhat,10)
  # ghat2=round(gamma%*%ans2$bhat,10)
  dfhat2=diff(fhat2)
  outtb=NULL
  
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    pvalue=2
  }else{
    md=max(which(dfhat2>0))
    if (!is.unsorted(fhat2[1:(md-1)]) ) {
      pvalue=2
    } else {
      cdf1=delta%*%ans1$bhat
      for(i in 2:2000){
        cdf1[i]=cdf1[i-1]+cdf1[i]
      }
      cdf1=cdf1-min(cdf1)
      cdf1=cdf1/cdf1[2000]
      outtb <- foreach(t=1:B,.combine = 'rbind') %dopar%{
        yb=sapply(1:n,function(o){u=runif(1);id=min(which(u<cdf1));alp=(cdf1[id]-u)/(cdf1[id]-cdf1[id-1]);alp*yp[id-1]+(1-alp)*yp[id]})
        yb=yb+rnorm(n,0,h)#;hist(yb,breaks=30)
        deconmodet(yb,tk,d,h,hmat,b0,wmat,D,delta,lam=lam,eps=eps,od)
      }
      # outtb <- rep(0,B)
      # for(t in 1:B){
      #   yb=sapply(1:n,function(o){u=runif(1);id=min(which(u<cdf1));alp=(cdf1[id]-u)/(cdf1[id]-cdf1[id-1]);alp*yp[id-1]+(1-alp)*yp[id]})
      #   yb=yb+rnorm(n,0,h)#;hist(yb,breaks=30)
      #   outtb[t]=deconmodet(yb,tk,d,h,hmat,b0,wmat,D,delta,lam=lam,eps=eps,od)
      # }
      pvalue=sum(outtb>t1)/B
    }
  }
  ans=new.env()
  ans$yp=yp
  ans$fhat1=delta%*%ans1$bhat
  ans$fhat2=delta%*%ans2$bhat
  ans$ghat1=gamma%*%ans1$bhat
  ans$ghat2=gamma%*%ans2$bhat
  ans$statistic=t1
  ans$tb=outtb
  ans$pvalue=pvalue
  ans$lam=c(ans1$lam,ans2$lam)
  ans$kn=tk
  ans$crit=c(ans1$crit,ans2$crit)
  ans$kurtosis=kurtosis(y)
  ans$skewness=skewness(y)
  ans$truncated=(nraw!=n)
  ans
}


##############
deconumfit=function(y,tk,hmat,cvec,b0,wmat,D,lam=NULL,eps,od){	
  n <- length(y)
  m=length(tk)-3
  k1=max(min(which(tk>quantile(y,.1)))-1,4)
  k2=min(max(which(tk<quantile(y,.9)))+1,m-1)   
  amatl1=list()
  for(i in k1:k2){
    tmpm=matrix(0,m+1,m)
    tmpm[1,1]=1
    tmpm[m+1,m]=1
    for(j in 2:m){
      if(j<(i-1)){
        tmpm[j,j-1]=-1
        tmpm[j,j]=1
      }
      else{
        tmpm[j,j-1]=1
        tmpm[j,j]=-1
      }
    }
    epsvec=c(rep(0,k1-3),rep(eps/n^(od*2)/diff(range(tk))^2,i-k1+1),0,0,rep(eps/n^(od*2)/diff(range(tk))^2,k2-i+1),rep(0,m-k2))
    amatl1[[i-k1+1]]=list(amat=tmpm, epsvec=epsvec)
    
  }
  
  zvec=t(wmat)%*%(cvec-hmat%*%b0-lam*t(D)%*%D%*%b0)
  qmat=t(wmat)%*%(hmat+lam*t(D)%*%D)%*%wmat 
  crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x[[1]]%*%wmat),x[[2]]-x[[1]]%*%b0);ans$value})
  amat1 <- amatl1[[which.min(crit)]][[1]]
  epsvec <- amatl1[[which.min(crit)]][[2]]
  ans1=solve.QP(qmat,zvec,t(amat1%*%wmat),epsvec-amat1%*%b0)
  alphahat1=ans1$solution
  bhat1=wmat%*%alphahat1+b0
  cr1=t(bhat1)%*%(hmat+lam*t(D)%*%D)%*%bhat1-2*sum(cvec*bhat1)
  
  ans1=new.env()
  ans1$bhat=bhat1
  ans1$lam=lam
  ans1$crit=cr1
  ans1
}
############################################################################
deconbmfit=function(y,tk,hmat,cvec,b0,wmat,D,lam=NULL){
  n <- length(y)
  m=length(tk)-3
  m1=max(min(which(tk>quantile(y,.1)))-1,3)
  m2=min(max(which(tk<quantile(y,.9)))+1,m)
  if((m2-m1)<5){m1=max(3,m1-1);m2=min(m2+1,m)}
  trips=matrix(0,nrow=choose((m2-m1+1),3),ncol=3)
  nr=0   
  for(i in m1:(m2-4)){
    for(j in (i+2):(m2-2)){
      for(k in (j+2):(m2)){
        nr=nr+1
        trips[nr,]=c(i,j,k)
      }
    }
  }
  trips=trips[1:nr,]
  slopes=matrix(0,m+1,m)
  for(i in 1:m){slopes[i,i]=1;slopes[i+1,i]=-1}
  amat2=list()
  for(i in 1:nr){
    amat=matrix(0,m+2,m)
    amat[1:(m+1),1:m]=slopes
    amat[(trips[i,1]+1):trips[i,2],]=-slopes[(trips[i,1]+1):trips[i,2],]
    amat[(trips[i,3]):(m+1),]=-slopes[(trips[i,3]):(m+1),]
    amat[m+2,trips[i,2]]=1
    amat2[[i]]=amat
  }
  zvec=t(wmat)%*%(cvec-hmat%*%b0-lam*t(D)%*%D%*%b0)
  qmat=t(wmat)%*%(hmat+lam*t(D)%*%D)%*%wmat 
  crit<-lapply(amat2, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x%*%wmat),-x%*%b0);ans$value})
  amat2 <- amat2[[which.min(crit)]]
  ans2=solve.QP(qmat,zvec,t(amat2%*%wmat),-amat2%*%b0)
  alphahat2=ans2$solution
  bhat2=wmat%*%alphahat2+b0
  cr2=t(bhat2)%*%(hmat+lam*t(D)%*%D)%*%bhat2-2*sum(cvec*bhat2)
  ans2=new.env()
  ans2$bhat=bhat2
  ans2$lam=lam
  ans2$crit=cr2
  ans2
}
##################################
deconmodet <- function(yb,tk,d,h,hmat,b0,wmat,D,delta,lam,eps,od){
  n=length(yb)
  m=length(tk)-3
  capk=length(tk) 
  s1=min(tk)
  s2=max(tk)
  cb=makec(yb,tk,d,h)
  fit1=deconumfit(yb,tk,hmat,cb,b0,wmat,D,lam,eps=eps,od)
  fit2=deconbmfit(yb,tk,hmat,cb,b0,wmat,D,lam)
  fhat2=round(delta%*%fit2$bhat,10)
  dfhat2=diff(fhat2)
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    return (-1e+3)
  }else{md=max(which(dfhat2>0))
  if (!is.unsorted(fhat2[1:(md-1)]) ){
    return (-1e+3)
  }else{return(as.numeric(fit1$crit-fit2$crit))}
  }
}




library(quadprog)
library(moments)
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores - 2)
registerDoParallel(cl)
clusterEvalQ(cl, library(quadprog))  # Load the quadprog library in each worker
clusterExport(cl, list("deconmodetest", "deconumfit", "deconbmfit", "deconmodet", "makec", "makegamma","FUN1","FUN2","FUN3","d1","d2","d3"))
#stopCluster(cl)

n=800
h=1
mu2=4
uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
y=x+rnorm(n,0,h)
ans=deconmodetest(y,h,v=0.001,B=800)
ans$pvalue
ans$lam
ans$crit


par(mfrow=c(1,2))
hist(x,breaks=30,freq=FALSE)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
lines(ans$yp,0.6*dnorm(ans$yp,0,1)+0.4*dnorm(ans$yp,mu2,1))

hist(y,breaks=30,freq=FALSE)
lines(ans$yp,ans$ghat1,col=2)
lines(ans$yp,ans$ghat2,col=3)
f0 <- function(x, y, h) {
  (0.6 * dnorm(x, mean = 0, sd = 1) + 0.4 * dnorm(x, mean = mu2, sd = 1))*dnorm(y-x,0,h)
}
density_with_error <- sapply(ans$yp, function(y) {
  integrate(f0, lower = -Inf, upper = Inf, y=y, h=h)$value
})
lines(ans$yp,density_with_error)


#################################
### simulation for parameters ###
### v = 0.01, v = 0.001, v = 0.0001
n=500
h=1
mu2=4
pvmat <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,B=600)
  ans2=deconmodetest(y,h,v=0.001,B=600)
  ans3=deconmodetest(y,h,v=0.0001,B=600)
  pvmat[irow,] = c(ans1$pvalue, ans2$pvalue, ans3$pvalue)
}
pvmat
#              [,1]        [,2]        [,3]
#  [1,] 0.048333333 0.068333333 0.125000000
#  [2,] 0.005000000 0.020000000 0.071666667
#  [3,] 0.016666667 0.028333333 0.048333333
#  [4,] 0.035000000 0.153333333 0.255000000
#  [5,] 0.026666667 0.066666667 0.068333333
#  [6,] 0.023333333 0.081666667 0.150000000
#  [7,] 0.013333333 0.010000000 0.025000000
#  [8,] 0.081666667 0.126666667 0.136666667
#  [9,] 0.025000000 0.028333333 0.040000000
# [10,] 0.003333333 0.001666667 0.000000000
# [11,] 0.000000000 0.005000000 0.020000000
# [12,] 0.023333333 0.076666667 0.165000000
# [13,] 0.000000000 0.005000000 0.023333333
# [14,] 0.050000000 0.138333333 0.228333333
# [15,] 0.005000000 0.003333333 0.008333333
# [16,] 0.006666667 0.010000000 0.041666667
# [17,] 0.008333333 0.026666667 0.053333333
# [18,] 0.001666667 0.005000000 0.025000000
# [19,] 0.000000000 0.001666667 0.001666667
# [20,] 0.010000000 0.030000000 0.021666667
# [21,] 0.000000000 0.000000000 0.000000000
# [22,] 0.006666667 0.033333333 0.066666667
# [23,] 0.020000000 0.036666667 0.070000000
# [24,] 0.001666667 0.001666667 0.003333333
# [25,] 0.003333333 0.003333333 0.006666667
# [26,] 0.081666667 0.216666667 0.296666667
# [27,] 0.000000000 0.000000000 0.000000000
# [28,] 0.200000000 0.385000000 0.645000000
# [29,] 0.015000000 0.038333333 0.025000000
# [30,] 0.013333333 0.035000000 0.046666667
# [31,] 0.006666667 0.006666667 0.001666667
# [32,] 0.018333333 0.008333333 0.001666667
# [33,] 0.000000000 0.011666667 0.025000000
# [34,] 0.000000000 0.000000000 0.000000000
# [35,] 0.000000000 0.000000000 0.000000000
# [36,] 0.025000000 0.045000000 0.086666667
# [37,] 0.000000000 0.000000000 0.000000000
# [38,] 0.000000000 0.000000000 0.003333333
# [39,] 0.006666667 0.008333333 0.001666667
# [40,] 0.000000000 0.000000000 0.000000000
# [41,] 0.223333333 0.358333333 0.558333333
# [42,] 0.003333333 0.010000000 0.020000000
# [43,] 0.786666667 0.496666667 0.218333333
# [44,] 0.005000000 0.016666667 0.041666667
# [45,] 0.056666667 0.133333333 0.236666667
# [46,] 0.000000000 0.000000000 0.000000000
# [47,] 0.061666667 0.115000000 0.153333333
# [48,] 0.046666667 0.130000000 0.246666667
# [49,] 0.036666667 0.108333333 0.111666667
# [50,] 0.001666667 0.011666667 0.020000000
apply(pvmat,2,function(o){sum(o<=0.05)})
# 43 35 30
n=500
h=1
mu2=0
pvmat2 <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,B=600)
  ans2=deconmodetest(y,h,v=0.001,B=600)
  ans3=deconmodetest(y,h,v=0.0001,B=600)
  pvmat2[irow,] = c(ans1$pvalue, ans2$pvalue, ans3$pvalue)
}
pvmat2
# [,1]       [,2]       [,3]
# [1,]    2 2.00000000 2.00000000
# [2,]    2 0.30833333 0.23000000
# [3,]    2 2.00000000 0.63166667
# [4,]    2 2.00000000 2.00000000
# [5,]    2 2.00000000 0.54000000
# [6,]    2 0.14666667 0.13333333
# [7,]    2 0.27000000 0.21666667
# [8,]    2 2.00000000 0.71000000
# [9,]    2 2.00000000 0.30166667
# [10,]    2 2.00000000 0.76166667
# [11,]    2 0.08333333 0.04833333
# [12,]    2 2.00000000 2.00000000
# [13,]    2 0.57166667 0.47500000
# [14,]    2 2.00000000 2.00000000
# [15,]    2 2.00000000 0.40333333
# [16,]    2 2.00000000 0.56333333
# [17,]    2 2.00000000 0.57500000
# [18,]    2 0.53166667 0.48666667
# [19,]    2 0.24000000 0.25500000
# [20,]    2 2.00000000 0.63333333
# [21,]    2 2.00000000 0.52333333
# [22,]    2 2.00000000 0.64166667
# [23,]    2 2.00000000 0.63166667
# [24,]    2 0.22833333 0.12666667
# [25,]    2 2.00000000 0.48666667
# [26,]    2 2.00000000 0.20666667
# [27,]    2 2.00000000 2.00000000
# [28,]    2 0.44333333 0.32166667
# [29,]    2 2.00000000 0.63333333
# [30,]    2 0.27000000 0.33666667
# [31,]    2 2.00000000 2.00000000
# [32,]    2 0.28333333 0.20833333
# [33,]    2 2.00000000 0.18500000
# [34,]    2 2.00000000 0.73500000
# [35,]    2 2.00000000 0.81666667
# [36,]    2 0.36833333 0.39833333
# [37,]    2 2.00000000 0.45166667
# [38,]    2 0.17500000 0.16166667
# [39,]    2 0.60500000 0.57833333
# [40,]    2 0.35333333 0.42166667
# [41,]    2 0.57333333 0.58833333
# [42,]    2 0.44333333 0.59666667
# [43,]    2 2.00000000 0.45833333
# [44,]    2 0.52000000 0.48833333
# [45,]    2 0.51500000 0.70500000
# [46,]    2 2.00000000 0.36000000
# [47,]    2 2.00000000 0.75333333
# [48,]    2 2.00000000 2.00000000
# [49,]    2 0.18000000 0.08000000
# [50,]    2 2.00000000 0.79166667
# 0 0 1

n=500
h=1
mu2=4.5
pvmat3 <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,B=600)
  ans2=deconmodetest(y,h,v=0.001,B=600)
  #ans3=deconmodetest(y,h,v=0.0001,B=600)
  pvmat3[irow,] = c(ans1$pvalue, ans2$pvalue, 0)
}
pvmat3
# 49 48

### eps = 1, 2 ; eps = 10. 20; eps = 0.1, 0.2
n=500
h=1
mu2=4
pvmat4 <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,eps1=2,eps2=1,B=600)
  ans2=deconmodetest(y,h,v=0.01,eps1=1,eps2=0.5,B=600)
  ans3=deconmodetest(y,h,v=0.01,eps1=0.5,eps2=0.25,B=600)
  pvmat4[irow,] = c(ans1$pvalue, ans2$pvalue, ans3$pvalue)
}
pvmat4
# 37 36 36
n=500
h=1
mu2=0
pvmat5 <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,eps1=2,eps2=1,B=600)
  ans2=deconmodetest(y,h,v=0.01,eps1=1,eps2=0.5,B=600)
  ans3=deconmodetest(y,h,v=0.01,eps1=0.5,eps2=0.25,B=600)
  pvmat5[irow,] = c(ans1$pvalue, ans2$pvalue, ans3$pvalue)
}
pvmat5
# 0 0 0

n=500
h=1
mu2=4
pvmat6 <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,B=700)
  ans2=deconmodetest(y,h,v=0.1,B=700)
  ans3=deconmodetest(y,h,v=1,B=700)
  pvmat6[irow,] = c(ans1$pvalue, ans2$pvalue, ans3$pvalue)
}
pvmat6
apply(pvmat6,2,function(o){sum(o<=0.05)})
# 34 36 0

n=500
h=1
mu2=0
pvmat7 <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,B=700)
  ans2=deconmodetest(y,h,v=0.1,B=700)
  ans3=deconmodetest(y,h,v=1,B=700)
  pvmat7[irow,] = c(ans1$pvalue, ans2$pvalue, ans3$pvalue)
}
pvmat7
apply(pvmat7,2,function(o){sum(o<=0.05)})
# 0 0 0

n=500
h=1
mu2=4
pvmat8 <- matrix(0,50,3)
for(irow in 1:50){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),mu2,1)
  y=x+rnorm(n,0,h)
  ans1=deconmodetest(y,h,v=0.01,B=700)
  ans2=deconmodetest(y,h,v=0.03,B=700)
  ans3=deconmodetest(y,h,v=0.06,B=700)
  pvmat8[irow,] = c(ans1$pvalue, ans2$pvalue, ans3$pvalue)
}
pvmat8
apply(pvmat8,2,function(o){sum(o<=0.05)})
# 0 0 0


t1=Sys.time()
ans1=deconmodetest(y,h,v=0.01,B=600)
t2=Sys.time()
t2-t1


n=500
mu2=4
uu=runif(n);x=rnorm(n,0,1);x[uu<0.5]=rnorm(sum(uu<0.5),mu2,1)
y=x+rnorm(n,0,h)
yp=seq(-5,9, length.out = 1000)
plot(yp,0.5*dnorm(yp,0,1)+0.5*dnorm(yp,mu2,1),type="l")
yp=seq(1.99,2.01, length.out = 1000)
h=1.732051
f0 <- function(x, y, h) {
  (0.5 * dnorm(x, mean = 0, sd = 1) + 0.5 * dnorm(x, mean = mu2, sd = 1))*dnorm(y-x,0,h)
}
density_with_error <- sapply(yp, function(y) {
  integrate(f0, lower = -Inf, upper = Inf, y=y, h=h)$value
})
plot(yp, density_with_error, type="l")


n=500
mu2=4


h=0.5 #100/100
h=1 #95/100
h = 1.5#58/200  #26/77 v=0.0001
h=2 #8/200
pvec = rep(0,1000)
pvec3 = rep(0,100)
t1=Sys.time()
for(i in 1:100){
  uu=runif(n);x=rnorm(n,0,1);x[uu<0.5]=rnorm(sum(uu<0.5),mu2,1)
  y=x+rnorm(n,0,h)
  ans=deconmodetest(y,h,v=0.0001,B=500)
  pvec3[i]=ans$pvalue
}
t2=Sys.time()
t2-t1
sum(pvec3[1:77]<=0.05)

h=2
uu=runif(n);x=rnorm(n,0,1);x[uu<0.5]=rnorm(sum(uu<0.5),mu2,1)
y=x+rnorm(n,0,h)
ans=deconmodetest(y,h,v=0.00001,B=1)
par(mfrow=c(1,2))
hist(x,breaks=30,freq=FALSE)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
lines(ans$yp,0.5*dnorm(ans$yp,0,1)+0.5*dnorm(ans$yp,mu2,1))
rug(ans$kn)
hist(y,breaks=30,freq=FALSE)
lines(ans$yp,ans$ghat1,col=2)
lines(ans$yp,ans$ghat2,col=3)
f0 <- function(x, y, h) {
  (0.5 * dnorm(x, mean = 0, sd = 1) + 0.5 * dnorm(x, mean = mu2, sd = 1))*dnorm(y-x,0,h)
}
density_with_error <- sapply(ans$yp, function(y) {
  integrate(f0, lower = -Inf, upper = Inf, y=y, h=h)$value
})
lines(ans$yp,density_with_error)
rug(ans$kn)
#################################

plot(yp,delta[,1],type = "l")
for(i in 2:(ncol(delta))){
  lines(yp,delta[,i])
}

plot(yp,gamma[,1],type = "l")
for(i in 2:ncol(gamma)){
  lines(yp,gamma[,i])
}