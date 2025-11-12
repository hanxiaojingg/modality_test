
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

t0=d
t1=t0+d
t2=t1+d
t3=t2+d
makegamma(y_values,t0,d)

makec=function(y,tk,d,h){
  nc=length(y)
  nm=length(tk)-3
  cvec=1:nm
  for(j in 1:nm){
    cvec[j]=sum(makegamma(y,tk[j],d,h))/nc
  }
  cvec
}



nraw=length(y)
# q1=quantile(y,.01)
# q2=quantile(y,.99)
# rng=q2-q1
# s1=as.numeric(q1-.4*rng)
# s2=as.numeric(q2+.4*rng)
# l = round((10*n^(1/9))*h/(s2-s1))
# l = 2
# d = h/l 
# ran=s2-s1
# s1=s1-(ceiling(ran/d)*d-ran)/2
# s2=s2+(ceiling(ran/d)*d-ran)/2
s1=min(x)
s2=max(x)
tk = seq(s1,s2,length.out = round((10*n^(1/9))))
d = tk[2]-tk[1]
m = length(tk)-3
yp = seq(s1-3*h,s2+3*h,length.out = 2000); np = length(yp)
yraw = y
y = y[y>=s1 & y<=s2]
n = length(y)
  K = kurtosis(y)
  if(K<2){
    lam = 10^2*n^(-1/7)
  } else if(K>2 & K<5){
    lam = 10^(4-K)*n^(-1/7)
  } else if(K>5 & K<9){
    lam = 10^(3/2-K/2)*n^(-1/7)
  } else{
    lam = 10^(-3)*n^(-1/7)
  }
  v=0.01
  lam = v*lam
  eps1=2;eps2=1
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
od=2/9
ans1=deconumfit(y,tk,hmat,cvec,b0,wmat,D,lam,eps=eps,od)
ans2=deconbmfit(y,tk,hmat,cvec,b0,wmat,D,lam)
t1=as.numeric((ans1$crit-ans2$crit))
#t2=as.numeric((ans1$crit-ans2$crit)/abs(ans1$crit))
fhat1=round(delta%*%ans1$bhat,10)
fhat2=round(delta%*%ans2$bhat,10)
ghat1=round(gamma%*%ans1$bhat,10)
ghat2=round(gamma%*%ans2$bhat,10)
dfhat2=diff(fhat2)
outtb=NULL
par(mfrow=c(1,2))
hist(x,breaks=30,freq=FALSE)
lines(yp,fhat1,col=2)
lines(yp,fhat2,col=3)
lines(yp,0.6*dnorm(yp,0,1)+0.4*dnorm(yp,3,1))

hist(y,breaks=30,freq=FALSE)
lines(yp,ghat1,col=2)
lines(yp,ghat2,col=3)
f0 <- function(x, y, h) {
 (0.6 * dnorm(x, mean = 0, sd = 1) + 0.4 * dnorm(x, mean = 3, sd = 1))*dnorm(y-x,0,h)
}
density_with_error <- sapply(yp, function(y) {
  integrate(f0, lower = -Inf, upper = Inf, y=y, h=h)$value
})
lines(yp,density_with_error)


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
        yb=yb+runif(n,-h,h)#;hist(yb,breaks=30)
        deconmodet(yb,tk,d,h,hmat,b0,wmat,D,delta,lam=lam,eps=eps,od)
      } 
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
    epsvec=c(rep(0,k1-3),rep(eps/n^od/diff(range(tk))^2,i-k1+1),0,0,rep(eps/n^od/diff(range(tk))^2,k2-i+1),rep(0,m-k2))
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
n=800
h=1
uu=runif(n);x=rnorm(n,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),3,1)
y=x+rnorm(n,0,h)
ans=deconmodetest(y,h)
ans$pvalue
hist(y,breaks=40,freq=FALSE)
lines(ans$yp,ans$ghat1,col=2)
lines(ans$yp,ans$ghat2,col=3)
lines(ans$yp,(0.6*pnorm(ans$yp+h,0,1)+0.4*pnorm(ans$yp+h,4,1))/2/h-(0.6*pnorm(ans$yp-h,0,1)+0.4*pnorm(ans$yp-h,4,1))/2/h)
hist(x,breaks=40,freq=FALSE)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
lines(ans$yp,0.6*dnorm(ans$yp,0,1)+0.4*dnorm(ans$yp,4,1))


plot(yp,delta[,1],type = "l")
for(i in 2:(ncol(delta))){
  lines(yp,delta[,i])
}

plot(yp,gamma[,1],type = "l")
for(i in 2:ncol(gamma)){
  lines(yp,gamma[,i])
}


