library("splines2")
library("quadprog")
library("doParallel")
library("diptest")
library("moments")
library("multimode")
bmodetest <- function(y,lower = NULL, upper = NULL,B=1000,lam=NULL,eps1=5,eps2=2,eps=NULL){
  nraw=length(y)
  s1=lower
  s2=upper
  q1=quantile(y,.01)
  q2=quantile(y,.99)
  rng=q2-q1
  if(is.null(lower)){s1=as.numeric(q1-.4*rng)}
  if(is.null(upper)){s2=as.numeric(q2+.4*rng)}
  capk=round(nraw^(1/7)*12)
  qy=as.numeric(quantile(y,1:capk/(capk+1)))
  ####add more knots on both sides
  if(s1<qy[1]){
    k1=min(floor(log((qy[1]-s1)/3/(qy[2]-qy[1])+1, 3/2) -1 ), round(capk/8))
    if(k1<2){
      qy = c(qy,s1)
      if((qy[1]-s1)/(qy[2]-qy[1])>=2) qy = c(qy,qy[1]-(qy[1]-s1)/2)
    }else{
      q1 = ((qy[1]-s1)/(qy[2]-qy[1]))^(1/(k1+1))
      qy = c(qy, qy[1] - (qy[2]-qy[1])*q1^(1:(k1)), s1)
    }
  }
  if(s2>qy[capk]){
    k2=min(floor(log((s2-qy[capk])/3/(qy[capk]-qy[capk-1])+1, 3/2) - 1), round(capk/8))
    if(k2<2){
      qy=c(qy,s2)
      if((s2-qy[capk])/(qy[capk]-qy[capk-1])>=2) qy = c(qy,(s2-qy[capk])/2+qy[capk])
    }else{
      q2 = ((s2-qy[capk])/(qy[capk]-qy[capk-1]))^(1/(k2+1))
      qy = c(qy, (qy[capk]-qy[capk-1])*q2^(1:(k2))+qy[capk], s2)
    }
  }
  qy = sort(qy)
  ####add more knots between the two modes
  d=min(diff(qy))
  dd=floor(diff(qy)/min(diff(qy)))
  m1=min(which(dd==1))
  m2=max(which(dd==1))
  if(sum(dd[m1:m2]>2)>0){ for(i in which(dd[m1:m2]>2)){qy=c(qy,(qy[m1+i-1]+qy[m1+i])/2)}
  }
  kn=sort(qy)
  
  s1=min(kn)
  s2=max(kn)
  yraw = y
  y = y[y>=s1 & y<=s2]
  n = length(y)
  if(is.null(lam)){
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
  }
  if(is.null(eps)){
    eps=ifelse(abs(skewness(y))>0.7,eps2,eps1)
  }else{eps1=eps2=eps}
  m=length(kn)+1
  bspl=bSpline(y,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  yp=0:4000/4000*(s2-s1)+s1
  s1=min(s1,min(yp))
  s2=max(s2,max(yp))
  bp=bSpline(yp,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  slopes=bSpline(kn,degree=2,derivs=1,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  D2=bSpline(kn,degree=2,derivs=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  d=(s2-s1)/m
  ###  make all basis functions integrate to one
  dp=yp[2]-yp[1]
  avec=1:m
  for(i in 1:m){avec[i]=sum(bp[,i])*dp}
  for(i in 1:m){
    bspl[,i]=bspl[,i]/avec[i]
    bp[,i]=bp[,i]/avec[i]
    slopes[,i]=slopes[,i]/avec[i]
    D2[,i]=D2[,i]/avec[i]
  }
  av1=avec
  avec=rep(1,m)
  hmat=matrix(nrow=m,ncol=m)
  cvec=1:m
  for(i in 1:m){
    for(j in i:m){
      pr=bp[,i]*bp[,j]
      hmat[i,j]=sum(pr)*dp
      hmat[j,i]=hmat[i,j]
    }
    cvec[i]=sum(bspl[,i])/n
  }
  b0=rep(1/m,m)
  
  wmat=matrix(0,nrow=m,ncol=m-1)
  for(i in 1:(m-1)){wmat[i,i]=-1;wmat[i+1,i]=1}
  
  D1=matrix(0,m-2,m-1)
  for(i in 1:(m-2)){D1[i,i]=-1;D1[i,i+1]=1}
  D=D1%*%D2*d^(5/2)
  ##  get unimodal
  ans1=umfit(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam,eps=eps)
  ans2=bmfit(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam)
  t1=as.numeric((ans1$crit-ans2$crit))
  #t2=as.numeric((ans1$crit-ans2$crit)/abs(ans1$crit))
  fhat2=round(bp%*%ans2$bhat,10)
  dfhat2=diff(fhat2)
  outtb=NULL
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    pvalue=2
  }else{
    md=max(which(dfhat2>0))
    if (!is.unsorted(fhat2[1:(md-1)]) ) {
      pvalue=2
    } else {
      cdf1=bp%*%ans1$bhat
      for(i in 2:4001){
        cdf1[i]=cdf1[i-1]+cdf1[i]
      }
      cdf1=cdf1-min(cdf1)
      cdf1=cdf1/cdf1[4001]
      outtb <- foreach(t=1:B,.combine = 'rbind') %dopar%{
        yb=sapply(1:n,function(o){u=runif(1);id=min(which(u<cdf1));alp=(cdf1[id]-u)/(cdf1[id]-cdf1[id-1]);alp*yp[id-1]+(1-alp)*yp[id]})
        modet(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam=lam,eps=eps)
      } 
      pvalue=sum(outtb>t1)/B
    }
  }
  ans=new.env()
  ans$yp=yp
  ans$fhat1=bp%*%ans1$bhat
  ans$fhat2=bp%*%ans2$bhat
  ans$statistic=t1
  ans$tb=outtb
  ans$pvalue=pvalue
  ans$lam=c(ans1$lam,ans2$lam)
  ans$kn=kn
  ans$crit=c(ans1$crit,ans2$crit)
  ans$kurtosis=kurtosis(y)
  ans$skewness=skewness(y)
  ans$truncated=(nraw!=n)
  ans
}
##############
umfit=function(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam=NULL,eps){	
  n <- length(y)
  m=length(kn)+1
  k1=max(min(which(kn>quantile(y,.2)))-1,3)
  k2=min(max(which(kn<quantile(y,.8)))+1,m-4)         
  amatl1=list()
  for(k in k1:k2){
    amat=matrix(0,nrow=m+1,ncol=m)
    amat[1:k,]=slopes[1:k,]
    amat[(k+1):(m-1),]=-slopes[(k+1):(m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    epsvec=c(rep(0,k1-2),rep(eps/n^(2/7)/diff(range(kn))^2,k-k1+1),0,0,rep(eps/n^(2/7)/diff(range(kn))^2,k2-k+1),rep(0,m-k2-3),0,0)
    amatl1[[k-k1+1]]=list(amat=amat, epsvec=epsvec)
  }
  lamt=lam
  zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*t(D)%*%D%*%b0)
  qmat=t(wmat)%*%(hmat+lamt*t(D)%*%D)%*%wmat 
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
bmfit=function(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam=NULL){
  n=length(y)
  m=dim(slopes)[2]
  m1=max(min(which(kn>quantile(y,.1)))-1,2)
  m2=min(max(which(kn<quantile(y,.9)))+1,m-2)
  if((m2-m1)<5){m1=max(2,m1-1);m2=min(m2+1,m-2)}
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
  amat2=list()
  for(i in 1:nr){
    amat=matrix(0,m+1,m)
    amat[1:(m-1),1:m]=slopes
    amat[(trips[i,1]+1):trips[i,2],]=-slopes[(trips[i,1]+1):trips[i,2],]
    amat[(trips[i,3]):(m-1),]=-slopes[(trips[i,3]):(m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    amat2[[i]]=amat
  }
  lamt=lam
  zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*t(D)%*%D%*%b0)
  qmat=t(wmat)%*%(hmat+lamt*t(D)%*%D)%*%wmat 
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
modet <- function(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam,eps){
  n=length(yb)
  m=dim(slopes)[2]
  capk=length(kn) 
  s1=min(kn)
  s2=max(kn)
  bb=bSpline(yb,degree=2,knots=kn[2:(capk-1)],Boundary.knots=c(s1,s2),intercept=TRUE)
  m=dim(bb)[2]
  for(i in 1:m){bb[,i]=bb[,i]/av1[i]}
  cb=1:m
  for(i in 1:m){
    cb[i]=sum(bb[,i])/n
  }
  fit1=umfit(yb,kn,hmat,cb,slopes,b0,wmat,D,bspl,lam,eps=eps)
  fit2=bmfit(yb,kn,hmat,cb,slopes,b0,wmat,D,bspl,lam)
  fhat2=round(bp%*%fit2$bhat,10)
  dfhat2=diff(fhat2)
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    return (-1e+3)
  }else{md=max(which(dfhat2>0))
  if (!is.unsorted(fhat2[1:(md-1)]) ){
    return (-1e+3)
  }else{return(as.numeric(fit1$crit-fit2$crit))}
  }
}








registerDoParallel(3)
n=1000
y=rnorm(n)
ans=bmodetest(y,lam1=n^(-5/9)*10^(5-kurtosis(y)/1))
hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
stopImplicitCluster()
