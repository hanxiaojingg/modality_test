##slope<=eps/n^od/range^2
bmodetest <- function(y,lower = NULL, upper = NULL,B=1000,lam=NULL,v=1,eps1=5,eps2=2,eps=NULL,od=2/7,nc=12){
  nraw=length(y)
  s1=lower
  s2=upper
  q1=quantile(y,.01)
  q2=quantile(y,.99)
  rng=q2-q1
  if(is.null(lower)){s1=as.numeric(q1-.4*rng)}
  if(is.null(upper)){s2=as.numeric(q2+.4*rng)}
  capk=round(nraw^(1/7)*nc)
  qy=as.numeric(quantile(y,1:capk/(capk+1)))
  #dt=max(qy[2]-qy[1],qy[capk]-qy[capk-1])
  #k1=max(floor((qy[1]-s1)/dt/2),2)
  #k2=max(floor((s2-qy[capk])/dt/2),2)
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
    lam = v*lam
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
  ans1=umfit(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam,eps=eps,od)
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
             modet(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam=lam,eps=eps,od)
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
umfit=function(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam=NULL,eps,od){	
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
    #epsvec=c(0,0,rep(eps/n^(1/7)/diff(range(kn))^2,k-3),0,0,rep(eps/n^(1/7)/diff(range(kn))^2,m-4-k),0,0,0,0)
    epsvec=c(rep(0,k1-2),rep(eps/n^od/diff(range(kn))^2,k-k1+1),0,0,rep(eps/n^od/diff(range(kn))^2,k2-k+1),rep(0,m-k2-3),0,0)
    amatl1[[k-k1+1]]=list(amat=amat, epsvec=epsvec)
  }
  lamt=lam
  #if(is.null(lamt)){lamt=n^(-5/9)*10^(5-kurtosis(y)/1)}
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
  #if(is.null(lamt)){lamt=n^(-5/9)*10^(5-kurtosis(y)/1)}
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
modet <- function(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam,eps,od){
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
  fit1=umfit(yb,kn,hmat,cb,slopes,b0,wmat,D,bspl,lam,eps=eps,od)
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

install.packages("quadprog",repos="http://cran.us.r-project.org")
install.packages("doParallel",repos="http://cran.us.r-project.org")
install.packages("diptest",repos="http://cran.us.r-project.org")
install.packages("splines2",repos="http://cran.us.r-project.org")
install.packages("moments",repos="http://cran.us.r-project.org")
library("splines2")
library("quadprog")
library("doParallel")
library("diptest")
library("moments")
registerDoParallel(64)

#stopImplicitCluster()



set.seed(1234)
n=200
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.8]=rnorm(sum(uu>0.8),0,3);
x[uu<0.4]=rnorm(sum(uu<0.4),0,1)
x=matrix(x,1000,n)
out <- matrix(0,1000,6)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,],B=1000,eps1=5,eps2=2,od=2/7)$pvalue
}

set.seed(1234)
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.8]=rnorm(sum(uu>0.8),0,3);
x[uu<0.4]=rnorm(sum(uu<0.4),1,1)
x=matrix(x,1000,n)
for(i in 1:1000){
  out[i,2]=bmodetest(x[i,],B=1000,eps1=5,eps2=2,od=2/7)$pvalue
}

set.seed(1234)
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.8]=rnorm(sum(uu>0.8),0,3);
x[uu<0.4]=rnorm(sum(uu<0.4),1.2,1)
x=matrix(x,1000,n)
for(i in 1:1000){
  out[i,3]=bmodetest(x[i,],B=1000,eps1=5,eps2=2,od=2/7)$pvalue
}

set.seed(1234)
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.8]=rnorm(sum(uu>0.8),0,3);
x[uu<0.4]=rnorm(sum(uu<0.4),1.4,1)
x=matrix(x,1000,n)
for(i in 1:1000){
  out[i,4]=bmodetest(x[i,],B=1000,eps1=5,eps2=2,od=2/7)$pvalue
}

set.seed(1234)
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.8]=rnorm(sum(uu>0.8),0,3);
x[uu<0.4]=rnorm(sum(uu<0.4),1.6,1)
x=matrix(x,1000,n)
for(i in 1:1000){
  out[i,5]=bmodetest(x[i,],B=1000,eps1=5,eps2=2,od=2/7)$pvalue
}

set.seed(1234)
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.8]=rnorm(sum(uu>0.8),0,3);
x[uu<0.4]=rnorm(sum(uu<0.4),1.8,1)
x=matrix(x,1000,n)
for(i in 1:1000){
  out[i,6]=bmodetest(x[i,],B=1000,eps1=5,eps2=2,od=2/7)$pvalue
}
write(c(sum(out[,1]<=0.05),sum(out[,2]<=0.05),sum(out[,3]<=0.05),sum(out[,4]<=0.05),sum(out[,5]<=0.05),sum(out[,6]<=0.05)),file="/home/.colostate.edu/hanxiao/Mode/test.txt",sep=",")



set.seed(1234)
n=200
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.55]=rnorm(sum(uu>0.55),2.5,1);x[uu<0.25]=rnorm(sum(uu<0.25),-2.5,1)
x=matrix(x,1000,n)
out <- matrix(0,1000,2)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,])$pvalue
  out[i,2]=dip.test(x[i,])$p.value
  }
write(c(sum(out[,1]<=0.05),sum(out[,2]<=0.05)),file="/home/.colostate.edu/hanxiao/Mode/bptest1000.txt",sep=",")

yy=0.25*dnorm(xp,-3,1)+0.35*dnorm(xp,0,1)+0.4*dnorm(xp,2.5,1)
set.seed(1234)
n=200
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.6]=rnorm(sum(uu>0.6),2.5,1);x[uu<0.25]=rnorm(sum(uu<0.25),-3,1)
x=matrix(x,1000,n)
out <- matrix(0,1000,2)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,])$pvalue
  out[i,2]=dip.test(x[i,])$p.value
}
write(c(sum(out[,1]<=0.05),sum(out[,2]<=0.05)),file="/home/.colostate.edu/hanxiao/Mode/bptest1000.txt",sep=",")

set.seed(1234)
n=2000
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.6]=rnorm(sum(uu>0.6),0,4)
x=matrix(x,1000,n)
out <- matrix(0,1000,3)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=NULL,v=1)$pvalue
  out[i,2]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=NULL,v=2)$pvalue
}
write(c(sum(out[,1]<=0.05),sum(out[,2]<=0.05)),file="/home/.colostate.edu/hanxiao/Mode/bptest2000.txt",sep=",")



set.seed(1234)
n=5000
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.6]=rnorm(sum(uu>0.6),0,4)
x=matrix(x,1000,n)
out <- matrix(0,1000,3)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=NULL,v=1)$pvalue
  out[i,2]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=NULL,v=2)$pvalue
}
write(c(sum(out[,1]<=0.05),sum(out[,2]<=0.05)),file="/home/.colostate.edu/hanxiao/Mode/bptest5000.txt",sep=",")






set.seed(1234)
n=800
uu=runif(n*1000);x=rnorm(n*1000,-1.5,1);x[uu>0.6]=rnorm(sum(uu>0.6),1.5,1)
x=matrix(x,1000,n)
out <- matrix(0,1000,3)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=n^(-5/9)*10^(5-kurtosis(x[i,]))/10)$pvalue
  out[i,2]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=n^(-5/9)*10^(5-kurtosis(x[i,]))/50)$pvalue
  out[i,3]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=n^(-5/9)*10^(5-kurtosis(x[i,]))/100)$pvalue
}
write(c(sum(out[,1]<=0.05)/1000,sum(out[,2]<=0.05)/1000,sum(out[,3]<=500)),file="/home/.colostate.edu/hanxiao/Mode/bmodepower500.txt",sep=",")

set.seed(1234)
n=500
uu=runif(n*1000);x=rnorm(n*1000,-3,1);x[uu>0.7]=rnorm(sum(uu>0.7),3,1);x[uu<0.3]=rnorm(sum(uu<0.3),0,1)
x=matrix(x,1000,n)
out <- matrix(0,1000,2)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,],B=1000,lam1=NULL,lam2=NULL,v=2)$pvalue
  out[i,2]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=NULL,v=2)$pvalue
}
write(c(sum(out[,1]<=0.05)/1000,sum(out[,2]<=0.05)/1000),file="/home/.colostate.edu/hanxiao/Mode/btpower500.txt",sep=",")






set.seed(1234)
n=1000
uu=runif(n*1000);x=rnorm(n*1000,-1,1);x[uu>0.6]=rnorm(sum(uu>0.6),3,6)
x=matrix(x,1000,n)
out <- matrix(0,1000,3)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=n^(-5/9)*10^(5-kurtosis(x[i,]))/10)$pvalue
  out[i,2]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=n^(-5/9)*10^(5-kurtosis(x[i,]))/50)$pvalue
  out[i,3]=bmodetest(x[i,],B=1000,lam1=n^(-5/9)*10^(5-kurtosis(x[i,])),lam2=n^(-5/9)*10^(5-kurtosis(x[i,]))/100)$pvalue
}
write(c(sum(out[,1]<=0.05),sum(out[,2]<=0.05),sum(out[,3]<=0.05)),file="/home/.colostate.edu/hanxiao/Mode/pptest1000.txt",sep=",")





makec=function(y,tk,d,h){
  nc=length(y)
  nm=length(tk)-3
  cvec=1:nm
  if(h==0){
    for(j in 1:nm){
      gt=rep(0,nc)
      gt[y>tk[j]&y<=tk[j+1]] = 2*(y[y>tk[j]&y<=tk[j+1]]-tk[j])^2/d^2/3
      gt[y>=tk[j+1]&y<=tk[j+2]] = 1-4*(y[y>=tk[j+1]&y<=tk[j+2]]-tk[j+1]-d/2)^2/d^2/3
      gt[y>tk[j+2]&y<=tk[j+3]] = 2*(y[y>tk[j+2]&y<=tk[j+3]]-tk[j+3])^2/d^2/3
      cvec[j]=sum(gt)/nc
    }
  }else{
    for(j in 1:nm){
      gt=rep(0,nc)
      gt[y>tk[j]-h&y<tk[j+1]-h]=(y[y>tk[j]-h&y<tk[j+1]-h]+h-tk[j])^3/9/h/d^2
      gt[y>=tk[j+1]-h&y<tk[j+2]-h]=-2*(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1]-d/2)^3/9/h/d^2+(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1])/2/h+d/h/12
      gt[y>=tk[j+2]-h&y<tk[j+3]-h]=(y[y>=tk[j+2]-h&y<tk[j+3]-h]+h-tk[j+3])^3/9/d/h^2+2*d/3/h
      gt[y>=tk[j+3]-h&y<=tk[j]+h]=2*d/3/h
      gt[y>tk[j]+h&y<=tk[j+1]+h]= -(y[y>tk[j]+h&y<=tk[j+1]+h]-h-tk[j])^3/9/d/h^2+2*d/3/h
      gt[y>tk[j+1]+h&y<=tk[j+2]+h]=2*(y[y>tk[j+1]+h&y<=tk[j+2]+h]-h-tk[j+1]-d/2)^3/9/h/d^2+(-y[y>tk[j+1]+h&y<=tk[j+2]+h]+tk[j+2]+h)/2/h+d/12/h
      gt[y>tk[j+2]+h&y<tk[j+3]+h]=(-y[y>tk[j+2]+h&y<tk[j+3]+h]+tk[j+3]+h)^3/9/h/d^2
      cvec[j]=sum(gt)/nc
    }
  }
  cvec
}


##slope<=eps/n^od/range^2
deconmodetest <- function(y,h,lower = NULL, upper = NULL,d=NULL,B=1000,lam=NULL,v=1,eps1=5,eps2=2,eps=NULL,od=2/9){
  nraw=length(y)
  s1=lower
  s2=upper
  q1=quantile(y,.01)
  q2=quantile(y,.99)
  rng=q2-q1
  if(is.null(lower)){s1=as.numeric(q1-.4*rng)}
  if(is.null(upper)){s2=as.numeric(q2+.4*rng)}
  if(!is.null(d)) {
    l = h/d
  }
  if(is.null(d)){
    l = round((22*n^(1/9))*h/(s2-s1))
    d = h/l 
  }
  if(l<2) {
    l = 2
    d = h/l
  }
  ran=s2-s1
  s1=s1-(ceiling(ran/d)*d-ran)/2
  s2=s2+(ceiling(ran/d)*d-ran)/2
  tk = seq(s1,s2,d)
  m = length(tk)-3
  yp = seq(s1-h,s2+h,length.out = 2000); np = length(yp)
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
    lam = v*lam
  }
  if(is.null(eps)){
    eps=ifelse(abs(skewness(y))>0.7,eps2,eps1)
  }else{eps1=eps2=eps}
  delta  = matrix(0,nrow=np,ncol=m)
  for(i in 1:m){
    delta[yp>tk[i]&yp<=tk[i+1],i] = 2*(yp[yp>tk[i]&yp<=tk[i+1]]-tk[i])^2/d^2/3
    delta[yp>=tk[i+1]&yp<=tk[i+2],i] = 1-4*(yp[yp>=tk[i+1]&yp<=tk[i+2]]-tk[i+1]-d/2)^2/d^2/3
    delta[yp>tk[i+2]&yp<=tk[i+3],i] = 2*(yp[yp>tk[i+2]&yp<=tk[i+3]]-tk[i+3])^2/d^2/3
  }
  gamma = matrix(0,nrow=np,ncol=m)
  for(j in 1:m){
    gamma[yp>=tk[j]-h&yp<tk[j+1]-h,j] = (yp[yp>=tk[j]-h&yp<tk[j+1]-h]+h-tk[j])^3/9/h/d^2
    gamma[yp>=tk[j+1]-h&yp<tk[j+2]-h,j] = -2*(yp[yp>=tk[j+1]-h&yp<tk[j+2]-h]+h-tk[j+1]-d/2)^3/9/h/d^2+(yp[yp>=tk[j+1]-h&yp<tk[j+2]-h]+h-tk[j+1])/2/h+d/h/12
    gamma[yp>=tk[j+2]-h&yp<tk[j+3]-h,j] = (yp[yp>=tk[j+2]-h&yp<tk[j+3]-h]+h-tk[j+3])^3/9/h/d^2+2*d/3/h
    gamma[yp>=tk[j+3]-h&yp<=tk[j]+h,j] = 2*d/3/h
    gamma[yp>tk[j]+h&yp<=tk[j+1]+h,j] = -(yp[yp>tk[j]+h&yp<=tk[j+1]+h]-h-tk[j])^3/9/h/d^2+2*d/3/h
    gamma[yp>tk[j+1]+h&yp<=tk[j+2]+h,j] = 2*(yp[yp>tk[j+1]+h&yp<=tk[j+2]+h]-h-tk[j+1]-d/2)^3/9/h/d^2+(-yp[yp>tk[j+1]+h&yp<=tk[j+2]+h]+tk[j+2]+h)/2/h+d/12/h
    gamma[yp>tk[j+2]+h&yp<=tk[j+3]+h,j] = (-yp[yp>tk[j+2]+h&yp<=tk[j+3]+h]+tk[j+3]+h)^3/9/h/d^2
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
  fhat2=round(delta%*%ans2$bhat,10)
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



install.packages("quadprog",repos="http://cran.us.r-project.org")
install.packages("doParallel",repos="http://cran.us.r-project.org")
install.packages("moments",repos="http://cran.us.r-project.org")
library("quadprog")
library("doParallel")
library("moments")
registerDoParallel(64)

#stopImplicitCluster()



set.seed(1234)
n=800
h=1
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu<0.4]=rnorm(sum(uu<0.4),4,1);x=x+runif(n*1000,-h,h)
x=matrix(x,1000,n)
out <- matrix(0,1000,3)
for(i in 1:1000){
  out[i,1]=deconmodetest(x[i,],h,v=1)$pvalue
  out[i,2]=deconmodetest(x[i,],h,v=0.1)$pvalue
  out[i,3]=deconmodetest(x[i,],h,v=0.01)$pvalue
}

write(c(sum(out[,1]<=0.05),sum(out[,2]<=0.05),sum(out[,3]<=0.05)),file="/home/.colostate.edu/hanxiao/Mode/decon.txt",sep=",")

