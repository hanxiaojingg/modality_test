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
  D1=matrix(0,m-2,m-1)
  for(i in 1:(m-2)){D1[i,i]=-1;D1[i,i+1]=1}
  D=D1%*%D2*d^(5/2)
  ##  get unimodal
  ans1=umfit(y,kn,hmat,cvec,slopes,D,bspl,lam,eps=eps,od)
  ans2=bmfit(y,kn,hmat,cvec,slopes,D,bspl,lam)
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
        modet(yb,kn,hmat,slopes,D,bspl,av1,bp,lam=lam,eps=eps,od)
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
umfit=function(y,kn,hmat,cvec,slopes,D,bspl,lam=NULL,eps,od){	
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
    amat=rbind(rep(1,m),amat)
    epsvec=c(1,epsvec)
    amatl1[[k-k1+1]]=list(amat=amat, epsvec=epsvec)
  }
  lamt=lam
  #if(is.null(lamt)){lamt=n^(-5/9)*10^(5-kurtosis(y)/1)}
  qmat=hmat+lamt*t(D)%*%D
  crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,cvec,t(x[[1]]),x[[2]],meq=1);ans$value})
  amat1 <- amatl1[[which.min(crit)]][[1]]
  epsvec <- amatl1[[which.min(crit)]][[2]]
  ans1=solve.QP(qmat,cvec,t(amat1),epsvec,meq=1)
  bhat1=ans1$solution
  cr1=t(bhat1)%*%(hmat+lam*t(D)%*%D)%*%bhat1-2*sum(cvec*bhat1)
  
  ans1=new.env()
  ans1$bhat=bhat1
  ans1$lam=lam
  ans1$crit=cr1
  ans1
}
############################################################################
bmfit=function(y,kn,hmat,cvec,slopes,D,bspl,lam=NULL){
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
    amat=rbind(rep(1,m),amat)
    amat2[[i]]=amat
  }
  lamt=lam
  #if(is.null(lamt)){lamt=n^(-5/9)*10^(5-kurtosis(y)/1)}
  qmat=hmat+lamt*t(D)%*%D
  crit<-lapply(amat2, function(x){ans <- quadprog::solve.QP(qmat,cvec,t(x),c(1,rep(0,m+1)),meq=1);ans$value})
  amat2 <- amat2[[which.min(crit)]]
  ans2=solve.QP(qmat,cvec,t(amat2),c(1,rep(0,m+1)),meq=1)
  bhat2=ans2$solution
  cr2=t(bhat2)%*%(hmat+lam*t(D)%*%D)%*%bhat2-2*sum(cvec*bhat2)
  ans2=new.env()
  ans2$bhat=bhat2
  ans2$lam=lam
  ans2$crit=cr2
  ans2
}
##################################
modet <- function(yb,kn,hmat,slopes,D,bspl,av1,bp,lam,eps,od){
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
  fit1=umfit(yb,kn,hmat,cb,slopes,D,bspl,lam,eps=eps,od)
  fit2=bmfit(yb,kn,hmat,cb,slopes,D,bspl,lam)
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

set.seed(1234)
n=400
y=rnorm(n,0,1)
ans2=bmodetest(y,B=1)
ans1$fhat1
ans1$statistic
ans2$statistic
ans2$fhat1-ans1$fhat1[2500:3000]
plot(ans1$yp,ans1$fhat2,type="l")
lines(ans2$yp,ans2$fhat2,col=2)
xgrid=seq(0,30,0.1)
plot(xgrid,0.6*dchisq(xgrid,5)+0.4*dchisq(xgrid,24),type="l")




Qmat <- function(m,l){
  mmat=matrix(0,m+2*l-1,m)
  for(j in 1:m){
    mmat[j:(j+2*l-1),j]=rep(1,2*l)
  }
  wmat=matrix(0,m,m-1)
  for(j in 1:(m-1)){
    wmat[j,j]=-1
    wmat[j+1,j]=1
  }
  t(wmat)%*%t(mmat)%*%mmat%*%wmat
}
Qmat(13,3)


par(bty='n')
par(mar=c(0,0,0,0))
par(mfrow=c(1,1))

plot(c(0,0.7),c(-.1,1),pch='',xaxt='n',yaxt='n')
lines(c(-1,2),c(.1,.1))

m=4;ell=2
mtil=m+2*ell-1
d=.07

for(i in -1:(m+1)){
  lines(c(.2+i*d,.2+i*d),c(.09,.11))
}
text(.2+1*d,.05,expression(t[1]))
text(.2+2*d,.05,expression(t[2]))
text(.2+3*d,.05,expression(t[3]))

lines(c(-1,.2+2*d),c(.1,.1),col='darkred')
lines(c(.2+2*d,.2+2*d),c(.1,.8),lty=3,col='darkred')
lines(c(.2+2*d,.2+3*d),c(.8,.8),lwd=3,col='darkred')
lines(c(.2+3*d,.2+3*d),c(.1,.8),lty=3,col='darkred')
lines(c(.2+3*d,2),c(.1,.1),col='darkred')
lines(c(.2+0*d,.2+1*d),c(.1,.3),lwd=2,col='darkred')
lines(c(.2+1*d,.2+4*d),c(.3,.3),lwd=2,col='darkred')
lines(c(.2+4*d,.2+5*d),c(.3,.1),lwd=2,col='darkred')


lines(c(-1,.2+1*d),c(.1,.1),col='steelblue')
lines(c(.2+1*d,.2+1*d),c(.1,.8),lty=3,col='steelblue')
lines(c(.2+1*d,.2+2*d),c(.8,.8),lwd=3,col='steelblue')
lines(c(.2+2*d,.2+2*d),c(.1,.8),lty=3,col='steelblue')
lines(c(.2+2*d,2),c(.1,.1),col='steelblue')
lines(c(.2-1*d,.2+0*d),c(.1,.3),lwd=2,col='steelblue')
lines(c(.2+0*d,.2+3*d),c(.3,.3),lwd=2,col='steelblue')
lines(c(.2+3*d,.2+4*d),c(.3,.1),lwd=2,col='steelblue')

text(.314,.87,expression(f[1]),col='steelblue')
text(.380,.87,expression(f[2]),col='darkred')
text(.15,.26,expression(g[1]),col='steelblue')
text(.52,.26,expression(g[2]),col='darkred')




quaddecon <- function(y,h,support=NULL,d=NULL,plot=TRUE,lam=NULL,mode0=TRUE){
  n = length(y)
  if(is.null(support)){
    support=c(0,0)
    dis=quantile(y,0.99)-quantile(y,0.01)
    support[1]=quantile(y,0.01)-dis/3
    support[2]=quantile(y,0.99)+dis/3
  }
  if(!is.null(d)) {
    l = h/d
    if(d>h|round(h/d)!=h/d) d=NULL
  }
  if(is.null(d)){
    l = round((22*n^(1/9))*h/diff(support))
    d = h/l 
  }
  if(l<2) {
    l = 2
    d = h/l
  }
  if(mode0){
    support[1] = floor(support[1]/d)*d
    support[2] = ceiling(support[2]/d)*d
  }else{
    ran=diff(support)
    support[1]=support[1]-(ceiling(ran/d)*d-ran)/2
    support[2]=support[2]+(ceiling(ran/d)*d-ran)/2
  }
  tk = seq(support[1],support[2],d)
  m = length(tk)-3
  yp = seq(support[1]-h,support[2]+h,length.out = 2000); np = length(yp)
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
  if(is.null(lam)){
    gammay = matrix(0,nrow=n,ncol=m)
    for(j in 1:m){
      gammay[y>=tk[j]-h&y<tk[j+1]-h,j] = (y[y>=tk[j]-h&y<tk[j+1]-h]+h-tk[j])^3/9/h/d^2
      gammay[y>=tk[j+1]-h&y<tk[j+2]-h,j] = -2*(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1]-d/2)^3/9/h/d^2+(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1])/2/h+d/h/12
      gammay[y>=tk[j+2]-h&y<tk[j+3]-h,j] = (y[y>=tk[j+2]-h&y<tk[j+3]-h]+h-tk[j+3])^3/9/h/d^2+2*d/3/h
      gammay[y>=tk[j+3]-h&y<=tk[j]+h,j] = 2*d/3/h
      gammay[y>tk[j]+h&y<=tk[j+1]+h,j] = -(y[y>tk[j]+h&y<=tk[j+1]+h]-h-tk[j])^3/9/h/d^2+2*d/3/h
      gammay[y>tk[j+1]+h&y<=tk[j+2]+h,j] = 2*(y[y>tk[j+1]+h&y<=tk[j+2]+h]-h-tk[j+1]-d/2)^3/9/h/d^2+(-y[y>tk[j+1]+h&y<=tk[j+2]+h]+tk[j+2]+h)/2/h+d/12/h
      gammay[y>tk[j+2]+h&y<=tk[j+3]+h,j] = (-y[y>tk[j+2]+h&y<=tk[j+3]+h]+tk[j+3]+h)^3/9/h/d^2
    }
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
  if(mode0){
    amat = matrix(0,m,m)
    tm = which(tk==0)
    amat[1,1] = 1
    for(i in 2:(tm-2)){
      amat[i,i-1] = -1
      amat[i,i] = 1
    }
    for(i in (tm-1):(m-1)){
      amat[i,i+1] = -1
      amat[i,i] = 1
      
    }
    amat[m,m] = 1
    amat=rbind(c(rep(0,tm-3),1,-1,rep(0,m+1-tm)),amat)
  }else{
    amat=list()
    for(i in 1:m){
      tmpm=matrix(0,m+1,m)
      tmpm[1,1]=1
      tmpm[m+1,m]=1
      for(j in 2:m){
        if(j<(i+1)){
          tmpm[j,j-1]=-1
          tmpm[j,j]=1
        }
        else{
          tmpm[j,j-1]=1
          tmpm[j,j]=-1
        }
      }
      amat[[i]]=tmpm
    }
    c0vec=t(wmat)%*%(cvec-hmat%*%b0-30*n^(-5/9)*t(D)%*%D%*%b0)
    q0mat=t(wmat)%*%(hmat+30*n^(-5/9)*t(D)%*%D)%*%wmat
    crit=lapply(amat,function(x){ans <- quadprog::solve.QP(q0mat,c0vec,t(x%*%wmat),-x%*%b0);ans$value})
    ind <- which.min(crit)
    amat <- amat[[ind]]
  }
  c0vec=t(wmat)%*%(cvec-hmat%*%b0)
  q0mat=t(wmat)%*%(hmat)%*%wmat
  nv <- ceiling(n/10)
  fold <- rep(1:10,nv)[1:n]
  ##if(nv<n/10) {fold <- c(fold,1:(n-10*nv))}
  nf <- 1:10
  for(i in 1:10){nf[i] <- sum(fold==i)}
  fold <- sample(fold,n)
  risk=NULL
  id=NULL
  if(is.null(lam)){
    penvals <- 2^(1:20)/2^20*n^(-5/9)*100
    # seq(1,diff(range(y)),length.out = 20)* n^(-11/18)
    if(mode0==FALSE){
      risk <- sapply(penvals,function(penval){
        c0vec=t(wmat)%*%(cvec-hmat%*%b0-penval*t(D)%*%D%*%b0)
        q0mat=t(wmat)%*%(hmat+penval*t(D)%*%D)%*%wmat
        ans <- quadprog::solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0)
        alphahat=ans$solution
        bhat=wmat%*%alphahat+b0
        tmpr=sum(sapply(1:10, function(j){xcv <- y[fold!=j];nm <- n-nf[j];cvecf <- makec(xcv,tk,d,h);
        solf <- quadprog::solve.QP(q0mat, t(wmat)%*%(cvecf-hmat%*%b0-penval*t(D)%*%D%*%b0), t(amat%*%wmat),-amat%*%b0)
        alphahatf <- solf$solution;bhatf <- wmat%*%alphahatf + b0;(gammay%*%bhatf)[fold==j]}))
        t(bhat)%*%hmat%*%bhat-2/n*tmpr
      })
    }else{  risk <- sapply(penvals,function(penval){
      c0vec=t(wmat)%*%(cvec-hmat%*%b0-penval*t(D)%*%D%*%b0)
      q0mat=t(wmat)%*%(hmat+penval*t(D)%*%D)%*%wmat
      ans <- quadprog::solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0,meq=1)
      alphahat=ans$solution
      bhat=wmat%*%alphahat+b0
      tmpr=sum(sapply(1:10, function(j){xcv <- y[fold!=j];nm <- n-nf[j];cvecf <- makec(xcv,tk,d,h);
      solf <- quadprog::solve.QP(q0mat, t(wmat)%*%(cvecf-hmat%*%b0-penval*t(D)%*%D%*%b0), t(amat%*%wmat),-amat%*%b0,meq=1)
      alphahatf <- solf$solution;bhatf <- wmat%*%alphahatf + b0;(gammay%*%bhatf)[fold==j]}))
      t(bhat)%*%hmat%*%bhat-2/n*tmpr
    })}
    rm=1:20*0
    if(risk[1]<risk[2]){rm[1]=1}
    for(i in 2:19){
      if(risk[i-1]>=risk[i]&risk[i+1]>=risk[i]){rm[i]=1}
    }
    if(risk[20]<risk[19]){rm[20]=1}
    #id <- max(which(rm==1))
    riskrange = diff(range(risk))
    id=1
    for(i in sort(which(rm==1),decreasing = TRUE)){
      if( (risk[i]-min(risk)) < riskrange/2) {id=i;break}
    }
    lam <- penvals[id]
    #risk <- risk[id]
  } 
  c0vec=t(wmat)%*%(cvec-hmat%*%b0-lam*t(D)%*%D%*%b0)
  q0mat=t(wmat)%*%(hmat+lam*t(D)%*%D)%*%wmat  
  if(mode0==TRUE){ans=solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0,meq=1)}else{
    ans=solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0)}
  crit=ans$value
  alphahat=ans$solution
  bhat=wmat%*%alphahat+b0
  #risk <- t(bhat)%*%hmat%*%bhat-2/n*sum(t(gamma)%*%bhatf)
  #alphatil=solve(q0mat)%*%c0vec
  #btil=wmat%*%alphatil+b0
  
  
  if(plot){
    par(mfrow=c(1,2))
    hist(y,freq=F,br=40,main="y density")
    lines(yp,gamma%*%bhat,col=4,lwd=3)
    plot(support,c(0,max(delta%*%bhat)*1.1),type="n",xlab="x",ylab = "Density",main = "x density")
    lines(yp,delta%*%bhat,col=4,lwd=3)
    rug(tk)
  }
  ans=new.env()
  ans$yp = yp
  ans$knots = tk
  ans$fhat = delta%*%bhat
  ans$ghat = gamma%*%bhat
  ans$bhat=bhat
  ans$support = support
  ans$d = d
  ans$amat=amat
  ans$mode=ifelse(mode0,0,yp[which.max(delta%*%bhat)])
  ans$lam=lam
  ans$risk=risk
  ans$id=id
  ans$crit=crit
  ans
  
}
n=800
x=rgamma(n,30,1)-29
y=x+runif(n,-10,10)
ans0=quaddecon(y,h=10,lam=0)
ans1=quaddecon(y,h=10,lam=0.001)
ans2=quaddecon(y,h=10,lam=0.01)
ans3=quaddecon(y,h=10,lam=0.1)
ans4=quaddecon(y,h=10,lam=1)

hist(x,breaks=20,freq=FALSE,xlab="",ylab="",ylim=c(0,0.095),main="unobserved (x)")
lines(ans0$yp,ans0$fhat,lwd=2.5,lty=2,col="darkred")
lines(ans1$yp,ans1$fhat,lwd=2.5,lty=6,col="steelblue")
lines(ans2$yp,ans2$fhat,lwd=2.5,lty=4,col="darkgreen")
lines(ans3$yp,ans3$fhat,lwd=2.5,lty=5,col=5)
lines(ans$yp,dgamma(ans$yp+29,30,1),lty=1,lwd=1.5)
legend("topright",legend = c(expression(paste(lambda, " = ", 0)),
                             expression(paste(lambda, " = ", .001)),
                             expression(paste(lambda, " = ", .01)),
                             expression(paste(lambda, " = ", .1)),'true'),cex=1,lwd=2.5,lty=c(2,6,4,5,1),col=c('darkred','steelblue','darkgreen',5,1))


hist(y,breaks=20,freq = FALSE, xlab="",ylab="",ylim=c(0,0.06),main="observed (y)")
lines(ans0$yp,ans0$ghat,lwd=2.5,lty=2,col="darkred")
lines(ans1$yp,ans1$ghat,lwd=2.5,lty=6,col="steelblue")
lines(ans2$yp,ans2$ghat,lwd=2.5,lty=4,col="darkgreen")
lines(ans3$yp,ans3$ghat,lwd=2.5,lty=5,col=5)
lines(ans$yp,(pgamma(ans$yp+29+10,30,1)-pgamma(ans$yp+29-10,30,1))/2/10,lty=1,lwd=1.5,col="steelblue")
legend("topright",legend = c(expression(paste(lambda, " = ", 0)),
                             expression(paste(lambda, " = ", .001)),
                             expression(paste(lambda, " = ", .01)),
                             expression(paste(lambda, " = ", .1)),'true'),cex=1,lwd=2.5,lty=c(2,6,4,5,1),col=c('darkred','steelblue','darkgreen',5,1))



n=800
x=rgamma(n,30,1)-29
y=x+runif(n,-10,10)
ans0=quaddecon(y,h=10,lam=0)
ans1=quaddecon(y,h=10,lam=0.001)
ans2=quaddecon(y,h=10,lam=0.01)
ans3=quaddecon(y,h=10,lam=0.1)
ans4=quaddecon(y,h=10,lam=1)

hist(x,breaks=20,freq=FALSE,xlab="",ylab="",ylim=c(0,0.095),main="unobserved (x)")
lines(ans0$yp,ans0$fhat,lwd=2.5,lty=2,col="darkred")
lines(ans1$yp,ans1$fhat,lwd=2.5,lty=6,col="steelblue")
lines(ans2$yp,ans2$fhat,lwd=2.5,lty=4,col="darkgreen")
lines(ans3$yp,ans3$fhat,lwd=2.5,lty=5,col=5)
lines(ans$yp,dgamma(ans$yp+29,30,1),lty=1,lwd=1.5)
legend("topright",legend = c(expression(paste(lambda, " = ", 0)),
                             expression(paste(lambda, " = ", .001)),
                             expression(paste(lambda, " = ", .01)),
                             expression(paste(lambda, " = ", .1)),'true'),cex=1,lwd=2.5,lty=c(2,6,4,5,1),col=c('darkred','steelblue','darkgreen',5,1))


hist(y,breaks=20,freq = FALSE, xlab="",ylab="",ylim=c(0,0.06),main="observed (y)")
lines(ans0$yp,ans0$ghat,lwd=2.5,lty=2,col="darkred")
lines(ans1$yp,ans1$ghat,lwd=2.5,lty=6,col="steelblue")
lines(ans2$yp,ans2$ghat,lwd=2.5,lty=4,col="darkgreen")
lines(ans3$yp,ans3$ghat,lwd=2.5,lty=5,col=5)
lines(ans$yp,(pgamma(ans$yp+29+10,30,1)-pgamma(ans$yp+29-10,30,1))/2/10,lty=1,lwd=1.5,col="steelblue")
legend("topright",legend = c(expression(paste(lambda, " = ", 0)),
                             expression(paste(lambda, " = ", .001)),
                             expression(paste(lambda, " = ", .01)),
                             expression(paste(lambda, " = ", .1)),'true'),cex=1,lwd=2.5,lty=c(2,6,4,5,1),col=c('darkred','steelblue','darkgreen',5,1))



n=1200
x=rnorm(n,3,1)
u=runif(n)
x[u<0.6]=rnorm(sum(u<0.6),-1,1)
y=x+runif(n,-8,8)
ans=quaddeconbm(y,h=8,lam=NULL)

hist(x,breaks=20,freq=FALSE,xlab="",ylab="",ylim=c(0,0.26),xlim=c(-5,7),main="unobserved (x)")
lines(ans$yp,ans$fhat,lwd=2.5,col="darkred")
lines(ans$yp,0.6*dnorm(ans$yp,-1,1)+0.4*dnorm(ans$yp,3,1),lty=2,lwd=2.5,col="steelblue")
f1=ker_amise(y,8)
f1$fhat[f1$fhat<0]=0
lines(f1$xt,f1$fhat,col='darkgreen',lwd=2.5,lty=4)
legend("topright",c('spline','kernel','true'),cex=1,lwd=2.5,lty=c(1,4,2),col=c('darkred','darkgreen','steelblue'))


hist(y,breaks=20,freq = FALSE, xlab="",ylab="",ylim=c(0,0.08),main="observed (y)")
lines(ans$yp,ans$ghat,lwd=2.5,col="darkred")
lines(ans$yp,(pnorm(ans$yp+8,-1,1)*0.6-pnorm(ans$yp-8,-1,1)*0.6+0.4*pnorm(ans$yp+8,3,1)-0.4*pnorm(ans$yp-8,3,1))/2/8,lty=2,lwd=2.5,col="steelblue")
legend("topright",c('spline','true'),cex=1,lwd=2.5,lty=c(1,2),col=c('darkred','steelblue'))


{
uu=runif(500);y=rnorm(500,0,1);y[uu<0.4]=rnorm(sum(uu<0.4),3,1);
nraw=length(y)
q1=quantile(y,.01)
q2=quantile(y,.99)
rng=q2-q1
s1=as.numeric(q1-.4*rng)
s2=as.numeric(q2+.4*rng)
capk=round(nraw^(1/7)*12)
qy=as.numeric(quantile(y,1:capk/(capk+1)))
qy1=qy
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
hist(y,breaks=30,freq=FALSE,xlim=c(s1,s2),main="")
lines(yp,0.6*dnorm(yp,0,1)+0.4*dnorm(yp,3,1),lwd=1.8,lty=2)
rug(kn,col="darkred",lwd=1.8)
rug(qy1,col="steelblue",lwd=1.8)
legend("topright",c('true density f','quantiles knots','additional knots'),lty=c(2,1,1),col=c('black','darkred','steelblue'),lwd=1.5)

}

s1=min(kn)
s2=max(kn)
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
  eps=ifelse(abs(skewness(y))>0.7,eps2,eps1)
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
D1=matrix(0,m-2,m-1)
for(i in 1:(m-2)){D1[i,i]=-1;D1[i,i+1]=1}
D=D1%*%D2*d^(5/2)
##  get unimodal
ans2=bmfit(y,kn,hmat,cvec,slopes,D,bspl,lam=0.9)
#t2=as.numeric((ans1$crit-ans2$crit)/abs(ans1$crit))
fhat2=round(bp%*%ans2$bhat,10)
hist(y,breaks=30,xlab="",ylab="",main="",freq=FALSE,ylim=c(0,0.35))
lines(yp,fhat2,lwd=2.5,col="darkred")
lines(yp,.7*dnorm(yp,0,1)+0.3*dnorm(yp,4,1),col='steelblue',lwd=2.5,lty=2)
legend("topright",c('penalized estimate','true'),cex=1,lwd=2.5,lty=c(1,2),col=c('darkred','steelblue'))


