library(quadprog)
library(doParallel)
#############################################################
#' @param y A vector of observed contaminated samples
#' @param h the measurement error size 
#' @param lower lower bound of the support of X
#' @param upper upper bound of the support of X
#' @param B the number of replicates used in the test
#' @param lam the smoothing penalty parameter, if lam=NULL, it is chosen via cross-validation
#' @param eps customized penalty parameter for 'flat spot' slope, the actual value used is eps*n^{-5/9}
#' @param parallel default = TRUE, whether use parallel for sampling
##############################################################
deconmodetest <- function(y,h,lower = NULL, upper = NULL,B=500,lam=NULL, eps=0.001, parallel=TRUE){
  n = length(y)
  sd = sd(x)
  y = sort(y)/sd
  h = h/sd
  if(is.null(lower)){
    if(n > 5){
      s1 = min(y) - max(diff(y[1:6]))
    } else {
      s1 = min(y)
    }
  } else {
    s1 = lower
  }
  
  if(is.null(upper)){
    if(n > 5){
      s2 = max(y) + max(diff(y[(n-5):n]))
    } else {
      s2 = max(y)
    }
  } else {
    s2 = upper
  }
  ran=s2-s1
  l = max(round((20*n^(1/9))*h/ran),2)
  d = h/l
  s1=min(s1-(ceiling(ran/d)*d-ran)/2, min(y)+h)
  s2=max(s2+(ceiling(ran/d)*d-ran)/2, max(y)-h)
  tk = seq(s1,s2,d)
  m = length(tk)-3
  yp = seq(s1-h,s2+h,length.out = 4000); np = length(yp)
  delta  = matrix(0,nrow=np,ncol=m)
  for(i in 1:m){
    delta[yp>tk[i]&yp<=tk[i+1],i] = 2*(yp[yp>tk[i]&yp<=tk[i+1]]-tk[i])^2/d^2/3
    delta[yp>=tk[i+1]&yp<=tk[i+2],i] = 1-4*(yp[yp>=tk[i+1]&yp<=tk[i+2]]-tk[i+1]-d/2)^2/d^2/3
    delta[yp>tk[i+2]&yp<=tk[i+3],i] = 2*(yp[yp>tk[i+2]&yp<=tk[i+3]]-tk[i+3])^2/d^2/3
  }
  slopes=matrix(0,m+1,m)
  for(i in 1:m){slopes[i,i]=1;slopes[i+1,i]=-1}
  slopes = 4/3/d*slopes
  bk = matrix(0,m+1,m)
  for(i in 1:m)  bk[i,i] = bk[i+1,i] = 2/3
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
  
  ##unimodal amat list
  k1=max(min(which(tk>quantile(y,.1)))-1,4)
  k2=min(max(which(tk<quantile(y,.9)))+1,m-1)   
  amatl1=list()
  ## when the mode is between t_i,t_i+1
  for(i in k1:k2){
    amat=matrix(0,nrow=m+3,ncol=m)
    amat[1:(i-2),]=slopes[1:(i-2),]-eps/n^(2/9)*bk[1:(i-2),]
    amat[c(i-1,i),]=slopes[c(i-1,i),]
    amat[c(i+1,i+2),]=-slopes[c(i+1,i+2),]
    amat[(i+3):(m+1),]=-slopes[(i+3):(m+1),]-eps/n^(2/9)*bk[(i+3):(m+1),]
    amat[m+2,1]=1
    amat[m+3,m]=1
    amatl1[[i-k1+1]]=list(amatw=t(amat%*%wmat), bvec= -amat %*% b0)
  }  
  ##bimodal amat list
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

  amatl2=list()
  for(i in 1:nr){
    amat=matrix(0,m+5,m)
    amat[1:(m+1),1:m]=slopes
    amat[(trips[i,1]+1):trips[i,2],]=-slopes[(trips[i,1]+1):trips[i,2],]
    amat[(trips[i,3]):(m+1),]=-slopes[(trips[i,3]):(m+1),]
    amat[m+2,]=bk[trips[i,2],]
    amat[m+3,]=bk[trips[i,2]+1,]
    amat[m+4,1]=1
    amat[m+5,m]=1
    amatl2[[i]]=list(amatw = t(amat %*% wmat), bvec = -amat %*% b0)
  }
  

  ##  get unimodal
  DtD = crossprod(D)
  ans2=deconbmfit(y,tk,amatl2,hmat,cvec,b0,wmat,DtD,lam=NULL, d, h)
  ans1=deconumfit(y,tk,amatl1,hmat,cvec,b0,wmat,DtD,lam=ans2$lam, eps=eps)
  t1=as.numeric((ans1$crit-ans2$crit))  
  fhat2=round(delta%*%ans2$bhat,10)
  fhat1=round(delta%*%ans1$bhat,10)
  dfhat2=diff(fhat2)
  outtb=NULL
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    pvalue=2
  }else{
    md=max(which(dfhat2>0))
    if (!is.unsorted(fhat2[1:(md-1)]) ) {
      pvalue=2
    } else {
      cdf1=cumsum(fhat1)
      cdf1=cdf1-min(cdf1)
      cdf1=cdf1/cdf1[4000]
      one_boot <- function(t) {
        u = runif(n)
        id = findInterval(u,cdf1) + 1
        id = pmax(2, pmin(id, length(cdf1)))
        alp = (cdf1[id]-u)/(cdf1[id]-cdf1[id-1])
        yb = alp*yp[id-1]+(1-alp)*yp[id] + runif(n, -h, h)
        deconmodet(yb,tk,amatl1,amatl2, d,h,hmat,b0,wmat,DtD,delta,lam=ans2$lam,eps=eps)#yb,tk,d,h,amatl1,amatl2, hmat,b0,wmat,DtD,delta,lam=ans$lam,eps=eps)
      }
      if(parallel){
        outtb <- foreach(t=1:B,.combine='c') %dopar% {
          one_boot(t)
        }
      } else {
        outtb <- sapply(1:B, one_boot)
      }
      pvalue=sum(outtb>t1)/B
    }
  }
  ans=new.env()
  ans$yp=yp*sd
  ans$fhat1=delta%*%ans1$bhat/sd
  ans$fhat2=delta%*%ans2$bhat/sd
  ans$ghat1=gamma%*%ans1$bhat/sd
  ans$ghat2=gamma%*%ans2$bhat/sd
  ans$statistic=t1
  ans$tb=outtb
  ans$pvalue=pvalue
  ans$lam=ans2$lam
  ans$kn=tk*sd
  ans$crit=c(ans1$crit,ans2$crit)
  ans
}
##############
deconumfit=function(y,tk,amatl1,hmat,cvec,b0,wmat,DtD,lam=NULL,eps){	
  n <- length(y)
  lamt = lam*n^{-5/9}
  zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*DtD%*%b0)
  qmat=t(wmat)%*%(hmat+lamt*DtD)%*%wmat 
  crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,x$amatw,x$bvec);ans$value})
  amatw1 <- amatl1[[which.min(crit)]][[1]]
  bvec <- amatl1[[which.min(crit)]][[2]]
  ans1=solve.QP(qmat,zvec,amatw1,bvec)
  alphahat1=ans1$solution
  bhat1=wmat%*%alphahat1+b0
  cr1=t(bhat1)%*%(hmat+lamt*DtD)%*%bhat1-2*sum(cvec*bhat1)
  
  ans1=new.env()
  ans1$bhat=bhat1
  ans1$lam=lamt
  ans1$crit=cr1
  ans1
}
############################################################################
deconbmfit=function(y,tk,amatl2,hmat,cvec,b0,wmat,DtD,lam=NULL,d ,h){
  n <- length(y)
  m=length(tk)-3
  if(is.null(lam)){
    zvec=t(wmat)%*%(cvec-hmat%*%b0-0.1*n^(-5/9)*DtD%*%b0)
    qmat=t(wmat)%*%(hmat+0.1*n^(-5/9)*DtD)%*%wmat 
    crit<-lapply(amatl2, function(x){ans <- quadprog::solve.QP(qmat,zvec,x$amatw,x$bvec);ans$value})
    amatw2 <- amatl2[[which.min(crit)]][[1]]
    bvec <- amatl2[[which.min(crit)]][[2]]
    L = 10
    lam_set = c(0.0001,0.001,0.01,0.1,1,10,100)
    fold_assignments <- sample(rep(1:L, length.out=n))
    crit_cv = NULL
    for( lamt in lam_set){
      err = rep(0,L)
      for (l in 1:L) {
        n_train = sum(fold_assignments != l)
        n_val   = sum(fold_assignments == l)
        cvec_train = makec(y[fold_assignments != l], tk, d, h)
        cvec_val = makec(y[fold_assignments == l], tk, d, h)
        zvec=t(wmat)%*%(cvec_train-hmat%*%b0-lamt*n_train^(-5/9)*DtD%*%b0)
        qmat=t(wmat)%*%(hmat+lamt*n_train^(-5/9)*DtD)%*%wmat 
        ans2=solve.QP(qmat,zvec,amatw2,bvec)
        alphahat2=ans2$solution
        bhat2=wmat%*%alphahat2+b0
        err[l]=t(bhat2)%*%hmat%*%bhat2-2*sum(cvec_val*bhat2)
      }
      crit_cv = c(crit_cv,mean(err))
    }
    lamt=lam_set[which(crit_cv==min(crit_cv))]
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-5/9)*DtD%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-5/9)*DtD)%*%wmat 
    crit<-lapply(amatl2, function(x){ans <- quadprog::solve.QP(qmat,zvec,x$amatw,x$bvec);ans$value})
    amatw2 <- amatl2[[which.min(crit)]][[1]]
    bvec <- amatl2[[which.min(crit)]][[2]]
    ans2=solve.QP(qmat,zvec,amatw2,bvec)
    alphahat2=ans2$solution
    bhat2=wmat%*%alphahat2+b0
    cr2=t(bhat2)%*%(hmat+lamt*n^(-5/9)*DtD)%*%bhat2-2*sum(cvec*bhat2)
  }else{
    lamt=lam*n^{-5/9}
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*DtD%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*DtD)%*%wmat 
    crit<-lapply(amatl2, function(x){ans <- quadprog::solve.QP(qmat,zvec,x$amatw,x$bvec);ans$value})
    amatw2 <- amatl2[[which.min(crit)]][[1]]
    bvec <- amatl2[[which.min(crit)]][[2]]
    ans2=solve.QP(qmat,zvec,amatw2,bvec)
    alphahat2=ans2$solution
    bhat2=wmat%*%alphahat2+b0
    cr2=t(bhat2)%*%(hmat+lamt*DtD)%*%bhat2-2*sum(cvec*bhat2)    
  }

  ans2=new.env()
  ans2$bhat=bhat2
  ans2$lam=lamt
  ans2$crit=cr2
  ans2
}
##################################
deconmodet <- function(yb,tk,amatl1,amatl2,d,h,hmat,b0,wmat,DtD,delta,lam,eps){
  n=length(yb)
  m=length(tk)-3
  capk=length(tk) 
  s1=min(tk)
  s2=max(tk)
  cb=makec(yb,tk,d,h)
  fit1=deconumfit(yb,tk,amatl1, hmat,cb,b0,wmat,DtD,lam,eps=eps)
  fit2=deconbmfit(yb,tk,amatl2, hmat,cb,b0,wmat,DtD,lam,d,h)
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
############################
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



library(quadprog)
library(doParallel)
ncores <- detectCores()-1
registerDoParallel(ncores)
stopImplicitCluster()

par(mfrow=c(1,2))
n=600
h=2
mu2=4
p=0.5
uu=runif(n);x=rnorm(n,0,1);x[uu<p]=rnorm(sum(uu<p),mu2,1)
y=x+runif(n,-h,h)
t1=Sys.time()
ans=deconmodetest(y,h, B=500,eps = 0.001)
t2=Sys.time()
t2-t1
ans$pvalue
hist(y,xlim = range(ans$yp), breaks=40,freq=FALSE)
lines(ans$yp,ans$ghat1,col=2)
lines(ans$yp,ans$ghat2,col=3)
lines(ans$yp,((1-p)*pnorm(ans$yp+h,0,1)+p*pnorm(ans$yp+h,mu2,1))/2/h-((1-p)*pnorm(ans$yp-h,0,1)+p*pnorm(ans$yp-h,mu2,1))/2/h)
hist(x,xlim = range(ans$yp),breaks=40,freq=FALSE)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
lines(ans$yp,(1-p)*dnorm(ans$yp,0,1)+p*dnorm(ans$yp,mu2,1))
rug(ans$kn)
ans$lam


par(mfrow=c(1,2))
n=600
h=2
mu2=3.5
p=0.5
uu=runif(n);x=rnorm(n,0,1);x[uu<p]=rnorm(sum(uu<p),mu2,1)
y=x+rnorm(n,0,h/sqrt(3))
t1=Sys.time()
ans=deconmodetest(y,h, B=500,eps = 0.001)
t2=Sys.time()
t2-t1
ans$pvalue
hist(y,breaks=40,freq=FALSE)
lines(ans$yp,ans$ghat1,col=2)
lines(ans$yp,ans$ghat2,col=3)
lines(ans$yp,((1-p)*pnorm(ans$yp+h,0,1)+p*pnorm(ans$yp+h,mu2,1))/2/h-((1-p)*pnorm(ans$yp-h,0,1)+p*pnorm(ans$yp-h,mu2,1))/2/h)
hist(x,breaks=40,freq=FALSE)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
lines(ans$yp,(1-p)*dnorm(ans$yp,0,1)+p*dnorm(ans$yp,mu2,1))
rug(ans$kn)
ans$lam




n=600
h=0.5 # 500 out of 500
pv1=rep(0,100)
t1=Sys.time()
for(rep in 1:100){
  uu=runif(n);x=rnorm(n,0,1);x[uu<p]=rnorm(sum(uu<p),mu2,1)
  y=x+runif(n,-h,h)
  ans=deconmodetest(y,h,v=0.01, B=1000)
  pv1[rep]=ans$pvalue
}
t2=Sys.time()
t2-t1
sum(pv1<=0.05)


n=600
h=1 # 1000 out of 1000
pv1=rep(0,200)
t1=Sys.time()
for(rep in 1:100){
  uu=runif(n);x=rnorm(n,0,1);x[uu<p]=rnorm(sum(uu<p),mu2,1)
  y=x+runif(n,-h,h)
  ans=deconmodetest(y,h,v=0.01, B=1000)
  pv1[rep]=ans$pvalue
}
t2=Sys.time()
t2-t1
sum(pv1<=0.05)

n=600
h=1.5 # v= 0.01 89 out of 100 v=0.1 494 out of 500 246 out of 250 247 out of 250
pv2=rep(0,300)
t1=Sys.time()
for(rep in 1:300){
  uu=runif(n);x=rnorm(n,0,1);x[uu<p]=rnorm(sum(uu<p),mu2,1)
  y=x+runif(n,-h,h)
  ans=deconmodetest(y,h,v=0.1, B=1000)
  pv2[rep]=ans$pvalue
} 
t2=Sys.time()
t2-t1
sum(pv2<=0.05)

n=600
h=2 # 365 out of 600, 0.1 193 out of 300; 227 out of 300 ;130 out of 200 179 out of 250 166, 171 out of 250
pv3=matrix(0,250,2)
t1=Sys.time()
for(rep in 90:250){
  uu=runif(n);x=rnorm(n,0,1);x[uu<p]=rnorm(sum(uu<p),mu2,1)
  y=x+runif(n,-h,h)
  ans1=deconmodetest(y,h,v=0.1, B=1000)#v=0.1 0.75 v=1 0.14
 #ans2=deconmodetest(y,h,v=0.09, B=1000)#v=0.001 0.52 v=0.01,0.05,0.08 0.65
  pv3[rep,1]=ans1$pvalue
  #pv3[rep,2]=ans2$pvalue
} 
t2=Sys.time()
t2-t1
sum(pv3[,1]<=0.05)


n=600
h=2.5 # v=0.1 0.66 0.56 122 out of 200 126 out of 200 v=0.05 0.68 0.52 v=0.09 0.63 0.63 123 out of 200 184,193 out of 300 178 out of 300 177 out of 300 v=0.08 0.65
pv4=matrix(0,300,2)
t1=Sys.time()
for(rep in 1:300){
  uu=runif(n);x=rnorm(n,0,1);x[uu<p]=rnorm(sum(uu<p),mu2,1)
  y=x+runif(n,-h,h)
  ans1=deconmodetest(y,h,v=0.09, B=1000)
  #ans2=deconmodetest(y,h,v=0.08, B=800)
  pv4[rep,1]=ans1$pvalue
 # pv4[rep,2]=ans2$pvalue
} 
t2=Sys.time()
t2-t1
sum(pv4[,1]<=0.05)
#sum(pv4[,2]<=0.05)
stopCluster(cl)

cl <- makeCluster(no_cores - 2)
registerDoParallel(cl)
clusterEvalQ(cl, library(quadprog))  # Load the quadprog library in each worker
clusterExport(cl, list("deconmodetest", "deconumfit", "deconbmfit", "deconmodet", "makec"))


