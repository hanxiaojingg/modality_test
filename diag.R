bmodetest <- function(x ,lower = NULL, upper = NULL,B=500,lam=NULL,eps=0.1,epsrate=2/7,cv=TRUE,parallel=FALSE,sample_from_true=FALSE){
  sd = sd(x)
  y = x/sd
  #sd=1
  #y=x
  n = length(y)
  y = sort(y)
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
  capk=round(n^(1/7)*12)
  qy=as.numeric(quantile(y,1:capk/(capk+1)))
  ####add more knots on both sides
  if(s1<qy[1]){
    k1=min(floor(log((qy[1]-s1)/3/(qy[2]-qy[1])+1, 3/2) ), round(capk/8))
    if(k1<2){
      qy = c(qy,s1)
      # if((qy[1]-s1)/(qy[2]-qy[1])>=2) qy = c(qy,qy[1]-(qy[1]-s1)/2)
    }else{
      gaps1 = (qy[2]-qy[1])*1.5^(1:(k1-1))
      qy = c(qy, qy[1] - cumsum(gaps1), s1)
    }
  }
  qy = sort(qy)
  capk = length(qy)
  if(s2>qy[capk]){
    k2=min(floor(log((s2-qy[capk])/3/(qy[capk]-qy[capk-1])+1, 3/2) ), round(capk/8))
    if(k2<2){
      qy=c(qy,s2)
      # if((s2-qy[capk])/(qy[capk]-qy[capk-1])>=2) qy = c(qy,(s2-qy[capk])/2+qy[capk])
    }else{
      # q2 = ((s2-qy[capk])/(qy[capk]-qy[capk-1]))^(1/(k2+1))
      gaps2 = (qy[capk]-qy[capk-1])*1.5^(1:(k2-1))
      qy = c(qy, cumsum(gaps2) + qy[capk], s2)
    }
  }
  qy = sort(qy)
  ####add more knots between the two modes
  d=min(diff(qy))
  dd=floor(diff(qy)/min(diff(qy)))
  m1=min(which(dd<=2))
  m2=max(which(dd<=2))
  if(sum(dd[m1:m2]>10)>0){ for(i in which(dd[m1:m2]>10)){qy=c(qy,(qy[m1+i-1]+qy[m1+i])/2)}
  }
  kn=sort(qy)
  yp=0:4000/4000*(s2-s1)+s1
  ## do no delete, might cause rounding issue
  s1=min(s1,min(yp))
  s2=max(s2,max(yp))
  m=length(kn)+1
  bspl=bSpline(y,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  bp=bSpline(yp,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  bk=bSpline(kn,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  slopes=bSpline(kn,degree=2,derivs=1,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  D2=bSpline(kn,degree=2,derivs=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  ###  make all basis functions integrate to one
  dp=yp[2]-yp[1]
  avec=1:m
  for(i in 1:m){avec[i]=sum(bp[,i])*dp}
  for(i in 1:m){
    bspl[,i]=bspl[,i]/avec[i]
    bp[,i]=bp[,i]/avec[i]
    bk[,i]=bk[,i]/avec[i]
    slopes[,i]=slopes[,i]/avec[i]
    D2[,i]=D2[,i]/avec[i]
  }
  av1=avec
  avec=rep(1,m)
  hmat=crossprod(bp) * dp
  b0=rep(1/m,m)
  cvec=colMeans(bspl)
  wmat=matrix(0,nrow=m,ncol=m-1)
  for(i in 1:(m-1)){wmat[i,i]=-1;wmat[i+1,i]=1}
  
  D1=matrix(0,m-2,m-1)
  for(i in 1:(m-2)){D1[i,i]=-1;D1[i,i+1]=1}
  D=D1%*%D2
  
  ## unimodal amat list
  amatl1=list()
  ## mode between t_k and t_{k+1} where k=1,...,m-2
  for(k in 1:2){
    amat=matrix(0,nrow=m+1,ncol=m)
    amat[1:k,]=slopes[1:k,]
    amat[(k+1):(k+2),]=-slopes[(k+1):(k+2),]
    amat[(k+3):(m-1),]=-slopes[(k+3):(m-1),]-eps/n^epsrate*bk[(k+3):(m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    #epsvec=c(rep(eps/n^(2/7)/diff(range(kn))^2,k-1),0,0,rep(eps/n^(2/7)/diff(range(kn))^2,m-k-2),0,0)
    #amatl1[[k]]=list(amatw=t(amat%*%wmat), epsbvec=epsvec-amat%*%b0)
    amatl1[[k]]=list(amatw=t(amat%*%wmat), bvec= -amat %*% b0)
  }
  for(k in 3:(m-4)){
    amat=matrix(0,nrow=m+1,ncol=m)
    amat[1:(k-2),]=slopes[1:(k-2),]-eps/n^epsrate*bk[1:(k-2),]
    amat[c(k-1,k),]=slopes[c(k-1,k),]
    amat[c(k+1,k+2),]=-slopes[c(k+1,k+2),]
    amat[(k+3):(m-1),]=-slopes[(k+3):(m-1),]-eps/n^epsrate*bk[(k+3):(m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    #epsvec=c(rep(eps/n^(2/7)/diff(range(kn))^2,k-1),0,0,rep(eps/n^(2/7)/diff(range(kn))^2,m-k-2),0,0)
    #amatl1[[k]]=list(amatw=t(amat%*%wmat), epsbvec=epsvec-amat%*%b0)
    amatl1[[k]]=list(amatw=t(amat%*%wmat), bvec= -amat %*% b0)
  }
  for(k in (m-3):(m-2)){
    amat=matrix(0,nrow=m+1,ncol=m)
    amat[1:(k-2),]=slopes[1:(k-2),]-eps/n^epsrate*bk[1:(k-2),]
    amat[c(k-1,k),]=slopes[c(k-1,k),]
    amat[c(k+1,m-1),]=-slopes[c(k+1,m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    #epsvec=c(rep(eps/n^(2/7)/diff(range(kn))^2,k-1),0,0,rep(eps/n^(2/7)/diff(range(kn))^2,m-k-2),0,0)
    #amatl1[[k]]=list(amatw=t(amat%*%wmat), epsbvec=epsvec-amat%*%b0)
    amatl1[[k]]=list(amatw=t(amat%*%wmat), bvec= -amat %*% b0)
  }
  ## mode at t_1 and t_{m-1}
  ## for amat_{m-1} and amat_m
  amatt1 = rbind(-slopes,c(rep(0,m-1),1))
  amatt1[3:(m-1),] = -slopes[3:(m-1),]-eps/n^epsrate*bk[3:(m-1),]
  amattm1 = rbind(slopes,c(1, rep(0,m-1)))
  amattm1[1:(m-3),]=slopes[1:(m-3),]-eps/n^epsrate*bk[1:(m-3),]
  amatl1[[m-1]]=list(amatw = t(amatt1%*%wmat), bvec = -amatt1 %*% b0 )
  amatl1[[m]]=list(amatw = t(amattm1%*%wmat), bvec = -amattm1 %*% b0)
  ## bimodal amat list
  ## triplets (i,j,k) where i<j<k 
  trips=matrix(0,nrow=choose(m,3),ncol=3)
  nr=0
  amatl2=list()
  ## modes and antimodes between t_i and t_{i+1}, t_j and t_{j+1}, t_k and t_{k+1}, k<=m-2
  for(i in 1:(m-6)){
    for(j in (i+2):(m-4)){
      for(k in (j+2):(m-2)){
        nr=nr+1
        trips[nr,]=c(i,j,k)
        amat=matrix(0,m+3,m)
        amat[1:(m-1),1:m]=slopes
        amat[(i+1):j,]=-slopes[(i+1):j,]
        amat[k:(m-1),]=-slopes[k:(m-1),]
        amat[m,1]=1
        amat[m+1,m]=1
        amat[m+2,]=bk[j,]
        amat[m+3,]=bk[j+1,]
        amatl2[[nr]]=list(amatw = t(amat %*% wmat), bvec = -amat %*% b0)
      }
    }
  }
  ## modes at the lower endpoint, i.e. i=1 and mode is at t_1
  for(j in 3:(m-4)){
    for(k in (j+2):(m-2)){
      nr=nr+1
      trips[nr,]=c(1,j,k)
      amat=matrix(0,m,m)
      amat[1:(m-1),1:m]=slopes
      amat[1:j,]=-slopes[1:j,]
      amat[(k+1):(m-1),]=-slopes[(k+1):(m-1),]
      amat[m,m]=1
      amatl2[[nr]]=list(amatw = t(amat %*% wmat), bvec = -amat %*% b0)
    }
  }
  ## modes at upper endpoint, i.e. k=m-2 and mode is at t_{m-1}
  for(i in 1:(m-6)){
    for(j in (i+2):(m-4)){
      nr=nr+1
      trips[nr,]=c(i,j,m-2)
      amat=matrix(0,m,m)
      amat[1:(m-1),1:m]=slopes
      amat[(i+1):j,]=-slopes[(i+1):j,]
      amat[m,1]=1
      amatl2[[nr]]=list(amatw = t(amat %*% wmat), bvec = -amat %*% b0)
    }
  }
  ## both modes at the endpoints
  for(j in 3:(m-4)){
    nr=nr+1
    trips[nr,]=c(1,j,m-2)
    amat=matrix(0,m-1,m)
    amat[1:(m-1),1:m]=slopes
    amat[1:j,]=-slopes[1:j,]
    amatl2[[nr]]=list(amatw = t(amat %*% wmat), bvec = -amat %*% b0)
  }
  trips = trips[1:nr,,drop=FALSE]
  
  
  ##  get unimodal
  DtD = crossprod(D)
  ans2=bmfit(y,kn,amatl2,hmat,cvec,slopes,b0,wmat,DtD,bspl,lam,cv=cv)
  ans1=umfit(y,kn,amatl1,hmat,cvec,slopes,b0,wmat,DtD,bspl,lam = ans2$lam,eps=eps, cv=FALSE)
  t1=as.numeric((ans1$crit-ans2$crit))
  #t2=as.numeric((ans1$crit-ans2$crit)/abs(ans1$crit))
  fhat2=round(bp%*%ans2$bhat,10)
  fhat1=round(bp%*%ans1$bhat,10)
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
      cdf1=cdf1/cdf1[4001]
      one_boot <- function(t) {
        #yb=sapply(1:n,function(o){u=runif(1);id=min(which(u<cdf1));alp=(cdf1[id]-u)/(cdf1[id]-cdf1[id-1]);alp*yp[id-1]+(1-alp)*yp[id]})
        u = runif(n)
        id = findInterval(u,cdf1) + 1
        id = pmax(2, pmin(id, length(cdf1)))
        alp = (cdf1[id]-u)/(cdf1[id]-cdf1[id-1])
        if(sample_from_true){
          yb <- rtrue(n) / sd
        }else{
          yb = alp*yp[id-1]+(1-alp)*yp[id] 
        }
        modet(yb,kn,amatl1,amatl2,hmat,slopes,b0,wmat,DtD,bspl,av1,bp,lam=ans2$lam,eps=eps)
      }
      if(parallel){
        outtb <- foreach(t=1:B,.combine='c') %dopar% {
          one_boot(t)
        }
      } else {
        outtb <- sapply(1:B, one_boot)
      }
      pvalue=sum(outtb>=t1)/B
    }
  }
  ### fhat2 correction
  fhat2 = bp%*%ans2$bhat
  #trapz <- function(x,y){sum(diff(x)*(head(y,-1)+tail(y,-1)))/2}
  #int_original <- trapz(yp,fhat2)
  #int_neg <- trapz(yp,pmin(fhat2,0))
  fhat2_corrected <- pmax(fhat2,0)
  int_corrected <- sum(diff(yp)*(head(fhat2_corrected,-1)+tail(fhat2_corrected,-1)))/2
  fhat2_corrected <- fhat2_corrected/int_corrected
  ans=new.env()
  ans$xp=yp*sd
  ans$yp=yp
  ans$fhat1=fhat1/sd
  ans$fhat2=fhat2_corrected/sd
  ans$statistic=t1
  ans$tb=outtb
  ans$pvalue=pvalue
  ans$lam=ans2$lam
  ans$kn=kn*sd
  ans$crit=c(ans1$crit,ans2$crit)
  ans
}


check_grid <- function(ans) {
  kn <- sort(ans$kn)
  xp <- ans$xp
  
  knot_gaps <- diff(kn)
  grid_step <- median(diff(xp))
  
  data.frame(
    support_width = diff(range(xp)),
    grid_step = grid_step,
    min_knot_gap = min(knot_gaps),
    median_knot_gap = median(knot_gaps),
    max_knot_gap = max(knot_gaps),
    min_gap_grid_points = min(knot_gaps) / grid_step,
    intervals_under_3_points =
      sum(knot_gaps / grid_step < 3),
    intervals_under_5_points =
      sum(knot_gaps / grid_step < 5)
  )
}

# After running bmodetest:
check_grid(ans)

n=200
z <- runif(n)
x = rnorm(n,0,9)
x[z < 0.4] <- rnorm(sum(z<0.4),0,1)
x[z > 0.6] <- rnorm(sum(z>0.6),2,1)
ans=bmodetest(x,B=1,parallel = FALSE)
check_grid(ans)



# True density: third component has SD = 9
ftrue <- function(x) {
  0.4 * dnorm(x, 0, 1) +
    0.4 * dnorm(x, 2, 1) +
    0.2 * dnorm(x, 0, 9)
}

plot(ans$xp, ans$fhat1,
     type = "l", col = "red", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "Three-normal mixture, n = 200",
     ylim = range(c(ans$fhat1, ans$fhat2,
                    ftrue(ans$xp))))

lines(ans$xp, ans$fhat2, col = "blue", lwd = 1)
lines(ans$xp, ftrue(ans$xp),
      col = "black", lwd = 1, lty = 2)
rug(ans$kn)
legend("topright",
       legend = c("Unimodal fit", "Bimodal fit",
                  "True density"),
       col = c("red", "blue", "black"),
       lty = c(1, 1, 2), lwd = 2)



set.seed(123)
ans_small <- bmodetest(
  x, B = 1, eps = 0.1, epsrate = 1/7,
  parallel = FALSE
)

set.seed(123)
ans_large <- bmodetest(
  x, B = 1, eps = 0.34, epsrate = 1/7,
  parallel = FALSE
)

plot(ans_small$xp, ans_small$fhat1,
     type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "Effect of epsilon on unimodal fit")

lines(ans_large$xp, ans_large$fhat1,
      col = "red", lwd = 1)

lines(ans_small$xp, ftrue(ans_small$xp),
      col = "black", lty = 2, lwd = 1)

legend("topright",
       legend = c("Epsilon = 0.1",
                  "Epsilon = 0.34",
                  "True density"),
       col = c("blue", "red", "black"),
       lty = c(1, 1, 2), lwd = 2)

data.frame(
  epsilon = c(0.1, 0.34),
  lambda = c(ans_small$lam, ans_large$lam),
  statistic = c(ans_small$statistic,
                ans_large$statistic)
)


ncores = max(1, detectCores() - 1)
registerDoParallel(ncores)
set.seed(123)
ans_small <- bmodetest(
  x, B = 200,
  eps = 0.1,
  epsrate = 1/7,
  parallel = TRUE
)

set.seed(123)
ans_large <- bmodetest(
  x, B = 200,
  eps = 0.34,
  epsrate = 1/7,
  parallel = TRUE
)

# Compare observed statistics and bootstrap distributions
summary_table <- data.frame(
  epsilon = c(0.1, 0.34),
  statistic = c(ans_small$statistic,
                ans_large$statistic),
  bootstrap_mean = c(mean(ans_small$tb),
                     mean(ans_large$tb)),
  bootstrap_q95 = c(
    unname(quantile(ans_small$tb, 0.95)),
    unname(quantile(ans_large$tb, 0.95))
  ),
  pvalue = c(ans_small$pvalue,
             ans_large$pvalue)
)

print(summary_table)
par(mfrow=c(1,2))
# Plot bootstrap statistics
hist(ans_small$tb[ans_small$tb>-1000],
     breaks = 25,
     freq = FALSE,
     main = "Bootstrap distribution, epsilon = 0.1",
     xlab = "Test statistic")
abline(v = ans_small$statistic, col = "red", lwd = 2)

hist(ans_large$tb[ans_large$tb>-1000],
     breaks = 25,
     freq = FALSE,
     main = "Bootstrap distribution, epsilon = 0.34",
     xlab = "Test statistic")
abline(v = ans_large$statistic, col = "red", lwd = 2)


rtrue <- function(n) {
  z <- runif(n)
  x <- rnorm(n, 0, 9)
  
  x[z < 0.4] <- rnorm(sum(z < 0.4), 0, 1)
  x[z > 0.6] <- rnorm(sum(z > 0.6), 2, 1)
  
  x
}
par(mfrow=c(1,1))
ncores = max(1, detectCores() - 1)
registerDoParallel(ncores)
#stopImplicitCluster()
datay=read.csv("/Users/sheronjing/Library/CloudStorage/Dropbox/Hanxiao/modality/revision/goofyrej.txt")[,1]
ans=bmodetest(datay,B=200,parallel = TRUE)
hist(datay,xlim=range(ans$xp),freq=FALSE,breaks=30)
lines(ans$xp,ans$fhat1,col=2)
lines(ans$xp,ans$fhat2,col=3)
rug(ans$kn)
ans$pvalue