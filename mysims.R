
n=200
nloop=100
pv=1:nloop
d=2
xp=-400:400/100
lamch=1:nloop
for(iloop in 1:nloop){
#	y=rnorm(n)
#	u=runif(n)
#	y[u<.4]=y[u<.4]+d
	y1=rnorm(n,d,1)
	y2=rnorm(n,0,9)
	y=rnorm(n)
	u=runif(n)
	y[u>.4]=y1[u>.4]
	y[u>.8]=y2[u>.8]
	my=mean(y);sy=sd(y)
	ys=(y-mean(y))/sd(y)
	hist(ys,freq=F,br=40,main='')
	lines(xp,.4*dnorm(xp*sy+my)*sy+.4*dnorm(xp*sy+my,d,1)*sy+.2*dnorm(xp*sy+my,0,9)*sy,col=3)
	atest=testunim(y,200,.16)
	pv[iloop]=atest$pval
	lines(atest$xp,atest$fhat1,lwd=2)
	lines(atest$xp,atest$fhat2,col=2,lwd=2)
	rug(atest$kn)
	print(iloop+sum(pv[1:iloop]<.05)/iloop)
	if(pv[iloop]>.05){ysave=y}
	lamch[iloop]=atest$lam
}
sum(pv<.05)/nloop




n=200
nloop=100
pv=1:nloop
d=2
xp=-400:400/100
lamch=1:nloop
for(iloop in 1:nloop){
  #	y=rnorm(n)
  #	u=runif(n)
  #	y[u<.4]=y[u<.4]+d
  y1=rnorm(n,d,1)
  y2=rnorm(n,0,9)
  y=rnorm(n)
  u=runif(n)
  y[u>.4]=y1[u>.4]
  y[u>.8]=y2[u>.8]
  my=mean(y);sy=sd(y)
  ys=(y-mean(y))/sd(y)
  hist(ys,freq=F,br=40,main='')
  lines(xp,.4*dnorm(xp*sy+my)*sy+.4*dnorm(xp*sy+my,d,1)*sy+.2*dnorm(xp*sy+my,0,9)*sy,col=3)
  atest=testunim(y,200,.16)
  pv[iloop]=atest$pval
  lines(atest$xp,atest$fhat1,lwd=2)
  lines(atest$xp,atest$fhat2,col=2,lwd=2)
  rug(atest$kn)
  print(iloop+sum(pv[1:iloop]<.05)/iloop)
  if(pv[iloop]>.05){ysave=y}
  lamch[iloop]=atest$lam
}
sum(pv<.05)/nloop


nsim=200
n=200
ncores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
t1 <- Sys.time()
sim_result <- foreach(
  r = seq_len(nsim),
  .combine = "rbind",
  .packages = c(
    "coneproj", "quadprog","splines2")
  # ),
  # .options.RNG = 1234 + i_n
) %dorng% {
  # Generate x once for this replicate
  y1=rnorm(n,d,1)
  y2=rnorm(n,0,9)
  y=rnorm(n)
  u=runif(n)
  y[u>.4]=y1[u>.4]
  y[u>.8]=y2[u>.8]
  my=mean(y);sy=sd(y)
  ys=(y-mean(y))/sd(y)
  atest=testunim(y,200,.16)
  atest$pval
}
elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
power = sum(sim_result<0.05)/nsim
power
###### close
parallel::stopCluster(cl)
foreach::registerDoSEQ()
