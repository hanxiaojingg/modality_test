
n=200
nloop=10000
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
	lines(xp,.6*dnorm(xp*sy+my)*sy+.4*dnorm(xp*sy+my,d,1)*sy,col=3)
	atest=testunim(y,100,.16)
	pv[iloop]=atest$pval
	lines(atest$xp,atest$fhat1,lwd=2)
	lines(atest$xp,atest$fhat2,col=2,lwd=2)
	print(iloop+sum(pv[1:iloop]<.05)/iloop)
	if(pv[iloop]>.05){ysave=y}
	lamch[iloop]=atest$lam
}
sum(pv<.05)/nloop


