### penden unimodal & bimodal

testunim=function(ys,nsim,eps){
	n=length(ys)
	nkn=round(8*n^(1/7))
### standardize for invariance to scale
	y=(ys-mean(ys))/sd(ys)
	yst=sort(y)
### define support
	maxd=0
	for(i in 2:n){
		if(yst[i]-yst[i-1]>maxd){maxd=yst[i]-yst[i-1]}
	}
	maxd=max(maxd,(max(y)-min(y))/8)
	sp=c(yst[1]-maxd,yst[n]+maxd)
	kn=c(sp[1],quantile(y,0:(nkn-3)/(nkn-3)),sp[2])
	kn=as.numeric(kn)
	stop=FALSE
### add in some extra knots if needed
	while(!stop){
		newk=c()
		for(i in 2:nkn){
			gap=kn[i]-kn[i-1]
			if(gap>maxd+.01){
				newk=c(newk,i)
			}
		}
		if(length(newk)==0){
			stop=TRUE
		}else{
			for(i in 1:length(newk)){
				kn=c(kn[1:(newk[i]-1)],(kn[newk[i]]+kn[newk[i]-1])/2,kn[newk[i]:nkn])
				nkn=nkn+1
			}
		}
	}
## splines
	bspl=bSpline(y,degree=2,knots=kn[2:(nkn-1)],Boundary.knots=c(kn[1],kn[nkn]),intercept=TRUE)
	bspl0=bSpline(kn,degree=2,knots=kn[2:(nkn-1)],Boundary.knots=c(kn[1],kn[nkn]),intercept=TRUE)
	bspl1=bSpline(kn,degree=2,knots=kn[2:(nkn-1)],Boundary.knots=c(kn[1],kn[nkn]),intercept=TRUE,derivs=1)
	m=dim(bspl)[2]
### penalty
	xm=1:(nkn-1)   ## midpoints of knots
	for(i in 1:(nkn-1)){xm[i]=(kn[i]+kn[i+1])/2}
	d2spl=bSpline(xm,degree=2,knots=kn[2:(nkn-1)],Boundary.knots=c(kn[1],kn[nkn]),intercept=TRUE,derivs=2)
	difmat=matrix(0,nrow=nkn-2,ncol=nkn-1)
	for(i in 1:(nkn-2)){difmat[i,i]=-1;difmat[i,i+1]=1}
	pmat=t(d2spl)%*%t(difmat)%*%difmat%*%d2spl
	pmat=pmat/abs(max(pmat))
### make hmat
	xp=0:1000/1000*(sp[2]-sp[1])+sp[1];np=length(xp)
	dspl=bSpline(xp,degree=2,knots=kn[2:(nkn-1)],Boundary.knots=c(kn[1]-1e-8,kn[nkn]+1e-8),intercept=TRUE)
	dx=xp[4]-xp[3]
	hmat=matrix(0,nrow=m,ncol=m)
	for(i in 1:m){
		for(j in 1:m){
			hmat[i,j]=sum(dspl[,i]*dspl[,j])*dx
		}
	}
### make cvec
	cvec=1:m
	for(i in 1:m){cvec[i]=sum(bspl[,i])/n}
## find modes using a small penalty parameter
	lam=n^(-2/7)*.01	
	mfit1=findmodes(y,kn,bspl,pmat,hmat,cvec,lam)
	modes=mfit1$bestmodes
### do cross-validation to find penalty parameter lam
	nt=round(n/10)
	fhatt=matrix(nrow=10,ncol=nt)
	nf=rep(1:10,nt)
	lams=n^(-2/7)*10^(-4:0)
	nlam=length(lams)
	critlam=1:nlam
	cvecx=1:m
	for(ilam in 1:nlam){
		ans=dens2mode(y,kn,bspl,modes,pmat,hmat,cvec,lam=lams[ilam])
		fp=dspl%*%ans$bhat
		intden=sum(fp)^2*dx
		nfs=sample(nf,n)
		for(ifold in 1:10){
			xs=y[nfs!=ifold]
			xt=y[nfs==ifold]
			for(i in 1:m){cvecx[i]=sum(bspl[nfs!=ifold,i])/n}
			ansx=dens2mode(xs,kn,bspl,modes,pmat,hmat,cvecx,lam=lams[ilam]) 
			for(i in 1:nt){
				fhatxp=dspl%*%ansx$bhat
				ii=which(abs(xp-xt[i])==min(abs(xp-xt[i])))
				fhatt[ifold,i]=fhatxp[ii]
			}
		}
		critlam[ilam]=intden-2*sum(fhatt)/n
	}
	ilammin=min(which(critlam==min(critlam)))
	lamch=lams[ilammin]
## find best mode for unimodal
	mincrit=1e10
	for(i in 2:(nkn-2)){
		ans=dens1mode(y,kn,bspl,i,pmat,hmat,cvec,lam=lamch,eps,bspl0,bspl1)
		crval=t(ans$bhat)%*%(hmat+lamch*pmat)%*%ans$bhat-2*t(cvec)%*%ans$bhat
		if(crval<mincrit){
			mincrit=crval
			bhat1=ans$bhat
			bestmode=i
		}
	}	
## re-find best modes for bimodal using chosen penalty
	mfit=findmodes(y,kn,bspl,pmat,hmat,cvec,lamch)
	modes2=mfit$bestmodes
	bhat2=mfit$bhat
### check if bimodal fit is unimodal
	numsw=0;f2=round(dspl%*%bhat2,7)
	ss=1:(np-1)
	ss[1]=sign(f2[2]-f2[1])
	for(i in 2:(np-1)){
		ss[i]=sign(f2[i+1]-f2[i])
		if(abs(ss[i]-ss[i-1])==2){numsw=numsw+1}
	}
	if(numsw>1){
### compute the value of the test statistic
		cvsim=1:m;ii=1:n;tsim=1:nsim
		t2=t(bhat2)%*%(hmat+lamch*pmat)%*%bhat2-2*t(cvec)%*%bhat2
		t1=t(bhat1)%*%(hmat+lamch*pmat)%*%bhat1-2*t(cvec)%*%bhat1
		tstat=as.numeric(t1-t2)
		f1=dspl%*%bhat1
		cdf1p=f1
		for(i in 2:np){
			cdf1p[i]=cdf1p[i-1]+f1[i]*dx
		}
		cdf1p=cdf1p/cdf1p[np]
		ysim=1:n
		for(isim in 1:nsim){
			u=runif(n)
			for(i in 1:n){
				ii[i]=which(abs(cdf1p-u[i])==min(abs(cdf1p-u[i])))
				ysim[i]=xp[ii[i]]
			}
			for(i in 1:m){cvsim[i]=sum(dspl[ii,i])/n}
## find best mode for unimodal
			mincrit=1e10
			for(i in 2:(nkn-2)){
				ans=dens1mode(ysim,kn,bspl,i,pmat,hmat,cvsim,lam=lamch,eps,bspl0,bspl1)
				crval=t(ans$bhat)%*%(hmat+lamch*pmat)%*%ans$bhat-2*t(cvsim)%*%ans$bhat
				if(crval<mincrit){
					mincrit=crval
					bhat1s=ans$bhat
					bestmodesim=i
				}
			}	
## find best modes for bimodal
			mfits=findmodes(ysim,kn,bspl,pmat,hmat,cvsim,lamch)
			modes2s=mfits$bestmodes
			bhat2s=mfits$bhat
### check if bimodal fit is unimodal
			numsw=0;f2s=round(dspl%*%bhat2s,7)
			ss=1:(np-1)
			ss[1]=sign(f2s[2]-f2s[1])
			for(i in 2:(np-1)){
				ss[i]=sign(f2s[i+1]-f2s[i])
				if(abs(ss[i]-ss[i-1])==2){numsw=numsw+1}
			}
			if(numsw>1){
				t2sim=t(bhat2s)%*%(hmat+lamch*pmat)%*%bhat2s-2*t(cvsim)%*%bhat2s
				t1sim=t(bhat1s)%*%(hmat+lamch*pmat)%*%bhat1s-2*t(cvsim)%*%bhat1s
				tsim[isim]=t1sim-t2sim
			}else{
				tsim[isim]=-1
			}
		}
		pval=sum(tsim>tstat)/nsim
	}else{
		pval=2
	}
	atest=new.env()
	atest$pval=pval
	atest$xp=xp
	atest$fhat1=dspl%*%bhat1
	atest$dspl=atest$dspl
	atest$fhat2=f2
	atest$bhat1=bhat1
	atest$bhat2=bhat2
	atest$kn=kn
	atest$lam=lamch
	atest

}


### find best modes given penalty parameter
findmodes=function(y,kn,bspl,pmat,hmat,cvec,lam,eps){
	n=length(y);m=dim(bspl)[2];nkn=length(kn)	
	bestmodes=rep(0,3)
	mincrit=1e10
	for(i in 1:(nkn-5)){
		for(j in (i+2):(nkn-3)){
			for(k in (j+2):(nkn-1)){
				modes=c(i,j,k)
				ansm=dens2mode(y,kn,bspl,modes,pmat,hmat,cvec,lam)
				crval=t(ansm$bhat)%*%(hmat+lam*pmat)%*%ansm$bhat-2*t(cvec)%*%ansm$bhat
				if(crval<mincrit){
					bestmodes=modes
					mincrit=crval
					bhat=ansm$bhat
				}
			}
		}
	}
	bans=new.env()
	bans$bestmodes=bestmodes
	bans$bhat=bhat
	bans
}



### find the one-mode density estimate given mode & everything
dens1mode=function(y,kn,bspl,mode1,pmat,hmat,cvec,lam,eps,bspl0,bspl1){
	n=length(y);m=dim(bspl)[2];nkn=length(kn);aadd=1:m
	amat=matrix(0,m,m)
#	for(i in 1:(mode1-1)){
#		amat[i,i]=-1;amat[i,i+1]=1
#	}
#	for(i in (mode1+1):nkn){
#		amat[i-1,i]=1;amat[i-1,i+1]=-1
#	}
	amat[1,1]=1;amat[2,m]=1
	amat[3,mode1-2]=-1;amat[3,mode1-1]=1
	amat[4,mode1+1]=1;amat[4,mode1+2]=-1
	e1=max(1,mode1-1)
	for(i in 1:e1){
		for(j in 1:m){
			aadd[j]=bspl1[i,j]-eps*bspl0[i,j]
		}		
		amat=rbind(amat,aadd)
	}
	e2=min(mode1+1,nkn)
	for(i in e2:nkn){
		for(j in 1:m){
			aadd[j]=-eps*bspl0[i,j]-bspl1[i,j]
		}
		amat=rbind(amat,aadd)
	}
	
	umat=chol(hmat+lam*pmat)
	uinv=solve(umat)
	zvec=t(uinv)%*%cvec
	atil=amat%*%uinv
	
	ans=coneA(zvec,atil)
	bhat=uinv%*%ans$thetahat
	fit=new.env()
	fit$bhat=bhat
	fit
}


### find the two-mode density spline coefficient vector given everything we need
dens2mode=function(y,kn,bspl,modes,pmat,hmat,cvec,lam){
	n=length(y);m=dim(bspl)[2];nkn=length(kn)
	amat=matrix(0,m,m)
	for(i in 1:(modes[1]-1)){
		amat[i,i]=-1;amat[i,i+1]=1
	}
	for(i in (modes[1]+1):(modes[2]-1)){
		amat[i-1,i]=1;amat[i-1,i+1]=-1
	}
	for(i in (modes[2]+1):(modes[3]-1)){
		amat[i-2,i]=-1;amat[i-2,i+1]=1
	}
	for(i in (modes[3]+1):nkn){
		amat[i-3,i]=1;amat[i-3,i+1]=-1
	}
	amat[nkn-2,1]=1;amat[nkn-1,m]=1;amat[nkn,modes[2]-1]=1;amat[nkn+1,modes[2]+1]=1
	
	umat=chol(hmat+lam*pmat)
	uinv=solve(umat)
	zvec=t(uinv)%*%cvec
	atil=amat%*%uinv
	
	ans=coneA(zvec,atil)
	bhat=uinv%*%ans$thetahat
	fit=new.env()
	fit$bhat=bhat
	fit
}