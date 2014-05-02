dyn.load("~/Documents/nuismerGrant/traitPaper/rCode/VCV_Iterate_forR.so")


predictWholeTree<-function(tt, z0, psi, Sm, eps, theta, rootBranch) {
	bb<-sort(branching.times(tt), decreasing=T)
	
	liveTips<-names(bb)[1]
	
	currV<-as.matrix(0)
	currM<-z0
	
	r1<-.C("competitionInterval", mm=as.double(currM), VCV=as.double(currV), S=as.double(Sm), Psi=as.double(psi), Eps=as.double(eps), theta=as.double(theta), Nspecies=as.integer(1.0), Gmax=as.integer(rootBranch))

	currV<-matrix(r1$VCV, nrow=1, ncol=1)
	currM<-r1$mm


	for(i in 1:(length(bb)-1)) {
		speciatingLineage<-names(bb)[i]
		daughterRows<-which(tt$edge[,1]==speciatingLineage)
		daughters<-tt$edge[daughterRows,2]
		replaceMe<-which(liveTips== speciatingLineage)
		liveTips<-c(liveTips[-replaceMe], daughters)
		
		newV1<-currV[-replaceMe,-replaceMe]
		newV2<-cbind(currV[-replaceMe, replaceMe], currV[-replaceMe, replaceMe])
		newV3<-rbind(currV[replaceMe, -replaceMe], currV[replaceMe, -replaceMe])
		newV4<-matrix(currV[replaceMe, replaceMe], nrow=2, ncol=2)
		
		newV<-rbind(cbind(newV1, newV2), cbind(newV3, newV4))
		
		newM<-c(currM[-replaceMe], currM[replaceMe],currM[replaceMe])

		branchLength<-round(bb[i]-bb[i+1])
		
		ns<-dim(newV)[1]
		res<-.C("competitionInterval", mm=as.double(newM), VCV=as.double(newV), S=as.double(Sm), Psi=as.double(psi), Eps=as.double(eps), theta=as.double(theta), Nspecies=as.integer(ns), Gmax=as.integer(branchLength))
		
		currV<-matrix(res$VCV, nrow=ns, ncol=ns)
		currM<-res$mm

	}
	
	i<-length(bb)
	speciatingLineage<-names(bb)[i]
		daughterRows<-which(tt$edge[,1]==speciatingLineage)
		daughters<-tt$edge[daughterRows,2]
		replaceMe<-which(liveTips== speciatingLineage)
		liveTips<-c(liveTips[-replaceMe], daughters)
		newV1<-currV[-replaceMe,-replaceMe]
		newV2<-cbind(currV[-replaceMe, replaceMe], currV[-replaceMe, replaceMe])
		newV3<-rbind(currV[replaceMe, -replaceMe], currV[replaceMe, -replaceMe])
		newV4<-matrix(currV[replaceMe, replaceMe], nrow=2, ncol=2)
		
		newV<-rbind(cbind(newV1, newV2), cbind(newV3, newV4))
		
		newM<-c(currM[-replaceMe], currM[replaceMe],currM[replaceMe])

		branchLength<-round(bb[i])
		
		ns<-dim(newV)[1]
		res<-.C("competitionInterval", mm=as.double(newM), VCV=as.double(newV), S=as.double(Sm), Psi=as.double(psi), Eps=as.double(eps), theta=as.double(theta), Nspecies=as.integer(ns), Gmax=as.integer(branchLength))
		
		currV<-matrix(res$VCV, nrow=ns, ncol=ns)
		currM<-res$mm
	
	
	rownames(currV)<-tt$tip.label[as.numeric(liveTips)]
	colnames(currV)<-tt$tip.label[as.numeric(liveTips)]
	names(currM)<-tt$tip.label[as.numeric(liveTips)]
	return(list(M=currM, V=currV))
}


vcvPlot<-function(vcv, nn, gMax=NULL, pp=diverge_hcl(100, c = c(100, 0), l = c(70, 0), h=c(250, 405)), mult=1) {
	plot("", xlim=c(0.5,nn+0.5), ylim=c(0.5,nn+0.5), axes=F)

	if(is.null(gMax)) gMax<-max(abs(vcv))
	vcv<- vcv/gMax * mult
	vcv[vcv > 1] <- 1
	vcv[vcv < -1] <- -1
	#palette(heatmap(100))
	palette(pp)

	#palette(diverge_hcl(100, c = 200, l = c(40, 90)))
	dd<-dim(vcv)[1]
	
	for(i in 1:dd)
		for(j in 1:dd)
			symbols(x=i, y=j, squares=1.0, fg= vcv[i,j]*50+50, pch=15, bg= vcv[i,j]*50+50, pch=15, cex=10, add=T, inches=F)
	
}


interPlot<-function(vcv, nn, gMax=0.025, pp=diverge_hcl(100, c = c(100, 0), l = c(70, 0), h=c(250, 405)), mult=1) {
	plot("", xlim=c(0.5,nn+0.5), ylim=c(0.5,nn+0.5), axes=F)

	vcv<- vcv/gMax

	vcv[vcv > 1] <- 1
	vcv[vcv < 0] <- 0
	#palette(heatmap(100))
	palette(pp)

	#palette(diverge_hcl(100, c = 200, l = c(40, 90)))
	dd<-dim(vcv)[1]
	
	for(i in 1:dd)
		for(j in 1:dd)
			symbols(x=i, y=j, squares=1.0, fg= vcv[i,j]*100, pch=15, bg= vcv[i,j]*100, pch=15, cex=10, add=T, inches=F)
	
}

vcvplotWholeTree<-function(tt, z0, psi, Sm, eps, theta, rootBranch, pp=diverge_hcl(100, c = c(100, 0), l = c(70, 0), h=c(250, 405)), mult=1) {
	bb<-sort(branching.times(tt), decreasing=T)
	
	nn<-length(tt$tip.label)
	
	liveTips<-names(bb)[1]
	
	currV<-as.matrix(0)
	currM<-z0
	
	r1<-.C("competitionInterval", mm=as.double(currM), VCV=as.double(currV), S=as.double(Sm), Psi=as.double(psi), Eps=as.double(eps), theta=as.double(theta), Nspecies=as.integer(1.0), Gmax=as.integer(rootBranch))

	currV<-matrix(r1$VCV, nrow=1, ncol=1)
	currM<-r1$mm
	#vcvPlot(currV, nn)

	allVCV<-list()

	for(i in 1:(length(bb)-1)) {
		speciatingLineage<-names(bb)[i]
		daughterRows<-which(tt$edge[,1]==speciatingLineage)
		daughters<-tt$edge[daughterRows,2]
		replaceMe<-which(liveTips== speciatingLineage)
		liveTips<-c(liveTips[-replaceMe], daughters)
		
		newV1<-currV[-replaceMe,-replaceMe]
		newV2<-cbind(currV[-replaceMe, replaceMe], currV[-replaceMe, replaceMe])
		newV3<-rbind(currV[replaceMe, -replaceMe], currV[replaceMe, -replaceMe])
		newV4<-matrix(currV[replaceMe, replaceMe], nrow=2, ncol=2)
		
		newV<-rbind(cbind(newV1, newV2), cbind(newV3, newV4))
		
		newM<-c(currM[-replaceMe], currM[replaceMe],currM[replaceMe])

		branchLength<-round(bb[i]-bb[i+1])
		
		ns<-dim(newV)[1]
		res<-.C("competitionInterval", mm=as.double(newM), VCV=as.double(newV), S=as.double(Sm), Psi=as.double(psi), Eps=as.double(eps), theta=as.double(theta), Nspecies=as.integer(ns), Gmax=as.integer(branchLength))
		
		currV<-matrix(res$VCV, nrow=ns, ncol=ns)
		currM<-res$mm
		allVCV[[i]]<-currV

	}
	
	i<-length(bb)
	speciatingLineage<-names(bb)[i]
		daughterRows<-which(tt$edge[,1]==speciatingLineage)
		daughters<-tt$edge[daughterRows,2]
		replaceMe<-which(liveTips== speciatingLineage)
		liveTips<-c(liveTips[-replaceMe], daughters)
		newV1<-currV[-replaceMe,-replaceMe]
		newV2<-cbind(currV[-replaceMe, replaceMe], currV[-replaceMe, replaceMe])
		newV3<-rbind(currV[replaceMe, -replaceMe], currV[replaceMe, -replaceMe])
		newV4<-matrix(currV[replaceMe, replaceMe], nrow=2, ncol=2)
		
		newV<-rbind(cbind(newV1, newV2), cbind(newV3, newV4))
		
		newM<-c(currM[-replaceMe], currM[replaceMe],currM[replaceMe])

		branchLength<-round(bb[i])
		
		ns<-dim(newV)[1]
		res<-.C("competitionInterval", mm=as.double(newM), VCV=as.double(newV), S=as.double(Sm), Psi=as.double(psi), Eps=as.double(eps), theta=as.double(theta), Nspecies=as.integer(ns), Gmax=as.integer(branchLength))
		
		currV<-matrix(res$VCV, nrow=ns, ncol=ns)
		currM<-res$mm
	
		allVCV[[i]]<-currV
	
	globalMax<-0
	for(i in 1:length(bb)) {
		mm<-max(abs(allVCV[[i]]))
		if(mm>globalMax) globalMax<-mm		
	}
	for(i in 1:length(bb)) {
		vcvPlot(allVCV[[i]], nn=nn, gMax=globalMax, pp=pp, mult=mult)
	}
}


lnlCompetition<-function(phy, data, z0, psi, Sm, eps, theta, rootBranch) {
	mvnPred<-predictWholeTree(phy, z0, psi, Sm, eps, theta, rootBranch)
	mm<-match(rownames(mvnPred$V),names(data))
	dd<-data[mm]
	lnl<-dmvnorm(dd, mean=mvnPred$M, sigma=mvnPred$V, log=T)
	return(lnl)
}

simpleSimulate<-function(tt, z0, psi, Sm, eps, theta, rootBranch) {
	mvModel<-predictWholeTree(tt, z0, psi, Sm, eps, theta, rootBranch)
	dd<-rmvnorm(1, mean=mvModel$M, sigma=mvModel$V)[1,]
	return(dd)
}

interpretOutput<-function(pp) {
	z0<-pp[1]
	psi<-exp(pp[2])
	Sm<-psi*pp[3]
	eps <- exp(pp[4])
	theta<-pp[5]
	return(list(z0=z0, psi=psi, Sm=Sm, eps=eps, theta=theta))
}

maxLikelihood<-function(phy, data) {

	foo<-function(x) {
		z0<-x[1]
		psi<-exp(x[2])
		Sm<-psi*x[3]
		eps <- exp(x[4])
		theta<-x[5]
	
		res<-lnlCompetition(phy, data, z0, psi, Sm, eps, theta, rootBranch)

		return(-res)
	}
	
	o<-optim(par=c(10, log(0.001), 0.5, log(0.00001), 10), foo, method="L-BFGS-B",lower=c(-1,-10,-0.999,-20,-10), upper=c(100,-2,0.999,3,20))

	return(o)
}


