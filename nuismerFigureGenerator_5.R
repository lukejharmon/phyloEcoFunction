library(geiger)
library(colorspace)
library(RColorBrewer)

pairwiseInteractions<-function(phy, alpha=0.0001, G=1, N=1000, s2=1) {
	ntips<-length(phy$tip.label)

	t<-vcv.phylo(phy) # this only works for ultrametric trees

	f<-rep(1/ntips, ntips)
	s2z<-rep(s2, ntips)

	Tau<-max(branching.times(phy))

	P<-0
	for(i in 1:ntips)
		for(j in 1:ntips) {
			P<-P+f[i]*f[j]*(1 - 2*alpha*G/N*(Tau-t[i,j]) - alpha*(s2z[i]+s2z[j]))
			#print(P)
			}
		
	P
}


pairwiseInteractionsOU<-function(phy, alpha=0.0001, G=1, N=1000, s2=1, gamma=0.0001) {
	ntips<-length(phy$tip.label)

	t<-vcv.phylo(phy) # this only works for ultrametric trees

	Tau<-max(branching.times(phy))

	f<-rep(1/ntips, ntips)
	s2z<-rep(s2, ntips)

	P<-0
	for(i in 1:ntips)
		for(j in 1:ntips) {
			f1<-alpha*(1-exp(-4*gamma*G*(Tau-t[i,j])))/(2*N*gamma)
			f2<-alpha*(s2z[i]+s2z[j])
			P<-P+f[i]*f[j]*(1 - f1 - f2)
			#print(P)
		}
	P
}

pairwiseInteractionsMatrix<-function(mm, s2=1, alpha=0.0001) {
	ntips<-dim(mm)[1]

	f<-rep(1/ntips, ntips)
	s2z<-rep(s2, ntips)

	P<-0
	for(i in 1:ntips)
		for(j in 1:ntips) {
			P<-P+f[i]*f[j]*(1 - alpha*(mm[i,i] + mm[j,j] - 2*mm[i,j]) - alpha*(s2z[i]+s2z[j]))
			#print(P)
			}
		
	P
}

wholePairwiseInteractionMatrix<-function(mm, s2=1, alpha=0.0001) {
	ntips<-dim(mm)[1]

	f<-rep(1/ntips, ntips)
	s2z<-rep(s2, ntips)

	P<-mm
	P[,]<-0

	for(i in 1:ntips)
		for(j in 1:ntips) {
			P[i,j]<-f[i]*f[j]*(1 - alpha*(mm[i,i] + mm[j,j] - 2*mm[i,j]) - alpha*(s2z[i]+s2z[j]))
			cat(alpha*(mm[i,i] + mm[j,j] - 2*mm[i,j]), "\n")
			}
		
	P
}



treeB<-read.tree(text="((A:1,B:1):1,(C:1, D:1):1);")

makeTreeA<-function(t1Ratio, tTot) {
 	t1<-tTot * t1Ratio
	treeAbase<-read.tree(text="(((A:1,B:1):1,C:2):1,D:3);")
	newBranches<-c(t1, t1, tTot-2*t1, tTot-2*t1, tTot-t1, tTot)
	treeAbase$edge.length<-newBranches
	treeAbase
}

makeTreeB<-function(t1Ratio, tTot) {
 	t1<-tTot * t1Ratio
	treeBbase<-read.tree(text="((A:1,B:1):1,(C:1, D:1):1);")
	newBranches<-c(t1, tTot-t1, tTot-t1, t1, tTot-t1, tTot-t1)
	treeBbase$edge.length<-newBranches
	treeBbase
}


# figure one

pdf("~/Documents/nuismerGrant/traitPaper/paper/figures/nuisHar_fig1_treeA.pdf", width=12, height=12)

tt<-makeTreeA(t1Ratio=0.2, tTot=1000)
plot(tt)

# difference model
ppMatrix<-matrix(nrow=100, ncol=100)
ppMatrix[]<-(0.25*0.25/2)*16
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))


# matching model
totalInt<-seq(from=10000, to=1000000, length.out=100)
frac<-seq(from=0.001, to=0.499, length.out=100)
ppMatrix<-matrix(nrow=100, ncol=100)

for(i in 1:100)
	for(j in 1:100) {
 		tt<-makeTreeA(t1Ratio=frac[i], tTot=totalInt[j])
		ppMatrix[i,j]<-pairwiseInteractions(tt, alpha=0.0005)
	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))

for(i in 1:100)
	for(j in 1:100) {
		tt<-makeTreeA(t1Ratio=frac[i], tTot=totalInt[j])
		ppMatrix[i,j]<-pairwiseInteractionsOU(tt, alpha=0.0005, gamma=0.000001)
	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))




dev.off()



pdf("~/Documents/nuismerGrant/traitPaper/paper/figures/nuisHar_fig1_treeB.pdf", width=12, height=12)

tt<-makeTreeB(t1Ratio=0.2, tTot=1000)
plot(tt)

# difference model
ppMatrix<-matrix(nrow=100, ncol=100)
ppMatrix[]<-(0.25*0.25/2)*16
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))


# matching model
totalInt<-seq(from=10000, to=1000000, length.out=100)
frac<-seq(from=0.001, to=0.5, length.out=100)
ppMatrix<-matrix(nrow=100, ncol=100)

for(i in 1:100)
	for(j in 1:100) {
 		tt<-makeTreeB(t1Ratio=frac[i], tTot=totalInt[j])
		ppMatrix[i,j]<-pairwiseInteractions(tt, alpha=0.0005)
	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))

for(i in 1:100)
	for(j in 1:100) {
		tt<-makeTreeB(t1Ratio=frac[i], tTot=totalInt[j])
		ppMatrix[i,j]<-pairwiseInteractionsOU(tt, alpha=0.0005, gamma=0.000001)
	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))




dev.off()





pdf("~/Documents/nuismerGrant/traitPaper/paper/figures/nuisHar_newBoxes.pdf", width=30, height=12)
pp<-diverge_hcl(100)



layout(matrix(1:10, nrow=2, byrow=T))

Sm<- 0
G<-1
gamma <- 0
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

treeA<-read.tree(text="(((A:10000.0,B:10000.0):10000.0,C:20000.0):10000.0,D:30000.0);")
treeB<-read.tree(text="((A:10000.0,B:10000.0):10000.0,(C:10000.0, D:10000.0):10000.0);")



ppA1<-predictWholeTree(treeA, z0, psi, Sm, eps, theta, rootBranch=10000)
piA1<-wholePairwiseInteractionMatrix(ppA1$V, alpha=0.01)

interPlot(piA1[4:1,], nn=nrow(piA1), pp=pp, gMax=1/16)

Sm<- 0
G<-1
gamma <- 0.00001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

ppA2<-predictWholeTree(treeA, z0, psi, Sm, eps, theta, rootBranch=10000)
piA2<-wholePairwiseInteractionMatrix(ppA2$V, alpha=0.01)

interPlot(piA2[4:1,], nn=nrow(piA2), pp=pp, gMax=0.07)


Sm<- -0.00003

G<-1
gamma <- 0.00001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

ppA3<-predictWholeTree(treeA, z0, psi, Sm, eps, theta, rootBranch=10000)
piA3<-wholePairwiseInteractionMatrix(ppA3$V, alpha=0.01)

interPlot(piA3[4:1,], nn=nrow(piA3), pp=pp, gMax=0.07)

Sm<- 0.00003
G<-1
gamma <- 0.00001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

ppA4<-predictWholeTree(treeA, z0, psi, Sm, eps, theta, rootBranch=10000)
piA4<-wholePairwiseInteractionMatrix(ppA4$V, alpha=0.01)

interPlot(piA4[4:1,], nn=nrow(piA3), pp=pp, gMax=0.07)


Sm<- 0
G<-1
gamma <- 0
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

plot(treeB)
ppB1<-predictWholeTree(treeB, z0, psi, Sm, eps, theta, rootBranch=10000)
piB1<-wholePairwiseInteractionMatrix(ppB1$V, alpha=0.01)

interPlot(piB1[4:1,], nn=nrow(piB1), pp=pp, gMax=0.07)

Sm<- 0
G<-1
gamma <- 0.00001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

ppB2<-predictWholeTree(treeB, z0, psi, Sm, eps, theta, rootBranch=10000)
piB2<-wholePairwiseInteractionMatrix(ppB2$V, alpha=0.01)

interPlot(piB2[4:1,], nn=nrow(piB2), pp=pp, gMax=0.07)


Sm<- -0.00003
G<-1
gamma <- 0.00001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

ppB3<-predictWholeTree(treeB, z0, psi, Sm, eps, theta, rootBranch=10000)
piB3<-wholePairwiseInteractionMatrix(ppB3$V, alpha=0.01)

interPlot(piB3[4:1,], nn=nrow(piB3), pp=pp, gMax=0.07)

Sm<- 0.00003
G<-1
gamma <- 0.00001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

ppB4<-predictWholeTree(treeB, z0, psi, Sm, eps, theta, rootBranch=10000)
piB4<-wholePairwiseInteractionMatrix(ppB4$V, alpha=0.01)

interPlot(piB4[4:1,], nn=nrow(piB4), pp=pp, gMax=0.07)


dev.off()

pdf("~/Documents/nuismerGrant/traitPaper/paper/figures/nuisHar_palette.pdf", width=6, height=6)

plot("", xlim=c(0.5,10+0.5), ylim=c(0.5,1+0.5), axes=F)
for(i in 1:10)
		symbols(x=i, y=1, squares=1.0, fg= i*10, pch=15, bg= i*10, pch=15, cex=10, add=T, inches=F)
dev.off()




# competition and mutualism models

pdf("~/Documents/nuismerGrant/traitPaper/paper/figures/nuisHar_fig4_raw.pdf")


Sm<- -0.000001
G<-1
gamma <- 0.000001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

totalInt<-seq(from=10000, to=1000000, length.out=100)
frac<-seq(from=0.001, to=0.499, length.out=100)
ppMatrix1<-matrix(nrow=100, ncol=100)
ppMatrix2<-matrix(nrow=100, ncol=100)

for(i in 1:100)
	for(j in 1:100) {
 		tt<-makeTreeA(t1Ratio=frac[i], tTot=totalInt[j])
 		vv<-predictWholeTree(tt, z0, psi, Sm, eps, theta, rootBranch=0)
 		oo<-order(rownames(vv$V))
 
 
		ppMatrix[i,j]<-pairwiseInteractionsMatrix(vv$V[oo,oo], alpha=0.0005)

	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))


Sm<- -0.000001
G<-1
gamma <- 0.000001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

totalInt<-seq(from=10000, to=1000000, length.out=100)
frac<-seq(from=0.001, to=0.499, length.out=100)
ppMatrix1<-matrix(nrow=100, ncol=100)
ppMatrix2<-matrix(nrow=100, ncol=100)

for(i in 1:100)
	for(j in 1:100) {
 		tt<-makeTreeB(t1Ratio=frac[i], tTot=totalInt[j])
 		vv<-predictWholeTree(tt, z0, psi, Sm, eps, theta, rootBranch=0)
 		oo<-order(rownames(vv$V))
 
 
		ppMatrix[i,j]<-pairwiseInteractionsMatrix(vv$V[oo,oo], alpha=0.0005)

	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))



Sm<- 0.000001
G<-1
gamma <- 0.000001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

totalInt<-seq(from=10000, to=1000000, length.out=100)
frac<-seq(from=0.001, to=0.499, length.out=100)
ppMatrix1<-matrix(nrow=100, ncol=100)
ppMatrix2<-matrix(nrow=100, ncol=100)

for(i in 1:100)
	for(j in 1:100) {
 		tt<-makeTreeA(t1Ratio=frac[i], tTot=totalInt[j])
 		vv<-predictWholeTree(tt, z0, psi, Sm, eps, theta, rootBranch=0)
 		oo<-order(rownames(vv$V))
 
 
		ppMatrix[i,j]<-pairwiseInteractionsMatrix(vv$V[oo,oo], alpha=0.0005)

	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))


Sm<- 0.000001
G<-1
gamma <- 0.000001
N<-1000
eps<-G/N
theta<-10
z0<-10
psi<- 2 * gamma * G

totalInt<-seq(from=10000, to=1000000, length.out=100)
frac<-seq(from=0.001, to=0.499, length.out=100)
ppMatrix1<-matrix(nrow=100, ncol=100)
ppMatrix2<-matrix(nrow=100, ncol=100)

for(i in 1:100)
	for(j in 1:100) {
 		tt<-makeTreeB(t1Ratio=frac[i], tTot=totalInt[j])
 		vv<-predictWholeTree(tt, z0, psi, Sm, eps, theta, rootBranch=0)
 		oo<-order(rownames(vv$V))
 
 
		ppMatrix[i,j]<-pairwiseInteractionsMatrix(vv$V[oo,oo], alpha=0.0005)

	}
ppMatrix[ppMatrix<0]<-0
#contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, method="simple")
filled.contour(x=frac, y=totalInt, z= ppMatrix, levels=1:10/10, col=diverge_hcl(10))

dev.off()

