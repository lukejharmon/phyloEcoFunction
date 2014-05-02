bmdd<-read.csv("~/Documents/nuismerGrant/traitPaper/paper/figures/traitPlot/Fossils_BM_Match_Figures.csv")

makeTreePlots<-function(file) {

dd<-read.csv(file)	
times<-as.numeric(as.character(dd[1:10001,1]))
species<-matrix(nrow=10001, ncol=12)

for(i in 1:12) species[,i]<-as.numeric(as.character(dd[1:10001, i+1]))

palette(gray(1:12/12))
plot(times, species[,1], type="l", lwd=2, ylim=c(-20,20))
for(i in 2:12) {
	lines(times, species[,i], col=i)
}

times2<-as.numeric(as.character(dd[10003:20003,1]))
species2<-matrix(nrow=10001, ncol=12)

for(i in 1:12) species2[,i]<-as.numeric(as.character(dd[10003:20003, i+1]))

palette(gray(1:12/12))
plot(times2, species2[,1], type="l", lwd=2, ylim=c(-20,20))
for(i in 2:12) {
	lines(times2, species2[,i], col=i)
}

}

pdf("~/Documents/nuismerGrant/traitPaper/paper/figures/traitPlot/traitPlot_v2.pdf", width=30, height=12)
layout(matrix(c(9, 10, 1:8), nrow=2, ncol=5))
makeTreePlots("~/Documents/nuismerGrant/traitPaper/paper/figures/traitPlot/Fossils_BM_Match_Figures.csv")
makeTreePlots("~/Documents/nuismerGrant/traitPaper/paper/figures/traitPlot/Fossils_OU_Match_Figures.csv")
makeTreePlots("~/Documents/nuismerGrant/traitPaper/paper/figures/traitPlot/Fossils_NH_Match_Comp_Figures.csv")
makeTreePlots("~/Documents/nuismerGrant/traitPaper/paper/figures/traitPlot/Fossils_NH_Match_Mut_Figures.csv")

dev.off()