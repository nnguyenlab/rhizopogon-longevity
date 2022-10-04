##Commands written by Kabir Peay, 2020 

##import the data
setwd("~/Documents/Projects/Spore Longevity/SporeData_2020")
sporeData<-read.csv("spore.longevity.060120.csv",as.is=TRUE)
##function below creates a logistic regression graph that 
##plots each year of a species individually

#temporarily rounding all the non-integer data
sporeData$colonized<-round(sporeData$colonized)

library("MASS")

logisticGraph<-function(SPECIES, dataSet, TITLE) {

#create subsets of the data for species of interest
subsetSp<-subset(dataSet,dataSet$species==SPECIES)
percent.colon<-subsetSp$colonized/(subsetSp$colonized+subsetSp$uncolonized)
spores<-subsetSp$spores.ml
time<-sort(as.numeric(unique(subsetSp$year)))

P50<-rep(NA,7)
P50.SE<-rep(NA,7)

##Set up the initial plot
G1<-plot(log10(spores), percent.colon, pch="", ylab="% of seedlings of colonized", xlab=expression(log[10]~spores~per~ml), main=TITLE, cex.main=1.4, cex.lab=1.4, cex.axis=1.4)
#COLOR<-c("#DDDC32","#75D054","#47C16E","#27AD81","#1F9A8A","#24868E","#2C728E")
COLOR<-c("#E3E418","#AADC32","#5DC863","#27AD81","#21908C","#2C728E","#443A83")#viridis
#COLOR<-c("#f1eef6","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#034e7b")#blues
#COLOR<-c("#FDE275","#8FD744","#35B779","#21908C","#31688E","#443A83","#440154")
#COLOR<-heat.colors(length(time), alpha=0.8)

##Create a for loop that plots each year
P50<-NA

for(i in 1:length(time)){
  
#create the variables
temp<-subset(subsetSp,subsetSp$year==time[i])
spores2<-temp$spores[which(temp$year==time[i])]
colon2<-temp$colonized[which(temp$year==time[i])]
uncolon2<-temp$uncolonized[which(temp$year==time[i])]
prop.colon2<-cbind(colon2, uncolon2)

##build the model
temp.model<-glm(prop.colon2~log10(spores2), binomial)

#extract 50% colonization values
LD50 <- dose.p(temp.model, p = 0.5)
LD01 <- dose.p(temp.model, p = 0.01)
P50[i]<-LD50[[1]]
#P01[i]<-LD01[[1]]
P01<-0 #LD01[[1]] gives a negative number near zero but R could not subset it so I made it 0.
P50.SE[i]<-attributes(LD50)$SE[,1]

##make the graph
p<-colon2/(colon2+uncolon2)
xv<-seq(min(log10(spores2)),max(log10(spores2)),0.1)
xAxis<-data.frame(xv=log10(spores2))
points(log10(spores2),p, pch=21,cex=1.5,lwd=1.3,bg=COLOR[i],col=COLOR[i])
predicted<-predict(temp.model,xAxis, type="resp")
lines(xAxis[,1], predicted, col=COLOR[i],lwd=2.5,lty="dashed")

points(x=LD50[[1]],y=0.5,pch=21,cex=3,col="black",bg="black")
text(x=LD50[[1]],y=0.5,labels=time[i],cex=1,col="white")

legend.text=paste("Yr ", time, sep="")
#legend(3.25,0.3, legend=c(legend.text), pch=21, bty="n",pt.bg=c(heat.colors(7,alpha=0.8)), col=c(heat.colors(7)),cex=1, pt.cex=1, pt.lwd=1.4) 
legend(3.25,0.3, legend=c(legend.text), pch=21, bty="n", pt.bg=COLOR, col=COLOR, cex=1.1, pt.cex=1, pt.lwd=1.4)

	}

P50<-data.frame(time,P50,P50.SE)

return(P50)

}

##Calls the function
par(mfrow=c(1,3))
RV<-logisticGraph("RV",sporeData, substitute(paste(italic('Rhizopogon vulgaris')))); RV$species<-rep("RV",length(RV$time))#; write.csv(RV,"RV.P50.csv")
RO<-logisticGraph("RO",sporeData, substitute(paste(italic('Rhizopogon occidentalis')))); RO$species<-rep("RO",length(RO$time))#; write.csv(RO,"RO.P50.csv")
RS<-logisticGraph("RS",sporeData, substitute(paste(italic('Rhizopogon salebrosus')))); RS$species<-rep("RS",length(RS$time))#; write.csv(RS,"RS.P50.csv")

##Finding significance between the P50 of each species over
combined.pogon<-rbind(RV,RO,RS)
combined.pogon$time2<-combined.pogon$time^2
combined.model<-with(combined.pogon,lm(P50~time+time2+species));summary(combined.model)

#note - I eliminated higher order terms 1 at a time and none of them were significant.
---------
par(mfrow=c(1,1))
#-#-# Rhizopogon vulgaris
RV$time2<-RV$time^2
rv.dd.model<-with(RV,lm(P50~time+time2)); summary(rv.dd.model)
plot(RV$time,RV$P50,pch=21,col="#d95f02",bg="#d95f02",cex=2, xlab="Years buried", ylab=expression(log[10]~spores~required~to~colonize~"50%"~seedlings~(C[50])), ylim=c(-0.5,3),xlim=c(0,20))
#xlab=expression(log[10]~spores~per~ml)
#lines(RV$time,RV$P50,col="red",lty="dashed")
curve((2.386023-0.990677)-0.274695*x+0.016789*x^2,0,100, col="#d95f02",lty="dashed",add=TRUE,lwd=2)

#-#-# Rhizopogon salebrosus
RS$time2<-RS$time^2
rs.dd.model<-with(RS,lm(P50~time+time2)); summary(rs.dd.model)
points(RS$time,RS$P50,pch=21,col="#7570b3",bg="#7570b3",cex=2)
#curve(2.35382-0.31277*x+0.01772*x^2,0,15, col="gold",lty="dashed",add=TRUE,lwd=2)
curve(2.386023-0.175395-0.274695*x+0.016789*x^2,0,40, col="#7570b3",lty="dashed",add=TRUE,lwd=2) #combined model

#-#-# Rhizopogon occidentalis
RO$time2<-RO$time^2
ro.dd.model<-with(RO,lm(P50~time+time2)); summary(ro.dd.model)
points(RO$time,RO$P50,pch=21,col="#1b9e77",bg="#1b9e77",cex=2)
#curve(2.187238-0.207737*x+0.014108*x^2,0,15, col="green",lty="dashed",add=TRUE,lwd=2)
curve(2.386023-0.274695*x+0.016789*x^2,0,40, col="#1b9e77",lty="dashed",add=TRUE,lwd=2) #combined model

legend("bottomright",legend=c(substitute(paste(italic('R. vulgaris'))), 
                              substitute(paste(italic('R. salebrosus'))), 
                              substitute(paste(italic('R. occidentalis')))), 
       pch=21, col=c("#d95f02","#7570b3","#1b9e77"), pt.bg=c("#d95f02","#7570b3","#1b9e77"), bty="n",cex=0.9)

#fit one model for all three species; this was not used in the manuscript.
#curve(1.450938-0.303578*x+0.018541*x^2,0,15, col="red",lty="dashed",add=TRUE,lwd=2)