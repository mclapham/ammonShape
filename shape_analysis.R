
#read ammonoid measurement data
ammon<-read.csv("http://files.figshare.com/1846468/ammon.csv")

#read time interval data
time_periods<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=1&min_ma=66&max_ma=425")
stages<-subset(time_periods,time_periods$level==5) #extracts stage level data
stages$interval_name<-factor(stages$interval_name) #resets factor levels for stage names
periods<-subset(time_periods,time_periods$level==3) #extracts period level data

#rearranges data so that each specimen is in a single row
ammonD<-subset(ammon,ammon$measurement=="diameter") #extract shell diameter measurements
ammonW<-subset(ammon,ammon$measurement=="width") #extract whorl width measurements
ammonH<-subset(ammon,ammon$measurement=="height") #extract whorl height measurements
ammonU<-subset(ammon,ammon$measurement=="inflation") #extract umbilical diameter measurements

specimen_ID<-levels(ammon$specimen.ID) #procedure uses specimen ID numbers as identifier to match W, D, H, and U values

D<-ammonD$mean[match(specimen_ID,ammonD$specimen.ID)] #shell diameter matched to specimen ID
W<-ammonW$mean[match(specimen_ID,ammonW$specimen.ID)] #whorl width matched to specimen ID
H<-ammonH$mean[match(specimen_ID,ammonH$specimen.ID)] #whorl height matched to specimen ID
U<-ammonU$mean[match(specimen_ID,ammonU$specimen.ID)] #umbilical diameter matched to specimen ID

stage<-ammon$stage[match(specimen_ID,ammon$specimen.ID)] #geological stage matched to specimen ID
species<-ammon$species[match(specimen_ID,ammon$specimen.ID)] #species name matched to specimen ID
coll_no<-ammon$collection_no[match(specimen_ID,ammon$specimen.ID)] #PBDB collection number matched to specimen ID

ammon_shape<-data.frame(specimen_ID,coll_no,species,stage,D,W,H,U) #combine values into data frame

ammon_shape<-na.omit(ammon_shape) #remove specimens with missing measurements

ammon_shape$WD<-ammon_shape$W/ammon_shape$D #calculate W/D ratio
ammon_shape$HD<-ammon_shape$H/ammon_shape$D #calculate H/D ratio
ammon_shape$UD<-ammon_shape$U/ammon_shape$D #calculate U/D ratio

ammon_shape$size<-(0.5*ammon_shape$D)^2*ammon_shape$W #shell volume assuming cylinder shape


#principal components analysis
ammon_pca<-prcomp(ammon_shape[,9:11])

pca_res<-data.frame(species=ammon_shape$species,stage=ammon_shape$stage,ammon_pca$x)

pca_res$D<-ammon_shape$D
pca_res$WD<-ammon_shape$WD
pca_res$HD<-ammon_shape$HD
pca_res$UD<-ammon_shape$UD
pca_res$size<-ammon_shape$size
pca_res$stage_ma<-rowMeans(cbind(stages$early_age[match(pca_res$stage,stages$interval_name)],stages$late_age[match(pca_res$stage,stages$interval_name)]))


#reads species taxonomy
ammon_species<-read.csv("http://paleobiodb.org/data1.1/taxa/list.txt?name=Ammonoidea&rel=all_children&rank=species&show=phylo&limit=99999")
pca_res$order<-ammon_species$order[match(pca_res$species,ammon_species$taxon_name)]
pca_res$family<-ammon_species$family[match(pca_res$species,ammon_species$taxon_name)]

#adds geological timescale colors to each
pca_res$color<-stages$color[match(pca_res$stage,stages$interval_name)]


specimen_ct<-sapply(split(ammon_shape,ammon_shape$stage),nrow)
species_ct<-sapply(split(ammon_shape$species,ammon_shape$stage),function(x) length(unique(x)))

#plot all specimens
par(mgp=c(2,0.75,0))
plot(pca_res$PC1,pca_res$PC2,xlab="PC1",ylab="PC2",pch=16,col=paste(pca_res$color))
arrows(0,0,ammon_pca$rotation[1,1]/1.3,ammon_pca$rotation[1,2]/1.3,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[1,1]/1.25,ammon_pca$rotation[1,2]/1.25,rownames(ammon_pca$rotation)[1])
arrows(0,0,ammon_pca$rotation[2,1]/1.3,ammon_pca$rotation[2,2]/1.3,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[2,1]/1.25,ammon_pca$rotation[2,2]/1.25,rownames(ammon_pca$rotation)[2])
arrows(0,0,ammon_pca$rotation[3,1]/1.3,ammon_pca$rotation[3,2]/1.3,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[3,1]/1.25,ammon_pca$rotation[3,2]/1.25,rownames(ammon_pca$rotation)[3])


#PCA results per stage
par(mfrow=c(1,2))
par(mar=c(3.5,3,1,0.5))
plot(pca_res$PC1,pca_res$PC2,xlab="PC1",ylab="PC2",col="gray20")
points(subset(pca_res$PC1,pca_res$stage=="Artinskian"),subset(pca_res$PC2,pca_res$stage=="Artinskian"),pch=16,col=paste(subset(pca_res$color,pca_res$stage=="Artinskian")))
text(0.6,0.4,"Artinskian")
text(0.6,0.35,paste(species_ct[which(names(species_ct)=="Artinskian")],"species",sep=" "))
text(0.6,0.3,paste(specimen_ct[which(names(specimen_ct)=="Artinskian")],"specimens",sep=" "))
arrows(0,0,ammon_pca$rotation[1,1]/1.35,ammon_pca$rotation[1,2]/1.35,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[1,1]/1.3,ammon_pca$rotation[1,2]/1.3,rownames(ammon_pca$rotation)[1])
arrows(0,0,ammon_pca$rotation[2,1]/1.35,ammon_pca$rotation[2,2]/1.35,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[2,1]/1.3,ammon_pca$rotation[2,2]/1.3,rownames(ammon_pca$rotation)[2])
arrows(0,0,ammon_pca$rotation[3,1]/1.35,ammon_pca$rotation[3,2]/1.35,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[3,1]/1.3,ammon_pca$rotation[3,2]/1.3,rownames(ammon_pca$rotation)[3])

plot(pca_res$PC1,pca_res$PC2,xlab="PC1",ylab="PC2",col="gray20")
points(subset(pca_res$PC1,pca_res$stage=="Toarcian"),subset(pca_res$PC2,pca_res$stage=="Toarcian"),pch=16,col=paste(subset(pca_res$color,pca_res$stage=="Toarcian")))
text(0.6,0.4,"Toarcian")
text(0.6,0.35,paste(species_ct[which(names(species_ct)=="Toarcian")],"species",sep=" "))
text(0.6,0.3,paste(specimen_ct[which(names(specimen_ct)=="Toarcian")],"specimens",sep=" "))
arrows(0,0,ammon_pca$rotation[1,1]/1.35,ammon_pca$rotation[1,2]/1.35,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[1,1]/1.3,ammon_pca$rotation[1,2]/1.3,rownames(ammon_pca$rotation)[1])
arrows(0,0,ammon_pca$rotation[2,1]/1.35,ammon_pca$rotation[2,2]/1.35,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[2,1]/1.3,ammon_pca$rotation[2,2]/1.3,rownames(ammon_pca$rotation)[2])
arrows(0,0,ammon_pca$rotation[3,1]/1.35,ammon_pca$rotation[3,2]/1.35,col="red",lwd=2,length=0.1)
text(ammon_pca$rotation[3,1]/1.3,ammon_pca$rotation[3,2]/1.3,rownames(ammon_pca$rotation)[3])



#PCA position of species by order
ammon_orders<-c("Goniatitida","Ceratitida","Ammonitida")
ammon_colors<-c("mediumpurple2","dark red","steelblue2")

par(mfrow=c(1,3))
par(mar=c(0.5,0.25,1.5,0.1))

for (i in 1:length(ammon_orders)) {
  order_data<-subset(pca_res,pca_res$order==ammon_orders[i])
  order_data$species<-as.factor(order_data$species)
  order_morph<-data.frame(PC1=sapply(split(order_data$PC1,order_data$species),mean),PC2=sapply(split(order_data$PC2,order_data$species),mean))
  
  plot(pca_res$PC1,pca_res$PC2,pch=16,xlab="",ylab="",xaxt="n",yaxt="n")
  points(order_morph$PC1,order_morph$PC2,col=ammon_colors[i],pch=16)
  mtext(ammon_orders[i],3,adj=0)
}



#variance based on single average shape per species (stage level)
var_sp<-sapply(split(pca_res,pca_res$stage),function(z) sum(apply(sapply(split(z[,3:5],z$species),function(y) apply(y,2,mean)),1,function(x) var(x,na.rm=T))))

sp_ct<-sapply(split(pca_res,pca_res$stage),function(x) length(unique(x$species)))

var_sp_res<-data.frame(interval=levels(pca_res$stage),age=pca_res$stage_ma[match(levels(pca_res$stage),pca_res$stage)],var=var_sp,ct=sp_ct)

var_sp_res<-subset(var_sp_res,var_sp_res$age>0)

var_sp_res<-var_sp_res[order(var_sp_res$age),]

par(mfrow=c(1,1))
par(mar=c(4,4,1,2))

plot(var_sp_res$age,var_sp_res$var,xlab="age (Ma)",ylab="disparity (sum of variances)",xlim=rev(range(var_sp_res$age)),ylim=c(min(var_sp_res$var)-0.005,max(var_sp_res$var)),type="o",pch=16,col="gray")

points(subset(var_sp_res$age,var_sp_res$ct>19),subset(var_sp_res$var,var_sp_res$ct>19),pch=16,cex=1.25,col="red")

rect(periods$early_age,min(var_sp_res$var)-0.008,periods$late_age,min(var_sp_res$var)-0.005,col=paste(periods$color))
text(rowMeans(cbind(periods$early_age,periods$late_age)),min(var_sp_res$var)-0.0065,periods$interval_name,cex=0.75)


#Foote 1993 partial disparity
par.disp<-function(pca_res) { #function to calculate partial disparity
  species_means<-sapply(split(pca_res,pca_res$species),function(x) c(mean(x$PC1),mean(x$PC2),mean(x$PC3)))
  centroid_pt<-apply(species_means,1,function(x) mean(x,na.rm=T))
  apply(species_means,2,function(x) sqrt(sum((x-centroid_pt)^2)))
}

#calculates partial disparity contribution of each species for each stage
distance_res<-sapply(split(pca_res,pca_res$stage_ma),par.disp)

#finds values greater than zero
distance_res<-distance_res[which(apply(distance_res,1,function(x) sum(x,na.rm=T))>0),]

#creates data frame with species name, stage, and partial disparity contribution
distance_frame<-data.frame(species=rep(rownames(distance_res),ncol(distance_res)),stage_ma=sort(as.numeric(rep(colnames(distance_res),nrow(distance_res)))),sp_dist=as.vector(distance_res))
distance_frame<-na.omit(distance_frame)
distance_frame$order<-pca_res$order[match(distance_frame$species,pca_res$species)]
distance_frame$family<-pca_res$family[match(distance_frame$species,pca_res$species)]

#sums partial disparity contribution of each order
part_disp_order<-sapply(split(distance_frame,distance_frame$stage_ma),function(x) sapply(split(x,x$order),function(y) sum(y$sp_dist)/(length(na.omit(x$sp_dist))-1)))
part_disp_order[which(part_disp_order==0)]<-NA

disparity_overall<-sapply(split(distance_frame,distance_frame$stage_ma),function(x) sum(x$sp_dist)/(length(na.omit(x$sp_dist))-1))

#plots overall partial disparity
plot(as.numeric(names(disparity_overall)),disparity_overall,type="o",pch=16,xlim=rev(range(as.numeric(names(disparity_overall)))),ylim=c(0,max(disparity_overall)),xlab="age (Ma)",ylab="partial disparity")

#adds lines for each order
lines(as.numeric(colnames(part_disp_order)),part_disp_order[1,],lty=2,col="red")
lines(as.numeric(colnames(part_disp_order)),part_disp_order[2,],lty=2,col="steelblue2")
lines(as.numeric(colnames(part_disp_order)),part_disp_order[4,],lty=2,col="dark red")
lines(as.numeric(colnames(part_disp_order)),part_disp_order[5,],lty=2,col="orange")
lines(as.numeric(colnames(part_disp_order)),part_disp_order[6,],lty=2,col="mediumpurple2")
lines(as.numeric(colnames(part_disp_order)),part_disp_order[7,],lty=2,col="blue")
lines(as.numeric(colnames(part_disp_order)),part_disp_order[9,],lty=2,col="seagreen3")

points(as.numeric(colnames(part_disp_order)),part_disp_order[1,],pch=16,col="red")
points(as.numeric(colnames(part_disp_order)),part_disp_order[2,],pch=16,col="steelblue2")
points(as.numeric(colnames(part_disp_order)),part_disp_order[4,],pch=16,col="dark red")
points(as.numeric(colnames(part_disp_order)),part_disp_order[5,],pch=16,col="orange")
points(as.numeric(colnames(part_disp_order)),part_disp_order[6,],pch=16,col="mediumpurple2")
points(as.numeric(colnames(part_disp_order)),part_disp_order[7,],pch=16,col="blue")
points(as.numeric(colnames(part_disp_order)),part_disp_order[9,],pch=16,col="seagreen3")

legend(125,0.32,c("Ammonitida","Phylloceratida","Ceratitida","Prolecanitida","Goniatitida","Clymeniida","Agoniatitida"),pch=16,lty=2,col=c("steelblue2","blue","dark red","seagreen3","mediumpurple2","orange","red"),cex=0.9,bty="n")

rect(periods$early_age,-0.01,periods$late_age,0,col=paste(periods$color))
text(rowMeans(cbind(periods$early_age,periods$late_age)),-0.005,periods$interval_name,cex=0.75)



#body size (stage level)
species_size<-sapply(split(pca_res,pca_res$stage_ma),function(x) sapply(split(x$size,x$species),mean)) #calculates mean size per species in each stage

par(mfrow=c(1,1))
par(mar=c(5,4,4,2))
plot(pca_res$stage_ma,pca_res$size,log="y",type="n",ylim=c(min(species_size,na.rm=T),max(species_size,na.rm=T)),xlim=rev(range(pca_res$stage_ma,na.rm=T)),xlab="age (Ma)",ylab="volume")
for (i in 1:ncol(species_size)) {
  points(rep(colnames(species_size)[i],nrow(species_size)),species_size[,i])
}
stage_mean<-apply(species_size,2,function(x) exp(mean(log(x),na.rm=T))) #mean volume of species in each stage
points(names(stage_mean),stage_mean,col="red",pch=15,cex=1.2)


#histogram of goniatite and ammonite shell volumes
par(mfrow=c(1,1))
par(mar=c(5,4,4,2))
hist(log10(subset(pca_res$size,pca_res$order=="Ammonitida")),col=rgb(92,172,238,0.5*255,maxColorValue=255),xlim=c(-1,9),main="",xlab="log(volume)",ylab="specimens")
abline(v=median(log10(subset(pca_res$size,pca_res$order=="Ammonitida"))),lwd=3,col="steelblue2")
hist(log10(subset(pca_res$size,pca_res$order=="Goniatitida")),col=rgb(159,121,238,0.5*255,maxColorValue=255),add=T)
abline(v=median(log10(subset(pca_res$size,pca_res$order=="Goniatitida"))),lwd=3,col="mediumpurple2")
legend(-1,600,c("Ammonitida","Goniatitida"),col=c(rgb(92,172,238,0.5*255,maxColorValue=255),rgb(159,121,238,0.5*255,maxColorValue=255)),pch=15,pt.cex=2,bty="n")



#relationship between shell volume and W/D ratio
par(mfrow=c(1,2))
par(mar=c(4,3,1,0))
plot(pca_res$size,pca_res$WD,log="xy",col="gray20",xlab="log(volume)",ylab="log(W/D)")
points(subset(pca_res$size,pca_res$order=="Goniatitida"),subset(pca_res$WD,pca_res$order=="Goniatitida"),pch=16,col="mediumpurple2")
text(2,1.4,"Goniatitida")
par(mar=c(4,1,1,2))
plot(pca_res$size,pca_res$WD,log="xy",col="gray20",xlab="log(volume)",ylab="",yaxt="n")
points(subset(pca_res$size,pca_res$order=="Ammonitida"),subset(pca_res$WD,pca_res$order=="Ammonitida"),pch=16,col="steelblue2")
text(2,1.4,"Ammonitida")
