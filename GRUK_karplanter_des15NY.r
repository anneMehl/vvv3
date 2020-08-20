## Script for trinnvise analyser av sammenhengen mellom fragmentering
## (isolasjon og patchstørrelse) og artsrikdom, artssammensetning og trekk
## for karplanter i åpen grunnlendt kalkmark
## 2015 12 18

setwd("R/Landskaps-SIS")

## laster biblioteker
library(vegan)
library(MASS)        
library(lme4)
library(nlme)

#Hele datasettet - KARPLANTER
#Inkluderer også alle antatt viktige forklaringsvariabler
community.full <- read.table("plants.full.txt", header=T, dec=".")
# head(community.full)
# str(community.full)
community.full$Subset <- as.factor(community.full$Subset)
community.full$Storrute <- as.factor(community.full$Storrute)
community.full$Polygon <- as.factor(community.full$Polygon)
community.full$PolygonID <- as.factor(community.full$PolygonID)
community.full$Rute<- as.factor(community.full$Rute)
str(community.full)
ncol(community.full)
# fjerner ruta som ble slått
community.full <- community.full[-108,]

# splitter datasettet på arter og miljøvariabler
species.full <- community.full[,20:177]
head(species.full)
nrow(species.full)
# sum(colSums(species.full) == 0) #ingen tomme kolonner, alle arter er til stede 
envir.full <- community.full[,1:19]
head(envir.full)

# kobler artsdata med info om spesialister og invasives
groups <- read.table("plants.groups.txt", header=T, dec=",") #hvilke arter som er spesialister/generalister
names(groups)
colnames(species.full) == groups$Art #Stemmer! 
# definerer spesialister, generalister og invasives
species.generalists <- species.full[groups$Habitat == "generalist"] #generalistene
# ncol(species.generalists)
# species.generalists
species.specialists <- species.full[groups$Habitat == "spesialist"] #spesialistene 
nrow(species.specialists)
# ncol(species.specialists)
# species.specialists
invasives <- species.full[groups$Inv. == "invasive"] #svartelistearter 
# ncol(invasives)
# invasives

################################################################################################
### SPECIES RICHNESS ###
## beregner artsrikdom pr. rute
#Alle arter
total.richness <- rep(0, times=nrow(species.full))
for(i in 1:nrow(species.full)){
  total.richness[i] <- sum(species.full[i,] > 0)
}
#Spesialister
specialist.richness <- rep(0, times=nrow(species.specialists))
for(i in 1:nrow(species.specialists)){
  specialist.richness[i] <- sum(species.specialists[i,] > 0)
}
#Generalister
generalist.richness <- rep(0, times=nrow(species.generalists))
for(i in 1:nrow(species.generalists)){
  generalist.richness[i] <- sum(species.generalists[i,] > 0)
}
#### TOTAL ABUNDANCE ####
## beregner total abundanse pr. rute
#Alle arter
total.abu <- rep(0, times=nrow(species.full))
for(i in 1:nrow(species.full)){
  total.abu[i] <- sum(species.full[i,])
}
#Spesialister
specialist.abu <- rep(0, times=nrow(species.specialists))
for(i in 1:nrow(species.specialists)){
  specialist.abu[i] <- sum(species.specialists[i,])
}
#Generalister
generalist.abu <- rep(0, times=nrow(species.generalists))
for(i in 1:nrow(species.generalists)){
  generalist.abu[i] <- sum(species.generalists[i,])
}
## invasive arter
invasive.abu <- rep(0, times=nrow(invasives))
for(i in 1:nrow(invasives)){
  invasive.abu[i] <- sum(invasives[i,])
}

## slår sammen alle disse til én fil
richness.data <- cbind(total.richness,specialist.richness,generalist.richness,total.abu,specialist.abu,invasive.abu,generalist.abu,envir.full)
str(richness.data)
richness.data$Proportion <- richness.data$specialist.abu/richness.data$total.abu
richness.data$Proportion.gen <- richness.data$generalist.abu/richness.data$total.abu
## beregner gjennomsnittsverdier til Tabell 2
# #artsrikdom spesialister
# min(richness.data$specialist.richness)
# max(richness.data$specialist.richness)
# mean(richness.data$specialist.richness)
# #artsrikdcom generalister
# min(richness.data$generalist.richness)
# max(richness.data$generalist.richness)
# mean(richness.data$generalist.richness)
# #relativ abundanse spesialister
# min(richness.data$Proportion)
# max(richness.data$Proportion)
# mean(richness.data$Proportion)
# #relativ abundanse generalister
# min(richness.data$Proportion.gen)
# max(richness.data$Proportion.gen)
# mean(richness.data$Proportion.gen)


################################################################################################
## miljøvariabler - planteanalysene
summary(lme(log2(Konn)~log2(Areal),random=~1|PolygonID,data=envir.full))
summary(lme(log2(Stone.gravel+1)~log2(Areal),random=~1|PolygonID,data=envir.full))
summary(lme(log2(Busk+1)~log2(Areal),random=~1|PolygonID,data=envir.full))
summary(lme(log2(invasive.abu+1)~log2(Areal),random=~1|PolygonID,data=envir.full))
summary(lme(log2(Stone.gravel+1)~log2(Konn),random=~1|PolygonID,data=envir.full))
summary(lme(log2(Busk+1)~log2(Konn),random=~1|PolygonID,data=envir.full))
summary(lme(log2(invasive.abu+1)~log2(Konn),random=~1|PolygonID,data=envir.full))

## tester effekter av patch size og connectivity på artsrikdom og -abundanse
## model simplification with log-likelihood testing
library(lme4)
library(nlme)
## alle arter, artsrikdom
sr1 <- glmer(total.richness~log2(Areal)*log2(Konn)+(1|PolygonID),data=richness.data,family = "poisson")
sr2 <- update(sr1,~.-log2(Areal):log2(Konn))
anova(sr1,sr2)
sr3 <- update(sr2,~.-log2(Konn))
anova(sr2,sr3)
sr4 <- update(sr3,~.-log2(Areal))
anova(sr3,sr4)
summary(sr1)
# ingenting signifikant

## spesialister, artsrikdom
ssr1 <- glmer(specialist.richness~log2(Areal)*log2(Konn)+(1|PolygonID),data=richness.data,family = "poisson")
 ssr2 <- update(ssr1,~.-log2(Areal):log2(Konn))
 anova(ssr1,ssr2)
summary(ssr1)
plot(ssr1)#ser ok ut


###############################################################################################
## GNMDS ##
library(MASS)
library(vegan)
# ## Kjører GNMDS på nytt, slik at vi kan lagre envfit-resultatene 
# #Alle arter
# dist.mat.species.full <- vegdist(species.full,method="bray")
# mds<-NULL
# for(i in 1:100){mds[[i]] <- isoMDS(dist.mat.species.full, initMDS(dist.mat.species.full), k=2,maxit=200, tol=1e-7)} #k=antall dimensjoner
# mds.stress<-unlist(lapply(mds,function(v){v[[2]]})) 
# mds.stress
# ordered <- order(mds.stress)
# mds.stress[ordered[1]]
# mds.stress[ordered[2]]
# mds.best<-postMDS(mds[[ordered[1]]],dist.mat.species.full,halfchange=T) #Skalering til HC-units
# mds.best
# mds.secbest<-postMDS(mds[[ordered[2]]],dist.mat.species.full,halfchange=T)
# mds.secbest
# #Procrustes comparisons
# procrustes(mds.best,mds.secbest,permutations=999)
# protest(mds.best,mds.secbest,permutations=999) #Hvis signifikant, er de to løsningene like og vi har funnet den beste løsningen
# #Plot
# gnmds1<-mds.best$points[,1]
# gnmds2<-mds.best$points[,2]
# plot(mds.best$points[,1],mds.best$points[,2],xlab="GNMDS1",ylab="GNMDS2",type="n") #ylim=c(-0.3,0.3), xlim=c(-0.4,0.7))
# points(gnmds1,gnmds2,cex=1.5) #svart, firkant
# species.scores1 <- wascores(gnmds1, species.full)
# species.scores2 <- wascores(gnmds2, species.full)
# text(species.scores1[groups$Habitat == "generalist"],species.scores2[groups$Habitat == "generalist"], labels=groups$Art[groups$Habitat == "generalist"], col="black")
# text(species.scores1[groups$Habitat == "spesialist"],species.scores2[groups$Habitat == "spesialist"], labels=groups$Art[groups$Habitat == "spesialist"], col="red")
# abline(h=0, lty=3)
# abline(v=0, lty=3)
# gnmds.full <- mds.best
# save(gnmds.full,file="PLANTE.GNMDS.RData")
# 
# gmds1.data <- as.data.frame(gnmds1)
# gmds2.data <- as.data.frame(gnmds2)
# y <- cbind(gmds1.data,gmds2.data,invasive.abu)
# 
### lager envfit-funksjoner av hver av ordinasjonene
# env.plants <- envfit(gnmds.full~log2(Stone.gravel+1)+log2(Busk+1)+log2(invasive.abu+1)+log2(Areal)+log2(Konn),permutations=999,data=envir.full)
# env.plants
# save(env.plants, file="ENV.PLANTS.RData")

# #Utregning av gjennomsnittlig aksescore (til senere) - finnes i egen fil, coordinates.txt
# uID <- unique(envir.full$PolygonID)
# #Akse1
# mean.gnmds1 <- rep(0, times=length(uID))
# for(i in 1:length(uID)) {
#   mean.gnmds1[i] <- mean(gnmds1[envir.full$PolygonID == (uID[i])], na.rm=T)
# }
# mean.gnmds1
# mean.gnmds1.data <- as.data.frame(mean.gnmds1)
# mean.gnmds1.data <- cbind(mean.gnmds1.data, uID)
# #Akse 2
# mean.gnmds2 <- rep(0, times=length(uID))
# for(i in 1:length(uID)) {
#   mean.gnmds2[i] <- mean(gnmds2[envir.full$PolygonID == (uID[i])], na.rm=T)
# }
# mean.gnmds2
# mean.gnmds2.data <- as.data.frame(mean.gnmds2)
# mean.gnmds2.data <- cbind(mean.gnmds2.data, uID)
# mean.gnmds1.data
# mean.gnmds2.data
# 

#Spesialister
str(species.specialists)
# fjerner ruter uten arter
rowSums(species.specialists)==0
species.specialists.r <- species.specialists[-61,]
#To ruter er helt identiske: frag.spp = 1 - endrer den ene til frag spp. = 2 så jeg får kjørt analysen: 
species.specialists.r[130,]
species.specialists.r[152,]
species.specialists.r[152,12] <- 2

# dist.mat.species.specialists <- vegdist(species.specialists.r,method="bray")
# mds.s<-NULL
# for(i in 1:100){mds.s[[i]] <- isoMDS(dist.mat.species.specialists, initMDS(dist.mat.species.specialists), k=2,maxit=200, tol=1e-7)} #k=antall dimensjoner
# mds.stress.s<-unlist(lapply(mds.s,function(v){v[[2]]})) 
# mds.stress.s
# ordered.s <- order(mds.stress.s)
# mds.stress.s[ordered.s[1]]
# mds.stress.s[ordered.s[2]]
# mds.best.s<-postMDS(mds.s[[ordered.s[1]]],dist.mat.species.specialists,halfchange=T) #Skalering til HC-units
# mds.best.s
# mds.secbest.s<-postMDS(mds.s[[ordered.s[2]]],dist.mat.species.specialists,halfchange=T)
# mds.secbest.s
# 
# #Procrustes comparisons
# procrustes(mds.best.s,mds.secbest.s,permutations=999)
# protest(mds.best.s,mds.secbest.s,permutations=999) #Hvis signifikant, er de to løsningene like og vi har funnet den beste løsningen
# 
# #Plot
# gnmds1.s<-mds.best.s$points[,1]
# gnmds2.s<-mds.best.s$points[,2]
# plot(mds.best.s$points[,1],mds.best.s$points[,2],xlab="GNMDS1",ylab="GNMDS2",type="n") # ylim=c(-0.6,0.6), xlim=c(-0.7,0.7))
# points(gnmds1.s,gnmds2.s,cex=1.5) #svart, firkant
# species.scores1.s <- wascores(gnmds1.s, species.specialists.r)
# species.scores2.s <- wascores(gnmds2.s, species.specialists.r)
# text(species.scores1.s,species.scores2.s, labels=rownames(species.scores1.s), col="red")
# abline(h=0, lty=3)
# abline(v=0, lty=3)
# 
# gnmds.spes <- mds.best.s
# save(gnmds.spes,file="PLANTE.GNMDS.SPES.RData")
# 
# gmds1s.data <- as.data.frame(gnmds1.s)
# gmds2s.data <- as.data.frame(gnmds2.s)
# y <- cbind(gmds1s.data,gmds2s.data)
# y

# #Utregning av gjennomsnittlig aksescore (til senere) - finnes i egen fil, coordinates.txt
# uID <- unique(envir.full$PolygonID)
# #Akse1
# mean.gnmds1s <- rep(0, times=length(uID))
# for(i in 1:length(uID)) {
#   mean.gnmds1s[i] <- mean(gnmds1.s[envir.full$PolygonID == (uID[i])], na.rm=T)
# }
# mean.gnmds1s
# mean.gnmds1s.data <- as.data.frame(mean.gnmds1s)
# #Akse 2
# mean.gnmds2s <- rep(0, times=length(uID))
# for(i in 1:length(uID)) {
#   mean.gnmds2s[i] <- mean(gnmds2.s[envir.full$PolygonID == (uID[i])], na.rm=T)
# }
# mean.gnmds2s
# mean.gnmds2s.data <- as.data.frame(mean.gnmds2s.data)
# mean.gnmds2s.data <- cbind(mean.gnmds1s.data,mean.gnmds2s.data, uID)
# mean.gnmds2s.data


## har lastet over alle koordinater - både polygon-means og ruteverdier - for
## alle arter, generalister og spesialister i egne filer
coord <- read.delim("coord.plants.txt",dec=".")
str(coord)
coord$PolygonID <- as.factor(coord$PolygonID)
coord$Storrute <- as.factor(coord$Storrute)

## laster inn gnmds-resultater og envfit-resultater
load("PLANTE.GNMDS.RData")
load("PLANTE.GNMDS.SPES.RData")
load("PLANTE.GNMDS.META.RData")
load("ENV.PLANTS.RData")
load("ENV.PLANTS.SPES.RData")
load("ENV.SCALED.PLANTS.RData")
## plotter GNMDS-resultatene - ALLE ARTER

plot(coord$gnmds1,coord$gnmds2,xlab="GNMDS1",ylab="GNMDS2",type="n") #ylim=c(-0.3,0.3), xlim=c(-0.4,0.7))
points(coord$gnmds1,coord$gnmds2,pch=16,cex=0.6,col=8) #liten, grå prikk
species.scores1 <- wascores(coord$gnmds1, species.full)
species.scores2 <- wascores(coord$gnmds2, species.full)
text(species.scores1[groups$Habitat == "generalist"],species.scores2[groups$Habitat == "generalist"], labels=groups$Art[groups$Habitat == "generalist"], col="black")
text(species.scores1[groups$Habitat == "spesialist"],species.scores2[groups$Habitat == "spesialist"], labels=groups$Art[groups$Habitat == "spesialist"], col="red")
abline(h=0, lty=3)
abline(v=0, lty=3)
plot(env.plants)

## et enklere plot
## SIRI: kan vi gjøre dette bedre?
plot(coord$gnmds1,coord$gnmds2,xlab="GNMDS1",ylab="GNMDS2",type="n") #ylim=c(-0.3,0.3), xlim=c(-0.4,0.7))
#points(coord$gnmds1,coord$gnmds2,pch=16,cex=0.6,col=8) #liten, grå prikk
abline(h=0, lty=3)
abline(v=0, lty=3)
# species.scores1 <- wascores(coord$gnmds1, species.full)
# species.scores2 <- wascores(coord$gnmds2, species.full)
idx <- colSums(species.full)>100
#text(species.scores1[idx],species.scores2[idx], labels=row.names(species.scores1)[idx], col="black")
text(species.scores1[idx&groups$Habitat=="generalist"],species.scores2[idx&groups$Habitat=="generalist"], labels=row.names(species.scores1)[idx&groups$Habitat=="generalist"], col="black")
text(species.scores1[idx&groups$Habitat=="spesialist"],species.scores2[idx&groups$Habitat=="spesialist"], labels=row.names(species.scores1)[idx&groups$Habitat=="spesialist"], col="red")
plot(env.plants,col="blue",add=T,cex=0.75, labels=c("Cover of gravel", "Shrub cover", "Cover of invasives", "Patch size", "Patch connectivity"))



#korrelasjonstester mellom akseskår og miljøvariabler
#akse 1
cor.test(coord$gnmds1,log2(coord$Areal),method="k")
cor.test(coord$gnmds1,log2(coord$Konn),method="k")
cor.test(coord$gnmds1,log2(coord$Busk+1),method="k")
cor.test(coord$gnmds1,log2(coord$Stone.gravel+1),method="k")
cor.test(coord$gnmds1,log2(coord$invasive.abu+1),method="k")
#akse 2
cor.test(coord$gnmds2,log2(coord$Areal),method="k")
cor.test(coord$gnmds2,log2(coord$Konn),method="k")
cor.test(coord$gnmds2,log2(coord$Busk+1),method="k")
cor.test(coord$gnmds2,log2(coord$Stone.gravel+1),method="k")
cor.test(coord$gnmds2,log2(coord$invasive.abu+1),method="k")


## plotter GNMDS-resultatene - HABITATSPESIALISTENE
# env.plants.spes <- envfit(gnmds.spes~log2(Stone.gravel+1)+log2(Busk+1)+log2(invasive.abu+1)+log2(Areal)*log2(Konn),data=coord[-61,])
# save(env.plants.spes, file="ENV.PLANTS.SPES.RData")
load("ENV.PLANTS.SPES.RData")

plot(coord$gnmds1.s,coord$gnmds2.s,xlab="GNMDS1",ylab="GNMDS2",type="n") #ylim=c(-0.3,0.3), xlim=c(-0.4,0.7))
points(coord$gnmds1.s,coord$gnmds2.s,pch=16,cex=0.6,col=8) #liten, grå prikk
species.scores1s <- wascores(coord$gnmds1.s, species.specialists)
species.scores2s <- wascores(coord$gnmds2.s, species.specialists)
text(species.scores1s,species.scores2s, labels=rownames(species.scores1s), col="red")
abline(h=0, lty=3)
abline(v=0, lty=3)
plot(env.plants.spes)

## ryddigere
## SIRI: må det fikses, eller er det ok?
plot(coord$gnmds1.s,coord$gnmds2.s,xlab="GNMDS1",ylab="GNMDS2",type="n",xlim=c(-1.5,1.5)) #ylim=c(-0.3,0.3), xlim=c(-0.4,0.7))
#points(coord$gnmds1.s,coord$gnmds2.s,pch=16,cex=0.6,col=8) #liten, grå prikk
abline(h=0, lty=3)
abline(v=0, lty=3)
idx <- colSums(species.specialists) > 40 & !colSums(species.specialists) == 334
text(species.scores1s[idx],species.scores2s[idx], labels= rownames(species.scores1s)[idx], col="red")
plot(env.plants.spes,col="blue",add=T,cex=0.9, labels=c("Cover of gravel", "Shrub cover", "Cover of invasives", "Patch size", "Patch connectivity"))


#korrelasjonstester mellom akseskår og miljøvariabler
#akse 1
cor.test(coord$gnmds1.s,log2(coord$Areal),method="k")
cor.test(coord$gnmds1.s,log2(coord$Konn),method="k")
cor.test(coord$gnmds1.s,log2(coord$Busk+1),method="k")
cor.test(coord$gnmds1.s,log2(coord$Stone.gravel+1),method="k")
cor.test(coord$gnmds1.s,log2(coord$invasive.abu+1),method="k")
#akse 2
cor.test(coord$gnmds2.s,log2(coord$Areal),method="k")
cor.test(coord$gnmds2.s,log2(coord$Konn),method="k")
cor.test(coord$gnmds2.s,log2(coord$Busk+1),method="k")
cor.test(coord$gnmds2.s,log2(coord$Stone.gravel+1),method="k")
cor.test(coord$gnmds2.s,log2(coord$invasive.abu+1),method="k")

## korrelasjonstest mellom total artssammensetning og spesialistsammensetning
cor.test(coord$gnmds1,coord$gnmds1.s,method="k")
cor.test(coord$gnmds1,coord$gnmds2.s,method="k")
cor.test(coord$gnmds2,coord$gnmds1.s,method="k")
cor.test(coord$gnmds2,coord$gnmds2.s,method="k")



##################################################################################################################################
############# BETINGET ORDINASJON - HVOR MYE FORKLARER PATCH SIZE AND ISOLATION AV VARIASJON I ART-RUTE-MATRISEN? ################
## constrained ordination - species-plot matrix by patch size and isolation  

## constrained ordination - habitat specialist species-plot matrix by patch size and isolation  
## tar ut ruten uten spesialistarter
m1s <- cca(species.specialists.r~log2(Areal),data=envir.full[-61,])
m1s # 2,5 % av variasjonen i artssammensetning forklart av areal
permutest(m1s,permutations = 999)
# det er signifikant. 
m2s <- cca(species.specialists.r~log2(Konn),data=envir.full[-61,])
m2s # 3,7 % av variasjonen i artssammensetning forklart av isolasjon
permutest(m2s,permutations = 999)# sign

#bidrar areal utover konnektivitet?
m3s <- cca(species.specialists.r~log2(Areal)+Condition(log2(Konn)),data=envir.full[-61,])
permutest(m3s,permutations = 999)
m3s
#JA, 1,7 % ekstra
#bidrar interaksjonen mellom areal og konnektivitet?
m4s <- cca(species.specialists.r~log2(Konn):log2(Areal)+Condition(log2(Areal)+log2(Konn)),data=envir.full[-61,])
m4s
permutest(m4s,permutations = 999)
# TJA, 0,7 % ekstra, p er ca. 0,05

m5s <- cca(species.specialists.r~log2(Konn)*log2(Areal),data=envir.full[-61,])
m5s
permutest(m5s,permutations = 999)

#Forklarer til sammen 6,2 % av variasjonen i art-rute-matrisen
plot(m5s,col=3,cex=0.75, display=c("sp","bp"))

###############################################################################################
#### VARIATION IN TRAIT COMPOSITION
## traits-plot-matrisen for habitatspesialistene
traits.full <- read.table("plants.traits.txt", header=T, dec=".") 

traits.full$Vit.navn == names(species.full) #stemmer 
dim(traits.full)
head(traits.full)
identical(as.character(traits.full$Vit.navn), names(species.full))  #OK

## calculating min, max and mean values for each trait - for generalists and specialists
## for table 2
table(traits.full$no.clonality[groups$Habitat=="generalist"])
table(traits.full$no.clonality[groups$Habitat=="spesialist"])
table(traits.full$short.clonality[groups$Habitat=="generalist"])
table(traits.full$short.clonality[groups$Habitat=="spesialist"])
table(traits.full$far.clonality[groups$Habitat=="generalist"])
table(traits.full$far.clonality[groups$Habitat=="spesialist"])

summary(traits.full$no.clonality[groups$Habitat=="generalist"])
summary(traits.full$short.clonality[groups$Habitat=="generalist"])
summary(traits.full$far.clonality[groups$Habitat=="generalist"])

traits.specialists <- traits.full[groups$Habitat == "spesialist",]
dim(traits.specialists)

weighted.specialists <- matrix(nrow=nrow(species.specialists), ncol=13) 
colnames(weighted.specialists) <- colnames(traits.full[c(7:8, 14, 16:17, 18:19, 21, 24:25, 30:32)]) 
colnames(weighted.specialists)
dim(weighted.specialists)
head(weighted.specialists)

for (i in 1:nrow(species.specialists)){
  
  #weighted.specialists[i,1] <- weighted.mean(traits.specialists$herb, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,2] <- weighted.mean(traits.specialists$gram, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,3] <- weighted.mean(traits.specialists$woody, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,1] <- weighted.mean(traits.specialists$annual, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,2] <- weighted.mean(traits.specialists$biennial, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,3] <- weighted.mean(traits.specialists$perennial, species.specialists[i,], na.rm=T)
  weighted.specialists[i,1] <- weighted.mean(traits.specialists$shortlived, species.specialists[i,], na.rm=T)
  weighted.specialists[i,2] <- weighted.mean(traits.specialists$longlived, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,7] <- weighted.mean(traits.specialists$monocarpic, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,8] <- weighted.mean(traits.specialists$polycarpic, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,9] <- weighted.mean(traits.specialists$rosette, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,10] <- weighted.mean(traits.specialists$half.rosette, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,11] <- weighted.mean(traits.specialists$no.rosette, species.specialists[i,], na.rm=T)
  weighted.specialists[i,3] <- weighted.mean(traits.specialists$no.clonality, species.specialists[i,], na.rm=T)
  weighted.specialists[i,4] <- weighted.mean(traits.specialists$short.clonality, species.specialists[i,], na.rm=T)
  weighted.specialists[i,5] <- weighted.mean(traits.specialists$far.clonality, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,4] <- weighted.mean(traits.specialists$clonality, species.specialists[i,], na.rm=T)
  weighted.specialists[i,6] <- weighted.mean(traits.specialists$canopy.height, species.specialists[i,], na.rm=T)
  weighted.specialists[i,7] <- weighted.mean(traits.specialists$SLA, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,17] <- weighted.mean(traits.specialists$LDMC, species.specialists[i,], na.rm=T)
  weighted.specialists[i,8] <- weighted.mean(traits.specialists$Leaf.size, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,19] <- weighted.mean(traits.specialists$leaf.mass, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,9] <- weighted.mean(traits.specialists$sl.index, species.specialists[i,], na.rm=T)
  weighted.specialists[i,9] <- weighted.mean(traits.specialists$seed.mass, species.specialists[i,], na.rm=T)
  weighted.specialists[i,10] <- weighted.mean(traits.specialists$seed.number, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,12] <- weighted.mean(traits.specialists$terminal.velocity, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,13] <- weighted.mean(traits.specialists$insect, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,14] <- weighted.mean(traits.specialists$windpoll, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,15] <- weighted.mean(traits.specialists$selfing, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,17] <- weighted.mean(traits.specialists$unassisted, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,18] <- weighted.mean(traits.specialists$wind, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,19] <- weighted.mean(traits.specialists$endozoochory, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,20] <- weighted.mean(traits.specialists$exozoochory, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,21] <- weighted.mean(traits.specialists$water, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,22] <- weighted.mean(traits.specialists$ants, species.specialists[i,], na.rm=T)
  #weighted.specialists[i,23] <- weighted.mean(traits.specialists$ballistically, species.specialists[i,], na.rm=T)
  weighted.specialists[i,11] <- weighted.mean(traits.specialists$unassist.disp, species.specialists[i,], na.rm=T)
  weighted.specialists[i,12] <- weighted.mean(traits.specialists$wind.disp, species.specialists[i,], na.rm=T)
  weighted.specialists[i,13] <- weighted.mean(traits.specialists$animal.disp, species.specialists[i,], na.rm=T)
  #weighted.means[i,34] <- weighted.mean(traits.full$light, species.full[i,], na.rm=T)
  #weighted.means[i,35] <- weighted.mean(traits.full$moisture, species.full[i,], na.rm=T)
  #weighted.means[i,36] <- weighted.mean(traits.full$pH, species.full[i,], na.rm=T)
  #weighted.means[i,37] <- weighted.mean(traits.full$nitrogen, species.full[i,], na.rm=T)
  #weighted.means[i,38] <- weighted.mean(traits.full$hab.spes, species.full[i,], na.rm=T)
  #weighted.means[i,39] <- weighted.mean(traits.full$kalkkrevende, species.full[i,], na.rm=T)  
}

head(weighted.specialists)
weighted.specialists.scaled <- scale(weighted.specialists) #Sentrerer og skalerer
summary(weighted.specialists)
# tar ut rute 61, som ikke har noen spesialister
weighted.specialists <- weighted.specialists[-61,]
weighted.specialists.scaled <- scale(weighted.specialists) #Sentrerer og skalerer

t1 <- rda(weighted.specialists.scaled~log2(Areal),data=envir.full[-61,])
t1 # 0,9 % av variasjonen i trekk-rute-matrisen blir forklart av areal
permutest(t1,permutations = 999)#ikke-sign

t2 <- rda(weighted.specialists.scaled~log2(Konn),data=envir.full[-61,])
t2 # 2,0 % av variasjonen i trekk-rute-matrisen blir forklart av areal
permutest(t2,permutations = 999)#sign
plot(t2,col=3,cex=0.75, display=c("sp","bp"))
plot(t2)

## forklarer areal noe i tillegg til konnektivitet?
t3 <- rda(weighted.specialists.scaled~log2(Areal)+Condition(log2(Konn)),data=envir.full[-61,])
permutest(t3,permutations = 999)#ja
t3 # 1,3 % ekstra

## forklarer interaksjonen noe i tillegg til konnektivitet og areal?
t4 <- rda(weighted.specialists.scaled~log2(Konn):log2(Areal)+Condition(log2(Konn)+log2(Areal)),data=envir.full[-61,])
permutest(t4,permutations = 999)#ja
t4 # 2,3 % ekstra

t5 <- rda(weighted.specialists.scaled~log2(Konn)*log2(Areal),data=envir.full[-61,])
permutest(t5,permutations = 999)#ja
t5 # til sammen 5,7 %


## SIRI: må fikses slik at vi har "ordentlige" pilnavn, jf. GNMDS-figurene. 
rda.t51 <- scores(t5,display="species",origin=TRUE)[,1]
rda.t52 <- scores(t5,display="species",origin=TRUE)[,2]
plot(rda.t51,rda.t52,type="n",xlab="RDA1", ylab="RDA2",xlim = c(-1,1))
text(rda.t51,rda.t52,labels=names(rda.t51),cex=0.6)
abline(h=0,lty=3)
abline(v=0,lty=3)


