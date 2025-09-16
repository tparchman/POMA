library("dplyr")
library("stringr")

# reading in the geno probabilities for k=2
# k4 best supported by DIC but k2 takes slots 2-4 for next best
datk2 <- read.csv("gprob2.txt")
ids <- read.csv("POMA_PopID.csv") #to bind so we know indivs
head(datk2[,c(1:6)])

k2.id <- cbind(ids,datk2)

all_pca <- prcomp(k2.id[,c(5:33801)], center=TRUE, scale=FALSE)
summary(all_pca)

head(pca.ids[,c(1:6)])
pca.ids <- cbind(ids,all_pca$x)
summary(all_pca)

#replaces the MV with WV in the Pop IDs
pca.ids <- pca.ids %>%
  mutate(Pop = str_replace(Pop, "^MV", "WV"))

#replaces the MV with WV in the individual IDs
pca.ids <- pca.ids %>%
  mutate(ID = str_replace(ID, "^PM_MV", "PM_WV"))

unique(pca.ids$Pop)

#subsetting the pops
hp <- pca.ids[pca.ids$Pop == "HP",]
ja <- pca.ids[pca.ids$Pop == "JA",]
lb <- pca.ids[pca.ids$Pop == "LB",]
me <- pca.ids[pca.ids$Pop == "ME",]
mt <- pca.ids[pca.ids$Pop == "MT",]
wv <- pca.ids[pca.ids$Pop == "WV",]
pp <- pca.ids[pca.ids$Pop == "PP",]

par(mar=c(5,5,1,1))
plot(all_pca$x[,1], all_pca$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)

points(hp$PC1,hp$PC2,pch=20,col="#3E4A1C",cex=1.5)
points(ja$PC1,ja$PC2,pch=20,col="#3F2438",cex=1.5)
points(lb$PC1,lb$PC2,pch=20,col="#F7DEFF",cex=1.5)
points(me$PC1,me$PC2,pch=20,col="#D4C588",cex=1.5) 
points(mt$PC1,mt$PC2,pch=20,col="#E4D7F5",cex=1.5) 
points(wv$PC1,wv$PC2,pch=20,col="#787737",cex=1.5) 
points(pp$PC1,pp$PC2,pch=20,col="#D8B37D",cex=1.5) 


####### K4 just to compare ########
datk4 <- read.csv("gprob4.txt")
k4.id <- cbind(ids,datk4)

all_pca4 <- prcomp(k4.id[,c(5:33801)], center=TRUE, scale=FALSE)
summary(all_pca4)
summary(all_pca) # to compare

pca.ids4 <- cbind(ids,all_pca4$x)


#replaces the MV with WV in the Pop IDs
pca.ids4 <- pca.ids4 %>%
  mutate(Pop = str_replace(Pop, "^MV", "WV"))

#replaces the MV with WV in the individual IDs
pca.ids4 <- pca.ids4 %>%
  mutate(ID = str_replace(ID, "^PM_MV", "PM_WV"))

unique(pca.ids4$Pop)

#subsetting the pops
hp4 <- pca.ids4[pca.ids4$Pop == "HP",]
ja4 <- pca.ids4[pca.ids4$Pop == "JA",]
lb4 <- pca.ids4[pca.ids4$Pop == "LB",]
me4 <- pca.ids4[pca.ids4$Pop == "ME",]
mt4 <- pca.ids4[pca.ids4$Pop == "MT",]
wv4 <- pca.ids4[pca.ids4$Pop == "WV",]
pp4 <- pca.ids4[pca.ids4$Pop == "PP",]

par(mar=c(5,5,1,1))
plot(all_pca4$x[,1], all_pca4$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)

points(hp4$PC1,hp4$PC2,pch=20,col="#3E4A1C",cex=1.5) #OR
points(ja4$PC1,ja4$PC2,pch=20,col="#3F2438",cex=1.5) #WA
points(lb4$PC1,lb4$PC2,pch=20,col="#F7DEFF",cex=1.5) #CA
points(me4$PC1,me4$PC2,pch=20,col="#D4C588",cex=1.5) #WA
points(mt4$PC1,mt4$PC2,pch=20,col="#E4D7F5",cex=1.5) #CA
points(wv4$PC1,wv4$PC2,pch=20,col="#787737",cex=1.5) #OR
points(pp4$PC1,pp4$PC2,pch=20,col="#D8B37D",cex=1.5) #WA

################
### just WA PCA ###
#### k=2 ######
###############
ja.gen2 <- subset(k2.id, Pop=="JA")
me.gen2 <- subset(k2.id, Pop=="ME")
pp.gen2 <- subset(k2.id, Pop=="PP")

#combining
step <- rbind(ja.gen2,me.gen2)
wa.genos2 <- rbind(step, pp.gen2)
dim(wa.genos2) # should be 32 33801

head(wa.genos2[,c(1:6)])
wa_pca2 <- prcomp(wa.genos2[,c(5:33801)], center=TRUE, scale=FALSE)
summary(wa_pca2)

wa.pca.ids2 <- cbind(wa.genos2[,c(1,3)],wa_pca2$x)
head(wa.pca.ids2[,c(1:6)])

#subsetting from PC scores
ja.wa2 <- wa.pca.ids2[wa.pca.ids2$Pop == "JA",]
me.wa2 <- wa.pca.ids2[wa.pca.ids2$Pop == "ME",]
pp.wa2 <- wa.pca.ids2[wa.pca.ids2$Pop == "PP",]

par(mar=c(5,5,1,1))
plot(wa_pca2$x[,1], wa_pca2$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)

points(ja.wa2$PC1,ja.wa2$PC2,pch=20,col="#3F2438",cex=1.5) #WA
points(me.wa2$PC1,me.wa2$PC2,pch=20,col="#D4C588",cex=1.5) #WA
points(pp.wa2$PC1,pp.wa2$PC2,pch=20,col="#D8B37D",cex=1.5) #WA


#################
### just OR PCA ###
###   k=2 ####
###################

k2.id <- k2.id %>%
  mutate(Pop = str_replace(Pop, "^MV", "WV"))

#replaces the MV with WV in the individual IDs
k2.id <- k2.id %>%
  mutate(ID = str_replace(ID, "^PM_MV", "PM_WV"))

hp.gen2 <- subset(k2.id, Pop=="HP")
wv.gen2 <- subset(k2.id, Pop=="WV")

#combining
or.genos2 <- rbind(hp.gen2, wv.gen2)
dim(or.genos2) # should be 25 33801

#PCA
or_pca2 <- prcomp(or.genos2[,c(5:33801)], center=TRUE, scale=FALSE)
summary(or_pca2)

or.pca.ids2 <- cbind(or.genos2[,c(1,3)],or_pca2$x)
head(or.pca.ids2[,c(1:6)])

#subsetting from PC scores
hp.or2 <- or.pca.ids2[or.pca.ids2$Pop == "HP",]
wv.or2 <- or.pca.ids2[or.pca.ids2$Pop == "WV",]

par(mar=c(5,5,1,1))
plot(or_pca2$x[,1], or_pca2$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)

points(hp.or2$PC1,hp.or2$PC2,pch=20,col="#3E4A1C",cex=1.5)
points(wv.or2$PC1,wv.or2$PC2,pch=20,col="#787737",cex=1.5)


#################
#### CA PCA
##### K - 2
#################

lb.gen2 <- subset(k2.id, Pop=="LB")
mt.gen2 <- subset(k2.id, Pop=="MT")

#combining
ca.genos2 <- rbind(lb.gen2, mt.gen2)
dim(ca.genos2) # should be 32 33801


#PCA
ca_pca2 <- prcomp(ca.genos2[,c(5:33801)], center=TRUE, scale=FALSE)
summary(ca_pca2)

ca.pca.ids2 <- cbind(ca.genos2[,c(1,3)],ca_pca2$x)
head(ca.pca.ids2[,c(1:6)])

#subsetting from PC scores
lb.ca2 <- ca.pca.ids2[ca.pca.ids2$Pop == "LB",]
mt.ca2 <- ca.pca.ids2[ca.pca.ids2$Pop == "MT",]

par(mar=c(5,5,1,1))
plot(ca_pca2$x[,1], ca_pca2$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)

points(lb.ca2$PC1,lb.ca2$PC2,pch=20,col="#F7DEFF",cex=1.5)
points(mt.ca2$PC1,mt.ca2$PC2,pch=20,col="#E4D7F5",cex=1.5)


############### 
#### CA + OR for fun
###############
south.genos2 <- rbind(ca.genos2,or.genos2)
dim(south.genos2)

so_pca2 <- prcomp(south.genos2[,c(5:33801)], center=TRUE, scale=FALSE)
summary(so_pca2) #PC1: 0.1621, PC2: 0.8233

so.pca.ids2 <- cbind(south.genos2[,c(1,3)],so_pca2$x)
head(so.pca.ids2[,c(1:6)])

lb.so2 <- so.pca.ids2[so.pca.ids2$Pop == "LB",]
mt.so2 <- so.pca.ids2[so.pca.ids2$Pop == "MT",]
hp.so2 <- so.pca.ids2[so.pca.ids2$Pop == "HP",]
wv.so2 <- so.pca.ids2[so.pca.ids2$Pop == "WV",]

dev.off()
plot.new()
par(mar=c(5,5,1,1))
plot(so_pca2$x[,1], so_pca2$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)

points(lb.so2$PC1,lb.so2$PC2,pch=20,col="#F7DEFF",cex=1.5) #ca
points(mt.so2$PC1,mt.so2$PC2,pch=20,col="#D4C588",cex=1.5) #ca
points(hp.so2$PC1,hp.so2$PC2,pch=20,col="#3E4A1C",cex=1.5) #or
points(wv.so2$PC1,wv.so2$PC2,pch=20,col="#787737",cex=1.5) #or

####################################################
####################################################
###########putting em all together##################
####################################################
####################################################
layout <- layout(matrix(1:4, ncol=2))

#all at k=2
par(mar=c(5,5,1,1))
plot(all_pca$x[,1], all_pca$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)
points(hp$PC1,hp$PC2,pch=20,col="#3E4A1C",cex=1.5) #dark green
points(ja$PC1,ja$PC2,pch=20,col="#3F2438",cex=1.5) #park purple
points(lb$PC1,lb$PC2,pch=20,col="#F7DEFF",cex=1.5) #pink
points(me$PC1,me$PC2,pch=20,col="#E4DF75",cex=1.5) #light purple
points(mt$PC1,mt$PC2,pch=20,col="#D4C588",cex=1.5) #yellow beige
points(wv$PC1,wv$PC2,pch=20,col="#787737",cex=1.5) #light green
points(pp$PC1,pp$PC2,pch=20,col="#D8B37D",cex=1.5)  #orange beige

#OR
par(mar=c(5,5,1,1))
plot(or_pca2$x[,1], or_pca2$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)
points(hp.or2$PC1,hp.or2$PC2,pch=20,col="#3E4A1C",cex=1.5)
points(wv.or2$PC1,wv.or2$PC2,pch=20,col="#787737",cex=1.5)


#WA
par(mar=c(5,5,1,1))
plot(wa_pca2$x[,1], wa_pca2$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)
points(ja.wa2$PC1,ja.wa2$PC2,pch=20,col="#3F2438",cex=1.5) #WA
points(me.wa2$PC1,me.wa2$PC2,pch=20,col="#E4DF75",cex=1.5) #WA
points(pp.wa2$PC1,pp.wa2$PC2,pch=20,col="#D8B37D",cex=1.5) #WA

#CA
par(mar=c(5,5,1,1))
plot(ca_pca2$x[,1], ca_pca2$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1, ylim=c(-50,50), xlim=c(-50,50))
points(lb.ca2$PC1,lb.ca2$PC2,pch=20,col="#F7DEFF",cex=1.5)
points(mt.ca2$PC1,mt.ca2$PC2,pch=20,col="#D4C588",cex=1.5)

summary(all_pca)# PC1: 0.139  PC2: 0.0919
summary(wa_pca2)# PC1: 0.2351 PC2: 0.04818
summary(or_pca2)# PC1: 0.2206 PC2: 0.0400
summary(ca_pca2)# PC1: 0.1129 PC2: 0.04932


#####
#### fuckin around
#####
dev.off()
plot.new()

str(wa_pca2)
