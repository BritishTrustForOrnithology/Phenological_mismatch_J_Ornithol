############################################################################
# Produce the results for the chapter:                                     #
# Results/Temporal distribution of survey visits                           #
############################################################################

# Define the study periods (eyears=period 1, lyears=period 2)
eyears <- 1994:1998
lyears <- 2013:2017

# Specify study region
region <- "rda7"

# Specify format of figures ("screen" to show figures on screen)
out <- "screen"

setwd("N:/Phenology and trends")

fyear <- eyears[1]
lyear <- lyears[length(lyears)]
fyr <- substr(eyears[1],3,4)
lyr <- substr(lyears[length(lyears)],3,4)

# Minimum requirements for the species to be included
minn <- 0

# Maximum population change for the species to be shown on the graphs
maxtr <- 100	# Set to 100-fold increase for the paper

# Get the list of BBS species
suffix <- ifelse(region=="uk", "ubbssws", "ubbs")
sprep <- read.csv(paste0("X:/CensusUnit/Trends",lyear+1,"/BBS/BBS_unsm_result_",fyr,lyr,"/index_",region,"_",suffix,".csv"),
	stringsAsFactors=F)
threshold <- ifelse(region=="uk", 40, 30)
sprep <- subset(sprep, (nsqus>=threshold & !species%in%c("BH","CM","LB","HG","GB","FF"))) #| species%in%c("N.","MN")))
splist <- subset(sprep, select=c(species, spname))
splist <- splist[order(splist$spname),]
splist$spname[splist$spname=="Great Spotted Woodpecker"] <- "G. Spotted Woodpecker"
splist$id <- 1:nrow(splist)
species <- splist$species

# Pull together all population indices
for (i in 1:length(species)){
	sp <- species[i]
	filename <- paste0("results/models/species/",sp,fyr,lyr,region,".csv")
	datsp <- read.csv(filename)
	colnames(datsp)[-1] <- paste0(sp,colnames(datsp)[-1])
	if(i==1) tab <- datsp else tab <- merge(tab, datsp, by="year")
}
write.csv(tab, paste0("results/models/allspecies",fyr,lyr,region,".csv"), row.names=F)

# Calculate trends between early years and late  years
classiccols <- grep("classicunsm",colnames(tab))
rvisitcols <- grep("rvisitunsm",colnames(tab))
gammdaycols <- grep("gammdayunsm",colnames(tab))
classic <- tab[,c(1,classiccols)]
rvisit <- tab[,c(1,rvisitcols)]
gammday <- tab[,c(1,gammdaycols)]
classicchg <- colMeans(classic[classic$year%in%lyears,])/colMeans(classic[classic$year%in%eyears,])
rvisitchg <- colMeans(rvisit[rvisit$year%in%lyears,])/colMeans(rvisit[rvisit$year%in%eyears,])
gammdaychg <- colMeans(gammday[gammday$year%in%lyears,])/colMeans(gammday[gammday$year%in%eyears,])
rvisitchg1 <- rvisitchg[-1]
gammdaychg1 <- gammdaychg[-1]

# Prepare table to export
BTOsplist <- read.csv("https://app.bto.org/general/taxonomy/species_list.csv",
			    head=F, stringsAsFactors = F)
BTOsplist <- subset(BTOsplist, V3!="", select=c("V1","V2","V3"))
colnames(BTOsplist) <- c("engname","scname","species")
BTOsplist$engname[BTOsplist$species=="PW"] <- "Pied wagtail"
BTOsplist$order <- 1:nrow(BTOsplist)
dfexport <- data.frame(substr(names(rvisitchg1),1,2), rvisitchg1,
			     gammdaychg1)
dfexport$log.rvisitchg1 <- log(dfexport$rvisitchg1)
dfexport$log.gammdaychg1 <- log(dfexport$gammdaychg1)
colnames(dfexport) <- c("species","Visit timing not accounted for",
				"Visit timing accounted for",
				"Visit timing not accounted for (log)",
				"Visit timing accounted for (log)")
dfexport1 <- merge(dfexport, BTOsplist, by="species", all.x=T)
intnames <- read.csv("data/IOC_Names_File_Plus-10.1.csv")
intnames <- subset(intnames, select=c(English.name, Scientific.Name))
dfexport2 <- merge(dfexport1, intnames, by.x="scname", by.y="Scientific.Name", all.x=T)

dfexport2$English.name[dfexport2$species=="LI"] <- "Common Linnet"
dfexport2$scname[dfexport2$species=="LI"] <- "Linaria cannabina"
dfexport2$English.name[dfexport2$species=="JD"] <- "Western Jackdaw"
dfexport2$scname[dfexport2$species=="JD"] <- "Coloeus monedula"
dfexport2$English.name[dfexport2$species=="GA"] <- "Gadwall"
dfexport2$scname[dfexport2$species=="GA"] <- "Mareca strepera"
dfexport2$English.name[dfexport2$species=="LR"] <- "Lesser Redpoll"
dfexport2$scname[dfexport2$species=="LR"] <- "Acanthis cabaret"
dfexport2$English.name[dfexport2$species=="SK"] <- "Eurasian Siskin"
dfexport2$scname[dfexport2$species=="SK"] <- "Spinus spinus"
dfexport2$English.name[dfexport2$species=="WT"] <- "Willow Tit"
dfexport2$scname[dfexport2$species=="WT"] <- "Poecile montanus"

dfexport2 <- dfexport2[order(dfexport2$scname),]

# Export table
write.csv(dfexport2,
	paste0("results/models/Trend comparison_",region,".csv"), row.names=F)

# Exclude species with too large increase from the graph
toolarge <- which(gammdaychg1>maxtr | rvisitchg1>maxtr)
rvisitchg2 <- rvisitchg1[-toolarge]
gammdaychg2 <- gammdaychg1[-toolarge]

# Plot gammdaychg2 against rvisitchg2 (figure 6 in the main text)
if (out=="screen") windows()
if (out=="png") png(paste0("results/graphs/",region,"/_trend comparison.png"),
				 width=800, height=800, pointsize=20)
if(out=="eps") {
	setEPS(width=12, height=12, pointsize=20)
	postscript(paste0("results/graphs/",region,"/_trend comparison.eps"))
} 
# Exclude species with too large increase from the graph
toolarge <- which(gammdaychg1>maxtr | rvisitchg1>maxtr)
rvisitchg2 <- rvisitchg1[-toolarge]
gammdaychg2 <- gammdaychg1[-toolarge]
log.rvisitchg2 <- log(rvisitchg2)
log.gammdaychg2 <- log(gammdaychg2)
xlim <- ylim <- range(log.rvisitchg2, log.gammdaychg2)
# Plot
plot(log.rvisitchg2, log.gammdaychg2, type="n", xlim=xlim, ylim=ylim,
     xlab="Trend estimate not accounting for visit timing",
     ylab="Trend estimate accounting for visit timing", axes=F)
axis(1, pos=ylim[1]-0.01)
axis(2, pos=xlim[1]-0.01)
segments(xlim[1]-0.01, ylim[1]-0.01, xlim[2]+0.01, ylim[2]-0.01, col="grey")
segments(xlim[1]-0.01, ylim[1]-0.01, xlim[2]+0.01, ylim[1]-0.01)
segments(xlim[1]-0.01, ylim[1]-0.01, xlim[1]-0.01, ylim[2]+0.01)
points(log.rvisitchg2, log.gammdaychg2, pch=16)
idsp <- c("RI","CG","MT","SH","KT")
log.rvisitchg2.id <- log.rvisitchg2[substr(names(rvisitchg2),1,2) %in% idsp]

if(out!="screen") dev.off()

########
# Assess potential changes in red-listing if models accounted for changes
# in visit timing. Bit added in response to referee.

tcomp <- subset(dfexport2, select=c("engname","scname",
						"Visit timing not accounted for",
						"Visit timing accounted for"))
colnames(tcomp) <- c("engname","scname","naf","af")

# Estimate red listing based on:
# 1) visit timing not accounted for
tcomp$naf25 <- (tcomp$naf)^(25/(mean(lyears)-mean(eyears)))
tcomp$naf.redl <- ifelse(tcomp$naf<=0.5, "red",
				 ifelse(tcomp$naf>0.5 & tcomp$naf <=0.75, "amber",
				 	 "green"))
tcomp$naf25.redl <- ifelse(tcomp$naf25<=0.5, "red",
				 ifelse(tcomp$naf25>0.5 & tcomp$naf25 <=0.75, "amber",
				 	 "green"))
# 2) visit timing accounted for
tcomp$af25 <- (tcomp$af)^(25/(mean(lyears)-mean(eyears)))
tcomp$af.redl <- ifelse(tcomp$af<=0.5, "red",
				 ifelse(tcomp$af>0.5 & tcomp$af <=0.75, "amber",
				 	 "green"))
tcomp$af25.redl <- ifelse(tcomp$af25<=0.5, "red",
				 ifelse(tcomp$af25>0.5 & tcomp$af25 <=0.75, "amber",
				 	 "green"))

# Real red list
bocc <- read.csv("data/bocc.csv")
names(bocc)[names(bocc)=="spname"] <- "engname"

# Merge our estimated red list to the real one
tcomp1 <- merge(tcomp, bocc, by="engname", all.x=T)

# Which species change classification if accounting for visit timing
change <- subset(tcomp1, naf.redl!=af.redl | naf25.redl!=af25.redl)

# Export table
write.csv(tcomp1, paste0("results/models/trend comparison with red listing_",
				 region,".csv"), row.names=F)


########
# Additional bit to add standard errors to the table with trend comparisons.
# This was added in response to a referee.

# Pull together all population indices
for (i in 1:length(species)){
	sp <- species[i]
	filename.boot <- paste0("results/models/species/",sp,fyr,lyr,region,"_boot.csv")
	datsp.boot <- read.csv(filename.boot)
	colnames(datsp.boot)[-c(1,ncol(datsp.boot))] <- paste0(sp,colnames(datsp.boot)[-c(1,ncol(datsp.boot))])
	if(i==1) tab.boot <- datsp.boot else tab.boot <- merge(tab.boot, datsp.boot, by=c("year","boot"))
}
tab.boot$boot <- as.numeric(as.character(tab.boot$boot))
write.csv(tab.boot, paste0("results/models/allspecies",fyr,lyr,region,"_boot.csv"), row.names=F)

# Read point trends
#tab <- read.csv(paste0("results/models/allspecies",fyr,lyr,region,".csv"))

# Calculate trend between early years and late  years
classiccols <- grep("classicunsm",colnames(tab))
rvisitcols <- grep("rvisitunsm",colnames(tab))
gammdaycols <- grep("gammdayunsm",colnames(tab))
classiccols.boot <- grep("classicunsm",colnames(tab.boot))
rvisitcols.boot <- grep("rvisitunsm",colnames(tab.boot))
gammdaycols.boot <- grep("gammdayunsm",colnames(tab.boot))
classic <- tab[,c(1,classiccols)]
rvisit <- tab[,c(1,rvisitcols)]
gammday <- tab[,c(1,gammdaycols)]
classic.boot <- tab.boot[,c(1,2,classiccols.boot)]
rvisit.boot <- tab.boot[,c(1,2,rvisitcols.boot)]
gammday.boot <- tab.boot[,c(1,2,gammdaycols.boot)]
classicchg <- colMeans(classic[classic$year%in%lyears,])/colMeans(classic[classic$year%in%eyears,])
rvisitchg <- colMeans(rvisit[rvisit$year%in%lyears,])/colMeans(rvisit[rvisit$year%in%eyears,])
gammdaychg <- colMeans(gammday[gammday$year%in%lyears,])/colMeans(gammday[gammday$year%in%eyears,])
classic.boot.lyears <- classic.boot[classic.boot$year%in%lyears,]
classic.boot.eyears <- classic.boot[classic.boot$year%in%eyears,]
rvisit.boot.lyears <- rvisit.boot[rvisit.boot$year%in%lyears,]
rvisit.boot.eyears <- rvisit.boot[rvisit.boot$year%in%eyears,]
gammday.boot.lyears <- gammday.boot[gammday.boot$year%in%lyears,]
gammday.boot.eyears <- gammday.boot[gammday.boot$year%in%eyears,]
classic.boot.lyears.av <- aggregate(.~boot, data=classic.boot.lyears, mean)
classic.boot.eyears.av <- aggregate(.~boot, data=classic.boot.eyears, mean)
rvisit.boot.lyears.av <- aggregate(.~boot, data=rvisit.boot.lyears, mean)
rvisit.boot.eyears.av <- aggregate(.~boot, data=rvisit.boot.eyears, mean)
gammday.boot.lyears.av <- aggregate(.~boot, data=gammday.boot.lyears, mean)
gammday.boot.eyears.av <- aggregate(.~boot, data=gammday.boot.eyears, mean)
classicchg.boot <- classic.boot.lyears.av/classic.boot.eyears.av
rvisitchg.boot <- rvisit.boot.lyears.av/rvisit.boot.eyears.av
gammdaychg.boot <- gammday.boot.lyears.av/gammday.boot.eyears.av
rvisitchg1 <- rvisitchg[-1]
gammdaychg1 <- gammdaychg[-1]

# Calculate standard errors and transpose
classicchg.boot.se <- apply(classicchg.boot, 2, sd)
rvisitchg.boot.se <- apply(rvisitchg.boot, 2, sd)
gammdaychg.boot.se <- apply(gammdaychg.boot, 2, sd)
rvisitchg.boot.seT <- data.frame(substr(names(rvisitchg.boot.se),1,2), rvisitchg.boot.se)
gammdaychg.boot.seT <- data.frame(substr(names(gammdaychg.boot.se),1,2), gammdaychg.boot.se)
colnames(rvisitchg.boot.seT)[1] <- colnames(gammdaychg.boot.seT)[1] <- "species"
boot.seT <- merge(rvisitchg.boot.seT, gammdaychg.boot.seT, by="species")
boot.seT <- subset(boot.seT, !species%in%c("bo","ye"))

# Prepare table to export
BTOsplist <- read.csv("https://app.bto.org/general/taxonomy/species_list.csv",
			    head=F, stringsAsFactors = F)
BTOsplist <- subset(BTOsplist, V3!="", select=c("V1","V2","V3"))
colnames(BTOsplist) <- c("engname","scname","species")
BTOsplist$engname[BTOsplist$species=="PW"] <- "Pied wagtail"
BTOsplist$order <- 1:nrow(BTOsplist)
dfexport <- data.frame(species=substr(names(rvisitchg1),1,2), rvisitchg1,
			     gammdaychg1)
dfexport.se <- merge(dfexport, boot.seT, by="species")
colnames(dfexport.se) <- c("species","Visit timing not accounted for",
				"Visit timing accounted for",
				"SE of visit timing not accounted for",
				"SE of visit timing accounted for")
dfexport1 <- merge(dfexport.se, BTOsplist, by="species", all.x=T)
intnames <- read.csv("data/IOC_Names_File_Plus-10.1.csv")
intnames <- subset(intnames, select=c(English.name, Scientific.Name))
dfexport2 <- merge(dfexport1, intnames, by.x="scname", by.y="Scientific.Name", all.x=T)

# Correct scientific names
dfexport2$English.name[dfexport2$species=="LI"] <- "Common Linnet"
dfexport2$scname[dfexport2$species=="LI"] <- "Linaria cannabina"
dfexport2$English.name[dfexport2$species=="JD"] <- "Western Jackdaw"
dfexport2$scname[dfexport2$species=="JD"] <- "Coloeus monedula"
dfexport2$English.name[dfexport2$species=="GA"] <- "Gadwall"
dfexport2$scname[dfexport2$species=="GA"] <- "Mareca strepera"
dfexport2$English.name[dfexport2$species=="LR"] <- "Lesser Redpoll"
dfexport2$scname[dfexport2$species=="LR"] <- "Acanthis cabaret"
dfexport2$English.name[dfexport2$species=="SK"] <- "Eurasian Siskin"
dfexport2$scname[dfexport2$species=="SK"] <- "Spinus spinus"
dfexport2$English.name[dfexport2$species=="WT"] <- "Willow Tit"
dfexport2$scname[dfexport2$species=="WT"] <- "Poecile montanus"

# Order by scientific name
dfexport2 <- dfexport2[order(dfexport2$scname),]

write.csv(dfexport2,
	paste0("results/models/Trend comparison_",region,"_with_SEs.csv"),
	row.names=F)

