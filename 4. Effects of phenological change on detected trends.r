############################################################################
# Estimate the biases in the estimated trends caused by changes in         #
# species' phenology and visit timing for the chapter:                     #
# Results/Effects of phenological change on detected trends                          #
############################################################################

# Define the study periods (eyears=period 1, lyears=period 2)
eyears <- 1994:1998
lyears <- 2013:2017

# Specify format of figures
out <- "eps"

# Specify study region
region <- "rda7"

# Specify days since 1th April to be used in the models. Visits outside the
# permitted season (April-June) are excluded
intmodel <- 1:91

k <- 10		# Dimension of the basis for the smooth term

# Graphical parameters:
cex <- 0.5
cex.main <- 0.8
cex.axis <- 0.8
cex.lab <- 0.8
lwd=1.5

setwd("N:/Phenology and trends")

library(mgcv)
library(ggplot2)
library(ggrepel)

options(stringsAsFactors=F)

fyr <- substr(eyears[1],3,4)
lyr <- substr(lyears[length(lyears)],3,4)

# Define name of the BBS weighting file in Rdata format
filename <- paste0("data/weights ",fyr,lyr,"/bbswgt-",fyr,lyr,"-",region,
			 "_surveyed_only.rdata")

# Load weighting file in Rdata format if it already exists, otherwise read it
# from the BBS archive and save it into Rdata format
if(file.exists(filename)) load(filename) else {
	datfilename <- paste0("N:/Visit timing/data/weights ",fyr,lyr,"/bbswgt-",
				    fyr,lyr,"-",region,".dat")
	dat <- read.table(datfilename, col.names=c("site","sp","count","year",
								 "SASdate","eorl","weight",
								 "region","special", "class"),
				fill=T)
	# For this analysis I am not interested in non-surveyed squares
	dat <- subset(dat, SASdate!=".")
	save(dat, file=filename)
}

# Create date variable day1A expressed as number of days since 1st April
# (1st April = 1)
dat$SASdate <- as.numeric(as.character(dat$SASdate))
dat$count <- as.numeric(as.character(dat$count))
dat$date <- as.Date(dat$SASdate, origin="1960-01-01")
dat$dayy <- as.numeric(strftime(dat$date, format="%j"))
dat$day1A <- dat$dayy-90

# Get rid of surveys outside the BBS season
dat1 <- subset(dat, day1A %in% intmodel)

# Get the list of BBS species
suffix <- ifelse(region=="uk", "ubbssws", "ubbs")
sprep <- read.csv(paste0("X:/CensusUnit/Trends",
				 lyears[length(lyears)]+1,"/BBS/BBS_unsm_result_",fyr,
				 lyr,"/index_",region,"_",suffix,".csv"),
	stringsAsFactors=F)
threshold <- ifelse(region=="uk", 40, 30)
sprep <- subset(sprep, (nsqus>=threshold &
					!species%in%c("BH","CM","LB","HG","GB","FF")))
splist <- subset(sprep, select=c(species, spname))
splist <- splist[order(splist$spname),]
splist$spname[splist$spname=="Great Spotted Woodpecker"] <- "G. Spotted Woodpecker"
splist$id <- 1:nrow(splist)

# Loop over all species
for (i in 1:nrow(splist)) {
	
	# Get BTO 2-letter code of the species being processed
	sp <- splist$species[i]
	print(paste("Processing",sp))
	
	# Subset dataset by the species being processed
	datsp <- dat1[dat1$sp==sp,]

	# Prepare the plotting window or file
	if(out=="screen") print("Can't show all single species graphs")
	if(out=="png") png(paste0("results/graphs/",region,"/",sp,".png"),
				 width=960, height=960, pointsize=20)
	if(out!="screen") par(mfrow=c(2,2))
	
	# Get distribution of visits per day. Only use squares where the
	# species has been recorded at least once.
	allvisitss <- unique(subset(datsp, select=-c(sp, count)))
	allvisitss.t0 <- subset(allvisitss, year%in% eyears)
	allvisitss.t1 <- subset(allvisitss, year%in% lyears)
	vispdays.t0 <- aggregate(site~day1A, data=allvisitss.t0, length)
	vispdays.t1 <- aggregate(site~day1A, data=allvisitss.t1, length)
	ymax <- max(vispdays.t0$site, vispdays.t1$site)
	if(out!="screen"){
		plot(vispdays.t0$day1A, vispdays.t0$site, ylim=c(0,ymax),
		     main=c(region, sp), xlab="Days since 1st April",
		     ylab="Number of visits", pch=16, col="lightblue", cex=cex,
		     cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab)
		points(vispdays.t1$day1A, vispdays.t1$site, pch=16, col="blue",
			 cex=cex)
	}
	
	# Smooth the distribution of visits per day
	gams.t0 <- gam(site~s(day1A), data=vispdays.t0, family="poisson")
	gams.t1 <- gam(site~s(day1A), data=vispdays.t1, family="poisson")
	newdats <- data.frame(day1A=intmodel)
	survs <- data.frame(newdats, surv0=predict(gams.t0, newdats, "response"),
				  surv1=predict(gams.t1, newdats, "response"))
	if(out!="screen") {
		lines(survs$day1A, survs$surv0, col="lightblue", lwd=lwd)
		lines(survs$day1A, survs$surv1, col="blue", lwd=lwd)
	}

	# Divide smoothed visits by their sum so that sum=1
	survs$surv0n <- survs$surv0/sum(survs$surv0)
	survs$surv1n <- survs$surv1/sum(survs$surv1)

	# Keep only years of interest
	datsp.t0pre <- subset(datsp, year %in% eyears)
	datsp.t1pre <- subset(datsp, year %in% lyears)
	datsp.t0 <- datsp.t0pre
	datsp.t1 <- datsp.t1pre

	# Get distribution of sightings per day
	sppday.t0 <- aggregate(count~day1A, data=datsp.t0, mean)
	temp <- aggregate(count~day1A, data=datsp.t0, length)
	sppday.t1 <- aggregate(count~day1A, data=datsp.t1, mean)
	
	# Smooth the distribution of sightings per day
	gamsp.t0 <- gam(count~s(day1A, k=k), data=datsp.t0, family="quasipoisson")
	gamsp.t1 <- gam(count~s(day1A, k=k), data=datsp.t1, family="quasipoisson")
	newdatsp <- data.frame(day1A=intmodel)
	spec <- data.frame(newdatsp, spec0=predict(gamsp.t0, newdatsp, "response"),
	 			  spec1=predict(gamsp.t1, newdatsp, "response"))

	# Plot distribution of sightings
	ymax <- max(sppday.t0$count, sppday.t1$count)
	if(out!="screen"){
		plot(sppday.t0$day1A, sppday.t0$count, ylim=c(0, ymax),
		     main=c(region, sp), xlab="Days since 1st April",
		     ylab="Mean count", pch=16, col="chartreuse3", cex=cex,
		     cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab)
		points(sppday.t1$day1A, sppday.t1$count, pch=16, cex=cex,
			 col="green4")
		lines(spec$day1A, spec$spec0, col="chartreuse3", lwd=lwd)
		lines(spec$day1A, spec$spec1, col="green4", lwd=lwd)
	}

	# Divide smoothed sightings by their sum so that sum=1
	spec$spec0n <- spec$spec0/sum(spec$spec0)
	spec$spec1n <- spec$spec1/sum(spec$spec1)

	# Plot distribution of visits and of sightings 
	ymax <- max(survs$surv0n, survs$surv1n, spec$spec0n, spec$spec1n)
	if(out!="screen"){
		plot(survs$day1A, survs$surv0n, type="l", xlim=range(intmodel),
		     ylim=c(0,ymax), main=c(region, sp),
		     xlab="Days since 1st April",
		     ylab="Density of visits and sightings", cex=cex,
		     cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab,
		     col="lightblue", lwd=lwd)
		lines(survs$day1A, survs$surv1n, col="blue", lwd=lwd)
		lines(spec$day1A, spec$spec0n, col="chartreuse3", lwd=lwd)
		lines(spec$day1A, spec$spec1n, col="green4", lwd=lwd)
		#abline(v=c(0,91))
		dev.off()
	}
	
	# If surveyors had not changed their own phenology, What would be the
	# bias in trends caused by changes in the species' phenology?
	survx <- subset(survs, day1A %in% intmodel)
	specx <- subset(spec, day1A %in% intmodel)
	survNspecY <- sum(survx$surv0n*specx$spec1n)/
		sum(survx$surv0n*specx$spec0n)
	
	# What is the bias in trends caused by combined changes in surveyors' AND
	# species' phenology
	survYspecY <- sum(survx$surv1n*specx$spec1n)/
		sum(survx$surv0n*specx$spec0n)

	# If the species had not changed its phenology, What would be the
	# bias in trends caused by changes in the surveyors' phenology?
	survYspecN <- sum(survx$surv1n*specx$spec0n)/
		sum(survx$surv0n*specx$spec0n)
	
	temp <- data.frame(sp=sp, survNspecY=survNspecY, survYspecY=survYspecY,
		survYspecN=survYspecN,
		n0=nrow(datsp.t0), n1=nrow(datsp.t1),
		nhil0=nrow(datsp.t0pre)-nrow(datsp.t0),
		nhil1=nrow(datsp.t1pre) - nrow(datsp.t1),
		meancount0=mean(datsp.t0$count), meancount1=mean(datsp.t1$count),
		medianE0=median(allvisitss.t0$day1A[allvisitss.t0$eorl=="E"]),
		medianE1=median(allvisitss.t1$day1A[allvisitss.t1$eorl=="E"]),
		medianL0=median(allvisitss.t0$day1A[allvisitss.t0$eorl=="L"]),
		medianL1=median(allvisitss.t1$day1A[allvisitss.t1$eorl=="L"]),
		mediansp0=median(rep(datsp.t0$day1A, time=datsp.t0$count)),
		mediansp1=median(rep(datsp.t1$day1A, time=datsp.t1$count)))
	
	if(i==1) tab <- temp else tab <- rbind(tab, temp)

}

# Read species' names
BTOsplist <- read.csv("https://app.bto.org/general/taxonomy/species_list.csv",
			    head=F, stringsAsFactors = F)
BTOsplist <- subset(BTOsplist, V3!="", select=c("V1","V2","V3"))
colnames(BTOsplist) <- c("engname","scname","species")
BTOsplist$engname[BTOsplist$species=="PW"] <- "Pied wagtail"
BTOsplist$order <- 1:nrow(BTOsplist)

# Merge big table to species' names
tab.exp <- merge(tab, BTOsplist, by.x="sp", by.y="species", all.x=T)
tab.exp <- tab.exp[order(tab.exp$scname),]

# Read international species' names
intnames <- read.csv("data/IOC_Names_File_Plus-10.1.csv")
intnames <- subset(intnames, select=c(English.name, Scientific.Name))

# Merge internatioinal species' names
tab.exp1 <- merge(tab.exp, intnames, by.x="scname", by.y="Scientific.Name",
		     all.x=T)
tab.exp1 <- subset(tab.exp1, select=-engname)
tab.exp1$English.name[tab.exp1$sp=="LI"] <- "Common Linnet"
tab.exp1$scname[tab.exp1$sp=="LI"] <- "Linaria cannabina"
tab.exp1$English.name[tab.exp1$sp=="JD"] <- "Western Jackdaw"
tab.exp1$scname[tab.exp1$sp=="JD"] <- "Coloeus monedula"
tab.exp1$English.name[tab.exp1$sp=="GA"] <- "Gadwall"
tab.exp1$scname[tab.exp1$sp=="GA"] <- "Mareca strepera"
tab.exp1$English.name[tab.exp1$sp=="LR"] <- "Lesser Redpoll"
tab.exp1$scname[tab.exp1$sp=="LR"] <- "Acanthis cabaret"
tab.exp1$English.name[tab.exp1$sp=="SK"] <- "Eurasian Siskin"
tab.exp1$scname[tab.exp1$sp=="SK"] <- "Spinus spinus"
tab.exp1$English.name[tab.exp1$sp=="WT"] <- "Willow Tit"
tab.exp1$scname[tab.exp1$sp=="WT"] <- "Poecile montanus"

# Order by scientific name
tab.exp1 <- tab.exp1[order(tab.exp1$scname),]

# Export table for the Supplementary Material
write.csv(tab.exp1, paste0("results/",region,"_apparent trends.csv"),
	    row.names=F)
	
hist(tab$survNspecY)
hist(tab$survYspecY)

# Quantiles to report in the paper
quantile(tab$survNspecY)
quantile(tab$survYspecY)

# Prepare plot for figure 5
if (out=="screen") windows()
tab$log.survNspecY <- log(tab$survNspecY)
tab$log.survYspecY <- log(tab$survYspecY)
xlim <- ylim <- range(tab$log.survNspecY, tab$log.survYspecY)
# Identify species with large bias so that they can be labelled on the graph
large.bias.survNspecY <- subset(tab, abs(log.survNspecY) > 0.045)
large.bias.survYspecY <- subset(tab, abs(log.survYspecY) > 0.045)
idsp <- unique(c(large.bias.survNspecY$sp, large.bias.survYspecY$sp))
ggplot(tab, aes(log.survNspecY, log.survYspecY, label=sp)) +
	theme(panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		axis.line = element_line(colour = "black")) +
	labs(x=expression("log(C"[s]*") - Bias produced by changes in species phenology only"),
	     y=expression("log(C"[sv]*") - Bias produced by changes in both species and surveyor phenology")) +
	xlim(xlim) + ylim(ylim) +
	geom_abline(slope=1, intercept=0, col="grey") +
	geom_vline(xintercept=0, col="grey") +
	geom_hline(yintercept=0, col="grey") + geom_point() +
	geom_text_repel(data=subset(tab, sp %in% idsp), min.segment.length=0.1)
if (out=="png") ggsave(paste0("results/graphs/",region,"/_scatterplot.png"),
			     width=28.2, height=28.2, units="cm", scale=0.6)
if (out=="eps") ggsave(paste0("results/graphs/",region,"/_scatterplot.eps"),
			     width=28.2, height=28.2, units="cm", scale=0.6)

