############################################################################
# Produce the results for the chapter:                                     #
# Results/Estimating population trends by accounting for changes in visit  #
# timing                                                                   #
############################################################################

# Define the study periods (eyears=period 1, lyears=period 2)
eyears <- 1994:1998
lyears <- 2013:2017

# Specify format of figures
out <- "eps"

# Specify study region
region <- "rda7"

# Specify 4 species for the example figure in the main text (all species will
# be shown in the supplementary material)
whichspecies <- c("WW","NH","CK","SL")

# Specify days since 1th April to be used in the models. Visits outside the
# permitted season (April-June) are excluded
intmodel <- 1:91

setwd("N:/Phenology and trends")

library(mgcv)

fyr <- substr(eyears[1],3,4)
lyr <- substr(lyears[length(lyears)],3,4)

# Define name of the BBS weighting file in Rdata format
filename <- paste0("data/weights ",fyr,lyr,"/bbswgt-",fyr,lyr,"-",region,
			 "_surveyed_only.rdata")

# Load weighting file in Rdata format if it already exists, otherwise read it
# from the BBS archive and save it into Rdata format
if(file.exists(filename)) load(filename) else {
	datfilename <- paste0("N:/Visit timing/data/weights ",fyr,lyr,
				    "/bbswgt-",fyr,lyr,"-",region,".dat")
	dat <- read.table(datfilename,
		col.names=c("site","sp","count","year","SASdate","eorl","weight",
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
sprep <- read.csv(paste0("X:/CensusUnit/Trends",lyears[length(lyears)]+1,
				 "/BBS/BBS_unsm_result_",fyr,lyr,"/index_",region,
				 "_",suffix,".csv"), stringsAsFactors=F)
threshold <- ifelse(region=="uk", 40, 30)
sprep <- subset(sprep, (nsqus>=threshold &
					!species%in%c("BH","CM","LB","HG","GB","FF")))
splistpre <- subset(sprep, select=c(species, spname))
splistpre$spname[splistpre$spname=="Great Spotted Woodpecker"] <-
	"G. Spotted Woodpecker"

# Read species' names
BTOsplist <- read.csv("https://app.bto.org/general/taxonomy/species_list.csv",
			    head=F, stringsAsFactors = F)
BTOsplist <- subset(BTOsplist, V3!="", select=c("V1","V2","V3"))
colnames(BTOsplist) <- c("engname","scname","species")
BTOsplist$engname[BTOsplist$species=="PW"] <- "Pied wagtail"
BTOsplist$order <- 1:nrow(BTOsplist)

# Read international species' names
intnames <- read.csv("data/IOC_Names_File_Plus-10.1.csv")
intnames <- subset(intnames, select=c(English.name, Scientific.Name))

# Merge international species' names
splist0 <- merge(splistpre, BTOsplist, by="species", all.x=T)
splist <- merge(splist0, intnames, by.x="scname", by.y="Scientific.Name",
		     all.x=T)
splist$English.name[splist$species=="LI"] <- "Common Linnet"
splist$scname[splist$species=="LI"] <- "Linaria cannabina"
splist$English.name[splist$species=="JD"] <- "Western Jackdaw"
splist$scname[splist$species=="JD"] <- "Coloeus monedula"
splist$English.name[splist$species=="GA"] <- "Gadwall"
splist$scname[splist$species=="GA"] <- "Mareca strepera"
splist$English.name[splist$species=="LR"] <- "Lesser Redpoll"
splist$scname[splist$species=="LR"] <- "Acanthis cabaret"
splist$English.name[splist$species=="SK"] <- "Eurasian Siskin"
splist$scname[splist$species=="SK"] <- "Spinus spinus"
splist$English.name[splist$species=="WT"] <- "Willow Tit"
splist$scname[splist$species=="WT"] <- "Poecile montanus"

# Export the splist list
write.csv(splist, paste0("results/species list ",region,".csv"), row.names=F)

########
# 1. Plot the 4 species for figure 3 in the main text 
#

splist.paper <- subset(splist, species %in% whichspecies)
splist.paper <- splist.paper[order(splist.paper$order),]

# Prepare the plotting window or file
if(out=="screen") windows()
filename <- paste0("results/graphs/", region, "/species only/",
				  paste0(whichspecies, collapse=""))
if(out=="png") png(paste0(filename, ".png"), width=1600, height=900,
			 pointsize=30)
if(out=="wmf") win.metafile(paste0(filename, ".wmf"), width=6, height=6)

if(out=="eps") {
	setEPS(width=12, height=8, pointsize=20)
	postscript(paste0(filename, ".eps"))
} 

par(mfrow=c(2,2), mar=c(3,4,3,1))

for (i in 1:nrow(splist.paper)) {
	
	sp <- splist.paper$species[i]
	engname <- splist.paper$spname[i]
	scname <- splist.paper$scname[i]
	print(paste(sp,engname,scname))
	
	datsp <- dat1[dat1$sp==sp,]

	# Keep only years of interest
	datsp.t0pre <- subset(datsp, year %in% eyears)
	datsp.t1pre <- subset(datsp, year %in% lyears)
	datsp.t0 <- datsp.t0pre
	datsp.t1 <- datsp.t1pre

	sppday.t0 <- aggregate(count~day1A, data=datsp.t0, mean)
	sppday.t1 <- aggregate(count~day1A, data=datsp.t1, mean)

	# Calculate median
	median.t0 <- median(rep(sppday.t0$day1A, sppday.t0$count*1000))
	median.t1 <- median(rep(sppday.t1$day1A, sppday.t1$count*1000))

	# Plot distribution of sightings
	ymax <- max(sppday.t0$count, sppday.t1$count)
	ylab <- if (i %in% c(1,3)) "Mean count" else ""
	plot(sppday.t1$day1A, sppday.t1$count, xlim=range(intmodel),
	     ylim=c(0, ymax), pch=16, axes=F, xlab="", ylab=ylab)
	mtext(paste0(letters[i],")"), side=3, adj=-0.1, line=1.5, font=2)
	mtext(scname, side=3, line=1.5, font=3)
	points(sppday.t0$day1A, sppday.t0$count, pch=1)
	axis(1, at=c(0,30,61), labels=c("1/4", "1/5", "1/6"), pos=0)
	segments(1,0,91,0)
	segments(0,0,0,ymax)
	axis(2, pos=0)
	abline(v=median.t0, lty=2)
	abline(v=median.t1, lty=1)
}

if(out!="screen") dev.off()

########
# 2. Plot all species for the supplementary material
# Also calculate the median shift by species

mfrow <- c(5,2)

gpp <- mfrow[1]*mfrow[2]

splist <- splist[order(splist$scname),]

# Prepare the plotting window or file
nwin <- ceiling(nrow(splist)/gpp)

for (j in 1:nwin) {
			
	if(out=="screen") windows(2.5*mfrow[2], 2.5*mfrow[1])
	if(out=="png") png(paste0("results/graphs/", region,
	 				  "/species only/all species - batch ",j,".png"),
	 			 width=640*mfrow[2], height=360*mfrow[1],
				 pointsize=30)
	par(mfrow=mfrow, mar=c(3,4,3,1))

	for (i in (j*gpp-gpp+1):(j*gpp)) {
	
		sp <- splist$species[i]
		engname <- splist$English.name[i]
		if(is.na(sp)) break()
		scname <- splist$scname[i]
		print(paste(sp,engname,scname))
	
		datsp <- dat1[dat1$sp==sp,]
	
		# Keep only years of interest
		datsp.t0 <- subset(datsp, year %in% eyears)
		datsp.t1 <- subset(datsp, year %in% lyears)

		sppday.t0 <- aggregate(count~day1A, data=datsp.t0, mean)
		sppday.t1 <- aggregate(count~day1A, data=datsp.t1, mean)

		# Calculate and save median
		median.t0 <- median(rep(sppday.t0$day1A, sppday.t0$count*1000))
		median.t1 <- median(rep(sppday.t1$day1A, sppday.t1$count*1000))
		temp <- data.frame(English_name=engname, Scientific_name=scname,
					 median_P1=median.t0, median_P2=median.t1)
		medians <- if(i==1 & j==1) temp else rbind(medians, temp)

		print(paste("Medians:",median.t0,median.t1))
	
		# Plot distribution of sightings
		ymax <- max(sppday.t0$count, sppday.t1$count)
		ylab <- if (i%%mfrow[2]==1) "Mean count" else ""
		plot(sppday.t1$day1A, sppday.t1$count, xlim=range(intmodel),
		     ylim=c(0, ymax), pch=16, axes=F, xlab="", ylab=ylab,
		     main=scname, font.main=3)
		points(sppday.t0$day1A, sppday.t0$count, pch=1)
		axis(1, at=c(0,30,61), labels=c("1/4", "1/5", "1/6"), pos=0)
		segments(1,0,91,0)
		segments(0,0,0,ymax)
		axis(2, pos=0)
		abline(v=median.t0, lty=2)
		abline(v=median.t1, lty=1)
		
		# Bootstrap over sites
		nboot <- 2000
		sites.t0 <- unique(datsp.t0$site)
		sites.t1 <- unique(datsp.t1$site)
		for (b in 1:nboot) {
			bt0 <- sort(sample(sites.t0, length(sites.t0), replace=T))
			bt1 <- sort(sample(sites.t1, length(sites.t1), replace=T))
			bootsites.t0 <- data.frame(site=bt0)
			bootsites.t1 <- data.frame(site=bt1)
			bootsp.t0 <- merge(bootsites.t0, datsp.t0, by="site",
						 all.x=T)
			bootsp.t1 <- merge(bootsites.t1, datsp.t1, by="site",
						 all.x=T)
			bootday.t0 <- aggregate(count~day1A, data=bootsp.t0, mean)
			bootday.t1 <- aggregate(count~day1A, data=bootsp.t1, mean)
			# Remove outliers (Tukey's fences)
			q0 <- quantile(bootday.t0$count, probs=c(0.25, 0.75))
			tukey.t0 <- c(q0[1]-3*(q0[2]-q0[1]), q0[1]+3*(q0[2]-q0[1]))
			bootday.t0 <- subset(bootday.t0, count>=tukey.t0[1] &
					  	count<=tukey.t0[2])
			q1 <- quantile(bootday.t1$count, probs=c(0.25, 0.75))
			tukey.t1 <- c(q1[1]-3*(q1[2]-q1[1]), q1[1]+3*(q1[2]-q1[1]))
			bootday.t1 <- subset(bootday.t1, count>=tukey.t1[1] &
					  	count<=tukey.t1[2])
			# Calculate and save median
			bootmedian.t0 <- median(rep(bootday.t0$day1A,
							    bootday.t0$count*1000))
			bootmedian.t1 <- median(rep(bootday.t1$day1A,
							    bootday.t1$count*1000))
			boottemp <- data.frame(English_name=engname,
						     Scientific_name=scname,
						     boot <- b, median_P1=bootmedian.t0,
						     median_P2=bootmedian.t1)
			bootmedians <- if(i==1 & j==1 & b==1) boottemp else
				rbind(bootmedians, boottemp)
		}
	}
	
	if(out!="screen") dev.off()
}

# Calculate difference in medians
medians$difference <- medians$median_P2 - medians$median_P1
bootmedians$difference <- bootmedians$median_P2 - bootmedians$median_P1

# Average across species
mean(medians$difference, na.rm=T)
bootmeans <- aggregate(difference~boot....b, data=bootmedians, mean)
quantile(bootmeans$difference, probs=c(0.025, 0.975))

# Calculate lower and upper confidence limits for each species
medians.lcl <- aggregate(cbind(median_P1,median_P2,difference)~English_name+
				 	Scientific_name,
				 data=bootmedians, quantile, probs=0.025)
medians.ucl <- aggregate(cbind(median_P1,median_P2,difference)~English_name+
				 	Scientific_name,
				 data=bootmedians, quantile, probs=0.975)
colnames(medians.lcl)[3:5] <- paste0(colnames(medians.lcl)[3:5],"_lcl")
colnames(medians.ucl)[3:5] <- paste0(colnames(medians.ucl)[3:5],"_ucl")

# Merge confidence limits to central estimates
medians.cis <- merge(medians.lcl, medians.ucl,
			   by=c("English_name","Scientific_name"))
medians.wcis <- merge(medians, medians.cis,
			    by=c("English_name","Scientific_name"), all.x=T)

# Export for table S2 in Supplementary Material
write.csv(medians, paste0("results/", region, "_median detection day.csv"),
	    row.names = F)
write.csv(medians.wcis, paste0("results/", region,
				  "_median detection day with Cis.csv"), row.names = F)

# order by difference to get in groups
medians1 <- medians.wcis[order(medians.wcis$difference),]
# Get a counter value for all species with same difference
medians1$order <- ave(medians1$difference, medians1$difference,
			    FUN = seq_along)
# Subtract 0.5 so the value can be used for plotting in the right place
medians1$order <- medians1$order - 0.5

# Add BTO 2-letter codes
medians2 <- merge(medians1, subset(splist, select=-order),
			by.x="Scientific_name", by.y="scname", all.x=T)

# Calculate which species had significant changes
medians2$sig <- sign(medians2$difference_lcl*medians2$difference_ucl)
medians.sig <- subset(medians2, sig>0)

# Plot histogram with 2-letter codes (figure 4 in the paper)
if(out=="screen") windows(9,5)
if(out=="png") png(paste0("results/graphs/", region,
			"/_histogram of species shifts with 2-letter codes.png"),
			width=1350, height=750, pointsize=24)
if(out=="eps") {
	setEPS(width=14.4, height=8, pointsize=20)
	postscript(paste0("results/graphs/", region,
			"/_histogram of species shifts with 2-letter codes.eps"))
} 
par(mar=c(4.5,4,0,0))
rangem <- range(medians$difference, na.rm=T)
hist(medians$difference, breaks=seq(rangem[1]-0.5, rangem[2]+0.5, by=1),
     main="", xlab="Phenological shift", ylab="Number of species",
     col="grey92", border="grey50")
rect(medians.sig$difference-0.5, medians.sig$order-0.5,
     medians.sig$difference+0.5, medians.sig$order+0.5, col="grey76",
     border=NA)
hist(medians$difference, breaks=seq(rangem[1]-0.5, rangem[2]+0.5, by=1),
     main="", xlab="Phenological shift", ylab="Number of species",
     col=NA, border="grey50", add=T)
text(medians2$difference, medians2$order, medians2$species, cex=0.6,
     col="grey20")
text(medians.sig$difference, medians.sig$order, medians.sig$species, cex=0.6,
     font=2)

if(out!="screen") dev.off()

########
# Additional bit to calculate changes sepatately for long-distance migrants
# and residents or short-distance migrants.
# This was added in response to a referee.

# Merge in information on species phenology
ph <- read.csv("results/species list with phenology.csv")
ph <- subset(ph, select=c(English.name, Long.distance.migrant))
medians.p <- merge(medians, ph, by.x="English_name", by.y="English.name")
bootmedians.p <- merge(bootmedians, ph, by.x="English_name", by.y="English.name")

medians.p$difference <- medians.p$median_P2 - medians.p$median_P1
bootmedians.p$difference <- bootmedians.p$median_P2 - bootmedians.p$median_P1

# Average across species
mean(medians.p$difference, na.rm=T)
bootmeans.p <- aggregate(difference~boot....b+Long.distance.migrant,
				 data=bootmedians.p, mean)
quantile(bootmeans.p$difference, probs=c(0.025, 0.975))

# Average across residents/short distance migrants only
mean(medians.p$difference[medians.p$Long.distance.migrant=="n"], na.rm=T)
mean(medians.p$difference[medians.p$Long.distance.migrant=="y"], na.rm=T)
quantile(bootmeans.p$difference[bootmeans.p$Long.distance.migrant=="n"], probs=c(0.025, 0.975))
quantile(bootmeans.p$difference[bootmeans.p$Long.distance.migrant=="y"], probs=c(0.025, 0.975))



