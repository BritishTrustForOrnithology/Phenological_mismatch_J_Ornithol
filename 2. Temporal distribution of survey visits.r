############################################################################
# Produce the results for the chapter:                                     #
# Results/Temporal distribution of survey visits'                          #
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
	datfilename <- paste0("N:/Visit timing/data/weights ",fyr,lyr,"/bbswgt-",
				    fyr,lyr,"-",region,".dat")
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
dat$date <- as.Date(dat$SASdate, origin="1960-01-01") # correct, checked
dat$dayy <- as.numeric(strftime(dat$date, format="%j"))
dat$day1A <- dat$dayy-90

# Get rid of surveys outside the BBS season
dat1 <- subset(dat, day1A %in% intmodel)

########
# Analyse the distribution of visits per day

# Get distribution of visits per day
dat2 <- subset (dat1, select=-c(sp, count, weight, region, special, class))
allvisits <- unique(dat2)
allvisits.t0 <- subset(allvisits, year%in% eyears)
allvisits.t1 <- subset(allvisits, year%in% lyears)

# Calculate medians
t0.E <- allvisits.t0$day1A[allvisits.t0$eorl=="E"]
t1.E <- allvisits.t1$day1A[allvisits.t1$eorl=="E"]
t0.L <- allvisits.t0$day1A[allvisits.t0$eorl=="L"]
t1.L <- allvisits.t1$day1A[allvisits.t1$eorl=="L"]
medianE.t0 <- median(t0.E) + 90
medianE.t1 <- median(t1.E) + 90
medianL.t0 <- median(t0.L) + 90
medianL.t1 <- median(t1.L) + 90
medians <- data.frame(median.t0=c(medianE.t0, medianL.t0),
	     median.t1=c(medianE.t1, medianL.t1))
medians$median.t0 <- as.Date(medians$median.t0, origin="1994-01-01")
medians$median.t1 <- as.Date(medians$median.t1, origin="1994-01-01")

# Print medians, to be reported in the paper
print(medians)

# Number of squares used, to be reported in the paper
length(unique(allvisits.t0$site))
length(unique(allvisits.t1$site))

# Bootstrap over sites
nboot <- 2000
sites.t0 <- unique(allvisits.t0$site)
sites.t1 <- unique(allvisits.t1$site)
for (b in 1:nboot) {
	if(b%%500==0) print(b)
	bt0 <- sort(sample(sites.t0, length(sites.t0), replace=T))
	bt1 <- sort(sample(sites.t1, length(sites.t1), replace=T))
	bootsites.t0 <- data.frame(site=bt0)
	bootsites.t1 <- data.frame(site=bt1)
	boot.t0 <- merge(bootsites.t0, allvisits.t0, by="site", all.x=T)
	boot.t1 <- merge(bootsites.t1, allvisits.t1, by="site", all.x=T)
	boot.t0.E <- boot.t0$day1A[boot.t0$eorl=="E"]
	boot.t1.E <- boot.t1$day1A[boot.t1$eorl=="E"]
	boot.t0.L <- boot.t0$day1A[boot.t0$eorl=="L"]
	boot.t1.L <- boot.t1$day1A[boot.t1$eorl=="L"]
	medianE.t0.boot <- median(boot.t0.E) + 90
	medianE.t1.boot <- median(boot.t1.E) + 90
	medianL.t0.boot <- median(boot.t0.L) + 90
	medianL.t1.boot <- median(boot.t1.L) + 90
	temp.boot <- data.frame(medianE.t0.boot, medianL.t0.boot,
					medianE.t1.boot, medianL.t1.boot)
	medians.boot <- if (b==1) temp.boot else rbind(medians.boot, temp.boot)
}
# Calculate difference of median dates for early and late visits
medians.boot$differenceE <- medians.boot$medianE.t1.boot-
	medians.boot$medianE.t0.boot
medians.boot$differenceL <- medians.boot$medianL.t1.boot-
	medians.boot$medianL.t0.boot
# Calculate bootstrap CIs
medians.boot.lcl <- apply(medians.boot, 2, quantile, probs=0.025)
medians.boot.ucl <- apply(medians.boot, 2, quantile, probs=0.975)
medians.boot.cis <- rbind(medians.boot.lcl, medians.boot.ucl)
# Transform into normal dates
temp <- data.frame(matrix(nrow=2, ncol=4))
for (j in 1:4) temp[,j] <- as.Date(medians.boot.cis[,j],
								 origin="1994-01-01")
medians.boot.cis2 <- cbind(temp, medians.boot.cis[,5:6])
colnames(medians.boot.cis2) <- colnames(medians.boot.cis)
# print Cis and p-values to be reported in the paper
print(medians.boot.cis2)
min(sum(medians.boot$differenceE<0),sum(medians.boot$differenceE>0))*2/nboot
min(sum(medians.boot$differenceL<0),sum(medians.boot$differenceL>0))*2/nboot

# Produce plot - Figure 2 of the paper
vispday.t0 <- aggregate(site~day1A, data=allvisits.t0, length)
vispday.t1 <- aggregate(site~day1A, data=allvisits.t1, length)
vispday.t0EL <- aggregate(site~day1A+eorl, data=allvisits.t0, length)
vispday.t1EL <- aggregate(site~day1A+eorl, data=allvisits.t1, length)
ymax <- max(vispday.t0$site, vispday.t1$site)

if(out=="screen") windows()
if(out=="png") png(paste0("results/graphs/", region, "/_survey dates.png"),
			 width=1600, height=900, pointsize=36)
if(out=="wmf") win.metafile(paste0("results/graphs/", region,
					     "/_survey dates.wmf"), width=6, height=4)
if(out=="eps") {
	setEPS(width=6, height=4, pointsize=10)
	postscript(paste0("results/graphs/", region, "/_survey dates.eps"))
} 

plot(vispday.t0$day1A, vispday.t0$site, xlim=range(intmodel),
     ylim=c(0,ymax), pch=1, axes=F, xlab="", ylab="Number of visits")
axis(1, at=c(0,30,61), labels=c("1/4", "1/5", "1/6"), pos=0)
segments(1,0,91,0)
axis(2, pos=0)
points(vispday.t1$day1A, vispday.t1$site, pch=16)
abline(v=medianE.t0-90, lty=2)
abline(v=medianL.t0-90, lty=2)
abline(v=medianE.t1-90, lty=1)
abline(v=medianL.t1-90, lty=1)

if(out!="screen") dev.off()

########
# Repeat everything but now with the subset of squares surveyed all years
sqsyrs <- unique(subset(allvisits, year%in%eyears | year%in%lyears,
				select=c(site, year)))
Nyrs <- aggregate(year~site, data=sqsyrs, FUN=length)
site.allyrs <- as.character(subset(Nyrs, year==10)$site)

allvisits.site.allyrs <- subset(allvisits, site %in% site.allyrs)
allvisits.site.allyrs.t0 <- subset(allvisits.site.allyrs, year%in% eyears)
allvisits.site.allyrs.t1 <- subset(allvisits.site.allyrs, year%in% lyears)
t0.E.sites <- allvisits.site.allyrs.t0$day1A[allvisits.site.allyrs.t0$eorl=="E"]
t1.E.sites <- allvisits.site.allyrs.t1$day1A[allvisits.site.allyrs.t1$eorl=="E"]
t0.L.sites <- allvisits.site.allyrs.t0$day1A[allvisits.site.allyrs.t0$eorl=="L"]
t1.L.sites <- allvisits.site.allyrs.t1$day1A[allvisits.site.allyrs.t1$eorl=="L"]
medianE.site.allyrs.t0 <- median(t0.E.sites) + 90
medianE.site.allyrs.t1 <- median(t1.E.sites) + 90
medianL.site.allyrs.t0 <- median(t0.L.sites) + 90
medianL.site.allyrs.t1 <- median(t1.L.sites) + 90
medians.site.allyrs <- data.frame(median.t0=c(medianE.site.allyrs.t0, medianL.site.allyrs.t0),
	     median.t1=c(medianE.site.allyrs.t1, medianL.site.allyrs.t1))
medians.site.allyrs$median.t0 <- as.Date(medians.site.allyrs$median.t0, origin="1994-01-01")
medians.site.allyrs$median.t1 <- as.Date(medians.site.allyrs$median.t1, origin="1994-01-01")

# Print medians to report in the paper
print(medians.site.allyrs)

# Number of squares used (report in the paper)
length(unique(allvisits.site.allyrs.t0$site))
length(unique(allvisits.site.allyrs.t1$site))

# Bootstrap over sites
nboot <- 2000
sites.allyears.t0 <- unique(allvisits.site.allyrs.t0$site)
sites.allyears.t1 <- unique(allvisits.site.allyrs.t1$site)
for (b in 1:nboot) {
	if(b%%500==0) print(b)
	bt0 <- sort(sample(sites.allyears.t0, length(sites.allyears.t0),
				 replace=T))
	bt1 <- sort(sample(sites.allyears.t1, length(sites.allyears.t1),
				 replace=T))
	bootsites.t0 <- data.frame(site=bt0)
	bootsites.t1 <- data.frame(site=bt1)
	boot.t0 <- merge(bootsites.t0, allvisits.site.allyrs.t0, by="site", all.x=T)
	boot.t1 <- merge(bootsites.t1, allvisits.site.allyrs.t1, by="site", all.x=T)
	boot.t0.E <- boot.t0$day1A[boot.t0$eorl=="E"]
	boot.t1.E <- boot.t1$day1A[boot.t1$eorl=="E"]
	boot.t0.L <- boot.t0$day1A[boot.t0$eorl=="L"]
	boot.t1.L <- boot.t1$day1A[boot.t1$eorl=="L"]
	medianE.t0.boot <- median(boot.t0.E) + 90
	medianE.t1.boot <- median(boot.t1.E) + 90
	medianL.t0.boot <- median(boot.t0.L) + 90
	medianL.t1.boot <- median(boot.t1.L) + 90
	temp.boot <- data.frame(medianE.t0.boot, medianL.t0.boot,
					medianE.t1.boot, medianL.t1.boot)
	medians.site.allyrs.boot <- if (b==1) temp.boot else
		rbind(medians.site.allyrs.boot, temp.boot)
}
# Calculate difference of median dates for early and late visits
medians.site.allyrs.boot$differenceE <- medians.site.allyrs.boot$medianE.t1.boot-
	medians.site.allyrs.boot$medianE.t0.boot
medians.site.allyrs.boot$differenceL <- medians.site.allyrs.boot$medianL.t1.boot-
	medians.site.allyrs.boot$medianL.t0.boot

# Calculate bootstrap CIs
medians.site.allyrs.boot.lcl <- apply(medians.site.allyrs.boot, 2, quantile, probs=0.025)
medians.site.allyrs.boot.ucl <- apply(medians.site.allyrs.boot, 2, quantile, probs=0.975)
medians.site.allyrs.boot.cis <- rbind(medians.site.allyrs.boot.lcl, medians.site.allyrs.boot.ucl)
# Transform into normal dates
temp <- data.frame(matrix(nrow=2, ncol=4))
for (j in 1:4) temp[,j] <- as.Date(medians.site.allyrs.boot.cis[,j],
								 origin="1994-01-01")
medians.site.allyrs.boot.cis2 <- cbind(temp, medians.site.allyrs.boot.cis[,5:6])
colnames(medians.site.allyrs.boot.cis2) <- colnames(medians.site.allyrs.boot.cis)
# print Cis and p-values to be reported in the paper
print(medians.site.allyrs.boot.cis2)
min(sum(medians.site.allyrs.boot$differenceE<0),sum(medians.site.allyrs.boot$differenceE>0))*2/nboot
min(sum(medians.site.allyrs.boot$differenceL<0),sum(medians.site.allyrs.boot$differenceL>0))*2/nboot

# Mood test
# (see https://stat.ethz.ch/pipermail/r-help/2010-April/234387.html)
median.test<-function(x,y){
    z<-c(x,y)
    g <- rep(1:2, c(length(x),length(y)))
    m<-median(z)
    fisher.test(z<m,g)
}
median.test(t0.E, t1.E)
median.test(t0.L, t1.L)
median.test(t0.E.sites, t1.E.sites)
median.test(t0.L.sites, t1.L.sites)


 
