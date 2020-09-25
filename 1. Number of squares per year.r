############################################################################
# Calculate number of squares per year in rda7 (South-East England)        #
# This is needed for the chapter Methods/Analyses                          #
############################################################################

# Read bird data
dat <- readLines("X:/CensusUnit/Trends2019/BBS_birds/data/ebird_all")

# Prepare data frame with unique combinations of year and site
year <- substr(dat,1,4)
gridref <- substr(dat,9,14)
df <- data.frame(year, gridref)
dfu <- unique(df)

# Prepare data frame with list
rdasqs <- readLines("Z:/bbs/rdasqs")
gridref <- substr(rdasqs,1,6)
#country <- substr(rdasqs,16,16)
rda <- substr(rdasqs,32,33)
#rdadf <- data.frame(gridref, country, rda)
rdadf <- data.frame(gridref, rda)

# Merge data frames together 
dfurda <- merge(dfu, rdadf, by="gridref")

# Remove possible trailing space in the rda variable
dfurda$rda <- sub(" ", "", dfurda$rda)

# Only keep rda==7 (South-East England)
rda7 <- subset(dfurda, rda==7)

# Calculate number of squares surveyed per year
rda7.agg <- aggregate(gridref~year, rda7, length)

# Exclude 2001 as usual as BBS was not properly run (foot-and-mouth disease)
rda7.agg <- subset(rda7.agg, year!=2001)

# Range and mean to be reported in the paper
range(rda7.agg$gridref)
mean(rda7.agg$gridref)



