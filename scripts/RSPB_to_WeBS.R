### Wrangle the RSPB data
source("scripts/helpers.R")

RSPBfiles <- list.files("Data/RSPB_data_filtered/",full.names = TRUE)
spp <- "eider"
site <- "Outer_Dornoch_Firth"
month <- "DEC"
year <- 2019
predict.to <- c(site = site,
                month = month,
                year = year)


rspb.data <- wrangle.RSPB.data(spp,save.output=T)
rspb.m <- merge.RSPB.to.WeBS(spp,rspb.data)
rspb.merged <- rspb.m$merged
WeBS.dat <- rspb.m$WeBS
rspb.v.webs.eider.jan <- RSPB.v.WeBS.plot("JAN",rspb.merged)
rspb.v.webs.eider.dec <- RSPB.v.WeBS.plot("DEC",rspb.merged)


rspb.merged$YEAR <- as.numeric(rspb.merged$YEAR)
mod1 <- with(rspb.merged,
             glm(rspb.rate ~ mean.rate,family="gaussian")
             )
summary(mod1)
plot(mod1)

predict.rspb(WeBS.dat,rspb.merged,predict.to,month,site)









      