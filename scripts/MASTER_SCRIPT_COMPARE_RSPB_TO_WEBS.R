##########################
### Compare historic RSPB counts to WeBS counts
##########################
## List all the RSPB files that were filtered for the three regions
## 

RSPBfiles <- list.files("Data/RSPB_data_filtered/",full.names = TRUE)
spp <- "eider"


rspb.data <- wrangle.RSPB.data(spp,save.output=T)
rspb.m <- merge.RSPB.to.WeBS(spp,rspb.data)
rspb.merged <- rspb.m$merged
WeBS.dat <- rspb.m$WeBS
rspb.v.webs.jan <- RSPB.v.WeBS.plot("JAN",rspb.merged)
rspb.v.webs.dec <- RSPB.v.WeBS.plot("DEC",rspb.merged)

## Run these to explore plots
rspb.v.webs.dec
rspb.v.webs.jan

## After you have run these above plots, see if any of them have a 
## r value > 0.6  If they do, then type in the region below
## as well as the month of interest.  
## Change year to 2019 if the month = DEC
## Change year to 2020 if the month = JAN
###############################################

## Nairn_Culbin_Bars, Outer_Dornoch_Firth, or Inverness_Beauly_Firth
site <- "Nairn_Culbin_Bars"
month <- "JAN"
year <- 2020
predict.to <- c(site = site,
                month = month,
                year = year)
######################################################
# This will generate predictions for 2020 for RSPB data
predict.rspb(WeBS.dat,rspb.merged,predict.to,month,site)

