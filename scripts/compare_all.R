
spp <- "Common scoter"

iww.pgons <- outputs$locs.dens

compare.all.data(iww.pgons,spp)

compare.all.data <- function(iww.pgons,spp){
  month <- "01"  
  WeBS.data <- readRDS("Data/WeBS/WeBS_wrangled_20May2020.rds")
  input.name <- paste0("outputs/",spp,"_month",month,".rds")  
  hd.raster <<- readRDS(input.name)
  pgons <- iww.pgons
  IWW.df <- pgons@data
  
  ## A list of sites to extract from WeBS to match RSPB
  sites.to.extract <- c("Rosemarkie to Avoch","Avoch Bay",
                        "Easterton - Fort George","Castle Stuart - Westerton",
                        "Craigiehowe to Craigton", "South Kessock",
                        "Culbin Bar","Nairn Bar","Embo","Golspie")
  
  
  ### First, get the HiDef data (predicted models) summarized in the IWW polygons..
  extracted <- raster::extract(hd.raster,pgons)
  dfout <- data.frame(hd.data=unlist(lapply(extracted,function(x) mean(x,na.rm=T))),
                      iww.data=pgons@data$dens,
                      iww.rate=pgons@data$Number/pgons@data$Time_obs,
                      iww.site=pgons@data$sector,
                      iww.area=pgons@data$area.x)
  
  MEANDUR.IWW <- mean(pgons@data$Time_obs)
  
  WeBS <- WeBS.data %>% dplyr::filter(SECTOR_NAME %in% sites.to.extract)
  
  WeBS$SECTOR_NAME_IWW <- plyr::revalue(WeBS$SECTOR_NAME,
                                        c("South Kessock"="Beauly Firth"))
  
  WeBS$SECTOR_NAME_RSPB <- plyr::revalue(WeBS$SECTOR_NAME_IWW,
                                         c("Rosemarkie to Avoch" = "Inverness_Beauly_Firth",
                                           "Avoch Bay" = "Inverness_Beauly_Firth",
                                           "Easterton - Fort George" = "Inverness_Beauly_Firth",
                                           "Castle Stuart - Westerton" = "Inverness_Beauly_Firth",
                                           "Craigiehowe to Craigton" = "Inverness_Beauly_Firth",
                                           "Beauly Firth" = "Inverness_Beauly_Firth",
                                           "Nairn Bar" = "Nairn_Culbin_Bars",
                                           "Culbin Bar" = "Nairn_Culbin_Bars",
                                           "Embo" = "Outer_Dornoch_Firth",
                                           "Golspie" = "Outer_Dornoch_Firth"))
  
  ## We can only compare to January here as there are no data in March
  WeBS.df <- WeBS %>% 
    unnest(data) %>%
    dplyr::filter(tolower(SPECIES) == tolower(spp),
                  YEAR == 2020, MONTH == "01")
  
  ## Using the mean survey duration for imputation of effort where not available
  WeBS.df$DURATION_HRS[is.na(WeBS.df$DURATION_HRS)] <- mean(WeBS.df$DURATION_HRS,na.rm=T)
  
  MEANDUR.WEBS <- mean(WeBS.df$DURATION_HRS,na.rm=T)
  
  WeBS.df$BIRD_RATE <- WeBS.df$BIRD_COUNT/WeBS.df$DURATION_HRS
  
  WeBS.summarized <- WeBS.df %>% 
    group_by(SECTOR_NAME_RSPB) %>%
    dplyr::summarise(WeBS.rate = mean(BIRD_RATE))
  
  ## Now extract IWW
  IWW.out <- dfout %>% dplyr::filter(iww.site %in% sites.to.extract)
  
  
  
  sites.to.extract <- c("Rosemarkie to Avoch","Avoch Bay",
                        "Easterton - Fort George","Castle Stuart - Westerton",
                        "Craigiehowe to Craigton", "South Kessock",
                        "Culbin Bar","Embo","Golspie to Brora")
  
  
  
  IWW.out$SECTOR_NAME_RSPB <- plyr::revalue(IWW.out$iww.site,
                                            c("Rosemarkie to Avoch" = "Inverness_Beauly_Firth",
                                              "Easterton - Fort George" = "Inverness_Beauly_Firth",
                                              "Castle Stuart - Westerton" = "Inverness_Beauly_Firth",
                                              "Craigiehowe to Craigton" = "Inverness_Beauly_Firth",
                                              "Culbin Bar" = "Nairn_Culbin_Bars",
                                              "Embo" = "Outer_Dornoch_Firth",
                                              "Golspie to Brora" = "Outer_Dornoch_Firth"))
  
  
  IWW.hd.sum <- IWW.out %>% 
    group_by(SECTOR_NAME_RSPB) %>%
    dplyr::summarise(hd.mean = mean(hd.data),
                     iww.mean.rate = mean(iww.rate),
                     iww.area = sum(iww.area))
  
  merged <- merge.data.frame(IWW.hd.sum,WeBS.summarized,by="SECTOR_NAME_RSPB")
  return(merged)
  
}




