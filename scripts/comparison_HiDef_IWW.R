### COMPARISON BETWEEN IWW DENSITIES AND HIDEF DENSITIES

setwd("C:/Workspace/Nextcloud/Github/HiDef/SNH_Moray_Firth/data")

geogr<-CRS("+init=epsg:4326")
main.crs = "+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"


# 1. Load data ------------------------------------------------------------

  spp <- "Common eider"  # Species subset ("Common eider", "Common scoter", "Red-throated diver")
  survey.date <- "2020-01-19"  # 1st aerial survey: "2020-01-19", 2nd aerial survey: "2020-03-08"

  # IWW DATA
  
    input.name <- paste0(spp,"_",survey.date,".RData") # IWW results were named as "spp_date.RData" 
    IWWData <- load(input.name)
  
    plot(buffer_ras)  # we are going to use the raster output


 # HiDef data
    spp.hd <- "Eider"  #In the HiDef database, "Common Eider" is called just "Eider"
    input.name <- paste0(spp.hd,"_survey01.rds")  #output was named "spp_survey01.rds" for results from the first survey
    hd.raster <- readRDS(input.name)
    plot(hd.raster)   # we check the raster



# 2. Make both rasters comparable -----------------------------------------



    crs(hd.raster) <- crs(buffer_ras)  # give same crs


    iww.res <- resample(buffer_ras, hd.raster)  # resample buffer_ras to get same extent and resolution than hd.raster
    
    # check resolution and extensions
      res(iww.res)  
      extent(iww.res)
      plot(iww.res)
      res(hd.raster)
      extent(hd.raster)

  # Mask raster hd.raster to get same area than iww.res
    hd.crop<- mask(hd.raster, iww.res)
 
    

# 3. Measure correlation --------------------------------------------------

    
  # extract values
    r1 <- values(hd.crop)
    r2 <- values(iww.res)
     
  # scatterplot
    plot(r1, r2)
    
  # Pearson?s correlation
    cor(r1, r2, use= "pairwise.complete.obs")

  # Linear model
     lm1 <- lm(r1 ~ r2)
     summary(lm1)  # significance of r2
     plot(lm1)
 