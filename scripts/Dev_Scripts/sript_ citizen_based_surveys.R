
setwd("C:/Workspace/Nextcloud/Github/HiDef/SNH_Moray_Firth/data")
library(rgdal)
library(readxl)
library(raster)
library(sf)
library(plyr)

geogr<-CRS("+init=epsg:4326")
main.crs = "+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"

clean_units <- function(x){
  attr(x,"units") <- NULL
  class(x) <- setdiff(class(x),"units")
  x
}

# 1. Load data ------------------------------------------------------------

  # accesory layers
    mask <- readOGR(dsn = "Data/Shapefile", layer = "Moray_Firth_Area_WGS84")
    mask<-spTransform(mask, main.crs)

  # We make a raster mask...
  # not really needed, but it might be useful (...who knows)
    extension <- extent(mask)
    raster.mask <- raster()
    extent(raster.mask) <- extension
    crs(raster.mask) <- main.crs
    res(raster.mask) <- c(1,1)
    raster.mask <- rasterize(mask, raster.mask,1)
    raster.mask[is.na(raster.mask[])] <- 0 

  # sighting locations

    locs <- read_excel("Data/VPs_sectors_coordinates_rv.xlsx",sheet=1)

      # Aggregate WeBS sectors inside the Beauly Firth (overlapping areas between both sides of the firth)

        locs$sector <- locs$WeBS_sector
        locs$sector<- revalue(locs$WeBS_sector, c("Lentran to Bunchrew"="Beauly Firth", "Beauly Firth North"="Beauly Firth",
                                                  "South Kessock"="Beauly Firth", "Ness Mouth"="Beauly Firth"))
    
        coordinates(locs) <- c("Lon","Lat")
        crs(locs) <- geogr
        locs<-spTransform(locs, main.crs)
        plot(locs)


  # Polygons with the WeBS areas
    polygs<- readOGR(dsn = "Data/Shapefile", layer = "Banff_to_Helmsdale_WeBS")
    polygs<-spTransform(polygs, main.crs)
    names(polygs)
    sort(unique(polygs$NAME))


      # merge polygons within the Beauly Firth
        polygs$NAME<- revalue(polygs$NAME, c("Lentran"="Beauly Firth","Lentran to Bunchrew"="Beauly Firth", "Beauly Firth North"="Beauly Firth",
                                                "South Kessock"="Beauly Firth", "Ness Mouth"="Beauly Firth"))

        polygs <- aggregate(polygs, by= "NAME", dissolve=T)

        plot(polygs)
        plot(locs, add=TRUE, pch=20, col="red")
    
        spplot(polygs, ycol="NAME")  # looks good



  # iww monitoring project observations

    iww <- read_excel("Data/IWW_data/IWW_2020_rv.xlsx", sheet = 2)
    names(iww)
    head(iww)
    unique(iww$taxonName)
    sort(unique(iww$VP_code))
    sort(unique(locs$VP_code))

    setdiff(iww$VP_code,locs$VP_code) # all codes in iww are in locs
    setdiff(locs$VP_code,iww$VP_code)  # not all codes in locs are in iww, but that is ok
    

# 2. Subset IWW data by date ----------------------------------------------

    # available dates: "2020-01-19", "2020-01-21", "2020-03-08", "2020-03-09"
    # 1st aerial survey: "2020-01-19", 2nd aerial survey: "2020-03-08"
    
    survey.date <- "2020-01-19"
  
    iww$Date <- as.character(iww$Date)
    iww <- iww[iww$Date == survey.date, ]
    iww <- iww[ c(3,4,8,13,20, 21) ]
    names(iww)

  # subset only locs with IWW data

    uniq.iww.vpcodes <- unique(iww$VP_code) # only for sectors with IWW VPs
    locs.iww <- subset(locs, VP_code %in% uniq.iww.vpcodes )


# 3. Calculate buffer 2km and mask ----------------------------------------

    locs.sf = st_as_sf(locs.iww)
    mask.sf <- st_as_sf(mask)
    
    locs.buffer = st_buffer(locs.sf, 2)
    locs.sea = st_intersection(locs.buffer, mask.sf)
    locs.sea$area = st_area(locs.sea)
    locs.sea <- as(locs.sea,"Spatial")


# 4. Merge datasets -------------------------------------------------------


  # Assign coordinates to IWW locations based on VP code
    iww.locs <- merge(iww, locs.sea, by.x="VP_code",by.y=c("VP_code"))
    names(iww.locs)
    unique(iww.locs$VP_code)
    
  # Unique list of dates and species
    unique(iww.locs$taxonName)
    unique(iww.locs$Date)

  # subset WeBS sectors with IWW observations inside

    uniq.iww.sectors <- unique(iww.locs$sector)
    polygs.iww <- subset(polygs, NAME %in% uniq.iww.sectors )
    
  # Plot
    plot(locs.sea, col="blue")  # buffers
    plot(mask, add=TRUE)   # polygon mask Moray Firth
    plot(locs.iww, add=TRUE, col="red", pch=20) # VPs in the IWW data
    plot(polygs.iww, add=TRUE, col="transparent", border="orange") # WebS sectors in the IWW data

  # clean tables
    locs.sea$area <- clean_units(locs.sea$area)
    locs.sea <- locs.sea[ -c(2,3,4,5,6,9,10) ]
    names(locs.sea)

  # first we take the first line of each VP code (several groups of the same spp can be sighted by VP)
  # we do it before subsetting the species to include effort in VPs with zero observations
    iww.aggr <- ddply(iww.locs, "sector", head, 1) 
    names
    iww.aggr <- iww.aggr[ c(4,13) ]  # TIME and UNIQUE VPs (with or without observation) 

  # Species subset ("Common eider", "Common scoter", "Red-throated diver")
    spp <- "Common eider"

    iww.sub <- iww.locs[iww.locs$taxonName == spp, ]
    iww.sub <- iww.sub[complete.cases(iww.sub[ , 1]),] # I do not know why the previous line generates an NA, we remove it

  # aggregate average of animals by sector, spp and date

    counts <- aggregate(Number~sector, iww.sub, mean)  # mean of individuals

  # Area

    locs.sector <- aggregate(locs.sea, by= "sector", dissolve=T)
    locs.sector.sf = st_as_sf(locs.sector)
    locs.sector.sf$area = st_area(locs.sector.sf)
    locs.sector <- as(locs.sector.sf,"Spatial")
    locs.sector$area <- clean_units(locs.sector$area)


  # MERGE all (counts, areas and time)
    # counts are in data.frame counts (only for sectors with positive counts)
    # areas are in locs.sector (for all sectors)
    # time is in iww.aggr  (for all sectors)

    iww.all <- merge(iww.aggr, locs.sector, by="sector")
    iww.all <- merge(iww.all, counts, by="sector",all.x=TRUE)
    iww.all$Number[is.na(iww.all$Number)] <- 0

  # calculate density by sector (counts/ hour* sq.km  by sector)
    
    iww.all$dens <- with(iww.all, Number/ (Time_obs * area))
    

  # merge with polygons
    polygs.dens <- merge(polygs, iww.all, by.x="NAME", by.y="sector", all.x=TRUE) 
    spplot(polygs.dens, zcol= "dens", colorkey=FALSE) 

  # merge with buffers
    locs.dens <- merge(locs.sector, iww.all, by="sector", all.x=TRUE) 
    spplot(locs.dens, zcol= "dens", colorkey=FALSE) 




# 5. Outputs --------------------------------------------------------------


  # convert buffers with density to raster

    r <- raster(ncol=1000, nrow=1000)
    extent(r) <- extent(locs.dens)

    buffer_ras<- rasterize(locs.dens, r, "dens", fun="mean")
    plot(buffer_ras)
    
  # plot
    plot(mask, add=TRUE)   # polygon mask Moray Firth
    plot(locs.iww, add=TRUE, col="red", pch=20) # VPs in the IWW data
    plot(polygs.iww, add=TRUE, col="transparent", border="orange") # WebS sectors in the IWW data
    
  # save outputs
    output.name <- paste0(spp,"_",survey.date,".RData")
    save(polygs.dens, locs.dens, buffer_ras, file= output.name)
