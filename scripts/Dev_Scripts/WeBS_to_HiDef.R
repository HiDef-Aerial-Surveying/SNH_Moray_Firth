spp <- "common scoter"
month <- "01"
Webspolys<- readOGR(dsn = "Data/Shapefile", layer = "Banff_to_Helmsdale_WeBS")
Webspolys<-spTransform(Webspolys, UTM30)
#polygs$NAME<- plyr::revalue(polygs$NAME, c("Lentran"="Beauly Firth","Lentran to Bunchrew"="Beauly Firth", "Beauly Firth North"="Beauly Firth",
#                                           "South Kessock"="Beauly Firth", "Ness Mouth"="Beauly Firth"))
#polygs<<- aggregate(polygs, by= "NAME", dissolve=T)

WeBS.data <- readRDS("Data/WeBS/WeBS_wrangled_20May2020.rds")

tt <- merge.data.frame(Webspolys@data,WeBS.data,by.x="NAME",by.y = "SECTOR_NAME")

tt <- tt %>% 
  unnest(data) %>%
  dplyr::filter(tolower(SPECIES) == tolower(spp),
                YEAR == 2020, MONTH == "01")


ts <- merge.data.frame(tt,Webspolys@data,by="NAME",all.y=T)
Webspolys@data$BIRD_COUNT <- ts$BIRD_COUNT
Webspolys@data$BIRD_DENSITY <- ts$BIRD_DENSITY



input.name <- paste0("outputs/",spp,"_month",month,".rds")  
hd.raster <<- readRDS(input.name)


extracted <- raster::extract(hd.raster,Webspolys)
dfout <- data.frame(hd.data=unlist(lapply(extracted,function(x) mean(x,na.rm=T))),
                    WeBS.count=Webspolys@data$BIRD_COUNT,
                    WeBS.dens=Webspolys@data$BIRD_DENSITY)



