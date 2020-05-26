##############################################
### Create survey area maps ###
### Project HC0049
### Grant Humphries & Raul Vilela
### May 2020
### HiDef Aerial Surveying, BioConsult SH
##############################################

# Libraries ---------------------------------------------------------------
source("load_libs.R")
LIBS <- c("tidyverse","foreach","rgeos","sp","maptools","rgdal","raster","Hmisc","ggthemes",
          "Cairo","doSNOW","readxl","INLA","mapproj","inlabru")
Load_Libs(LIBS)


# Species select ----------------------------------------------------------
source('scripts/species_list.R')
SppName <- 'Common Scoter'

if(length(names(which(Species==SppName)))>0){
  Spp <- names(which(Species==SppName))
}else{
  print("WARNING - CHECK SPELLING OF SppName")
}



# Spatial Projections -------------------------------------------------------------
WGS84<-CRS("+init=epsg:4326")
UTM30<-CRS("+init=epsg:32630")


# Load Shapefiles ---------------------------------------------------------

boundary.shape <- readOGR(dsn="./Data/Shapefile/Moray_Firth_Area_WGS84.shp",
                          layer="Moray_Firth_Area_WGS84")
boundaryUTM30 <- spTransform(boundary.shape,UTM30)
boundary_DF <- fortify(boundaryUTM30)

Transects <- readOGR(dsn="./Data/Aerial_Surveys/2020 - Month 01 - Survey 01/Output/Zone87_M01_S01_20_Output/Zone87_M01_S01_20_Output-Day1-Transects.shp",
                     layer="Zone87_M01_S01_20_Output-Day1-Transects")
TransectsUTM30 <- spTransform(Transects,UTM30)
transdf <- data.frame(coordinates(TransectsUTM30))

load(file = 'D:/GIS_DATA/UK_coastline.rda')
coastline_sp_utm <- spTransform(tt, UTM30)
coastData_DF <- fortify(coastline_sp_utm)


SppDatajan <- readOGR(dsn="./Data/Aerial_Surveys/2020 - Month 01 - Survey 01/Output/Zone87_M01_S01_20_Output/Zone87_M01_S01_20_Output-Day1-CentCount.shp",
                   layer="Zone87_M01_S01_20_Output-Day1-CentCount")
SppDatajanUTM30 <- spTransform(SppDatajan,UTM30)
Sppdfjan <- data.frame(coordinates(SppDatajanUTM30))

SppDatamar <- readOGR(dsn="./Data/2020 - Month 03 - Survey 01/Output/Zone87_M03_S01_20_Output/Zone87_M03_S01_20_Output-Day1-CentCount.shp",
                      layer="Zone87_M03_S01_20_Output-Day1-CentCount")
SppDatamarUTM30 <- spTransform(SppDatamar,UTM30)
Sppdfmar <- data.frame(coordinates(SppDatamarUTM30))


# Generate transect map ---------------------------------------------------


xlims <- c(min(transdf$coords.x1)*0.99, max(transdf$coords.x1)*1.01)
ylims <- c(min(transdf$coords.x2)*0.999, max(transdf$coords.x2)*1.001)

G <- ggplot(transdf,aes(x=coords.x1,y=coords.x2)) +
  geom_polygon(aes(x = long, y = lat, group = group), data = coastData_DF, fill = "darkolivegreen3", 
               col = "grey35", alpha = 0.7) + 
  geom_polygon(aes(x=long,y=lat,group=group,fill="SPA"),data=boundary_DF,
               col="black",alpha=0.6)+
  scale_fill_manual(name="",values=c("SPA" = "lightblue"))+
  geom_point(size=0.8,pch=21,fill="orange")+
  scale_x_continuous(breaks=seq(425000,525000,25000),label=seq(4.25,5.25,.25))+
  scale_y_continuous(breaks=seq(6380000,6440000,20000),label=seq(6.38,6.44,.02))+
  xlab(expression(X(meters~x~10^{"6"})))+
  ylab(expression(Y(meters~x~10^{"7"})))+
  coord_equal(xlim = xlims, ylim = ylims)+
  
  ggsn::scalebar(boundary_DF,location = "bottomright",st.size = 5,height = 0.04,dist = 20,dist_unit = "km",transform = FALSE)+
  ggsn::north(boundary_DF,scale = 0.15,location = "topleft")+
  ggtitle("Transects")+
  theme_gdocs()+
  theme(
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black"),
    plot.background = element_blank()
  )

G

ggsave(plot=G,filename="Transects_Map.jpeg",device = "jpeg",width=10,height=10)



# Generate species map ----------------------------------------------------


source('scripts/species_list.R')
SppName <- 'Long-tailed Duck'

if(length(names(which(Species==SppName)))>0){
  Spp <- names(which(Species==SppName))
}else{
  print("WARNING - CHECK SPELLING OF SppName")
}



Sppdfjan$spp <- as.vector(SppDatajanUTM30@data[,Spp]) 
xlims <- c(min(transdf$coords.x1)*0.99, max(transdf$coords.x1)*1.01)
ylims <- c(min(transdf$coords.x2)*0.999, max(transdf$coords.x2)*1.001)

G <- ggplot(Sppdfjan,aes(x=coords.x1,y=coords.x2)) +
  #geom_polygon(aes(x = long, y = lat, group = group), data = coastData_DF, fill = "darkolivegreen3", 
  #             col = "grey35", alpha = 0.7) +
  geom_polygon(aes(x=long,y=lat,group=group,fill="SPA"),data=boundary_DF,
               col="black",alpha=0.6)+
  scale_fill_manual(name="",values=c("SPA" = "lightblue"))+
  geom_point(aes(size=spp),pch=21,fill="orange")+
  scale_size_continuous(range=c(0.5,8),name="Count")+
  
  scale_x_continuous(breaks=seq(425000,525000,25000),label=seq(4.25,5.25,.25))+
  scale_y_continuous(breaks=seq(6380000,6440000,20000),label=seq(6.38,6.44,.02))+
  xlab(expression(X(meters~x~10^{"6"})))+
  ylab(expression(Y(meters~x~10^{"7"})))+
  coord_equal(xlim = xlims, ylim = ylims)+
  ggsn::scalebar(boundary_DF,location = "bottomright",st.size = 5,height = 0.04,dist = 20,dist_unit = "km",transform = FALSE)+
  ggsn::north(boundary_DF,scale = 0.15,location = "topleft")+
  ggtitle(SppName)+
  theme_gdocs()+
  theme(
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black"),
    plot.background = element_blank()
  )


G

ggsave(plot=G,filename="Output.jpeg",device = "jpeg",width=10,height=10)



