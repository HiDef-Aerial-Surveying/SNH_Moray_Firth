### Code to open and wrangle WeBS datasets ###
### Project HC0049
### Grant Humphries & Raul Vilela
### May 2020
### HiDef Aerial Surveying, BioConsult SH


# Libraries ---------------------------------------------------------------
source("scripts/load_libs.R")
LIBS <- c("tidyverse","foreach","readxl","stringi",
          "rgdal","ggthemes","raster")
Load_Libs(LIBS)



# Functions ---------------------------------------------------------------


load_webs_sheet <- function(x,y){
  df <- vroom::vroom(x)
  df$Species <- plyr::revalue(df$Species,c("Eider (except Shetland)"="Eider"))
  df <- df[df$Species %in% Species,]
  df <- df[df$Count > 0,]
  if(nrow(df)>0){
    df$VisitDate <- NA
    df$VisitDate[df$Visit=="OCT-2019"] <- "2019-10-01"
    df$VisitDate[df$Visit=="NOV-2019"] <- "2019-11-01"
    df$VisitDate[df$Visit=="DEC-2019"] <- "2019-12-01"
    df$VisitDate[df$Visit=="JAN-2020"] <- "2020-01-01"
    df$VisitDate[df$Visit=="FEB-2020"] <- "2020-02-01"
    df$VisitDate[df$Visit=="MAR-2020"] <- "2020-03-01"
    
    df$VisitDate <- as.POSIXct(df$VisitDate)
    
    df <- df %>% arrange(VisitDate)
    
    WeB_sheet <- All_WeBS %>% dplyr::filter(SECTOR_NAME == y[stri_detect_fixed(x,y)])
    
    if(nrow(WeB_sheet)==0){
      print(x)
      print("uh oh... no data!")
    }
    
    outdf <- data.frame(
      YEAR = format(as.Date(df$VisitDate,format="%Y-%m-%d"),"%Y"),
      SECTOR_CODE = rep(WeB_sheet[1,]$SECTOR_CODE,nrow(df)),
      SECTOR_NAME = rep(WeB_sheet[1,]$SECTOR_NAME,nrow(df)),
      SITE_CODE = rep(WeB_sheet[1,]$SITE_CODE,nrow(df)),
      SITE_NAME = rep(WeB_sheet[1,]$SITE_NAME,nrow(df)),
      SECTORCENTREGRID = rep(WeB_sheet[1,]$SECTORCENTREGRID,nrow(df)),
      COUNT_DATE = df$VisitDate,
      WEBS_METHOD = rep(WeB_sheet[1,]$WEBS_METHOD,nrow(df)),
      TAXONOMIC_ORDER = rep(99999,nrow(df)),
      SPECIES = df$Species,
      BIRD_COUNT = df$Count
    )
    return(outdf)
  }
}


checkGrid <- function(x){
  unique(x)[which(!is.na(unique(x)))]
}

# Load WeBS dataset -----------------------------------------------------------
source("scripts/species_list.R")
All_WeBS <- readxl::read_xlsx("Data/WeBS/Marine pSPA Seaduck WeBS second extract.xlsx")
All_WeBS <- All_WeBS[All_WeBS$SPECIES %in% Species,]


WeBS_shape <- readOGR("./Data/Shapefile/Banff_to_Helmsdale_WeBS.shp",
                      layer="Banff_to_Helmsdale_WeBS",verbose = F)



# Load up and append sheets -----------------------------------------------

fs <- list.files(pattern="*-1YrsTo2019.csv",path="./Data/Webs/",full.names = T)

y <- unique(All_WeBS$SECTOR_NAME)

out <- foreach(x=fs,.combine='rbind')%do%{
  tt <- load_webs_sheet(x,y)
  print(tt)
}
        

finalmerge <- rbind(All_WeBS,out)
finalmerge <- finalmerge %>% arrange(SECTOR_NAME,YEAR)
saveRDS(finalmerge,"WeBS_merged_to_2020.rds")

# Expand merged data to include all species -------------------------------

finalmerge <- readRDS("WeBS_merged_to_2020.rds")
spp_in_finalmerge <- unique(finalmerge$SPECIES)


expanded <- finalmerge %>% group_by(SECTOR_NAME,COUNT_DATE,SECTOR_CODE,SECTORCENTREGRID)%>%
  complete(SPECIES = spp_in_finalmerge) %>%
  mutate(
    YEAR = format(as.Date(COUNT_DATE,format="%Y-%m-%d"),"%Y"),
    MONTH = format(as.Date(COUNT_DATE,format="%Y-%m-%d"),"%m")#,
    ) %>% 
  dplyr::select(SECTOR_NAME,SECTOR_CODE,SECTORCENTREGRID,COUNT_DATE,YEAR,MONTH,SPECIES,BIRD_COUNT)

expanded$BIRD_COUNT[is.na(expanded$BIRD_COUNT)] <- 0  
  

saveRDS(expanded,"WeBS_expanded_data.rds")


# Assess WeBS data --------------------------------------------------------
expanded <- readRDS("WeBS_expanded_data.rds")

assess_exp <- expanded %>% 
  nest(data = c(COUNT_DATE,YEAR,MONTH,SPECIES,BIRD_COUNT)) %>%
  mutate(
    TOTAL_YRS =map_dbl(data,~length(unique(.x$YEAR)))
  )
names(assess_exp)[grep("SECTORCENTREGRID",names(assess_exp))] <- "GRIDREF"
assess_exp <- assess_exp %>% ungroup()
WScopy <- WeBS_shape


tt <- WScopy@data
tomerge <- assess_exp %>%dplyr::select(SECTOR_NAME,GRIDREF,TOTAL_YRS)

outtt <- tomerge %>% right_join(tt,copy = T) %>% dplyr::filter(!is.na(NAME))

WScopy@data$TOTAL_YRS <- outtt$TOTAL_YRS



# Mapping WeBS data -------------------------------------------------------
WGS84<-CRS("+init=epsg:4326")
UTM30<-CRS("+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

Transects <- readOGR(dsn="./Data/Aerial_Surveys/2020 - Month 01 - Survey 01/Output/Zone87_M01_S01_20_Output/Zone87_M01_S01_20_Output-Day1-Transects.shp",
                     layer="Zone87_M01_S01_20_Output-Day1-Transects")
TransectsUTM30 <- spTransform(Transects,UTM30)
transdf <- data.frame(coordinates(TransectsUTM30))

WeBS.poly.UTM30 <- spTransform(WScopy,UTM30)
WeBS.poly.UTM30@data$id <- rownames(WeBS.poly.UTM30@data)
WeBS.fort <- fortify(WeBS.poly.UTM30,region="id")


WeBS.fort.df <- merge(WeBS.fort,WeBS.poly.UTM30@data,by="id")


xlims <- c(min(transdf$coords.x1)*0.99, max(transdf$coords.x1)*1.01)
ylims <- c(min(transdf$coords.x2)*0.999, max(transdf$coords.x2)*1.001)

G <- ggplot(transdf,aes(x=coords.x1,y=coords.x2)) +
  geom_polygon(aes(x=long,y=lat,group=group,fill=TOTAL_YRS),data=WeBS.fort.df,
               col="black",alpha=0.6)+
  scale_fill_gradient2(high="red",
                      na.value="white")+
  #scale_fill_manual(name="",values=c("SPA" = "lightblue"))+
  #geom_point(size=0.8,pch=21,fill="orange")+
  scale_x_continuous(breaks=seq(425000,525000,25000),label=seq(4.25,5.25,.25))+
  scale_y_continuous(breaks=seq(6380000,6440000,20000),label=seq(6.38,6.44,.02))+
  xlab(expression(X(meters~x~10^{"6"})))+
  ylab(expression(Y(meters~x~10^{"7"})))+
  coord_equal(xlim = xlims, ylim = ylims)+
  #annotation_sca
  #ggsn::scalebar(WeBS.fort,location = "bottomright",st.size = 5,height = 0.04,dist = 20,dist_unit = "km",transform = FALSE)+
  #ggsn::north(WeBS.fort,scale = 0.15,location = "topleft")+
  ggtitle("Total number of years in each WeBS sector")+
  theme_gdocs()+
  theme(
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black"),
    plot.background = element_blank()
  )

  

# Merge effort data -------------------------------------------------------

  
expanded <- readRDS("WeBS_expanded_data.rds")

effort <- readxl::read_xlsx("Data/WeBS/Marine pSPA WeBS duration - Effort data.xlsx")  
names(effort)[2] <- "SECTOR_CODE"

effort$MONTH <- format(effort$NOMINALDATE,"%m")
effort$YEAR <- format(effort$NOMINALDATE,"%Y")
effort$id <- paste0(effort$MONTH,effort$YEAR,as.character(effort$SECTOR_CODE))
expanded$id <- paste0(expanded$MONTH,expanded$YEAR,as.character(expanded$SECTOR_CODE))

joined_dat <- effort %>% right_join(expanded,by="id",copy=T) %>%
  dplyr::select(id,SECTOR_NAME,SECTOR_CODE.y,YEAR.y,MONTH.y,DURATION_HRS,
                SPECIES,BIRD_COUNT)
names(joined_dat)[3:5] <- c("SECTOR_CODE","YEAR","MONTH")

## Polygon WeBS area 
area_sqkm <- raster::area(WeBS_shape)/1000000
df.to.merge <- tbl_df(data.frame(SECTOR_CODE=WeBS_shape@data$LOC_LABEL,AREA = area_sqkm))


joined_area <- df.to.merge %>% right_join(joined_dat,by="SECTOR_CODE")


joined_area$BIRD_DENSITY <- joined_area$BIRD_COUNT/joined_area$AREA



assess_exp <- 
  joined_area %>% 
  nest(data = c(id,YEAR,MONTH,DURATION_HRS,SPECIES,BIRD_COUNT,BIRD_DENSITY))

saveRDS(assess_exp,"WeBS_wrangled_20May2020.rds")




