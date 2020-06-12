source("scripts/helpers.R")
source("scripts/species_list.R")
#View(Species)
spp <- "common scoter"
#month <- "01"
month <- "03"
#survey.date <- "2020-01-19"  
survey.date <- "2020-03-08"  

# Load data ------------------------------------------------------------
data.loaded <- load_month_data(month,spp)



# Create a mesh (fishnet) -------------------------------------------------
MESH <- create.mesh(mask = mask,
                    MaxEdge = 4,
                    track=data.loaded$track,
                    obs.model = data.loaded$obs.model)

# Prepare the model -------------------------------------------------------
model.data <- model.prep(obj.list = data.loaded,
                         jitter.points = TRUE, 
                         MESH = MESH,
                         priorrange = c(15,0.5),
                         priorsigma = c(1,0.5))

# Run the spde model ------------------------------------------------------
lgcp.model <- lgcp.spde(modprep.list = model.data)

# Create model predictions in space ---------------------------------------
predictions <- model.predictions(MESH = MESH,lgcp.fit = lgcp.model)

# Save the predictions to the output folder
save.predictions(predicted=predictions$predicted,month = month,spp=spp)

##
#prediction.raster <- readRDS("outputs/Common scoter_month01.rds")

# Plot the predicted distribution -----------------------------------------
plot.abundance(predicted=predictions$predicted,track=data.loaded$track,
               plot.title = "Density of Common Scoter (January 2020)",
               plot.boundary = T,plot.coastline = T,round.scale = 0)

# Plot population estimate ------------------------------------------------
plot.pop.priors(pred.list = predictions)

# Get population estimate numerically -------------------------------------
predictions$marginals

# Wrangling the IWW data --------------------------------------------------

spp <- "common scoter"
month <- "01"
#month <- "03"
survey.date <- "2020-01-19"
#survey.date <- "2020-03-08"


load_VP_data()
locs.iww <- subset_IWW_date(survey.date = survey.date)
locs.sea <- calculate.buffer(locs.iww = locs.iww,
                             mask = mask,
                             buffer.size=2)

iww$taxonName <- plyr::revalue(iww$taxonName, c("Common eider"="Eider"))
merge.iww_to_WeBS_polys(iww = iww,
                        locs.sea = locs.sea,
                        polygs = polygs)
plot.basic.view(locs.sea = locs.sea,
                mask = mask,
                locs.iww = locs.iww,
                polygs.iww = polygs.iww)

counts <- subset_species(iww.locs = iww.locs,spp = spp)
outputs <- create.polygons(locs.sea,
                           iww.aggr,
                           counts,
                           save.output = T,
                           spp = spp,
                           survey.date = survey.date)

plot.title = "Red-Throated Diver IWW counts (March 2020)"
plot.iww.pgons(pgons=outputs$locs.dens,track=data.loaded$track,
               plot.title=plot.title)



# Comparing IWW data to HiDef flights -------------------------------------
load.rasters(spp,survey.date,month)

# Comparing hidef raster to vantage point raster --------------------------
valuesout <- compare.rasters(hd.raster,outputs$buffer_ras,log.scale=F)

# Comparing Hidef raster to vantage point polygons ------------------------
psout <- compare.polygons(pgons = outputs$locs.dens,
                          hd.raster,
                          log.scale = F)


ggsave(plot=psout$Gplt,filename="outputs/Common_scoter_iww_v_Hidef_Mar.jpeg",device="jpeg",
       width=8,height=6)
       

# Comparing a subset of polygons to the raster ----------------------------
print(outputs$locs.dens@data$sector)

site.list <- c("Embo","Culbin Bar","Burghead to Hopeman",
               "Lossie Mouth to Spey Mouth","Dornoch Firth",
               "Tarbat Ness to Portmahomack","Wilkhaven to Rockfield")
pgonsout <- extract.polygons(outputs$locs.dens,site.list)
subsetcompare <- compare.polygons(pgons = pgonsout,
                                  hd.raster = hd.raster,
                                  log.scale = T)

# Compare IWW to WeBS -----------------------------------------------------
# Can only do this for January... 
WeBS.To.IWW <- webs.to.iww(spp,log.scale=F)

WeBS.To.IWW
ggsave(plot = WeBS.To.IWW,
       filename="outputs/scoter_IWWvWeBS.jpeg",
       device="jpeg",width=8,height=6)

# Compare RSPB to WeBS ----------------------------------------------------
RSPBfiles <- list.files("Data/RSPB_data_filtered/",full.names = TRUE)
spp <- "eider"

## Nairn_Culbin_Bars, Outer_Dornoch_Firth, or Inverness_Beauly_Firth
site <- "Nairn_Culbin_Bars"
month <- "JAN"
year <- 2020
predict.to <- c(site = site,
                month = month,
                year = year)


rspb.data <- wrangle.RSPB.data(spp,save.output=T)
rspb.m <- merge.RSPB.to.WeBS(spp,rspb.data)
rspb.merged <- rspb.m$merged
WeBS.dat <- rspb.m$WeBS
rspb.v.webs.jan <- RSPB.v.WeBS.plot("JAN",rspb.merged)
rspb.v.webs.dec <- RSPB.v.WeBS.plot("DEC",rspb.merged)


rspb.merged$YEAR <- as.numeric(rspb.merged$YEAR)
mod1 <- with(rspb.merged,
             glm(rspb.rate ~ mean.rate,family="gaussian")
)
summary(mod1)
plot(mod1)

predict.rspb(WeBS.dat,rspb.merged,predict.to,month,site)



# Compare all data --------------------------------------------------------

spp <- "Common eider"

iww.pgons <- outputs$locs.dens

compare.all.data(iww.pgons,spp)







