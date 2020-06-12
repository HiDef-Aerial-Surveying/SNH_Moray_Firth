### Moray Firth Modelling
### Mapping the distribution of seabirds from digital aerial data
### Code by Raul Vilela and Grant Humphries
### Written for Scottish Natural Heritage
### May, 2020


source("scripts/helpers.R")
# Select Species ----------------------------------------------------------
## Load the species list and view it
source("scripts/species_list.R")
View(Species)

## Exchange "Common scoter" with your species of choice
spp <- "Common scoter"

# Select month ------------------------------------------------------------

## Select the month you wish to model in MM format as a CHARACTER (string)
## I.E. Jan = "01", Feb = "02", ... Dec = "12"
month <- "01"

# Select the date of the survey -------------------------------------------

# 1st aerial survey: "2020-01-19", 2nd aerial survey: "2020-03-08"
survey.date <- "2020-01-19"  


# Load data ------------------------------------------------------------
## This will load up the track data
# month = the month of interest as a character string ("01" or "03")
# spp = the common name of the species. Refer to species list loaded above. 
data.loaded <- load_month_data(month,spp)



# Create a mesh (fishnet) -------------------------------------------------

# This mesh is used to integrate information on survey effort and observations
# to generate the voronoi polygons. The final predictions are also made in this fishnet/mesh
# You only have to run this one time and save it out to load in where needed. 
# e.g.,   MESH <- create.mesh(...)
#         saveRDS('mesh1.rds')
#         mesh1 <- readRDS('mesh1.rds)
# mask = the boundary of the study area - this is loaded in "helpers.R"
# MaxEdge = about 1/5 of the maximum spatial correlation we can expect.
# track = the transect track (containing survey effort in km2)
# obs.model = the observations as a shape, generated in the load_month_data function
MESH <- create.mesh(mask = mask,
                    MaxEdge = 4,
                    track=data.loaded$track,
                    obs.model = data.loaded$obs.model)

# Prepare the model -------------------------------------------------------

# This creates the matern surface which is used to generate the model. It also adds some noise
# to the data for species where data distributions are very non-randomly distributed
# E.G., common scoter has only a few points where they are located, but in large numbers

# obj.list = loaded month data (this is the output from the load_month_data function)
# If you get an error when running this, try adding jitter to the data by going jitter.points=TRUE
# MESH = the fishnet/mesh that was generated from create.mesh()

model.data <- model.prep(obj.list = data.loaded,
                         jitter.points = TRUE, 
                         MESH = MESH,
                         priorrange = c(15,0.5),
                         priorsigma = c(1,0.5))

# Run the spde model ------------------------------------------------------

# This will run the model and may take a few minutes to complete.
# This is parameterized with only the spatial structure of the underlying aerial survey data
# modprep.list = the output from the model.prep function
lgcp.model <- lgcp.spde(modprep.list = model.data)


# Create model predictions in space ---------------------------------------

# This generates the predictions
# MESH = the mesh/fishnet generated from the create.mesh() function
# lgcp.fit = the output from lcpe.spde()
predictions <- model.predictions(MESH = MESH,lgcp.fit = lgcp.model)

# Save the predictions to the output folder
save.predictions(predicted=predictions$predicted,month = month,spp=spp)

# Plot the predicted distribution -----------------------------------------

# This generates the output plot and can be saved using ggsave()
# predicted = output from model.predictions()function
# track = the survey track data from the load_month_data()function
# plot.title = the name of the plot (could be the name of the species)
# plot.boundary = T or F - to plot the outline of the SPA
# plot.coastline = T or F - to plot the coastline
# round.scale = this is for rounding the values in the density scale

plot.abundance(predicted=predictions$predicted,track=data.loaded$track,plot.title = "test",
                    plot.boundary = T,plot.coastline = F,round.scale = 0)



# Plot population estimate ------------------------------------------------

# This plots the predicted population size. Red line is the mean, blue lines are 
# 2.5% and 97.5% limits.
# pred.list = output from model.predictions

plot.pop.priors(pred.list = predictions)


# Get population estimate numerically -------------------------------------

# To get the estimated population go: 
predictions$marginals


# Wrangling the IWW data --------------------------------------------------

## This command loads the VP data into the environment
load_VP_data()

## Subsets the IWW data by the survey date above
locs.iww <- subset_IWW_date(survey.date = survey.date)

## Creates the buffer around the VPs for comparing against HD data
## This is set to 2 (km) at the minute - see buffer.size=2
## But could be changed to something else
## We are recommending this stay at 2km as it's the maximum distance at which 
## birds can be reliably identified. 

locs.sea <- calculate.buffer(locs.iww = locs.iww,
                             mask = mask,
                             buffer.size=2)

## Merges IWW vantage points with the WeBS polygon for a naming convention. 

#Match up any differing species names in HiDef/IWW data
iww$taxonName <- plyr::revalue(iww$taxonName, c("Common eider"="Eider"))

merge.iww_to_WeBS_polys(iww = iww,
                        locs.sea = locs.sea,
                        polygs = polygs)

## displays a simple plot (optional)
plot.basic.view(locs.sea = locs.sea,
                mask = mask,
                locs.iww = locs.iww,
                polygs.iww = polygs.iww)

## Subset species and get mean counts for selected locations (i.e. for when adjacent areas are merged)
counts <- subset_species(iww.locs = iww.locs,spp = spp)

## This outputs a list of the vantage point polygons (locs.dens),
## the WeBS polygons (polygs.dens), and a raster representation
## of the buffers. 

outputs <- create.polygons(locs.sea,
                           iww.aggr,
                           counts,
                           save.output = F,
                           spp = spp,
                           survey.date = survey.date)

### If you want to have a quick look at the outputs, try these commands
spplot(outputs$locs.dens, zcol= "dens")
spplot(outputs$polygs.dens, zcol= "dens", colorkey=FALSE) 

plot(mask)   # polygon mask Moray Firth
plot(locs.iww, add=TRUE, col="red", pch=20) # VPs in the IWW data
plot(polygs.iww, add=TRUE, col="transparent", border="orange")


# Comparing IWW data to HiDef flights -------------------------------------

### This will load up the HiDef modelled data that was saved using the 
### save.predictions function
load.rasters(spp,survey.date,month)


# Comparing hidef raster to vantage point raster --------------------------

## Vantage point polygons were converted to rasters in the create.polygons function
## we can compare those rasters at face value (grid cell to grid cell)
valuesout <- compare.rasters(hd.raster,outputs$buffer_ras,log.scale=T)



# Comparing Hidef raster to vantage point polygons ------------------------

## Using this function we can look at the mean density of birds/hr in the 
## IWW polygons (2km buffers) and compare against the mean of the HiDef 
## modelled output. outputs$locs.dens is the polygon shapefile for IWW

psout <- compare.polygons(pgons = outputs$locs.dens,
                          hd.raster,
                          log.scale = T)


# Comparing a subset of polygons to the raster ----------------------------

## In some cases you might want to explore a subset of the IWW polygons. 
## First, look at the list of names and then put the ones you want to examine
## into the list. 

## NOTE the log.scale=T function. This is used to make the comparison
## on a log scale, which will give a better idea of if data are correlated.

print(outputs$locs.dens@data$sector)


site.list <- c("Embo","Culbin Bar","Burghead to Hopeman",
               "Lossie Mouth to Spey Mouth","Dornoch Firth",
               "Tarbat Ness to Portmahomack","Wilkhaven to Rockfield")

## Use the extract.polygons function to get the desired polygons
pgonsout <- extract.polygons(outputs$locs.dens,site.list)

## Compare the polygons to the HD raster using the compare polygons function
subsetcompare <- compare.polygons(pgons = pgonsout,
                                  hd.raster = hd.raster,
                                  log.scale = T)



# Compare IWW to WeBS -----------------------------------------------------
## Compares the species of choice to the IWW data. Use the species list
## i.e. Species,  to view the spellings

WeBS.To.IWW <- webs.to.iww(spp,log.scale=F)
#Plot by using: 
WeBS.To.IWW




# Compare RSPB to WeBS ----------------------------------------------------


RSPBfiles <- list.files("Data/RSPB_data_filtered/",full.names = TRUE)
spp <- "eider"

## Nairn_Culbin_Bars, Outer_Dornoch_Firth, or Inverness_Beauly_Firth
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










