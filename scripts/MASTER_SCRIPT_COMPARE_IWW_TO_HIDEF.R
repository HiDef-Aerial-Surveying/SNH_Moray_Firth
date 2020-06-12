########################################
### This script compares IWW to HiDef models
#######################################

source("scripts/helpers.R")
source("scripts/species_list.R")
#View(Species)  ## View this if you need to
spp <- "common scoter"
#month <- "01"
month <- "03"
#survey.date <- "2020-01-19"  
survey.date <- "2020-03-08"  

## Set the title of the output plot here
plot.title = "Common scoter IWW counts (March 2020)"

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

### This creates a nice looking plot of the buffers
plot.iww.pgons(pgons=outputs$locs.dens,track=data.loaded$track,
               plot.title=plot.title)


