###############
## Run this to compare all the contemporary data
#################
source("scripts/helpers.R")
source("scripts/species_list.R")
#View(Species)  ## View this if you need to
spp <- "common scoter"
month <- "01"  #For this analysis, only use January
## Set the title of the output plot here


# Wrangling the IWW data --------------------------------------------------

## This command loads the VP data into the environment
load_VP_data()

## Subsets the IWW data by the survey date - set to the january survey
## because this is the only 
locs.iww <- subset_IWW_date(survey.date = "2020-01-19" )

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


### This will compare the data in the three RSPB sub-regions
iww.pgons <- outputs$locs.dens
compare.all.data(iww.pgons,spp)


