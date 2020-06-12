
source("scripts/helpers.R")

spp <- "Common scoter"  # Species subset ("Common eider", "Common scoter", "Red-throated diver")
survey.date <- "2020-01-19"  # 1st aerial survey: "2020-01-19", 2nd aerial survey: "2020-03-08"
month <- "01"

load_VP_data()
locs.iww <- subset_IWW_date(survey.date = survey.date)

locs.sea <- calculate.buffer(locs.iww = locs.iww,
                             mask = mask,
                             buffer.size=2)


merge.iww_to_WeBS_polys(iww = iww,
                        locs.sea = locs.sea,
                        polygs = polygs)
## displays the sites and WeBS sectors that have been extracted
plot.basic.view(locs.sea = locs.sea,
                mask = mask,
                locs.iww = locs.iww,
                polygs.iww = polygs.iww)

## Subset species and get mean counts for selected locations (i.e. for when adjacent areas are merged)
counts <- subset_species(iww.locs = iww.locs,spp = spp)

outputs <- create.polygons(locs.sea,
                           iww.aggr,
                           counts,
                           save.output = T,
                           spp = spp,
                           survey.date = survey.date)

### If you want to have a quick look at the outputs, try these commands
spplot(outputs$locs.dens, zcol= "dens")
spplot(outputs$polygs.dens, zcol= "dens", colorkey=FALSE) 
plot(mask)   # polygon mask Moray Firth
plot(locs.iww, add=TRUE, col="red", pch=20) # VPs in the IWW data
plot(polygs.iww, add=TRUE, col="transparent", border="orange")



# Comparing IWW to aerial survey data -------------------------------------


load.rasters(spp,survey.date,month)

valuesout <- compare.rasters(hd.raster,outputs$buffer_ras,log.scale=T)


site.list <- c("Embo","Culbin Bar","Burghead to Hopeman",
               "Lossie Mouth to Spey Mouth","Dornoch Firth",
               "Tarbat Ness to Portmahomack","Wilkhaven to Rockfield")
pgonsout <- extract.polygons(outputs$locs.dens,site.list)

psout <- compare.polygons(outputs$locs.dens,hd.raster,log.scale = F)


subsetcompare <- compare.polygons(pgons = pgonsout,
                 hd.raster = hd.raster,
                 log.scale = F)


psout <- compare.polygons(polygs.dens=pgonsout,hd.raster,log.scale = F)


