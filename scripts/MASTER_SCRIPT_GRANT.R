source("scripts/helpers.R")
source("scripts/species_list.R")
spp <- "Common scoter"
month <- "01"
data.loaded <- load_month_data(month,spp)

MESH <- create.mesh(mask = mask,
                    MaxEdge = 4,
                    track=data.loaded$track,
                    obs.model = data.loaded$obs.model)

model.data <- model.prep(obj.list = data.loaded,
                         jitter.points = TRUE, 
                         MESH = MESH,
                         priorrange = c(15,0.5),
                         priorsigma = c(1,0.5))

lgcp.model <- lgcp.spde(modprep.list = model.data)

predictions <- model.predictions(MESH = MESH,lgcp.fit = lgcp.model)

plot.abundance(predicted=predictions$predicted,track=data.loaded$track,plot.title = "test",
               plot.boundary = F,plot.coastline = F,round.scale = 0)

plot.pop.priors(pred.list = predictions)

predictions$marginals

save.predictions(predicted=predictions$predicted,
                 spp = spp,
                 month=month)

# Wrangle the IWW data ----------------------------------------------------


SurveyDate <- "2020-01-19"
spp <- "Common scoter"

load_VP_data()
locs.iww <- subset_IWW_date(survey.date = SurveyDate)

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
                           survey.date = SurveyDate)




