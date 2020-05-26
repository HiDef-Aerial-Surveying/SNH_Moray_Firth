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
                         priorrange = c(15,0.95),
                         priorsigma = c(1,0.5))

lgcp.model <- lgcp.spde(modprep.list = model.data)

predictions <- model.predictions(MESH = MESH,lgcp.fit = lgcp.model)

plot.abundance(predicted=predictions$predicted,track=data.loaded$track,plot.title = "test",
               plot.boundary = F,plot.coastline = F,round.scale = 0)

plot.pop.priors(pred.list = predictions)

predictions$marginals
