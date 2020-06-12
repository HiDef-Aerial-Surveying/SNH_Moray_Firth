################################
## Spatial modeling code
## This is the master script for running the digital aerial survey models
################################
source("scripts/helpers.R")
source("scripts/species_list.R")
#View(Species)  ## View this if you need to
spp <- "common scoter"
#month <- "01"
month <- "03"
#survey.date <- "2020-01-19"  
survey.date <- "2020-03-08"  

## Set the title of the abundance plot here
plot.title = "Density of Common Scoter (March 2020)"

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

#### You can read a raster this way and plot it if 
#### you would like to explore it a little bit
#### Change the name of the file if need be
#prediction.raster <- readRDS("outputs/Common scoter_month01.rds")

# Plot the predicted distribution -----------------------------------------
plot.abundance(predicted=predictions$predicted,track=data.loaded$track,
               plot.title = plot.title,
               plot.boundary = T,plot.coastline = T,round.scale = 0)

# Plot population estimate ------------------------------------------------
plot.pop.priors(pred.list = predictions)

# Get population estimate numerically -------------------------------------
predictions$marginals
