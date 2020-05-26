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



# Plot the predicted distribution -----------------------------------------

# This generates the output plot and can be saved using ggsave()
# predicted = output from model.predictions()function
# track = the survey track data from the load_month_data()function
# plot.title = the name of the plot (could be the name of the species)
# plot.boundary = T or F - to plot the outline of the SPA
# plot.coastline = T or F - to plot the coastline
# round.scale = this is for rounding the values in the density scale

plot.abundance(predicted=predictions,track=data.loaded$track,plot.title = "test",
                    plot.boundary = T,plot.coastline = F,round.scale = 0)



# Plot population estimate ------------------------------------------------

# This plots the predicted population size. Red line is the mean, blue lines are 
# 2.5% and 97.5% limits.
# pred.list = output from model.predictions

plot.pop.priors(pred.list = predictions)


# Get population estimate numerically -------------------------------------

# To get the estimated population go: 
predictions$marginals



