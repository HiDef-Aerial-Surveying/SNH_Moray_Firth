## Modeling functions for Bayesian (INLA) spatial models
## Written for SNH Moray firth research project, May 2020
## Grant Humphries and Raul Vilela

# Load libraries ----------------------------------------------------------
source("scripts/load_libs.R")
LIBS <- c("tidyverse","foreach","rgeos","sp","maptools","rgdal","raster","Hmisc","ggthemes",
          "Cairo","doSNOW","readxl","INLA","mapproj","inlabru","ggnewscale")
Load_Libs(LIBS)

# Spatial Projections -------------------------------------------------------------
WGS84<-CRS("+init=epsg:4326")
## We use the longer form to define the CRS (instead of just +init) so that
## we are speaking in KM instead of M - this is to make sure the MESH function doesn't
## crash when using a value of 4
UTM30<-CRS("+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Species list ------------------------------------------------------------
source("scripts/species_list.R")


# Load study area ---------------------------------------------------------------
# This is the study area, which is used to define the boundary of the final maps

mask <- readOGR(dsn = "data/Shapefile", layer = "Moray_Firth_Area_WGS84")
mask <- spTransform(mask, UTM30)  # Transforms data to UTM30 spatial reference


# Functions ---------------------------------------------------------------

## Function to first load up the data filtered by month and species
load_month_data <- function(month,spp){
  ## Function to load up the transect and observation data
  ## month = month of the year as an integer (1 = 12)
  ## spp = common name of species to model
  spp <- tolower(spp)
  
  trckpath <- paste0("Data/Aerial_Surveys/2020 - Month ",
                     month,
                     " - Survey 01/Output/Zone87_M",
                     month,
                     "_S01_20_Output")
  
  trcklyr <- paste0("Zone87_M",
                    month,
                    "_S01_20_Output-Day1-Transects")
  
  obstabpath <- paste0("Data/Aerial_Surveys/2020 - Month ",
                       month,
                       " - Survey 01/Observations/Zone87_M",
                       month,
                       "_S01_20_Observations_MSS.xlsx")
  ## Loads the transect shapefile which contains effort (i.e. area surveyed at each point)
  track <- readOGR(dsn=trckpath,layer=trcklyr)
  track<-spTransform(track, UTM30) 
  ## Read the observation spreadsheet
  obs <- read_excel(obstabpath, sheet = 1)
  ## Removes capitalization so it's not a problem anymore
  obs$Species <- tolower(obs$Species)
  ## Filter the data to only get your chosen species
  obs <- subset(obs,obs$Species==spp)
  
  ## Turn obs into a vector (point) shapefile (i.e. a spatial points data frame)
  coordinates(obs) <- c("Longitude","Latitude")
  crs(obs) <- WGS84
  obs.model <- spTransform(obs,UTM30)
  print(paste0("Loaded ",spp," for month ",month))
  return(list(obs.model=obs.model,track=track))
}


integrate_effort = function(mesh,effort,columns){
  
  res <- INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, points2mesh = coordinates(effort))
  tri <- res$p2m.t
  
  data <- list()
  for (k in 1:length(columns)) {
    cn <- columns[k]
    nw <- effort@data[, columns] * res$p2m.b
    w.by <- by(as.vector(nw), as.vector(mesh$graph$tv[tri, ]), sum, simplify = TRUE)
    data[[cn]] <- as.vector(w.by)
  }
  
  data <- data.frame(data)
  coords <- mesh$loc[as.numeric(names(w.by)), c(1, 2)]
  data$vertex <- as.numeric(names(w.by))
  
  fit.ips <- SpatialPointsDataFrame(coords,
                                    proj4string = CRS(proj4string(effort)),
                                    data = data,
                                    match.ID = FALSE
  )
  coordnames(fit.ips) <- coordnames(effort) 
  
  return(fit.ips)
  
}



jitter_sightings <- function(sights){
  ## add random noise to points
  obs.cs <- as.data.frame(sights, xy=TRUE)
  obs.cs$Longitude.j <-jitter(obs.cs$Longitude, factor=1, amount = NULL)
  obs.cs$Latitude.j <-jitter(obs.cs$Latitude, factor=1, amount = NULL)
  coordinates(obs.cs) <- c("Longitude.j", "Latitude.j")
  crs(obs.cs) <- UTM30
  return(obs.cs)
  
}


prediction_fun <- function(nx,ny, mesh)  {
  
  x <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = nx)
  y <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = ny)
  
  lattice <- INLA::inla.mesh.lattice(x = x, y = y)
  pixels <- data.frame(x = lattice$loc[, 1], y = lattice$loc[, 2])
  coordinates(pixels) <- c("x", "y")
  pxl <- SpatialPixels(pixels)
  pxl <- pixels[is.inside(mesh, coordinates(pixels))]
  crs(pxl) <- UTM30
  pxl.sub <-  pxl[mask,] 
  return(pxl.sub)
  
}

create.mesh <- function(mask,MaxEdge=4,track,obs.model){
  # CREATE MESH
  # MaxEdge is about 1/5 of the maximum spatial correlation we can expect.
  # maximum spatial correlation would be the distance at which the variogram 
  # would reach the sill
  # The mesh is our prediction "grid", but instead of a regular mesh, 
  # we use an irregular mesh with higher node density in heavier sampled areas
  # In adittion, we want to use a double mesh that extends far from the 
  # prediction area to avoid boundary effects
  
  print('creating mesh...')
  mesh1 <- inla.mesh.2d(boundary=mask, max.edge=c(1,5)*MaxEdge, offset=c(1,5)) 
  # the first value (1) sets the edge size of the mesh in the inner zone, and the second value (5) sets the edge size of the outer mesh
  print('mesh created...')
  #we plot the mesh and the observations together, bear in mind that inlabru can plot objects directly using the gg() function
  M <- ggplot() + gg(track) + gg(mesh1) + gg(obs.model, col="red") + gg(mask) + coord_fixed(ratio = 1)
  return(list(mesh = mesh1, M = M))
}



model.prep <- function(obj.list,jitter.points=FALSE,MESH,
                       priorrange=c(15,0.5),priorsigma=c(1,0.5)){
  mesh1 <- MESH$mesh
  obs.model <- obj.list$obs.model
  track <- obj.list$track
  
  # CREATE MATERN MODEL
  # based in the mesh and priors. prior sigma is related with the smoothness of the 
  # matern surface and prior.range with the maximum expected reach of the spatial correlation (in km)
  # the 0.5 is our confidence in the sigma and range provided as priors, we do not want to make 
  # any assumption, so we give a very low confidence: 0.5
  
  print("creating matern surface")
  # non-informative priors
  matern2D = inla.spde2.pcmatern(mesh1,  prior.sigma = priorsigma, prior.range = priorrange)  
  
  # INLA options (just asking for WAIC and DIC for model validation)
  c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)     
  
  # project effort to the 2D mesh vertices 
  # This is done for computational convenience. basically que aggregate the effort in the mesh nodes
  #not need to understand the code, I will wrap it into a function.
  # Important!, the area observed at each point must be called "weight", so we have to change the column names before
  
  names(track) <- c("length", "weight", "transect" )
  columns = "weight"
  
  print("integrating effort into mesh")
  fit.ips <- integrate_effort(mesh=mesh1, effort=track,columns)
  
  # IMPORTANT! for common scoter many sightings have the same coordinates and INLA will crash
  # one solution (better) is to assum the number of sightings at each location as a mark and run a marked point process model
  # However this will complicate the model, so we will go for the "dirty" solution:
  # WE will add some random noise to the sights, just few meters. 
  # We assume that moving points few meters away so they do not overlap do not have any impact at the scale of the model
  
  if(jitter.points == TRUE){
    print("Adding Jitter to points...")
    # check number of unique positions and compare with total number of sights
    obs.df <- data.frame(coordinates(obs.model))
    obs.df$count <- 1
    
    obs.cs <- aggregate(.~Longitude+Latitude, obs.df, sum)
    nrow(obs.cs)  # number of unique pair of coordinates
    nrow(obs.df)  # total number of sightings
    
    
    #for Common scoter it is 55 unique pair of coordinates and 1694 total sightings...we need to jitter sightings
    #Jittering sightings (only if needed!)
    obs.model <- jitter_sightings(obs.model)
  }
  
  print("Ready")
  return(list(obs.model=obs.model,matern2D=matern2D,c.c=c.c,fit.ips=fit.ips))
  
}

lgcp.spde <- function(modprep.list){
  ## loaddata.list = output from the load_month_data function
  ## modprep.list = output from the model.prep function
  
  obs.model <- modprep.list$obs.model
  matern2D <- modprep.list$matern2D
  c.c <- modprep.list$c.c
  fit.ips <- modprep.list$fit.ips
  # model components
  # We define the model here, we model the point locations (coordinates) using a random gaussin field model we defined before (matern2D) and an intercept
  print("running spatial spde")
  cmp = coordinates ~ spatial_spde(map = coordinates, model = matern2D) + Intercept   # just spatial random effect, no covariates
  
  print("fitting the point process model")
  #fit the point process model
  lgcp.fit= lgcp(components = cmp,               # model formula
                 data = obs.model,                 # sighting locations
                 ips = fit.ips,                     # integrated effort 
                 options=list(control.compute=c.c))  # INLA options (defined before)
  
  
  print(INLA:::summary.inla(lgcp.fit))
  print(lgcp.fit$dic$dic)  # for model validation and compare models
  print(lgcp.fit$waic$waic)
  return(lgcp.fit)
}



model.predictions <- function(MESH,lgcp.fit){
  
  mesh1 <- MESH$mesh
  # set a PREDICTION AREA
  pxl.sub <- prediction_fun(nx=150, ny=150, mesh= mesh1) # nx and ny set the spatial resoltion output
  # full prediction
  # predict model lgcp.fir on pxl.sub using the formula ~ exp(spatial_spde + Intercept
  sp.intensity = predict(lgcp.fit, pxl.sub, ~ exp(spatial_spde + Intercept)) 
  predicted <- as.data.frame(sp.intensity, xy=TRUE)
  
  # estimated posterior number of individuals + IC
  
  min.range <-0 # from (min expected range of N)
  max.range <- 50000  # to (max expected range of N)
  
  # for a full bayesian estimation of the mean and Interval confidence, 
  # we set min.range and max.range as the limits of the prediction interval
  
  abund <- predict(lgcp.fit, ipoints(mask, mesh1), 
                   ~ data.frame(N = min.range:max.range,
                                dpois(min.range:max.range,
                                      lambda = sum(weight * exp(spatial_spde + Intercept)))))
  
  # get Q0.025, median and Q0.975
  marginals <- inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x=abund$N, y = abund$mean)); marginals 
  
  #get mean estimated population  
  marg.mean <- inla.emarginal(identity, marginal = list(x=abund$N, y = abund$mean)) ; marg.mean  
  
  return(list(predicted=predicted,abund=abund,marginals=marginals,marg.mean=marg.mean))
}


plot.abundance <- function(predicted,track,plot.title="",plot.boundary=F,
                           plot.coastline=F,round.scale=0){

  breaks<-c(seq(min(predicted$mean),
                      max(predicted$mean),
                      ceiling(max(predicted$mean)/3)),max(predicted$mean))
  labels <- round(breaks,round.scale)
  
  transdf <- data.frame(coordinates(track))
  
  xlims <- c(min(transdf$coords.x1)*0.99, max(transdf$coords.x1)*1.01)
  ylims <- c(min(transdf$coords.x2)*0.999, max(transdf$coords.x2)*1.001)
  
  G <- ggplot(predicted, aes(x=x, y=y)) + geom_raster(aes(fill= mean)) +
    scale_fill_gradient(trans = 'sqrt', low = "blue", high = "red",
                         breaks=breaks,
                         labels=labels,
                         name=expression(Density~(Ind/km^{"2"}))) +
    scale_x_continuous(breaks=seq(425,525,25),label=seq(4.25,5.25,.25))+
    scale_y_continuous(breaks=seq(6380,6440,20),label=seq(6.38,6.44,.02))+
    xlab(expression(X(meters~x~10^{"6"})))+
    ylab(expression(Y(meters~x~10^{"7"})))+
    ggtitle(plot.title)+
    theme_gdocs()+
    theme(
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black"),
      plot.background = element_blank()
    )+
    coord_equal(xlim = xlims, ylim = ylims)
    
  
  if(plot.coastline==TRUE){
    print("loading high resolution coastline.. this may take a moment..")
    load(file = 'D:/GIS_DATA/UK_coastline.rda')
    coastline_sp_utm <- spTransform(tt, UTM30)
    coastData_DF <- fortify(coastline_sp_utm)
    
    G <- G + 
    geom_polygon(aes(x = long, y = lat, group = group), data = coastData_DF, fill = "darkolivegreen3", 
                 col = "grey35", alpha = 0.7) 
  }
  if(plot.boundary==TRUE){
    boundary.shape <- readOGR(dsn="./Data/Shapefile/Moray_Firth_Area_WGS84.shp",
                              layer="Moray_Firth_Area_WGS84")
    boundaryUTM30 <- spTransform(boundary.shape,UTM30)
    boundary_DF <- fortify(boundaryUTM30)
    G <- G + 
      geom_polygon(aes(x=long,y=lat,group=group),data=boundary_DF,
                   col="black",alpha=0)
  }
  
  

  
  return(G)
}



plot.pop.priors <- function(pred.list){
  
  predicted<-pred.list$predicted
  abund<-pred.list$abund
  marginals<-pred.list$marginals
  marg.mean<-pred.list$marg.mean
  
  
  lower.lim <- marginals[1] - 1000  # adjust later
  upper.lim <- marginals[3] + 1000  # adjust later
  
  #plot posterior mean and IC
  G <- ggplot(data = abund, aes(x = N, y = mean)) +
    geom_area(fill="orange",color="black") +
    geom_vline(xintercept=marginals[1], linetype="dashed", color = "blue") +
    geom_vline(xintercept=marginals[2], color = "red") +
    geom_vline(xintercept=marginals[3], linetype="dashed", color = "blue") +
    ggtitle("Estimated number of Individuals") +
    xlab("Population size")+
    ylab("Mean")+
    scale_y_continuous(expand=c(0,0))+
    theme_gdocs()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black"),
      plot.background = element_blank()
    )+
    #theme_classic(base_size = 16) +
    xlim(lower.lim, upper.lim)
  
  # check range
  # this is the posterior calculated by the model for the range 
  # (that we very roughly estimated for the non-informative prior)
  #int.plot <- plot(lgcp.fit, "Intercept")
  #spde.range <- spde.posterior(lgcp.fit, "spatial_spde", what = "range")
  #plot(spde.range)
  return(G) 
}



