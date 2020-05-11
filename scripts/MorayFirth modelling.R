
library(readxl)
library(rgdal)
library(rgeos)
library(INLA)
library(mapproj)
library(raster)
library(inlabru)

setwd("C:/Workspace/Nextcloud/Github/HiDef/SNH_Moray_Firth")
geogr<-CRS("+init=epsg:4326")
# we are going to use the UTM30N as main reference system 
# (note that we set crs units to km...this will ensure that we are always taling about km...i.e. predicted density will be Individuals/sq. km, range will be in km, etc)
main.crs = "+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"


# 1. Load data ------------------------------------------------------------


    # accesory layers
      mask <- readOGR(dsn = "data/Shapefile", layer = "Moray_Firth_Area_WGS84")
      mask<-spTransform(mask, main.crs)

    # EFFORT
      track01 <- readOGR(dsn = "data/2020 - Month 01 - Survey 01/Output/Zone87_M01_S01_20_Output", layer = "Zone87_M01_S01_20_Output-Day1-Transects")
      plot(track01)
      head(track01)
      track01<-spTransform(track01, main.crs)

    # SIGHTINGS
      table <- "data/2020 - Month 01 - Survey 01/Observations/Zone87_M01_S01_20_Observations_MSS.xlsx"

      obs01 <- read_excel(table, sheet = 1)
      unique(obs01$Species)
      obs01 <- subset(obs01, Species != "N/A")
 
      coordinates(obs01) <- c("Longitude", "Latitude")
      crs(obs01) <- geogr

    # Select species
      unique(obs01$Species)
      table(obs01$Species)
      
      spp <- "Common scoter"  #"Herring gull" ; "Razorbill" ;"Common scoter"
      
      obs01.model <- subset(obs01, Species == spp)
      obs01.model<-spTransform(obs01.model, main.crs)
      
      #plot
      plot(mask, col="grey80")
      plot(track01, add=TRUE, pch=20)
      plot(obs01.model, col="red", add=TRUE, pch=19)
      

# 2. Model preparation ----------------------------------------------------
      
      # CREATE MESH
      MaxEdge    <- 4                     # MaxEdge is about 1/5 of the maximum spatial correaltion we can expect.
                                          # maximum spatial correlation would be the distance at which the variogram would reach the sill
      
      # The mesh is our prediction "grid", but instead of a regular mesh, we use an irregular mesh with higher node density in heavier sampled areas
      # In adittion, we want to use a double mesh that extends far from the prediction area to avoid boundary effects
      mesh1 <- inla.mesh.2d(boundary=mask, max.edge=c(1,5)*MaxEdge, offset=c(1,5)) 
      # the first value (1) sets the edge size of the mesh in the inner zone, and the second value (5) sets the edge size of the outer mesh

      #we plot the mesh and the observations together, bear in mind that inlabru can plot objects directly using the gg() function
      ggplot() + gg(track01) + gg(mesh1) + gg(obs01.model, col="red") + gg(mask) + coord_fixed(ratio = 1)

      mesh1$n # this is the total number of mesh nodes (more nodes -> more processing time)

      # CREATE MATERN MODEL
      # based in the mesh and priors. prior sigma is related with the smoothness of the matern surface and prior.range with the maximum expected reach of the spatial correlation (in km)
      # the 0.5 is our confidence in the sigma and range provided as priors, we do not want to make any assumption, so we give a very low confidence: 0.5
      matern2D = inla.spde2.pcmatern(mesh1,  prior.sigma = c(1, 0.5), prior.range = c(15, 0.5))  # non-informative priors
      
      # INLA options (just asking for WAIC and DIC for model validation)
      c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)     
      
      
      # project effort to the 2D mesh vertices 
      # This is done for computational convenience. basically que aggregate the effort in the mesh nodes
      #not need to understand the code, I will wrap it into a function.
      # Important!, the area observed at each point must be called "weight", so we have to change the column names before
      
      names(track01) <- c("length", "weight", "transect" )
      
      columns = "weight"
 
      source("./scripts/tools.R")
      
      fit.ips <- integrate_effort(mesh=mesh1, effort=track01)
        
        
        #IMPORTANT! for common scoter many sightings have the same coordinates and INLA will crash
        # one solution (better) is to assum the number of sightings at each location as a mark and run a marked point process model
        # However this will complicate the model, so we will go for the "dirty" solution:
        # WE will add some random noise to the sights, just few meters. 
        # We assume that moving points few meters away so they do not overlap do not have any impact at the scale of the model
      
            # check number of unique positions and compare with total number of sights
                obs.df <- data.frame(coordinates(obs01.model))
                obs.df$count <- 1
                
                obs.cs <- aggregate(.~Longitude+Latitude, obs.df, sum)
                nrow(obs.cs)  # number of unique pair of coordinates
                nrow(obs.df)  # total number of sightings
            
            
            
            # for Common scoter it is 55 unique pair of coordinates and 1694 total sightings...we need to jitter sightings
            
            #Jittering sightings (only if needed!)
            
                obs01.model <- jitter_sightings(obs01.model)
      
      
      
      plot(mask, col="grey80")
      plot(track01, add=TRUE, pch=20)
      plot(obs01.model, col="red", add=TRUE, pch=19)

        
# 3. LGCP-SPDE model ------------------------------------------------------

      
      # model components
          # We define the model here, we model the point locations (coordinates) using a random gaussin field model we defined before (matern2D) and an intercept
      cmp = coordinates ~ spatial_spde(map = coordinates, model = matern2D) + Intercept   # just spatial random effect, no covariates
      
      #fit the point process model
      lgcp.fit= lgcp(components = cmp,               # model formula
                 data = obs01.model,                 # sighting locations
                 ips = fit.ips,                     # integrated effort 
                 options=list(control.compute=c.c))  # INLA options (defined before)
      
      
      INLA:::summary.inla(lgcp.fit)
      lgcp.fit$dic$dic  # for model validation and compare models
      lgcp.fit$waic$waic

      

# 4. prediction -----------------------------------------------------------

      
      # set a PREDICTION AREA
      
     pxl.sub <- prediction_fun(nx=150, ny=150, mesh= mesh1) # nx and ny set the spatial resoltion output
      
      plot(pxl.sub)

      
      # full prediction
      sp.intensity = predict(lgcp.fit, pxl.sub, ~ exp(spatial_spde + Intercept)) # predict model lgcp.fir on pxl.sub using the formula ~ exp(spatial_spde + Intercept
      
      predicted <- as.data.frame(sp.intensity, xy=TRUE)
      
      #plot
      ggplot(predicted, aes(x=x, y=y, fill= mean)) + geom_raster() +
         scale_fill_gradient( trans = 'sqrt', low = "blue", high = "red",name="Density [Ind/sq. km]" ) +
       # scale_fill_gradient(low = "blue", high = "red",name="Density [Ind/sq. km]")+
       # scale_fill_viridis(name = "Intensity") +
        theme_bw(base_size=16) +
        coord_equal()
      
      
      # estimated posterior number of individuals + IC
      
      min.range <-0 # from (min expected range of N)
      max.range <- 50000  # to (max expected range of N)
      
      # for a full bayesian estimation of the mean and Interval confidence, we set min.range and max.range as the limits of the prediction interval
      abund <- predict(lgcp.fit, ipoints(mask, mesh1), 
                      ~ data.frame(N = min.range:max.range,
                                   dpois(min.range:max.range,
                                         lambda = sum(weight * exp(spatial_spde + Intercept)))))
      
      marginals <- inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x=abund$N, y = abund$mean)); marginals # get Q0.025, median and Q0.975
      marg.mean <- inla.emarginal(identity, marginal = list(x=abund$N, y = abund$mean)) ; marg.mean  #get mean estimated population
      
      lower.lim <- marginals[1] - 1000  # adjust later
      upper.lim <- marginals[3] + 1000  # adjust later
      
      #plot posterior mean and IC
      ggplot(data = abund, aes(x = N, y = mean)) +
        geom_area(fill="orange") +
        geom_vline(xintercept=marginals[1], linetype="dashed", color = "blue") +
        geom_vline(xintercept=marginals[2], color = "red") +
        geom_vline(xintercept=marginals[3], linetype="dashed", color = "blue") +
        ggtitle("Estimated number of Individuals") +
        theme_classic(base_size = 16) +
        xlim(lower.lim, upper.lim)
      
      
      
      # check range
          # this is the posterior calculated by the model for the range (that we very roughly estimated for the non-informative prior)
      int.plot <- plot(lgcp.fit, "Intercept")
      spde.range <- spde.posterior(lgcp.fit, "spatial_spde", what = "range")
      plot(spde.range)

      