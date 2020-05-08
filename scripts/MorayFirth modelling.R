
library(readxl)
library(rgdal)
library(rgeos)
library(INLA)
library(mapproj)
library(viridis)
library(raster)

#we need the development version for a couple of functions
#devtools::install_github("fbachl/inlabru")
library(inlabru)

setwd("C:/Workspace/Nextcloud/Bioconsult/MorayFirth")
geogr<-CRS("+init=epsg:4326")
main.crs = "+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"



# 1. Load data ------------------------------------------------------------


    # accesory layers
      area <- readOGR(dsn = "Shapefile", layer = "Moray_Firth_Area_WGS84")
      area<-spTransform(area, main.crs)

    # EFFORT
      track01 <- readOGR(dsn = "2020 - Month 01 - Survey 01/Output/Zone87_M01_S01_20_Output", layer = "Zone87_M01_S01_20_Output-Day1-Transects")
      plot(track01)
      head(track01)
      track01<-spTransform(track01, main.crs)

    # SIGHTINGS
      table <- "2020 - Month 01 - Survey 01/Observations/Zone87_M01_S01_20_Observations_MSS.xlsx"

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
      plot(area, col="grey80")
      plot(track01, add=TRUE, pch=20)
      plot(obs01.model, col="red", add=TRUE, pch=19)
      

# 2. Model preparation ----------------------------------------------------
      
      # Mesh
      MaxEdge    <- 4
      mesh1 <- inla.mesh.2d(boundary=area, max.edge=c(1,5)*MaxEdge, offset=c(1,5))

      
      ggplot() + gg(track01) + gg(mesh1) + gg(obs01.model, col="red") + gg(area) + coord_fixed(ratio = 1)

      mesh1$n # total number of mesh nodes

      # Create Matern 
      matern2D = inla.spde2.pcmatern(mesh1,  prior.sigma = c(1, 0.5), prior.range = c(15, 0.5))  # 10, 0.1
      # INLA options
      c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)     
      
      # project effort to the 2D mesh vertices (for computational convenience)
      
      names(track01) <- c("length", "weight", "transect" )
      
      columns = "weight"
        res <- INLA::inla.fmesher.smorg(mesh1$loc, mesh1$graph$tv, points2mesh = coordinates(track01))
        tri <- res$p2m.t
        
        data <- list()
        for (k in 1:length(columns)) {
          cn <- columns[k]
          nw <- track01@data[, columns] * res$p2m.b
          w.by <- by(as.vector(nw), as.vector(mesh1$graph$tv[tri, ]), sum, simplify = TRUE)
          data[[cn]] <- as.vector(w.by)
        }
        
        data <- data.frame(data)
        coords <- mesh1$loc[as.numeric(names(w.by)), c(1, 2)]
        data$vertex <- as.numeric(names(w.by))
        
        fit.ips <- SpatialPointsDataFrame(coords,
                                      proj4string = CRS(proj4string(track01)),
                                      data = data,
                                      match.ID = FALSE
        )
        coordnames(fit.ips) <- coordnames(track01)
      
      

# 3. LGCP-SPDE model ------------------------------------------------------

      
      # model components
      cmp = coordinates ~ spatial_spde(map = coordinates, model = matern2D) + Intercept   # just spatial random effect
      
      #fit the point process model
      lgcp.fit= lgcp(components = cmp,               # model formula
                 data = obs01.model,                 # sighting locations
                 ips = fit.ips,                     # integrated effort
                 options=list(control.compute=c.c))  # options
      
      
      INLA:::summary.inla(lgcp.fit)
      lgcp.fit$dic$dic
      lgcp.fit$waic$waic

      

# 4. prediction -----------------------------------------------------------

      
      # set a PREDICTION AREA
      
      
      nx= 150
      ny= 150
      x <- seq(min(mesh1$loc[, 1]), max(mesh1$loc[, 1]), length = nx)
      y <- seq(min(mesh1$loc[, 2]), max(mesh1$loc[, 2]), length = ny)
      
      lattice <- INLA::inla.mesh.lattice(x = x, y = y)
      pixels <- data.frame(x = lattice$loc[, 1], y = lattice$loc[, 2])
      coordinates(pixels) <- c("x", "y")
      pxl <- SpatialPixels(pixels)
      pxl <- pixels[is.inside(mesh1, coordinates(pixels))]
      crs(pxl) <- main.crs
      pxl.sub <-  pxl[area,] 
      
      plot(pxl.sub)

      
      # full prediction
      sp.intensity = predict(lgcp.fit, pxl.sub, ~ exp(spatial_spde + Intercept))
      
      predicted <- as.data.frame(sp.intensity, xy=TRUE)
      
      #plot
      ggplot(predicted, aes(x=x, y=y, fill=mean)) + geom_raster() +
        scale_fill_gradient(low = "blue", high = "red",name="Density [Ind/sq. km]")+
       # scale_fill_viridis(name = "Intensity") +
        theme_bw(base_size=16) +
        coord_equal()
      
      
      # estimated posterior number of individuals + IC
      
      min.range <-100 # from (min expected range of N)
      max.range <- 3000  # to (max expected range of N)
      abund <- predict(lgcp.fit, ipoints(area, mesh1), 
                      ~ data.frame(N = min.range:max.range,
                                   dpois(min.range:max.range,
                                         lambda = sum(weight * exp(spatial_spde + Intercept)))))
      
      marginals <- inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x=abund$N, y = abund$mean)); marginals
      marg.mean <- inla.emarginal(identity, marginal = list(x=abund$N, y = abund$mean)) ; marg.mean
      
      
      #plot posterior mean and IC
      ggplot(data = abund, aes(x = N, y = mean)) +
        geom_area(fill="orange") +
        geom_vline(xintercept=marginals[1], linetype="dashed", color = "blue") +
        geom_vline(xintercept=marginals[2], color = "red") +
        geom_vline(xintercept=marginals[3], linetype="dashed", color = "blue") +
        ggtitle("Estimated number of Individuals") +
        theme_classic(base_size = 16) +
        xlim(600, 1500)
      
      
      
      # check range
      int.plot <- plot(lgcp.fit, "Intercept")
      spde.range <- spde.posterior(lgcp.fit, "spatial_spde", what = "range")
      plot(spde.range)

      