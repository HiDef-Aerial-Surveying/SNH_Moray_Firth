
integrate_effort = function(mesh,effort){
  
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
  coords <- mesh1$loc[as.numeric(names(w.by)), c(1, 2)]
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
  crs(obs.cs) <- main.crs
  return(obs.cs)
  
}


prediction_fun <- function(nx,ny, mesh)  {

  x <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = nx)
  y <- seq(min(mesh1$loc[, 2]), max(mesh1$loc[, 2]), length = ny)
  
  lattice <- INLA::inla.mesh.lattice(x = x, y = y)
  pixels <- data.frame(x = lattice$loc[, 1], y = lattice$loc[, 2])
  coordinates(pixels) <- c("x", "y")
  pxl <- SpatialPixels(pixels)
  pxl <- pixels[is.inside(mesh, coordinates(pixels))]
  crs(pxl) <- main.crs
  pxl.sub <-  pxl[mask,] 
  return(pxl.sub)
 
}

