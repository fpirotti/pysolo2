library(terra)
library(sf)
library(tidyverse)
library(nabor)
library(dbscan)

cities <- sf::read_sf("data/cities.gpkg") |> sf::st_transform(3035)


biorefineries <- sf::read_sf("data/Biorefineries.shp")|> sf::st_transform(3035)
companies <- sf::read_sf("data/Companies.shp")|> sf::st_transform(3035)

maxDist <- 100000
biorefDist <- nabor::knn(sf::st_coordinates(biorefineries), sf::st_coordinates(cities), k=1 )  # fast radius neighbors
companiesDist <- nabor::knn(sf::st_coordinates(companies), sf::st_coordinates(cities), k=1)  # fast radius neighbors

cities$score.compDists <-  (1 / companiesDist$nn.dists^1)
cities$score.biorefDists <-  (1 / biorefDist$nn.dists^1)


pts <- list()
load("pts.rda")
countries <- unique(cities$country)

for(countryn in countries){

  if(countryn!="Spain") next

  if(is.null(pts[[countryn]])){
    r <- terra::rast(sprintf("data/AgriSuitabilitPysolo3035_%s.tif", countryn))
    r[r==0] <- NA
    r.df <- terra::as.data.frame(r, xy=T)
    r.df.sf <- sf::st_as_sf(r.df,coords = c("x", "y"))
    pts[[countryn]] <- r.df.sf
    save(pts, file="pts.rda")
  }
  city <- cities |> filter(country==countryn)

  city$cluster <- as.factor(km$cluster)
  fields <- sf::st_coordinates(pts[[countryn]])

  res <- frNN(fields, eps = maxDist, query=sf::st_coordinates(city), sort = FALSE)  # fast radius neighbors

  suit <- lapply( 1:nrow(city), function(id){
    w <- res$dist[[id]]
    w[ w == 0] <- 1e6
    w[] <- 1 / (w^1)
    sum(pts[[countryn]]$mean[ res$id[[id]] ] * w)
  })
  city$score.suitability <- unlist(suit)
  sf::write_sf(city, sprintf("%sSuitableCities.gpkg", countryn))


  km <- kmeans(sf::st_coordinates(city), centers = 5)
  for(k in levels(city$cluster)){
    candidateCities <- city |> filter(cluster==k)





  }



  plot(city["cluster"], pch = 19)

}

#plot(r.df[,c("x","y")], col=r.df["cluster"], pch = 19)

plot(r)


