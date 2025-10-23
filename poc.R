library(terra)
library(sf)
library(tidyverse)
library(nabor)
library(dbscan)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)

cities <- sf::read_sf("data/cities.gpkg") |> sf::st_transform(3035)


biorefineries <- sf::read_sf("data/Biorefineries.shp")|> sf::st_transform(3035)
companies <- sf::read_sf("data/Companies.shp")|> sf::st_transform(3035)
ports.es <- sf::read_sf("data/Ports_ES.shp")|> sf::st_transform(3035)
ports.it <- sf::read_sf("data/Ports_IT.shp")|> sf::st_transform(3035)
ports.gr <- sf::read_sf("data/Ports_GR.shp")|> sf::st_transform(3035)
ports<- do.call(rbind, list(ports.es, ports.it, ports.gr))

maxDist <- 100000
biorefDist <- nabor::knn(sf::st_coordinates(biorefineries), sf::st_coordinates(cities), k=1 )  # fast radius neighbors
companiesDist <- nabor::knn(sf::st_coordinates(companies), sf::st_coordinates(cities), k=1)  # fast radius neighbors
portsDist <- nabor::knn(sf::st_coordinates(ports), sf::st_coordinates(cities), k=1)  # fast radius neighbors

cities$score.compDists <-  companiesDist$nn.dists[,1]
cities$score.biorefDists <- biorefDist$nn.dists[,1]
cities$score.portsDists <-  portsDist$nn.dists[,1]


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
    sf::st_crs(r.df.sf) <- terra::crs(r)
    pts[[countryn]] <- r.df.sf
    save(pts, file="pts.rda")
  }

  city <- cities |> filter(country==countryn) |> select(3:7)
  city$cityId <- 1:nrow(city)

  vor <- sf::st_voronoi(st_union(city))
  vor_sf <- st_collection_extract(vor)
  vor_sf <- st_sf(geometry = vor_sf)

  vor_sf <- st_join(vor_sf, city, join = st_intersects)

  pts_with_poly <- st_join(pts[[countryn]], vor_sf, join = st_intersects)

  suit <- lapply( city$cityId, function(id){
    ct <- city |> filter(cityId==id)
    pt <- pts_with_poly |> filter(cityId==id)
    w <- sf::st_distance(pt, ct)
    w <- as.numeric(w)
    w[ w == 0] <- 1e6
    w <- 1 / (w^1)
    sum(pt$mean * w)
  })

  city$score.suitability <- unlist(suit)
  sf::write_sf(city, sprintf("%sSuitableCities2.gpkg", countryn))

  n <- length(city$score.suitability)
  weights <- city$score.suitability
  costsO <- (city$score.compDists + city$score.biorefDists + city$score.portsDists)
  #normalizzo
  costsO <- costsO / max(costsO)
  # hist(costsO)
  weights <- weights / max(weights)
  #hist(weights)
  output <- list(count=c(), totSuitability=c(), totCost=c() )
  for(cost in 1:20/2){

    costs <- costsO * cost

    model <- MIPModel() %>%
      add_variable(x[i], i = 1:n, type = "binary") %>%
      set_objective(sum_expr(weights[i]*x[i] - costs[i]*x[i], i = 1:n), "max") %>%
      # optionally: add_constraint(sum_expr(costs[i]*x[i], i = 1:n) <= B)
      solve_model(with_ROI(solver = "glpk"))

    solution <- get_solution(model, x[i]) #%>% filter(value > 0.5)
    city$selectedCities  <- solution$value
    plot(city[,"selectedCities"] )
    output$count <- c(output$count, sum(city$selectedCities) )
    output$totSuitability <- c(output$totSuitability, sum(weights*city$selectedCities) )
    output$totCost <- c(output$totCost, sum(costs) )
  }

  sf::write_sf(city, sprintf("%sSuitableCities.gpkg", countryn))


}




km <- kmeans(sf::st_coordinates(city), centers = 5)
city$cluster <- as.factor(km$cluster)
for(k in levels(city$cluster)){
  candidateCities <- city |> filter(cluster==k)





}



