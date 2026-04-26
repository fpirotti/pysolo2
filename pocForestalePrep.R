library(terra)
library(sf)
library(tidyverse)
library(nabor)
library(dbscan)
library(ompr)
library(ompr.roi)
library(tidyterra)
library(s2)
library(magick)
library(parallel)

library(ROI.plugin.glpk)

cities <- sf::read_sf("data/cities3.gpkg") #|> select(3) #|> sf::st_transform(3035)

# Download city center  ----
## points for Spain, Italy, and Greece - uncomment if needed
# library(giscoR)
# countries <- c("ES", "IT", "EL") # Note: Greece is EL in Eurostat
# lau_data <- gisco_get_communes(country = c("ES", "IT", "EL"), spatialtype="LB")
# sf::write_sf(lau_data, "data/cities3.gpkg") #|> select(3) #|> sf::st_transform(3035)

# nuts <- sf::read_sf("/archivio/shared/geodati/vector/NUTS_2024_all_4326v2.gpkg")

biorefineries <- sf::read_sf("data/forest/Main Refineries_FINAL.shp") #|> sf::st_transform(3035)

ports <- sf::read_sf("data/forest/Official Ports_FINAL.shp")#|> sf::st_transform(3035)

maxDist <- c(50000, 100000, 150000)

calculate_weight <- function(distance, max_dist ) {
  weight <- 1 - (distance / max_dist)
  return(pmax(0, weight)) # pmax is the vectorized version of max()
}

biorefDist <- nabor::knn(sf::st_coordinates(sf::st_transform(biorefineries, 3035)),
                         sf::st_coordinates(sf::st_transform(cities, 3035)), k=1 )  # fast radius neighbors
# companiesDist <- nabor::knn(sf::st_coordinates(companies), sf::st_coordinates(cities), k=1)  # fast radius neighbors
portsDist <-  nabor::knn(sf::st_coordinates(sf::st_transform(ports, 3035)),
                            sf::st_coordinates(sf::st_transform(cities, 3035)), k=1 )

# cities$score.compDists <-  companiesDist$nn.dists[,1]
cities$score.biorefDistsT1 <- pmin(1,  biorefDist$nn.dists[,1] / maxDist[[1]])
cities$score.biorefDistsT2 <- pmin(1,  biorefDist$nn.dists[,1] / maxDist[[2]])
cities$score.biorefDistsT3 <- pmin(1,  biorefDist$nn.dists[,1] / maxDist[[3]])

cities$score.portsDistsT1 <-  pmin(1,  portsDist$nn.dists[,1] / maxDist[[1]])
cities$score.portsDistsT2 <-  pmin(1,  portsDist$nn.dists[,1] / maxDist[[2]])
cities$score.portsDistsT3 <-  pmin(1,  portsDist$nn.dists[,1] / maxDist[[3]])

cities$costTotT1 <- mean(c(cities$score.portsDistsT1, cities$score.biorefDistsT1))
cities$costTotT2 <- mean(c(cities$score.portsDistsT2, cities$score.biorefDistsT2))
cities$costTotT3 <- mean(c(cities$score.portsDistsT3, cities$score.biorefDistsT3))

cities$sumBiomassSuitabilityT1 <- 0
cities$sumBiomassSuitabilityT2 <- 0
cities$sumBiomassSuitabilityT3 <- 0

cities$name<-cities$COMM_NAME
cities$dni <- 0
pts <- list()
# load("pts.rda")


df <- list()
ccode <- list("Spain"="ES", "Greece"="EL", "Italy"="IT")
cname <- names(ccode)
names(cname) <- ccode
countries <- cname[ unique(cities$CNTR_CODE) ]

out_dir <- "images"
for(countryn in countries){
  dirout <- sprintf("%s/%s", out_dir, countryn)
  if(!dir.exists(dirout)){
    dir.create(dirout)
  }
  message(countryn)
  city <- cities |> filter(CNTR_ID== ccode[[countryn]])
  city$cityId <- 1:nrow(city)
  message(nrow(city), " cities in ", countryn)
  borders <- nuts |> filter(LEVL_CODE==0)
  provs <- nuts |> filter(LEVL_CODE==2) |>
    filter(CNTR_CODE==ccode[[countryn]])

  r <- terra::rast(sprintf("data/forest/suitability1kmCluster/suitability10-5-40-5-40_%s.tif", countryn))
  if(st_crs(r) != st_crs(provs) ){
    r <- terra::project(r, provs )
  }
  # r <- terra::rast(sprintf("data/AgriSuitabilitPysolo_%s.tif", countryn))
  r[r==0] <- NA


  if(is.null(pts[[countryn]])){
    r.df <- terra::as.data.frame(r, xy=T)
    r.df.sf <- sf::st_as_sf(r.df,coords = c("x", "y"))
    sf::st_crs(r.df.sf) <- terra::crs(r)
    pts[[countryn]] <- r.df.sf
    save(pts, file="pts.rda")
  }

  # plot(r)
  r_extent <- as.polygons(ext(r), crs = crs(r)) |> st_as_sf()
  provs <- provs[st_intersects(provs, r_extent, sparse = FALSE), ]
  city2 <- city[!st_intersects(city, st_union(provs), sparse = FALSE), ]

  cities <- cities %>%
    filter(!COMM_ID %in% city2$COMM_ID)

  rDNI <- terra::rast(sprintf("data/DNI_%s.tif", countryn))
  dniCity <- terra::extract(rDNI, city, ID=F)
  city$dni <- dniCity[,1]
  cities$dni[match(city$COMM_ID, cities$COMM_ID)]  <- city$dni

  sf::write_sf(cities, "data/cities3.gpkg")
  message(nrow(city2), " cities in ", countryn)
  message(nrow(cities), " TOT cities ")



  if(!file.exists(sprintf("%s_cityScore3.gpkg", countryn))){
    city <- city |> filter(dni>1400)
    message(nrow(city))

    # r[r==0] <- NA
    message(countryn)

    # vor <- sf::st_voronoi( st_union(city) |> sf::st_transform(3035)  ) |> sf::st_transform(4326)
    # vor_sf <- st_collection_extract(vor)
    # vor_sf <- st_sf(geometry = vor_sf)
    #
    # vor_sf <- st_join(vor_sf, city, join = st_intersects)

    #sf::write_sf(vor_sf, "vor.sf.gpkg")

    message("Intersection")
    # pts_with_poly <- st_join(pts[[countryn]], vor_sf[,c( "cityId")], join = st_intersects)
    # # pts_with_poly
    sf_use_s2(F)

    message("start lapply")
    ncities <- length(city$cityId)
    nnn <- 0
    suit <- mclapply( city$cityId, function(id){
      nnn <<- nnn+1
      #if(as.integer(id)%% as.integer((ncities/100))==0){
        print(sprintf("%d / %d", nnn,  ncities))
      #}
      ct <- city |> filter(cityId==id)
      pt <- pts[[countryn]] #pts_with_poly |> filter(cityId==id)
      dist<- sf::st_distance(pt, ct)
      dist <- as.numeric(dist)
      w1 <- calculate_weight(dist, maxDist[[1]])
      w2 <- calculate_weight(dist, maxDist[[2]])
      w3 <- calculate_weight(dist, maxDist[[3]])
      c(sumBiomassSuitabilityT1=sum(pt$sum * w1),
           sumBiomassSuitabilityT2=sum(pt$sum * w2),
           sumBiomassSuitabilityT3=sum(pt$sum * w3)
           )
    }  , mc.cores=120)

# FOREST SUITABILITY SUM AT CITY -----
    df <- do.call(rbind,suit)

    scores <- cbind(city, df)

    cities$score.biomassSuitabilityT1[match(scores$COMM_ID, cities$COMM_ID)]  <-
      scores$sumBiomassSuitabilityT1
    cities$score.biomassSuitabilityT2[match(score.biomassSuitability$COMM_ID, cities$COMM_ID)]  <-
      score.biomassSuitability$sumBiomassSuitabilityT2
    cities$score.biomassSuitabilityT3[match(score.biomassSuitability$COMM_ID, cities$COMM_ID)]  <-
      score.biomassSuitability$sumBiomassSuitabilityT3

    sf::write_sf(city, sprintf("%s_cityScore.gpkg", countryn))
  } else {
    city <- sf::read_sf(sprintf("%s_cityScore.gpkg", countryn))
  }


  next


}

dfFinal <- data.table::rbindlist(df, idcol = "Country")
dfFinal$ratio <- dfFinal$Suitability/dfFinal$Cost

dfFinal <- dfFinal |> dplyr::group_by(Country) |>
  dplyr::mutate(UnitCost = row_number(),
                NcitiesNorm = Ncities/max(Ncities)*100,
                CostNorm = Cost/max(Cost)*100,
                SuitabilityNorm = Suitability/max(Suitability)*100,
                title = sprintf("%s (%d cities near biomass)", Country, max(Ncities)))

# dfFinal$UnitCost <- 0:(nrow(dfFinal)-1)%%43

# save(dfFinal, file="dfFinal.rda")

# Build the plot
p1 <- ggplot(dfFinal ) +
  # Map color to a constant with a label
  # geom_line(aes(y = Suitability, color = "Total Suitability"),  linewidth = 1) +
  # geom_line(aes(x=UnitCost, y = ratio, color = "Total Cost"),  linewidth = 1) +
  geom_col(width = 0.7, fill = "grey30", aes(x=UnitCost, y=NcitiesNorm)) +
  facet_wrap(vars(title), nrow=3) +

  labs(
    y = "% of total cities",
    x = "Increase in Unit Cost"
  ) +
  coord_cartesian(xlim = c(0, 30)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )


ggsave(
  filename = sprintf("Cities.png"),
  plot = p1,
  width = 7,
  height = 7,
  dpi = 300,
  bg = "white"   # ensures no transparency
)


p2 <- ggplot(dfFinal ) +

  geom_col(width = 0.7, fill = "grey30", alpha = 0.6, aes(x=UnitCost, y=NcitiesNorm, fill="% of cities")) +
  geom_line(aes(x=UnitCost, y = SuitabilityNorm , color = "Total Gain (normalized)"),  linewidth = 1) +
  geom_line(aes(x=UnitCost, y = CostNorm, color = "Total Loss (normalized)"),  linewidth = 1) +
  # geom_line(aes(x=UnitCost, y = UnitCost/30*100, color = "UnitCost Loss"),  linewidth = 1) +
  # geom_line(aes(x=UnitCost, y = ratio*10, color = "Ratio"),  linewidth = 1) +
  facet_wrap(vars(title), nrow=3, scales ="free_y") +

  labs(
    y = "Normalized Value",
    x = "Unit Cost (UCpp)"
  ) +
  coord_cartesian(xlim = c(0, 20)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
p2
ggsave(
  filename = sprintf("Ratio.png"),
  plot = p2,
  width = 6,
  height = 8,
  dpi = 300,
  bg = "white"   # ensures no transparency
)



# km <- kmeans(sf::st_coordinates(city), centers = 5)
# city$cluster <- as.factor(km$cluster)

# system("ffmpeg -framerate 1 -i images/frame_%03d.png -pix_fmt yuv420p output.mp4")


