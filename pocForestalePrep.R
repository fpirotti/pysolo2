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

cities <- sf::read_sf("data/cities.gpkg") #|> select(3) #|> sf::st_transform(3035)
cities <- cities |> select(!c(LAT, LON, TRUE_FLAG, COMM_NAME, NSI_CODE, NAME_NSI))
# Download city center  ----
## points for Spain, Italy, and Greece - uncomment if needed
# library(giscoR)
# countries <- c("ES", "IT", "EL") # Note: Greece is EL in Eurostat
# lau_data <- gisco_get_communes(country = c("ES", "IT", "EL"), spatialtype="LB")
# sf::write_sf(lau_data, "data/cities.gpkg") #|> select(3) #|> sf::st_transform(3035)

 nuts <- sf::read_sf("/archivio/shared/geodati/vector/NUTS_2024_all_4326v2.gpkg")
rk <- matrix(rep(1, 81), 9, 9)
biorefineries <- sf::read_sf("data/forest/Main Refineries_FINAL.shp") #|> sf::st_transform(3035)

ports <- sf::read_sf("data/forest/Official Ports_FINAL.shp")#|> sf::st_transform(3035)

maxDist <- c(50000, 100000, 150000)


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

cities$costTotT1 <-  (cities$score.portsDistsT1+cities$score.biorefDistsT1)/2
cities$costTotT2 <-  (cities$score.portsDistsT2+cities$score.biorefDistsT2)/2
cities$costTotT3 <-  (cities$score.portsDistsT3+cities$score.biorefDistsT3)/2

biom <- terra::rast("/archivio/shared/geodati/raster/Biomass/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2022-fv6.0.nc")

cities$name<-cities$COMM_NAME
cities$dni <- 0
pts <- list()
# load("pts.rda")


df <- list()
ccode <- list("Greece"="EL", "Italy"="IT", "Spain"="ES")
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
    r.df <- terra::as.data.frame(r, xy=T, cells=F)
    r.df.sf <- sf::st_as_sf(r.df,coords = c("x", "y"))
    message("Start extraction..." , Sys.time())

    nr <- terra::crop( biom , r)
    nr2agb <- terra::aggregate(nr[[1]], 11, fun="sum")
    nr2sd <- terra::aggregate(nr[[2]], 11, fun=function(x){ sum(x^2)^0.5 }, cores=20)
    biomassNew <- c(nr2agb, nr2sd)
    names(biomassNew) <- c("agb", "agb_sd")
    agbs <- terra::extract(biomassNew, r.df.sf)

    r.df.sf2 <- cbind(r.df.sf, agbs)

    sf::st_crs(r.df.sf2) <- terra::crs(r)
    pts[[countryn]] <- r.df.sf2
    save(pts, file="pts.rda")
  }
  next
  # plot(r)
  r_extent <- as.polygons(ext(r), crs = crs(r)) |> st_as_sf()
  provs <- provs[st_intersects(provs, r_extent, sparse = FALSE), ]
  city2 <- city[!st_intersects(city, st_union(provs), sparse = FALSE), ]
  city  <- city[st_intersects(city, st_union(provs), sparse = FALSE), ]

  cities <- cities %>%
    filter(!COMM_ID %in% city2$COMM_ID)

  rDNI <- terra::rast(sprintf("data/DNI_%s.tif", countryn))
  dniCity <- terra::extract(rDNI, city, ID=F)
  city$dni <- dniCity[,1]
  mm1<- which(cities$COMM_ID%in%city$COMM_ID)
  cities$dni[mm1]  <- city$dni

  sf::write_sf(cities, "data/cities3.gpkg")
  message(nrow(city2), " cities in ", countryn)
  message(nrow(cities), " TOT cities ")

  if(!file.exists(sprintf("%s_cityScore3.gpkg", countryn))){

    message(countryn)

    sf_use_s2(F)

    message("start lapply")
    ncities <- length(city$cityId)
    nnn <- 0
    suit <- mclapply( city$cityId, distanceMatrix , mc.cores=120)

# FOREST SUITABILITY SUM AT CITY -----
    dfc <- do.call(rbind,suit)
    scores <- cbind(city, dfc)
    # mm <- which(cities$COMM_ID %in%scores$COMM_ID)
    # cities$sumBiomassSuitabilityT1[mm]  <-  scores$sumBiomassSuitabilityT1
    # cities$sumBiomassSuitabilityT2[mm]  <-  scores$sumBiomassSuitabilityT2
    # cities$sumBiomassSuitabilityT3[mm]  <-  scores$sumBiomassSuitabilityT3

    # sf::write_sf(cities, "data/cities3.gpkg")
    scores$sumBiomassSuitabilityT1n <- scores$sumBiomassSuitabilityT1 /
      max(scores$sumBiomassSuitabilityT1)
    scores$sumBiomassSuitabilityT2n <- scores$sumBiomassSuitabilityT2 /
      max(scores$sumBiomassSuitabilityT2)
    scores$sumBiomassSuitabilityT3n <- scores$sumBiomassSuitabilityT3 /
      max(scores$sumBiomassSuitabilityT3)
    sf::write_sf(scores, sprintf("%s_cityScore3.gpkg", countryn))
  } else {
    message(countryn, " already done!")
    scores <-  sf::read_sf(sprintf("%s_cityScore3.gpkg", countryn))
    scores$sumBiomassSuitabilityT1n <- scores$sumBiomassSuitabilityT1 /
      max(scores$sumBiomassSuitabilityT1)
    scores$sumBiomassSuitabilityT2n <- scores$sumBiomassSuitabilityT2 /
      max(scores$sumBiomassSuitabilityT2)
    scores$sumBiomassSuitabilityT3n <- scores$sumBiomassSuitabilityT3 /
      max(scores$sumBiomassSuitabilityT3)

    sf::write_sf(scores, sprintf("%s_cityScore3.gpkg", countryn))
    # city <- sf::read_sf(sprintf("%s_cityScore3.gpkg", countryn))
  }
  df[[countryn]] <-  sf::read_sf(sprintf("%s_cityScore3.gpkg", countryn))

}

dfFinal <- data.table::rbindlist(df, idcol = "Country")
sf::write_sf(dfFinal, "data/citiesFinal.gpkg")
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


