library(mapgl)
library(viridisLite)
library(htmlwidgets)
library(tidyterra)
library(sf)
library(ggplot2)
library(viridis)
library(classInt)
library(scales)
library(ompr)
library(ggtext)


pkgs <- c(  "patchwork", "tidyterra", "ggplot2",    "this.path")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    warning("Package \033[1m", p, "\033[0m not found... installing it.")
    install.packages(p)
  }
  message("Package \033[1m", p, "\033[0m  found... loading it.")
  library(p, character.only = TRUE, quietly = TRUE)
}


nuts <- sf::read_sf("/archivio/shared/geodati/vector/NUTS_2024_all_4326v2.gpkg")
maxDist <- c(50000, 100000, 150000)


calculate_weight <- function(distance, max_dist ) {
  weight <- 1 - (distance / max_dist)
  return(pmax(0, weight)) # pmax is the vectorized version of max()
}

distanceMatrix <- function(id, countryn){
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
}


plotCluster <- function(dat, cluster_bounds, ret=FALSE, ccode){

  vals <- na.omit(dat$sumBiomassSuitabilityT1n)
  breaks_fisher <- classIntervals(vals, n = 12, style = "fisher")$brks
  cluster_centers <- st_centroid(cluster_bounds)
  provs <- nuts |> filter(LEVL_CODE==2) |>
    filter(CNTR_CODE==ccode ) |>
    st_filter(dat, .predicate = st_intersects)

  p <- ggplot() +
    geom_sf(data =  dat,
            aes(color = round(k, 2) ) ) +
    # geom_sf(data = cluster_bounds, fill = NA, color = "black", lwd=2, linetype = "dashed") +
    # geom_sf(data = cluster_centers, color = "red", size = 6, shape = 18) +

    geom_sf(data =  provs,    fill = NA) +
    labs(subtitle = "Diamonds represent K-means centers; Dashed lines are cluster bounds") +
    scale_color_stepsn(
      colors = viridis_pal(option = "turbo")(12),
      # breaks = round(breaks_fisher,2),
      # values = round(rescale(breaks_fisher), 2), # Rescales breaks to 0-1 range
      name = "Suitability"
    ) +
    coord_sf() +
    theme_minimal()+
    # transition_manual(iter) +
    # shadow_null() +
    labs(title = "Iteration: {current_frame}/10000")

  if(!ret){
    png("plotCluster", res=100, width=1000, height=1000)
     print(p)
    dev.off()
  }
 p
}




plotCitiesBiomass <- function(dat, data, idBest, breaks_fisher,text,  lowres=T){


  provs <- nuts |> filter(LEVL_CODE==2) |>
    filter(CNTR_CODE==data$code ) |>
    st_filter(dat, .predicate = st_intersects)

  p1 <- ggplot() +
    geom_sf(data =  dat,
            aes(color = tBiomInit ), alpha=0.7 ) +

    geom_sf(data =  provs,    fill = NA) +
    # geom_sf(data = cluster_bounds, fill = NA, color = "black", lwd=2, linetype = "dashed") +
    geom_sf(data = dat[idBest,],  shape = 21,
            size = 10,
            fill = "#ff000000",
            color = "red",
            stroke = 1) +
    geom_sf(data = dat[idBest[[length(idBest)]],],  shape = 21,
            size = 6,
            fill = "#ff000000",
            color = "red",
            stroke = 1) +

    scale_color_stepsn(

      colors = viridis_pal(option = "viridis")(6),
       breaks = round(breaks_fisher,2),
       values = round(rescale(breaks_fisher), 2), # Rescales breaks to 0-1 range
      name = "Suitability"
    ) +
    coord_sf() +
    theme_minimal(  ) #+
    # transition_manual(iter) +
    # shadow_null() +
    #labs(title =  "City") )

  vals <- na.omit(dat$tBiom)
  breaks_fisher2 <- classIntervals(vals, n = 6, style = "pretty")$brks

  p4 <- ggplot() +
    geom_sf(data =  dat,
            aes(color = tBiom ), alpha=0.7 ) +

    geom_sf(data =  provs,    fill = NA) +
    geom_sf(data = dat[idBest[[length(idBest)]],],  shape = 21,
            size = 10,
            fill = "#ff000000",
            color = "red",
            stroke = 1) +
    geom_sf(data = dat[idBest[[length(idBest)]],],  shape = 21,
            size = 6,
            fill = "#ff000000",
            color = "red",
            stroke = 1) +
    scale_color_stepsn(

      colors = viridis_pal(option = "viridis")(6),
      breaks = round(breaks_fisher2,2),
      values = round(rescale(breaks_fisher2), 2), # Rescales breaks to 0-1 range
      name = "Suitability"
    ) +
    coord_sf() +
    theme_minimal(  )
  # transition_manual(iter) +
  # shadow_null() +
  # labs(title =  "Best City Choice",
  #      subtitle = paste(text,collapse=" ") )

  final_dashboard <- p1  / p4 +
    plot_annotation(
      title = sprintf("Optimal Pyrolysis Plant Location <span style='font-family:serif;'>&tau;</span> = %s km ", data$kmThreshold),
      subtitle = gsub("\n","<br>", text),
      theme = theme(plot.title =  element_markdown(),
                    plot.subtitle =  element_markdown() )
    )

    outd <- file.path("plots", data$country, paste0(data$kmThreshold, "km" ))
    dir.create(outd,recursive = T,showWarnings = F)
    # png(file.path(outd, sprintf("%04d_plotCities%s_v2.png", length(idBest), ccode)),
    #     res=ifelse(lowres, 96, 300),
    #     width=ifelse(lowres, 700, 2000),
    #     height=ifelse(lowres,700, 2000) )
    #   print(p)
      ggsave(
        filename = file.path(outd, sprintf("%04d_plotCities%s.png",
                                           length(idBest), data$code)),
        plot = final_dashboard,
        width = ifelse(lowres, 1350, 1350*2),
        height = ifelse(lowres, 1080*2, 1080*4),
        units = "px",
        # width = ifelse(lowres, 7, 20),
        # height =ifelse(lowres, 7, 20),
        dpi = ifelse(lowres, 200, 300),
        bg = "white"   # ensures no transparency
      )

  # p
}


plotHTML <- function(result, tit="Spain"){
  # result <- candidateCities
  # Get selected facility locations
  selected <- result$facilities |>
    filter(.selected) |>
    mutate(id = as.character(COMM_NAME))

  # Color demand points by their assigned facility
  demand_colored <- result$demand |>
    mutate(.facility = as.character(.facility))

  # Map the results
  mg <- maplibre(bounds = candidateCities) |>
    add_circle_layer(
      id = "Candidate Cities",
      source = candidateCities,
      circle_color = "black",
 circle_radius = interpolate(
      column = "zoom",            # This refers to the map's current scale
      values = c(1, 18),          # From zoom level 5 to 18
      stops = c(1, 3000)            # Size scales from 2px to 30px
    ),
      circle_opacity = 0.7
    ) |>
    add_circle_layer(
      id = "Biomass",
      source = pts$Spain ,
      circle_color = interpolate("sum", values=(1:5)/5,
                                 stops = viridisLite::viridis(5) ),
      circle_radius = 3,
      circle_opacity = 0.7,
      cluster_options = cluster_options(max_zoom = 7, circle_opacity=0.5)
    )
  htmlwidgets::saveWidget(mg, file=sprintf("%s.html", tit) )
}
