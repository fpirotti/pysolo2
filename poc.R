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

cities <- sf::read_sf("data/cities2.gpkg") #|> select(3) #|> sf::st_transform(3035)

nuts <- sf::read_sf("/archivio/shared/geodati/vector/NUTS_2024_all_4326v2.gpkg")

biorefineries <- sf::read_sf("data/Biorefineries.shp") #|> sf::st_transform(3035)
companies <- sf::read_sf("data/Companies.shp")#|> sf::st_transform(3035)
ports <- sf::read_sf("data/Ports_ES.shp")#|> sf::st_transform(3035)

maxDist <- 100000
biorefDist <- nabor::knn(sf::st_coordinates(biorefineries), sf::st_coordinates(cities), k=1 )  # fast radius neighbors
companiesDist <- nabor::knn(sf::st_coordinates(companies), sf::st_coordinates(cities), k=1)  # fast radius neighbors
portsDist <- nabor::knn(sf::st_coordinates(ports), sf::st_coordinates(cities), k=1)  # fast radius neighbors

cities$score.compDists <-  companiesDist$nn.dists[,1]
cities$score.biorefDists <- biorefDist$nn.dists[,1]
cities$score.portsDists <-  portsDist$nn.dists[,1]

cities$name<-cities$name_2

pts <- list()
load("pts.rda")
countries <- unique(cities$Country)

df <- list()
ccode <- list("Spain"="ES", "Greece"="EL", "Italy"="IT")

out_dir <- "images"
for(countryn in countries){
  dirout <- sprintf("%s/%s", out_dir, countryn)
  if(!dir.exists(dirout)){
    dir.create(dirout)
  }
  message(countryn)
  # Read all PNGs in order
  # frames <- list.files(path = dirout, pattern = "frame.*\\.png$", full.names = TRUE)
  # imgs <- image_read(frames)
  # imgs <- image_background(imgs, "white")
  # imgs <- image_quantize(imgs, max = 32, colorspace = "rgb")
  #
  # # # Animate at 10 frames per second
  # animation <- image_animate(imgs, fps = 4,optimize = T, dispose ="background")
  #
  # # Save the GIF
  # image_write(animation, sprintf("animation2_%s.gif", countryn) )
  #
  # next
  # if(countryn!="Spain") next

  city <- cities |> filter(Country==countryn)
  city$cityId <- 1:nrow(city)
  message(nrow(city))
  borders <- nuts |> filter(LEVL_CODE==0)
  provs <- nuts |> filter(LEVL_CODE==2) |>
    filter(CNTR_CODE==ccode[[countryn]])

  r <- terra::rast(sprintf("data/AgriSuitabilitPysolo_%s.tif", countryn))
  r[r==0] <- NA
  r_extent <- as.polygons(ext(r), crs = crs(r)) |> st_as_sf()
  provs <- provs[st_intersects(provs, r_extent, sparse = FALSE), ]
  if(!file.exists(sprintf("%s_cityScore.gpkg", countryn))){

    rDNI <- terra::rast(sprintf("data/DNI_%s.tif", countryn))
    dniCity <- terra::extract(rDNI, city, ID=F)
    city$dni <- dniCity[,1]
    city <- city |> filter(dni>1000)
    message(nrow(city))


    r[r==0] <- NA
    message(countryn)





    if(is.null(pts[[countryn]])){

      r.df <- terra::as.data.frame(r, xy=T)
      r.df.sf <- sf::st_as_sf(r.df,coords = c("x", "y"))
      sf::st_crs(r.df.sf) <- terra::crs(r)
      pts[[countryn]] <- r.df.sf
      save(pts, file="pts.rda")
    }



    vor <- sf::st_voronoi( st_union(city) |> sf::st_transform(3035)  ) |> sf::st_transform(4326)
    vor_sf <- st_collection_extract(vor)
    vor_sf <- st_sf(geometry = vor_sf)

    vor_sf <- st_join(vor_sf, city, join = st_intersects)

    #sf::write_sf(vor_sf, "vor.sf.gpkg")

    message("Intersection")
    pts_with_poly <- st_join(pts[[countryn]], vor_sf, join = st_intersects)
    # pts_with_poly
    sf_use_s2(F)

    message("start lapply")
    ncities <- length(city$cityId)

    suit <- lapply( city$cityId, function(id){
      # browser()
      if(id%%(ncities/100)==0){
        print(sprintf("%d / %d", id,  ncities))
      }
      ct <- city |> filter(cityId==id)
      pt <- pts_with_poly |> filter(cityId==id)
      dist<- sf::st_distance(pt, ct)
      dist <- as.numeric(dist)
      if(sum(dist == 0)>0){
        dist[ dist == 0] <- 1e6
      }
      w <- 1 / (dist^1)
      sum(pt$mean * w)
    } ) #, mc.cores=100)

    city$score.suitability <- unlist(suit)
    sf::write_sf(city, sprintf("%s_cityScore.gpkg", countryn))
  } else {
    city <- sf::read_sf(sprintf("%s_cityScore.gpkg", countryn))
  }

  n <- length(city$score.suitability)
  weights <- city$score.suitability
  costsO <- (city$score.compDists + city$score.biorefDists + city$score.portsDists)
  #normalizzo
  costsO <- costsO / max(costsO)
  weights <- weights / max(weights)
  city$gains <- weights
  city$losses <- costsO
  history <- list()
  iter <- 1
  # score_fun4SANN <- function(cityid) {
  #   cc<-city |> filter(cityId%in%cityid)
  #
  #   d <- as.matrix(dist(st_coordinates(cc)))
  #   d[upper.tri(d, diag=TRUE)] <- NA
  #   penaltyd <- sum(d/1000 < 100, na.rm=TRUE)
  #   score <- sum(cc$gains)
  #   penalty <- sum(cc$losses)
  #   fina <- score*10  - (penalty+penaltyd/max(penaltyd))
  #   cc$score <- fina
  #   history[[iter]] <<- cc[,c("cityId","score")]
  #   iter <<- iter + 1
  #   # message(nrow(cc))
  #   # message(fina)
  #   fina
  # }

  # hist(costsO)
  #hist(weights)

  output <- list(count=c(), totSuitability=c(), totCost=c() )
  count <- 0

  n <- nrow(city)  # number of candidate rows
  gain <- city$gains
  loss <- city$losses
  lambda <- 1
  # k <- 300
  # init <- sample(1:nrow(city), k)

  init <-  rep(T, nrow(city))
  ncities <- NROW(city)
  prop <- function(x) {
    r <- runif(1)
    n <-    x == 1
    nc0 <- which(!n)
    nc1 <- which(n)

    ## force to add a city if only one city left
    if(length(nc1)<2){
      i0 <- sample(nc0, 1)
      x[i0] <- 1
      return(x)
    }
    ## force to remove a city if all cities
    if(length(nc1) >= ncities){
      i1 <- sample(nc1, 1)
      x[i1] <- 0
      return(x)
    }

    i0 <- sample(nc0, 1)
    i1 <- sample(nc1, 1)

    if (r < 0.33 && length(x) > 1) {
      # remove one city
      x[i1] <- 0
    } else if (r < 0.66) {
      # add a city
      x[i0] <- 1
    } else {
      # replace one city
      x[i0] <- 1
      x[i1] <- 0
    }

    x
  }


  trace_values <- numeric(50000)
  trace_ncities <- numeric(50000)
  ## lambda is the unit cost of building a plant ... can be changed
  objective <- function(idx, gain, loss, lambda=1) {

    idx <- as.logical(idx)
    # idx <- pmax(1, pmin(n, round(idx)))   # clip & round to 1..n
    total_gain <- sum(gain[idx])
    total_loss <- sum(loss[idx])


    trace_ncities[[iter]] <<- sum(idx)
    score <- total_gain*100 - total_loss*20 - lambda*trace_ncities[[iter]]

    trace_values[[iter]] <<- score
    if(iter%%100==0){
      cc <- city[idx,]
      cc$score <- score
      cc$ncities <- length(idx)

      history[[as.integer(iter/100)]] <<- cc[,c("cityId","score", "ncities" )]
    }
    iter <<- iter + 1
    # print(score)
    return(-score)  # because optim minimizes by default
  }
  # history <- list()
  iter <- 1
  n <- 3000

  history <- vector("list", n)
  # ## SANN simulated annhealing -------
  res <- optim(
    par = init,
    fn = function(x) objective(x, gain, loss, lambda),
    gr = prop,
    method = "SANN",
    control = list(maxit = 50000, temp = 100, trace = TRUE)
  )

  scale_factor <- max(trace_ncities) / max(trace_values)
  offset <- min(trace_ncities) - min(trace_values) * scale_factor
  p <- ggplot() +
    geom_col(aes(x = 1:length(trace_ncities), y = trace_ncities, fill = "Cities"),
             alpha = 0.5) +
    geom_line(aes(x = 1:length(trace_ncities),  y = trace_values * scale_factor/1.5,
                  color = "Score" ),
              linewidth = 0.31) +
    scale_y_continuous(
      name = "Number of cities",
      sec.axis = sec_axis(~ . / scale_factor, name = "Total score")
    ) +
    scale_fill_manual(
      name = "Cities",
      values = c("Cities" = "steelblue")
    ) +

    scale_color_manual(
      name = "Score",
      values = c("Score" = "red")
    )  + guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    xlab("Iteration") +
    theme_minimal()


  ggsave(
    sprintf("spainSANN.png", i),
    plot = p,
    width = 8, height = 5, dpi = 120
  )
  # best_x <- as.logical(round(res$par))
  # best_rows <- city[res$par, ]

  # res <- optim(
  #   par = init,
  #   fn = function(x) -score_fun4SANN(round(x)),
  #   gr = prop,
  #   method = "SANN",
  #   control = list(maxit = 10000, temp=50)
  # )

  ########## plot ----
  frames <- data.table::rbindlist(history,idcol = "iter")

  library(ggplot2)
  library(gganimate)
  frames <- sf::st_set_geometry(frames, value = frames$geom)
  frames$iter <- as.numeric(frames$iter)
  frames$score <- as.numeric(frames$iter)

  frames$scoreTot <- (frames$score - min(frames$score))/ diff(range(frames$score))
  library(viridis)
  p <- ggplot() +
     # geom_spatraster(data = r, na.rm = T) +
    # scale_fill_whitebox_c(
    #   palette = "viridi",
    #   n.breaks = 12,
    #   guide = guide_legend(reverse = TRUE)
    # ) +
    geom_sf(data =  provs,    fill = NA) +
    geom_sf(data = frames ,
               # aes(x=x, y=y),
            aes(color = scoreTot),
               size=2) +
    scale_color_viridis(option = "turbo", direction=-1) +
    coord_sf() +
    theme_minimal()+
    transition_manual(iter) +
    shadow_null() +
    labs(title = "Iteration: {current_frame}/10000")

  # print(p)
  animate(p, fps=10,
          width = 800,
          height = 600,
          renderer = ffmpeg_renderer("optimization.mp4") )

  opt_coords <- matrix(res$par, ncol=2, byrow=TRUE)

  opt_pts <- vect(opt_coords, type="points", crs=crs(r))

  for(cost in (1:50)){
    count <- count+1
    costs <- costsO * cost / 10

    model <- MIPModel() %>%
      add_variable(x[i], i = 1:n, type = "binary") %>%
      set_objective(sum_expr(weights[i]*x[i] - costs[i]*x[i], i = 1:n), "max") %>%
      # optionally: add_constraint(sum_expr(costs[i]*x[i], i = 1:n) <= B)
      solve_model(with_ROI(solver = "glpk"))

    message("solution loop")
    solution <- get_solution(model, x[i]) #%>% filter(value > 0.5)
    city$selectedCities  <- solution$value
    city[[sprintf("restrictionLevel%02d", count)]]  <- solution$value
    p <- ggplot() +
      geom_sf(data=provs, fill = NA) +
      geom_spatraster(data = r, na.rm = T) +
      scale_fill_whitebox_c(
        palette = "viridi",
        n.breaks = 12,
        guide = guide_legend(reverse = TRUE)
      ) +
      geom_sf(data=city[city$selectedCities==0,"selectedCities"], size = 0.5, color="#000000") +
      geom_sf(data=city[city$selectedCities==1,"selectedCities"], fill=NA, color = "#ff000099", size = 3) +
      ggtitle(sprintf("Unit cost=%02d - N. Selected Cities=%d", count-1, sum(city$selectedCities) )) +
      theme_minimal() +
      theme(legend.position = "none")

    ggsave(
      filename = sprintf("%s/frame_%02d.png", dirout, count),
      plot = p,
      width = 6, height = 5, dpi = 100
    )
    output$count <- c(output$count, sum(city$selectedCities) )
    output$totSuitability <- c(output$totSuitability, sum(weights*city$selectedCities) )
    output$totCost <- c(output$totCost, sum(costs*city$selectedCities) )
    message(sum(city$selectedCities))
  }

  sf::write_sf(city, sprintf("%sSuitableCities.gpkg", countryn))
  ## numero di città candidate


  # Create a data frame with both series
  df[[countryn]] <- data.frame(
    Suitability = output$totSuitability,
    Cost = output$totCost,
    Ncities=output$count
  )


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


