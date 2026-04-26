

cities <- sf::read_sf("data/cities3.gpkg") #|> select(3) #|> sf::st_transform(3035)
df <- list()
ccode <- list( "Greece"="EL", "Italy"="IT", "Spain"="ES")
cname <- names(ccode)
names(cname) <- ccode
countries <- cname[ unique(cities$CNTR_CODE) ]

out_dir <- "images"
for(countryn in names(ccode)){
  # dirout <- sprintf("%s/%s", out_dir, countryn)
  # if(!dir.exists(dirout)){
  #   dir.create(dirout)
  # }
  message(countryn)
  candidateCities <- cities |> filter(CNTR_ID==ccode[[countryn]])
  message(nrow(demand_pts), " cities ")
  result_pmedian <-   p_median(
    demand = pts[[countryn]][1:100,],
    facilities = candidateCities,
    n_facilities = 44,
    weight_col = "sum"
  )

  result_pmedian <-  cflp(
      demand = pts[[countryn]][1:100,],
      facilities = candidateCities,
      n_facilities = 44,  # Let solver determine optimal number
      weight_col = "sum" ,
      capacity_col = "capacity",
      facility_cost_col = "costTotT1"
    )

}


# km <- kmeans(sf::st_coordinates(city), centers = 5)
# city$cluster <- as.factor(km$cluster)

# system("ffmpeg -framerate 1 -i images/frame_%03d.png -pix_fmt yuv420p output.mp4")


