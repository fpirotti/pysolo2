library(sf)
library(parallelDist)
library(tidyverse)
library(fields)
library(tictoc)
library(this.path)
setwd(this.path::this.dir())
library(classInt)
source("pocForestaleFunctions.R")
ccode <- list( "Greece"="EL", "Italy"="IT", "Spain"="ES")
cname <- names(ccode)
names(cname) <- ccode
load("pts.rda")


bestCities<-list()
for(cn in names(ccode)){
  # dirout <- sprintf("%s/%s", out_dir, cn)
  # if(!dir.exists(dirout)){
  #   dir.create(dirout)
  # }
  cat("## ", cn, "\t\n\n")
  # if(cn!="Spain") next
  candidateCitiest <- sf::read_sf(sprintf("%s_cityScore3.gpkg", cn))


  bestCities[[cn]] <- list()
  candidateCities <- candidateCitiest |> filter(dni > 1400)
  cat("### ", nrow(candidateCities), "/",nrow(candidateCitiest), " candidate cities after DNI filter at 1400\n\n")

  for(i in 1:3){
    ct <- candidateCities |> sf::st_transform(3035)
    pt <- pts[[cn]] |> sf::st_transform(3035)
    # system.time({ dist<- sf::st_distance(pt, candidateCities , tolerance=thresholds[[i]] ) })
     biomass = sf::st_coordinates(pt)
     cities = sf::st_coordinates(ct)

    # system.time({
    dists <- fields::rdist(biomass, cities, compact=T)
    storage.mode(dists) <- "integer"
    idx <- which(dists > maxDist[[i]])
    dists[idx] <- NA
    wm <- (1 - dists/maxDist[[i]])*pt$sum

    cat("## First run we assign all neighbour biomass  to all cities.\n
Other iterations remove the picked city and destroy all the biomass that it takes.\n
It is faster as it shows the new assigned biomass after the added city grabs the biomass around it.\n
As we add cities, more biomass is taken and the higher the cost, so the lower the chance that the \n
next best city will be worth adding.")
# })
    candidateCities$tBiomInit <- colSums( wm ,na.rm = T)
    nCitiesTested <- 20
    bestC<- c()
    costWeight <- c(20, 40)
    breaks_fisher <- NULL

    if(!is.null(bestCities[[cn]][[ paste(as.character(maxDist[[i]]/1000), "km") ]])) next
    bestCities[[cn]][[ paste(as.character(maxDist[[i]]/1000), "km") ]] <-  vector("list", nCitiesTested)
    for(gg in 1:nCitiesTested){
        message("\n====>>", gg)

        tic("## Sum and get max")
        candidateCities$tBiom <- colSums( wm ,na.rm = T)
        if( sum(candidateCities$tBiom, na.rm=T) < 1){
          message("No more biomass!\n\n==========\n\n")
          toc()
          break
        }
        if(is.null(breaks_fisher)){
          vals <- na.omit(candidateCities$tBiomInit)
          breaks_fisher <- classIntervals(vals, n = 6, style = "pretty")$brks
        }
        wmax <- which.max(candidateCities$tBiom)
        if(length(wmax)==0){
          warning("Reached no more cities with max values!")
          break
        }
        toc()

        tic("## Calculate Marginal Costs")
        ## picked city removes biomass so the wm matrix has to be updated

        bestC <- c(bestC, wmax)

        marCost1 <- candidateCities$tBiom[wmax] - candidateCities[[sprintf("costTotT%d", i)]][wmax]*costWeight[[1]]
        cumCost1 <- sum(candidateCities$tBiom[bestC] - candidateCities[[sprintf("costTotT%d", i)]][bestC]*costWeight[[1]])


        marCost2 <- candidateCities$tBiom[wmax] - candidateCities[[sprintf("costTotT%d", i)]][wmax]*costWeight[[2]]
        cumCost2 <- sum(candidateCities$tBiom[bestC] - candidateCities[[sprintf("costTotT%d", i)]][bestC]*costWeight[[2]])


        num.free.cells <- !is.na( wm[ , wmax])
        dat <- list(
          country=cn,
          iteration=gg,
          code=ccode[[cn]],
          city=candidateCities$NAME_LATN[[wmax]],
          cityId=candidateCities$COMM_ID[[wmax]],
          suitabilityUnits = sum(wm[ , wmax], na.rm=T),
          n1kmsqPixelsWithAvailableBiomass = sum( num.free.cells ),
          totBiomass  = sum(pt$agb[ num.free.cells ]),
          totBiomassSD  = (sum(pt$agb_sd[ num.free.cells ]^2))^0.5,
          avgBiomassXha  = sum(pt$agb[ num.free.cells ]) / sum( num.free.cells ) / 100,
          avgBiomassSD_Xha  = (sum(pt$agb_sd[ num.free.cells ]^2))^0.5 / sum( num.free.cells )^0.5 / 10,
          cumGain = cumCost2,
          marGain  = marCost2,
          cumGain1 = cumCost1,
          marGain1 = marCost1,
          cumGain2 = cumCost2,
          marGain2 = marCost2,
          kmThreshold = as.integer(maxDist[[i]]/1000)
        )

        # str <- c("Adding **",  dat$city, "** which harvests ",
        #          round(dat$suitabilityUnits,2) , " biomass units (u)  (",dat$npixels," pixels)\n",
        #          "Marginal gain: ", round(dat$marCost,2), " u, cumulative gain:", paste0("<span style='color:",
        #                                                                                  ifelse(dat$cumCost < 0 , "red", "blue"),
        #                                                                                  ";'>",
        #                                                                                         round(dat$cumCost,2), "</span> u" )
        #          )

        bestCities[[cn]][[ paste(as.character(maxDist[[i]]/1000), "km") ]][[gg]] <- dat

        toc()

        # tic("## Plot image")
        # na.omit(candidateCities$tBiom)
        # plotCitiesBiomass( candidateCities, dat, bestC,
        #                    breaks_fisher, paste(str, collapse="") )
        #
        # toc()
        # cat(str, "\n\n" )
        tic("## Recalculate matrix")
        r2rem <- which(!is.na(wm[, wmax]))
        wm[r2rem,] <- NA
        wm[  , wmax] <- NA

        toc()

    }

  }
  save(bestCities, file="bestcities.rda")

  # candidateCities |> filter()
  # rr <- terra::rast(
  #   sprintf("data/forest/suitabilityStackEffective/suitabilityStackEffective_%s.tif", cn)
  # )

}

load(file="bestcities.rda")



library(dplyr)
library(purrr)

final_table <- bestCities |>
  list_flatten() |> # Removes the Country level
  list_flatten() |> # Removes the Distance level
  map_dfr(as.data.frame)

# names(final_table)[7:12] <- names(dat)[7:12]


writexl::write_xlsx(final_table, "pysolo2.xlsx" )

library(ggplot2)
ggplot(final_table) +
   geom_line(aes(x=iteration, y=cumGain1,  col="Cost1x") )  +
  geom_line(aes(x=iteration, y=cumGain2 ,  col="Cost2x" ) )  +

  xlab("Number of progressive best cities") +
  ylab("Cum. Gain Units (Gain - Costs)") +
  geom_hline(yintercept = 0 ) +
  theme_minimal() +
  scale_color_manual(name = "Scenario",
                     values = c("Cost1x" = "blue", "Cost2x" = "red")) +
  ggtitle("Optimal City Identification via Progressive Local Maximum")+
  facet_grid( as.factor(kmThreshold) ~ country, scales = "free", space = "free")

ggplot(final_table) +
  geom_line(aes(x=harvestUnits    , y=cumGain1,  col="Cost1x") )  +
  geom_line(aes(x=harvestUnits    , y=cumGain2 ,  col="Cost2x" ) )  +

  xlab("Number of progressive best cities") +
  ylab("Cum. Gain Units (Gain - Costs)") +
  geom_hline(yintercept = 0 ) +
  theme_minimal() +
  scale_color_manual(name = "Scenario",
                     values = c("Cost1x" = "blue", "Cost2x" = "red")) +
  ggtitle("Optimal City Identification via Progressive Local Maximum")+
  facet_grid( as.factor(kmThreshold) ~ country, scales = "free", space = "free")


# install.packages("flextable")
 # install.packages("officer")
library(flextable)
library(tidyverse)

# Prepare data: must be sorted by the grouping variable
grouped_df <- final_table |> select(country, city, code, iteration, kmThreshold,
                                    MargGain=cumGain, areakm2=npixels) |>
  arrange(country) |>
  as_grouped_data(groups = c("country", "kmThreshold") )

grouped_df$MargGain <- as.integer(grouped_df$MargGain)

ina <- is.na(grouped_df$city)
grouped_df$city <- paste0(grouped_df$city, " (", grouped_df$code,")")
grouped_df$city[ina] <- NA
grouped_df$code <- NULL

# grouped_df$MargGainXarea <- as.integer(grouped_df$MarginalGain) / grouped_df$areakm2
# Create a table where the country name sits in its own row above the data
ft_grouped <- flextable(grouped_df) |>
  theme_vanilla()  |>
  fontsize(i = 1, size = 9, part = "header") |>
  fontsize(size = 8, part = "all") |>
  # Make the country rows look distinct
  bold(i = ~ !is.na(country), part = "body") |>
  bg(i = ~ !is.na(country), bg = "#eeeeee", part = "body")  |>
  autofit()

print(ft_grouped, preview = "docx")


library(av)

# Get a list of your images
images <- list.files(path = "plots/Spain/50km", pattern = "*.png", full.names  = TRUE)

# Sort them to ensure iteration 1 comes before iteration 10
# (Standard R sorting might put 10 after 1, so mixedsort is helpful)
images <- gtools::mixedsort(images)

# Create the video
av::av_encode_video(
  input = images,
  output = "optimization_results.mp4",
  framerate = 1 # 2 images per second (adjust as needed)
)
