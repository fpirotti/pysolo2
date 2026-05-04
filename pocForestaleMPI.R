library(deldir)   # Voronoi tessellation
library(sf)       # Spatial operations
library(ompr)     # Optimization modeling
library(ompr.roi) # Solver interface
library(ROI.plugin.glpk)  # GLPK solver
library(Matrix)
library(RcppAnnoy)
# install.packages("ROI.plugin.highs")
library(ROI.plugin.highs)
# 1. Genera un set di punti con costi casuali
set.seed(42)


# 2. Funzione per calcolare la matrice di copertura
compute_cover_matrix <- function(citiesPts, biomassSuitPts, radius) {
  browser()
  dist_matrix <- fields::rdist(citiesPts, biomassSuitPts, compact=T)
  # storage.mode(dists) <- "integer"
  # idx <- which(dists > maxDist[[i]])
  # dists[idx] <- NA
  # wm <- (1 - dists/maxDist[[i]])*pt$sum
  # dist_matrix <- (st_distance(citiesPts, biomassSuitPts))
  # class(dist_matrix)<-"numeric"
  nearest_idx <- apply(dist_matrix, 2, which.min)  # Trova il punto primario più vicino

  cover_matrix <- matrix(0, nrow = nrow(citiesPts), ncol = nrow(biomassSuitPts))
  cover_matrix2 <- ifelse(dist_matrix <= radius, 1, 0)  # Coverage within given radius


  for (j in 1:nrow(biomassSuitPts)) {
    cover_matrix[nearest_idx[j], j] <- 1*cover_matrix2[nearest_idx[j], j]
  }
  # browser()
  return(cover_matrix)
}

# 3. Funzione per ottimizzare il set primario minimizzando il costo
optimize_primary_points <- function(citiesPts, biomassSuitPts, lambda = 0.99, radius=100000) {

  citiesPts <- ct
  biomassSuitPts <-pt
  # system.time({ dist<- sf::st_distance(pt, candidateCities , tolerance=thresholds[[i]] ) })
  biomass = sf::st_coordinates(pt)
  cities = sf::st_coordinates(ct)

  cover_matrix <- compute_cover_matrix(cities, biomass, radius)

  costs <- citiesPts$costTotT2*10
  weights <- biomassSuitPts$sum


  model <- MIPModel() %>%
    add_variable(x[i], i = 1:nrow(citiesPts), type = "binary") %>%
    add_variable(y[j], j = 1:nrow(biomassSuitPts), type = "binary") %>%

    # Objective: Minimize primary set cost and maximize secondary set weights
    set_objective(
      sum_expr(lambda * costs[i] * x[i], i = 1:nrow(citiesPts)) -
        sum_expr((1 - lambda) * weights[j] * y[j], j = 1:nrow(biomassSuitPts)),
      sense = "min"
    ) %>%

    # Constraint: A secondary point is covered if at least one primary point covering it is selected
    add_constraint(
      y[j] <= sum_expr(cover_matrix[i, j] * x[i], i = 1:nrow(citiesPts)),
      j = 1:nrow(biomassSuitPts)
    )

  system.time({ result.glpk <- solve_model(model, with_ROI(solver = "glpk")) })
  system.time({ result.highs <- solve_model(model, with_ROI(solver = "highs", control = list(threads = 16, parallel = "on") )) })


  if (result.glpk$status == "success") {
    selected_primary <- which(result.glpk$solution[1:nrow(citiesPts)] > 0.5)
    selected_secondary <- which(result.glpk$solution[(nrow(citiesPts) + 1):(nrow(citiesPts) + nrow(biomassSuitPts))] > 0.5)

    df1<- list(
      citiesPts = citiesPts[selected_primary, ],
      biomassSuitPts = biomassSuitPts[selected_secondary, ]
    )
    selected_primary <- which(result.highs$solution[1:nrow(citiesPts)] > 0.5)
    selected_secondary <- which(result.highs$solution[(nrow(citiesPts) + 1):(nrow(citiesPts) + nrow(biomassSuitPts))] > 0.5)

    df2<- list(
      citiesPts = citiesPts[selected_primary, ],
      biomassSuitPts = biomassSuitPts[selected_secondary, ]
    )

    return(list(glpk = df1, highs=df2))
  } else {
    stop("Solution not found!")
  }
}


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
    citiesPts <- candidateCities |> sf::st_transform(3035)
    biomassSuitPts <- pts[[cn]] |> sf::st_transform(3035)

    # 4. Esegui l'ottimizzazione
    bestCities[[cn]][[ paste(as.character(maxDist[[i]]/1000), "km") ]]  <- optimize_primary_points(citiesPts, biomassSuitPts)
  }
}
# # 5. Visualizza i risultati
plot(citiesPts$geom, col = "gray", pch = 19)
plot(optimal_citiesPts$highs$biomassSuitPts$geom, col = "blue", pch = ".", add=T)
plot(optimal_citiesPts$highs$citiesPts$geom, col = "red", pch = 8, cex = 1.5, add=T)
sum(optimal_citiesPts$glpk$biomassSuitPts$sum) -
  sum(optimal_citiesPts$glpk$citiesPts$costTotT2*40)
# # Stampa il costo totale
# cat("Costo totale selezionato:", sum(optimal_citiesPts$cost), "\n")



cat('Distance Calculation: It uses fields::rdist to calculate the straight-line distance between every city and every biomass point.

The "Nearest Neighbor" Rule: The line apply(dist_matrix, 2, which.min) is a very strict constraint. It assumes a biomass point can only be assigned to its single closest city.

The Radius Filter: Even if a city is the closest, if it is further away than the radius (e.g., 100km), it cannot "cover" that biomass.

Result: You get a matrix where a 1 means "City A can harvest Biomass B," and a 0 means it cannot.')
