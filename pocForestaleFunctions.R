library(mapgl)
library(viridisLite)
library(htmlwidgets)

plotHTML <- function(result, tit){

  # Get selected facility locations
  selected <- result$facilities |>
    filter(.selected) |>
    mutate(id = as.character(COMM_NAME))

  # Color demand points by their assigned facility
  demand_colored <- result$demand |>
    mutate(.facility = as.character(.facility))

  # Map the results
  mg <- maplibre(bounds = selected) |>

    add_circle_layer(
      id = "demand",
      source = demand_colored,
      circle_color = match_expr(
        column = ".facility",
        values = selected$id,
        stops = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")
      ),
      circle_radius = 4,
      circle_opacity = 0.7
    ) |>
    add_circle_layer(
      id = "facilities",
      source = selected,
      circle_color = match_expr(
        column = "id",
        values = selected$id,
        stops = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")
      ),
      circle_radius = 10,
      circle_stroke_color = "white",
      circle_stroke_width = 2
    )

  htmlwidgets::saveWidget(mg, file=sprintf("%s.html", tit) )
}
