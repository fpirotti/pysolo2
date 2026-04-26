## must have a frames table with iter column
library(parallel)
library(ggplot2)
library(SpatialKDE)
library(viridis)

df <- frames
times <- unique(frames$iter)
n <- length(times)
qq <- quantile(r$mean[], c(0.1,0.9), na.rm=T)
r2 <- terra::aggregate(r, 12, fun=mean, na.rm=TRUE)
score_range <- range(df$scoreTot, na.rm = TRUE)

makeHeat <- function(df.iter){
  radius <- 100000
  df.iter_ <- sf::st_transform(df.iter, 3035)
  bb<-provs |> sf::st_transform(3035)  |> st_union() |> st_bbox() |>
    st_as_sfc() |> st_buffer(radius) |>
    st_as_sf()
  r <- create_raster(bb, cell_size = 2000)

  kde_raster <- kde(
    points = df.iter_,
    band_width = radius,
    weights = df.iter_$scoreTot, scaled = T,
    grid = r #, quiet=T
  )
  # plot(kde_raster)
  tr <- terra::rast(kde_raster)
  tr[tr==0] <- NA
  tr
}

res <- mclapply((1:50), #seq_along(times),
                function(i) {
  df.iter <- df[df$iter == times[i], ]
  rheat <- suppressMessages(makeHeat(df.iter))

  p <- ggplot(df, aes(color = scoreTot)) +
     geom_spatraster(data =  rheat , na.rm = T ) +
    # scale_fill_whitebox_c( #breaks = round(seq(qq[[1]], qq[[2]], length.out=12)),
    #                       palette = "viridi") +

    scale_fill_viridis_c(option = "inferno",
                         na.value = tail(inferno(256),1),
                         direction = -1, limits= score_range) +
    geom_sf(data =  provs,  color="grey" ,   fill = NA) +
    # geom_sf(data = df.iter ,
    #         # aes(x=x, y=y),
    #
    #         size=2) +
    # scale_color_viridis(option = "inferno", direction=-1,   limits = score_range) +
    coord_sf() +
    theme_minimal()+
    labs(title = sprintf("Iteration: %02d/%d", i, n) )


  ggsave(
    sprintf("frames/frame_%04d.png", i),
    plot = p,
    width = 6, height = 5, dpi = 100
  )

}, mc.cores = 20)
# ffmpeg -start_number 1 -i frames/frame_%04d.png -c:v libx264 out.mp4
system("ffmpeg -y -framerate 5   -i frames/frame_%04d.png -pix_fmt yuv420p output.mp4")
browseURL("output.mp4")
