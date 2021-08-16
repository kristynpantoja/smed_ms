grid12 = list("numx" = numx,
              "x_seq" = x_seq,
              "x_grid" = x_grid,
              "null_cov1d" = null_cov1d,
              "null_mean1d" = null_mean1d,
              "null_cov2d" = null_cov2d,
              "null_mean2d" = null_mean2d)
saveRDS(grid12, paste(home, "/grid_gpvs12.rds", sep = ""))
