test_that("example basins run and have the expected structure", {
  res <- run_basin_model(brabant)
  expect_s3_class(res, "basin_model")
  expect_contains(names(res), c("profile", "temp_grid", "ages", "phases"))
  expect_equal(nrow(res$temp_grid), length(res$ages))
  expect_equal(ncol(res$temp_grid), length(res$depth_grid))
})

test_that("Brabant burial reproduces the published thicknesses", {
  res <- run_basin_model(brabant)
  expect_equal(res$total_thickness[res$ages == 420], 9600)
  expect_equal(res$total_thickness[res$ages == 300], 8600)
  expect_equal(res$total_thickness[res$ages == 0], 6100)
})

test_that("hiatuses are derived automatically", {
  res <- run_basin_model(condroz)
  types <- character(length(res$phases))
  for (k in seq_along(res$phases)) {
    types[k] <- res$phases[[k]]$type
  }
  expect_contains(types, "hiatus")
})

test_that("temperature increases with depth in every time step", {
  res <- run_basin_model(condroz)
  prof <- res$profile
  expect_gte(min(prof$T_base - prof$T_top), 0)
})

test_that("compaction is irreversible", {
  res <- run_basin_model(condroz)
  prof <- res$profile
  for (u in unique(prof$unit)) {
    d <- prof[prof$unit == u & prof$age <= res$units$end[u], ]
    d <- d[order(-d$age), ]
    if (nrow(d) > 1) {
      expect_lte(max(diff(d$h)), 1e-6)
    }
  }
})

test_that("constant boundary conditions are accepted", {
  basin <- list(
    units = data.frame(name = "A", start = 100, end = 50, thickness = 1,
                       phi0 = 0.5, ck = 0.3, K = 3, A = 1),
    erosions = NULL,
    water_depth = 0, heat_flow = 60, surface_temp = 15
  )
  res <- run_basin_model(basin)
  expect_equal(max(res$total_thickness), 1000, tolerance = 1e-6)
})

test_that("invalid input gives an error", {
  bad <- brabant
  bad$units$phi0[1] <- 1.2
  expect_snapshot(run_basin_model(bad), error = TRUE)
  expect_snapshot(run_basin_model(), error = TRUE)
})

test_that("legacy and cumulative radiogenic schemes differ", {
  a <- run_basin_model(brabant)
  b <- run_basin_model(brabant, radiogenic = "legacy")
  expect_gt(max(abs(a$profile$T_base - b$profile$T_base)), 0)
})
