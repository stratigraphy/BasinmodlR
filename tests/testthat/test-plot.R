test_that("plot_basin_model() draws and leaves par() unchanged", {
  res <- run_basin_model(condroz)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  op <- graphics::par(no.readonly = TRUE)
  expect_silent(plot_basin_model(res))
  expect_equal(graphics::par(no.readonly = TRUE), op)
})
