test_that("print() summarises the model", {
  res <- run_basin_model(brabant)
  expect_output(print(res), "Brabant")
})

test_that("plot() draws the composite figure", {
  res <- run_basin_model(condroz)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(res))
})
