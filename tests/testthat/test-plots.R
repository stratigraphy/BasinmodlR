test_that("plot functions run on a file device", {
  res <- run_basin_model(condroz)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot_basin_model(res))
  expect_silent(plot(res))
  ev <- plot_event_chart(res, source_units = c(2, 4), reservoir_units = 5,
                         seal_units = 5,
                         trap_formation = data.frame(start = 320, end = 295),
                         critical_moment = 310)
  expect_type(ev, "list")
  expect_error(plot_event_chart(res, source_units = "no such unit"),
               "unknown unit")
})

test_that("print method works", {
  res <- run_basin_model(brabant)
  expect_output(print(res), "Brabant")
})
