test_that("plot_event_chart() returns the event intervals", {
  res <- run_basin_model(condroz)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  ev <- plot_event_chart(res, source_units = c(2, 4), reservoir_units = 5,
                         seal_units = 5,
                         trap_formation = data.frame(start = 320, end = 295),
                         critical_moment = 310)
  expect_type(ev, "list")
  expect_contains(names(ev), c("Source Rock", "Overburden", "Preservation"))
})

test_that("plot_event_chart() leaves par() unchanged", {
  res <- run_basin_model(condroz)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  op <- graphics::par(no.readonly = TRUE)
  plot_event_chart(res, source_units = 2)
  expect_equal(graphics::par(no.readonly = TRUE), op)
})

test_that("plot_event_chart() rejects unknown units", {
  res <- run_basin_model(condroz)
  expect_snapshot(plot_event_chart(res, source_units = "no such unit"),
                  error = TRUE)
})
