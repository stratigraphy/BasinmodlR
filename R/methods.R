#' @title Print and plot methods for basin model results
#'
#' @description Methods for objects of class \code{"basin_model"}, the result of
#' \code{\link{run_basin_model}}. \code{print} gives a short summary of the
#' modelled burial and temperature history and \code{plot} draws the composite
#' figure of \code{\link{plot_basin_model}}.
#'
#' @param x An object of class \code{"basin_model"}.
#' @param ... For \code{plot}, further arguments passed on to
#' \code{\link{plot_basin_model}}. Not used by \code{print}.
#'
#' @return \code{x}, invisibly.
#'
#' @author Michiel Arts
#'
#' @seealso \code{\link{run_basin_model}}, \code{\link{plot_basin_model}}
#'
#' @examples
#' res <- run_basin_model(brabant)
#' print(res)
#' plot(res)
#' @name basin_model_methods
NULL

#' @rdname basin_model_methods
#' @export
print.basin_model <- function(x, ...) {

  cat("One-dimensional basin model")
  if (!is.na(x$site_name)) cat(": ", x$site_name, sep = "")
  cat("\n")

  cat("  time span            : ", max(x$ages) + 1, " to ", min(x$ages),
      " Ma (", length(x$ages), " time steps of 1 Myr)\n", sep = "")
  cat("  stratigraphic units  : ", nrow(x$units), "\n", sep = "")

  types <- vapply(x$phases, function(p) p$type, "")
  cat("  phases               : ", sum(types == "deposition"), " deposition, ",
      sum(types == "erosion"), " erosion, ", sum(types == "hiatus"),
      " hiatus\n", sep = "")

  i_max <- which.max(x$total_thickness)
  cat("  maximum burial       : ", round(max(x$total_thickness)),
      " m of sediment at ", x$ages[i_max], " Ma\n", sep = "")
  cat("  preserved today      : ", round(x$total_thickness[length(x$ages)]),
      " m of sediment\n", sep = "")
  cat("  maximum temperature  : ", round(max(x$profile$T_base), 1),
      " degrees Celsius\n", sep = "")
  cat("  radiogenic heat      : ", x$radiogenic, "\n", sep = "")

  invisible(x)
}

#' @rdname basin_model_methods
#' @export
plot.basin_model <- function(x, ...) {
  plot_basin_model(x, ...)
}
