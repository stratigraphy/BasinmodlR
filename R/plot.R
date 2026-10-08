#' @title Composite burial and temperature figure
#'
#' @description Draws the standard figure of a basin model run: a geological
#' period bar, an event bar (deposition, erosion and hiatus phases), a
#' temperature raster with the water column, the burial path of every
#' stratigraphic unit (horizon lines) and labelled isotherms, a colour key, and
#' one generation-window panel for every selected unit.
#'
#' @details
#' \strong{Burial and temperature panel.} The horizontal axis is the age (Ma)
#' running from the oldest age on the left to the present on the right. The
#' vertical axis is the depth below sea level (m). The coloured raster is the
#' temperature of the column, the light blue area at the top is the water
#' column, and every black line is the base of one stratigraphic unit. Where a
#' line moves up the unit has been uplifted and eroded, and where it moves down
#' it has been buried deeper. The thin lines are isotherms.
#'
#' \strong{Generation windows.} The lower panels show for a selected unit when
#' it is inside a temperature window (for example the oil window) \emph{for
#' the first time}. Because organic matter that has been heated cannot go back
#' to an immature state, generation can only continue when the unit reaches a
#' temperature that is higher than any temperature it has had before. A unit is
#' therefore drawn as "gen" (generating) at an age when its temperature lies
#' inside \code{window} and is equal to the highest temperature that the unit has
#' reached up to that age (Tissot and Welte, 1984; Allen and Allen, 2013). The
#' default window of 95 to 135 degrees Celsius is a nominal choice and should be
#' adjusted to the kerogen type and the maturity model of the study at hand.
#'
#' \strong{Continuity check.} With \code{verbose = TRUE} the function prints
#' a table that compares the change in layer depth and temperature over the one
#' million year step that straddles every phase boundary with the background
#' rate of change on either side of it. Ratios close to 1 show that no artificial
#' step is introduced at the boundary between two phases.
#'
#' The geological periods and their colours follow the International Chronostratigraphic
#' Chart (Gradstein et al., 2020). The colour scale of the temperature raster
#' is a rainbow-like palette (the 'turbo' colour scheme).
#'
#' @param res An object of class \code{"basin_model"}, the result of
#' \code{\link{run_basin_model}}.
#' @param window Numeric vector of length two, the temperature window in degrees
#' Celsius that is used for the generation panels. Default is \code{c(95, 135)}.
#' @param window_units Which units get a generation panel, given as row numbers
#' of \code{res$units} or as unit names. \code{NULL} uses
#' \code{res$source_units} when the basin description provides it and draws no
#' panel otherwise.
#' @param n_levels Number of colour levels in the temperature raster. Default is
#' 250.
#' @param scale How the colours are spread over the temperature range, either
#' \code{"quantile"} (equal areas get equal numbers of colours, the default) or
#' \code{"linear"} (equal temperature steps get equal colours).
#' @param contour_n Approximate number of isotherms. Default is 20.
#' @param xmax Oldest age shown in Ma. \code{NULL} (default) uses the start of the
#' model.
#' @param tick Spacing of the age axis labels in Myr. Default is 25.
#' @param new_device Open a new graphics window before plotting. Leave this
#' \code{FALSE} (default) when writing to a file device or to the plot pane.
#' @param use The temperature that is used for the generation window test, one of
#' \code{"T_mid"} (the middle of the layer, the default), \code{"T_top"} or
#' \code{"T_base"}.
#' @param span Half width, in Myr, of the interval on either side of a phase
#' boundary that is used to measure the background rate for the continuity
#' check. Default is 5.
#' @param verbose Print the continuity check table. Default is \code{FALSE}.
#'
#' @return \code{res}, invisibly. The function is called for its side effect, the
#' figure.
#'
#' @author
#' Michiel Arts
#'
#' @references
#' Allen, P.A., and Allen, J.R. (2013). Basin Analysis: Principles and
#' Applications, 3rd ed. Blackwell Publishing, Oxford, 619 p. \cr
#'
#' Gradstein, F.M., Ogg, J.G., Schmitz, M.D., and Ogg, G.M. (eds) (2020). Geologic
#' Time Scale 2020. Elsevier, Amsterdam, 1357 p.
#' \doi{10.1016/C2020-1-02369-3} \cr
#'
#' Tissot, B.P., and Welte, D.H. (1984). Petroleum Formation and Occurrence,
#' 2nd ed. Springer, Berlin, 699 p. \doi{10.1007/978-3-642-87813-8}
#'
#' @seealso \code{\link{run_basin_model}}, \code{\link{plot_event_chart}}
#'
#' @examples
#' # Burial and temperature history of the Huy section (Condroz inlier).
#' # The generation panels are drawn for the units listed in condroz$source_units.
#' res <- run_basin_model(condroz)
#' plot_basin_model(res)
#'
#' # Show the generation window of the Mousty equivalent only, using a
#' # narrower window, and print the continuity check.
#' plot_basin_model(res,
#'                  window = c(100, 130),
#'                  window_units = "Mousty eq. (unit 2)",
#'                  verbose = TRUE)
#' @export
#' @importFrom graphics abline arrows axis box lines mtext par rect
#' @importFrom graphics text polygon image contour layout
#' @importFrom grDevices colorRampPalette dev.new
#' @importFrom stats quantile median

plot_basin_model <- function(res,
                             window = c(95, 135),
                             window_units = NULL,
                             n_levels = 250,
                             scale = c("quantile", "linear"),
                             contour_n = 20,
                             xmax = NULL,
                             tick = 25,
                             new_device = FALSE,
                             use = c("T_mid", "T_top", "T_base"),
                             span = 5,
                             verbose = FALSE) {

  scale <- match.arg(scale)
  use <- match.arg(use)

  # Colour palette of the temperature raster (the 'turbo' colour scheme).
  turbo_pal <- grDevices::colorRampPalette(c(
    "#30123B", "#4145AB", "#4675ED", "#39A2FC", "#1BCFD4", "#24ECA6", "#61FC6C",
    "#A4FC3B", "#D1E834", "#F3C63A", "#FE9B2D", "#F36315", "#D93806", "#B11901",
    "#7A0403"
  ))

  # Geological periods (ages in Ma) and their colours for the period bar.
  ics_periods <- data.frame(
    name = c("Quaternary", "Neogene", "Paleogene", "Cretaceous", "Jurassic",
             "Triassic", "Permian", "Carboniferous", "Devonian", "Silurian",
             "Ordovician", "Cambrian"),
    abbr = c(NA, "Neo.", "Pale.", "Cret.", "Jura.", "Tria.", "Perm.", "Carb.",
             "Devo.", "Silu.", "Ordo.", "Camb."),
    base = c(2.58, 23.03, 66, 143.1, 201.4, 251.9, 298.9, 358.9, 419.2, 443.8,
             486.9, 538.8),
    top = c(0, 2.58, 23.03, 66, 143.1, 201.4, 251.9, 298.9, 358.9, 419.2, 443.8,
            486.9),
    col = c("#F9F97F", "#FFE619", "#FD9A52", "#7FC64E", "#34B2C9", "#812B92",
            "#F04028", "#67A599", "#CB8C37", "#B3E1B6", "#009270", "#7FA056"),
    stringsAsFactors = FALSE
  )

  # Continuity check at the phase boundaries
  # Change in base depth (dz) and base temperature (dT) of every unit from one
  # age to the next.
  d_cc <- do.call(rbind, lapply(split(res$profile, res$profile$unit), function(x) {
    x <- x[order(-x$age), ]
    if (nrow(x) < 2) return(NULL)
    data.frame(age = x$age[-1], unit = x$unit[-1],
               dz = abs(diff(x$z_base)), dT = abs(diff(x$T_base)))
  }))

  # Ages of the boundaries between two phases.
  bnd <- vapply(res$phases[-length(res$phases)], function(p) p$end, 0)

  out_cc <- do.call(rbind, lapply(bnd, function(b) {
    at <- d_cc[d_cc$age == b - 1, ]
    nb <- d_cc[abs(d_cc$age - (b - 1)) <= span & d_cc$age != b - 1, ]
    data.frame(boundary_Ma = b,
               jump_z = if (nrow(at)) max(at$dz) else NA_real_,
               bg_z = if (nrow(nb)) stats::median(nb$dz) else NA_real_,
               jump_T = if (nrow(at)) max(at$dT) else NA_real_,
               bg_T = if (nrow(nb)) stats::median(nb$dT) else NA_real_)
  }))

  if (verbose && !is.null(out_cc) && nrow(out_cc) > 0) {
    out_cc$ratio_z <- round(out_cc$jump_z / pmax(out_cc$bg_z, 1e-9), 2)
    out_cc$ratio_T <- round(out_cc$jump_T / pmax(out_cc$bg_T, 1e-9), 2)
    cat("depth and temperature step at each phase boundary versus the background rate\n")
    print(out_cc, row.names = FALSE, digits = 3)
  }

  # Prepare the data for plotting
  # The x axis is flipped so that the oldest age is on the left: x = xmax - age.
  if (is.null(xmax)) xmax <- max(res$ages) + 1

  xf <- function(a) xmax - a
  o <- order(-res$ages)
  ages <- res$ages[o]
  xx <- xf(ages)
  z <- res$temp_grid[o, , drop = FALSE]
  wd <- res$water_depth[o]
  dg <- res$depth_grid
  xtick <- seq(0, floor(xmax / tick) * tick, by = tick)

  if (is.null(window_units)) window_units <- res$source_units
  if (is.character(window_units)) {
    window_units <- match(window_units, res$units$name)
  }
  window_units <- window_units[!is.na(window_units)]
  nw <- length(window_units)

  # Colour breaks: equal numbers of grid cells per colour (quantile) or equal
  # temperature steps per colour (linear).
  if (scale == "quantile") {
    br <- unique(stats::quantile(
      z,
      probs = seq(0, 1, length.out = n_levels + 1),
      na.rm = TRUE
    ))
  } else {
    br <- seq(
      min(z, na.rm = TRUE),
      max(z, na.rm = TRUE),
      length.out = n_levels + 1
    )
  }

  cols <- turbo_pal(length(br) - 1)

  if (new_device) grDevices::dev.new(width = 10, height = 9)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  # Page layout
  # Row 1: period bar, row 2: event bar, row 3: temperature raster and colour
  # key, followed by one row for every generation panel.
  m <- rbind(c(1, 0), c(2, 0), c(3, 4))
  if (nw > 0) {
    for (i in seq_len(nw)) {
      m <- rbind(m, c(4 + i, 0))
    }
  }

  graphics::layout(
    m,
    heights = c(1, 1, 10, rep(1, nw)),
    widths = c(5, 1)
  )

  # Panel 1: geological periods
  graphics::par(mar = c(0.5, 4, 0.5, 0.5))
  plot(
    NA,
    xlim = c(0, xmax),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )

  graphics::mtext("Period", side = 2, line = 1, las = 2, cex = 0.75)

  for (i in seq_len(nrow(ics_periods))) {
    b <- min(ics_periods$base[i], xmax)
    tp <- ics_periods$top[i]

    if (tp >= xmax) next

    graphics::rect(
      xf(b), 0,
      xf(tp), 1,
      col = ics_periods$col[i]
    )

    if (!is.na(ics_periods$abbr[i]) && (b - tp) > xmax / 40) {
      graphics::text(
        xf((b + tp) / 2),
        0.5,
        ics_periods$abbr[i],
        cex = 0.8
      )
    }
  }

  graphics::box()

  plot(
    NA,
    xlim = c(0, xmax),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )

  graphics::mtext("Event", side = 2, line = 1, las = 2, cex = 0.75)

  # Panel 2: deposition, erosion and hiatus phases (numbered)
  ecol <- c(
    deposition = "gold",
    erosion = "grey80",
    hiatus = "white"
  )

  for (i in seq_along(res$phases)) {
    p <- res$phases[[i]]

    graphics::rect(
      xf(p$start), 0,
      xf(p$end), 1,
      col = ecol[[p$type]]
    )

    graphics::text(
      xf((p$start + p$end) / 2),
      0.5,
      i,
      cex = 0.8
    )
  }

  graphics::box()

  # Panel 3: temperature raster, water column, horizons and isotherms
  graphics::par(mar = c(4, 4, 1, 0.5))

  max_plot_depth <- max(c(dg, res$profile$z_base + res$profile$water_depth), na.rm = TRUE)

  graphics::image(
    x = xx,
    y = dg,
    z = z,
    col = cols,
    breaks = br,
    useRaster = TRUE,
    xlim = c(0, xmax),
    ylim = c(max_plot_depth, 0),
    xaxt = "n",
    xlab = "Age (Ma)",
    ylab = "Depth below sea level (m)"
  )

  graphics::axis(
    1,
    at = xf(xtick),
    labels = xtick,
    las = 2
  )

  # Water column.
  graphics::polygon(
    c(xx, rev(xx)),
    c(rep(0, length(xx)), rev(wd)),
    col = "lightblue2",
    border = NA
  )

  graphics::lines(xx, wd, lwd = 2)

  # Base of every unit (burial path), plotted below sea level.
  for (u in sort(unique(res$profile$unit))) {
    d <- res$profile[res$profile$unit == u, ]
    d <- d[order(-d$age), ]

    graphics::lines(
      xf(d$age),
      d$z_base + d$water_depth,
      lwd = 2
    )
  }

  # Isotherms.
  graphics::contour(
    x = xx,
    y = dg,
    z = z,
    add = TRUE,
    drawlabels = TRUE,
    labcex = 0.6,
    method = "flattest",
    col = "grey15",
    levels = pretty(range(z, na.rm = TRUE), contour_n)
  )

  graphics::box()

  # Colour key
  graphics::par(mar = c(4, 1, 1, 4))

  nb <- length(br) - 1

  graphics::image(
    x = 1,
    y = seq_len(nb),
    z = matrix(seq_len(nb), nrow = 1),
    col = cols,
    useRaster = TRUE,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = ""
  )

  at <- round(seq(1, nb, length.out = 6))

  graphics::axis(
    4,
    at = at,
    labels = round(br[at + 1]),
    las = 2
  )

  graphics::mtext(
    "Temperature (\u00B0C)",
    side = 4,
    line = 2.5,
    cex = 0.75
  )

  graphics::box()

  # Generation panels
  # A unit generates when its temperature is inside the window and is the
  # highest temperature it has experienced so far (cummax).
  if (nw > 0) {
    mw_ages <- res$ages
    present <- sort(unique(res$profile$unit))
    mw <- matrix(0, nrow = length(mw_ages), ncol = length(present),
                 dimnames = list(NULL, res$units$name[present]))

    for (j in seq_along(present)) {
      d_mw <- res$profile[res$profile$unit == present[j], ]
      d_mw <- d_mw[order(-d_mw$age), ]
      tt <- d_mw[[use]]
      peak <- cummax(tt)
      hot <- tt >= peak - 1e-9 & tt > window[1] & tt < window[2]
      mw[match(d_mw$age, mw_ages), j] <- as.numeric(hot)
    }

    graphics::par(mar = c(0.5, 4, 0.5, 0.5))

    for (u in window_units) {
      plot(
        NA,
        xlim = c(0, xmax),
        ylim = c(-0.1, 1.2),
        xaxs = "i",
        yaxs = "i",
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n"
      )

      j <- which(present == u)

      if (length(j)) {
        graphics::lines(
          xf(res$ages),
          mw[, j]
        )
      }

      graphics::axis(2, at = c(0, 1), labels = NA)

      graphics::mtext(
        c("no-gen", "gen"),
        side = 2,
        at = c(0, 1),
        line = 0.5,
        las = 2,
        cex = 0.7
      )

      graphics::text(
        xf(xmax / 3),
        0.5,
        res$units$name[u],
        cex = 0.9
      )

      graphics::box()
    }
  }

  invisible(res)
}

