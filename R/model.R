#' @title Run a one-dimensional burial and thermal history model
#'
#' @description Forward models the burial, compaction and temperature history
#' of a one-dimensional sedimentary column (a single "well" or stratigraphic
#' section) through geological time. The basin is described by one row per
#' stratigraphic unit and one row per erosion event. Deposition phases, hiatuses,
#' layer stacks, maximum-burial windows and grain (solid) thicknesses are all
#' derived from that description, so nothing has to be re-specified for each
#' time interval.
#'
#' @details
#' \strong{Time and depth.} Time runs from the oldest age of the oldest unit to
#' the youngest age of the last phase in steps of 1 Myr. Depths are measured
#' from the sediment-water interface, so the water column is not part of the
#' burial depth. A changing water depth therefore does not compact or decompact
#' the column; it only shifts the column relative to sea level.
#'
#' \strong{Phases.} The timeline is divided into deposition phases (one for each
#' unit in \code{units}), erosion phases (one for each row of \code{erosions})
#' and hiatuses. Every part of the timeline that is covered by neither a unit nor
#' an erosion event is automatically treated as a hiatus (no deposition, no
#' erosion). Phases may not overlap.
#'
#' \strong{Compaction.} Porosity decreases exponentially with burial depth
#' (Sclater and Christie, 1980; Allen and Allen, 2013):
#' \deqn{\phi(z) = \phi_0 e^{-c z}}{phi(z) = phi0 * exp(-c * z)}
#' where \eqn{\phi_0} is the depositional porosity (\code{phi0}), \eqn{c} is the
#' porosity-depth coefficient (\code{ck}) and \eqn{z} the burial depth. The volume
#' of solid grains in a layer does not change during burial, so the present
#' thickness of a layer is found by "decompaction" (the reverse of backstripping,
#' Steckler and Watts, 1978; Allen and Allen, 2013). The decompaction equation is
#' solved by fixed-point iteration for every layer at every time step, rather than
#' on a search grid, so layer boundaries are not quantised.
#'
#' \strong{Irreversible compaction.} Compaction is treated as irreversible: the
#' thickness of a layer is never larger than the thickness it had at its deepest
#' burial. An uplifted (exhumed) layer therefore keeps the porosity it acquired
#' at maximum burial instead of re-expanding.
#'
#' \strong{Deposition.} During deposition grain material is added at a constant
#' rate. A unit therefore reaches exactly its maximum-burial thickness at the end
#' of its own deposition phase and hands over to the next phase without a step.
#'
#' \strong{Erosion.} Erosion removes the stated amount (km) from the top of the
#' column. It may cut through several units, and the layer that is partly eroded
#' keeps the porosity of its remaining part.
#'
#' \strong{Temperature.} Temperatures are the steady-state conductive solution
#' for a layered column with radiogenic heat production (Allen and Allen, 2013;
#' Turcotte and Schubert, 2014). The bulk thermal conductivity of a layer is the
#' porosity-weighted arithmetic mean of the conductivity of the grains
#' (\code{K}) and of the pore water (\code{k_water}). Within a layer of thickness
#' \eqn{h}, the temperature at depth \eqn{d} below the top of the layer is
#' \deqn{T(d) = T_{top} + q_{top} d / K - A d^2 / (2 K)}{T(d) = T_top + q_top * d / K - A * d^2 / (2 * K)}
#' where \eqn{q_{top}} is the heat flow entering the top of the layer and
#' \eqn{A} is its radiogenic heat production.
#' With \code{radiogenic = "cumulative"} the surface heat flow is the basal heat
#' flow plus the radiogenic production of the whole column, and the heat flow
#' entering each layer is reduced by the production of the layers above it. With
#' \code{radiogenic = "legacy"} every layer instead receives the basal heat flow
#' plus its own production over the full column thickness; use this option only
#' to reproduce results that were calibrated against that older scheme.
#'
#' \strong{Time-dependent boundary conditions.} \code{water_depth},
#' \code{heat_flow} and \code{surface_temp} can be given as a table of
#' (age, value) knots. Values between the knots are linearly interpolated and
#' the curves are held constant beyond the first and last knot. A single number
#' is used as a constant for the whole run. A path to a csv file with the age in
#' the first and the value in the second column is also accepted.
#'
#' @param site A basin description: a list with the elements \code{units},
#' \code{erosions}, \code{water_depth}, \code{heat_flow} and \code{surface_temp},
#' and optionally \code{name} (used by the print method and the figures) and
#' \code{source_units} (used by \code{\link{plot_basin_model}}). See
#' \code{\link{brabant}}, \code{\link{condroz}} and \code{\link{ardennes}} for
#' worked examples. Can be left \code{NULL} when \code{units},
#' \code{erosions}, \code{water_depth}, \code{heat_flow} and
#' \code{surface_temp} are supplied separately.
#' @param units Data frame of stratigraphic units, one row per unit, oldest
#' first. It needs the following columns: \code{name} (character label),
#' \code{start} and \code{end} (deposition interval in Ma, \code{start} being the
#' older age), \code{thickness} (km, the thickness at maximum burial, that is
#' as deposited and before any later erosion of that unit), \code{phi0}
#' (depositional porosity, between 0 and 1), \code{ck} (porosity-depth
#' coefficient in 1/km), \code{K} (thermal conductivity of the grains in W/m/K)
#' and \code{A} (radiogenic heat production in uW/m3). The optional column
#' \code{rho} (grain density in kg/m3) is carried through but is not used in the
#' calculation yet. Typical values of \code{phi0} and \code{ck} for common
#' lithologies are given by Sclater and Christie (1980) and Allen and Allen
#' (2013). Defaults to \code{site$units}.
#' @param erosions Data frame of erosion events with the columns \code{start}
#' and \code{end} (Ma, \code{start} being the older age) and \code{amount} (km of
#' section removed from the top of the column), or \code{NULL} when no erosion
#' occurred. Defaults to \code{site$erosions}.
#' @param water_depth Water depth in metres. Either a two column matrix or data
#' frame with age (Ma) and water depth (m), the path of a csv file in the same
#' layout, or a single number for a constant water depth. Defaults to
#' \code{site$water_depth}.
#' @param heat_flow Heat flow into the base of the column in mW/m2, given in the
#' same formats as \code{water_depth}. Defaults to \code{site$heat_flow}.
#' @param surface_temp Temperature at the sediment-water interface in degrees
#' Celsius, given in the same formats as \code{water_depth}. Defaults to
#' \code{site$surface_temp}.
#' @param dz Depth step (m) of the temperature raster that is returned in
#' \code{temp_grid}. Default is 10.
#' @param k_water Thermal conductivity of the pore water in W/m/K. Default is
#' 0.6.
#' @param radiogenic How radiogenic heat is carried through the column, either
#' \code{"cumulative"} (default) or \code{"legacy"}. See Details.
#' @param verbose Print a one line summary of the run. Default is
#' \code{FALSE}.
#' @param genplot Plot the result with \code{\link{plot_basin_model}}. Default is
#' \code{FALSE}.
#'
#' @return
#' An object of class \code{"basin_model"}. This is a list with the elements:
#' \describe{
#'   \item{\code{units}}{the unit table, sorted from old to young.}
#'   \item{\code{phases}}{list of all deposition, erosion and hiatus phases.
#'   Every phase has a \code{type}, a \code{start} and \code{end} age (Ma), the
#'   \code{unit} (row of \code{units}, \code{NA} when not a deposition phase)
#'   and the eroded \code{amount} (km).}
#'   \item{\code{ages}}{the ages (Ma) of the time steps, from old to young.}
#'   \item{\code{water_depth}}{water depth (m) at each age.}
#'   \item{\code{heat_flow}}{basal heat flow (mW/m2) at each age.}
#'   \item{\code{surface_temp}}{surface temperature (degrees Celsius) at each
#'   age.}
#'   \item{\code{total_thickness}}{thickness (m) of the preserved column at
#'   each age.}
#'   \item{\code{profile}}{a data frame with one row for every layer at every
#'   age: \code{age}, \code{unit}, \code{name}, \code{z_top}, \code{z_base} and
#'   \code{z_mid} (burial depth in m below the sediment-water interface),
#'   \code{h} (thickness in m), \code{porosity}, \code{K} (bulk thermal
#'   conductivity in W/m/K), \code{A} (heat production in W/m3), \code{q_top}
#'   (heat flow into the top of the layer in W/m2), \code{T_top}, \code{T_base}
#'   and \code{T_mid} (degrees Celsius) and \code{water_depth} (m).}
#'   \item{\code{depth_grid}}{depths (m below sea level) of the temperature
#'   raster.}
#'   \item{\code{temp_grid}}{temperature raster (degrees Celsius) with one row per
#'   age and one column per element of \code{depth_grid}.}
#'   \item{\code{radiogenic}}{the radiogenic heat scheme that was used.}
#'   \item{\code{site_name}}{the name of the basin.}
#'   \item{\code{source_units}}{the \code{source_units} of \code{site}, if any.}
#' }
#'
#' @author
#' The porosity-depth law, the decompaction approach and the heat conduction
#' equation follow the textbook treatment of Allen and Allen (2013).
#'
#' @references
#' Allen, P.A., and Allen, J.R. (2013). Basin Analysis: Principles and
#' Applications, 3rd ed. Blackwell Publishing, Oxford, 619 p. \cr
#'
#' Sclater, J.G., and Christie, P.A.F. (1980). Continental stretching: an
#' explanation of the post-mid-Cretaceous subsidence of the central North Sea
#' Basin. Journal of Geophysical Research, 85, 3711--3739.
#' \doi{10.1029/JB085iB07p03711} \cr
#'
#' Steckler, M.S., and Watts, A.B. (1978). Subsidence of the Atlantic-type
#' continental margin off New York. Earth and Planetary Science Letters, 41,
#' 1--13. \doi{10.1016/0012-821X(78)90036-5} \cr
#'
#' Turcotte, D.L., and Schubert, G. (2014). Geodynamics, 3rd ed. Cambridge
#' University Press, Cambridge.
#'
#' @seealso \code{\link{plot_basin_model}}, \code{\link{plot_event_chart}},
#' \code{\link{BasinmodlR_datasets}}
#'
#' @examples
#' # Burial and thermal history of the Huy section (Condroz inlier, Belgium)
#' # using the example basin that is shipped with the package.
#' res <- run_basin_model(condroz)
#' res
#'
#' # The layer table has one row for every layer at every age.
#' # Show the burial depth and temperature of the Mousty equivalent at 300 Ma.
#' prof <- res$profile
#' prof[prof$age == 300 & prof$unit == 2, c("age", "name", "z_top", "z_base", "T_mid")]
#'
#' # Build your own basin: three units, one erosion event and constant boundary
#' # conditions. Thicknesses and erosion amounts are in km, ages in Ma.
#' my_basin <- list(
#'   name = "Synthetic foreland basin",
#'   units = data.frame(
#'     name = c("Sandstone", "Shale", "Shaly sand"),
#'     start = c(300, 250, 200),   # older age of the deposition interval (Ma)
#'     end = c(250, 200, 150),     # younger age of the deposition interval (Ma)
#'     thickness = c(1.2, 0.8, 0.5),
#'     phi0 = c(0.49, 0.63, 0.56),  # depositional porosity
#'     ck = c(0.27, 0.51, 0.39),    # porosity-depth coefficient (1/km)
#'     K = c(3.0, 1.7, 2.4),        # grain thermal conductivity (W/m/K)
#'     A = c(1.0, 1.8, 1.4)         # radiogenic heat production (uW/m3)
#'   ),
#'   # 0.4 km of section is removed between 100 and 0 Ma. The interval between
#'   # 150 and 100 Ma is not covered by a unit or an erosion event and
#'   # becomes a hiatus automatically.
#'   erosions = data.frame(start = 100, end = 0, amount = 0.4),
#'   water_depth = 0,     # constant water depth (m)
#'   heat_flow = 60,      # constant basal heat flow (mW/m2)
#'   surface_temp = 15    # constant surface temperature (degrees Celsius)
#' )
#' my_res <- run_basin_model(my_basin)
#' max(my_res$profile$T_base)   # maximum temperature reached by any layer
#' @export
#' @importFrom stats approx
#' @importFrom utils read.csv

run_basin_model <- function(site = NULL,
                            units = site$units,
                            erosions = site$erosions,
                            water_depth = site$water_depth,
                            heat_flow = site$heat_flow,
                            surface_temp = site$surface_temp,
                            dz = 10,
                            k_water = 0.6,
                            radiogenic = c("cumulative", "legacy"),
                            verbose = FALSE,
                            genplot = FALSE) {

  radiogenic <- match.arg(radiogenic)

  # 1. Check the input

  if (is.null(units)) {
    stop("no stratigraphic units supplied (use 'site' or 'units')",
         call. = FALSE)
  }

  need <- c("name", "start", "end", "thickness", "phi0", "ck", "K", "A")
  miss <- setdiff(need, names(units))
  if (length(miss)) {
    stop("the unit table is missing the column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (any(units$start <= units$end)) {
    stop("every unit needs start larger than end", call. = FALSE)
  }
  if (any(units$thickness <= 0)) {
    stop("thickness must be positive", call. = FALSE)
  }
  if (any(units$phi0 <= 0 | units$phi0 >= 1)) {
    stop("phi0 must lie between 0 and 1", call. = FALSE)
  }
  if (any(units$ck <= 0)) {
    stop("ck must be positive", call. = FALSE)
  }
  if (any(units$K <= 0)) {
    stop("K must be positive", call. = FALSE)
  }

  # Oldest unit first.
  units <- units[order(-units$start), , drop = FALSE]
  rownames(units) <- NULL
  n <- nrow(units)

  # Convert to SI-based working units: thickness in m, porosity coefficient in
  # 1/m and heat production in W/m3.
  thk  <- units$thickness * 1000
  phi0 <- units$phi0
  cc   <- units$ck / 1000
  Kgr  <- units$K
  Arad <- units$A * 1e-6

  # 2. Build the list of phases (deposition, erosion and hiatus)

  ph <- list()
  for (i in seq_len(nrow(units))) {
    ph[[length(ph) + 1L]] <- list(type = "deposition", start = units$start[i],
                                  end = units$end[i], unit = i, amount = 0)
  }
  if (!is.null(erosions) && nrow(erosions) > 0) {
    for (i in seq_len(nrow(erosions))) {
      ph[[length(ph) + 1L]] <- list(type = "erosion", start = erosions$start[i],
                                    end = erosions$end[i], unit = NA_integer_,
                                    amount = erosions$amount[i])
    }
  }
  # Sort all phases from old to young.
  ph <- ph[order(-vapply(ph, function(p) p$start, 0))]

  # Walk through the sorted phases and fill every gap with a hiatus.
  phases <- list()
  for (k in seq_along(ph)) {
    p <- ph[[k]]
    if (p$start <= p$end) {
      stop("phase start must be the older age", call. = FALSE)
    }
    if (k > 1) {
      prev <- phases[[length(phases)]]$end
      if (p$start > prev) {
        stop("phases overlap", call. = FALSE)
      }
      if (p$start < prev) {
        phases[[length(phases) + 1L]] <- list(type = "hiatus", start = prev,
                                              end = p$start, unit = NA_integer_,
                                              amount = 0)
      }
    }
    phases[[length(phases) + 1L]] <- p
  }

  # The model starts one Myr after the oldest age and ends at the youngest age.
  t_base <- phases[[1]]$start
  t_top  <- phases[[length(phases)]]$end
  ages   <- seq(t_base - 1, t_top)
  nA     <- length(ages)

  # Position of an age in the arrays that hold one value per time step.
  idx_of <- function(a) t_base - a

  # 3. Boundary conditions on the model time grid

  if (is.null(water_depth)) {
    stop("water_depth is missing", call. = FALSE)
  }
  if (is.null(heat_flow)) {
    stop("heat_flow is missing", call. = FALSE)
  }
  if (is.null(surface_temp)) {
    stop("surface_temp is missing", call. = FALSE)
  }

  wd    <- curve_on_ages(water_depth, ages, scale = 1, fallback = 20,
                         label = "water_depth")
  qbase <- curve_on_ages(heat_flow, ages, scale = 1e-3, fallback = 20e-3,
                         label = "heat_flow")   # mW/m2 -> W/m2
  tsurf <- curve_on_ages(surface_temp, ages, scale = 1, fallback = 20,
                         label = "surface_temp")

  # A negative water depth (land) is treated as zero.
  wd[wd < 0] <- 0

  # 4. First pass: find the deepest position of every unit
  # Compaction is irreversible, so the state of a unit at its maximum burial
  # fixes its grain (solid) thickness. In this pass the column is stacked
  # without compaction (every unit keeps what is left of its thickness at
  # maximum burial) to find the window of burial depths [zt0, zb0] of each unit
  # at the time that it is buried the deepest.

  rem <- rep(0, n)            # remaining thickness of each unit (m)
  stack <- integer(0)         # units in the column, youngest (top) first
  zt0 <- rep(NA_real_, n)     # top of the unit at maximum burial
  zb0 <- rep(NA_real_, n)     # base of the unit at maximum burial

  for (p in phases) {
    if (p$type == "deposition") {
      stack <- c(p$unit, stack)
      rem[p$unit] <- thk[p$unit]
    }

    if (p$type == "erosion") {
      # Strip the eroded amount from the top, layer after layer.
      left <- p$amount * 1000

      while (left > 1e-9 && length(stack) > 0) {
        u <- stack[1]
        take <- min(left, rem[u])
        rem[u] <- rem[u] - take
        left <- left - take

        if (rem[u] <= 1e-9) {
          stack <- stack[-1]
        }
      }
    }

    # Record the deepest position that each unit has reached so far.
    z <- 0
    for (u in stack) {
      z1 <- z
      z <- z + rem[u]

      if (is.na(zb0[u]) || z > zb0[u]) {
        zb0[u] <- z
        zt0[u] <- z1
      }
    }
  }

  # Volume of pore space between depths z1 and z2 for one square metre of
  # column: the integral of phi0 * exp(-c * z) from z1 to z2.
  pore_volume <- function(phi0, cc, z1, z2) {
    (phi0 / cc) * (exp(-cc * z1) - exp(-cc * z2))
  }

  # Grain (solid) thickness of every unit: the thickness at maximum burial
  # minus the pore space. This stays the same during the whole burial history.
  solid_full <- (zb0 - zt0) - pore_volume(phi0, cc, zt0, zb0)

  if (any(is.na(solid_full)) || any(solid_full <= 0)) {
    stop("zero or negative grain thickness detected", call. = FALSE)
  }

  # 5. Second pass: burial, compaction and temperature at every time step

  hmin   <- rep(Inf, n)       # smallest thickness reached so far (irreversible compaction)
  zt     <- rep(NA_real_, n)  # top of the unit at its deepest burial so far
  zb     <- rep(NA_real_, n)  # base of the unit at its deepest burial so far
  solid  <- rep(0, n)         # grain thickness present in the column
  active <- rep(FALSE, n)     # has the unit started to be deposited?
  stack  <- integer(0)

  layers <- vector("list", nA)
  tot    <- numeric(nA)

  for (p in phases) {
    steps <- seq(p$start - 1, p$end)
    N <- length(steps)

    for (j in seq_len(N)) {
      age <- steps[j]
      k <- idx_of(age)

      if (p$type == "deposition") {
        u <- p$unit
        if (!active[u]) {
          active[u] <- TRUE
          stack <- c(u, stack)
        }
        # Constant grain supply: the unit is complete at the end of its phase.
        solid[u] <- (j / N) * solid_full[u]
        hmin[u] <- Inf

      } else if (p$type == "erosion") {
        # Spread the eroded amount evenly over the time steps of the phase.
        left <- (p$amount * 1000) / N

        while (left > 1e-9 && length(stack) > 0) {
          u <- stack[1]
          if (!is.finite(hmin[u])) {
            break
          }
          take <- min(left, hmin[u])
          zt[u] <- zt[u] + take
          hmin[u] <- zb[u] - zt[u]
          solid[u] <- hmin[u] - pore_volume(phi0[u], cc[u], zt[u], zb[u])
          left <- left - take

          if (hmin[u] <= 1e-6) {
            active[u] <- FALSE
            stack <- stack[-1]
          }
        }
      }

      ns <- length(stack)
      if (ns == 0) {
        tot[k] <- 0
        next
      }

      # Decompaction: thickness, depth and porosity of every layer
      z_top <- numeric(ns)
      z_bas <- numeric(ns)
      hh <- numeric(ns)
      por <- numeric(ns)
      z <- 0

      for (i in seq_len(ns)) {
        u <- stack[i]
        a <- (phi0[u] / cc[u]) * exp(-cc[u] * z)
        h <- solid[u]

        # Solve h = solid + a * (1 - exp(-cc * h)) by fixed-point iteration.
        for (iter in seq_len(500)) {
          hn <- solid[u] + a * (1 - exp(-cc[u] * h))
          if (abs(hn - h) < 1e-10) {
            h <- hn
            break
          }
          h <- hn
        }

        # A layer cannot become thicker than it was at its deepest burial.
        hc <- h
        if (hc < hmin[u]) {
          hmin[u] <- hc
          zt[u] <- z
          zb[u] <- z + hc
        }

        hh[i] <- hmin[u]
        z_top[i] <- z
        z <- z + hh[i]
        z_bas[i] <- z

        # Mean porosity of the layer, using the depth window of its deepest
        # burial.
        por[i] <- if (hh[i] > 0) {
          pore_volume(phi0[u], cc[u], zt[u], zb[u]) / hh[i]
        } else {
          0
        }
      }

      tot[k] <- z

      # Temperature: steady-state conduction through the layer stack
      us <- stack
      Kb <- (1 - por) * Kgr[us] + por * k_water   # bulk conductivity
      Ai <- Arad[us]

      Tt <- numeric(ns)
      Tb <- numeric(ns)
      Tm <- numeric(ns)
      qt <- numeric(ns)

      # Heat flow at the surface: basal heat flow plus the heat produced inside
      # the column.
      qs <- qbase[k] + sum(Ai * hh)
      Tp <- tsurf[k]

      for (i in seq_len(ns)) {
        qi <- if (radiogenic == "cumulative") qs else qbase[k] + Ai[i] * z
        qt[i] <- qi
        Tt[i] <- Tp
        Tb[i] <- Tp + qi * hh[i] / Kb[i] - Ai[i] * hh[i]^2 / (2 * Kb[i])
        Tm[i] <- Tp + qi * (hh[i] / 2) / Kb[i] - Ai[i] * (hh[i] / 2)^2 / (2 * Kb[i])

        if (radiogenic == "cumulative") {
          # The heat produced in this layer no longer flows through the layers
          # below it.
          qs <- qs - Ai[i] * hh[i]
        }
        Tp <- Tb[i]
      }

      layers[[k]] <- data.frame(
        age = age,
        unit = us,
        name = units$name[us],
        z_top = z_top,
        z_base = z_bas,
        z_mid = z_top + hh / 2,
        h = hh,
        porosity = por,
        K = Kb,
        A = Ai,
        q_top = qt,
        T_top = Tt,
        T_base = Tb,
        T_mid = Tm,
        water_depth = wd[k],
        stringsAsFactors = FALSE
      )
    }
  }

  profile <- do.call(rbind, layers[!vapply(layers, is.null, TRUE)])

  # 6. Temperature raster on a common age-depth grid
  # The raster is used by plot_basin_model(). Depths are measured below sea
  # level, so the water column is added on top of the sediment column.

  zmax <- max(tot + wd)
  depth_grid <- seq(0, ceiling(zmax / dz) * dz, by = dz)
  temp_grid <- matrix(NA_real_, nrow = nA, ncol = length(depth_grid))

  for (k in seq_len(nA)) {
    L <- layers[[k]]
    v <- rep(NA_real_, length(depth_grid))
    # Inside the water column the temperature equals the surface temperature.
    v[depth_grid <= wd[k]] <- tsurf[k]

    if (!is.null(L)) {
      zz <- depth_grid - wd[k]
      for (i in seq_len(nrow(L))) {
        sel <- zz > L$z_top[i] & zz <= L$z_base[i]
        if (any(sel)) {
          d <- zz[sel] - L$z_top[i]
          v[sel] <- L$T_top[i] + L$q_top[i] * d / L$K[i] - L$A[i] * d^2 / (2 * L$K[i])
        }
      }
    }
    temp_grid[k, ] <- v
  }

  # 7. Collect the output

  res <- list(
    units = units,
    phases = phases,
    ages = ages,
    water_depth = wd,
    heat_flow = qbase * 1000,
    surface_temp = tsurf,
    total_thickness = tot,
    profile = profile,
    depth_grid = depth_grid,
    temp_grid = temp_grid,
    radiogenic = radiogenic,
    site_name = if (is.null(site$name)) NA_character_ else site$name,
    source_units = site$source_units
  )

  class(res) <- "basin_model"

  if (verbose) {
    cat("1-D basin model",
        if (!is.na(res$site_name)) paste0(": ", res$site_name) else "",
        "\n", sep = "")
  }

  if (genplot) {
    plot_basin_model(res)
  }

  invisible(res)
}


# Internal helper: bring a boundary condition curve onto the model time grid.
#
# 'x' is a two column matrix or data frame (age in Ma, value), the path of a
# csv file with the same two columns, or a single number. The values are
# multiplied by 'scale' (for unit conversion) and linearly interpolated onto
# 'ages'. Ages outside the range of the knots get the value of the nearest
# knot. When a csv path does not exist the 'fallback' is used with a warning.

curve_on_ages <- function(x, ages, scale = 1, fallback, label) {

  if (is.character(x)) {
    if (!file.exists(x)) {
      warning("file for '", label, "' not found, using a constant value of ",
              fallback / scale, call. = FALSE)
      return(rep(fallback, length(ages)))
    }
    x <- as.matrix(utils::read.csv(x))
  } else if (is.numeric(x) && is.null(dim(x)) && length(x) == 1L) {
    return(rep(x * scale, length(ages)))
  } else {
    x <- as.matrix(x)
  }

  # Only the first two columns (age and value) are used.
  if (ncol(x) > 2) {
    x <- x[, 1:2, drop = FALSE]
  }
  x <- x[order(x[, 1]), , drop = FALSE]

  stats::approx(x[, 1], x[, 2] * scale, xout = ages, rule = 2)$y
}
