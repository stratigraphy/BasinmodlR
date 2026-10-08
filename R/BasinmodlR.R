#' @name BasinmodlR
#' @aliases BasinmodlR-package
#' @title One-Dimensional Burial and Thermal History Modelling of Sedimentary
#' Basins
#' @description BasinmodlR is an R package for one-dimensional (1D) burial,
#' compaction and thermal history modelling of sedimentary basins. It is aimed
#' at geologists and stratigraphers who want to reconstruct how deep a
#' stratigraphic unit was buried, how much it was compacted and how hot it
#' became through geological time, starting from nothing more than a
#' stratigraphic column. The whole basin history is described by two small
#' tables: one row per stratigraphic unit and one row per erosion event. From
#' these tables the package derives the deposition phases, hiatuses, layer
#' stacks, the maximum-burial depth of each unit and its grain (solid)
#' thickness. Compaction follows an exponential porosity-depth law and is treated
#' as irreversible, so that uplifted layers keep the porosity they acquired at
#' their deepest burial. Temperatures are the steady-state conductive
#' geotherm of the layered column including radiogenic heat production. The
#' package includes functions to draw a composite figure of the burial path,
#' temperature and hydrocarbon generation windows and to draw a petroleum system
#' event chart. The methods follow the textbook treatment of Allen and Allen
#' (2013).
#'
#' @section Describing a basin:
#' A basin is a list with five elements (see \code{\link{brabant}},
#' \code{\link{condroz}} and \code{\link{ardennes}} for worked examples):
#' \describe{
#'   \item{\code{units}}{one row per stratigraphic unit, oldest first, with the
#'   deposition interval (\code{start} and \code{end} in Ma, \code{start} being
#'   the older age), the thickness at maximum burial (\code{thickness}, km), the
#'   depositional porosity (\code{phi0}), the porosity-depth coefficient
#'   (\code{ck}, 1/km), the thermal conductivity of the grains (\code{K}, W/m/K),
#'   the radiogenic heat production (\code{A}, uW/m3) and optionally the grain
#'   density (\code{rho}, kg/m3).}
#'   \item{\code{erosions}}{one row per erosion event, with \code{start},
#'   \code{end} (Ma) and the amount of section removed (\code{amount}, km).
#'   Erosion is stripped off the top of the column and may cut through several
#'   units.}
#'   \item{\code{water_depth}}{water depth (m) as (age, value) knots.}
#'   \item{\code{heat_flow}}{heat flow into the base of the column (mW/m2) as
#'   (age, value) knots.}
#'   \item{\code{surface_temp}}{temperature at the sediment-water interface
#'   (degrees Celsius) as (age, value) knots.}
#' }
#' Deposition phases follow from the unit intervals, erosion phases from the
#' erosion table, and every remaining gap in the timeline automatically becomes
#' a hiatus. Nothing has to be repeated for each phase: layer stacks, layer
#' identifiers, maximum-burial windows and the amount of section previously
#' removed are all worked out by \code{\link{run_basin_model}}.
#'
#' @section How the model works:
#' Porosity decreases exponentially with depth below the sediment-water
#' interface: \eqn{\phi(z) = \phi_0 e^{-c z}}{phi(z) = phi0 * exp(-c * z)}. For
#' every unit the model first determines the deepest position that the unit ever
#' occupies and from the thickness at that position it computes the volume of
#' grains (the solid thickness), which does not change during burial. At every
#' time step the thickness of a layer whose top is at depth \eqn{z_1} is then
#' found by solving
#' \deqn{h = h_{solid} + (\phi_0/c)(e^{-c z_1} - e^{-c (z_1 + h)}),}{h = h_solid + (phi0/c) * (exp(-c * z1) - exp(-c * (z1 + h))),}
#' which is iterated to convergence. Compaction is irreversible, so the
#' thickness that is used is the smaller of this solution and the thickness the
#' layer had at its deepest burial; an exhumed layer keeps its porosity. During
#' deposition grain material is added at a constant rate, so a unit reaches
#' exactly its maximum-burial thickness at the end of its own deposition phase
#' and hands over to the next phase without a step.
#'
#' Temperatures are the steady-state conductive solution for the layered
#' column. The bulk conductivity of a layer is the porosity-weighted mean of the
#' conductivity of the grains and of the pore water. The surface heat flow is the
#' basal heat flow plus the radiogenic production of the whole column, and the
#' heat flow entering each layer is reduced by the production of the layers
#' above it, so that within a layer of thickness \eqn{h}
#' \deqn{T_{base} = T_{top} + q h / K - A h^2 / (2 K).}{T_base = T_top + q * h / K - A * h^2 / (2 * K).}
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{run_basin_model}}}{runs the model and returns a table with
#'   one row for every layer at every age plus a temperature raster on a common
#'   age-depth grid.}
#'   \item{\code{\link{plot_basin_model}}}{composite figure: geological periods,
#'   event bar, temperature raster with burial paths and isotherms, colour key
#'   and hydrocarbon generation-window panels. Also available as
#'   \code{plot(<basin_model>)}.}
#'   \item{\code{\link{plot_event_chart}}}{petroleum system event chart.}
#' }
#'
#' @section Units:
#' Input thicknesses and erosion amounts are in km, porosity coefficients in
#' 1/km, radiogenic heat production in uW/m3, heat flow in mW/m2 and water depth
#' in m. The results are returned in metres, degrees Celsius, W/m/K and Ma.
#'
#' @references
#' The 'BasinmodlR' package builds upon the textbook treatment of basin
#' analysis by Allen and Allen (2013). The following list of articles and books
#' is relevant for the 'BasinmodlR' R package and its functions. Individual
#' references are also cited in the descriptions of the functions when relevant.
#' The list can be grouped in three subjects: (1) basin analysis, compaction and
#' backstripping, (2) heat flow and thermal modelling and petroleum systems and
#' (3) the example basins. \cr
#'
#' #1. Basin analysis, compaction and backstripping \cr
#'
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
#' #2. Thermal modelling and petroleum systems \cr
#'
#' Hantschel, T., and Kauerauf, A.I. (2009). Fundamentals of Basin and Petroleum
#' Systems Modeling. Springer, Berlin, 476 p.
#' \doi{10.1007/978-3-540-72318-9} \cr
#'
#' Magoon, L.B., and Dow, W.G. (1994). The petroleum system. In: Magoon, L.B.
#' and Dow, W.G. (eds), The Petroleum System - From Source to Trap. American
#' Association of Petroleum Geologists Memoir 60, 3--24.
#' \doi{10.1306/M60585C1} \cr
#'
#' Tissot, B.P., and Welte, D.H. (1984). Petroleum Formation and Occurrence,
#' 2nd ed. Springer, Berlin, 699 p. \doi{10.1007/978-3-642-87813-8} \cr
#'
#' Turcotte, D.L., and Schubert, G. (2014). Geodynamics, 3rd ed. Cambridge
#' University Press, Cambridge. \cr
#'
#' #3. Example basins \cr
#'
#' Arts, M.C.M., Djouder, H., Devleeschouwer, X., and Da Silva, A.-C.
#' (2026). Solid bitumen in the Silurian Bonne Esperance Formation at Huy,
#' Belgium: petroleum system elements reconstructed from integrated maturity
#' indicators and 1D petroleum system modelling. Manuscript submitted to
#' Geologica Belgica. \cr
#'
#' Gradstein, F.M., Ogg, J.G., Schmitz, M.D., and Ogg, G.M. (eds) (2020).
#' Geologic Time Scale 2020. Elsevier, Amsterdam, 1357 p.
#' \doi{10.1016/C2020-1-02369-3} \cr
#'
#' Li, X., Hu, Y., Guo, J., Lan, J., Lin, Q., Bao, X., Yuan, S., Wei, M., Li, Z.,
#' Man, K., Yin, Z., Han, J., Zhang, J., Zhu, C., Zhao, Z., Liu, Y., Yang, J.,
#' and Nie, J. (2022). A high-resolution climate simulation dataset for the past
#' 540 million years. Scientific Data, 9, 371.
#' \doi{10.1038/s41597-022-01490-4}
#'
#' @details Package: 'BasinmodlR'
#'
#' Type: R package
#'
#' Version: 0.1.0 (4th quarter 2026)
#'
#' License: GPL-3
#'
#' @note
#' If you want to use this package for publication or research
#' purposes, please cite:
#'
#' Arts, M.C.M. (2026).
#' BasinmodlR: One-Dimensional Burial and Thermal History Modelling of
#' Sedimentary Basins.
#' https://CRAN.R-project.org/package=BasinmodlR
#'
#' @examples
#' # Burial and thermal history of the Huy section (Condroz inlier, Belgium)
#' res <- run_basin_model(condroz)
#' res
#' head(res$profile)
#'
#' @author Michiel Arts
#'
#' Maintainer: Michiel Arts \email{michiel.arts@stratigraphy.eu}
NULL
