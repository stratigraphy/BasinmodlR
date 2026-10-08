#' @title Example basins for the 'BasinmodlR' package
#' @name BasinmodlR_datasets
#' @description Three worked basin descriptions for the Lower Palaeozoic of
#' Belgium, each ready to be handed to \code{\link{run_basin_model}}:\cr
#' \cr
#' The \code{brabant} data set is a highly generalised 1D model of the Brabant
#' Massif (525-0 Ma).\cr
#' \cr
#' The \code{condroz} data set is the 1D model of the Rue Bonne Espérance
#' section at Huy in the Condroz Inlier (520-0 Ma).\cr
#' \cr
#' The \code{ardennes} data set is a highly generalised 1D model of the
#' Cambro-Ordovician inliers of the Ardennes, such as the Stavelot Massif
#' (525-0 Ma).\cr
#' \cr
#'
#' @details
#' All three data sets have the same structure: a list with the elements
#' \describe{
#'   \item{\code{name}}{character, the label that is used by \code{print} and by
#'   the figures.}
#'   \item{\code{units}}{data frame, one row per stratigraphic unit, oldest
#'   first, with the columns \code{name} (character), \code{start} and
#'   \code{end} (deposition interval in Ma, \code{start} being the older age),
#'   \code{thickness} (km at maximum burial, that is as deposited and before any
#'   later erosion of that unit), \code{phi0} (depositional porosity),
#'   \code{ck} (porosity-depth coefficient, 1/km), \code{K} (thermal
#'   conductivity of the grains, W/m/K), \code{A} (radiogenic heat production,
#'   uW/m3) and \code{rho} (grain density, kg/m3, carried through but not used
#'   in the calculation yet).}
#'   \item{\code{erosions}}{data frame with \code{start}, \code{end} (Ma) and
#'   \code{amount} (km of section that is removed from the top of the column).}
#'   \item{\code{water_depth}}{two column matrix with the age (Ma) and the water
#'   depth (m).}
#'   \item{\code{heat_flow}}{two column matrix with the age (Ma) and the heat
#'   flow into the base of the column (mW/m2).}
#'   \item{\code{surface_temp}}{data frame with the age (Ma) and the
#'   temperature at the sediment-water interface (degrees Celsius), at 10 Myr
#'   resolution.}
#'   \item{\code{source_units}}{integer, the rows of \code{units} that get a
#'   generation panel in \code{\link{plot_basin_model}}.}
#' }
#' Ages between the knots of \code{water_depth}, \code{heat_flow} and
#' \code{surface_temp} are filled in by linear interpolation, and the curves are
#' held constant beyond their first and last knot. Gaps in the timeline, such as
#' 480-467 Ma in the Brabant Massif and at Huy, become hiatuses automatically.
#'
#' Stratigraphic units are grouped where possible, so that one modelled unit
#' can contain several formations or groups. Unit names such as "layer 1" are
#' generalised units without a formal lithostratigraphic name. Thermal
#' conductivity, heat production, depositional porosity and the porosity-depth
#' coefficient are reference values for generalised lithologies (Allen and Allen,
#' 2013), and the basal heat flow was assigned to each time interval according
#' to the prevailing tectonic regime. The heat flow and the eroded thicknesses
#' are therefore scenario assumptions and not independently constrained
#' quantities, and the models should be regarded as scenario-based
#' reconstructions (Hantschel and Kauerauf, 2009). The palaeo-surface
#' temperatures were extracted for the palaeo-position of Belgium from the
#' climate simulations of Li et al. (2022) and the same curve is used for all
#' three basins.
#'
#' The maximum-burial depth window of every unit, which has to be listed by
#' hand for every phase in layer-by-layer scripts, follows from the units and
#' the erosions and is not part of the input. For the Brabant Massif this gives
#' 9600 m of section at 420 Ma, 8600 m at 300 Ma and 6100 m preserved today.
#'
#'
#' @format Each data set is a list of seven elements, described under Details.
#'
#' @source
#' The stratigraphic framework (unit thickness, depositional environment and the
#' position of hiatuses) was compiled from the published stratigraphy of the
#' Brabant Massif, the Condroz Inlier and the Ardennes and is described in Arts
#' et al. (2026). The surface temperature curve is derived from Li et al.
#' (2022).
#'
#' @references
#' Allen, P.A., and Allen, J.R. (2013). Basin Analysis: Principles and
#' Applications, 3rd ed. Blackwell Publishing, Oxford, 619 p. \cr
#'
#' Arts, M.C.M., Djouder, H., Devleeschouwer, X., and Da Silva, A.-C.
#' (2026). Solid bitumen in the Silurian Bonne Espérance Formation at Huy,
#' Belgium: petroleum system elements reconstructed from integrated maturity
#' indicators and 1D petroleum system modelling. Manuscript submitted to
#' Geologica Belgica. \cr
#'
#' Hantschel, T., and Kauerauf, A.I. (2009). Fundamentals of Basin and Petroleum
#' Systems Modeling. Springer, Berlin, 476 p.
#' \doi{10.1007/978-3-540-72318-9} \cr
#'
#' Li, X., Hu, Y., Guo, J., Lan, J., Lin, Q., Bao, X., Yuan, S., Wei, M., Li, Z.,
#' Man, K., Yin, Z., Han, J., Zhang, J., Zhu, C., Zhao, Z., Liu, Y., Yang, J.,
#' and Nie, J. (2022). A high-resolution climate simulation dataset for the past
#' 540 million years. Scientific Data, 9, 371.
#' \doi{10.1038/s41597-022-01490-4}
#'
#' @examples
#' # Structure of an example basin
#' str(brabant, max.level = 1)
#'
#' # The stratigraphic units of the Brabant Massif
#' brabant$units
#'
#' # The erosion events at Huy
#' condroz$erosions
#'
#' # The water depth curve of the Ardennes (age in Ma, water depth in m)
#' ardennes$water_depth
#'
#' @keywords datasets
NULL

#' @rdname BasinmodlR_datasets
#' @format NULL
"brabant"

#' @rdname BasinmodlR_datasets
#' @format NULL
"condroz"

#' @rdname BasinmodlR_datasets
#' @format NULL
"ardennes"
