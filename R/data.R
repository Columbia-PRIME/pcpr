#' Daily chemical concentrations of 26 PM2.5 species from Queens, NY (2001-2021)
#'
#' A dataset containing the chemical concentrations (in µg/m^3) of 26 PM2.5
#' species measured every three to six days from 04/04/2001 through 12/30/2021
#' in Queens, New York City. Data obtained from the U.S. Environmental
#' Protection Agency's Air Quality System data mart (site ID: `36-081-0124`).
#'
#' @format A tibble with 2443 rows and 27 variables:
#' * `Date`: The date the PM2.5 measurements were made
#' * `...`: The remaining 26 variables are the 26 PM2.5 species:
#'   Al, NH4, As, Ba, Br, Cd, Ca, Cl, Cr, Cu, EC, Fe, Pb, Mg, Mn, Ni, OC, K, Se,
#'   Si, Na, S, Ti, NO3, V, Zn
#' @source <https://epa.maps.arcgis.com/apps/webappviewer/index.html?id=5f239fd3e72f424f98ef3d5def547eb5>
#' @references US Environmental Protection Agency. Air Quality System Data Mart
#'   internet database available via
#'   <https://www.epa.gov/outdoor-air-quality-data>. Accessed July 15, 2022.
"queens"
