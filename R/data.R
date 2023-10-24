#' Daily chemical concentrations of 26 PM2.5 species from Queens, NY (2001-2021).
#'
#' A dataset containing the daily chemical concentrations (in ug/m^3) of 26 PM2.5 species from the Queens,
#' NY air monitor in the EPA's AQS data mart (AQS Site ID: 36-081-0124). The dataset spans 2001-2021.
#'
#' @format A tibble with 2443 rows and 27 variables:
#' \describe{
#'   \item{Date}{The date of the PM2.5 measurement}
#'   \item{...}{The remaining 26 variables are the 26 PM2.5 species: Al, NH4, As, Ba, Br, Cd, Ca, Cl, Cr, Cu, EC, Fe, Pb, Mg, Mn, Ni, OC, K, Se, Si, Na, S, Ti, NO3, V, Zn}
#' }
#' @source \url{https://epa.maps.arcgis.com/apps/webappviewer/index.html?id=5f239fd3e72f424f98ef3d5def547eb5}
"queens"
