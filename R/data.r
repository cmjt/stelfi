#' Serial killers of the UK
#'
#' A dataset containing the names and number of recorded
#' murders commited by UK serial killers
#'
#' @format A dataframe with 62 rows and 8 variables
#' \describe{
#' \item{Number.of.kills}{approx number of murders committed}
#' \item{Years}{The years of operation}
#' \item{Name}{Name of convicted serial killer}
#' \item{AKA}{Some serial killers were given nicknames}
#' \item{Year.Start}{Year the murders began}
#' \item{Year.End}{Year the murders ended}
#' \item{Date.of.first.kill}{Date, if known, first murder was committed}
#' \item{population..M.}{Est popn in millions at time of first murder}
#' }
#' @source \url{http://www.murderuk.com/}
#' @name serial_uk
NULL
#' Murders of NZ
#' A dataset of recorded murder cases in New Zealand 2004--2019
#'
#' @format A dataframe with 967 rows and 13 variables
#' \describe{
#' \item{Latitude}{approximate latitude of murder}
#' \item{Longitude}{approximate longitude of murder}
#' \item{Sex}{biological sex of victim}
#' \item{Age}{age of victim (years)}
#' \item{Date}{month and day of murder}
#' \item{Year}{year}
#' \item{Cause}{cause of death}
#' \item{Killer}{killer}
#' \item{Name}{name of victim}
#' \item{Full_date}{date object of observation on single days}
#' \item{Month}{month name of observation}
#' \item{Cause_cat}{cause of death as category}
#' \item{Region}{NZ region}
#' }
#' @references \url{https://interactives.stuff.co.nz/2019/the-homicide-report/}
#' @source Data scraped and cleaned by Charlie Timmings,
#' honours student at the University of Auckland.
#' @name murders_nz
NULL
#' NIWA's viral leopard seal Tweet
#'
#' A dataset of retweet times of NIWA's viral leopard seal tweet
#' on the 5th Feb 2019
#'
#' @format A vector of length 4890
#' \describe{
#' \item{}{date and time of retweet}
#' }
#' @source \url{https://twitter.com/niwa_nz/status/1092610541401587712}
#' @name retweets_niwa
NULL
#' Canterbury, NZ earthquakes
#'
#' A dataset taken containing earthquake information in Canterbury,
#' New Zealand 16-Jan-2010--24-Dec-2014
#'
#' @format A data frame containing 3824 rows and 5 variables:
#' \describe{
#' \item{origintime}{The UTC time of the event's occurrence}
#' \item{longitude}{longitude location}
#' \item{latitude}{latitude location}
#' \item{magnitude}{The magnitude of the earthquake}
#' \item{depth}{The focal depth of the event (km)}
#' }
#' @source \url{http://quakesearch.geonet.org.nz/}
#' @name earthquakes
NULL
#' Terrorism in Iraq 2013--2017
#'
#' A dataset containing information of terrorism activity carried out
#' by the Islamic State of Iraq and the Levant (ISIL) in Iraq,
#' 2013--2017
#'
#' @format A data frame with 4208 rows and 16 variables:
#' \describe{
#' \item{iyear}{numeric year 2013--2017}
#' \item{imonth}{numeric month index 1--12}
#' \item{iday}{numeric day 1--31 (zeros are a non-entry)}
#' \item{country}{country (IRAQ)}
#' \item{latitude}{latitude location}
#' \item{longitude}{longitude location}
#' \item{utm.x}{x-coord location UTM}
#' \item{utm.y}{y-coord location UTM}
#' \item{success}{\code{logical} wheather attack was fatal or not, \code{TRUE} = fatal}
#' \item{nkill}{number of fatalities per attack}
#' \item{specificity}{spatial accuracy of event: 1 = most accurate, 5 = worst}
#' \item{gname}{character name of attack perpetrators (ISIL)}
#' \item{x.coord}{x coordinate from location projected onto a sphere}
#' \item{y.coord}{y coordinate from location projected onto a sphere}
#' \item{z.coord}{z coordinate from location projected onto a sphere}
#' \item{popdensity}{scaled: number of people per kilometer squared}
#' \item{luminosity}{scaled: luminosity}
#' \item{tt}{scaled: time to nearest city in minutes}
#' }
#' @source \url{http://www.start.umd.edu/gtd/}
#' @name terrorism
NULL
#' New Zealand polygon
#'
#' Spatial Polygon of NZ
#'
#' @format  A \code{SpatialPolygonsDataFrame} of New Zeland in NZTM
#' @name nz
NULL
#' Example marked point pattern data set
#'
#' @format A data frame with 159 rows and 5 variables:
#' \describe{
#' \item{x}{x coordinate}
#' \item{y}{y coordinate}
#' \item{m1}{mark, Gaussian distributed}
#' \item{m2}{mark, Bernoulli distributed}
#' \item{m3}{mark, Gamma distributed}
#' }
#' @name marked
NULL
#' A self-exciting spatiotemporal point pattern of class \code{stppp}
#'
#' @format A \code{stppp} object with 653 observations
#' \describe{
#' \item{window}{domain of the point pattern of class \code{owin}}
#' \item{n}{number of observations, 653}
#' \item{x}{x coordinate}
#' \item{y}{y coordinate}
#' \item{markformat}{\code{none}}
#' \item{t}{timestamp of points}
#' \item{tlim}{time frame \code{0 2}}
#' }
#' @name xyt
NULL
#' UFO sightings in the USA
#'
#' Subset of data sourced from the UFO Sightings Map
#' (\url{https://www.arcgis.com/apps/webappviewer/index.html?id=ddda71d5211f47e782b12f3f8d06246e}).
#' 
#' @format A data frame with 97027 observations and 12 variables
#' \describe{
#' \item{city}{city where sighting was reported}
#' \item{state}{state code where sighting was reported}
#' \item{date_time}{\code{POSIXct}, time and date of sighting}
#' \item{city_latitude}{latitude of sighting}
#' \item{city_longitude}{longitude of sighting}
#' }
#' @source \url{https://data.world/timothyrenner/ufo-sightings#}
#' @name ufo
NULL
#' Example \code{inla.mesh}
#'
#' @format A \code{inla.mesh} based on the outline of a horse
#' @name horse_mesh
NULL
