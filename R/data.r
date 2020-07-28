#' A dataset of recorded serial killers in the UK
#' constructed from http://www.murderuk.com/
#' @name serial_uk
#' @format A dataframe with 62 rows and 8 variables
#' \describe{
#' \item{Number.of.kills}{approximate number of murders committed by each person}
#' \item{Years}{The years of operation}
#' \item{Name}{Name of convicted serial killer}
#' \item{AKA}{Some serial killers were given nicknames}
#' \item{Year.Start}{Year the murders began}
#' \item{Year.End}{Year the murders ended}
#' \item{Date.of.first.kill}{Date, if known, first murder was committed}
#' \item{population..M.}{Estimated population of UK in millions at time of first murder}
#' }
#' @docType data
#' @usage data(serial_uk)
NULL
#' A dataset of recorded murder cases in New Zealand 2004--2019
#' scraped from the online database https://interactives.stuff.co.nz/2019/the-homicide-report/
#' @name murders_nz
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
#' @docType data
#' @usage data(murders_nz)
NULL
#' A dataset of retweets of NIWA's "viral" leopard seal tweet https://twitter.com/niwa_nz/status/1092610541401587712
#' on the 5th Feb 2019
#' @name retweets_niwa
#' @format A vector of length 4890
#' \describe{
#' \item{}{date and time of retweet}
#' }
NULL
#' A dataset taken from the GeoNet Quake (http://quakesearch.geonet.org.nz/)
#' containing earthquake information in Canterbury, New Zealand 16-Jan-2010--24-Dec-2014
#' @name earthquakes
#' @format A data frame containing 3824 rows and 5 variables:
#' \describe{
#' \item{origintime}{The UTC time of the event's occurrence}
#' \item{longitude}{longitude location}
#' \item{latitude}{latitude location}
#' \item{magnitude}{The magnitude of the earthquake}
#' \item{depth}{The focal depth of the event (km)}
#' }
#' @docType data
#' @usage data(earthquakes)
NULL
#' A dataset taken from the global terrorism database (GTD) (http://www.start.umd.edu/gtd/)
#' containing information of terrorism activity 2010--2017
#' @name terrorism
#' @format A data frame with 72366 rows and 16 variables:
#' \describe{
#' \item{iyear}{numeric year 2010--2017}
#' \item{imonth}{numeric month index 1--12}
#' \item{iday}{numeric day 1--31 (zeros are a non-entry)}
#' \item{country}{country}
#' \item{latitude}{latitude location}
#' \item{longitude}{longitude location}
#' \item{sucess}{Logical fatal or not. TRUE = fatal}
#' \item{nkill}{number of fatalities per attack}
#' \item{specificity}{factor which represents the apatial accuracy of the evets: 1 = most accurate, 5 = worst}
#' \item{gname}{character name of attack perpetrators}
#' \item{x.coord}{x coordinate from location projected onto a sphere}
#' \item{y.coord}{y coordinate from location projected onto a sphere}
#' \item{z.coord}{z coordinate from location projected onto a sphere}
#' \item{popdensity}{scaled: number of people per kilometer squared}
#' \item{luminosity}{scaled: luminosity}
#' \item{tt}{scaled: time to nearest city in minutes}
#' }
#' @docType data
#' @usage data(terrorism)
NULL
#' A dataset of IS attacks in Iraq in 2017 taken from the global terrorism database
#' (GTD) (http://www.start.umd.edu/gtd/)
#' @name iraq
#' @format a list with two elements occ (occurances of IS attacks) and covs.mesh
#' (a list of monthly covariates at set locations in over Iraq, only used as an illustration of
#' a model in Python et al. (TODO)). The third element is a spatial polygon of Iraq.
#' @docType data
#' @usage data(iraq)
