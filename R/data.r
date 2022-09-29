#' Serial killers of the UK, 1828 - 2015
#'
#' A dataset containing the names and number of recorded
#' murders commited by some infamous UK serial killers 1828 - 2015.
#'
#' @format A dataframe with 62 rows and 8 variables:
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
#' Murders of NZ, 2004 - 2019
#' 
#' A dataset of recorded murder cases in New Zealand 2004 - 2019.
#'
#' @format A simple features dataframe, of type \code{POINT}, with 967 observations and 11 fields:
#' \describe{
#' \item{Sex}{Biological sex of victim}
#' \item{Age}{Age of victim (years)}
#' \item{Date}{Month and day of murder}
#' \item{Year}{Year}
#' \item{Cause}{Cause of death}
#' \item{Killer}{Killer}
#' \item{Name}{Name of victim}
#' \item{Full_date}{Date object of observation on single days}
#' \item{Month}{Month name of observation}
#' \item{Cause_cat}{Cause of death as category}
#' \item{Region}{NZ region}
#' }
#' @source Data scraped and cleaned by Charlie Timmings,
#' honours student at the University of Aucklandfrom the website
#' \url{https://interactives.stuff.co.nz/2019/the-homicide-report/}
#' @name murders_nz
NULL
#' Retweets of NIWA's viral leopard seal Tweet
#'
#' A dataset of retweet times of NIWA's viral leopard seal tweet
#' on the 5th Feb 2019 (\url{https://twitter.com/niwa_nz/status/1092610541401587712}).
#'
#' @format A vector of length 4890
#' \describe{
#' \item{}{date and time of retweet}
#' }
#' @source \url{https://twitter.com/niwa_nz/status/1092610541401587712}
#' @name retweets_niwa
NULL
#' Earthquakes in Canterbury, NZ, 2010 - 2014
#'
#' Earthquake data from Canterbury,
#' New Zealand 16-Jan-2010--24-Dec-2014.
#'
#' @format A simple features dataframe, of type \code{POINT}, with 3824 observations and 3 fields:
#' \describe{
#' \item{origintime}{The UTC time of the event's occurrence}
#' \item{magnitude}{The magnitude of the earthquake}
#' \item{depth}{The focal depth of the event (km)}
#' }
#' @source \url{http://quakesearch.geonet.org.nz/}
#' @name earthquakes
NULL
#' Terrorism in Iraq, 2013 - 2017
#'
#' A dataset containing information of terrorism activity carried out
#' by the Islamic State of Iraq and the Levant (ISIL) in Iraq,
#' 2013 - 2017.
#'
#' @format  A simple features dataframe, of type \code{POINT}, with 4208 observations and 16 fields:
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
#' Self-exciting point pattern
#'
#' Simulated self-exciting spatiotemporal point pattern of class \code{stppp}
#'
#' @format A \code{stppp} object with 653 observations
#' \describe{
#' \item{window}{Domain of the point pattern of class \code{owin}}
#' \item{n}{Number of observations, 653}
#' \item{x}{x coordinate}
#' \item{y}{y coordinate}
#' \item{markformat}{\code{none}}
#' \item{t}{Timestamp of points}
#' \item{tlim}{Time frame \code{0 2}}
#' }
#' @name xyt
NULL
#' UFO sightings in the USA, 2006 - 2021
#'
#' Subset of data sourced from the UFO Sightings Map
#' (\url{https://www.arcgis.com/apps/webappviewer/index.html?id=ddda71d5211f47e782b12f3f8d06246e}).
#' 
#' @format A simple features dataframe, of type \code{POINT}, with 65603 observations and 3 fields:
#' \describe{
#' \item{city}{city where sighting was reported}
#' \item{state}{state code where sighting was reported}
#' \item{date_time}{\code{POSIXct}, time and date of sighting}
#' }
#' @source \url{https://data.world/timothyrenner/ufo-sightings#}
#' @name ufo
NULL
#' Bigfoot (Sasquatch) sightings in the USA, 2000 - 2005
#'
#' Subset of data sourced from the Bigfoot Field Researchers Organization (BFRO)
#' (\url{https://data.world/timothyrenner/bfro-sightings-data}).
#' 
#' @format A simple features dataframe, of type \code{POINT}, with 972 observations and 27 fields:
#' \describe{
#' \item{observed}{Text observation summary}
#' \item{location_details}{Text location summary}
#' \item{county}{County}
#' \item{state}{State}
#' \item{season}{Season}
#' \item{title}{Report title}
#' \item{date}{Date}
#' \item{number}{Report number}
#' \item{classification}{Report classification}
#' \item{geohash}{Geohash code}
#' \item{temperature_high}{Repotred weather measure}
#' \item{temperature_mid}{Repotred weather measure}
#' \item{temperature_low}{Repotred weather measure}
#' \item{dew_point}{Repotred weather measure}
#' \item{humidity}{Repotred weather measure}
#' \item{cloud_cover}{Repotred weather measure}
#' \item{moon_phase}{Repotred measure}
#' \item{precip_intensity}{Repotred weather measure}
#' \item{precip_probability}{Repotred weather measure}
#' \item{precip_type}{Repotred weather measure}
#' \item{pressure}{Repotred weather measure}
#' \item{summary}{Text weather summary}
#' \item{uv_index}{Repotred weather measure}
#' \item{visibility}{Repotred weather measure}
#' \item{wind_bearing}{Repotred weather measure}
#' \item{wind_speed}{Repotred weather measure}
#' \item{year}{Year}
#' }
#' @source \url{https://data.world/timothyrenner/bfro-sightings-data}
#' @name bigfoot
NULL
#' Example \code{inla.mesh}
#'
#' @format A \code{inla.mesh} based on the outline of a horse
#' @name horse_mesh
NULL
