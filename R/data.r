#' Serial killers of the UK, 1828 - 2015
#'
#' A dataset containing the names and number of recorded
#' murders committed by some of the infamous UK serial killers
#' from 1828 to 2015.
#'
#' @format A dataframe with 62 rows and 8 variables:
#' \describe{
#' \item{number_of_kills}{approx number of murders committed}
#' \item{years}{The years of operation}
#' \item{name}{Name of convicted serial killer}
#' \item{aka}{Some serial killers were given nicknames}
#' \item{year_start}{Year the murders began}
#' \item{tear_end}{Year the murders ended}
#' \item{date_of_first_kill}{Date, if known, first murder was committed}
#' \item{population_million}{Est popn in millions at time of first murder}
#' }
#' @source \url{https://www.murderuk.com/}
#' @name uk_serial
NULL
#' Murders of NZ, 2004 - 2019
#' 
#' A dataset of recorded murder cases in New Zealand 2004 - 2019.
#'
#' @format A simple features dataframe, of type \code{POINT}, with 967 observations and 11 fields:
#' \describe{
#' \item{sex}{Biological sex of victim}
#' \item{age}{Age of victim (years)}
#' \item{date}{Month and day of murder}
#' \item{year}{Year}
#' \item{cause}{Cause of death}
#' \item{killer}{Killer}
#' \item{name}{Name of victim}
#' \item{full_date}{Date object of observation on single days}
#' \item{month}{Month name of observation}
#' \item{cause_cat}{Cause of death as category}
#' \item{region}{NZ region}
#' }
#' @source Data scraped and cleaned by Charlie Timmings,
#' honours student at the University of Auckland from the website
#' \url{https://interactives.stuff.co.nz/2019/the-homicide-report/}
#' @name nz_murders
NULL
#' Retweets of NIWA's viral leopard seal Tweet
#'
#' A dataset of retweet times of NIWA's viral leopard seal tweet
#' on the 5th Feb 2019 (\url{https://twitter.com/niwa_nz/status/1092610541401587712}).
#'
#' @format A vector of length 4890 specifying the date and time of retweet (UTC)
#' 
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
#' @source \url{https://quakesearch.geonet.org.nz/}
#' @name nz_earthquakes
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
#' \item{utm_x}{x-coord location UTM}
#' \item{utm_y}{y-coord location UTM}
#' \item{success}{logical, was fatal? \code{TRUE} = fatal}
#' \item{nkill}{number of fatalities per attack}
#' \item{specificity}{spatial accuracy of event: 1 = most accurate, 5 = worst}
#' \item{gname}{character name of attack perpetrators (ISIL)}
#' \item{x_coord}{x coordinate from location projected onto a sphere}
#' \item{y_coord}{y coordinate from location projected onto a sphere}
#' \item{z_coord}{z coordinate from location projected onto a sphere}
#' \item{popdensity}{scaled: number of people per kilometer squared}
#' \item{luminosity}{scaled: luminosity}
#' \item{tt}{scaled: time to nearest city in minutes}
#' }
#' @source \url{https://www.start.umd.edu/data-tools/GTD}
#' @name iraq_terrorism
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
#' Sasquatch (bigfoot) sightings in the USA, 2000 - 2005
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
#' \item{temperature_high}{Reported weather measure}
#' \item{temperature_mid}{Reported weather measure}
#' \item{temperature_low}{Reported weather measure}
#' \item{dew_point}{Reported weather measure}
#' \item{humidity}{Reported weather measure}
#' \item{cloud_cover}{Reported weather measure}
#' \item{moon_phase}{Reported measure}
#' \item{precip_intensity}{Reported weather measure}
#' \item{precip_probability}{Reported weather measure}
#' \item{precip_type}{Reported weather measure}
#' \item{pressure}{Reported weather measure}
#' \item{summary}{Text weather summary}
#' \item{uv_index}{Reported weather measure}
#' \item{visibility}{Reported weather measure}
#' \item{wind_bearing}{Reported weather measure}
#' \item{wind_speed}{Reported weather measure}
#' \item{year}{Year}
#' }
#' @source \url{https://data.world/timothyrenner/bfro-sightings-data}
#' @name sasquatch
NULL
#' Example Delaunay triangulation
#'
#' @format A \code{fmesher::fm_mesh_2d()} based on the outline of a horse
#' @name horse_mesh
NULL
#' Example \code{sf} \code{POLYGON}
#'
#' @format A \code{sf} \code{POLYGON} of a horse outline
#' @name horse_sf
NULL
#' Example multivariate Hawkes dataset
#'
#' Two-stream multivariate Hawkes data with \eqn{mu_1 = \mu_2 = 0.2},
#' \eqn{beta_1 = \beta_2 = 0.7}, and \eqn{\boldsymbol{\alpha}} =
#' \code{matrix(c(0.5,0.1,0.1,0.5),byrow = TRUE,nrow = 2)}.
#'
#' @format A data frame with 213 observations and 2 variables:
#' \describe{
#' \item{times}{ordered time stamp of observation}
#' \item{stream}{character giving stream ID (i.e., Stream 1 or
#' Stream 2) of observation}
#' }
#' @name multi_hawkes
NULL
