#' A dataset of recorded serial killers in the UK
#' constructed from http://www.murderuk.com/
#' @name serialUK
#' @format A dataframe with 62 rows and 8 variables
#' \describe{
#' \item {Number.of.kills}{approximate number of murders committed by each person}
#' \item {Years}{The years of "operation"}
#' \item {Name}{Name of convicted serial killer}
#' \item {AKA}{Some serial killers were given nicknames}
#' \item {Year.Start}{Year the murders began}
#' \item {Year.End}{Year the murders ended}
#' \item {Date.of,first.kill}{Date, if known, first murder was committed}
#' \item {population..M.}{Estimated population of UK in millions at time of first murder}
#' }
#' @docType data
#' @usage data(serialUK)
NULL
#' A dataset of recorded murder cases in New Zealand 2004--2019
#' scraped from the online database https://interactives.stuff.co.nz/2019/the-homicide-report/
#' @name murdersNZ
#' @format A dataframe with 1019 rows and 9 variables
#' \describe{
#' \item{Latitude}{approximate latitude of murder}
#' \item{Longitude} approximate longitude of murder}
#' \item{Sex}{biological sex of victim}
#' \item{Age}{age of victim (years)}
#' \item{Date}{month and day of murder}
#' \item{Year}{year}
#' \item{Cause}{cause of death}
#' \item{Killer}{killer}
#' \item{Name}{name of victim}
#' }
#' @docType data
#' @usage data(murdersNZ)
