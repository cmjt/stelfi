To install the `R` package `stelfi` run
`devtools::install_github("cmjt/stelfi")`, Use
`devtools::install_github("cmjt/stelfi",build_vignettes = TRUE)` if you
want to build vignettes at the same time.

    library(stelfi)

ETAS model for earthquake data
------------------------------

    data(earthquakes)
    times <- earthquakes$origintime
    ##  2012
    times <- subset(times, lubridate::year(times) == "2012")

![Daily rate of earthquakes in Canterbury in
2012](self-exciting_files/figure-markdown_strict/hist_earth-1.png)

    t.idx <- sort(as.numeric(difftime(times, min(times), unit = "days")))
    params <- c(mu = 10,alpha = 13,beta = 20)
    fit <- fit_hawkes(times = t.idx,parameters = params)

    fit

    ## sdreport(.) result
    ##        Estimate Std. Error
    ## mu    0.3104913  0.0616009
    ## alpha 0.6232689  0.2134116
    ## beta  1.0915660  0.4834086
    ## Maximum gradient component: 0.1576034

    par <- fit$par.fixed
    show_hawkes(t.idx, mu = par[1], alpha = par[2],beta = par[3])

![Fitted model for daily occurence of earthquakes in Canterbury,
2012](self-exciting_files/figure-markdown_strict/inference_earth-1.png)

Self-exciting models for terrorism attacks perpetrated by ISIL in Iraq, 2013--2017
----------------------------------------------------------------------------------

    data(terrorism)
    ## make date vector
    times <- lubridate::dmy(paste(terrorism$iday, terrorism$imonth, terrorism$iyear, sep = "-"),tz = "Asia/Baghdad")

![Daily attack rates by ISIL in Iraq
2013--2017](self-exciting_files/figure-markdown_strict/hist_iq-1.png)

![Daily attack rates by ISIL in Iraq
2016](self-exciting_files/figure-markdown_strict/hist2_iq-1.png)
