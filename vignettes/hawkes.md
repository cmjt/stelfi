Fitting a Hawkes process to Twitter data
----------------------------------------

A [NIWA](https://niwa.co.nz/) scientist [found a working USB in the scat
of a leopard
seal](https://www.nzherald.co.nz/nz/news/article.cfm?c_id=1&objectid=12201147),
they then [tweeted about
it](https://twitter.com/niwa_nz/status/1092610541401587712) in the hopes
of finding its owner.

Below we use the `stelfi` package to model the times of retweets as a
self-exciting (Hawkes) process.

    library(stelfi)
    data(retweets_niwa)
    head(retweets_niwa)

    ## numeric time stamps
    times <- sort(as.numeric(difftime(retweets_niwa ,min(retweets_niwa),units = "mins")))

    # NIWA seal first tweet according to Twitter
    start <- lubridate::ymd_hms("2019-02-05 02:27:07")
    ## NIWA seal found owner tweet according to Twitter
    end <- lubridate::ymd_hms("2019-02-07 06:50:08")
    ## hist
    hist(retweets_niwa, breaks = "hours", axes = FALSE, 
         xlab = "", ylab = "",main = "", col = "grey",border = "grey",freq = TRUE)
    axis(2,at = c(0,250))
    mtext(2, line = 0, text = "Number of retweets per hour",cex = 0.8)
    mtext(1,line = 1, at = c(start,end),text = c("NIWA \n tweeted","USB owner \n found"),cex = 0.9)

![Observed counts of retweet
times.](hawkes_files/figure-markdown_strict/plot%20hist-1.png)

### Model fitting using `TMB`

Before using the `TMB` templates in `stelfi` you should use
`compile.stelfi()` to compile them.

    compile.stelfi()

    params <- c(mu = 9,alpha = 3,beta = 10)
    ## must have compiled TMB templates first use compile.stelfi()
    fit <- fit_hawkes(times = times,parameters = params) 

    ## print out estimated parameters
    fit

    ## sdreport(.) result
    ##         Estimate  Std. Error
    ## mu    0.06461416 0.017781393
    ## alpha 0.08026173 0.007992007
    ## beta  0.08355835 0.008334237
    ## Maximum gradient component: 3.84318

    par <- fit$par.fixed
    par

    ##         mu      alpha       beta 
    ## 0.06461416 0.08026173 0.08355835

![Fitted intensity (top plot) and observed counts of retweets
(bottom).](hawkes_files/figure-markdown_strict/plot-1.png)
