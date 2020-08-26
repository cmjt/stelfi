## ---- libs
library(stelfi)

## ---- break
##########################################################
##########################################################
##########################################################

## ---- earthquakes
data(earthquakes)
times <- earthquakes$origintime
##  2012
times <- subset(times, lubridate::year(times) == "2012")

## ---- hist_earth
h <- hist(times, breaks = "days",axes = FALSE,, 
          xlab = "", ylab = "",main = "", col = "grey",freq = TRUE)
axis(2,at = c(0,25))
mtext(2, line = 1, text = "Number of earthquakes per day")
abline(v =  h$mids[seq(1,365,by = 31)], col = "grey", lwd = 2, lty = 2)
p <- c(h$mids[seq(1,365,by = 31)],h$mids[length(h$mids)])
roll.mean <- numeric(12)
for(i in 1:12){roll.mean[i] <- mean(c(p[i],p[(i+1)]))}
text(x = roll.mean , y = rep(25,12),labels = names(table(lubridate::month(times,abbr = TRUE,label = TRUE))),
     col = "grey")

## ---- fit_earth
t.idx <- sort(as.numeric(difftime(times, min(times), unit = "days")))
params <- c(mu = 10,alpha = 13,beta = 20)
fit <- fit_hawkes(times = t.idx,parameters = params)

## ---- inference_earth
fit
par <- fit$par.fixed
show_hawkes(t.idx, mu = par[1], alpha = par[2],beta = par[3])

## ---- break
##########################################################
##########################################################
##########################################################

## ---- data_iq
data(terrorism)
## make date vector
times <- lubridate::dmy(paste(terrorism$iday, terrorism$imonth, terrorism$iyear, sep = "-"),tz = "Asia/Baghdad")

## ---- hist_iq
h <- hist(times, breaks = "days",axes = FALSE,, 
          xlab = "", ylab = "",main = "", col = "grey",freq = TRUE)
axis(2,at = c(0,30))
mtext(2, line = 1, text = "Number of attacks per day")
abline(v =  h$mids[seq(1,2920,by = 365)], col = "grey", lwd = 2, lty = 2)
p <- c(h$mids[seq(1,length(h$mids),by = 365)],h$mids[length(h$mids)])
roll.mean <- numeric(length(p) - 1)
for(i in 1:length(roll.mean)){roll.mean[i] <- mean(c(p[i],p[(i+1)]))}
text(x = roll.mean , y = rep(30,8),labels =  names(table(lubridate::year(times))), col = "grey")

## --- 2016
times <- subset(times, lubridate::year(times) == "2016")


## ---- hist2_iq
h <- hist(times, breaks = "days",axes = FALSE,, 
          xlab = "", ylab = "",main = "", col = "grey",freq = TRUE)
axis(2,at = c(0,20))
mtext(2, line = 1, text = "Number of attacks per day")
abline(v =  h$mids[seq(1,365,by = 31)], col = "grey", lwd = 2, lty = 2)
p <- c(h$mids[seq(1,365,by = 31)],h$mids[length(h$mids)])
roll.mean <- numeric(12)
for(i in 1:12){roll.mean[i] <- mean(c(p[i],p[(i+1)]))}
text(x = roll.mean , y = rep(20,12),labels = names(table(lubridate::month(times,abbr = TRUE,label = TRUE))),
     col = "grey")


t.idx <- lubridate::yday(times)

## ---- fit
## t.idx <- as.numeric(difftime(times, min(times), unit = "weeks"))
## t.idx <- sort(unique(t.idx)) ## unique so dealing with aggregate days ignoring "number"
## t.idx <- t.idx/max(t.idx) ## hmm need to think about this
## jitter as can't have exactly the same timestamps
set.seed(123)
jitter <- sort(runif(length(t.idx)))
t.idx <- sort(t.idx + jitter)
params <- c(mu = 0.1,alpha = 0.3,beta = 0.5)
## must have compiled TMB templates first use compile.stelfi()
fit <- fit_hawkes(times = t.idx,parameters = params)#, lower = c(0,0.1,0.2),
                  upper = c(30,15,30),
                  method = "L-BFGS-B")
fit
## ---- inference
par <- fit$par.fixed
show_hawkes(t.idx, mu = par[1], alpha = par[2],beta = par[3])




