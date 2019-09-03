TMB::compile("inst/src/hawkes.cpp")

dyn.load(TMB::dynlib("inst/src/hawkes"))


source("~/Git/selfE/R/sim_hawkes.r")
params <- c(mu = 0.57,alpha = 3,beta = 4)
simt <- rhawkes(params[1], params[2], params[3], n = 10,plot = TRUE)
simo <- simulate_hawkes(params[1], params[2], params[3],maxT = 10,plot = TRUE)
data <- list(times = simt)


obj <- TMB::MakeADFun(data = data, parameters = params,DLL = "hawkes")
opt <- stats::nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max = 1000, iter.max = 1000),silent = TRUE)
rep <- TMB::sdreport(obj)

obj$fn(obj$par)


loglik <- function(params, arrivals){
    alpha_i <- params[1]
    beta_i <- params[2]
    mu_i <- params[3]
    n <- length(arrivals)
    term_1 <- -mu_i*arrivals[n]
    term_2 <- sum(alpha_i/beta_i*(exp( -beta_i * (arrivals[n] - arrivals)) - 1))
    Ai <- c(0, sapply(2:n, function(z) {
        sum(exp( -beta_i * (arrivals[z]- arrivals[1:(z - 1)])))
    }))
    term_3 <- sum(log( mu_i + alpha_i * Ai))
    return(-term_1- term_2 -term_3)
}
loglik(c(params[2],params[3],params[1]),simt)

library(lgcpSPDE)

data(package = "lgcpSPDE")
terrorism <- subset(terrorism, terrorism$iyear == 2016)
terrorism <- subset(terrorism, terrorism$country == "Iraq")
t <- lubridate::ymd(paste(terrorism$iyear,terrorism$imonth, terrorism$iday,sep = ":"))
t <- na.omit(t)
t <- as.numeric(difftime(t,min(t),units = "days"))
data <- list(times = t)
params <- c(mu = 15,alpha = 5,beta = 17)
obj <- TMB::MakeADFun(data = data, parameters = params,DLL = "hawkes")
opt <- stats::nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max = 1000, iter.max = 1000),silent = TRUE)
rep <- TMB::sdreport(obj)
summary(rep)

plot.hawkes(times = t,mu = rep$par.fixed[1],alpha = rep$par.fixed[2],beta = rep$par.fixed[3])


plot(1:length(t),sapply(t, function(x) length(t[t < x]))) ## Hawkes process realisation...
abline(a = 0, b = 1,lty = 2, col = 2)

plot(table(cut(t,breaks = ceiling(length(t))))) ## compensator?

