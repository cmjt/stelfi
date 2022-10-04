#' Plot Hawkes intensity
#'
#' Plots a Hawkes intensity function, options to extend to 
#' non-homogeneous background intensity.
#' 
#' @param obj Either object returned by \code{fit_hawkes()}/\code{fit_hawkes_cbf()} 
#' or a named list with elements \code{times} and \code{params}. If the latter, 
#' then \code{times} should be a numeric vector of observed time points, and 
#' \code{params} must contain, \code{alpha} (intensity jump after an event occurrence) 
#' and \code{beta} (exponential intensity decay). In addition, should contain 
#' either \code{mu} (base rate of the Hawkes process) or \code{background_parameters}
#' (parameter(s) for the assumed non-homogeneous background function;
#' could be a list of multiple values). May also contain \code{marks} (a vector of numerical marks).
#' 
#' @seealso [hawkes_intensity()] [fit_hawkes_cbf()]
#' @examples \dontrun{
#' data(retweets_niwa, package = "stelfi")
#' times <- unique(sort(as.numeric(difftime(retweets_niwa, min(retweets_niwa),units = "mins"))))
#' params <- c(mu = 9, alpha = 3, beta = 10)
#' show_hawkes(list(times = times, params = params))
#' ## must have compiled TMB templates first use compile_stelfi()
#' fit <- fit_hawkes(times = times, parameters = params)
#' show_hawkes(fit)
#' }
#' @export
show_hawkes <-  function(obj) {
    if (!"times" %in% names(obj)) {
        times <- obj$env$data$times
        marks <- obj$env$data$marks
        pars <- get_coefs(obj)
        if ("background_parameters" %in% names(obj)) {
            alpha <- pars[1,1]
            beta <- pars[2,1]
            background_parameters = obj$background_parameters
        } else {
            mu <- pars[1,1]
            alpha <- pars[2,1]
            beta <- pars[3,1]
        }
    } else {
        times <- obj$times
        if(!is.null(obj$marks)) {
            marks <- obj$marks
        } else {
            marks <- rep(1, length(times))
        }
        pars <- obj$params
        alpha <- pars[["alpha"]]
        beta <- pars[["beta"]]
        if ("background_parameters" %in% names(pars)) {
            background_parameters <- pars$background_parameters
        } else {
            mu <- pars[["mu"]]
        }  
    }
    n <- length(times)
    max <- max(times)
    p <- seq(0, max, length.out = 500)
    lam.p <- hawkes_intensity(mu = mu, alpha = alpha, beta = beta, times = times,
                             p = p, marks = marks, background_parameters = background_parameters)
    ylab <- expression(lambda(t))
    col <- 1
    lmax <- max(lam.p)
    lmin <- min(lam.p)
    data <-  data.frame(xp = p, yp = lam.p)
    line <- ggplot2::ggplot(data = data,
                            ggplot2::aes(x = .data$xp, y = .data$yp)) +
        ggplot2::xlab("") +
        ggplot2::ylab(expression(lambda(t))) + 
        ggplot2::geom_line() +  ggplot2::theme_minimal()
    data <- data.frame(times = times)
    hist <-  ggplot2::ggplot(data = data,  ggplot2::aes(x = .data$times)) +
        ggplot2::geom_histogram() +  ggplot2::theme_minimal() +
        ggplot2::xlab("times") +  ggplot2::ylab("Number of events")
    gridExtra::grid.arrange(line, hist, ncol = 1)
}
#' Compensator and other goodness-of-fit metrics for a Hawkes process
#'
#' Plots a number of goodness-of-fit plots for a fitted
#' Hawkes process. Includes 1) a comparison of the  compensator and observed 
#' events, 2) a histogram of transformed interarrival times, 3) a Q-Q plot of 
#' transformed interarrival times, and 4) the CDF of consecutive interarrival 
#' times, In addition, results of a Kolmogorov-Smirnov and
#' Ljung-Box hypothesis test for the interarrival times are printed.
#' 
#' @param plot Logical, whether to plot  goodness-of-fit plots. Default \code{TRUE}.
#' @param return_values Logical, whether to return GOF values. Default \code{FALSE}.
#' @examples \dontrun{
#' data(retweets_niwa, package = "stelfi")
#' times <- unique(sort(as.numeric(difftime(retweets_niwa, min(retweets_niwa),units = "mins"))))
#' params <- c(mu = 9, alpha = 3, beta = 10)
#' show_hawkes_GOF(list(times = times, params = params))
#' ## must have compiled TMB templates first use compile_stelfi()
#' fit <- fit_hawkes(times = times, parameters = params)
#' show_hawkes_GOF(fit)
#' }
#' @export
#' @rdname show_hawkes
show_hawkes_GOF <-  function(obj, plot = TRUE, return_values = FALSE) {
    if (!"times" %in% names(obj)) {
        times <- obj$env$data$times
        marks <- obj$env$data$marks
        pars <- get_coefs(obj)
        if ("background_parameters" %in% names(obj)) {
            alpha <- pars[1,1]
            beta <- pars[2,1]
            background_parameters = obj$background_parameters
        } else {
            mu <- pars[1,1]
            alpha <- pars[2,1]
            beta <- pars[3,1]
        }
    } else {
        times <- obj$times
        if(!is.null(obj$marks)) {
            marks <- obj$marks
        } else {
            marks <- rep(1, length(times))
        }
        pars <- obj$params
        alpha <- pars[["alpha"]]
        beta <- pars[["beta"]]
        if ("background_parameters" %in% names(pars)) {
            background_parameters <- pars$background_parameters
        } else {
            mu <- pars[["mu"]]
        }  
    }
    A <- numeric(length = length(times))
    for(i in 2:length(times)) {
        A[i] <- exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1])
    }
    compensator <- numeric(length = length(times))
    if (!inherits(mu, "function")) {
        for(i in 1:length(times)) {
            compensator[i] <- (mu * times[i]) - ((alpha/beta)*A[i]) +
                ((alpha / beta) * (sum(marks[1:i])-marks[i]))
        }
    } else {
        for(i in 1:length(times)) {
            compensator[i] <- mu(background_parameters, times[i]) - ((alpha/beta) * A[i]) +
                ((alpha / beta) * (sum(marks[1:i]) - marks[i]))
        }
        compensator <- compensator - mu(background_parameters,0) ## Subtract integral at zero
    }
    interarrivals <- compensator[2:length(compensator)] - compensator[1:(length(compensator)-1)]
    ## Kolmogorov-Smirnov Test
    print(stats::ks.test(interarrivals, "pexp"))
    ## Ljung-Box Test
    print(stats::Box.test(interarrivals, type = "Ljung"))
    if (plot) {
        data <- data.frame(xs = times, observed = 1:length(times), compensator = compensator)
        ## Plot of compensator and actual events
        data <- reshape(data, direction = "long", idvar = "xs",
                        varying = c("observed", "compensator"), v.names = "val",
                        times = c("observed", "compensator"),
                        new.row.names = 1:(2*length(times)))
        lineplot <- ggplot2::ggplot(data = data,
                                    ggplot2::aes(x = .data$xs, y = .data$val, colour = .data$time)) +
            ggplot2::xlab("Time") +
            ggplot2::ylab("Events") +
            ggplot2::geom_line() +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position=c(0.8,0.2)) +
            ggplot2::ggtitle("Actual Events and Compensator")
        ## Histogram of transformed interarrival times
        data <- data.frame(data = interarrivals)
        hist <-  ggplot2::ggplot(data = data,  ggplot2::aes(x = .data$data)) +
            ggplot2::geom_histogram(stat = "density") +  ggplot2::theme_minimal() +
            ggplot2::xlab("Interarrival times") +  ggplot2::ylab("Density") +
            ggplot2::stat_function(fun = "dexp", args = (mean = 1), color = "red") +
            ggplot2::ggtitle("Transformed Interarrival Times")
        p <- ppoints(100)    ## 100 equally spaced points on (0,1), excluding endpoints
        q <- quantile(interarrivals,p = p) ## percentiles of the sample distribution
        ## Q-Q plot of transformed interarrival times
        data <- data.frame(x = qexp(p), y = q)
        qqplot <- ggplot2::ggplot(data =  data,
                                  ggplot2::aes(x = .data$x, y = .data$y)) +
            ggplot2::xlab("Theoretical Quantiles") +
            ggplot2::ylab("Observed Quantiles") + 
            ggplot2::geom_point() +  ggplot2::theme_minimal() + 
            ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
            ggplot2::ggtitle("Transformed Interarrival Times")
        U <- numeric(length=length(interarrivals))
        U <- 1 - exp(-compensator[2:length(compensator)] + compensator[1:(length(compensator) - 1)])
        ## Scatterplot of the CDF of consecutive interarrival times
        data <-  data.frame(x = U[1:(length(U)-1)], y = U[2:length(U)])
        scatter <- ggplot2::ggplot(data = data,
                                   ggplot2::aes(x = .data$x, y = .data$y)) +
            ggplot2::xlab("F(Interarrival time k)") +
            ggplot2::ylab("F(Interarrival time k+1)") + 
            ggplot2::geom_point() +  ggplot2::theme_minimal() +
            ggplot2::ggtitle("Consecutive Interarrival Times")
        gridExtra::grid.arrange(lineplot, qqplot, hist, scatter, ncol = 2)
    }
    if(return_values) {
        return(list(interarrivals = interarrivals))
    }
}
#' Plot the estimated random field(s) of a fitted LGCP
#'
#' Plots the values of \code{x} at each node of \code{smesh}, with
#' optional control over resolutions using \code{dims}.
#' 
#' @param x A vector of values for each \code{smesh} node.
#' @param sf Optional, \code{sf} of type \code{POLYGON} specifying the region
#' of the domain.
#' @inheritParams show_lambda
#' @seealso \code{\link{show_lambda}} and \code{\link{get_fields}}
#' @examples \dontrun{
#' data(xyt, package = "stelfi")
#' domain <- sf::st_as_sf(xyt$window)
#' bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- INLA::inla.mesh.2d(boundary = bnd, 
#' max.edge = 0.75, cutoff = 0.3)
#' parameters <- c(beta = 1, log_tau = log(1), log_kappa = log(1))
#' simdata <- simulate_lgcp(parameters = parameters, sf = domain, smesh = smesh)
#' show_field(simdata$x, smesh = smesh, sf = domain)
#' }
#' @export
show_field <- function(x, smesh, sf, dims = c(500,500)) {
    nx <- dims[1]
    ny <- dims[2]
    xs <- seq(min(smesh$loc[, 1]), max(smesh$loc[, 1]), length = nx)
    ys <- seq(min(smesh$loc[, 2]), max(smesh$loc[, 2]), length = ny)
    data <- expand.grid(xs = xs, ys = ys)
    pxl <- sf::st_multipoint(as.matrix(data))
    A <- INLA::inla.spde.make.A(smesh, pxl)
    data$colz <-  as.vector(A %*% x)
    plt <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = data, ggplot2::aes(x = .data$xs, y = .data$ys, 
                                                     fill = .data$colz)) +
        ggplot2::labs(fill = "") +
        ggplot2::xlab("") + ggplot2::ylab("") + 
        ggplot2::scale_fill_viridis_c(option = "D") +
        ggplot2::coord_equal() 
    if (!missing(sf)) {
        plt <- plt +
            ggplot2::geom_sf(data = sf, fill = NA, size = 2)
    }
    plt

}
#' Plot the estimated intensity from a fitted LGCP model
#'
#' Plots the estimated spatial intensity from
#' a fitted log-Gaussian Cox process model. If \code{obj} is a
#' spatiotemporal model then \code{timestamp} provides control
#' over which temporal index to plot the estimated spatial intensity.
#' 
#' @param obj A fitted LGCP model object for, for example, \code{fit_lgcp()}.
#' @param dims A numeric vector of length 2 specifying
#' the spatial pixel resolution. Default \code{c(500,500)}.
#' @param timestamp The index of time stamp to plot. Default \code{1}.
#' @inheritParams fit_lgcp
#' @seealso \code{\link{show_field}} and \code{\link{get_fields}}
#' @examples \dontrun{
#' data(xyt, package = "stelfi")
#' domain <- sf::st_as_sf(xyt$window)
#' locs <- data.frame(x = xyt$x, y = xyt$y)
#' bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- INLA::inla.mesh.2d(boundary = bnd, max.edge = 0.75, cutoff = 0.3)
#' fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
#' parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1)))
#' show_lambda(fit, smesh = smesh, sf = domain)
#' }
#' @export
show_lambda <- function(obj, smesh, sf, tmesh,
                        covariates, dims = c(500,500),
                        timestamp = 1) {
    if(!missing(tmesh)) {
        if(!missing(covariates)) {
            designmat <- cbind(1, covariates)
        } else {
            designmat <- matrix(rep(1, smesh$n*tmesh$n), ncol = 1)
        }
        
        res <- TMB::sdreport(obj)
        field <- res$par.random
        beta <- res$value["beta"]
        beta <- as.matrix(beta, ncol = 1)
        lambda <- exp(field + designmat%*%beta)
        ind <- rep(seq(tmesh$n), each = smesh$n)
        x <- split(lambda, ind)
        plt <- list()
        for(i in seq(tmesh$n)) {
            if (missing(sf)) {
                plt[[i]] <- show_field(x = x[[i]], smesh = smesh, dims = dims)
            } else {
                plt[[i]] <- show_field(x = x[[i]], smesh = smesh, sf = sf, dims = dims)
            }
        }
        plt[[timestamp]]
    }else{
        if(!missing(covariates)) {
            designmat <- cbind(1, covariates)
        } else {
            designmat <- matrix(rep(1, smesh$n), ncol = 1)
        }
        
        res <- TMB::sdreport(obj)
        field <- res$par.random
        beta <- res$value["beta"]
        beta <- as.matrix(beta, ncol = 1)
        lambda <- exp(field + designmat%*%beta)
        if (missing(sf)) {
            show_field(x = lambda, smesh = smesh,dims = dims)
        } else {
            show_field(x = lambda, smesh = smesh, sf = sf, dims = dims)
        }
    }           
}




