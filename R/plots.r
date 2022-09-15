#' Plots the Hawkes intensity function with decay historical dependence
#' 
#' @inheritParams sim_hawkes
#' @inheritParams fit_hawkes
#' @inheritParams fit_hawkes_cbf
#' @param times A numeric vector of observed time points.
#' @param obj An object returned by fit_hawkes() or fit_hawkes_cbf()
#' If provided, the only other parameters required are \code{mu} and \code{background_parameters}, and only if \code{mu} is not constant. 
#' @examples \dontrun{
#' data(retweets_niwa, package = "stelfi")
#' times <- unique(sort(as.numeric(difftime(retweets_niwa ,min(retweets_niwa),units = "mins"))))
#' params <- c(mu = 9, alpha = 3, beta = 10)
#' ## must have compiled TMB templates first use compile_stelfi()
#' fit <- fit_hawkes(times = times, parameters = params)
#' pars <- get_coefs(fit)
#' show_hawkes(times = times, mu = pars[1,1], alpha = pars[2,1], beta = pars[3,1])
#' }
#' @export
show_hawkes <-  function(times = NULL, mu, alpha, beta, marks = rep(1, length(times)), 
                         obj = NULL, background_parameters) {
    if (!is.null(obj)) {
        times = obj$env$data$times
        marks = obj$env$data$marks
        pars = get_coefs(obj)
        if ("background_parameters" %in% names(obj)) {
            alpha = pars[1,1]
            beta = pars[2,1]
            background_parameters = obj$background_parameters
        } else {
            mu = pars[1,1]
            alpha = pars[2,1]
            beta = pars[3,1]
        }
    }
    n <- length(times)
    max <- max(times)
    p <- seq(0, max, length.out = 500)
    lam.p <- hawke_intensity(mu = mu, alpha = alpha, beta = beta, times = times,
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
#' Plot the compensator and other associated goodness of fit metrics for a Hawkes process
#' 
#' @inheritParams show_hawkes
#' @param plot Logical, whether to plot GOF plots. Default \code{TRUE}.
#' @param return_values Logical, whether to return GOF values. Default \code{TRUE}.
#' @examples \dontrun{
#' data(retweets_niwa, package = "stelfi")
#' times <- unique(sort(as.numeric(difftime(retweets_niwa,
#' min(retweets_niwa),units = "mins"))))
#' params <- c(mu = 9, alpha = 3, beta = 10)
#' ## must have compiled TMB templates first use compile_stelfi()
#' fit <- fit_hawkes(times = times, parameters = params)
#' pars <- get_coefs(fit)
#' show_hawkes_GOF(times = times, mu = pars[1,1], alpha = pars[2,1],
#' beta = pars[3,1], return_values = FALSE)
#' }
#' @export
show_hawkes_GOF <-  function(times, mu, alpha, beta, marks, obj, background_parameters, plot, return_values) {
    
    if (!is.null(obj)) {
        times = obj$env$data$times
        marks = obj$env$data$marks
        pars = get_coefs(obj)
        if ("background_parameters" %in% names(obj)) {
            alpha = pars[1,1]
            beta = pars[2,1]
            background_parameters = obj$background_parameters
        } else {
            mu = pars[1,1]
            alpha = pars[2,1]
            beta = pars[3,1]
        }
    }
    A <- numeric(length = length(times))
    for(i in 2:length(times)) {
        A[i] <- exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1])
    }
    compensator <- numeric(length = length(times))
    if (class(mu) == "numeric") {
        for(i in 1:length(times)) {
            compensator[i] <- (mu * times[i]) - ((alpha/beta)*A[i]) +
                ((alpha / beta) * (sum(marks[1:i])-marks[i]))
        }
    } else {
        for(i in 1:length(times)) {
            compensator[i] <- mu(background_parameters,times[i]) - ((alpha/beta)*A[i]) +
                ((alpha / beta) * (sum(marks[1:i])-marks[i]))
        }
        compensator <- compensator - mu(background_parameters,0) ## Subtract integral at zero
        
    }
    interarrivals <- compensator[2:length(compensator)] - compensator[1:(length(compensator)-1)]
    data <- data.frame(xs = times, observed = 1:length(times), compensator = compensator)
    ## Plot of compensator and actual events
    data <- reshape2::melt(data, id.vars = "xs",, value.name = "val")
    lineplot <- ggplot2::ggplot(data = data,
                                ggplot2::aes(x = .data$xs, y = .data$val, colour = .data$variable)) +
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
    
    ## Kolmogorov-Smirnov Test
    KS <- stats::ks.test(interarrivals, "pexp")$p.value
    ## Ljung-Box Test
    LBQ <- stats::Box.test(interarrivals, type = "Ljung")$p.value
    
    title <- paste("LBQ p-value = ", round(LBQ,3), "KS p-value = ", round(KS,3))
    if (plot) {
        gridExtra::grid.arrange(lineplot, qqplot, hist, scatter, ncol = 2,
                                bottom = grid::textGrob(
                                                   title,
                                                   gp = grid::gpar(fontface = 3, fontsize = 9),
                                                   hjust = 1,
                                                   x = 1
                                               ))
    }
    if(return_values) {
        return(list(interarrivals = interarrivals, KS = KS, LBQ = LBQ))
    }
}
#' Plot estimated random field(s) of a fitted LGCP
#' 
#' @param x A vector of values for each \code{smesh} node.
#' @param dims A numeric vector of length 2 specifying
#' the spatial pixel resolution. By default \code{ = c(500,500)}.
#' @param sf optional, an \code{sf} of type \code{POLYGON} specifying the region
#' of the domain.
#' @inheritParams get_fields
#' @inheritParams fit_lgcp
#' @seealso \code{\link{plot_lambda}} and \code{\link{get_fields}}
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
show_field <- function(x, smesh, dims = c(500,500), sf) {
    nx <- dims[1]
    ny <- dims[2]
    xs <- seq(min(smesh$loc[, 1]), max(smesh$loc[, 1]), length = nx)
    ys <- seq(min(smesh$loc[, 2]), max(smesh$loc[, 2]), length = ny)
    data <- expand.grid(xs = xs, ys = ys)
    pxl <- sf::st_multipoint(as.matrix(data))
    A <- INLA::inla.spde.make.A(smesh, pxl)
    data$colz <-  as.vector(A %*% x)
    plt <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = data, ggplot2::aes(x = .data$xs, y = .data$ys, fill = .data$colz)) +
        ggplot2::labs(fill = "") + 
        ggplot2::scale_fill_viridis_c(option = "D") +
        ggplot2::coord_equal() 
    
    if (!missing(sf)) {
        plt <- plt +
            ggplot2::geom_sf(data = sf, fill = NA, size = 2)
    }
    plt

}
#' Plot the lambda estimate from a fitted LGCP model
#' 
#' @inheritParams get_fields
#' @inheritParams fit_lgcp
#' @param fit a fitted model object for, for example, \code{fit_lgcp()}.
#' @param dims 2D spatial resolution to plot, default \code{c(500,500)}.
#' @param timestamp index of time stamp to plot, default 1.
#' @seealso \code{\link{show_field}} and \code{\link{get_fields}}
#' @examples \dontrun{
#' data(xyt, package = "stelfi")
#' domain <- sf::st_as_sf(xyt$window)
#' locs <- data.frame(x = xyt$x, y = xyt$y)
#' bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- INLA::inla.mesh.2d(boundary = bnd, 
#' max.edge = 0.75, cutoff = 0.3)
#' fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
#' parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1)))
#' plot_lambda(fit, smesh = smesh, sf = domain)
#' }
#' @export
plot_lambda <- function(fit, covariates, smesh, tmesh,
                        dims = c(500,500),
                        sf,
                        timestamp = 1) {
    
    if(!missing(tmesh)) {
        if(!missing(covariates)) {
            designmat <- cbind(1, covariates)
        } else {
            designmat <- matrix(rep(1, smesh$n*tmesh$n), ncol = 1)
        }
        
        res <- TMB::sdreport(fit)
        field <- res$par.random
        beta <- res$value["beta"]
        beta <- as.matrix(beta, ncol = 1)
        lambda <- exp(field + designmat%*%beta)
        ind <- rep(seq(tmesh$n), each = smesh$n)
        x <- split(lambda, ind)
        plt <- list()
        for(i in seq(tmesh$n)) {
            if (missing(sf)) {
                plt[[i]] <- show_field(x = x[[i]], smesh = smesh,dims = dims)
            } else {
                plt[[i]] <- show_field(x = x[[i]], smesh, dims, sf)
            }
        }
        plt[[timestamp]]
    }else{
        if(!missing(covariates)) {
            designmat <- cbind(1, covariates)
        } else {
            designmat <- matrix(rep(1, smesh$n), ncol = 1)
        }
        
        res <- TMB::sdreport(fit)
        field <- res$par.random
        beta <- res$value["beta"]
        beta <- as.matrix(beta, ncol = 1)
        lambda <- exp(field + designmat%*%beta)
        if (missing(sf)) {
            show_field(x = lambda, smesh = smesh,dims = dims)
        } else {
            show_field(lambda, smesh, dims, sf)
        }
    }           
}




