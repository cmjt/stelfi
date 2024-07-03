#' Hawkes intensity function
#' 
#' Calculate the Hawkes intensity given a vector
#' of time points and parameter values
#' 
#' @inheritParams sim_hawkes
#' @param p An optional vector of pseudo times at which to calculate the intensity.
#' @inheritParams fit_hawkes_cbf
#' @noRd
hawkes_intensity <- function(times, mu, alpha, beta,
                            p, marks = rep(1, length(times)), 
                            background_parameters) {
    if (missing(p)) p <- times
    lam <- function(p, mu) {
      A <- exp(-beta * (p - times))[times < p]
      A <- append(A, rep(0, length(times) - length(A)))
      mu + alpha * sum(marks * A)
    }
    lam_p <- rep(0, length(p))
    if (inherits(mu, "function")) {
        mus <- mu(background_parameters, p)
    } else {
        mus <- rep(mu, length(p))
    }
    for (i in seq_along(p)) {
        lam_p[i] <- lam(p[i], mus[i])
    }
    return(lam_p)
}
#' Multivariate Hawkes intensity function
#' @inheritParams fit_mhawkes
#' @param p An optional vector of pseudo times at which to calculate the intensity.
#' @inheritParams fit_hawkes_cbf
#' @noRd
multi_hawkes_intensity <- function(times, mu, alpha, beta,
                           p, stream = rep("1", length(times))) {
    if (missing(p)) p <- times
    n_stream <- length(table(stream))
    if( n_stream != length(beta))
        stop("beta should be a vector of length equal to the number of streams")
    if( n_stream != length(mu))
        stop("mu should be a vector of length equal to the number of streams")   
    ## error checks for multivariate
    if(n_stream > 1){
        if( n_stream != nrow(alpha) |  n_stream != ncol(alpha))
            stop("alpha should be a matrix with dims matching the number of streams (> 1)")
    }
    ## split times by streams
    splt <- split(times, stream)
    streams <-  rep(NA, length(p))
    idx <- findInterval(times, p)
    stream_idx <- match(stream,  names(table(stream)))
    streams[idx] <- stream_idx
    tmp <- streams |> as.data.frame()
    streams <- tidyr::fill(tmp, streams)$streams
    lam_p <- se <- list()
    for(k in 1:n_stream){
        lam_p[[k]] <-  rep(0, length(p))
        se[[k]] <- list()
        ## self-exciting components
        se[[k]] <- lapply(1:n_stream,
                          function(n) stelfi:::hawkes_intensity(splt[[n]], 0, alpha[k, n], beta[n], p))
        lam_p[[k]] <- mu[k] + Reduce('+', se[[k]])
    }
    return(lam_p)
}
#' Extract reported parameter estimates
#'
#' Return parameter estimates from a fitted model.
#' 
#' @param obj A fitted model as returned by one of  \code{\link{fit_hawkes}}, \code{\link{fit_hawkes_cbf}}, \code{\link{fit_lgcp}},  \code{\link{fit_mlgcp}}, or \code{\link{fit_stelfi}}.
#' @return A matrix of estimated parameters and standard errors returned by
#' \code{TMB::sdreport} (\code{"report"}).
#' @seealso  \code{\link{fit_hawkes}},  \code{\link{fit_hawkes_cbf}}, \code{\link{fit_lgcp}},  \code{\link{fit_mlgcp}}, and \code{\link{fit_stelfi}}
#' @examples
#' ## Hawkes
#' data(retweets_niwa, package = "stelfi")
#' times <- unique(sort(as.numeric(difftime(retweets_niwa, min(retweets_niwa),units = "mins"))))
#' params <- c(mu = 9, alpha = 3, beta = 10)
#' fit <- fit_hawkes(times = times, parameters = params)
#' get_coefs(fit)
#' ## LGCP
#' if(requireNamespace("fmesher")) {
#' data(xyt, package = "stelfi")
#' domain <- sf::st_as_sf(xyt$window)
#' locs <- data.frame(x = xyt$x, y = xyt$y)
#' bnd <- fmesher::fm_as_segm(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = 0.75, cutoff = 0.3)
#' fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
#' parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1)))
#' get_coefs(fit)
#' }
#' @export
get_coefs <- function(obj) {
    table <- summary(TMB::sdreport(obj), "report")  
    ## The code below is for fit_hawkes_cbf()
    if("background_parameters" %in% names(obj)) {
        for (j in 1:length(obj$background_parameters)) {
            table <- rbind(table, c(obj$background_parameters[j], NA))
            row.names(table)[j + 2] <- paste("BP", j)
        }
    }
    return(table)
}

#' Extract the compensator differences
#'
#' Extract the compensator differences from a fitted Hawkes
#' model.
#' @inheritParams show_hawkes_GOF
#' @export
#'

#' Estimated random field(s)
#' 
#' Extract the estimated mean, or standard deviation, of the 
#' values of the Gaussian Markov random field for a fitted log-Gaussian
#' Cox process model at each node of \code{smesh}.
#' @seealso \code{\link{show_hawkes_GOF}}
#' @export
#' @rdname show_hawkes
compensator_differences <- function(obj){
    x <- show_hawkes_GOF(obj,plot = FALSE,
                         return_values = TRUE, tests = FALSE)$compensator_differences
    return(x)
}
#' 
#' @param obj A fitted model object returned by \code{\link{fit_lgcp}}.
#' @param plot Logical, if \code{TRUE} then the returned values are plotted.
#' Default \code{FALSE}.
#' @param sd Logical, if \code{TRUE} then standard errors returned.
#' Default \code{FALSE}.
#' @inheritParams fit_lgcp
#' @return A \code{numeric} vector or a \code{list} of returned values at each \code{smesh} node.
#' @examples \donttest{
#' if(requireNamespace("fmesher")) {
#' data(xyt, package = "stelfi")
#' domain <- sf::st_as_sf(xyt$window)
#' locs <- data.frame(x = xyt$x, y = xyt$y)
#' bnd <- fmesher::fm_as_segm(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = 0.75, cutoff = 0.3)
#' fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
#' parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1)))
#' get_fields(fit, smesh, plot = TRUE)
#' }
#' }
#' @seealso \code{\link{fit_lgcp}} and \code{\link{fit_mlgcp}}
#' @export
get_fields <- function(obj, smesh, tmesh, plot = FALSE, sd = FALSE) {
    if(missing(smesh)) stop("argument smesh missing")
    idx <- ifelse(sd, 2, 1)
    x <- cbind(summary(TMB::sdreport(obj),"random")[,idx])
    if(!missing(tmesh)) {
        ind <- rep(seq(tmesh$n), each = smesh$n)
        x <- split(x, ind)
        x <- do.call('cbind', x)
        pre <- ifelse(sd, "gmrf_sd_", "gmrf_mean_")
        colnames(x) <- paste(pre, "time_index_", seq(tmesh$n), sep = "")
        if(plot) {
            for(i in seq(tmesh$n)) {
                dev.new()
                print(show_field(x[, i], smesh))
            }
        }
    }else{
        if("betapp" %in% names(obj$par)){
            num_f <-  as.numeric(table(names(obj$par))["log_kappa"])
            ind <- rep(seq(num_f), each = smesh$n)
            x <- split(x, ind)
            x <- do.call('cbind', x)
            pre <- ifelse(sd, "sd_", "mean_")
            colnames(x) <- paste(pre, "field_", seq(num_f), sep = "")
            if(plot) {
                for(i in seq(num_f)) {
                    dev.new()
                    print(show_field(x[, i], smesh))
                }
            }
        }else{    
            colnames(x) <- ifelse(sd, "gmrf_sd", "gmrf_mean")
            if(plot) print(show_field(x, smesh))
        }}
    return(x)
}
#' Mesh weights
#' 
#' Calculate the  areas (weights) around the mesh nodes that 
#' are within the specified spatial polygon \code{sf} of the domain.
#'
#' @param mesh A spatial mesh of class \code{fmesher::fm_mesh_2d()} 
#' @param sf An \code{sf} of type \code{POLYGON} specifying the region
#' of the domain.
#' @param plot Logical, whether to plot the calculated \code{mesh} weights. 
#' Default, \code{FALSE}.
#' @return Either a simple features, \code{sf}, object or values returned by \code{geom_sf}.
#' @examples 
#' data(horse_mesh, package = "stelfi")
#' data(horse_sf, package = "stelfi")
#' get_weights(horse_mesh, horse_sf, plot = TRUE)
#' @seealso \url{https://becarioprecario.bitbucket.io/spde-gitbook/}.
#' @export
get_weights <- function(mesh, sf, plot = FALSE) {
    dmesh <- dual_mesh(mesh)
    sf::st_crs(sf) <- NA
    if (!all(sf::st_is_valid(sf))) {
        ## check for invalid geometries
        sf <- sf::st_make_valid(sf)
    }
    w <- sapply(1:nrow(dmesh), function(i) {
        coord <- sf::st_coordinates(dmesh[i, ])[, 1:2]
        coord <- sf::st_polygon(list(coord))
        if (any(sf::st_intersects(coord, sf, sparse = FALSE)))
            return(sf::st_area(sf::st_intersection(coord, sf)))
        else return(0)
    })
    weights <- data.frame(ID = 1:nrow(dmesh), weights = unlist(w))
    res <- dplyr::left_join(dmesh, weights, by = "ID")
    if(plot) {
        p <-  ggplot2::ggplot(res) +
            ggplot2::geom_sf(ggplot2::aes(fill = weights)) +
            ggplot2::geom_sf(data = mesh_2_sf(mesh),
                             col = "white", fill = NA) +
            ggplot2::geom_sf(data = sf, fill = NA) +
            ggplot2::theme_void()
        print(p)
    }else{
        return(res)
    }
}

#' Internal function that takes in a list of points
#' (and optionally a weight or covariate for each point)
#' and a mesh of polygons.
#' @param xy A data frame of locations.
#' @param dmesh An object returned by \code{\link{dual_mesh}}.
#' @param weights A vector of weights/covariates for each point.
#' @noRd
points_in_mesh <- function(xy, dmesh, weights) {
  xy <- sf::st_as_sf(xy, coords = c("x", "y"))
  xy <- sf::st_geometry(xy)
  if (missing(weights)) {
    sapply(1:nrow(dmesh), function(i) {
      coord <- sf::st_coordinates(sf::st_as_sf(dmesh[i, ]))[, 1:2]
      coord <- sf::st_polygon(list(coord))
      sum(sf::st_contains(coord, xy, sparse = FALSE) > 0)
    })
  }
  else {
    result <- rep(0, nrow(dmesh))
    for (i in 1:nrow(dmesh)) {
      coord <- sf::st_coordinates(sf::st_as_sf(dmesh[i, ]))[, 1:2]
      coord <- sf::st_polygon(list(coord))
      temp_result <- sf::st_contains(coord, xy, sparse = FALSE)
      temp_result <- temp_result * weights
      result[i] <- sum(temp_result)
    }
    return(result)
  }
}
#' Internal function to construct the dual mesh
#'
#' @param mesh A spatial mesh of class \code{fmesher::fm_mesh_2d()}.
#' @return An\ simple features, code{sf}, object of the Voronoi tessellation
#' centered at each \code{mesh} node.
#' @source \url{http://www.r-inla.org/spde-book},  \url{https://becarioprecario.bitbucket.io/spde-gitbook/}
#' @noRd
dual_mesh <- function(mesh) {
    if (mesh$manifold == 'R2') {
        ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
            colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
        pls <- lapply(1:mesh$n, function(i) {
            p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                j <- which(mesh$graph$tv[, k] == i)
                if (length(j) > 0)
                    return(rbind(ce[j, , drop = FALSE],
                                 cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 1],
                                       mesh$loc[mesh$graph$tv[j, k], 2] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 2]) / 2))
                else return(ce[j, , drop = FALSE])
            })))
            j1 <- which(mesh$segm$bnd$idx[, 1] == i)
            j2 <- which(mesh$segm$bnd$idx[, 2] == i)
            if ((length(j1) > 0) | (length(j2) > 0)) {
                p <- unique(rbind(mesh$loc[i, 1:2], p,
                                  mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2] / 2 +
                                  mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2] / 2,
                                  mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2] / 2 +
                                  mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2] / 2))
                yy <- p[, 2] - mean(p[, 2]) / 2 - mesh$loc[i, 2] / 2
                xx <- p[, 1] - mean(p[, 1]) / 2 - mesh$loc[i, 1] / 2
            }
            else {
                yy <- p[, 2] - mesh$loc[i, 2]
                xx <- p[, 1] - mesh$loc[i, 1]
            }
            p <- p[order(atan2(yy, xx)), ]
            ## close polygons
            p <- rbind(p, p[1, ])
            sf::st_polygon(list(p))
        })
        geometry <- sf::st_sfc(pls)
        dat <- data.frame(ID = 1:length(geometry))
        return(sf::st_sf(dat,
                         geometry = geometry))
    }
    else stop("It only works for R2!")
}
