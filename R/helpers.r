setClassUnion("numeric_or_missing", c("numeric", "missing"))
#' Hawkes intensity function
#' @inheritParams sim_hawkes
#' @inheritParams show_hawkes
#' @param p An optional vector of pseudo times at which to calculate the intensity.
#' @export
setGeneric("hawke_intensity",
           function(times, mu, alpha, beta, p, marks, background_param) {
           })

setMethod("hawke_intensity",
          c(times = "vector", alpha = "numeric",
            beta  = "numeric", p = "numeric_or_missing"),
          function(times, mu, alpha, beta,
                   p, marks, background_param) {
              if (missing(p)) p <- times
              lam <- function(p, mark, mu) {
                  mu + mark * alpha * sum(exp(-beta * (p - times))[times < p])
              }
              lam_p <- rep(0, length(p))
              if (class(mu) == "function") {
                mus <- mu(background_param, p)
              } else {
                mus <- rep(mu, length(p))
              }
              for (i in seq_along(p)) {
                  lam_p[i] <- lam(p[i], marks[i], mus[i])
              }
              return(lam_p)
          })
#' Reported parameter estimates
#'
#' \code{get_coefs} returns the parameter estimates for the fitted model.
#' @param object result of a call to \code{\link{fit_hawkes}} or \code{\link{fit_lgcp}}
#' @export
setGeneric("get_coefs",
           function(object) {
               standardGeneric("get_coefs")
           }
           )
setMethod("get_coefs",
          signature(object = "list"),
          function(object) {
              table <- summary(TMB::sdreport(object), "report")  
              ## The code below is for fit_hawkes_cbf()
              if("background_parameters" %in% names(object)) {
                for (j in 1:length(object$background_parameters)) {
                  table <- rbind(table, c(object$background_parameters[j], NA))
                  row.names(table)[j + 2] <- paste("BP", j)
              }
              }
          return(table)
          }
          )
#' Estimated random field(s)
#'
#' \code{get_fields} returns
#' @param object The result of a call to \code{\link{fit_lgcp}}.
#' @param plot Logical, if \code{TRUE} then the returned values are plotted.
#' @param sd Logical, if \code{TRUE} then standard errors of field are returned.
#' @export
get_fields <- function(object, smesh, tmesh, plot = FALSE, sd = FALSE) {
    idx <- ifelse(sd, 2, 1)
    x <- summary(TMB::sdreport(object),"random")[,idx]
    if(!missing(tmesh)) {
        ind <- rep(seq(tmesh$n), each = smesh$n)
        x <- split(x, ind)
        if(plot) {
            for(i in seq(tmesh$n)) {
                dev.new()
                print(show_field(x[[i]], smesh))
            }
        }
    }else{
        if(plot) print(show_field(x, smesh))
    }
    return(x)
}
      
#' Function to find areas (weights) around the mesh nodes which are
#' within the specified spatial polygon.
#' Relies on the internal \code{\link{dual_mesh}} function.
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
    if(missing(plot)) plot <- FALSE
    if(plot) {
        p <-  ggplot2::ggplot(res) +
            ggplot2::geom_sf(ggplot2::aes(fill = weights)) +
            ggplot2::geom_sf(data = stelfi:::mesh_2_sf(mesh),
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
#' @param weights A vector of weights/covariates for each point
#' @return Either the number of points in each polygon or sum of the weights.
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
#' Internal function to construct the `dual` mesh
#'
#' @return a \code{sf} object of the Voronoi tesselation
#' centered at each \code{mesh} node.
#' @seealso \url{http://www.r-inla.org/spde-book} and \url{https://becarioprecario.bitbucket.io/spde-gitbook/}
#' @source \url{http://www.r-inla.org/spde-book}
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
