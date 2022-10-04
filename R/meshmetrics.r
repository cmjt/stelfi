#' Calculate the distance between two points
#' @param A A numeric vector of length 2 specifying (first) location.
#' @param B A numeric vector of length 2 specifying (second) location.
#' @examples 
#' dist(c(1, 1), c(2, 1))
#' @noRd
dist <- function(A, B) {
  sqrt((A[1] - B[1])^2 + (A[2] - B[2])^2)
}
#' Cosine triangle function to calculate angle C (in degrees)
#' given edge lengths a, b, and c
#' 
#' @param a A numeric triangle edgelength (opposite angle A).
#' @param b A numeric triangle edgelength (opposite angle B).
#' @param c A numeric triangle edgelength (opposite angle C).
#' @examples
#' ang(1, 2, 1.5)
#' @noRd
ang <- function(a, b, c) {
  cosC <- (a^2 + b^2 - c^2) / (2 * a * b)
  angC <- acos(cosC)
  return(angC * (180 / pi))
}
#' Sin function to calculate the area of a triangle given
#' edge lengths a, b and angle C.
#' 
#' @param a A numeric triangle edgelength (opposite angle A).
#' @param b A numeric triangle edgelength (opposite angle B).
#' @param angC angle (in degrees)
#' @examples 
#' A <- c(1, 1)
#' B <- c(2, 1)
#' C <- c(2.5, 1.5)
#' a <- dist(B, C)
#' b <- dist(A, C)
#' c <- dist(A, B)
#' angC <- ang( a, b, c)
#' tri_area(a, b, angC)
#' @noRd
tri_area <- function(a, b, angC) {
    (1 / 2) * a * b * sin(angC / (180 / pi))
}
#' Calculate cicumcircle radius given vertex A, B, C
#' locations of a triangle
#' 
#' @inheritParams incircle_r
#' @examples
#' circum_R(c(1, 1), c(2, 1), c(3, 2.5))
#' @noRd
circum_R <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  abc <- a * b * c
  d1 <- a + b + c
  d2 <- b + c - a
  d3 <- c + a - b
  d4 <- a + b - c
  return((abc) / (sqrt(d1 * d2 * d3 * d4)))
}
#' Calculate cicumcircle centroid given vertex A, B, C
#' locations of a triangle
#' 
#' @inheritParams incircle_r
#' @examples
#' circum_O(c(1, 1), c(2, 1), c(3, 2.5))
#' @noRd
circum_O <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  angA <- ang(b, c, a) * (pi / 180)
  angB <- ang(a, c, b) * (pi / 180)
  angC <- ang(b, a, c) * (pi / 180)
  sumsins <- sin(2 * angA) + sin(2 * angB) + sin(2 * angC)
  xo <- (A[1] * sin(2 * angA) + B[1] * sin(2 * angB) + C[1] * sin(2 * angC)) / sumsins
  yo <- (A[2] * sin(2 * angA) + B[2] * sin(2 * angB) + C[2] * sin(2 * angC)) / sumsins
  return(c(xo, yo))
}
#' Calculate incumcircle radius given vertex A, B, C
#' locations of a triangle
#' 
#' @inheritParams incircle_r
#' @examples 
#' incircle_O(c(1, 1), c(2, 1), c(3, 2.5))
#' @noRd
incircle_O <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  abc <- a + b + c
  xc <- (a * A[1] + b * B[1] + c *  C[1]) / abc
  yc <- (a * A[2] + b * B[2] + c *  C[2]) / abc
  return(c(xc, yc))
}
#' Calculate incumcircle centroid given vertex A, B, C
#' locations of a triangle
#' 
#' @param A A numeric vector of length 2 specifying vertex location "A".
#' @param B A numeric vector of length 2 specifying vertex location "B".
#' @param C A numeric vector of length 2 specifying vertex location "C".
#' @examples 
#' incircle_r(c(1, 1), c(2, 1), c(3, 2.5))
#' @noRd
incircle_r <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  s <- (a + b + c) / 2
  sqrt(((s - a) * (s - b) * (s - c)) / s)
}
#' Extract a dataframe of mesh triangle segments (start and end locations)
#' 
#' @inheritParams meshmetrics
#' @examples 
#' data(horse_mesh, package = "stelfi")
#' segments(horse_mesh)
#' @noRd
segments <- function(mesh) {
  df <- rbind(data.frame(a = mesh$loc[mesh$graph$tv[, 1], c(1, 2)],
                         b = mesh$loc[mesh$graph$tv[, 2], c(1, 2)]),
              data.frame(a = mesh$loc[mesh$graph$tv[, 2], c(1, 2)],
                         b = mesh$loc[mesh$graph$tv[, 3], c(1, 2)]),
              data.frame(a = mesh$loc[mesh$graph$tv[, 1], c(1, 2)],
                         b = mesh$loc[mesh$graph$tv[, 3], c(1, 2)]))
  colnames(df) <- c("x", "y", "xend", "yend")
  df$length <- c(unlist(dist(df[, 1:2], df[, 3:4])))
  return(df)
}
#' Calculate all interior Delaunay triangulation triangle angles
#' 
#' @inheritParams meshmetrics
#' @examples 
#' data(horse_mesh, package = "stelfi")
#' mesh_ang(horse_mesh)
#' @noRd
mesh_ang <- function(mesh) {
  tv <- mesh$graph$tv
  angs <- matrix(numeric(3 * nrow(tv)), ncol = 3)
  for(i in 1:nrow(tv)) {
    ## the three verts of one triangle
    vs <- data.frame(x = rep(mesh$loc[tv[i,], 1], each = 2),
                     y = rep(mesh$loc[tv[i,], 2], each = 2))
    vs$xend <- vs$x[c(3, 5, 1, 5, 1, 3)];vs$yend <- vs$y[c(3, 5, 1, 5, 1, 3)]
    vs <- vs[c(1,2,4), ]
    dists <- c(unlist(dist(vs[, 1:2], vs[, 3:4])))
    angs[i, 1] <- ang(dists[1], dists[2], dists[3])
    angs[i, 2] <- ang(dists[2], dists[3], dists[1])
    angs[i, 3] <- ang(dists[3], dists[1], dists[2])
  }
  angs <- as.data.frame(angs)
  names(angs) <- c("angleA", "angleB", "angleC")
  angs$ID <- 1:nrow(angs)
  return(angs)
}
#' Calculate minimum edge length for each triangle in a given Delaunay triangulation
#' 
#' @inheritParams meshmetrics
#' @examples
#' data(horse_mesh, package = "stelfi")
#' lmin(horse_mesh)
#' @noRd
lmin <- function(mesh) {
  tv <- mesh$graph$tv
  lmin <- numeric(nrow(tv))
  for (i in 1:nrow(tv)) {
    A <- mesh$loc[tv[i, 1], 1:2]
    B <- mesh$loc[tv[i, 2], 1:2]
    C <- mesh$loc[tv[i, 3], 1:2]
    a <- dist(B, C)
    b <- dist(A, C)
    c <- dist(A, B)
    lmin[i] <- min(a, b, c)
  }
  return(lmin)
}
#' Transform a \code{INLA::inla.mesh.2d} into a \code{sf} object
#' 
#' @inheritParams meshmetrics
#' @source Modified from \code{sp} based function suggested by Finn in the
#' R-inla discussion Google Group
#' \url{https://groups.google.com/g/r-inla-discussion-group/c/z1n1exlZrKM}.
#' @return A simple features, \code{sf}, object.
#' @seealso \code{\link{meshmetrics}}
#' @examples
#' data(horse_mesh, package = "stelfi")
#' sf <- mesh_2_sf(horse_mesh)
#' if(require("ggplot2")) {
#' ggplot(sf) + geom_sf()
#' }
#' @export
mesh_2_sf <- function(mesh) {
    st <- sf::st_sfc(lapply(
                  1:nrow(mesh$graph$tv),
                  function(x) {
                      tv <- mesh$graph$tv[x, , drop = TRUE]
                      sf::st_polygon(list(mesh$loc[tv[c(1, 3, 2, 1)],
                                                   1:2,
                                                   drop = FALSE]))
                  }
              )
              )
    dat <- as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE])
    dat$ID <- 1:length(st)
    res <- sf::st_sf(dat,
                     geometry = st)
    return(res)
}
#' Calculate a number of different geometric attributes of a Delaunay triangulation
#' 
#' Calculates a number of geometric attributes for a given
#' Delaunay triangulation based on the circumscribed and inscribed circle of each triangle.
#' 
#' @param mesh A \code{INLA::inla.mesh.2d()} object.
#' @return An object of class \code{sf} with the following data for each triangle in the
#' triangulation
#' \itemize{
#' \item \code{V1}, \code{V2}, and \code{V3} corresponding vertices
#' of \code{mesh} matches \code{mesh$graph$tv};
#' \item \code{ID}, numeric triangle id;
#' \item \code{angleA}, \code{angleB}, and \code{angleC}, the
#' interior angles;
#' \item circumcircle radius, circumradius, \code{circumcircle_R} (\eqn{R});
#' \item incircle radius \code{incircle_r} (\eqn{r});
#' \item centroid locations of the circumcircle, circumcenter, (\code{c_Ox, cOy});
#' \item centroid locations of the incircle, incenter, (\code{i_Ox, iOy});
#' \item the radius-edge ratio \code{radius_edge} \eqn{\frac{R}{l_{min}}},
#' where \eqn{l_{min}} is the minimum edge length;
#' \item the radius ratio \code{radius_ratio} \eqn{\frac{r}{R}};
#' \item \code{area}, area (\eqn{A});
#' \item \code{quality} a measure of "quality" defined as
#' \eqn{\frac{4\sqrt{3}|A|}{\Sigma_{i = 1}^3 L_i^2}},
#' where \eqn{L_i} is the length of edge \eqn{i}.
#' }
#' @details A triangle's circumcircle (circumscribed circle) is the unique circle that passes
#' through each of its three vertices. A triangle's incircle (inscribed circle) is the
#' largest circle that can be contained within it (i.e., touches it's three edges).
#' @examples
#' data(horse_mesh, package = "stelfi")
#' metrics <- meshmetrics(horse_mesh)
#' if(require("ggplot2")) {
#' ggplot(metrics) + geom_sf(aes(fill = radius_ratio))
#' }
#' @export
meshmetrics <- function(mesh) {
  angles <- mesh_ang(mesh = mesh)
  tv <- mesh$graph$tv
  c_R <- i_R <- area <- quality <- numeric(nrow(tv))
  c_O <- i_O <- matrix(rep(0, 2 * nrow(tv)), ncol = 2)
  for (i in 1:nrow(tv)) {
    A <- mesh$loc[tv[i, 1], 1:2]
    B <- mesh$loc[tv[i, 2], 1:2]
    C <- mesh$loc[tv[i, 3], 1:2]
    a <-  dist(B, C)
    b <- dist(A, C)
    c <- dist(A, B)
    c_R[i] <- circum_R(A, B, C)
    i_R[i] <- incircle_r(A, B, C)
    c_O[i, ] <- circum_O(A, B, C)
    i_O[i, ] <- incircle_O(A, B, C)
    area[i] <- abs(tri_area(a = a, b = b,
                            angC = ang(a, b, c)))
    quality[i] <- (4 * sqrt(3) * abs(area[i])) / (sum(c(a, b, c)^2))
  }
  mn <- lmin(mesh)
  sf <- mesh_2_sf(mesh)
  df <- data.frame(ID = 1:nrow(sf), incircle_r = i_R,
                   circumcircle_R = c_R,
                   c_Ox = c_O[, 1],  c_Oy = c_O[, 2],
                   i_Ox = i_O[, 1],  i_Oy = i_O[, 2],
                   radius_edge = c_R / mn, radius_ratio = i_R / c_R,
                   area = area, quality = quality)
  sf <- dplyr::left_join(sf, angles, by = "ID")
  sf <- dplyr::left_join(sf, df, by = "ID")
  return(sf)
}
