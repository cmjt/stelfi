#' A package to fit Hawkes and Log-Gaussian Cox Point Process models
#' using Template Model Builder
#'
#' Fit Hawkes and log-Gaussian Cox process models with extensions.
#' Introduced in Hawkes (1971) a Hawkes process is a
#' self-exciting temporal point process where the occurrence of an event
#' immediately increases the chance of another. We extend this to consider
#' self-inhibiting process and a non-homogeneous background rate.
#' A log-Gaussian Cox process is a Poisson point process where the
#' log-intensity is given by a Gaussian random field. We extend this
#' to a joint likelihood formulation fitting a marked log-Gaussian Cox model.
#' In addition, the package offers functionality to fit self-exciting
#' spatiotemporal point processes. Models are fitted via maximum likelihood
#' using 'TMB' (Template Model Builder) (Kristensen, Nielsen, Berg, Skaug,
#' and Bell, 2016). Where included 1) random fields are assumed to be
#' Gaussian and are integrated over using the Laplace approximation and
#' 2) a stochastic partial differential equation model, introduced by
#' Lindgren, Rue, and Lindström. (2011), is defined for the field(s).
#'
#' @section Model fitting:
#'
#' \itemize{
#' \item The functions \code{\link{fit_hawkes}} and \code{\link{fit_hawkes_cbf}}
#' fit self-exciting Hawkes (Hawkes AG., 1971) processes to temporal point pattern data.
#' \item The function \code{\link{fit_lgcp}} fit a log-Gaussian Cox process to
#' either spatial or spatiotemporal point pattern data. If a spatiotemporal
#' model is fitted a AR1 process is assumed for the temporal progression.
#' \item The function \code{\link{fit_mlgcp}} fits a joint likelihood model between
#' the point locations and the mark(s).
#' \item The function \code{\link{fit_stelfi}} fits self-exciting spatiotemporal
#' Hawkes models to point pattern data. The self-excitement is Gaussian in space
#' and exponentially decaying over time. In addition, a GMRF can be included
#' to account for latent spatial dependency.
#' }
#'
#' @references Hawkes, AG. (1971) Spectra of some self-exciting and
#' mutually exciting point processes. \emph{Biometrika}, \strong{58}: 83--90.
#'
#' @references Lindgren, F., Rue, H., and Lindström, J. (2011)
#' An explicit link between {G}aussian fields and {G}aussian {M}arkov random
#' fields: the stochastic partial differential equation approach.
#' \emph{Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology)}, \strong{73}: 423--498.
#'
#' @references Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., and
#'  Bell B. M. (2016). {TMB}: Automatic Differentiation and Laplace
#' Approximation. \emph{Journal of Statistical Software}, \strong{70}: 1--21.
#' @keywords internal
"_PACKAGE"
NULL
