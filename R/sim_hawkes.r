#' Simulate a self-exciting Hawkes process
#'
#' Simulates a self-exciting Hawkes process. 
#' 
#' @details Option of two methods to simulate a Hawkes process: 
#' if \code{method = "1"} then a univariate Hawkes process as 
#' algorithm 2 in Ogata, Y. (1981) is simulated,
#' if \code{method = "2"} then an accept/reject
#' framework is used).
#' @param mu A numeric specifying base rate of the Hawkes process.
#' @param alpha A numeric specifying intensity jump after an event
#' occurrence.
#' @param beta A numeric specifying exponential intensity decay
#' @param n A numeric depending on method: if \code{method = "1"} specifies
#' end of the time line within which to simulate the process, 
#' if \code{method = "2"} specifies the number of observations to simulate.
#' Default, \code{100}.
#' @param plot Logical, if \code{TRUE} data plotted along with the intensity.
#' Default, \code{FALSE}.
#' @param seed The seed. Default, \code{123}
#' @param method A character "1" or "2" specifying the method (see details)
#' to simulate Hawkes process. Default,\code{"1"}.
#' @return A \code{numeric} vector of simulated event times.
#' 
#' @examples
#' sim_hawkes(10.2, 3.1, 8.9)
#' sim_hawkes(10.2, 3.1, 8.9, method = "2")
#' 
#' @seealso \code{\link{fit_hawkes}}
#' @export
sim_hawkes <- function(mu, alpha, beta, n = 100, plot = FALSE, 
                       seed = 123, method = "1") {
    if (alpha >= beta) {
        stop("alpha must be less then beta")
    }
    if (!method %in% c("1", "2")) {
        stop("method should be either \"1\" or \"2\"")
    }
    set.seed(seed)
    times <- numeric()
    if (method == "1") {
        s <- 0
        t <- 0
        lambda_star <- mu
        s <- s - log(runif(1)) / lambda_star
        t <- s
        dlambda <- alpha
        times <- c(times, t)
        while (s < n) {
            U <- runif(1)
            s <- s - log(U) / lambda_star
            u <- stats::runif(1)
            if (u <= (mu + dlambda * exp(-beta * (s-t))) / lambda_star) {
                dlambda <- alpha + dlambda * exp(-beta * (s - t))
                lambda_star <- lambda_star + alpha
                t <- s
                times <- c(times, t)
            }
        }
        if (plot) show_hawkes(list(times = times,
                                   params = c(mu = mu, alpha = alpha, beta = beta)))
        return(times)
    }else{
        if (method == "2"){
            eps <- 1e-6
            S <- function(k) {
                if (k == 1) return(1)
                return(1 + S(k - 1) * exp(-beta * (times[k] - times[k - 1])))
            }
            fu <- function(u, U, k) {
                return(log(U) + mu * (u - times[k]) +
                       alpha / beta * S(k) * (1 - exp(-beta * (u - times[k]))))
            }
            fu_prime <- function(u, U, k) {
                return(mu + alpha * S(k) * (exp(-beta * (u - times[k]))))
            }
            solve_u <- function(U, k) {
                u_prev <- times[k] - log(U) / mu
                u_next <- u_prev - fu(u_prev, U, k) / fu_prime(u_prev, U, k)
                while (abs(u_next - u_prev) > eps) {
                    u_prev <- u_next
                    u_next <- u_prev - fu(u_prev, U, k) / fu_prime(u_prev, U, k)
                }
                return(0.5 * (u_prev + u_next))
            }
            t1 <- -log(runif(1)) / mu
            times <- c(times, t1)
            k <- length(times)
            while (k < n) {
                U <- stats::runif(1)
                t_next <- solve_u(U, k)
                times <- c(times, t_next)
                k <- length(times)
            }
            if (plot) show_hawkes(list(times = times,
                                       params =  c(mu = mu, alpha = alpha, beta = beta)))
            return(times)
        }
    }
}
          
