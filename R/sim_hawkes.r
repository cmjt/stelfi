#' Function to simulate a simple univariate Hawkes process
#'
#' Code to simulate a realisation from a Hawkes process (code from R  package 'hawkes')
#' 
#'
#' @param mu numeric integer specifying base rate of the hawkes process
#' @param alpha numeric integer specifying intensity jump after an event occurence
#' @param beta numeric integer specifying exponential intensity decay
#' @param maxT numeric specifying end of the time line within which to simulat the proces, by default = 10
#' @param plot.proc logical if TRUE the process is plotted along with the intensity
#' @param seed by default = 1
#'
#' @export
 
simulate_hawkes <- function(mu = NULL,alpha = NULL,beta = NULL,maxT = 10, plot = FALSE,seed = 123){
    if(alpha >= beta){stop("alpha must be less then beta to ensure intensity decreases quicker than new events increase it")}
    set.seed(seed)
    times <- numeric()
    s <- 0
    t <- 0
    lambda_star <- mu
    s <- s -log(runif(1))/lambda_star
    t <- s
    dlambda <- alpha
    times <- c(times, t)
    while(s < maxT) {
        U <- runif(1)
        s <- s -log(U)/lambda_star
        u <- runif(1)
        if(u <= (mu + dlambda*exp(-beta*(s-t)))/lambda_star){
            dlambda <- alpha + dlambda*exp(-beta*(s-t))
            lambda_star <- lambda_star + alpha
            t <- s
            times <- c(times, t)
        }
    }
    if(plot){
        plot.hawkes(times, mu, alpha, beta)
    }
    return(times)
}

#' Using function from https://radhakrishna.typepad.com/mle-of-hawkes-self-exciting-process.pdf
#' @inheritParams simulate_hawkes
#' @param n number of timepoints to simulate
#' @export

rhawkes <- function(n = NULL,mu = NULL, alpha = NULL, beta = NULL,plot = FALSE, seed = 123){
    set.seed(seed)
    times <- numeric()
    eps <- 1e-6
    ## recurrence function
    S <- function(k){
        if(k==1) return(1)
        return( 1 + S(k-1)* exp(-beta * (times[k]- times[k-1])) )
    }
    ## the function to be solved to obtain, u
    fu <- function(u, U, k){
        return(log(U) + mu*(u - times[k]) +
               alpha/beta * S(k) * (1 - exp(-beta * (u - times[k]))))
    }
    fu_prime <- function(u, U , k){
        return( mu + alpha * S(k) * (exp(-beta * (u - times[k]))))
    }
    ## iterative procedure to solve f(u)
    solve_u <- function(U,k) {
        u_prev <- times[k] - log(U)/mu
        u_next <- u_prev - fu(u_prev, U, k)/ fu_prime(u_prev, U, k)
        while( abs(u_next - u_prev) > eps){
            u_prev <- u_next
            u_next <- u_prev - fu(u_prev, U, k)/fu_prime(u_prev, U, k)
        }
        return(0.5 * (u_prev + u_next))
    }
    t1 <- -log(runif(1))/mu
    times <- c( times, t1)
    k <- length(times)
    while(k < n) {
        U <- runif(1)
        t_next <- solve_u( U, k)
        times <- c(times, t_next)
        k <- length(times)
    }
    if(plot){
        plot.hawkes(times, mu, alpha, beta)
    }
    return(times)
}
#' Plotting intensity of a hawkes process
#' @inheritParams simulate_hawkes
#' @export

plot.hawkes <- function(times = NULL, mu = NULL, alpha = NULL, beta = NULL){
    n <- length(times)
    max <- max(times)
    p <- seq(0,max,length.out = 500)
    lam.p <- hawke.intensity(mu = mu, alpha = alpha, beta = beta, times = times, p = p)
    ylab <- expression(lambda(t))
    col <- 1
    lmax <- max(lam.p)
    lmin <- min(lam.p)
    plot(times,rep(lmin-1,n),ylim = c(lmin-2,lmax),xlab="time",ylab = ylab,col = col,pch=20)
    lines(p,lam.p,col="grey")
}
#' Hawkes intensty function with decay historical dependence 
#' @inheritParams simulate_hawkes
#' @export
hawke.intensity <- function(mu, alpha, beta,times,p = NULL){
    if(is.null(p)) p <- times
    lam <- function(p){
        mu + alpha*sum(exp(-beta*(p-times))[times<p])
    }
    lam.p <- rep(0,length(p))
    for(i in 1:length(p)){
        lam.p[i] <- lam(p[i])
    }
    lam.p
}

