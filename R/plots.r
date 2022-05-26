setClassUnion("missing_or_numeric", c("numeric", "missing"))
#' Plots the Hawkes intensty function with decay historical dependence
#' @docType methods
#' @rdname show_hawkes
#' @inheritParams sim_hawkes
#' @param times a vector of numeric observation time points
#' @export
setGeneric("show_hawkes",
           function(times, mu, alpha, beta, marks=c(rep(1,length(times))), background_param=NULL){
               standardGeneric("show_hawkes")
           })

setMethod("show_hawkes",
          c(times = "vector", alpha = "numeric", beta  = "numeric"),
          function(times, mu, alpha, beta, marks, background_param){
              n = length(times)
              max = max(times)
              p = seq(0,max,length.out = 500)
              lam.p = hawke_intensity(mu = mu, alpha = alpha, beta = beta, times = times, 
                                      p = p, marks= marks, background_param = background_param)
              ylab = expression(lambda(t))
              col = 1
              lmax = max(lam.p)
              lmin = min(lam.p)
              data = data.frame(x = p, y = lam.p)
              line = ggplot2::ggplot(data = data, ggplot2::aes(x = .data$x, y = .data$y)) +
                  ggplot2::xlab("") +
                  ggplot2::ylab(expression(lambda(t))) + 
                  ggplot2::geom_line() +  ggplot2::theme_minimal()
              hist =  ggplot2::ggplot(data = data.frame(times = times),  ggplot2::aes(x = times)) +
                  ggplot2::geom_histogram() +  ggplot2::theme_minimal() +
                  ggplot2::xlab("times") +  ggplot2::ylab("Number of events")
              gridExtra::grid.arrange(line, hist, ncol = 1)
          })
setGeneric("show_hawkes_GOF", # only for constant mu at this stage
           function(times, mu, alpha, beta, marks=c(rep(1,length(times))),
                    background_param){
                   standardGeneric("show_hawkes_GOF")
           })

setMethod("show_hawkes_GOF",
          c(times = "vector", alpha = "numeric", beta  = "numeric"),
          function(times, mu, alpha, beta, marks, background_param){
                  
                  A <- numeric(length=length(times))
                  for(i in 2:length(times)){
                          A[i] <- exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1])
                  }
                  
                  compensator = numeric(length=length(times))
                  if (class(mu) == "numeric"){
                        for(i in 1:length(times)){
                                compensator[i] <- (mu * times[i]) - ((alpha/beta)*A[i])+ ((alpha / beta) * (i - 1))
                        }
                  } else {
                          for(i in 1:length(times)){
                                  compensator[i] <- mu(background_param,times[i]) - ((alpha/beta)*A[i])+ ((alpha / beta) * (i - 1))
                          }
                          compensator = compensator - mu(background_param,0) # Subtract integral at zero
                          
                  }
                  
                  interarrivals <- compensator[2:length(compensator)] - compensator[1:(length(compensator)-1)]

                  data <- data.frame(x = times, observed = 1:length(times), compensator = compensator)
                  melted <- reshape2::melt(data, id.vars = "x")
                  
                  # Plot of compensator and actual events
                  lineplot = ggplot2::ggplot(data = melted, ggplot2::aes(x = x, y = value, colour = variable)) +
                          ggplot2::xlab("Time") +
                          ggplot2::ylab("Events") +
                          ggplot2::geom_line() +
                          ggplot2::theme_minimal() +
                          ggplot2::theme(legend.position=c(0.8,0.2)) +
                          ggplot2::ggtitle("Actual Events and Compensator")
                  
                  
                  # Histogram of transformed interarrival times
                  hist =  ggplot2::ggplot(data = data.frame(data = interarrivals),  ggplot2::aes(x = data)) +
                          ggplot2::geom_histogram(ggplot2::aes(y=..density..)) +  ggplot2::theme_minimal() +
                          ggplot2::xlab("Interarrival times") +  ggplot2::ylab("Density") +
                          ggplot2::stat_function(fun="dexp", args=(mean=1), color = "red") +
                          ggplot2::ggtitle("Transformed Interarrival Times")
                  
                  p <- ppoints(100)    # 100 equally spaced points on (0,1), excluding endpoints
                  q <- quantile(interarrivals,p=p) # percentiles of the sample distribution
                  data <- data.frame(x = qexp(p), y = q)
                  # Q-Q plot of transformed interarrival times
                  qqplot = ggplot2::ggplot(data = data, ggplot2::aes(x = .data$x, y = .data$y)) +
                          ggplot2::xlab("Expected Quantiles") +
                          ggplot2::ylab("Observed Quantiles") + 
                          ggplot2::geom_point() +  ggplot2::theme_minimal() + 
                          ggplot2::geom_abline(intercept = 0,slope = 1, color = "red") +
                          ggplot2::ggtitle("Transformed Interarrival Times")

                  U <- numeric(length=length(interarrivals))
                  U <- 1 - exp(-compensator[2:length(compensator)] + compensator[1:(length(compensator)-1)])
                  
                  data <- data.frame(x = U[1:(length(U)-1)], y = U[2:length(U)])
                  # Scatterplot of consecutive interarrival times
                  scatter = ggplot2::ggplot(data = data, ggplot2::aes(x = .data$x, y = .data$y)) +
                          ggplot2::xlab("Interarrival time k") +
                          ggplot2::ylab("Interarrival time k+1") + 
                          ggplot2::geom_point() +  ggplot2::theme_minimal() +
                          ggplot2::ggtitle("Consecutive Interarrival Times")
                
                  # Kolmogorov-Smirnov Test
                  KS <- stats::ks.test(interarrivals, "pexp")$p.value
                  # Ljung-Box Test
                  LBQ <- stats::Box.test(interarrivals, type = "Ljung")$p.value
                  
                  title = paste("LBQ p-value = ", round(LBQ,3), "KS p-value = ", round(KS,3))
                  gridExtra::grid.arrange(lineplot, qqplot, hist, scatter, ncol = 2, bottom = grid::textGrob(
                          title,
                          gp = grid::gpar(fontface = 3, fontsize = 9),
                          hjust = 1,
                          x = 1
                  ))
                  return(list(interarrivals=interarrivals, KS = KS))
                  })
#' Plots estimated random field(s) of a LGCP
#' @docType methods
#' @rdname show_field
#' @param x a vector of values for each mesh node
#' @param dims by default \code{c(500,500)}, vector of length 2 specifying
#' spatial pixel resolution
#' @inheritParams get_fields
#' @export
show_field <- function(x,smesh, dims = c(500,500), border, colour_option="D",
                       title="Plot"){
     nx <- dims[1]
     ny <- dims[2]
     px <- inlabru::pixels(smesh, nx = nx, ny = ny)
     A <- INLA::inla.spde.make.A(smesh, px)
     px$color <- as.vector(A %*% x)
     # Convert border data from Spatial Polygons to dataframe
     if (!missing(border)) {
        border_f<-ggplot2::fortify(border)
        }

     plt <- ggplot2::ggplot(as.data.frame(px), ggplot2::aes(x, y)) +
         ggplot2::geom_tile(ggplot2::aes(fill = color)) +
         ggplot2::scale_fill_viridis_c(option=colour_option) +
         ggplot2::coord_equal()+
         ggplot2::ggtitle(title)
     
     if (!missing(border)){
     plt<- plt+ggplot2::geom_polygon(data=border_f,ggplot2::aes(x=long,y=lat,group=group),fill=NA,color='red')
     }
     plt

 }
#' Plot the lambda estimate from a spatial model

plot_lambda <-function(fit,covariates, smesh,dims=c(500,500), border, colour_option="D",
                       title="Plot"){
        if(!missing(covariates)) {
                designmat <- cbind(1, covariates)
        } else {
                designmat <- matrix(rep(1, smesh$n), ncol = 1)
        }
        
        res = TMB::sdreport(fit)
        field = res$par.random
        beta = res$value[1:ncol(designmat)]
        beta = as.matrix(beta,ncol=1)
        lambda = field + designmat%*%beta
        if (missing(border)){
                stelfi::show_field(x=lambda,smesh=smesh,dims=dims,
                                  colour_option=colour_option, title=title)
        } else {
                stelfi::show_field(lambda,smesh,dims,border,colour_option,title)
        }
        
}
        
        