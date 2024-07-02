#' Sim study
#'
#' @param n 
#' @param times 
#' @param locs 
#' @param param 
#' @param domain 
#' @param smesh 
#' @param time_independent 
#'
#' @return
#' @export
#'
#' @examples
sim_study <- function (n, times, locs, param, domain, smesh, time_independent=F) {
  message(paste0("Simulation study with n=", n))
  
  fit <- fit_stelfi(times = times, locs = locs, sf = domain, smesh = smesh, parameters = param, time_independent = time_independent)
  
  # Setup parallel computing stuff
  num_cores <- floor(parallel::detectCores() * 0.5)
  message(paste0("Running simulation study using ", num_cores, " cores"))
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`
  
  est_values <- foreach::foreach(i = 1:n, .packages = 'stelfi', .errorhandling = 'pass') %dopar% {
  # for (i in 1:n) {
    simulated_data <- fit$simulate()
    
    # Convert simulated data to correct format
    sim_locs <- data.frame(x = simulated_data$locs[,1], 
                           y = simulated_data$locs[,2])
    
    simulated_fit <- fit_stelfi(times = simulated_data$times, locs = sim_locs, 
                                sf = domain, smesh = smesh, parameters = param,
                                time_independent = time_independent)
    
    fit_coefs <- get_coefs(simulated_fit)
    
    return(list(fit_coefs))
  }
  
  parallel::stopCluster(cl)
  
  message("Finished sim study successfully.")
  
  saveRDS(est_values, paste0("stelfi/n_", n))
  
  return(est_values)
}
