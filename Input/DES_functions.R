# This script specifies several general functions relating to (discrete event) simulation and which are used across
# the different PRIMCAT models. Even though they are defined in this separate script, they are not necessarily 
# stand alone functions are refer to variables in the global environment. Therefor, this script is to be sourced 
# from within the respective functions that define the PRIMCAT models.
#
# The following functions are defined:
# - draw            randomly samples an event/index from a vector of probabilities
# - rcustom         convenience function to perform a random draw from any distribution
# - rmixture        custom function to sample from a mixture distribution
# - rgompertz       link to the flexsurv rgompertz function
# - rweibullPH      link to the flexsurv weibullPH function
# - priorityQueue   defines a event queue as backbone for a discrete event simulation
# - schedule        adds an event to the event queue
# - seize           seizes a resource by recording the start and end time
# - runOne          runs the discrete event simulation for one individual
# - runAll          runs the complete discrete event simulation for all individuals
#
# Version history:
# - May 2021        Koen Degeling       R v4.0.3    Initial version
# - November 2021   Koen Degeling       R v4.0.3    Updated functions and commenting

# The draw function randomly samples an event/index from a vector of probabilities. Note that there is no check to 
# ensure all probabilities end up to 1, nor is it possible to return anything else than the index, such as the names.
# It is programmed this way for speed because the base R sample function is slow for sampling 1 value. Example of 
# the 'probs' argument (names not necessary): 
# c(p_L1_SG = 0.6, p_L2_RT = 0.4)
draw <- function(probs) {
  if(is.null(probs)) stop("The 'probs'' provided to draw() are NULL")
  if(abs(sum(probs) - 1) > 1*10^(-15)) stop("The sum of the 'probs' provided to draw() is not 1")
  return(sum(cumsum(probs) < runif(1)) + 1)
}

# The rcustom function is a convenience function to perform a random draw from any distribution that is supported. 
# There are no checks other than that the distribution should be supported, which can be specified through the name
# in R (e.g. 'gompertz') or by its capitalized name (e.g. 'Gompertz'). The 'distlist' argument should be a list 
# object including a 'dist' object that specifies the name of the distribution and a 'pars' object that is a named 
# list specifying the parameter values. Note that the parameters need to be on real scale, i.e. the scale used in 
# function calls. If the dist is a standard parametric distribution (e.g., Weibull), the corresponding parameters 
# need to be provided (e.g., shape and scale for Weibull distributions). The order of the parameters is important,
# since the parameters are selected using an index rather than the parameter names. If dist is a mixture, the 
# following parameters have to be provided in the dist list:
# - d1_prob   the mixing proportion for the first distribution in the mixture
# - d1        the name of the first distribution, which needs to be a standard parametric distribution
# - d1_par1   the first parameter for the first distribution
# - d1_par2   the second parameter for the first distribution
# - d2        the name of the second distribution, which needs to be a standard parametric distribution
# - d2_par1   the first parameter for the second distribution
# - d2_par2   the second parameter for the second distribution
# Example of the 'distlist' argument:
# list(
#   dist = 'weibull',
#   pars = list(shape = 1.2, scale = 10)
# )
rcustom <- function(distlist, n = 1) {
  
  t <- Inf
  
  while(is.infinite(t)) { 
    
    t <- switch(
      distlist$dist,
      'gompertz'  = rgompertz (n = n, shape = distlist$pars[[1]], rate = distlist$pars[[2]]),
      'weibullPH' = rweibullPH(n = n, shape = distlist$pars[[1]], scale = distlist$pars[[2]]),
      'mixture'   = rmixture  (n = n, d1_prob = distlist$pars$d1_prob, 
                               d1_distlist = list(dist = distlist$pars$d1, pars = list(distlist$pars$d1_par1, distlist$pars$d1_par2)),
                               d2_distlist = list(dist = distlist$pars$d2, pars = list(distlist$pars$d2_par1, distlist$pars$d2_par2))),
      stop(paste0("Distribution type '", distlist$dist, "' not supported by the rcustom() function!"))
    )
    
    if(is.infinite(t)) warning('Infinite value sampled from distribution: ', distlist)
    
  }
  
  return(t)
  
}

# The rmixture function performs a draw from a mixture of two distributions. It uses a random number to first select
# which distribution is used based on the probability for the first distribution that is specified through the 
# 'd1_prob' argument. After selecting the distribution, the corresponding distlist (i.e., 'd1_distlist' or
# 'd2_distlist') is used to call the rcustom function to perform a draw from that distribution. See the information
# for the rcustom function about how this distlist needs to be specified.
rmixture   <- function(n, d1_prob, d1_distlist, d2_distlist) ifelse(runif(n) < d1_prob, rcustom(n = n, distlist = d1_distlist), rcustom(n = n, distlist = d2_distlist))

# The rgompertz and rweibullPH functions are simply links from the flexsurv package. This is done so that no
# packages are required to run the simulation. These functions could be replaced by hardcoded functions if the 
# runtime needs to be reduced, as the flexsurv functions are quite slow.
rgompertz  <- flexsurv::rgompertz
rweibullPH <- flexsurv::rweibullPH

# The priorityQueue function defines all the elements and functions required for using an event queue in the discrete
# event simulation. It defines the queue through a list of 'events' and a vector of 'times'. The functions defined 
# are:
# - push      add an event to the queue
# - pop       select and remove the first event from the queue
# - empty     check whether the queue is empty
# - clear     clear the queue
# - remove    remove a specific event from the queue based on a predicate
priorityQueue <- function() {
  
  times  <- numeric()
  events <- list()
  
  push <- function(event, time) {
    insert <-  findInterval(time, times)
    times  <<- append(times,  time,  insert)
    events <<- append(events, event, insert)
  }
  
  pop <- function() {
    head   <-  structure(events[[1]], time = times[1])
    times  <<- times[-1]
    events <<- events[-1]
    return(head)
  }
  
  empty <- function() {
    return(length(times) == 0)
  }
  
  clear <- function() {
    times  <<- numeric()
    events <<- list()
  }
  
  remove <- function(predicate, ...) {
    i <- sapply(events, predicate, ...)
    stopifnot(is.logical(i))
    i[is.na(i)] <- TRUE
    times  <<- times[!i]
    events <<- events[!i]
  }
  
  return(list(
    push   = push, 
    pop    = pop, 
    empty  = empty, 
    clear  = clear, 
    remove = remove
  ))
  
}

# The schedule functions adds an event to the event queue based on the name of the event that is specified through
# the 'event' argument and the time at which the event is to occur through the numerical 'time' argument. Note that
# the function is hardcoded to use an event queue that is named 'queue' as well as a variable 'currentTime' that are
# both to be specified in the parent environment.
schedule <- function(event, time) {
  attr(event, 'time') <- time
  attr(event, 'sendingTime') <- currentTime
  queue$push(event, time)
}

# The seize function records the use of a 'resource' by entering the start and end date to the 'resourceLog' array 
# based on the 'start' argument (default is 'currentTime' and specified 'duration'. Note that the 'resourceLog' 
# array is a 3-dimensional array with dimensions equal to: 1) the number of simulated individuals, 2) the number of
# resources, and 3) the start and end time. The second dimension is named according to the resources specified in 
# the parent environment and, hence, only prespecified resources can be seized. Additionally, the function uses 
# 'until', which specifies for how long the simulation runs, to ensure resources are not seized beyond the time 
# horizon of the simulation. The 'resourceLog', 'currentTime', 'id, and 'until' objects have to be specified in the
# parent environment before the simulation can be run.
seize <- function(resource, start = currentTime, duration = NA_real_) resourceLog[id, resource, ] <<- c(start, min(start + duration, until))

# The runOne function runs the discrete event simulation for one individual by first initializing the simulation for
# the individual by calling the init function that is to be specified and, subsequently, processing all events in 
# the event queue one by one until all events are processed. Note that it does not look at the maximum simulation
# time, so the maximum runtime specified by 'until' is used in the init function to schedule an event 'End' at time
# 'until' at which the queue is emptied. In processing the events, the function repeats the following steps:
# - select the next event using the pop function in queue
# - update the time of simulation to the time of the event
# - log the event to the eventLog
# - process the event, which may include scheduling new events
# Note that events are recorded in the 'eventLog' matrix with the number of rows equal to the number of simulated 
# individuals and the columns equal to the number of events. The possible events have to be prespecified in the 
# parent environment to define the eventLog matrix. The 'queue', 'currentTime', and 'previousEventTime' obviously 
# also have to be prespecified in the parent environment.
runOne = function() {
  init()
  while(!queue$empty()) {
    event <- queue$pop()
    currentTime <<- attr(event, 'time')
    eventLog[id, event] <<- currentTime
    handleEvent(event)
  }
}

# The runAll functions loops over the number of individuals that are to be simulated specified through the 
# 'nIndividuals' variable. It tracks and prints the progress of the simulation and updates the 'id' variable that 
# should be present in the parent environment and which is used to track the index of the individual in matrices and
# arrays etcetera. It finally returns the 'resourceLog' and 'eventLog'. At the beginning it sets the random number 
# seed based on the 'seed' variable that has to be specified in the parent environment.
runAll = function() {
  
  # Print progress
  tic <- Sys.time()
  cat(paste0('Simulation of ', format(nIndividuals, big.mark = ','), ' individuals started on ', tic, '\n\nSimulation progress:\n'))
  
  set.seed(seed)
  
  for(iIndividual in 1:nIndividuals) {
    
    # Run the DES for an individual
    id <<- iIndividual
    runOne()
    
    # Print progress
    if(id == 1 | (id %% ((floor(nIndividuals/20)*20)/20)) == 0 | id == nIndividuals) {
      progress <- id / nIndividuals * 100
      cat(sprintf('\r[%-50s] %d%%', paste(rep('=', progress/2), collapse = ''), floor(progress)), ' runtime:', round(difftime(Sys.time(), tic, units = 'mins'), 1), 'minutes')
      if(id == nIndividuals) cat("\n")
    }
    
  }
  
  # Create output and remove any other objects to handle memory pressure
  out <- list(eventLog = eventLog, resourceLog = resourceLog)
  rm(list = ls()[ls() != 'out'])
  
  return(out)
  
}



