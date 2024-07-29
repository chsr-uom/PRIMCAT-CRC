#-----------------------------------------------------------------------------------------#
# The PRIMCAT-CRC discrete event simulation model was built using the established         #
# modelling framework and extracted parameters to populate it.                            #
# The framework dictates the flow and discrete events patients can experience over time.  #
#                                                                                         #
# PRIMCAT-CRC is comprised of four functions:                                             #  
# 1)	Parameter function:                                                                 #
#   loads extracted parameters and returns them in a list that can be used to run         #
#   the simulation. The list of parameters can be defined for each time point at which    #
#   the parameters are to be updated, with the option to update parameter values over     #
#   time by passing an "updates" argument that contains a list of lists.                  #
#                                                                                         #
# 2)	Population function:                                                                #
#   defines the patient population to be simulated, including the number of               #
#   individuals (incidence) and their disease stage at diagnosis, for each year           #
#   that is to be simulated. The resulting data frame can be used to run the simulation.  #
#                                                                                         #
# 3)	Simulation function:                                                                #
#   corresponds to the actual DES that simulates individual patient trajectories.         #
#   It takes in a data frame or matrix that specifies the number of individuals and       #
#   their start time and first event, a list of parameters used in the simulation,        #
#   and a numeric value defining the time until which the simulation is to be run.        #
#                                                                                         #
# 4)	PRIMCAT-CRC:                                                                        #
#   wrapper function that calls the above functions multiple times and in parallel        #
#   to perform multiple runs of the discrete event simulation and average the outputs.    #
#-----------------------------------------------------------------------------------------#
#
#
# More details regarding the different functions are provided with the functions themselves.
# The PRIMCAT-CRC model structure is described in the manuscript and available 
# as supplementary material at the following link: https://doi.org/10.1016/j.jval.2024.06.006
#
# Version history:
# - September 2021  Koen Degeling       R v4.0.3    Initial version
# - November 2021   Koen Degeling       R v4.0.3    Checked code and improved commenting
# - February 2022   Fanny Franchini     R v4.0.5    Code review and further testing
# - January 2023    Fanny Franchini     R v4.1.2    Rerun on new R version, added comments for publication
# - June 2024       Fanny Franchini     R v4.3.1    Rerun on new R version, revamped for GitHub repository
#
#
#-----------------------------------------------------------------------------------------#
#
#
# The PRIMCAT_CRC_Parameters function loads the parameters that have been exported from the 1_parameters.R script
# and returns them in a list that can be used to run the simulation. A list of parameters is defined for each time
# point at which the parameters are to be defined/updated. The deterministic basecase values are used to define ALL
# parameters at the start time, which is specified through the "starttime" argument. These values can be updated 
# through the list that can be passed to the "updates" argument. This way, parameter value changes over time can be
# incorporated. The list that is returned will be used in the simulation to define events at the times defined by
# the names of the pars list, at which the parameters will be loaded/updated. The "updates" argument can be used, 
# for example, to define the scenarios for the horizon scanning. The "updates" arguments requires a list of list, 
# each of which should be named according to the time at which the updates are to be performed, with the list
# itself being a named list of parameters to be updated and their values (which can again be a list, for example to
# define a distribution). Note that parameters are not changed for a period of time, but from the moment at which
# they are updated for the remainder of the simulation or until the are updated again. an example of the "update"
# list is as follows:
#   list(
#   '1' = list(
#     p_C1_SG = 0.5,
#     p_C1_NEW = 0.5,
#     p_C1_NEW_recur = 0.05,
#   ),
#   '2' = list(
#     p_C1_SG = 0,
#     p_C1_NEW = 1
#   )
#   )
# This example updates the three parameters at time 1, and two of the three parameters again at time 2. The
# p_C1_NEW_recur parameter that is updated only once, remains at it's updated value of 0.05.
PRIMCAT_CRC_Parameters <- function(starttime = 0, updates = NULL) {
  
  # Load deterministic base-case parameter values as the start values
  pars <- list(readRDS(file = 'Input/pars_deterministic.RDS'))
  names(pars) <- as.character(starttime)
  
  # Adding the updates to be performed
  if(!is.null(updates)) {
    for(i in 1:length(updates)) {
      
      # Extract the name and time of the update
      updatename <- names(updates)[i]
      updatetime <- as.numeric(updatename)
      
      # Check whether the update is after the starttime
      if(updatetime < starttime) stop('Times of updates should be at or after the starttime')
      
      # Add the update to the pars list
      pars[[i+1]]        <- updates[[i]]
      names(pars)[[i+1]] <- updatename
      
    }
  }
  
  # Return the list
  return(pars)
  
}


# The PRIMCAT_CRC_Population function defines the patient population that is to be simulated. For each year that is
# to be simulated, it defines the number of individuals and their disease stage at diagnosis and at what time they
# enter the simulation. The resulting data.frame can be used to run the simulation. All the information about the
# population characteristics for the years that are to be simulated are provided through the 'm_years' argument.
# This function input m_years should be a matrix or data.frame with a row for each year for which a population is 
# to be generated.
#
# If function input by_total_incidence is TRUE, then the population is defined by applying a stage distribution 
# based on probability to the total incidence, for which the following columns are required:
# - t       the simulation time at which the year should start
# - n       the total incidence of that year
# - p_C1    the probability of being C1 in that year
# - P_C2    the probability of being C2 in that year
# - p_C3    the probability of being C3 in that year
# - p_R1    the probability of being R1 in that year
# - P_R2    the probability of being R2 in that year
# - p_R3    the probability of being R3 in that year
# - p_CR4   the probability of being CR4 in that year
#
# If function input by_total_incidence is FALSE, then the population is defined by stage-specific incidences using
# the following columns:
# - t       the simulation time at which the year should start
# - n_C1    the incidence of C1 in that year
# - n_C2    the incidence of C2 in that year
# - n_C3    the incidence of C3 in that year
# - n_R1    the incidence of R1 in that year
# - n_R2    the incidence of R2 in that year
# - n_R3    the incidence of R3 in that year
# - n_CR4   the incidence of CR4 in that year
PRIMCAT_CRC_Population <- function(m_years, by_total_incidence = TRUE, seed = 1) {
  
  # Random seed for reproducibility
  set.seed(seed)
  
  # Names of the stages that are considered
  stages <- c('C1', 'C2', 'C3', 'R1', 'R2', 'R3', 'CR4')
  
  # Loop over the years (i.e., rows) and determine the number of individuals in each stage, resulting in a list of
  # year-specific data.frames
  ls_df_start <- apply(m_years, 1, function(year) {
    
    # If the numbers are to be determined by total incidence, multiply the incidence by the stage proportions and
    # round to ensure only 'whole' patients are simulated
    if(by_total_incidence) {
      
      n_stages <- round(x = year['n'] * year[paste0('p_', stages)], digits = 0)
      
      # If the numbers per stage are already provided, these simply can be extracted and rounded to be sure that only
      # 'whole' patients are simulated
    } else {
      
      n_stages <- round(x = year[paste0('n_', stages)], digits = 0)
      
    }
    
    # Sample a random time within the year from a Uniform distribution and save the disease stage as first event
    data.frame(
      startTime  = runif(n = sum(n_stages), min = year['t'], max = year['t'] + 1),
      firstEvent = rep(x = stages, times = n_stages)
    )
    
  })
  
  # Combine the list of all the year-specific data.frames into one data.frame and remove any remaining objects to
  # handle memory efficiently
  df_start <- do.call('rbind', ls_df_start)
  rm(list = ls()[ls() != 'df_start'])
  
  return(df_start)
  
}

# The PRIMCAT_CRC_Simulation function is the actual DES that simulated the individual patient trajectories. Time in
# the simulation is defined in years It has the following input arguments:
# - start     data.frame or matrix that specifies the number of individuals (rows) and must contains a startTime 
#             column that specifies the time at which the individual enters the simulation and a firstEvent column 
#             that specifies the first event that is to be simulated for the individual. Such a data.frame is 
#             returned by the PRIMCAT_CRC_Population function
# - pars      list containing all the parameters used in the simulation, which is returned from the 
#             PRIMCAT_CRC_Parameters function
# - until     numeric defining the time until which the simulation is to be run, with the default being Inf 
#             resulting in a simulation until all events have occurred
# - seed      random number seed on for reproducibility
# 
# The function has three main parts: 
#   1. INITIALIZATION     loading the general DES functions and defining all the objects required in the simulation
#   2. EVENT HANDLING     processing all the events that define the simulation model structure
#   3. SIMULATION         running the simulation
#
# It returns the output of the runAll() function, which is a list containing both the eventLog and resourceLog, 
# which contain a individual-level record of the events occurred and resources used, respectively.
PRIMCAT_CRC_Simulation <- function(start, pars, until = Inf, seed = 1) {
  
  ## 1. INITIALIZATION ----
  
  # Save the enviroment in which the model parameters are to be loaded
  env_sim <- environment()
  
  # Loading the general DES functions:
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
  source(file = '../Input/DES_functions.R', local = TRUE)
  
  # Specifying objects (variables, functions, etc.) that have to be available from within
  # all sections/functions within the function
  queue <- priorityQueue()          # the event queue
  nIndividuals <- nrow(start)       # number of individuals to be simulated
  
  # Specifying the resources that can be seized and, subsequently, the resourceLog object
  resources <- c(
    'C1_Tx_SG', 'C1_Tx_NEW',
    'C2_Tx_SG', 'C2_Tx_SGadj', 'C2_Tx_NEW',
    'C3_Tx_SG', 'C3_Tx_SGadj', 'C3_Tx_NEW',
    'R1_Tx_SG', 'R1_Tx_neoSG', 'R1_Tx_SGadj', 'R1_Tx_neoSGadj', 'R1_Tx_NEW',
    'R2_Tx_SG', 'R2_Tx_neoSG', 'R2_Tx_SGadj', 'R2_Tx_neoSGadj', 'R2_Tx_NEW',
    'R3_Tx_SG', 'R3_Tx_neoSG', 'R3_Tx_SGadj', 'R3_Tx_neoSGadj', 'R3_Tx_NEW',
    'LRR_Tx_SG', 'LRR_Tx_neoSG', 'LRR_Tx_SGadj', 'LRR_Tx_neoSGadj', 'LRR_Tx_systemic', 'LRR_Tx_NEW',
    'CR4_testRAS',
    'RASwt_L0_SGmets', 'RASwt_L0_SGprim', 'RASwt_L0_RT', 
    'RASwt_L1_SGmets', 'RASwt_L1_SGprim', 'RASwt_L1_RT', 
    'RASwt_L2_SGmets', 'RASwt_L2_SGprim', 'RASwt_L2_RT', 
    'RASwt_L3_SGmets', 'RASwt_L3_SGprim', 'RASwt_L3_RT', 
    'RASwt_L4_SGmets', 'RASwt_L4_SGprim', 'RASwt_L4_RT', 
    'RASwt_L1_chemo', 'RASwt_L1_EGFR', 'RASwt_L1_chemoEGFR', 'RASwt_L1_chemoVEGF', 'RASwt_L1_NEW',  
    'RASwt_L2_chemo', 'RASwt_L2_EGFR', 'RASwt_L2_chemoEGFR', 'RASwt_L2_chemoVEGF', 'RASwt_L2_NEW',
    'RASwt_L3_chemo', 'RASwt_L3_NEW',
    'RASwt_L4_chemo', 'RASwt_L4_NEW',
    'RASmt_L0_SGmets', 'RASmt_L0_SGprim', 'RASmt_L0_RT', 
    'RASmt_L1_SGmets', 'RASmt_L1_SGprim', 'RASmt_L1_RT', 
    'RASmt_L2_SGmets', 'RASmt_L2_SGprim', 'RASmt_L2_RT', 
    'RASmt_L3_SGmets', 'RASmt_L3_SGprim', 'RASmt_L3_RT', 
    'RASmt_L4_SGmets', 'RASmt_L4_SGprim', 'RASmt_L4_RT',
    'RASmt_L1_chemo', 'RASmt_L1_chemoVEGF', 'RASmt_L1_NEW',
    'RASmt_L2_chemo', 'RASmt_L2_chemoVEGF', 'RASmt_L2_NEW',  
    'RASmt_L3_chemo', 'RASmt_L3_NEW',
    'RASmt_L4_chemo', 'RASmt_L4_NEW'
  )
  resourceLog <- array(dim = c(nIndividuals, length(resources), 2), dimnames = list(NULL, resources, c('start', 'end')))
  
  # Specifying the events that can occur and, hence, are tracked in the eventLog object
  events <- c(
    'C1', 'C1_Tx',
    'C2', 'C2_Tx',
    'C3', 'C3_Tx',
    'R1', 'R1_Tx',
    'R2', 'R2_Tx',
    'R3', 'R3_Tx',
    'LRR', 'LRR_Tx',
    'CR4_progr', 'CR4', 
    'RASwt', 'RASwt_L0', 'RASwt_L1', 'RASwt_L2', 'RASwt_L3', 'RASwt_L4', 
    'RASmt', 'RASmt_L0', 'RASmt_L1', 'RASmt_L2', 'RASmt_L3', 'RASmt_L4', 
    'ParameterUpdate', 'Until', 'End')
  eventLog <- array(dim = c(nIndividuals, length(events)), dimnames = list(NULL, events))
  
  # Individual level characteristics/history
  id <- NA_real_                    # numeric for the index of the individual
  currentTime <- NA_real_           # current time of the simulation
  L1 <- NA_real_                    # first-line systemic treatment
  
  
  ## 2. EVENT HANDLING ----
  
  # The goal of the init function is to ensure all global variable are reset or set to the value specified for 
  # that individual by the start input matrix/data.frame, as well as to schedule the first event and the 'Until'
  # event based on the until argument
  init = function() {
    
    # Ensure the event queue is empty
    queue$clear()
    
    # (Re)set the global variables
    #currentTime <<- as.numeric(start[id, 'startTime'])
    currentTime <<- NA_real_
    
    # Schedule the events that define the parameter 
    for(i in 1:length(pars)) schedule('ParameterUpdate', as.numeric(names(pars)[i]))
    
    # Schedule the first event based on the disease stage and the 'End' event
    schedule(start[id, 'firstEvent'], as.numeric(start[id, 'startTime']))
    schedule('Until', until)
    
  }
  
  # The handleEvent function defines what actions have to be performed when a particular event occurs. It is
  # ordered according to the different disease stages and treatment lines. A switch statement is used to identify
  # the event that needs to be processed, because this is much faster than a large if-elseif-elseif-elseif...
  # statement. The names of the events correspond to the names of the treatment lines etc. in the conceptual model.
  # Note that extensive comments are provided for the first events and for subsequent events when relevant. To
  # reduce the length of the code, comments are not repeated too much.
  handleEvent = function(event) {
    switch(EXPR = event,
           
           ## 2.1 Stage I colon (C1) ----
           
           C1 = {
             
             # Determine whether any treatment will be provided and, if so, when that will be. Note that some individuals
             # will have no time-to-treatment (TTT), i.e. a TTT of zero, in which case the treatment is scheduled for the
             # currentTime. If there is no treatment, the End event is scheduled at the currentTime and the individual
             # leaves the simulation.
             if(runif(1) < p_C1_anyTx) {
               nextEvent <- 'C1_Tx'
               nextTime  <- if(runif(1) < p_C1_noTTT) {currentTime} else {currentTime + rcustom(t_C1_TTT)}
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           C1_Tx = {
             
             # Sample the treatment to be received. For this also a switch statement is used because this is very 
             # efficient. The draw function is used, which returns an number/index that corresponds to the place of the
             # corresponding probability in the vector of probabilities that is provided as input. Based on the treatment
             # the probability of recurrence and distribution for TTR is saved so that this can be used later to determine
             # the next event.
             switch(
               EXPR = draw(c(p_C1_SG, p_C1_NEW)),
               # 1 = SG
               {
                 p_recur <- p_C1_SG_recur
                 t_recur <- t_C1_SG_TTR
                 seize('C1_Tx_SG')
               }, 
               # 2 = NEW
               {
                 p_recur <- p_C1_NEW_recur
                 t_recur <- t_C1_NEW_TTR
                 seize('C1_Tx_NEW')
               },
               stop('Index in selecting treatment for C1_Tx out of range')
             )
             
             # Determine whether, which, and when recurrence occurs based on the saved probability of recurrence and 
             # distribution for TTR. Also determine whether the recurrence will be a locoregional or distant recurrence.
             # If there is no recurrence, the End event is scheduled at the currentTime and the individual leaves the
             # simulation.
             if(runif(1) < p_recur) {
               nextEvent <- if(runif(1) < p_C1_recur_distant) {'CR4_progr'} else {'LRR'}
               nextTime  <- currentTime + rcustom(t_recur)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.2 Stage II colon (C2) ----
           
           C2 = {
             
             if(runif(1) < p_C2_anyTx) {
               nextEvent <- 'C2_Tx'
               nextTime  <- if(runif(1) < p_C2_noTTT) {currentTime} else {currentTime + rcustom(t_C2_TTT)}
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           C2_Tx = {
             
             # Sample the treatment to be received. Note that there are now options, so the draw function returns a 1 or
             # a 2 corresponding to the order in which the probabilities are provided.
             switch(
               EXPR = draw(c(p_C2_SG, p_C2_SGadj, p_C2_NEW)),
               # 1 = SG
               {
                 p_recur <- p_C2_SG_recur
                 t_recur <- t_C2_SG_TTR
                 seize('C2_Tx_SG')
               }, 
               # 2 = SGadj
               {
                 p_recur <- p_C2_SGadj_recur
                 t_recur <- t_C2_SGadj_TTR
                 seize('C2_Tx_SGadj')
               }, 
               # 3 = NEW
               {
                 p_recur <- p_C2_NEW_recur
                 t_recur <- t_C2_NEW_TTR
                 seize('C2_Tx_NEW')
               }, 
               stop('Index in selecting treatment for C2_Tx out of range')
             )
             
             # Determine whether, which, and when recurrence occurs
             if(runif(1) < p_recur) {
               nextEvent <- if(runif(1) < p_C2_recur_distant) {'CR4_progr'} else {'LRR'}
               nextTime  <- currentTime + rcustom(t_recur)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.3 Stage III colon (C3) ----
           
           C3 = {
             
             if(runif(1) < p_C3_anyTx) {
               nextEvent <- 'C3_Tx'
               nextTime  <- if(runif(1) < p_C3_noTTT) {currentTime} else {currentTime + rcustom(t_C3_TTT)}
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           C3_Tx = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_C3_SG, p_C3_SGadj, p_C3_NEW)),
               # 1 = SG
               {
                 p_recur <- p_C3_SG_recur
                 t_recur <- t_C3_SG_TTR
                 seize('C3_Tx_SG')
               }, 
               # 2 = SGadj
               {
                 p_recur <- p_C3_SGadj_recur
                 t_recur <- t_C3_SGadj_TTR
                 seize('C3_Tx_SGadj')
               }, 
               # 3 = NEW
               {
                 p_recur <- p_C3_NEW_recur
                 t_recur <- t_C3_NEW_TTR
                 seize('C3_Tx_NEW')
               }, 
               stop('Index in selecting treatment for C3_Tx out of range')
             )
             
             # Determine whether, which, and when recurrence occurs
             if(runif(1) < p_recur) {
               nextEvent <- if(runif(1) < p_C3_recur_distant) {'CR4_progr'} else {'LRR'}
               nextTime  <- currentTime + rcustom(t_recur)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.4 Stage I rectal (R1) ----
           
           R1 = {
             
             if(runif(1) < p_R1_anyTx) {
               nextEvent <- 'R1_Tx'
               nextTime  <- if(runif(1) < p_R1_noTTT) {currentTime} else {currentTime + rcustom(t_R1_TTT)}
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           R1_Tx = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_R1_SG, p_R1_neoSG, p_R1_SGadj, p_R1_neoSGadj, p_R1_NEW)),
               # 1 = SG
               {
                 p_recur <- p_R1_SG_recur
                 t_recur <- t_R1_SG_TTR
                 seize('R1_Tx_SG')
               }, 
               # 2 = neoSG
               {
                 p_recur <- p_R1_neoSG_recur
                 t_recur <- t_R1_neoSG_TTR
                 seize('R1_Tx_neoSG')
               }, 
               # 3 = SGadj
               {
                 p_recur <- p_R1_SGadj_recur
                 t_recur <- t_R1_SGadj_TTR
                 seize('R1_Tx_SGadj')
               }, 
               # 4 = neoSGadj
               {
                 p_recur <- p_R1_neoSGadj_recur
                 t_recur <- t_R1_neoSGadj_TTR
                 seize('R1_Tx_neoSGadj')
               }, 
               # 5 = NEW
               {
                 p_recur <- p_R1_NEW_recur
                 t_recur <- t_R1_NEW_TTR
                 seize('R1_Tx_NEW')
               }, 
               stop('Index in selecting treatment for R1_Tx out of range')
             )
             
             # Determine whether, which, and when recurrence occurs
             if(runif(1) < p_recur) {
               nextEvent <- if(runif(1) < p_R1_recur_distant) {'CR4_progr'} else {'LRR'}
               nextTime  <- currentTime + rcustom(t_recur)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.5 Stage II rectal (R2) ----
           
           R2 = {
             
             if(runif(1) < p_R2_anyTx) {
               nextEvent <- 'R2_Tx'
               nextTime  <- if(runif(1) < p_R2_noTTT) {currentTime} else {currentTime + rcustom(t_R2_TTT)}
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           R2_Tx = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_R2_SG, p_R2_neoSG, p_R2_SGadj, p_R2_neoSGadj, p_R2_NEW)),
               # 1 = SG
               {
                 p_recur <- p_R2_SG_recur
                 t_recur <- t_R2_SG_TTR
                 seize('R2_Tx_SG')
               }, 
               # 2 = neoSG
               {
                 p_recur <- p_R2_neoSG_recur
                 t_recur <- t_R2_neoSG_TTR
                 seize('R2_Tx_neoSG')
               }, 
               # 3 = SGadj
               {
                 p_recur <- p_R2_SGadj_recur
                 t_recur <- t_R2_SGadj_TTR
                 seize('R2_Tx_SGadj')
               }, 
               # 4 = neoSGadj
               {
                 p_recur <- p_R2_neoSGadj_recur
                 t_recur <- t_R2_neoSGadj_TTR
                 seize('R2_Tx_neoSGadj')
               }, 
               # 5 = NEW
               {
                 p_recur <- p_R2_NEW_recur
                 t_recur <- t_R2_NEW_TTR
                 seize('R2_Tx_NEW')
               }, 
               stop('Index in selecting treatment for R2_Tx out of range')
             )
             
             # Determine whether, which, and when recurrence occurs
             if(runif(1) < p_recur) {
               nextEvent <- if(runif(1) < p_R2_recur_distant) {'CR4_progr'} else {'LRR'}
               nextTime  <- currentTime + rcustom(t_recur)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.6 Stage III rectal (R3) ----
           
           R3 = {
             
             if(runif(1) < p_R3_anyTx) {
               nextEvent <- 'R3_Tx'
               nextTime  <- if(runif(1) < p_R3_noTTT) {currentTime} else {currentTime + rcustom(t_R3_TTT)}
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           R3_Tx = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_R3_SG, p_R3_neoSG, p_R3_SGadj, p_R3_neoSGadj, p_R3_NEW)),
               # 1 = SG
               {
                 p_recur <- p_R3_SG_recur
                 t_recur <- t_R3_SG_TTR
                 seize('R3_Tx_SG')
               }, 
               # 2 = neoSG
               {
                 p_recur <- p_R3_neoSG_recur
                 t_recur <- t_R3_neoSG_TTR
                 seize('R3_Tx_neoSG')
               }, 
               # 3 = SGadj
               {
                 p_recur <- p_R3_SGadj_recur
                 t_recur <- t_R3_SGadj_TTR
                 seize('R3_Tx_SGadj')
               }, 
               # 4 = neoSGadj
               {
                 p_recur <- p_R3_neoSGadj_recur
                 t_recur <- t_R3_neoSGadj_TTR
                 seize('R3_Tx_neoSGadj')
               }, 
               # 5 = NEW
               {
                 p_recur <- p_R3_NEW_recur
                 t_recur <- t_R3_NEW_TTR
                 seize('R3_Tx_NEW')
               }, 
               stop('Index in selecting treatment for R3_Tx out of range')
             )
             
             # Determine whether, which, and when recurrence occurs
             if(runif(1) < p_recur) {
               nextEvent <- if(runif(1) < p_R3_recur_distant) {'CR4_progr'} else {'LRR'}
               nextTime  <- currentTime + rcustom(t_recur)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.7 Locoregional recurrence (LRR) ----
           
           LRR = {
             
             if(runif(1) < p_LRR_anyTx) {
               nextEvent <- 'LRR_Tx'
               nextTime  <- if(runif(1) < p_LRR_noTTT) {currentTime} else {currentTime + rcustom(t_LRR_TTT)}
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           LRR_Tx = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_LRR_SG, p_LRR_neoSG, p_LRR_SGadj, p_LRR_neoSGadj, p_LRR_systemic, p_LRR_NEW)),
               # 1 = SG
               {
                 p_recur <- p_LRR_SG_recur
                 t_recur <- t_LRR_SG_TTR
                 seize('LRR_Tx_SG')
               }, 
               # 2 = neoSG
               {
                 p_recur <- p_LRR_neoSG_recur
                 t_recur <- t_LRR_neoSG_TTR
                 seize('LRR_Tx_neoSG')
               }, 
               # 3 = SGadj
               {
                 p_recur <- p_LRR_SGadj_recur
                 t_recur <- t_LRR_SGadj_TTR
                 seize('LRR_Tx_SGadj')
               }, 
               # 4 = neoSGadj
               {
                 p_recur <- p_LRR_neoSGadj_recur
                 t_recur <- t_LRR_neoSGadj_TTR
                 seize('LRR_Tx_neoSGadj')
               }, 
               # 5 = systemic
               {
                 p_recur <- p_LRR_systemic_recur
                 t_recur <- t_LRR_systemic_TTR
                 seize('LRR_Tx_systemic')
               }, 
               # 6 = NEW
               {
                 p_recur <- p_LRR_NEW_recur
                 t_recur <- t_LRR_NEW_TTR
                 seize('LRR_Tx_NEW')
               }, 
               stop('Index in selecting treatment for LRR_Tx out of range')
             )
             
             # Determine whether, which, and when recurrence occurs
             if(runif(1) < p_recur) {
               nextEvent <- if(runif(1) < p_LRR_recur_distant) {'CR4_progr'} else {'LRR'}
               nextTime  <- currentTime + rcustom(t_recur)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.8 Stage 4 colorectal (CR4) ----
           
           CR4_progr = {
             
             # The CR4_progr event is used so that at a later stage it can be separated out how many CR4 events 
             # corresponded to de novo diagnosis or patients who progressed from earlier stages.
             schedule('CR4', currentTime)
             
           },
           
           CR4 = {
             
             # For those who receive treatment, it is determined whether they are tested for their RAS status and, if so,
             # wheterh they are RASwt or RASmt, based on which the first treatment event is scheduled.
             if(runif(1) < p_CR4_anyTx) {
               nextEvent <- if((runif(1) < p_CR4_testRAS) & (runif(1) < p_CR4_RASwt)) {'RASwt'} else {'RASmt'} 
               nextTime  <- currentTime
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           
           ## 2.9 RASwt ----
           
           RASwt = {
             
             # Determine whether the individual will go to L0 or L1
             if(runif(1) < p_RASwt_L0) {
               nextEvent <- 'RASwt_L0' 
               nextTime  <- if(runif(1) < p_RASwt_L0_noTTT) {currentTime} else {currentTime + rcustom(t_RASwt_L0_TTT)} 
               schedule(nextEvent, nextTime)
             } else {
               nextEvent <- 'RASwt_L1' 
               nextTime  <- if(runif(1) < p_RASwt_L1_noTTT) {currentTime} else {currentTime + rcustom(t_RASwt_L1_TTT)} 
               schedule(nextEvent, nextTime)
             }
             
           },
           
           RASwt_L0 = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_RASwt_L0_SGmets, p_RASwt_L0_SGprim, p_RASwt_L0_RT, p_RASwt_L0_SGmetsSGprim, p_RASwt_L0_SGmetsRT)),
               # 1 = SGmets
               {
                 p_progr <- p_RASwt_L0_progr
                 t_progr <- t_RASwt_L0_TTP
                 seize('RASwt_L0_SGmets')
               }, 
               # 2 = SGprim
               {
                 p_progr <- p_RASwt_L0_progr
                 t_progr <- t_RASwt_L0_TTP
                 seize('RASwt_L0_SGprim')
               }, 
               # 3 = RT
               {
                 p_progr <- p_RASwt_L0_progr
                 t_progr <- t_RASwt_L0_TTP
                 seize('RASwt_L0_RT')
               }, 
               # 4 = SGmets + SGPrim
               {
                 p_progr <- p_RASwt_L0_progr
                 t_progr <- t_RASwt_L0_TTP
                 seize('RASwt_L0_SGmets')
                 seize('RASwt_L0_SGprim')
               },           
               # 5 = SGmets + RT
               {
                 p_progr <- p_RASwt_L0_progr
                 t_progr <- t_RASwt_L0_TTP
                 seize('RASwt_L0_SGmets')
                 seize('RASwt_L0_RT')
               }, 
               stop('Index in selecting treatment for RASwt_L0 out of range')
             )
             
             # Determine whether and when progression occurs
             # IMPORTANT! Note that at this point in the simulation, 12 weeks need to be added to the sampled times 
             # because 12 weeks were subtracted from the TTP to account for a plateau at the beginning of the curve that
             # came to existence due to the classification algorithm used to define whether surgery should be considered
             # part of first-line systemic treatment or not.
             if(runif(1) < p_progr) {
               nextEvent <- 'RASwt_L1' 
               nextTime  <- currentTime + 12/52 + rcustom(t_progr) # add 12w 
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           RASwt_L1 = {
             
             # Sample the treatment to be received. Note that we same the type of treatment in the L1 variable because
             # this is used in second-line treatment to select the treatment.
             switch(
               EXPR = draw(c(p_RASwt_L1_chemo, p_RASwt_L1_EGFR, p_RASwt_L1_chemoEGFR, p_RASwt_L1_chemoVEGF, p_RASwt_L1_NEW)),
               # 1 = chemo
               {
                 L1 <<- 'chemo'
                 p_progr <- p_RASwt_L1_chemo_progr
                 t_progr <- t_RASwt_L1_chemo_TTP
                 seize('RASwt_L1_chemo')
               }, 
               # 2 = EGFR
               {
                 L1 <<- 'EGFR'
                 p_progr <- p_RASwt_L1_EGFR_progr
                 t_progr <- t_RASwt_L1_EGFR_TTP
                 seize('RASwt_L1_EGFR')
               }, 
               # 3 = chemoEGFR
               {
                 L1 <<- 'chemoEGFR'
                 p_progr <- p_RASwt_L1_chemoEGFR_progr
                 t_progr <- t_RASwt_L1_chemoEGFR_TTP
                 seize('RASwt_L1_chemoEGFR')
               }, 
               # 4 = chemoVEGF
               {
                 L1 <<- 'chemoVEGF'
                 p_progr <- p_RASwt_L1_chemoVEGF_progr
                 t_progr <- t_RASwt_L1_chemoVEGF_TTP
                 seize('RASwt_L1_chemoVEGF')
               }, 
               # 5 = NEW
               {
                 L1 <<- 'NEW'
                 p_progr <- p_RASwt_L1_NEW_progr
                 t_progr <- t_RASwt_L1_NEW_TTP
                 seize('RASwt_L1_NEW')
               }, 
               stop('Index in selecting treatment for RASwt_L1 out of range')
             )
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASwt_L1_SGmets) seize('RASwt_L1_SGmets')
             if(runif(1) < p_RASwt_L1_SGprim) seize('RASwt_L1_SGprim')
             if(runif(1) < p_RASwt_L1_RT)     seize('RASwt_L1_RT')
             
             # Determine whether and when progression occurs
             if(runif(1) < p_progr) {
               nextEvent <- 'RASwt_L2' 
               nextTime  <- currentTime + rcustom(t_progr)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           RASwt_L2 = {
             
             # Sample the treatment to be received. Note that this is conditional on the type of treatment received as
             # first-line systemic treatment.
             if(L1 %in% c('chemo')) {
               switch(
                 EXPR = draw(c(p_RASwt_L1chemo_L2_chemo, p_RASwt_L1chemo_L2_EGFR, p_RASwt_L1chemo_L2_chemoEGFR, p_RASwt_L1chemo_L2_chemoVEGF, p_RASwt_L1chemo_L2_NEW)),
                 # 1 = chemo
                 {
                   p_progr <- p_RASwt_L2_chemo_progr
                   t_progr <- t_RASwt_L2_chemo_TTP
                   seize('RASwt_L2_chemo')
                 }, 
                 # 2 = EGFR
                 {
                   p_progr <- p_RASwt_L2_EGFR_progr
                   t_progr <- t_RASwt_L2_EGFR_TTP
                   seize('RASwt_L2_EGFR')
                 }, 
                 # 3 = chemoEGFR
                 {
                   p_progr <- p_RASwt_L2_chemoEGFR_progr
                   t_progr <- t_RASwt_L2_chemoEGFR_TTP
                   seize('RASwt_L2_chemoEGFR')
                 }, 
                 # 4 = chemoVEGF
                 {
                   L2 <<- 'chemoVEGF'
                   p_progr <- p_RASwt_L2_chemoVEGF_progr
                   t_progr <- t_RASwt_L2_chemoVEGF_TTP
                   seize('RASwt_L2_chemoVEGF')
                 }, 
                 # 5 = NEW
                 {
                   L2 <<- 'NEW'
                   p_progr <- p_RASwt_L2_NEW_progr
                   t_progr <- t_RASwt_L2_NEW_TTP
                   seize('RASwt_L2_NEW')
                 }, 
                 stop('Index in selecting treatment for RASwt_L2 out of range')
               )
             } else if(L1 %in% c('EGFR', 'chemoEGFR')) {
               switch(
                 EXPR = draw(c(p_RASwt_L1EGFR_L2_chemo, p_RASwt_L1EGFR_L2_chemoVEGF, p_RASwt_L1EGFR_L2_NEW)),
                 # 1 = chemo
                 {
                   p_progr <- p_RASwt_L2_chemo_progr
                   t_progr <- t_RASwt_L2_chemo_TTP
                   seize('RASwt_L2_chemo')
                 }, 
                 # 2 = chemoVEGF
                 {
                   p_progr <- p_RASwt_L2_chemoVEGF_progr
                   t_progr <- t_RASwt_L2_chemoVEGF_TTP
                   seize('RASwt_L2_chemoVEGF')
                 }, 
                 # 3 = NEW
                 {
                   p_progr <- p_RASwt_L2_NEW_progr
                   t_progr <- t_RASwt_L2_NEW_TTP
                   seize('RASwt_L2_NEW')
                 }, 
                 stop('Index in selecting treatment for RASwt_L2 out of range')
               )
             } else if(L1 %in% c('chemoVEGF')) {
               switch(
                 EXPR = draw(c(p_RASwt_L1VEGF_L2_chemo, p_RASwt_L1VEGF_L2_EGFR, p_RASwt_L1VEGF_L2_chemoEGFR, p_RASwt_L1VEGF_L2_NEW)),
                 # 1 = chemo
                 {
                   p_progr <- p_RASwt_L2_chemo_progr
                   t_progr <- t_RASwt_L2_chemo_TTP
                   seize('RASwt_L2_chemo')
                 }, 
                 # 2 = EGFR
                 {
                   p_progr <- p_RASwt_L2_EGFR_progr
                   t_progr <- t_RASwt_L2_EGFR_TTP
                   seize('RASwt_L2_EGFR')
                 }, 
                 # 3 = chemoEGFR
                 {
                   p_progr <- p_RASwt_L2_chemoEGFR_progr
                   t_progr <- t_RASwt_L2_chemoEGFR_TTP
                   seize('RASwt_L2_chemoEGFR')
                 }, 
                 # 4 = NEW
                 {
                   p_progr <- p_RASwt_L2_NEW_progr
                   t_progr <- t_RASwt_L2_NEW_TTP
                   seize('RASwt_L2_NEW')
                 }, 
                 stop('Index in selecting treatment for RASwt_L2 out of range')
               )
             } else if(L1 %in% c('NEW')) {
               switch(
                 EXPR = draw(c(p_RASwt_L1NEW_L2_chemo, p_RASwt_L1NEW_L2_EGFR, p_RASwt_L1NEW_L2_chemoEGFR, p_RASwt_L1NEW_L2_chemoVEGF, p_RASwt_L1NEW_L2_NEW)),
                 # 1 = chemo
                 {
                   p_progr <- p_RASwt_L2_chemo_progr
                   t_progr <- t_RASwt_L2_chemo_TTP
                   seize('RASwt_L2_chemo')
                 }, 
                 # 2 = EGFR
                 {
                   p_progr <- p_RASwt_L2_EGFR_progr
                   t_progr <- t_RASwt_L2_EGFR_TTP
                   seize('RASwt_L2_EGFR')
                 }, 
                 # 3 = chemoEGFR
                 {
                   p_progr <- p_RASwt_L2_chemoEGFR_progr
                   t_progr <- t_RASwt_L2_chemoEGFR_TTP
                   seize('RASwt_L2_chemoEGFR')
                 }, 
                 # 4 = chemoVEGF
                 {
                   L2 <<- 'chemoVEGF'
                   p_progr <- p_RASwt_L2_chemoVEGF_progr
                   t_progr <- t_RASwt_L2_chemoVEGF_TTP
                   seize('RASwt_L2_chemoVEGF')
                 }, 
                 # 5 = NEW
                 {
                   L2 <<- 'NEW'
                   p_progr <- p_RASwt_L2_NEW_progr
                   t_progr <- t_RASwt_L2_NEW_TTP
                   seize('RASwt_L2_NEW')
                 }, 
                 stop('Index in selecting treatment for RASwt_L2 out of range')
               )
             } else {
               stop('Recorded value for L1 not supported')
             }
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASwt_L2_SGmets) seize('RASwt_L2_SGmets')
             if(runif(1) < p_RASwt_L2_SGprim) seize('RASwt_L2_SGprim')
             if(runif(1) < p_RASwt_L2_RT)     seize('RASwt_L2_RT')
             
             # Determine whether and when progression occurs
             if(runif(1) < p_progr) {
               nextEvent <- 'RASwt_L3' 
               nextTime  <- currentTime + rcustom(t_progr)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
             
           },
           
           RASwt_L3 = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_RASwt_L3_chemo, p_RASwt_L3_NEW)),
               # 1 = chemo
               {
                 p_progr <- p_RASwt_L3_chemo_progr
                 t_progr <- t_RASwt_L3_chemo_TTP
                 seize('RASwt_L3_chemo')
               },
               # 2 = NEW
               {
                 p_progr <- p_RASwt_L3_NEW_progr
                 t_progr <- t_RASwt_L3_NEW_TTP
                 seize('RASwt_L3_NEW')
               },
               stop('Index in selecting treatment for RASwt_L3 out of range')
             )
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASwt_L3_RT) seize('RASwt_L3_RT')
             
             # Determine whether and when progression occurs
             if(runif(1) < p_progr) {
               nextEvent <- 'RASwt_L4' 
               nextTime  <- currentTime + rcustom(t_progr)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           RASwt_L4 = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_RASwt_L4_chemo, p_RASwt_L4_NEW)),
               # 1 = chemo
               {
                 seize('RASwt_L4_chemo')
               },
               # 2 = NEW
               {
                 seize('RASwt_L4_NEW')
               },
               stop('Index in selecting treatment for RASwt_L4 out of range')
             )
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASwt_L4_RT) seize('RASwt_L4_RT')
             
             # No further treatments are considered, so the End event is scheduled
             schedule('End', currentTime)
             
           },
           
           
           ## 2.10 RASmt ----
           
           RASmt = {
             
             if(runif(1) < p_RASmt_L0) {
               nextEvent <- 'RASmt_L0' 
               nextTime  <- if(runif(1) < p_RASmt_L0_noTTT) {currentTime} else {currentTime + rcustom(t_RASmt_L0_TTT)} 
               schedule(nextEvent, nextTime)
             } else {
               nextEvent <- 'RASmt_L1' 
               nextTime  <- if(runif(1) < p_RASmt_L1_noTTT) {currentTime} else {currentTime + rcustom(t_RASmt_L1_TTT)} 
               schedule(nextEvent, nextTime)
             }
             
           },
           
           RASmt_L0 = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_RASmt_L0_SGmets, p_RASmt_L0_SGprim, p_RASmt_L0_RT, p_RASmt_L0_SGmetsSGprim, p_RASmt_L0_SGmetsRT)),
               # 1 = SGmets
               {
                 p_progr <- p_RASmt_L0_progr
                 t_progr <- t_RASmt_L0_TTP
                 seize('RASmt_L0_SGmets')
               }, 
               # 2 = SGprim
               {
                 p_progr <- p_RASmt_L0_progr
                 t_progr <- t_RASmt_L0_TTP
                 seize('RASmt_L0_SGprim')
               }, 
               # 3 = RT
               {
                 p_progr <- p_RASmt_L0_progr
                 t_progr <- t_RASmt_L0_TTP
                 seize('RASmt_L0_RT')
               }, 
               # 4 = SGmets + SGPrim
               {
                 p_progr <- p_RASmt_L0_progr
                 t_progr <- t_RASmt_L0_TTP
                 seize('RASmt_L0_SGmets')
                 seize('RASmt_L0_SGprim')
               },           
               # 5 = SGmets + RT
               {
                 p_progr <- p_RASmt_L0_progr
                 t_progr <- t_RASmt_L0_TTP
                 seize('RASmt_L0_SGmets')
                 seize('RASmt_L0_RT')
               }, 
               stop('Index in selecting treatment for RASmt_L0 out of range')
             )
             
             # Determine whether and when progression occurs
             # IMPORTANT! Note that at this point in the simulation, 12 weeks need to be added to the sampled times 
             # because 12 weeks were subtracted from the TTP to account for a plateau at the beginning of the curve that
             # came to existence due to the classification algorithm used to define whether surgery should be considered
             # part of first-line systemic treatment or not.
             if(runif(1) < p_progr) {
               nextEvent <- 'RASmt_L1' 
               nextTime  <- currentTime + 12/52 + rcustom(t_progr) # add 12w 
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           RASmt_L1 = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_RASmt_L1_chemo, p_RASmt_L1_chemoVEGF, p_RASmt_L1_NEW)),
               # 1 = chemo
               {
                 L1 <<- 'chemo'
                 p_progr <- p_RASmt_L1_chemo_progr
                 t_progr <- t_RASmt_L1_chemo_TTP
                 seize('RASmt_L1_chemo')
               }, 
               # 2 = chemoVEGF
               {
                 L1 <<- 'chemoVEGF'
                 p_progr <- p_RASmt_L1_chemoVEGF_progr
                 t_progr <- t_RASmt_L1_chemoVEGF_TTP
                 seize('RASmt_L1_chemoVEGF')
               }, 
               # 3 = NEW
               {
                 L1 <<- 'NEW'
                 p_progr <- p_RASmt_L1_NEW_progr
                 t_progr <- t_RASmt_L1_NEW_TTP
                 seize('RASmt_L1_NEW')
               }, 
               stop('Index in selecting treatment for RASmt_L1 out of range')
             )
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASmt_L1_SGmets) seize('RASmt_L1_SGmets')
             if(runif(1) < p_RASmt_L1_SGprim) seize('RASmt_L1_SGprim')
             if(runif(1) < p_RASmt_L1_RT)     seize('RASmt_L1_RT')
             
             # Determine whether and when progression occurs
             if(runif(1) < p_progr) {
               nextEvent <- 'RASmt_L2' 
               nextTime  <- currentTime + rcustom(t_progr)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           RASmt_L2 = {
             
             # Sample the treatment to be received
             if(L1 %in% c('chemo')) {
               switch(
                 EXPR = draw(c(p_RASmt_L1chemo_L2_chemo, p_RASmt_L1chemo_L2_chemoVEGF, p_RASmt_L1chemo_L2_NEW)),
                 # 1 = chemo
                 {
                   p_progr <- p_RASmt_L2_chemo_progr
                   t_progr <- t_RASmt_L2_chemo_TTP
                   seize('RASmt_L2_chemo')
                 }, 
                 # 2 = chemoVEGF
                 {
                   L2 <<- 'chemoVEGF'
                   p_progr <- p_RASmt_L2_chemoVEGF_progr
                   t_progr <- t_RASmt_L2_chemoVEGF_TTP
                   seize('RASmt_L2_chemoVEGF')
                 }, 
                 # 3 = NEW
                 {
                   L2 <<- 'NEW'
                   p_progr <- p_RASmt_L2_NEW_progr
                   t_progr <- t_RASmt_L2_NEW_TTP
                   seize('RASmt_L2_NEW')
                 }, 
                 stop('Index in selecting treatment for RASmt_L2 out of range')
               )
             } else if(L1 %in% c('chemoVEGF')) {
               switch(
                 EXPR = draw(c(p_RASmt_L1VEGF_L2_chemo, p_RASmt_L1VEGF_L2_NEW)),
                 # 1 = chemo
                 {
                   p_progr <- p_RASmt_L2_chemo_progr
                   t_progr <- t_RASmt_L2_chemo_TTP
                   seize('RASmt_L2_chemo')
                 }, 
                 # 2 = NEW
                 {
                   p_progr <- p_RASmt_L2_NEW_progr
                   t_progr <- t_RASmt_L2_NEW_TTP
                   seize('RASmt_L2_NEW')
                 }, 
                 stop('Index in selecting treatment for RASmt_L2 out of range')
               )
             } else if(L1 %in% c('NEW')) {
               switch(
                 EXPR = draw(c(p_RASmt_L1NEW_L2_chemo, p_RASmt_L1NEW_L2_chemoVEGF, p_RASmt_L1NEW_L2_NEW)),
                 # 1 = chemo
                 {
                   p_progr <- p_RASmt_L2_chemo_progr
                   t_progr <- t_RASmt_L2_chemo_TTP
                   seize('RASmt_L2_chemo')
                 }, 
                 # 2 = chemoVEGF
                 {
                   L2 <<- 'chemoVEGF'
                   p_progr <- p_RASmt_L2_chemoVEGF_progr
                   t_progr <- t_RASmt_L2_chemoVEGF_TTP
                   seize('RASmt_L2_chemoVEGF')
                 }, 
                 # 3 = NEW
                 {
                   L2 <<- 'NEW'
                   p_progr <- p_RASmt_L2_NEW_progr
                   t_progr <- t_RASmt_L2_NEW_TTP
                   seize('RASmt_L2_NEW')
                 }, 
                 stop('Index in selecting treatment for RASmt_L2 out of range')
               )
             } else {
               stop('Recorded value for L1 not supported')
             }
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASmt_L2_SGmets) seize('RASmt_L2_SGmets')
             if(runif(1) < p_RASmt_L2_SGprim) seize('RASmt_L2_SGprim')
             if(runif(1) < p_RASmt_L2_RT)     seize('RASmt_L2_RT')
             
             # Determine whether and when progression occurs
             if(runif(1) < p_progr) {
               nextEvent <- 'RASmt_L3' 
               nextTime  <- currentTime + rcustom(t_progr)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
             
           },
           
           RASmt_L3 = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_RASmt_L3_chemo, p_RASmt_L3_NEW)),
               # 1 = chemo
               {
                 p_progr <- p_RASmt_L3_chemo_progr
                 t_progr <- t_RASmt_L3_chemo_TTP
                 seize('RASmt_L3_chemo')
               },
               # 2 = NEW
               {
                 p_progr <- p_RASmt_L3_NEW_progr
                 t_progr <- t_RASmt_L3_NEW_TTP
                 seize('RASmt_L3_NEW')
               },
               stop('Index in selecting treatment for RASmt_L3 out of range')
             )
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASmt_L3_RT) seize('RASmt_L3_RT')
             
             # Determine whether and when progression occurs
             if(runif(1) < p_progr) {
               nextEvent <- 'RASmt_L4' 
               nextTime  <- currentTime + rcustom(t_progr)
               schedule(nextEvent, nextTime)
             } else {
               schedule('End', currentTime)
             }
             
           },
           
           RASmt_L4 = {
             
             # Sample the treatment to be received
             switch(
               EXPR = draw(c(p_RASmt_L4_chemo, p_RASmt_L4_NEW)),
               # 1 = chemo
               {
                 seize('RASmt_L4_chemo')
               },
               # 2 = NEW
               {
                 seize('RASmt_L4_NEW')
               },
               stop('Index in selecting treatment for RASmt_L4 out of range')
             )
             
             # Sample the local treatment(s) to be received
             if(runif(1) < p_RASmt_L4_RT) seize('RASmt_L4_RT')
             
             # Determine whether and when progression occurs
             schedule('End', currentTime)
             
           },
           
           
           ## 2.11 Update the parameter values ----
           
           ParameterUpdate = {
             
             list2env(pars[[as.character(currentTime)]], envir = env_sim)
             
           },
           
           
           ## 2.12 Until or End ----
           
           # Remove all remaining events from the queue if the simulation is to end
           End = ,
           Until = {
             queue$clear()
           },
           
           stop('Unmatched event in handleEvent function: ', event)
           
    )
    
  } # end of handleEvent
  
  
  ## 3. SIMULATION ----
  
  # After loading the parameters and functions, defining all the global and individual-level variables, and
  # specifying how the events should be handled, the simulation can be started and the results returned.
  out <- runAll()
  
  return(out)
  
}


# The PRIMCAT_CRC function runs the simulation multiple times and summarizes and averages the number of events and
# resource utilization in certain time periods. It has the following input arguments:
# - n_runs                the number of runs to be performed, each using a different random number seed
# - n_nodes               the number of CPU cores to be used to perform the runs of the simulation in parallel
# - m_years               see the PRIMCAT_CRC_Population function for more information
# - timegrid              numerical vector defining the breaks of the time periods for which the events and use of
#                         resources is to be summarized. If NULL (default), year quarters are used corresponding to
#                         four time periods per year
# - gridnames             vector of names for the time periods specified through the timegrid argument. If NULL
#                         (default), the quarters are named as Q1, Q2 etc for each year
# - by_total_incidence    see the PRIMCAT_CRC_Population function for more information
# - until                 see the PRIMCAT_CRC_Simulation function for more information
#
# Note that this function uses a Fork-cluster to perform the runs in parallel, which does not work on Windows
# machines. If you are using a Windows machine, you will need to adapt the code to use a SOCK-cluster instead.
PRIMCAT_CRC <- function(n_runs, n_nodes, m_years, starttime = 0, updates = NULL, timegrid = NULL, gridnames = NULL, by_total_incidence = TRUE, until = Inf) {
  
  # Defining a standard timegrid and gridnames if not specified
  if(is.null(timegrid))  timegrid  <- seq(from = min(m_years[, 't']), to = max(m_years[, 't']) + 1, by = 0.25)
  if(is.null(gridnames)) gridnames <- paste(rep(x = c(min(m_years[, 't']):max(m_years[, 't'])), each = 4), rep(x = c('Q1', 'Q2', 'Q3', 'Q4'), times = nrow(m_years)))
  
  # Setting up the Fork-cluster to perform the runs in parallel
  library(parallel)
  cl <- makeForkCluster(nnodes = n_nodes)
  
  ls_n_events <- parLapply(cl, 1:n_runs, function(i_run) {
    
    # Get the required objects and run the simulation, with the random seed set to the run number, and removing any
    # objects that are not longer needed to control memory usage
    pars   <- PRIMCAT_CRC_Parameters(starttime = starttime, updates = updates)
    start  <- PRIMCAT_CRC_Population(m_years = m_years, by_total_incidence = by_total_incidence, seed = i_run)
    sim    <- PRIMCAT_CRC_Simulation(start = start, pars = pars, until = until, seed = i_run)
    rm(list = ls()[ls() != 'sim'])
    
    # Combine the eventLog and resourceLog into one matrix, using only the start time for the resourceLog given 
    # that duration of resource use is not included in the current model anyway, and removing any objects that are
    # not longer needed to control memory usage
    events <- cbind(sim$eventLog, sim$resourceLog[ , , 'start'])
    rm(list = ls()[ls() != 'events'])
    
    # Counting the number of events and resource uses in each time period, and removing any objects that are not 
    # longer needed to control memory usage
    n_events <- sapply(colnames(events), function(col) {
      bin <- cut(x = events[, col], breaks = timegrid, labels = gridnames)
      table(bin)
    })
    rm(list = ls()[ls() != 'n_events'])
    
    # Returning the count of the events and resources for each time period
    n_events
    
  })
  
  # Stop the Fork-cluster to control memory usage
  stopCluster(cl)
  
  # Averaging the list of matrices with the count of events and resource for each run into one matrix
  n_events <- Reduce("+", ls_n_events) / length(ls_n_events)
  
  # Creating the output object with the names of the time periods
  df_events <- data.frame(time_period = rownames(n_events))
  df_events <- cbind(df_events, n_events)
  rownames(df_events) <- NULL
  
  # removing any objects that are not longer needed to control memory usage
  rm(list = ls()[ls() != 'df_events'])
  
  return(df_events)
  
}



