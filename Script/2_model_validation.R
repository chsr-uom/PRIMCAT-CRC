#-----------------------------------------------------------------------------------------#
# The DES model was subjected to a validation process whereby the input parameters were   #
# compared against the outcomes, events, and resource usage of the simulation.            #
#                                                                                         #
# The correctness of the implementation of time-dependent parameters and                  #
# new treatment parameters was also verified, with separate verification conducted        #
# for each disease stage.                                                                 #
#                                                                                         #
# To confirm the accuracy of the simulation scenario, the new treatment was not used      #
# in the first year, followed by some use in the second year, and subsequently being the  #
# sole treatment used from the third year onward.                                         #
#                                                                                         #
# The verification process was conducted to ensure that the new treatment was             #
# correctly implemented and aligned with the intended simulation scenario.                #
#-----------------------------------------------------------------------------------------#
#
#
# This script perform the verification of the model by comparing the simulation model 
# input parameters with the simulated events, outcomes and resource use. 
#
# The PRIMCAT-CRC model structure is described in the manuscript and available 
# as supplementary material at the following link: https://doi.org/10.1016/j.jval.2024.06.006
#
# Version history:
# - November 2021   Koen Degeling       R v4.0.3    Initial version
# - February 2022   Fanny Franchini     R v4.0.5    Code review, Rerun verification with updated R
# - June 2024       Fanny Franchini     R v4.3.1    Rerun on new R, revamped for GitHub repository
#
#
#-----------------------------------------------------------------------------------------#
#
#
## 1. INITIALIZATION ----

# Clearing the global environment and console
rm(list = ls()); gc();

# Libraries
library(flexsurv)   # v2.0
library(ggplot2)    # v3.3.3

# Load the PRIMCAT-CRC model
source(file = 'Script/1_PRIMCAT_CRC_functions.R')

# Path to the plots that were exported as part of the data analysis
path_plots <- '../../WP1_treatment_patterns/PRIMCAT_CRC/plots/'

# Loading the parameters
pars  <- PRIMCAT_CRC_Parameters()
pars0 <- pars[['0']]


## 3. COLON I ----

n_sim <- 2*10^5
start_colonI <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'C1'
)

sim_colonI <- PRIMCAT_CRC_Simulation(start = start_colonI, pars = pars)

events_colonI <- as.data.frame(sim_colonI$eventLog)
resources_colonI <- as.data.frame(sim_colonI$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_colonI, 2, function(x) sum(!is.na(x)))


# Check time to treatment (TTT)
pars0$p_C1_noTTT
mean(events_colonI$C1_Tx == 0)

fit_C1_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_colonI.RDS'))
km_C1_TTT  <- survfit(formula = Surv(C1_Tx) ~ 1, data = events_colonI, subset = (events_colonI$C1_Tx > 0))

fit_C1_TTT + 
  geom_line(mapping = aes(x = km_C1_TTT$time, y = km_C1_TTT$surv), color = 'red')


# Check treatment utilization
mean(!is.na(resources_colonI$C1_Tx_SG))


# Check time to outcome (TTO)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
pars0$t_C1_SG_TTR

C1_TTO <- ifelse(!is.na(events_colonI$LRR), events_colonI$LRR - events_colonI$C1_Tx, events_colonI$CR4_progr - events_colonI$C1_Tx)
plot(density(na.omit(C1_TTO)), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.300364, scale = 0.1258425))


# Check probability of locoregional recurrence
(p_C1_LRR <- pars0$p_C1_SG_recur * (1 - pars0$p_C1_recur_distant))
mean(!is.na(events_colonI$LRR))


# Check probability of distant recurrence, potentially also through LRR
p_C1_CR4 <- pars0$p_C1_SG_recur * pars0$p_C1_recur_distant
p_C1_LRR_CR4 <- p_C1_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SG_recur) * pars0$p_LRR_recur_distant
p_C1_CR4 + p_C1_LRR_CR4
mean(!is.na(events_colonI$CR4_progr))




## 3. COLON II ----

n_sim <- 2*10^5
start_colonII <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'C2'
)

sim_colonII <- PRIMCAT_CRC_Simulation(start = start_colonII, pars = pars)

events_colonII <- as.data.frame(sim_colonII$eventLog)
resources_colonII <- as.data.frame(sim_colonII$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_colonII, 2, function(x) sum(!is.na(x)))


# Check time to treatment (TTT)
pars0$p_C2_noTTT
mean(events_colonII$C2_Tx == 0)

fit_C2_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_colonII.RDS'))
km_C2_TTT  <- survfit(formula = Surv(C2_Tx) ~ 1, data = events_colonII, subset = (events_colonII$C2_Tx > 0))

fit_C2_TTT + 
  geom_line(mapping = aes(x = km_C2_TTT$time, y = km_C2_TTT$surv), color = 'red')


# Check treatment utilization
pars0$p_C2_SG
mean(!is.na(resources_colonII$C2_Tx_SG))

pars0$p_C2_SGadj
mean(!is.na(resources_colonII$C2_Tx_SGadj))


# Check time to outcome (TTO)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
pars0$t_C2_SG_TTR

C2_TTO <- ifelse(!is.na(events_colonII$LRR), events_colonII$LRR - events_colonII$C2_Tx, events_colonII$CR4_progr - events_colonII$C2_Tx)
plot(density(na.omit(C2_TTO)), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.127459, scale = 0.24097))


# Check probability of locoregional recurrence
(p_C2_LRR <- pars0$p_C2_SG_recur * (1 - pars0$p_C2_recur_distant))
mean(!is.na(events_colonII$LRR))


# Check probability of distant recurrence, potentially also through LRR
p_C2_CR4 <- pars0$p_C2_SG_recur * pars0$p_C2_recur_distant
p_C2_LRR_CR4 <- p_C2_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SG_recur) * pars0$p_LRR_recur_distant
p_C2_CR4 + p_C2_LRR_CR4
mean(!is.na(events_colonII$CR4_progr))




## 4. COLON III ----

n_sim <- 2*10^5
start_colonIII <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'C3'
)

sim_colonIII <- PRIMCAT_CRC_Simulation(start = start_colonIII, pars = pars)

events_colonIII <- as.data.frame(sim_colonIII$eventLog)
resources_colonIII <- as.data.frame(sim_colonIII$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_colonIII, 2, function(x) sum(!is.na(x)))


# Check time to treatment (TTT)
pars0$p_C3_noTTT
mean(events_colonIII$C3_Tx == 0)

fit_C3_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_colonIII.RDS'))
km_C3_TTT  <- survfit(formula = Surv(C3_Tx) ~ 1, data = events_colonIII, subset = (events_colonIII$C3_Tx > 0))

fit_C3_TTT + 
  geom_line(mapping = aes(x = km_C3_TTT$time, y = km_C3_TTT$surv), color = 'red')


# Check treatment utilization
pars0$p_C3_SG
mean(!is.na(resources_colonIII$C3_Tx_SG))

pars0$p_C3_SGadj
mean(!is.na(resources_colonIII$C3_Tx_SGadj))


# Check time to outcome (TTO)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
C3_TTO <- ifelse(!is.na(events_colonIII$LRR), events_colonIII$LRR - events_colonIII$C3_Tx, events_colonIII$CR4_progr - events_colonIII$C3_Tx)

pars0$t_C3_SG_TTR
plot(density(na.omit(C3_TTO[!is.na(resources_colonIII$C3_Tx_SG)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.358801, scale = 0.5355751))

pars0$t_C3_SGadj_TTR
plot(density(na.omit(C3_TTO[!is.na(resources_colonIII$C3_Tx_SGadj)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.358801, scale = 0.3309068))


# Check probability of locoregional recurrence
(p_C3_LRR <- pars0$p_C3_SG_recur * (1 - pars0$p_C3_recur_distant))
mean(!is.na(events_colonIII$LRR[!is.na(resources_colonIII$C3_Tx_SG)]))

(p_C3_LRR <- pars0$p_C3_SGadj_recur * (1 - pars0$p_C3_recur_distant))
mean(!is.na(events_colonIII$LRR[!is.na(resources_colonIII$C3_Tx_SGadj)]))


# Check probability of distant recurrence, potentially also through LRR
p_C3_CR4 <- pars0$p_C3_SG_recur * pars0$p_C3_recur_distant
p_C3_LRR_CR4 <- p_C3_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SG_recur) * pars0$p_LRR_recur_distant
p_C3_CR4 + p_C3_LRR_CR4
mean(!is.na(events_colonIII$CR4_progr[!is.na(resources_colonIII$C3_Tx_SG)]))

p_C3_CR4 <- pars0$p_C3_SGadj_recur * pars0$p_C3_recur_distant
p_C3_LRR_CR4 <- p_C3_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SGadj_recur) * pars0$p_LRR_recur_distant
p_C3_CR4 + p_C3_LRR_CR4
mean(!is.na(events_colonIII$CR4_progr[!is.na(resources_colonIII$C3_Tx_SGadj)]))




## 5. RECTAL I ----

n_sim <- 2*10^5
start_rectalI <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'R1'
)

sim_rectalI <- PRIMCAT_CRC_Simulation(start = start_rectalI, pars = pars)

events_rectalI <- as.data.frame(sim_rectalI$eventLog)
resources_rectalI <- as.data.frame(sim_rectalI$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_rectalI, 2, function(x) sum(!is.na(x)))


# Check time to treatment (TTT)
pars0$p_R1_noTTT
mean(events_rectalI$R1_Tx == 0)

fit_R1_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_rectalI.RDS'))
km_R1_TTT  <- survfit(formula = Surv(R1_Tx) ~ 1, data = events_rectalI, subset = (events_rectalI$R1_Tx > 0))

fit_R1_TTT + 
  geom_line(mapping = aes(x = km_R1_TTT$time, y = km_R1_TTT$surv), color = 'red')


# Check treatment utilization
pars0$p_R1_neoSG
mean(!is.na(resources_rectalI$R1_Tx_neoSG))

pars0$p_R1_SG
mean(!is.na(resources_rectalI$R1_Tx_SG))

pars0$p_R1_SGadj
mean(!is.na(resources_rectalI$R1_Tx_SGadj))

pars0$p_R1_neoSGadj
mean(!is.na(resources_rectalI$R1_Tx_neoSGadj))


# Check time to outcome (TTO)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
pars0$t_R1_SG_TTR

R1_TTO <- ifelse(!is.na(events_rectalI$LRR), events_rectalI$LRR - events_rectalI$R1_Tx, events_rectalI$CR4_progr - events_rectalI$R1_Tx)
plot(density(na.omit(R1_TTO)), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.493804, scale = 0.1905196))


# Check probability of locoregional recurrence
(p_R1_LRR <- pars0$p_R1_SG_recur * (1 - pars0$p_R1_recur_distant))
mean(!is.na(events_rectalI$LRR))


# Check probability of distant recurrence, potentially also through LRR
p_R1_CR4 <- pars0$p_R1_SG_recur * pars0$p_R1_recur_distant
p_R1_LRR_CR4 <- p_R1_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SG_recur) * pars0$p_LRR_recur_distant
p_R1_CR4 + p_R1_LRR_CR4
mean(!is.na(events_rectalI$CR4_progr))




## 6. RECTAL II ----

n_sim <- 2*10^5
start_rectalII <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'R2'
)

sim_rectalII <- PRIMCAT_CRC_Simulation(start = start_rectalII, pars = pars)

events_rectalII <- as.data.frame(sim_rectalII$eventLog)
resources_rectalII <- as.data.frame(sim_rectalII$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_rectalII, 2, function(x) sum(!is.na(x)))


# Check time to treatment (TTT)
pars0$p_R2_noTTT
mean(events_rectalII$R2_Tx == 0)

fit_R2_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_rectalII.RDS'))
km_R2_TTT  <- survfit(formula = Surv(R2_Tx) ~ 1, data = events_rectalII, subset = (events_rectalII$R2_Tx > 0))

fit_R2_TTT + 
  geom_line(mapping = aes(x = km_R2_TTT$time, y = km_R2_TTT$surv), color = 'red')


# Check treatment utilization
pars0$p_R2_neoSG
mean(!is.na(resources_rectalII$R2_Tx_neoSG))

pars0$p_R2_SG
mean(!is.na(resources_rectalII$R2_Tx_SG))

pars0$p_R2_SGadj
mean(!is.na(resources_rectalII$R2_Tx_SGadj))

pars0$p_R2_neoSGadj
mean(!is.na(resources_rectalII$R2_Tx_neoSGadj))


# Check time to outcome (TTO)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
pars0$t_R2_SG_TTR

R2_TTO <- ifelse(!is.na(events_rectalII$LRR), events_rectalII$LRR - events_rectalII$R2_Tx, events_rectalII$CR4_progr - events_rectalII$R2_Tx)
plot(density(na.omit(R2_TTO)), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.403666, scale = 0.2798643))


# Check probability of locoregional recurrence
(p_R2_LRR <- pars0$p_R2_SG_recur * (1 - pars0$p_R2_recur_distant))
mean(!is.na(events_rectalII$LRR))


# Check probability of distant recurrence, potentially also through LRR
p_R2_CR4 <- pars0$p_R2_SG_recur * pars0$p_R2_recur_distant
p_R2_LRR_CR4 <- p_R2_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SG_recur) * pars0$p_LRR_recur_distant
p_R2_CR4 + p_R2_LRR_CR4
mean(!is.na(events_rectalII$CR4_progr))




## 7. RECTAL III ----

n_sim <- 2*10^5
start_rectalIII <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'R3'
)

sim_rectalIII <- PRIMCAT_CRC_Simulation(start = start_rectalIII, pars = pars)

events_rectalIII <- as.data.frame(sim_rectalIII$eventLog)
resources_rectalIII <- as.data.frame(sim_rectalIII$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_rectalIII, 2, function(x) sum(!is.na(x)))


# Check time to treatment (TTT)
pars0$p_R3_noTTT
mean(events_rectalIII$R3_Tx == 0)

fit_R3_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_rectalIII.RDS'))
km_R3_TTT  <- survfit(formula = Surv(R3_Tx) ~ 1, data = events_rectalIII, subset = (events_rectalIII$R3_Tx > 0))

fit_R3_TTT + 
  geom_line(mapping = aes(x = km_R3_TTT$time, y = km_R3_TTT$surv), color = 'red')


# Check treatment utilization
pars0$p_R3_neoSG
mean(!is.na(resources_rectalIII$R3_Tx_neoSG))

pars0$p_R3_SG
mean(!is.na(resources_rectalIII$R3_Tx_SG))

pars0$p_R3_SGadj
mean(!is.na(resources_rectalIII$R3_Tx_SGadj))

pars0$p_R3_neoSGadj
mean(!is.na(resources_rectalIII$R3_Tx_neoSGadj))


# Check time to outcome (TTO)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
R3_TTO <- ifelse(!is.na(events_rectalIII$LRR), events_rectalIII$LRR - events_rectalIII$R3_Tx, events_rectalIII$CR4_progr - events_rectalIII$R3_Tx)

pars0$t_R3_neoSG_TTR
plot(density(na.omit(R3_TTO[!is.na(resources_rectalIII$R3_Tx_neoSG)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.561969, scale = 0.3884443))

pars0$t_R3_SG_TTR
plot(density(na.omit(R3_TTO[!is.na(resources_rectalIII$R3_Tx_SG)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.561969, scale = 0.3397029))

pars0$t_R3_SGadj_TTR
plot(density(na.omit(R3_TTO[!is.na(resources_rectalIII$R3_Tx_SGadj)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.561969, scale = 0.2639929))

pars0$t_R3_neoSGadj_TTR
plot(density(na.omit(R3_TTO[!is.na(resources_rectalIII$R3_Tx_neoSGadj)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.561969, scale = 0.3958142))


# Check probability of locoregional recurrence
(p_R3_LRR <- pars0$p_R3_neoSG_recur * (1 - pars0$p_R3_recur_distant))
mean(!is.na(events_rectalIII$LRR[!is.na(resources_rectalIII$R3_Tx_neoSG)]))

(p_R3_LRR <- pars0$p_R3_SG_recur * (1 - pars0$p_R3_recur_distant))
mean(!is.na(events_rectalIII$LRR[!is.na(resources_rectalIII$R3_Tx_SG)]))

(p_R3_LRR <- pars0$p_R3_SGadj_recur * (1 - pars0$p_R3_recur_distant))
mean(!is.na(events_rectalIII$LRR[!is.na(resources_rectalIII$R3_Tx_SGadj)]))

(p_R3_LRR <- pars0$p_R3_neoSGadj_recur * (1 - pars0$p_R3_recur_distant))
mean(!is.na(events_rectalIII$LRR[!is.na(resources_rectalIII$R3_Tx_neoSGadj)]))


# Check probability of distant recurrence, potentially also through LRR
p_R3_CR4 <- pars0$p_R3_neoSG_recur * pars0$p_R3_recur_distant
p_R3_LRR_CR4 <- p_R3_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SG_recur) * pars0$p_LRR_recur_distant
p_R3_CR4 + p_R3_LRR_CR4
mean(!is.na(events_rectalIII$CR4_progr[!is.na(resources_rectalIII$R3_Tx_neoSG)]))

p_R3_CR4 <- pars0$p_R3_SG_recur * pars0$p_R3_recur_distant
p_R3_LRR_CR4 <- p_R3_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SG_recur) * pars0$p_LRR_recur_distant
p_R3_CR4 + p_R3_LRR_CR4
mean(!is.na(events_rectalIII$CR4_progr[!is.na(resources_rectalIII$R3_Tx_SG)]))

p_R3_CR4 <- pars0$p_R3_SGadj_recur * pars0$p_R3_recur_distant
p_R3_LRR_CR4 <- p_R3_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SGadj_recur) * pars0$p_LRR_recur_distant
p_R3_CR4 + p_R3_LRR_CR4
mean(!is.na(events_rectalIII$CR4_progr[!is.na(resources_rectalIII$R3_Tx_SGadj)]))

p_R3_CR4 <- pars0$p_R3_neoSGadj_recur * pars0$p_R3_recur_distant
p_R3_LRR_CR4 <- p_R3_LRR * (pars0$p_LRR_systemic *  pars0$p_LRR_systemic_recur + (1 - pars0$p_LRR_systemic) * pars0$p_LRR_SGadj_recur) * pars0$p_LRR_recur_distant
p_R3_CR4 + p_R3_LRR_CR4
mean(!is.na(events_rectalIII$CR4_progr[!is.na(resources_rectalIII$R3_Tx_neoSGadj)]))




## 8. LOCOREG. RECURRENCE ----

n_sim <- 2*10^5
start_locreg <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'LRR'
)

sim_locreg <- PRIMCAT_CRC_Simulation(start = start_locreg, pars = pars)

events_locreg <- as.data.frame(sim_locreg$eventLog)
resources_locreg <- as.data.frame(sim_locreg$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_locreg, 2, function(x) sum(!is.na(x)))


# Check time to treatment (TTT)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
pars0$p_LRR_noTTT
mean(events_locreg$LRR_Tx == 0, na.rm = TRUE)

LRR_TTT <- events_locreg$LRR_Tx - events_locreg$LRR

pars0$t_LRR_TTT
plot(density(na.omit(LRR_TTT)), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dgompertz(x = seq(0, 30, 0.1), shape = -0.9299574, rate = 6.143221))


# Check treatment utilization
pars0$p_LRR_neoSG
mean(!is.na(resources_locreg$LRR_Tx_neoSG)[!is.na(events_locreg$LRR_Tx)])

pars0$p_LRR_SG
mean(!is.na(resources_locreg$LRR_Tx_SG)[!is.na(events_locreg$LRR_Tx)])

pars0$p_LRR_SGadj
mean(!is.na(resources_locreg$LRR_Tx_SGadj)[!is.na(events_locreg$LRR_Tx)])

pars0$p_LRR_neoSGadj
mean(!is.na(resources_locreg$LRR_Tx_neoSGadj)[!is.na(events_locreg$LRR_Tx)])

pars0$p_LRR_systemic
mean(!is.na(resources_locreg$LRR_Tx_systemic)[!is.na(events_locreg$LRR_Tx)])


# Check time to outcome (TTO)
# - Cannot be checked through Kaplan-Meier, as the competing event of death is not simulated
LRR_TTO <- events_locreg$CR4_progr - events_locreg$LRR_Tx

pars0$t_LRR_SG_TTR
plot(density(na.omit(LRR_TTO[!is.na(resources_locreg$LRR_Tx_SG)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.091366, scale = 0.1250659))


# Check probability of distant recurrence
pars0$p_LRR_SG_recur
mean(!is.na(events_locreg$CR4_progr[!is.na(resources_locreg$LRR_Tx_SG)]))

pars0$p_LRR_systemic_recur
mean(!is.na(events_locreg$CR4_progr[!is.na(resources_locreg$LRR_Tx_systemic)]))




## 9. CRC IV ----

n_sim <- 2*10^4
start_CR4 <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'CR4'
)

sim_CR4 <- PRIMCAT_CRC_Simulation(start = start_CR4, pars = pars)

events_CR4 <- as.data.frame(sim_CR4$eventLog)
resources_CR4 <- as.data.frame(sim_CR4$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_CR4, 2, function(x) sum(!is.na(x)))


# Check proportion RASwt
pars0$p_CR4_testRAS * pars0$p_CR4_RASwt
mean(!is.na(events_CR4$RASwt))

(1 - pars0$p_CR4_testRAS) + pars0$p_CR4_testRAS * (1 - pars0$p_CR4_RASwt)
mean(!is.na(events_CR4$RASmt))




## 10. RASwt ----

n_sim <- 2*10^5
start_RASwt <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'RASwt'
)

sim_RASwt <- PRIMCAT_CRC_Simulation(start = start_RASwt, pars = pars)

events_RASwt <- as.data.frame(sim_RASwt$eventLog)
resources_RASwt <- as.data.frame(sim_RASwt$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_RASwt, 2, function(x) sum(!is.na(x)))


### 10.1 Dx ----

# Proportion L0 vs L1
pars0$p_RASwt_L0
mean(!is.na(events_RASwt$RASwt_L0))


# Time to L0
pars0$p_RASwt_L0_noTTT
mean(na.omit(events_RASwt$RASwt_L0) == 0)

pars0$t_RASwt_L0_TTT
fit_RASwt_L0_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_DxL0.RDS'))
km_RASwt_L0_TTT  <- survfit(formula = Surv(RASwt_L0) ~ 1, data = events_RASwt, subset = (events_RASwt$RASwt_L0 > 0))

fit_RASwt_L0_TTT + 
  geom_line(mapping = aes(x = km_RASwt_L0_TTT$time, y = km_RASwt_L0_TTT$surv), color = 'red')


# Time to L1
pars0$p_RASwt_L1_noTTT
mean(na.omit(events_RASwt$RASwt_L1) == 0)

pars0$t_RASwt_L1_TTT
fit_RASwt_L1_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_DxL1.RDS'))
km_RASwt_L1_TTT  <- survfit(formula = Surv(RASwt_L1) ~ 1, data = events_RASwt, subset = (events_RASwt$RASwt_L1 > 0 & is.na(events_RASwt$RASwt_L0)))

fit_RASwt_L1_TTT + 
  geom_line(mapping = aes(x = km_RASwt_L1_TTT$time, y = km_RASwt_L1_TTT$surv), color = 'red')


### 10.2 L0 -----

events_RASwt_L0 <- events_RASwt[!is.na(events_RASwt$RASwt_L0), ]
resources_RASwt_L0 <- resources_RASwt[!is.na(events_RASwt$RASwt_L0), ]

# Treatment utilization
pars0$p_RASwt_L0_SGmets
mean(!is.na(resources_RASwt_L0$RASwt_L0_SGmets) & is.na(resources_RASwt_L0$RASwt_L0_SGprim) & is.na(resources_RASwt_L0$RASwt_L0_RT))

pars0$p_RASwt_L0_SGprim
mean(is.na(resources_RASwt_L0$RASwt_L0_SGmets) & !is.na(resources_RASwt_L0$RASwt_L0_SGprim) & is.na(resources_RASwt_L0$RASwt_L0_RT))

pars0$p_RASwt_L0_RT
mean(is.na(resources_RASwt_L0$RASwt_L0_SGmets) & is.na(resources_RASwt_L0$RASwt_L0_SGprim) & !is.na(resources_RASwt_L0$RASwt_L0_RT))

pars0$p_RASwt_L0_SGmetsSGprim
mean(!is.na(resources_RASwt_L0$RASwt_L0_SGmets) & !is.na(resources_RASwt_L0$RASwt_L0_SGprim) & is.na(resources_RASwt_L0$RASwt_L0_RT))

pars0$p_RASwt_L0_SGmetsRT
mean(!is.na(resources_RASwt_L0$RASwt_L0_SGmets) & is.na(resources_RASwt_L0$RASwt_L0_SGprim) & !is.na(resources_RASwt_L0$RASwt_L0_RT))


# Time to outcome
# - note that 12 weeks have to be added to the fitted model
RASwt_L0_TTO <- na.omit(events_RASwt_L0$RASwt_L1 - events_RASwt_L0$RASwt_L0)

pars0$t_RASwt_L0_TTP
plot(density(RASwt_L0_TTO), col = 'red', lwd = 2, lty = 2, xlim = c(0, 10))
lines(density(rweibullPH(n = 10^5, shape = 0.5263444, scale = 0.6644272) + (12/52)))


# Probability of progression
pars0$p_RASwt_L0_progr
length(RASwt_L0_TTO) / nrow(events_RASwt_L0)


### 10.3 L1 -----

events_RASwt_L1 <- events_RASwt[!is.na(events_RASwt$RASwt_L1), ]
resources_RASwt_L1 <- resources_RASwt[!is.na(events_RASwt$RASwt_L1), ]

# Treatment utilization
pars0$p_RASwt_L1_chemo
mean(!is.na(resources_RASwt_L1$RASwt_L1_chemo))

pars0$p_RASwt_L1_EGFR
mean(!is.na(resources_RASwt_L1$RASwt_L1_EGFR))

pars0$p_RASwt_L1_chemoEGFR
mean(!is.na(resources_RASwt_L1$RASwt_L1_chemoEGFR))

pars0$p_RASwt_L1_chemoVEGF
mean(!is.na(resources_RASwt_L1$RASwt_L1_chemoVEGF))


# Adjunct treatment utilization
pars0$p_RASwt_L1_SGmets
mean(!is.na(resources_RASwt_L1$RASwt_L1_SGmets))

pars0$p_RASwt_L1_SGprim
mean(!is.na(resources_RASwt_L1$RASwt_L1_SGprim))

pars0$p_RASwt_L1_RT
mean(!is.na(resources_RASwt_L1$RASwt_L1_RT))


# Time to outcome
RASwt_L1_TTO <- events_RASwt_L1$RASwt_L2 - events_RASwt_L1$RASwt_L1

pars0$t_RASwt_L1_chemo_TTP
plot(density(na.omit(RASwt_L1_TTO[!is.na(resources_RASwt_L1$RASwt_L1_chemo)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 0.9404586, scale = 0.5117598))

pars0$t_RASwt_L1_EGFR_TTP
plot(density(na.omit(RASwt_L1_TTO[!is.na(resources_RASwt_L1$RASwt_L1_chemoEGFR)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 0.9404586, scale = 0.4284408))


# Probability of progression
pars0$p_RASwt_L1_chemo_progr
nrow(events_RASwt_L1[!is.na(resources_RASwt_L1$RASwt_L1_chemo) & !is.na(events_RASwt_L1$RASwt_L2) , ]) / nrow(events_RASwt_L1[!is.na(resources_RASwt_L1$RASwt_L1_chemo) , ])

pars0$p_RASwt_L1_chemoEGFR_progr
nrow(events_RASwt_L1[!is.na(resources_RASwt_L1$RASwt_L1_chemoEGFR) & !is.na(events_RASwt_L1$RASwt_L2) , ]) / nrow(events_RASwt_L1[!is.na(resources_RASwt_L1$RASwt_L1_chemoEGFR) , ])


### 10.4 L2 -----

events_RASwt_L2 <- events_RASwt[!is.na(events_RASwt$RASwt_L2), ]
resources_RASwt_L2 <- resources_RASwt[!is.na(events_RASwt$RASwt_L2), ]

# Treatment utilization - L1 chemo
pars0$p_RASwt_L1chemo_L2_chemo
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemo) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_chemo))

pars0$p_RASwt_L1chemo_L2_EGFR
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemo) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_EGFR))

pars0$p_RASwt_L1chemo_L2_chemoEGFR
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemo) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_chemoEGFR))

pars0$p_RASwt_L1chemo_L2_chemoVEGF
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemo) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_chemoVEGF))


# Treatment utilization - L1 chemoEGFR
pars0$p_RASwt_L1EGFR_L2_chemo
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemoEGFR) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_chemo))

pars0$p_RASwt_L1EGFR_L2_chemoVEGF
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemoEGFR) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_chemoVEGF))


# Treatment utilization - L1 chemoVEGF
pars0$p_RASwt_L1VEGF_L2_chemo
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemoVEGF) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_chemo))

pars0$p_RASwt_L1VEGF_L2_EGFR
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemoVEGF) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_EGFR))

pars0$p_RASwt_L1VEGF_L2_chemoEGFR
mean(!is.na(resources_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L1_chemoVEGF) & !is.na(events_RASwt_L2$RASwt_L2), ]$RASwt_L2_chemoEGFR))


# Adjunct treatment utilization
pars0$p_RASwt_L2_SGmets
mean(!is.na(resources_RASwt_L2$RASwt_L2_SGmets))

pars0$p_RASwt_L2_SGprim
mean(!is.na(resources_RASwt_L2$RASwt_L2_SGprim))

pars0$p_RASwt_L2_RT
mean(!is.na(resources_RASwt_L2$RASwt_L2_RT))


# Time to outcome
RASwt_L2_TTO <- events_RASwt_L2$RASwt_L3 - events_RASwt_L2$RASwt_L2

pars0$t_RASwt_L2_chemo_TTP
plot(density(na.omit(RASwt_L2_TTO[!is.na(resources_RASwt_L2$RASwt_L2_chemo)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 5, 0.1), y = dweibullPH(x = seq(0, 5, 0.1), shape = 1.478937, scale = 1.572553))

pars0$t_RASwt_L2_chemoVEGF_TTP
plot(density(na.omit(RASwt_L2_TTO[!is.na(resources_RASwt_L2$RASwt_L2_chemoVEGF)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 5, 0.1), y = dweibullPH(x = seq(0, 5, 0.1), shape = 1.478937, scale = 1.031103))


# Probability of progression
pars0$p_RASwt_L2_chemo_progr
nrow(events_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L2_chemo) & !is.na(events_RASwt_L2$RASwt_L3) , ]) / nrow(events_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L2_chemo) , ])

pars0$p_RASwt_L2_chemoVEGF_progr
nrow(events_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L2_chemoVEGF) & !is.na(events_RASwt_L2$RASwt_L3) , ]) / nrow(events_RASwt_L2[!is.na(resources_RASwt_L2$RASwt_L2_chemoVEGF) , ])




### 10.5 L3 -----

events_RASwt_L3 <- events_RASwt[!is.na(events_RASwt$RASwt_L3), ]
resources_RASwt_L3 <- resources_RASwt[!is.na(events_RASwt$RASwt_L3), ]

# Treatment utilization
pars0$p_RASwt_L3_chemo
mean(!is.na(resources_RASwt_L3$RASwt_L3_chemo))


# Adjunct treatment utilization
pars0$p_RASwt_L3_RT
mean(!is.na(resources_RASwt_L3$RASwt_L3_RT))


# Time to outcome
RASwt_L3_TTO <- events_RASwt_L3$RASwt_L4 - events_RASwt_L3$RASwt_L3

pars0$t_RASwt_L3_chemo_TTP
plot(density(na.omit(RASwt_L3_TTO)), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.387365, scale = 1.213713))


# Probability of progression
pars0$p_RASwt_L3_chemo_progr
nrow(events_RASwt_L3[!is.na(resources_RASwt_L3$RASwt_L3_chemo) & !is.na(events_RASwt_L3$RASwt_L4) , ]) / nrow(events_RASwt_L3[!is.na(resources_RASwt_L3$RASwt_L3_chemo) , ])




### 10.6 L4 -----

events_RASwt_L4 <- events_RASwt[!is.na(events_RASwt$RASwt_L4), ]
resources_RASwt_L4 <- resources_RASwt[!is.na(events_RASwt$RASwt_L4), ]

# Treatment utilization
pars0$p_RASwt_L4_chemo
mean(!is.na(resources_RASwt_L4$RASwt_L4_chemo))


# Adjunct treatment utilization
pars0$p_RASwt_L4_RT
mean(!is.na(resources_RASwt_L4$RASwt_L4_RT))




## 11. RASmt ----

n_sim <- 2*10^5
start_RASmt <- data.frame(
  startTime  = rep(0, n_sim),
  firstEvent = 'RASmt'
)

sim_RASmt <- PRIMCAT_CRC_Simulation(start = start_RASmt, pars = pars)

events_RASmt <- as.data.frame(sim_RASmt$eventLog)
resources_RASmt <- as.data.frame(sim_RASmt$resourceLog[ , , 'start'])


# Check what events occur to assess whether there are events that should not happen
apply(events_RASmt, 2, function(x) sum(!is.na(x)))


### 11.1 Dx ----

# Proportion L0 vs L1
pars0$p_RASmt_L0
mean(!is.na(events_RASmt$RASmt_L0))


# Time to L0
pars0$p_RASmt_L0_noTTT
mean(na.omit(events_RASmt$RASmt_L0) == 0)

pars0$t_RASmt_L0_TTT
fit_RASmt_L0_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_DxL0.RDS'))
km_RASmt_L0_TTT  <- survfit(formula = Surv(RASmt_L0) ~ 1, data = events_RASmt, subset = (events_RASmt$RASmt_L0 > 0))

fit_RASmt_L0_TTT + 
  geom_line(mapping = aes(x = km_RASmt_L0_TTT$time, y = km_RASmt_L0_TTT$surv), color = 'red')


# Time to L1
pars0$p_RASmt_L1_noTTT
mean(na.omit(events_RASmt$RASmt_L1) == 0)

pars0$t_RASmt_L1_TTT
fit_RASmt_L1_TTT <- readRDS(file = paste0(path_plots, 'plot_TTT_DxL1.RDS'))
km_RASmt_L1_TTT  <- survfit(formula = Surv(RASmt_L1) ~ 1, data = events_RASmt, subset = (events_RASmt$RASmt_L1 > 0 & is.na(events_RASmt$RASmt_L0)))

fit_RASmt_L1_TTT + 
  geom_line(mapping = aes(x = km_RASmt_L1_TTT$time, y = km_RASmt_L1_TTT$surv), color = 'red')


### 11.2 L0 -----

events_RASmt_L0 <- events_RASmt[!is.na(events_RASmt$RASmt_L0), ]
resources_RASmt_L0 <- resources_RASmt[!is.na(events_RASmt$RASmt_L0), ]

# Treatment utilization
pars0$p_RASmt_L0_SGmets
mean(!is.na(resources_RASmt_L0$RASmt_L0_SGmets) & is.na(resources_RASmt_L0$RASmt_L0_SGprim) & is.na(resources_RASmt_L0$RASmt_L0_RT))

pars0$p_RASmt_L0_SGprim
mean(is.na(resources_RASmt_L0$RASmt_L0_SGmets) & !is.na(resources_RASmt_L0$RASmt_L0_SGprim) & is.na(resources_RASmt_L0$RASmt_L0_RT))

pars0$p_RASmt_L0_RT
mean(is.na(resources_RASmt_L0$RASmt_L0_SGmets) & is.na(resources_RASmt_L0$RASmt_L0_SGprim) & !is.na(resources_RASmt_L0$RASmt_L0_RT))

pars0$p_RASmt_L0_SGmetsSGprim
mean(!is.na(resources_RASmt_L0$RASmt_L0_SGmets) & !is.na(resources_RASmt_L0$RASmt_L0_SGprim) & is.na(resources_RASmt_L0$RASmt_L0_RT))

pars0$p_RASmt_L0_SGmetsRT
mean(!is.na(resources_RASmt_L0$RASmt_L0_SGmets) & is.na(resources_RASmt_L0$RASmt_L0_SGprim) & !is.na(resources_RASmt_L0$RASmt_L0_RT))


# Time to outcome
# - note that 12 weeks have to be added to the fitted model
RASmt_L0_TTO <- na.omit(events_RASmt_L0$RASmt_L1 - events_RASmt_L0$RASmt_L0)

pars0$t_RASmt_L0_TTP
plot(density(RASmt_L0_TTO), col = 'red', lwd = 2, lty = 2, xlim = c(0, 10))
lines(density(rweibullPH(n = 10^5, shape = 0.5263444, scale = 0.491539) + (12/52)))


# Probability of progression
pars0$p_RASmt_L0_progr
length(RASmt_L0_TTO) / nrow(events_RASmt_L0)


### 11.3 L1 -----

events_RASmt_L1 <- events_RASmt[!is.na(events_RASmt$RASmt_L1), ]
resources_RASmt_L1 <- resources_RASmt[!is.na(events_RASmt$RASmt_L1), ]

# Treatment utilization
pars0$p_RASmt_L1_chemo
mean(!is.na(resources_RASmt_L1$RASmt_L1_chemo))

pars0$p_RASmt_L1_chemoVEGF
mean(!is.na(resources_RASmt_L1$RASmt_L1_chemoVEGF))


# Adjunct treatment utilization
pars0$p_RASmt_L1_SGmets
mean(!is.na(resources_RASmt_L1$RASmt_L1_SGmets))

pars0$p_RASmt_L1_SGprim
mean(!is.na(resources_RASmt_L1$RASmt_L1_SGprim))

pars0$p_RASmt_L1_RT
mean(!is.na(resources_RASmt_L1$RASmt_L1_RT))


# Time to outcome
RASmt_L1_TTO <- events_RASmt_L1$RASmt_L2 - events_RASmt_L1$RASmt_L1

pars0$t_RASmt_L1_chemo_TTP
plot(density(na.omit(RASmt_L1_TTO[!is.na(resources_RASmt_L1$RASmt_L1_chemo)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.555574, scale = 0.7350749))

pars0$t_RASmt_L1_chemoVEGF_TTP
plot(density(na.omit(RASmt_L1_TTO[!is.na(resources_RASmt_L1$RASmt_L1_chemoVEGF)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.555574, scale = 0.6392009))


# Probability of progression
pars0$p_RASmt_L1_chemo_progr
nrow(events_RASmt_L1[!is.na(resources_RASmt_L1$RASmt_L1_chemo) & !is.na(events_RASmt_L1$RASmt_L2) , ]) / nrow(events_RASmt_L1[!is.na(resources_RASmt_L1$RASmt_L1_chemo) , ])

pars0$p_RASmt_L1_chemoVEGF_progr
nrow(events_RASmt_L1[!is.na(resources_RASmt_L1$RASmt_L1_chemoVEGF) & !is.na(events_RASmt_L1$RASmt_L2) , ]) / nrow(events_RASmt_L1[!is.na(resources_RASmt_L1$RASmt_L1_chemoVEGF) , ])


### 11.4 L2 -----

events_RASmt_L2 <- events_RASmt[!is.na(events_RASmt$RASmt_L2), ]
resources_RASmt_L2 <- resources_RASmt[!is.na(events_RASmt$RASmt_L2), ]

# Treatment utilization - L1 chemo
pars0$p_RASmt_L1chemo_L2_chemo
mean(!is.na(resources_RASmt_L2[!is.na(resources_RASmt_L2$RASmt_L1_chemo) & !is.na(events_RASmt_L2$RASmt_L2), ]$RASmt_L2_chemo))

pars0$p_RASmt_L1chemo_L2_chemoVEGF
mean(!is.na(resources_RASmt_L2[!is.na(resources_RASmt_L2$RASmt_L1_chemo) & !is.na(events_RASmt_L2$RASmt_L2), ]$RASmt_L2_chemoVEGF))


# Treatment utilization - L1 chemoVEGF
pars0$p_RASmt_L1VEGF_L2_chemo
mean(!is.na(resources_RASmt_L2[!is.na(resources_RASmt_L2$RASmt_L1_chemoVEGF) & !is.na(events_RASmt_L2$RASmt_L2), ]$RASmt_L2_chemo))


# Adjunct treatment utilization
pars0$p_RASmt_L2_SGmets
mean(!is.na(resources_RASmt_L2$RASmt_L2_SGmets))

pars0$p_RASmt_L2_SGprim
mean(!is.na(resources_RASmt_L2$RASmt_L2_SGprim))

pars0$p_RASmt_L2_RT
mean(!is.na(resources_RASmt_L2$RASmt_L2_RT))


# Time to outcome
RASmt_L2_TTO <- events_RASmt_L2$RASmt_L3 - events_RASmt_L2$RASmt_L2

pars0$t_RASmt_L2_chemo_TTP
plot(density(na.omit(RASmt_L2_TTO[!is.na(resources_RASmt_L2$RASmt_L2_chemo)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 5, 0.1), y = dweibullPH(x = seq(0, 5, 0.1), shape = 1.478937, scale = 1.572553))

pars0$t_RASmt_L2_chemoVEGF_TTP
plot(density(na.omit(RASmt_L2_TTO[!is.na(resources_RASmt_L2$RASmt_L2_chemoVEGF)])), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 5, 0.1), y = dweibullPH(x = seq(0, 5, 0.1), shape = 1.478937, scale = 1.031103))


# Probability of progression
pars0$p_RASmt_L2_chemo_progr
nrow(events_RASmt_L2[!is.na(resources_RASmt_L2$RASmt_L2_chemo) & !is.na(events_RASmt_L2$RASmt_L3) , ]) / nrow(events_RASmt_L2[!is.na(resources_RASmt_L2$RASmt_L2_chemo) , ])

pars0$p_RASmt_L2_chemoVEGF_progr
nrow(events_RASmt_L2[!is.na(resources_RASmt_L2$RASmt_L2_chemoVEGF) & !is.na(events_RASmt_L2$RASmt_L3) , ]) / nrow(events_RASmt_L2[!is.na(resources_RASmt_L2$RASmt_L2_chemoVEGF) , ])




### 11.5 L3 -----

events_RASmt_L3 <- events_RASmt[!is.na(events_RASmt$RASmt_L3), ]
resources_RASmt_L3 <- resources_RASmt[!is.na(events_RASmt$RASmt_L3), ]

# Treatment utilization
pars0$p_RASmt_L3_chemo
mean(!is.na(resources_RASmt_L3$RASmt_L3_chemo))


# Adjunct treatment utilization
pars0$p_RASmt_L3_RT
mean(!is.na(resources_RASmt_L3$RASmt_L3_RT))


# Time to outcome
RASmt_L3_TTO <- events_RASmt_L3$RASmt_L4 - events_RASmt_L3$RASmt_L3

pars0$t_RASmt_L3_chemo_TTP
plot(density(na.omit(RASmt_L3_TTO)), col = 'red', lwd = 2, lty = 2)
lines(x = seq(0, 30, 0.1), y = dweibullPH(x = seq(0, 30, 0.1), shape = 1.387365, scale = 1.213713))


# Probability of progression
pars0$p_RASmt_L3_chemo_progr
nrow(events_RASmt_L3[!is.na(resources_RASmt_L3$RASmt_L3_chemo) & !is.na(events_RASmt_L3$RASmt_L4) , ]) / nrow(events_RASmt_L3[!is.na(resources_RASmt_L3$RASmt_L3_chemo) , ])




### 11.6 L4 -----

events_RASmt_L4 <- events_RASmt[!is.na(events_RASmt$RASmt_L4), ]
resources_RASmt_L4 <- resources_RASmt[!is.na(events_RASmt$RASmt_L4), ]

# Treatment utilization
pars0$p_RASmt_L4_chemo
mean(!is.na(resources_RASmt_L4$RASmt_L4_chemo))


# Adjunct treatment utilization
pars0$p_RASmt_L4_RT
mean(!is.na(resources_RASmt_L4$RASmt_L4_RT))




## 12. TIME-DEPENDENT NEW TREATMENTS ----

# This section verifies whether the time-dependent parameters and parameters for the NEW treatment have been
# implemented correctly. This is done for each disease stage separately. The new treatment is to be simulated to
# be not utilized in year 0, then in year 1 there is some use, and from year 2 it is the only treatment used.

library(data.table)   # v1.14.0


### 12.1 COLON I ----

n_sim <- 2*10^5
start_colonI_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'C1'
)

pars0$p_C1_SG_recur
pars0$t_C1_SG_TTR

pars_colonI_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_C1_SG  = 0.5,
    p_C1_NEW = 0.5,
    p_C1_NEW_recur = 0.05,
    t_C1_NEW_TTR   = pars0$t_C1_SG_TTR
  ),
  '2' = list(
    p_C1_SG  = 0,
    p_C1_NEW = 1
  )
))

sim_colonI_NEW <- PRIMCAT_CRC_Simulation(start = start_colonI_NEW, pars = pars_colonI_NEW)

events_colonI_NEW <- as.data.frame(sim_colonI_NEW$eventLog)
resources_colonI_NEW <- as.data.frame(sim_colonI_NEW$resourceLog[ , , 'start'])

df_colonI_NEW <- as.data.table(cbind(events_colonI_NEW, resources_colonI_NEW))

# Treatment utilization
df_colonI_NEW[ , .(n = .N, p_SG  = mean(!is.na(C1_Tx_SG))),  by = floor(C1_Tx)][order(floor)]
df_colonI_NEW[ , .(n = .N, p_NEW = mean(!is.na(C1_Tx_NEW))), by = floor(C1_Tx)][order(floor)]

# Probability of recurrence
df_colonI_NEW[!is.na(C1_Tx_SG),  .(n = .N, p_SG_recur  = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C1_Tx)][order(floor)]
df_colonI_NEW[!is.na(C1_Tx_NEW), .(n = .N, p_NEW_recur = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C1_Tx)][order(floor)]


### 12.2 COLON II ----

n_sim <- 2*10^5
start_colonII_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'C2'
)

pars0$p_C2_SG
pars0$p_C2_SGadj
pars0$p_C2_SG_recur
pars0$t_C2_SG_TTR

pars_colonII_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_C2_SG    = 0.5*pars0$p_C2_SG,
    p_C2_SGadj = 0.5*pars0$p_C2_SGadj,
    p_C2_NEW   = 0.5,
    p_C2_NEW_recur = 0.10,
    t_C2_NEW_TTR   = pars0$t_C2_SG_TTR
  ),
  '2' = list(
    p_C2_SG = 0,
    p_C2_SGadj = 0,
    p_C2_NEW = 1
  )
))

sim_colonII_NEW <- PRIMCAT_CRC_Simulation(start = start_colonII_NEW, pars = pars_colonII_NEW)

events_colonII_NEW <- as.data.frame(sim_colonII_NEW$eventLog)
resources_colonII_NEW <- as.data.frame(sim_colonII_NEW$resourceLog[ , , 'start'])

df_colonII_NEW <- as.data.table(cbind(events_colonII_NEW, resources_colonII_NEW))

# Treatment utilization
df_colonII_NEW[ , .(n = .N, p_SG    = mean(!is.na(C2_Tx_SG))),    by = floor(C2_Tx)][order(floor)]
df_colonII_NEW[ , .(n = .N, p_SGadj = mean(!is.na(C2_Tx_SGadj))), by = floor(C2_Tx)][order(floor)]
df_colonII_NEW[ , .(n = .N, p_NEW   = mean(!is.na(C2_Tx_NEW))),   by = floor(C2_Tx)][order(floor)]

# Probability of recurrence
df_colonII_NEW[!is.na(C2_Tx_SG),    .(n = .N, p_SG_recur   = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C2_Tx)][order(floor)]
df_colonII_NEW[!is.na(C2_Tx_SGadj), .(n = .N, p_SGadjrecur = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C2_Tx)][order(floor)]
df_colonII_NEW[!is.na(C2_Tx_NEW),   .(n = .N, p_NEW_recur  = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C2_Tx)][order(floor)]


### 12.3 COLON III ----

n_sim <- 2*10^5
start_colonIII_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'C3'
)

pars0$p_C3_SG
pars0$p_C3_SGadj
pars0$p_C3_SG_recur
pars0$t_C3_SG_TTR

pars_colonIII_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_C3_SG    = 0.5*pars0$p_C3_SG,
    p_C3_SGadj = 0.5*pars0$p_C3_SGadj,
    p_C3_NEW   = 0.5,
    p_C3_NEW_recur = 0.20,
    t_C3_NEW_TTR   = pars0$t_C3_SG_TTR
  ),
  '2' = list(
    p_C3_SG = 0,
    p_C3_SGadj = 0,
    p_C3_NEW = 1
  )
))

sim_colonIII_NEW <- PRIMCAT_CRC_Simulation(start = start_colonIII_NEW, pars = pars_colonIII_NEW)

events_colonIII_NEW <- as.data.frame(sim_colonIII_NEW$eventLog)
resources_colonIII_NEW <- as.data.frame(sim_colonIII_NEW$resourceLog[ , , 'start'])

df_colonIII_NEW <- as.data.table(cbind(events_colonIII_NEW, resources_colonIII_NEW))

# Treatment utilization
df_colonIII_NEW[ , .(n = .N, p_SG    = mean(!is.na(C3_Tx_SG))),    by = floor(C3_Tx)][order(floor)]
df_colonIII_NEW[ , .(n = .N, p_SGadj = mean(!is.na(C3_Tx_SGadj))), by = floor(C3_Tx)][order(floor)]
df_colonIII_NEW[ , .(n = .N, p_NEW   = mean(!is.na(C3_Tx_NEW))),   by = floor(C3_Tx)][order(floor)]

# Probability of recurrence
df_colonIII_NEW[!is.na(C3_Tx_SG),    .(n = .N, p_SG_recur    = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C3_Tx)][order(floor)]
df_colonIII_NEW[!is.na(C3_Tx_SGadj), .(n = .N, p_SGadj_recur = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C3_Tx)][order(floor)]
df_colonIII_NEW[!is.na(C3_Tx_NEW),   .(n = .N, p_NEW_recur   = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(C3_Tx)][order(floor)]


### 12.4 RECTAL I ----

n_sim <- 2*10^5
start_rectalI_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'R1'
)

pars0$p_R1_SG
pars0$p_R1_neoSG
pars0$p_R1_SGadj
pars0$p_R1_neoSGadj
pars0$p_R1_SG_recur
pars0$t_R1_SG_TTR

pars_rectalI_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_R1_SG       = 0.5*pars0$p_R1_SG,
    p_R1_neoSG    = 0.5*pars0$p_R1_neoSG,
    p_R1_SGadj    = 0.5*pars0$p_R1_SGadj,
    p_R1_neoSGadj = 0.5*pars0$p_R1_neoSGadj,
    p_R1_NEW      = 0.5,
    p_R1_NEW_recur = 0.05,
    t_R1_NEW_TTR   = pars0$t_R1_SG_TTR
  ),
  '2' = list(
    p_R1_SG       = 0,
    p_R1_neoSG    = 0,
    p_R1_SGadj    = 0,
    p_R1_neoSGadj = 0,
    p_R1_NEW      = 1
  )
))

sim_rectalI_NEW <- PRIMCAT_CRC_Simulation(start = start_rectalI_NEW, pars = pars_rectalI_NEW)

events_rectalI_NEW <- as.data.frame(sim_rectalI_NEW$eventLog)
resources_rectalI_NEW <- as.data.frame(sim_rectalI_NEW$resourceLog[ , , 'start'])

df_rectalI_NEW <- as.data.table(cbind(events_rectalI_NEW, resources_rectalI_NEW))

# Treatment utilization
df_rectalI_NEW[ , .(n = .N, p_SG       = mean(!is.na(R1_Tx_SG))),       by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[ , .(n = .N, p_neoSG    = mean(!is.na(R1_Tx_neoSG))),    by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[ , .(n = .N, p_SGadj    = mean(!is.na(R1_Tx_SGadj))),    by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[ , .(n = .N, p_neoSGadj = mean(!is.na(R1_Tx_neoSGadj))), by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[ , .(n = .N, p_NEW      = mean(!is.na(R1_Tx_NEW))),      by = floor(R1_Tx)][order(floor)]

# Probability of recurrence
df_rectalI_NEW[!is.na(R1_Tx_SG),       .(n = .N, p_SG_recur       = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[!is.na(R1_Tx_neoSG),    .(n = .N, p_neoSG_recur    = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[!is.na(R1_Tx_SGadj),    .(n = .N, p_SGadj_recur    = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[!is.na(R1_Tx_neoSGadj), .(n = .N, p_SGneoadj_recur = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R1_Tx)][order(floor)]
df_rectalI_NEW[!is.na(R1_Tx_NEW),      .(n = .N, p_NEW_recur      = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R1_Tx)][order(floor)]


### 12.5 RECTAL II ----

n_sim <- 2*10^5
start_rectalII_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'R2'
)

pars0$p_R2_SG
pars0$p_R2_neoSG
pars0$p_R2_SGadj
pars0$p_R2_neoSGadj
pars0$p_R2_SG_recur
pars0$t_R2_SG_TTR

pars_rectalII_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_R2_SG       = 0.5*pars0$p_R2_SG,
    p_R2_neoSG    = 0.5*pars0$p_R2_neoSG,
    p_R2_SGadj    = 0.5*pars0$p_R2_SGadj,
    p_R2_neoSGadj = 0.5*pars0$p_R2_neoSGadj,
    p_R2_NEW      = 0.5,
    p_R2_NEW_recur = 0.10,
    t_R2_NEW_TTR   = pars0$t_R2_SG_TTR
  ),
  '2' = list(
    p_R2_SG       = 0,
    p_R2_neoSG    = 0,
    p_R2_SGadj    = 0,
    p_R2_neoSGadj = 0,
    p_R2_NEW      = 1
  )
))

sim_rectalII_NEW <- PRIMCAT_CRC_Simulation(start = start_rectalII_NEW, pars = pars_rectalII_NEW)

events_rectalII_NEW <- as.data.frame(sim_rectalII_NEW$eventLog)
resources_rectalII_NEW <- as.data.frame(sim_rectalII_NEW$resourceLog[ , , 'start'])

df_rectalII_NEW <- as.data.table(cbind(events_rectalII_NEW, resources_rectalII_NEW))

# Treatment utilization
df_rectalII_NEW[ , .(n = .N, p_SG       = mean(!is.na(R2_Tx_SG))),       by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[ , .(n = .N, p_neoSG    = mean(!is.na(R2_Tx_neoSG))),    by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[ , .(n = .N, p_SGadj    = mean(!is.na(R2_Tx_SGadj))),    by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[ , .(n = .N, p_neoSGadj = mean(!is.na(R2_Tx_neoSGadj))), by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[ , .(n = .N, p_NEW      = mean(!is.na(R2_Tx_NEW))),      by = floor(R2_Tx)][order(floor)]

# Probability of recurrence
df_rectalII_NEW[!is.na(R2_Tx_SG),       .(n = .N, p_SG_recur       = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[!is.na(R2_Tx_neoSG),    .(n = .N, p_neoSG_recur    = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[!is.na(R2_Tx_SGadj),    .(n = .N, p_SGadj_recur    = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[!is.na(R2_Tx_neoSGadj), .(n = .N, p_neoSGadj_recur = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R2_Tx)][order(floor)]
df_rectalII_NEW[!is.na(R2_Tx_NEW),      .(n = .N, p_NEW_recur      = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R2_Tx)][order(floor)]


### 12.6 RECTAL III ----

n_sim <- 2*10^5
start_rectalIII_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'R3'
)

pars0$p_R3_SG
pars0$p_R3_neoSG
pars0$p_R3_SGadj
pars0$p_R3_neoSGadj
pars0$p_R3_SG_recur
pars0$t_R3_SG_TTR

pars_rectalIII_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_R3_SG       = 0.5*pars0$p_R3_SG,
    p_R3_neoSG    = 0.5*pars0$p_R3_neoSG,
    p_R3_SGadj    = 0.5*pars0$p_R3_SGadj,
    p_R3_neoSGadj = 0.5*pars0$p_R3_neoSGadj,
    p_R3_NEW      = 0.5,
    p_R3_NEW_recur = 0.20,
    t_R3_NEW_TTR   = pars0$t_R3_SG_TTR
  ),
  '2' = list(
    p_R3_SG       = 0,
    p_R3_neoSG    = 0,
    p_R3_SGadj    = 0,
    p_R3_neoSGadj = 0,
    p_R3_NEW      = 1
  )
))

sim_rectalIII_NEW <- PRIMCAT_CRC_Simulation(start = start_rectalIII_NEW, pars = pars_rectalIII_NEW)

events_rectalIII_NEW <- as.data.frame(sim_rectalIII_NEW$eventLog)
resources_rectalIII_NEW <- as.data.frame(sim_rectalIII_NEW$resourceLog[ , , 'start'])

df_rectalIII_NEW <- as.data.table(cbind(events_rectalIII_NEW, resources_rectalIII_NEW))

# Treatment utilization
df_rectalIII_NEW[ , .(n = .N, p_SG       = mean(!is.na(R3_Tx_SG))),       by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[ , .(n = .N, p_neoSG    = mean(!is.na(R3_Tx_neoSG))),    by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[ , .(n = .N, p_SGadj    = mean(!is.na(R3_Tx_SGadj))),    by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[ , .(n = .N, p_neoSGadj = mean(!is.na(R3_Tx_neoSGadj))), by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[ , .(n = .N, p_NEW      = mean(!is.na(R3_Tx_NEW))),      by = floor(R3_Tx)][order(floor)]

# Probability of recurrence
df_rectalIII_NEW[!is.na(R3_Tx_SG),       .(n = .N, p_SG_recur       = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[!is.na(R3_Tx_neoSG),    .(n = .N, p_neoSG_recur    = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[!is.na(R3_Tx_SGadj),    .(n = .N, p_SGadj_recur    = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[!is.na(R3_Tx_neoSGadj), .(n = .N, p_neoSGadj_recur = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R3_Tx)][order(floor)]
df_rectalIII_NEW[!is.na(R3_Tx_NEW),      .(n = .N, p_NEW_recur      = mean(!is.na(LRR) | !is.na(CR4_progr))), by = floor(R3_Tx)][order(floor)]


### 12.7 LOCOREG. RECURRENCE ----

n_sim <- 2*10^5
start_locreg_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'LRR'
)

pars0$p_LRR_SG
pars0$p_LRR_neoSG
pars0$p_LRR_SGadj
pars0$p_LRR_neoSGadj
pars0$p_LRR_systemic
pars0$p_LRR_SG_recur
pars0$t_LRR_SG_TTR

pars_locreg_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_LRR_SG       = 0.5*pars0$p_LRR_SG,
    p_LRR_neoSG    = 0.5*pars0$p_LRR_neoSG,
    p_LRR_SGadj    = 0.5*pars0$p_LRR_SGadj,
    p_LRR_neoSGadj = 0.5*pars0$p_LRR_neoSGadj,
    p_LRR_systemic = 0.5*pars0$p_LRR_systemic,
    p_LRR_NEW      = 0.5,
    p_LRR_NEW_recur = 0.30,
    t_LRR_NEW_TTR   = pars0$t_LRR_SG_TTR
  ),
  '2' = list(
    p_LRR_SG       = 0,
    p_LRR_neoSG    = 0,
    p_LRR_SGadj    = 0,
    p_LRR_neoSGadj = 0,
    p_LRR_systemic = 0,
    p_LRR_NEW      = 1
  )
))

sim_locreg_NEW <- PRIMCAT_CRC_Simulation(start = start_locreg_NEW, pars = pars_locreg_NEW)

events_locreg_NEW <- as.data.frame(sim_locreg_NEW$eventLog)
resources_locreg_NEW <- as.data.frame(sim_locreg_NEW$resourceLog[ , , 'start'])

df_locreg_NEW <- as.data.table(cbind(events_locreg_NEW, resources_locreg_NEW))

# Treatment utilization
df_locreg_NEW[ , .(n = .N, p_SG       = mean(!is.na(LRR_Tx_SG))),       by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[ , .(n = .N, p_neoSG    = mean(!is.na(LRR_Tx_neoSG))),    by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[ , .(n = .N, p_SGadj    = mean(!is.na(LRR_Tx_SGadj))),    by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[ , .(n = .N, p_neoSGadj = mean(!is.na(LRR_Tx_neoSGadj))), by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[ , .(n = .N, p_systemic = mean(!is.na(LRR_Tx_systemic))), by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[ , .(n = .N, p_NEW      = mean(!is.na(LRR_Tx_NEW))),      by = floor(LRR_Tx)][order(floor)]

# Probability of recurrence
df_locreg_NEW[!is.na(LRR_Tx_SG),       .(n = .N, p_SG_recur       = mean(!is.na(CR4_progr))), by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[!is.na(LRR_Tx_neoSG),    .(n = .N, p_neoSG_recur    = mean(!is.na(CR4_progr))), by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[!is.na(LRR_Tx_SGadj),    .(n = .N, p_SGadj_recur    = mean(!is.na(CR4_progr))), by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[!is.na(LRR_Tx_neoSGadj), .(n = .N, p_neoSGadj_recur = mean(!is.na(CR4_progr))), by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[!is.na(LRR_Tx_systemic), .(n = .N, p_systemic_recur = mean(!is.na(CR4_progr))), by = floor(LRR_Tx)][order(floor)]
df_locreg_NEW[!is.na(LRR_Tx_NEW),      .(n = .N, p_NEW_recur      = mean(!is.na(CR4_progr))), by = floor(LRR_Tx)][order(floor)]


### 12.8 RASwt L1 ----

n_sim <- 2*10^5
start_RASwtL1_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'RASwt_L1'
)

pars0$p_RASwt_L1_chemo
pars0$p_RASwt_L1_EGFR
pars0$p_RASwt_L1_chemoEGFR
pars0$p_RASwt_L1_chemoVEGF

pars0$p_RASwt_L1_chemo_progr
pars0$t_RASwt_L1_chemo_TTP

pars_RASwtL1_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_RASwt_L1_chemo     = 0.5*pars0$p_RASwt_L1_chemo,
    p_RASwt_L1_EGFR      = 0.5*pars0$p_RASwt_L1_EGFR,
    p_RASwt_L1_chemoEGFR = 0.5*pars0$p_RASwt_L1_chemoEGFR,
    p_RASwt_L1_chemoVEGF = 0.5*pars0$p_RASwt_L1_chemoVEGF,
    p_RASwt_L1_NEW       = 0.5,
    p_RASwt_L1NEW_L2_chemo     = 0.4,
    p_RASwt_L1NEW_L2_EGFR      = 0.3,
    p_RASwt_L1NEW_L2_chemoEGFR = 0.2,
    p_RASwt_L1NEW_L2_chemoVEGF = 0.1,
    p_RASwt_L1NEW_L2_NEW       = 0.0,
    p_RASwt_L1_NEW_progr = 0.80,
    t_RASwt_L1_NEW_TTP   = pars0$t_RASwt_L1_chemo_TTP
  ),
  '2' = list(
    p_RASwt_L1_chemo     = 0,
    p_RASwt_L1_EGFR      = 0,
    p_RASwt_L1_chemoEGFR = 0,
    p_RASwt_L1_chemoVEGF = 0,
    p_RASwt_L1_NEW       = 1
  )
))

sim_RASwtL1_NEW <- PRIMCAT_CRC_Simulation(start = start_RASwtL1_NEW, pars = pars_RASwtL1_NEW)

events_RASwtL1_NEW <- as.data.frame(sim_RASwtL1_NEW$eventLog)
resources_RASwtL1_NEW <- as.data.frame(sim_RASwtL1_NEW$resourceLog[ , , 'start'])

df_RASwtL1_NEW <- as.data.table(cbind(events_RASwtL1_NEW, resources_RASwtL1_NEW))

# Treatment utilization
df_RASwtL1_NEW[ , .(n = .N, p_chemo     = mean(!is.na(RASwt_L1_chemo))),     by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[ , .(n = .N, p_EGFR      = mean(!is.na(RASwt_L1_EGFR))),      by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[ , .(n = .N, p_chemoEGFR = mean(!is.na(RASwt_L1_chemoEGFR))), by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[ , .(n = .N, p_chemoVEGF = mean(!is.na(RASwt_L1_chemoVEGF))), by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[ , .(n = .N, p_NEW       = mean(!is.na(RASwt_L1_NEW))),       by = floor(RASwt_L1)][order(floor)]

# Probability of progression
df_RASwtL1_NEW[!is.na(RASwt_L1_chemo),     .(n = .N, p_chemo_progr     = mean(!is.na(RASwt_L2))), by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_EGFR),      .(n = .N, p_EGFR_progr      = mean(!is.na(RASwt_L2))), by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_chemoEGFR), .(n = .N, p_chemoEGFR_progr = mean(!is.na(RASwt_L2))), by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_chemoVEGF), .(n = .N, p_chemoVEGF_progr = mean(!is.na(RASwt_L2))), by = floor(RASwt_L1)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_NEW),       .(n = .N, p_NEW_progr       = mean(!is.na(RASwt_L2))), by = floor(RASwt_L1)][order(floor)]

# Treatment utilization (L2)
df_RASwtL1_NEW[!is.na(RASwt_L1_NEW), .(n = .N, p_chemo     = mean(!is.na(RASwt_L2_chemo))),     by = floor(RASwt_L2)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_NEW), .(n = .N, p_EGFR      = mean(!is.na(RASwt_L2_EGFR))),      by = floor(RASwt_L2)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_NEW), .(n = .N, p_chemoEGFR = mean(!is.na(RASwt_L2_chemoEGFR))), by = floor(RASwt_L2)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_NEW), .(n = .N, p_chemoVEGF = mean(!is.na(RASwt_L2_chemoVEGF))), by = floor(RASwt_L2)][order(floor)]
df_RASwtL1_NEW[!is.na(RASwt_L1_NEW), .(n = .N, p_NEW       = mean(!is.na(RASwt_L2_NEW))),       by = floor(RASwt_L2)][order(floor)]


### 12.9 RASwt L3 ----

n_sim <- 2*10^5
start_RASwtL3_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'RASwt_L3'
)

pars0$p_RASwt_L3_chemo

pars0$p_RASwt_L3_chemo_progr
pars0$t_RASwt_L3_chemo_TTP

pars_RASwtL3_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_RASwt_L3_chemo = 0.5,
    p_RASwt_L3_NEW   = 0.5,
    p_RASwt_L3_NEW_progr = 0.80,
    t_RASwt_L3_NEW_TTP   = pars0$t_RASwt_L3_chemo_TTP
  ),
  '2' = list(
    p_RASwt_L3_chemo = 0,
    p_RASwt_L3_NEW   = 1
  )
))

sim_RASwtL3_NEW <- PRIMCAT_CRC_Simulation(start = start_RASwtL3_NEW, pars = pars_RASwtL3_NEW)

events_RASwtL3_NEW <- as.data.frame(sim_RASwtL3_NEW$eventLog)
resources_RASwtL3_NEW <- as.data.frame(sim_RASwtL3_NEW$resourceLog[ , , 'start'])

df_RASwtL3_NEW <- as.data.table(cbind(events_RASwtL3_NEW, resources_RASwtL3_NEW))

# Treatment utilization
df_RASwtL3_NEW[ , .(n = .N, p_chemo = mean(!is.na(RASwt_L3_chemo))), by = floor(RASwt_L3)][order(floor)]
df_RASwtL3_NEW[ , .(n = .N, p_NEW   = mean(!is.na(RASwt_L3_NEW))),   by = floor(RASwt_L3)][order(floor)]

# Probability of progression
df_RASwtL3_NEW[!is.na(RASwt_L3_chemo), .(n = .N, p_chemo_progr = mean(!is.na(RASwt_L4))), by = floor(RASwt_L3)][order(floor)]
df_RASwtL3_NEW[!is.na(RASwt_L3_NEW),   .(n = .N, p_NEW_progr   = mean(!is.na(RASwt_L4))), by = floor(RASwt_L3)][order(floor)]


### 12.10 RASmt L1 ----

n_sim <- 2*10^5
start_RASmtL1_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'RASmt_L1'
)

pars0$p_RASmt_L1_chemo
pars0$p_RASmt_L1_chemoVEGF

pars0$p_RASmt_L1_chemo_progr
pars0$t_RASmt_L1_chemo_TTP

pars_RASmtL1_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_RASmt_L1_chemo     = 0.5*pars0$p_RASmt_L1_chemo,
    p_RASmt_L1_chemoVEGF = 0.5*pars0$p_RASmt_L1_chemoVEGF,
    p_RASmt_L1_NEW       = 0.5,
    p_RASmt_L1NEW_L2_chemo     = 0.6,
    p_RASmt_L1NEW_L2_chemoVEGF = 0.4,
    p_RASmt_L1NEW_L2_NEW       = 0.0,
    p_RASmt_L1_NEW_progr = 0.80,
    t_RASmt_L1_NEW_TTP  = pars0$t_RASmt_L1_chemo_TTP
  ),
  '2' = list(
    p_RASmt_L1_chemo     = 0,
    p_RASmt_L1_chemoVEGF = 0,
    p_RASmt_L1_NEW       = 1
  )
))

sim_RASmtL1_NEW <- PRIMCAT_CRC_Simulation(start = start_RASmtL1_NEW, pars = pars_RASmtL1_NEW)

events_RASmtL1_NEW <- as.data.frame(sim_RASmtL1_NEW$eventLog)
resources_RASmtL1_NEW <- as.data.frame(sim_RASmtL1_NEW$resourceLog[ , , 'start'])

df_RASmtL1_NEW <- as.data.table(cbind(events_RASmtL1_NEW, resources_RASmtL1_NEW))

# Treatment utilization
df_RASmtL1_NEW[ , .(n = .N, p_chemo     = mean(!is.na(RASmt_L1_chemo))),     by = floor(RASmt_L1)][order(floor)]
df_RASmtL1_NEW[ , .(n = .N, p_chemoVEGF = mean(!is.na(RASmt_L1_chemoVEGF))), by = floor(RASmt_L1)][order(floor)]
df_RASmtL1_NEW[ , .(n = .N, p_NEW       = mean(!is.na(RASmt_L1_NEW))),       by = floor(RASmt_L1)][order(floor)]

# Probability of progression
df_RASmtL1_NEW[!is.na(RASmt_L1_chemo),     .(n = .N, p_chemo_progr     = mean(!is.na(RASmt_L2))), by = floor(RASmt_L1)][order(floor)]
df_RASmtL1_NEW[!is.na(RASmt_L1_chemoVEGF), .(n = .N, p_chemoVEGF_progr = mean(!is.na(RASmt_L2))), by = floor(RASmt_L1)][order(floor)]
df_RASmtL1_NEW[!is.na(RASmt_L1_NEW),       .(n = .N, p_NEW_progr       = mean(!is.na(RASmt_L2))), by = floor(RASmt_L1)][order(floor)]

# Treatment utilization (L2)
df_RASmtL1_NEW[!is.na(RASmt_L1_NEW), .(n = .N, p_chemo     = mean(!is.na(RASmt_L2_chemo))),     by = floor(RASmt_L2)][order(floor)]
df_RASmtL1_NEW[!is.na(RASmt_L1_NEW), .(n = .N, p_chemoVEGF = mean(!is.na(RASmt_L2_chemoVEGF))), by = floor(RASmt_L2)][order(floor)]
df_RASmtL1_NEW[!is.na(RASmt_L1_NEW), .(n = .N, p_NEW       = mean(!is.na(RASmt_L2_NEW))),       by = floor(RASmt_L2)][order(floor)]


### 12.11 RASmt L3 ----

n_sim <- 2*10^5
start_RASmtL3_NEW <- data.frame(
  startTime  = sample(x = 0:2, size = n_sim, replace = TRUE),
  firstEvent = 'RASmt_L3'
)

pars0$p_RASmt_L3_chemo

pars0$p_RASmt_L3_chemo_progr
pars0$t_RASmt_L3_chemo_TTP

pars_RASmtL3_NEW <- PRIMCAT_CRC_Parameters(starttime = 0, updates = list(
  '1' = list(
    p_RASmt_L3_chemo = 0.5,
    p_RASmt_L3_NEW   = 0.5,
    p_RASmt_L3_NEW_progr = 0.80,
    t_RASmt_L3_NEW_TTP   = pars0$t_RASmt_L3_chemo_TTP
  ),
  '2' = list(
    p_RASmt_L3_chemo = 0,
    p_RASmt_L3_NEW   = 1
  )
))

sim_RASmtL3_NEW <- PRIMCAT_CRC_Simulation(start = start_RASmtL3_NEW, pars = pars_RASmtL3_NEW)

events_RASmtL3_NEW <- as.data.frame(sim_RASmtL3_NEW$eventLog)
resources_RASmtL3_NEW <- as.data.frame(sim_RASmtL3_NEW$resourceLog[ , , 'start'])

df_RASmtL3_NEW <- as.data.table(cbind(events_RASmtL3_NEW, resources_RASmtL3_NEW))

# Treatment utilization
df_RASmtL3_NEW[ , .(n = .N, p_chemo = mean(!is.na(RASmt_L3_chemo))), by = floor(RASmt_L3)][order(floor)]
df_RASmtL3_NEW[ , .(n = .N, p_NEW   = mean(!is.na(RASmt_L3_NEW))),   by = floor(RASmt_L3)][order(floor)]

# Probability of progression
df_RASmtL3_NEW[!is.na(RASmt_L3_chemo), .(n = .N, p_chemo_progr = mean(!is.na(RASmt_L4))), by = floor(RASmt_L3)][order(floor)]
df_RASmtL3_NEW[!is.na(RASmt_L3_NEW),   .(n = .N, p_NEW_progr   = mean(!is.na(RASmt_L4))), by = floor(RASmt_L3)][order(floor)]



