#-----------------------------------------------------------------------------------------#
# The Scenario Analysis capability of PRIMCAT-CRC is illustrated in a case study.         #
#                                                                                         #
# To demonstrate the use of PRIMCAT-CRC, the 5-year impact (2022-2026) of introducing     #
# pembrolizumab as first-line treatment for dMMR metastatic CRC                           #
# in Australia was estimated. For context, pembrolizumab received taxpayer reimbursement  #
# approval in August 2021, post our study’s data cut-off.                                 #
#                                                                                         #
# The dMMR prevalence among patients diagnosed with stage IV was estimated from the       #
# literature and the TRACC registry, yielding estimates of 15%32 and 6.9%, respectively.  #
#                                                                                         #
# The uptake rate of dMMR testing was modeled on historical RAS testing trends from the   #
# TRACC registry, which initially rose from 42% to 82%. Anticipating a similar pattern    #
# for dMMR, we project an increase from 42% in 2022 to 82% by 2023.                       #
#                                                                                         #
# André et al. provided key evidence on pembrolizumab’s outcomes, with a hazard ratio     #
# of 0.6 and a 19% mortality reduction. We used these data to calculate relative risks    #
# (RRs) for disease progression in the context of RAS status, comparing with standard     #
# chemotherapy and biological regimens. The derived RRs were 1.07 for RASwt and           #
# 1.11 for RASmt.                                                                         #
#                                                                                         #
# Four scenarios were considered:                                                         #
# • S1: dMMR prevalence 15%, testing uptake 100%.                                         #
# • S2: as per S1, considering HR 0.6 and RR progression of 1.07 (RASwt) and 1.11 (RASmt).#
# • S3: dMMR prevalence 6.9%, testing uptake 42% (2022), 84% from 2023.                   #
# • S4: as per S3, considering HR 0.6 and RR progression of 1.07 (RASwt) and 1.11 (RASmt).#
#                                                                                         #
#-----------------------------------------------------------------------------------------#
#
#
# This script perform the DES model for the case study analysis.
#
# The PRIMCAT-CRC model structure is described in the manuscript and available 
# as supplementary material at the following link: https://doi.org/10.1016/j.jval.2024.06.006
#
#
#-----------------------------------------------------------------------------------------#
#
#
### 1. INITIALIZATION ----

rm(list = ls()); gc(); 

library(nnet)
library(dplyr)
library(tidyr)
library(ggplot2)

source(file = 'Script/1_PRIMCAT_CRC_functions.R')


### 1. SIMULATION ----

# Pembrolizumab scenario analysis
# - First-line treatment of unresectable or metastatic mismatch repair deficient colorectal cancer
# - Estimated prevalence of mismatch repair deficiency cited in NEJM paper is 15%
# - Estimated prevalence of mismatch repair deficiency in Australian patients with metastatic colorectal is 6.9% based on TRACC registry 
# - Assumed uptake based on uptake of RAS testing in TRACC:
#   - 42% tested in 2022 and 84% tested from 2023 onward
# - Comparator is doublet + targeted (EGFR for RASwt, VEGF for RASmt)
# - Further progression assumed to be in line with chemo + cetuximab/panitumumab (RASwt) or chemo + bevacizumab (RASmt)
# - Scenario analysis based on 69 (control) vs 56 (pembro) deaths in NEJM paper
# - Further treatments assumed to be in line with base case first-line treatment utilization
# - Hazard ratio on time to progression: 0.60 based on NEJM paper
# - 28 parameters to be defined based on length(unlist(update_pembrolizumab))

# Scenarios:
#   1. Base case - utilization only: all tested and 15% uptake, other impact in line with chemoEGFR/chemoVEGF
#   2. Base case - utilization + progression: Scenario 1 + hazard ratio 0.6 on TTP + 19% reduction in death
#   3. Australian - utilization only: 6.9% uptake, and 42% tested in 2022 and 84% tested from 2023 onward
#   4. Australian - utilization + progression: Scenario 3 + hazard ratio 0.6 on TTP + 19% reduction in death

pars_basecase <- PRIMCAT_CRC_Parameters()$`0`

updates_pembrolizumab_scenario1 <- list(
  '2022' = list(
    # RASwt
    p_RASwt_L1_chemo     = pars_basecase$p_RASwt_L1_chemo * (1 - 0.15),
    p_RASwt_L1_EGFR      = pars_basecase$p_RASwt_L1_EGFR * (1 - 0.15),
    p_RASwt_L1_chemoEGFR = pars_basecase$p_RASwt_L1_chemoEGFR * (1 - 0.15),
    p_RASwt_L1_chemoVEGF = pars_basecase$p_RASwt_L1_chemoVEGF * (1 - 0.15),
    p_RASwt_L1_NEW       = 0.15,
    p_RASwt_L1_NEW_progr = pars_basecase$p_RASwt_L1_chemoEGFR_progr,
    t_RASwt_L1_NEW_TTP   = pars_basecase$t_RASwt_L1_chemoEGFR_TTP,
    p_RASwt_L1NEW_L2_chemo     = pars_basecase$p_RASwt_L1_chemo,
    p_RASwt_L1NEW_L2_EGFR      = pars_basecase$p_RASwt_L1_EGFR,
    p_RASwt_L1NEW_L2_chemoEGFR = pars_basecase$p_RASwt_L1_chemoEGFR,
    p_RASwt_L1NEW_L2_chemoVEGF = pars_basecase$p_RASwt_L1_chemoVEGF,
    p_RASwt_L1NEW_L2_NEW       = 0,
    # RASmt
    p_RASmt_L1_chemo     = pars_basecase$p_RASmt_L1_chemo * (1 - 0.15),
    p_RASmt_L1_chemoVEGF = pars_basecase$p_RASmt_L1_chemoVEGF * (1 - 0.15),
    p_RASmt_L1_NEW       = 0.15,
    p_RASmt_L1_NEW_progr = pars_basecase$p_RASmt_L1_chemoVEGF_progr,
    t_RASmt_L1_NEW_TTP   = pars_basecase$t_RASmt_L1_chemoVEGF_TTP,
    p_RASmt_L1NEW_L2_chemo     = pars_basecase$p_RASmt_L1_chemo,
    p_RASmt_L1NEW_L2_chemoVEGF = pars_basecase$p_RASmt_L1_chemoVEGF,
    p_RASmt_L1NEW_L2_NEW       = 0
  )
)

updates_pembrolizumab_scenario2 <- updates_pembrolizumab_scenario1
updates_pembrolizumab_scenario2$`2022`$p_RASwt_L1_NEW_progr          <- 1 - (1 - 0.19) * (1 - updates_pembrolizumab_scenario1$`2022`$p_RASwt_L1_NEW_progr)
updates_pembrolizumab_scenario2$`2022`$t_RASwt_L1_NEW_TTP$pars$scale <- 0.60 * updates_pembrolizumab_scenario1$`2022`$t_RASwt_L1_NEW_TTP$pars$scale
updates_pembrolizumab_scenario2$`2022`$p_RASmt_L1_NEW_progr          <- 1 - (1 - 0.19) * (1 - updates_pembrolizumab_scenario1$`2022`$p_RASmt_L1_NEW_progr)
updates_pembrolizumab_scenario2$`2022`$t_RASmt_L1_NEW_TTP$pars$scale <- 0.60 * updates_pembrolizumab_scenario1$`2022`$t_RASmt_L1_NEW_TTP$pars$scale

updates_pembrolizumab_scenario3 <- list(
  '2022' = list(
    # RASwt
    p_RASwt_L1_chemo     = pars_basecase$p_RASwt_L1_chemo * (1 - 0.42*0.069),
    p_RASwt_L1_EGFR      = pars_basecase$p_RASwt_L1_EGFR * (1 - 0.42*0.069),
    p_RASwt_L1_chemoEGFR = pars_basecase$p_RASwt_L1_chemoEGFR * (1 - 0.42*0.069),
    p_RASwt_L1_chemoVEGF = pars_basecase$p_RASwt_L1_chemoVEGF * (1 - 0.42*0.069),
    p_RASwt_L1_NEW       = 0.42*0.069,
    p_RASwt_L1_NEW_progr = pars_basecase$p_RASwt_L1_chemoEGFR_progr,
    t_RASwt_L1_NEW_TTP   = pars_basecase$t_RASwt_L1_chemoEGFR_TTP,
    p_RASwt_L1NEW_L2_chemo     = pars_basecase$p_RASwt_L1_chemo,
    p_RASwt_L1NEW_L2_EGFR      = pars_basecase$p_RASwt_L1_EGFR,
    p_RASwt_L1NEW_L2_chemoEGFR = pars_basecase$p_RASwt_L1_chemoEGFR,
    p_RASwt_L1NEW_L2_chemoVEGF = pars_basecase$p_RASwt_L1_chemoVEGF,
    p_RASwt_L1NEW_L2_NEW       = 0,
    # RASmt
    p_RASmt_L1_chemo     = pars_basecase$p_RASmt_L1_chemo * (1 - 0.42*0.069),
    p_RASmt_L1_chemoVEGF = pars_basecase$p_RASmt_L1_chemoVEGF * (1 - 0.42*0.069),
    p_RASmt_L1_NEW       = 0.42*0.069,
    p_RASmt_L1_NEW_progr = pars_basecase$p_RASmt_L1_chemoVEGF_progr,
    t_RASmt_L1_NEW_TTP   = pars_basecase$t_RASmt_L1_chemoVEGF_TTP,
    p_RASmt_L1NEW_L2_chemo     = pars_basecase$p_RASmt_L1_chemo,
    p_RASmt_L1NEW_L2_chemoVEGF = pars_basecase$p_RASmt_L1_chemoVEGF,
    p_RASmt_L1NEW_L2_NEW       = 0
  ),
  '2023' = list(
    # RASwt
    p_RASwt_L1_chemo     = pars_basecase$p_RASwt_L1_chemo * (1 - 0.84*0.069),
    p_RASwt_L1_EGFR      = pars_basecase$p_RASwt_L1_EGFR * (1 - 0.84*0.069),
    p_RASwt_L1_chemoEGFR = pars_basecase$p_RASwt_L1_chemoEGFR * (1 - 0.84*0.069),
    p_RASwt_L1_chemoVEGF = pars_basecase$p_RASwt_L1_chemoVEGF * (1 - 0.84*0.069),
    p_RASwt_L1_NEW       = 0.84*0.069,
    # RASmt
    p_RASmt_L1_chemo     = pars_basecase$p_RASmt_L1_chemo * (1 - 0.84*0.069),
    p_RASmt_L1_chemoVEGF = pars_basecase$p_RASmt_L1_chemoVEGF * (1 - 0.84*0.069),
    p_RASmt_L1_NEW       = 0.84*0.069
  )
)

updates_pembrolizumab_scenario4 <- updates_pembrolizumab_scenario3
updates_pembrolizumab_scenario4$`2022`$p_RASwt_L1_NEW_progr          <- 1 - (1 - 0.19) * (1 - updates_pembrolizumab_scenario3$`2022`$p_RASwt_L1_NEW_progr)
updates_pembrolizumab_scenario4$`2022`$t_RASwt_L1_NEW_TTP$pars$scale <- 0.60 * updates_pembrolizumab_scenario3$`2022`$t_RASwt_L1_NEW_TTP$pars$scale
updates_pembrolizumab_scenario4$`2022`$p_RASmt_L1_NEW_progr          <- 1 - (1 - 0.19) * (1 - updates_pembrolizumab_scenario3$`2022`$p_RASmt_L1_NEW_progr)
updates_pembrolizumab_scenario4$`2022`$t_RASmt_L1_NEW_TTP$pars$scale <- 0.60 * updates_pembrolizumab_scenario3$`2022`$t_RASmt_L1_NEW_TTP$pars$scale

system.time({
  sim_pembrolizumab_scenario1 <- PRIMCAT_CRC(updates = updates_pembrolizumab_scenario1, n_runs = 50, n_nodes = 10, m_years = m_years, by_total_incidence = FALSE)
  saveRDS(object = sim_pembrolizumab_scenario1, file = 'Output/sim_pembrolizumab_scenario1 50runs 20211213.RDS')
})

system.time({
  sim_pembrolizumab_scenario2 <- PRIMCAT_CRC(updates = updates_pembrolizumab_scenario2, n_runs = 50, n_nodes = 10, m_years = m_years, by_total_incidence = FALSE)
  saveRDS(object = sim_pembrolizumab_scenario2, file = 'Output/sim_pembrolizumab_scenario2 50runs 20211213.RDS')
})

system.time({
  sim_pembrolizumab_scenario3 <- PRIMCAT_CRC(updates = updates_pembrolizumab_scenario3, n_runs = 50, n_nodes = 10, m_years = m_years, by_total_incidence = FALSE)
  saveRDS(object = sim_pembrolizumab_scenario3, file = 'Output/sim_pembrolizumab_scenario3 50runs 20211213.RDS')
})

system.time({
  sim_pembrolizumab_scenario4 <- PRIMCAT_CRC(updates = updates_pembrolizumab_scenario4, n_runs = 50, n_nodes = 10, m_years = m_years, by_total_incidence = FALSE)
  saveRDS(object = sim_pembrolizumab_scenario4, file = 'Output/sim_pembrolizumab_scenario4 50runs 20211213.RDS')
})


### 2. RESULTS ANALYSIS -----

df <- bind_rows(
  readRDS(file = 'Output/sim_basecase 50runs 20211213.RDS') %>% mutate(scenario = 'Base case'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario1 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 1'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario2 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 2'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario3 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 3'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario4 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 4')
)

df <- df %>%
  group_by(scenario) %>% 
  mutate(t_index = 1:n()) %>%
  ungroup() %>% 
  rename(t_label = time_period) %>%
  relocate(t_index, t_label, .before = C1)




