#-----------------------------------------------------------------------------------------#
# The Base Case scenario uses the default model parameters informed by ACCORD and TRACC   #
# analyses to reflect the current treatment pathways in colorectal cancer.                #
#                                                                                         #
# The base-case analysis forecasts the number of CRC patients treated under current       #
# standard care in Australia from 2022 to 2026, serving as a baseline for evaluating      #
# the impact of new treatments.                                                           #
#                                                                                         #
# There are multiple steps:                                                               #
# 1. Modelling the incidence of patients with colorectal cancer                           #
# 2. Modelling the stage distribution at diagnosis                                        #
# 3. Simulating the number of patients as per the model framework established             #
#    3.1 Colon cancer stage I (C1), stage II (C2), stage III (C3)                         #
#    3.2 Rectal cancer stage I (R1), stage II (R2), stage III (R3)                        #
#    3.3 Colorectal cancer stage IV (CR4)                                                 #

#-----------------------------------------------------------------------------------------#
#
#
# This script perform the DES model for the base case scenario analysis of CRC, not 
# accounting for any new treatment options as per 2022 standard of care in Australia. 
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

AIHW_CRC_Incidence         <- readxl::read_xlsx(path = 'Input/AIHW_CRC_Incidence.xlsx', sheet = 'Output')
CVDL_CRC_StageDistribution <- readxl::read_xlsx(path = 'Input/CVDL_CRC_StageDistribution.xlsx', sheet = 'Output')


### 2. MODELING INCIDENCE AND STAGE DISTRIBUTION ----

## 2.1 Incidence ----

df_incidence <- data.frame(
  row.names = NULL,
  year   = as.integer(gsub(pattern = 'year_', replacement = '', x = colnames(AIHW_CRC_Incidence[, -1]))),
  colon  = unlist(AIHW_CRC_Incidence[1, -1]),
  rectal = unlist(AIHW_CRC_Incidence[2, -1])
) %>% 
  mutate(
    year_indexed = as.integer(year - 2010)
  )

lm_colon  <- lm(formula = colon ~ year_indexed, data = df_incidence)
lm_rectal <- lm(formula = rectal ~ year_indexed, data = df_incidence)

x_real <- 2010:2026
x_pred <- 0:16

extrapolated_incidence_colon  <- predict(lm_colon, data.frame(year_indexed = x_pred))
extrapolated_incidence_rectal <- predict(lm_rectal, data.frame(year_indexed = x_pred))

ggplot_incidence <- ggplot() + 
  geom_point(data = filter(df_incidence, year <= 2017), mapping = aes(x = year, y = colon, shape = 'AIHW reported', color = 'Colon')) +  
  geom_point(data = filter(df_incidence, year <= 2017), mapping = aes(x = year, y = rectal, shape = 'AIHW reported', color = 'Rectal')) + 
  geom_point(data = filter(df_incidence, year > 2017), mapping = aes(x = year, y = colon, shape = 'AIHW predicted', color = 'Colon')) +  
  geom_point(data = filter(df_incidence, year > 2017), mapping = aes(x = year, y = rectal, shape = 'AIHW predicted', color = 'Rectal')) + 
  geom_line(mapping = aes(x = x_real, y = extrapolated_incidence_colon, linetype = 'Extrapolation', color = 'Colon')) +
  geom_line(mapping = aes(x = x_real, y = extrapolated_incidence_rectal, linetype = 'Extrapolation', color = 'Rectal')) +
  theme(legend.title = element_blank()) +
  labs(x = 'Year', y = 'Incidence') + 
  ylim(0, 12000)

# Get incidences for simulation based on data and extrapolation
incidence_colon <- extrapolated_incidence_colon
incidence_colon[x_real <= 2021] <- df_incidence$colon[df_incidence$year <= 2021]

incidence_rectal <- extrapolated_incidence_rectal
incidence_rectal[x_real <= 2021] <- df_incidence$rectal[df_incidence$year <= 2021]


## 2.2 Stage distribution ----  

df_stage_colon <- data.frame(
  row.names = NULL,
  year      = as.integer(gsub(pattern = 'year_', replacement = '', x = colnames(CVDL_CRC_StageDistribution[, 2:11]))),
  stageI    = unlist(CVDL_CRC_StageDistribution[1, 2:11]),
  stageII   = unlist(CVDL_CRC_StageDistribution[2, 2:11]),
  stageIII  = unlist(CVDL_CRC_StageDistribution[3, 2:11]),
  stageIV   = unlist(CVDL_CRC_StageDistribution[4, 2:11])
) %>% 
  mutate(
    year_indexed = as.integer(year - 2010)
  )

df_stage_rectal <- data.frame(
  row.names = NULL,
  year      = as.integer(gsub(pattern = 'year_', replacement = '', x = colnames(CVDL_CRC_StageDistribution[, 2:11]))),
  stageI    = unlist(CVDL_CRC_StageDistribution[5, 2:11]),
  stageII   = unlist(CVDL_CRC_StageDistribution[6, 2:11]),
  stageIII  = unlist(CVDL_CRC_StageDistribution[7, 2:11]),
  stageIV   = unlist(CVDL_CRC_StageDistribution[8, 2:11])
) %>% 
  mutate(
    year_indexed = as.integer(year - 2010)
  )

mnm_colon  <- multinom(formula = as.matrix(df_stage_colon[ , 2:5]) ~ year_indexed, data = df_stage_colon)
mnm_rectal <- multinom(formula = as.matrix(df_stage_rectal[ , 2:5]) ~ year_indexed, data = df_stage_colon)

extrapolated_stage_colon  <- as.data.frame(predict(mnm_colon, data.frame(year_indexed = x_pred), 'probs'))
extrapolated_stage_rectal <- as.data.frame(predict(mnm_rectal, data.frame(year_indexed = x_pred), 'probs'))

ggplot_stage_colon <- ggplot() + 
  geom_point(data = df_stage_colon, mapping = aes(x = year, y = stageI, color = 'Stage I')) +  
  geom_point(data = df_stage_colon, mapping = aes(x = year, y = stageII, color = 'Stage II')) +  
  geom_point(data = df_stage_colon, mapping = aes(x = year, y = stageIII, color = 'Stage III')) +  
  geom_point(data = df_stage_colon, mapping = aes(x = year, y = stageIV, color = 'Stage IV')) +  
  geom_line(data = extrapolated_stage_colon, mapping = aes(x = x_real, y = stageI, color = 'Stage I')) +
  geom_line(data = extrapolated_stage_colon, mapping = aes(x = x_real, y = stageII, color = 'Stage II')) +
  geom_line(data = extrapolated_stage_colon, mapping = aes(x = x_real, y = stageIII, color = 'Stage III')) +
  geom_line(data = extrapolated_stage_colon, mapping = aes(x = x_real, y = stageIV, color = 'Stage IV')) +
  labs(x = 'Year', y = 'Stage Proportion', color = 'Colon') + 
  ylim(0, 0.5)

ggplot_stage_rectal <- ggplot() + 
  geom_point(data = df_stage_rectal, mapping = aes(x = year, y = stageI, color = 'Stage I')) +  
  geom_point(data = df_stage_rectal, mapping = aes(x = year, y = stageII, color = 'Stage II')) +  
  geom_point(data = df_stage_rectal, mapping = aes(x = year, y = stageIII, color = 'Stage III')) +  
  geom_point(data = df_stage_rectal, mapping = aes(x = year, y = stageIV, color = 'Stage IV')) +  
  geom_line(data = extrapolated_stage_rectal, mapping = aes(x = x_real, y = stageI, color = 'Stage I')) +
  geom_line(data = extrapolated_stage_rectal, mapping = aes(x = x_real, y = stageII, color = 'Stage II')) +
  geom_line(data = extrapolated_stage_rectal, mapping = aes(x = x_real, y = stageIII, color = 'Stage III')) +
  geom_line(data = extrapolated_stage_rectal, mapping = aes(x = x_real, y = stageIV, color = 'Stage IV')) +
  labs(x = 'Year', y = 'Stage Proportion', color = 'Rectal') + 
  ylim(0, 0.5)

# Get stage distributions for simulation based on data and extrapolation
stage_colon <- extrapolated_stage_colon
stage_colon[x_real <= 2019, ] <- df_stage_colon[df_stage_colon$year <= 2019, c('stageI', 'stageII', 'stageIII', 'stageIV')]

stage_rectal <- extrapolated_stage_rectal
stage_rectal[x_real <= 2019, ] <- df_stage_rectal[df_stage_rectal$year <= 2019, c('stageI', 'stageII', 'stageIII', 'stageIV')]


### 3. SIMULATION ----

# - t       the simulation time at which the year should start
# - n_C1    the incidence of C1 in that year
# - n_C2    the incidence of C2 in that year
# - n_C3    the incidence of C3 in that year
# - n_R1    the incidence of R1 in that year
# - n_R2    the incidence of R2 in that year
# - n_R3    the incidence of R3 in that year
# - n_CR4   the incidence of CR4 in that year
m_years <- cbind(
  t = x_real,
  n_C1  = round(incidence_colon * stage_colon$stageI, digits = 0),
  n_C2  = round(incidence_colon * stage_colon$stageII, digits = 0),
  n_C3  = round(incidence_colon * stage_colon$stageIII, digits = 0),
  n_R1  = round(incidence_rectal * stage_rectal$stageI, digits = 0),
  n_R2  = round(incidence_rectal * stage_rectal$stageII, digits = 0),
  n_R3  = round(incidence_rectal * stage_rectal$stageIII, digits = 0),
  n_CR4 = round(incidence_colon * stage_colon$stageIV + incidence_rectal * stage_rectal$stageIV, digits = 0)
)

m_years

# Simulation of 262,346 individuals: 1.7 minute
sim_start <- PRIMCAT_CRC_Population(m_years = m_years, by_total_incidence = FALSE)
sim_pars  <- PRIMCAT_CRC_Parameters()
sim_test  <- PRIMCAT_CRC_Simulation(start = sim_start, pars = sim_pars)

# 10 runs of 262,346 individuals on 10 nodes: 4 minutes
# 10 runs of 262,346 individuals on 5 nodes:  5 minutes
# 50 runs of 262,346 individuals on 10 nodes: 21 minutes
system.time({
  sim_basecase <- PRIMCAT_CRC(n_runs = 50, n_nodes = 10, m_years = m_years, by_total_incidence = FALSE)
  saveRDS(object = sim_basecase, file = 'Output/sim_basecase 50runs 20211213.RDS')
})



### 4. RESULTS ANALYSIS -----

df <- readRDS(file = 'Output/sim_basecase 50runs 20211213.RDS') %>% mutate(scenario = 'Base case')

df <- df %>%
  group_by(scenario) %>% 
  mutate(t_index = 1:n()) %>%
  ungroup() %>% 
  rename(t_label = time_period) %>%
  relocate(t_index, t_label, .before = C1)

# Incidence by disease stage
ggplot_sim_incidence <- df %>%
  select(t_index, t_label, scenario, C1, C2, C3, R1, R2, R3, LRR, CR4) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(name = factor(x = name, levels = c('C1', 'R1', 'C2', 'R2', 'C3', 'R3', 'LRR', 'CR4'))) %>%
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = name, linetype = scenario)) +
  scale_x_continuous(breaks = seq(0, 28*4, 8), labels = df$t_label[seq(0, 28*4, 8)+1]) +
  labs(
    title = 'Different Disease Types and Stages',
    y     = 'Number of Individuals'
  ) +
  theme(
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

# Treatment utilization for colon cancer
ggplot_sim_treatment_colon <- df %>% 
  filter(t_index >= 49) %>% 
  filter(scenario == 'Base case') %>% 
  select(t_index, t_label, scenario, contains('C1_Tx_'), contains('C2_Tx_'), contains('C3_Tx_')) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(
    stage = case_when(
      grepl(pattern = '^C1', x = name) ~ 'stage I',
      grepl(pattern = '^C2', x = name) ~ 'stage II',
      grepl(pattern = '^C3', x = name) ~ 'stage III'
    ),
    treatment = factor(
      x = gsub(pattern = 'C._Tx_', replacement = '', x = name),
      levels = c('SG', 'SGadj', 'NEW')
    )
  ) %>% 
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment, linetype = stage)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  labs(
    title = 'Treatment utilization for stage I-III colon cancer',
    y     = 'Number of Individuals'
  ) +
  theme(
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

# Treatment utilization for rectal cancer
ggplot_sim_treatment_rectal <- df %>% 
  filter(t_index >= 49) %>% 
  filter(scenario == 'Base case') %>% 
  select(t_index, t_label, scenario, contains('R1_Tx_'), contains('R2_Tx_'), contains('R3_Tx_')) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(
    stage = case_when(
      grepl(pattern = '^R1', x = name) ~ 'stage I',
      grepl(pattern = '^R2', x = name) ~ 'stage II',
      grepl(pattern = '^R3', x = name) ~ 'stage III'
    ),
    treatment = factor(
      x = gsub(pattern = 'R._Tx_', replacement = '', x = name),
      levels = c('SG', 'neoSG', 'SGadj', 'neoSGadj', 'NEW')
    )
  ) %>% 
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment, linetype = stage)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  labs(
    title = 'Treatment utilization for stage I-III rectal cancer',
    y     = 'Number of Individuals'
  ) +
  theme(
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

# Treatment utilization for locoregional recurrence
ggplot_sim_treatment_locreg <- df %>% 
  filter(t_index >= 49) %>% 
  filter(scenario == 'Base case') %>% 
  select(t_index, t_label, scenario, contains('LRR_Tx_')) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(
    treatment = factor(
      x = gsub(pattern = 'LRR_Tx_', replacement = '', x = name),
      levels = c('SG', 'neoSG', 'SGadj', 'neoSGadj', 'systemic', 'NEW')
    )
  ) %>% 
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  labs(
    title = 'Treatment utilization for locoregional recurrences',
    y     = 'Number of Individuals'
  ) +
  theme(
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

# Breakdown of CR4 diagnosis between de novo and progressed
ggplot_sim_CR4_breakdown <- df %>%
  filter(t_index >= 49) %>% 
  filter(scenario == 'Base case') %>% 
  mutate(CR4_novo = CR4 - CR4_progr) %>%
  select(t_index, t_label, CR4_novo, CR4_progr) %>%
  pivot_longer(cols = c(-t_index, -t_label)) %>%
  ggplot() +
  geom_area(mapping = aes(x = t_index, y = value, fill = name)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  labs(
    title = 'Breakdown of CR4',
    y     = 'Number of Individuals'
  ) +
  theme(
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

# Breakdown of CR4 by treatment lines and RAS status
ggplot_sim_CR4_treatmentlines <- df %>%
  filter(t_index >= 49) %>% 
  filter(scenario == 'Base case') %>% 
  select(t_index, t_label, RASwt_L1, RASwt_L2, RASwt_L3, RASwt_L4, RASmt_L1, RASmt_L2, RASmt_L3, RASmt_L4) %>%
  pivot_longer(cols = c(-t_index, -t_label)) %>%
  mutate(
    line = substr(name, 7, 8),
    type = substr(name, 1, 5)
  ) %>%
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = line, linetype = type)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  labs(
    title = 'Systemic Treatment Lines for stage IV colorectal cancer',
    y     = 'Number of Individuals'
  ) +
  theme(
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

ggplot_sim_CR4_systemic <- df %>%
  filter(t_index >= 49) %>% 
  select(
    t_index, t_label, scenario, 
    contains('_L1_'), contains('_L2_'), contains('_L3_'), contains('_L4_'),
    -contains('_SG'), -contains('_RT')
  ) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(
    treatmentline = factor(
      x = substr(name, 7, 8),
      levels = c('L1', 'L2', 'L3', 'L4'),
      labels = c('1st line systemic', '2nd line systemic', '3rd line systemic', '4th line systemic')
    ),
    treatment = factor(
      x = gsub(pattern = 'RAS.t_L._', replacement = '', x = name),
      levels = c('chemo', 'EGFR', 'chemoEGFR', 'chemoVEGF', 'NEW')
    )
  ) %>%
  group_by(t_index, t_label, scenario, treatmentline, treatment) %>% 
  summarize(value = sum(value)) %>% 
  ungroup %>% 
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment, linetype = scenario)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  labs(
    title = 'Systemic treatment lines for stage IV colorectal cancer',
    y     = 'Number of Individuals'
  ) +
  theme(
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  ) +
  facet_grid(~treatmentline) + 
  ylim(0, 600)


tbl_out <- df %>% 
  filter(t_index >= 49) %>% 
  relocate(scenario, .before = t_index)



### 5. FINAL ACTIONS ----

saveRDS(object = ggplot_incidence,    file = 'Output/ggplot_incidence 20211026.RDS')
saveRDS(object = ggplot_stage_colon,  file = 'Output/ggplot_stage_colon 20211026.RDS')
saveRDS(object = ggplot_stage_rectal, file = 'Output/ggplot_stage_rectal 20211026.RDS')

saveRDS(object = ggplot_sim_incidence,          file = 'Output/ggplot_sim_incidence 20211213.RDS')
saveRDS(object = ggplot_sim_treatment_colon,    file = 'Output/ggplot_sim_treatment_colon 20211213.RDS')
saveRDS(object = ggplot_sim_treatment_rectal,   file = 'Output/ggplot_sim_treatment_rectal 20211213.RDS')
saveRDS(object = ggplot_sim_treatment_locreg,   file = 'Output/ggplot_sim_treatment_locreg 20211213.RDS')
saveRDS(object = ggplot_sim_CR4_breakdown,      file = 'Output/ggplot_sim_CR4_breakdown 20211213.RDS')
saveRDS(object = ggplot_sim_CR4_treatmentlines, file = 'Output/ggplot_sim_CR4_treatmentlines 20211213.RDS')
saveRDS(object = ggplot_sim_CR4_systemic,       file = 'Output/ggplot_sim_CR4_systemic 20211213.RDS')

saveRDS(object = tbl_out, file = 'Output/tbl_out 20211213.RDS')
write.csv(x = tbl_out, file = 'Output/tbl_out 20211213.csv')


