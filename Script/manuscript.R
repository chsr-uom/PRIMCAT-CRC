#-----------------------------------------------------------------------------------------#
# This script contains the R code used to generate the tables and figures included        #
# in the manuscript.                                                                      #
#-----------------------------------------------------------------------------------------#
#
#
### Libraries and set up ----
library(tidyverse)
library(reshape2)
library(ggtext)
library(patchwork)


# colour palette blind-friendly 
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7","#F0E442", "#0072B2", "#D55E00", "#AB1113")


### Load results from analyses ----
df <- bind_rows(
  readRDS(file = 'Output/sim_basecase 50runs 20211213.RDS') %>% mutate(scenario = 'Base case'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario1 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 1'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario2 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 2'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario3 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 3'),
  readRDS(file = 'Output/sim_pembrolizumab_scenario4 50runs 20211213.RDS') %>% mutate(scenario = 'Scenario 4')) %>%
  group_by(scenario) %>% 
  mutate(t_index = 1:n()) %>%
  ungroup() %>% 
  rename(t_label = time_period) %>%
  relocate(t_index, t_label, .before = C1)


### Base case scenario, i.e., standard of care Australia ----

#### Incidence and prevalence from 2010 onwards ----
incidence <- df %>% 
  select(t_index, t_label, scenario, C1, C2, C3, R1, R2, R3, LRR, CR4) %>%
  filter(scenario == "Base case") %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(name = factor(x = name, levels = c('CR4', 'C1', 'C3', 'C2', 'R1', 'R3', 'R2', 'LRR'),
                       labels = c('Colorectal IV','Colon I', 'Colon III', 'Colon II',
                                  'Rectal I', 'Rectal III','Rectal II', 'Locreg. Recur'))) %>%
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = name)) +
  scale_x_continuous(breaks = seq(0, 28*4, 8), labels = df$t_label[seq(0, 28*4, 8)+1]) +
  scale_color_manual(values = cbp1) +
  # labs(title = 'Number of Patients by Disease Type and Stage',
  #   y     = 'Number of Individuals') +
  labs(y = 'Number of Individuals') +
  theme(axis.text.x  = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) +
  geom_vline(mapping = aes(xintercept = 40), linetype = 2, color = 'darkred', linewidth = 0.7) + 
  geom_text(mapping = aes(x = 55, y = 1050, label = 'End of warm-up to simulate\nprevalent population'), size = 4, color = 'darkred')


#### Overall treatment utilisation colon and rectal cancers ----
combined_data <- df %>%
  filter(t_index >= 49, scenario == 'Base case') %>% 
  select(t_index, t_label, scenario, starts_with('C'), starts_with('R')) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(cancer_type = ifelse(str_detect(name, "^C"), "Colon", "Rectal"),
         stage = case_when(
           str_detect(name, '1') ~ 'stage I',
           str_detect(name, '2') ~ 'stage II',
           str_detect(name, '3') ~ 'stage III'),
         treatment = factor(
           x = gsub(pattern = '[CR]._Tx_', replacement = '', x = name),
           levels = c('SG', 'neoSG', 'SGadj', 'neoSGadj'))) %>%
  filter(!is.na(treatment)) %>%
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment, linetype = cancer_type)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  scale_color_manual(values = cbp1) +
  labs(y = 'Number of treatments per quarter') +
  theme(axis.text.x  = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15))  +
  facet_grid(~stage)


#### Treatment lines colorectal IV ----
ggplot_crc <- df %>% 
  filter(t_index >= 49) %>% 
  relocate(scenario, .before = t_index) %>%
  select(t_index, t_label, scenario, 
         contains('_L1_'), contains('_L2_'), contains('_L3_'), contains('_L4_'),
         -contains('_SG'), -contains('_RT')) %>%
  filter(scenario == "Base case") %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(treatmentline = factor(
    x = substr(name, 7, 8),
    levels = c('L1', 'L2', 'L3', 'L4'),
    labels = c('First line', 'Second line', 'Third line', 'Fourth line')),
    treatment = factor(
      x = gsub(pattern = 'RAS.t_L._', replacement = '', x = name),
      levels = c('chemo', 'EGFR', 'chemoEGFR', 'chemoVEGF'))) %>%
  group_by(t_index, t_label, scenario, treatmentline, treatment) %>% 
  summarise(value = sum(value)) %>% 
  ungroup %>%
  filter(!is.na(treatment)) %>%
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  scale_color_manual(values = cbp1) +
  labs(y = 'Number of treatments per quarter') +
  theme(axis.text.x  = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) +
  facet_grid(~treatmentline) + 
  ylim(0, 600)

## Figure generated
figure_03 <- combined_data / ggplot_crc + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 15))



#### CRC IV de novo and progressed ----
ggplot_crc_breakdown <- df %>%
  filter(t_index >= 49) %>% 
  filter(scenario == 'Base case') %>% 
  mutate(CR4_novo = CR4 - CR4_progr) %>%
  select(t_index, t_label, CR4_novo, CR4_progr) %>%
  pivot_longer(cols = c(-t_index, -t_label)) %>%
  ggplot() +
  geom_area(mapping = aes(x = t_index, y = value, fill = name)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  labs(y = 'Number of Individuals') +
  theme(axis.text.x  = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) 


#### Treatment type numbers over the forecast period ----
total_counts <- df %>% 
  filter(t_index >= 49, scenario == 'Base case') %>% 
  select(t_index, t_label, scenario, contains('1_Tx_'), contains('2_Tx_'), contains('3_Tx_'), 
         contains('LRR_Tx_'), contains('_L1_'), contains('_L2_'), 
         contains('_L3_'), contains('_L4_')) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(cancer_type = case_when(
    str_detect(name, "LRR") ~ "LocReg Rec", 
    str_detect(name, "RAS") ~ "Colorectal",          
    str_detect(name, "^C") ~ "Colon", 
    str_detect(name, "R") ~ "Rectal"), 
    stage = case_when(
      str_detect(name, '1') & cancer_type %in% c('Colon', 'Rectal') ~ 'stage I',
      str_detect(name, '2') & cancer_type %in% c('Colon', 'Rectal') ~ 'stage II',
      str_detect(name, '3') & cancer_type %in% c('Colon', 'Rectal') ~ 'stage III',
      cancer_type %in% c('LocReg Rec') ~ 'LocReg Rec',
      cancer_type %in% c('Colorectal') ~ 'stage IV'),
    treatmentline = factor(
      x = substr(name, 7, 8),
      levels = c('L1', 'L2', 'L3', 'L4'),
      labels = c('First line', 'Second line', 'Third line', 'Fourth line'))) %>%
  mutate(treatment = ifelse(cancer_type %in% c('Colon', 'Rectal'), 
                            str_extract(name, "(neoSGadj|SGadj|neoSG|SG)"),
                            NA),
         treatment = ifelse(cancer_type %in% c('LocReg Rec'), 
                            str_extract(name, "(neoSGadj|SGadj|neoSG|SG|systemic)"),
                            treatment),
         treatment = ifelse(cancer_type %in% c('Colorectal'), 
                            str_extract(name, "(chemo(?!EGFR|VEGF)|chemoEGFR|chemoVEGF|EGFR|SGmets|SGprim|RT)"),
                            treatment)) %>%
  filter(!is.na(treatment)) %>%
  filter(!treatment %in% c("SGmets", "SGprim", "RT")) %>% # if don't want local tx for stage IV
  mutate(treatment = factor(treatment, levels = c('SG', 'neoSG', 'SGadj', 'neoSGadj',
                                                  'systemic', 'chemo', 'EGFR', 'chemoEGFR', 'chemoVEGF')))

# Create a new column 'Year'
total_counts$Year <- ifelse(between(total_counts$t_index, 49, 52), 2022,
                            ifelse(between(total_counts$t_index, 53, 56), 2023,
                                   ifelse(between(total_counts$t_index, 57, 60), 2024,
                                          ifelse(between(total_counts$t_index, 61, 64), 2025,
                                                 ifelse(between(total_counts$t_index, 65, 68), 2026, NA)))))


## Summary high level table
total_counts_by_year_type <- total_counts %>%
  group_by(Year, cancer_type) %>%
  summarise(total_patients = sum(value, na.rm = TRUE))

total_counts_wide <- total_counts_by_year_type %>%
  pivot_wider(names_from = cancer_type, values_from = total_patients) %>%
  mutate(Overall = Colon + Rectal + Colorectal) 

print(total_counts_wide) 


#### Adjunct treatment numbers in stage IV ----
total_counts_local <- df %>%
  filter(t_index >= 49, scenario == 'Base case') %>% 
  select(t_index, t_label, scenario, contains('1_Tx_'), contains('2_Tx_'), contains('3_Tx_'), 
         contains('LRR_Tx_'), contains('_L1_'), contains('_L2_'), 
         contains('_L3_'), contains('_L4_')) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(cancer_type = case_when(
    str_detect(name, "LRR") ~ "LocReg Rec", 
    str_detect(name, "RAS") ~ "Colorectal",          
    str_detect(name, "^C") ~ "Colon", 
    str_detect(name, "R") ~ "Rectal"), 
    stage = case_when(
      str_detect(name, '1') & cancer_type %in% c('Colon', 'Rectal') ~ 'stage I',
      str_detect(name, '2') & cancer_type %in% c('Colon', 'Rectal') ~ 'stage II',
      str_detect(name, '3') & cancer_type %in% c('Colon', 'Rectal') ~ 'stage III',
      cancer_type %in% c('LocReg Rec') ~ 'LocReg Rec',
      cancer_type %in% c('Colorectal') ~ 'stage IV'),
    treatmentline = factor(
      x = substr(name, 7, 8),
      levels = c('L1', 'L2', 'L3', 'L4'),
      labels = c('First line', 'Second line', 'Third line', 'Fourth line'))) %>%
  mutate(treatment = ifelse(cancer_type %in% c('Colon', 'Rectal'), 
                            str_extract(name, "(neoSGadj|SGadj|neoSG|SG)"),
                            NA),
         treatment = ifelse(cancer_type %in% c('LocReg Rec'), 
                            str_extract(name, "(neoSGadj|SGadj|neoSG|SG|systemic)"),
                            treatment),
         treatment = ifelse(cancer_type %in% c('Colorectal'), 
                            str_extract(name, "(chemo(?!EGFR|VEGF)|chemoEGFR|chemoVEGF|EGFR|SGmets|SGprim|RT)"),
                            treatment)) %>%
  filter(!is.na(treatment)) %>%
  filter(treatment %in% c("SGmets", "SGprim", "RT")) %>% # want local tx for stage IV
  mutate(treatment = factor(treatment, levels = c("SGmets", "SGprim", "RT")))

# Create a new column 'Year'
total_counts_local$Year <- ifelse(between(total_counts_local$t_index, 49, 52), 2022,
                            ifelse(between(total_counts_local$t_index, 53, 56), 2023,
                                   ifelse(between(total_counts_local$t_index, 57, 60), 2024,
                                          ifelse(between(total_counts_local$t_index, 61, 64), 2025,
                                                 ifelse(between(total_counts_local$t_index, 65, 68), 2026, NA)))))

## Appendix 3, surgery vs other treatments
# First, summarise total patients by year and cancer_type
surgery_only <- total_counts %>%
  filter(!stage %in% c("stage IV", 'LocReg Rec')) %>%
  filter(treatment %in% c("SG")) %>%
  group_by(Year, cancer_type, stage, treatment) %>%
  summarise(total_patients = round(sum(value, na.rm = TRUE))) %>%
  ungroup() %>%
  group_by(Year, cancer_type, stage) %>%
  mutate(Frequency = total_patients / sum(total_patients) * 100) 
write.csv(surgery_only, file = "Output/surgery_only_appendix3.csv")

surgery_adj <- total_counts %>%
  filter(!stage %in% c("stage IV", 'LocReg Rec')) %>%
  filter(!treatment %in% c("SG")) %>%
  group_by(Year, cancer_type, stage) %>%
  summarise(total_patients = round(sum(value, na.rm = TRUE))) %>%
  ungroup() %>%
  group_by(Year, cancer_type, stage) %>%
  mutate(Frequency = total_patients / sum(total_patients) * 100) 
write.csv(surgery_adj, file = "Output/surgery_adj_appendix3.csv")

tx_line_iv <- total_counts %>%
  filter(stage %in% c("stage IV")) %>%
  group_by(Year, treatmentline) %>%
  summarise(total_patients = round(sum(value, na.rm = TRUE))) %>%
  filter(!total_patients %in% c(0)) %>%
  ungroup() %>%
  group_by(Year) %>%
  mutate(Frequency = total_patients / sum(total_patients) * 100) 
write.csv(tx_line_iv, file = "Output/tx_line_iv_appendix3.csv")

tx_stage_iv <- total_counts %>%
  filter(stage %in% c("stage IV")) %>%
  group_by(Year, treatmentline, treatment) %>%
  summarise(total_patients = round(sum(value, na.rm = TRUE))) %>%
  ungroup() %>%
  group_by(Year) %>%
  mutate(Frequency = total_patients / sum(total_patients) * 100) 
write.csv(tx_stage_iv, file = "Output/tx_stage_iv_appendix3.csv")


summary_counts <- total_counts %>%
  group_by(Year, stage, treatment) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  group_by(stage) %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100) 

# factor stage so that it appears in the right order in the graph below
summary_counts$stage <- factor(summary_counts$stage, 
                               levels = c('stage I', "stage II", "stage III", "LocReg Rec", "stage IV"),
                               labels = c('Stage I', "Stage II", "Stage III", "LocReg Rec", "Stage IV")) 
write.csv(summary_counts, file = "Output/summary_forecast.csv")


ggplot(summary_counts, aes(fill=treatment, y=TotalCount, x=stage)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Number of treatments for 2022-2026", x = "Disease Stage", fill = "Treatment") +
  scale_fill_manual(values = cbp1) +
  theme(axis.text.x  = element_text(angle = 60),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) +
  ggtitle("Distribution of treatments across disease stages") 



#### Appendix ----
# for all treatments 
total_counts_appendix3 <- total_counts %>%
  group_by(Year, name, scenario, cancer_type, stage, treatmentline) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(total_counts_appendix3, file = "Output/appendix3_total_baseline.csv")

# for local therapies in stage IV
total_counts_local_appendix3 <- total_counts_local %>%
  group_by(Year, name, scenario, cancer_type, stage, treatmentline) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(total_counts_local_appendix3, file = "Output/appendix3_total_local_baseline.csv")

summary_counts_all <- total_counts %>%
  group_by(stage, treatment) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(summary_counts_all, file = "Output/summary_all_forecast.csv")

summary_counts_cancer <- total_counts %>%
  group_by(cancer_type, stage, treatment) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  group_by(stage) %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(summary_counts_cancer, file = "Output/summary_cancertype_forecast.csv")

summary_counts_cancer_all <- total_counts %>%
  group_by(cancer_type, stage, treatment) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(summary_counts_cancer_all, file = "Output/summary_cancertype_all_forecast.csv")

summary_tx_IV <- total_counts %>%
  filter(stage == "stage IV") %>%
  group_by(treatmentline, treatment) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  group_by(treatmentline) %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(summary_tx_IV, file = "Output/summary_txIV_forecast.csv")

summary_line_IV <- total_counts %>%
  filter(stage == "stage IV") %>%
  group_by(treatmentline) %>%
  summarise(TotalCount = round(sum(value))) %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(summary_line_IV, file = "Output/summary_lineIV_forecast.csv")


### Pembrolizumab case study ----

#### Treatment utilisation colorectal stage IV ----
t_labels <- c("2022 Q1", "2023 Q1", "2024 Q1", "2025 Q1", "2026 Q1")
pembro_df <- df %>%
  filter(t_index >= 49) %>% 
  select(t_index, t_label, scenario, 
         contains('_L1_'), contains('_L2_'), contains('_L3_'), contains('_L4_'),
         -contains('_SG'), -contains('_RT')) %>%
  filter(scenario %in% c("Base case","Scenario 2","Scenario 4")) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario)) %>%
  mutate(name = ifelse(str_detect(name, "NEW"), str_replace(name, "NEW", "pembrolizumab"), name)) %>%
  mutate(treatmentline = factor(
    x = substr(name, 7, 8),
    levels = c('L1', 'L2', 'L3', 'L4'),
    labels = c('First line', 'Second line', 'Third line', 'Fourth line')),
    treatment = factor(
      x = gsub(pattern = 'RAS.t_L._', replacement = '', x = name),
      levels = c('chemo', 'EGFR', 'chemoEGFR', 'chemoVEGF', 'pembrolizumab'))) %>%
  group_by(t_index, t_label, scenario, treatmentline, treatment) %>% 
  summarise(value = sum(value)) %>% 
  ungroup %>%
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment, linetype = scenario)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = t_labels) +
  scale_color_manual(values = cbp1) +
  labs(y = 'Number of Individuals') +
  theme(axis.text.x  = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) +
  facet_grid(~treatmentline) + 
  ylim(0, 600)

mutate(treatmentline = factor(
    x = substr(name, 7, 8),
    levels = c('L1', 'L2', 'L3', 'L4'),
    labels = c('First line', 'Second line', 'Third line', 'Fourth line')),
    treatment = factor(
      x = gsub(pattern = 'RAS.t_L._', replacement = '', x = name),
      levels = c('chemo', 'EGFR', 'chemoEGFR', 'chemoVEGF'))) %>%
  group_by(t_index, t_label, scenario, treatmentline, treatment) %>% 
  summarise(value = sum(value)) %>% 
  ungroup %>%
  filter(!is.na(treatment)) %>%
  ggplot() +
  geom_line(mapping = aes(x = t_index, y = value, color = treatment)) +
  scale_x_continuous(breaks = seq(49, 68, 4), labels = df$t_label[seq(49, 68, 4)]) +
  scale_color_manual(values = cbp1) +
  labs(y = 'Number of Individuals') +
  theme(axis.text.x  = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) +
  facet_grid(~treatmentline) + 
  ylim(0, 600)


## full table for pembrolizumab case study
total_pembro_appendix3 <- df %>%
  filter(t_index >= 49) %>% 
  select(t_index, t_label, scenario, L2, L3, L4, Pembro, Year,
         contains('_L1_'), contains('_L2_'), contains('_L3_'), contains('_L4_'),
         -contains('_SG'), -contains('_RT')) %>%
  pivot_longer(cols = c(-t_index, -t_label, -scenario, -L2, -L3, -L4, -Pembro, -Year)) %>%
  mutate(name = ifelse(str_detect(name, "NEW"), str_replace(name, "NEW", "pembrolizumab"), name)) %>%
  group_by(Year, name, scenario) %>%
  summarise(TotalCount = round(sum(value))) %>%
  ungroup() %>%
  mutate(Frequency = TotalCount / sum(TotalCount) * 100)
write.csv(total_pembro_appendix3, file = "Output/appendix3_total_pembro.csv")


progression_appendix3 <- df %>%
  group_by(scenario, Year) %>% 
  summarise(n_stageIV_progr = sum(CR4_progr),
            n_stageIV_total = sum(CR4)) %>% 
  mutate(n_stageIV_novo  = n_stageIV_total - n_stageIV_progr,
         p_stageIV_progr = n_stageIV_progr / n_stageIV_total,
         p_stageIV_novo  = n_stageIV_novo / n_stageIV_total) %>%
  arrange(Year)
write.csv(progression_appendix3, file = "Output/progression_appendix3.csv")

