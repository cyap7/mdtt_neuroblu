
#============================================================
# 
# PREDICTORS of depression from the manic/mixed episode - ALL bipolar
#
#============================================================

source("Setup_background.r")
source("fastPhenotype.r")
library(tidyverse)
library(forcats)
library(lme4)
library(viridis)
library(ggrepel)
library(broom)
library(survival)
library(tidycmprsk)
library(ggsurvfit)
library(sandwich)
library(lmtest) # for robust SEs for repeated measures
library(gee)

source("read_in_routine.r")

#======================================================
# Read-in data
#======================================================

meds_psych_bipolar_dx <- getDf("meds_psych_bipolar_dx")

bipolar_mm_7d_hosp_psychmed <- getDf("bipolar_mm_7d_hosp_psychmed")
bipolar_mm_7d_hosp_psychmeds_count <- getDf("bipolar_mm_7d_hosp_psychmeds_count")
bipolar_mm_7d_hosp_psychmeds_count.long <- getDf("bipolar_mm_7d_hosp_psychmeds_count.long")

bipolar_mm_7d_hosp_medsind_count <- getDf("bipolar_mm_7d_hosp_medsind_count")
bipolar_mm_7d_hosp_medsind_count.long <- getDf("bipolar_mm_7d_hosp_medsind_count.long")

bipolar_mm_7d_hosp_medsind_count <- getDf("bipolar_mm_7d_hosp_medsind_count")
bipolar_mm_7d_hosp_medsind_count.long <- getDf("bipolar_mm_7d_hosp_medsind_count.long")

bipolar_epall_hospital_overlap <- getDf("bipolar_epall_hospital_overlap")

psych_med_df <- read.csv("psych_med_df.csv", header = T)
psych_med_df <- rbind(psych_med_df, data.frame(drug = "trospium", broad_class = "anticholinergic", narrow_class = ""))

#============================================================
# Initial munge
#============================================================

#--------------------------------
# All individuals, hospital only
#--------------------------------
# - exclude when the number of people with a class is <= 100
cols.tmp <- unique(bipolar_mm_7d_hosp_psychmed$broad_class)
cols.tmp <- bipolar_mm_7d_hosp_psychmeds_count[c(cols.tmp, "mm_split")] %>%
  group_by(mm_split) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
  filter(mm_split != "both") %>%
  pivot_longer(!mm_split, names_to = "drug", values_to = "n") %>%
  rename(broadclass = drug)
count_broadclass_exc <- cols.tmp %>% filter(n<=200)
count_broadclass_keep <- cols.tmp %>% filter(n>200)

# - ensure at least 50 for each drug
# - exclude when the number of people with a class is <= 50
count_min_medsind <- bipolar_mm_7d_hosp_medsind_count %>%
  dplyr::select(-c("person_id", "pmmd01", "pmmd01", "diagnosis_date", "visit_start_date", "visit_end_date", "mm_split", "dep_date", "hospital_los", "polypharmacy")) %>%
  colSums() %>%
  data.frame() %>% rownames_to_column("term") %>% rename(n = ".")

count_min_medsind_mania <- bipolar_mm_7d_hosp_medsind_count %>%
  filter(mm_split == "mania") %>%
  dplyr::select(-c("person_id", "pmmd01", "diagnosis_date", "visit_start_date", "visit_end_date", "mm_split", "dep_date", "hospital_los", "polypharmacy")) %>%
  colSums() %>%
  data.frame() %>% rownames_to_column("meds") %>% rename(n = ".") %>%
  mutate(mm_split = "mania")
count_min_medsind_mixed <- bipolar_mm_7d_hosp_medsind_count %>%
  filter(mm_split == "mixed") %>%
  dplyr::select(-c("person_id", "pmmd01", "diagnosis_date", "visit_start_date", "visit_end_date", "mm_split", "dep_date", "hospital_los", "polypharmacy")) %>%
  colSums() %>%
  data.frame() %>% rownames_to_column("meds") %>% rename(n = ".") %>%
  mutate(mm_split = "mixed")
count_min_medsind <- rbind(count_min_medsind_mania, count_min_medsind_mixed) %>% arrange(mm_split, -n)
saveDf(count_min_medsind, "count_min_medsind", overwrite = T)
count_medsind_exc <- count_min_medsind %>% filter(n<=200)
count_medsind_keep <- count_min_medsind %>% filter(n>200) %>%
  left_join(., psych_med_df %>% rename(meds = drug)) %>%
  inner_join(., count_broadclass_keep %>% rename(broad_class = broadclass) %>% dplyr::select(-n))

# UpSet plot for coprescription
bipolar_mm_7d_hosp_psychmeds_count_forupset_mania <- bipolar_mm_7d_hosp_psychmeds_count %>% 
  filter(mm_split == "mania") %>%
  dplyr::select(person_id, antidepressant, antiepileptic, lithium, antipsychotic, antihistamine, benzodiazepine) %>%
  group_by(antidepressant, antiepileptic, lithium, antipsychotic, antihistamine, benzodiazepine) %>%
  summarise(breakdown=n())
saveDf(bipolar_mm_7d_hosp_psychmeds_count_forupset_mania[1:50,], "bipolar_mm_7d_hosp_psychmeds_count_forupset_mania_part1", overwrite = T)
saveDf(bipolar_mm_7d_hosp_psychmeds_count_forupset_mania[50:nrow(bipolar_mm_7d_hosp_psychmeds_count_forupset_mania),], "bipolar_mm_7d_hosp_psychmeds_count_forupset_mania_part2", overwrite = T)

bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed <- bipolar_mm_7d_hosp_psychmeds_count %>% 
  filter(mm_split == "mixed") %>%
  dplyr::select(person_id, antidepressant, antiepileptic, lithium, antipsychotic, antihistamine, benzodiazepine) %>%
  group_by(antidepressant, antiepileptic, lithium, antipsychotic, antihistamine, benzodiazepine) %>%
  summarise(breakdown=n())
saveDf(bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed[1:50,], "bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed_part1", overwrite = T)
saveDf(bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed[50:nrow(bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed),], "bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed_part2", overwrite = T)

# GO TO LOCAL
library(UpSetR)
# - mania
toupset1 <- read.csv("~/Downloads/bipolar_mm_7d_hosp_psychmeds_count_forupset_mania_part1")
toupset2 <- read.csv("~/Downloads/bipolar_mm_7d_hosp_psychmeds_count_forupset_mania_part2")
toupset <- full_join(toupset1, toupset2)
upset_data <- toupset %>%
  # Expand each row according to the breakdown count
  uncount(breakdown)
pdf("~/OneDrive - Nexus365/Research/Projects/Postmania_depression/EHR_Neuroblu/Write-up/figures/bipolar_mm_7d_hosp_psychmeds_count_forupset_mania.pdf", onefile=FALSE)
upset(upset_data, 
      sets = c("antidepressant", "antiepileptic", "lithium", "antipsychotic", "antihistamine", "benzodiazepine"),
      order.by = "freq")
dev.off()
# - mixed
toupset1 <- read.csv("~/Downloads/bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed_part1")
toupset2 <- read.csv("~/Downloads/bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed_part2")
toupset <- full_join(toupset1, toupset2)
upset_data <- toupset %>%
  # Expand each row according to the breakdown count
  uncount(breakdown)
pdf("~/OneDrive - Nexus365/Research/Projects/Postmania_depression/EHR_Neuroblu/Write-up/figures/bipolar_mm_7d_hosp_psychmeds_count_forupset_mixed.pdf", onefile=FALSE)
upset(upset_data, 
      sets = c("antidepressant", "antiepileptic", "lithium", "antipsychotic", "antihistamine", "benzodiazepine"),
      order.by = "freq")
dev.off()

#--------------------------------
# Munge input dataframe
#--------------------------------

pmmd_ep <- getDf("pmmd_ep")
pmmd_ep.tmp <- pmmd_ep %>% 
  dplyr::select(person_id, mm_episode, dep_episode, diagnosis_date, diagnosis_date_end, dep_date) %>%
  rename(dep_episode_update = dep_episode)
mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
mds1m_group <- getDf("mds1m_group")
observation_period <- getDf("observation_period")

bipolar_dx7d_hosp_mmsplit_cov_csh_180d.prep0 <- mania_dep_dates_toepisodes_qc %>%
  #left_join(., firstmm_end) %>%
  # take the min obs time in case of multiple 
  left_join(., observation_period %>% dplyr::select(person_id, observation_period_end_date)) %>%
  mutate(obs_time_diff = as.numeric(difftime(observation_period_end_date, diagnosis_date, units = "days"))) %>% 
  group_by(across(c(-obs_time_diff))) %>% slice_min(abs(obs_time_diff)) %>% dplyr::select(-obs_time_diff) %>% ungroup %>%
  # add dep from additional case-finding
  left_join(., pmmd_ep.tmp) %>%
  mutate(dep_episode = ifelse(!is.na(dep_episode_update), dep_episode_update, dep_episode)) %>% dplyr::select(-dep_episode_update) %>%
  # censor past 180 days post-MM
  mutate(postmm_180d_date = if_else(state == "mm", as.Date(diagnosis_date_end) + days(180), NA)) %>%
  # calculate start and end times
  group_by(person_id) %>%
  arrange(mm_episode, dep_episode) %>%
  mutate(start_csh = 0,
  end_csh = as.numeric(difftime(
    pmin(
      if_else(dep_date > diagnosis_date_end, dep_date, as.Date(NA)),
      if_else(observation_period_end_date > diagnosis_date_end, observation_period_end_date, as.Date(NA)),
      if_else(postmm_180d_date > diagnosis_date_end, postmm_180d_date, as.Date(NA)),
      if_else(lead(diagnosis_date) > diagnosis_date_end, lead(diagnosis_date), as.Date(NA)),
      na.rm = TRUE),
    diagnosis_date_end,
    units = "days"
  ))) %>%
  # add status/outcome information per episode. Censor to 180d
  mutate(next_dx = ifelse(lead(state) == "mm", lead(diagnosis_date), NA)) %>%
  mutate(min_date = pmin(dep_date, observation_period_end_date, postmm_180d_date, next_dx, na.rm = TRUE)) %>%
  mutate(status.tmp = case_when(
    dep_date == min_date ~ "dep",
    observation_period_end_date == min_date ~ "censor",
    postmm_180d_date == min_date ~ "180d",
    next_dx == min_date ~ "next_dx")) %>%
  mutate(status0 = case_when(status.tmp == "dep" ~ 1,
                             status.tmp == "next_dx" ~ 2,
                             .default = 0)) %>%
  mutate(status0 = factor(status0, levels = c(0, 1, 2))) %>%
  ungroup() %>%
  # filter for mm episodes
  filter(state == "mm") %>%
  # add hospital data, rapid switch group, filter for meds that are sufficiently commonly prescribed
  left_join(., bipolar_epall_hospital_overlap %>% filter(state == "mm") %>%
                dplyr::select(person_id, diagnosis_date, diagnosis_date_end) %>%
                mutate(hospital = 1)) %>%
  mutate(hospital01 = ifelse(is.na(hospital), 0, hospital)) %>% dplyr::select(-hospital) %>%
  left_join(., mds1m_group) %>%
  # add covariates
  left_join(., mm_split) %>%
  inner_join(., demog_cov %>% dplyr::select(person_id, birth_year, starts_with("cov_"))) %>%
  mutate(cov_age_ep = as.numeric(difftime(diagnosis_date, as.Date(birth_year), units = "days"))/365) %>% dplyr::select(-birth_year) %>%
  filter(cov_gender != "No matching concept") %>%
  mutate(cov_edu2 = case_when(cov_edu %in% c("8th grade or less", "Some high school, but did not graduate", "High school diploma or GED") ~ "School",
                              cov_edu %in% c("Some college", "Associate degree (e.g., AA, AS)", "Bachelor's degree", "Graduate degree") ~ "College",
                              .default = cov_edu)) %>%
  mutate(cov_employ2 = case_when(cov_employ %in% c("Full-time employment", "Part-time employment") ~ "Employed",
                                cov_employ %in% c("Not in labor force", "Unemployed", "Retired", "Permanently unable to perform work activities due to medical condition", "On work experience") ~ "NotEmployed",
                                .default = cov_employ)) %>%
  mutate(cov_race2 = ifelse(cov_race %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander"), "AmIndian_Alaskan_PacIslander", cov_race))

cov_state_n_meds.csh180d <- bipolar_dx7d_hosp_mmsplit_cov_csh_180d.prep0 %>% 
  left_join(., mm_split) %>%
  inner_join(., bipolar_mm_7d_hosp_medsind_count.long) %>%
  group_by(person_id, mm_split, cov_state) %>% summarise(n=n()) %>% 
  group_by(mm_split, cov_state) %>% summarise(n=n()) %>%
  mutate(cov_state2 = ifelse(n <= 50, "grouped", cov_state)) %>% ungroup

# Broad class
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep <- bipolar_dx7d_hosp_mmsplit_cov_csh_180d.prep0 %>%
  left_join(., cov_state_n_meds.csh180d %>% dplyr::select(cov_state, cov_state2, mm_split)) %>%
  rename(cov_state0 = cov_state, cov_state = cov_state2) %>%
  # add medication data
  left_join(., bipolar_mm_7d_hosp_psychmeds_count.long) %>%
  inner_join(., count_broadclass_keep) %>%
  filter(!broadclass %in% c("vitamin", "analgesia"))
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep", overwrite = T)

# - wide format for joint analysis
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  dplyr::select(-c("n")) %>%
  pivot_wider(names_from = broadclass, values_from = drug01)
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep", overwrite = T)

# Individual medications
bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep <- bipolar_dx7d_hosp_mmsplit_cov_csh_180d.prep0 %>%
  left_join(., cov_state_n_meds.csh180d %>% dplyr::select(cov_state, cov_state2, mm_split)) %>%
  rename(cov_state0 = cov_state, cov_state = cov_state2) %>%
  # add medication data
  left_join(., bipolar_mm_7d_hosp_medsind_count.long) %>%
  inner_join(., count_medsind_keep) %>%
  left_join(., days_since_hosp)
saveDf(bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep, "bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep", overwrite = T)

# - wide format for joint analysis
bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d_wide.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  dplyr::select(-c("n", "broad_class", "narrow_class")) %>%
  pivot_wider(names_from = meds, values_from = drug01)
saveDf(bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d_wide.prep, "bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d_wide.prep", overwrite = T)

#============================================================
# EDA plots
#============================================================

# Medications
# 1. polypharmacy 
# 2. number of mm episodes
# 3. demographics
# combine 1+3 in a summary plot

# Polypharmacy
bipolar_mm_7d_hosp_psychmeds_count_polypharmacy <- bipolar_mm_7d_hosp_psychmeds_count.long %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end) %>%
  mutate(polypharmacy_class = sum(drug01), broadclass = broadclass) %>% ungroup() %>%
  group_by(broadclass, polypharmacy_class, drug01) %>%
  summarise(n=n()) %>% ungroup() %>%
  filter(drug01 == 1) %>%
  group_by(broadclass) %>%
  mutate(n_total = sum(n)) %>%
  mutate(prop = n/n_total) %>%
  mutate(cumulative_prop = cumsum(prop)) %>% ungroup() 

bipolar_mm_7d_hosp_psychmeds_count_polypharmacy.plot <- ggplot(bipolar_mm_7d_hosp_psychmeds_count_polypharmacy, aes(x = polypharmacy_class, y = cumulative_prop, colour = broadclass)) +
  geom_point() + geom_line() + 
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_colour_manual(values = c(millerstone)) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  xlab("Number of co-prescribed drug classes") +
  ylab("Cumulative prop")

# Number of episodes at first prescription coinciding with first 7 days of episode
bipolar_dx7d_dx_medsbroadclass_mmsplit_cov_csh_180d.prep <- getDf("bipolar_dx7d_dx_medsbroadclass_mmsplit_cov_csh_180d.prep")
bipolar_dx7d_dx_medsbroadclass_mmsplit_cov_csh_180d_epnumber <- bipolar_dx7d_dx_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(person_id, broadclass, drug01) %>%
  summarise(mm_episode=min(mm_episode)) %>% ungroup() %>%
  group_by(broadclass, mm_episode, drug01) %>%
  summarise(n=n()) %>% ungroup() %>%
  filter(drug01 == 1) %>%
  group_by(broadclass) %>%
  mutate(n_total = sum(n)) %>%
  mutate(prop = n/n_total) %>%
  mutate(cumulative_prop = cumsum(prop)) %>% ungroup() 

bipolar_dx7d_dx_medsbroadclass_mmsplit_cov_csh_180d_epnumber.plot <- ggplot(bipolar_dx7d_dx_medsbroadclass_mmsplit_cov_csh_180d_epnumber, aes(x = mm_episode, y = cumulative_prop, colour = broadclass)) +
  geom_point() + geom_line() + 
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c(millerstone)) +
  coord_cartesian(xlim = c(1, 3)) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  xlab("Nth manic/mixed episodes at 1st record of prescription") +
  ylab("Cumulative prop of all manic/mixed episodes")

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep <- getDf("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep")

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_demo_plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"))

# - Age at ep
meds_hosp_age_ep.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_demo_plot %>% 
  filter(drug01 == 1) %>%
  ggplot(., aes(x = cov_age_ep, colour = broadclass)) +
  geom_density() +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_colour_manual(values = c(millerstone))

# - Years EHR
meds_hosp_years_ehr.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_demo_plot %>% 
  filter(drug01 == 1) %>%
  ggplot(., aes(x = cov_years_ehr, colour = broadclass)) +
  geom_density() +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_colour_manual(values = c(millerstone))

# - Sex
meds_hosp_sex <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(broadclass, drug01, cov_gender) %>%
  summarise(n=n()) %>% filter(drug01 == 1) %>%
  ungroup() %>%
  group_by(broadclass) %>%
  mutate(drug1sum = sum(n)) %>%
  mutate(prop = n/drug1sum)

meds_hosp_sex.plot <- meds_hosp_sex %>% 
  ggplot(., aes(x = cov_gender, y = prop, fill = broadclass)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_fill_manual(values = c(millerstone))
ggsave("tmp.png")

# - cov_race2
meds_hosp_cov_race <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(broadclass, drug01, cov_race2) %>%
  summarise(n=n()) %>% filter(drug01 == 1) %>%
  ungroup() %>%
  group_by(broadclass) %>%
  mutate(drug1sum = sum(n)) %>%
  mutate(prop = n/drug1sum)

meds_hosp_cov_race.plot <- meds_hosp_cov_race %>% filter(drug01 == 1) %>%
  ggplot(., aes(x = cov_race2, y = prop, fill = broadclass)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values = c(millerstone))
  #facet_wrap(~broadclass, scale = "free_y")

# - cov_edu2
meds_hosp_cov_edu <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(broadclass, drug01, cov_edu2) %>%
  summarise(n=n()) %>% filter(drug01 == 1) %>%
  ungroup() %>%
  group_by(broadclass) %>%
  mutate(drug1sum = sum(n)) %>%
  mutate(prop = n/drug1sum)

meds_hosp_cov_edu.plot <- meds_hosp_cov_edu %>% filter(drug01 == 1) %>%
  ggplot(., aes(x = cov_edu2, y = prop, fill = broadclass)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_fill_manual(values = c(millerstone))

# - cov_employ2
meds_hosp_cov_employ <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(broadclass, drug01, cov_employ2) %>%
  summarise(n=n()) %>% filter(drug01 == 1) %>%
  ungroup() %>%
  group_by(broadclass) %>%
  mutate(drug1sum = sum(n)) %>%
  mutate(prop = n/drug1sum)

meds_hosp_cov_employ.plot <- meds_hosp_cov_employ %>% filter(drug01 == 1) %>%
  ggplot(., aes(x = cov_employ2, y = prop, fill = broadclass)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_fill_manual(values = c(millerstone))

# - cov_rel
meds_hosp_cov_rel <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(broadclass, drug01, cov_rel) %>%
  summarise(n=n()) %>% filter(drug01 == 1) %>%
  ungroup() %>%
  group_by(broadclass) %>%
  mutate(drug1sum = sum(n)) %>%
  mutate(prop = n/drug1sum)

meds_hosp_cov_rel.plot <- meds_hosp_cov_rel %>% filter(drug01 == 1) %>%
  ggplot(., aes(x = cov_rel, y = prop, fill = broadclass)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_fill_manual(values = c(millerstone))

# - cov_state
meds_hosp_cov_state <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  group_by(broadclass, drug01, cov_state) %>%
  summarise(n=n()) %>% filter(drug01 == 1) %>%
  ungroup() %>%
  group_by(broadclass) %>%
  mutate(drug1sum = sum(n)) %>%
  mutate(prop = n/drug1sum)

meds_hosp_cov_state.plot <- meds_hosp_cov_state %>% filter(drug01 == 1) %>%
  ggplot(., aes(x = cov_state, y = prop, fill = broadclass)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  #theme(legend.position = "bottom") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values = c(millerstone))

# AGGREGATE plots
library(patchwork)
bipolar_mm_7d_hosp_psychmeds_count_polypharmacy.plot /
(meds_hosp_age_ep.plot + meds_hosp_years_ehr.plot) /
(meds_hosp_sex.plot + meds_hosp_cov_race.plot) /
(meds_hosp_cov_edu.plot + meds_hosp_cov_employ.plot) /
(meds_hosp_cov_rel.plot + meds_hosp_cov_state.plot) +
  #plot_layout(heights = c(c(1.25, 1.25))) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect")
ggsave("meds_hosp_demog_cov.png", height = 13, width = 11)

#------------------------------------------------------------
# Calculate days until follow-up in outpatient after hospital
#------------------------------------------------------------

bipolar_visit <- getDf("bipolar_visit") %>% dplyr::select(-visit_occurrence_id) %>% distinct

# Get outpatient visits for people in the MM hospital analysis (ie. no "Pharmacy visit" or "Laboratory Visit") - do this later as sensitivity
bipolar_visit_mm_nonhosp <- bipolar_visit %>% 
  #filter(visit_type == "Outpatient") %>%
  filter(person_id %in% bipolar_mm_7d_hosp_psychmeds_count.long$person_id) %>%
  rename(visit_start_date_all = visit_start_date, visit_end_date_all = visit_end_date)

# Join with hospital df ##### decide or not whether to keep this days_since_hosp
bipolar_mm_7d_hosp_psychmeds_count_visit <- bipolar_mm_7d_hosp_psychmeds_count.long %>%
  full_join(., bipolar_visit_mm_nonhosp, relationship = "many-to-many") %>%
  mutate(days_since_hosp = case_when(dep_date <= visit_end_date ~ 0, # if the dep date comes before the end of the initial hospital visit, follow-up time = 0
                                    visit_start_date_all > visit_end_date ~ as.numeric(difftime(visit_start_date_all, visit_end_date, units = "days"))
                                    )) %>%
  mutate(next_visit = case_when(dep_date <= visit_end_date ~ "during_hosp", # if the dep date comes before the end of the initial hospital visit, follow-up time = 0
                                    visit_start_date_all > visit_end_date ~ "after_hosp"
                                    )) %>%
  filter(days_since_hosp >= 0)

# Look at spread of earliest visit date per medication class
bipolar_mm_7d_hosp_psychmeds_count_visit_followup_first <- bipolar_mm_7d_hosp_psychmeds_count_visit %>%
  #filter(days_since_hosp <= 180) %>%
  dplyr::select(-c("visit_start_date_all", "visit_end_date_all", "visit_type")) %>%
  group_by(across(c(-days_since_hosp))) %>%
  slice_min(days_since_hosp) %>% ungroup

hosp_followup_first_per_meds <- bipolar_mm_7d_hosp_psychmeds_count_visit_followup_first %>%
  group_by(mm_split, broadclass, drug01) %>%
  summarise(days_since_hosp_median = median(days_since_hosp),
        days_since_hosp_mean = mean(days_since_hosp)) %>%
  filter(broadclass %in% c("antipsychotic", "antiepileptic", "lithium", "antihistamine", "benzodiazepine", "antidepressant"))

bipolar_mm_7d_hosp_psychmeds_count_visit_followup_first.gg <- bipolar_mm_7d_hosp_psychmeds_count_visit_followup_first %>% 
  filter(broadclass %in% c("antipsychotic", "antiepileptic", "lithium", "antihistamine", "benzodiazepine", "antidepressant")) %>%
  mutate(drug_prescribed_during_mania_mixed = ifelse(drug01 == 1, "yes", "no")) %>%
  filter(days_since_hosp > 0) %>%
  ggplot(., aes(x = broadclass, y = days_since_hosp, fill = drug_prescribed_during_mania_mixed)) +
  #geom_violin() + 
  geom_boxplot() + #geom_jitter(alpha = 0.5, colour = "gray50", width = 0.25)  +
  theme_bw() + theme(legend.position = "bottom") +
  scale_fill_manual(values = millerstone[c(1,3,7,9,8,2)]) +
  coord_cartesian(ylim=c(0,100)) +  
  ylab("Days to first follow-up (exclude switches in hospital)") +
  facet_wrap(~ mm_split, nrow = 2)

# Look at sum of visits within 180 days
bipolar_mm_7d_hosp_psychmeds_count_visit_followup_visitspp <- bipolar_mm_7d_hosp_psychmeds_count_visit %>%
  #filter(days_since_hosp <= 180) %>%
  dplyr::select(-c("visit_start_date_all", "visit_end_date_all", "visit_type")) %>%
  group_by(across(c(-days_since_hosp))) %>%
  summarise(visits_sum = n())

hosp_followup_visitspp_per_meds <- bipolar_mm_7d_hosp_psychmeds_count_visit_followup_visitspp %>%
  group_by(mm_split, broadclass, drug01) %>%
  summarise(visits_sum_median = median(visits_sum),
        visits_sum_mean = mean(visits_sum)) %>%
  filter(broadclass %in% c("antipsychotic", "antiepileptic", "lithium", "antihistamine", "benzodiazepine", "antidepressant"))

bipolar_mm_7d_hosp_psychmeds_count_visit_followup_visitspp.gg <- bipolar_mm_7d_hosp_psychmeds_count_visit_followup_visitspp %>% 
  filter(broadclass %in% c("antipsychotic", "antiepileptic", "lithium", "antihistamine", "benzodiazepine", "antidepressant")) %>%
  mutate(drug_prescribed_during_mania_mixed = ifelse(drug01 == 1, "yes", "no")) %>%
  ggplot(., aes(x = broadclass, y = visits_sum, fill = drug_prescribed_during_mania_mixed)) +
  geom_boxplot() + #geom_jitter(alpha = 0.5, colour = "gray50", width = 0.25)  +
  theme_bw() + theme(legend.position = "bottom") +
  scale_fill_manual(values = millerstone[c(1,3,7,9,8,2)]) +
  ylab("Number of visits within 180 days after hospital visit") +
  coord_cartesian(ylim=c(0,180)) +  
  facet_wrap(~ mm_split, nrow = 2)

# Plot
library(patchwork)
bipolar_mm_7d_hosp_psychmeds_count_visit_followup_first.gg / bipolar_mm_7d_hosp_psychmeds_count_visit_followup_visitspp.gg +
  plot_annotation(tag_levels = 'a')
ggsave("bipolar_mm_7d_hosp_psychmeds_count_visit_followup_first_visitspp.png", height = 11, width = 7)

# Get days_since_hosp covariate
days_since_hosp <- bipolar_mm_7d_hosp_psychmeds_count_visit_followup_first %>%
  ungroup() %>%
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, days_since_hosp) %>%
  distinct()
saveDf(days_since_hosp, "days_since_hosp", overwrite = T)

missing <- anti_join(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep, days_since_hosp) # remaining do not have follow-up evidence

#============================================================
# Modelling
#============================================================

#------------------------------------------------------------
# Cox with repeated risk intervals ANALYSIS: impact of drug classes on time until depression (instead of 0/1 within 6 months)
#------------------------------------------------------------

#------------------------------------------------------------
# Check PH assumptions
#------------------------------------------------------------

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_coxph.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
 mutate(statuscoxph = case_when(status0 == 0 ~ 1,
                             status0 == 1 ~ 2,
                             status0 == 2 ~ 1))

ph_jointmeds_mania_hosp <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_coxph.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania", 0, antidepressant)) %>%
  filter(mm_split == "mania") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  coxph(Surv(end_csh, statuscoxph) ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)
data.frame(cox.zph(ph_jointmeds_mania_hosp)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mania") %>% mutate(sample = "bipolar_hosp")
 
ph_jointmeds_mixed_hosp <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_coxph.prep %>%
  filter(mm_split == "mixed") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("PA", "CA", "FL", "CO", "DE", "FL", "IL", "SC", "MI", "IL"), "grouped", cov_state)) %>% # due to collinearity issue
  coxph(Surv(end_csh, statuscoxph) ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)
data.frame(cox.zph(ph_jointmeds_mixed_hosp)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mixed") %>% mutate(sample = "bipolar_hosp")

cox_zph_hospbipolar_meds_summary <- rbind(
  data.frame(cox.zph(ph_jointmeds_mania_hosp)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mania"),
  data.frame(cox.zph(ph_jointmeds_mixed_hosp)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mixed")) %>% 
  mutate(sample = "bipolar_hosp")
saveDf(cox_zph_hospbipolar_meds_summary, "cox_zph_hospbipolar_meds_summary", overwrite = T)

#============================================================
# Cox with repeated risk intervals
#============================================================

#---------------------------
# Broad class: 
#---------------------------

bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broadclass == "antidepressant", 0, drug01)) %>%
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  nest(.by = c("broadclass", "mm_split")) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status0 == 1) ~ drug01 + strata(MDS1m_ever) + cluster(person_id), data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% mutate(p_adj = p.value*8) %>% ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("HR=", sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA))
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh", overwrite = T)

# - covariates (NOT censored for 180d)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broadclass == "antidepressant", 0, drug01)) %>%
  left_join(., days_since_hosp) %>%
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  nest(.by = c("broadclass", "mm_split")) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status0 == 1) ~ drug01 + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term == "drug01", p.value*8, NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("HR=", sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA))
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh", overwrite = T)

# - joint model with broad_class lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine (as sig in the cov analysis)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania", 0, antidepressant)) %>%
  left_join(., days_since_hosp) %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status0 == 1) ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("HR=", sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA))
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh", overwrite = T)

# - joint model + account for time to follow-up, with broad_class lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine (as sig in the cov analysis)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania", 0, antidepressant)) %>%
  left_join(., days_since_hosp) %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status0 == 1) ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + cluster(person_id) + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("HR=", sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA))
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp", overwrite = T)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp_part1 <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp[1:50,]
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp_part2 <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp[51:nrow(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh),]
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp_part1, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp_part1", overwrite = T)
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp_part2, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp_part2", overwrite = T)

# PLOT
# - include sample size per drug
class_counts_cov <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  group_by(mm_split, broadclass) %>%
  summarise(drug_count = sum(drug01)) %>%
  mutate(drug_name_count = paste0(broadclass, "(", drug_count, ")"))

class_counts_joint <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  group_by(mm_split, broadclass) %>%
  summarise(drug_count = sum(drug01)) %>%
  mutate(drug_name_count = paste0(broadclass, "(", drug_count, ")")) %>%
  rename(term = broadclass)

class_counts_joint_days_since_hosp <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  inner_join(., days_since_hosp) %>%
  group_by(mm_split, broadclass) %>%
  summarise(drug_count = sum(drug01)) %>%
  mutate(drug_name_count = paste0(broadclass, "(", drug_count, ")")) %>%
  rename(term = broadclass)

# Plot - NO COVARIATES
bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh %>%
  inner_join(., class_counts_cov) %>% filter(!is.na(estimate_HR)) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(term == "drug01") %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_HR), y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("HR estimate (no covariates)") +
  scale_y_log10() + coord_flip(ylim = c(0.6,1.2)) +
  facet_wrap(~ mm_split, scales = "free")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh.png", bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh.plot, width = 7, height = 3)

# Plot - COVARIATE model
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh %>%
  left_join(., class_counts_cov) %>% filter(!is.na(estimate_HR)) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(term == "drug01") %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_HR), y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("HR estimate (adj. for covariates)") +
  scale_y_log10() + coord_flip(ylim = c(0.6,1.2)) +
  facet_wrap(~ mm_split, scales = "free")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh.png", bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh.plot, width = 7, height = 3)

# plot - JOINT model
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh %>%
  left_join(., class_counts_joint) %>% filter(!is.na(estimate_HR)) %>%
  mutate(name = "broad class") %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_HR), y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("HR estimate (joint model, adj. for covariates)") +
  scale_y_log10() + coord_flip(ylim = c(0.6,1.2)) +
  facet_wrap(~ mm_split, scales = "free")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh.png", bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh.plot, width = 7, height = 3)

# plot - JOINT model + days_since_hosp
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp %>%
  left_join(., class_counts_joint_days_since_hosp) %>% filter(!is.na(estimate_HR)) %>%
  mutate(name = "broad class") %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_HR), y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("HR estimate (joint model, adj. for covariates + follow-up)") +
  scale_y_log10() + coord_flip(ylim = c(0.6,1.2)) +
  facet_wrap(~ mm_split, scales = "free")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp.png", bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp.plot, width = 7, height = 3)

library(patchwork)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_csh.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_cov_csh_dayssincehosp.plot + 
  plot_annotation(tag_levels = 'a')
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_csh.png", height = 10, width = 8.5)

#---------------------------
# Individual medications: 
#---------------------------

# - here want to "explode" the classes into constituent medications
bipolar_dx7d_hosp_medsind_mmsplit_nocov_csh <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broad_class == "antidepressant", 0, drug01)) %>%
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  mutate(meds = ifelse(broad_class %in% c("antiepileptic", "lithium", "antihistamine"), meds, broad_class)) %>%
  dplyr::select(-c("broad_class", "narrow_class", "n")) %>%
  group_by(across(c(-drug01))) %>% # group_by all except drug01 --> want to have drug01 == 1 if any drug within that class is given
  summarise(drug01 = max(drug01)) %>% 
  ungroup() %>%
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  nest(.by = c("meds", "mm_split")) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status0 == 1) ~ drug01 + cluster(person_id), data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>%
  mutate(p_adj = ifelse(term == "drug01", p.value*6, NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste("HR=", sprintf('%.2f', estimate_HR), ", ", 
    "SE=", sprintf('%.2f', std.error), ", ", 
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA))
saveDf(bipolar_dx7d_hosp_medsind_mmsplit_nocov_csh, "bipolar_dx7d_hosp_medsind_mmsplit_nocov_csh", overwrite = T)

# - covariates (censored for 180d)
bipolar_dx7d_hosp_medsind_mmsplit_cov_csh <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broad_class == "antidepressant", 0, drug01)) %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  mutate(meds = ifelse(broad_class %in% c("antiepileptic", "lithium", "antihistamine"), meds, broad_class)) %>%
  dplyr::select(-c("broad_class", "narrow_class", "n")) %>%
  group_by(across(c(-drug01))) %>% # group_by all except drug01 --> want to have drug01 == 1 if any drug within that class is given
  summarise(drug01 = max(drug01)) %>% 
  ungroup() %>%
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  nest(.by = c("meds", "mm_split")) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status0 == 1) ~ drug01 + cluster(person_id) + hospital_los + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term == "drug01", p.value*6, NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  #mutate(p_FDR = ifelse(term == "drug01", p.adjust(p.value, "fdr"), NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste("HR=", sprintf('%.2f', estimate_HR), ", ", 
    "SE=", sprintf('%.2f', std.error), ", ", 
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA))
saveDf(bipolar_dx7d_hosp_medsind_mmsplit_cov_csh, "bipolar_dx7d_hosp_medsind_mmsplit_cov_csh", overwrite = T)

# JOINT
bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d_wide_ms_ah.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broad_class == "antidepressant", 0, drug01)) %>%
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  mutate(meds = ifelse(broad_class %in% c("antiepileptic", "lithium", "antihistamine"), meds, broad_class)) %>%
  dplyr::select(-c("broad_class", "narrow_class", "n")) %>%
  group_by(across(c(-drug01))) %>% # group_by all except drug01 --> want to have drug01 == 1 if any drug within that class is given
  summarise(drug01 = max(drug01)) %>% 
  ungroup() %>%
  pivot_wider(names_from = "meds", values_from = "drug01")

bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d_wide_ms_ah.prep %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  filter(mm_split == "mania") %>%
  coxph(Surv(start_csh, end_csh, status0 == 1) ~ lithium + valproate + lamotrigine + hydroxyzine + diphenhydramine + antipsychotic + antidepressant + benzodiazepine + cluster(person_id) + cov_age_ep + cov_age_entry + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .) %>% 
  #coxph(Surv(start_csh, end_csh, status0 == 1) ~ lithium + valproate + lamotrigine + haloperidol + aripiprazole + olanzapine + quetiapine + risperidone + ziprasidone + chlorpromazine + paliperidone + sertraline + trazodone + strata(MDS1m_ever) + cluster(person_id) + cov_age_ep + cov_age_entry + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .) %>% 
  tidy() %>% 
  mutate(mm_split = "mania") %>%
  mutate(p_adj = ifelse(term %in% c("lithium", "valproate", "lamotrigine", "hydroxyzine", "diphenhydramine"), p.value*5, NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  #mutate(p_FDR = ifelse(term %in% c("lithium", "valproate", "lamotrigine", "hydroxyzine", "diphenhydramine"), p.adjust(p.value, "fdr"), NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("HR=", sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA))

bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d_wide_ms_ah.prep %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "CA", "FL"), "grouped", cov_state)) %>% # due to collinearity issue
  mutate(end_csh = ifelse(end_csh == start_csh, end_csh + 1, end_csh)) %>%
  filter(mm_split == "mixed") %>%
  coxph(Surv(start_csh, end_csh, status0 == 1) ~ lithium + valproate + lamotrigine + hydroxyzine + diphenhydramine + antipsychotic + antidepressant + benzodiazepine + cluster(person_id) + cov_age_ep + cov_age_entry + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .) %>% 
  #coxph(Surv(start_csh, end_csh, status0 == 1) ~ lithium + valproate + lamotrigine + haloperidol + aripiprazole + olanzapine + quetiapine + risperidone + ziprasidone + chlorpromazine + paliperidone + sertraline + trazodone + strata(MDS1m_ever) + cluster(person_id) + cov_age_ep + cov_age_entry + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .) %>% 
  tidy() %>% 
  mutate(mm_split = "mixed") %>%
  mutate(p_adj = ifelse(term %in% c("lithium", "valproate", "lamotrigine", "hydroxyzine", "diphenhydramine"), p.value*5, NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  #mutate(p_FDR = ifelse(term %in% c("lithium", "valproate", "lamotrigine", "hydroxyzine", "diphenhydramine"), p.adjust(p.value, "fdr"), NA)) %>% ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)), 
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("HR=", sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA))

bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh <- rbind(bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh_mania, bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh_mixed)
saveDf(bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh, "bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh", overwrite = T)

# Plot - NO COVARIATES
bipolar_dx7d_hosp_medsind_mmsplit_nocov_csh.plot <- bipolar_dx7d_hosp_medsind_mmsplit_nocov_csh %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(term == "drug01") %>%
  left_join(., psych_med_df %>% rename(meds = drug)) %>%
  mutate(broad_class = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  #filter(broad_class %in% c("moodstabiliser", "antihistamine")) %>%
  mutate(broad_class = factor(broad_class, levels = c("moodstabiliser", "antihistamine"))) %>%  # Apply custom ordering
  ggplot(aes(x = fct_reorder(meds, estimate_HR), y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("HR estimate (unadj. for covariates)") +
  coord_flip() +
  facet_grid(broad_class ~ mm_split, scales = "free_y", space = "free_y")
ggsave("bipolar_dx7d_hosp_medsind_mmsplit_nocov_csh.png", bipolar_dx7d_hosp_medsind_mmsplit_nocov_csh.plot, width = 7, height = 3)

# Plot - COVARIATE model
bipolar_dx7d_hosp_medsind_mmsplit_cov_csh.plot <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(term == "drug01") %>%
  left_join(., psych_med_df %>% rename(meds = drug)) %>%
  mutate(broad_class = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  #filter(broad_class %in% c("moodstabiliser", "antihistamine")) %>%
  mutate(broad_class = factor(broad_class, levels = c("moodstabiliser", "antihistamine"))) %>%  # Apply custom ordering
  ggplot(aes(x = fct_reorder(meds, estimate_HR), y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("HR estimate (adj. for covariates)") +
  coord_flip() +
  facet_grid(broad_class ~ mm_split, scales = "free_y", space = "free_y")
ggsave("bipolar_dx7d_hosp_medsind_mmsplit_cov_csh.png", bipolar_dx7d_hosp_medsind_mmsplit_cov_csh.plot, width = 7, height = 3)

# - JOINT model
bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh.plot <- bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh %>%
  mutate(name = "broad class") %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  left_join(., psych_med_df %>% rename(term = drug)) %>%
  mutate(broad_class = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  filter(broad_class %in% c("moodstabiliser", "antipsychotic", "antidepressant", "antihistamine")) %>%
  mutate(broad_class = factor(broad_class, levels = c("moodstabiliser", "antipsychotic", "antidepressant", "antihistamine"))) %>%  # Apply custom ordering
  ggplot(aes(x = fct_reorder(term, estimate_HR), y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("HR estimate (joint model, adj. for covariates)") +
  coord_flip() +
  facet_grid(broad_class ~ mm_split, scales = "free_y", space = "free_y")
ggsave("bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh.png", bipolar_dx7d_hosp_medsind_mmsplit_joint_cov_csh.plot, width = 7, height = 3)

#============================================================
# Modified poisson regression models for RR
#============================================================

#---------------------------
# RR, bipolar group, broad class: 
#---------------------------

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d.prep %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) #%>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) # due to collinearity issue 

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep <- getDf("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep")
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) #%>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) # due to collinearity issue 

# Modified poisson regression models for RR
# - NOCOV
bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_glmrr <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d.prep %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broadclass == "antidepressant", 0, drug01)) %>%  
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  nest(.by = c("broadclass", "mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ drug01, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term == "drug01", p.value*8, NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "no covariates")
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_glmrr, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_glmrr", overwrite = T)

# - COV
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d.prep %>%
  filter(broadclass %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broadclass == "antidepressant", 0, drug01)) %>%  
  nest(.by = c("broadclass", "mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ drug01 + hospital_los + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term == "drug01", p.value*8, NA)) %>% 
  #mutate(p_FDR = ifelse(term == "drug01", p.adjust(p.value, "fdr"), NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "covariates")
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr", overwrite = T)

# - JOINT, not adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>%  
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% 
  #mutate(p_FDR = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.adjust(p.value, "bonf"), NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint")
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr", overwrite = T)

# - JOINT, adjust for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>%  
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + days_since_hosp + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% 
  #mutate(p_FDR = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.adjust(p.value, "bonf"), NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint")
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup", overwrite = T)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_part1 <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup[1:50,]
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_part2 <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup[51:nrow(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup),]
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_part1, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_part1", overwrite = T)
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_part2, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_part2", overwrite = T)

# Plot
# - include sample size per drug
class_counts_cov <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d.prep %>%
  group_by(mm_split, broadclass) %>%
  summarise(drug_count = sum(drug01)) %>%
  mutate(drug_name_count = paste0(broadclass, "(", drug_count, ")"))

class_counts_joint <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d.prep %>%
  group_by(mm_split, broadclass) %>%
  summarise(drug_count = sum(drug01)) %>%
  mutate(drug_name_count = paste0(broadclass, "(", drug_count, ")")) %>%
  rename(term = broadclass)

class_counts_joint_days_since_hosp <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d.prep %>%
  inner_join(., days_since_hosp) %>%
  group_by(mm_split, broadclass) %>%
  summarise(drug_count = sum(drug01)) %>%
  mutate(drug_name_count = paste0(broadclass, "(", drug_count, ")")) %>%
  rename(term = broadclass)

# - NO COVARIATES
bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_glmrr.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_glmrr %>%
  filter(term == "drug01") %>%
  left_join(., class_counts_cov) %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate for Dep switch in 1m (unadj. for covariates)") +
  scale_y_log10() + coord_flip(ylim = c(0.7,1.1)) +
  facet_wrap(~ mm_split, scales = "free")

# - COVARIATE model
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr %>%
  filter(term == "drug01") %>%
  left_join(., class_counts_cov) %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate for Dep switch in 1m (adj. for covariates)") +
  scale_y_log10() + coord_flip(ylim = c(0.7,1.1)) +
  facet_wrap(~ mm_split, scales = "free")

# - JOINT, not adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  left_join(., class_counts_joint) %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for covariates)") +
  scale_y_log10() + coord_flip(ylim = c(0.7,1.1)) +
  facet_wrap(~ mm_split, scales = "free")

# - JOINT, adjust for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  left_join(., class_counts_joint_days_since_hosp) %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for cov + follow-up)") +
  scale_y_log10() + coord_flip(ylim = c(0.7,1.1)) +
  facet_wrap(~ mm_split, scales = "free")

library(patchwork)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_nocov_glmrr.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup.plot + 
  plot_annotation(tag_levels = 'a')
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_glmrr.png", height = 10, width = 8.5)

# FINAL main text figure
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  left_join(., class_counts_joint_days_since_hosp) %>%
  ggplot(aes(x = fct_reorder(drug_name_count, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for cov + follow-up)") +
  scale_y_log10() + coord_flip(ylim = c(0.75,1.1)) +
  facet_wrap(~ mm_split, scales = "free")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup.svg", bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup.plot, width = 9, height = 3)
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup.png", bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup.plot, width = 9, height = 3)

#---------------------------
# RR, bipolar group, individual drugs
#---------------------------

bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) %>%
  left_join(., days_since_hosp)

# - mania
bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_ms_ap_mania.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  mutate(meds = ifelse(broad_class %in% c("antiepileptic", "lithium", "antipsychotic", "benzodiazepine"), meds, broad_class)) %>%
  dplyr::select(-c("broad_class", "narrow_class", "n")) %>%
  group_by(across(c(-drug01))) %>% # group_by all except drug01 --> want to have drug01 == 1 if any drug within that class is given
  summarise(drug01 = max(drug01)) %>% 
  ungroup()

bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_wide_ms_ap_mania.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_ms_ap_mania.prep %>%
  pivot_wider(names_from = "meds", values_from = "drug01") %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) %>%
  filter(cov_gender != "No matching concept")

meds_list_mania.tmp <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "benzodiazepine")) %>%
  filter(drug01 == 1) %>%
  group_by(meds, mm_split) %>%
  summarise(n=n()) %>% arrange(mm_split) %>% filter(mm_split == "mania")
meds_list_mania <- meds_list_mania.tmp %>% pull(meds) %>% unique

# - mixed
bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_ms_ap_mixed.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  mutate(meds = ifelse(broad_class %in% c("antiepileptic", "lithium", "antipsychotic"), meds, broad_class)) %>%
  dplyr::select(-c("broad_class", "narrow_class", "n")) %>%
  group_by(across(c(-drug01))) %>% # group_by all except drug01 --> want to have drug01 == 1 if any drug within that class is given
  summarise(drug01 = max(drug01)) %>% 
  ungroup()

bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_wide_ms_ap_mixed.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_ms_ap_mixed.prep %>%
  pivot_wider(names_from = "meds", values_from = "drug01") %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) %>%
  filter(cov_gender != "No matching concept")

meds_list_mixed.tmp <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic")) %>%
  filter(drug01 == 1) %>%
  group_by(meds, mm_split) %>%
  summarise(n=n()) %>% arrange(mm_split) %>% filter(mm_split == "mixed")
meds_list_mixed <- meds_list_mixed.tmp %>% pull(meds) %>% unique

# - JOINT
# - mania
bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_wide_ms_ap_mania.prep %>%
  left_join(., days_since_hosp) %>%
  filter(mm_split == "mania") %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + valproate + lamotrigine + aripiprazole + chlorpromazine + haloperidol + olanzapine + paliperidone + quetiapine + risperidone + ziprasidone + alprazolam + chlordiazepoxide + clonazepam + diazepam + lorazepam + midazolam + antihistamine + days_since_hosp + hospital_los + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  #mutate(p_adj = ifelse(term %in% meds_list_mania, p.value*5, NA)) %>% 
  #mutate(p_FDR = ifelse(!term %in% c("hospital_los", "cov_age_ep", "cov_age_ehr", "cov_years_ehr", "cov_gender", "cov_state", "cov_rel", "cov_edu2", "cov_employ2"), p.adjust(p.value, "fdr"), NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint")

# - mixed
bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d_wide_ms_ap_mixed.prep %>%
  left_join(., days_since_hosp) %>%
  filter(mm_split == "mixed") %>%
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + valproate + lamotrigine + aripiprazole + haloperidol + lurasidone + olanzapine + quetiapine + risperidone + ziprasidone + antidepressant + antihistamine + benzodiazepine + days_since_hosp + hospital_los + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  #mutate(p_adj = ifelse(term %in% meds_list_mixed, p.value*5, NA)) %>% 
  #mutate(p_FDR = ifelse(!term %in% c("hospital_los", "cov_age_ep", "cov_age_ehr", "cov_years_ehr", "cov_gender", "cov_state", "cov_rel", "cov_edu2", "cov_employ2"), p.adjust(p.value, "fdr"), NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint")

# - Plot JOINT model
# - include per drug count
drug_counts <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  group_by(mm_split, meds) %>%
  summarise(drug_count = sum(drug01)) %>%
  mutate(drug_name_count = paste0(meds, "(", drug_count, ")")) %>%
  rename(term = meds)

# - munge for plot
bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr <- full_join(bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mania, bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mixed) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(!grepl("cov|Intercept|hospital_los|polypharmacy", term)) %>%
  left_join(., psych_med_df %>% rename(term = drug)) %>%
  mutate(broad_class = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  filter(broad_class %in% c("moodstabiliser", "antipsychotic", "benzodiazepine")) %>%
  mutate(broad_class = factor(broad_class, levels = c("moodstabiliser", "antipsychotic", "antidepressant", "benzodiazepine", "antihistamine"))) %>% # %>%  # Apply custom ordering
  left_join(., drug_counts)

bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mania.plot <- bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr %>%
  filter(mm_split == "mania") %>%
  ggplot(., aes(x = fct_reorder(drug_name_count, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  #geom_label_repel(size = 3.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  #expand_limits(x= c(0, 6)) +  
  xlab("Variable") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for covariates)") +
  coord_flip() + scale_y_log10(breaks = c(0.7, 0.80, 1, 1.3)) +
  facet_grid(broad_class ~ mm_split, scales = "free_y", space = "free_y")

bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mixed.plot <- bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr %>%
  filter(mm_split == "mixed") %>%
  ggplot(., aes(x = fct_reorder(drug_name_count, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  #geom_label_repel(size = 3.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  #expand_limits(x= c(0, 6)) +  
  xlab("Variable") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for covariates)") +
  coord_flip() + scale_y_log10(breaks = c(0.7, 0.80, 1, 1.3)) +
  facet_grid(broad_class ~ mm_split, scales = "free_y", space = "free_y")

library(patchwork)
bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mania.plot + bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr_mixed.plot +
  plot_annotation(tag_levels = 'a') + plot_layout(axes = 'collect')
ggsave("bipolar_dx7d_hosp_medsind_mmsplit_joint_glmrr.png", width = 8.5, height = 7)

#====================================================
# Sensitivity analyses
#====================================================

#---------------------------
# Narrow ICD depression definition
#---------------------------

mm0_mmenddepstart0_dep0_toexclude <- getDf("mm0_mmenddepstart0_dep0_toexclude")
bipolar_epall_rm_qc <- getDf("bipolar_epall_rm_qc")

# Identify which episodes should switch from MDS1m == 1 to MDS1m == 0
mania_dep_dates_toepisodes_qc_defbroad_noad <- getDf("mania_dep_dates_toepisodes_qc_defbroad_noad")
mania_dep_dates_toepisodes_qc_defbroad_MDS1m_noad <- mania_dep_dates_toepisodes_qc_defbroad_noad %>%
  mutate(MDS1m_noad = ifelse(timetodep_months <= 1, 1, 0)) %>%
  mutate(MDS1m_noad = ifelse(is.na(MDS1m_noad), 0, MDS1m_noad)) %>%
  dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, mm_episode, dep_episode, MDS1m_noad)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide_MDS1m_noad.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  left_join(., mania_dep_dates_toepisodes_qc_defbroad_MDS1m_noad) %>%
  mutate(MDS1m_noad = ifelse(is.na(MDS1m_noad), 0, MDS1m_noad))

# - JOINT, adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_MDS1m_noad_glmrr <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide_MDS1m_noad.prep %>%
  left_join(., days_since_hosp) %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, NA, antidepressant)) %>% 
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m_noad ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender, #cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2
    family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint")

# - JOINT, adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_MDS1m_noad_glmrr.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_MDS1m_noad_glmrr %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(term, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint, no antidepressant definition)") +
  coord_flip() +
  facet_grid(name ~ mm_split, scales = "free_x")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_MDS1m_noad_glmrr.png", bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_MDS1m_noad_glmrr.plot, width = 7, height = 3, dpi = 600)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide_MDS1m_noad.prep %>% group_by(mm_split, lithium, MDS1m) %>% summarise(n=n())

#---------------------------
# State prescribing differences
#---------------------------

cov_state_keep_analysis <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  group_by(mm_split, cov_state) %>%
  summarise(cov_state_n=n()) %>%
  filter(cov_state_n>=20) %>%
  mutate(cov_state = ifelse(cov_state == "No matching concept", "N/A", cov_state)) %>%
  mutate(cov_state_new = paste0(cov_state, "(", cov_state_n, ")"))

# JOINT model stratifying by state, not adjusting for follow-up (excluded covariates cov_employ2, cov_race2, cov_edu2)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>% 
  inner_join(., cov_state_keep_analysis) %>% # ensure that there are enough observations per comparison
  nest(.by = c("mm_split", "cov_state_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_rel, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% 
  #mutate(p_FDR = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.adjust(p.value, "bonf"), NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint")

# Plots
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state_mania.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state %>%
  filter(mm_split == "mania") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(cov_state_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for cov)") +
  scale_y_log10() +
  coord_flip() +
  facet_grid(term ~ mm_split)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state_mixed.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state %>%
  filter(mm_split == "mixed") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(cov_state_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for cov)") +
  scale_y_log10() + coord_flip() +
  facet_grid(term ~ mm_split)

library(patchwork)
(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state_mania.plot | bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state_mixed.plot) +
  plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_state.png", height = 11, width = 10)

#---------------------------
# Race prescribing differences
#---------------------------

cov_race_keep_analysis <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  group_by(mm_split, cov_race2) %>%
  summarise(cov_race_n=n()) %>%
  filter(cov_race_n>=20) %>%
  mutate(cov_race = ifelse(cov_race2 == "No matching concept", "N/A", cov_race2)) %>%
  mutate(cov_race_new = paste0(cov_race, "(", cov_race_n, ")")) %>% dplyr::select(-cov_race)

# JOINT model stratifying by state, not adjusting for follow-up (excluded covariates cov_employ2, cov_race2, cov_edu2)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>% 
  inner_join(., cov_race_keep_analysis) %>% # ensure that there are enough observations per comparison
  nest(.by = c("mm_split", "cov_race_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + days_since_hosp + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_rel + cov_state, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) %>% 
  mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% 
  #mutate(p_FDR = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.adjust(p.value, "bonf"), NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint")

# Plots
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race_mania.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race %>%
  filter(mm_split == "mania") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(cov_race_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for cov)") +
  scale_y_log10() + coord_flip(ylim = c(0.5, 2)) +
  facet_grid(term ~ mm_split)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race_mixed.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race %>%
  filter(mm_split == "mixed") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(cov_race_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for cov)") +
  scale_y_log10() + coord_flip(ylim = c(0.5, 2)) +
  facet_grid(term ~ mm_split)

library(patchwork)
(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race_mania.plot | bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race_mixed.plot) +
  plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_cov_race.png", height = 8, width = 10)

#---------------------------
# Medication start date #####
#---------------------------
# Check how many manic/mixed episodes have had lithium prescription in the preceding 30 days

meds_psych_bipolar_dx <- getDf("meds_psych_bipolar_dx")

# BROAD CLASS
# Get continuous distribution of days before diagnosis the medication was last prescribed
timesince_lastprescription_broadclass0 <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  filter(drug01 == 1) %>% filter(!meds == "topiramate") %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  #mutate(days_lastmed_from_diagnosis7d = case_when(start_date < diagnosis_date & end_date > diagnosis_date ~ 0,
  #                            .default = as.numeric(difftime(end_date, diagnosis_date, units = "days")))) %>% # relative to end date: -ive means med end_date came before
  mutate(med_check_end = diagnosis_date-7) %>% # consider the window up to 7 days BEFORE the manic/mixed episode code
  mutate(days_lastmed_from_diagnosis7d = case_when(start_date < med_check_end & end_date > med_check_end ~ 0, # when the medication prescription period overlaps with the mm diagnosis date, make 0 ie. they were prescribed the medication contemporaneously
                              .default = as.numeric(difftime(end_date, med_check_end, units = "days")))) %>% # when was the medication last prescribed (end_date) relative to the relative to the mm diagnosis (and up to 7 days before): -ive means med end_date came before
  mutate(days_lastmed_from_diagnosis7d = ifelse(days_lastmed_from_diagnosis7d > 0, -99999, days_lastmed_from_diagnosis7d)) %>% # assign extreme value and then crop out
  #filter(days_lastmed_from_diagnosis7d <= 0) %>% # want medications to have been prescribed before diagnosis
  mutate(days_fromdiagnosis7d_to_lastmed = days_lastmed_from_diagnosis7d * -1)

timesince_lastprescription_broadclass_separateantiepileptic <- timesince_lastprescription_broadclass0 %>%
  mutate(broad_class = ifelse(broad_class == "antiepileptic", meds, broad_class)) %>% # separate out lamotrigine and valproate
  group_by(person_id, mm_episode, dep_episode, diagnosis_date, diagnosis_date_end, mm_split, broad_class) %>%
  summarise(days_fromdiagnosis7d_to_lastmed = min(days_fromdiagnosis7d_to_lastmed)) %>% # ie. the least-negative value
  ungroup() %>%
  mutate(drug01_no90dbefore = ifelse(days_fromdiagnosis7d_to_lastmed > 90, 1, 2))

timesince_lastprescription_broadclass_separateantiepileptic.plot <- ggplot(timesince_lastprescription_broadclass_separateantiepileptic %>% mutate(drug = broad_class), aes(x = days_fromdiagnosis7d_to_lastmed, colour = drug)) +
  stat_ecdf() +
  theme_bw() + theme(legend.position = "bottom") +
  scale_colour_manual(values = millerstone) +
  xlab("Days before MM start period since last prescription") + # start period = MM diagnosis -7 days
  ylab("Cumulative proportion") +
  coord_cartesian(xlim = c(0,90)) +
  facet_wrap(~mm_split, scales = "free_x")

# Sensitivity analysis excluding people from joint model with prior prescription in the past 90 days
# - munge to get table per episode of who has had that broad_class drug in the previous 90 days --> will exclude
timesince_lastprescription_broadclass <- timesince_lastprescription_broadclass0 %>%
  mutate(drug01 = ifelse(mm_split == "mania" & broad_class == "antidepressant", 0, drug01)) %>%
  group_by(person_id, mm_episode, dep_episode, diagnosis_date, diagnosis_date_end, mm_split, broad_class) %>%
  summarise(days_fromdiagnosis7d_to_lastmed = min(days_fromdiagnosis7d_to_lastmed)) %>% # ie. the least-negative value
  ungroup() %>%
  mutate(drug01_no90dbefore = ifelse(days_fromdiagnosis7d_to_lastmed > 90, "start", "90d")) %>%
  rename(broadclass = broad_class)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_drug01_no90dbefore.prep <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>% # first remake the drug01 column dividing up the antiepileptics
  left_join(., days_since_hosp) %>%  
  mutate(broad_class = ifelse(broad_class == "antiepileptic", meds, broad_class)) %>%
  dplyr::select(-c("narrow_class", "n")) %>%
  group_by(across(c(-meds, -drug01))) %>% 
  summarise(drug01 = max(drug01)) %>%
  ungroup() %>%
  filter(broad_class %in% c("lithium", "valproate", "lamotrigine", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  left_join(., timesince_lastprescription_broadclass_separateantiepileptic) %>%
  mutate(drug01_no90dbefore = case_when(drug01 == 1 & days_fromdiagnosis7d_to_lastmed <= 90 ~ "90d",
                                        drug01 == 1 & days_fromdiagnosis7d_to_lastmed > 90 ~ "start",
                                        TRUE ~ "0")) %>%
  dplyr::select(-c("days_fromdiagnosis7d_to_lastmed", "drug01")) %>%
  pivot_wider(names_from = broad_class, values_from = drug01_no90dbefore) %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) #%>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) # due to collinearity issue 

# - model
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr_mania <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_drug01_no90dbefore.prep %>%
  filter(mm_split == "mania") %>%
  left_join(., days_since_hosp) %>%
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + valproate + lamotrigine + antipsychotic + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) 

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr_mixed <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_drug01_no90dbefore.prep %>%
  filter(mm_split == "mixed") %>%
  left_join(., days_since_hosp) %>% 
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ lithium + valproate + lamotrigine + antipsychotic + antidepressant + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) 

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr <- full_join(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr_mania, bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr_mixed) %>%
  group_by(mm_split) %>% 
  #mutate(p_adj = ifelse(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"), p.value*6, NA)) %>% 
  ungroup() %>% # treat mania and mixed as separate tests
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value <= 0.05, paste("RR=", sprintf('%.2f', estimate_RR), ",", 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 1), sep = ""), NA)) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "joint") 

# - JOINT, not adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr %>%
  #mutate(name = "broad_class") %>%
  filter(grepl("start", term) | grepl("90d", term)) %>%
  mutate(name = ifelse(grepl("start", term), "from MM start period", "within 90d before MM start period")) %>%
  ggplot(aes(x = term, estimate_RR, y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, adj. for covariates)") +
  coord_flip() +
  facet_grid(name ~ mm_split, scales = "free")

# PRETTY PLOTS
library(patchwork)
timesince_lastprescription_broadclass_separateantiepileptic.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_drugno90dbefore_glmrr.plot +
  plot_annotation(tag_levels = 'a') + plot_layout(heights = unit(c(0.5, 1), 'null'))
ggsave("drugno90dbefore.png", height = 10, width = 8)

#-----------------------
# - plot single cumulative incidence curves
#-----------------------
# - for aggregated plot
lithium_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mania" & meds == "lithium") %>%
  mutate(drug01 = ifelse(drug01 == 1, "lithium", "no lithium")) %>%
  cuminc(Surv(end_csh, status0) ~ drug01, data = .) %>%
  ggcuminc() +
  xlim(0,180) + ylim(0,1) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[8], "gray20")) +
  scale_fill_manual(values = c(millerstone[8], "gray20")) +
  theme(legend.position = "none") +
  ggtitle("Mania-to-Dep ~ Lithium") +
  add_censor_mark() +
  add_confidence_interval()
valproate_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mania" & meds == "valproate") %>%
  mutate(drug01 = ifelse(drug01 == 1, "valproate", "no valproate")) %>%
  cuminc(Surv(end_csh, status0) ~ drug01, data = .) %>%
  ggcuminc() +
  xlim(0,180) + ylim(0,1) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c("gray20", millerstone[10])) +
  scale_fill_manual(values = c("gray20", millerstone[10])) +
  theme(legend.position = "none") +
  ggtitle("Mania-to-Dep ~ Valproate") +
  add_censor_mark() +
  add_confidence_interval()
lamotrigine_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mania" & meds == "lamotrigine") %>%
  mutate(drug01 = ifelse(drug01 == 1, "lamotrigine", "no lamotrigine")) %>%
  cuminc(Surv(end_csh, status0) ~ drug01, data = .) %>%
  ggcuminc() +
  xlim(0,180) + ylim(0,1) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[10], "gray20")) +
  scale_fill_manual(values = c(millerstone[10], "gray20")) +
  theme(legend.position = "none") +
  ggtitle("Mania-to-Dep ~ Lamotrigine") +
  add_censor_mark() +
  add_confidence_interval()

lithium_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mixed" & meds == "lithium") %>%
  mutate(drug01 = ifelse(drug01 == 1, "lithium", "no lithium")) %>%
  cuminc(Surv(end_csh, status0) ~ drug01, data = .) %>%
  ggcuminc() +
  xlim(0,180) + ylim(0,1) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[8], "gray20")) +
  scale_fill_manual(values = c(millerstone[8], "gray20")) +
  theme(legend.position = "none") +
  ggtitle("Mixed-to-Dep ~ Lithium") +
  add_censor_mark() +
  add_confidence_interval()
valproate_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mixed" & meds == "valproate") %>%
  mutate(drug01 = ifelse(drug01 == 1, "valproate", "no valproate")) %>%
  cuminc(Surv(end_csh, status0) ~ drug01, data = .) %>%
  ggcuminc() +
  xlim(0,180) + ylim(0,1) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c("gray20", millerstone[8])) +
  scale_fill_manual(values = c("gray20", millerstone[8])) +
  theme(legend.position = "none") +
  ggtitle("Mixed-to-Dep ~ Valproate") +
  add_censor_mark() +
  add_confidence_interval()
lamotrigine_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mixed" & meds == "lamotrigine") %>%
  mutate(drug01 = ifelse(drug01 == 1, "lamotrigine", "no lamotrigine")) %>%
  cuminc(Surv(end_csh, status0) ~ drug01, data = .) %>%
  ggcuminc() +
  xlim(0,180) + ylim(0,1) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[2], "gray20")) +
  scale_fill_manual(values = c(millerstone[2], "gray20")) +
  theme(legend.position = "none") +
  ggtitle("Mixed-to-Dep ~ Lamotrigine") +
  add_censor_mark() +
  add_confidence_interval()

ggarrange(lithium_mania, valproate_mania, lamotrigine_mania,
  lithium_mixed, valproate_mixed, lamotrigine_mixed,
  nrow = 2, ncol = 3)
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_maniamixed_allsig_cuminc_180d.png", width = 10, height = 10)

#-----------------------
# - plot joint cumulative incidence curves
#-----------------------

ms_ap_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mania") %>%
  mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  ##filter(mm_split == "mania" & broad_class %in% c("lithium", "antiepileptic", "antipsychotic")) %>%
  mutate(broad_class1 = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, broad_class1) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "broad_class1", values_from = "sum") %>%
  mutate(med_status = case_when(moodstabiliser == 0 & antipsychotic == 0 ~ "neither",
         moodstabiliser >= 1 & antipsychotic == 0 ~ "MS-only",
         moodstabiliser == 0 & antipsychotic >= 1 ~ "1.AP",
         moodstabiliser >= 1 & antipsychotic >= 1 ~ "MS+AP"))
ms_ap_mania.cuminc <- ms_ap_mania %>%
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  theme(legend.text=element_text(size=20)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mania-to-Dep: MS (any)") +
  add_censor_mark() +
  add_confidence_interval()

ms_ap_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mixed") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
#  filter(mm_split == "mixed" & broad_class %in% c("lithium", "antiepileptic", "antipsychotic")) %>%
  mutate(broad_class1 = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, broad_class1) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "broad_class1", values_from = "sum") %>%
  mutate(med_status = case_when(moodstabiliser == 0 & antipsychotic == 0 ~ "neither",
         moodstabiliser >= 1 & antipsychotic == 0 ~ "MS-only",
         moodstabiliser == 0 & antipsychotic >= 1 ~ "1.AP",
         moodstabiliser >= 1 & antipsychotic >= 1 ~ "MS+AP"))
ms_ap_mixed.cuminc <- ms_ap_mixed %>% 
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mixed-to-Dep: MS (any)") +
  add_censor_mark() +
  add_confidence_interval()

lithium_ap_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mania") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
#  filter(mm_split == "mania" & broad_class %in% c("lithium", "antipsychotic")) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, broad_class) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "broad_class", values_from = "sum") %>%
  filter(antiepileptic == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lithium == 0 & antipsychotic == 0 ~ "neither",
         lithium >= 1 & antipsychotic == 0 ~ "MS-only",
         lithium == 0 & antipsychotic >= 1 ~ "1.AP",
         lithium >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lithium_ap_mania.cuminc <- lithium_ap_mania %>%
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mania-to-Dep: LITHIUM") +
  add_censor_mark() +
  add_confidence_interval()

lithium_ap_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mixed") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
#  filter(mm_split == "mixed" & broad_class %in% c("lithium", "antipsychotic")) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, broad_class) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "broad_class", values_from = "sum") %>%
  filter(antiepileptic == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lithium == 0 & antipsychotic == 0 ~ "neither",
         lithium >= 1 & antipsychotic == 0 ~ "MS-only",
         lithium == 0 & antipsychotic >= 1 ~ "1.AP",
         lithium >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lithium_ap_mixed.cuminc <- lithium_ap_mixed %>%
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mixed-to-Dep: LITHIUM") +
  add_censor_mark() +
  add_confidence_interval()

valproate_ap_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mania") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
#  filter((mm_split == "mania" & meds == "valproate") | (mm_split == "mania" & broad_class == "antipsychotic")) %>%
  mutate(cov_state = ifelse(cov_state %in% c("FL", "PA", "IL", "TN", "SC", "VA"), "grouped", cov_state)) %>% # due to singularities with model
  mutate(meds = case_when(meds == "valproate" ~ "valproate", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & lamotrigine == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(valproate == 0 & antipsychotic == 0 ~ "neither",
         valproate >= 1 & antipsychotic == 0 ~ "MS-only",
         valproate == 0 & antipsychotic >= 1 ~ "1.AP",
         valproate >= 1 & antipsychotic >= 1 ~ "MS+AP"))
valproate_ap_mania.cuminc <- valproate_ap_mania %>%
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mania-to-Dep: VALP") +
  add_censor_mark() +
  add_confidence_interval()

valproate_ap_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mixed") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
#  filter((mm_split == "mixed" & meds == "valproate") | (mm_split == "mixed" & broad_class == "antipsychotic")) %>%
  mutate(meds = case_when(meds == "valproate" ~ "valproate", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & lamotrigine == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(valproate == 0 & antipsychotic == 0 ~ "neither",
         valproate >= 1 & antipsychotic == 0 ~ "MS-only",
         valproate == 0 & antipsychotic >= 1 ~ "1.AP",
         valproate >= 1 & antipsychotic >= 1 ~ "MS+AP"))
valproate_ap_mixed.cuminc <- valproate_ap_mixed %>%
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mixed-to-Dep: VALP") +
  add_censor_mark() +
  add_confidence_interval()

lamotrigine_ap_mania <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mania") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  mutate(meds = case_when(meds == "lamotrigine" ~ "lamotrigine", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & valproate == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lamotrigine == 0 & antipsychotic == 0 ~ "neither",
         lamotrigine >= 1 & antipsychotic == 0 ~ "MS-only",
         lamotrigine == 0 & antipsychotic >= 1 ~ "1.AP",
         lamotrigine >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lamotrigine_ap_mania.cuminc <- lamotrigine_ap_mania %>%
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mania-to-Dep: LAM") +
  add_censor_mark() +
  add_confidence_interval()

lamotrigine_ap_mixed <- bipolar_dx7d_hosp_medsind_mmsplit_cov_csh_180d.prep %>%
  filter(mm_split == "mixed") %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  mutate(meds = case_when(meds == "lamotrigine" ~ "lamotrigine", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, start_csh, end_csh, status0, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  ungroup() %>% mutate(end_csh = ifelse(end_csh == start_csh, end_csh+1, end_csh)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & valproate == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lamotrigine == 0 & antipsychotic == 0 ~ "neither",
         lamotrigine >= 1 & antipsychotic == 0 ~ "MS-only",
         lamotrigine == 0 & antipsychotic >= 1 ~ "1.AP",
         lamotrigine >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lamotrigine_ap_mixed.cuminc <- lamotrigine_ap_mixed %>%
  cuminc(Surv(end_csh, status0) ~ med_status, data = .) %>%
  ggcuminc() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,1)) +
  labs(x = "Days", y = "Cumulative Incidence") +
  scale_color_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  scale_fill_manual(values = c(millerstone[c(8,1,7)], "gray20")) +
  ggtitle("Mixed-to-Dep: LAM") +
  add_censor_mark() +
  add_confidence_interval()

# Plots of cumulative incidence curves
library(ggpubr)

# - for main text
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_maniamixed_maintext_cuminc_180d.plot <- ggarrange(
  ms_ap_mania.cuminc, lithium_ap_mania.cuminc, valproate_ap_mania.cuminc, lamotrigine_ap_mania.cuminc, 
  ms_ap_mixed.cuminc, lithium_ap_mixed.cuminc, valproate_ap_mixed.cuminc, lamotrigine_ap_mixed.cuminc,
  nrow = 2, ncol = 4, common.legend = TRUE)
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_maniamixed_maintext_cuminc_180d.png", width = 14.5, height = 7.5, bg="white")

# put all dfs together --> barplot to show sample size included in MS/AP combo analysis
ms_ap_aggregate_mania_mixed <- rbind(ms_ap_mania %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mania", analysis = "MS = Any"),
  ms_ap_mixed %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mixed", analysis = "MS = Any"),
  lithium_ap_mania %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mania", analysis = "MS = Lithium"),
  lithium_ap_mixed %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mixed", analysis = "MS = Lithium"),
  valproate_ap_mania %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mania", analysis = "MS = Valproate"),
  valproate_ap_mixed %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mixed", analysis = "MS = Valproate"),
  lamotrigine_ap_mania %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mania", analysis = "MS = Lamotrigine"),
  lamotrigine_ap_mixed %>% ungroup %>% dplyr::select(med_status) %>% mutate(mm_split = "mixed", analysis = "MS = Lamotrigine"))

ms_ap_aggregate_mania_mixed_counts <- ms_ap_aggregate_mania_mixed %>% 
  group_by(med_status, analysis, mm_split) %>% summarise(n=n()) 

ms_ap_aggregate_mania_mixed_counts %>%
  pivot_wider(id_cols = "analysis", names_from = "med_status", values_from = "n")

ggplot(ms_ap_aggregate_mania_mixed_counts, aes(x = med_status, y = n, fill = med_status)) +
  geom_col(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(millerstone[8], millerstone[7], millerstone[1], "gray40")) +
  facet_grid(analysis ~ mm_split)
ggsave("ms_ap_aggregate_mania_mixed_counts.png")

#---------------------
# forest plots to compare cumulative incidence curves
#---------------------

# CSH approach
# - Mania: ms vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_ms_ap.csh <-  coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = ms_ap_mania)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_ms_ap.csh_coef <- cbind(mm_split = "Mania-to-Dep", analysis = "Any Mood Stabiliser", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_ms_ap.csh))
# - Mixed: ms vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_ms_ap.csh <- coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = ms_ap_mania)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_ms_ap.csh_coef <- cbind(mm_split = "Mixed-to-Dep", analysis = "Any Mood Stabiliser", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_ms_ap.csh))

# - Mania: lithium vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lithium_ap.csh <- coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = lithium_ap_mania)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lithium_ap.csh_coef <- cbind(mm_split = "Mania-to-Dep", analysis = "Lithium", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lithium_ap.csh))
# - Mixed: lithium vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lithium_ap.csh <- coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = lithium_ap_mixed %>% mutate(cov_state = ifelse(cov_state %in% c("CO", "DE", "MI"), "grouped", cov_state)))
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lithium_ap.csh_coef <- cbind(mm_split = "Mixed-to-Dep", analysis = "Lithium", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lithium_ap.csh))

# - Mania: valproate vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_valproate_ap.csh <- coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = valproate_ap_mania)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_valproate_ap.csh_coef <- cbind(mm_split = "Mania-to-Dep", analysis = "Valproate", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_valproate_ap.csh))
# - Mixed: valproate vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_valproate_ap.csh <- coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = valproate_ap_mixed %>% mutate(cov_state = ifelse(cov_state %in% c("CO", "DE", "MI"), "grouped", cov_state)))
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_valproate_ap.csh_coef <- cbind(mm_split = "Mixed-to-Dep", analysis = "Valproate", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_valproate_ap.csh))

# - Mania: lamotrigine vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lamotrigine_ap.csh <- coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_race2, data = lamotrigine_ap_mania %>% mutate(cov_state = ifelse(cov_state %in% c("MO", "IN", "AL", "NC", "NY"), "grouped", cov_state)))
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lamotrigine_ap.csh_coef <- cbind(mm_split = "Mania-to-Dep", analysis = "Lamotrigine", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lamotrigine_ap.csh))
# - Mixed: lamotrigine vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lamotrigine_ap.csh <- coxph(Surv(start_csh, end_csh, status0 == 1) ~ med_status + cluster(person_id) + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_race2, data = lamotrigine_ap_mixed %>% mutate(cov_state = ifelse(cov_state %in% c("CO", "DE", "MI"), "grouped", cov_state)))
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lamotrigine_ap.csh_coef <- cbind(mm_split = "Mixed-to-Dep", analysis = "Lamotrigine", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lamotrigine_ap.csh))

# Forest plots
csh_coef <- rbind(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_ms_ap.csh_coef,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_ms_ap.csh_coef,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lithium_ap.csh_coef,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lithium_ap.csh_coef,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_valproate_ap.csh_coef,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_valproate_ap.csh_coef,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lamotrigine_ap.csh_coef,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lamotrigine_ap.csh_coef) %>%
  mutate(HR = exp(estimate), lower95 = exp(estimate - (2*robust.se)), upper95 = exp(estimate + (2*robust.se)))
saveDf(csh_coef, "bipolar_csh_coef", overwrite = T)

csh_coef_download <- csh_coef %>% filter(grepl("med_status", term))
saveDf(csh_coef_download, "bipolar_csh_coef_download", overwrite = T)

cox_zph_hospbipolar_ms_ap_combo_summary_mania <- rbind(data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_ms_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "MS_any", mm_split = "mania"),
  data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lithium_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "lithium", mm_split = "mania"),
  data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_valproate_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "valproate", mm_split = "mania"),
  data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mania_lamotrigine_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "lamotrigine", mm_split = "mania")) %>% 
  mutate(sample = "bipolar_hosp")
saveDf(cox_zph_hospbipolar_ms_ap_combo_summary_mania, "cox_zph_hospbipolar_ms_ap_combo_summary_mania", overwrite = T)

cox_zph_hospbipolar_ms_ap_combo_summary_mixed <- rbind(data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_ms_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "MS_any", mm_split = "mixed"),
  data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lithium_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "lithium", mm_split = "mixed"),
  data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_valproate_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "valproate", mm_split = "mixed"),
  data.frame(cox.zph(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_crr_180d_mixed_lamotrigine_ap.csh)$table) %>% rownames_to_column(var = "term") %>% mutate(analysis = "lamotrigine", mm_split = "mixed")) %>% 
  mutate(sample = "bipolar_hosp")
saveDf(cox_zph_hospbipolar_ms_ap_combo_summary_mixed, "cox_zph_hospbipolar_ms_ap_combo_summary_mixed", overwrite = T)

# PLOT
term_uniq <- unique(csh_coef$term)
term_uniq1 <- term_uniq[which(grepl("med_status", term_uniq))]
term_uniq2 <- substring(term_uniq1, 11, nchar(term_uniq1))
term_uniq2 <- c("MS+AP", "MS-only", "neither")
analysis_uniq <- c("Any Mood Stabiliser", "Lithium", "Valproate", "Lamotrigine")

headtohead_ms_ap_forest_csh.plot <- csh_coef %>% 
  mutate(col = case_when(term == "med_statusneither" ~ "gray40",
                        term == "med_statusMS-only" ~ millerstone[7],
                        term == "med_statusMS+AP" ~ millerstone[1])) %>%
  mutate(sig = ifelse(p.value <= 0.05, col, "white")) %>%
  filter(grepl("med_status", term)) %>%
  mutate(term0 = substring(term, 11, nchar(term))) %>%
  mutate(term = factor(term0, levels=rev(term_uniq2))) %>%
  mutate(analysis = factor(analysis, levels=analysis_uniq)) %>%
  ggplot(., aes(x = term, y = exp(estimate), colour = term)) +
  geom_pointrange(aes(ymin = lower95, ymax = upper95, fill = sig, colour = term), shape = 21) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  theme_bw() + theme(legend.position = "none") +
  xlab("Drugs during hospital for mania/mixed") +
  ylab("Hazard ratio vs Antipsychotic only (adj. for covariates)") +
  scale_y_log10(breaks=c(0.7, 1, 1.5, 2.0)) + 
  coord_flip() + 
  scale_colour_manual(values = c("gray40", millerstone[7], millerstone[1])) +
  scale_fill_identity() +
  facet_grid(mm_split ~ analysis)
ggsave("bipolar_hosp_headtohead_ms_ap_forest_plot_csh.png", headtohead_ms_ap_forest_csh.plot, height = 3.5, width = 8)

# GLM RR APPROACH
ms_ap_mania_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mania") %>%
  mutate(broad_class1 = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, broad_class1) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "broad_class1", values_from = "sum") %>%
  mutate(med_status = case_when(moodstabiliser == 0 & antipsychotic == 0 ~ "neither",
         moodstabiliser >= 1 & antipsychotic == 0 ~ "MS-only",
         moodstabiliser == 0 & antipsychotic >= 1 ~ "1.AP",
         moodstabiliser >= 1 & antipsychotic >= 1 ~ "MS+AP"))
ms_ap_mixed_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mixed") %>%
  mutate(broad_class1 = ifelse(broad_class %in% c("lithium", "antiepileptic"), "moodstabiliser", broad_class)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, broad_class1) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "broad_class1", values_from = "sum") %>%
  mutate(med_status = case_when(moodstabiliser == 0 & antipsychotic == 0 ~ "neither",
         moodstabiliser >= 1 & antipsychotic == 0 ~ "MS-only",
         moodstabiliser == 0 & antipsychotic >= 1 ~ "1.AP",
         moodstabiliser >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lithium_ap_mania_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mania") %>%
  mutate(meds = case_when(meds == "lithium" ~ "lithium", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(valproate == 0 & lamotrigine == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lithium == 0 & antipsychotic == 0 ~ "neither",
         lithium >= 1 & antipsychotic == 0 ~ "MS-only",
         lithium == 0 & antipsychotic >= 1 ~ "1.AP",
         lithium >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lithium_ap_mixed_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mixed") %>%
  mutate(meds = case_when(meds == "lithium" ~ "lithium", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(valproate == 0 & lamotrigine == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lithium == 0 & antipsychotic == 0 ~ "neither",
         lithium >= 1 & antipsychotic == 0 ~ "MS-only",
         lithium == 0 & antipsychotic >= 1 ~ "1.AP",
         lithium >= 1 & antipsychotic >= 1 ~ "MS+AP"))
valproate_ap_mania_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mania") %>%
  mutate(meds = case_when(meds == "valproate" ~ "valproate", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & lamotrigine == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(valproate == 0 & antipsychotic == 0 ~ "neither",
         valproate >= 1 & antipsychotic == 0 ~ "MS-only",
         valproate == 0 & antipsychotic >= 1 ~ "1.AP",
         valproate >= 1 & antipsychotic >= 1 ~ "MS+AP"))
valproate_ap_mixed_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mixed") %>%
  mutate(meds = case_when(meds == "valproate" ~ "valproate", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & lamotrigine == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(valproate == 0 & antipsychotic == 0 ~ "neither",
         valproate >= 1 & antipsychotic == 0 ~ "MS-only",
         valproate == 0 & antipsychotic >= 1 ~ "1.AP",
         valproate >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lamotrigine_ap_mania_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mania") %>%
  mutate(meds = case_when(meds == "lamotrigine" ~ "lamotrigine", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & valproate == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lamotrigine == 0 & antipsychotic == 0 ~ "neither",
         lamotrigine >= 1 & antipsychotic == 0 ~ "MS-only",
         lamotrigine == 0 & antipsychotic >= 1 ~ "1.AP",
         lamotrigine >= 1 & antipsychotic >= 1 ~ "MS+AP"))
lamotrigine_ap_mixed_glmrr <- bipolar_dx7d_hosp_medsind_mmsplit_cov_glmrr_180d.prep %>%
  filter(mm_split == "mixed") %>%
  mutate(meds = case_when(meds == "lamotrigine" ~ "lamotrigine", 
                          broad_class == "antipsychotic" ~ "antipsychotic",
                          .default = meds)) %>%
  group_by(person_id, MDS1m, pmmd01, mm_split, hospital_los, days_since_hosp, cov_age_ep, cov_age_ehr, cov_years_ehr, cov_gender, cov_state, cov_rel, cov_edu2, cov_employ2, cov_race2, meds) %>%
  summarise(sum=sum(drug01)) %>%
  pivot_wider(names_from = "meds", values_from = "sum") %>%
  filter(lithium == 0 & valproate == 0) %>% # clean analysis with no other MS
  mutate(med_status = case_when(lamotrigine == 0 & antipsychotic == 0 ~ "neither",
         lamotrigine >= 1 & antipsychotic == 0 ~ "MS-only",
         lamotrigine == 0 & antipsychotic >= 1 ~ "1.AP",
         lamotrigine >= 1 & antipsychotic >= 1 ~ "MS+AP"))

# - Mania: ms vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_ms_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = ms_ap_mania_glmrr), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_ms_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mania-to-Dep", analysis = "Any Mood Stabiliser", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_ms_ap.glmrr_days_since_hosp))
# - Mixed: ms vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_ms_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = ms_ap_mixed_glmrr), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_ms_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mixed-to-Dep", analysis = "Any Mood Stabiliser", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_ms_ap.glmrr_days_since_hosp))

# - Mania: lithium vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lithium_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = lithium_ap_mania_glmrr %>% left_join(., days_since_hosp %>% distinct)), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lithium_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mania-to-Dep", analysis = "Lithium", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lithium_ap.glmrr_days_since_hosp))
# - Mixed: lithium vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lithium_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = lithium_ap_mixed_glmrr %>% left_join(., days_since_hosp)), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lithium_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mixed-to-Dep", analysis = "Lithium", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lithium_ap.glmrr_days_since_hosp))

# - Mania: valproate vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_valproate_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = valproate_ap_mania_glmrr %>% left_join(., days_since_hosp)), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_valproate_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mania-to-Dep", analysis = "Valproate", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_valproate_ap.glmrr_days_since_hosp))
# - Mixed: valproate vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_valproate_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = valproate_ap_mixed_glmrr %>% left_join(., days_since_hosp)), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_valproate_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mixed-to-Dep", analysis = "Valproate", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_valproate_ap.glmrr_days_since_hosp))

# - Mania: lamotrigine vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lamotrigine_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = lamotrigine_ap_mania_glmrr %>% left_join(., days_since_hosp)), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lamotrigine_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mania-to-Dep", analysis = "Lamotrigine", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lamotrigine_ap.glmrr_days_since_hosp))
# - Mixed: lamotrigine vs ap
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lamotrigine_ap.glmrr_days_since_hosp <- coeftest(glm(MDS1m ~ med_status + hospital_los + days_since_hosp + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = lamotrigine_ap_mixed_glmrr %>% left_join(., days_since_hosp)), vcov = vcovCL, cluster = ~person_id)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lamotrigine_ap.glmrr_coef_days_since_hosp <- cbind(mm_split = "Mixed-to-Dep", analysis = "Lamotrigine", tidy(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lamotrigine_ap.glmrr_days_since_hosp))

# Forest plots
glmrr_coef_days_since_hosp <- rbind(bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_ms_ap.glmrr_coef_days_since_hosp,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_ms_ap.glmrr_coef_days_since_hosp,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lithium_ap.glmrr_coef_days_since_hosp,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lithium_ap.glmrr_coef_days_since_hosp,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_valproate_ap.glmrr_coef_days_since_hosp,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_valproate_ap.glmrr_coef_days_since_hosp,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mania_lamotrigine_ap.glmrr_coef_days_since_hosp,
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_mixed_lamotrigine_ap.glmrr_coef_days_since_hosp) %>%
  mutate(RR = exp(estimate), lower95 = exp(estimate - (2*std.error)), upper95 = exp(estimate + (2*std.error)))
saveDf(glmrr_coef_days_since_hosp, "glmrr_coef_days_since_hosp", overwrite = T)

bipolar_glmrr_coef_days_since_hosp_download <- bipolar_glmrr_coef_days_since_hosp %>% filter(grepl("med_status", term))
saveDf(bipolar_glmrr_coef_days_since_hosp_download, "bipolar_glmrr_coef_days_since_hosp_download", overwrite = T)

# PLOT
term_uniq <- unique(glmrr_coef_days_since_hosp$term)
term_uniq1 <- term_uniq[which(grepl("med_status", term_uniq))]
term_uniq2 <- substring(term_uniq1, 11, nchar(term_uniq1))
term_uniq2 <- c("MS+AP", "MS-only", "neither")
analysis_uniq <- c("Any Mood Stabiliser", "Lithium", "Valproate", "Lamotrigine")

bipolar_hosp_headtohead_ms_ap_forest_plot_glmrr_days_since_hosp.plot <- glmrr_coef_days_since_hosp %>% 
  mutate(col = case_when(term == "med_statusneither" ~ "gray40",
                        term == "med_statusMS-only" ~ millerstone[7],
                        term == "med_statusMS+AP" ~ millerstone[1])) %>%
  mutate(sig = ifelse(p.value <= 0.05, col, "white")) %>%
  filter(grepl("med_status", term)) %>%
  mutate(term0 = substring(term, 11, nchar(term))) %>%
  mutate(term = factor(term0, levels=rev(term_uniq2))) %>%
  mutate(analysis = factor(analysis, levels=analysis_uniq)) %>%
  ggplot(., aes(x = term, y = exp(estimate), colour = term)) +
  geom_pointrange(aes(ymin = lower95, ymax = upper95, fill = sig, colour = term), shape = 21) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  theme_bw() + theme(legend.position = "none") +
  xlab("Drugs during hospital for mania/mixed") +
  ylab("Risk ratio vs Antipsychotic only (adj. for covariates)") +
  scale_y_log10() +
#  scale_y_log10(breaks=c(0.6, 0.8, 1, 1.2, 1.4)) + 
  coord_flip() +
  scale_colour_manual(values = c("gray40", millerstone[7], millerstone[1])) +
  scale_fill_identity() +
  facet_grid(mm_split ~ analysis)
ggsave("bipolar_hosp_headtohead_ms_ap_forest_plot_glmrr_days_since_hosp.svg", bipolar_hosp_headtohead_ms_ap_forest_plot_glmrr_days_since_hosp.plot, height = 3.5, width = 8)

#---------------------------------------------------
# Compare to TTE results
#---------------------------------------------------

tte_output <- read_csv("PostMania_Discharge_Landmark_additional_estimands_30d.csv")

tte_output.plot <- tte_output %>%
  filter(treatment != "AP_only") %>% filter(estimand == "RR") %>%
  mutate(sig = ifelse(p_value < 0.05, millerstone[1], "white")) %>%
  mutate(analysis = case_when(outcome == "depression_composite" ~ "Depression",
                              outcome == "urti_negative_control" ~ "URTI (neg control)",
                              outcome == "allergic_rhinitis_negative_control" ~ "Allergic rhinitis (neg control)")) %>%
  mutate(analysis = factor(analysis, levels=rev(c("Depression", "URTI (neg control)", "Allergic rhinitis (neg control)")))) %>%  
  mutate(treatment = gsub("_plus_", "+", treatment)) %>%                
  ggplot(., aes(x = analysis, y = estimate, colour = sig)) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, fill = sig, colour = sig), shape = 21) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  theme_bw() + theme(legend.position = "none") +
  xlab("Outcome (within 30d of hospital visit end)") +
  ylab("Risk ratio vs Antipsychotic only (target trial emulation)") +
  scale_y_log10() +
  coord_flip(ylim = c(0.1, 7)) +
  scale_colour_manual(values = c(millerstone[1], millerstone[1])) +
  scale_fill_identity() +
  facet_wrap(~treatment)
ggsave("tte_output.png", tte_output.plot, height = 3.5, width = 8, dpi = 600)

