
#============================================================
# 
# PREDICTORS of depression from the manic/mixed episode
# ##### RESTORED from 22/12/25 version ... need to restore back on 30/12/25
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

#------------------------------------------------------------
# Episode-specific
# - hospitalisation
# - duration of mania/mixed episode
# - CGI-S maximum
# - mania/mixed psychotic
#------------------------------------------------------------

#============================================================
# Munge
#============================================================

#------------------------------------------------------------
# Severity proxies
#------------------------------------------------------------

icdseverity_split <- getDf("icdseverity_split")
pmmdgroup_cgis <- getDf("pmmdgroup_cgis")
nopmmdgroup_cgis <- getDf("nopmmdgroup_cgis")
cgis_mania <- rbind(pmmdgroup_cgis, nopmmdgroup_cgis)
cgis_mania$person_id <- as.character(cgis_mania$person_id)

# prior script in timetodepression_outcome.r: source('remove_spurious_episodes.r'), so have QC-ed already
bipolar_epall_any_overlap <- getDf("bipolar_epall_any_overlap") 
bipolar_epall_hospital_overlap <- getDf("bipolar_epall_hospital_overlap")

# - munge 
  # mm_split already there
bipolar_mm_clinical.tmp <- bipolar_epall_any_overlap %>%
  mutate(reltoPMMD1 = ifelse(is.na(reltoPMMD1), "noPMMD", reltoPMMD1)) %>%
  filter(state == "mm") %>%
  # add hospital episodes
  left_join(., bipolar_epall_hospital_overlap %>% 
                dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end) %>%
                mutate(hospital = 1)) %>%
  mutate(hospital01 = ifelse(is.na(hospital), 0, hospital)) %>%
  mutate(MDS = ifelse((is.na(PMMD) | PMMD == "nonPMMD"), 0, 1)) %>%
  # convert los to weeks
  mutate(los_weeks = los/7) %>%
  # add mania/mixed duration
  mutate(mm_duration_weeks = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days")/7)) %>%
  # add psychotic features
  left_join(., icdseverity_split) %>% mutate(icdseverity = ifelse(icdseverity == "2mild", "0mild", icdseverity)) %>%
  left_join(., demog_cov) %>%
  mutate(cov_age_ep = as.numeric(difftime(diagnosis_date, birth_year))/365) %>%
  # make sure the "no matching concept" variable is the reference factor level for regression
  mutate(cov_employ = ifelse(cov_employ == "No matching concept", "0No matching concept", cov_employ)) %>%
  mutate(cov_edu = ifelse(cov_edu == "No matching concept", "0No matching concept", cov_edu)) %>%
  mutate(cov_rel = ifelse(cov_rel == "No matching concept", "0No matching concept", cov_rel)) %>%
  mutate(cov_state = ifelse(cov_state == "No matching concept", "0No matching concept", cov_state)) %>%
  # because there are only n=2 for self-employed, absorb within "0No matching concept"
  mutate(cov_employ = ifelse(cov_employ == "Self-employed", "0No matching concept", cov_employ)) %>%
  mutate(los_per_visit = los/record_count_visit) %>%
  # add short MDS field (<1m after)
  left_join(pmmd_ep %>% filter(state == "mm") %>% dplyr::select(person_id, mm_episode, dep_episode, state, diagnosis_date, dep_date) %>% distinct()) %>%
  mutate(timetodep = as.numeric(difftime(dep_date, diagnosis_date_end, units = "days"))) %>%
  filter(cov_gender != "No matching concept") %>%
  mutate(cov_edu2 = case_when(cov_edu %in% c("8th grade or less", "Some high school, but did not graduate", "High school diploma or GED") ~ "School",
                              cov_edu %in% c("Some college", "Associate degree (e.g., AA, AS)", "Bachelor's degree", "Graduate degree") ~ "College",
                              .default = cov_edu)) %>%
  mutate(cov_employ2 = case_when(cov_employ %in% c("Full-time employment", "Part-time employment", "Self-employed") ~ "Employed",
                                cov_employ %in% c("Not in labor force", "Unemployed", "Retired", "Permanently unable to perform work activities due to medical condition", "On work experience") ~ "NotEmployed",
                                .default = cov_employ)) %>%
  mutate(cov_race2 = ifelse(cov_race %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander"), "AmIndian_Alaskan_PacIslander", cov_race))

# - add on prevMDS
tmp_prev_mds.rr <- bipolar_mm_clinical.tmp %>%
  filter(timetodep <= 30) %>%
  group_by(person_id) %>%
  mutate(first_mm_episode_mds = min(mm_episode)) %>%
  dplyr::select(person_id, first_mm_episode_mds) %>%
  distinct()
# - apply minimum cov_state group size of 1000 (0.5% of 200K)
cov_state_n.csh180d <- bipolar_mm_clinical.tmp %>% group_by(person_id, cov_state) %>% summarise(n=n()) %>% group_by(cov_state) %>% summarise(n=n()) %>%
  mutate(cov_state2 = ifelse(n <= 1000, "grouped", cov_state)) %>% ungroup
# - apply minimum hospital cov_state size (as some sites have barely any hospitalisation)
cov_state_hospital01_n <- bipolar_mm_clinical.tmp %>% 
  group_by(cov_state, hospital01, mm_split) %>% summarise(n=n()) %>%
  mutate(hospital01 = ifelse(hospital01 == 0, "nohosp", "hosp")) %>%
  pivot_wider(names_from = hospital01, values_from = n) %>%
  mutate(hosp = ifelse(is.na(hosp), 0, hosp), nohosp = ifelse(is.na(nohosp), 0, nohosp)) %>%
  mutate(prop_hosp = hosp/(hosp+nohosp))
cov_state_hospital01_n_keep <- cov_state_hospital01_n %>% 
  mutate(hospital_include = ifelse(prop_hosp >= 0.01, 1, 0)) %>%
  dplyr::select(cov_state, mm_split, hospital_include)

bipolar_mm_clinical <- bipolar_mm_clinical.tmp %>%
  left_join(., tmp_prev_mds.rr, relationship = "many-to-many") %>%
  left_join(., cov_state_n.csh180d) %>%
  mutate(prev_MDS = case_when(reltoPMMD1 == "after" & first_mm_episode_mds < mm_episode ~ "prevMDS", 
                              mm_episode == 1 ~ "noprevMM",
                              .default = "0prevMM")) %>%
  mutate(MDS1m = ifelse(timetodep <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) %>%
  left_join(., cov_state_hospital01_n_keep)

saveDf(bipolar_mm_clinical, "bipolar_mm_clinical", overwrite = T)

# add on peak CGI-S mania per mm episodes (with a window of 0-7 days) - done separately as only available for a subset
bipolar_mm_cgis_maxd7 <- bipolar_epall_any_overlap %>%
  filter(state == "mm") %>%
  left_join(., cgis_mania %>% distinct, relationship = "many-to-many") %>%
  mutate(time_since_diagnosis = as.numeric(difftime(measurement_date, diagnosis_date, units = "days"))) %>%
  filter(time_since_diagnosis >= 0 & time_since_diagnosis <= 7) %>% # filter for CGI-S within 7 days of diagnosis
  group_by(person_id, diagnosis_date) %>% # get per diagnosis date, get the maximum CGI-S
  arrange(desc(score), time_since_diagnosis) %>%
  slice_head(n=1) %>%
  rename(severitycgis_max_d07 = score)

# - apply minimum cov_state group size of 100 (0.5% of 2000)
cov_state_n_cgis.csh180d <- bipolar_mm_clinical.tmp %>% 
  left_join(., bipolar_mm_cgis_maxd7) %>%
  filter(!is.na(severitycgis_max_d07)) %>% 
  group_by(person_id, cov_state) %>% summarise(n=n()) %>% group_by(cov_state) %>% summarise(n=n()) %>%
  mutate(cov_state2 = ifelse(n <= 100, "grouped", cov_state)) %>% ungroup

bipolar_mm_clinical_cgis <- bipolar_mm_clinical.tmp %>%
  left_join(., tmp_prev_mds.rr, relationship = "many-to-many") %>%
  mutate(prev_MDS = case_when(reltoPMMD1 == "after" & first_mm_episode_mds < mm_episode ~ "prevMDS", 
                              mm_episode == 1 ~ "noprevMM",
                              .default = "0prevMM")) %>%
  mutate(MDS1m = ifelse(timetodep <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) %>%
  left_join(., cov_state_n_cgis.csh180d) %>%
  left_join(., bipolar_mm_cgis_maxd7) %>%
  filter(!is.na(severitycgis_max_d07)) %>%
  left_join(., cov_state_hospital01_n_keep)

saveDf(bipolar_mm_clinical_cgis, "bipolar_mm_clinical_cgis", overwrite = T)

# Compare CGI-S in 7d to icdseverity_rank
bipolar_mm_clinical_cgis %>% 
  mutate(icdseverity = ifelse(icdseverity == "0mild", "2mild", icdseverity)) %>%
  filter(hospital_include == 1) %>%
  mutate(hospital01 = ifelse(hospital01 == 1, "hospital", "no hospital")) %>%
  ggplot(., aes(x = icdseverity, y = severitycgis_max_d07, fill = icdseverity)) +
  geom_boxplot() +
  scale_fill_manual(values = millerstone) +
  theme_bw() + theme(legend.position = "none") +
  ylab("CGI-S maximum score within first 7 days of diagnosis") + xlab("ICD severity category") +
  facet_grid(hospital01 ~ mm_split)
ggsave("cgis_icdseverity_boxplot.png", width = 10, height = 6)

# CHECKS:
# quantify N
bipolar_mm_clinical %>% dim
length(unique(bipolar_mm_clinical$person_id))
bipolar_mm_clinical_cgis %>% filter(!is.na(severitycgis_max_d07)) %>% dim
bipolar_mm_clinical_cgis %>% filter(!is.na(severitycgis_max_d07)) %>% group_by(person_id) %>% summarise(n=n()) %>% dim
bipolar_mm_clinical %>% filter(mm_episode >= 2) %>% dim
bipolar_mm_clinical %>% filter(mm_episode >= 2) %>% group_by(person_id) %>% summarise(n=n()) %>% dim
bipolar_mm_clinical %>% filter(hospital_include == 1) %>% dim
bipolar_mm_clinical %>% filter(hospital_include == 1) %>% pull(person_id) %>% unique %>% length
bipolar_mm_clinical %>% filter(hospital_include == 1) %>% group_by(hospital01) %>% summarise(n=n())
bipolar_mm_clinical %>% filter(hospital_include == 1) %>% group_by(person_id, hospital01) %>% summarise(n=n()) %>% group_by(hospital01) %>% summarise(n=n())
# absolute rates in groups
bipolar_mm_clinical %>% filter(mm_episode >= 2) %>% group_by(mm_split, prev_MDS, MDS1m) %>% summarise(n=n())
bipolar_mm_clinical %>% group_by(mm_split, hospital, MDS1m) %>% summarise(n=n())

# Check how max CGI-S time maps to mm diagnosis date end (ie. does it happen before or after the end)
bipolar_mm_cgis_maxd7_mm_time <- bipolar_mm_cgis_maxd7 %>% 
  filter(!is.na(severitycgis_max_d07)) %>% 
  mutate(time_mm_start_cgis = as.numeric(difftime(measurement_date, diagnosis_date, units = "days")),
         time_mm_end_cgis = as.numeric(difftime(measurement_date, diagnosis_date_end, units = "days")))

bipolar_mm_cgis_maxd7_mm_time %>% pull(time_mm_start_cgis) %>% summary
bipolar_mm_cgis_maxd7_mm_time %>% pull(time_mm_end_cgis) %>% summary
# CONFIRM that almost all are on the same date as the diagnosis start date (1562/1793)

#======================================================
# Further confounding by indication plots and data exploration
#======================================================

source("predictdep_mm_medications_all_novisitconstraint.r")
# - gives the df bipolar_dx7d_dx_medsbroadclass_mmsplit_joint_cov_crr 
# - and plot bipolar_dx7d_dx_medsbroadclass_mmsplit_joint_cov_crr.plot

#======================================================
#
# Bipolar group (primary analysis)
#
#======================================================
#----------------------
# Explore relationships between variables relating to the mm episode

bipolar_mm_clinical %>% group_by(person_id, MDS) %>% summarise(n=n()) %>% group_by(MDS) %>% summarise(n=n())

bipolar_mm_clinical %>% group_by(hospital01, icdseverity) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = "hospital01", values_from = "n") %>%
  mutate(sum = `0` + `1`, prop1 = `1`/sum)

bipolar_mm_clinical %>% group_by(hospital01, mm_split) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = "hospital01", values_from = "n") %>%
  mutate(sum = `0` + `1`, prop1 = `1`/sum)

bipolar_mm_clinical_cgis %>% group_by(hospital01, severitycgis_max_d07) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = "hospital01", values_from = "n") %>%
  mutate(sum = `0` + `1`, prop1 = `1`/sum) # - CGI-S subset

# - check association with MDS status

bipolar_mm_clinical %>%
  mutate(MDS = ifelse(MDS == 1, "MDS", "noMDS")) %>%
  ggplot(., aes(x = MDS, y = icdseverity_rank, fill = as.factor(hospital01))) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = millerstone)
ggsave("tmp.png")

bipolar_mm_clinical_cgis %>%
  mutate(MDS = ifelse(MDS == 1, "MDS", "noMDS")) %>%
  ggplot(., aes(x = MDS, y = severitycgis_max_d07, fill = as.factor(hospital01))) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = millerstone)
ggsave("tmp.png")

bipolar_mm_clinical %>% group_by(MDS, mm_split, hospital01) %>% summarise(n=n()) %>% 
  mutate(MDS = ifelse(MDS == 1, "MDS", "noMDS")) %>%
  pivot_wider(names_from = "MDS", values_from = "n") %>%
  mutate(noMDS = ifelse(is.na(noMDS), 0, noMDS)) %>%
  mutate(sum = MDS + noMDS, prop1 = MDS/sum)
bipolar_mm_clinical_icdseverity_hosp_barplot <- bipolar_mm_clinical %>% group_by(MDS, icdseverity, mm_split, hospital01) %>% summarise(n=n()) %>% 
  mutate(MDS = ifelse(MDS == 1, "MDS", "noMDS")) %>%
  pivot_wider(names_from = "MDS", values_from = "n") %>%
  mutate(noMDS = ifelse(is.na(noMDS), 0, noMDS)) %>%
  mutate(sum = MDS + noMDS, prop1 = MDS/sum)
bipolar_mm_clinical %>% group_by(MDS, icdseverity) %>% summarise(n=n()) %>% 
  mutate(MDS = ifelse(MDS == 1, "MDS", "noMDS")) %>%
  pivot_wider(names_from = "MDS", values_from = "n") %>%
  mutate(sum = MDS + noMDS, prop1 = MDS/sum)

bipolar_mm_clinical_cgis_barplot <- bipolar_mm_clinical_cgis %>% group_by(MDS1m, severitycgis_max_d07, mm_split) %>% summarise(n=n()) %>% 
  mutate(MDS1m = ifelse(MDS1m == 1, "MDS", "noMDS")) %>%
  pivot_wider(names_from = "MDS1m", values_from = "n") %>%
  mutate(noMDS = ifelse(is.na(noMDS), 0, noMDS)) %>%
  mutate(sum = MDS + noMDS, prop1 = MDS/sum) %>%
  arrange(severitycgis_max_d07)

# - how many have MDS with prevMDS
prevmds_tmp <- bipolar_mm_clinical %>% filter(mm_episode >= 2) %>% filter(prev_MDS == "prevMDS") %>% group_by(mm_split, prev_MDS, MDS1m) %>% summarise(n=n())
prevmds_tmp %>% group_by(mm_split) %>% mutate(sum = sum(n)) %>% ungroup %>% mutate(prop = n/sum)

# - how many have MDS when CGI-S < 5
cgis_less5 <- bipolar_mm_clinical_cgis_barplot %>% filter(severitycgis_max_d07 %in% 1:4)
sum(cgis_less5$MDS)/(sum(cgis_less5$MDS) + sum(cgis_less5$noMDS))
cgis_more5 <- bipolar_mm_clinical_cgis_barplot %>% filter(severitycgis_max_d07 %in% 5:7)
sum(cgis_more5$MDS)/(sum(cgis_more5$MDS) + sum(cgis_more5$noMDS))

# - how many have MDS with hospitalisation
hospital_tmp <- bipolar_mm_clinical %>% group_by(mm_split, hospital, MDS1m) %>% summarise(n=n())

# empirical CDFs (ecdf)
# - prevMDS
prevMDS_ttd_ecdf_plot <- bipolar_mm_clinical %>% 
  filter(mm_episode >= 2) %>% 
  mutate(prev_MDS = ifelse(prev_MDS == "0prevMM", "noprevMDS", "prevMDS")) %>%
  ggplot(., aes(x = timetodep, colour = prev_MDS)) +
  stat_ecdf() +
  theme_bw() +
  scale_colour_manual(values = millerstone) +
  xlab("Depression transition time (days)") +
  ylab("Cumulative proportion") +
  coord_cartesian(xlim = c(0,180)) +
  facet_wrap(~ mm_split, nrow = 1)

# - CGI-S
cgis_ttd_ecdf_plot <- bipolar_mm_clinical_cgis %>% 
  filter(!is.na(severitycgis_max_d07)) %>%
  mutate(CGIS_grouped = case_when(severitycgis_max_d07 %in% 0:3 ~ "0-3",
                                  severitycgis_max_d07 %in% 6:7 ~ "6-7",
                                  .default = as.character(severitycgis_max_d07))) %>% # ensure at least 50 or so each in category
  ggplot(., aes(x = timetodep, colour = CGIS_grouped)) +
  stat_ecdf() +
  theme_bw() +
  scale_colour_manual(values = millerstone) +
  xlab("Depression transition time (days)") +
  ylab("Cumulative proportion") +
  coord_cartesian(xlim = c(0,180)) +
  facet_wrap(~ mm_split, nrow = 1)

# - hospital
hospital_ttd_ecdf_plot <- bipolar_mm_clinical %>% 
  filter(hospital_include == 1) %>%
  mutate(hospital = ifelse(is.na(hospital), "nohosp", "hosp")) %>%
  ggplot(., aes(x = timetodep, colour = hospital)) +
  stat_ecdf() +
  theme_bw() +
  scale_colour_manual(values = millerstone) +
  xlab("Depression transition time (days)") +
  ylab("Cumulative proportion") +
  coord_cartesian(xlim = c(0,180)) +
  facet_wrap(~ mm_split, nrow = 1)

library(patchwork)
prevMDS_ttd_ecdf_plot / cgis_ttd_ecdf_plot / hospital_ttd_ecdf_plot +
plot_annotation(tag_levels = 'a')
ggsave("clinical_ttd_ecdf_plot.png", height = 10, width = 9)

# Plot MDS against increasing MM severity
bipolar_mm_clinical_cgis_barplot_gg <- ggplot(bipolar_mm_clinical_cgis_barplot, aes(x = severitycgis_max_d07, y = prop1)) +
  geom_col(fill = millerstone[2]) +
  ylab("Proportion of MM episodes followed by Dep") +
  xlab("CGI-S score") +
  theme_bw() +
  facet_wrap(~mm_split)
ggsave("bipolar_mm_clinical_cgis_barplot_gg.png", bipolar_mm_clinical_cgis_barplot_gg, height = 5, width = 7)

# Plot MDS against increasing ICD severity, stratified by hospitalisation
bipolar_mm_clinical_icdseverity_hosp_barplot_gg <- ggplot(bipolar_mm_clinical_icdseverity_hosp_barplot %>% mutate(hospital01 = as.factor(hospital01)), aes(x = icdseverity, y = prop1, fill = hospital01)) +
  geom_col(position = "dodge") +
  ylab("Proportion of MM episodes followed by Dep") +
  xlab("ICD severity") +
  scale_fill_manual(values = millerstone[c(3,9)]) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~mm_split)
ggsave("bipolar_mm_clinical_icdseverity_hosp_barplot.png", bipolar_mm_clinical_icdseverity_hosp_barplot_gg, height = 5, width = 7)

#============================================================
# Modified poisson regression models for RR
#============================================================

#-------------------
# RR: Bipolar group
#-------------------

# Modified poisson regression models for RR
# - NOCOV
bipolar_mm_icdseverity_glmrr <- bipolar_mm_clinical %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ icdseverity, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

bipolar_mm_hospital_glmrr <- bipolar_mm_clinical %>%
  filter(hospital_include == 1) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ hospital01, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

bipolar_mm_cgis_glmrr <- bipolar_mm_clinical_cgis %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ severitycgis_max_d07, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

bipolar_mm_prevMDS_glmrr <- bipolar_mm_clinical %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ prev_MDS, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

bipolar_mm_clinical_glmrr_nocov <- rbind(bipolar_mm_prevMDS_glmrr,
    bipolar_mm_cgis_glmrr,
    bipolar_mm_hospital_glmrr,
    bipolar_mm_icdseverity_glmrr
    ) %>%
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(p_adj = p.adjust(p.value*3)) %>% # adjust for prevMDS, severity, hospital
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "no covariates")

# - COV
bipolar_mm_icdseverity_cov_glmrr <- bipolar_mm_clinical %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ icdseverity + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_hospital_cov_glmrr <- bipolar_mm_clinical %>%
  filter(hospital_include == 1) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_cgis_cov_glmrr <- bipolar_mm_clinical_cgis %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_prevMDS_cov_glmrr <- bipolar_mm_clinical %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_glmrr_cov <- rbind(bipolar_mm_clinical_prevMDS_cov_glmrr %>% mutate(primary_var = "prevMDS"),
  bipolar_mm_cgis_cov_glmrr %>% mutate(primary_var = "CGI-S severity"),
  bipolar_mm_hospital_cov_glmrr %>% mutate(primary_var = "hospital01"),
  bipolar_mm_icdseverity_cov_glmrr %>% mutate(primary_var = "ICD severity")) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(p_adj = p.value*3) %>% # 3 tests for each of prevMDS, severity, hospital
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "covariates")

#------------------------------------------------------------
# Summary plot
#------------------------------------------------------------

#-------------------
# Bipolar
#-------------------

bipolar_mm_clinical_glmrr_plot <- full_join(bipolar_mm_clinical_glmrr_nocov,
  bipolar_mm_clinical_glmrr_cov) # do away with bipolar_mm_clinical_glmrr_joint)
saveDf(bipolar_mm_clinical_glmrr_plot, "bipolar_mm_clinical_glmrr_plot", overwrite = T)

bipolar_mm_clinical_glmrr_plot.tmp <- bipolar_mm_clinical_glmrr_plot %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(term = case_when(term == "severitycgis_max_d07" ~ "severity: CGI-S score +1",
    term == "prev_MDSprevMDS" ~ "prev short MDTT (<1m)",
    #term == "mm_splitmixed" ~ "mixed ep (vs mania)",
    term == "hospital01" ~ "hospital for mania/mixed",
    .default = term))

bipolar_mm_clinical_glmrr_plot.tmp %>%
  filter(facet != "joint model") %>%
  mutate(facet = factor(facet, levels = c("no covariates", "covariates"))) %>%
#  mutate(term2 = factor(term, levels = rev(c("previous MDS", "severity: CGI-S score +1", "mixed ep (vs mania", "hospital for mania/mixed", "0icdseveritymild", "icdseverity1hypomania", "icdseverity3moderate", "icdseverity4severe", "icdseverity5psychosis", "icdseverityunspecified"))) %>%
  mutate(term2 = factor(term, levels = rev(unique(bipolar_mm_clinical_glmrr_plot.tmp$term)))) %>%
  ggplot(., aes(x = term2, y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 2.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate for mania/mixed-depression switch") +
  coord_flip(ylim = c(0, 5)) +
  facet_grid(facet~mm_split)
ggsave("bipolar_mm_clinical_glmrr_plot.png", width = 7, height = 5)

# For final figure
bipolar_mm_clinical_glmrr_plot_covfinal <- bipolar_mm_clinical_glmrr_plot.tmp %>%
  filter(!grepl("icdseverity", term)) %>%
  filter(facet == "covariates") %>%
  mutate(title = "Risk Ratio of Switch in 1m") %>%
#bipolar_mm_clinical_all_glmrr_plot_gg <- bipolar_mm_clinical_glmrr_cov.plot %>%
#  full_join(., bipolar_dx7d_dx_medsbroadclass_mmsplit_joint_cov_glmrr.plot) %>%
  mutate(term = ifelse(term == "previous Dep switch (1m)", "prev short MDTT (<1m)", term)) %>%
  mutate(term2 = factor(term, levels = rev(c("prev short MDTT (<1m)", "severity: CGI-S score +1", "hospital for mania/mixed", "antidepressant", "antiepileptic", "antipsychotic", "benzodiazepine", "lithium", "antihistamine")))) %>%
  ggplot(., aes(x = term2, y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR), shape = 17) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_label_repel(size = 2.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  #scale_y_continuous(limits = c(0.4, 2.5)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate, full bipolar disorder sample (adj. covariates)") +
  coord_flip() +
  facet_grid(facet ~ mm_split, scales = "free",  space = "free_y")
ggsave("bipolar_mm_clinical_glmrr_plot_covfinal.svg", bipolar_mm_clinical_glmrr_plot_covfinal, width = 7, height = 3)

#============================================================
# Bipolar group: Time to depression analysis (CSH)
#============================================================

#-------------------
# Prep data
#-------------------

mania_dep_dates_toepisodes_qc_defbroad <- getDf("mania_dep_dates_toepisodes_qc_defbroad") 

# - survival analysis of time to depression
nonpmmdep_censor_epall_visit <- mania_dep_dates_toepisodes_qc_defbroad %>%
  filter(is.na(dep_date)) %>%
  dplyr::select(c(person_id, diagnosis_date)) %>%
  rename(index_date = diagnosis_date)

nonpmmdep_nextep_epall_visit <- nonpmmdep_censor_epall_visit %>%
  full_join(., mania_dep_dates_toepisodes_qc_defbroad, relationship = "many-to-many") %>%
  filter(is.na(dep_date)) %>%
  mutate(timediff = as.numeric(difftime(diagnosis_date, index_date, units = "days"))) %>%
  filter(timediff > 0) %>%
  arrange(timediff) %>%
  group_by(person_id, index_date) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  rename(nextep = diagnosis_date) %>%
  dplyr::select(person_id, index_date, state, nextep, timediff)

observation_period <- getDf("observation_period")
observation_period.tmp <- observation_period %>%
  mutate(person_id = as.character(person_id)) %>%
  filter(person_id %in% nonpmmdep_censor_epall_visit$person_id)
nonpmmdep_obsperiod_epall_visit <- nonpmmdep_censor_epall_visit %>%
  left_join(., observation_period.tmp, relationship = "many-to-many") %>%
  rename(obsend = observation_period_end_date) %>%
  dplyr::select(person_id, index_date, obsend)

nonpmmdep_censor_date_epall_visit <- full_join(nonpmmdep_obsperiod_epall_visit, nonpmmdep_nextep_epall_visit) %>%
  mutate(end_date0 = pmin(obsend, nextep, na.rm = T)) %>%
  mutate(end_time0 = as.numeric(difftime(end_date0, index_date, units = "days"))) %>%
  mutate(status0 = case_when(is.na(state) ~ 0,
                            state == "dep" ~ 1, 
                            state == "mm" ~ 2)) %>%
  mutate(diagnosis_date = index_date) %>%
  dplyr::select(person_id, diagnosis_date, end_date0, end_time0, status0)
saveDf(nonpmmdep_censor_date_epall_visit, "nonpmmdep_censor_date_epall_visit", overwrite = T)

# - Clinical variable competing risk
mania_dep_dates_toepisodes_qc_defbroad.tmp <- mania_dep_dates_toepisodes_qc_defbroad %>%
  filter(state == "mm") %>%
  dplyr::select(person_id, diagnosis_date, dep_date, mm_episode, state)

bipolar_mm_clinical_csh.prep <- bipolar_mm_clinical %>%
  left_join(., nonpmmdep_censor_date_epall_visit) %>%
  left_join(., mania_dep_dates_toepisodes_qc_defbroad.tmp) %>%
  left_join(., bipolar_mm_cgis_maxd7) %>%  
  mutate(status = ifelse(is.na(PMMD), status0, 1)) %>%
  mutate(end_date = coalesce(dep_date, end_date0)) %>%
  mutate(end_time = ifelse(!is.na(dep_date), as.numeric(difftime(dep_date, diagnosis_date_end, units = "days")), end_time0)) %>%
  dplyr::select(-c(end_date0, end_time0, status0, end_date)) %>%
  filter(mm_split != "both") %>%
  mutate(cov_age_ep = as.numeric(difftime(diagnosis_date, as.Date(birth_year), units = "days"))/365) %>% dplyr::select(-birth_year) %>%
  mutate(time_mm_dep = as.numeric(difftime(dep_date, diagnosis_date, units = "days"))) %>%
  #dplyr::select(-c("diagnosis_date", "visit_start_date", "visit_end_date", "dep_date")) %>%
  mutate(censor = ifelse(end_time >= 180, 1, 0)) %>%
  mutate(end_time = ifelse(censor == 1, 180, end_time),
         status = ifelse(censor == 1, 0, status)) %>%
  mutate(start_csh = 0) %>%
  mutate(end_csh = ifelse(end_time == start_csh, end_time + 1, end_time)) %>%
  mutate(status = as.factor(status))

saveDf(bipolar_mm_clinical_csh.prep, "bipolar_mm_clinical_csh.prep", overwrite = T)

bipolar_mm_clinical_csh.prep <- getDf("bipolar_mm_clinical_csh.prep")

#-------------------
# CSH: all bipolar
#-------------------
# - nocov
bipolar_mm_clinical_hospital_nocov.csh <- bipolar_mm_clinical_csh.prep %>%
  filter(hospital_include == 1) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ hospital01 + cluster(person_id), data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

bipolar_mm_clinical_icdseverity_nocov.csh <- bipolar_mm_clinical_csh.prep %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ icdseverity + cluster(person_id), data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

bipolar_mm_clinical_cgis_nocov.csh <- bipolar_mm_clinical_csh.prep %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ severitycgis_max_d07 + cluster(person_id), data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

bipolar_mm_clinical_prevMDS_nocov.csh <- bipolar_mm_clinical_csh.prep %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ prev_MDS + cluster(person_id), data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "univariable")

#bipolar_mm_clinical_hospital_nocov.crr <- tidy(crr(Surv(end_time, status) ~ hospital01, data = bipolar_mm_clinical_crr.prep))
#bipolar_mm_clinical_mmsplit_nocov.crr <- tidy(crr(Surv(end_time, status) ~ mm_split, data = bipolar_mm_clinical_crr.prep))
#bipolar_mm_clinical_icdseverity_nocov.crr <- tidy(crr(Surv(end_time, status) ~ icdseverity, data = bipolar_mm_clinical_crr.prep))
#bipolar_mm_clinical_cgis_nocov.crr <- tidy(crr(Surv(end_time, status) ~ severitycgis_max_d07, data = bipolar_mm_clinical_crr.prep))
#bipolar_mm_clinical_prevMDS_nocov.crr <- tidy(crr(Surv(end_time, status) ~ prev_MDS, data = bipolar_mm_clinical_crr.prep %>% filter(mm_episode >= 2)))
#bipolar_mm_clinical_prevMDS_nocov.crr <- tidy(crr(Surv(end_time, status) ~ timetodep1, data = bipolar_mm_clinical_crr.prep %>% filter(mm_episode >= 2)))

bipolar_mm_clinical_nocov_csh <- rbind(bipolar_mm_clinical_hospital_nocov.csh, 
  bipolar_mm_clinical_icdseverity_nocov.csh,
  bipolar_mm_clinical_cgis_nocov.csh,
  bipolar_mm_clinical_prevMDS_nocov.csh) %>%
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)),
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(p_adj = p.value*3) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80"))
saveDf(bipolar_mm_clinical_nocov_csh, "bipolar_mm_clinical_nocov_csh", overwrite = T)

# - cov
bipolar_mm_clinical_hospital_cov.csh_coef <- bipolar_mm_clinical_csh.prep %>%
  filter(hospital_include == 1) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ hospital01 + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  mutate(analysis = "cov_hospital01")

bipolar_mm_clinical_icdseverity_cov.csh_coef <- bipolar_mm_clinical_csh.prep %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ icdseverity + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  mutate(analysis = "cov_ICD")

bipolar_mm_clinical_cgis_cov.csh_coef <- bipolar_mm_clinical_csh.prep %>% # LEAVE OUT cov_employ and cov_state due to singularities
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ severitycgis_max_d07 + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  mutate(analysis = "cov_CGI-S")

state_tmp <- bipolar_mm_clinical_csh.prep %>% filter(mm_episode >= 2) %>% group_by(prev_MDS, cov_state2) %>% summarise(n=n()) %>%
  mutate(cov_state2 = ifelse(n <= 10, "grouped", cov_state2)) # - do this due to small N
bipolar_mm_clinical_prevMDS_cov.csh_coef <- bipolar_mm_clinical_csh.prep %>%
  filter(mm_episode >= 2) %>%
  left_join(., state_tmp) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coxph(Surv(start_csh, end_csh, status == 1) ~ prev_MDS + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  mutate(analysis = "cov_prevMDS")

bipolar_mm_clinical_cov_csh <- rbind(bipolar_mm_clinical_hospital_cov.csh_coef, 
  bipolar_mm_clinical_icdseverity_cov.csh_coef,
  bipolar_mm_clinical_cgis_cov.csh_coef,
  bipolar_mm_clinical_prevMDS_cov.csh_coef) %>%
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_HR = exp(estimate),
        lowerHR = exp(estimate - (std.error * 1.96)),
        higherHR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_HR), 
    " [", sprintf('%.2f', lowerHR), "-", sprintf('%.2f', higherHR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  group_by(mm_split) %>%
  mutate(p_adj = ifelse(!grepl("cov|Intercept", term), p.value*3, NA)) %>%
  ungroup() %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80"))
saveDf(bipolar_mm_clinical_cov_csh, "bipolar_mm_clinical_cov_csh", overwrite = T)

#------------------------------------------------------------
# Bipolar: Summary plot of CSH results
#------------------------------------------------------------

bipolar_mm_clinical_csh_plot <- bipolar_mm_clinical_nocov_csh %>% mutate(analysis = "univariable", facet = "no covariates") %>%
  full_join(., bipolar_mm_clinical_cov_csh %>% mutate(facet = "covariates")) %>%
  filter(!grepl("cov|Intercept", term))
saveDf(bipolar_mm_clinical_csh_plot, "bipolar_mm_clinical_csh_plot", overwrite = T)

bipolar_mm_clinical_csh_plot.tmp <- bipolar_mm_clinical_csh_plot %>%
  filter(!grepl("cov", term)) %>%
  filter(!grepl("(Intercept)", term)) %>%
  filter(!(analysis == "joint_ICDprevMDS" & term != "prev_MDSprevMDS")) %>%
  filter(!(analysis == "joint_CGI-S" & term != "severitycgis_max_d07")) %>%
  mutate(term = case_when(term == "severitycgis_max_d07" ~ "severity: CGI-S score +1",
    term == "prev_MDSprevMDS" ~ "prev short MDTT (<1m)",
    #term == "mm_splitmixed" ~ "mixed ep (vs mania)",
    term == "hospital01" ~ "hospital for mania/mixed", 
    .default = term))

bipolar_mm_clinical_csh_plot.tmp %>%
  filter(facet != "joint model") %>%
  mutate(facet = factor(facet, levels = c("no covariates", "covariates", "joint model"))) %>%
  mutate(term2 = factor(term, levels = unique(bipolar_mm_clinical_csh_plot.tmp$term))) %>%
  ggplot(., aes(x = term2, y = estimate_HR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerHR, ymax = higherHR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 2.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("HR estimate for time to depression switch") +
  coord_flip(ylim = c(0,5)) +
  facet_grid(facet~mm_split)
ggsave("bipolar_mm_clinical_csh_plot.png", width = 7, height = 6)

#------------------------------------------------------------
# Test PH assumptions - VIOLATED
#------------------------------------------------------------
ph_hospitalmania <- bipolar_mm_clinical_csh.prep %>%
  filter(hospital_include == 1) %>%
  filter(mm_split == "mania") %>%
  coxph(Surv(start_csh, end_csh, status == 1) ~ hospital01 + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)
cox.zph(ph_hospitalmania) # <2e-16
ph_hospitalmixed <- bipolar_mm_clinical_csh.prep %>%
  filter(hospital_include == 1) %>%
  filter(mm_split == "mixed") %>%
  coxph(Surv(start_csh, end_csh, status == 1) ~ hospital01 + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)
cox.zph(ph_hospitalmixed) # <2e-16

ph_prevMDSmania <- bipolar_mm_clinical_csh.prep %>%
  filter(mm_episode >= 2) %>%
  filter(mm_split == "mania") %>%
  coxph(Surv(start_csh, end_csh, status == 1) ~ prev_MDS + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)
cox.zph(ph_prevMDSmania) # 7.1e-11
ph_prevMDSmixed <- bipolar_mm_clinical_csh.prep %>%
  filter(mm_episode >= 2) %>%
  filter(mm_split == "mixed") %>%
  coxph(Surv(start_csh, end_csh, status == 1) ~ prev_MDS + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .)
cox.zph(ph_prevMDSmixed) # 2.3e-12

ph_cgismania <- bipolar_mm_clinical_csh.prep %>%
  filter(mm_split == "mania") %>%
  dplyr::select(-c("cov_state2", "n")) %>%
  left_join(., state_tmp, relationship = "many-to-many") %>%
  coxph(Surv(start_csh, end_csh, status == 1) ~ severitycgis_max_d07 + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_rel + cov_edu2 + cov_employ2 + cov_race2, data = .) # singularity so remove cov_state2
cox.zph(ph_cgismania) # 1.5e-08
ph_cgismixed <- bipolar_mm_clinical_csh.prep %>%
  filter(mm_split == "mixed") %>%
  dplyr::select(-c("cov_state2", "n")) %>%
  left_join(., state_tmp, relationship = "many-to-many") %>%
  coxph(Surv(start_csh, end_csh, status == 1) ~ severitycgis_max_d07 + cluster(person_id) + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender, data = .) # singularity so remove cov_state2 + cov_edu2 + cov_employ2 + cov_rel +  cov_race2
cox.zph(ph_cgismixed) # < 2e-16

cox_zph_allbipolar_mania_clinical_summary <- rbind(data.frame(cox.zph(ph_hospitalmania)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mania", analysis = "hospital01"),
  data.frame(cox.zph(ph_prevMDSmania)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mania", analysis = "prevMDS"),
  data.frame(cox.zph(ph_cgismania)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mania", analysis = "cgis")) %>%
  filter(analysis != "GLOBAL") %>%
  mutate(mm_split = "mania", sample = "bipolar_all")
cox_zph_allbipolar_mixed_clinical_summary <- rbind(data.frame(cox.zph(ph_hospitalmixed)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mixed", analysis = "hospital01"),
  data.frame(cox.zph(ph_prevMDSmixed)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mixed", analysis = "prevMDS"),
  data.frame(cox.zph(ph_cgismixed)$table) %>% rownames_to_column(var = "term") %>% mutate(mm_split = "mixed", analysis = "cgis")) %>%
  filter(analysis != "GLOBAL") %>%
  mutate(mm_split = "mixed", sample = "bipolar_all")
saveDf(cox_zph_allbipolar_mania_clinical_summary, "cox_zph_allbipolar_mania_clinical_summary", overwrite = T)
saveDf(cox_zph_allbipolar_mixed_clinical_summary, "cox_zph_allbipolar_mixed_clinical_summary", overwrite = T)

#======================================================
# Sensitivity analyses
#======================================================

#---------------------------
# Analysis with narrow ICD-only definition
#---------------------------

# Identify which episodes should switch from MDS1m == 1 to MDS1m == 0
mania_dep_dates_toepisodes_qc_defbroad <- getDf("mania_dep_dates_toepisodes_qc_defbroad")
mania_dep_dates_toepisodes_qc_defbroad_MDS1m_dx <- mania_dep_dates_toepisodes_qc_defbroad %>%
  mutate(months_timetodep_dx = as.numeric(difftime(dep_date_dx, diagnosis_date_end, units = "days")/30)) %>%
  mutate(MDS1m_dx = ifelse(months_timetodep_dx <= 1, 1, 0)) %>%
  mutate(MDS1m_dx = ifelse(is.na(MDS1m_dx), 0, MDS1m_dx)) %>%
  filter(MDS1m_dx == 1) %>%
  dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, mm_episode, dep_episode, MDS1m_dx)

# Redo clinical features analysis
bipolar_mm_clinical <- getDf("bipolar_mm_clinical")
bipolar_mm_clinical_MDS1m_dx <- bipolar_mm_clinical %>% 
  left_join(., mania_dep_dates_toepisodes_qc_defbroad_MDS1m_dx) %>%
  mutate(MDS1m_dx = ifelse(is.na(MDS1m_dx), 0, MDS1m_dx))
bipolar_mm_clinical_cgis_MDS1m_dx <- bipolar_mm_clinical_cgis %>% 
  left_join(., mania_dep_dates_toepisodes_qc_defbroad_MDS1m_dx) %>%
  mutate(MDS1m_dx = ifelse(is.na(MDS1m_dx), 0, MDS1m_dx))

# - COV
bipolar_mm_hospital_cov_MDS1m_dx_glmrr <- bipolar_mm_clinical_MDS1m_dx %>%
  filter(hospital_include == 1) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m_dx ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_cgis_cov_MDS1m_dx_glmrr <- bipolar_mm_clinical_cgis_MDS1m_dx %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m_dx ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_prevMDS_cov_MDS1m_dx_glmrr <- bipolar_mm_clinical_MDS1m_dx %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m_dx ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_glmrr_cov_MDS1m_dx <- rbind(bipolar_mm_clinical_prevMDS_cov_MDS1m_dx_glmrr,
  bipolar_mm_cgis_cov_MDS1m_dx_glmrr,
  bipolar_mm_hospital_cov_MDS1m_dx_glmrr) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(p_adj = p.value*3) %>% # 3 tests for each of prevMDS, severity, hospital
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(facet = "covariates")

# PLOT
bipolar_mm_clinical_glmrr_cov_MDS1m_dx.tmp <- bipolar_mm_clinical_glmrr_cov_MDS1m_dx %>%
  mutate(term = case_when(term == "severitycgis_max_d07" ~ "severity: CGI-S score +1",
    term == "prev_MDSprevMDS" ~ "prev short MDTT (<1m)",
    #term == "mm_splitmixed" ~ "mixed ep (vs mania)",
    term == "hospital01" ~ "hospital for mania/mixed",
    .default = term))

# For final figure
bipolar_mm_clinical_glmrr_plot_covfinal_MDS1m_dx <- bipolar_mm_clinical_glmrr_cov_MDS1m_dx.tmp %>%
  filter(!grepl("icdseverity", term)) %>%
  filter(facet == "covariates") %>%
  mutate(title = "Risk Ratio of Switch in 1m") %>%
#bipolar_mm_clinical_all_glmrr_plot_gg <- bipolar_mm_clinical_glmrr_cov.plot %>%
#  full_join(., bipolar_dx7d_dx_medsbroadclass_mmsplit_joint_cov_glmrr.plot) %>%
  mutate(term = ifelse(term == "previous Dep switch (1m)", "prev short MDTT (<1m)", term)) %>%
  mutate(term2 = factor(term, levels = rev(c("prev short MDTT (<1m)", "severity: CGI-S score +1", "hospital for mania/mixed", "antidepressant", "antiepileptic", "antipsychotic", "benzodiazepine", "lithium", "antihistamine")))) %>%
  ggplot(., aes(x = term2, y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR), shape = 17) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_label_repel(size = 2.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  #scale_y_continuous(limits = c(0.4, 2.5)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate, ICD depression only (adj. covariates)") +
  coord_flip() +
  facet_grid(facet ~ mm_split, scales = "free",  space = "free_y")
ggsave("bipolar_mm_clinical_glmrr_plot_covfinal_MDS1m_dx.png", bipolar_mm_clinical_glmrr_plot_covfinal_MDS1m_dx, width = 7, height = 3, dpi = 600)

#####
# Identify which episodes should switch from MDS1m == 1 to MDS1m == 0
mania_dep_dates_toepisodes_qc_defbroad_noad <- getDf("mania_dep_dates_toepisodes_qc_defbroad_noad")
mania_dep_dates_toepisodes_qc_defbroad_MDS1m_noad <- mania_dep_dates_toepisodes_qc_defbroad_noad %>%
  mutate(MDS1m_noad = ifelse(timetodep_months <= 1, 1, 0)) %>%
  mutate(MDS1m_noad = ifelse(is.na(MDS1m_noad), 0, MDS1m_noad)) %>%
  filter(MDS1m_noad == 1) %>%
  dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, mm_episode, dep_episode, MDS1m_noad)

# Redo clinical features analysis
bipolar_mm_clinical <- getDf("bipolar_mm_clinical")
bipolar_mm_clinical_MDS1m_noad <- bipolar_mm_clinical %>% 
  left_join(., mania_dep_dates_toepisodes_qc_defbroad_MDS1m_noad) %>%
  mutate(MDS1m_noad = ifelse(is.na(MDS1m_noad), 0, MDS1m_noad))
bipolar_mm_clinical_cgis_MDS1m_noad <- bipolar_mm_clinical_cgis %>% 
  left_join(., mania_dep_dates_toepisodes_qc_defbroad_MDS1m_noad) %>%
  mutate(MDS1m_noad = ifelse(is.na(MDS1m_noad), 0, MDS1m_noad))

# - COV
bipolar_mm_hospital_cov_MDS1m_noad_glmrr <- bipolar_mm_clinical_MDS1m_noad %>%
  filter(hospital_include == 1) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m_noad ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_cgis_cov_MDS1m_noad_glmrr <- bipolar_mm_clinical_cgis_MDS1m_noad %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m_noad ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_prevMDS_cov_MDS1m_noad_glmrr <- bipolar_mm_clinical_MDS1m_noad %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m_noad ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_glmrr_cov_MDS1m_noad <- rbind(bipolar_mm_clinical_prevMDS_cov_MDS1m_noad_glmrr,
  bipolar_mm_cgis_cov_MDS1m_noad_glmrr,
  bipolar_mm_hospital_cov_MDS1m_noad_glmrr) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(p_adj = p.value*3) %>% # 3 tests for each of prevMDS, severity, hospital
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(facet = "covariates")

# PLOT
bipolar_mm_clinical_glmrr_cov_MDS1m_noad.tmp <- bipolar_mm_clinical_glmrr_cov_MDS1m_noad %>%
  mutate(term = case_when(term == "severitycgis_max_d07" ~ "severity: CGI-S score +1",
    term == "prev_MDSprevMDS" ~ "prev short MDTT (<1m)",
    #term == "mm_splitmixed" ~ "mixed ep (vs mania)",
    term == "hospital01" ~ "hospital for mania/mixed",
    .default = term))

# For final figure
bipolar_mm_clinical_glmrr_plot_covfinal_MDS1m_noad <- bipolar_mm_clinical_glmrr_cov_MDS1m_noad.tmp %>%
  filter(!grepl("icdseverity", term)) %>%
  filter(facet == "covariates") %>%
  mutate(title = "Risk Ratio of Switch in 1m") %>%
#bipolar_mm_clinical_all_glmrr_plot_gg <- bipolar_mm_clinical_glmrr_cov.plot %>%
#  full_join(., bipolar_dx7d_dx_medsbroadclass_mmsplit_joint_cov_glmrr.plot) %>%
  mutate(term = ifelse(term == "previous Dep switch (1m)", "prev short MDTT (<1m)", term)) %>%
  mutate(term2 = factor(term, levels = rev(c("prev short MDTT (<1m)", "severity: CGI-S score +1", "hospital for mania/mixed", "antidepressant", "antiepileptic", "antipsychotic", "benzodiazepine", "lithium", "antihistamine")))) %>%
  ggplot(., aes(x = term2, y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR), shape = 17) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_label_repel(size = 2.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  #scale_y_continuous(limits = c(0.4, 2.5)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate (adj. cov, antidepressant sensitivity)") +
  coord_flip() +
  facet_grid(facet ~ mm_split, scales = "free",  space = "free_y")
ggsave("bipolar_mm_clinical_glmrr_plot_covfinal_MDS1m_noad.png", bipolar_mm_clinical_glmrr_plot_covfinal_MDS1m_noad, width = 7, height = 3, dpi = 600)

#---------------------------
# State stratified differences
#---------------------------

#cov_state_keep_analysis_clinical0 <- bipolar_mm_clinical %>%
#  group_by(mm_split, cov_state) %>%
#  summarise(cov_state_n=n()) %>%
#  mutate(cov_state3 = ifelse(cov_state == "0No matching concept", "N/A", cov_state)) %>%
#  mutate(cov_state3 = ifelse(cov_state_n < 1000, "group", cov_state3)) %>%
#  ungroup() 

#cov_state_keep_analysis_clinical <- cov_state_keep_analysis_clinical0 %>%
#  group_by(mm_split, cov_state3) %>%
#  summarise(cov_state_n2=sum(cov_state_n)) %>%
#  mutate(cov_state_new = paste0(cov_state3, "(", cov_state_n2, ")")) %>%
#  left_join(., cov_state_keep_analysis_clinical0 %>% dplyr::select("mm_split", "cov_state", "cov_state3"), relationship = "many-to-many") %>% 
#  distinct()

cov_state_keep_analysis_clinical <- bipolar_mm_clinical %>%
  group_by(mm_split, cov_state2) %>%
  summarise(cov_state_n2=n()) %>%
  mutate(cov_state_new = paste0(cov_state2, "(", cov_state_n2, ")")) %>%
  #left_join(., cov_state_keep_analysis_clinical0 %>% dplyr::select("mm_split", "cov_state", "cov_state3"), relationship = "many-to-many") %>% 
  distinct()

bipolar_mm_clinical_state <- bipolar_mm_clinical %>% left_join(., cov_state_keep_analysis_clinical)

cov_state_keep_analysis_clinical_cgis0 <- bipolar_mm_clinical_cgis %>%
  filter(!is.na(severitycgis_max_d07)) %>%
  group_by(mm_split, cov_state) %>%
  summarise(cov_state_n=n()) %>%
  mutate(cov_state2 = ifelse(cov_state == "0No matching concept", "N/A", cov_state)) %>%
  mutate(cov_state2 = ifelse(cov_state_n < 100, "group", cov_state2)) %>%
  ungroup() 

cov_state_keep_analysis_clinical_cgis <- cov_state_keep_analysis_clinical_cgis0 %>%
  group_by(mm_split, cov_state2) %>%
  summarise(cov_state_n2=sum(cov_state_n)) %>%
  mutate(cov_state_new = paste0(cov_state2, "(", cov_state_n2, ")")) %>%
  left_join(., cov_state_keep_analysis_clinical_cgis0 %>% dplyr::select("mm_split", "cov_state", "cov_state2"), relationship = "many-to-many") %>% 
  distinct()

bipolar_mm_clinical_cgis_state <- bipolar_mm_clinical_cgis %>% inner_join(., cov_state_keep_analysis_clinical_cgis)

# - COV
bipolar_mm_hospital_cov_glmrr_state <- bipolar_mm_clinical_state %>%
  filter(hospital_include == 1) %>%
  nest(.by = c("mm_split", "cov_state_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "Hospital")

bipolar_mm_cgis_cov_glmrr_state <- bipolar_mm_clinical_cgis_state %>%
  nest(.by = c("mm_split", "cov_state_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "CGI-S")

bipolar_mm_clinical_prevMDS_cov_glmrr_state <- bipolar_mm_clinical_state %>%
  filter(mm_episode >= 2) %>%
  nest(.by = c("mm_split", "cov_state_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "Prev MDTT (<1m)")

bipolar_mm_clinical_glmrr_cov_state <- rbind(bipolar_mm_clinical_prevMDS_cov_glmrr_state,
  bipolar_mm_cgis_cov_glmrr_state,
  bipolar_mm_hospital_cov_glmrr_state) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(p_adj = p.value*3) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(facet = "covariates")

# Plots
bipolar_mm_clinical_glmrr_cov_state_mania.plot <- bipolar_mm_clinical_glmrr_cov_state %>%
  filter(!is.na(cov_state_new)) %>%
  filter(mm_split == "mania") %>% filter(!term %in% c("cov_age_ep", "cov_age_ehr", "cov_years_ehr", "cov_gender")) %>%
  ggplot(aes(x = fct_reorder(cov_state_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (adj. for cov)") +
  coord_flip(ylim = c(0,6)) +
  facet_grid(analysis ~ mm_split, scales = "free_y", space = "free_y")

bipolar_mm_clinical_glmrr_cov_state_mixed.plot <- bipolar_mm_clinical_glmrr_cov_state %>%
  filter(!is.na(cov_state_new)) %>%
  filter(mm_split == "mixed") %>% filter(!term %in% c("cov_age_ep", "cov_age_ehr", "cov_years_ehr", "cov_gender")) %>%
  ggplot(aes(x = fct_reorder(cov_state_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (adj. for cov)") +
  coord_flip(ylim = c(0,6)) +
  facet_grid(analysis ~ mm_split, scales = "free_y", space = "free_y")

library(patchwork)
(bipolar_mm_clinical_glmrr_cov_state_mania.plot | bipolar_mm_clinical_glmrr_cov_state_mixed.plot) +
  plot_annotation(tag_levels = 'a') + plot_layout(axes = "collect")
ggsave("bipolar_mm_clinical_glmrr_cov_state.png", height = 11, width = 10)

#---------------------------
# Race stratified differences
#---------------------------

cov_race_keep_analysis_clinical0 <- bipolar_mm_clinical %>%
  group_by(mm_split, cov_race2) %>%
  summarise(cov_race2_n=n()) %>%
  mutate(cov_race_new = paste0(cov_race2, "(", cov_race2_n, ")"))

bipolar_mm_clinical_race <- bipolar_mm_clinical %>% left_join(., cov_race_keep_analysis_clinical0)

cov_race_keep_analysis_clinical_cgis0 <- bipolar_mm_clinical_cgis %>%
  filter(!is.na(severitycgis_max_d07)) %>%
  group_by(mm_split, cov_race2) %>%
  summarise(cov_race2_n=n()) %>%
  mutate(cov_race_new = paste0(cov_race2, "(", cov_race2_n, ")"))

bipolar_mm_clinical_cgis_race <- bipolar_mm_clinical_cgis %>% inner_join(., cov_race_keep_analysis_clinical_cgis0)

# - COV
bipolar_mm_hospital_cov_glmrr_race <- bipolar_mm_clinical_race %>%
  nest(.by = c("mm_split", "cov_race_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "Hospital")

bipolar_mm_cgis_cov_glmrr_race <- bipolar_mm_clinical_cgis_race %>%
  nest(.by = c("mm_split", "cov_race_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "CGI-S severity")

bipolar_mm_clinical_prevMDS_cov_glmrr_race <- bipolar_mm_clinical_race %>%
  filter(mm_episode >= 2) %>%
  nest(.by = c("mm_split", "cov_race_new")) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "Prev MDTT (<1m)")

bipolar_mm_clinical_glmrr_cov_race <- rbind(bipolar_mm_clinical_prevMDS_cov_glmrr_race,
  bipolar_mm_cgis_cov_glmrr_race,
  bipolar_mm_hospital_cov_glmrr_race) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  mutate(p_adj = p.value*3) %>%
  mutate(sig = case_when(#p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(facet = "covariates")

# Plots
bipolar_mm_clinical_glmrr_cov_race_mania.plot <- bipolar_mm_clinical_glmrr_cov_race %>%
  filter(!is.na(cov_race_new)) %>%
  filter(mm_split == "mania") %>% filter(!term %in% c("cov_age_ep", "cov_age_ehr", "cov_years_ehr", "cov_gender")) %>%
  ggplot(aes(x = fct_reorder(cov_race_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (adj. for cov)") +
  coord_flip(ylim = c(0,6)) +
  facet_grid(analysis ~ mm_split, scales = "free_y", space = "free_y")

bipolar_mm_clinical_glmrr_cov_race_mixed.plot <- bipolar_mm_clinical_glmrr_cov_race %>%
  filter(!is.na(cov_race_new)) %>%
  filter(mm_split == "mixed") %>% filter(!term %in% c("cov_age_ep", "cov_age_ehr", "cov_years_ehr", "cov_gender")) %>%
  ggplot(aes(x = fct_reorder(cov_race_new, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (adj. for cov)") +
  coord_flip(ylim = c(0,6)) +
  facet_grid(analysis ~ mm_split, scales = "free_y", space = "free_y")

library(patchwork)
(bipolar_mm_clinical_glmrr_cov_race_mania.plot | bipolar_mm_clinical_glmrr_cov_race_mixed.plot) +
  plot_annotation(tag_levels = 'a') + plot_layout(axes = "collect")
ggsave("bipolar_mm_clinical_glmrr_cov_race.png", height = 7, width = 10)

#======================================================
#
# MM severity --> Dep severity?
#
#======================================================

mania_dep_dates_toepisodes_qc_broaddef_cgis_mmend <- getDf("mania_dep_dates_toepisodes_qc_broaddef_cgis_mmend")
mm_cgis_dep_cgis <- mania_dep_dates_toepisodes_qc_broaddef_cgis_mmend %>%
  inner_join(., bipolar_mm_clinical_cgis %>% dplyr::select(person_id, mm_episode, dep_episode, severitycgis_max_d07)) %>%
  filter(!is.na(severitycgis_max_d07) & !is.na(value)) %>%
  filter(months_postmm > 0)

ggplot(mm_cgis_dep_cgis %>% mutate(episode_type = mm_split), aes(x = as.factor(severitycgis_max_d07), y = value, fill = episode_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(colour = "gray50", alpha = 0.7) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "CGI-S mania/mixed (maximum in first 7 days)", y = "CGI-S depression (maximum)") +
  scale_fill_manual(values = millerstone[c(7,4)])
ggsave("mm_cgis_dep_cgis.png", height = 4, width = 7)

cgis_mm_dep_cor <- cor.test(mm_cgis_dep_cgis$value, mm_cgis_dep_cgis$severitycgis_max_d07)
print(cgis_mm_dep_cor)
print(cgis_mm_dep_cor$p.value)

#======================================================
#
# MOVED from predictdep_mm_medications_pmmdgroup.r
# [tidy up for this script OR put into code_graveyard.r]
#======================================================

# CGI-S differences between mania/mixed groups
pmmdgroup_hosp_cgis <- cgis_bipolar %>% filter(person_id %in% pmmdgroup_mm_7d_hosp_psychmeds_count$person_id)

pmmdgroup_mm_7d_hosp_psychmeds_count_cgis <- pmmdgroup_hosp_cgis %>%
  mutate(person_id = as.character(person_id)) %>%
  rename(cgis = value) %>%
  left_join(pmmdgroup_mm_7d_hosp_psychmeds_count, ., relationship = "many-to-many") %>%
  mutate(time_dx_cgis = as.numeric(difftime(diagnosis_date, measurement_date, units = "days"))) %>%
  filter(time_dx_cgis %in% -7:28) %>%
  arrange(time_dx_cgis) %>%
  group_by(person_id, visit_start_date, measurement_date) %>%
  slice_head(n=1) %>%
  group_by(person_id, visit_start_date) %>%
  mutate(timefromobs1 = as.numeric(difftime(measurement_date, min(measurement_date), units = "days")))

# - check how many observations per person per episode, median = 7!
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis %>% 
  group_by(person_id, visit_start_date) %>% 
  summarise(n_obs = n()) %>%
  pull(n_obs) %>% summary

# - derive first score, max score, instability score (only have for n=342)
tRMSSD <- function(times,values){
  n=length(values)
  dtime = times[2:n]-times[1:n-1]
  dvals = values[2:n]-values[1:n-1]
  res = sqrt(1/n*sum((dvals/dtime)^2))
  return(res)
}

pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep <- pmmdgroup_mm_7d_hosp_psychmeds_count_cgis %>%
  arrange(measurement_date) %>%
  group_by(person_id, visit_start_date) %>%
  mutate(max = max(cgis), first = cgis[1], tRMSSD = tRMSSD(timefromobs1, cgis), n_cgis = n()) %>%
  dplyr::select(-c("measurement_date", "cgis", "time_dx_cgis", "timefromobs1")) %>%
  distinct()

# - check differences in each measure between CGI-S and mm_split: 
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>% 
  filter(mm_split != "both") %>%
  pivot_longer(cols = c("max", "first", "tRMSSD", "n_cgis"), names_to = "metric", values_to = "cgis_measure") %>%
  mutate(PMMDep = ifelse(pmmd01 == 0, "nonPMMDep", "PMMDep")) %>%
ggplot(., aes(x = mm_split, y = cgis_measure, fill = mm_split)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, colour = "gray30", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(1,4)]) +
  ylab("CGI-S metric") +
  theme(legend.position = "bottom") +
  facet_wrap(~ metric, scales = "free_y")
ggsave("pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_mmsplit.png")

# - check differences in each measure between CGI-S and PMMDep vs nonPMMDep: 
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>% 
  filter(mm_split != "both") %>%
  pivot_longer(cols = c("max", "first", "tRMSSD", "n_cgis"), names_to = "metric", values_to = "cgis_measure") %>%
  mutate(PMMDep = ifelse(pmmd01 == 0, "nonPMMDep", "PMMDep")) %>%
ggplot(., aes(x = mm_split, y = cgis_measure, fill = PMMDep)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, colour = "gray30", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  ylab("CGI-S metric") +
  theme(legend.position = "bottom") +
  facet_wrap(~ metric, scales = "free_y")
ggsave("pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_pmmd01.png")

# - check differences in each measure between CGI-S and prescribing: 
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>% 
  filter(mm_split != "both") %>%
  pivot_longer(cols = c("max", "first", "tRMSSD", "n_cgis"), names_to = "metric", values_to = "cgis_measure") %>%
  pivot_longer(cols = c("lithium", "antiepileptic", "antipsychotic", "antidepressant"), names_to = "broad_class", values_to = "prescribed") %>%
  mutate(prescribed = as.character(prescribed)) %>%
  mutate(PMMDep = ifelse(pmmd01 == 0, "nonPMMDep", "PMMDep")) %>%
ggplot(., aes(x = broad_class, y = cgis_measure, fill = prescribed)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, size = 1, colour = "gray30", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(2,3)]) +
  ylab("CGI-S metric") +
  theme(legend.position = "bottom") +
  facet_grid(metric ~ mm_split, scales = "free_y")
ggsave("pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_broadclassmmsplit.png", height = 8, width = 8)

# - check differences in each measure between PMMD0/1 and mm_split: 
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>% 
  filter(mm_split != "both") %>%
  pivot_longer(cols = c("max", "first", "tRMSSD", "n_cgis"), names_to = "metric", values_to = "cgis_measure") %>%
  mutate(PMMDep = ifelse(pmmd01 == 0, "nonPMMDep", "PMMDep")) %>%
ggplot(., aes(x = mm_split, y = cgis_measure, fill = PMMDep)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, colour = "gray30", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  ylab("CGI-S metric") +
  theme(legend.position = "bottom") +
  facet_wrap(~ metric, scales = "free_y")
ggsave("pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_pmmd01mmsplit.png")

# - check differences in each measure between PMMD01 and medication types: 
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>% 
  filter(mm_split != "both") %>%
  pivot_longer(cols = c("max", "first", "tRMSSD", "n_cgis"), names_to = "metric", values_to = "cgis_measure") %>%
  pivot_longer(cols = c("lithium", "antiepileptic", "antipsychotic", "antidepressant"), names_to = "broad_class", values_to = "prescribed") %>%
  mutate(prescribed = as.character(prescribed)) %>%
  mutate(PMMDep = ifelse(pmmd01 == 0, "nonPMMDep", "PMMDep")) %>%
ggplot(., aes(x = prescribed, y = cgis_measure, fill = PMMDep)) +
  geom_boxplot(position = "dodge") +
  geom_point(alpha = 0.2, size = 1, colour = "gray10", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +  
  ylab("CGI-S metric") +
  theme(legend.position = "bottom") +
  facet_grid(metric ~ broad_class, scales = "free")
ggsave("pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_pmmd01broadclass.png")

# - check differences in each measure between mm_split and medication types: 
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>% 
  filter(mm_split != "both") %>%
  pivot_longer(cols = c("max", "first", "tRMSSD", "n_cgis"), names_to = "metric", values_to = "cgis_measure") %>%
  pivot_longer(cols = c("lithium", "antiepileptic", "antipsychotic", "antidepressant"), names_to = "broad_class", values_to = "prescribed") %>%
  mutate(prescribed = as.factor(prescribed)) %>%
ggplot(., aes(x = prescribed, y = cgis_measure, fill = mm_split)) +
  geom_boxplot(position = "dodge") +
  geom_point(alpha = 0.2, size = 1, colour = "gray10", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(1,4)]) +  
  ylab("CGI-S metric") +
  theme(legend.position = "bottom") +
  facet_grid(metric ~ broad_class, scales = "free")
ggsave("pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_mmsplitbroadclass.png")

# - logistic regression for first CGI-S ~ antidepressants, lithium, antiepileptic (based on plots above)
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_first.lm <- pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>%
  filter(mm_split != "both") %>%
  ungroup() %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~lm(first ~ lithium + antidepressant + antiepileptic, data = .)), tidied = map(model, tidy)) %>% 
  unnest(tidied)
print(pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_first.lm)

pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_tRMSSD.lm <- pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>%
  filter(mm_split != "both") %>%
  ungroup() %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~lm(tRMSSD ~ lithium + antidepressant + antiepileptic, data = .)), tidied = map(model, tidy)) %>% 
  unnest(tidied)
print(pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep_tRMSSD.lm)

# - logistic regression on PMMD01 ~ CGI-S measures
pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep.glm <- pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep %>%
  filter(mm_split != "both") %>%
  ungroup() %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~glm(pmmd01 ~ first + max + tRMSSD, family = 'binomial', data = .)), tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  mutate(estimate_OR = exp(estimate),
        lowerOR = exp(estimate - (std.error * 1.96)), 
        higherOR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste("OR=", sprintf('%.2f', estimate_OR), ", ", 
                                      "SE=", sprintf('%.2f', std.error), ", ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA))

pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep.plot <- pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep.glm %>%
  filter(term != "(Intercept)") %>%
  #left_join(., colours.tmp, relationship = "many-to-many") %>%
  ggplot(aes(x = term, y = estimate_OR, colour = term, label = lab)) +
  geom_pointrange(aes(ymin = lowerOR, ymax = higherOR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_label_repel(size = 3.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_manual(values = millerstone[c(9,10,11)]) +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("OR estimate") +
  coord_flip() +
  facet_wrap(~mm_split)
ggsave("pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep.png", pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep.plot, width = 6, height = 4)

print(pmmdgroup_mm_7d_hosp_psychmeds_count_cgis_perep.glm)
