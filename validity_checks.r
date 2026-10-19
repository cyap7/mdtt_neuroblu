#======================================================
# Libraries and functions
#======================================================

source("Setup_background.r")
source("fastPhenotype.r")
source("read_in_routine.r")

library(tidyverse)
library(forcats)
library(lme4)
library(viridis)
library(ggrepel)

library(MASS)
library(broom)
library(zoo) # for loca()
library(pracma) # for trapz()
library(gee)

library(gtsummary)

# Palettes
millerstone <- c("#4F6980FF", "#849DB1FF", "#A2CEAAFF", "#638B66FF", "#BFBB60FF", "#F47942FF", "#FBB04EFF", "#B66353FF", "#D7CE9FFF", "#B9AA97FF", "#7E756DFF", "gray70")
antique <- c("#855C75FF", "#D9AF6BFF", "#AF6458FF", "#736F4CFF", "#526A83FF", "#625377FF", "#68855CFF", "#9C9C5EFF", "#A06177FF", "#8C785DFF", "#467378FF", "#7C7C7CFF")

mm_split <- getDf("mm_split")
mm_split_fullep <- getDf("mm_split_fullep")

pmmdgroup_visit <- getDf("pmmdgroup_visit")
nopmmdgroup_visit <- getDf("nopmmdgroup_visit")
mmonlygroup_visit <- getDf("mmonlygroup_visit")

#============================================================
#
# Group differences: demographics
#
#============================================================

group_summary_table = NULL

bipolar_epall <- getDf("bipolar_epall")

##id_group <- getDf("id_group") %>% filter(is.na(samedayswitch_anyorder)) # EXCLUDE samedayswitch to avoid confusion
##id_group_dx <- getDf("id_group_dx")
# this id_group version excludes samedayswitch_anyorder and also episodes with mm and dep on the same day
bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")
id_group_clean <- bipolar_epall_count_mm_to_depafter_inf %>% 
  mutate(MTTD_1m_01 = ifelse(mm_to_depafter <= 1, 1, 0)) %>%
  group_by(person_id) %>%
  summarise(MTTD_1m.tmp = max(MTTD_1m_01)) %>%
  mutate(MTTD_1m.tmp = ifelse(is.na(MTTD_1m.tmp), 0, MTTD_1m.tmp)) %>%
  mutate(group = ifelse(MTTD_1m.tmp == 1, "MTTD_less1m",  "MTTD_more1m")) %>%
  dplyr::select(person_id, group)

# id_group_dx is equivalent to switching all of those with addpmmd == "PMMDadd" to MTTD_more1m
# - have checked this against mania_dep_dates_toepisodes_qc_defbroad
id_group_full_pmmdadd <- getDf("id_group") %>% filter(addpmmd == "PMMDadd")
id_group_clean_dx <- id_group_clean %>%
  mutate(group = ifelse(person_id %in% id_group_full_pmmdadd$person_id, "MTTD_more1m", group))

n_patients <- id_group_clean %>% 
  #filter(is.na(samedayswitch_anyorder)) %>%
  group_by(group) %>% summarise(n=n()) %>% 
  column_to_rownames(., var = "group") %>%
  mutate(n = as.character(n)) %>%
  t() %>% as.data.frame

#------------------------------------------------------------
# Compare number of people in each group
#------------------------------------------------------------



#------------------------------------------------------------
# Compare demographics
#------------------------------------------------------------

demog_cov <- getDf("demog_cov")

age_firstmoodep.tmp <- bipolar_epall %>% 
  group_by(person_id) %>%
  summarise(firstmoodep_date = min(diagnosis_date))

demog_compare <- demog_cov %>%
  #left_join(., observation_period) %>%
  left_join(., age_firstmoodep.tmp) %>%
  #mutate(observation_period = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days"))/365) %>%
  mutate(birth_date = as.Date(paste(birth_year, "-01-01", sep = ""))) %>%
  mutate(age_entry = as.numeric(difftime(observation_period_start_date, birth_date, units = "days"))/365) %>%
  mutate(age_firstmoodep = as.numeric(difftime(firstmoodep_date, birth_date, units = "days"))/365) %>%
  mutate(age = as.numeric(difftime("2025-01-01", birth_date, units = "days"))/365) %>%
  rename(gender = cov_gender) %>%
  left_join(id_group_clean_dx %>% dplyr::select(person_id, group) %>% rename(group_dx = group)) %>%
  left_join(id_group_clean) %>%
  filter(!is.na(group)) 

#---------------
# - plots for birth date, age_firstmoodep, sex, observation years
yob_plot <- ggplot(demog_compare, aes(x = birth_date)) +
  geom_histogram(fill = millerstone[9], alpha = 0.5) +
  theme_bw() +
  xlab("Year of birth")
agefirstmoodep <- ggplot(demog_compare, aes(x = age_firstmoodep)) +
  geom_histogram(fill = millerstone[4], alpha = 0.5) +
  theme_bw() +
  xlab("Age 1st recorded mood ep")
sex_plot <- ggplot(demog_compare, aes(x = gender, fill = gender)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(2,4,8)]) +
  xlab("Gender") + theme(legend.position = "none")
obstime_plot <- ggplot(demog_compare, aes(x = observation_years)) +
  geom_histogram(fill = millerstone[2], alpha = 0.5) +
  theme_bw() +
  xlab("Observation years")
library(patchwork)
cohort_demog <- (yob_plot + agefirstmoodep) / (sex_plot + obstime_plot) + 
            patchwork::plot_annotation(tag_levels = "a")
ggsave("cohort_demog.png", cohort_demog, height = 6, width = 8)
#---------------

# - add to summary table (compare when PMMD is defined via diagnosis only via broader definition)
age_sex <- demog_compare %>%
  group_by(group) %>%
  summarise(Age_entry_mean = paste(signif(mean(age_entry),2), " (", signif(sd(age_entry),2), ")", sep = ""),
            Age_1stmooddx_mean = paste(signif(mean(age_firstmoodep),2), " (", signif(sd(age_firstmoodep),2), ")", sep = ""),
            Age_mean = paste(signif(mean(age),2), " (", signif(sd(age),2), ")", sep = ""),
            Sex_F_pc = paste(n(), " (", round(length(which(gender == "FEMALE"))/n()*100,1), ")", sep = "")) %>%
  column_to_rownames(., var = "group") %>%
  t() %>% data.frame
age_sex_dx <- demog_compare %>%
  group_by(group_dx) %>%
  summarise(Age_entry_mean = paste(signif(mean(age_entry),2), " (", signif(sd(age_entry),2), ")", sep = ""),
            Age_1stmooddx_mean = paste(signif(mean(age_firstmoodep),2), " (", signif(sd(age_firstmoodep),2), ")", sep = ""),
            Age_mean = paste(signif(mean(age),2), " (", signif(sd(age),2), ")", sep = ""),
            Sex_F_pc = paste(n(), " (", round(length(which(gender == "FEMALE"))/n()*100,1), ")", sep = "")) %>%
  column_to_rownames(., var = "group_dx") %>%
  t() %>% data.frame 

age_entry_test <- tidy(aov(age_entry ~ group, data = demog_compare %>% filter(group %in% c("MTTD_less1m", "MTTD_more1m"))))
age_entry_test_out <- data.frame(variable = "Age at EHR entry", p.value = age_entry_test %>% filter(term == "group") %>% pull(p.value) %>% as.vector)
age_entry_test_dx <- tidy(aov(age_entry ~ group_dx, data = demog_compare %>% filter(group_dx %in% c("MTTD_less1m", "MTTD_more1m"))))
age_entry_test_out_dx <- data.frame(variable = "Age at EHR entry", p.value = age_entry_test_dx %>% filter(term == "group_dx") %>% pull(p.value) %>% as.vector)

age_1stmooddx_test <- tidy(aov(age_firstmoodep ~ group, data = demog_compare %>% filter(group %in% c("MTTD_less1m", "MTTD_more1m"))))
age_1stmooddx_test_out <- data.frame(variable = "Age at first mood diagnosis", p.value = age_1stmooddx_test %>% filter(term == "group") %>% pull(p.value) %>% as.vector)
age_1stmooddx_test_dx <- tidy(aov(age_firstmoodep ~ group_dx, data = demog_compare %>% filter(group_dx %in% c("MTTD_less1m", "MTTD_more1m"))))
age_1stmooddx_test_out_dx <- data.frame(variable = "Age at EHR entry", p.value = age_1stmooddx_test_dx %>% filter(term == "group_dx") %>% pull(p.value) %>% as.vector)

age_test <- tidy(aov(age ~ group, data = demog_compare %>% filter(group %in% c("MTTD_less1m", "MTTD_more1m"))))
age_test_lm <- tidy(lm(age ~ group, data = demog_compare %>% filter(group %in% c("MTTD_less1m", "MTTD_more1m"))))
age_test_out <- data.frame(variable = "Age", p.value = age_test %>% filter(term == "group") %>% pull(p.value) %>% as.vector)
age_test_dx <- tidy(aov(age ~ group_dx, data = demog_compare %>% filter(group_dx %in% c("MTTD_less1m", "MTTD_more1m"))))
age_test_out_dx <- data.frame(variable = "Age", p.value = age_test_dx %>% filter(term == "group_dx") %>% pull(p.value) %>% as.vector)

sex_test <- demog_compare %>% group_by(group) %>% summarise(M = length(which(gender == "MALE")), F = length(which(gender == "FEMALE"))) %>% filter(group %in% c("MTTD_less1m", "MTTD_more1m")) %>% column_to_rownames(., var = "group") %>% as.matrix
sex_test_out <- data.frame(variable = "Sex", p.value = tidy(chisq.test(sex_test)) %>% pull(p.value) %>% as.vector)
sex_test_dx <- demog_compare %>% group_by(group_dx) %>% summarise(M = length(which(gender == "MALE")), F = length(which(gender == "FEMALE"))) %>% filter(group_dx %in% c("MTTD_less1m", "MTTD_more1m")) %>% column_to_rownames(., var = "group_dx") %>% as.matrix
sex_test_out_dx <- data.frame(variable = "Sex", p.value = tidy(chisq.test(sex_test_dx)) %>% pull(p.value) %>% as.vector)

# - output
age_sex$p_MTTD_less1m <- NA
age_sex[which(rownames(age_sex) == "Age_entry_mean"),]$p_MTTD_less1m <- age_entry_test_out$p.value
age_sex[which(rownames(age_sex) == "Age_1stmooddx_mean"),]$p_MTTD_less1m <- age_1stmooddx_test_out$p.value
age_sex[which(rownames(age_sex) == "Age_mean"),]$p_MTTD_less1m <- age_test_out$p.value
age_sex[which(rownames(age_sex) == "Sex_F_pc"),]$p_MTTD_less1m <- sex_test_out$p.value
age_sex$comment <- ""
age_sex[which(rownames(age_sex) == "Age_entry_mean"),]$comment <- "Age of entry into EHR"
age_sex[which(rownames(age_sex) == "Age_1stmooddx_mean"),]$comment <- "Age of 1st recorded mood diagnosis"
age_sex[which(rownames(age_sex) == "Age_mean"),]$comment <- "Age current"
age_sex[which(rownames(age_sex) == "Sex_F_pc"),]$comment <- ""

age_sex_dx$p_MTTD_less1m <- NA
age_sex_dx[which(rownames(age_sex_dx) == "Age_entry_mean"),]$p_MTTD_less1m <- age_entry_test_out_dx$p.value
age_sex_dx[which(rownames(age_sex_dx) == "Age_mean"),]$p_MTTD_less1m <- age_test_out_dx$p.value
age_sex_dx[which(rownames(age_sex_dx) == "Sex_F_pc"),]$p_MTTD_less1m <- sex_test_out_dx$p.value
age_sex_dx$comment <- ""
age_sex_dx[which(rownames(age_sex_dx) == "Age_entry_mean"),]$comment <- "Age of entry into EHR"
age_sex_dx[which(rownames(age_sex_dx) == "Age_mean"),]$comment <- "Age current"
age_sex_dx[which(rownames(age_sex_dx) == "Sex_F_pc"),]$comment <- ""

#------------------------------------------------------------
# Compare observation period
#------------------------------------------------------------

demog_cov <- getDf("demog_cov") 
observation_period <- demog_cov %>% dplyr::select(person_id, observation_years)

obs_period <- id_group_clean %>% 
  left_join(., observation_period) %>%
  dplyr::select(person_id, group, observation_years) %>%
  group_by(group) %>% 
  summarise(obs_years_median = paste0(round(median(observation_years), 2), " [", 
                                      round(quantile(observation_years, 0.25), 2), "-", 
                                      round(quantile(observation_years, 0.75), 2), "]")) %>%
  pivot_wider(names_from = "group", values_from = "obs_years_median") %>%
  mutate(var = "observationduration_years (median [IQR])")

obs_period_dx <- id_group_clean_dx %>% 
  left_join(., observation_period) %>%
  dplyr::select(person_id, group, observation_years) %>%
  group_by(group) %>% 
  summarise(obs_years_median = paste0(round(median(observation_years), 2), " [", 
                                      round(quantile(observation_years, 0.25), 2), "-", 
                                      round(quantile(observation_years, 0.75), 2), "]")) %>%
  pivot_wider(names_from = "group", values_from = "obs_years_median") %>%
  mutate(var = "observationduration_years (median [IQR])")

#------------------------------------------------------------
# Compare associations with state (PMMD group status, PMMD episode rate, PMMD definition, PManicD rate, % u30 at FEM, prescription availability)
# NO LONGER INCLUDING IN SUPP - DISREGARD
#------------------------------------------------------------

# munge meds data
meds_psych_bipolar_dx <- getDf("meds_psych_bipolar_dx")
bipolar_epall_meds_psych <- bipolar_epall %>%
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMD", "PMMD")) %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  mutate(dxstart_medstart = as.numeric(difftime(diagnosis_date, start_date, units = "days")),
        dxstart_medend = as.numeric(difftime(diagnosis_date, end_date, units = "days")),
        dxend_medstart = as.numeric(difftime(diagnosis_date_end, start_date, units = "days")),
        dxend_medend = as.numeric(difftime(diagnosis_date_end, end_date, units = "days"))) %>%
  filter((dxstart_medstart <= 0 & dxend_medend >= 0) |
        (dxstart_medstart >= 0 & dxstart_medend <= 0) |
        (dxend_medstart >= 0 & dxend_medend <= 0) |
        (dxstart_medstart >= 0 & dxend_medend <= 0))

# initial stats
bipolar_epall %>% group_by(person_id, group) %>% summarise(n=n()) %>% group_by(group) %>% summarise(n=n())
pmmd_ep_temporal %>% dplyr::select(person_id, diagnosis_date, mm_episode, dep_episode, PMMD, group) %>% distinct() %>% group_by(group, PMMD) %>% summarise(n=n())

bipolar_epall_tmp <- bipolar_epall %>% dplyr::select(person_id, diagnosis_date, mm_episode, dep_episode, PMMD, group) %>%
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMD", "PMMD"))
drugs <- bipolar_epall_meds_psych %>% dplyr::select(person_id, diagnosis_date, mm_episode, dep_episode, PMMD, group) %>% distinct()
tmp2 <- full_join(bipolar_epall_tmp, drugs) %>% left_join(., demog_cov) %>% mutate(set = "full")
tmp3 <- anti_join(bipolar_epall_tmp, drugs) %>% left_join(., demog_cov) %>% mutate(set = "no_drugs")
tmp4 <- inner_join(bipolar_epall_tmp, drugs) %>% left_join(., demog_cov) %>% mutate(set = "drugs")

# table of state representation with drug info 
state_by_drugs <- rbind(tmp2, tmp3, tmp4) %>%
  group_by(set, cov_state) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = set, values_from = n) %>%
  mutate(state_pc = full/nrow(tmp2), no_drugs_pc = no_drugs/full, drugs_pc = drugs/full) %>%
  arrange(-state_pc)

# table of state by PMMD episode rate
state_by_pmmdprop <- bipolar_epall_tmp %>% left_join(., demog_cov) %>%
  group_by(PMMD, cov_state) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = PMMD, values_from = n) %>%
  mutate(sum = PMMD+nonPMMD, state_pc = (PMMD+nonPMMD)/nrow(bipolar_epall), PMMD_pc = PMMD/(PMMD+nonPMMD)) %>%
  arrange(-sum)

# - combine state x drug info + state
state_by_drugs_pmmdprop <- left_join(state_by_drugs, state_by_pmmdprop)
state_by_drugs_pmmdprop_gg <- ggplot(state_by_drugs_pmmdprop, aes(x = state_pc, y = drugs_pc)) + 
  geom_point(aes(size = state_pc, colour = cov_state)) + geom_smooth(method = "lm") + 
  theme_bw() + theme(legend.position = "none") +
  xlab("Proportion of full sample in state") + ylab("Proportion of N per state with prescription data")
ggsave("state_by_drugs_pmmdprop.png", state_by_drugs_pmmdprop_gg)

state_by_drugs_pmmdprop2 <- rbind(tmp3, tmp4) %>%
  group_by(set, cov_state, PMMD) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = c("set", "PMMD"), values_from = n) %>%
  arrange(-drugs_PMMD) %>%
  mutate(PMMDprop_fordrugs = drugs_PMMD/(drugs_PMMD+drugs_nonPMMD), PMMDprop_forno_drugs = no_drugs_PMMD/(no_drugs_PMMD+no_drugs_nonPMMD)) %>%
  dplyr::select(cov_state, PMMDprop_fordrugs, PMMDprop_forno_drugs)

# table of state by PMMDgroup proportion
state_by_group <- bipolar_epall %>% dplyr::select(person_id, group) %>% distinct() %>% left_join(., demog_cov) %>%
  group_by(group, cov_state) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = group, values_from = n) %>%
  mutate(sum = PMMD+noPMMD, state_pc = (PMMD+noPMMD)/nrow(bipolar_epall), PMMD_pc = PMMD/(PMMD+noPMMD)) %>%
  arrange(-sum)
state_by_group_gg <- ggplot(state_by_group, aes(x = state_pc, y = PMMD_pc)) + 
  geom_point(aes(size = state_pc, colour = cov_state)) + geom_smooth(method = "lm") + 
  theme_bw() + theme(legend.position = "none") +
  xlab("Proportion of full sample in state") + ylab("Proportion of N per state in MDSgroup from broad definition")
ggsave("state_by_group.png", state_by_group_gg)

# table of state by PMMDgroup_dx definition proportion
pmmd_ep_dx_temporal <- getDf("pmmd_ep_dx_temporal")
nopmmd_epall_dx <- getDf("nopmmd_epall_dx")
bipolar_epall_dx <- pmmd_ep_dx_temporal %>% mutate(PMMD = ifelse(is.na(PMMD), "nonPMMD", "PMMD")) %>%
  full_join(., nopmmd_epall_dx %>% mutate(PMMD = "nonPMMD", group = "noPMMD"))
state_by_groupdx <- bipolar_epall_dx %>% dplyr::select(person_id, group) %>% distinct() %>% left_join(., demog_cov) %>%
  group_by(group, cov_state) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = group, values_from = n) %>%
  mutate(sum = PMMD+noPMMD, state_pc = (PMMD+noPMMD)/nrow(bipolar_epall), PMMDdx_pc = PMMD/(PMMD+noPMMD)) %>%
  arrange(-sum)

# plot of PMMDgroup vs PMMDgroup_dx (1 point per state)
compare_pmmdrate_definition_bystate <- state_by_group %>% dplyr::select(sum, cov_state, PMMD_pc) %>%
  left_join(., state_by_groupdx %>% dplyr::select(cov_state, PMMDdx_pc))
compare_pmmdrate_definition_bystate_gg <- ggplot(compare_pmmdrate_definition_bystate, aes(x = PMMDdx_pc, y = PMMD_pc)) + 
  geom_point(aes(size = sum, colour = cov_state)) + geom_smooth(method = "lm") + 
  theme_bw() + theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, color= "black") + xlim(c(0.2, 0.8)) + ylim(c(0.2, 0.8)) +
  xlab("Proportion of N per state in MDSgroup from ICD definition") + ylab("Proportion of N per state in MDSgroup from broad definition")
ggsave("stateprop_pmmdbroad_vs_pmmddx.png", compare_pmmdrate_definition_bystate_gg)

# table of state by PManicD proportion
state_by_pmmd_mmsplit <- bipolar_epall %>% left_join(mm_split) %>% filter(!is.na(mm_split)) %>% left_join(., demog_cov) %>%
  group_by(mm_split, cov_state) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = mm_split, values_from = n) %>%
  mutate(PMMD_pc = (mania+mixed)/nrow(pmmd_ep), PMMD_mania_pc = mania/(mania+mixed)) %>%
  arrange(-PMMD_pc)
state_by_pmmd_mmsplit_gg <- ggplot(state_by_pmmd_mmsplit, aes(x = PMMD_pc, y = PMMD_mania_pc)) + 
  geom_point(aes(size = PMMD_pc, colour = cov_state)) + geom_smooth(method = "lm") + 
  theme_bw() + theme(legend.position = "none") +
  xlab("Proportion of MDS episodes per state") + ylab("Proportion of ManicDS episodes per state")
ggsave("state_by_pmmd_mmsplit.png", state_by_pmmd_mmsplit_gg)

# table of state by number of people with first MM episode <30
state_by_mmage30 <- bipolar_epall %>% left_join(., demog_cov) %>%
  filter(state == "mm" & mm_episode == 1) %>%
  mutate(age_FEM = as.numeric(difftime(diagnosis_date, birth_year, units = "days"))/365) %>%
  mutate(age_FEM_u30 = ifelse(age_FEM <= 30, "u30", "o30")) %>%
  group_by(age_FEM_u30, cov_state) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = age_FEM_u30, values_from = n) %>%
  mutate(state_pc = (u30+o30)/nrow(pmmd_ep), PMMD_u30_pc = u30/(u30+o30)) %>%
  arrange(-state_pc)
state_by_mmage30_gg <- ggplot(state_by_mmage30, aes(x = state_pc, y = PMMD_u30_pc)) + 
  geom_point(aes(size = state_pc, colour = cov_state)) + geom_smooth(method = "lm") + 
  theme_bw() + theme(legend.position = "none") +
  xlab("Proportion of full sample per state") + ylab("Proportion of first recorded MM episodes at <30y per state")
ggsave("state_by_mmage30.png", state_by_mmage30_gg)

library(patchwork)

patchwork::wrap_plots(list(state_by_group_gg, compare_pmmdrate_definition_bystate_gg, state_by_pmmd_mmsplit_gg, state_by_mmage30_gg, state_by_drugs_pmmdprop_gg), ncol = 2) + 
  patchwork::plot_annotation(tag_levels = 'a')
ggsave("cov_state_eda_plots.png", height = 15, width = 10)

#------------------------------------------------------------
# Compare relationship status
#------------------------------------------------------------

rel_status <- demog_compare %>%
  mutate(cov_rel = ifelse(cov_rel %in% c("HUNK", "No matching concept"), "Unknown", cov_rel)) %>%
  #filter(group %in% c("MMonly", "PMMD", "noPMMD")) %>%
  group_by(group, cov_rel) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_rel, values_from = n) %>%
  column_to_rownames(., var = "group") %>%
  mutate(sum = rowSums(.)) %>%
  mutate(Rel_Married_pc = paste0(Married, " (", round(Married/sum*100, digits = 1), ")"), 
         Rel_Single_pc = paste0(Single, " (", round(Single/sum*100, digits = 1), ")"),
         Rel_Other_pc = paste0(Other, " (", round(Other/sum*100, digits = 1), ")"),
         Rel_Unknown_pc = paste0(Unknown, " (", round(Unknown/sum*100, digits = 1), ")")) %>%
  #mutate(across(where(is.numeric), ~ round(., digits = 1))) %>%
  dplyr::select(Rel_Married_pc, Rel_Single_pc, Rel_Other_pc, Rel_Unknown_pc) %>%
  #mutate(Rel_Married_Single_Unk_pc = paste(signif(Rel_Married_pc, 2), signif(Rel_Single_pc, 2), signif(Rel_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Rel_Married_Single_Unk_pc) %>%
  t() %>% data.frame 

rel_status_dx <- demog_compare %>%
  mutate(cov_rel = ifelse(cov_rel %in% c("HUNK", "No matching concept"), "Unknown", cov_rel)) %>%
  #filter(group_dx %in% c("MMonly", "PMMD", "noPMMD")) %>%
  group_by(group_dx, cov_rel) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_rel, values_from = n) %>%
  column_to_rownames(., var = "group_dx") %>%
  mutate(sum = rowSums(.)) %>%
  mutate(Rel_Married_pc = paste0(Married, " (", round(Married/sum*100, digits = 1), ")"), 
         Rel_Single_pc = paste0(Single, " (", round(Single/sum*100, digits = 1), ")"),
         Rel_Other_pc = paste0(Other, " (", round(Other/sum*100, digits = 1), ")"),
         Rel_Unknown_pc = paste0(Unknown, " (", round(Unknown/sum*100, digits = 1), ")")) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 1))) %>%
  dplyr::select(Rel_Married_pc, Rel_Single_pc, Rel_Other_pc, Rel_Unknown_pc) %>%
  #mutate(Rel_Married_Single_Unk_pc = paste(signif(Rel_Married_pc, 2), signif(Rel_Single_pc, 2), signif(Rel_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Rel_Married_Single_Unk_pc) %>%
  t() %>% data.frame 

# - choose to restrict statistical test to those with recorded data
rel_status_known <- demog_compare %>%
  group_by(group, cov_rel) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_rel, values_from = n) %>%
  column_to_rownames(., var = "group") %>%
  mutate(sum = Married + Single) %>%
  mutate(Rel_Married_pc = Married, 
         Rel_Single_pc = Single) %>%
  dplyr::select(c(Rel_Married_pc, Rel_Single_pc)) %>%
  t() %>% data.frame 
rel_status_test_out <- tidy(chisq.test(rel_status_known))

rel_status_known_dx <- demog_compare %>%
  group_by(group_dx, cov_rel) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_rel, values_from = n) %>%
  column_to_rownames(., var = "group_dx") %>%
  mutate(sum = Married + Single) %>%
  mutate(Rel_Married_pc = Married, 
         Rel_Single_pc = Single) %>%
  dplyr::select(c(Rel_Married_pc, Rel_Single_pc)) %>%
  t() %>% data.frame 
rel_status_test_out_dx <- tidy(chisq.test(rel_status_known_dx))

# - output
rel_status$p_MTTD_less1m <- ""
rel_status$p_MTTD_less1m[1] <- rel_status_test_out$p.value
rel_status$comment <- ""
rel_status$comment[1] <- "Married vs Single (exclude Unknown)"
rel_status_dx$p_MTTD_less1m <- ""
rel_status_dx$p_MTTD_less1m[1] <- rel_status_test_out_dx$p.value
rel_status_dx$comment <- ""
rel_status_dx$comment[1] <- "Married vs Single (exclude Unknown)"

#------------------------------------------------------------
# Compare employment status
#------------------------------------------------------------

employ_status <- demog_compare %>%
  group_by(group, cov_employ) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_employ, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Emp_Employed_any = (`Full-time employment` + `On work experience` + `Part-time employment` + `Self-employed`),
         Emp_Not_working = (`Not in labor force` + `Permanently unable to perform work activities due to medical condition` + `Retired` + `Unemployed`),
         Emp_Unknown = `No matching concept`) %>%
  mutate(Emp_Employed_any_pc = paste0(Emp_Employed_any, " (", round(Emp_Employed_any / sum * 100, digits = 1), ")"),
         Emp_Not_working_pc = paste0(Emp_Not_working, " (", round(Emp_Not_working / sum * 100, digits = 1), ")"),
         Emp_Unknown_pc = paste0(Emp_Unknown, " (", round(Emp_Unknown / sum * 100, digits = 1), ")")) %>%
  #mutate(across(where(is.numeric), ~ round(., digits = 1))) %>%
  dplyr::select(Emp_Employed_any_pc, Emp_Not_working_pc, Emp_Unknown_pc) %>%
  #mutate(Emp_Employed_Notworking_Unk_pc = paste(signif(Emp_Employed_any_pc, 2), signif(Emp_Not_working_pc, 2), signif(Emp_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Emp_Employed_Notworking_Unk_pc) %>%
  t() %>% data.frame 

employ_status_dx <- demog_compare %>%
  group_by(group_dx, cov_employ) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_employ, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group_dx") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Emp_Employed_any = (`Full-time employment` + `On work experience` + `Part-time employment` + `Self-employed`),
         Emp_Not_working = (`Not in labor force` + `Permanently unable to perform work activities due to medical condition` + `Retired` + `Unemployed`),
         Emp_Unknown = `No matching concept`) %>%
  mutate(Emp_Employed_any_pc = paste0(Emp_Employed_any, " (", round(Emp_Employed_any / sum * 100, digits = 1), ")"),
         Emp_Not_working_pc = paste0(Emp_Not_working, " (", round(Emp_Not_working / sum * 100, digits = 1), ")"),
         Emp_Unknown_pc = paste0(Emp_Unknown, " (", round(Emp_Unknown / sum * 100, digits = 1), ")")) %>%
  dplyr::select(Emp_Employed_any_pc, Emp_Not_working_pc, Emp_Unknown_pc) %>%
  #mutate(Emp_Employed_Notworking_Unk_pc = paste(signif(Emp_Employed_any_pc, 2), signif(Emp_Not_working_pc, 2), signif(Emp_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Emp_Employed_Notworking_Unk_pc) %>%
  t() %>% data.frame 

# - choose to restrict statistical test to those with recorded data
employ_status_known <- demog_compare %>%
  group_by(group, cov_employ) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_employ, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Emp_Employed_any_pc = (`Full-time employment` + `On work experience` + `Part-time employment` + `Self-employed`),
         Emp_Not_working_pc = (`Not in labor force` + `Permanently unable to perform work activities due to medical condition` + `Retired` + `Unemployed`)) %>%
  dplyr::select(c(Emp_Employed_any_pc, Emp_Not_working_pc)) %>%
  t() %>% data.frame 
employ_status_test_out <- tidy(chisq.test(employ_status_known))

employ_status_known_dx <- demog_compare %>%
  group_by(group_dx, cov_employ) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_employ, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group_dx") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Emp_Employed_any_pc = (`Full-time employment` + `On work experience` + `Part-time employment` + `Self-employed`),
         Emp_Not_working_pc = (`Not in labor force` + `Permanently unable to perform work activities due to medical condition` + `Retired` + `Unemployed`)) %>%
  dplyr::select(c(Emp_Employed_any_pc, Emp_Not_working_pc)) %>%
  t() %>% data.frame 
employ_status_test_out_dx <- tidy(chisq.test(employ_status_known_dx))

# - output
employ_status$p_MTTD_less1m <- ""
employ_status$p_MTTD_less1m[1] <- employ_status_test_out$p.value
employ_status$comment <- ""
employ_status$comment[1] <- "Employed vs Not Working (exclude Unknown)"

employ_status_dx$p_MTTD_less1m <- ""
employ_status_dx$p_MTTD_less1m[1] <- employ_status_test_out_dx$p.value
employ_status_dx$comment <- ""
employ_status_dx$comment[1] <- "Employed vs Not Working (exclude Unknown)"

#------------------------------------------------------------
# Compare education level
#------------------------------------------------------------

edu_level <- demog_compare %>%
  group_by(group, cov_edu) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_edu, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Edu_Highschoolupto = (`8th grade or less` + `Some high school, but did not graduate` + `High school diploma or GED`),
         Edu_Tertiaryany = (`Associate degree (e.g., AA, AS)` + `Bachelor's degree` + `Graduate degree` + `Some college`),
         Edu_Unknown = `No matching concept`) %>%
  mutate(Edu_Highschoolupto_pc = paste0(Edu_Highschoolupto, " (", round(Edu_Highschoolupto / sum * 100, digits = 1), ")"),
         Edu_Tertiaryany_pc = paste0(Edu_Tertiaryany, " (", round(Edu_Tertiaryany / sum * 100, digits = 1), ")"),
         Edu_Unknown_pc = paste0(Edu_Unknown, " (", round(Edu_Unknown / sum * 100, digits = 1), ")")) %>%
  #mutate(across(where(is.numeric), ~ round(., digits = 1))) %>%
  dplyr::select(Edu_Highschoolupto_pc, Edu_Tertiaryany_pc, Edu_Unknown_pc) %>%
  #mutate(Edu_Highschoolupto_Tertiaryany_Unk_pc = paste(signif(Edu_Highschoolupto_pc, 2), signif(Edu_Tertiaryany_pc, 2), signif(Edu_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Edu_Highschoolupto_Tertiaryany_Unk_pc) %>%
  t() %>% data.frame

edu_level_dx <- demog_compare %>%
  group_by(group_dx, cov_edu) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_edu, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group_dx") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Edu_Highschoolupto = (`8th grade or less` + `Some high school, but did not graduate` + `High school diploma or GED`),
         Edu_Tertiaryany = (`Associate degree (e.g., AA, AS)` + `Bachelor's degree` + `Graduate degree` + `Some college`),
         Edu_Unknown = `No matching concept`) %>%
  mutate(Edu_Highschoolupto_pc = paste0(Edu_Highschoolupto, " (", round(Edu_Highschoolupto / sum * 100, digits = 1), ")"),
         Edu_Tertiaryany_pc = paste0(Edu_Tertiaryany, " (", round(Edu_Tertiaryany / sum * 100, digits = 1), ")"),
         Edu_Unknown_pc = paste0(Edu_Unknown, " (", round(Edu_Unknown / sum * 100, digits = 1), ")")) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 1))) %>%
  dplyr::select(Edu_Highschoolupto_pc, Edu_Tertiaryany_pc, Edu_Unknown_pc) %>%
  #mutate(Edu_Highschoolupto_Tertiaryany_Unk_pc = paste(signif(Edu_Highschoolupto_pc, 2), signif(Edu_Tertiaryany_pc, 2), signif(Edu_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Edu_Highschoolupto_Tertiaryany_Unk_pc) %>%
  t() %>% data.frame 

edu_level_known <- demog_compare %>%
  left_join(., id_group_clean) %>%
  group_by(group, cov_edu) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_edu, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group") %>%
  dplyr::select(-c(`No matching concept`)) %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Edu_Highschoolupto_pc = (`8th grade or less` + `Some high school, but did not graduate` + `High school diploma or GED`),
         Edu_Tertiaryany_pc = (`Associate degree (e.g., AA, AS)` + `Bachelor's degree` + `Graduate degree` + `Some college`)) %>%
  dplyr::select(Edu_Highschoolupto_pc, Edu_Tertiaryany_pc) %>%
  t() %>% data.frame 
edu_level_test_out <- tidy(chisq.test(edu_level_known))

edu_level_known_dx <- demog_compare %>%
  group_by(group_dx, cov_edu) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_edu, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group_dx") %>%
  dplyr::select(-c(`No matching concept`)) %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Edu_Highschoolupto_pc = (`8th grade or less` + `Some high school, but did not graduate` + `High school diploma or GED`),
         Edu_Tertiaryany_pc = (`Associate degree (e.g., AA, AS)` + `Bachelor's degree` + `Graduate degree` + `Some college`)) %>%
  dplyr::select(Edu_Highschoolupto_pc, Edu_Tertiaryany_pc) %>%
  t() %>% data.frame 
edu_level_test_out_dx <- tidy(chisq.test(edu_level_known_dx))

# - output
edu_level$p_MTTD_less1m <- ""
edu_level$p_MTTD_less1m[1] <- edu_level_test_out$p.value
edu_level$comment <- ""
edu_level$comment[1] <- "High school vs Tertiary (exclude Unknown)"

edu_level_dx$p_MTTD_less1m <- ""
edu_level_dx$p_MTTD_less1m[1] <- edu_level_test_out_dx$p.value
edu_level_dx$comment <- ""
edu_level_dx$comment[1] <- "High school vs Tertiary (exclude Unknown)"

#------------------------------------------------------------
# Compare race
#------------------------------------------------------------

race <- demog_compare %>%
  group_by(group, cov_race) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_race, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Race_White_pc = paste0(White, " (", round(White / sum * 100, digits = 1), ")"),
         Race_BlackAA_pc = paste0(`Black or African American`, " (", round(`Black or African American` / sum * 100, digits = 1), ")"),
         Race_Asian_pc = paste0(Asian, " (", round(Asian / sum * 100, digits = 1), ")"),
         Race_AmIndian_Alaskan_pc = paste0(`American Indian or Alaska Native`, " (", round(`American Indian or Alaska Native` / sum * 100, digits = 1), ")"),
         Race_PacificIslander_pc = paste0(`Native Hawaiian or Other Pacific Islander`, " (", round(`Native Hawaiian or Other Pacific Islander` / sum * 100, digits = 1), ")"),
         Race_Other_pc = paste0(`Other Race`, " (", round(`Other Race` / sum * 100, digits = 1), ")"),
         Race_Unknown_pc = paste0(`No matching concept`, " (", round(`No matching concept` / sum * 100, digits = 1), ")")) %>%
  #mutate(across(where(is.numeric), ~ round(., digits = 1))) %>%
  dplyr::select(Race_White_pc, Race_BlackAA_pc, Race_AmIndian_Alaskan_pc, Race_PacificIslander_pc, Race_Other_pc, Race_Unknown_pc) %>%
  #mutate(Race_White_BlackAA_Asian_AmIndAlaskan_Pacific_Other_Unk_pc = paste(signif(Race_White_pc, 2), signif(Race_BlackAA_pc, 2), signif(Race_Asian_pc, 2), signif(Race_AmIndian_Alaskan_pc, 2), signif(Race_PacificIslander_pc, 2), signif(Race_Other_pc, 2), signif(Race_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Race_White_BlackAA_Asian_AmIndAlaskan_Pacific_Other_Unk_pc) %>%
  t() %>% data.frame 

race_dx <- demog_compare %>%
  group_by(group_dx, cov_race) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = cov_race, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group_dx") %>%
  mutate(sum = rowSums(., na.rm = T)) %>%
  mutate(Race_White_pc = paste0(White, " (", round(White / sum * 100, digits = 1), ")"),
         Race_BlackAA_pc = paste0(`Black or African American`, " (", round(`Black or African American` / sum * 100, digits = 1), ")"),
         Race_Asian_pc = paste0(Asian, " (", round(Asian / sum * 100, digits = 1), ")"),
         Race_AmIndian_Alaskan_pc = paste0(`American Indian or Alaska Native`, " (", round(`American Indian or Alaska Native` / sum * 100, digits = 1), ")"),
         Race_PacificIslander_pc = paste0(`Native Hawaiian or Other Pacific Islander`, " (", round(`Native Hawaiian or Other Pacific Islander` / sum * 100, digits = 1), ")"),
         Race_Other_pc = paste0(`Other Race`, " (", round(`Other Race` / sum * 100, digits = 1), ")"),
         Race_Unknown_pc = paste0(`No matching concept`, " (", round(`No matching concept` / sum * 100, digits = 1), ")")) %>%
  #mutate(across(where(is.numeric), ~ round(., digits = 1))) %>%
  dplyr::select(Race_White_pc, Race_BlackAA_pc, Race_AmIndian_Alaskan_pc, Race_PacificIslander_pc, Race_Other_pc, Race_Unknown_pc) %>%
  #mutate(Race_White_BlackAA_Asian_AmIndAlaskan_Pacific_Other_Unk_pc = paste(signif(Race_White_pc, 2), signif(Race_BlackAA_pc, 2), signif(Race_Asian_pc, 2), signif(Race_AmIndian_Alaskan_pc, 2), signif(Race_PacificIslander_pc, 2), signif(Race_Other_pc, 2), signif(Race_Unknown_pc, 2), sep = "/")) %>%
  #dplyr::select(Race_White_BlackAA_Asian_AmIndAlaskan_Pacific_Other_Unk_pc) %>%
  t() %>% data.frame 

race_known <- demog_compare %>% 
  group_by(group, cov_race) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "cov_race", values_from = "n") %>% 
  replace(is.na(.), 0) %>% 
  mutate(Race_White = White,
         Race_BlackAA = `Black or African American`,
         Race_Asian = `Asian`,
         Race_AmericanIndian_Alaskan = `American Indian or Alaska Native`,
         Race_PacificIslander = `Native Hawaiian or Other Pacific Islander`) %>%
  dplyr::select(group, Race_White, Race_BlackAA, Race_Asian, Race_AmericanIndian_Alaskan, Race_PacificIslander) %>%
  column_to_rownames(., var = "group") %>% as.matrix
race_test_out <- data.frame(variable = "Race", p.value = tidy(chisq.test(race_known)) %>% pull(p.value) %>% as.vector)

race_known_dx <- demog_compare %>% 
  group_by(group_dx, cov_race) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "cov_race", values_from = "n") %>% 
  replace(is.na(.), 0) %>% 
  mutate(Race_White = White,
         Race_BlackAA = `Black or African American`,
         Race_Asian = `Asian`,
         Race_AmericanIndian_Alaskan = `American Indian or Alaska Native`,
         Race_PacificIslander = `Native Hawaiian or Other Pacific Islander`) %>%
  dplyr::select(group_dx, Race_White, Race_BlackAA, Race_Asian, Race_AmericanIndian_Alaskan, Race_PacificIslander) %>%
  column_to_rownames(., var = "group_dx") %>% as.matrix
race_test_out_dx <- data.frame(variable = "Race", p.value = tidy(chisq.test(race_known_dx)) %>% pull(p.value) %>% as.vector)

# - output
race$p_MTTD_less1m <- ""
race$p_MTTD_less1m[1] <- race_test_out$p.value
race$comment <- ""
race$comment[1] <- "Exclude Other and Unknown"

race_dx$p_MTTD_less1m <- ""
race_dx$p_MTTD_less1m[1] <- race_test_out_dx$p.value
race_dx$comment <- ""
race_dx$comment[1] <- "Exclude Other and Unknown"

#------------------------------------------------------------
# Compare state
#------------------------------------------------------------

state <- demog_compare %>%
  group_by(group, cov_state) %>%
  summarise(n=n()) %>% 
  mutate(cov_state = paste("state", cov_state, sep = "_")) %>%
  pivot_wider(names_from = cov_state, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group")
tmpsum <- colSums(state)
tmpsum_keep <- names(tmpsum[which(tmpsum > 100)]) # approximate 1%
tmpsum_agg <- names(tmpsum[which(tmpsum <= 100)])
state2 <- state %>%
  mutate(state_Other = rowSums(across(all_of(tmpsum_agg)), na.rm = TRUE)) %>%
  dplyr::select(all_of(tmpsum_keep), state_Other) %>%
  mutate(total = rowSums(across(everything()), na.rm = TRUE),
         across(-total, ~ paste0(., " (", round(. / total * 100, 1), ")"))) %>%
  dplyr::select(-total) %>% t() %>% data.frame()

state_dx <- demog_compare %>%
  group_by(group, cov_state) %>%
  summarise(n=n()) %>% 
  mutate(cov_state = paste("state", cov_state, sep = "_")) %>%
  pivot_wider(names_from = cov_state, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var = "group")
tmpsum <- colSums(state_dx)
tmpsum_keep <- names(tmpsum[which(tmpsum > 100)]) # approximate 1%
tmpsum_agg <- names(tmpsum[which(tmpsum <= 100)])
state2_dx <- state_dx %>%
  mutate(state_Other = rowSums(across(all_of(tmpsum_agg)), na.rm = TRUE)) %>%
  dplyr::select(all_of(tmpsum_keep), state_Other) %>%
  mutate(total = rowSums(across(everything()), na.rm = TRUE),
         across(-total, ~ paste0(., " (", round(. / total * 100, 1), ")"))) %>%
  dplyr::select(-total) %>% t() %>% data.frame()

state_test0 <- demog_compare %>% 
  group_by(group, cov_state) %>% 
  summarise(n = n()) %>% 
  mutate(cov_state = paste("cov_state", cov_state, sep = "_")) %>%
  pivot_wider(names_from = "cov_state", values_from = "n") %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames(., var = "group")

tmpsum2 <- colSums(state_test0)
tmpsum2_keep <- names(tmpsum2[which(tmpsum2 > 100)]) # approximate 1%
tmpsum2_agg <- names(tmpsum2[which(tmpsum2 <= 100)])

state_test0_dx <- demog_compare %>% 
  group_by(group_dx, cov_state) %>% 
  summarise(n = n()) %>% 
  mutate(cov_state = paste("cov_state", cov_state, sep = "_")) %>%
  pivot_wider(names_from = "cov_state", values_from = "n") %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames(., var = "group_dx")

tmpsum2_dx <- colSums(state_test0_dx)
tmpsum2_keep_dx <- names(tmpsum2[which(tmpsum2_dx > 100)]) # approximate 1%
tmpsum2_agg_dx <- names(tmpsum2[which(tmpsum2_dx <= 100)])

state_test <- state_test0 %>%
  mutate(state_Other = rowSums(across(all_of(tmpsum2_agg)))) %>%
  dplyr::select(all_of(tmpsum2_keep), state_Other) %>% as.matrix
state_test_out <- data.frame(variable = "State", p.value = tidy(chisq.test(state_test)) %>% pull(p.value) %>% as.vector)

state_test_dx <- state_test0_dx %>%
  mutate(state_Other = rowSums(across(all_of(tmpsum2_agg_dx)))) %>%
  dplyr::select(all_of(tmpsum2_keep_dx), state_Other) %>% as.matrix
state_test_out_dx <- data.frame(variable = "State", p.value = tidy(chisq.test(state_test_dx)) %>% pull(p.value) %>% as.vector)

# - output
state2$p_MTTD_less1m <- ""
state2$p_MTTD_less1m[1] <- state_test_out$p.value
state2$comment <- ""
state2$comment[1] <- "Include only states with >= 100 total counts"

state2_dx$p_MTTD_less1m <- ""
state2_dx$p_MTTD_less1m[1] <- state_test_out_dx$p.value
state2_dx$comment <- ""
state2_dx$comment[1] <- "Include only states with >= 100 total counts"

#------------------------------------------------------------
# Compare mixed/manic episodes
#------------------------------------------------------------

mm_split <- getDf("mm_split")

# Determine how many people have had a mixed episode at some point in EHR
mm_split %>% 
  group_by(person_id, mm_split) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = "mm_split", values_from = "n") %>% 
  mutate(mixed_ever = ifelse(!is.na(mixed), "yes", "no")) %>%
  group_by(mixed_ever) %>% summarise(n=n())
# Determine how many episodes were mixed episodes
mm_split %>% group_by(mm_split) %>% summarise(n=n())

# Add to table
mm_split_table <- mm_split %>% 
  inner_join(., id_group_clean) %>%
  group_by(person_id, group, mm_split) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = "mm_split", values_from = "n") %>% 
  mutate(mixed_ever = ifelse(!is.na(mixed), "mixed_ever", "no_mixed_ever")) %>%
  group_by(group, mixed_ever) %>% summarise(n=n()) %>%
  pivot_wider(names_from = "mixed_ever", values_from = "n")

mm_split_test_out <- tidy(chisq.test(mm_split_table[,2:3]))

mm_split_table2 <- mm_split_table %>%
  mutate(mixed_ever_pc = mixed_ever/(mixed_ever+no_mixed_ever)*100) %>%
  pivot_longer(cols = c("mixed_ever_pc"), names_to = "demographic", values_to = "mixed_ever_pc") %>%
  dplyr::select(-c("no_mixed_ever", "mixed_ever")) %>%
  pivot_wider(names_from = "group", values_from = "mixed_ever_pc") %>%
  mutate(p_MTTD_less1m = as.character(mm_split_test_out$p.value)) %>%
  mutate(comment = "Chi-squared test on count") %>%
  mutate(var = "mixedever_pc", demographic = "Mania_v_Mixed") %>%
  mutate(MTTD_less1m = as.character(signif(MTTD_less1m,2)), MTTD_more1m = as.character(signif(MTTD_more1m,2)))

# Diagnosis definition
mm_split_table_dx <- mm_split %>% 
  inner_join(., id_group_clean_dx) %>%
  group_by(person_id, group, mm_split) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = "mm_split", values_from = "n") %>% 
  mutate(mixed_ever = ifelse(!is.na(mixed), "mixed_ever", "no_mixed_ever")) %>%
  group_by(group, mixed_ever) %>% summarise(n=n()) %>%
  pivot_wider(names_from = "mixed_ever", values_from = "n")

mm_split_test_out_dx <- tidy(chisq.test(mm_split_table_dx[,2:3]))

mm_split_table2_dx <- mm_split_table_dx %>%
  mutate(mixed_ever_pc = mixed_ever/(mixed_ever+no_mixed_ever)*100) %>%
  pivot_longer(cols = c("mixed_ever_pc"), names_to = "demographic", values_to = "mixed_ever_pc") %>%
  dplyr::select(-c("no_mixed_ever", "mixed_ever")) %>%
  pivot_wider(names_from = "group", values_from = "mixed_ever_pc") %>%
  mutate(p_MTTD_less1m = as.character(mm_split_test_out_dx$p.value)) %>%
  mutate(comment = "Chi-squared test on count") %>%
  mutate(var = "mixedever_pc", demographic = "Mania_v_Mixed") %>%
  mutate(MTTD_less1m = as.character(signif(MTTD_less1m,2)), MTTD_more1m = as.character(signif(MTTD_more1m,2)))

#------------------------------------------------------------
# Compare number of episodes
#------------------------------------------------------------

# - identify who is in short switch group
# - look at number of episodes overall
# - look at number of episodes per person

# - exclude same day mania/mixed episodes from these counts
mm0_mmenddepstart0_dep0_toexclude <- getDf("mm0_mmenddepstart0_dep0_toexclude")
mm0_mmenddepstart0_dep0_toexclude_date <- mm0_mmenddepstart0_dep0_toexclude %>%
  ungroup %>%
  dplyr::select(person_id, diagnosis_date_mm, diagnosis_date_end_mm) %>%
  rename(diagnosis_date = diagnosis_date_mm, diagnosis_date_end = diagnosis_date_end_mm)

bipolar_epall_excludemmdep0 <- bipolar_epall %>% 
  anti_join(., mm0_mmenddepstart0_dep0_toexclude_date)

# episode count
# all of the "NA" group (representing IDs missing from id_group_clean but present in bipolar_epall) 
# represent those with unresolvable mm/dep sequence (which were excluded from the file used to derive id_group_clean)
# -> do not include as unable to tell which group they truly belong to
episodes_number <- bipolar_epall_excludemmdep0 %>%
  dplyr::select(-group) %>%
  left_join(id_group_clean) %>%
  filter(!is.na(group)) %>%
  left_join(mm_split) %>%
  mutate(state_mmdep = ifelse(state == "mm", mm_split, state)) %>%
  group_by(group, state_mmdep) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = "group", values_from = "n") %>%
  mutate(state_mmdep = paste0("episodecount_", state_mmdep)) %>%
  rename(var = state_mmdep)

# per-person episodes
episodes_number_pp <- bipolar_epall_excludemmdep0 %>%
  dplyr::select(-group) %>%
  left_join(id_group_clean) %>%
  filter(!is.na(group)) %>%
  left_join(mm_split) %>%
  mutate(state_mmdep = ifelse(state == "mm", mm_split, state)) %>%
  group_by(person_id, group, state_mmdep) %>%
  summarise(n=n()) %>%
  group_by(group, state_mmdep) %>%
  summarise(median=paste0(median(n), "[", quantile(n, 0.25), "-", q3 = quantile(n, 0.75), "]")) %>%
  pivot_wider(names_from = "group", values_from = "median") %>%
  mutate(state_mmdep = paste0("episodecount_", state_mmdep)) %>%
  rename(var = state_mmdep)

#------------------------------------------------------------
# Compare number of visits
#------------------------------------------------------------

bipolar_visit <- getDf("bipolar_visit")

visits_number <- bipolar_visit %>%
  group_by(person_id, visit_type) %>% summarise(n=n()) %>%
  left_join(., id_group_clean) %>%
  filter(!is.na(group)) %>%
  dplyr::select(person_id, n, group, visit_type) %>%
  group_by(group, visit_type) %>% summarise(visit_count=sum(n)) %>%
  pivot_wider(names_from = "group", values_from = "visit_count") %>%
  mutate(var = paste0("visitcount_", visit_type)) %>%
  dplyr::select(-visit_type)

visits_number_pp <- bipolar_visit %>%
  group_by(person_id, visit_type) %>% summarise(n=n()) %>%
  left_join(., id_group_clean) %>%
  filter(!is.na(group)) %>%
  dplyr::select(person_id, n, group, visit_type) %>%
  group_by(person_id, group, visit_type) %>% 
  summarise(n=sum(n)) %>%
  group_by(group, visit_type) %>% 
  summarise(median=paste0(median(n), "[", quantile(n, 0.25), "-", q3 = quantile(n, 0.75), "]")) %>%
  pivot_wider(names_from = "group", values_from = "median") %>%
  mutate(var = paste0("visitcount_", visit_type, "_pp")) %>%
  dplyr::select(-visit_type)

#------------------------------------------------------------
# Aggregate
#------------------------------------------------------------

demog_compare_table <- full_join(n_patients %>% rownames_to_column(., "var"), 
    rbind(age_sex, rel_status, employ_status, edu_level, race, state2) %>% rownames_to_column(., "var")) %>%
  full_join(., obs_period) %>%
  full_join(., mm_split_table2) %>%
  full_join(., episodes_number %>% mutate(MTTD_less1m = as.character(MTTD_less1m), MTTD_more1m = as.character(MTTD_more1m))) %>%
  full_join(., episodes_number_pp %>% mutate(MTTD_less1m = as.character(MTTD_less1m), MTTD_more1m = as.character(MTTD_more1m))) %>%
  #full_join(., episodes_state_number_pp %>% mutate(MTTD_less1m = as.character(MTTD_less1m), MTTD_more1m = as.character(MTTD_more1m))) %>%
  full_join(., visits_number %>% mutate(MTTD_less1m = as.character(MTTD_less1m), MTTD_more1m = as.character(MTTD_more1m))) %>%
  full_join(., visits_number_pp %>% mutate(MTTD_less1m = as.character(MTTD_less1m), MTTD_more1m = as.character(MTTD_more1m))) %>%
    mutate(demographic = unlist(lapply(strsplit(var, "_"), function(x) x[1]))) %>%
  relocate(demographic, .before = var) %>%
  #dplyr::select(-"NA") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, ""))) %>%
  mutate(across(where(is.character), ~replace_na(.x, ""))) 

demog_compare_table_dx <- full_join(n_patients %>% rownames_to_column(., "var"), 
    rbind(age_sex_dx, rel_status_dx, employ_status_dx, edu_level_dx, race_dx, state2_dx) %>% rownames_to_column(., "var")) %>%
  full_join(., obs_period_dx) %>%
  full_join(., mm_split_table2_dx) %>%
  mutate(demographic = unlist(lapply(strsplit(var, "_"), function(x) x[1]))) %>%
  relocate(demographic, .before = var) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, ""))) %>%
  mutate(across(where(is.character), ~replace_na(.x, ""))) 

saveDf(demog_compare_table[1:50,], "demog_compare_table_part1", overwrite = T)
saveDf(demog_compare_table[51:nrow(demog_compare_table),], "demog_compare_table_part2", overwrite = T)

saveDf(demog_compare_table_dx[1:50,], "demog_compare_table_dx_part1", overwrite = T)
saveDf(demog_compare_table_dx[51:nrow(demog_compare_table_dx),], "demog_compare_table_dx_part2", overwrite = T)

demog_compare_table_Table1 <- demog_compare_table %>%
  #filter(demographic != "state") %>%
  mutate(var = ifelse(demographic == "state", paste0(var, "_pc"), var)) %>%
  filter(demographic != "mixedever") %>%
  #dplyr::select(-"MMonly") %>%
  mutate(var = sub(".*?_","", var)) %>%
  mutate(demographic = gsub("Rel", "Relationship status", demographic)) %>%
  mutate(demographic = ifelse(demographic == "n" & var == "n", "N", demographic)) %>%
  mutate(var = ifelse(demographic == "N" & var == "n", "individuals", var)) %>%
  mutate(var = ifelse(demographic == "Age" & var == "entry_mean", "entry mean (SD)", var)) %>%
  mutate(var = ifelse(demographic == "Age" & var == "mean", "current mean (SD)", var)) %>%
  mutate(var = gsub("_pc", "%", var)) %>%
  mutate(var = gsub("pc", "%", var)) %>%
  mutate(var = gsub("_pp", " per-person (median [IQR])", var)) %>%
  #mutate(var = gsub("nonMDS", "TTD>6m", var)) %>%
  #mutate(var = gsub("MDS", "TTD<6m", var)) %>%
  mutate(var = gsub("_", " ", var)) %>%
  mutate(demographic = ifelse(demographic == "Rel", "Relationship status", demographic)) %>%
  mutate(demographic = ifelse(demographic == "Emp", "Employment status", demographic)) %>%
  mutate(demographic = ifelse(demographic == "Edu", "Educational level (highest)", demographic)) %>%
  mutate(demographic = ifelse(demographic == "state", "State", demographic)) %>%
  mutate(demographic = ifelse(demographic == "observationduration", "Observation time", demographic)) %>%
  #mutate(demographic = ifelse(demographic == "mixedever", "Mixed episodes ever", demographic)) %>%
  mutate(demographic = ifelse(demographic == "episodecount", "Count of episodes", demographic)) %>%
  mutate(demographic = ifelse(demographic == "visitcount", "Count of visits", demographic)) %>%
  mutate(var = case_when(var == "dep" ~ "Depression",
                         var == "mania" ~ "Mania",
                         var == "mixed" ~ "Mixed",
                          .default = var))
colnames(demog_compare_table_Table1) <- c("Demographic","Detail","Lifetime MDTT<1m","Never MDTT>1m", "p-value", "comment")
saveDf(demog_compare_table_Table1[1:50,], "demog_compare_table_Table1_part1", overwrite = T)
saveDf(demog_compare_table_Table1[51:nrow(demog_compare_table_Table1),], "demog_compare_table_Table1_part2", overwrite = T)

#============================================================
# 
# Check MMonly group
#
#============================================================

id_group <- getDf("id_group")
all_mania_dep <- getDf("all_mania_dep")
demog_cov <- getDf("demog_cov")
bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")

mmonly_id <- id_group %>% filter(group == "MMonly")
bipolar_id <- bipolar_epall_count_mm_to_depafter_inf %>% dplyr::select(person_id) %>% distinct

all_mania_dep_mmonly <- all_mania_dep %>%
  filter(person_id %in% mmonly_id$person_id)
all_mania_dep_bipolar <- all_mania_dep %>%
  filter(person_id %in% bipolar_id$person_id)

# How many mm ICD diagnoses per person?
all_mania_dep_mmonly_ppcount <- all_mania_dep_mmonly %>%
  group_by(person_id) %>%
  summarise(n=n()) %>%
  mutate(group = "MMonly")
all_mania_dep_bipolar_ppcount <- all_mania_dep_bipolar %>% # compare to mm subset to be fair
  filter(state == "mm") %>%
  group_by(person_id) %>%
  summarise(n=n()) %>%
  mutate(group = "MM_Dep")
mmicdcodes_ecdf <- full_join(all_mania_dep_mmonly_ppcount, all_mania_dep_bipolar_ppcount)

all_mania_dep_mmonly_ppcount %>% group_by(group, n) %>% summarise(m=n()) %>% mutate(prop = m/sum(m))
all_mania_dep_bipolar_ppcount %>% group_by(group, n) %>% summarise(m=n()) %>% mutate(prop = m/sum(m))

mmicdcodes_ecdf.plot <- ggplot(mmicdcodes_ecdf, aes(x = n, colour = group, fill = group)) +
  stat_ecdf() +
  theme_bw() +
  #theme(legend.position = "none") +
  xlab("Number of ICD diagnosis (mania/mixed)") +
  scale_fill_manual(values = millerstone[c(2,3,9)]) +
  scale_colour_manual(values = millerstone[c(2,3,9)]) +
  coord_cartesian(xlim = c(1, 50))
ggsave("mmicdcodes_ecdf.png", width = 7.5, height = 3.5)

# Do most of these people only have 1 (mm) entry?
#mmicdcodes_n <- mmicdcodes_hist %>%
#  group_by(group) %>%
#  mutate(n_indiv = n()) %>%
#  group_by(group, n_indiv, n) %>%
#  summarise(m=n()) %>%
#  mutate(prop_with_n = m/n_indiv)
#mmicdcodes_n %>% filter(n == 1)

# What are the specific codes and how do they differ between groups?
mmonly_codes_bar.tmp <- all_mania_dep_mmonly %>% 
  group_by(diagnosis, icd_code) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(group = "MMonly")
bipolar_codes_bar.tmp <- all_mania_dep_bipolar %>%
  filter(state == "mm") %>%
  group_by(diagnosis, icd_code) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(group = "MM_Dep")
mmicdcodes_dx_barchart <- full_join(mmonly_codes_bar.tmp, bipolar_codes_bar.tmp) %>%
  mutate(icd_version = case_when(substring(icd_code, 1,1) == "2" ~ "icd9",
                             substring(icd_code, 1,1) == "F" ~ "icd10",
                            .default = "other")) #%>%
#  mutate(diagnosis = str_wrap(diagnosis, width = 50))
#  mutate(recency = case_when(grepl("remission", diagnosis) ~ "remission",
#                            grepl("recent", diagnosis) ~ "recent",
#                             grepl("current", diagnosis) ~ "current",
#                            .default = "other"))

mmicdcodes_dx_barchart.plot <- ggplot(mmicdcodes_dx_barchart, aes(x = diagnosis, y = prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_bw() +
  #theme(axis.text.x = element_text(size = 24)) +
  xlab("ICD diagnosis (mania/mixed)") +
  ylab("Proportion of group") +
  scale_fill_manual(values = millerstone[c(2,3,9)]) +
  facet_wrap(~ icd_version, nrow = 1)
ggsave("mmicdcodes_dx_barchart.png", width = 15, height = 15)

# What is follow-up time like for these people?
observation_period <- getDf("observation_period")

observation_period_mmonly <- observation_period %>%
  filter(person_id %in% mmonly_id$person_id) %>%
   mutate(observation_period_years = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days"))/365)
observation_period_bipolar <- observation_period %>%
  filter(person_id %in% bipolar_id$person_id) %>%
  mutate(observation_period_years = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days"))/365)

observation_period_mmonly_ppcount <- observation_period_mmonly %>%
  group_by(person_id) %>%
  summarise(obs_period_years_sum=sum(observation_period_years)) %>%
  mutate(group = "MMonly")
observation_period_bipolar_ppcount <- observation_period_bipolar %>%
  group_by(person_id) %>%
  summarise(obs_period_years_sum=sum(observation_period_years)) %>%
  mutate(group = "MM_Dep")

observation_period_mmonly_bipolar_ppcount <- full_join(observation_period_mmonly_ppcount, observation_period_bipolar_ppcount)
observation_period_mmonly_bipolar_ppcount.plot <- ggplot(observation_period_mmonly_bipolar_ppcount, aes(x = obs_period_years_sum, colour = group)) +
  stat_ecdf() +
  theme_bw() +
  xlab("Observation years per individual") +
  scale_fill_manual(values = millerstone[c(2,3,9)]) +
  scale_colour_manual(values = millerstone[c(2,3,9)])
ggsave("observation_period_mmonly_bipolar_ppcount_ecdf.png", width = 7.5, height = 3.5)

# EHR entry date
ehr_start_date_mmonly_bipolar <- full_join(observation_period_mmonly %>% dplyr::select(person_id, observation_period_start_date) %>% mutate(group = "MMonly"),
  observation_period_bipolar %>% dplyr::select(person_id, observation_period_start_date) %>% mutate(group = "MM_Dep"))

ehr_start_date_mmonly_bipolar.plot <- ggplot(ehr_start_date_mmonly_bipolar, aes(x = observation_period_start_date, colour = group)) +
  stat_ecdf() +
  theme_bw() +
  xlab("Observation start date") +
  scale_fill_manual(values = millerstone[c(2,3,9)]) +
  scale_colour_manual(values = millerstone[c(2,3,9)])
ggsave("ehr_start_date_mmonly_bipolar_ecdf.png", width = 7.5, height = 3.5)

# Relationship of group to state
demog_cov <- getDf("demog_cov")
demog_cov_state <- demog_cov %>% dplyr::select("person_id", "cov_state")

demog_cov_state_mmonly <- demog_cov_state %>%
  filter(person_id %in% mmonly_id$person_id) %>%
  mutate(group = "MMonly")
demog_cov_state_bipolar <- demog_cov_state %>%
  filter(person_id %in% bipolar_id$person_id) %>%
  mutate(group = "MM_Dep")

demog_cov_state_mmonly_bipolar <- full_join(demog_cov_state_mmonly, demog_cov_state_bipolar) %>%
  group_by(group, cov_state) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)

demog_cov_state_mmonly_bipolar_barchart.plot <- ggplot(demog_cov_state_mmonly_bipolar, aes(x = fct_reorder(cov_state, prop), y = prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_bw() +
  xlab("State (of residence)") +
  ylab("Proportion of group") +
  scale_fill_manual(values = millerstone[c(2,3,9)])
ggsave("demog_cov_state_mmonly_bipolar_barchart.png", width = 15, height = 7)

# Check other systematic differences
mmonly_demog_cov <- inner_join(mmonly_id, demog_cov) 
bipolar_demog_cov <- inner_join(bipolar_id, demog_cov) %>% mutate(group = "bipolar")

# - race
mmonly_cov_race <- mmonly_demog_cov %>%
  group_by(group, cov_race) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
bipolar_cov_race <- bipolar_demog_cov %>%
  group_by(group, cov_race) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
all_cov_race <- full_join(mmonly_cov_race, bipolar_cov_race) %>%
  mutate(cov_race = ifelse(cov_race == "Native Hawaiian or Other Pacific Islander", "Native Hawaiian/Pacific Islander", cov_race)) %>%
  mutate(cov_race = ifelse(cov_race == "No matching concept", "N/A", cov_race))

demog_cov_race_mmonly_bipolar_barchart.plot <- ggplot(all_cov_race, aes(x = fct_reorder(cov_race, prop), y = prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_bw() +
  xlab("Race") +
  ylab("Proportion") +
  scale_fill_manual(values = millerstone[c(2,3,9)])

# - cov_rel
mmonly_cov_rel <- mmonly_demog_cov %>%
  group_by(group, cov_rel) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
bipolar_cov_rel <- bipolar_demog_cov %>%
  group_by(group, cov_rel) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
all_cov_rel <- full_join(mmonly_cov_rel, bipolar_cov_rel) %>%
  mutate(cov_rel = ifelse(cov_rel == "No matching concept", "N/A", cov_rel))

demog_cov_rel_mmonly_bipolar_barchart.plot <- ggplot(all_cov_rel, aes(x = fct_reorder(cov_rel, prop), y = prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_bw() +
  xlab("Relationship status") +
  ylab("Proportion") +
  scale_fill_manual(values = millerstone[c(2,3,9)])

# - cov_employ
mmonly_cov_employ <- mmonly_demog_cov %>%
  group_by(group, cov_employ) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
bipolar_cov_employ <- bipolar_demog_cov %>%
  group_by(group, cov_employ) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
all_cov_employ <- full_join(mmonly_cov_employ, bipolar_cov_employ) %>%
  mutate(cov_employ = ifelse(cov_employ == "Permanently unable to perform work activities due to medical condition", "Permanently unable to work", cov_employ)) %>%
  mutate(cov_employ = ifelse(cov_employ == "No matching concept", "N/A", cov_employ))

demog_cov_employ_mmonly_bipolar_barchart.plot <- ggplot(all_cov_employ, aes(x = fct_reorder(cov_employ, prop), y = prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_bw() +
  xlab("Employment status") +
  ylab("Proportion") +
  scale_fill_manual(values = millerstone[c(2,3,9)])

# - cov_edu
mmonly_cov_edu <- mmonly_demog_cov %>%
  group_by(group, cov_edu) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
bipolar_cov_edu <- bipolar_demog_cov %>%
  group_by(group, cov_edu) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop = n/total_n)
all_cov_edu <- full_join(mmonly_cov_edu, bipolar_cov_edu) %>%
  mutate(cov_edu = ifelse(cov_edu == "Some high school, but did not graduate", "High school, did not grad", cov_edu)) %>%
  mutate(cov_edu = ifelse(cov_edu == "Associate degree (e.g., AA, AS)", "Associate degree", cov_edu)) %>%
  mutate(cov_edu = ifelse(cov_edu == "No matching concept", "N/A", cov_edu))

demog_cov_edu_mmonly_bipolar_barchart.plot <- ggplot(all_cov_edu, aes(x = fct_reorder(cov_edu, prop), y = prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_bw() +
  xlab("Education") +
  ylab("Proportion") +
  scale_fill_manual(values = millerstone[c(2,3,9)])

# Combine into 1 plot
library(patchwork)

p_a <- mmicdcodes_ecdf.plot + theme(legend.position = "none")
p_b <- observation_period_mmonly_bipolar_ppcount.plot + theme(legend.position = "none")

p_c <- mmicdcodes_dx_barchart.plot +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

p_d <- ehr_start_date_mmonly_bipolar.plot + theme(legend.position = "none")

p_e <- demog_cov_state_mmonly_bipolar_barchart.plot +
  theme(
    axis.text.y = element_text(size = 5),
    legend.position = "none"
  )

p_f <- demog_cov_race_mmonly_bipolar_barchart.plot  + theme(legend.position = "none")
p_g <- demog_cov_rel_mmonly_bipolar_barchart.plot  + theme(legend.position = "none")
p_h <- demog_cov_employ_mmonly_bipolar_barchart.plot  + theme(legend.position = "none")
p_i <- demog_cov_edu_mmonly_bipolar_barchart.plot + theme(legend.position = "none")

mmonly_checks <-
  wrap_plots(p_a, p_b, nrow = 1) /
  free(p_c) /
  wrap_plots(p_d, p_e, nrow = 1) /
  wrap_plots(p_f, p_g, p_h, p_i, nrow = 1) +
  plot_layout(
    heights = c(1, 2.4, 1.4, 1)
  ) +
  plot_annotation(tag_levels = "a")

ggsave("mmonly_checks.png", mmonly_checks, height = 18, width = 14, units = "in")

#======================================================
#
# INSTABILITY analysis
#
#======================================================

pmmdgroup_cgis <- getCGIS(unique(as.numeric(pmmd_group$person_id)))
saveDf(pmmdgroup_cgis, "pmmdgroup_cgis")
nopmmdgroup_cgis <- getCGIS(unique(as.numeric(nopmmd_group$person_id)))
saveDf(nopmmdgroup_cgis, "nopmmdgroup_cgis")

#======================================================
#
# EMERGENCE analysis
#
#======================================================

bipolar_epall <- rbind(pmmd_epall %>% mutate(group = "PMMD"), nopmmd_epall %>% mutate(PMMD = NA, group = "noPMMD"))

observation_period <- getDf("observation_period")
observation_period$person_id <- as.character(observation_period$person_id)

# MM all
observation_period_all <- observation_period %>% 
  dplyr::select(person_id, observation_period_start_date, observation_period_end_date) %>% 
  mutate(observation_duration = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days"))/365) %>%
  inner_join(id_group, .) 

#id_group_dx <- getDf("id_group_dx") # hash out as have new id_group_dx
observation_period_all_dx <- observation_period %>% 
  dplyr::select(person_id, observation_period_start_date, observation_period_end_date) %>% 
  mutate(observation_duration = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days"))/365) %>%
  inner_join(id_group_dx, .) 

observation_period_all_compare <- rbind(data.frame(observation_period_all %>% filter(!is.na(group)) %>% dplyr::select(-c("addpmmd", "samedayswitch_anyorder", "rapidcycler")), def = "additional case-finding"), 
    data.frame(observation_period_all_dx, def = "dx only")) %>%
  mutate(group = gsub("PMMD", "MDS", group))

# Compare date of observation_period start
ggplot(observation_period_all_compare, aes(x = observation_period_start_date, colour = group)) +
  geom_density() +
  theme_bw() +
  scale_colour_manual(values = millerstone) +
  facet_wrap(~def)
ggsave("group_observation_start_density.png")  

# Compare years of observation_period total
ggplot(observation_period_all_compare, aes(x = observation_duration, colour = group)) +
  geom_density() +
  theme_bw() +
  scale_colour_manual(values = millerstone) +
  facet_wrap(~def)
ggsave("group_observation_years_density.png")

ggplot(observation_period_all_compare, aes(x = observation_duration, colour = group)) +
  stat_ecdf(geom = "step") +  
  theme_bw() +
  scale_colour_manual(values = millerstone) +
  facet_wrap(~def)
ggsave("group_observation_years_ecdf.png")  

# Statistical test for difference
tmp <- observation_period_all %>% filter(group %in% c("PMMD", "noPMMD"))
summary(aov(observation_duration ~ group, data = tmp))
tidy(lm(observation_duration ~ group, data = tmp))
tmp_dx <- observation_period_all_dx %>% filter(group %in% c("PMMD", "noPMMD"))
summary(aov(observation_duration ~ group, data = tmp_dx))
tidy(lm(observation_duration ~ group, data = tmp_dx))

#-----------------
# Articulate the case that need to subset to u30 for emergence
#-----------------

# BIPOLAR ONLY
observation_period_bipolar <- observation_period %>% 
  dplyr::select(person_id, observation_period_start_date, observation_period_end_date) %>% 
  inner_join(bipolar_epall, .) %>%
  inner_join(demog_cov %>% dplyr::select(person_id, birth_year)) %>%
  mutate(birth_year = as.Date(paste(birth_year, "-01-01", sep = "")))

# Compare age at observation_period start
observation_period_bipolar_firstep <- observation_period_bipolar %>%
  group_by(person_id) %>%
  slice_min(diagnosis_date, n=1, with_ties = FALSE) %>%
  mutate(age_obsstart = as.numeric(difftime(observation_period_start_date, birth_year, units = "days"))/365) %>%
  mutate(age_firstep = as.numeric(difftime(diagnosis_date, birth_year, units = "days"))/365) %>%
  mutate(years_obsstart_firstep = age_obsstart - age_firstep) 

observation_period_bipolar_firstep.long <- observation_period_bipolar_firstep %>%
  pivot_longer(cols = c("age_obsstart", "age_firstep"), names_to = "age_at", values_to = "age")

observation_period_bipolar_firstep.long %>% group_by(group, age_at) %>%
  summarise(n=n(), mean = mean(age), sd = sd(age), median = median(age))

ggplot(observation_period_bipolar_firstep.long, aes(x = group, y = age, fill = group)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(alpha=0.5) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  theme(legend.position = "bottom") +
  facet_wrap(~age_at, nrow = 1)
ggsave("observation_period_bipolar_pmmdvsnopmmd_ageobsstart_agefirstep_violin.png", height = 4, width = 4)

#------------------------------------------------------
# Time from 1st episode of depressive to 1st mm episode
#------------------------------------------------------

pmmd_prepmmddep <- pmmd_ep_temporal %>% filter(state == "dep" & reltoPMMD1 == "before")
pmmd_firstmm <- pmmd_epall %>% 
  filter(state == "mm") %>%
  group_by(person_id) %>%
  filter(mm_episode == 1) %>%
  filter(dep_episode == min(dep_episode)) %>%
  dplyr::select(person_id, diagnosis_date, dep_episode) %>%
  rename(diagnosis_date_mm = diagnosis_date) %>%
  rename(mm1_depepcount = dep_episode)
pmmd_firstdep <- pmmd_epall %>% 
  filter(state == "dep") %>%
  group_by(person_id) %>%
  filter(dep_episode == 1) %>%
  filter(mm_episode == min(mm_episode)) %>%
  dplyr::select(person_id, diagnosis_date, mm_episode) %>%
  rename(diagnosis_date_dep = diagnosis_date) %>%
  rename(dep1_mmepcount = mm_episode)

mm_split.tmp <- mm_split %>% dplyr::select(person_id, diagnosis_date, mm_split) %>% rename(diagnosis_date_mm = diagnosis_date)
pmmd_firstep <- inner_join(pmmd_firstmm, pmmd_firstdep) %>%
  mutate(time_dep_pre_mm = as.numeric(difftime(diagnosis_date_mm, diagnosis_date_dep, units = "days"))) %>%
  left_join(., id_group) %>%
  left_join(., mm_split.tmp)

nopmmd_prepmmddep <- nopmmd_epall %>% filter(state == "dep" & mm_episode == 0)
nopmmd_firstmm <- nopmmd_epall %>% 
  filter(state == "mm") %>%
  group_by(person_id) %>%
  filter(mm_episode == 1) %>%
  filter(dep_episode == min(dep_episode)) %>%
  dplyr::select(person_id, diagnosis_date, dep_episode) %>%
  rename(diagnosis_date_mm = diagnosis_date) %>%
  rename(mm1_depepcount = dep_episode)
nopmmd_firstdep <- nopmmd_epall %>% 
  filter(state == "dep") %>%
  group_by(person_id) %>%
  filter(dep_episode == 1) %>%
  filter(mm_episode == min(mm_episode)) %>%
  dplyr::select(person_id, diagnosis_date, mm_episode) %>%
  rename(diagnosis_date_dep = diagnosis_date) %>%
  rename(dep1_mmepcount = mm_episode)

mm_split.tmp <- mm_split %>% dplyr::select(person_id, diagnosis_date, mm_split) %>% rename(diagnosis_date_mm = diagnosis_date)

nopmmd_firstep <- inner_join(nopmmd_firstmm, nopmmd_firstdep) %>%
  mutate(time_dep_pre_mm = as.numeric(difftime(diagnosis_date_mm, diagnosis_date_dep, units = "days"))) %>%
  left_join(., id_group) %>%
  left_join(., mm_split.tmp)

bipolar_firstep <- full_join(pmmd_firstep, nopmmd_firstep) %>%
  mutate(mm_split = ifelse(rapidcycler == "rapidcycler", "rapidcycler", mm_split)) %>%
  mutate(mm_split = factor(mm_split, levels = c("mania", "mixed", "rapidcycler"))) %>%
  left_join(., demog_cov %>% dplyr::select(person_id, birth_year)) %>%
  mutate(age_dep1 = as.numeric(difftime(diagnosis_date_dep, birth_year, units = "days"))/365) %>%
  mutate(age_mm1 = as.numeric(difftime(diagnosis_date_mm, birth_year, units = "days"))/365) %>%
  mutate(first_ep = case_when(diagnosis_date_mm < diagnosis_date_dep ~ "mm",
                              diagnosis_date_mm > diagnosis_date_dep ~ "dep",
                              .default = NA)) %>%
  mutate(pmmdgroup01 = case_when(group == "PMMD" ~ 1,
                            group == "noPMMD" ~ 0, 
                            .default = NA))

bipolar_firstep %>% 
  filter(!is.na(age_dep1)) %>%
  group_by(group) %>% 
  summarise(mean = mean(age_dep1), sd = sd(age_dep1), med = median(age_dep1), q1 = quantile(age_dep1, 0.25), q3 = quantile(age_dep1, 0.75))

# What is the time difference from dep1 to mm1?
# - use this one for figure purposes
ggplot(bipolar_firstep %>% filter(time_dep_pre_mm > 0) %>% filter(!is.na(pmmdgroup01)), aes(x = mm_split, y = time_dep_pre_mm/365, fill = group)) +
  geom_boxplot() + 
  geom_point(alpha = 0.05, aes(colour = group), position = position_jitterdodge()) +
  theme_bw() +
  theme(legend.text = element_text(size = 6), legend.title = element_text(size = 6), legend.margin = margin(0, 0, 0, 0)) +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  ylim(0,10) +
  ylab("Years: 1st dep to 1st MM episode") +
  xlab("Category") +
  theme(legend.position = "right")
ggsave("bipolar_firstep_time_dep_pre_mm_3650d_subsetdepfirst.png", height = 3.3, width = 4)

ggplot(bipolar_firstep, aes(x = mm_split, y = time_dep_pre_mm/365, fill = group)) +
  geom_boxplot() + 
  geom_point(alpha = 0.01, aes(colour = group), position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  ylim(-10,10) +
  ylab("Years from 1st dep to 1st MM episode") +
  xlab("Category") +
  theme(legend.position = "bottom")
ggsave("bipolar_firstep_time_dep_pre_mm_3650d.png", height = 3.3, width = 4.5)

bipolar_firstep %>% 
  filter(time_dep_pre_mm > 0) %>% 
  filter(age_dep1 <= 30) %>%
ggplot(., aes(x = mm_split, y = time_dep_pre_mm/365, fill = group)) +
  geom_boxplot() + 
  geom_point(alpha = 0.05, aes(colour = group), position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  ylim(0,10) +
  ylab("Years from 1st depressive episode to 1st manic/mixed episode recorded") +
  xlab("Category") +
  theme(legend.position = "bottom")
ggsave("bipolar_firstep_time_dep_pre_mm_3650d_subsetdepfirst_ageu30.png", height = 3.3, width = 4.5)

bipolar_firstep %>% 
  filter(time_dep_pre_mm > 0) %>% 
  filter(age_dep1 <= 20) %>%
ggplot(., aes(x = mm_split, y = time_dep_pre_mm/365, fill = group)) +
  geom_boxplot() + 
  geom_point(alpha = 0.05, aes(colour = group), position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  ylim(0,10) +
  ylab("Years from 1st depressive episode to 1st manic/mixed episode recorded") +
  xlab("Category") +
  theme(legend.position = "bottom")
ggsave("bipolar_firstep_time_dep_pre_mm_3650d_subsetdepfirst_ageu20.png")

# Linear models
summary(lm(time_dep_pre_mm ~ group + mm_split, data = bipolar_firstep))
summary(lm(time_dep_pre_mm ~ group + mm_split, data = bipolar_firstep %>% filter(time_dep_pre_mm > 0)))
summary(lm(time_dep_pre_mm ~ group + mm_split, data = bipolar_firstep %>% filter(time_dep_pre_mm > 0 & age_dep1 <= 30)))
summary(lm(time_dep_pre_mm ~ group + mm_split, data = bipolar_firstep %>% filter(time_dep_pre_mm > 0 & age_dep1 <= 20)))

#------------------------------------------------------
# Number of mm episodes until first PMMD
#------------------------------------------------------

# people to include
bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")
mds1m_id <- bipolar_epall_count_mm_to_depafter_inf %>% filter(mm_to_depafter <= 1) %>% pull(person_id)

# only look at diagnoses when people are aged <30y
prepmmd1_ep_number <- bipolar_epall_count_mm_to_depafter_inf %>% # used to be pmmd_ep_temporal
  filter(person_id %in% mds1m_id) %>% # only include individuals with MDS at some point
  #mutate(mds1m = ifelse(mm_to_depafter <= 1, "mds1m", "no")) %>%
  filter(mm_to_depafter <= 1) %>% # filter for MDS episodes among these individuals
  group_by(person_id) %>% 
  summarise(min_mmep = min(mm_episode), diagnosis_date_mm = min(diagnosis_date_mm)) %>% # find the mm index for the first MDS episode
  left_join(., demog_cov) %>%
  mutate(cov_age = as.numeric(difftime(diagnosis_date_mm, birth_year, units = "days"))/365) %>%
  mutate(ageu30 = ifelse(cov_age <= 30, "1st Mania: <30y", NA)) %>%
  mutate(ageu25 = ifelse(cov_age <= 25, "1st Mania: <25y", NA)) %>%
  mutate(ageu20 = ifelse(cov_age <= 20, "1st Mania: <20y", NA)) %>%
  filter(!is.na(ageu30)) %>%
  pivot_longer(cols = c("ageu30", "ageu25", "ageu20"), names_to = "age_category", values_to = "age") %>%
  filter(!is.na(age)) %>%
  group_by(age, min_mmep) %>%
  summarise(n=n()) %>% ungroup %>% group_by(age) %>%
  mutate(percent = n/sum(n)*100)
saveDf(prepmmd1_ep_number, "prepmmd1_ep_number")

prepmmd1_ep_number.plot <- ggplot(prepmmd1_ep_number, aes(x = min_mmep, y = percent, fill = age)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("MM episode at 1st Mania-Dep Switch") +
  ylab("Percentage of group (%)") +
  scale_fill_manual(values = millerstone[c(2,3,9)], name = "Age threshold") +
  facet_wrap(~ age, ncol = 3)
ggsave("preMDS1_mmep_number_percent.png", height = 3.3, width = 4)

#------------------------------------------------------
# Compare first mood episodes
#------------------------------------------------------

# Identify individuals with MDS
bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")
mds1m <- bipolar_epall_count_mm_to_depafter_inf %>% filter(mm_to_depafter <= 1) %>% dplyr::select(person_id) %>% distinct() %>% mutate(group_MDS1m = "yes")
bipolar_epall_mm_excl0 <- getDf("bipolar_epall_mm_excl0") %>% dplyr::select(person_id, mm_episode, state) %>% mutate(keep = "keep")

# Find 1st mood episode ever
bipolar_epall_mds1m_moodep1 <- bipolar_epall %>% left_join(., bipolar_epall_mm_excl0) %>% filter(!(is.na(keep) & state == "mm")) %>%
  filter(mm_episode == 0 | dep_episode == 0) %>%
  filter(mm_episode == 1 | dep_episode == 1) %>%
  left_join(., mds1m) %>%  mutate(group_MDS1m = ifelse(is.na(group_MDS1m), "no", group_MDS1m)) %>%
  inner_join(., demog_cov) %>%
  mutate(cov_age = as.numeric(difftime(diagnosis_date, birth_year, units = "days"))/365) %>%
  mutate(ageu30 = ifelse(cov_age <= 30, "1st Mania: <30y", NA)) %>%
  mutate(ageu25 = ifelse(cov_age <= 25, "1st Mania: <25y", NA)) %>%
  mutate(ageu20 = ifelse(cov_age <= 20, "1st Mania: <20y", NA)) %>%
  filter(!is.na(ageu30)) %>%
  pivot_longer(cols = c("ageu30", "ageu25", "ageu20"), names_to = "age_category", values_to = "age")

bipolar_firstep_prop <- bipolar_epall_mds1m_moodep1 %>% 
  group_by(age, group_MDS1m, state) %>%
  summarise(n=n()) %>%
  ungroup %>% group_by(age, group_MDS1m) %>%
  mutate(percent = n/sum(n)*100) %>%
  filter(!is.na(age))
saveDf(bipolar_firstep_prop, "bipolar_firstep_prop")

bipolar_firstep_prop_mm.plot <- ggplot(bipolar_firstep_prop %>% filter(state == "mm"), aes(x = group_MDS1m, y = percent, fill = age)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = millerstone[c(2,3,9)], name = "Age threshold") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("% with Mania as 1st episode") + xlab("Ever Mania-Dep Switch?") +
  facet_wrap(~ age, ncol = 3)
ggsave("bipolar_firstep_prop.png", height = 3.3, width = 4)

library(patchwork)
prepmmd1_ep_number.plot + bipolar_firstep_prop_mm.plot + plot_annotation(tag_levels = 'a')
ggsave("bipolar_natural_history.png", height = 3.3, width = 8)

#======================================================
# Rate of episodes (episodes / year)
#======================================================

# Get observation window
observation_period <- getDf("observation_period")
observation_period$person_id <- as.character(observation_period$person_id)
bipolar <- id_group %>% filter(group %in% c("PMMD", "noPMMD"))
observation_period_bipolar <- observation_period %>% filter(person_id %in% bipolar$person_id) %>%
  dplyr::select(person_id, observation_period_start_date, observation_period_end_date) %>%
  mutate(obs_period = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days")))

# Get episodes
pmmd_epall <- getDf("pmmd_epall")
nopmmd_epall <- getDf("nopmmd_epall")

bipolar_epall <- pmmd_epall %>% 
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMDep", "PMMDep"), 
         group = "PMMD") %>%
  full_join(., nopmmd_epall %>% mutate(group = "noPMMD")) %>%
  mutate(PMMD = ifelse(is.na(PMMD), "noPMMD", PMMD)) 

bipolar_epall_count <- bipolar_epall %>%
  group_by(person_id, state, group) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = state, values_from = n)

bipolar_epall_count_obsperiod <- left_join(bipolar_epall_count, observation_period_bipolar) %>%
  mutate(ep_annual_rate_mm = mm/(obs_period/365)) %>% 
  mutate(ep_annual_rate_dep = dep/(obs_period/365)) %>%
  full_join(., id_group) %>%
  filter(obs_period > 365) %>%
  pivot_longer(cols = c("ep_annual_rate_mm", "ep_annual_rate_dep"), names_to = "state", values_to = "ep_rate") %>%
  mutate(state = ifelse(state == "ep_annual_rate_mm", "Annual mm episode rate", "Annual dep episode rate"))

ggplot(bipolar_epall_count_obsperiod, aes(x = rapidcycler, y = ep_rate, fill = group)) +
  geom_boxplot() + 
  geom_point(alpha = 0.01, colour = "gray30", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  facet_wrap(~ state)
ggsave("bipolar_epall_count_obsperiod.png")

ggplot(bipolar_epall_count_obsperiod, aes(x = rapidcycler, y = ep_rate, fill = group)) +
  geom_boxplot() + 
  geom_point(alpha = 0.01, colour = "gray30", position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  ylim(0,3) +
  facet_wrap(~ state)
ggsave("bipolar_epall_count_obsperiod_eprate3.png")

#------------------------------------------------------
# Find PMMD onset peak
#------------------------------------------------------

# Find preceding and following depressive episode
bipolar_epall_mm <- bipolar_epall %>% filter(state == "mm") %>% 
  left_join(., pmmd_ep %>% dplyr::select(person_id, diagnosis_date, timetodep)) %>%
  left_join(., mm_split %>% dplyr::select(person_id, diagnosis_date, diagnosis_date_end, mm_split)) %>%
  rename(diagnosis_date_mm = diagnosis_date, diagnosis_date_end_mm = diagnosis_date_end) %>%
  dplyr::select(person_id, mm_episode, dep_episode, state, mm_split, diagnosis_date_mm, diagnosis_date_end_mm, group)
bipolar_epall_dep <- bipolar_epall %>% filter(state == "dep") %>% 
  rename(diagnosis_date_dep = diagnosis_date, diagnosis_date_end_dep = diagnosis_date_end) %>%
  dplyr::select(person_id, dep_episode, diagnosis_date_dep, diagnosis_date_end_dep, group)

# - after mm episode
bipolar_epall_count_mm_to_depafter <- bipolar_epall_mm %>%
  mutate(dep_episode_after = dep_episode + 1) %>%
  dplyr::select(-dep_episode) %>%
  rename(dep_episode = dep_episode_after) %>%
  left_join(., bipolar_epall_dep) %>%
  mutate(mm_to_depafter = as.numeric(difftime(diagnosis_date_dep, diagnosis_date_mm, units = "days"))/30) %>%
  filter(dep_episode >= 1) %>% # exclude if no subsequent depressive episodes
  group_by(person_id, dep_episode) %>% # for each dep episode, take the least time since mania
  filter(mm_to_depafter == min(mm_to_depafter))

# - before mm episode
bipolar_epall_count_mm_to_depbefore <- bipolar_epall_mm %>%
  left_join(., bipolar_epall_dep) %>% # no need to do dep episode index change as for the after (as already linked)
  mutate(mm_to_depbefore = as.numeric(difftime(diagnosis_date_dep, diagnosis_date_mm, units = "days"))/30) %>%
  filter(dep_episode >= 1) %>% # exclude if no previous dep episodes
  group_by(person_id, dep_episode) %>% # for each dep episode, take the least-negative time since mania
  filter(mm_to_depbefore == max(mm_to_depbefore))

bipolar_epall_count_mm_to_depafter.tmp <- bipolar_epall_count_mm_to_depafter %>% 
  rename(timetodep = mm_to_depafter) %>% 
  mutate(type = "Depression after")
bipolar_epall_count_mm_to_depbefore.tmp <- bipolar_epall_count_mm_to_depbefore %>% 
  mutate(timetodep = mm_to_depbefore) %>% 
  mutate(type = "Depression before")

bipolar_epall_count_mm_to_deprecent <- full_join(bipolar_epall_count_mm_to_depbefore.tmp, bipolar_epall_count_mm_to_depafter.tmp) %>%
  mutate(type = factor(type, levels = c("Depression before", "Depression after")))

ggplot(bipolar_epall_count_mm_to_deprecent, aes(x = timetodep, colour = mm_split)) +
  geom_density() +
  theme_bw() +
  scale_colour_manual(values = millerstone[c(1,4,2)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 3)) + 
  facet_grid(.~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_deprecent.png", height = 2, width = 6)

#ggplot(bipolar_epall_count_mm_to_deprecent %>% filter(mm_split == "mania") %>% filter(timetodep %in% -24:24), aes(x = timetodep, colour = type)) + ##### too many with timedep == 0?
ggplot(bipolar_epall_count_mm_to_deprecent %>% filter(mm_split == "mania") %>% filter(timetodep != 0 & timetodep %in% -24:24), aes(x = timetodep, colour = type)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 3)) + 
  facet_grid(.~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_deprecent_mania.png", height = 2.2, width = 4)

ggplot(bipolar_epall_count_mm_to_deprecent %>% filter(mm_split == "mixed") %>% filter(timetodep != 0 & timetodep %in% -24:24), aes(x = timetodep, colour = type)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 3)) + 
  facet_grid(.~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_deprecent_mixed.png", height = 2.2, width = 4)

ggplot(bipolar_epall_count_mm_to_deprecent %>% filter(timetodep != 0 & timetodep %in% -24:24), aes(x = timetodep, colour = type)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 3)) + 
  facet_grid(mm_split~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_deprecent_maniamixed.png", height = 2.2, width = 4)


bipolar_epall_count_mm_to_depafter_survfit <- survfit(Surv(mm_to_depafter, dep_episode_after_ind) ~ group, data = bipolar_epall_count_mm_to_depafter)
ggsurvplot(bipolar_epall_count_mm_to_depafter_survfit, data = bipolar_epall_count_mm_to_depafter,
          surv.median.line = "hv",
          legend.title = "Years: mania/mixed episode to depression AFTER",
          pval = TRUE, conf.int = TRUE,
          palette = millerstone[c(7,8)], ggtheme = theme_bw())
ggsave("bipolar_epall_count_mm_to_depafter_survfit.png") # doesn't make sense to display as the separation just shows ascertainment ...

#------------------------------------------------------
# Compare depression rate between groups
#------------------------------------------------------

# - overall depression rate
bipolar_epall_count_obsperiod


#============================================================
# Does antipsychotic prescribing explain depression?
#============================================================

# Read in data
pmmdgroup_mm_7d_anyvisit_psychmed <- getDf("pmmdgroup_mm_7d_anyvisit_psychmed")
nopmmdgroup_mm_7d_anyvisit_psychmed <- getDf("nopmmdgroup_mm_7d_anyvisit_psychmed")

bipolar_mm_7d_anyvisit_psychmed <- full_join(pmmdgroup_mm_7d_anyvisit_psychmed, nopmmdgroup_mm_7d_anyvisit_psychmed) %>%
  left_join(., demog_cov) %>%
  mutate(cov_age_visit = as.numeric(difftime(visit_start_date, birth_year, units = "days")/365))

bipolar_mm_7d_anyvisit_psychmed_ap <- bipolar_mm_7d_anyvisit_psychmed %>%
  group_by(person_id, diagnosis_date, PMMD, mm_split, cov_age_visit, cov_gender, cov_state, cov_race, cov_rel, cov_employ, cov_edu, broad_class) %>%
  summarise(n=n()) %>%
  group_by(person_id, diagnosis_date, PMMD, mm_split, cov_age_visit, cov_gender, cov_state, cov_race, cov_rel, cov_employ, cov_edu, broad_class) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = broad_class, values_from = n) %>%
  mutate_if(is.numeric, replace_na, replace = 0) %>%
  mutate(pmmd01 = ifelse(PMMD == "PMMD", 1, 0))

pmmd_mm_7d_anyvisit_psychmed_ap <- pmmdgroup_mm_7d_anyvisit_psychmed %>%
  left_join(., demog_cov) %>%
  mutate(cov_age_visit = as.numeric(difftime(visit_start_date, birth_year, units = "days")/365)) %>%
  group_by(person_id, diagnosis_date, PMMD, mm_split, cov_age_visit, cov_gender, cov_state, cov_race, cov_rel, cov_employ, cov_edu, broad_class) %>%
  summarise(n=n()) %>%
  group_by(person_id, diagnosis_date, PMMD, mm_split, cov_age_visit, cov_gender, cov_state, cov_race, cov_rel, cov_employ, cov_edu, broad_class) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = broad_class, values_from = n) %>%
  mutate_if(is.numeric, replace_na, replace = 0) %>%
  mutate(pmmd01 = ifelse(PMMD == "PMMD", 1, 0))

# pmmd01 ~ antipsychotic use during mania
bipolar_pmmd01_antipsychotic_gee <- gee(pmmd01 ~ antipsychotic + mm_split + cov_age_visit + cov_gender,
               data = bipolar_mm_7d_anyvisit_psychmed_ap, 
               id = person_id, 
               family = binomial,
               corstr = "exchangeable")
bipolar_pmmd01_antipsychotic_gee_coef <- summary(bipolar_pmmd01_antipsychotic_gee)$coefficients %>% data.frame %>% rownames_to_column("term") %>%
  mutate(EstimateOR = exp(Estimate)) %>%
  mutate(p = ifelse(Robust.z <= 0, 2*pnorm(q = -1*Robust.z, lower.tail=FALSE), 2*pnorm(q = Robust.z, lower.tail=FALSE)))
saveDf(bipolar_pmmd01_antipsychotic_gee_coef, "bipolar_pmmd01_antipsychotic_gee_coef")

pmmd_pmmd01_antipsychotic_gee <- gee(pmmd01 ~ antipsychotic + mm_split + cov_age_visit + cov_gender,
               data = pmmd_mm_7d_anyvisit_psychmed_ap, 
               id = person_id, 
               family = binomial,
               corstr = "exchangeable")
pmmd_pmmd01_antipsychotic_gee_coef <- summary(pmmd_pmmd01_antipsychotic_gee)$coefficients %>% data.frame %>% rownames_to_column("term") %>%
  mutate(EstimateOR = exp(Estimate)) %>%
  mutate(p = ifelse(Robust.z <= 0, 2*pnorm(q = -1*Robust.z, lower.tail=FALSE), 2*pnorm(q = Robust.z, lower.tail=FALSE)))
saveDf(pmmd_pmmd01_antipsychotic_gee_coef, "pmmd_pmmd01_antipsychotic_gee_coef")

# Plot results
# - all bipolar
bipolar_mm_7d_anyvisit_psychmed_ap_barplot <- bipolar_mm_7d_anyvisit_psychmed_ap %>%
  group_by(pmmd01, antipsychotic) %>%
  summarise(n=n()) %>% 
  group_by(antipsychotic) %>%
  rename(episode = pmmd01) %>%
  mutate(episode = ifelse(episode == 0, "nonPMMD", "PMMD")) %>%
  mutate(antipsychotic = ifelse(antipsychotic == 0, "no", "yes")) %>%
  pivot_wider(names_from = antipsychotic, values_from = n) %>%
  mutate(sum = no+yes) %>%
  mutate(no = no/sum, yes = yes/sum) %>%
  pivot_longer(cols = c("no", "yes"), names_to = "antipsychotic", values_to = "prop") %>%
  mutate(comparison = "PMMDgroup + noPMMDgroup")

# - PMMD group only
pmmd_mm_7d_anyvisit_psychmed_ap_barplot <- pmmd_mm_7d_anyvisit_psychmed_ap %>%
  group_by(pmmd01, antipsychotic) %>%
  summarise(n=n()) %>% 
  group_by(antipsychotic) %>%
  rename(episode = pmmd01) %>%
  mutate(episode = ifelse(episode == 0, "nonPMMD", "PMMD")) %>%
  mutate(antipsychotic = ifelse(antipsychotic == 0, "no", "yes")) %>%
  pivot_wider(names_from = antipsychotic, values_from = n) %>%
  mutate(sum = no+yes) %>%
  mutate(no = no/sum, yes = yes/sum) %>%
  pivot_longer(cols = c("no", "yes"), names_to = "antipsychotic", values_to = "prop") %>%
  mutate(comparison = "within PMMDgroup")

# - combine
rbind(bipolar_mm_7d_anyvisit_psychmed_ap_barplot, pmmd_mm_7d_anyvisit_psychmed_ap_barplot) %>%
  ggplot(., aes(x = antipsychotic, y = prop, fill = episode)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  theme_bw() +
  ylab("Proportion of MM episodes") + xlab("Antipsychotic during MM episode") +
  facet_wrap(~ comparison)
ggsave("bipolar_mm_7d_anyvisit_psychmed_ap_barplot.png", height = 3.3, width = 6)

#============================================================
# 
# Diagnosis code flipping between mixed and depressive ep
#
#============================================================

pmmd_ep_temporal <- getDf("pmmd_ep_temporal")
mm_split <- getDf("mm_split")

ep_sum <- pmmd_ep_temporal %>%
  group_by(person_id, state) %>%
  summarise(total_ep=n())

pmmd_ep_temporal_mm_split <- left_join(pmmd_ep_temporal, mm_split) %>%
  left_join(demog_cov) %>%
  mutate(cov_age_ep = as.numeric(difftime(diagnosis_date, birth_year, units = "days")/365)) %>%
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  left_join(., id_group) %>%
  left_join(., ep_sum)

library(gee)
#pmmd_vs_mixeddx_gee <- gee(pmmd01 ~ reltoPMMD1 + rapidcycler + cov_age_ep + cov_gender + cov_state + cov_race + cov_rel + cov_employ + cov_edu,

bipolar_epall_mmsplit <- bipolar_epall %>% 
  filter(state == "mm") %>%
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  left_join(., mm_split)

bipolar_epall_mmsplit %>% filter(pmmd01 == 1) %>% group_by(mm_split) %>% summarise(n=n())

glm(pmmd01 ~ mm_split, data = bipolar_epall_mmsplit, family=poisson(link='log')) %>%
  coeftest(., vcov = vcovCL, cluster = ~person_id) %>%
  tidy %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)), 
        higherRR = exp(estimate + (std.error * 1.96)))

pmmd_vs_mixeddx_gee <- gee(pmmd01 ~ mm_split + rapidcycler + total_ep + cov_age_ep + cov_gender,
               data = pmmd_ep_temporal_mm_split, 
               id = person_id, 
               family = binomial,
               corstr = "exchangeable")
pmmd_vs_mixeddx_gee_coef <- summary(pmmd_vs_mixeddx_gee)$coefficients %>% data.frame %>% rownames_to_column("term") %>%
  mutate(EstimateOR = exp(Estimate)) %>%
  mutate(p = ifelse(Robust.z <= 0, 2*pnorm(q = -1*Robust.z, lower.tail=FALSE), 2*pnorm(q = Robust.z, lower.tail=FALSE)))
saveDf(pmmd_vs_mixeddx_gee_coef, "pmmd_vs_mixeddx_gee_coef")

pmmd_vs_mixeddx <- pmmd_ep_temporal_mm_split %>%
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMD", "PMMD")) %>%
  filter(state == "mm") %>%
  group_by(PMMD, mm_split) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = mm_split, values_from = n) %>%
  mutate(sum = mania+mixed) %>%
  mutate(mania_prop = mania/sum, mixed_prop = mixed/sum)

print(pmmd_vs_mixeddx)
fisher.test(as.matrix(pmmd_vs_mixeddx[,2:3]))

#============================================================
#
# Distribution of PMMD episodes diagnosis date
#
#============================================================

# - broader PMMD definition
bipolar_epall <- pmmd_ep_temporal %>%
  mutate(PMMDep = ifelse(is.na(PMMD), "nonPMMD", "PMMD")) %>%
  full_join(., nopmmd_epall %>% mutate(group = "noPMMD", PMMDep = "nonPMMD") %>% mutate(duration = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days")))) %>%
  mutate(group = gsub("PMMD", "MDS", group), PMMDep = gsub("PMMD", "MDS", PMMDep)) %>%
  mutate(category = factor(paste("group", group, "_", "ep", PMMDep, sep = ""), levels = c("groupMDS_epMDS", "groupMDS_epnonMDS", "groupnoMDS_epnonMDS"))) %>%
  mutate(definition = "Broad definition")

# - dx only PMMD definition
pmmd_epall_dx <- getDf("pmmd_epall_dx")
nopmmd_epall_dx <- getDf("nopmmd_epall_dx")
bipolar_epall_dx <- pmmd_epall_dx %>%
  mutate(group = "PMMD", PMMDep = ifelse(is.na(PMMD), "nonPMMD", "PMMD")) %>%
  full_join(., nopmmd_epall_dx %>% mutate(group = "noPMMD", PMMDep = "nonPMMD")) %>%
  mutate(group = gsub("PMMD", "MDS", group), PMMDep = gsub("PMMD", "MDS", PMMDep)) %>%
  mutate(category = factor(paste("group", group, "_", "ep", PMMDep, sep = ""), levels = c("groupMDS_epMDS", "groupMDS_epnonMDS", "groupnoMDS_epnonMDS"))) %>%
  mutate(definition = "ICD diagnosis-only definition") %>%
  dplyr::select(-duration)

ggplot(full_join(bipolar_epall, bipolar_epall_dx) %>% filter(state == "mm"), aes(x = diagnosis_date, colour = category)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = millerstone[c(8,2,7)]) +
  facet_wrap(~ definition, nrow = 2)
ggsave("bipolar_mm_epall_dxdate_density.png", height = 8, width = 9)

#============================================================
#
# Natural history differences
#
#============================================================

#============================================================
# Number of each episode pole
#============================================================

bipolar_episodes_count <- getDf("bipolar_episodes_count")

bipolar_episodes_count_numbermm_dep <- bipolar_episodes_count %>%
  left_join(., id_group_clean %>% dplyr::select(person_id, rapidcycler)) %>%
  pivot_longer(cols = c("mm", "dep"), names_to = "polarity", values_to = "episodes")

bipolar_episodes_count %>%
  left_join(., id_group %>% dplyr::select(person_id, rapidcycler)) %>%
  filter(mm>1) %>%
  group_by(group, rapidcycler) %>%
  summarise(n=n())

ggplot(bipolar_episodes_count_numbermm_dep, aes(x = polarity, y = episodes, fill = group)) +
  geom_boxplot() + 
  geom_point(alpha = 0.05, aes(colour = group), position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  ylab("Number of episodes") +
  theme(legend.position = "bottom") +
  facet_wrap(~ rapidcycler, scale = "free_y")
ggsave("bipolar_episodes_count_numbermm_dep.png")

#============================================================
# Ratio of episode poles
#============================================================

bipolar_episodes_count_mmdepratio <- bipolar_episodes_count %>%
  mutate(mm_dep_ratio = mm/dep) %>%
  left_join(., id_group %>% dplyr::select(person_id, rapidcycler))

bipolar_episodes_count_mmdepratio %>%
  group_by(group) %>%
  summarise(n=n(), mean = mean(mm_dep_ratio), sd = sd(mm_dep_ratio), med = median(mm_dep_ratio), q1 = quantile(mm_dep_ratio, 0.25), q3 = quantile(mm_dep_ratio, 0.75))

ggplot(bipolar_episodes_count_mmdepratio, aes(x = rapidcycler, y = mm_dep_ratio, fill = group)) +
  geom_boxplot() + 
  #geom_point(alpha = 0.05, aes(colour = group), position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(7,8)]) +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  ylim(0,10) +
  ylab("Ratio of manic/mixed:depressed episodes") +
  theme(legend.position = "bottom")  
ggsave("bipolar_episodes_count_mmdepratio.png")

#============================================================
# Emergence
#============================================================

pmmd_id <- id_group %>% filter(group == "PMMD") %>% pull(person_id)
pmmd_alldx <- getDiagnosis(as.numeric(pmmd_id))
pmmd_alldx <- pmmd_alldx %>% mutate(person_id = as.character(person_id))
nopmmd_id <- id_group %>% filter(group == "noPMMD") %>% pull(person_id)
nopmmd_alldx <- getDiagnosis(as.numeric(nopmmd_id))
nopmmd_alldx <- nopmmd_alldx %>% mutate(person_id = as.character(person_id))

demog <- getDf("demog")
age <- demog %>% mutate(person_id = as.character(person_id)) %>% dplyr::select(person_id, birth_year)

# First ever diagnosis
pmmd_alldx_first <- pmmd_alldx %>% 
  distinct() %>%
  arrange(diagnosis_date) %>% 
  group_by(person_id) %>% 
  slice_min(diagnosis_date) %>%
  ungroup() %>% 
  left_join(., age) %>%
  mutate(birth_year = as.Date(paste(birth_year, "-01-01", sep = ""))) %>%
  mutate(age = as.numeric(difftime(diagnosis_date, birth_year, units = "days")/365)) %>%
  filter(age > 0)

summary(pmmd_alldx_first$age)

pmmd_alldx_first_summary <- pmmd_alldx_first %>%
  group_by(diagnosis) %>% 
  summarise(n=n(), prop = n/length(nopmmd_id), age_mean = mean(age), age_sd = sd(age)) %>%
  arrange(-prop)
pmmd_alldx_first_summary
saveDf(pmmd_alldx_first_summary, "pmmd_alldx_first_summary")

nopmmd_alldx_first <- nopmmd_alldx %>% 
  distinct() %>%
  arrange(diagnosis_date) %>% 
  group_by(person_id) %>% 
  slice_min(diagnosis_date) %>%
  ungroup() %>% 
  left_join(., age) %>%
  mutate(birth_year = as.Date(paste(birth_year, "-01-01", sep = ""))) %>%
  mutate(age = as.numeric(difftime(diagnosis_date, birth_year, units = "days")/365)) %>%
  filter(age > 0)

summary(nopmmd_alldx_first$age)

nopmmd_alldx_first_summary <- nopmmd_alldx_first %>%
  group_by(diagnosis) %>% 
  summarise(n=n(), prop = n/length(nopmmd_id), age_mean = mean(age), age_sd = sd(age)) %>%
  arrange(-prop)
nopmmd_alldx_first_summary

saveDf(nopmmd_alldx_first_summary, "nopmmd_alldx_first_summary", overwrite = T)

pmmd_alldx_first_summary.tmp <- pmmd_alldx_first_summary
colnames(pmmd_alldx_first_summary.tmp) <- c("diagnosis", "n_pmmd", "prop_pmmd", "age_mean_pmmd", "age_sd_pmmd")
nopmmd_alldx_first_summary.tmp <- nopmmd_alldx_first_summary
colnames(nopmmd_alldx_first_summary.tmp) <- c("diagnosis", "n_nopmmd", "prop_nopmmd", "age_mean_nopmmd", "age_sd_nopmmd")
pmmd_nopmmd_alldx_first_summary <- left_join(pmmd_alldx_first_summary.tmp, nopmmd_alldx_first_summary.tmp)
saveDf(pmmd_nopmmd_alldx_first_summary, "pmmd_nopmmd_alldx_first_summary", overwrite = T)

# First mood diagnosis
pmmd_nopmmd_mooddx_first <- full_join(pmmd_epall %>% mutate(group = "PMMD"), nopmmd_epall %>% mutate(group = "noPMMD")) %>%
  mutate(total_n = ifelse(group == "PMMD", nrow(pmmd_group), nrow(nopmmd_group))) %>%
  distinct() %>%
  arrange(diagnosis_date) %>% 
  group_by(person_id) %>% 
  slice_min(diagnosis_date) %>%
  ungroup() %>% 
  left_join(., age) %>%
  mutate(birth_year = as.Date(paste(birth_year, "-01-01", sep = ""))) %>%
  mutate(age = as.numeric(difftime(diagnosis_date, birth_year, units = "days")/365)) %>%
  filter(age > 0)

summary(pmmd_mooddx_first$age)

pmmd_nopmmd_mooddx_first_summary <- pmmd_nopmmd_mooddx_first %>%
  group_by(group, state) %>% 
  summarise(n=n(), age_mean = mean(age), age_sd = sd(age)) %>%
  mutate(prop = ifelse(group == "PMMD", n/4856, n/7873))
pmmd_nopmmd_mooddx_first_summary
saveDf(pmmd_nopmmd_mooddx_first_summary, "pmmd_nopmmd_mooddx_first_summary")

nopmmd_alldx_first <- nopmmd_alldx %>% 
  distinct() %>%
  arrange(diagnosis_date) %>% 
  group_by(person_id) %>% 
  slice_min(diagnosis_date) %>%
  ungroup() %>% 
  left_join(., age) %>%
  mutate(birth_year = as.Date(paste(birth_year, "-01-01", sep = ""))) %>%
  mutate(age = as.numeric(difftime(diagnosis_date, birth_year, units = "days")/365)) %>%
  filter(age > 0)

summary(nopmmd_alldx_first$age)

nopmmd_alldx_first_summary <- nopmmd_alldx_first %>%
  group_by(diagnosis) %>% 
  summarise(n=n(), prop = n/length(nopmmd_id), age_mean = mean(age), age_sd = sd(age)) %>%
  arrange(-prop)
nopmmd_alldx_first_summary

saveDf(nopmmd_alldx_first_summary, "nopmmd_alldx_first_summary", overwrite = T)

pmmd_alldx_first_summary.tmp <- pmmd_alldx_first_summary
colnames(pmmd_alldx_first_summary.tmp) <- c("diagnosis", "n_pmmd", "prop_pmmd", "age_mean_pmmd", "age_sd_pmmd")
nopmmd_alldx_first_summary.tmp <- nopmmd_alldx_first_summary
colnames(nopmmd_alldx_first_summary.tmp) <- c("diagnosis", "n_nopmmd", "prop_nopmmd", "age_mean_nopmmd", "age_sd_nopmmd")
pmmd_nopmmd_alldx_first_summary <- left_join(pmmd_alldx_first_summary.tmp, nopmmd_alldx_first_summary.tmp)
saveDf(pmmd_nopmmd_alldx_first_summary, "pmmd_nopmmd_alldx_first_summary", overwrite = T)

###### TO PUT SOMEWHERE ELSE

#============================================================
# MOVE ELSEWHERE: Time to PMMD switch
#============================================================
# - episode-level statistics vs patient-level median statistics
summary(pmmd_ep$timetodep)
pmmd_ep %>% filter(timetodep > 0) %>% pull(timetodep) %>% summary()
timetodep_dx_person <- pmmd_ep %>% filter(state == "mm") %>% group_by(person_id) %>% summarise(timetodep_dx_med = median(timetodep))
summary(timetodep_dx_person$timetodep_dx_med)

timetodep_dx_rm0_person <- pmmd_ep %>% filter(state == "mm") %>% filter(timetodep > 0) %>% group_by(person_id) %>% summarise(timetodep_dx_med = median(timetodep))
summary(timetodep_dx_rm0_person$timetodep_dx_med)

pmmd_ep %>% filter(timetodep > 0) %>% group_by(person_id) %>% summarise(timetodep_dx_med = median(timetodep)) %>% pull(timetodep_dx_med) %>% summary()

pmmd_ep %>% filter(state == "mm") %>% group_by(person_id) %>% summarise(timetodep_dx_med = median(timetodep)) %>% filter(timetodep_dx_med <= 30) %>% nrow()
pmmd_ep %>% filter(state == "mm") %>% group_by(person_id) %>% summarise(timetodep_dx_med = median(timetodep)) %>% filter(timetodep_dx_med <= 90) %>% nrow()

timetodep_dx_person <- ggplot(timetodep_dx_person, aes(x = timetodep_dx_med)) + 
  geom_density() + 
  geom_vline(xintercept = c(0, 12, 55.5), colour = "gray50", linetype = "dashed") + # numbers come from summary(timetodep_dx_person$timetodep_dx_med)
  theme_bw()
ggsave("timetodep_dx_person.png", timetodep_dx_person)

timetodep_dx_rm0_person <- ggplot(timetodep_dx_rm0_person, aes(x = timetodep_dx_med)) + 
  geom_density() + 
  geom_vline(xintercept = c(12.5, 41, 86), colour = "gray50", linetype = "dashed") +  # numbers come from summary(timetodep_dx_rm0_person$timetodep_dx_med)
  theme_bw()
ggsave("timetodep_dx_rm0_person.png", timetodep_dx_rm0_person)

timetodep_dx_episode <- ggplot(pmmd_ep %>% filter(state == "mm"), aes(x = timetodep_dx)) + 
  geom_density() + 
  geom_vline(xintercept = c(4, 29, 110), colour = "gray50", linetype = "dashed") + 
  theme_bw()
ggsave("timetodep_dx_episode.png", timetodep_dx_episode)

timetodep_dx_rm0_episode <- ggplot(pmmd_ep %>% filter(state == "mm") %>% filter(timetodep_dx > 0), aes(x = timetodep_dx)) + 
  geom_density() + 
  geom_vline(xintercept = c(13, 44, 133), colour = "gray50", linetype = "dashed") + 
  theme_bw()
ggsave("timetodep_dx_rm0_episode.png", timetodep_dx_rm0_episode)

#============================================================
# MOVE ELSEWHERE: Time to PMMD switch
#============================================================

