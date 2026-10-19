# RESPIRATORY INFECTIONS
rti <- c(addICDdecimal('B34'), addICDdecimal('J00'), addICDdecimal('J01'), addICDdecimal('J02'), 
                             addICDdecimal('J04'), addICDdecimal('J05'), addICDdecimal('J06'), 
                             addICDdecimal('J09'), addICDdecimal('J10'), addICDdecimal('J11'), addICDdecimal('J12'))
icdCodes_rti <- list(rti = rti)
allPatients_rti <- getPatients(icdCodes_rti) # n=324049 > 346932

# Overlap with bipolar 
bipolar_mm_clinical <- getDf("bipolar_mm_clinical")
allPatients_rti_bipolar <- allPatients_rti %>% filter(person_id %in% bipolar_mm_clinical$person_id)
bipolar_rti_dx <- getDiagnosis(allPatients_rti_bipolar$person_id, icd_codes = c(rti)) %>% distinct(diagnosis_date, icd_code, .keep_all = TRUE) %>%
  mutate(person_id = as.character(person_id)) %>%
  mutate(rti = "rti") %>%
  rename(rti_date = diagnosis_date) #8050

# ALLERGIC RHINITIS
allergic_rhinitis_codes <- c(
  # ICD-9-CM
  "477",
  "477.0",
  "477.1",
  "477.2",
  "477.8",
  "477.9",
  # ICD-10-CM
  "J30.1",
  "J30.2",
  "J30.5",
  "J30.8",
  "J30.81",
  "J30.89",
  "J30.9"
)

icdCodes_rhinitis <- list(rhinitis = allergic_rhinitis_codes)
allPatients_rhinitis <- getPatients(icdCodes_rhinitis) # n=324049 > 346932

# Overlap with bipolar 
bipolar_mm_clinical <- getDf("bipolar_mm_clinical")
allPatients_rhinitis_bipolar <- allPatients_rhinitis %>% filter(person_id %in% bipolar_mm_clinical$person_id)
bipolar_rhinitis_dx <- getDiagnosis(allPatients_rhinitis_bipolar$person_id, icd_codes = c(allergic_rhinitis_codes)) %>% distinct(diagnosis_date, icd_code, .keep_all = TRUE) %>%
  mutate(person_id = as.character(person_id)) %>%
  mutate(rhinitis = "rhinitis") %>%
  rename(rhinitis_date = diagnosis_date) #14475

# COMPOSITE
icdCodes_composite <- list(rhinitis = allergic_rhinitis_codes, rti = rti)
allPatients_composite <- getPatients(icdCodes_composite) # n=12478115

# Overlap with bipolar 
bipolar_mm_clinical <- getDf("bipolar_mm_clinical")
allPatients_composite_bipolar <- allPatients_composite %>% filter(person_id %in% bipolar_mm_clinical$person_id)
bipolar_composite_dx <- getDiagnosis(allPatients_composite_bipolar$person_id, icd_codes = c(allergic_rhinitis_codes, rti)) %>% distinct(diagnosis_date, icd_code, .keep_all = TRUE) %>%
  mutate(person_id = as.character(person_id)) %>%
  mutate(composite = "composite") %>%
  rename(composite_date = diagnosis_date) #26465

#===================================================================
# Negative control: Diagnosis rate after mania/mixed episodes_per_id
#===================================================================

bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")

# RTI
# Join rti diagnosis date to df
bipolar_epall_count_mm_to_depafter_inf_rti <- bipolar_epall_count_mm_to_depafter_inf %>%
  left_join(., bipolar_rti_dx %>% dplyr::select(person_id, rti, rti_date), relationship = "many-to-many") %>%
  filter(rti_date > diagnosis_date_end_mm | is.na(rti_date)) %>%
  mutate(months_mmend_to_rti = as.numeric(difftime(rti_date, diagnosis_date_end_mm, units = "days"))/30) %>%
  group_by(person_id, diagnosis_date_end_mm) %>%
  slice_min(months_mmend_to_rti) %>% ungroup() %>% distinct()
saveDf(bipolar_epall_count_mm_to_depafter_inf_rti, "bipolar_epall_count_mm_to_depafter_inf_rti", overwrite = T)

# Plot
totalpercent_mds_rti <- bipolar_epall_count_mm_to_depafter_inf_rti %>%
  mutate(ever_rti = ifelse(!is.na(months_mmend_to_rti), "rti", "norti")) %>%
  group_by(mm_split, ever_rti) %>%
  summarise(rti_transition_ep = n()) %>%
  pivot_wider(names_from = "ever_rti", values_from = "rti_transition_ep") %>%
  mutate(total_ep = rti + norti,
         prop_rti = rti/total_ep)

bipolar_epall_count_mm_to_depafter_inf_rti.plot <- bipolar_epall_count_mm_to_depafter_inf_rti %>% 
  left_join(totalpercent_mds_rti %>% dplyr::select(mm_split, prop_rti)) %>%
  rename(monthstorti = months_mmend_to_rti) %>% 
  mutate(monthstorti = ifelse(is.na(monthstorti), 100, monthstorti)) %>%
 # rename(monthstodep = mm_to_depafter) %>% 
 # mutate(monthstodep = ifelse(is.na(monthstodep), 100, monthstodep)) %>%
ggplot(., aes(x = monthstorti, colour = mm_split)) +
  stat_ecdf() +
  theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  xlab("Time to URTI (months)") +
  ylab("Cumulative proportion") +
  geom_hline(aes(yintercept = prop_rti), linetype = "dashed") +
  coord_cartesian(xlim = c(0,24)) + scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  facet_wrap(~mm_split, scales = "free_x")

# Join RHINITIS diagnosis date to df
bipolar_epall_count_mm_to_depafter_inf_rhinitis <- bipolar_epall_count_mm_to_depafter_inf %>%
  left_join(., bipolar_rhinitis_dx %>% dplyr::select(person_id, rhinitis, rhinitis_date), relationship = "many-to-many") %>%
  filter(rhinitis_date > diagnosis_date_end_mm | is.na(rhinitis_date)) %>%
  mutate(months_mmend_to_rhinitis = as.numeric(difftime(rhinitis_date, diagnosis_date_end_mm, units = "days"))/30) %>%
  group_by(person_id, diagnosis_date_end_mm) %>%
  slice_min(months_mmend_to_rhinitis) %>% ungroup() %>% distinct()
saveDf(bipolar_epall_count_mm_to_depafter_inf_rhinitis, "bipolar_epall_count_mm_to_depafter_inf_rhinitis", overwrite = T)

# Plot
totalpercent_mds_rhinitis <- bipolar_epall_count_mm_to_depafter_inf_rhinitis %>%
  mutate(ever_rhinitis = ifelse(!is.na(months_mmend_to_rhinitis), "rhinitis", "norhinitis")) %>%
  group_by(mm_split, ever_rhinitis) %>%
  summarise(rhinitis_transition_ep = n()) %>%
  pivot_wider(names_from = "ever_rhinitis", values_from = "rhinitis_transition_ep") %>%
  mutate(total_ep = rhinitis + norhinitis,
         prop_rhinitis = rhinitis/total_ep)

bipolar_epall_count_mm_to_depafter_inf_rhinitis.plot <- bipolar_epall_count_mm_to_depafter_inf_rhinitis %>% 
  left_join(totalpercent_mds_rhinitis %>% dplyr::select(mm_split, prop_rhinitis)) %>%
  rename(monthstorhinitis = months_mmend_to_rhinitis) %>% 
  mutate(monthstorhinitis = ifelse(is.na(monthstorhinitis), 100, monthstorhinitis)) %>%
 # rename(monthstodep = mm_to_depafter) %>% 
 # mutate(monthstodep = ifelse(is.na(monthstodep), 100, monthstodep)) %>%
ggplot(., aes(x = monthstorhinitis, colour = mm_split)) +
  stat_ecdf() +
  theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  xlab("Time to rhinitis (months)") +
  ylab("Cumulative proportion") +
  geom_hline(aes(yintercept = prop_rhinitis), linetype = "dashed") +
  coord_cartesian(xlim = c(0,24)) + scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  facet_wrap(~mm_split, scales = "free_x")

# COMPOSITE
# Join composite diagnosis date to df
bipolar_epall_count_mm_to_depafter_inf_composite <- bipolar_epall_count_mm_to_depafter_inf %>%
  left_join(., bipolar_composite_dx %>% dplyr::select(person_id, composite, composite_date), relationship = "many-to-many") %>%
  filter(composite_date > diagnosis_date_end_mm | is.na(composite_date)) %>%
  mutate(months_mmend_to_composite = as.numeric(difftime(composite_date, diagnosis_date_end_mm, units = "days"))/30) %>%
  group_by(person_id, diagnosis_date_end_mm) %>%
  slice_min(months_mmend_to_composite) %>% ungroup() %>% distinct()
saveDf(bipolar_epall_count_mm_to_depafter_inf_composite, "bipolar_epall_count_mm_to_depafter_inf_composite", overwrite = T)

# Plot
totalpercent_mds_composite <- bipolar_epall_count_mm_to_depafter_inf_composite %>%
  mutate(ever_composite = ifelse(!is.na(months_mmend_to_composite), "composite", "nocomposite")) %>%
  group_by(mm_split, ever_composite) %>%
  summarise(composite_transition_ep = n()) %>%
  pivot_wider(names_from = "ever_composite", values_from = "composite_transition_ep") %>%
  mutate(total_ep = composite + nocomposite,
         prop_composite = composite/total_ep)

bipolar_epall_count_mm_to_depafter_inf_composite.plot <- bipolar_epall_count_mm_to_depafter_inf_composite %>% 
  left_join(totalpercent_mds_composite %>% dplyr::select(mm_split, prop_composite)) %>%
  rename(monthstocomposite = months_mmend_to_composite) %>% 
  mutate(monthstocomposite = ifelse(is.na(monthstocomposite), 100, monthstocomposite)) %>%
 # rename(monthstodep = mm_to_depafter) %>% 
 # mutate(monthstodep = ifelse(is.na(monthstodep), 100, monthstodep)) %>%
ggplot(., aes(x = monthstocomposite, colour = mm_split)) +
  stat_ecdf() +
  theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  xlab("Time to URTI + Rhinitis (months)") +
  ylab("Cumulative proportion") +
  geom_hline(aes(yintercept = prop_composite), linetype = "dashed") +
  coord_cartesian(xlim = c(0,24)) + scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  facet_wrap(~mm_split, scales = "free_x")
ggsave("bipolar_epall_count_mm_to_depafter_inf_maniamixed_composite_negcontrol.png", bipolar_epall_count_mm_to_depafter_inf_composite.plot, height = 3, width = 7)

library(patchwork)
bipolar_epall_count_mm_to_depafter_inf_rti.plot / bipolar_epall_count_mm_to_depafter_inf_rhinitis.plot +
plot_annotation(tag_levels = 'a')
ggsave("bipolar_epall_count_mm_to_depafter_inf_maniamixed_negcontrol.png", height = 6, width = 7)

#===================================================================
# Negative control: Per-month rate
#===================================================================
#----------------------------------
# Baseline depression rate
#----------------------------------

# - get observation period - SAME as before in this section
observation_period <- getDf("observation_period")
mm_split <- getDf("mm_split")
bipolar_epall_mm_excl0 <- getDf("bipolar_epall_mm_excl0")
bipolar_epall_count_mm_to_deprecent <- getDf("bipolar_epall_count_mm_to_deprecent")
observation_period$person_id <- as.character(observation_period$person_id)
observation_period_bipolar_time_perind <- observation_period %>% 
  dplyr::select(person_id, observation_period_start_date, observation_period_end_date) %>% 
  inner_join(id_group %>% filter(group %in% c("PMMD", "noPMMD")), .) %>%
  mutate(observation_period_years = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days"))/365) %>%
  group_by(person_id) %>%
  summarise(observation_period_years_sum = sum(observation_period_years)) %>%
  left_join(., mm_split %>% dplyr::select("person_id", "mm_split")) %>% distinct

observation_period_bipolar <- observation_period %>% filter(person_id %in% bipolar_epall_mm_excl0$person_id)
mmep_denominator <- bipolar_epall_count_mm_to_depafter_inf %>% # have checked no overlap with bipolar_epall_mm_excl0
  left_join(., observation_period_bipolar) %>%
  filter(observation_period_start_date <= diagnosis_date_dep | is.na(diagnosis_date_dep)) %>%
  mutate(obs_monthssincemm = as.numeric(difftime(observation_period_end_date, diagnosis_date_end_mm, units = "days"))/30) %>%
  mutate(obs_monthssincemm = ifelse(obs_monthssincemm == 0, obs_monthssincemm+0.1, obs_monthssincemm)) %>% # include 0 in bin=1
  mutate(obs_monthssincemm_bin = ceiling(obs_monthssincemm)+1) %>% # round up and +1 so that drop-outs are excluded in next bracket (this is so that same-day-switches+dropout in the numerator make sense)
  filter(obs_monthssincemm_bin > 0) %>%
  group_by(obs_monthssincemm_bin, mm_split) %>%
  summarise(n=n()) %>% ungroup() %>%
  pivot_wider(., id_cols = obs_monthssincemm_bin, names_from = mm_split, values_from = n) %>%
  mutate(mania = ifelse(is.na(mania), 0, mania), mixed = ifelse(is.na(mixed), 0, mixed)) %>%
  mutate(mania_cumsum = cumsum(mania), mixed_cumsum = cumsum(mixed)) %>%
  mutate(mania_epleft = max(mania_cumsum) - mania_cumsum, mixed_epleft = max(mixed_cumsum) - mixed_cumsum)

month1_mm_denominator <- data.frame(obs_monthssincemm_bin = 1, mania_epleft = max(mmep_denominator$mania_cumsum), mixed_epleft = max(mmep_denominator$mixed_cumsum))

ttd_incidence_mm_denominator <- mmep_denominator %>% 
  dplyr::select(obs_monthssincemm_bin, mania_epleft, mixed_epleft) %>%
  full_join(month1_mm_denominator, .) %>%
  pivot_longer(cols = c("mania_epleft", "mixed_epleft"), names_to = "mm_split", values_to = "n_mm_ep") %>%
  mutate(mm_split = ifelse(mm_split == "mania_epleft", "mania", "mixed")) %>%
  rename(month = obs_monthssincemm_bin) %>%
  mutate(postmm_personyears = n_mm_ep/12)

# baseline - URTI
all_rti_ep_count <- bipolar_rti_dx %>% filter(rti == "rti") %>% nrow
all_mm_ep_count <- bipolar_epall_mm_excl0 %>% nrow
all_obstime_years <- sum(observation_period_bipolar_time_perind$observation_period_years_sum)

ttrti_month <- bipolar_epall_count_mm_to_depafter_inf_rti %>% 
  filter(months_mmend_to_rti >= 0) %>%
  mutate(month = ceiling(months_mmend_to_rti)) %>% 
  mutate(month = ifelse(month == 0, 1, month)) %>%
  group_by(month, mm_split) %>% summarise(n_rti_transition = n()) %>%
  inner_join(., ttd_incidence_mm_denominator) %>% mutate(rate = n_rti_transition/postmm_personyears) %>% 
  mutate(all_rti_ep_count = all_rti_ep_count, all_obstime_years = all_obstime_years) %>%
  mutate(base_rate = (all_rti_ep_count/all_obstime_years)) %>%
  mutate(vs_base = rate / base_rate) %>%
  mutate(rti_pc_permonth = n_rti_transition / all_rti_ep_count)
saveDf(ttrti_month, "ttrti_month", overwrite = T)

# - plot
ttrti_month.plot <- ggplot(ttrti_month, aes(x = month, y = vs_base, colour = mm_split)) +
  geom_point() +
  theme_bw() + scale_colour_manual(values = millerstone[c(7,4)]) +
  theme(legend.position = "none") +
  xlab("Month after mania/mixed episode") +
  ylab("Monthly / Overall URTI rate") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_x_continuous(limits = c(1, 12), breaks = seq(1, 12, by = 1)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
  facet_wrap(~ mm_split, ncol = 2)

# baseline - RHINITIS
all_rhinitis_ep_count <- bipolar_rhinitis_dx %>% filter(rhinitis == "rhinitis") %>% nrow
all_mm_ep_count <- bipolar_epall_mm_excl0 %>% nrow
all_obstime_years <- sum(observation_period_bipolar_time_perind$observation_period_years_sum)

ttrhinitis_month <- bipolar_epall_count_mm_to_depafter_inf_rhinitis %>% 
  filter(months_mmend_to_rhinitis >= 0) %>%
  mutate(month = ceiling(months_mmend_to_rhinitis)) %>% 
  mutate(month = ifelse(month == 0, 1, month)) %>%
  group_by(month, mm_split) %>% summarise(n_rhinitis_transition = n()) %>%
  inner_join(., ttd_incidence_mm_denominator) %>% mutate(rate = n_rhinitis_transition/postmm_personyears) %>% 
  mutate(all_rhinitis_ep_count = all_rhinitis_ep_count, all_obstime_years = all_obstime_years) %>%
  mutate(base_rate = (all_rhinitis_ep_count/all_obstime_years)) %>%
  mutate(vs_base = rate / base_rate) %>%
  mutate(rhinitis_pc_permonth = n_rhinitis_transition / all_rhinitis_ep_count)
saveDf(ttrhinitis_month, "ttrhinitis_month", overwrite = T)

# - plot
ttrhinitis_month.plot <- ggplot(ttrhinitis_month, aes(x = month, y = vs_base, colour = mm_split)) +
  geom_point() +
  theme_bw() + scale_colour_manual(values = millerstone[c(7,4)]) +
  theme(legend.position = "none") +
  xlab("Month after mania/mixed episode") +
  ylab("Monthly / Overall rhinitis rate") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_x_continuous(limits = c(1, 12), breaks = seq(1, 12, by = 1)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
  facet_wrap(~ mm_split, ncol = 2)

library(patchwork)
ttrti_month.plot / ttrhinitis_month.plot +
plot_annotation(tag_levels = 'a')
ggsave("ttnegcontrol_month.png", height = 6, width = 7)

# baseline - COMPOSITE
all_composite_ep_count <- bipolar_composite_dx %>% filter(composite == "composite") %>% nrow
all_mm_ep_count <- bipolar_epall_mm_excl0 %>% nrow
all_obstime_years <- sum(observation_period_bipolar_time_perind$observation_period_years_sum)

ttcomposite_month <- bipolar_epall_count_mm_to_depafter_inf_composite %>% 
  filter(months_mmend_to_composite >= 0) %>%
  mutate(month = ceiling(months_mmend_to_composite)) %>% 
  mutate(month = ifelse(month == 0, 1, month)) %>%
  group_by(month, mm_split) %>% summarise(n_composite_transition = n()) %>%
  inner_join(., ttd_incidence_mm_denominator) %>% mutate(rate = n_composite_transition/postmm_personyears) %>% 
  mutate(all_composite_ep_count = all_composite_ep_count, all_obstime_years = all_obstime_years) %>%
  mutate(base_rate = (all_composite_ep_count/all_obstime_years)) %>%
  mutate(vs_base = rate / base_rate) %>%
  mutate(composite_pc_permonth = n_composite_transition / all_composite_ep_count)
saveDf(ttcomposite_month, "ttcomposite_month", overwrite = T)

# - plot
ttcomposite_month.plot <- ggplot(ttcomposite_month, aes(x = month, y = vs_base, colour = mm_split)) +
  geom_point() +
  theme_bw() + scale_colour_manual(values = millerstone[c(7,4)]) +
  theme(legend.position = "none") +
  xlab("Month after mania/mixed episode") +
  ylab("Monthly / Overall URTI + rhinitis rate") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_x_continuous(limits = c(1, 12), breaks = seq(1, 12, by = 1)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
  facet_wrap(~ mm_split, ncol = 2)
ggsave("ttcomposite_month.plot.png", ttcomposite_month.plot, height = 3, width = 7)

#===================================================================
# Clinical features analysis
#===================================================================

bipolar_mm_clinical_cgis <- getDf("bipolar_mm_clinical_cgis")

bipolar_rti_1m <- bipolar_epall_count_mm_to_depafter_inf_rti %>% 
  mutate(rti1m = ifelse(months_mmend_to_rti <= 1, 1, 0)) %>%
  mutate(rti1m = ifelse(is.na(months_mmend_to_rti), 0, rti1m)) %>%
  rename(diagnosis_date = diagnosis_date_mm, diagnosis_date_end = diagnosis_date_end_mm) %>%
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, rti1m, months_mmend_to_rti) %>%
  distinct()

bipolar_rhinitis_1m <- bipolar_epall_count_mm_to_depafter_inf_rhinitis %>% 
  mutate(rhinitis1m = ifelse(months_mmend_to_rhinitis <= 1, 1, 0)) %>%
  mutate(rhinitis1m = ifelse(is.na(months_mmend_to_rhinitis), 0, rhinitis1m)) %>%
  rename(diagnosis_date = diagnosis_date_mm, diagnosis_date_end = diagnosis_date_end_mm) %>%
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, rhinitis1m, months_mmend_to_rhinitis) %>%
  distinct()

bipolar_mm_clinical_negcontrol <- bipolar_mm_clinical %>%
  left_join(., bipolar_rti_1m) %>% left_join(., bipolar_rhinitis_1m)

bipolar_mm_clinical_cgis_negcontrol <- bipolar_mm_clinical_cgis %>%
  left_join(., bipolar_rti_1m) %>% left_join(., bipolar_rhinitis_1m) 

# - COV with URTI negative control
bipolar_mm_icdseverity_cov_glmrr_rti <- bipolar_mm_clinical_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rti1m ~ icdseverity + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_hospital_cov_glmrr_rti <- bipolar_mm_clinical_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rti1m ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_cgis_cov_glmrr_rti <- bipolar_mm_clinical_cgis_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rti1m ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_prevMDS_cov_glmrr_rti <- bipolar_mm_clinical_negcontrol %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rti1m ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_glmrr_cov_rti <- rbind(bipolar_mm_clinical_prevMDS_cov_glmrr_rti,
  bipolar_mm_cgis_cov_glmrr_rti,
  bipolar_mm_hospital_cov_glmrr_rti) %>%
#  bipolar_mm_icdseverity_cov_glmrr) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(p_adj = p.adjust(p.value, "bonf")) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "URTI")
saveDf(bipolar_mm_clinical_glmrr_cov_rti, "bipolar_mm_clinical_glmrr_cov_rti", overwrite = T)

# - COV with allergic rhinitis negative control
bipolar_mm_icdseverity_cov_glmrr_rhinitis <- bipolar_mm_clinical_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rhinitis1m ~ icdseverity + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_hospital_cov_glmrr_rhinitis <- bipolar_mm_clinical_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rhinitis1m ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_cgis_cov_glmrr_rhinitis <- bipolar_mm_clinical_cgis_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rhinitis1m ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_prevMDS_cov_glmrr_rhinitis <- bipolar_mm_clinical_negcontrol %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(rhinitis1m ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_glmrr_cov_rhinitis <- rbind(bipolar_mm_clinical_prevMDS_cov_glmrr_rhinitis,
  bipolar_mm_cgis_cov_glmrr_rhinitis,
  bipolar_mm_hospital_cov_glmrr_rhinitis) %>%
#  bipolar_mm_icdseverity_cov_glmrr) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(p_adj = p.adjust(p.value, "bonf")) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "rhinitis")
saveDf(bipolar_mm_clinical_glmrr_cov_rhinitis, "bipolar_mm_clinical_glmrr_cov_rhinitis", overwrite = T)

# Plot
bipolar_mm_clinical_glmrr_plot_negcontrol <- rbind(
  bipolar_mm_clinical_glmrr_cov_rti,
  bipolar_mm_clinical_glmrr_cov_rhinitis)
saveDf(bipolar_mm_clinical_glmrr_plot_negcontrol, "bipolar_mm_clinical_glmrr_plot_negcontrol", overwrite = T)

bipolar_mm_clinical_glmrr_plot_negcontrol.tmp <- bipolar_mm_clinical_glmrr_plot_negcontrol %>%
  mutate(term = case_when(term == "severitycgis_max_d07" ~ "severity: CGI-S score +1",
    term == "prev_MDSprevMDS" ~ "prev short MDTT (<1m)",
    #term == "mm_splitmixed" ~ "mixed ep (vs mania)",
    term == "hospital01" ~ "hospital for mania/mixed",
    .default = term))

# For final figure
bipolar_mm_clinical_glmrr_plot_covfinal_negcontrol <- bipolar_mm_clinical_glmrr_plot_negcontrol.tmp %>%
  mutate(title = "Risk Ratio of Switch in 1m") %>%
  mutate(term = ifelse(term == "previous Dep switch (1m)", "prev short MDTT (<1m)", term)) %>%
  mutate(term2 = factor(term, levels = rev(c("prev short MDTT (<1m)", "severity: CGI-S score +1", "hospital for mania/mixed", "antidepressant", "antiepileptic", "antipsychotic", "benzodiazepine", "lithium", "antihistamine")))) %>%
  ggplot(., aes(x = term2, y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR), shape = 17) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_label_repel(size = 2.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  scale_y_log10(limits = c(0.1, 7)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate, full bipolar disorder sample (adj. covariates)") +
  coord_flip() +
  facet_grid(facet ~ mm_split, scales = "free", space = "free_y")
ggsave("bipolar_mm_clinical_glmrr_plot_covfinal_negcontrol.png", bipolar_mm_clinical_glmrr_plot_covfinal_negcontrol, width = 7, height = 3, dpi = 600)

# COMPOSITE
bipolar_composite_1m <- bipolar_epall_count_mm_to_depafter_inf_composite %>% 
  mutate(composite1m = ifelse(months_mmend_to_composite <= 1, 1, 0)) %>%
  mutate(composite1m = ifelse(is.na(months_mmend_to_composite), 0, composite1m)) %>%
  rename(diagnosis_date = diagnosis_date_mm, diagnosis_date_end = diagnosis_date_end_mm) %>%
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, composite1m, months_mmend_to_composite) %>%
  distinct()

bipolar_mm_clinical_negcontrol <- bipolar_mm_clinical %>%
  left_join(., bipolar_composite_1m)

bipolar_mm_clinical_cgis_negcontrol <- bipolar_mm_clinical_cgis %>%
  left_join(., bipolar_composite_1m)

# - COV with Composite negative control
bipolar_mm_icdseverity_cov_glmrr_composite <- bipolar_mm_clinical_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(composite1m ~ icdseverity + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_hospital_cov_glmrr_composite <- bipolar_mm_clinical_negcontrol %>%
  filter(hospital_include == 1) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(composite1m ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_cgis_cov_glmrr_composite <- bipolar_mm_clinical_cgis_negcontrol %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(composite1m ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_prevMDS_cov_glmrr_composite <- bipolar_mm_clinical_negcontrol %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(composite1m ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .)), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_glmrr_cov_composite <- rbind(bipolar_mm_clinical_prevMDS_cov_glmrr_composite,
  bipolar_mm_cgis_cov_glmrr_composite,
  bipolar_mm_hospital_cov_glmrr_composite) %>%
#  bipolar_mm_icdseverity_cov_glmrr) %>% 
  dplyr::select(-c(data, model)) %>%
  mutate(estimate_RR = exp(estimate),
        lowerRR = exp(estimate - (std.error * 1.96)),
        higherRR = exp(estimate + (std.error * 1.96))) %>%
  mutate(lab = ifelse(p.value < 0.05, paste(sprintf('%.2f', estimate_RR), 
    " [", sprintf('%.2f', lowerRR), "-", sprintf('%.2f', higherRR), "], ",
    "p=", formatC(p.value, format = "e", digits = 2), sep = ""), NA)) %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(p_adj = p.adjust(p.value, "bonf")) %>%
  mutate(sig = case_when(p_adj < 0.05 ~ millerstone[8], 
                         p.value <= 0.05 ~ millerstone[2],
                         p.value > 0.05 ~ "gray80")) %>%
  mutate(facet = "URTI + Rhinitis")
saveDf(bipolar_mm_clinical_glmrr_cov_composite, "bipolar_mm_clinical_glmrr_cov_composite", overwrite = T)

# Plot COMPOSITE
bipolar_mm_clinical_glmrr_cov_composite.tmp <- bipolar_mm_clinical_glmrr_cov_composite %>%
  mutate(term = case_when(term == "severitycgis_max_d07" ~ "severity: CGI-S score +1",
    term == "prev_MDSprevMDS" ~ "prev short MDTT (<1m)",
    #term == "mm_splitmixed" ~ "mixed ep (vs mania)",
    term == "hospital01" ~ "hospital for mania/mixed",
    .default = term))

# For final figure
bipolar_mm_clinical_glmrr_cov_composite_plot <- bipolar_mm_clinical_glmrr_cov_composite.tmp %>%
  mutate(title = "Risk Ratio of Switch in 1m") %>%
  mutate(term = ifelse(term == "previous Dep switch (1m)", "prev short MDTT (<1m)", term)) %>%
  mutate(term2 = factor(term, levels = rev(c("prev short MDTT (<1m)", "severity: CGI-S score +1", "hospital for mania/mixed", "antidepressant", "antiepileptic", "antipsychotic", "benzodiazepine", "lithium", "antihistamine")))) %>%
  ggplot(., aes(x = term2, y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR), shape = 17) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_label_repel(size = 2.5, nudge_x = 0.15, force = 5, colour = "black", min.segment.length = 0) +
  scale_colour_identity() +
  scale_y_log10(limits = c(0.1, 7)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("Variable") +
  ylab("RR estimate, full bipolar disorder sample (adj. covariates)") +
  coord_flip() +
  facet_grid(facet ~ mm_split, scales = "free", space = "free_y")
ggsave("bipolar_mm_clinical_glmrr_cov_composite_plot.png", bipolar_mm_clinical_glmrr_cov_composite_plot, width = 7, height = 3, dpi = 600)

#===================================================================
# Medication analysis
#===================================================================

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep <- getDf("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep")
days_since_hosp <- getDf("days_since_hosp")

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_negcontrol.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
  left_join(., bipolar_rti_1m) %>%
  left_join(., bipolar_rhinitis_1m)
  
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rti <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_negcontrol.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>%  
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(rti1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
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
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rti, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rti", overwrite = T)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rhinitis <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_negcontrol.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>%  
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(rhinitis1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
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
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rhinitis, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rhinitis")

# Plot
# - JOINT, not adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rti.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rti %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(term, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  scale_y_log10() + coord_flip(ylim = c(0.1, 7)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for URTI in 1m (joint model, adj. for covariates)") +
  coord_flip() +
  facet_grid(name ~ mm_split, scales = "free_x")

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rhinitis.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rhinitis %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(term, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  scale_y_log10() + coord_flip(ylim = c(0.1, 7)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Allergic Rhinitis in 1m (joint model, adj. for covariates)") +
  coord_flip() +
  facet_grid(name ~ mm_split, scales = "free_x")

library(patchwork)
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rti.plot / bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_rhinitis.plot +
plot_annotation(tag_levels = 'a')
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_negcontrol.png", height = 6, width = 7)

# CHECK
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_negcontrol.prep %>% 
  filter(!is.na(cellulitis1m)) %>% group_by(lithium, cellulitis1m, mm_split) %>% summarise(n=n())

# COMPOSITE OUTCOME
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_negcontrol.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
  left_join(., bipolar_composite_1m) %>%
  left_join(., days_since_hosp)
  
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_composite <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide_negcontrol.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>%  
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(glm(composite1m ~ hospital_los + days_since_hosp + lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id), tidied = map(model, tidy)), 
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
saveDf(bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_composite, "bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_composite", overwrite = T)

# Plot
# - JOINT, not adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_composite.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_composite %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(term, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  scale_y_log10() + coord_flip(ylim = c(0.1, 7)) +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for URTI + Rhinitis in 1m (joint model, adj. for covariates)") +
  coord_flip() +
  facet_grid(name ~ mm_split, scales = "free_x")

ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_negcontrol_composite.png", height = 3, width = 7)
