#==========================================
# Mediation analysis
#==========================================

bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep <- getDf("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep")
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) #%>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) # due to collinearity issue 

# Overall effect
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mania_prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  mutate(antidepressant = ifelse(mm_split == "mania" & antidepressant == 1, 0, antidepressant)) %>% 
  left_join(., days_since_hosp) %>%  
  filter(mm_split == "mania") %>%
  mutate(log_days_since_hosp = log1p(days_since_hosp))
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mixed_prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide.prep %>%
  left_join(., days_since_hosp) %>%  
  filter(mm_split == "mixed") %>%
  mutate(log_days_since_hosp = log1p(days_since_hosp))

# No follow-up adjustment
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mania <- glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + 
  hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, 
  family='binomial', 
#  family=poisson(link='log'), 
  data = bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mania_prep)#, 
#  vcov = vcovCL, cluster = ~person_id)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mixed <- glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + 
  hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_edu2 + cov_employ2 + cov_race2, 
  family='binomial', 
#  family=poisson(link='log'), 
  data = bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mixed_prep)#, 

# Effect of medications on days_since_hosp
days_since_hosp_outcome_mania <- lm(log_days_since_hosp ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + 
  hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_rel + cov_edu2 + cov_employ2 + cov_race2, 
  #family = quasipoisson(link = "log"),
  data = bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mania_prep)#, 
#  vcov = vcovCL, cluster = ~person_id)

days_since_hosp_outcome_mixed <- lm(log_days_since_hosp ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + 
  hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_rel + cov_edu2 + cov_employ2 + cov_race2, 
  #family = quasipoisson(link = "log"),
  data = bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mixed_prep)#, 

# Effect of days_since_hosp on MDS1m
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mania <- glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + 
  log_days_since_hosp + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_rel + cov_edu2 + cov_employ2 + cov_race2, 
  family='binomial', 
  data = bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mania_prep)#, 
#  vcov = vcovCL, cluster = ~person_id)

bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mixed <- glm(MDS1m ~ lithium + antiepileptic + antipsychotic + antidepressant + antihistamine + benzodiazepine + 
  log_days_since_hosp + hospital_los + cov_age_ep + cov_years_ehr + cov_gender + cov_rel + cov_edu2 + cov_employ2 + cov_race2, 
  family='binomial', 
  data = bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_mixed_prep)#, 

source('mediate.r')
tidy_mediate <- function(x, average_only = FALSE) {
  s <- summary(x)
  out <- data.frame(
    effect = c("ACME (control)", "ACME (treated)", "ADE (control)", "ADE (treated)", "Total Effect", "Prop. Mediated (control)", "Prop. Mediated (treated)", "ACME (average)", "ADE (average)", "Prop. Mediated (average)"),
    estimate = c(s$d0, s$d1, s$z0, s$z1, s$tau.coef, s$n0, s$n1, s$d.avg, s$z.avg, s$n.avg),
    conf.low = c(s$d0.ci[1], s$d1.ci[1], s$z0.ci[1], s$z1.ci[1], s$tau.ci[1], s$n0.ci[1], s$n1.ci[1], s$d.avg.ci[1], s$z.avg.ci[1], s$n.avg.ci[1]),
    conf.high = c(s$d0.ci[2], s$d1.ci[2], s$z0.ci[2], s$z1.ci[2], s$tau.ci[2], s$n0.ci[2], s$n1.ci[2], s$d.avg.ci[2], s$z.avg.ci[2], s$n.avg.ci[2]),
    p.value = c(s$d0.p, s$d1.p, s$z0.p, s$z1.p, s$tau.p, s$n0.p, s$n1.p, s$d.avg.p, s$z.avg.p, s$n.avg.p)
  )
  if (average_only) out <- out[c(8, 9, 5, 10), ]
  rownames(out) <- NULL
  out
}

# limitations of analysis: no coeftest (so no clustered standard errors) - so just doing Poisson regression essentially. 
results_lithium_mania <- mediate(days_since_hosp_outcome_mania, 
  bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mania, 
  treat = 'lithium', mediator = 'log_days_since_hosp', boot = TRUE, sims = 500)
results_lithium_mania_tidy <- tidy_mediate(results_lithium_mania, average_only = TRUE) %>% mutate(mm_split = "mania", drug = "lithium")

results_antiepileptic_mania <- mediate(days_since_hosp_outcome_mania, 
  bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mania, 
  treat = 'antiepileptic', mediator = 'log_days_since_hosp', boot = TRUE, sims = 500)
results_antiepileptic_mania_tidy <- tidy_mediate(results_antiepileptic_mania, average_only = TRUE) %>% mutate(mm_split = "mania", drug = "antiepileptic")

results_antipsychotic_mania <- mediate(days_since_hosp_outcome_mania, 
  bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mania, 
  treat = 'antipsychotic', mediator = 'log_days_since_hosp', boot = TRUE, sims = 500)
results_antipsychotic_mania_tidy <- tidy_mediate(results_antipsychotic_mania, average_only = TRUE) %>% mutate(mm_split = "mania", drug = "antipsychotic")

# limitations of analysis: no coeftest (so no clustered standard errors) - so just doing Poisson regression essentially. 
results_lithium_mixed <- mediate(days_since_hosp_outcome_mixed, 
  bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mixed, 
  treat = 'lithium', mediator = 'log_days_since_hosp', boot = TRUE, sims = 500)
results_lithium_mixed_tidy <- tidy_mediate(results_lithium_mixed, average_only = TRUE) %>% mutate(mm_split = "mixed", drug = "lithium")

results_antiepileptic_mixed <- mediate(days_since_hosp_outcome_mixed, 
  bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mixed, 
  treat = 'antiepileptic', mediator = 'log_days_since_hosp', boot = TRUE, sims = 500)
results_antiepileptic_mixed_tidy <- tidy_mediate(results_antiepileptic_mixed, average_only = TRUE) %>% mutate(mm_split = "mixed", drug = "antiepileptic")

results_antipsychotic_mixed <- mediate(days_since_hosp_outcome_mixed, 
  bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_sens_followup_mixed, 
  treat = 'antipsychotic', mediator = 'log_days_since_hosp', boot = TRUE, sims = 500)
results_antipsychotic_mixed_tidy <- tidy_mediate(results_antipsychotic_mixed, average_only = TRUE) %>% mutate(mm_split = "mixed", drug = "antipsychotic")

mediation_results <- rbind(results_lithium_mania_tidy, results_antiepileptic_mania_tidy, results_antipsychotic_mania_tidy, 
                           results_lithium_mixed_tidy, results_antiepileptic_mixed_tidy, results_antipsychotic_mixed_tidy) %>%
  dplyr::select(mm_split, drug, effect, estimate, conf.low, conf.high, p.value)

saveDf(mediation_results, "mediation_results", overwrite = T)

#==========================================
# Availability of EHR data across outcomes per state
#==========================================

bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")

meds_psych_bipolar_dx <- getDf("meds_psych_bipolar_dx")
meds_bipolar <- meds_psych_bipolar_dx %>% filter(person_id %in% bipolar_epall_count_mm_to_depafter_inf$person_id) %>% dplyr::select(person_id) %>% distinct() %>% mutate(meds = 1)

pmmdgroup_cgis <- getDf("pmmdgroup_cgis")
nopmmdgroup_cgis <- getDf("nopmmdgroup_cgis")
cgis_mania <- rbind(pmmdgroup_cgis, nopmmdgroup_cgis)
cgis_mania$person_id <- as.character(cgis_mania$person_id)
cgis_bipolar <- cgis_mania %>% filter(person_id %in% bipolar_epall_count_mm_to_depafter_inf$person_id) %>% dplyr::select(person_id) %>% distinct() %>% mutate(cgis = 1)

phq9_mania <- getDf("phq9_mm_sql")
phq9_mania$person_id <- as.character(phq9_mania$person_id)
phq9_bipolar <- phq9_mania %>% filter(person_id %in% bipolar_epall_count_mm_to_depafter_inf$person_id) %>% dplyr::select(person_id) %>% distinct() %>% mutate(phq9 = 1)

madrscgis_mania <- getDf("madrs_mm_sql")
madrscgis_mania$person_id <- as.character(madrscgis_mania$person_id)
madrscgis_bipolar <- madrscgis_mania %>% filter(person_id %in% bipolar_epall_count_mm_to_depafter_inf$person_id) %>% dplyr::select(person_id) %>% distinct() %>% mutate(madrs = 1)

pmmdgroup_visit <- getDf("pmmdgroup_visit")
nopmmdgroup_visit <- getDf("nopmmdgroup_visit")
bipolar_visit <- full_join(pmmdgroup_visit, nopmmdgroup_visit)
bipolar_visit$person_id <- as.character(bipolar_visit$person_id)
visit_bipolar <- bipolar_visit %>% filter(person_id %in% bipolar_epall_count_mm_to_depafter_inf$person_id) %>% dplyr::select(person_id) %>% distinct() %>% mutate(visit = 1)

demog_cov <- getDf("demog_cov")

demog_cov_state_outcomes <- demog_cov %>%
  mutate(cov_state = ifelse(cov_state == "No matching concept", "N/A", cov_state)) %>%
  dplyr::select(person_id, cov_state, cov_race) %>%
  filter(person_id %in% bipolar_epall_count_mm_to_depafter_inf$person_id) %>%
  left_join(., meds_bipolar) %>%
  left_join(., phq9_bipolar) %>%
  left_join(., cgis_bipolar) %>%
  left_join(., madrscgis_bipolar) %>%
  left_join(., visit_bipolar) %>%
  replace(., is.na(.), 0)

state_summary <- demog_cov_state_outcomes %>% 
  dplyr::select(-cov_race) %>%
  pivot_longer(!c("person_id", "cov_state"), names_to = "data", values_to = "available") %>%
  group_by(cov_state, data) %>%
  summarise(prop_dataavailable = sum(available)/n(), total_n = n()) 

state_race_summary <- demog_cov_state_outcomes %>%
  dplyr::select(person_id, cov_state, cov_race) %>%
  pivot_longer(!c("person_id", "cov_state"), names_to = "var", values_to = "race") %>%
  group_by(cov_state, race) %>%
  summarise(n=n()) %>% ungroup %>%
  group_by(cov_state) %>% mutate(total_perstate = sum(n)) %>% ungroup %>%
  mutate(race_prop = n/total_perstate)

state_n <- state_summary %>% distinct(cov_state, total_n)

# Order states by decreasing sample size
state_order <- state_n %>%
  arrange(desc(total_n)) %>%
  pull(cov_state)

# Get one total_n per state
state_n <- state_summary %>%
  distinct(cov_state, total_n)

# Order states by decreasing sample size
state_order <- state_n %>%
  arrange(desc(total_n)) %>%
  pull(cov_state)

data_cols <- c(
  "Medications" = millerstone[2],
  "PHQ-9"       = millerstone[4],
  "CGI-S"       = millerstone[6],
  "MADRS (derived)"       = millerstone[8],
  "Visits"       = millerstone[10]
)

# Recode data labels and set facet order
state_summary <- state_summary %>%
  mutate(
    cov_state = factor(cov_state, levels = state_order),
    data = recode(
      data,
      "meds"  = "Medications",
      "phq9"  = "PHQ-9",
      "cgis"  = "CGI-S",
      "madrs" = "MADRS (derived)",
      "visit" = "Visits"
    ),
    data = factor(
      data,
      levels = c("Visits", "Medications", "PHQ-9", "CGI-S", "MADRS (derived)")
    )
  )

state_n <- state_n %>%
  mutate(
    cov_state = factor(cov_state, levels = state_order)
  )

# - add race barplot
# Make state ordering match existing plots
state_race_summary <- state_race_summary %>%
  mutate(
    cov_state = factor(cov_state, levels = state_order)
  )

# Proportion of data available
p1 <- ggplot(
  state_summary,
  aes(x = cov_state, y = prop_dataavailable, fill = data)) +
  geom_col(width = 0.7) +
  facet_grid(data ~ .) +
  scale_fill_manual(values = data_cols) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent) +
  labs(x = NULL, y = "% of individuals with data available") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(),
    legend.position = "none")

# Race composition panel
p_race <- ggplot(
  state_race_summary,
  aes(x = cov_state, y = race_prop, fill = race)
) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Race"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.position = "bottom"                  
  )

# Sample size, normal linear scale
p2 <- ggplot(state_n, aes(x = cov_state, y = total_n)) +
  geom_col(width = 0.7, fill = millerstone[12]) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))) +
  labs(x = "State", y = "Sample size (n)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5),
    panel.grid.major.x = element_blank())


# Combine
library(patchwork)
p1 / p_race / p2 + plot_layout(heights = c(5, 1, 1))

ggsave("perstate_dataavailability.png", height = 10, width = 7)

#==========================================
# Confirming depression "provenance"
# Analysis checking antidepressant for anxiety
#==========================================

# How many total depressive episodes in the cohort?
bipolar_epall <- getDf("bipolar_epall")
source('diagnosis_definitions.r')

table(bipolar_epall$state) # 315526 dep episodes

# How many within 6m of mania?
bipolar_epall_count_mm_to_depafter_inf %>% filter(mm_to_depafter <= 6 & mm_to_depafter >= 0) %>% dim

# How many are due to the antidepressant definition?
# - get the updated mm_episode/dep_episode numbering file
mania_dep_dates_toepisodes_qc_broaddef <- getDf("mania_dep_dates_toepisodes_qc_broaddef")
# - add it to the df with the supersede/addpmmd info
mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest <- getDf("mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest")
mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest_update <- mania_dep_dates_toepisodes_qc_broaddef %>%
  dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, mm_episode, dep_episode) %>%
  inner_join(., mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest %>% dplyr::select(-c("mm_episode", "dep_episode")))

checkdefinition_dep <- bipolar_epall_count_mm_to_depafter_inf %>% # start with what is used for Fig 1
  filter(mm_to_depafter >= 0) %>%
  filter(mm_to_depafter <= 6) %>% # take dep under 6 months
  dplyr::select(person_id, mm_episode) %>%
  inner_join(., mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest_update) # join with df indicating where dep definition came from

addpmmd_n_useddata <- checkdefinition_dep %>%  # previously addpmmd_n
  group_by(supersede, addpmmd, method) %>%
  summarise(n=n()) %>%
  mutate(supersede = ifelse(is.na(supersede), "no", supersede)) %>%
  mutate(addpmmd = ifelse(is.na(addpmmd), "no", addpmmd)) %>%
  rename(dep_detection_method = method)
saveDf(addpmmd_n_useddata, "addpmmd_n_useddata", overwrite = T)

# How many of the antidepressant defined group have anxiety diagnosis within the previous year
icdCodes_anxiety <- list(anxiety=anxiety)
allPatients_anxiety <- getPatients(icdCodes_anxiety) # n=324049 > 346932
allPatients_anxiety_bipolar <- allPatients_anxiety %>% filter(person_id %in% (id_group %>% filter(group != "MMonly") %>% pull(person_id)))
anxiety_bipolar_diagnosis <- getDiagnosis(allPatients_anxiety_bipolar$person_id, icd_codes = c(anxiety)) %>% distinct(diagnosis_date, icd_code, .keep_all = TRUE) %>%
  mutate(person_id = as.character(person_id))
saveDf(anxiety_bipolar_diagnosis, "anxiety_bipolar_diagnosis")

checkdefinition_dep_anxiety <- checkdefinition_dep %>%
  left_join(., anxiety_bipolar_diagnosis %>% rename(diagnosis_date_anxiety = diagnosis_date) %>% dplyr::select(person_id, diagnosis_date_anxiety) %>% distinct, relationship = "many-to-many") %>%
  mutate(months_anxiety_before_dep = as.numeric(difftime(dep_date, diagnosis_date_anxiety, unit = "days"))/30) %>% # +ive means depression comes after anxiety
  mutate(months_anxiety_before_dep = ifelse(months_anxiety_before_dep < 0, NA, months_anxiety_before_dep)) %>% # ignore where anxiety comes after depression
  arrange(months_anxiety_before_dep) %>% # as positive values here, put the shortest time first
  group_by(person_id, mm_episode, dep_episode, state, diagnosis_date, diagnosis_date_end, dep_date, addpmmd, supersede, method) %>%
  summarise(months_anxiety_before_dep = min(months_anxiety_before_dep)) %>%
  mutate(anxiety_within1yprev = ifelse(months_anxiety_before_dep <= 12, "yes_anxiety1y", "noevidence_anxiety1y")) %>%
  mutate(anxiety_within1yprev = ifelse(is.na(anxiety_within1yprev), "noevidence_anxiety1y", anxiety_within1yprev))
saveDf(checkdefinition_dep_anxiety, "checkdefinition_dep_anxiety")

# - look at how many per depression definition group had anxiety within the past year
checkdefinition_dep_anxiety_summary <- checkdefinition_dep_anxiety %>% group_by(method, anxiety_within1yprev) %>% 
  summarise(n=n()) %>% 
  pivot_wider(., names_from = "anxiety_within1yprev", values_from = "n") %>%
  mutate(prop_yes_anxiety1y = yes_anxiety1y/(yes_anxiety1y+noevidence_anxiety1y)) %>%
  rename(dep_detection_method = method)
saveDf(checkdefinition_dep_anxiety_summary, "checkdefinition_dep_anxiety_summary", overwrite = T)

# - episodes to exclude
exclude_noadanxiety <- checkdefinition_dep_anxiety %>%
  filter(method == "ad" & anxiety_within1yprev == "yes_anxiety1y") %>%
  dplyr::select(person_id, mm_episode, state, diagnosis_date, diagnosis_date_end)

# Clinical features result
bipolar_mm_clinical <- getDf("bipolar_mm_clinical")
bipolar_mm_clinical_cgis <- getDf("bipolar_mm_clinical_cgis")
bipolar_mm_clinical_noadanxiety <- bipolar_mm_clinical %>% anti_join(., exclude_noadanxiety)
bipolar_mm_clinical_cgis_noadanxiety <- bipolar_mm_clinical_cgis %>% anti_join(., exclude_noadanxiety)

# - COV
bipolar_mm_hospital_cov_glmrr_noadanxiety <- bipolar_mm_clinical_noadanxiety %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ hospital01 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_cgis_cov_glmrr_noadanxiety <- bipolar_mm_clinical_cgis_noadanxiety %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ severitycgis_max_d07 + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_prevMDS_cov_glmrr_noadanxiety <- bipolar_mm_clinical_noadanxiety %>%
  filter(mm_episode >= 2) %>%
  nest(.by = mm_split) %>%
  mutate(model = map(data, ~coeftest(glm(MDS1m ~ prev_MDS + cov_age_ep + cov_age_ehr + cov_years_ehr + cov_gender + cov_state2 + cov_rel + cov_edu2 + cov_employ2 + cov_race2, family=poisson(link='log'), data = .), vcov = vcovCL, cluster = ~person_id)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(analysis = "cov")

bipolar_mm_clinical_glmrr_cov_noadanxiety <- rbind(bipolar_mm_clinical_prevMDS_cov_glmrr_noadanxiety %>% mutate(primary_var = "prevMDS"),
  bipolar_mm_cgis_cov_glmrr_noadanxiety %>% mutate(primary_var = "CGI-S severity"),
  bipolar_mm_hospital_cov_glmrr_noadanxiety %>% mutate(primary_var = "hospital01")) %>% 
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

# - plot
bipolar_mm_clinical_glmrr_cov_noadanxiety.tmp <- bipolar_mm_clinical_glmrr_cov_noadanxiety %>%
  filter(!grepl("cov|Intercept", term)) %>%
  mutate(term = case_when(term == "severitycgis_max_d07" ~ "severity: CGI-S score +1",
    term == "prev_MDSprevMDS" ~ "prev short MDTT (<1m)",
    #term == "mm_splitmixed" ~ "mixed ep (vs mania)",
    term == "hospital01" ~ "hospital for mania/mixed",
    .default = term))

bipolar_mm_clinical_glmrr_cov_noadanxiety.plot <- bipolar_mm_clinical_glmrr_cov_noadanxiety.tmp %>%
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
ggsave("bipolar_mm_clinical_glmrr_cov_noad.png", bipolar_mm_clinical_glmrr_cov_noadanxiety.plot, width = 7, height = 3, dpi = 600)

# Medications analysis
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep <- getDf("bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep")
bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide_noadanxiety.prep <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_csh_180d_wide.prep %>%
  mutate(MDS1m = ifelse(end_csh <= 30, 1, 0)) %>% mutate(MDS1m = ifelse(is.na(MDS1m), 0, MDS1m)) %>%
  #mutate(cov_state = ifelse(cov_state %in% c("TN", "PA", "VA", "IL", "FL", "SC", "MA", "OR", "WV"), "grouped", cov_state)) %>% # due to collinearity issue 
  anti_join(., exclude_noadanxiety)

# - JOINT, not adjusting for follow-up
days_since_hosp <- getDf("days_since_hosp")
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_noadanxiety <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_cov_glmrr_180d_wide_noadanxiety.prep %>%
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

# - JOINT, not adjusting for follow-up
bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_noadanxiety.plot <- bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_noadanxiety %>%
  mutate(name = "broad class") %>%
  filter(term %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine")) %>%
  ggplot(aes(x = fct_reorder(term, estimate_RR), y = estimate_RR, colour = sig, label = lab)) +
  geom_pointrange(aes(ymin = lowerRR, ymax = higherRR)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", colour = "navy") +
  geom_text_repel(size = 3, nudge_x = 0.20, force = 10, colour = "black", min.segment.length = 0.5) +
  scale_colour_identity() +
  theme_bw() + theme(legend.position = "none") +
  xlab("") +
  ylab("RR estimate for Dep switch in 1m (joint model, antidepressant sensitivity)") +
  coord_flip() +
  facet_grid(name ~ mm_split, scales = "free_x")
ggsave("bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_noadanxiety.png", bipolar_dx7d_hosp_medsbroadclass_mmsplit_joint_glmrr_noadanxiety.plot, width = 7, height = 3, dpi = 600)


#==========================================
# Mixed features - are they depression with mixed features vs mania with mixed features?
#==========================================
bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")
bipolar_id <- bipolar_epall_count_mm_to_depafter_inf %>% dplyr::select(person_id) %>% distinct

categories <- c("Mood", "Affect")

bipolar_depression_nlp <- getDf("bipolar_depression_nlp") %>% mutate(person_id = as.character(person_id))

bipolar_epall_count_mm_to_depafter_inf_depevidence <- left_join(bipolar_epall_count_mm_to_depafter_inf, bipolar_depression_nlp, relationship = "many-to-many") %>%
  filter(measurement_date >= diagnosis_date_mm & measurement_date <= diagnosis_date_end_mm) %>%
  group_by(person_id, diagnosis_date_mm, diagnosis_date_end_mm) %>%
  arrange(measurement_date) %>% slice_head() %>% ungroup()

bipolar_epall_count_mm_to_depafter_inf_depevidence %>% 
  group_by(mm_split, value_as_string) %>%
  summarise(n=n())

bipolar_epall_count_mm_to_depafter_inf_depevidence %>% 
  mutate(mm_duration_days = as.numeric(difftime(diagnosis_date_end_mm, diagnosis_date_mm, units = "days"))) %>%
  filter(mm_split == "mixed" & value_as_string == "Evidence of") %>%
  pull(mm_duration_days) %>% summary
  summary(mm_duration_days)
