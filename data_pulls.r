#=============================================================
# Demographics
#=============================================================

demog <- getDf("demog_cov_sql")

#observation_period_sql <- getDf("observation_period_sql")
#all_mania_mixed_dep_count <- getDf("all_mania_mixed_dep_count")
#observation_period_sql$person_id <- as.character(observation_period_sql$person_id)
#observation_period <- observation_period_sql %>% filter(person_id %in% all_mania_mixed_dep_count$person_id)
#saveDf(observation_period, "observation_period")

demog_cov <- demog %>% 
  mutate(person_id = as.character(person_id)) %>%
  dplyr::select(person_id, year_of_birth, gender, state, race, marital_status, occupation, education_level, observation_period_start_date, observation_period_end_date, observation_years) %>%
  rename(cov_gender = gender, cov_state = state, cov_race = race, cov_rel = marital_status, cov_employ = occupation, cov_edu = education_level) %>%
  mutate(birth_year = as.Date(paste(year_of_birth, "-01-01", sep = ""))) %>%
  #left_join(., observation_period) %>%
  mutate(cov_years_ehr = as.numeric(difftime(as.Date("2025-01-01"), observation_period_start_date, units = "days")/365)) %>%
  mutate(cov_age_ehr = as.numeric(difftime(as.Date("2025-01-01"), as.Date(birth_year), units = "days"))/365) %>%
  mutate(cov_age_entry = as.numeric(difftime(observation_period_start_date, birth_year, units = "days"))/365) %>%
  mutate(cov_state = ifelse(is.na(cov_state), "No matching concept", cov_state)) %>%
  mutate(cov_rel = ifelse(is.na(cov_rel), "No matching concept", cov_rel)) %>%
  mutate(cov_rel = ifelse(cov_rel == "HUNK", "Unknown", cov_rel))
saveDf(demog_cov, "demog_cov", overwrite = T)

#============================================================
# MDS1m assignment code
#============================================================

id_group <- getDf("id_group")

mds1m_group <- pmmd_ep %>% 
  filter(!is.na(dep_date)) %>%
  group_by(person_id) %>%
  summarise(timebwdep_min = min(difftime(dep_date, diagnosis_date_end, units = "days"))) %>%
  mutate(MDS1m_ever = ifelse(timebwdep_min <= 30, 1, 0)) %>% dplyr::select(-timebwdep_min) %>%
  ungroup %>%
  full_join(., id_group %>% filter(group %in% c("PMMD", "noPMMD")) %>% dplyr::select(person_id)) %>%
  mutate(MDS1m_ever = ifelse(is.na(MDS1m_ever), 0, MDS1m_ever))
saveDf(mds1m_group, "mds1m_group", overwrite = T)

#============================================================
# mm_split assignment code
#============================================================

source("diagnosis_definitions.r")
all_mania_dep <- getDf("all_mania_dep")
mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
#mania_dep_dates_toepisodes_qc_broaddef <- getDf("mania_dep_dates_toepisodes_qc_broaddef")
pmmd_ep <- getDf("pmmd_ep")

# Take all_mania_dep
# - for each episode in pmmd_ep_dx, take mm episodes within this time
# - classify into mania_only

# MOVED this code to postmania_depression_epidem.r or the code becomes circular
#all_mania_dep.tmp <- all_mania_dep %>% mutate(diagnosis_date0 = diagnosis_date)
#mm_split <- mania_dep_dates_toepisodes_qc %>%
#  filter(state == "mm") %>%
#  full_join(., all_mania_dep.tmp, relationship = "many-to-many") %>%
#  mutate(state0 = ifelse(icd_code %in% c(mania1, mania2), "mania", "mixed")) %>%
#  filter(diagnosis_date0 >= diagnosis_date & diagnosis_date0 <= diagnosis_date_end) %>%
#  group_by(person_id, diagnosis_date, diagnosis_date_end, state) %>%
#  summarise(mania = length(which(state0 == "mania")), mixed = length(which(state0 == "mixed"))) %>%
#  mutate(mm_split = case_when(mania >= 1 & mixed == 0 ~ "mania",
#                                mania >= 1 & mixed >= 1 & (mania > mixed) ~ "mania",
#                                mania >= 1 & mixed >= 1 & (mania <= mixed) ~ "mixed",
#                                mania == 0 & mixed >= 1 ~ "mixed"))
#saveDf(mm_split, "mm_split", overwrite = T)
mm_split <- getDf("mm_split")

mm_split_fullep.tmp <- pmmd_ep %>% 
  dplyr::select(c(person_id, diagnosis_date, diagnosis_date_end, dep_date, PMMD, group, diagnosis_date_PMMD1, reltoPMMD1, timebwdep, timetodep)) %>%
  left_join(., mm_split) %>%
  rename(mm_date = diagnosis_date, mm_date_end = diagnosis_date_end)

mm_split_dep <- mm_split_fullep.tmp %>%
  filter(!is.na(dep_date)) %>%
  mutate(state = "dep") %>%
  dplyr::select(person_id, dep_date, state, mm_split) %>%
  rename(diagnosis_date = dep_date) %>% distinct %>%
  left_join(., pmmd_ep %>% dplyr::select(person_id, diagnosis_date, diagnosis_date_end, state) %>% distinct)
saveDf(mm_split_dep, "mm_split_dep", overwrite = T)

# - per-episode mm_split
mm_split_mmdep.tmp <- full_join(mm_split, mm_split_dep) %>%
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, state, mm_split) %>%
  left_join(mania_dep_dates_toepisodes_qc_broaddef, .) %>%
  mutate(mm_split = ifelse(is.na(mm_split), "dep_alone", mm_split)) %>%
  dplyr::select(person_id, mm_episode, dep_episode, state, diagnosis_date, diagnosis_date_end, mm_split)
# - per-individual mm_split
mm_split_perind <- mm_split_mmdep.tmp %>% 
  filter(state == "mm") %>% 
  group_by(person_id, mm_split) %>% 
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "mm_split", values_from = "n") %>%
  mutate_if(is.numeric, list(~replace_na(., 0))) %>%
  mutate(prop = mania/(mania+mixed)) %>%
  mutate(mm_split_group = case_when(prop == 1 & mania >= 1 ~ "mania",
                                    prop == 0 & mixed >= 1 ~ "mixed",
                                    .default = "both"))
saveDf(mm_split_perind, "mm_split_perind", overwrite = T)
# - combine together
mm_split_mmdep <- left_join(mm_split_mmdep.tmp, mm_split_perind)
saveDf(mm_split_mmdep, "mm_split_mmdep", overwrite = T)

#----------------------------------
# ICD-code PMMD definition
#----------------------------------

mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
pmmd_ep_dx <- getDf("pmmd_ep_dx")

all_mania_dep.tmp <- all_mania_dep %>% mutate(diagnosis_date0 = diagnosis_date)
mm_split_dx <- mania_dep_dates_toepisodes_qc %>%
  filter(state == "mm") %>%
  full_join(., all_mania_dep.tmp, relationship = "many-to-many") %>%
  mutate(state0 = ifelse(icd_code %in% c(mania1, mania2), "mania", "mixed")) %>%
  filter(diagnosis_date0 >= diagnosis_date & diagnosis_date0 <= diagnosis_date_end) %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end, state) %>%
  summarise(mania = length(which(state0 == "mania")), mixed = length(which(state0 == "mixed"))) %>%
  mutate(mm_split = case_when(mania >= 1 & mixed == 0 ~ "mania",
                                mania >= 1 & mixed >= 1 & (mania > mixed) ~ "mania",
                                mania >= 1 & mixed >= 1 & (mania <= mixed) ~ "mixed",
                                mania == 0 & mixed >= 1 ~ "mixed"))
saveDf(mm_split_dx, "mm_split_dx", overwrite = T)

mm_split_fullep_dx.tmp <- pmmd_ep_dx %>% 
  filter(state == "mm") %>%
  mutate(dep_date = diagnosis_date_end + days(timetodep_dx)) %>%
  dplyr::select(c(person_id, state, diagnosis_date, diagnosis_date_end, dep_date)) %>%
  mutate(group = "PMMD") %>% mutate(PMMD = "PMMD") %>%
  left_join(., mm_split_dx) %>%
  rename(mm_date = diagnosis_date, mm_date_end = diagnosis_date_end)

mm_split_dep_dx <- mm_split_fullep_dx.tmp %>%
  filter(!is.na(dep_date)) %>%
  mutate(state = "dep") %>%
  dplyr::select(person_id, dep_date, state, mm_split) %>%
  rename(diagnosis_date = dep_date) %>%
  left_join(., pmmd_ep_dx %>% dplyr::select(person_id, diagnosis_date, diagnosis_date_end, state))
saveDf(mm_split_dep_dx, "mm_split_dep_dx", overwrite = T)

# - per-episode mm_split
mm_split_mmdep_dx.tmp <- full_join(mm_split_dx, mm_split_dep_dx) %>%
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, state, mm_split) %>%
  left_join(mania_dep_dates_toepisodes_qc, .) %>%
  mutate(mm_split = ifelse(is.na(mm_split), "dep_alone", mm_split)) %>%
  dplyr::select(person_id, mm_episode, dep_episode, state, diagnosis_date, diagnosis_date_end, mm_split)
# - per-individual mm_split
mm_split_perind_dx <- mm_split_mmdep_dx.tmp %>% 
  filter(state == "mm") %>% 
  group_by(person_id, mm_split) %>% 
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "mm_split", values_from = "n") %>%
  mutate_if(is.numeric, list(~replace_na(., 0))) %>%
  mutate(prop = mania/(mania+mixed)) %>%
  mutate(mm_split_group = case_when(prop == 1 & mania >= 1 ~ "mania",
                                    prop == 0 & mixed >= 1 ~ "mixed",
                                    .default = "both"))
saveDf(mm_split_perind_dx, "mm_split_perind_dx", overwrite = T)
# - combine together
mm_split_mmdep_dx <- left_join(mm_split_mmdep_dx.tmp, mm_split_perind_dx)
saveDf(mm_split_mmdep_dx, "mm_split_mmdep_dx", overwrite = T)

#============================================================
# psychosis assignment code
#============================================================

all_mania_dep <- getDf("all_mania_dep")
#mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
mania_dep_dates_toepisodes_qc_broaddef <- getDf("mania_dep_dates_toepisodes_qc_broaddef")
pmmd_ep <- getDf("pmmd_ep")
source("diagnosis_definitions.r")

# Take all_mania_dep
# - for each episode in pmmd_ep_dx, take mm episodes within this time
# - classify into mania_only

all_mania_dep.tmp <- all_mania_dep %>% rename(diagnosis_date0 = diagnosis_date) %>%
  filter(state == "mm") %>%
  mutate(state0 = ifelse(icd_code %in% c(psychosis1, psychosis2), "psychosis", "no_psychosis"))

psychosis_split <- mania_dep_dates_toepisodes_qc_broaddef %>%
  filter(state == "mm") %>%
  left_join(., all_mania_dep.tmp, relationship = "many-to-many") %>%
  filter(diagnosis_date0 >= diagnosis_date & diagnosis_date0 <= diagnosis_date_end) %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end, state) %>%
  summarise(psychosis = length(which(state0 == "psychosis")), no_psychosis = length(which(state0 == "no_psychosis"))) %>%
  mutate(psychosis_split = ifelse(psychosis >= 1, "psychosis", "no_psychosis")) %>%
  dplyr::select(-c("psychosis", "no_psychosis"))
saveDf(psychosis_split, "psychosis_split", overwrite = T)

#============================================================
# severity (max for episode) assignment code
#============================================================

all_mania_dep <- getDf("all_mania_dep")
#mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
mania_dep_dates_toepisodes_qc_broaddef <- getDf("mania_dep_dates_toepisodes_qc_broaddef")
pmmd_ep <- getDf("pmmd_ep")
source("diagnosis_definitions.r")

all_mania_dep.tmp <- all_mania_dep %>% # directly join on diagnosis_date - then are getting mm at the start of the episode
  filter(state == "mm") %>%
  left_join(., icdseverity_score)

#all_mania_dep.tmp <- all_mania_dep %>% rename(diagnosis_date0 = diagnosis_date) %>%
#  filter(state == "mm") %>%
#  left_join(., icdseverity_score)
  
icdseverity_split <- mania_dep_dates_toepisodes_qc_broaddef %>%
  filter(state == "mm") %>%
  left_join(., all_mania_dep.tmp, relationship = "many-to-many") %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end) %>%
  summarise(icdseverity_rank = max(icdseverity_rank)) %>%
  left_join(., icdseverity_score %>% dplyr::select(icdseverity_rank, icdseverity) %>% distinct) # re-add back the mild/mod/severe rating
saveDf(icdseverity_split, "icdseverity_split", overwrite = T)

#=============================================================
# Clinical scales
#=============================================================

#madrs_raw <- getMADRS(as.numeric(id_group %>% filter(group == "PMMD") %>% pull(person_id)))
#saveDf(madrs_raw, "madrs_raw", overwrite = T)

#pmmdgroup_gaf <- getGAF(as.numeric(pmmd_group$person_id))
#nopmmdgroup_gaf <- getGAF(as.numeric(nopmmd_group$person_id))
#saveDf(pmmdgroup_gaf, "pmmdgroup_gaf")
#saveDf(nopmmdgroup_gaf, "nopmmdgroup_gaf")

#gaf_bipolar <- rbind(pmmdgroup_gaf, nopmmdgroup_gaf)
#saveDf(gaf_bipolar, "gaf_bipolar", overwrite = T)

cgis_mm_sql <- getDf("cgis_mm_sql")

pmmdgroup_cgis <- cgis_mm_sql %>% filter(person_id %in% pmmd_group$person_id)
nopmmdgroup_cgis <- cgis_mm_sql %>% filter(person_id %in% nopmmd_group$person_id)
saveDf(pmmdgroup_cgis, "pmmdgroup_cgis")
saveDf(nopmmdgroup_cgis, "nopmmdgroup_cgis")
pmmdgroup_cgis <- getDf("pmmdgroup_cgis")
nopmmdgroup_cgis <- getDf("nopmmdgroup_cgis")

#=============================================================
# Visits
#=============================================================

mm_visit_sql <- getDf("mm_visit_sql")

pmmdgroup_visit <- mm_visit_sql %>% filter(person_id %in% pmmd_group$person_id) %>%
  rename(visit_start_date = visit_date)
saveDf(pmmdgroup_visit, "pmmdgroup_visit", overwrite = T)

nopmmdgroup_visit <- mm_visit_sql %>% filter(person_id %in% nopmmd_group$person_id) %>%
  rename(visit_start_date = visit_date)
saveDf(nopmmdgroup_visit, "nopmmdgroup_visit", overwrite = T)

bipolar_visit <- full_join(pmmdgroup_visit, nopmmdgroup_visit)
bipolar_visit <- rbind(pmmdgroup_visit, nopmmdgroup_visit) %>% mutate(person_id = as.character(person_id)) %>%
  filter(visit_type %in% c("Emergency Room Visit", "Inpatient Visit", "Outpatient Visit", "Home Visit", "Case Management Visit", "Emergency Room and Inpatient Visit", "Telehealth Provided Other than in Patient's Home")) %>%
  mutate(visit_type = case_when(visit_type == "Emergency Room Visit" ~ "Emergency",
                                visit_type %in% c("Inpatient Visit", "Emergency Room and Inpatient Visit") ~ "Inpatient",
                                visit_type %in% c("Outpatient Visit", "Home Visit", "Case Management Visit", "Telehealth Provided Other than in Patient's Home") ~ "Outpatient"))
saveDf(bipolar_visit, "bipolar_visit", overwrite = T)

tmp <- id_group %>% filter(group == "MMonly") %>% filter(is.na(samedayswitch_anyorder)) %>% pull(person_id)
#mmonly_visits <- getVisit(as.numeric(tmp))
#saveDf(mmonly_visits, "mmonlygroup_visit", overwrite = T)

#------------------------------------------------------------
# Overlap episodes with types of visits: PMMD group
#------------------------------------------------------------

pmmdgroup_visit <- getDf("pmmdgroup_visit")
nopmmdgroup_visit <- getDf("nopmmdgroup_visit")

#----------------------
# HOSPITAL
#----------------------
# - take people with any hospital visit
# - take visits overlapping with mania diagnosis meeting following criteria ##### DOUBLE CHECK THIS
# - diagnosis start during first visit (or up to 7 days before) AND before end of visit (ie. start during hospital episode or up to 7 days preceding)
# - diagnosis end up to 7 days after start of visit AND before end of visit (ie. end of mania 7 days after hospital start but before end)
# - diagnosis start before start and end after end (ie. mania extending beyond hospital visit both ways)

# PMMD group
#--------------
pmmd_visit_hospital <- pmmdgroup_visit %>% mutate(person_id = as.character(person_id)) %>% 
  filter(visit_type %in% c("Inpatient Visit", "Emergency Room Visit", "Emergency Room and Inpatient Visit")) %>%
  mutate(hospital_los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")))

pmmd_epall_hospital_overlap.tmp <- inner_join(pmmd_epall, pmmd_visit_hospital, relationship = "many-to-many") %>%
  mutate(timediff_dxstart_vstart = diagnosis_date - visit_start_date,
         timediff_dxstart_vend = diagnosis_date - visit_end_date,
         timediff_dxend_vstart = diagnosis_date_end - visit_start_date,
         timediff_dxend_vend = diagnosis_date_end - visit_end_date) %>%
  filter((timediff_dxstart_vstart >= -7 & timediff_dxstart_vend <= 0) | (timediff_dxend_vstart >= 7 & timediff_dxend_vend <= 0) | (timediff_dxstart_vstart <= 0 & timediff_dxend_vend >= 0))
     
pmmd_epall_hospital_overlap.tmp2 <- pmmd_epall_hospital_overlap.tmp %>% 
  dplyr::select(person_id, mm_episode, dep_episode, state, diagnosis_date, visit_start_date, visit_end_date, hospital_los) %>% distinct() %>%
  inner_join(., pmmd_ep_temporal) %>%
  arrange(person_id, visit_start_date, visit_end_date) %>%
  group_by(person_id, state, mm_episode, dep_episode, ) %>%
  mutate(
     new_visit_episode = case_when(
      row_number() == 1 ~ TRUE,
      (visit_start_date - lag(visit_start_date) >= 7) ~ TRUE,
      # Otherwise, continue the current episode
      TRUE ~ FALSE
    ),
    # Create episode ID by cumulative sum of new_episode flags
    visit_episode = cumsum(new_visit_episode)
  )

pmmd_epall_hospital_overlap <- pmmd_epall_hospital_overlap.tmp2 %>%
  ungroup() %>%
  group_by(person_id, mm_episode, dep_episode, state, diagnosis_date, diagnosis_date_end, PMMD, diagnosis_date_PMMD1, reltoPMMD1) %>%
#  group_by(person_id, mania_episode, dep_episode, state, diagnosis_date, diagnosis_date_end, PMMD, diagnosis_date_PMMD1, reltoPMMD1, visit_episode) %>%
  summarise(
    record_count_visit = n(),
    visit_start_date = min(visit_start_date, na.rm = T),
    visit_end_date = max(visit_end_date, na.rm = T),
    hospital_los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")),
    .groups = "drop"
  ) %>%
  left_join(., mm_split) %>% dplyr::select(-c(mania, mixed))  

pmmd_phase_hosp_ppavg <- pmmd_epall_hospital_overlap %>% 
  group_by(person_id, PMMD, reltoPMMD1, state, mm_split) %>% 
  summarise(n=n()) %>% 
  ungroup() %>%
  group_by(PMMD, reltoPMMD1, state, mm_split) %>% 
  summarise(patients = n(), sum_hosp = sum(n), mean_pp_hosp_n=mean(n), sd_n=sd(n), med_pp_hosp_n=median(n)) %>% 
  ungroup() %>%
  mutate(reltoPMMD1 = fct_relevel(reltoPMMD1, c("before", "PMMD1", "after"))) %>%
  arrange(reltoPMMD1, PMMD, state)

# - check why number of individuals mm < dep
# seems to be because of hospitalisation filter
pmmd_epall_hospital_overlap %>% 
  group_by(person_id, PMMD, reltoPMMD1, state) %>% 
  summarise(n=n()) %>% 
  ungroup() %>%
  pivot_wider(names_from = "state", values_from = "n") %>%
  filter(PMMD == "PMMD" & reltoPMMD1 == "after" & is.na(mm))

saveDf(pmmd_epall_hospital_overlap, "pmmd_epall_hospital_overlap", overwrite = T)
saveDf(pmmd_phase_hosp_ppavg, "pmmd_phase_hosp_ppavg", overwrite = T)

#----------------------
# ANY VISIT
#----------------------

type <- c("Inpatient Visit", "Emergency Room Visit", "Emergency Room and Inpatient Visit", "Case Management Visit", "Home Visit", "Outpatient Visit", "Telehealth Provided Other than in Patient's Home")
pmmd_visit_any <- pmmdgroup_visit %>% mutate(person_id = as.character(person_id)) %>% 
  filter(visit_type %in% type) %>%
  mutate(los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")))
pmmd_epall_any_overlap.tmp <- inner_join(pmmd_epall, pmmd_visit_any, relationship = "many-to-many") %>%
  mutate(timediff_dxstart_vstart = diagnosis_date - visit_start_date,
         timediff_dxstart_vend = diagnosis_date - visit_end_date,
         timediff_dxend_vstart = diagnosis_date_end - visit_start_date,
         timediff_dxend_vend = diagnosis_date_end - visit_end_date) %>%
  filter((timediff_dxstart_vstart >= -7 & timediff_dxstart_vend <= 0) | (timediff_dxend_vstart >= 7 & timediff_dxend_vend <= 0) | (timediff_dxstart_vstart <= 0 & timediff_dxend_vend >= 0))
pmmd_epall_any_overlap.tmp2 <- pmmd_epall_any_overlap.tmp %>% 
  dplyr::select(person_id, state, diagnosis_date, mm_episode, dep_episode, visit_start_date, visit_end_date, los) %>% distinct() %>%
  inner_join(., pmmd_ep_temporal) %>%
  arrange(person_id, visit_start_date, visit_end_date) %>%
  group_by(person_id, state, diagnosis_date) %>%
  mutate(
     new_visit_episode = case_when(
      # First row of each person
      row_number() == 1 ~ TRUE,
      # Both current and previous are mania, but more than 1 month apart
      (visit_start_date - lag(visit_start_date) >= 7) ~ TRUE,
      # Otherwise, continue the current episode
      TRUE ~ FALSE
    ),
    # Create episode ID by cumulative sum of new_episode flags
    visit_episode = cumsum(new_visit_episode)
  )

pmmd_epall_any_overlap <- pmmd_epall_any_overlap.tmp2 %>%
  ungroup() %>%
  group_by(person_id, state, diagnosis_date, diagnosis_date_end, mm_episode, dep_episode, PMMD, diagnosis_date_PMMD1, reltoPMMD1) %>%
  summarise(
    record_count_visit = n(),
    visit_start_date = min(visit_start_date, na.rm = T),
    visit_end_date = max(visit_end_date, na.rm = T),
    los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")),
    .groups = "drop"
  ) %>%
  left_join(., mm_split) %>% dplyr::select(-c(mania, mixed))  

pmmd_phase_visitany_ppavg <- pmmd_epall_any_overlap %>% 
  group_by(person_id, PMMD, reltoPMMD1, state) %>% 
  summarise(n=n()) %>% 
  ungroup() %>%
  group_by(PMMD, reltoPMMD1, state) %>% 
  summarise(patients = n(), sum_any = sum(n), mean_pp_any_n=mean(n), sd_n=sd(n), med_pp_any_n=median(n)) %>% 
  ungroup() %>%
  mutate(reltoPMMD1 = fct_relevel(reltoPMMD1, c("before", "PMMD1", "after"))) %>%
  arrange(reltoPMMD1, PMMD, state)

saveDf(pmmd_phase_visitany_ppavg, "pmmd_phase_visitany_ppavg", overwrite = T)
saveDf(pmmd_epall_any_overlap, "pmmd_epall_any_overlap", overwrite = T)

#------------------------------------------------------------
# Overlap episodes with types of visits: noPMMD group
#------------------------------------------------------------

nopmmdgroup_visit <- getDf("nopmmdgroup_visit")

#----------------------
# HOSPITAL
#----------------------

# noPMMD group
#--------------
nopmmd_visit_hospital <- nopmmdgroup_visit %>% mutate(person_id = as.character(person_id)) %>% 
  filter(visit_type %in% c("Inpatient Visit", "Emergency Room Visit", "Emergency Room and Inpatient Visit")) %>%
  mutate(hospital_los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")))

nopmmd_epall_hospital_overlap.tmp <- inner_join(nopmmd_epall, nopmmd_visit_hospital, relationship = "many-to-many") %>%
  mutate(timediff_dxstart_vstart = diagnosis_date - visit_start_date,
         timediff_dxstart_vend = diagnosis_date - visit_end_date,
         timediff_dxend_vstart = diagnosis_date_end - visit_start_date,
         timediff_dxend_vend = diagnosis_date_end - visit_end_date) %>%
  filter((timediff_dxstart_vstart >= -7 & timediff_dxstart_vend <= 0) | (timediff_dxend_vstart >= 7 & timediff_dxend_vend <= 0) | (timediff_dxstart_vstart <= 0 & timediff_dxend_vend >= 0))
    
nopmmd_epall_hospital_overlap.tmp2 <- nopmmd_epall_hospital_overlap.tmp %>% 
  dplyr::select(person_id, mm_episode, dep_episode, state, diagnosis_date, visit_start_date, visit_end_date, hospital_los) %>% distinct() %>%
  inner_join(., nopmmd_epall) %>%
  arrange(person_id, visit_start_date, visit_end_date) %>%
  group_by(person_id, state, mm_episode, dep_episode, ) %>%
  mutate(
     new_visit_episode = case_when(
      row_number() == 1 ~ TRUE,
      (visit_start_date - lag(visit_start_date) >= 7) ~ TRUE,
      # Otherwise, continue the current episode
      TRUE ~ FALSE
    ),
    # Create episode ID by cumulative sum of new_episode flags
    visit_episode = cumsum(new_visit_episode)
  )

nopmmd_epall_hospital_overlap <- nopmmd_epall_hospital_overlap.tmp2 %>%
  ungroup() %>%
  group_by(person_id, mm_episode, dep_episode, state, diagnosis_date, diagnosis_date_end) %>%
#  group_by(person_id, mania_episode, dep_episode, state, diagnosis_date, diagnosis_date_end, PMMD, diagnosis_date_PMMD1, reltoPMMD1, visit_episode) %>%
  summarise(
    record_count_visit = n(),
    visit_start_date = min(visit_start_date, na.rm = T),
    visit_end_date = max(visit_end_date, na.rm = T),
    hospital_los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")),
    .groups = "drop"
  ) %>%
  left_join(., mm_split) %>% dplyr::select(-c(mania, mixed))  

nopmmd_phase_hosp_ppavg <- nopmmd_epall_hospital_overlap %>% 
  group_by(person_id, state, mm_split) %>% 
  summarise(n=n()) %>% 
  ungroup() %>%
  group_by(state, mm_split) %>% 
  summarise(patients = n(), sum_hosp = sum(n), mean_pp_hosp_n=mean(n), sd_n=sd(n), med_pp_hosp_n=median(n)) %>% 
  ungroup()

# - check why number of individuals mm < dep
# seems to be because of hospitalisation filter
nopmmd_epall_hospital_overlap %>% 
  group_by(person_id, state) %>% 
  summarise(n=n()) %>% 
  ungroup() %>%
  pivot_wider(names_from = "state", values_from = "n") %>%
  mutate(PMMD = "nonPMMD")

saveDf(nopmmd_epall_hospital_overlap, "nopmmd_epall_hospital_overlap", overwrite = T)
saveDf(nopmmd_phase_hosp_ppavg, "nopmmd_phase_hosp_ppavg", overwrite = T)

#----------------------
# ANY VISIT
#----------------------
type <- c("Inpatient Visit", "Emergency Room Visit", "Emergency Room and Inpatient Visit", "Case Management Visit", "Home Visit", "Outpatient Visit", "Telehealth Provided Other than in Patient's Home")
nopmmd_visit_any <- nopmmdgroup_visit %>% mutate(person_id = as.character(person_id)) %>% 
  filter(visit_type %in% type) %>%
  mutate(los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")))
nopmmd_epall_any_overlap.tmp <- inner_join(nopmmd_epall, nopmmd_visit_any, relationship = "many-to-many") %>%
  mutate(timediff_dxstart_vstart = diagnosis_date - visit_start_date,
         timediff_dxstart_vend = diagnosis_date - visit_end_date,
         timediff_dxend_vstart = diagnosis_date_end - visit_start_date,
         timediff_dxend_vend = diagnosis_date_end - visit_end_date) %>%
  filter((timediff_dxstart_vstart >= -7 & timediff_dxstart_vend <= 0) | (timediff_dxend_vstart >= 7 & timediff_dxend_vend <= 0) | (timediff_dxstart_vstart <= 0 & timediff_dxend_vend >= 0))
nopmmd_epall_any_overlap.tmp2 <- nopmmd_epall_any_overlap.tmp %>% 
  dplyr::select(person_id, state, diagnosis_date, mm_episode, dep_episode, visit_start_date, visit_end_date, los) %>% distinct() %>%
  inner_join(., nopmmd_epall) %>%
  arrange(person_id, visit_start_date, visit_end_date) %>%
  group_by(person_id, state, diagnosis_date) %>%
  mutate(
     new_visit_episode = case_when(
      # First row of each person
      row_number() == 1 ~ TRUE,
      # Both current and previous are mania, but more than 1 month apart
      (visit_start_date - lag(visit_start_date) >= 7) ~ TRUE,
      # Otherwise, continue the current episode
      TRUE ~ FALSE
    ),
    # Create episode ID by cumulative sum of new_episode flags
    visit_episode = cumsum(new_visit_episode)
  )

nopmmd_epall_any_overlap <- nopmmd_epall_any_overlap.tmp2 %>%
  ungroup() %>%
  group_by(person_id, state, diagnosis_date, diagnosis_date_end, mm_episode, dep_episode) %>%
  summarise(
    record_count_visit = n(),
    visit_start_date = min(visit_start_date, na.rm = T),
    visit_end_date = max(visit_end_date, na.rm = T),
    hospital_los = as.numeric(difftime(visit_end_date, visit_start_date, units = "days")),
    .groups = "drop"
  ) %>%
  left_join(., mm_split) %>% dplyr::select(-c(mania, mixed)) %>%
  mutate(PMMD = "nonPMMD")

nopmmd_phase_visitany_ppavg <- nopmmd_epall_any_overlap %>% 
  group_by(person_id, PMMD, state) %>% 
  summarise(n=n()) %>% 
  ungroup() %>%
  group_by(PMMD, state) %>% 
  summarise(patients = n(), sum_any = sum(n), mean_pp_any_n=mean(n), sd_n=sd(n), med_pp_any_n=median(n)) %>% 
  ungroup() %>%
  arrange(PMMD, state)

saveDf(nopmmd_phase_visitany_ppavg, "nopmmd_phase_visitany_ppavg", overwrite = T)
saveDf(nopmmd_epall_any_overlap, "nopmmd_epall_any_overlap", overwrite = T)

#=============================================================
# Hospitalisations associated with depression
#=============================================================

pmmdgroup_visit <- getDf("pmmdgroup_visit")
nopmmdgroup_visit <- getDf("nopmmdgroup_visit")
all_mania_dep <- getDf("all_mania_dep")

bipolar_visit <- full_join(pmmdgroup_visit, nopmmdgroup_visit) %>% mutate(person_id = as.character(person_id))
bipolar_visit_hospital <- bipolar_visit %>% filter(visit_type %in% c("Inpatient Visit", "Emergency Room Visit", "Emergency Room and Inpatient Visit"))
bipolar_visit_hospital_dx_7d <- bipolar_visit_hospital %>% 
  left_join(all_mania_dep, relationship = "many-to-many") %>%
  mutate(visit_dx_diff_days = as.numeric(difftime(diagnosis_date, visit_start_date, units = "days"))) %>%
  filter(abs(visit_dx_diff_days) <= 7) %>%
  group_by(person_id, visit_start_date, visit_end_date, visit_type) %>%
  slice_min(abs(visit_dx_diff_days)) %>% 
  ungroup() %>%
  group_by(person_id, visit_type, visit_start_date, visit_end_date, state) %>%
  summarise(dx_n=n())

saveDf(bipolar_visit_hospital_dx_7d, "bipolar_visit_hospital_dx_7d", overwrite = T)

#=============================================================
# Diagnoses (comorbidity)
#=============================================================

alldx_pmmd <- getDiagnosis(as.numeric(pmmd_group$person_id))
saveDf(alldx_pmmd, "alldx_pmmd2", overwrite = T)
alldx_nopmmd <- getDiagnosis(as.numeric(nopmmd_group$person_id))
saveDf(alldx_nopmmd, "alldx_nopmmd", overwrite = T)
alldx_bipolar <- rbind(data.frame(alldx_pmmd, group = "PMMD"), data.frame(alldx_nopmmd, group = "noPMMD"))
saveDf(alldx_bipolar, "alldx_bipolar", overwrite = T)

#------------------------------------------------------------
# Process ICD codes
#------------------------------------------------------------

library(icd.data)
# Split ICD-9 vs ICD-10
alldx_bipolar_icd <- alldx_bipolar %>% 
  mutate(ch1 = substring(icd_code, 1, 1)) %>%
  mutate(icd = ifelse(grepl("[A-Za-z]", ch1), "icd10", "icd9"))

icd9 <- alldx_bipolar_icd %>% filter(icd == "icd9") %>% mutate(person_id = as.character(person_id))
saveDf(icd9, "icd9", overwrite = T)
icd10 <- alldx_bipolar_icd %>% filter(icd == "icd10") %>% mutate(person_id = as.character(person_id))
saveDf(icd10, "icd10", overwrite = T)

# - get number of individuals with >=1 counts of each of ICD-10 broad codes
icd10_labels <- icd10 %>% 
  rename(code = icd_code) %>% mutate(code = gsub("\\.", "", code)) %>%
  inner_join(., icd10cm2016 %>% mutate(code = as.character(code))) %>%
  mutate(icd = "icd10")
saveDf(icd10_labels, "icd10_labels", overwrite = T)

icd9_labels <- icd9 %>% 
  rename(code = icd_code) %>% mutate(code = gsub("\\.", "", code)) %>%
  inner_join(., icd9cm_hierarchy %>% mutate(code = as.character(code))) %>%
  mutate(icd = "icd9")
saveDf(icd9_labels, "icd9_labels", overwrite = T)

#------------------------------------------------------------
# Identify first recorded mood episode for u30
#------------------------------------------------------------

icd10_firstepmood <- icd10 %>% 
  rename(code = icd_code) %>% mutate(code = gsub("\\.", "", code)) %>%
  filter(grepl("F3", code)) %>%
  group_by(person_id) %>%
  filter(diagnosis_date == min(diagnosis_date)) %>%
  distinct() %>%
  rename(diagnosis_date1 = diagnosis_date) %>%
  inner_join(., icd10cm2016 %>% mutate(code = as.character(code))) %>%
  mutate(icd = "icd10")

icd9_firstepmood <- icd9 %>%
  rename(code = icd_code) %>% mutate(code = gsub("\\.", "", code)) %>%
  inner_join(., icd9cm_hierarchy %>% mutate(code = as.character(code))) %>%
  distinct() %>%
  filter(three_digit %in% c(296, 311)) %>%
  group_by(person_id) %>%
  filter(diagnosis_date == min(diagnosis_date)) %>%
  distinct() %>%
  rename(diagnosis_date1 = diagnosis_date) %>%
  mutate(icd = "icd9")

icd_firstepmood <- rbind(icd10_firstepmood %>% dplyr::select(person_id, diagnosis_date1, icd),
                      icd9_firstepmood %>% dplyr::select(person_id, diagnosis_date1, icd)) %>%
  group_by(person_id) %>%
  filter(diagnosis_date1 == min(diagnosis_date1)) %>%
  ungroup() %>%
  dplyr::select(-icd) %>% distinct()
saveDf(icd_firstepmood, "icd_firstepmood", overwrite = T)

icd_firstepmood_u30 <- icd_firstepmood %>%
  left_join(., demog_cov) %>%
  mutate(age_firstep = as.numeric(difftime(diagnosis_date1, birth_year, units = "days")/365)) %>%
  filter(age_firstep <= 30) %>%
  distinct()
saveDf(icd_firstepmood_u30, "icd_firstepmood_u30", overwrite = T)

#=============================================================
# Medications: PMMD group
#=============================================================

#------------------------------------------------------------
# Overlap episodes with types of visits + medications
#------------------------------------------------------------

meds_psych_bipolar_dx <- getDf("meds_psych_bipolar_dx")
pmmd_epall_hospital_overlap <- getDf("pmmd_epall_hospital_overlap")
pmmd_epall_any_overlap <- getDf("pmmd_epall_any_overlap")
mm_split <- getDf("mm_split")

#----------------------
# HOSPITALISATIONS
#----------------------
# Define group
# - hospitalisations coinciding with diagnosis code (to be more robust)
# - medications window: <=7d before mania episode diagnosis and before depression switch
# Prep steps:
# - join hospitalisations and medications data
# - exclude medications given if there is an in-hospital PMMD switch

# Prep data for regression
# - here: define mania/mixed-related hospitalisation as start date within (before or after) 7d of mania/mixed dx
# - filter for:
# - 1) mm only
# - 2) diagnosis date within 7 days of hospital visit start date
# - 3) meds start_date up to 7 days before and 28 days after the hospital visit start date
# - 4) meds start_date preceding dep_date (ie. depression diagnosis)
# - med_episode variable is just to collapse multiple same-day medication records into 1 row
pmmdgroup_mm_7d_hosp_psychmed.tmp <- pmmd_epall_hospital_overlap %>% 
  filter(state == "mm") %>%
  filter((abs(diagnosis_date - visit_start_date) <= 7)) %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  mutate(visit_interval = interval(visit_start_date - 7, visit_start_date + 14)) %>%
  filter(start_date %within% visit_interval) %>%
  left_join(pmmd_ep %>% 
              filter(state == "mm") %>%
              distinct()) %>%
  filter(start_date < dep_date | is.na(dep_date)) %>%
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  group_by(person_id, diagnosis_date, drug) %>%
  arrange(person_id, drug, start_date) %>%
  mutate(
    new_med_episode = case_when(
      row_number() == 1 ~ TRUE,
      (start_date - lag(start_date) >= 7) ~ TRUE,
      TRUE ~ FALSE
    ),
    med_episode = cumsum(new_med_episode)
  )

pmmdgroup_mm_7d_hosp_psychmed <- pmmdgroup_mm_7d_hosp_psychmed.tmp %>%
  ungroup() %>%
  group_by(person_id, mm_episode, dep_episode, diagnosis_date, diagnosis_date_end, drug, PMMD, diagnosis_date_PMMD1, reltoPMMD1, visit_start_date, visit_end_date, hospital_los, broad_class, narrow_class, dep_date, pmmd01, mm_split) %>%
  summarise(
    record_count_med = n(),
    start_date = min(start_date, na.rm = T),
    end_date = max(end_date, na.rm = T),
    drug_window = end_date - start_date,
    .groups = "drop"
  )

# CHECK distribution of medication start date vs visit date
days_medstart_visitstart <- pmmdgroup_mm_7d_hosp_psychmed %>% 
  mutate(days_medstart_visitstart = as.numeric(difftime(start_date, visit_start_date, units = "days"))) %>%
  group_by(person_id, mm_episode, dep_episode, diagnosis_date, diagnosis_date_end, broad_class) %>%
  summarise(days_medstart_visitstart1st = min(days_medstart_visitstart)) %>% 
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"))

ggplot(days_medstart_visitstart, aes(x = days_medstart_visitstart1st)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  coord_cartesian(xlim = c(-7,14)) + scale_x_continuous(breaks = seq(-7, 14, by = 7)) +
  xlab("Days from hospital visit start to earliest prescription record within window") +
  facet_wrap(~broad_class)#, scale = "free_y")
ggsave("days_medstart_visitstart.png", height = 3.5, width = 7)

# CHECK distribution of medication start date vs manic/mixed date
days_medstart_mmstart <- pmmdgroup_mm_7d_hosp_psychmed %>% 
  mutate(days_medstart_mmstart = as.numeric(difftime(start_date, diagnosis_date, units = "days"))) %>%
  group_by(person_id, mm_episode, dep_episode, diagnosis_date, diagnosis_date_end, broad_class) %>%
  summarise(days_medstart_mmstart1st = min(days_medstart_mmstart)) %>% 
  filter(broad_class %in% c("lithium", "antiepileptic", "antipsychotic", "antidepressant", "antihistamine", "benzodiazepine"))

ggplot(days_medstart_mmstart, aes(x = days_medstart_mmstart1st)) +
  geom_histogram(binwidth=1) +
  theme_bw() +
#  coord_cartesian(xlim = c(-7,28)) + scale_x_continuous(breaks = seq(-7, 28, by = 7)) +
  xlab("Days from manic/mixed episode start to earliest prescription record within window") +
  facet_wrap(~broad_class)#, scale = "free_y")
ggsave("days_medstart_mmstart.png", height = 3.5, width = 7)

# CHECK that medications are prescribed before depression
days_medstart_to_depdate <- pmmdgroup_mm_7d_hosp_psychmed %>% 
  dplyr::select(person_id, diagnosis_date, dep_date, start_date, end_date, broad_class) %>% 
  distinct %>%
  mutate(days_medstart_to_depdate = as.numeric(difftime(dep_date, start_date, units = "days")))

ggplot(days_medstart_to_depdate, aes(x = days_medstart_to_depdate)) +
  geom_density() +
  theme_bw() +
  coord_cartesian(xlim = c(0, 180)) +
  facet_wrap(~broad_class)
ggsave("days_medstart_to_depdate.png")

days_medstart_to_depdate %>% 
  group_by(broad_class) %>%
  summarise(n=n(),
            mean = mean(days_medstart_to_depdate, na.rm = T), 
            sd = sd(days_medstart_to_depdate, na.rm = T), 
            median = median(days_medstart_to_depdate, na.rm = T), 
            q25 = quantile(days_medstart_to_depdate, 0.25, na.rm = T), 
            q75 = quantile(days_medstart_to_depdate, 0.75, na.rm = T))

# - pad out
pmmdgroup_mm_7d_hosp_psychmeds_count <- pmmdgroup_mm_7d_hosp_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  mutate(polypharmacy_time = polypharmacy / as.numeric(difftime(visit_end_date, visit_start_date, units = "days"))) %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, visit_start_date, visit_end_date, hospital_los, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  pivot_wider(., names_from = "broad_class", values_from = "present", values_fill = 0)
pmmdgroup_mm_7d_hosp_psychmeds_count

# - make into long format, padded out with 0s
broadclass_unique <- unique(pmmdgroup_mm_7d_hosp_psychmed$broad_class)
pmmdgroup_mm_7d_hosp_psychmeds_count.long <- pmmdgroup_mm_7d_hosp_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique), names_to = "broadclass", values_to = "drug01")

pmmdgroup_mm_7d_hosp_psychmeds_count_quant <- pmmdgroup_mm_7d_hosp_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, visit_start_date, visit_end_date, hospital_los, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  pivot_wider(., names_from = "broad_class", values_from = "n", values_fill = 0)

# - exclude when the number of people with a class is <= 20
cols.tmp <- unique(pmmdgroup_mm_7d_hosp_psychmed$broad_class)
count_broadclass_exc <- pmmdgroup_mm_7d_hosp_psychmeds_count[c(cols.tmp, "mm_split")] %>%
  group_by(mm_split) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
  filter(mm_split != "both") %>%
  pivot_longer(!mm_split, names_to = "drug", values_to = "n") %>%
  filter(n<=20) %>%
  pull(drug)
count_broadclass_keep <- cols.tmp[-which(cols.tmp %in% count_broadclass_exc)]
#count_broadclass <- colSums(pmmdgroup_mm_7d_hosp_psychmeds_count[,cols.tmp])
#count_broadclass_exc <- names(count_broadclass[which(count_broadclass <= 20)])

saveDf(pmmdgroup_mm_7d_hosp_psychmed, "pmmdgroup_mm_7d_hosp_psychmed", overwrite = T)
saveDf(pmmdgroup_mm_7d_hosp_psychmeds_count, "pmmdgroup_mm_7d_hosp_psychmeds_count", overwrite = T)
saveDf(pmmdgroup_mm_7d_hosp_psychmeds_count.long, "pmmdgroup_mm_7d_hosp_psychmeds_count.long", overwrite = T)

#-------------
# Specific medications
#-------------
pmmdgroup_mm_7d_hosp_medsind_count <- pmmdgroup_mm_7d_hosp_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, pmmd01, visit_start_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  ungroup() %>%
  #filter(broad_class %in% c("antiepileptic", "lithium", "antipsychotic", "antidepressant", "antihistamine")) %>%
  mutate(drug = case_when(broad_class == "lithium" ~ "lithium",
                          drug %in% c("aripiprazole", "aripiprazole lauroxil") ~ "aripiprazole",
                         .default = drug)) %>%
  group_by(person_id, visit_start_date) %>%
  #mutate(tmp_rm = ifelse(any(broad_class %in% c("antiepileptic", "lithium", "antipsychotic")), "rm_nomsap", "keep")) %>%
  ungroup() %>%
  #filter(!(drug == "no.MS.or.AP" & tmp_rm == "rm_nomsap")) %>%
  #dplyr::select(-tmp_rm) %>%
  group_by(person_id, pmmd01, diagnosis_date, visit_start_date, visit_end_date, mm_split, dep_date, hospital_los, drug, polypharmacy) %>%
  count(drug) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  distinct() %>%
  pivot_wider(., names_from = "drug", values_from = "present", values_fill = 0)
saveDf(pmmdgroup_mm_7d_hosp_medsind_count, "pmmdgroup_mm_7d_hosp_medsind_count", overwrite = T)

# - make into long format, padded out with 0s
meds_unique <- pmmdgroup_mm_7d_hosp_psychmed$drug
pmmdgroup_mm_7d_hosp_medsind_count.long <- pmmdgroup_mm_7d_hosp_medsind_count %>%
  pivot_longer(cols = any_of(meds_unique), names_to = "meds", values_to = "drug01") %>%
  left_join(., psych_med_df %>% rename(meds = drug))
saveDf(pmmdgroup_mm_7d_hosp_medsind_count.long, "pmmdgroup_mm_7d_hosp_medsind_count.long", overwrite = T)

#----------------------
# ANY EPISODE
#----------------------
pmmdgroup_mm_7d_anyvisit_psychmed.tmp <- pmmd_epall_any_overlap %>% 
  filter(state == "mm") %>%
  filter((abs(diagnosis_date - visit_start_date) <= 7)) %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  mutate(visit_interval = interval(visit_start_date - 7, visit_start_date + 14)) %>%
  filter(start_date %within% visit_interval) %>%
  left_join(pmmd_ep %>%
              filter(state == "mm") %>%
              distinct()) %>%
  filter(start_date < dep_date | is.na(dep_date)) %>%
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  arrange(person_id, drug, start_date) %>%
  group_by(person_id, diagnosis_date, drug) %>%
  mutate(new_med_episode = case_when(
      row_number() == 1 ~ TRUE,
      (start_date - lag(start_date) >= 7) ~ TRUE,
      TRUE ~ FALSE),
    med_episode = cumsum(new_med_episode)) %>%
  ungroup()

pmmdgroup_mm_7d_anyvisit_psychmed <- pmmdgroup_mm_7d_anyvisit_psychmed.tmp %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end, drug, PMMD, diagnosis_date_PMMD1, reltoPMMD1, visit_start_date, visit_end_date, broad_class, narrow_class, dep_date, pmmd01, mm_split) %>%
  summarise(
    record_count_med = n(),
    start_date = min(start_date, na.rm = T),
    end_date = max(end_date, na.rm = T),
    drug_window = end_date - start_date,
    .groups = "drop") %>%
  ungroup() %>%
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMD", "PMMD")) %>%
  left_join(., pmmdgroup_mm_7d_hosp_psychmed) %>% mutate(hospital_los = ifelse(is.na(hospital_los), 0, hospital_los))

pmmdgroup_mm_7d_anyvisit_psychmeds_count <- pmmdgroup_mm_7d_anyvisit_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>% ungroup %>%
  mutate(polypharmacy_time = polypharmacy / as.numeric(difftime(visit_end_date, visit_start_date, units = "days"))) %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, visit_start_date, visit_end_date, hospital_los, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  pivot_wider(., names_from = "broad_class", values_from = "present", values_fill = 0)

# - make into long format, padded out with 0s
broadclass_unique <- unique(pmmdgroup_mm_7d_anyvisit_psychmed$broad_class)
pmmdgroup_mm_7d_anyvisit_psychmeds_count.long <- pmmdgroup_mm_7d_anyvisit_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique), names_to = "broadclass", values_to = "drug01")

pmmdgroup_mm_7d_anyvisit_psychmeds_count_quant <- pmmdgroup_mm_7d_anyvisit_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, visit_start_date, visit_end_date, hospital_los, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  pivot_wider(., names_from = "broad_class", values_from = "n", values_fill = 0)

saveDf(pmmdgroup_mm_7d_anyvisit_psychmed, "pmmdgroup_mm_7d_anyvisit_psychmed", overwrite = T)
saveDf(pmmdgroup_mm_7d_anyvisit_psychmeds_count, "pmmdgroup_mm_7d_anyvisit_psychmeds_count", overwrite = T)
saveDf(pmmdgroup_mm_7d_anyvisit_psychmeds_count.long, "pmmdgroup_mm_7d_anyvisit_psychmeds_count.long", overwrite = T)

#----------------------
# NON-HOSPITAL EPISODES
#----------------------

pmmdgroup_mm_7d_nohosp_psychmed <- pmmdgroup_mm_7d_anyvisit_psychmed %>%
  anti_join(., pmmdgroup_mm_7d_hosp_psychmed, by = c("person_id", "diagnosis_date", "diagnosis_date_end"))

pmmdgroup_mm_7d_nohosp_psychmeds_count <- pmmdgroup_mm_7d_anyvisit_psychmeds_count %>%
  anti_join(., pmmdgroup_mm_7d_hosp_psychmeds_count, by = c("person_id", "diagnosis_date", "diagnosis_date_end"))

broadclass_unique <- unique(pmmdgroup_mm_7d_nohosp_psychmed$broad_class)
pmmdgroup_mm_7d_nohosp_psychmeds_count.long <- pmmdgroup_mm_7d_nohosp_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique), names_to = "broadclass", values_to = "drug01")

saveDf(pmmdgroup_mm_7d_nohosp_psychmed, "pmmdgroup_mm_7d_nohosp_psychmed", overwrite = T)
saveDf(pmmdgroup_mm_7d_nohosp_psychmeds_count, "pmmdgroup_mm_7d_nohosp_psychmeds_count", overwrite = T)
saveDf(pmmdgroup_mm_7d_nohosp_psychmeds_count.long, "pmmdgroup_mm_7d_nohosp_psychmeds_count.long", overwrite = T)

#----------------------
# NO VISIT CONSTRAINT
#----------------------

pmmdgroup_mm_7d_dx_psychmed.tmp <- pmmd_epall %>% left_join(pmmd_ep_temporal %>% dplyr::select(person_id, diagnosis_date, diagnosis_date_end, reltoPMMD1) %>% distinct) %>%
  left_join(., mm_split) %>%
  filter(state == "mm") %>% distinct %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  mutate(valid_interval = interval(diagnosis_date - 7, diagnosis_date + 14)) %>%
  filter(start_date %within% valid_interval) %>%
  left_join(pmmd_ep %>% filter(state == "mm") %>% distinct()) %>%
  filter(start_date < dep_date | is.na(dep_date)) %>% # ensure that medication is started before depression
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  arrange(person_id, drug, start_date) %>%
  group_by(person_id, diagnosis_date, drug) %>%
  mutate(new_med_episode = case_when(
      row_number() == 1 ~ TRUE,
      (start_date - lag(start_date) >= 7) ~ TRUE,
      TRUE ~ FALSE),
    med_episode = cumsum(new_med_episode)) %>%
  ungroup

pmmdgroup_mm_7d_dx_psychmed <- pmmdgroup_mm_7d_dx_psychmed.tmp %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end, drug, PMMD, diagnosis_date_PMMD1, reltoPMMD1, broad_class, narrow_class, dep_date, pmmd01, mm_split) %>%
  summarise(
    record_count_med = n(),
    start_date = min(start_date, na.rm = T),
    end_date = max(end_date, na.rm = T),
    drug_window = end_date - start_date,
    .groups = "drop") %>%
  ungroup %>%
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMD", "PMMD"))

pmmdgroup_mm_7d_dx_psychmeds_count <- pmmdgroup_mm_7d_dx_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>% ungroup() %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  pivot_wider(., names_from = "broad_class", values_from = "present", values_fill = 0)

# - make into long format, padded out with 0s
broadclass_unique <- unique(pmmdgroup_mm_7d_dx_psychmed$broad_class)
pmmdgroup_mm_7d_dx_psychmeds_count.long <- pmmdgroup_mm_7d_dx_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique), names_to = "broadclass", values_to = "drug01")

pmmdgroup_mm_7d_dx_psychmeds_count_quant <- pmmdgroup_mm_7d_dx_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>% ungroup %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  pivot_wider(., names_from = "broad_class", values_from = "n", values_fill = 0)

saveDf(pmmdgroup_mm_7d_dx_psychmed, "pmmdgroup_mm_7d_dx_psychmed", overwrite = T)
saveDf(pmmdgroup_mm_7d_dx_psychmeds_count, "pmmdgroup_mm_7d_dx_psychmeds_count", overwrite = T)
saveDf(pmmdgroup_mm_7d_dx_psychmeds_count.long, "pmmdgroup_mm_7d_dx_psychmeds_count.long", overwrite = T)

#-------------
# Specific medications
#-------------
psych_med_df <- read.csv("psych_med_df.csv")
pmmdgroup_mm_7d_dx_medsind_count <- pmmdgroup_mm_7d_dx_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  ungroup() %>%
  #filter(broad_class %in% c("antiepileptic", "lithium", "antipsychotic", "antidepressant", "antihistamine")) %>%
  mutate(drug = case_when(broad_class == "lithium" ~ "lithium",
                          drug %in% c("aripiprazole", "aripiprazole lauroxil") ~ "aripiprazole",
                         .default = drug)) %>%
  #group_by(person_id, diagnosis_date) %>%
  #mutate(tmp_rm = ifelse(any(broad_class %in% c("antiepileptic", "lithium", "antipsychotic")), "rm_nomsap", "keep")) %>%
  #ungroup() %>%
  #filter(!(drug == "no.MS.or.AP" & tmp_rm == "rm_nomsap")) %>%
  #dplyr::select(-tmp_rm) %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, drug, polypharmacy) %>%
  count(drug) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  distinct() %>%
  pivot_wider(., names_from = "drug", values_from = "present", values_fill = 0)

# - make into long format, padded out with 0s
meds_unique <- pmmdgroup_mm_7d_dx_psychmed$drug
pmmdgroup_mm_7d_dx_medsind_count.long <- pmmdgroup_mm_7d_dx_medsind_count %>%
  pivot_longer(cols = any_of(meds_unique), names_to = "meds", values_to = "drug01") %>%
  left_join(., psych_med_df %>% rename(meds = drug))

saveDf(pmmdgroup_mm_7d_dx_medsind_count, "pmmdgroup_mm_7d_dx_medsind_count", overwrite = T)
saveDf(pmmdgroup_mm_7d_dx_medsind_count.long, "pmmdgroup_mm_7d_dx_medsind_count.long", overwrite = T)

#=============================================================
# Medications: noPMMD group
#=============================================================

#------------------------------------------------------------
# Overlap episodes with types of visits + medications
#------------------------------------------------------------

meds_psych_bipolar_dx <- getDf("meds_psych_bipolar_dx")
nopmmd_epall_hospital_overlap <- getDf("nopmmd_epall_hospital_overlap")
nopmmd_epall_any_overlap <- getDf("nopmmd_epall_any_overlap")

#----------------------
# HOSPITALISATIONS
#----------------------

nopmmdgroup_mm_7d_hosp_psychmed.tmp <- nopmmd_epall_hospital_overlap %>% 
  filter(state == "mm") %>%
  filter((abs(diagnosis_date - visit_start_date) <= 7)) %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  mutate(visit_interval = interval(visit_start_date - 7, visit_start_date + 14)) %>%
  filter(start_date %within% visit_interval) %>%
  left_join(pmmd_ep %>% 
              filter(state == "mm") %>%
              distinct()) %>%
  filter(start_date < dep_date | is.na(dep_date)) %>%
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  arrange(person_id, drug, start_date) %>%
  group_by(person_id, diagnosis_date, drug) %>%
  mutate(
    new_med_episode = case_when(
      row_number() == 1 ~ TRUE,
      (start_date - lag(start_date) >= 7) ~ TRUE,
      TRUE ~ FALSE),
    med_episode = cumsum(new_med_episode)) %>%
  ungroup

nopmmdgroup_mm_7d_hosp_psychmed <- nopmmdgroup_mm_7d_hosp_psychmed.tmp %>%
  ungroup() %>%
  group_by(person_id, mm_episode, dep_episode, diagnosis_date, diagnosis_date_end, drug, PMMD, diagnosis_date_PMMD1, reltoPMMD1, visit_start_date, visit_end_date, hospital_los, broad_class, narrow_class, dep_date, pmmd01, mm_split) %>%
  summarise(
    record_count_med = n(),
    start_date = min(start_date, na.rm = T),
    end_date = max(end_date, na.rm = T),
    drug_window = end_date - start_date,
    .groups = "drop"
  ) %>% ungroup

# - pad out
nopmmdgroup_mm_7d_hosp_psychmeds_count <- nopmmdgroup_mm_7d_hosp_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  ungroup %>%
  mutate(polypharmacy_time = polypharmacy / as.numeric(difftime(visit_end_date, visit_start_date, units = "days"))) %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, visit_start_date, visit_end_date, hospital_los, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  pivot_wider(., names_from = "broad_class", values_from = "present", values_fill = 0)
nopmmdgroup_mm_7d_hosp_psychmeds_count

# - make into long format, padded out with 0s
broadclass_unique_nopmmd <- unique(nopmmdgroup_mm_7d_hosp_psychmed$broad_class)
nopmmdgroup_mm_7d_hosp_psychmeds_count.long <- nopmmdgroup_mm_7d_hosp_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique_nopmmd), names_to = "broadclass", values_to = "drug01")

nopmmdgroup_mm_7d_hosp_psychmeds_count_quant <- nopmmdgroup_mm_7d_hosp_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>% ungroup %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, visit_start_date, visit_end_date, hospital_los, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  pivot_wider(., names_from = "broad_class", values_from = "n", values_fill = 0)

# - exclude when the number of people with a class is <= 20
cols_nopmmd.tmp <- unique(nopmmdgroup_mm_7d_hosp_psychmed$broad_class)
count_broadclass_exc_nopmmd <- pmmdgroup_mm_7d_hosp_psychmeds_count[c(cols_nopmmd.tmp, "mm_split")] %>%
  group_by(mm_split) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
  ungroup %>%
  filter(mm_split != "both") %>%
  pivot_longer(!mm_split, names_to = "drug", values_to = "n") %>%
  filter(n<=20) %>%
  pull(drug)
count_broadclass_keep_nopmmd <- cols.tmp[-which(cols.tmp %in% count_broadclass_exc_nopmmd)]
#count_broadclass <- colSums(pmmdgroup_mm_7d_hosp_psychmeds_count[,cols.tmp])
#count_broadclass_exc <- names(count_broadclass[which(count_broadclass <= 20)])

saveDf(nopmmdgroup_mm_7d_hosp_psychmed, "nopmmdgroup_mm_7d_hosp_psychmed", overwrite = T)
saveDf(nopmmdgroup_mm_7d_hosp_psychmeds_count, "nopmmdgroup_mm_7d_hosp_psychmeds_count", overwrite = T)
saveDf(nopmmdgroup_mm_7d_hosp_psychmeds_count.long, "nopmmdgroup_mm_7d_hosp_psychmeds_count.long", overwrite = T)

#-------------
# Specific medications
#-------------
nopmmdgroup_mm_7d_hosp_medsind_count <- nopmmdgroup_mm_7d_hosp_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, pmmd01, visit_start_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  ungroup() %>%
  #filter(broad_class %in% c("antiepileptic", "lithium", "antipsychotic", "antidepressant", "antihistamine")) %>%
  mutate(drug = case_when(broad_class == "lithium" ~ "lithium",
                          drug %in% c("aripiprazole", "aripiprazole lauroxil") ~ "aripiprazole",
                         .default = drug)) %>%
  #group_by(person_id, visit_start_date) %>%
  #mutate(tmp_rm = ifelse(any(broad_class %in% c("antiepileptic", "lithium", "antipsychotic")), "rm_nomsap", "keep")) %>%
  #ungroup() %>%
  #filter(!(drug == "no.MS.or.AP" & tmp_rm == "rm_nomsap")) %>%
  #dplyr::select(-tmp_rm) %>%
  group_by(person_id, pmmd01, diagnosis_date, visit_start_date, visit_end_date, mm_split, dep_date, hospital_los, drug, polypharmacy) %>%
  count(drug) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  distinct() %>%
  pivot_wider(., names_from = "drug", values_from = "present", values_fill = 0)
saveDf(nopmmdgroup_mm_7d_hosp_medsind_count, "nopmmdgroup_mm_7d_hosp_medsind_count", overwrite = T)

# - make into long format, padded out with 0s
meds_unique <- nopmmdgroup_mm_7d_hosp_psychmed$drug
nopmmdgroup_mm_7d_hosp_medsind_count.long <- nopmmdgroup_mm_7d_hosp_medsind_count %>%
  pivot_longer(cols = any_of(meds_unique), names_to = "meds", values_to = "drug01") %>%
  left_join(., psych_med_df %>% rename(meds = drug))
saveDf(nopmmdgroup_mm_7d_hosp_medsind_count.long, "nopmmdgroup_mm_7d_hosp_medsind_count.long", overwrite = T)

#----------------------
# ANY EPISODE
#----------------------
nopmmdgroup_mm_7d_anyvisit_psychmed.tmp <- nopmmd_epall_any_overlap %>% 
  filter(state == "mm") %>%
  filter((abs(diagnosis_date - visit_start_date) <= 7)) %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  mutate(visit_interval = interval(visit_start_date - 7, visit_start_date + 14)) %>%
  filter(start_date %within% visit_interval) %>%
  left_join(nopmmd_epall %>% 
              filter(state == "mm") %>%
              distinct()) %>%
  #filter(start_date < dep_date | is.na(dep_date)) %>% # THIS STEP only for PMMD group
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  arrange(person_id, drug, start_date) %>%
  group_by(person_id, diagnosis_date, drug) %>%
  mutate(new_med_episode = case_when(
      row_number() == 1 ~ TRUE,
      (start_date - lag(start_date) >= 7) ~ TRUE,
      TRUE ~ FALSE),
    med_episode = cumsum(new_med_episode))

nopmmdgroup_mm_7d_anyvisit_psychmed <- nopmmdgroup_mm_7d_anyvisit_psychmed.tmp %>%
  ungroup() %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end, drug, PMMD, visit_start_date, visit_end_date, broad_class, narrow_class, pmmd01, mm_split) %>%
  summarise(
    record_count_med = n(),
    start_date = min(start_date, na.rm = T),
    end_date = max(end_date, na.rm = T),
    drug_window = end_date - start_date,
    .groups = "drop") %>% ungroup %>%
  left_join(., nopmmdgroup_mm_7d_hosp_psychmed) %>% mutate(hospital_los = ifelse(is.na(hospital_los), 0, hospital_los))

nopmmdgroup_mm_7d_anyvisit_psychmeds_count <- nopmmdgroup_mm_7d_anyvisit_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>% ungroup %>%
  mutate(polypharmacy_time = polypharmacy / as.numeric(difftime(visit_end_date, visit_start_date, units = "days"))) %>%
  group_by(person_id, pmmd01, PMMD, mm_split, diagnosis_date, diagnosis_date_end, visit_start_date, visit_end_date, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  pivot_wider(., names_from = "broad_class", values_from = "present", values_fill = 0)

# - make into long format, padded out with 0s
broadclass_unique <- nopmmdgroup_mm_7d_anyvisit_psychmed %>% filter(!is.na(broad_class)) %>% pull(broad_class) %>% unique
nopmmdgroup_mm_7d_anyvisit_psychmeds_count.long <- nopmmdgroup_mm_7d_anyvisit_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique), names_to = "broadclass", values_to = "drug01")

nopmmdgroup_mm_7d_anyvisit_psychmeds_count_quant <- nopmmdgroup_mm_7d_anyvisit_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  group_by(person_id, pmmd01, PMMD, mm_split, diagnosis_date, diagnosis_date_end, visit_start_date, visit_end_date, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  pivot_wider(., names_from = "broad_class", values_from = "n", values_fill = 0)

saveDf(nopmmdgroup_mm_7d_anyvisit_psychmed, "nopmmdgroup_mm_7d_anyvisit_psychmed", overwrite = T)
saveDf(nopmmdgroup_mm_7d_anyvisit_psychmeds_count, "nopmmdgroup_mm_7d_anyvisit_psychmeds_count", overwrite = T)
saveDf(nopmmdgroup_mm_7d_anyvisit_psychmeds_count.long, "nopmmdgroup_mm_7d_anyvisit_psychmeds_count.long", overwrite = T)

#----------------------
# NON-HOSPITAL EPISODES
#----------------------

nopmmdgroup_mm_7d_nohosp_psychmed <- nopmmdgroup_mm_7d_anyvisit_psychmed %>%
  anti_join(., nopmmdgroup_mm_7d_hosp_psychmed, by = c("person_id", "diagnosis_date", "diagnosis_date_end"))

nopmmdgroup_mm_7d_nohosp_psychmeds_count <- nopmmdgroup_mm_7d_anyvisit_psychmeds_count %>%
  anti_join(., nopmmdgroup_mm_7d_hosp_psychmeds_count, by = c("person_id", "diagnosis_date", "diagnosis_date_end"))

broadclass_unique <- nopmmdgroup_mm_7d_nohosp_psychmed %>% filter(!is.na(broad_class)) %>% pull(broad_class) %>% unique
nopmmdgroup_mm_7d_nohosp_psychmeds_count.long <- nopmmdgroup_mm_7d_nohosp_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique), names_to = "broadclass", values_to = "drug01")

saveDf(nopmmdgroup_mm_7d_nohosp_psychmed, "nopmmdgroup_mm_7d_nohosp_psychmed", overwrite = T)
saveDf(nopmmdgroup_mm_7d_nohosp_psychmeds_count, "nopmmdgroup_mm_7d_nohosp_psychmeds_count", overwrite = T)
saveDf(nopmmdgroup_mm_7d_nohosp_psychmeds_count.long, "nopmmdgroup_mm_7d_nohosp_psychmeds_count.long", overwrite = T)

#----------------------
# NO VISIT CONSTRAINT
#----------------------

nopmmdgroup_mm_7d_dx_psychmed.tmp <- nopmmd_epall %>% mutate(duration = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days"))) %>%
  left_join(., mm_split) %>%
  filter(state == "mm") %>% distinct %>%
  left_join(., meds_psych_bipolar_dx, relationship = "many-to-many") %>%
  mutate(valid_interval = interval(diagnosis_date - 7, diagnosis_date + 14)) %>%
  filter(start_date %within% valid_interval) %>%
  left_join(pmmd_ep %>%
              filter(state == "mm") %>%
              distinct()) %>%
  filter(start_date < dep_date | is.na(dep_date)) %>% # ensure that medication is started before depression
  mutate(pmmd01 = ifelse(is.na(PMMD), 0, 1)) %>%
  arrange(person_id, drug, start_date) %>%
  group_by(person_id, diagnosis_date, drug) %>%
  mutate(new_med_episode = case_when(
      row_number() == 1 ~ TRUE,
      (start_date - lag(start_date) >= 7) ~ TRUE,
      TRUE ~ FALSE),
    med_episode = cumsum(new_med_episode)) %>% ungroup

nopmmdgroup_mm_7d_dx_psychmed <- nopmmdgroup_mm_7d_dx_psychmed.tmp %>%
  group_by(person_id, diagnosis_date, diagnosis_date_end, drug, PMMD, diagnosis_date_PMMD1, reltoPMMD1, broad_class, narrow_class, dep_date, pmmd01, mm_split) %>%
  summarise(
    record_count_med = n(),
    start_date = min(start_date, na.rm = T),
    end_date = max(end_date, na.rm = T),
    drug_window = end_date - start_date,
    .groups = "drop") %>% ungroup %>%
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMD", "PMMD"))

nopmmdgroup_mm_7d_dx_psychmeds_count <- nopmmdgroup_mm_7d_dx_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>% ungroup %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  pivot_wider(., names_from = "broad_class", values_from = "present", values_fill = 0)

# - make into long format, padded out with 0s
broadclass_unique <- nopmmdgroup_mm_7d_dx_psychmed %>% filter(!is.na(broad_class)) %>% pull(broad_class) %>% unique
nopmmdgroup_mm_7d_dx_psychmeds_count.long <- nopmmdgroup_mm_7d_dx_psychmeds_count %>%
  pivot_longer(cols = all_of(broadclass_unique), names_to = "broadclass", values_to = "drug01")

nopmmdgroup_mm_7d_dx_psychmeds_count_quant <- nopmmdgroup_mm_7d_dx_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, PMMD, reltoPMMD1, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  group_by(person_id, pmmd01, PMMD, reltoPMMD1, mm_split, diagnosis_date, diagnosis_date_end, dep_date, broad_class, polypharmacy) %>%
  count(broad_class) %>% ungroup() %>%
  pivot_wider(., names_from = "broad_class", values_from = "n", values_fill = 0)

saveDf(nopmmdgroup_mm_7d_dx_psychmed, "nopmmdgroup_mm_7d_dx_psychmed", overwrite = T)
saveDf(nopmmdgroup_mm_7d_dx_psychmeds_count, "nopmmdgroup_mm_7d_dx_psychmeds_count", overwrite = T)
saveDf(nopmmdgroup_mm_7d_dx_psychmeds_count.long, "nopmmdgroup_mm_7d_dx_psychmeds_count.long", overwrite = T)

#-------------
# Specific medications
#-------------
nopmmdgroup_mm_7d_dx_medsind_count <- nopmmdgroup_mm_7d_dx_psychmed %>%
  left_join(., mm_split) %>%
  group_by(person_id, pmmd01, diagnosis_date) %>%
  mutate(polypharmacy = length(unique(drug))) %>%
  ungroup() %>%
  #filter(broad_class %in% c("antiepileptic", "lithium", "antipsychotic", "antidepressant", "antihistamine")) %>%
  mutate(drug = case_when(broad_class == "lithium" ~ "lithium",
                          drug %in% c("aripiprazole", "aripiprazole lauroxil") ~ "aripiprazole",
                         .default = drug)) %>%
  #group_by(person_id, diagnosis_date) %>%
  #mutate(tmp_rm = ifelse(any(broad_class %in% c("antiepileptic", "lithium", "antipsychotic")), "rm_nomsap", "keep")) %>%
  #ungroup() %>%
  #filter(!(drug == "no.MS.or.AP" & tmp_rm == "rm_nomsap")) %>%
  #dplyr::select(-tmp_rm) %>%
  group_by(person_id, pmmd01, mm_split, diagnosis_date, diagnosis_date_end, dep_date, drug, polypharmacy) %>%
  count(drug) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(-n) %>%
  distinct() %>%
  pivot_wider(., names_from = "drug", values_from = "present", values_fill = 0)

# - make into long format, padded out with 0s
meds_unique <- nopmmdgroup_mm_7d_dx_psychmed$drug
nopmmdgroup_mm_7d_dx_medsind_count.long <- nopmmdgroup_mm_7d_dx_medsind_count %>%
  pivot_longer(cols = any_of(meds_unique), names_to = "meds", values_to = "drug01") %>%
  left_join(., psych_med_df %>% rename(meds = drug))

saveDf(nopmmdgroup_mm_7d_dx_medsind_count, "nopmmdgroup_mm_7d_dx_medsind_count", overwrite = T)
saveDf(nopmmdgroup_mm_7d_dx_medsind_count.long, "nopmmdgroup_mm_7d_dx_medsind_count.long", overwrite = T)

