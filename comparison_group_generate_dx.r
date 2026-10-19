#======================================================
# Libraries and functions
#======================================================

source("Setup_background.r")
source("fastPhenotype.r")
library(tidyverse)
library(forcats)
library(lme4)
library(viridis)
library(ggrepel)

# Palettes
millerstone <- c("#4F6980FF", "#849DB1FF", "#A2CEAAFF", "#638B66FF", "#BFBB60FF", "#F47942FF", "#FBB04EFF", "#B66353FF", "#D7CE9FFF", "#B9AA97FF", "#7E756DFF", "gray70")
antique <- c("#855C75FF", "#D9AF6BFF", "#AF6458FF", "#736F4CFF", "#526A83FF", "#625377FF", "#68855CFF", "#9C9C5EFF", "#A06177FF", "#8C785DFF", "#467378FF", "#7C7C7CFF")

#=============================================================
# 
# 3 major groups to compare
# - PMMD / noPMMD / MMonly +/- depression_only
#
#=============================================================

# PMMD
pmmd_group_dx <- getDf("pmmd_group_dx")
pmmd_ep_dx <- getDf("pmmd_ep_dx")
pmmd_epall_dx <- getDf("pmmd_epall_dx")

# noPMMD
nopmmd_group_dx <- getDf("nopmmd_group_dx")
nopmmd_epall_dx <- getDf("nopmmd_epall_dx")

# MMonly
#all_mania_mixed_dep_count <- getDf("all_mania_mixed_dep_count")
#all_mania_dep <- getDf("all_mania_dep")
#mmonly_group <- all_mania_mixed_dep_count %>%
#  filter(is.na(dep)) %>%
#  replace(is.na(.), 0) %>%
#  mutate(sum = mania + mixed) 
#mmonly_epall_dates.tmp <- all_mania_dep %>% filter(person_id %in% mmonly_group$person_id)

#mmonly_dates_toepisodes.tmp <- mmonly_epall_dates.tmp %>%
#  arrange(person_id, diagnosis_date) %>%
#  group_by(person_id) %>%
#  mutate(
#    # Flag if this is a mania row
#    is_mania = (state == "mm"),
#    # Start a new group when:
#    # 1. State changes from non-mania to mania, OR
#    # 2. Current mania row is more than 1 month after previous mania row
#    new_mania_episode = case_when(
#      # First row of each person
#      row_number() == 1 & is_mania ~ TRUE,
#      # Current row is mania, previous row was not mania
#      is_mania & !lag(is_mania, default = FALSE) ~ TRUE,
#      # Both current and previous are mania, but more than 1 month apart
#      is_mania & lag(is_mania) & (diagnosis_date - lag(diagnosis_date) >= 28) ~ TRUE,
#      # Otherwise, continue the current episode
#      TRUE ~ FALSE
#    ),
#    # Create episode ID by cumulative sum of new_episode flags
#    mania_episode = cumsum(new_mania_episode))

#mmonly_epall <- mmonly_dates_toepisodes.tmp %>%
#  ungroup() %>%
#  filter(is_mania) %>%
#  group_by(person_id, mania_episode) %>%
#  summarise(
#    state = "mm",
#    record_count = n(),
#    diagnosis_date_start = min(diagnosis_date, na.rm = T),
#    diagnosis_date_end = max(diagnosis_date, na.rm = T),
#    duration = diagnosis_date_end - diagnosis_date_start,
#    .groups = "drop"
#  ) %>%
#  rename(diagnosis_date = diagnosis_date_start) %>%
#  arrange(person_id, diagnosis_date)

#mmonly_count_toepisodes <- mmonly_epall %>%
#  group_by(person_id) %>%
#  summarise(mm = n()) %>%
#  mutate(dep = 0) %>%
#  mutate(sum = mm) %>%
#  mutate(group = "MMonly")
#saveDf(mmonly_epall, "mmonly_epall", overwrite = T)
#saveDf(mmonly_count_toepisodes, "mmonly_count_toepisodes", overwrite = T)

# Make person_id > group df
id_group_dx <- data.frame(person_id = c(pmmd_group_dx$person_id, nopmmd_group_dx$person_id, mmonly_group$person_id), 
                       group = c(rep("PMMD", nrow(pmmd_group_dx)), 
                             rep("noPMMD", nrow(nopmmd_group_dx)),
                             rep("MMonly", nrow(mmonly_group))))
saveDf(id_group_dx, "id_group_dx", overwrite = TRUE)

# Identify manic vs mixed episodes (as not using antidepressants to define depression after mixed episodes)
all_mania_dep.tmp <- all_mania_dep %>% mutate(diagnosis_date0 = diagnosis_date)
mm_split <- mania_dep_dates_toepisodes_qc %>%
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
saveDf(mm_split, "mm_split", overwrite = T)
mm_split <- getDf("mm_split")

#------------------------------------------------------------
# Set-up for within-group comparisons
#------------------------------------------------------------

print("Set-up for within-group comparisons")
pmmd_ep_dx.tmp <- pmmd_ep_dx %>% 
  dplyr::select(person_id, diagnosis_date, state) %>% 
  mutate(group = "PMMD") %>% 
  group_by(person_id, state) %>%
  filter(diagnosis_date == min(diagnosis_date)) %>%
  rename(diagnosis_date_PMMD1 = diagnosis_date) %>% 
  ungroup() %>%
  left_join(pmmd_epall_dx, ., relationship = "many-to-many")

pmmd_ep_dx_temporal <- pmmd_ep_dx.tmp %>% 
  group_by(person_id) %>%
  mutate(reltoPMMD1 = case_when(diagnosis_date < diagnosis_date_PMMD1 ~ "before",
                                diagnosis_date == diagnosis_date_PMMD1 ~ "PMMD1",
                                diagnosis_date > diagnosis_date_PMMD1 ~ "after"))
saveDf(pmmd_ep_dx_temporal, "pmmd_ep_dx_temporal", overwrite = T)

#=============================================================
# Medications
#=============================================================

print("Medications")
bipolar_ids_dx <- id_group_dx %>% filter(group %in% c("PMMD", "noPMMD")) %>% pull(person_id)
#meds_bipolar_dx <- getIngredient(as.numeric(bipolar_ids_dx))
#saveDf(meds_bipolar_dx, "meds_bipolar_dx", overwrite = TRUE)

#meds_bipolar_dx <- getDf("meds_bipolar_dx")
#meds_bipolar_dx$person_id <- as.character(meds_bipolar_dx$person_id)

#meds_mania_dx <- getIngredient(as.numeric(id_group_dx$person_id))
#saveDf(meds_mania_dx, "meds_mania_dx", overwrite = TRUE)
#meds_mania_dx <- getDf("meds_mania_dx")
#meds_mania_dx$person_id <- as.character(meds_mania_dx$person_id)

# Extract the psychotropic medications from SQL query
#psych_meds <- getDf("psych_med_list")
#meds_psych_bipolar_dx <- meds_bipolar_dx %>% filter(drug %in% psych_meds$ingredient_concept_name) # 3M to 1M records

#meds_psych_mania_dx <- meds_mania_dx %>% filter(drug %in% psych_meds$ingredient_concept_name) # 6M to 2M records

# - Read in psychotropic drug annotations
print("Extract psychotropic medications")
psych_med_df <- read.csv("psych_med_df.csv", header = T)
psych_med_df <- rbind(psych_med_df, data.frame(drug = "trospium", broad_class = "anticholinergic", narrow_class = ""))

meds_psych_mm_dx_sql <- getDf("meds_psych_mm_dx_sql") %>%
  mutate(person_id = as.character(person_id)) %>%
  filter(start_date >= as.Date(1992-01-01))

meds_psych_bipolar_dx <- meds_psych_mm_dx_sql %>% filter(person_id %in% bipolar_ids_dx) %>% left_join(., psych_med_df, by = "drug")
saveDf(meds_psych_bipolar_dx, "meds_psych_bipolar_dx", overwrite = T)

meds_psych_mania_dx <- meds_psych_mm_dx_sql %>% filter(!person_id %in% bipolar_ids_dx) %>% left_join(., psych_med_df, by = "drug")
saveDf(meds_psych_mania_dx, "meds_psych_mania_dx", overwrite = T)

# Extract non-psychotropic medications
print("Extract non-psychotropic medications - LEFT OUT for code speed")
#meds_nonpsych_bipolar_dx <- meds_bipolar_dx %>% filter(!drug %in% psych_meds$ingredient_concept_name) # 3M to 1M records
#saveDf(meds_nonpsych_bipolar_dx, "meds_nonpsych_bipolar_dx", overwrite = TRUE)

#meds_nonpsych_mania_dx <- meds_mania_dx %>% filter(!drug %in% psych_meds$ingredient_concept_name) # 6M to 4M records
#saveDf(meds_nonpsych_mania_dx, "meds_nonpsych_mania_dx", overwrite = TRUE)

#tmp <- mergePrescriptions(meds_psych_bipolar_dx, gap = 7)

# - check what trazodone and bupropion prescribing trends look like in the depression-only cohort
# - doing this to exclude the possibility that the bipolar group is intentionally being prescribed trazodone and bupropion
# - bottom line is that trazodone and bupropion are the top-prescribed in the MDD datasets
# - helps to justify to not include trazodone and bupropion in capturing more PMMD cases
#icdCodes_dep <- list(dep=c(mdd1, mdd2))
#allPatients_dep <- getPatients(icdCodes_dep) # n=324049
#all_dep <- getDiagnosis(allPatients_dep$person_id, icd_codes = c(dep)) %>% distinct(diagnosis_date, icd_code, .keep_all = TRUE) %>%
#  mutate(person_id = as.character(person_id))
#
#meds_dep <- getIngredient(as.numeric(unique(all_dep$person_id)))
#meds_dep$person_id <- as.character(meds_dep$person_id)
#psych_med_df <- read.csv("psych_med_df.csv", header = T)
#meds_psych_dep_ad <- meds_dep %>% left_join(., psych_med_df, by = "drug") %>% filter(broad_class == "antidepressant")
#tmp <- meds_psych_dep_ad %>% group_by(person_id, drug) %>% summarise(n=n()) %>% ungroup() %>% group_by(drug) %>% summarise(m=n()) %>% arrange(-m)
