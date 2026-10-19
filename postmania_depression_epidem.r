#======================================================
# Libraries and functions
#======================================================

source("Setup_background.r")
source("fastPhenotype.r")
library(forcats)
library(viridis)

# Palettes
millerstone <- c("#4F6980FF", "#849DB1FF", "#A2CEAAFF", "#638B66FF", "#BFBB60FF", "#F47942FF", "#FBB04EFF", "#B66353FF", "#D7CE9FFF", "#B9AA97FF", "#7E756DFF", "gray70")
antique <- c("#855C75FF", "#D9AF6BFF", "#AF6458FF", "#736F4CFF", "#526A83FF", "#625377FF", "#68855CFF", "#9C9C5EFF", "#A06177FF", "#8C785DFF", "#467378FF", "#7C7C7CFF")

#======================================================
# Read in data
#======================================================

source("diagnosis_definitions.r")

# - helpful https://psychcentral.com/bipolar/dsm-5-diagnostic-codes-for-bipolar-disorder#patient-resources

# Diagnostic code inclusion
icdCodes_mania <- list(mm=mm)

# Load population of interest (those with a history of mania)
allPatients_mania <- getPatients(icdCodes_mania) # n=324049 > 346932

#======================================================
# Depression defined using diagnosis code change
#======================================================

# Find individuals with dx of mania and also of depression (n=288870 > 323701)
all_mania_dep <- getDiagnosis(allPatients_mania$person_id, icd_codes = c(mm, dep)) %>% 
  distinct(person_id, diagnosis_date, icd_code, .keep_all = TRUE) %>%
  mutate(person_id = as.character(person_id)) %>%
  mutate(state = ifelse(icd_code %in% mm, "mm", "dep"))
dim(all_mania_dep)

# Clean up cases with both ICD-9 and ICD-10 codes given on the same day
# - if ICD-10 F-code is present, then keep that. Otherwise, take ICD-9
all_mania_dep_rmmultiicd <- remove_icd_dup(all_mania_dep)

# Remove if multi-ICD (n=232527 > 251201)
all_mania_dep <- all_mania_dep %>%
  mutate(tmpid = paste0(person_id, diagnosis_date, icd_code)) %>%
  filter(!tmpid %in% all_mania_dep_rmmultiicd$tmpid) %>% 
  dplyr::select(-c("tmpid"))

# Summarise counts
all_mania_dep_count <- all_mania_dep %>% group_by(person_id, state) %>% summarise(n=n()) %>% pivot_wider(names_from = "state", values_from = "n")
mania_dep_count <- all_mania_dep_count %>% filter(!is.na(mm)) %>% filter(!is.na(dep)) # n=12711
mania_dep_dates <- all_mania_dep %>% filter(person_id %in% mania_dep_count$person_id) # n=107975

# - subdivide by mania vs mixed
all_mania_mixed_dep_count <- all_mania_dep %>% 
  mutate(state = case_when(icd_code %in% c(mania1, mania2) ~ "mania",
                          icd_code %in% c(mixed1, mixed2) ~ "mixed",
                          icd_code %in% dep ~ "dep")) %>%
  group_by(person_id, state) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = "state", values_from = "n")               
mania_mixed_dep_count <- all_mania_mixed_dep_count %>% filter(!is.na(mania) | !is.na(mixed)) %>% filter(!is.na(dep)) %>%
  mutate(mania = ifelse(is.na(mania), 0, mania),
         mixed = ifelse(is.na(mixed), 0, mixed),
         dep = ifelse(is.na(dep), 0, dep))
# - check how many have mania only
mania_mixed_dep_count %>% filter(mania != 0) %>% dim
mania_mixed_dep_count %>% filter(mixed != 0) %>% dim

mania_dep_dates$diagnosis_date <- as.Date(mania_dep_dates$diagnosis_date)

#======================================================
# Pts with depressive switch <1y after mania/mixed
#======================================================

# 1. DEFINE based on diagnosis code switch
# -------------------------------------

# Rearrange the row order where:
# - iterate per-row for an individual
# - find samedayswitch and put the depression entry first (just to standardise)
# - if samedayswitchdep == "dep1", look to the row above
# - switch the row order if the state of the row above is discordant
# - add a row number to maintain order

mania_dep_dates_samedayswitch_count <- mania_dep_dates %>%
  arrange(person_id, diagnosis_date) %>%
  group_by(person_id, diagnosis_date, state) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = "state", values_from = "count") %>%
  filter(!is.na(dep) & !is.na(mm)) %>%
  mutate(total_switch_ep = dep + mm) # n=3852 > 4055

# - collapse those with >2 switches into just 1 mania and 1 dep #n=107136 > 112491
mania_dep_dates_samedayswitch_count_multi <- mania_dep_dates_samedayswitch_count %>% filter(total_switch_ep >= 3)
mania_dep_dates_samedayswitch_multi_collapse <- mania_dep_dates %>% 
  arrange(person_id, diagnosis_date, state) %>%
  left_join(., mania_dep_dates_samedayswitch_count_multi) %>%
  group_by(person_id, diagnosis_date) %>%
  mutate(
    row_count = n(),
    row_num = row_number(),
    is_first = row_num == 1,
    is_last = row_num == row_count) %>%
  filter(
    # Keep all rows in groups with 2 or fewer rows
    row_count <= 2 |
    # OR keep only first and last rows in larger groups
    (row_count > 2 & (is_first | is_last))) %>%
  dplyr::select(-row_count, -row_num, -is_first, -is_last, -dep, -mm)

# - now that there are only a max of 2 switches per day, work with this # n=3852 > 4055
mania_dep_dates_samedayswitch2 <- mania_dep_dates_samedayswitch_count %>%
  filter(total_switch_ep >= 2) %>%
  dplyr::select(c("person_id", "diagnosis_date")) %>%
  mutate(samedayswitch = TRUE)

# - iterate to flip order of rows to minimise amount of potential switches (conservative)
mania_dep_dates_samedayswitch2.tmp <- mania_dep_dates_samedayswitch_multi_collapse %>%
  left_join(.,mania_dep_dates_samedayswitch2) %>%
  arrange(diagnosis_date, state) %>%
  group_by(person_id) %>%
  mutate(rowtmp = row_number()) %>%
  mutate(samedayswitchdep = case_when(samedayswitch == TRUE & state == "dep" ~ "dep1", 
                                      samedayswitch == TRUE & state == "mm" ~ "mm2",
                                      .default = NA)) %>%
  mutate(qswitch = case_when(is.na(samedayswitchdep) ~ "keep0",
                            samedayswitchdep == "dep1" & is.na(lag(state)) ~ "first_ep",
                            samedayswitchdep == "dep1" & lag(state) == "dep" ~ "keep",
                            samedayswitchdep == "dep1" & lag(state) == "mm" & lag(state, n = 2) == "mm" ~ "switch",
                            .default = "to_resolve")) %>%
  mutate(qswitch = case_when(samedayswitchdep == "mm2" & lag(qswitch) == "first_ep" ~ "first_ep",
                             samedayswitchdep == "mm2" & lag(qswitch) == "keep" ~ "keep",
                             samedayswitchdep == "mm2" & lag(qswitch) == "switch" ~ "switch",
                             samedayswitchdep == "mm2" & lag(qswitch) == "to_resolve" ~ "to_resolve",
                            .default = qswitch)) %>%
  ungroup() %>%
  mutate(rowtmp = case_when(qswitch == "switch" & samedayswitchdep == "dep1" ~ rowtmp + 1,
                            qswitch == "switch" & samedayswitchdep == "mm2" ~ rowtmp - 1,
                            .default = rowtmp)) %>%
  arrange(person_id, diagnosis_date, rowtmp) %>%
  mutate(samedayswitch = ifelse(qswitch == "keep0", "keep0", samedayswitch),
         samedayswitch = ifelse(qswitch == "first_ep", "first_ep", samedayswitch),
         samedayswitch = ifelse(qswitch == "keep", "keep", samedayswitch),
         samedayswitch = ifelse(qswitch == "switch", "switched", samedayswitch),
         qswitch = ifelse(qswitch == "keep0", NA, qswitch),
         qswitch = ifelse(qswitch == "first_ep", NA, qswitch),
         qswitch = ifelse(qswitch == "keep", NA, qswitch),
         qswitch = ifelse(qswitch == "switch", NA, qswitch)) %>%
  dplyr::select(-samedayswitchdep)

person_todo <- mania_dep_dates_samedayswitch2.tmp %>% filter(samedayswitch == TRUE) %>% pull(person_id) %>% unique
todo <- mania_dep_dates_samedayswitch2.tmp %>% filter(person_id %in% person_todo)
done <- mania_dep_dates_samedayswitch2.tmp %>% filter(!person_id %in% person_todo)

while (TRUE %in% todo$samedayswitch) {
  print(length(which(todo$samedayswitch == TRUE)))

  todo <- todo %>%
    group_by(person_id) %>% arrange(rowtmp) %>%
    mutate(samedayswitchdep = case_when(samedayswitch == TRUE & state == "dep" ~ "dep1", 
                                      samedayswitch == TRUE & state == "mm" ~ "mm2",
                                      .default = NA)) %>%
    mutate(qswitch = case_when(is.na(samedayswitchdep) ~ NA,
                            samedayswitchdep == "dep1" & lag(state) == "mm" & lag(samedayswitch, n=2) %in% c(NA, "keep0", "keep", "first_ep", "switched") ~ "switch",
                            samedayswitchdep == "dep1" & lag(state) == "dep" ~ "keep",
                            .default = "to_resolve")) %>%
    mutate(qswitch = case_when(samedayswitchdep == "mm2" & lag(qswitch) == "switch" ~ "switch",
                             samedayswitchdep == "mm2" & lag(qswitch) == "keep" ~ "keep",
                             samedayswitchdep == "mm2" & lag(qswitch) == "to_resolve" ~ "to_resolve",
                            .default = qswitch)) %>%
    mutate(rowtmp = case_when(qswitch == "switch" & samedayswitchdep == "dep1" ~ rowtmp + 1,
                            qswitch == "switch" & samedayswitchdep == "mm2" ~ rowtmp - 1,
                            #qswitch == "keep" & samedayswitchdep == "dep1" ~ rowtmp,
                            #qswitch == "keep" & samedayswitchdep == "mm2" ~ rowtmp,
                            .default = rowtmp)) %>%
    mutate(samedayswitch = case_when(qswitch == "switch" ~ "switched",
                                    qswitch == "keep" ~ "keep", 
                                    .default = samedayswitch),
         qswitch = case_when(qswitch == "switch" ~ NA, 
                             qswitch == "keep" ~ NA, 
                             .default = qswitch)) %>%
    dplyr::select(-samedayswitchdep) %>%
    arrange(person_id, diagnosis_date, rowtmp)

    # - find individuals where all samedayswitch != TRUE
  done_add <- todo %>% 
    group_by(person_id, qswitch) %>%
    summarise(n=n()) %>%
    pivot_wider(names_from = qswitch, values_from = n) 

  if ("to_resolve" %in% colnames(done_add)) {
    done_add <- done_add %>% filter(is.na("to_resolve"))
  } else {
    done_add <- done_add
  }
  
  # - separate the rows into those that are accounted for and those that need further iteration
  done <- rbind(done, todo %>% filter(person_id %in% done_add$person_id))
  todo <- todo %>% filter(!person_id %in% done_add$person_id)
  
  print(length(which(todo$samedayswitch == TRUE)))
}

mania_dep_dates_samedayswitch2.tmp <- done # n=2209091

mania_dep_dates_samedayswitch2.tmp %>%
  group_by(person_id) %>%
  summarise(first_ep = length(which(samedayswitch == "first_ep")), switched = length(which(samedayswitch == "switched")), keep = length(which(samedayswitch == "keep"))) %>%
  mutate(sum = first_ep + switched + keep) %>%
  arrange(-sum)

mania_dep_dates_samedayswitch2_clean <- mania_dep_dates_samedayswitch2.tmp %>%
  dplyr::select(-total_switch_ep, -samedayswitch, -qswitch) %>%
  rename(row_order = rowtmp) %>%
  arrange(row_order)

# 0) Collapse consecutive diagnoses if they are adjacent + the same state + occur within 3 months (90 days) of each other
mania_dep_dates_toepisodes.tmp <- mania_dep_dates_samedayswitch2_clean %>%
  arrange(row_order) %>%
  group_by(person_id) %>%
  mutate(
    # Flag if this is a mania row
    is_mania = (state == "mm"),
    # Start a new group when:
    # 1. State changes from non-mania to mania, OR
    # 2. Current mania row is more than 1 month after previous mania row
    new_mm_episode = case_when(
      # First row of each person
      row_number() == 1 & is_mania ~ TRUE,
      # Current row is mania, previous row was not mania
      is_mania & !lag(is_mania, default = FALSE) ~ TRUE,
      # Both current and previous are mania, but more than 1 month apart
      is_mania & lag(is_mania) & (diagnosis_date - lag(diagnosis_date) >= 90) ~ TRUE,
      # Otherwise, continue the current episode
      TRUE ~ FALSE
    ),
    # Create episode ID by cumulative sum of new_episode flags
    mm_episode = cumsum(new_mm_episode)
  ) %>%
  mutate(
    # Flag if this is a dep row
    is_dep = (state == "dep"),
    # Start a new group when:
    # 1. State changes from non-dep to dep, OR
    # 2. Current dep row is more than 1 month after previous dep row
    new_dep_episode = case_when(
      # First row of each person
      row_number() == 1 & is_dep ~ TRUE,
      # Current row is mania, previous row was not mania
      is_dep & !lag(is_dep, default = FALSE) ~ TRUE,
      # Both current and previous are mania, but more than 1 month apart
      is_dep & lag(is_dep) & (diagnosis_date - lag(diagnosis_date) >= 90) ~ TRUE,
      # Otherwise, continue the current episode
      TRUE ~ FALSE
    ),
    # Create episode ID by cumulative sum of new_episode flags
    dep_episode = cumsum(new_dep_episode)
  )

print("mania_collapsed")
mania_collapsed <- mania_dep_dates_toepisodes.tmp %>%
  ungroup() %>%
  filter(is_mania) %>%
  group_by(person_id, mm_episode, dep_episode) %>%
  summarise(
    state = "mm",
    record_count = n(),
    diagnosis_date_start = min(diagnosis_date, na.rm = T),
    diagnosis_date_end = max(diagnosis_date, na.rm = T),
    duration = diagnosis_date_end - diagnosis_date_start,
    .groups = "drop"
  ) %>%
  rename(diagnosis_date = diagnosis_date_start)

print("dep_collapsed")
dep_collapsed <- mania_dep_dates_toepisodes.tmp %>%
  ungroup() %>%
  filter(is_dep) %>%
  group_by(person_id, mm_episode, dep_episode) %>%
  summarise(
    state = "dep",
    record_count = n(),
    diagnosis_date_start = min(diagnosis_date),
    diagnosis_date_end = max(diagnosis_date),
    duration = diagnosis_date_end - diagnosis_date_start,
    .groups = "drop"
  ) %>%
  rename(diagnosis_date = diagnosis_date_start)

# Combine the collapsed mania rows with collapsed depression rows (n=61204 > 53264? relic)
# - given mania_collapsed is added before dep_collapsed, increased number of calling PMMD if going by lag/lead
mania_dep_dates_toepisodes <- bind_rows(mania_collapsed, dep_collapsed) %>% arrange(mm_episode, dep_episode)

# a) Find individuals where mania was followed by depression 
print("earlier_ep")
earlier_ep <- mania_dep_dates_toepisodes %>% 
  group_by(person_id) %>%
  arrange(diagnosis_date) %>%
  mutate(keep = ifelse(state == "mm" & lead(state) == "dep" | (state == "dep" & lag(state) == "mm"), "keep", "remove")) %>%
  filter(keep == "keep") %>%
  dplyr::select(-"keep") %>%
  arrange(person_id)

dim(earlier_ep) #n=20966 > 22168
length(unique(earlier_ep$person_id)) #n=7512 > 8013

# b) Find individuals where mania was followed by depression within 1 year
print("earlier_ep_1y_dx")
earlier_ep_1y_dx <- earlier_ep %>%
  arrange(diagnosis_date) %>%
  group_by(person_id) %>%
  mutate(timetodep_dx = difftime(diagnosis_date, lag(diagnosis_date_end), units = "days")) %>%
  mutate(episodeid = case_when(state == "dep" & timetodep_dx >= 0 & timetodep_dx <= 180 ~ paste0(person_id, "_", state, "_", diagnosis_date, "_", timetodep_dx), 
                              lead(state == "dep" & timetodep_dx >= 0 & timetodep_dx <= 180) ~ paste0(person_id, "_", lead(state), "_", lead(diagnosis_date), "_", lead(timetodep_dx)), 
                              .default = "remove")) %>%
  mutate(timetodep_dx = as.numeric(unlist(lapply(strsplit(episodeid, "_"), function(x) x[4])))) %>%
  filter(episodeid != "remove")

dim(earlier_ep_1y_dx) #n=12822 > 13490
length(unique(earlier_ep_1y_dx$person_id)) #n=4038 > 4283

# c) Check episodes where M-D episode is 0 days and D-M episode is 0 days (ie. samedayswitch)
# - have coded 0 as either meaning sameday switch back to mania, or if that is the end of the EHR
print("mania_dep_dates_toepisodes_pmmd")
mania_dep_dates_toepisodes_pmmd <- mania_dep_dates_toepisodes %>%
  left_join(., earlier_ep_1y_dx %>% mutate(PMMD = "PMMD")) %>%
  arrange(person_id, diagnosis_date) %>%
  group_by(person_id) %>%
  mutate(mm_time = ifelse(state == "dep", diagnosis_date - lag(diagnosis_date), NA),
         rown = row_number(),
         dep_time = case_when(state == "dep" & rown != max(rown) ~ lead(diagnosis_date) - diagnosis_date, 
                              state == "dep" & rown == max(rown) ~ as.difftime(-99, units = "days"),
                             .default = NA)) %>% dplyr::select(-rown)

earlier_ep_1y_dx_deptimethresh <- mania_dep_dates_toepisodes_pmmd %>% #filter(PMMD == "PMMD") %>%
  mutate(dep_time = as.numeric(dep_time))

# Exclusion/QC decisions
# - exclude: episodes/people where there is only 1mm + 1dep in the EHR record, and where the switch date is the same day (can't resolve which is first)
# - sensitivity analysis: excluding samedayswitch episodes
# - see May 6, 2025 milestone for working underlying this decision

# - exclude + re-filter datasets
samedayswitch_only.tmp <- earlier_ep_1y_dx_deptimethresh %>% 
  filter(state == "dep") %>% 
  group_by(person_id) %>% 
  summarise(n_ep = n(), sum_mm_time = sum(mm_time, na.rm = T)) %>%
  filter(sum_mm_time == 0)
samedayswitch_only <- earlier_ep_1y_dx_deptimethresh %>% 
  filter(state == "dep") %>% 
  mutate(sum_mm_time0 = ifelse(person_id %in% samedayswitch_only.tmp$person_id, "sum_mm_time0", "other"))
table(samedayswitch_only$sum_mm_time0) # people kept when sum PMMD mm_time = 0 excluded (n=204 excl): 6907 --> 6089 with 6m threshold --> 4744 when merging codes <90d apart > 4508
length(unique(samedayswitch_only %>% filter(sum_mm_time0 != "sum_mm_time0") %>% pull(person_id))) # 4378 --> 3797 with 6m threshold --> 3704 when merging codes <90d apart > 9092?
samedayswitch_only_sum_mm_time0 <- samedayswitch_only %>% filter(sum_mm_time0 == "sum_mm_time0")

mania_dep_count_1mm1dep_only <- mania_dep_count %>% filter(mm == 1 & dep == 1)
mania_dep_count_samedayswitch_only <- mania_dep_dates_toepisodes_pmmd %>% 
  filter(person_id %in% samedayswitch_only_sum_mm_time0$person_id) %>% 
  group_by(person_id, state) %>% count() %>% pivot_wider(names_from = "state", values_from = "n") %>% mutate(sum = mm + dep)
mania_dep_count_samedayswitch_only_1mm1dep <- mania_dep_count_samedayswitch_only %>% filter(person_id %in% mania_dep_count_1mm1dep_only$person_id)
dim(mania_dep_count_samedayswitch_only_1mm1dep) #n=557 > 2869
earlier_ep_1y_dx_deptimethresh_excl1mm1dep <- earlier_ep_1y_dx_deptimethresh %>% filter(state == "dep") %>% 
  mutate(ep1mm1dep = ifelse(person_id %in% mania_dep_count_samedayswitch_only_1mm1dep$person_id, "only1mm1dep", "multiple"))
# - exclude
excl1mm1dep_sameday <- earlier_ep_1y_dx_deptimethresh_excl1mm1dep %>% 
  filter(person_id %in% mania_dep_count_samedayswitch_only_1mm1dep$person_id) %>% 
  filter(dep_time == 0)

# - sensitivity
exclmm0_sensitivity <- earlier_ep_1y_dx_deptimethresh %>% filter(state == "dep") %>% filter(mm_time == 0)
dim(exclmm0_sensitivity) #260
length(unique(exclmm0_sensitivity$person_id)) #228

# Make dfs for the QC-ed groups
mania_dep_dates_qc <- mania_dep_dates %>% filter(!person_id %in% excl1mm1dep_sameday$person_id)
mania_dep_count_qc <- mania_dep_count %>% filter(!person_id %in% excl1mm1dep_sameday$person_id) %>% mutate(sum = mm+dep)

pmmd_group_dx <- mania_dep_count_qc %>% filter(person_id %in% earlier_ep_1y_dx$person_id)
nopmmd_group_dx <- mania_dep_count_qc %>% filter(!person_id %in% pmmd_group_dx$person_id)

mania_dep_dates_toepisodes_qc <- mania_dep_dates_toepisodes %>% filter(!person_id %in% excl1mm1dep_sameday$person_id)
mania_dep_count_toepisodes_qc <- mania_dep_dates_toepisodes_qc %>% 
  group_by(person_id, state) %>% 
  summarise(n = n()) %>% pivot_wider(names_from = state, values_from = n) %>%
  mutate(sum = mm+dep) %>%
  mutate(group = ifelse(person_id %in% pmmd_group_dx$person_id, "PMMD", "noPMMD"))

# ON NEW RUN: see some people with 300+ episodes total ... Exclude those with >50
excl_50ep <- mania_dep_count_toepisodes_qc %>% filter(sum >= 50)
mania_dep_count_toepisodes_qc <- mania_dep_count_toepisodes_qc %>% filter(!person_id %in% excl_50ep$person_id)
mania_dep_dates_toepisodes_qc <- mania_dep_dates_toepisodes_qc %>% filter(!person_id %in% excl_50ep$person_id)

#earlier_ep_1y_dx <- earlier_ep_1y_dx %>% filter(!person_id %in% excl_50ep$person_id)
pmmd_group_dx <- pmmd_group_dx %>% filter(!person_id %in% excl_50ep$person_id)
nopmmd_group_dx <- nopmmd_group_dx %>% filter(!person_id %in% excl_50ep$person_id)
pmmd_ep_dx <- earlier_ep_1y_dx %>% filter(!person_id %in% excl1mm1dep_sameday$person_id) %>% filter(!person_id %in% excl_50ep$person_id)
pmmd_ep_dx_mm <- pmmd_ep_dx %>% filter(state == "mm")
nopmmd_epall_dx <- mania_dep_dates_toepisodes_qc %>% filter(!person_id %in% pmmd_ep_dx$person_id)
pmmd_epall_dx <- mania_dep_dates_toepisodes_qc %>% 
  filter(person_id %in% pmmd_ep_dx$person_id) %>% 
  left_join(., pmmd_ep_dx %>% dplyr::select(person_id, diagnosis_date, state) %>% mutate(PMMD = "PMMD"))

# Save
print("saveDf")
saveDf(pmmd_group_dx, "pmmd_group_dx", overwrite = T)
saveDf(nopmmd_group_dx, "nopmmd_group_dx", overwrite = T)
saveDf(pmmd_ep_dx, "pmmd_ep_dx", overwrite = T)
saveDf(pmmd_ep_dx_mm, "pmmd_ep_dx_mm", overwrite = T)
saveDf(pmmd_epall_dx, "pmmd_epall_dx", overwrite = T)
saveDf(nopmmd_epall_dx, "nopmmd_epall_dx", overwrite = T)

saveDf(all_mania_dep, "all_mania_dep", overwrite = T)
saveDf(all_mania_mixed_dep_count, "all_mania_mixed_dep_count", overwrite = T)

saveDf(mania_dep_count_qc, "mania_dep_count_qc", overwrite = T)
saveDf(mania_dep_dates_qc, "mania_dep_dates_qc", overwrite = T)
saveDf(mania_dep_dates_toepisodes, "mania_dep_dates_toepisodes", overwrite = T)

saveDf(mania_dep_count_toepisodes_qc, "mania_dep_count_toepisodes_qc", overwrite = T)
saveDf(mania_dep_dates_toepisodes_qc, "mania_dep_dates_toepisodes_qc", overwrite = T)

# Generate comparison groups
source("comparison_group_generate_dx.r")

bipolar_ids_dx <- id_group_dx %>% filter(group %in% c("PMMD", "noPMMD")) %>% pull(person_id)

# 2. DEFINE based on PHQ-9 >= 10
# -------------------------------------

#phq9_mania <- getPHQ9(as.numeric(unique(all_mania_dep$person_id)))
#phq9_mania <- getPHQ9(as.numeric(unique(mania_dep_dates_toepisodes_qc$person_id)))
#saveDf(phq9_mania, "phq9_mania", overwrite = T)
phq9_mania <- getDf("phq9_mm_sql") %>% 
  mutate(person_id = as.character(person_id))

# Filter for depression based on PHQ-9 >= 10
mania_phq9dep <- phq9_mania %>% 
  dplyr::select(c("person_id", "measurement_date", "score")) %>%
  filter(score >= 10)

# Join with mania diagnosis date
mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
#all_mania <- all_mania_dep %>% filter(state == "mm")
mania_dx_phq9dep <- inner_join(mania_dep_dates_toepisodes_qc %>% filter(state == "mm"), mania_phq9dep, relationship = "many-to-many")

# Filter for measurements within 6m of mania diagnosis
# - remove rows if there is an intervening mania recurrence
# - for PHQ-9s with the same person_id and diagnosis_date, take the first PHQ-9 date only
# - for PHQ-9s mapping to multiple measurement dates, take the earliest measurement_date only
mania_dx_phq9dep_1y <- mania_dx_phq9dep %>%
  filter((measurement_date - diagnosis_date_end) %in% 0:180) %>%
  arrange(diagnosis_date) %>%
  group_by(person_id, diagnosis_date) %>%
  filter(measurement_date == min(measurement_date)) %>%
  mutate(timetodep_phq9 = measurement_date - diagnosis_date_end) %>%
  group_by(person_id, measurement_date) %>%
  filter(timetodep_phq9 == min(timetodep_phq9)) %>%
  ungroup() %>%
  mutate(phq9_date = measurement_date)

# - where do these new adds come from?
nrow(mania_dx_phq9dep_1y)
length(unique(mania_dx_phq9dep_1y$person_id))
mania_dx_phq9dep_1y.tmp <- mania_dx_phq9dep_1y %>% dplyr::select(person_id, diagnosis_date)
# - look within PMMD group
pmmd_ep_forantijoin.tmp <- pmmd_ep_dx %>% filter(state =="mm") %>% dplyr::select(person_id, diagnosis_date)
phq9_identified_pmmd <- anti_join(mania_dx_phq9dep_1y.tmp, pmmd_ep_forantijoin.tmp)
phq9_identified_pmmd_jointemporal <- inner_join(pmmd_ep_dx_temporal %>% filter(state == "mm"), phq9_identified_pmmd)
nrow(phq9_identified_pmmd_jointemporal)
length(unique(phq9_identified_pmmd_jointemporal$person_id))
table(phq9_identified_pmmd_jointemporal$reltoPMMD1)
# - look within noPMMD group
nopmmd_id_dx <- id_group_dx %>% filter(group == "noPMMD")
phq9_identified_nopmmd <- mania_dx_phq9dep_1y %>% filter(person_id %in% nopmmd_id_dx$person_id)
dim(phq9_identified_nopmmd)
length(unique(phq9_identified_nopmmd$person_id))

mania_dx_phq9dep_1y_group <- left_join(mania_dx_phq9dep_1y, id_group_dx) %>% 
  left_join(., pmmd_ep_dx_temporal %>% filter(state == "mm") %>% dplyr::select(person_id, diagnosis_date, PMMD, reltoPMMD1)) %>%
  mutate(addpmmd_phq9 = ifelse(is.na(PMMD), "PMMDadd", "PMMD"))
mania_dx_phq9dep_1y_group %>% group_by(group, PMMD, reltoPMMD1) %>% summarise(n=n())
mania_dx_phq9dep_1y_group %>% group_by(person_id, group, PMMD) %>% summarise(n=n()) %>% group_by(group, PMMD) %>% summarise(m=n())

# Find number of all people with mania/mixed diagnosis with PHQ-9 within 1 year of mania
phq9_mania_1y_total <- phq9_mania %>% 
  inner_join(mania_dep_dates_toepisodes_qc %>% filter(state == "mm"), relationship = "many-to-many") %>%
  filter(as.numeric(difftime(measurement_date, diagnosis_date, units = "days"))

# 3. DEFINE based on derived MADRS (imputed from CGI-S)
# -------------------------------------
#madrs_mania <- getDerivedMADRS(as.numeric(id_group_dx$person_id))
#saveDf(madrs_mania, "madrs_mania", overwrite = T)

madrs_mm_sql <- getDf("madrs_mm_sql")

# Filter for depression based on MADRS >= 20
mania_madrsdep <- madrs_mm_sql %>% 
  mutate(person_id = as.character(person_id)) %>%
  filter(score >= 20)

# Join with mania diagnosis date
mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
mania_dx_madrsdep <- inner_join(mania_dep_dates_toepisodes_qc %>% filter(state == "mm"), mania_madrsdep, relationship = "many-to-many")

# Filter for measurements within 1 year of mania diagnosis
# - remove rows if there is an intervening mania recurrence
# - for CGI-S with the same person_id and diagnosis_date, take the first CGI-S date only
# - for CGI-S mapping to multiple measurement dates, take the earliest measurement_date only
mania_dx_madrsdep_1y <- mania_dx_madrsdep %>%
  filter((measurement_date - diagnosis_date_end) %in% 0:180) %>%
  arrange(diagnosis_date) %>%
  group_by(person_id, diagnosis_date) %>%
  filter(measurement_date == min(measurement_date)) %>%
  mutate(timetodep_madrscgis = measurement_date - diagnosis_date_end) %>%
  group_by(person_id, measurement_date) %>%
  filter(timetodep_madrscgis == min(timetodep_madrscgis)) %>%
  ungroup() %>%
  group_by(person_id, diagnosis_date) %>%
  slice_max(score, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(madrscgis_date = measurement_date)

# - where do these new adds come from?
nrow(mania_dx_madrsdep_1y)
length(unique(mania_dx_madrsdep_1y$person_id))
mania_dx_madrsdep_1y.tmp <- mania_dx_madrsdep_1y %>% dplyr::select(person_id, diagnosis_date)
# - look within PMMD group
pmmd_ep_forantijoin.tmp <- pmmd_ep_dx %>% dplyr::select(person_id, diagnosis_date)
madrs_identified_pmmd <- anti_join(mania_dx_madrsdep_1y.tmp, pmmd_ep_forantijoin.tmp)
madrs_identified_pmmd_jointemporal <- inner_join(pmmd_ep_dx_temporal %>% filter(state == "mm"), madrs_identified_pmmd)
nrow(madrs_identified_pmmd_jointemporal)
length(unique(madrs_identified_pmmd_jointemporal$person_id))
table(madrs_identified_pmmd_jointemporal$reltoPMMD1)
# - look within noPMMD group
nopmmd_id_dx <- id_group_dx %>% filter(group == "noPMMD")
madrs_identified_nopmmd <- madrs_identified_pmmd %>% filter(person_id %in% nopmmd_id_dx$person_id)
dim(madrs_identified_nopmmd)
length(unique(madrs_identified_nopmmd$person_id))

mania_dx_madrsdep_1y_group <- left_join(mania_dx_madrsdep_1y, id_group_dx) %>% 
  left_join(., pmmd_ep_dx_temporal %>% filter(state == "mm") %>% dplyr::select(person_id, diagnosis_date, PMMD, reltoPMMD1)) %>%
  mutate(addpmmd_madrscgis = ifelse(is.na(PMMD), "PMMDadd", "PMMD"))
mania_dx_madrsdep_1y_group %>% group_by(group, PMMD, reltoPMMD1) %>% summarise(n=n())
mania_dx_madrsdep_1y_group %>% group_by(person_id, group, PMMD) %>% summarise(n=n()) %>% group_by(group, PMMD) %>% summarise(m=n())

# 4. DEFINE based on QIDS
# -------------------------------------
#qids_mania <- getQIDS(as.numeric(id_group_dx$person_id))
#saveDf(qids_mania, "qids_mania", overwrite = T)

qids_mm_sql <- getDf("qids_mm_sql")

# Filter for depression based on QIDS >= 11
mania_qidsdep <- qids_mm_sql %>% 
  mutate(person_id = as.character(person_id)) %>%
  filter(score >= 11)

# Join with mania diagnosis date
mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")
mania_dx_qidsdep <- inner_join(mania_dep_dates_toepisodes_qc %>% filter(state == "mm"), mania_qidsdep, relationship = "many-to-many")

# Filter for measurements within 1 year of mania diagnosis
# - remove rows if there is an intervening mania recurrence
# - for QIDS with the same person_id and diagnosis_date, take the first CGI-S date only
# - for QIDS mapping to multiple measurement dates, take the earliest measurement_date only
mania_dx_qidsdep_1y <- mania_dx_qidsdep %>%
  filter((measurement_date - diagnosis_date_end) %in% 0:180) %>%
  arrange(diagnosis_date) %>%
  group_by(person_id, diagnosis_date) %>%
  filter(measurement_date == min(measurement_date)) %>%
  mutate(timetodep_qids = measurement_date - diagnosis_date_end) %>%
  group_by(person_id, measurement_date) %>%
  filter(timetodep_qids == min(timetodep_qids)) %>%
  ungroup() %>%
  group_by(person_id, diagnosis_date) %>%
  slice_max(score, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(qids_date = measurement_date)

# - where do these new adds come from?
nrow(mania_dx_qidsdep_1y)
length(unique(mania_dx_qidsdep_1y$person_id))
mania_dx_qidsdep_1y.tmp <- mania_dx_qidsdep_1y %>% dplyr::select(person_id, diagnosis_date)
# - look within PMMD group
pmmd_ep_forantijoin.tmp <- pmmd_ep_dx %>% dplyr::select(person_id, diagnosis_date)
qids_identified_pmmd <- anti_join(mania_dx_qidsdep_1y.tmp, pmmd_ep_forantijoin.tmp)
qids_identified_pmmd_jointemporal <- inner_join(pmmd_ep_dx_temporal %>% filter(state == "mm"), qids_identified_pmmd)
nrow(qids_identified_pmmd_jointemporal)
length(unique(qids_identified_pmmd_jointemporal$person_id))
table(qids_identified_pmmd_jointemporal$reltoPMMD1)
# - look within noPMMD group
nopmmd_id_dx <- id_group_dx %>% filter(group == "noPMMD")
qids_identified_nopmmd <- qids_identified_pmmd %>% filter(person_id %in% nopmmd_id_dx$person_id)
dim(qids_identified_nopmmd)
length(unique(qids_identified_nopmmd$person_id))

mania_dx_qidsdep_1y_group <- left_join(mania_dx_qidsdep_1y, id_group_dx) %>% 
  left_join(., pmmd_ep_dx_temporal %>% filter(state == "mm") %>% dplyr::select(person_id, diagnosis_date, PMMD, reltoPMMD1)) %>%
  mutate(addpmmd_qids = ifelse(is.na(PMMD), "PMMDadd", "PMMD"))
mania_dx_qidsdep_1y_group %>% group_by(group, PMMD, reltoPMMD1) %>% summarise(n=n())
mania_dx_qidsdep_1y_group %>% group_by(person_id, group, PMMD) %>% summarise(n=n()) %>% group_by(group, PMMD) %>% summarise(m=n())

# 5. DEFINE based on antidepressant prescription (manic episodes only)
# -------------------------------------

meds_psych_mania_dx <- getDf("meds_psych_mania_dx")
mania_dep_dates_toepisodes_qc <- getDf("mania_dep_dates_toepisodes_qc")

meds_psych_mania_ad <- meds_psych_bipolar_dx %>% filter(broad_class == "antidepressant") %>% mutate(person_id = as.character(person_id)) # 500K
meds_psych_mania_ad_ep <- inner_join(mania_dep_dates_toepisodes_qc %>% filter(state == "mm"), 
  meds_psych_mania_ad, relationship = "many-to-many")

# - identify added episodes
meds_psych_mania_ad_ep_1y <- meds_psych_mania_ad_ep %>%
  left_join(., mm_split) %>% # EXCLUDE those with mixed episode (as antidepressant may be for this)
  filter(mm_split == "mania") %>%
  filter((start_date - diagnosis_date_end) %in% 0:180) %>%
  group_by(person_id, diagnosis_date) %>%
  filter(start_date == min(start_date)) %>% ungroup() %>%
  mutate(timetodep_ad = start_date - diagnosis_date_end) %>%
  group_by(person_id, start_date) %>%
  filter(timetodep_ad == min(timetodep_ad)) %>%
  ungroup()

# - check what meds are being prescribed for added episodes (trazodone 3x the next of bupropion)
meds_psych_mania_ad_ep_1y %>% group_by(drug) %>% summarise(n=n()) %>% arrange(-n)

# - exclude trazodone (have also checked that they aren't enriched for bipolar vs MDD only)
meds_psych_mania_ad_ep_1y <- meds_psych_mania_ad_ep_1y %>% filter(!drug %in% c("trazodone"))

# - where do these new adds come from?
nrow(meds_psych_mania_ad_ep_1y)
length(unique(meds_psych_mania_ad_ep_1y$person_id))
meds_psych_mania_ad_ep_1y.tmp <- meds_psych_mania_ad_ep_1y %>% dplyr::select(person_id, diagnosis_date)
# - look within PMMD group
pmmd_ep_forantijoin.tmp <- pmmd_ep_dx %>% dplyr::select(person_id, diagnosis_date)
ad_identified_pmmd <- anti_join(meds_psych_mania_ad_ep_1y.tmp, pmmd_ep_forantijoin.tmp)
ad_identified_pmmd_jointemporal <- inner_join(pmmd_ep_dx_temporal %>% filter(state == "mm"), ad_identified_pmmd)
nrow(ad_identified_pmmd_jointemporal)
length(unique(ad_identified_pmmd_jointemporal$person_id))
table(ad_identified_pmmd_jointemporal$reltoPMMD1)
# - look within noPMMD group
nopmmd_id_dx <- id_group_dx %>% filter(group == "noPMMD")
ad_identified_nopmmd <- ad_identified_pmmd %>% filter(person_id %in% nopmmd_id_dx$person_id)
dim(ad_identified_nopmmd)
length(unique(ad_identified_nopmmd$person_id))

mania_dx_addep_1y_group <- left_join(meds_psych_mania_ad_ep_1y, id_group_dx) %>% 
  left_join(., pmmd_ep_dx_temporal %>% filter(state == "mm") %>% dplyr::select(person_id, diagnosis_date, PMMD, reltoPMMD1)) %>%
  mutate(addpmmd_ad = ifelse(is.na(PMMD), "PMMDadd", "PMMD"))
mania_dx_addep_1y_group %>% group_by(group, PMMD, reltoPMMD1) %>% summarise(n=n())
mania_dx_addep_1y_group %>% group_by(person_id, group, PMMD) %>% summarise(n=n()) %>% group_by(group, PMMD) %>% summarise(m=n())

# JOIN together
# -------------------------------------

mania_dx_phq9dep_1y_group.tmp <- mania_dx_phq9dep_1y_group %>% dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, timetodep_phq9, group, phq9_date, addpmmd_phq9)
mania_dx_madrsdep_1y_group.tmp <- mania_dx_madrsdep_1y_group %>% dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, timetodep_madrscgis, group, madrscgis_date, addpmmd_madrscgis)
mania_dx_qidsdep_1y_group.tmp <- mania_dx_qidsdep_1y_group %>% dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, timetodep_qids, group, qids_date, addpmmd_qids)
mania_dx_addep_1y_group.tmp <- mania_dx_addep_1y_group %>% dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, timetodep_ad, group, start_date, addpmmd_ad)
pmmd_ep_dx.tmp <- pmmd_ep_dx %>% 
  filter(state == "mm") %>% 
  dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, episodeid) %>%
  mutate(dep_date_dx = as.Date(unlist(lapply(strsplit(episodeid, "_"), function(x) x[3])))) %>%
  dplyr::select(-episodeid) %>%
  mutate(timetodep_dx = dep_date_dx - diagnosis_date_end)
pmmd_ep_dx_temporal.tmp <- pmmd_ep_dx_temporal %>% dplyr::select(c(person_id, state, diagnosis_date, diagnosis_date_end, reltoPMMD1))
mania_dep_dates_toepisodes_qc.tmp <- mania_dep_dates_toepisodes_qc %>% 
  filter(state == "mm") %>% 
  left_join(., id_group_dx)

mania_dep_dates_toepisodes_qc_defbroad <- mania_dep_dates_toepisodes_qc.tmp %>%
  left_join(., pmmd_ep_dx.tmp) %>% left_join(., pmmd_ep_dx_temporal.tmp) %>%
  left_join(., mania_dx_phq9dep_1y_group.tmp) %>%
  left_join(., mania_dx_madrsdep_1y_group.tmp) %>%
  left_join(., mania_dx_qidsdep_1y_group.tmp) %>%
  left_join(., mania_dx_addep_1y_group.tmp %>% distinct) %>%
  distinct() %>%
  mutate(addpmmd = ifelse(is.na(dep_date_dx) & (!is.na(addpmmd_phq9) | !is.na(addpmmd_madrscgis) | !is.na(addpmmd_qids) | !is.na(addpmmd_ad)), "PMMDadd", NA)) %>%
  rowwise() %>%
  mutate(dep_date = pmin(dep_date_dx, phq9_date, madrscgis_date, qids_date, start_date, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(timetodep = difftime(dep_date, diagnosis_date_end, units="days")) %>%
  mutate(supersede = ifelse(dep_date < dep_date_dx, "yes", "no"))

mania_dep_dates_toepisodes_qc_defbroad_noad <- mania_dep_dates_toepisodes_qc.tmp %>%
  left_join(., pmmd_ep_dx.tmp) %>% left_join(., pmmd_ep_dx_temporal.tmp) %>%
  left_join(., mania_dx_phq9dep_1y_group.tmp) %>%
  left_join(., mania_dx_madrsdep_1y_group.tmp) %>%
  left_join(., mania_dx_qidsdep_1y_group.tmp) %>%
  distinct() %>%
  mutate(addpmmd = ifelse(is.na(dep_date_dx) & (!is.na(addpmmd_phq9) | !is.na(addpmmd_madrscgis) | !is.na(addpmmd_qids)), "PMMDadd", NA)) %>%
  rowwise() %>%
  mutate(dep_date = pmin(dep_date_dx, phq9_date, madrscgis_date, qids_date, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(timetodep_months = as.numeric(difftime(dep_date, diagnosis_date_end, units="days"))/30) %>%
  mutate(supersede = ifelse(dep_date < dep_date_dx, "yes", "no"))
saveDf(mania_dep_dates_toepisodes_qc_defbroad_noad, "mania_dep_dates_toepisodes_qc_defbroad_noad", overwrite = T)

# - added episodes total
mania_dep_dates_toepisodes_qc_defbroad %>% group_by(group, addpmmd, supersede) %>% summarise(n=n())
# - added individuals total
length(unique(mania_dep_dates_toepisodes_qc_defbroad %>% filter(group == "noPMMD") %>% filter(addpmmd == "PMMDadd") %>% pull(person_id)))
mania_dep_dates_toepisodes_qc_defbroad %>% group_by(person_id, group, addpmmd, supersede) %>% summarise(n=n()) %>% group_by(group, addpmmd, supersede) %>% summarise(m=n())

#summary(as.numeric(mania_dx_phq9dep_cgisdep_1y_group$timelagdx))
#mania_dx_phq9dep_cgisdep_1y_group %>% filter(timelagdx == 417) %>% data.frame
#mania_dx_phq9dep_cgisdep_1y_group %>% group_by(addpmmd_phq9, addpmmd_cgis) %>% summarise(n=n())

# - plot for differences in time to diagnosis
mania_dep_dates_toepisodes_qc_defbroad_timetodep <- mania_dep_dates_toepisodes_qc_defbroad %>%
  pivot_longer(cols = c("timetodep_dx", "timetodep_phq9", "timetodep_madrscgis", "timetodep_qids", "timetodep_ad"), names_to = "method", values_to = "time") %>%
  filter(!is.na(time)) %>%
  filter(!is.infinite(time)) %>%
  mutate(time = as.numeric(time, units="days")) %>%
  mutate(method = gsub("timetodep_", "", method)) %>%
  mutate(method = factor(method, levels = c("dx", "phq9", "madrscgis", "qids", "ad"))) %>%
  arrange(method)

ggplot(mania_dep_dates_toepisodes_qc_defbroad_timetodep, aes(x = time, colour = method)) +
  geom_density() +
  theme_bw() +
  scale_colour_manual(values = millerstone[c(4,1,6,3,2)]) +
  theme(legend.position = "bottom")
ggsave("timetodep_density_dxphq9madrsqidsad.png", height = 5, width = 5)

ggplot(mania_dep_dates_toepisodes_qc_defbroad_timetodep, aes(x = time, colour = method)) +
  stat_ecdf(geom = "step") +  
  theme_bw() +
  scale_colour_manual(values = millerstone[c(4,1,6,3,2)]) +
  theme(legend.position = "bottom")
ggsave("timetodep_ecdf_dxphq9madrsqidsad.png", height = 5, width = 5)

# Diagnosis delay vs ICD codes
mania_dep_dates_toepisodes_qc_defbroad_timelagdx <- mania_dep_dates_toepisodes_qc_defbroad %>%
  mutate(adstart_date = start_date) %>%
  pivot_longer(cols = c("phq9_date", "madrscgis_date", "qids_date", "adstart_date"), names_to = "method", values_to = "dep_date_other") %>%
  filter(!is.na(dep_date_other)) %>%
  mutate(icd_vs = gsub("_date", "", method)) %>%
  mutate(timelagdx = as.numeric(difftime(dep_date_dx, dep_date_other, units = "days"))) %>%
  mutate(icd_vs = factor(icd_vs, levels = c("phq9", "madrscgis", "qids", "adstart")))

mania_dep_dates_toepisodes_qc_defbroad_timelagdx_summary <- mania_dep_dates_toepisodes_qc_defbroad_timelagdx %>% 
  left_join(., mm_split) %>%
  group_by(icd_vs, mm_split) %>%
  summarise(n=length(which(!is.na(timelagdx))), 
              n_dx_delay = length(which(timelagdx > 0)), 
              prop_dx_delay = n_dx_delay/n, 
              med_delay = median(timelagdx, na.rm = T), 
              q25_75_delay = paste0("[", quantile(timelagdx, 0.25, na.rm = T), "-", quantile(timelagdx, 0.75, na.rm = T), "]"), 
              mean_delay = mean(timelagdx, na.rm = T), 
              sd = sd(timelagdx, na.rm = T),
              n_7d_delay = length(which(timelagdx >= 7)),
              n_14d_delay = length(which(timelagdx >= 14)),
              n_28d_delay = length(which(timelagdx >= 28))) %>%
  arrange(mm_split, icd_vs)
saveDf(mania_dep_dates_toepisodes_qc_defbroad_timelagdx_summary, "mania_dep_dates_toepisodes_qc_defbroad_timelagdx_summary", overwrite = T)

# - plot for diagnosis delay
ggplot(mania_dep_dates_toepisodes_qc_defbroad_timelagdx %>% filter(icd_vs != "qids") %>% left_join(., mm_split), aes(y = timelagdx, x = icd_vs, fill = icd_vs)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = millerstone[c(1,6,2)]) +
  theme(legend.position = "none") +
  ylab("Time lag of ICD diagnosis code") + xlab("Case-finding alternative") +
  facet_wrap(~mm_split)
ggsave("timelagdx_boxplot_phq9vmadrsvsqidsvsadstart.png", height = 4, width = 5)

ggplot(mania_dep_dates_toepisodes_qc_defbroad_timelagdx, aes(x = timelagdx, colour = icd_vs)) +
  stat_ecdf(geom = "step") +  
  theme_bw() +
  scale_colour_manual(values = millerstone[c(1,6,3,2)]) +
  theme(legend.position = "bottom")
ggsave("timelagdx_ecdf_phq9vmadrsvsqidsvsadstart.png", height = 4, width = 4)

saveDf(mania_dep_dates_toepisodes_qc_defbroad, "mania_dep_dates_toepisodes_qc_defbroad", overwrite = TRUE)
saveDf(mania_dep_dates_toepisodes_qc_defbroad_timelagdx, "mania_dep_dates_toepisodes_qc_defbroad_timelagdx", overwrite = T)

# Extract broadly-identified depression as as depression "diagnosis" (n=10441 episodes -> n=10969 episodes with icd code update)
mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest0 <- mania_dep_dates_toepisodes_qc_defbroad_timetodep %>%
  group_by(person_id, diagnosis_date) %>%
  arrange(dep_date) %>%
  slice_head(n=1) %>%
  #slice_min(dep_date, n=1, with_ties = FALSE) %>% ##### WITH THIS, LOSE n=448 ...
  ungroup() %>%
  mutate(method = as.character(method))

pmmd_broaddx <- mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest0 %>% 
  dplyr::select(c("person_id", "dep_date", "method", "group", "addpmmd", "supersede")) %>%
  filter(method != "dx") %>%
  mutate(state = "dep") %>%
  rename(diagnosis_date = dep_date) %>%
  mutate(diagnosis_date_end = diagnosis_date)

#------------------------------------------
# Update important files
#------------------------------------------
# Update id_group
add_pmmd_group <- pmmd_broaddx %>% filter(group == "noPMMD" & addpmmd == "PMMDadd") %>% dplyr::select("person_id", "addpmmd") %>% distinct()
id_group0 <- id_group_dx %>% 
  mutate(group = ifelse(person_id %in% add_pmmd_group$person_id, "PMMD", group)) %>%
  left_join(., add_pmmd_group)
saveDf(id_group0, "id_group0", overwrite = T)

# Update bipolar_episodes_count
add_pmmd_ep_count <- pmmd_broaddx %>% filter(addpmmd == "PMMDadd") %>% group_by(person_id) %>% summarise(addpmmd=n())
bipolar_episodes_count <- mania_dep_count_toepisodes_qc %>% 
  left_join(., add_pmmd_ep_count) %>%
  rename(dep0 = dep) %>%
  mutate(dep = ifelse(!is.na(addpmmd), dep0 + addpmmd, dep0)) %>%
  mutate(sum = mm + dep) %>%
  dplyr::select(person_id, mm, dep, sum, addpmmd) %>%
  left_join(., id_group0 %>% dplyr::select(-addpmmd))
saveDf(bipolar_episodes_count, "bipolar_episodes_count", overwrite = T)

# Plot number of counts per person
bipolar_episodes_count <- getDf("bipolar_episodes_count")
bipolar_episodes_count.long <- bipolar_episodes_count %>%
  arrange(mm, dep) %>%
  mutate(mm = mm*-1) %>%
  mutate(id = factor(person_id, levels = person_id)) %>%
  pivot_longer(., cols = c("mm", "dep"), values_to = "episode_count_pp", names_to = "ep_type") %>%
  mutate(ep_type = ifelse(ep_type == "mm", "mania_mixed", "depression"))

ggplot(bipolar_episodes_count.long, aes(x = episode_count_pp, y = id, fill = ep_type)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  scale_fill_manual(values = millerstone[c(2,4,8)]) +
  coord_cartesian(xlim = c(-32,30)) + scale_x_continuous(breaks = seq(-35, 30, by = 5))
ggsave("bipolar_episodes_countbarplot.png", height = 17, width = 10)

# Update pmmd_ep_temporal
# 0. 
# - identify all episodes occurring on the same day 
# - why: use when adding broaddef PMMD episodes which may have picked up dep/mm switched around
# - removes n=1584 episodes from 1489 individuals in mania_dep_dates_toepisodes_qc_defbroad_timetodep
samedayswitch_anyorder <- mania_dep_dates_toepisodes %>% 
  group_by(person_id, diagnosis_date) %>%
  summarise(n=n()) %>%
  filter(n>=2) %>% ungroup
saveDf(samedayswitch_anyorder, "samedayswitch_anyorder", overwrite = T)

# - add excluded samedayswitch_anyorder with unresolvable order to id_group0
id_group0_addsamedayswitch <- id_group0 %>% 
  mutate(samedayswitch_anyorder = ifelse(person_id %in% samedayswitch_anyorder$person_id, "yes", NA))

# Update pmmd_group and nopmmd_group with updated groupings
pmmd_group <- id_group0_addsamedayswitch %>% filter(group == "PMMD") %>% filter(is.na(samedayswitch_anyorder))
nopmmd_group <- id_group0_addsamedayswitch %>% filter(group == "noPMMD") %>% filter(is.na(samedayswitch_anyorder))

# 1. find records/people to supersede for people already in PMMD group
# - a) simple updating of existing record
mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest <- mania_dep_dates_toepisodes_qc_defbroad_timetodep %>%
  group_by(person_id, diagnosis_date) %>%
  #slice_head(n=1) %>%
  slice_min(time, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(method = as.character(method)) %>%
  anti_join(., samedayswitch_anyorder %>% dplyr::select(-n))
saveDf(mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest, "mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest", overwrite = T)

addpmmd_n <- mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest %>%  # previously addpmmd_n
  group_by(supersede, addpmmd, method) %>%
  summarise(n=n()) %>%
  mutate(supersede = ifelse(is.na(supersede), "no", supersede)) %>%
  mutate(addpmmd = ifelse(is.na(addpmmd), "no", addpmmd))
#supersede_n <- mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest %>% 
#  pivot_longer(cols = c("phq9_date", "madrscgis_date", "qids_date", "start_date"), names_to = "tmp", values_to = "tmp_date") %>%
#  filter(!is.na(tmp_date) & supersede == "yes") %>%
#  arrange(tmp_date) %>%
#  group_by(person_id, diagnosis_date) %>%
#  slice_head(n=1) %>%
#  mutate(method = case_when(tmp == "phq9_date" ~ "phq9",
#                         tmp == "madrscgis_date" ~ "madrscgis",
#                         tmp == "qids_date" ~ "qids",
#                         tmp == "start_date" ~ "ad")) %>%
#  group_by(method) %>%
#  summarise(n=n()) %>%
#  mutate(supersede = "yes", addpmmd = "no")

addpmmd_supersede <- addpmmd_n %>% pivot_wider(names_from = "method", values_from = "n")
saveDf(addpmmd_supersede, "addpmmd_supersede", overwrite = T)

supersede_update <- mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest %>% 
  filter(supersede == "yes") %>% 
  mutate(supersede = "supersede") %>%
  filter(is.na(addpmmd)) %>%
  rename(method1 = method) %>%
  rename(dep_date1 = dep_date) %>%
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, dep_date1, timetodep, method1, supersede)
   # extract the updates in this category

pmmd_ep_temporal_supersede_update.tmp <- pmmd_ep_dx_temporal %>% 
  inner_join(., supersede_update) %>%
  dplyr::select(person_id, mm_episode, dep_episode, dep_date1, timetodep, method1, supersede) %>%
  mutate(dep_episode = dep_episode + 1) # - update is currently relative to the mm episode, so then find the matching dep_episode

pmmd_ep_temporal_supersede_update <- pmmd_ep_dx_temporal %>%
  left_join(., pmmd_ep_temporal_supersede_update.tmp) %>%
  mutate(diagnosis_date = dplyr::if_else(!is.na(method1), as.Date(dep_date1), as.Date(diagnosis_date))) %>%
  mutate(duration = as.numeric(difftime(diagnosis_date_end, dep_date1, units = "days"))) %>%
  mutate(diagnosis_date_PMMD1 = if_else(reltoPMMD1 == "PMMD1" & !is.na(method1), dep_date1, diagnosis_date_PMMD1))

# - b) addPMMD within the PMMD group
# - this involves the following:
# - updating cols PMMD for the mm PMMD row
# - adding a corresponding dep PMMD row
# - changing reltoPMMD1, diagnosis_date_PMMD1

# - take the added PMMD episodes
addpmmd_updatemm <- mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest %>% 
  filter(group == "PMMD") %>%
  filter(!is.na(addpmmd)) %>%
  mutate(PMMD1 = "PMMD", group = "PMMD") %>%
  rename(method2 = method) %>%
  rename(dep_date2 = dep_date) %>%
  dplyr::select(person_id, diagnosis_date, PMMD1, dep_date2, addpmmd, method2)

addpmmd_updatedep <- addpmmd_updatemm %>%
  dplyr::select(person_id, dep_date2, method2) %>%
  mutate(diagnosis_date = dep_date2) %>%
  mutate(state = "dep") %>%
  mutate(PMMD = "PMMD") %>%
  mutate(group = "PMMD") %>%
  mutate(addpmmd = "PMMDadd")
  
# - join to pmmd_ep_dx_temporal
pmmd_ep_temporal_supersede_addpmmd_update <- pmmd_ep_temporal_supersede_update %>%
  left_join(., addpmmd_updatemm) %>%
  full_join(., addpmmd_updatedep)
# - how many episodes added? (Add n=1398)
nrow(pmmd_ep_temporal_supersede_addpmmd_update) - nrow(pmmd_ep_temporal_supersede_update)

# - clean up:
# - consolidate columns
pmmd_ep_temporal_supersede_addpmmd_update2 <- pmmd_ep_temporal_supersede_addpmmd_update %>%
  mutate(PMMD = ifelse(is.na(PMMD), PMMD1, PMMD)) %>%
  mutate(dep_date1 = ifelse(is.na(dep_date1), dep_date2, dep_date1)) %>% rename(dep_date = dep_date1) %>%
  mutate(method1 = ifelse(is.na(method1), method2, method1)) %>% rename(method = method1) %>%
  mutate(addpmmd = ifelse(is.na(addpmmd), supersede, addpmmd)) %>%
  dplyr::select(-c("PMMD1", "dep_date", "timetodep", "dep_date2", "method2", "supersede", "timetodep"))
# - update diagnosis_date_PMMD1
foo <- pmmd_ep_temporal_supersede_addpmmd_update2 %>% 
  filter(PMMD == "PMMD") %>%
  arrange(diagnosis_date) %>%
  group_by(person_id, state) %>%
  slice_min(diagnosis_date) %>%
  ungroup() %>%
  dplyr::select(person_id, state, diagnosis_date) %>%
  rename(diagnosis_date_PMMD1 = diagnosis_date) %>% distinct
# - add back to the original df and update relative timing
pmmd_ep_temporal_supersede_addpmmd_update3 <- pmmd_ep_temporal_supersede_addpmmd_update2 %>%
  dplyr::select(-c("diagnosis_date_PMMD1", "reltoPMMD1")) %>%
  left_join(., foo) %>%
  mutate(reltoPMMD1 = case_when(diagnosis_date < diagnosis_date_PMMD1 ~ "before",
                                diagnosis_date == diagnosis_date_PMMD1 ~ "PMMD1",
                                diagnosis_date > diagnosis_date_PMMD1 ~ "after"))
# - renumber episodes
pmmd_oldpmmd_ep_temporal <- pmmd_ep_temporal_supersede_addpmmd_update3 %>%
  arrange(diagnosis_date, mm_episode, dep_episode) %>%
  group_by(person_id) %>%
  mutate(is_mm = ifelse(state == "mm", 1, 0)) %>%
  mutate(is_dep = ifelse(state == "dep", 1, 0)) %>%
  mutate(mm_episode = cumsum(is_mm)) %>%
  mutate(dep_episode = cumsum(is_dep)) %>%
  ungroup() %>%
  dplyr::select(-c("is_mm", "is_dep")) %>%
  mutate(duration = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days")))

# - c) addPMMD within the noPMMD group
# - get the relevant episodes per reassigned id
nopmmd_addmm <- mania_dep_dates_toepisodes_qc_defbroad_timetodep_shortest %>% 
  filter(group == "noPMMD") %>%
  dplyr::select(person_id, state, diagnosis_date, diagnosis_date_end, dep_date, method) %>%
  mutate(PMMD = "PMMD", group = "PMMD") %>%
  mutate(addpmmd = "PMMDadd")

nopmmd_adddep <- nopmmd_addmm %>% 
  dplyr::select(-c("diagnosis_date", "diagnosis_date_end")) %>%
  mutate(state = "dep") %>%
  rename(diagnosis_date = dep_date) %>%
  mutate(addpmmd = "PMMDadd") %>%
  mutate(dep_date = diagnosis_date)

nopmmd_convert_epall <- nopmmd_epall_dx %>% filter(person_id %in% unique(nopmmd_addmm$person_id))

# - add nopmmd_add to nopmmd_convert_epall
pmmd_addnopmmd_ep_temporal.tmp <- nopmmd_convert_epall %>%
  left_join(., nopmmd_addmm) %>%
  full_join(., nopmmd_adddep) %>% 
  arrange(person_id, diagnosis_date) %>%
  mutate(group = "PMMD")

pmmd_addnopmmd_ep_temporal_pmmdep <- pmmd_addnopmmd_ep_temporal.tmp %>%
  dplyr::select(person_id, diagnosis_date, state, PMMD) %>%
  arrange(diagnosis_date) %>%
  mutate(group = "PMMD")

# - add diagnosis_date_PMMD1
pmmd_addnopmmd_ep_temporal_pmmd1 <- pmmd_addnopmmd_ep_temporal_pmmdep %>% 
  filter(PMMD == "PMMD") %>%
  group_by(person_id, state) %>%
  slice_min(diagnosis_date) %>%
  rename(diagnosis_date_PMMD1 = diagnosis_date) %>%
  ungroup() %>%
  dplyr::select(person_id, state, diagnosis_date_PMMD1) %>%
  left_join(pmmd_addnopmmd_ep_temporal.tmp, ., relationship = "many-to-many")

# - add timing
pmmd_addnopmmd_ep_temporal <- pmmd_addnopmmd_ep_temporal_pmmd1 %>% 
  #group_by(person_id) %>%
  mutate(reltoPMMD1 = case_when(diagnosis_date < diagnosis_date_PMMD1 ~ "before",
                                diagnosis_date == diagnosis_date_PMMD1 ~ "PMMD1",
                                diagnosis_date > diagnosis_date_PMMD1 ~ "after")) %>%
  #arrange(diagnosis_date, mm_episode, dep_episode) %>%
  mutate(is_mm = ifelse(state == "mm", 1, 0)) %>%
  mutate(is_dep = ifelse(state == "dep", 1, 0)) %>%
  group_by(person_id) %>%
  mutate(mm_episode = cumsum(is_mm)) %>%
  mutate(dep_episode = cumsum(is_dep)) %>%
  ungroup() %>%
  dplyr::select(-c("is_mm", "is_dep")) %>%
  mutate(duration = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days")))

saveDf(pmmd_addnopmmd_ep_temporal, "pmmd_addnopmmd_ep_temporal", overwrite = T)

# d) combine all episodes and then get a revised order
pmmd_ep_temporal <- full_join(pmmd_oldpmmd_ep_temporal, pmmd_addnopmmd_ep_temporal) %>%
  dplyr::select(-dep_date) %>%
  filter(!person_id %in% samedayswitch_anyorder$person_id)

# Update pmmd_ep
pmmd_ep <- pmmd_ep_temporal %>% 
  arrange(mm_episode, dep_episode) %>%
  filter(PMMD == "PMMD") %>%
  group_by(person_id, state) %>%
  mutate(indpmmd = row_number()) %>%
  ungroup() %>%
  group_by(person_id, indpmmd) %>%
  mutate(timebwdep = difftime(lead(diagnosis_date), diagnosis_date_end, units = "days")) %>%
  mutate(timetodep = difftime(lead(diagnosis_date), diagnosis_date, units = "days")) %>%
  mutate(dep_date = lead(diagnosis_date)) %>%
  ungroup()

# Update pmmd_epall
pmmd_epall <- pmmd_ep_temporal %>%
  dplyr::select(person_id, mm_episode, dep_episode, state, record_count, diagnosis_date, diagnosis_date_end, duration, PMMD)

# Update pmmd_ep_mm
pmmd_ep_mm <- pmmd_ep %>% filter(state == "mm")

# Update nopmmd_epall_dx (and exclude same day switches)
nopmmd_epall <- nopmmd_epall_dx %>% filter(!person_id %in% unique(pmmd_group$person_id)) %>% filter(!person_id %in% samedayswitch_anyorder$person_id) 

# Update mania_dep_dates_toepisodes_qc to include these newly-identified dep episodes
# - replace all of the new PMMD rows
mania_dep_dates_toepisodes_qc_broaddef <- mania_dep_dates_toepisodes_qc %>%
  filter(!person_id %in% unique(pmmd_ep_temporal$person_id)) %>%
  mutate(duration = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days"))) %>%
  full_join(pmmd_ep_temporal %>% mutate(duration = as.numeric(duration, units = "days"))) %>%
  dplyr::select(person_id, mm_episode, dep_episode, state, record_count, diagnosis_date, diagnosis_date_end, duration)
saveDf(mania_dep_dates_toepisodes_qc_broaddef, "mania_dep_dates_toepisodes_qc_broaddef")

# identify sensitivity analysis individuals 
switch_timeafter <- pmmd_ep_mm %>% arrange(timetodep)
length(which(switch_timeafter$timetodep == 0))
length(which(switch_timeafter$timetodep != 0))
sensitivity_pmmd_timetodep0 <- pmmd_ep_mm %>% filter(timetodep == 0)

# Compare other constructs (only possible after above) and finalise id_group with additional labels
source("compare_bipolar_constructs.r")
id_group <- id_group0_addsamedayswitch %>% 
  left_join(., rapidcycler_id %>% dplyr::select(person_id, category)) %>%
  rename(rapidcycler = category)

#-------------

saveDf(id_group, "id_group", overwrite = T)
saveDf(pmmd_group, "pmmd_group", overwrite = T)
saveDf(pmmd_ep_temporal, "pmmd_ep_temporal", overwrite = T)
saveDf(pmmd_ep, "pmmd_ep", overwrite = T)
saveDf(pmmd_epall, "pmmd_epall", overwrite = T)
saveDf(pmmd_ep_mm, "pmmd_ep_mm", overwrite = T)

saveDf(sensitivity_pmmd_timetodep0, "sensitivity_pmmd_timetodep0", overwrite = T)

saveDf(nopmmd_group, "nopmmd_group", overwrite = T)
saveDf(nopmmd_epall, "nopmmd_epall", overwrite = T)

saveDf(bipolar_episodes_count, "bipolar_episodes_count", overwrite = T)
#saveDf(mania_dep_dates_toepisodes_qc_broaddef, "mania_dep_dates_toepisodes_qc_broaddef", overwrite = T) # already saved above

#======================================================
# Quantify switches (overall + per individual)
#======================================================

# Find individuals with multiple episodes
episodes_per_id <- table(pmmd_ep_mm$person_id)
pmmd_ep_mm_multi <- pmmd_ep_mm %>% filter(person_id %in% (names(episodes_per_id[which(episodes_per_id > 1)])))

# Explore these individuals
print(paste("Total number of switches = ", sum(episodes_per_id), sep = "")) # number of episodes n=7679

print("Summary of switches per person:")
print(summary(as.vector(episodes_per_id))) # how many switches per person, max=50

print(paste("Unique individuals with a switch = ", length(unique(pmmd_ep_mm$person_id)), sep = "")) # n=4309 individuals

print(paste("Number of individuals with 1 switch = ", length(which(episodes_per_id == 1)), sep = ""))
  
print(paste("Number of individuals with >1 switch = ", length(which(episodes_per_id > 1)), sep = "")) # n=949 with >1 switch episodes

print(paste("Number of individuals with >1 mm episode = ", mania_dep_count_toepisodes_qc %>% filter(group == "PMMD") %>% filter(mm >= 2) %>% nrow(), sep = "")) # n=949 with >1 switch episodes

switchstats <- rbind.data.frame(
  c("Total number of switches", sum(episodes_per_id)),
  c("Unique individuals with a switch", length(unique(pmmd_ep_mm$person_id))),
  c("Number of individuals with 1 switch", length(which(episodes_per_id == 1))),
  c("Number of individuals with >1 switch", length(which(episodes_per_id > 1))),
  c("Number of individuals with >1 manic/mixed episode", mania_dep_count_toepisodes_qc %>% filter(group == "PMMD") %>% filter(mm >= 2) %>% nrow()))
colnames(switchstats) <- c("variable", "value")

png("hist2.png")
hist(as.vector(episodes_per_id[which(episodes_per_id > 1)]))
dev.off()

#======================================================
# Switch time (mania --> depression)
#======================================================

# Time between switch
switch_timeafter <- pmmd_ep_mm %>% arrange(timetodep)
length(which(switch_timeafter$timetodep == 0))
length(which(switch_timeafter$timetodep != 0))

# - find that the strict dx based PMMD definition has fewer
switch_timeafter_dx <- pmmd_ep_mm %>% arrange(timetodep) %>% filter(is.na(addpmmd))
length(which(switch_timeafter_dx$timetodep == 0))
length(which(switch_timeafter_dx$timetodep != 0))

print(paste("Mean time between switch (<1y), incl zero = ", mean(switch_timeafter$timetodep, na.rm = TRUE), sep = "")) # mean = 46 days
print(paste("Median time between switch (<1y), incl zero = ", median(switch_timeafter$timetodep, na.rm = TRUE), sep = "")) # mean = 14 days
switchstats <- rbind(switchstats, c(variable = "Mean (SD) time between switch (<1y), incl zero", value = paste0(mean(switch_timeafter$timetodep), "(", sd(switch_timeafter$timetodep), ")")))
switchstats <- rbind(switchstats, c(variable = "Median time between switch (<1y), incl zero", value = median(switch_timeafter$timetodep)))

print(paste("Mean time between switch (<1y), non-zero only = ", mean(switch_timeafter$timetodep[which(switch_timeafter$timetodep != 0)], na.rm = TRUE), sep = "")) # mean = 64 days
print(paste("Median time between switch (<1y), non-zero only = ", median(switch_timeafter$timetodep[which(switch_timeafter$timetodep != 0)], na.rm = TRUE), sep = "")) # median = 36 days
print(summary(as.numeric(switch_timeafter$timetodep)))
print(summary(as.numeric(switch_timeafter$timetodep[which(switch_timeafter$timetodep != 0)])))
switchstats <- rbind(switchstats, c(variable = "Mean (SD) time between switch (<1y), non-zero only", value = paste0(mean(switch_timeafter$timetodep[which(switch_timeafter$timetodep != 0)]), "(", sd(switch_timeafter$timetodep[which(switch_timeafter$timetodep != 0)]), ")")))
switchstats <- rbind(switchstats, c(variable = "Median time between switch (<1y), non-zero only", value = median(switch_timeafter$timetodep[which(switch_timeafter$timetodep != 0)])))

switch_dep_per_id <- switch_timeafter %>% group_by(person_id) %>% summarise(mean_switch_t = mean(timetodep), median_switch_t = median(timetodep), n = n()) %>% arrange(desc(n))
switchstats <- rbind(switchstats, c(variable = "Mean (SD) time of individual mean between switch (<1y)", value = paste0(mean(switch_dep_per_id$mean_switch_t), "(", sd(switch_dep_per_id$mean_switch_t), ")")))
switchstats <- rbind(switchstats, c(variable = "Median time of individual median between switch (<1y)", value = median(switch_dep_per_id$median_switch_t)))
switchstats <- rbind(switchstats, c(variable = "Mean (SD) time of individual mean between switch (<1y), non-zero only", value = paste0(mean(switch_dep_per_id %>% filter(mean_switch_t > 0) %>% pull(mean_switch_t)), "(", sd(switch_dep_per_id %>% filter(mean_switch_t > 0) %>% pull(mean_switch_t)), ")")))
switchstats <- rbind(switchstats, c(variable = "Median time of individual median between switch (<1y), non-zero only", value = median(switch_dep_per_id %>% filter(mean_switch_t > 0) %>% pull(median_switch_t))))
ggplot(switch_dep_per_id, aes(x = n, y = median_switch_t)) +
  geom_point()
ggsave("n_switch_vs_time_switch.png")

print(paste("Switch rate in 1 month = ", length(which(switch_timeafter$timetodep < 30))/nrow(switch_timeafter), sep = ""))
switch_timeafter_nonzero <- switch_timeafter %>% filter(timetodep != 0)
print(paste("Switch rate in 1 months, non-zero only = ", length(which(switch_timeafter_nonzero$timetodep < 30))/nrow(switch_timeafter_nonzero), sep = ""))
switchstats <- rbind(switchstats, c(variable = "Switch rate in 1 month", value = length(which(switch_timeafter$timetodep < 30))/nrow(switch_timeafter)))
switchstats <- rbind(switchstats, c(variable = "Switch rate in 1 month, non-zero only", value = length(which(switch_timeafter_nonzero$timetodep < 30))/nrow(switch_timeafter_nonzero)))

print(paste("Switch rate in 3 months = ", length(which(switch_timeafter$timetodep < 90))/nrow(switch_timeafter), sep = ""))
switch_timeafter_nonzero <- switch_timeafter %>% filter(timetodep != 0)
print(paste("Switch rate in 3 months, non-zero only = ", length(which(switch_timeafter_nonzero$timetodep < 90))/nrow(switch_timeafter_nonzero), sep = ""))
switchstats <- rbind(switchstats, c(variable = "Switch rate in 3 months", value = length(which(switch_timeafter$timetodep < 90))/nrow(switch_timeafter)))
switchstats <- rbind(switchstats, c(variable = "Switch rate in 3 months, non-zero only", value = length(which(switch_timeafter_nonzero$timetodep < 90))/nrow(switch_timeafter_nonzero)))

ggplot(switch_timeafter %>% mutate(time = as.numeric(timetodep)), aes(x = time)) +
  geom_density()
ggsave("switch_time_density.png")
ggplot(switch_timeafter %>% mutate(time = as.numeric(timetodep)) %>% filter(time != 0), aes(x = time)) +
  geom_density()
ggsave("switch_time_density_nonzerotime.png")

# - check switch time = 0 is not related to diagnosis_date (seems to peak 2020)
switch_zero <- switch_timeafter %>% mutate(samedayswitch = ifelse(timetodep == 0, TRUE, FALSE))
ggplot(switch_zero %>% mutate(time = as.numeric(timetodep)), aes(x = diagnosis_date, colour = samedayswitch)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("switch_time_date.png")

# - check switch time = 0 is not an individual-specific problem
zeroid <- switch_timeafter %>% filter(timetodep == 0) %>% dplyr::select(person_id) %>% unique() # n=2240
nonzeroid <- switch_timeafter %>% filter(timetodep != 0) %>% dplyr::select(person_id) %>% unique() #n=4608
print(length(which(zeroid$person_id %in% nonzeroid$person_id))) # 494 overlap

# Visualise trajectories (although note this is filtered for people with mania > depressive switch within 1y) 
#ggplot(earlier_ep_1y_multi %>% filter(person_id == 100001505534), aes(x = diagnosis_date, y = statenum, colour = person_id)) +
#  geom_line() +
#  theme(legend.position = "none")
#ggsave("trajectory_multi.png", width = 20, height = 7)

#======================================================
# Comparisons to broader population
#======================================================

# 1. To all people who experience mania/mixed
print(paste("All with mania/mixed = ", nrow(id_group_dx %>% filter(group == "MMonly")), sep = ""))
switchstats <- rbind(switchstats, c(variable = "All individuals with mania", value = length(unique(mmonly_count_toepisodes$person_id))))

# 2. To people who experience both poles (but no time limit)
print(paste("All with both poles = ", nrow(bipolar_episodes_count), sep = ""))
switchstats <- rbind(switchstats, c(variable = "All individuals with both poles", value = nrow(bipolar_episodes_count)))

#======================================================
# Time until remission of post-mania depressive episode
#======================================================

# 1. Diagnostic coding for remission
dep_remit_code0 <- getDf("dep_remit_code0")
dep_remit_code_rmmultiicd <- remove_icd_dup(dep_remit_code0)
dep_remit_code <- dep_remit_code0 %>%
  mutate(tmpid = paste0(person_id, diagnosis_date, icd_code)) %>%
  filter(!tmpid %in% dep_remit_code_rmmultiicd$tmpid) %>% 
  dplyr::select(-c("tmpid")) %>%
  rename(remit_date = diagnosis_date) %>%
  dplyr::select(-c("icd_code")) %>%
  mutate(person_id = as.character(person_id))
saveDf(dep_remit_code, "dep_remit_code", overwrite = T)

pmmd_ep_temporal_depremit <- pmmd_ep_temporal %>% 
  filter(state == "dep" & PMMD == "PMMD") %>%
  left_join(dep_remit_code, relationship = "many-to-many") %>%
  mutate(difftime.tmp = remit_date - diagnosis_date) %>%
  filter(difftime.tmp > 0) %>% # ensure remission happens after diagnosis_date_end
  group_by(person_id, diagnosis_date) %>% 
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get 1 per PMMD episode
  ungroup() %>%
  group_by(person_id, remit_date) %>%
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get one per remission code
  ungroup() %>%
  dplyr::select(-difftime.tmp)

# 2. Diagnostic coding for mania (ie. return to)
mm_code <- pmmd_ep_temporal %>% filter(state == "mm") %>%
  dplyr::select(c("person_id", "diagnosis_date")) %>%
  rename(recur_date = diagnosis_date)

pmmd_ep_temporal_mmrecur <- pmmd_ep_temporal %>% 
  filter(state == "dep" & PMMD == "PMMD") %>%
  left_join(mm_code, relationship = "many-to-many") %>%
  mutate(difftime.tmp = recur_date - diagnosis_date) %>%
  filter(difftime.tmp > 0) %>% # ensure recurrence happens after diagnosis_date_end
  group_by(person_id, diagnosis_date) %>% 
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get 1 per PMMD episode
  ungroup() %>%
  group_by(person_id, recur_date) %>%
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get one per recurrence code
  ungroup() %>%
  dplyr::select(-difftime.tmp)

# 3. PHQ-9
phq9_mania <- getDf("phq9_mm_sql") %>%  mutate(person_id = as.character(person_id))
phq9_mania_less10 <- phq9_mania %>% filter(score < 10) %>%
  dplyr::select(c("person_id", "measurement_date", "score")) %>%
  rename(phq9_date = measurement_date)

pmmd_ep_temporal_phq9 <- pmmd_ep_temporal %>% 
  filter(state == "dep" & PMMD == "PMMD") %>%
  left_join(phq9_mania_less10, relationship = "many-to-many") %>%
  mutate(difftime.tmp = phq9_date - diagnosis_date) %>%
  filter(difftime.tmp > 0) %>% # ensure remission happens after diagnosis_date_end
  group_by(person_id, diagnosis_date) %>% 
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get 1 per PMMD episode
  ungroup() %>%
  group_by(person_id, phq9_date) %>%
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get one per phq9 remission
  ungroup() %>%
  dplyr::select(-difftime.tmp)

# 4. CGI-S
madrscgis_mania <- getDf("madrs_mm_sql") %>% mutate(person_id = as.character(person_id))
madrscgis_mania_less20 <- madrscgis_mania %>% filter(score < 20) %>%
  dplyr::select(c("person_id", "measurement_date", "score")) %>%
  rename(madrscgis_date = measurement_date) %>%
  rename(madrscgis_value = score)

pmmd_ep_temporal_madrscgis <- pmmd_ep_temporal %>% 
  filter(state == "dep" & PMMD == "PMMD") %>%
  left_join(madrscgis_mania_less20, relationship = "many-to-many") %>%
  mutate(difftime.tmp = madrscgis_date - diagnosis_date) %>%
  filter(difftime.tmp > 0) %>% # ensure remission happens after diagnosis_date_end
  group_by(person_id, diagnosis_date) %>% 
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get 1 per PMMD episode
  ungroup() %>%
  group_by(person_id, madrscgis_date) %>%
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get one per phq9 remission
  ungroup() %>%
  dplyr::select(-difftime.tmp)

# 5. No visit in >1 year (assume discharge)
#pmmd_visit <- getVisit(as.numeric(pmmd_group$person_id))
#saveDf(pmmd_visit, "pmmd_visit")
pmmdgroup_visit <- getDf("pmmdgroup_visit")
# Find "discharge" proxy
# - filter out visits c("Failed encounter", "No matching concept") as no assessment is possible here
    # - keep "Pharmacy visit" and "Laboratory Visit" as these at least indicate that care is ongoing
# - group by person, order visits, take time between each visit (time_visit_diff)
# - define a careepisode_new if >365 days since last encounter
# - when time_visit_diff >1y, take the visit BEFORE as the discharge date 
  # - will align rows to have the timeto_prevvisit corresponding to the visit date before so filtering is easier
pmmd_visit_timebw <- pmmdgroup_visit %>% ungroup %>%
  mutate(visit_duration = as.numeric(difftime(visit_end_date, visit_start_date, units = "days"))) %>% 
  filter(visit_duration <= 365) %>% # as some visit_durations are >1000d
  dplyr::select(-visit_duration) %>%
  mutate(person_id = as.character(person_id)) %>%
  filter(!visit_type %in% c("Failed encounter", "No matching concept")) %>%
  arrange(visit_start_date) %>%
  group_by(person_id) %>%
  mutate(timeto_nextvisit = as.numeric(difftime(lead(visit_start_date), visit_start_date, units = "days"))) %>%
  ungroup %>%
  mutate(careepisode_new = ifelse(timeto_nextvisit >= 365, "yes", "no"))

pmmd_visit_careepisode <- pmmd_visit_timebw %>% 
  filter(careepisode_new == "yes") %>%
  dplyr::select(-c("careepisode_new"))

pmmd_ep_temporal_dc <- pmmd_ep_temporal %>% 
  filter(state == "dep" & PMMD == "PMMD") %>%
  left_join(pmmd_visit_careepisode, relationship = "many-to-many") %>%
  mutate(difftime.tmp = visit_end_date - diagnosis_date) %>%
  filter(difftime.tmp > 0) %>% # ensure remission happens after diagnosis_date_end
  group_by(person_id, diagnosis_date) %>% 
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get 1 per PMMD episode
  ungroup() %>%
  group_by(person_id, visit_end_date) %>%
  filter(difftime.tmp == min(difftime.tmp)) %>% # ensure you only get one per phq9 remission
  ungroup() %>%
  dplyr::select(-difftime.tmp)

# COMBINE 5 definitions
pmmd_ep_temporal_depend <- pmmd_ep_temporal_depremit %>% distinct %>%
  full_join(., pmmd_ep_temporal_mmrecur %>% distinct) %>%
  full_join(., pmmd_ep_temporal_phq9 %>% distinct) %>%
  full_join(., pmmd_ep_temporal_madrscgis %>% distinct) %>%
  full_join(., pmmd_ep_temporal_dc %>% distinct) %>%
  pivot_longer(cols = c("remit_date", "recur_date", "phq9_date", "madrscgis_date", "visit_end_date"), names_to = "end_def", values_to = "date_dep_end") %>%
  mutate(end_def = gsub("_date", "", end_def)) %>%
  mutate(time_dep_end = as.numeric(difftime(date_dep_end, diagnosis_date, units = "days"))) %>%
  group_by(person_id, diagnosis_date) %>%
  slice_min(time_dep_end, n=1, with_ties = FALSE) %>%
  ungroup()

table(pmmd_ep_temporal_depend$end_def)

ggplot(pmmd_ep_temporal_depend, aes(x = time_dep_end)) +
  stat_ecdf(geom = "step") +  
  theme_bw() +
  scale_colour_manual(values = millerstone[c(4,1,6,3)]) +
  theme(legend.position = "bottom")
ggsave("timedepend_ecdf.png", height = 4, width = 4)

ggplot(pmmd_ep_temporal_depend, aes(x = time_dep_end)) +
  stat_ecdf(geom = "step") +  
  theme_bw() +
  scale_colour_manual(values = millerstone[c(4,1,6,3)]) +
  xlim(0, 1000) +
  theme(legend.position = "bottom")
ggsave("timedepend_ecdf_lim1000d.png", height = 4, width = 4)

ggplot(pmmd_ep_temporal_depend, aes(x = time_dep_end, colour = end_def)) +
  stat_ecdf(geom = "step") +  
  theme_bw() +
  scale_colour_manual(values = millerstone[c(4,1,6,3,2)]) +
  theme(legend.position = "bottom")
ggsave("timedepend_ecdf_remitvsrecurvsphq9vsmadrscgisvsvisit.png", height = 5, width = 5)

ggplot(pmmd_ep_temporal_depend, aes(x = time_dep_end, colour = end_def)) +
  stat_ecdf(geom = "step") +  
  theme_bw() +
  scale_colour_manual(values = millerstone[c(4,1,6,3,2)]) +
  xlim(0, 1000) +
  theme(legend.position = "bottom")
ggsave("timedepend_ecdf_remitvsrecurvsphq9vsmadrscgisvsvisit_lim1000d.png", height = 5, width = 5)

# - summary stats of duration
summary(pmmd_ep_temporal_depend$time_dep_end)

# - % completeness
nrow(pmmd_ep_temporal_depend %>% filter(PMMD == "PMMD" & state == "dep"))
nrow(pmmd_ep_temporal %>% filter(PMMD == "PMMD" & state == "dep"))

saveDf(pmmd_ep_temporal_depend, "pmmd_ep_temporal_depend", overwrite = T)

# Join together for time course
# - get mm episodes for those with complete data
pmmd_ep_temporal_course <- pmmd_ep_temporal_depend %>%
  dplyr::select(person_id, mm_episode, dep_episode) %>%
  mutate(dep_episode = dep_episode - 1) %>%
  inner_join(., pmmd_ep_temporal) %>%
  left_join(., mm_split %>% dplyr::select(person_id, diagnosis_date, diagnosis_date_end, mm_split)) %>%
  rename(mm_diagnosis_date = diagnosis_date) %>%
  rename(mm_diagnosis_date_end.tmp = diagnosis_date_end) %>%
  mutate(dep_episode = dep_episode + 1) %>%
  dplyr::select(person_id, mm_episode, dep_episode, mm_diagnosis_date, mm_diagnosis_date_end.tmp, mm_split) %>%
  inner_join(., pmmd_ep_temporal_depend) %>%
  mutate(mm_diagnosis_date_end.tmp2 = ifelse(as.numeric(difftime(as.Date(mm_diagnosis_date_end.tmp), as.Date(diagnosis_date))) > 0, as.character(diagnosis_date), as.character(mm_diagnosis_date_end.tmp))) %>%
  mutate(mm_diagnosis_date_end = as.Date(mm_diagnosis_date_end.tmp2)) %>%
  dplyr::select(-c("mm_diagnosis_date_end.tmp", "mm_diagnosis_date_end.tmp2")) %>%
  mutate(time_mmstart_mmend = as.numeric(difftime(mm_diagnosis_date_end, mm_diagnosis_date, units = "days"))) %>%
  mutate(time_mmstart_depstart = as.numeric(difftime(diagnosis_date, mm_diagnosis_date, units = "days"))) %>%
  mutate(time_mmend_depstart = as.numeric(difftime(diagnosis_date, mm_diagnosis_date_end, units = "days"))) %>%
  mutate(time_depstart_depend = as.numeric(difftime(date_dep_end, diagnosis_date, units = "days"))) %>%
  #mutate(timetodep = as.numeric(difftime(diagnosis_date, mm_diagnosis_date, units = "days"))) %>%
  mutate(timetoend = time_mmstart_depstart + time_dep_end) %>%
  full_join(pmmd_ep_temporal) # join to complete dataset

saveDf(pmmd_ep_temporal_course, "pmmd_ep_temporal_course", overwrite = T)

#================================================
# Further data munging now that pheno is defined
#================================================

source("data.pulls.r")

#====================
#------------------>>>----------------------
# - some EDA of the characteristics of visits in this population
pmmd_visit %>% pull(person_id) %>% as.numeric(.) %>% table %>% summary
pmmd_visit_type_sum <- pmmd_visit %>% group_by(person_id, visit_type) %>% summarise(n=n()) 
pmmd_visit_type_sum_wide <- pmmd_visit_type_sum %>% pivot_wider(names_from = "visit_type", values_from = "n")
head(data.frame(pmmd_visit_type_sum_wide))
ggplot(pmmd_visit_type_sum, aes(x = n, fill = visit_type)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ visit_type, scale = "free")
ggsave("pmmd_visit_hist.png")

# - time between visit_start_date and visit_end_date
tmp <- pmmd_visit %>% 
  mutate(visit_duration = visit_end_date - visit_start_date)
tmp2 <- tmp %>% filter(visit_duration > 0)
unique(tmp2$visit_type) # find multi-day visits are "Inpatient" "Emergency Room Visit", "Outpatient Visit", "Emergency Room and Inpatient Visit"
ggplot(tmp2, aes(x = visit_duration, fill = visit_type)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ visit_type, scale = "free")
ggsave("pmmd_visit_duration_hist.png")
# - find that visit duration goes up to 1752 days and that these are "Outpatient Visit"
# - seems like a QC issue ...?
# - don't need to apply filter now as will largely just take visit_start_date anyway
# - accurate dates will be more important for inpatient times (as a predictor)
#-------------------<<<---------------------

#======================================================
# Save relevant files
#======================================================

saveDf(switchstats, "switchstats", overwrite = TRUE)

