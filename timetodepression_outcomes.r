#======================================================
#
# Association of time to depression with outcomes
#
#======================================================

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
library(MASS)

source("read_in_routine.r")

#======================================================
# Plot density for overview of the issue 
#======================================================

# Select for:
# - people with mania 
# - next episode is depression
# - episodes from 1999 (as this is when Neuroblu data official starts coming from

bipolar_epall.tmp <- pmmd_epall %>% 
  mutate(PMMD = ifelse(is.na(PMMD), "nonPMMDep", "PMMDep"), 
         group = "PMMD") %>%
  full_join(., nopmmd_epall %>% mutate(group = "noPMMD") %>% mutate(duration = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days")))) %>%
  mutate(PMMD = ifelse(is.na(PMMD), "noPMMD", PMMD))
# not yet saved
bipolar_epall_rm_qc <- bipolar_epall.tmp %>%
  mutate(ep_duration = as.numeric(difftime(diagnosis_date_end, diagnosis_date, units = "days"))) %>%
  filter(diagnosis_date <= as.Date("1999-01-01") | (state == "mm" & ep_duration >= 365)) %>% # exclude manic episodes >1 year in duration
  dplyr::select(person_id, diagnosis_date, diagnosis_date_end, mm_episode, dep_episode, state)
saveDf(bipolar_epall_rm_qc, "bipolar_epall_rm_qc", overwrite = T)

bipolar_epall <- bipolar_epall.tmp %>% anti_join(., bipolar_epall_rm_qc)
saveDf(bipolar_epall, "bipolar_epall", overwrite = T)

mm_split <- getDf("mm_split")

# Find preceding and following depressive episode
bipolar_epall_mm <- bipolar_epall %>% filter(state == "mm") %>% 
  #left_join(., pmmd_ep %>% dplyr::select(person_id, mm_episode, dep_episode, diagnosis_date)) %>%
  left_join(., mm_split %>% dplyr::select(person_id, diagnosis_date, diagnosis_date_end, mm_split)) %>%
  rename(diagnosis_date_mm = diagnosis_date, diagnosis_date_end_mm = diagnosis_date_end) %>%
  dplyr::select(person_id, mm_episode, dep_episode, state, mm_split, diagnosis_date_mm, diagnosis_date_end_mm, group) #227220
bipolar_epall_dep <- bipolar_epall %>% filter(state == "dep") %>% 
  rename(diagnosis_date_dep = diagnosis_date, diagnosis_date_end_dep = diagnosis_date_end) %>%
  dplyr::select(person_id, dep_episode, diagnosis_date_dep, diagnosis_date_end_dep, group) #315573

# ALL people so rate axis includes everyone
bipolar_epall_nextep_join <- bipolar_epall %>% 
  mutate(ep_ind_next = mm_episode + dep_episode) %>%
  rename(state_next = state,
         diagnosis_date_next = diagnosis_date,
         diagnosis_date_end_next = diagnosis_date_end)

# - identify the next episode
bipolar_epall_count_mm_to_nextep.tmp <- bipolar_epall_mm %>%
  mutate(ep_ind_next = mm_episode + dep_episode + 1) %>% # joining on the episode number index
  left_join(., bipolar_epall_nextep_join %>% dplyr::select(person_id, ep_ind_next, state_next, diagnosis_date_next, diagnosis_date_end_next))%>%
  mutate(state_next = ifelse(is.na(state_next), "no_next", state_next))

bipolar_epall_count_mm_to_depafter_inf0 <- bipolar_epall_count_mm_to_nextep.tmp %>%
  mutate(diagnosis_date_next = if_else(state_next == "dep", as.Date(diagnosis_date_next), as.Date(NA)),
         diagnosis_date_end_next = if_else(state_next == "dep", as.Date(diagnosis_date_end_next), as.Date(NA))) %>%
  rename(diagnosis_date_dep = diagnosis_date_next,
         diagnosis_date_end_dep = diagnosis_date_end_next) %>%
  mutate(mm_to_depafter = as.numeric(difftime(diagnosis_date_dep, diagnosis_date_end_mm, units = "days"))/30) %>%
  filter(mm_to_depafter >= 0 | is.na(mm_to_depafter)) 
# remove dep episodes coming before mm. This is an artefact of depression episodes 
# being identified (broad case-finding) 0 days after the manic episode, and then when
# the episodes are reordered they are put after that same manic episode
# leading to negative mm_to-depafter eg. 559880009.
# Given they were identified by case-finding anyway, have decided to exclude them

bipolar_epall_count_mm_to_depafter_windows <- bipolar_epall_count_mm_to_depafter_inf0 %>%
  filter(state_next == "dep") %>%
  mutate(months_mm = as.numeric(difftime(diagnosis_date_end_mm, diagnosis_date_mm, units = "days"))/30,
         months_mmend_depstart = as.numeric(difftime(diagnosis_date_dep, diagnosis_date_end_mm, units = "days"))/30,
         months_dep = as.numeric(difftime(diagnosis_date_end_dep, diagnosis_date_dep, units = "days"))/30,
         id = paste0(person_id, "_", mm_episode)) %>%
  arrange(months_mmend_depstart, months_mm, months_dep) %>%
  mutate(id = factor(id, levels = unique(id))) %>%
  pivot_longer(cols = c("months_mm", "months_mmend_depstart", "months_dep"), names_to = "window", values_to = "months") %>%
  mutate(window = factor(window, levels = c("months_mm", "months_mmend_depstart", "months_dep")))

# 1) compute windows and id
bipolar_epall_count_mm_to_depafter_windows0 <- bipolar_epall_count_mm_to_depafter_inf0 %>%
  mutate(months_mm = as.numeric(difftime(diagnosis_date_mm, diagnosis_date_end_mm, units = "days"))/30,
    months_mmend_depstart = as.numeric(difftime(diagnosis_date_dep, diagnosis_date_end_mm, units = "days"))/30,
    months_dep = as.numeric(difftime(diagnosis_date_end_dep, diagnosis_date_dep, units = "days"))/30,
    id = paste0(person_id, "_", mm_episode)) %>% distinct() %>%
  filter(!is.na(months_mmend_depstart)) %>%
  mutate(months_dep = ifelse(is.na(months_dep), 0, months_dep))

#--------------------
# GENERATE list of spurious associations to exclude ...
mm0_mmenddepstart0_dep0_toexclude <- bipolar_epall_count_mm_to_depafter_windows0 %>% 
  filter(months_mm == 0 & months_mmend_depstart == 0 & months_dep == 0)
saveDf(mm0_mmenddepstart0_dep0_toexclude, "mm0_mmenddepstart0_dep0_toexclude", overwrite = T)

# - For Results para1: Find the number of manic/mixed episodes excluding these spurious associations (calculate % with dep as next mood ep from this)
bipolar_epall_mm_excl0 <- bipolar_epall_mm %>% anti_join(., mm0_mmenddepstart0_dep0_toexclude)
saveDf(bipolar_epall_mm_excl0, "bipolar_epall_mm_excl0", overwrite = T)
# - total bipolar sample size
length(unique(bipolar_epall_mm_excl0$person_id)) # 100791
# - total mm mood episodes
dim(bipolar_epall_mm_excl0) # 208825

bipolar_epall_dep_excl0 <- bipolar_epall_dep %>% anti_join(., mm0_mmenddepstart0_dep0_toexclude)
saveDf(bipolar_epall_dep_excl0, "bipolar_epall_dep_excl0", overwrite = T)
# - total bipolar sample size
length(unique(bipolar_epall_dep_excl0$person_id)) # 108561
# - total mm mood episodes
dim(bipolar_epall_dep_excl0) # 315526

# For time to depression results, exclude same day switch
bipolar_epall_count_mm_to_nextep <- bipolar_epall_count_mm_to_nextep.tmp %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude)  
bipolar_epall_count_mm_to_depafter_inf <- bipolar_epall_count_mm_to_depafter_inf0 %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude)
bipolar_epall_count_mm_to_depafter_windows <- bipolar_epall_count_mm_to_depafter_windows0 %>% 
  anti_join(., mm0_mmenddepstart0_dep0_toexclude) 
saveDf(bipolar_epall_count_mm_to_depafter_inf, "bipolar_epall_count_mm_to_depafter_inf", overwrite = T)

# - For Results para1: Find the number of manic/mixed episodes excluding these spurious associations (calculate % with dep as next mood ep from this)
table(bipolar_epall_count_mm_to_depafter_inf$mm_split)
bipolar_epall_count_mm_to_depafter_inf %>% group_by(person_id, mm_split) %>% summarise(n=n()) %>% group_by(mm_split) %>% summarise(n=n())

# - For Results para1: Find the number of manic/mixed episodes transitioning to dep vs mm vs none
# -- quantify prop of each following episode
bipolar_epall_count_mm_to_nextep %>% group_by(state_next) %>% summarise(n=n(), prop = n()/nrow(bipolar_epall_count_mm_to_nextep))
# -- split by mm_split
bipolar_epall_count_mm_to_nextep %>% 
  group_by(mm_split, state_next) %>% summarise(n=n()) %>%
  group_by(mm_split) %>% mutate(mm_sum = sum(n)) %>%
  mutate(mm_prop = n/mm_sum)

source('remove_spurious_episodes.r')

#--------------------

# 2) create ordering vector by months_mmend_depstart (ascending)
order_ids <- bipolar_epall_count_mm_to_depafter_windows %>%
  arrange(-months_mmend_depstart, months_dep, months_mm) %>%    # small -> large; reverse if you want large->small
  pull(id)

# 3) pivot and set factor levels using that ordering
bipolar_epall_count_mm_to_depafter.plot <- bipolar_epall_count_mm_to_depafter_windows %>%
  filter(months_mmend_depstart >= 0) %>%
  pivot_longer(cols = c("months_mm", "months_mmend_depstart", "months_dep"),
    names_to = "window",
    values_to = "months") %>%
  mutate(window = factor(window, levels = rev(c("months_mm", "months_mmend_depstart", "months_dep"))),
    id = factor(id, levels = order_ids))   # IMPORTANT: set levels to the precomputed order

# 4) plot (make sure you reference the same dataframe plot_df)
ggplot(bipolar_epall_count_mm_to_depafter.plot, aes(x = months, y = id, fill = window)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  scale_fill_manual(values = millerstone[c(2,4,8)]) +
  coord_cartesian(xlim = c(-6, 6))
ggsave("mm_dep_windows.png", height = 20, width = 15)

# - quantify proportion of all episodes with transition
deptransition_cumulative_stats <- bipolar_epall_count_mm_to_depafter_inf %>%
  mutate(dep_transition_threshold = case_when(mm_to_depafter <= 1 ~ "upto01m",
                                              mm_to_depafter <= 3 ~ "upto03m",
                                              mm_to_depafter <= 6 ~ "upto06m",
                                              mm_to_depafter <= 12 ~ "upto12m",
                                              .default = "z12m")) %>%
  group_by(mm_split, dep_transition_threshold) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = "mm_split", values_from = "n") %>%
  mutate(mania_cumsum = cumsum(mania), mixed_cumsum = cumsum(mixed)) %>%
  mutate(mania_cumprop = mania_cumsum/sum(mania),
         mixed_cumprop = mixed_cumsum/sum(mixed))
saveDf(deptransition_cumulative_stats, "deptransition_cumulative_stats", overwrite = T)

# - quantify of all people
bipolar_epall_count_mm_to_depafter_inf %>% 
  group_by(person_id, mm_split) %>% summarise(n=n()) %>% group_by(mm_split) %>% summarise(n=n())

totalpercent_mds <- bipolar_epall_count_mm_to_depafter_inf %>%
  mutate(ever_deptransition = ifelse(!is.na(mm_to_depafter), "dep", "nodep")) %>%
  group_by(mm_split, ever_deptransition) %>%
  summarise(md_transition_ep = n()) %>%
  pivot_wider(names_from = "ever_deptransition", values_from = "md_transition_ep") %>%
  mutate(total_ep = dep + nodep,
         prop_dep = dep/total_ep)

bipolar_epall_count_mm_to_depafter_inf.plot <- bipolar_epall_count_mm_to_depafter_inf %>% 
  left_join(totalpercent_mds %>% dplyr::select(mm_split, prop_dep)) %>%
  rename(monthstodep = mm_to_depafter) %>% 
  mutate(monthstodep = ifelse(is.na(monthstodep), 100, monthstodep)) %>%
ggplot(., aes(x = monthstodep, colour = mm_split)) +
  stat_ecdf() +
  theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  xlab("Time to depression (months)") +
  ylab("Cumulative proportion") +
  geom_hline(aes(yintercept = prop_dep), linetype = "dashed") +
  coord_cartesian(xlim = c(0,24)) + scale_x_continuous(breaks = seq(0, 24, by = 3), minor_breaks = NULL) +
  facet_wrap(~mm_split, scales = "free_x")
ggsave("bipolar_epall_count_mm_to_depafter_inf_maniamixed.png", bipolar_epall_count_mm_to_depafter_inf.plot, height = 2.2, width = 4)

# Calculate gradient at per-month intervals on the cumulative density plot
cumprop <- bipolar_epall_count_mm_to_depafter_inf %>% 
  left_join(., totalpercent_mds %>% dplyr::select(mm_split, total_ep)) %>%
  filter(mm_to_depafter >= 0) %>%
  group_by(mm_split) %>%
  arrange(mm_to_depafter) %>%
  mutate(id = row_number()) %>%
  mutate(prop = id/total_ep)

ttd_cutoff <- 0:24

cumprop_mania <- tibble(ttd_cutoff = ttd_cutoff) %>%
  mutate(
    data = map(ttd_cutoff, ~ cumprop %>% filter(mm_to_depafter >= .x & mm_split == "mania")),
    model = map(data, ~ lm(prop ~ mm_to_depafter, data = .x)),
    tidied = map(model, tidy))
cumprop_gradient_mania <- cumprop_mania %>% dplyr::select(ttd_cutoff, tidied) %>% unnest(tidied) %>% filter(term == "mm_to_depafter")
cumprop_mixed <- tibble(ttd_cutoff = ttd_cutoff) %>%
  mutate(
    data = map(ttd_cutoff, ~ cumprop %>% filter(mm_to_depafter >= .x & mm_split == "mixed")),
    model = map(data, ~ lm(prop ~ mm_to_depafter, data = .x)),
    tidied = map(model, tidy))
cumprop_gradient_mixed <- cumprop_mixed %>% dplyr::select(ttd_cutoff, tidied) %>% unnest(tidied) %>% filter(term == "mm_to_depafter")

cumprop_gradient_mm <- rbind(data.frame(cumprop_gradient_mania %>% mutate(mm_split = "mania")),
                            data.frame(cumprop_gradient_mania %>% mutate(mm_split = "mixed")))
ggplot(cumprop_gradient_mm, aes(x = ttd_cutoff, y = estimate, colour = mm_split)) +
  geom_point() +
  theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  facet_wrap(~mm_split) +
  coord_cartesian(xlim = c(0,24)) + scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  ylab("Gradient of Time to depression ~ Cumulative prop") + xlab("Time to depression cut-off (months)")
ggsave("gradient_cumprop.png", height = 4, width = 7)

# - before mm episode
bipolar_epall_count_mm_to_depbefore <- bipolar_epall_mm_excl0 %>%
  left_join(., bipolar_epall_dep) %>% # no need to do dep episode index change as for the after (as already linked)
  mutate(mm_to_depbefore = as.numeric(difftime(diagnosis_date_dep, diagnosis_date_mm, units = "days"))/30) %>%
  filter(dep_episode >= 1) %>% # exclude if no previous dep episodes
  group_by(person_id, dep_episode) %>% # for each dep episode, take the least-negative time since mania
  filter(mm_to_depbefore == max(mm_to_depbefore))

bipolar_epall_count_mm_to_depafter.tmp <- bipolar_epall_count_mm_to_depafter_inf %>% 
  rename(monthstodep = mm_to_depafter) %>% 
  mutate(type = "Depression after")
bipolar_epall_count_mm_to_depbefore.tmp <- bipolar_epall_count_mm_to_depbefore %>% 
  rename(monthstodep = mm_to_depbefore) %>% 
  mutate(type = "Depression before")

bipolar_epall_count_mm_to_deprecent <- full_join(bipolar_epall_count_mm_to_depbefore.tmp, bipolar_epall_count_mm_to_depafter.tmp) %>%
  mutate(monthsmm = as.numeric(difftime(diagnosis_date_end_mm, diagnosis_date_mm, units = "days"))/30) %>%
  mutate(type = factor(type, levels = c("Depression before", "Depression after")))
saveDf(bipolar_epall_count_mm_to_deprecent, "bipolar_epall_count_mm_to_deprecent", overwrite = T)

#----------------------------------
# Density plot
#----------------------------------

ggplot(bipolar_epall_count_mm_to_deprecent, aes(x = monthstodep, colour = mm_split)) +
  geom_density() +
  theme_bw() +
  scale_colour_manual(values = millerstone[c(1,4,2)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 3)) + 
  facet_grid(.~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_deprecent.png", height = 2, width = 6)

# This plot removes monthstodep == 0
# - density scale is of people who have depression within the 24 months
bipolar_epall_count_mm_to_deprecent_maniamixed <- ggplot(bipolar_epall_count_mm_to_deprecent %>% filter(monthstodep != 0 & monthstodep >= -24 & monthstodep <= 24), aes(x = monthstodep, colour = type)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(7,8)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 6)) + 
  facet_grid(mm_split~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_deprecent_maniamixed.png", bipolar_epall_count_mm_to_deprecent_maniamixed, height = 2.2, width = 4)

# - keep monthstodep == 0
bipolar_epall_count_mm_to_deprecent_maniamixed_keepmonthstodep0 <- ggplot(bipolar_epall_count_mm_to_deprecent %>% filter(monthstodep >= -24 & monthstodep <= 24), aes(x = monthstodep, colour = type)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(9,8)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 6)) + 
  facet_grid(mm_split~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_deprecent_maniamixed_keepmonthstodep0.png", bipolar_epall_count_mm_to_deprecent_maniamixed_keepmonthstodep0, height = 2.2, width = 4)

bipolar_epall_count_mm_to_depafter_maniamixed_keepmonthstodep0 <- bipolar_epall_count_mm_to_deprecent %>% 
  filter(type == "Depression after") %>%
  filter(monthstodep >= 0 & monthstodep <= 24) %>%
  ggplot(., aes(x = monthstodep, colour = type)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[c(8)]) +
  xlab("Time to depression (months)") +
  scale_x_continuous(breaks = seq(-24, 24, by = 6)) + 
  facet_grid(mm_split~type, scales = "free_x", space = "free")
ggsave("bipolar_epall_count_mm_to_depafter_maniamixed_keepmonthstodep0.png", bipolar_epall_count_mm_to_depafter_maniamixed_keepmonthstodep0, height = 2.2, width = 4)

#----------------------------------
# Baseline depression rate
#----------------------------------

# - get observation period
observation_period <- getDf("observation_period")
mm_split <- getDf("mm_split")
bipolar_epall_mm_excl0 <- getDf("bipolar_epall_mm_excl0")
bipolar_epall_count_mm_to_deprecent <- getDf("bipolar_epall_count_mm_to_deprecent")

demog_cov <- getDf("demog_cov")
observation_period_bipolar_time_perind <- demog_cov %>% 
  dplyr::select(person_id, observation_years) %>%
  rename(observation_period_years_sum = observation_years) %>%
  filter(person_id %in% bipolar_epall_count_mm_to_deprecent$person_id) %>%
  left_join(., mm_split %>% dplyr::select("person_id", "mm_split")) %>% distinct # end up with more rows than total number - fine to double count in denominator as separate calculations

#observation_period_bipolar_time_perind <- observation_period %>% 
#  dplyr::select(person_id, observation_period_start_date, observation_period_end_date) %>% 
#  inner_join(id_group %>% filter(group %in% c("PMMD", "noPMMD")), .) %>%
#  mutate(observation_period_years = as.numeric(difftime(observation_period_end_date, observation_period_start_date, units = "days"))/365) %>%
#  group_by(person_id) %>%
#  summarise(observation_period_years_sum = sum(observation_period_years)) %>%
#  left_join(., mm_split %>% dplyr::select("person_id", "mm_split")) %>% distinct

observation_period_bipolar <- observation_period %>% filter(person_id %in% bipolar_epall_mm_excl0$person_id)
mmep_denominator <- bipolar_epall_count_mm_to_depafter_inf %>% # have checked no overlap with bipolar_epall_mm_excl0
  #filter(mm_to_depafter >= 0) %>%
  left_join(., observation_period_bipolar) %>%
  filter(observation_period_start_date <= diagnosis_date_dep | is.na(diagnosis_date_dep)) %>%
  mutate(obs_monthssincemm = as.numeric(difftime(observation_period_end_date, diagnosis_date_end_mm, units = "days"))/30) %>%
  mutate(obs_monthssincemm = ifelse(obs_monthssincemm <= 0, 0.1, obs_monthssincemm)) %>% # include 0 of negative obs time in bin=1
  mutate(obs_monthssincemm_bin = ceiling(obs_monthssincemm)+1) %>% # round up and +1 so that drop-outs are excluded in next bracket (this is so that same-day-switches+dropout in the numerator make sense)
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

# baseline
all_dep_ep_count <- bipolar_epall %>% filter(state == "dep") %>% nrow
all_mm_ep_count <- bipolar_epall_mm_excl0 %>% nrow
all_obstime_years <- sum(observation_period_bipolar_time_perind$observation_period_years_sum)
all_dep_ep_count/all_obstime_years
all_mm_ep_count/all_obstime_years

ttd_month <- bipolar_epall_count_mm_to_deprecent %>% 
  filter(type == "Depression after") %>% filter(monthstodep >= 0) %>%
  mutate(month = ceiling(monthstodep)) %>% 
  mutate(month = ifelse(month == 0, 1, month)) %>%
  group_by(month, mm_split) %>% summarise(n_dep_transition = n()) %>%
  ungroup %>%
  inner_join(., ttd_incidence_mm_denominator) %>% mutate(rate = n_dep_transition/postmm_personyears) %>% 
  mutate(all_dep_ep_count = all_dep_ep_count, all_obstime_years = all_obstime_years) %>%
  mutate(base_rate = (all_dep_ep_count/all_obstime_years)) %>%
  mutate(vs_base = rate / base_rate) %>%
  mutate(dep_pc_permonth = n_dep_transition / all_dep_ep_count)
saveDf(ttd_month, "ttd_month", overwrite = T)

# - calculate % of episodes switching within 6m
ttd_month %>% filter(mm_split == "mania") %>% filter(month <= 6) %>% pull(n_dep_transition) %>% sum
ttd_month %>% filter(mm_split == "mixed") %>% filter(month <= 6) %>% pull(n_dep_transition) %>% sum

# - plot
ttd_month.plot <- ggplot(ttd_month, aes(x = month, y = vs_base, colour = mm_split)) +
  geom_point() +
  theme_bw() + scale_colour_manual(values = millerstone[c(7,4)]) +
  theme(legend.position = "none") +
  xlab("Month after mania/mixed episode") +
  ylab("Monthly / Overall depression rate") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_x_continuous(limits = c(1, 12), breaks = seq(1, 12, by = 1), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(0, 12, by = 2), minor_breaks = NULL) +
  facet_wrap(~ mm_split, ncol = 2)
ggsave("ttd_month.png", ttd_month.plot, height = 3.3, width = 6)

# - proportion of all depressive episodes
ttd_month %>% filter(month == 1) %>% pull(n_dep_transition) %>% sum / nrow(bipolar_epall_dep)

# OPTION: to look at depression just before mania/mixed
ttd_month_before <- bipolar_epall_count_mm_to_deprecent %>% 
  filter(type == "Depression before") %>% filter(monthstodep <= 0) %>%
  mutate(month = ceiling(monthstodep)) %>% 
  mutate(month = ifelse(month == 0, -1, month)) %>%
  group_by(month, mm_split) %>% summarise(n_dep = n()) %>%
  inner_join(., ttd1m_count_total) %>% mutate(rate = n_dep/n_mm) %>% 
  mutate(vs_base = rate / (all_dep_ep_count/all_obstime_years))

ttd_month_before.plot <- ggplot(ttd_month_before, aes(x = month, y = vs_base, colour = mm_split)) +
  geom_point() +
  theme_bw() + scale_colour_manual(values = millerstone[c(7,4)]) +
  theme(legend.position = "none") +
  xlab("Month before mania/mixed episode") +
  ylab("Month dep rate / Base dep rate") +
  geom_hline(yintercept = 1) +
  scale_x_continuous(limits = c(-1, -36), breaks = seq(-1, -36, by = -3)) +
  facet_wrap(~ mm_split, ncol = 2)
ggsave("ttd_month_before.png", ttd_month_before.plot, height = 3.3, width = 6)

#======================================================
# Post-episode PHQ-9 maximum
#======================================================

mania_dep_dates_toepisodes_qc_broaddef <- getDf("mania_dep_dates_toepisodes_qc_broaddef")
bipolar_epall_count_mm_to_depafter_inf <- getDf("bipolar_epall_count_mm_to_depafter_inf")
bipolar_visit <- getDf("bipolar_visit")
phq9_mania <- getDf("phq9_mm_sql") %>% mutate(person_id = as.character(person_id)) %>% rename(phq9_total = score)
phq9_mania %>% filter(person_id %in% bipolar_epall_count_mm_to_depafter_inf$person_id) %>% pull(person_id) %>% unique %>% length
madrscgis_mania <- getDf("madrs_mm_sql") %>% mutate(person_id = as.character(person_id))

# TRIAL
pmmd_epall_hospital_overlap <- getDf("pmmd_epall_hospital_overlap")
nopmmd_epall_hospital_overlap <- getDf("nopmmd_epall_hospital_overlap")
bipolar_epall_hospital_overlap <- full_join(pmmd_epall_hospital_overlap, nopmmd_epall_hospital_overlap)
bipolar_epall_hospital_overlap_mm <- bipolar_epall_hospital_overlap %>% 
  filter(state == "mm") %>% dplyr::select(person_id, diagnosis_date, diagnosis_date_end) %>%
  mutate(hospital01 = "yes")

# - take maximum PHQ-9 within 365 days from mm episode END
mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend.tmp <- mania_dep_dates_toepisodes_qc_broaddef %>%
  inner_join(bipolar_epall_count_mm_to_depafter_inf %>% filter(!is.na(mm_to_depafter)) %>% dplyr::select(person_id, mm_episode, dep_episode, state)) %>%
  inner_join(., phq9_mania_total, relationship = "many-to-many") %>%
  left_join(., bipolar_epall_hospital_overlap_mm) %>% mutate(hospital01 = ifelse(is.na(hospital01), "no", "yes")) %>%
  mutate(months_postmm = as.numeric(difftime(measurement_date, diagnosis_date_end, units = "days"))/30) %>%
  filter(months_postmm >= 0 & months_postmm <= 12) %>%
  distinct()

# - have to have a PHQ-9 measurement in 6m
phq9_less6m <- mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend.tmp %>%
  group_by(person_id, diagnosis_date) %>%
  summarise(keep = ifelse(min(months_postmm) <= 6, "phq9_6m_keep", "rm")) %>%
  filter(keep == "phq9_6m_keep") %>% dplyr::select(-keep)

mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend <- mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend.tmp %>%
  inner_join(., phq9_less6m) %>% 
  group_by(person_id, measurement_date, phq9_total) %>% # filter for the closest matching PHQ-9 score to the episode
  filter(months_postmm == min(months_postmm)) %>% ungroup() %>%
  group_by(person_id, mm_episode, dep_episode) %>% # take the soonest score that is the max
  mutate(n_obs = n()) %>%
  arrange(desc(phq9_total), months_postmm, .by_group = TRUE) %>%
  slice_head(n=1) %>% ungroup() %>%
  distinct() %>%
  left_join(., demog_cov) %>%
  mutate(cov_age_ep = as.numeric(difftime(diagnosis_date, as.Date(birth_year), units = "days"))/365) %>%
  mutate(cov_age_ehr = as.numeric(difftime(as.Date(paste("2025-01-01")), as.Date(birth_year), units = "days"))/365) %>%
  dplyr::select(-birth_year) %>%
  left_join(., mm_split_mmdep %>% filter(state == "mm"))
saveDf(mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend, "mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend", overwrite = T)

phq9_max_mmend.long <- mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend %>% 
  pivot_longer(cols = "phq9_total", names_to = "scale", values_to = "score")
saveDf(phq9_max_mmend.long, "phq9_max_mmend.long", overwrite = T)

# Plot phq9_total vs time_postep (n=11860)
phq9_max_mmend.long %>% filter(n_obs >= 2) %>% filter(score >= 10) %>% dim #11669
phq9_max_mmend.long %>% filter(n_obs >= 2) %>% filter(score >= 10) %>% group_by(person_id, diagnosis_date) %>% summarise(n=n()) %>% dim
phq9_max_mmend.long %>% filter(n_obs >= 2) %>% filter(score >= 10) %>% group_by(person_id) %>% summarise(n=n()) %>% dim #9165
phq9_max_mmend.long %>% filter(n_obs >= 2) %>% filter(score >= 10) %>%
  ggplot(., aes(x = months_postmm, colour = mm_split)) +
  geom_density() +
  theme_bw() + theme(legend.position = "none") + scale_colour_manual(values = millerstone[c(7,4)]) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 24, by = 3)) +
  xlab("Months post-MM at max PHQ-9\u226510 (obs\u22652)") +
  facet_wrap(~ mm_split)
ggsave("phq9_max_mmend.png", height = 2.2, width = 4.4)

#----------------------------------
# Permutation
#----------------------------------

# - take maximum PHQ-9 within 365 days from mm episode END but RANDOMISE PHQ-9 values for permutation
base_df0 <- mania_dep_dates_toepisodes_qc_broaddef %>%
  inner_join(bipolar_epall_count_mm_to_depafter_inf %>% filter(!is.na(mm_to_depafter)) %>% dplyr::select(person_id, mm_episode, dep_episode, state)) %>%
  filter(state == "mm") %>%
  inner_join(phq9_mania_total, relationship = "many-to-many") %>% 
  mutate(months_postmm = as.numeric(difftime(measurement_date, diagnosis_date_end, units = "days"))/30) %>%
  filter(months_postmm >= 0 & months_postmm <= 12) %>%
  distinct() %>%
  group_by(person_id, measurement_date, phq9_total) %>%
  filter(months_postmm == min(months_postmm)) %>%
  ungroup()
obs2_phq910_filter <- base_df0 %>%
  group_by(person_id, diagnosis_date) %>%
  summarise(n=n(), phq9_max = max(phq9_total))
obs2_phq910_filter_keep <- obs2_phq910_filter %>% filter(n>=2 & phq9_max>=10) %>% dplyr::select(person_id, diagnosis_date)
base_df <- base_df0 %>% inner_join(., obs2_phq910_filter_keep)

# NOTE: this arrange makes a big difference: arrange(desc(PHQ9_Total_random), months_postmm)) vs arrange(desc(PHQ9_Total_random)))
# - because of there being duplication of PHQ-9 values in 610/3385 cases, so you get basically a "double go" at drawing a sooner date
# - base_df %>% group_by(person_id, mm_episode, dep_episode) %>% summarise(n = n(), n_unique = n_distinct(phq9_total)) %>% filter(n > n_unique) %>% nrow()
# - the point of this permutation is basically to get an empirical distribution as if sampling random measurement dates
# - so when there is a tie, it is appropriate to take one of the dates at random

permute_once <- function(df) {
  df %>%
    group_by(person_id, mm_episode, dep_episode) %>%
    mutate(n_obs = n(), PHQ9_Total_random = phq9_total[sample(n())], tiebreak = runif(n())) %>%
    arrange(desc(PHQ9_Total_random), tiebreak, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
}

# n=1000 permutations
nsim <- 1000
# reproducible seeds so the whole simulation is reproducible
set.seed(2025)
seeds <- sample.int(.Machine$integer.max, nsim)

permuted_all <- purrr::map_dfr(seq_len(nsim), function(i) {
  set.seed(seeds[i])
  permute_once(base_df) %>%
    mutate(iter = i)
})
saveDf(permuted_all, "permuted_all", overwrite = T)

# - get the distribution of median time to max depression 
permuted_phq9_date <- permuted_all %>% 
  left_join(., mm_split_mmdep %>% filter(state == "mm")) %>%
  filter(n_obs >= 2 & PHQ9_Total_random >= 10) %>%
  group_by(iter, mm_split) %>%
  summarise(median = median(months_postmm)) %>%
  pivot_wider(names_from = "mm_split", values_from = "median")
# - get the median for the real data
true_phq9_date_median <- mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend %>%
  filter(n_obs >= 2 & phq9_total >= 10) %>%
  group_by(mm_split) %>%
  summarise(median = median(months_postmm), q25 = quantile(months_postmm, 0.25), q75 = quantile(months_postmm, 0.75)) 
true_phq9_date_median.long <- true_phq9_date_median %>%
  pivot_longer(cols = c("median", "q25", "q75"), names_to = "quantile", values_to = "value")
# - calculate the rank
permuted_phq9_date %>% filter(mania < true_phq9_date_median$median[1]) %>% dim+1
permuted_phq9_date %>% filter(mixed < true_phq9_date_median$median[2]) %>% dim+1

permuted_all_gg <- permuted_all %>% 
  left_join(., mm_split_mmdep %>% filter(state == "mm")) %>%
  filter(n_obs >= 2 & PHQ9_Total_random >= 10) %>%
  dplyr::select(months_postmm, n_obs, mm_split, iter) %>% 
  rbind(., data.frame(mania_dep_dates_toepisodes_qc_broaddef_phq9_mmend %>% 
                      filter(n_obs >= 2 & phq9_total >= 10) %>%
                      dplyr::select(months_postmm, n_obs, mm_split) %>% mutate(iter = 0))) 
saveDf(permuted_all_gg, "permuted_all_gg", overwrite = T)

phq9_max_mmend_permutation.plot <- ggplot() +
  geom_density(
    data = permuted_all_gg %>% filter(iter != 0),
    aes(x = months_postmm, group = iter),
    color = "grey20", alpha = 0.08, size = 0.3) +
  geom_density(
    data = permuted_all_gg %>% filter(iter == 0),
    aes(x = months_postmm), color = "red", fill = NA, size = 1.2) +
  geom_vline(
    data = true_phq9_date_median.long %>% filter(quantile == "median"),
    aes(xintercept = value),
    #linetype = "dashed",
    color = millerstone[2],
    linewidth = 0.5) +
  geom_vline(
    data = true_phq9_date_median.long %>% filter(quantile != "median"),
    aes(xintercept = value),
    linetype = "dashed",
    color = millerstone[2],
    linewidth = 0.5) +
  theme_bw() + theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 24, by = 3)) +
  xlab("Months post-MM at max PHQ-9\u226510 (obs\u22652)") + 
  facet_wrap(~ mm_split, ncol = 2)
ggsave("phq9_max_mmend_permutation.png", phq9_max_mmend_permutation.plot, height = 2.2, width = 4.4)

#======================================================
# TTD on severity of depression
#======================================================

#==================================
# VISITS
#==================================

#----------------------------------
# Check visit timing relative to MDS
#----------------------------------

bipolar_visit <- getDf("bipolar_visit")
bipolar_visit$person_id <- as.character(bipolar_visit$person_id)

bipolar_epall_count_mm_to_depafter_inf_visitdistribution <- bipolar_epall_count_mm_to_depafter_inf %>% 
  left_join(bipolar_visit, relationship = "many-to-many") %>%
  filter(visit_start_date >= diagnosis_date_end_mm) %>% # visit has to come after or with end of mm episode
  mutate(months_from_mm_end = as.numeric(difftime(visit_start_date, diagnosis_date_end_mm, units = "days"))/30) %>%
  mutate(months_from_dep_start = as.numeric(difftime(visit_start_date, diagnosis_date_dep, units = "days"))/30)

# Calculation % of manic episodes with follow-up within 1 month
bipolar_epall_count_mm_to_depafter_inf_visitdistribution_visit1 <- bipolar_epall_count_mm_to_depafter_inf %>% 
  left_join(bipolar_visit, relationship = "many-to-many") %>%
  filter((visit_start_date >= diagnosis_date_end_mm & diagnosis_date_mm != diagnosis_date_end_mm) | (visit_start_date > diagnosis_date_end_mm)) %>% # visit has to come after mm episode if start date == end date, or can coincide with end date if start date != end date
  mutate(months_from_mm_end = as.numeric(difftime(visit_start_date, diagnosis_date_end_mm, units = "days"))/30) %>%
  mutate(months_from_dep_start = as.numeric(difftime(visit_start_date, diagnosis_date_dep, units = "days"))/30) %>%
  group_by(person_id, mm_episode, state, mm_split, diagnosis_date_mm, diagnosis_date_end_mm) %>%
  summarise(visit1_months = min(months_from_mm_end))
(bipolar_epall_count_mm_to_depafter_inf_visitdistribution_visit1 %>% filter(visit1_months <= 1) %>% nrow)/nrow(bipolar_epall_count_mm_to_depafter_inf_visitdistribution_visit1)

ggplot(bipolar_epall_count_mm_to_depafter_inf_visitdistribution %>% filter(months_from_mm_end >= 0 & months_from_mm_end <= 24), aes(x = months_from_mm_end, colour = visit_type)) +
  geom_density() +
  theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = millerstone[1:5]) +
  xlab("Months from depression start") +
  ylab("Visits following mania/mixed episode end") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 3)) +
  facet_grid(visit_type ~ mm_split, scales = "free_y")
ggsave("visits_monthssincemmend_density.png")

#----------------------------------
# TTD ~ all-cause hospital visits in 6 months following depression diagnosis
#----------------------------------

bipolar_visit <- getDf("bipolar_visit")

bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost <- bipolar_epall_count_mm_to_deprecent %>% 
  ungroup %>%
  filter(type == "Depression after") %>% filter(abs(monthstodep) <= 6) %>%
  left_join(distinct(bipolar_visit), relationship = "many-to-many") %>%
  filter(person_id %in% bipolar_visit$person_id) %>%
  filter(visit_start_date >= diagnosis_date_end_mm) %>% # visit has to come after mm episode
  mutate(months_from_dep_start = as.numeric(difftime(visit_start_date, diagnosis_date_dep, units = "days"))/30) %>%
  #filter(months_from_dep_start >= 0) %>% # want to keep this so you don't have total dropout of some individuals (still want to keep as 0)
  mutate(visit6m = ifelse(months_from_dep_start >= 0 & months_from_dep_start <= 6, 1, 0)) %>%
  mutate(visit6m = ifelse(is.na(months_from_dep_start), 0, visit6m)) %>% # keep people without further visits in the denominator
  mutate(dep_180d_date = as.Date(diagnosis_date_dep) + 180) %>%
  mutate(endpoint_180d = if_else(difftime(visit_end_date+1, dep_180d_date, units = "days") > 0, as.Date(dep_180d_date), as.Date(visit_end_date)+1)) %>% # add +1 as hospitalisation days starts forom 1
  mutate(visit_duration = as.numeric(difftime(endpoint_180d, visit_start_date, units = "days"))) %>% # quantify the time spent in hospital
  mutate(visit_duration = ifelse(is.na(visit_duration), 0, visit_duration)) %>%
  mutate(visit_duration = ifelse(visit6m == 0, 0, visit_duration))

# - separately pivot_wider for VisitsN and VisitsDays or get duplicate IDs
bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_count_VisitsN <- bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost %>%
  group_by(person_id, mm_episode, dep_episode, state, mm_split, diagnosis_date_mm, diagnosis_date_end_mm, group, diagnosis_date_dep, diagnosis_date_end_dep, monthstodep, monthsmm, visit_type) %>%
  summarise(n_visits_6m = sum(visit6m)) %>% #, visits_days = sum(visit_duration))
  ungroup %>%
  pivot_wider(names_from = "visit_type", values_from = "n_visits_6m") %>%
  mutate(OutpatientVisitsN = replace_na(Outpatient, 0), Inpatient = replace_na(Inpatient, 0), Emergency = replace_na(Emergency, 0)) %>%
  mutate(TotalVisitsN = Outpatient + Inpatient + Emergency) %>%
  mutate(HospitalVisitsN = Inpatient + Emergency) %>%
  dplyr::select(-c("Emergency", "Outpatient", "Inpatient")) 

bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_count_VisitsDays <- bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost %>%
  group_by(person_id, mm_episode, dep_episode, state, mm_split, diagnosis_date_mm, diagnosis_date_end_mm, group, diagnosis_date_dep, diagnosis_date_end_dep, monthstodep, monthsmm, visit_type) %>%
  summarise(visits_days = sum(visit_duration)) %>% ungroup %>%
  pivot_wider(names_from = "visit_type", values_from = "visits_days") %>%
  mutate(OutpatientVisitsDays = replace_na(Outpatient, 0), Inpatient = replace_na(Inpatient, 0), Emergency = replace_na(Emergency, 0)) %>%
  mutate(TotalVisitsDays = Outpatient + Inpatient + Emergency) %>%
  mutate(HospitalVisitsDays = Inpatient + Emergency) %>%
  dplyr::select(-c("Emergency", "Outpatient", "Inpatient")) %>%
  mutate(HospitalVisitsDays = ifelse(HospitalVisitsDays > 180, 180, HospitalVisitsDays))

bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_count <- bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_count_VisitsN %>%
  inner_join(., bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_count_VisitsDays) %>%
  filter(monthstodep >= 0)

bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_hosp_count <- bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_count %>%
  dplyr::select(-c("OutpatientVisitsN", "TotalVisitsN", "OutpatientVisitsDays", "TotalVisitsDays"))
bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_hosp_count.long <- bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_hosp_count %>%
  pivot_longer(cols = c("HospitalVisitsN", "HospitalVisitsDays"), names_to = "visit_type", values_to = "units")

bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_hosp_count.long %>%
  filter(visit_type == "HospitalVisitsDays") %>%
  ggplot(., aes(x = monthstodep, y = units)) +
  geom_jitter(width = 0.15, height = 0.1,
              size = 0.8, alpha = 0.1, color = "grey20") +
  theme_bw() + theme(legend.position = "none") +
  scale_fill_manual(values = millerstone[c(9,2,5)]) +
  xlab("Time to depression (months)") +
  ylab("Number of hospital days within 6 months of depression onset") +
  facet_grid(visit_type ~ mm_split, scales = "free")
ggsave("hospitaldays_monthstodep_scatter.png")

# Prep for regression
hospital_visit_ttd.glm_prep0 <- bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost_hosp_count %>%
  mutate(hospital01 = ifelse(HospitalVisitsN == 0, 0, 1)) %>%
  left_join(., demog_cov)
cov_state_n.tmp <- hospital_visit_ttd.glm_prep0 %>%
  group_by(person_id, cov_state) %>% summarise(n=n()) %>% group_by(cov_state) %>% summarise(n=n()) %>%
  mutate(cov_state2 = ifelse(n <= 20, "grouped", cov_state))
hospital_visit_ttd.glm_prep <- hospital_visit_ttd.glm_prep0 %>%
  left_join(., cov_state_n.tmp %>% dplyr::select(cov_state, cov_state2)) %>%
  rename(cov_state0 = cov_state, cov_state = cov_state2) %>%
  filter(cov_gender != "No matching concept") %>%
  mutate(cov_edu2 = case_when(cov_edu %in% c("8th grade or less", "Some high school, but did not graduate", "High school diploma or GED") ~ "School",
                              cov_edu %in% c("Some college", "Associate degree (e.g., AA, AS)", "Bachelor's degree", "Graduate degree") ~ "College",
                              .default = cov_edu)) %>%
  mutate(cov_employ2 = case_when(cov_employ %in% c("Full-time employment", "Part-time employment") ~ "Employed",
                                cov_employ %in% c("Not in labor force", "Unemployed", "Retired", "Permanently unable to perform work activities due to medical condition", "On work experience") ~ "NotEmployed",
                                .default = cov_employ)) %>%
  mutate(cov_age_ep = as.numeric(difftime(diagnosis_date_mm, as.Date(birth_year), units = "days"))/365) %>%
  mutate(cov_age_ehr = as.numeric(difftime(as.Date(paste("2025-01-01")), as.Date(birth_year), units = "days"))/365) %>%
  dplyr::select(-birth_year)

# GEE to model repeated measures in Poisson (within-person correlation accounted for with clustered SE approach)
library(geepack)
hospital_visitdays_ttd_glm_cov_mania.glmgee <- geeglm(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_employ2 + cov_edu2 + cov_rel,# + cov_state
  id = person_id,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  data = hospital_visit_ttd.glm_prep %>% filter(mm_split == "mania")) %>%
  tidy() %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospital_visitdays_ttd_glm_cov_mania.glmgee, "hospital_visitdays_ttd_glm_cov_mania.glmgee", overwrite = T)

hospital_visitdays_ttd_glm_cov_mixed.glmgee <- geeglm(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_employ2 + cov_edu2  + cov_rel,#cov_state
  id = person_id,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  data = hospital_visit_ttd.glm_prep %>% filter(mm_split == "mixed")) %>%
  tidy() %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospital_visitdays_ttd_glm_cov_mixed.glmgee, "hospital_visitdays_ttd_glm_cov_mixed.glmgee", overwrite = T)

# Sensitivity: Negative binomial regression for counts (very overdispersed so can't do Poisson)
hospital_visitdays_ttd_glm_cov_mania.glmnb <- glm.nb(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_rel + cov_employ2 + cov_edu2 + cov_state, data = hospital_visit_ttd.glm_prep %>% filter(mm_split == "mania")) 
hospital_visitdays_ttd_glm_cov_mania.glmnb_out <- tidy(hospital_visitdays_ttd_glm_cov_mania.glmnb) %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospital_visitdays_ttd_glm_cov_mania.glmnb_out, "hospital_visitdays_ttd_glm_cov_mania.glmnb_out", overwrite = T)
hospital_visitdays_ttd_glm_cov_mixed.glmnb <- glm.nb(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_rel + cov_employ2 + cov_edu2 + cov_state, data = hospital_visit_ttd.glm_prep %>% filter(mm_split == "mixed")) 
hospital_visitdays_ttd_glm_cov_mixed.glmnb_out <- tidy(hospital_visitdays_ttd_glm_cov_mixed.glmnb) %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospital_visitdays_ttd_glm_cov_mixed.glmnb_out, "hospital_visitdays_ttd_glm_cov_mixed.glmnb_out", overwrite = T)

# Aggregate sumstats
hospital_visitdays_ttd_glm_cov_mania.allmodels <- rbind(
  data.frame(hospital_visitdays_ttd_glm_cov_mania.glmgee, mm_split = "mania", model = "Poisson_generalised_estimating_equation"),
  data.frame(hospital_visitdays_ttd_glm_cov_mixed.glmgee, mm_split = "mixed", model = "Poisson_generalised_estimating_equation"),
  data.frame(hospital_visitdays_ttd_glm_cov_mania.glmnb_out, mm_split = "mania", model = "negative_binomial"),
  data.frame(hospital_visitdays_ttd_glm_cov_mixed.glmnb_out, mm_split = "mixed", model = "negative_binomial"))

#----------------------------------
# TTD ~ depression hospital visits in 6 months following depression diagnosis
#----------------------------------

# Filter for visits within this dataframe which reflect visits associated with dep diagnosis within 7 days
bipolar_visit_hospital_dx_7d <- getDf("bipolar_visit_hospital_dx_7d")
bipolar_visit_hospital_dx_7d_add <- bipolar_visit_hospital_dx_7d %>% 
  filter(state == "dep") %>% 
  mutate(visit_type = ifelse(visit_type == "Emergency Room Visit", "Emergency", "Inpatient"), 
        depassoc = 1) %>% 
  dplyr::select(-c("state", "dx_n")) %>% distinct
bipolar_visit <- getDf("bipolar_visit")

# Try starting from full visits dataset and add an indicator for those with depression-assoc hospitalisation
bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost.tmp <- bipolar_epall_count_mm_to_deprecent_visit_dep6m_6mpost %>%
  ungroup() %>%
  filter(type == "Depression after") %>% filter(abs(monthstodep) <= 6) %>% filter(monthstodep >= 0) %>%
  left_join(., bipolar_visit_hospital_dx_7d_add) 

# - have already excluded visits that come before the end of the manic/mixed episode 
# - CHECK distribution of start of mania to start of visit
days_moodep_to_visitstart_check_inpatient <- bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost.tmp %>% 
  filter(visit_type == "Inpatient") %>%
  mutate(days_mmend_to_visitstart = as.numeric(difftime(visit_start_date, diagnosis_date_end_mm, units = "days"))/30,
         days_depstart_to_visitstart = as.numeric(difftime(visit_start_date, diagnosis_date_dep, units = "days"))/30)
summary(days_moodep_to_visitstart_check_inpatient$days_mmend_to_visitstart)
summary(days_moodep_to_visitstart_check_inpatient$days_depstart_to_visitstart) 
days_moodep_to_visitstart_check_inpatient %>% filter(days_depstart_to_visitstart < 0) %>% dim #4053 out of 55768
# from these 2 checks, all come after mm_end and at most 6 days before dep start

bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost <- bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost.tmp %>%
  # sum visits and visit duration
  group_by(person_id, mm_episode, dep_episode, state, mm_split, diagnosis_date_mm, diagnosis_date_end_mm, group, diagnosis_date_dep, diagnosis_date_end_dep, monthstodep, visit_type, depassoc) %>%
  summarise(VisitsN = sum(visit6m), VisitsDays = sum(visit_duration)) %>% ungroup %>%
  mutate(VisitsDays = ifelse(VisitsDays > 180, 180, VisitsDays)) # Some IDs have overlapping hospital visit dates eg. 100000574806 so >180 sum

# Separately do VisitsN and VisitsDays
bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_VisitsN <- bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost %>%
  dplyr::select(-VisitsDays) %>%
  mutate(visit_type = paste0(visit_type, "VisitsN")) %>%
  pivot_wider(names_from = c("visit_type", "depassoc"), values_from = "VisitsN") %>%
  mutate(InpatientVisitsN_1 = ifelse(is.na(InpatientVisitsN_1), 0, InpatientVisitsN_1),
         EmergencyVisitsN_1 = ifelse(is.na(EmergencyVisitsN_1), 0, EmergencyVisitsN_1)) %>%
  mutate(HospitalVisitsN = InpatientVisitsN_1 + EmergencyVisitsN_1) %>%
  dplyr::select(-c("EmergencyVisitsN_NA", "OutpatientVisitsN_NA", "InpatientVisitsN_NA", "EmergencyVisitsN_1", "InpatientVisitsN_1"))

bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_VisitsDays <- bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost %>%
  dplyr::select(-VisitsN) %>%
  mutate(visit_type = paste0(visit_type, "VisitsDays")) %>%
  pivot_wider(names_from = c("visit_type", "depassoc"), values_from = "VisitsDays") %>%
  mutate(InpatientVisitsDays_1 = ifelse(is.na(InpatientVisitsDays_1), 0, InpatientVisitsDays_1),
         EmergencyVisitsDays_1 = ifelse(is.na(EmergencyVisitsDays_1), 0, EmergencyVisitsDays_1)) %>%
  mutate(HospitalVisitsDays = InpatientVisitsDays_1 + EmergencyVisitsDays_1) %>%
  dplyr::select(-c("EmergencyVisitsDays_NA", "OutpatientVisitsDays_NA", "InpatientVisitsDays_NA", "EmergencyVisitsDays_1", "InpatientVisitsDays_1"))

bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_count <- bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_VisitsN %>%
  inner_join(., bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_VisitsDays)

bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_count.long <- bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_count %>%
  pivot_longer(cols = c("HospitalVisitsN", "HospitalVisitsDays"), names_to = "visit_type", values_to = "units")

# Plot scatter
bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_count.long %>%
  filter(visit_type == "HospitalVisitsDays") %>%
  ggplot(., aes(x = monthstodep, y = units)) +
  geom_jitter(width = 0.15, height = 0.1,
              size = 0.8, alpha = 0.1, color = "grey20") +
  theme_bw() + theme(legend.position = "none") +
  scale_fill_manual(values = millerstone[c(9,2,5)]) +
  xlab("Time to depression (months)") +
  ylab("Number of hospital days for depression within 6 months of depression onset") +
  facet_grid(visit_type ~ mm_split, scales = "free_y")
ggsave("hospitaldaysfordep_monthstodep_scatter.png")

# Prep data for regression
hospitaldep_visit_ttd.glm_prep0 <- bipolar_epall_count_mm_to_deprecent_hospfordep_dep6m_6mpost_count %>%
  mutate(hospital01 = ifelse(HospitalVisitsN == 0, 0, 1)) %>%
  left_join(., demog_cov)
cov_state_n.tmp <- hospitaldep_visit_ttd.glm_prep0 %>%
  group_by(person_id, cov_state) %>% summarise(n=n()) %>% group_by(cov_state) %>% summarise(n=n()) %>%
  mutate(cov_state2 = ifelse(n <= 20, "grouped", cov_state))
hospitaldep_visit_ttd.glm_prep <- hospitaldep_visit_ttd.glm_prep0 %>%
  left_join(., cov_state_n.tmp %>% dplyr::select(cov_state, cov_state2)) %>%
  rename(cov_state0 = cov_state, cov_state = cov_state2) %>%
  filter(cov_gender != "No matching concept") %>%
  mutate(cov_edu2 = case_when(cov_edu %in% c("8th grade or less", "Some high school, but did not graduate", "High school diploma or GED") ~ "School",
                              cov_edu %in% c("Some college", "Associate degree (e.g., AA, AS)", "Bachelor's degree", "Graduate degree") ~ "College",
                              .default = cov_edu)) %>%
  mutate(cov_employ2 = case_when(cov_employ %in% c("Full-time employment", "Part-time employment") ~ "Employed",
                                cov_employ %in% c("Not in labor force", "Unemployed", "Retired", "Permanently unable to perform work activities due to medical condition", "On work experience") ~ "NotEmployed",
                                .default = cov_employ)) %>%
  mutate(cov_age_ep = as.numeric(difftime(diagnosis_date_mm, as.Date(birth_year), units = "days"))/365) %>%
  mutate(cov_age_ehr = as.numeric(difftime(as.Date(paste("2025-01-01")), as.Date(birth_year), units = "days"))/365) %>%
  dplyr::select(-birth_year) %>%
  filter(monthstodep >= 0)

# GEE to model repeated measures in Poisson (within-person correlation accounted for with clustered SE approach)
# - glmer.nb takes too long to run
library(geepack)
hospitaldep_visitdays_ttd_glm_cov_mania.glmgee <- geeglm(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_employ2 + cov_edu2 + cov_rel,# + cov_state,
  id = person_id,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  data = hospitaldep_visit_ttd.glm_prep %>% filter(mm_split == "mania")) %>%
  tidy() %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospitaldep_visitdays_ttd_glm_cov_mania.glmgee, "hospitaldep_visitdays_ttd_glm_cov_mania.glmgee", overwrite = T)

hospitaldep_visitdays_ttd_glm_cov_mixed.glmgee <- geeglm(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_employ2 + cov_edu2 + cov_rel,# + cov_state,
  id = person_id,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  data = hospitaldep_visit_ttd.glm_prep %>% filter(mm_split == "mixed")) %>%
  tidy() %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospitaldep_visitdays_ttd_glm_cov_mixed.glmgee, "hospitaldep_visitdays_ttd_glm_cov_mixed.glmgee", overwrite = T)

# Sensitivity: Negative binomial regression for counts (very overdispersed so can't do Poisson)
hospitaldep_visitdays_ttd_glm_cov_mania.glmnb <- glm.nb(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_rel + cov_employ2 + cov_edu2 + cov_state, data = hospitaldep_visit_ttd.glm_prep %>% filter(mm_split == "mania")) 
hospitaldep_visitdays_ttd_glm_cov_mania.glmnb_out <- tidy(hospitaldep_visitdays_ttd_glm_cov_mania.glmnb) %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospitaldep_visitdays_ttd_glm_cov_mania.glmnb_out, "hospitaldep_visitdays_ttd_glm_cov_mania.glmnb_out", overwrite = T)
hospitaldep_visitdays_ttd_glm_cov_mixed.glmnb <- glm.nb(HospitalVisitsDays ~ monthstodep + cov_gender + cov_age_ep + cov_years_ehr + cov_rel + cov_employ2 + cov_edu2 + cov_state, data = hospitaldep_visit_ttd.glm_prep %>% filter(mm_split == "mixed")) 
hospitaldep_visitdays_ttd_glm_cov_mixed.glmnb_out <- tidy(hospitaldep_visitdays_ttd_glm_cov_mixed.glmnb) %>%
  mutate(estimate_IRR = exp(estimate),
        lowerIRR = exp(estimate - (std.error * 1.96)), 
        higherIRR = exp(estimate + (std.error * 1.96)))
saveDf(hospitaldep_visitdays_ttd_glm_cov_mixed.glmnb_out, "hospitaldep_visitdays_ttd_glm_cov_mixed.glmnb_out", overwrite = T)

# Aggregate sumstats
hospitaldep_visitdays_ttd_glm_cov_mania.allmodels <- rbind(
  data.frame(hospitaldep_visitdays_ttd_glm_cov_mania.glmgee, mm_split = "mania", model = "Poisson_generalised_estimating_equation"),
  data.frame(hospitaldep_visitdays_ttd_glm_cov_mixed.glmgee, mm_split = "mixed", model = "Poisson_generalised_estimating_equation"),
  data.frame(hospitaldep_visitdays_ttd_glm_cov_mania.glmnb_out, mm_split = "mania", model = "negative_binomial"),
  data.frame(hospitaldep_visitdays_ttd_glm_cov_mixed.glmnb_out, mm_split = "mixed", model = "negative_binomial"))

# Aggregate BOTH analyses
hospitaldep_visitdays_ttd_glm_cov_mania.hospall_hospdep_allmodels <- rbind(
  data.frame(hospitaldep_visitdays_ttd_glm_cov_mania.allmodels, analysis = "Depression-related hospital LoS"), 
  data.frame(hospital_visitdays_ttd_glm_cov_mania.allmodels, analysis = "All-cause hospital LoS")) %>%
  mutate(p_bonf = ifelse(term == "monthstodep", p.value*3, NA)) %>% # Bonferroni correction for 3 tests incl PHQ-9
  filter(term == "monthstodep")
saveDf(hospitaldep_visitdays_ttd_glm_cov_mania.hospall_hospdep_allmodels, "hospitaldep_visitdays_ttd_glm_cov_mania.hospall_hospdep_allmodels", overwrite = T)

hospitaldep_visitdays_ttd_glm_cov_mania.hospall_hospdep_allmodels %>% filter(!is.na(p_bonf))

#----------------------------------
# TTD ~ PHQ-9 depression severity
#----------------------------------

# Plot
horiz <- ggplot(phq9_max_mmend.long %>% filter(n_obs >= 2) %>% filter(months_postmm <= 6), aes(x = months_postmm, colour = mm_split)) +
  geom_density(alpha = 0.6) +
  theme_bw() +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  theme(axis.text.x=element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs(y = NULL, x = NULL)
vert <- ggplot(phq9_max_mmend.long %>% filter(n_obs >= 2) %>% filter(months_postmm <= 6), aes(x = score, colour = mm_split)) +
  geom_density(alpha = 0.6) +
  theme_bw() +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  theme(axis.text.x=element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs(y = NULL, x = NULL) +
  coord_flip()
ttd_phq9_scatter <- ggplot(phq9_max_mmend.long %>% filter(n_obs >= 2) %>% filter(months_postmm <= 6), aes(x = months_postmm, y = score, colour = mm_split)) +
  geom_point() +
  geom_smooth() +
  scale_colour_manual(values = millerstone[c(7,4)]) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Time to depression (months)") + ylab("Max PHQ-9 score") +
  scale_x_continuous(limits = c(0, 6))

# Combine plot
(horiz + plot_spacer() + ttd_phq9_scatter + vert) +
  plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))
ggsave("ttd_phq9severity_scatterplot.png")

# Statistical test
phq9_max_mmend.lm <- phq9_max_mmend.long %>%
  filter(months_postmm <= 6 & months_postmm > 0) %>%
  mutate(cov_race = ifelse(cov_race %in% c("Native Hawaiian or Other Pacific Islander", "American Indian or Alaska Native", "Asian"), "grouped", cov_race)) %>%
  nest(.by = c("mm_split")) %>%
  mutate(model = map(data, ~coeftest(lm(score ~ months_postmm + cov_age_ep + cov_years_ehr + cov_gender + cov_state + cov_rel + cov_race, data = .), vcov = vcovCL, cluster = ~person_id)), 
                 tidied = map(model, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-c("data", "model")) %>%
  group_by(mm_split) 
phq9_max_mmend.lm %>% filter(term == "months_postmm")

#======================================================
# Figure 1 plot
#======================================================

library(patchwork)

empty_plot <- ggplot() + theme_void()
combined <- bipolar_epall_count_mm_to_depafter_inf.plot + ttd_month.plot +
            patchwork::plot_annotation(tag_levels = "a")
ggsave("timetodep_figure_rates.svg", combined, height = 3.5, width = 9)

