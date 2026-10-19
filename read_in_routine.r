# Commonly used dataframes
# PMMD
pmmd_group <- getDf("pmmd_group")
pmmd_ep <- getDf("pmmd_ep")
pmmd_epall <- getDf("pmmd_epall")
pmmd_ep_temporal <- getDf("pmmd_ep_temporal")

# noPMMD
nopmmd_group <- getDf("nopmmd_group")
nopmmd_epall <- getDf("nopmmd_epall")

# MMonly
mmonly_epall <- getDf("mmonly_epall")
mmonly_count_toepisodes <- getDf("mmonly_count_toepisodes")

# id_group
id_group <- getDf("id_group")
id_group <- id_group %>% distinct()

# episodes
mania_dep_dates_toepisodes_qc_broaddef <- getDf("mania_dep_dates_toepisodes_qc_broaddef")

# covariates
demog_cov <- getDf("demog_cov")

# get subset of those with episodes under 30 (ie. first ever mood episode)
#icd_firstepmood_u30 <- getDf("icd_firstepmood_u30")

# mm_split_mmdep
mm_split <- getDf("mm_split")
mm_split_mmdep <- getDf("mm_split_mmdep")