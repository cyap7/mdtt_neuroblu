#====================================================
# Combine pmmd and nopmmd dataframes into bipolar df
# Also ensure that QC filters are applied in 1 script
#====================================================

mm0_mmenddepstart0_dep0_toexclude <- getDf("mm0_mmenddepstart0_dep0_toexclude")
bipolar_epall_rm_qc <- getDf("bipolar_epall_rm_qc")

#--------------------------------
# Visit overlap
#--------------------------------
pmmd_epall_any_overlap <- getDf("pmmd_epall_any_overlap")
pmmd_epall_hospital_overlap <- getDf("pmmd_epall_hospital_overlap")
nopmmd_epall_any_overlap <- getDf("nopmmd_epall_any_overlap")
nopmmd_epall_hospital_overlap <- getDf("nopmmd_epall_hospital_overlap")

bipolar_epall_any_overlap <- full_join(pmmd_epall_any_overlap, nopmmd_epall_any_overlap) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_epall_hospital_overlap <- full_join(pmmd_epall_hospital_overlap, nopmmd_epall_hospital_overlap) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
saveDf(bipolar_epall_any_overlap, "bipolar_epall_any_overlap", overwrite = T)
saveDf(bipolar_epall_hospital_overlap, "bipolar_epall_hospital_overlap", overwrite = T)

# Medications
#--------------------------------
# All individuals, NO VISIT CONSTRAINT
#--------------------------------

pmmdgroup_mm_7d_dx_psychmed <- getDf("pmmdgroup_mm_7d_dx_psychmed")
pmmdgroup_mm_7d_dx_psychmeds_count <- getDf("pmmdgroup_mm_7d_dx_psychmeds_count")
pmmdgroup_mm_7d_dx_psychmeds_count.long <- getDf("pmmdgroup_mm_7d_dx_psychmeds_count.long")

pmmdgroup_mm_7d_dx_medsind_count <- getDf("pmmdgroup_mm_7d_dx_medsind_count")
pmmdgroup_mm_7d_dx_medsind_count.long <- getDf("pmmdgroup_mm_7d_dx_medsind_count.long")

nopmmdgroup_mm_7d_dx_psychmed <- getDf("nopmmdgroup_mm_7d_dx_psychmed")
nopmmdgroup_mm_7d_dx_psychmeds_count <- getDf("nopmmdgroup_mm_7d_dx_psychmeds_count")
nopmmdgroup_mm_7d_dx_psychmeds_count.long <- getDf("nopmmdgroup_mm_7d_dx_psychmeds_count.long")

nopmmdgroup_mm_7d_dx_medsind_count <- getDf("nopmmdgroup_mm_7d_dx_medsind_count")
nopmmdgroup_mm_7d_dx_medsind_count.long <- getDf("nopmmdgroup_mm_7d_dx_medsind_count.long")

# - broad classes
bipolar_mm_7d_dx_psychmed <- full_join(pmmdgroup_mm_7d_dx_psychmed, nopmmdgroup_mm_7d_dx_psychmed) %>% 
  filter(!is.na(broad_class)) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_dx_psychmeds_count <- full_join(pmmdgroup_mm_7d_dx_psychmeds_count, nopmmdgroup_mm_7d_dx_psychmeds_count) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_dx_psychmeds_count.long <- full_join(pmmdgroup_mm_7d_dx_psychmeds_count.long, nopmmdgroup_mm_7d_dx_psychmeds_count.long) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
saveDf(bipolar_mm_7d_dx_psychmed, "bipolar_mm_7d_dx_psychmed", overwrite = T)
saveDf(bipolar_mm_7d_dx_psychmeds_count, "bipolar_mm_7d_dx_psychmeds_count", overwrite = T)
saveDf(bipolar_mm_7d_dx_psychmeds_count.long, "bipolar_mm_7d_dx_psychmeds_count.long", overwrite = T)

# - individual meds
bipolar_mm_7d_dx_medsind_count <- full_join(pmmdgroup_mm_7d_dx_medsind_count, nopmmdgroup_mm_7d_dx_medsind_count) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_dx_medsind_count.long <- full_join(pmmdgroup_mm_7d_dx_medsind_count.long, nopmmdgroup_mm_7d_dx_medsind_count.long) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
saveDf(bipolar_mm_7d_dx_medsind_count, "bipolar_mm_7d_dx_medsind_count", overwrite = T)
saveDf(bipolar_mm_7d_dx_medsind_count.long, "bipolar_mm_7d_dx_medsind_count.long", overwrite = T)

#--------------------------------
# Hospitalisation
#--------------------------------

pmmdgroup_mm_7d_hosp_psychmed <- getDf("pmmdgroup_mm_7d_hosp_psychmed")
nopmmdgroup_mm_7d_hosp_psychmed <- getDf("nopmmdgroup_mm_7d_hosp_psychmed")

pmmdgroup_mm_7d_hosp_psychmeds_count <- getDf("pmmdgroup_mm_7d_hosp_psychmeds_count")
nopmmdgroup_mm_7d_hosp_psychmeds_count <- getDf("nopmmdgroup_mm_7d_hosp_psychmeds_count")

pmmdgroup_mm_7d_hosp_psychmeds_count.long <- getDf("pmmdgroup_mm_7d_hosp_psychmeds_count.long")
nopmmdgroup_mm_7d_hosp_psychmeds_count.long <- getDf("nopmmdgroup_mm_7d_hosp_psychmeds_count.long")

pmmdgroup_mm_7d_nohosp_psychmed <- getDf("pmmdgroup_mm_7d_nohosp_psychmed")
nopmmdgroup_mm_7d_nohosp_psychmed <- getDf("nopmmdgroup_mm_7d_nohosp_psychmed")

pmmdgroup_mm_7d_hosp_medsind_count <- getDf("pmmdgroup_mm_7d_hosp_medsind_count")
nopmmdgroup_mm_7d_hosp_medsind_count <- getDf("nopmmdgroup_mm_7d_hosp_medsind_count")
pmmdgroup_mm_7d_hosp_medsind_count.long <- getDf("pmmdgroup_mm_7d_hosp_medsind_count.long")
nopmmdgroup_mm_7d_hosp_medsind_count.long <- getDf("nopmmdgroup_mm_7d_hosp_medsind_count.long")

bipolar_mm_7d_hosp_psychmed <- full_join(pmmdgroup_mm_7d_hosp_psychmed, nopmmdgroup_mm_7d_hosp_psychmed) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_hosp_psychmeds_count <- full_join(pmmdgroup_mm_7d_hosp_psychmeds_count, nopmmdgroup_mm_7d_hosp_psychmeds_count) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_hosp_psychmeds_count.long <- full_join(pmmdgroup_mm_7d_hosp_psychmeds_count.long, nopmmdgroup_mm_7d_hosp_psychmeds_count.long) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_hosp_medsind_count <- full_join(pmmdgroup_mm_7d_hosp_medsind_count, nopmmdgroup_mm_7d_hosp_medsind_count) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_hosp_medsind_count.long <- full_join(pmmdgroup_mm_7d_hosp_medsind_count.long, nopmmdgroup_mm_7d_hosp_medsind_count.long) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)
bipolar_mm_7d_nohosp_psychmed <- full_join(pmmdgroup_mm_7d_nohosp_psychmed, nopmmdgroup_mm_7d_nohosp_psychmed) %>%
  anti_join(., mm0_mmenddepstart0_dep0_toexclude.tmp) %>%
  anti_join(., bipolar_epall_rm_qc)

saveDf(bipolar_mm_7d_hosp_psychmed, "bipolar_mm_7d_hosp_psychmed", overwrite = T)
saveDf(bipolar_mm_7d_hosp_psychmeds_count, "bipolar_mm_7d_hosp_psychmeds_count", overwrite = T)
saveDf(bipolar_mm_7d_hosp_psychmeds_count.long, "bipolar_mm_7d_hosp_psychmeds_count.long", overwrite = T)
saveDf(bipolar_mm_7d_hosp_medsind_count, "bipolar_mm_7d_hosp_medsind_count", overwrite = T)
saveDf(bipolar_mm_7d_hosp_medsind_count.long, "bipolar_mm_7d_hosp_medsind_count.long", overwrite = T)
saveDf(bipolar_mm_7d_nohosp_psychmed, "bipolar_mm_7d_nohosp_psychmed", overwrite = T)
