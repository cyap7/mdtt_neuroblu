############

## Objective: Load all necessary libraries and associated scripts

# Load all the packages
library('NeuroBlu')
library('dplyr')
library('plyr')
library('tidyr')
library('lubridate')
library('survival')
library('survminer')
library('tidyverse')

# Palettes
millerstone <- c("#4F6980FF", "#849DB1FF", "#A2CEAAFF", "#638B66FF", "#BFBB60FF", "#F47942FF", "#FBB04EFF", "#B66353FF", "#D7CE9FFF", "#B9AA97FF", "#7E756DFF", "gray70")

############

## Objective: Setup simple necessary functions

addICDdecimal=function(ICD){
  res=c()
  for (i in seq_along(ICD)){
    res=c(res,ICD[i])
    for (k in 0:9){
      ICDik=paste(ICD[i],".",as.character(k),sep='')
      res=c(res,ICDik)
    }  
  }
  return(res)
}


addICDdigit=function(ICD){
  res=c()
  for (i in seq_along(ICD)){
    res=c(res,ICD[i])
    for (k in 0:9){
      ICDik=paste(ICD[i],as.character(k),sep='')
      res=c(res,ICDik)
    }  
  }
  return(res)
}

getPatients=function(icdCodes){
  allPatients=data.frame(diagnosis=NULL,person_id=NULL)
  for (c in names(icdCodes)){
    pt=filterPatientsByICD(icdCodes[[c]])
    pt_dt=data.frame(diagnosis=c,person_id=pt)
    allPatients = allPatients %>% bind_rows(pt_dt)
  }
  allPatients$diagnosis=as.factor(allPatients$diagnosis)
  return(allPatients)
}

# Clean up cases with both ICD-9 and ICD-10 codes given on the same day
# - if ICD-10 F-code is present, then keep that. Otherwise, take ICD-9
# - input requires: columns with person_id, diagnosis_date, icd_code
remove_icd_dup <- function(x) {
  dx_id_count <- x %>%
    arrange(diagnosis_date) %>%
    group_by(person_id, diagnosis_date) %>%
    summarise(n=n()) %>% ungroup()
  dx_id_count_multi <- dx_id_count %>% filter(n >= 2) %>% 
    mutate(tmpid = paste0(person_id, diagnosis_date))
  # by removing duplicated ICD-9/ICD-10 codes for the same episode, cuts down from 88350 to 46417 rows
  rmmultiicd <- x %>% 
    mutate(tmp = paste0(person_id, diagnosis_date)) %>%
    filter(tmp %in% dx_id_count_multi$tmpid) %>%
    dplyr::select(-"tmp") %>%
    mutate(icd_num = ifelse(substring(icd_code, 1, 1) == "F", "ICD-10", "ICD-9")) %>%
    arrange(diagnosis_date) %>%
    group_by(person_id, diagnosis_date) %>%
    mutate(icd_rm = case_when(length(unique(icd_num)) == 2 & icd_num == "ICD-10" ~ "keep",
                       length(unique(icd_num)) == 2 & icd_num == "ICD-9" ~ "rm",
                       length(unique(icd_num)) == 1 ~ "keep")) %>%
    ungroup() %>%
    filter(icd_rm == "rm") %>%
    dplyr::select(-c("icd_num")) %>%
    mutate(tmpid = paste0(person_id, diagnosis_date, icd_code))
  return(rmmultiicd)
}

# Function to extract glm coefficients from glm.fit object
tidy_glmfit <- function(glm_fit_obj) {
  # Extract coefficients
  coefs <- glm_fit_obj$coefficients
  
  # Calculate covariance matrix
  R <- glm_fit_obj$qr$qr[1:glm_fit_obj$qr$rank, 1:glm_fit_obj$qr$rank, drop = FALSE]
  R[lower.tri(R)] <- 0
  covar <- chol2inv(R)
  
  # Standard errors
  se <- sqrt(diag(covar))
  
  # Test statistics and p-values
  z_stat <- coefs / se
  p_val <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
  
  tibble(
    term = names(coefs),
    estimate = coefs,
    std.error = se,
    statistic = z_stat,
    p.value = p_val
  )
}
