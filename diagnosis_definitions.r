#========================================
# Define diagnosis -> mood classification
#========================================

#Generate full vector of ICD codes eg addICDdecimal("F32") would return F32, F32.0, F32.1 etc by defining a "addIcDecimal" function
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

#Generate full vector of ICD codes as ICDdecimal but without decimal as some datasets record without decimal. F32 to F320, F321 etc.
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

mania1 <- c('296', '296.0', '296.00', '296.01', '296.02', '296.03', '296.04', '296.05', '296.4', '296.40', '296.41', '296.42', '296.43', '296.44', '296.45')
mania2 <- c('F30', 'F30.0', 'F30.10', 'F30.11', 'F30.12', 'F30.13', 'F30.8', 'F30.9', 'F30.1', 'F30.2', 'F31.0', 'F31.10', 'F31.1', 'F31.11', 'F31.12', 'F31.13', 'F31.2', 'F31.71', 'F31.73')
mixed1 <- c('296.6', '296.60', '296.61', '296.62', '296.63', '296.64' ,'296.65')
mixed2 <- c('F30.6', 'F31.60', 'F31.61', 'F31.62', 'F31.63', 'F31.64', 'F31.77')
mm <- c(mania1, mania2, mixed1, mixed2)

bpad_dep1 <- c('296.50', '296.51', '296.52', '296.53', '296.54', '296.55')
bpad_dep2 <- c('F31.30', 'F31.31', 'F31.32', 'F31.4', 'F31.5', 'F31.75')
mdd1 <- c('296.20', '296.21', '296.22', '296.23', '296.24', '296.25', '296.30', '296.31', '296.32', '296.33', '296.34', '296.35') # include in case of miscoding
mdd2 <- c('F32.0', 'F32.1', 'F32.2', 'F32.3', 'F32.4', 'F32.8', 'F32.9', 'F33.0', 'F33.1', 'F33.2', 'F33.3', 'F33.4', 'F33.8', 'F33.9')
dep <- c(bpad_dep1, bpad_dep2, mdd1, mdd2)
#dep <- c(bpad_dep2, mdd2)

# severity specifiers
mm_psychosis <- c("F31.2", "F31.5", 'F31.64', 'F30.2', "296.04", "296.14", "296.44", "296.64")
dep_psychosis <- c('F32.3', 'F33.3', "296.24", "296.34", '296.54')

mm_hypo <- c('F30.0', 'F31.0') # 'F31.71'
mm_mild <- c('F31.11', '296.01', '296.41', 'F31.61', '296.61')
#mm_mod <- c('F30.1', 'F31.12', 'F31.73', '296.02', '296.05', '296.42', 'F31.62', '296.62')
mm_mod <- c('F30.12', 'F31.12', '296.02', '296.42', 'F31.62', '296.62') # 'F31.73', '296.05', 
mm_severe <- c('F30.13', 'F31.13', '296.03', '296.43', 'F31.63', '296.63')
#mm_unspec <- c('F31.10', '296', '296.40', '296.45', 'F31.60', 'F31.77', '296.6', '296.60', '296.65') ##### was 295.45 = schizophreniform, so need to change to bpad!!
mm_unspec <- c('F30', 'F30.1', 'F30.10', 'F30.6', 'F30.8', 'F30.9', 'F31.10', '296', '296.0', '296.00', '296.4', '296.40', 'F31.60', '296.6', '296.60') # 'F31.77', '296.45', '296.65'
dep_unspec <- c('F31.75', 'F32.4', 'F32.8', 'F32.9', 'F33.4', 'F33.8', 'F33.9', '296.50', '296.55', '296.20', '296.25', '296.30', '296.35')
dep_mild <- c('F31.31', 'F32.0', '296.51', 'F33.0', '296.21', '296.31')
dep_mod <- c('F31.30', 'F31.32', 'F32.1', 'F33.1', '296.52', '296.22', '296.32')
dep_severe <- c('F31.5', 'F32.2', 'F33.2', '296.53', '296.23', '296.33')

icdseverity_score <- rbind(
data.frame(icd_code = mm_hypo, icdseverity = "1hypomania", icdseverity_rank = 1),
data.frame(icd_code = mm_mild, icdseverity = "2mild", icdseverity_rank = 2),
data.frame(icd_code = mm_mod, icdseverity = "3moderate", icdseverity_rank = 3),
data.frame(icd_code = mm_severe, icdseverity = "4severe", icdseverity_rank = 4),
data.frame(icd_code = mm_psychosis, icdseverity = "5psychosis", icdseverity_rank = 5), 
data.frame(icd_code = mm_unspec, icdseverity = "unspecified", icdseverity_rank = 0),
data.frame(icd_code = dep_mild, icdseverity = "1mild", icdseverity_rank = -1),
data.frame(icd_code = dep_mod, icdseverity = "2moderate", icdseverity_rank = -2),
data.frame(icd_code = dep_severe, icdseverity = "3severe", icdseverity_rank = -3),
data.frame(icd_code = dep_psychosis, icdseverity = "4psychosis", icdseverity_rank = -4),
data.frame(icd_code = dep_unspec, icdseverity = "unspecified", icdseverity_rank = 0))

mania_fullremission <- c('296.06', '296.46', '296.66', 'F31.72', 'F31.74', 'F31.78')

bpad_dep_fullremission <- c('296.56')
mdd_fullremission <- c('296.26', '296.36', 'F32.5', 'F33.42', 'F31.76')
dep_remission <- c(bpad_dep_fullremission, mdd_fullremission)

# Anxiety
anxiety_icd10 <- c(
  "F40.00", "F40.01", "F40.02",
  "F40.10", "F40.11",
  "F40.210", "F40.218", "F40.220", "F40.228",
  "F40.230", "F40.231", "F40.232", "F40.233",
  "F40.240", "F40.241", "F40.242", "F40.243", "F40.248",
  "F40.290", "F40.291", "F40.298",
  "F40.8", "F40.9",
  "F41.0", "F41.1", "F41.3", "F41.8", "F41.9")
anxiety_icd9 <- c(
  "300.00", "300.01", "300.02", "300.09",
  "300.20", "300.21", "300.22", "300.23", "300.29")

anxiety <- c(anxiety_icd10, anxiety_icd9)

# - helpful https://psychcentral.com/bipolar/dsm-5-diagnostic-codes-for-bipolar-disorder#patient-resources

# Suicide
#1) Diagnosis codes
suicide_codes <- c(
   # Suicidal ideation
   "R45.851",   # ICD-10-CM
   "V62.84",    # ICD-9-CM
   # Suicide attempt / intentional self-harm
   "T14.91",
   "X71", "X72", "X73", "X74", "X75", "X76", "X77",
   "X78", "X79", "X80", "X81", "X82", "X83",
   "E950", "E951", "E952", "E953", "E954",
   "E955", "E956", "E957", "E958")

# 2) NLP-based identification
suicide_nlp <- list(
          disorder_groups = c("Depressive Disorders"),
          features = c("suicidal thoughts", "suicide attempt", "suicidal behavior"),
          exclude_values = c("no evidence of", "2"))

# OLD 251209

#mania1 <- c('296', '296.01', '296.02', '296.03', '296.04', '296.05', '296.40', '296.41', '296.42', '296.43', '296.44', '295.45')
#mania2 <- c('F30.0', 'F30.1', 'F30.2', 'F31.0', 'F31.10', 'F31.11', 'F31.12', 'F31.13', 'F31.2', 'F31.71', 'F31.73')
#mixed1 <- c('296.6', '296.60', '296.61', '296.62', '296.63', '296.64' ,'296.65')
#mixed2 <- c('F31.60', 'F31.61', 'F31.62', 'F31.63', 'F31.64', 'F31.77')
#mm <- c(mania1, mania2, mixed1, mixed2)
#
#bpad_dep1 <- c('296.50', '296.51', '296.52', '296.53', '296.54', '296.55')
#bpad_dep2 <- c('F31.30', 'F31.31', 'F31.32', 'F31.4', 'F31.5', 'F31.75')
#mdd1 <- c('296.20', '296.21', '296.22', '296.23', '296.24', '296.25', '296.30', '296.31', '296.32', '296.33', '296.34', '296.35') # include in case of miscoding
#mdd2 <- c('F32.0', 'F32.1', 'F32.2', 'F32.3', 'F32.4', 'F32.8', 'F32.9', 'F33.0', 'F33.1', 'F33.2', 'F33.3', 'F33.4', 'F33.8', 'F33.9')
#dep <- c(bpad_dep1, bpad_dep2, mdd1, mdd2)
##dep <- c(bpad_dep2, mdd2)
#
#psychosis1 <- c("F31.2", "F31.5", "F32.3", 'F31.64')
#psychosis2 <- c("296.04", "296.14", "296.24", "296.34", "296.44", "296.54", "296.64")
#psychosis <- c(psychosis1, psychosis2)
#
#mania_fullremission <- c('296.06', '296.46', '296.66', 'F31.72', 'F31.74', 'F31.78')
#
#bpad_dep_fullremission <- c('296.56')
#mdd_fullremission <- c('296.26', '296.36', 'F32.5', 'F33.42', 'F31.76')
#dep_remission <- c(bpad_dep_fullremission, mdd_fullremission)
