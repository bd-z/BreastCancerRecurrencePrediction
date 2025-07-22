library(GEOquery)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
options(scipen = 999)


gset_7390 <- getGEO("GSE7390", GSEMatrix=TRUE)
expr_7390 <- exprs(gset_7390[[1]])
clinical_7390 <- pData(gset_7390[[1]])


gset_11121 <- getGEO("GSE11121", GSEMatrix=TRUE)
expr_11121 <- exprs(gset_11121[[1]])
clinical_11121 <- pData(gset_11121[[1]])


BiocManager::install("breastCancerUNT")
library(breastCancerUNT)
data(unt)
expr_oxfu <- exprs(unt)
clinical_oxfu <- pData(unt)



# 
col_selected_7390 <- c(
  "geo_accession",
  "t.dmfs:ch1",
  "e.dmfs:ch1",
  "grade:ch1",
  "er:ch1",
  "age:ch1",
  "size:ch1",
  "node:ch1"
)

#selects and cleans clinical columns, standardizes column names, converts values
# to numeric where possible, and transforms low-cardinality columns to factors.
clinical_cleaned_7390 <- clean_clinical_data(col_selected_7390, df_clinic_raw = clinical_7390)

colSums(is.na(clinical_cleaned_7390))

col_selected_11121 <- c(
  "geo_accession",
  "t.dmfs:ch1",
  "e.dmfs:ch1",
  "grade:ch1",
  #"er:ch1",
  #"age:ch1",
  "size_in_cm:ch1",
  "node:ch1"
)

clinical_cleaned_11121 <- clean_clinical_data(col_selected_11121, df_clinic_raw = clinical_11121)
colSums(is.na(clinical_cleaned_11121))


col_selected_oxfu <- c(
  "samplename",
  "t.dmfs",
  "e.dmfs",
  "grade",
  "er",
  "age",
  "size",
  "node"
)

clinical_cleaned_oxfu <- clean_clinical_data(col_selected_oxfu, df_clinic_raw = clinical_oxfu)
colSums(is.na(clinical_cleaned_oxfu))


# Box plot of continuos variable
continuous_variables <- c("size", "age", "t_dmfs")
generate_box_plots(clinical_cleaned_7390, continuous_variables)

selected_col <- c("grade", "age", "er", "size")
missed_col <- "grade"
clinical_cleaned_7390_imputed <- impute_missing_value(clinical_cleaned_7390, selected_col, missed_col)

selected_col <- c("grade", "age", "er", "size")
missed_col <- "er"
clinical_cleaned_oxfu_imputed <- impute_missing_value(clinical_cleaned_oxfu, selected_col, missed_col)

missed_col <- "grade"
clinical_cleaned_oxfu_imputed <- impute_missing_value(clinical_cleaned_oxfu_imputed, selected_col, missed_col)
