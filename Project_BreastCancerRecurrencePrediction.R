library(GEOquery)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(randomForestSRC)
library(pec)
options(scipen = 999)

## load dataset
#gset_7390 <- getGEO("GSE7390", GSEMatrix=TRUE)
expr_7390 <- exprs(gset_7390[[1]])
clinical_7390 <- pData(gset_7390[[1]])

#gset_11121 <- getGEO("GSE11121", GSEMatrix=TRUE)
expr_11121 <- exprs(gset_11121[[1]])
clinical_11121 <- pData(gset_11121[[1]])

#BiocManager::install("breastCancerUNT")
#library(breastCancerUNT)
data(unt)
expr_oxfu <- exprs(unt)
expr_oxfu_df <-  as.data.frame(expr_oxfu)
clinical_oxfu <- pData(unt)

gset_2990 <- getGEO("GSE2990", GSEMatrix=TRUE)
expr_2990 <- exprs(gset_2990[[1]])
clinical_2990_poor <- pData(gset_2990[[1]])

## Transform
#clean dataset 7390 
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

# clean dataset 11121
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

#clean dataset oxfu
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
clinical_cleaned_oxfu_naomt <- clinical_cleaned_oxfu %>%
  na.omit()

#clean dataset gse2990

# replace "?" with NA
clinical_2990_poor[clinical_2990_poor == "?"] <- NA

#load clinical_2990 from txt file with has more info.
file_suppl = "G:/clinical_data_science/final project/breast cancer/GSE2990_suppl_info.txt"
clinical_2990_good <- read.table(file_suppl, fileEncoding = "UTF-8", header = TRUE)

# fills missing (NA) values in clinical_2990_good$grade using corresponding
# non-missing values from clinical_2990_poor$"grade:ch1", matching rows by geo_accn
for (i in which(is.na(clinical_2990_good$grade))) {
  geo <- clinical_2990_good$geo_accn[i]
  if (geo %in% rownames(clinical_2990_poor)) {
    if (!is.na(clinical_2990_poor[geo, "grade:ch1"])) {
      clinical_2990_good$grade[i] <- clinical_2990_poor[geo, "grade:ch1"]
    }
  }
}

col_selected_2990 <- c(
  "geo_accn",
  "time.dmfs",
  "event.dmfs",
  "grade",
  "er",
  "age",
  "size",
  "node",
  "treatment"
)
clinical_cleaned_2990 <- clean_clinical_data(col_selected_2990, df_clinic_raw = clinical_2990_good)
colSums(is.na(clinical_cleaned_2990))
# since GSE7390 patients didn't receive szstemic therapy and lympha-note negative
clinical_cleaned_2990_naomt <- clinical_cleaned_2990 %>%
  na.omit() %>%
  filter(treatment == "none" & node == "0")


# Box plot of continuos variable
continuous_variables <- c("size", "age", "t_dmfs")
generate_box_plots(clinical_cleaned_7390, continuous_variables)

selected_col <- c("grade", "age", "er", "size")
missed_col <- "grade"
clinical_cleaned_7390_imputed <- impute_missing_value(clinical_cleaned_7390, selected_col, missed_col)
# 
# selected_col <- c("grade", "age", "er", "size")
# missed_col <- "er"
# clinical_cleaned_oxfu_imputed <- impute_missing_value(clinical_cleaned_oxfu, selected_col, missed_col)
# 
# missed_col <- "grade"
# clinical_cleaned_oxfu_imputed <- impute_missing_value(clinical_cleaned_oxfu_imputed, selected_col, missed_col)


#build and test a Cox proportional hazards model
predictors <- c("grade", "size", "age", "er")
results <- fit_cox_model(predictors, clinical_cleaned_7390_imputed)
results <- fit_cox_model(predictors, clinical_cleaned_oxfu_naomt)



predictors <- c("grade", "size")
results <- fit_cox_model(predictors, clinical_cleaned_11121)

# Access model or proportional hazards test if needed:
# results$model
# results$ph_test
# ggcoxzph(results$ph_test)

# validate the model
calculate_time_auc_cindex("Cox", fitted_model = results$model, df = clinical_cleaned_7390_imputed)
# Evalute with external data
calculate_time_auc_cindex("Cox", fitted_model = results$model, df = clinical_cleaned_oxfu_naomt)
calculate_time_auc_cindex("Cox", fitted_model = results$model, df = clinical_cleaned_2990_naomt)

# Plot KM curve by 'grade'
plot_km_by_group(clinical_cleaned_7390_imputed, group_var = "er")
plot_km_by_group(clinical_cleaned_oxfu_naomt, group_var = "er")
plot_km_by_group(clinical_cleaned_11121, group_var = "grade")


### randomforest
predictors <- c("grade", "size", "age", "er")
clinical_rsf <- clinical_cleaned_7390_imputed[, c("t_dmfs", "e_dmfs", predictors)]
# 构建 RSF 模型
result <- rsf_kfold_cv_best(clinical_rsf, K = 5)
# Best model
rsf_fit <- result$best_model
# Variable importance from retrained model with complete data
importance <- result$importance

# evaluate the model
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_7390_imputed)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_oxfu_naomt)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_2990_naomt)

# rebuild the model with non negative importance predictors
predictors <- c("size", "grade")
clinical_rsf <- clinical_cleaned_7390_imputed[, c("t_dmfs", "e_dmfs", predictors)]
# 构建 RSF 模型
rsf_fit <- rfsrc(Surv(t_dmfs, e_dmfs) ~ ., data = clinical_rsf, ntree = 1000, importance = TRUE)
rsf_fit$importance

# For RSF model
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_7390_imputed)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_oxfu_naomt)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_2990_naomt)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_11121)
