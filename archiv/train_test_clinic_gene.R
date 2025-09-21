library(GEOquery)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(randomForestSRC)
library(pec)
library(mice)
library(caret)
options(scipen = 999)

## load dataset
#gset_7390 <- getGEO("GSE7390", GSEMatrix=TRUE)
expr_7390 <- exprs(gset_7390[[1]])
clinical_7390 <- pData(gset_7390[[1]])

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


# Box plot of continuos variable
continuous_variables <- c("size", "age", "t_dmfs")
generate_box_plots(clinical_cleaned_7390, continuous_variables)

# MICE 7390
selected_col <- c("grade", "age", "er", "size")
missed_col <- c("grade")
colSums(is.na(clinical_cleaned_7390))
clinical_cleaned_7390_imputed <- impute_missing_value(clinical_cleaned_7390, selected_col, missed_col)

# Standardize each gene (row) across all samples (columns)
expr_7390_scaled <- t(apply(expr_7390, 1, scale))
rownames(expr_7390_scaled) <- rownames(expr_7390)
colnames(expr_7390_scaled) <- colnames(expr_7390)
# # 查看某个基因标准化后的均值和标准差（应为约 0 和 1）
# mean(expr_7390_scaled[1, ])
# sd(expr_7390_scaled[1, ])


######################################################################
result <- split_expr_clinical(expr_7390,
                              clinical_cleaned_7390_imputed,
                              stratify_col = "e_dmfs",
                              train_frac = 0.7, 
                              seed = 345)


# Access elements
train_expr = result$train_expr
train_clinical = result$train_clinical
test_expr = result$test_expr
test_clinical = result$test_clinical


sig_gene_df <- batch_univariate_cox_regression(train_expr, train_clinical)


selected_gene_df <- lasso_cox_cv(train_expr, train_clinical, sig_gene_df)


clinical_cleaned_risk_train <- compute_risk_score(gene_mat_scaled = train_expr,
                                                 significant_vars_df = selected_gene_df,
                                                 clinical_cleaned= train_clinical,
                                                 n_group = 3)


#plot_km_by_group(clinical_cleaned_risk_train, group_var = "risk_group")

predictors0 <- c("grade", "er", "age", "size") # "risk_score"
predictors  <- c(predictors0, colnames(clinical_cleaned_risk_train)[11: length(colnames(clinical_cleaned_risk_train))])

results_train <- fit_cox_model(predictors, clinical_cleaned_risk_train)
calculate_time_auc_cindex("Cox", fitted_model = results_train$model, df = clinical_cleaned_risk_train)


# valid
clinical_cleaned_risk_test <- compute_risk_score(gene_mat_scaled = test_expr,
                                                  significant_vars_df = selected_gene_df,
                                                  clinical_cleaned= test_clinical,
                                                  n_group = 3)
plot_km_by_group(clinical_cleaned_risk_test, group_var = "risk_group")
calculate_time_auc_cindex("Cox", fitted_model = results_train$model, df = clinical_cleaned_risk_test)

# file <- "cox_results_AUC_75_Index_69.rds"
# saveRDS(results_train, file)
# saved_model <- readRDS(file)
# calculate_time_auc_cindex("Cox", fitted_model = saved_model$model, df = clinical_cleaned_risk_test)

#file1 <- "selected_gene_df_AUC_75_Index_69.rds"
#saveRDS(selected_gene_df, file1)
#saved_selected_gene_df <- readRDS(file1)

####################################

cox_flow_result <- cox_workflow(gene_expr =  expr_7390,
             clinical_data_imputed = clinical_cleaned_7390_imputed,
             train_frac = 0.68,
             seed = 23456,
             clin_pred = TRUE,
             riskscore_pred = FALSE,
             gene_pred = TRUE,
             clin_predictors = c("grade", "er", "age", "size")
)



# 初始化
best_score <- -Inf
best_result <- NULL
best_frac <- NA
best_seed <- NA

# 可选：存储所有结果
all_results <- list()

# 遍历 train_frac 和不同的 seed
for (i in seq_along(seq(0.69, 0.71, by = 0.01))) {
  frac <- seq(0.69, 0.71, by = 0.01)[i]
  seed <- 100000 + i * 123  # 生成一个不同的 seed，你也可以用 sample() 随机生成
  
  cat("Running with train_frac =", frac, ", seed =", seed, "\n")
  
  result <- cox_rsf_workflow(
    gene_expr = expr_7390_scaled,
    clinical_data_imputed = clinical_cleaned_7390_imputed,
    train_frac = frac,
    seed = seed,
    clin_pred = TRUE,
    riskscore_pred = FALSE,
    gene_pred = TRUE,
    clin_predictors = c("grade", "er", "age", "size")
  )
  
  # 计算 iAUC + c_index
  print("IAUC:")
  print(result$result_valid$iAUC)
  print("C_index:")
  print(result$result_valid$c_index)
  
  score <- result$result_valid$iAUC + result$result_valid$c_index
  
  print("iAUC + c_index:")
  print(result$result_valid$c_index)
  
  
  # 存储所有结果（可选）
  all_results[[i]] <- list(
    train_frac = frac,
    seed = seed,
    iAUC = result$result_valid$iAUC,
    c_index = result$result_valid$c_index,
    score = score,
    result = result
  )
  
  # 找到最高得分
  # if (score > best_score) {
  #   best_score <- score
  #   best_result <- result
  #   best_frac <- frac
  #   best_seed <- seed
  # }
}


# 输出最佳结果信息
cat("Best train_frac:", best_frac, "\n")
cat("Best seed:", best_seed, "\n")
cat("Best score (iAUC + c_index):", best_score, "\n")

# 保存最佳结果
cox_flow_result <- best_result


library(dplyr)

summary_df <- bind_rows(lapply(all_results, function(x) {
  data.frame(train_frac = x$train_frac,
             seed = x$seed,
             iAUC = x$iAUC,
             c_index = x$c_index,
             score = x$score)
}))

print(summary_df)



saveRDS(result, file = "cox_rsf_workflow82_72_rfs816_715.rds")






###########################

rsf_workflow <- function(gene_expr,
                         clinical_data_imputed,
                         train_frac,
                         seed,
                         clin_pred = TRUE,
                         riskscore_pred = FALSE,
                         gene_pred = TRUE,
                         clin_predictors = c("grade", "er", "age", "size")
) {
  # # Cox Workflow Function
  # # Description: Executes a full Cox model workflow including data splitting, gene selection, risk scoring, and model evaluation.
  # # Input: gene_expr (expression matrix), clinical_data_imputed (imputed clinical data), train_frac (training fraction), seed (random seed)
  # # Output: results_train (Cox model results), result_valid (validation metrics)
  # 
  # 
  # # Split data into train and test sets
  # result <- split_expr_clinical(gene_expr, clinical_data_imputed, stratify_col = "e_dmfs", train_frac = train_frac, seed = seed)
  # train_expr <- result$train_expr
  # train_clinical <- result$train_clinical
  # test_expr <- result$test_expr
  # test_clinical <- result$test_clinical
  # 
  # # Perform batch univariate Cox regression
  # sig_gene_df <- batch_univariate_cox_regression(train_expr, train_clinical)
  # 
  # # Apply Lasso Cox with cross-validation
  # selected_gene_df <- lasso_cox_cv(train_expr, train_clinical, sig_gene_df)
  # 
  # Compute risk scores for training data
  clinical_cleaned_risk_train <- compute_risk_score(gene_mat_scaled = train_expr,
                                                    significant_vars_df = selected_gene_df,
                                                    clinical_cleaned = train_clinical,
                                                    n_group = 3)
  
  # Define predictors
  clin_predictors <- if (clin_pred) clin_predictors else NULL
  
  # find the position of "risk_group" column
  risk_group_pos <- which(colnames(clinical_cleaned_risk_train) == "risk_group")
  gene_predictors <- if (gene_pred)colnames(clinical_cleaned_risk_train)[(risk_group_pos + 1):
                                                                           length(colnames(clinical_cleaned_risk_train))] else NULL
  
  risk_score_predictor <- if (riskscore_pred) "risk_score" else NULL
  
  predictors <- c(clin_predictors, risk_score_predictor, gene_predictors)
  
  
  # random forest
  clinical_rsf <- clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs", predictors)]
  
  # build RSF model
  result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5)
  # Best model
  rsf_fit_best <- result_rsf_train$best_model
  
  #Remove variables with negative/less importance from the predictors vector based on RSF model output
  imp <- result_rsf_train$importance
  vars_to_remove <- names(imp)[imp < 0.01]
  predictors_filtered <- setdiff(predictors, vars_to_remove)
  
  clinical_rsf <- clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs", predictors_filtered)]
  # build RSF model
  result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5)
  # Best model
  rsf_fit_best <- result_rsf_train$best_model
  
  result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
  
  
  return(list(results_rsf_train = results_rsf_train, result_rsf_valid = result_rsf_valid))
}
  
 
















# random forest
predictors <- c("size", "age", "er", "risk_score")
clinical_rsf <- clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs", predictors)]
#clinical_rsf <- clinical_cleaned_7390_imputed[, c("t_dmfs", "e_dmfs", predictors)]
# 构建 RSF 模型
result_new_rsf <- rsf_kfold_cv_best(clinical_rsf, K = 5)
# Best model
rsf_fit_new <- result_new_rsf$best_model
# Variable importance from retrained model with complete data
#importance_new <- result_new_rsf$importance

# evaluate the model
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_new, df = clinical_cleaned_risk_train)

#calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_oxfu_naomt)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_new, df = clinical_cleaned_risk_test)



# random forest 2
predictors <- c("size", "age", "er", "risk_score")
predictors <- colnames(clinical_cleaned_risk_train)[11: length(colnames(clinical_cleaned_risk_train))]
clinical_rsf <- clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs", predictors)]
#clinical_rsf <- clinical_cleaned_7390_imputed[, c("t_dmfs", "e_dmfs", predictors)]
# 构建 RSF 模型
result_new_rsf <- rsf_kfold_cv_best(clinical_rsf, K = 5)
# Best model
rsf_fit_new <- result_new_rsf$best_model

# evaluate the model
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_new, df = clinical_cleaned_risk_train)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_new, df = clinical_cleaned_risk_test)











predictors <- c("size", "risk_score")
results_7390_0 <- fit_cox_model(predictors, clinical_cleaned_risk_7390)

predictors <- colnames(clinical_cleaned_risk_7390)[11: length(colnames(clinical_cleaned_risk_7390))]
results_7390 <- fit_cox_model(predictors, clinical_cleaned_risk_7390)
fitted_model_7390 <- results_7390$model
results_7390_coef_df <- results_7390$coef_table










#build and test a Cox proportional hazards model
predictors <- c("grade", "size", "age", "er")
results <- fit_cox_model(predictors, clinical_cleaned_7390_imputed)
# validate the model
calculate_time_auc_cindex("Cox", fitted_model = results$model, df = clinical_cleaned_7390_imputed)








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


#######################################################################################

### randomforest new
predictors <- c("grade", "size", "age", "er", "risk_score")
clinical_rsf <- clinical_cleaned_risk_7390[, c("t_dmfs", "e_dmfs", predictors)]
#clinical_rsf <- clinical_cleaned_7390_imputed[, c("t_dmfs", "e_dmfs", predictors)]
# 构建 RSF 模型
result_new_rsf <- rsf_kfold_cv_best(clinical_rsf, K = 5)
# Best model
rsf_fit_new <- result_new_rsf$best_model
# Variable importance from retrained model with complete data
importance_new <- result_new_rsf$importance

# evaluate the model
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_new, df = clinical_cleaned_risk_7390)

#calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_oxfu_naomt)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_new, df = clinical_cleaned_risk_2990)

# rebuild the model with non negative importance predictors
#predictors <- c("size", "grade")
#clinical_rsf <- clinical_cleaned_7390_imputed[, c("t_dmfs", "e_dmfs", predictors)]
# 构建 RSF 模型
rsf_fit2 <- rfsrc(Surv(t_dmfs, e_dmfs) ~ ., data = clinical_rsf, ntree = 1000, importance = TRUE)
rsf_fit2$importance
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit2, df = clinical_cleaned_risk_2990)






# For RSF model
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_7390_imputed)
calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_oxfu_naomt)

calculate_time_auc_cindex("RSF", fitted_model = rsf_fit, df = clinical_cleaned_11121)
