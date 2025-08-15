library(GEOquery)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(randomForestSRC)
library(pec)
library(mice)
library(caret)
library(riskRegression)
library(stringr)
options(scipen = 999)

## load dataset
#gset_7390 <- getGEO("GSE7390", GSEMatrix=TRUE)
#gset_39582 <- getGEO("GSE39582", GSEMatrix=TRUE)
expr_39582 <- exprs(gset_39582[[1]])
clinical_39582 <- pData(gset_39582[[1]])

col_selected_39582 <- c(
  "characteristics_ch1.2",
  "characteristics_ch1.3",
  "characteristics_ch1.4",
  "characteristics_ch1.5",
  "characteristics_ch1.6",
  "characteristics_ch1.7",
  "characteristics_ch1.8",
  "characteristics_ch1.9",
  "characteristics_ch1.11",
  "characteristics_ch1.12",
  "characteristics_ch1.13",
  "characteristics_ch1.14",
  "characteristics_ch1.15",
  "characteristics_ch1.16",
  "characteristics_ch1.17",
  "characteristics_ch1.18",
  "characteristics_ch1.22",
  "characteristics_ch1.26",
  "characteristics_ch1.30"
)



clinical_selected <- clinical_39582[col_selected_39582]





#创建一个函数来清理前缀

clean_prefix <- function(x, prefix) { str_replace(x, paste0("^", prefix, ": "), "") }

#数据清理

cleaned_data <- clinical_selected %>% mutate( 
  
  Sex = as.factor(clean_prefix(characteristics_ch1.2, "Sex")),
  # 清理诊断时年龄并转换为数值型
  Age_at_Diagnosis = as.numeric(clean_prefix(characteristics_ch1.3, "age.at.diagnosis \\(year\\)")),
  
  # 清理TNM分期并转换为有序因子
  TNM_Stage = factor(clean_prefix(characteristics_ch1.4, "tnm.stage"), 
                     levels = c("1", "2", "3", "4"), ordered = TRUE),
  
  # 清理T分期并转换为有序因子
  TNM_T = factor(clean_prefix(characteristics_ch1.5, "tnm.t"), 
                 levels = c("T1", "T2", "T3", "T4"), ordered = TRUE),
  
  # 清理N分期并转换为有序因子
  TNM_N = factor(clean_prefix(characteristics_ch1.6, "tnm.n"), 
                 levels = c("N0", "N1", "N2"), ordered = TRUE),
  
  # 清理M分期并转换为因子
  TNM_M = factor(clean_prefix(characteristics_ch1.7, "tnm.m"), 
                 levels = c("M0", "M1")),
  
  # 清理肿瘤位置并转换为因子
  Tumor_Location = factor(clean_prefix(characteristics_ch1.8, "tumor.location"), 
                          levels = c("proximal", "distal")),
  
  # 清理辅助化疗状态并转换为因子
  Chemotherapy_Adjuvant = factor(clean_prefix(characteristics_ch1.9, "chemotherapy.adjuvant"), 
                                 levels = c("N", "Y")),
  
  # # 清理无复发生存事件并转换为因子
  # RFS_Event = factor(clean_prefix(characteristics_ch1.11, "rfs.event"), 
  #                    levels = c("0", "1")),
  # 
  # # 清理无复发生存时间并转换为数值型
  # RFS_Delay = as.numeric(clean_prefix(characteristics_ch1.12, "rfs.delay")),
  
  # 清理总生存事件并转换为因子
  OS_Event = as.numeric(clean_prefix(characteristics_ch1.13, "os.event")), 
    # factor(clean_prefix(characteristics_ch1.13, "os.event"), 
    #                 levels = c("0", "1")),
    # 
  # 清理总生存时间并转换为数值型
  OS_Delay = as.numeric(clean_prefix(characteristics_ch1.14, "os.delay \\(months\\)")) * 30.42,
  
  # 清理MMR状态并转换为因子
  MMR_Status = factor(clean_prefix(characteristics_ch1.15, "mmr.status"), 
                      levels = c("pMMR", "dMMR")),
  
  # 清理CIMP状态并转换为因子
  CIMP_Status = factor(clean_prefix(characteristics_ch1.16, "cimp.status"), 
                       levels = c("-", "+")),
  
  # 清理CIN状态并转换为因子
  CIN_Status = factor(clean_prefix(characteristics_ch1.17, "cin.status"), 
                      levels = c("-", "+")),
  
  # 清理TP53突变状态并转换为因子
  TP53_Mutation = factor(clean_prefix(characteristics_ch1.18, "tp53.mutation"), 
                         levels = c("WT", "M")),
  
  # 清理KRAS突变状态并转换为因子
  KRAS_Mutation = factor(clean_prefix(characteristics_ch1.22, "kras.mutation"), 
                         levels = c("WT", "M")),
  
  # 清理BRAF突变状态并转换为因子
  BRAF_Mutation = factor(clean_prefix(characteristics_ch1.26, "braf.mutation"), 
                         levels = c("WT", "M")),
  
  # 清理分子亚型并转换为因子
  Molecular_Subtype = factor(clean_prefix(characteristics_ch1.30, "cit.molecularsubtype"), 
                             levels = c("C1", "C2", "C3", "C4", "C5", "C6"))) %>%
  # 选择清理后的列 #RFS_Event, RFS_Delay,
  select(Sex, Age_at_Diagnosis, TNM_T, TNM_N, TNM_M, Tumor_Location,
         Chemotherapy_Adjuvant,
         MMR_Status, KRAS_Mutation,
         OS_Event, OS_Delay) # Molecular_Subtype, TNM_Stage, # TP53_Mutation,214个NA BRAF_Mutation, 71N, CIN_Status, 102NA
#  CIMP_Status, NA70

# # %>%
# #   # 处理缺失值：数值型变量用中位数填充，分类变量保持不变
# #   mutate_if(is.numeric, ~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
# #   # 确保分类变量的因子水平一致
# #   mutate_if(is.factor, ~factor(., levels = levels(.)))
  
  
colSums(is.na(cleaned_data))

clinical_cleaned_39582 <- cleaned_data %>%
  filter(!is.na(OS_Delay),
         !is.na(OS_Event),
         rowSums(is.na(.)) < 6,
         OS_Delay > 0) 

colSums(is.na(clinical_cleaned_39582))
sum(is.na(clinical_cleaned_39582))



generate_box_plots(data=clinical_cleaned_39582, continuous_variables = c("Age_at_Diagnosis", "OS_Delay", "OS_Event")) 

# 查找包含 NA 的列名
na_columns <- colnames(clinical_cleaned_39582)[colSums(is.na(clinical_cleaned_39582)) > 0]
na_columns
selected_col <- colnames(clinical_cleaned_39582)[1:(ncol(clinical_cleaned_39582) - 2)]

clinical_cleaned_39582_imputed <- impute_missing_value(clinical_cleaned = clinical_cleaned_39582, selected_col, missed_col=na_columns)


col_mean <- mean(clinical_cleaned_39582_imputed$Age_at_Diagnosis, na.rm = TRUE)
clinical_cleaned_39582_imputed$Age_at_Diagnosis[is.na(clinical_cleaned_39582_imputed$Age_at_Diagnosis)] <- col_mean

#md.pattern(df_miss)


rownames_39582 <- rownames(clinical_cleaned_39582_imputed)
expr_39582_aligned <- expr_39582[, rownames_39582, drop = FALSE]


clinical_cleaned_39582_imputed <- clinical_cleaned_39582_imputed %>%
  rename(
    e_dmfs = OS_Event,
    t_dmfs = OS_Delay
  )


expr_mat=expr_39582_aligned

clinical_df = clinical_cleaned_39582_imputed %>%
  mutate(
    geo_accession = rownames(.)
  )


#############################################
# var_filter <- apply(expr_mat, 1, var)
# hist(var_filter, breaks = 100, main = "Gene Variance Distribution", xlab = "Variance")
# abline(v = quantile(var_filter, 0.5), col = "red", lwd = 2, lty = 2)
# low_var_genes <- expr_mat[var_filter <= quantile(var_filter, 0.25), ]
# summary(apply(low_var_genes, 1, sd))
# 
# # 2) thresholds
# q10 <- quantile(var_filter, 0.10, na.rm = TRUE)
# q25 <- quantile(var_filter, 0.25, na.rm = TRUE)
# 
# n_genes <- length(var_filter)
# n_drop10 <- sum(var_filter <= q10)
# n_drop25 <- sum(var_filter <= q25)
# 
# # 3) density plot + lines
# d <- density(var_filter, na.rm = TRUE)
# plot(d, main = sprintf("Gene variance density (q10=%.4f; drop=%d/%d, q25=%.4f; drop=%d/%d)",
#                        q10, n_drop10, n_genes, q25, n_drop25, n_genes),
#      xlab = "Per-gene variance", ylab = "Density")
# abline(v = q10, lty = 2)
# abline(v = q25, lty = 3)
# text(x = q10, y = max(d$y)*0.9, labels = "10% cutoff", srt = 90, pos = 4, cex = 0.8)
# text(x = q25, y = max(d$y)*0.8, labels = "25% cutoff", srt = 90, pos = 4, cex = 0.8)
# 
# 
# 
# expr_mat <- expr_mat[var_filter > quantile(var_filter, 0.25), ]
# 
# 








# 
# run_bootstrap_validation_safe <- function(expr_mat, clinical_df, 
#                                           B = 10, 
#                                           #max_vars = 30, 
#                                           seed = 90000,
#                                           min_epv = 2.5, 
#                                           coef_max = 10,
#                                           min_concord = 0.98) {
#   
#   n <- ncol(expr_mat)
#   perf_list <- list()
#   gene_list <- list()
#   rsf_predictor_list <- list()
#   
#   # 保存最优模型信息
#   best_model <- NULL
#   best_perf  <- -Inf
#   best_genes <- NULL
#   best_iter <- NULL
#   best_method <- NULL
#   train_indices_list <- list()
#   
#   for (b in 1:B) {
#     
#     # Bootstrap 抽样
#     current_seed <- seed + b 
#     set.seed(current_seed) 
#     message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
#     
#     train_idx <- sample(seq_len(n), size = n, replace = TRUE)
#     test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
#     
#     train_indices_list[[b]] <- train_idx
#     
#     if (length(test_idx) == 0) next
#     
#     train_expr     <- expr_mat[, train_idx]
#     test_expr      <- expr_mat[, test_idx]
#     train_clinical <- clinical_df[train_idx, ]
#     test_clinical  <- clinical_df[test_idx, ]
#     
#     # 事件数
#     events_train <- sum(train_clinical$e_dmfs)
#     message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
#     
#     # Median Absolute Deviation more robust than variance
#     mad_train_expr <- apply(train_expr, 1, mad)
#     cutoff <- quantile(mad_train_expr, 0.25)
#     #plot(density(mad_train_expr))
#     #abline(v = cutoff, col = "red", lty = 2)
#     
#     #hist(mad_train_expr, breaks = 50, main = "MAD Distribution",
#     #     xlab = "MAD", col = "lightblue", border = "white")
#     #abline(v = cutoff, col = "red", lwd = 2, lty = 2)
#     
#     keep_mad <- mad_train_expr > cutoff      # 过滤底部 20%（可调 10–30%）
#     # Filter the expression matrix
#     train_expr2 <- train_expr[keep_mad, ]
#     test_expr2 <- test_expr[keep_mad,]
#     
#     
#     
#     # Step 1: 单变量 Cox
#     sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = 0.01)
#     # sig_gene_df <- sig_gene_df[sig_gene_df$p.value < 0.05, ]
#     # if (nrow(sig_gene_df) == 0) next
#     
#     #sig_gene_df_39582 <- sig_gene_df
#     # 
#     # # 动态设置 p-value 阈值
#     # p_thresh <- if (events_train < 40) {
#     #   0.01
#     # } else {
#     #   0.05
#     # }
#     # 
#     # sig_gene_df <- sig_gene_df[sig_gene_df$p.value < p_thresh, ]
#     # if (nrow(sig_gene_df) == 0) {
#     #   message(sprintf("Bootstrap %d skipped: No genes passed p < %.2f (Events: %d)", b, p_thresh, events_train))
#     #   next
#     # }
#     # 
#     
#     
#     significant_gene = sig_gene_df$gene
#     
#     # scale train gene data
#     train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
#                                                 gene_mat_valid= train_expr2,
#                                                 significant_gene = significant_gene)
#     #dim(train_expr_scaled)
#     train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = 0.90)
#     #dim(train_expr_filtered)  # check new dimensions
#    
#     significant_gene2 = rownames(train_expr_filtered)
#     # scale test data
#     test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
#                                                gene_mat_valid= test_expr2,
#                                                significant_gene = significant_gene2)
#     
#     
#     # scale train clinical data
#     train_clinical_scaled <- standardize_with_train_clinical(train_clinical,
#                                                              train_clinical,
#                                                              scale_cols = c("Age_at_Diagnosis"))
#     
#     
#     # scale test data
#     test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
#                                                             test_clinical,
#                                                             scale_cols = c("Age_at_Diagnosis"))
#     
#     # Step 2: LASSO Cox
#     #selected_gene_df <- lasso_cox_cv(train_expr_filtered, train_clinical, sig_genes = significant_gene2)
#     
#     gene_freq_df <- repeat_cv_lasso_cox(train_expr = train_expr_filtered,
#                                          train_clinical,
#                                          significant_gene_vec= significant_gene2,
#                                          repeats = 5,
#                                          nfolds = 10,
#                                          alpha = 1)
#     
#     
#     #############################
#     # Step: 动态筛选变量，满足 EPV 要求
#     gene_freq_df_best <- gene_freq_df %>%
#       filter(freq >= 0.8)
#     
#     max_vars_allowed <- floor(events_train / min_epv)
#     
#     if (max_vars_allowed == 0 || nrow(gene_freq_df) == 0) {
#       message(sprintf("Bootstrap %d skipped: Too few events or no selected genes", b))
#       next
#     }
#     
#     # 限制变量数量（最多 max_vars_allowed 个）
#     
#     selected_gene_df <- gene_freq_df[1:min(nrow(gene_freq_df_best), max_vars_allowed), ] %>%
#       mutate(coef = mean_coef)
#     
#     # 打印信息
#     message(sprintf("Bootstrap %d: Events = %d, Vars = %d, EPV = %.2f", 
#                     b, events_train, nrow(selected_gene_df), events_train / nrow(selected_gene_df)))
#     
#     
#     
#     # Step 3: 风险分数
#     clinical_cleaned_risk_train <- compute_risk_score(
#       gene_mat_scaled = train_expr_filtered,
#       significant_vars_df = selected_gene_df,
#       clinical_cleaned = train_clinical_scaled,
#       n_group = 3
#     )
#     
#     predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
#                      "Chemotherapy_Adjuvant", "MMR_Status",
#                      "KRAS_Mutation")
#     #    "grade", "er", "age", "size"
#     
#     predictors  <- c(predictors0, colnames(clinical_cleaned_risk_train)[15:ncol(clinical_cleaned_risk_train)])
#     
#     # 相关性过滤
#     predictor_data <- clinical_cleaned_risk_train[, predictors]
#     numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
#     
#     filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.9)
#     
#     filtered_data <- cbind(
#       predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
#       filtered_numeric_vars
#     )
#     
#     predictors_filtered <- colnames(filtered_data)
#     df <- cbind(filtered_data, clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs")])
#     
#     # 拟合模型
#     results_train <- fit_cox_model(predictors_filtered, df)
#     
#     if (is.null(results_train$model)) {
#       message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
#       next
#     }
#     
#     
#     # 发散检测
#     if (any(abs(coef(results_train$model)) > coef_max)) {
#       message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
#       next
#     }
#     
#     # 检查完全分离
#     concordance_val <- tryCatch({
#       suppressWarnings(summary(results_train$model)$concordance[1])
#     }, error = function(e) NA)
#     if (!is.na(concordance_val) && concordance_val >= min_concord) {
#       message(
#         sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
#       next
#     }
#     
#     # Step 4: 测试集评估
#     clinical_cleaned_risk_test <- compute_risk_score(
#       gene_mat_scaled = test_expr_scaled,
#       significant_vars_df = selected_gene_df,
#       clinical_cleaned = test_clinical_scaled,
#       n_group = 3
#     )
#     
#     result_valid <- calculate_time_auc_cindex(
#       "Cox", 
#       fitted_model = results_train$model, 
#       df = clinical_cleaned_risk_test
#     )
#     
#     # perf_list[[b]] <- result_valid
#     # gene_list[[b]] <- selected_gene_df$gene
#     # 
#     # # 保存最优模型
#     # if (!is.na(result_valid$iAUC) && result_valid$iAUC > best_perf) {
#     #   best_perf  <- result_valid$iAUC
#     #   best_model <- results_train$model
#     #   best_genes <- selected_gene_df$gene
#     # }
#     
#     #
#     # random forest
#     clinical_rsf <- df
#     
#     # build RSF model
#     result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000) #, seed = current_seed
#     # Best model
#     rsf_fit_best <- result_rsf_train$best_model
#     
#     #Remove variables with negative/less importance from the predictors vector based on RSF model output
#     # imp <- result_rsf_train$importance
#     # vars_to_remove <- names(imp)[imp < 0.01]
#     # predictors_filtered_rsf <- setdiff(predictors_filtered, vars_to_remove)
#     # 
#     # clinical_rsf <- df[, c("t_dmfs", "e_dmfs", predictors_filtered_rsf)]
#     # # build RSF model with predictors_filtered_rsf
#     # result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5, ntree = 1000) #, seed = current_seed
#     # # Best model
#     # rsf_fit_best <- result_rsf_train$best_model
#     result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
#     
#     
#     # 把 Cox + RSF 一起保存
#     perf_list[[b]] <- list(
#       cox  = result_valid,
#       rsf  = result_rsf_valid
#     )
#     
#     # 保存基因
#     gene_list[[b]] <- selected_gene_df$gene
#     rsf_predictor_list[[b]] <- predictors_filtered#predictors_filtered_rsf
#     
#     
#     # 当前 Cox 和 RSF 的综合评分（iAUC + C-index）
#     cox_score <- if (!is.na(result_valid$iAUC) && !is.na(result_valid$c_index)) {
#       result_valid$iAUC + result_valid$c_index
#     } else {
#       -Inf
#     }
#     
#     rsf_score <- if (!is.na(result_rsf_valid$iAUC) && !is.na(result_rsf_valid$c_index)) {
#       result_rsf_valid$iAUC + result_rsf_valid$c_index
#     } else {
#       -Inf
#     }
#     cat(sprintf("Cox_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f",
#                 b,
#                 cox_score,
#                 result_valid$iAUC,
#                 result_valid$c_index
#                 ))
#     cat(sprintf("rsf_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f",
#                 b,
#                 rsf_score,
#                 result_rsf_valid$iAUC,
#                 result_rsf_valid$c_index))
#     
#     
#     # 判断是否更优（Cox）
#     if (cox_score > best_perf) {
#       best_perf    <- cox_score
#       best_model   <- results_train$model
#       best_genes   <- selected_gene_df$gene
#       best_iter    <- b
#       best_seed    <- current_seed
#       best_method  <- "Cox"
#     }
#     
#     # 判断是否更优（RSF）
#     if (rsf_score > best_perf) {
#       best_perf    <- rsf_score
#       best_model   <- rsf_fit_best
#       best_genes   <- predictors_filtered#predictors_filtered_rsf
#       best_iter    <- b
#       best_seed    <- current_seed
#       best_method  <- "RSF"
#     }
#     
#     
#     cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
#     cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
#   }
#   
#   # 汇总
#  # iAUCs   <- sapply(perf_list, function(x) x$iAUC)
#   #cindexs <- sapply(perf_list, function(x) x$c_index)
#   
#   
#   
#   ######################
#   # 过滤掉为 NULL 的轮次
#   perf_list_nz <- Filter(Negate(is.null), perf_list)
#   
#   # 取 Cox 指标
#   cox_iAUC  <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$cox) && !is.null(x$cox$iAUC)) x$cox$iAUC else NA_real_)
#   cox_cidx  <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$cox) && !is.null(x$cox$c_index)) x$cox$c_index else NA_real_)
#   
#   # 取 RSF 指标
#   rsf_iAUC  <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$rsf) && !is.null(x$rsf$iAUC)) x$rsf$iAUC else NA_real_)
#   rsf_cidx  <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$rsf) && !is.null(x$rsf$c_index)) x$rsf$c_index else NA_real_)
#   
#   # 汇总（去掉 NA）
#   mean_cox_iAUC <- mean(cox_iAUC, na.rm = TRUE)
#   mean_cox_cidx <- mean(cox_cidx, na.rm = TRUE)
#   mean_rsf_iAUC <- mean(rsf_iAUC, na.rm = TRUE)
#   mean_rsf_cidx <- mean(rsf_cidx, na.rm = TRUE)
#   
#   
#   ######################################################
#   
#   all_genes <- unlist(gene_list)
#   gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
#   
#   return (list(
#     #mean_iAUC   = mean(iAUCs, na.rm = TRUE),
#     #mean_cindex = mean(cindexs, na.rm = TRUE),
#     mean_cox_iAUC = mean_cox_iAUC,
#     mean_cox_cidx = mean_cox_cidx,
#     mean_rsf_cidx = mean_rsf_cidx,
#     mean_rsf_iAUC = mean_rsf_iAUC,
#     gene_frequency = gene_freq,
#     all_results = perf_list,
#     best_model  = best_model,   # ⬅ 最优模型对象
#     best_perf   = best_perf,    # ⬅ 最优 iAUC
#     best_predictors = best_genes,   # ⬅ 最优模型的基因列表或者RSF预测变量
#     best_iter = best_iter,
#     best_seed = best_seed,
#     best_method = best_method,
#     train_indices_list = train_indices_list
#   ))
# }
# 
# 
# res39582_g <- run_bootstrap_validation_safe(expr_mat,
#                                       clinical_df, 
#                                       B = 5,
#                                       seed = 90000,
#                                       # max_vars = 20,
#                                       min_epv = 4,
#                                       coef_max = 10,
#                                       min_concord = 0.97) 
# 




run_bootstrap_validation_safe <- function(expr_mat, clinical_df, 
                                          B = 10, 
                                          seed = 90000,
                                          min_epv = 2.5, 
                                          coef_max = 10,
                                          min_concord = 0.98) {
  
  n <- ncol(expr_mat)
  perf_list <- list()
  gene_list <- list()
  rsf_predictor_list <- list()
  
  # 保存最优模型信息
  best_model <- NULL
  best_perf  <- -Inf
  best_genes <- NULL
  best_iter <- NULL
  best_method <- NULL
  train_indices_list <- list()
  
  message("初始化完成：expr_mat 维度 = ", toString(dim(expr_mat)), 
          ", clinical_df 维度 = ", toString(dim(clinical_df)))
  
  for (b in 1:B) {
    # Bootstrap 采样
    current_seed <- seed + b 
    set.seed(current_seed) 
    message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
    
    train_idx <- sample(seq_len(n), size = n, replace = TRUE)
    test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
    message(sprintf("Bootstrap %d: train_idx length = %d, test_idx length = %d", 
                    b, length(train_idx), length(test_idx)))
    
    if (length(test_idx) == 0) {
      message(sprintf("Bootstrap %d skipped: No OOB samples", b))
      next
    }
    
    train_expr     <- expr_mat[, train_idx, drop = FALSE]
    test_expr      <- expr_mat[, test_idx, drop = FALSE]
    train_clinical <- clinical_df[train_idx, , drop = FALSE]
    test_clinical  <- clinical_df[test_idx, , drop = FALSE]
    
    message("Bootstrap 采样完成：train_expr 维度 = ", toString(dim(train_expr)), 
            ", test_expr 维度 = ", toString(dim(test_expr)),
            ", train_clinical 维度 = ", toString(dim(train_clinical)),
            ", test_clinical 维度 = ", toString(dim(test_clinical)))
    
    if (ncol(train_expr) == 0 || nrow(train_expr) == 0) {
      message(sprintf("Bootstrap %d skipped: Invalid train_expr dimensions", b))
      next
    }
    
    # 事件数
    events_train <- sum(train_clinical$e_dmfs)
    message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
    message("事件数计算完成")
    
    # MAD 过滤
    mad_train_expr <- apply(train_expr, 1, mad)
    cutoff <- quantile(mad_train_expr, 0.25)
    keep_mad <- mad_train_expr > cutoff
    message(sprintf("Bootstrap %d: MAD 过滤完成，保留基因数 = %d, cutoff = %.4f", 
                    b, sum(keep_mad), cutoff))
    
    if (sum(keep_mad) == 0) {
      message(sprintf("Bootstrap %d skipped: No genes passed MAD filter", b))
      next
    }
    
    train_expr2 <- train_expr[keep_mad, , drop = FALSE]
    test_expr2 <- test_expr[keep_mad, , drop = FALSE]
    message("MAD 过滤后矩阵生成完成：train_expr2 维度 = ", toString(dim(train_expr2)), 
            ", test_expr2 维度 = ", toString(dim(test_expr2)))
    
    # 单变量 Cox 回归
    sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = 0.01)
    message(sprintf("Bootstrap %d: 单变量 Cox 回归完成，显著基因数 = %d", 
                    b, if (is.null(sig_gene_df)) 0 else nrow(sig_gene_df)))
    
    if (is.null(sig_gene_df) || nrow(sig_gene_df) == 0) {
      message(sprintf("Bootstrap %d skipped: No significant genes from univariate Cox", b))
      next
    }
    
    significant_gene <- sig_gene_df$gene
    message(sprintf("Bootstrap %d: 提取显著基因完成，significant_gene 长度 = %d", 
                    b, length(significant_gene)))
    
    if (length(significant_gene) < 2) {
      message(sprintf("Bootstrap %d skipped: significant genes < 2", b))
      next
    }
    
    # 标准化训练集基因数据
    train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                                gene_mat_valid = train_expr2,
                                                significant_gene = significant_gene)
    message("标准化训练集基因数据完成：train_expr_scaled 维度 = ", 
            toString(dim(train_expr_scaled)))
    
    if (is.null(train_expr_scaled) || nrow(train_expr_scaled) == 0) {
      message(sprintf("Bootstrap %d skipped: Invalid standardized training data", b))
      next
    }
    
    if (nrow(train_expr_scaled) < 2) {
      message(sprintf("Bootstrap %d skipped: nrow(train_expr_scaled) < 2", b))
      next
    }
    
    # 相关性过滤
    train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = 0.90)
    message("相关性过滤完成：train_expr_filtered 维度 = ", 
            toString(dim(train_expr_filtered)))
    
    if (is.null(train_expr_filtered) || nrow(train_expr_filtered) < 2) {
      message(sprintf("Bootstrap %d skipped: No genes after correlation filtering", b))
      next
    }
    
    significant_gene2 <- rownames(train_expr_filtered)
    message(sprintf("Bootstrap %d: 提取过滤后基因完成，significant_gene2 长度 = %d", 
                    b, length(significant_gene2)))
    
    # 标准化测试集基因数据
    test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                               gene_mat_valid = test_expr2,
                                               significant_gene = significant_gene2)
    message("标准化测试集基因数据完成：test_expr_scaled 维度 = ", 
            toString(dim(test_expr_scaled)))
    
    # 标准化训练集临床数据
    train_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                             train_clinical,
                                                             scale_cols = c("Age_at_Diagnosis"))
    message("标准化训练集临床数据完成：train_clinical_scaled 维度 = ", 
            toString(dim(train_clinical_scaled)))
    
    # 标准化测试集临床数据
    test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                            test_clinical,
                                                            scale_cols = c("Age_at_Diagnosis"))
    message("标准化测试集临床数据完成：test_clinical_scaled 维度 = ", 
            toString(dim(test_clinical_scaled)))
    
    # LASSO Cox
    gene_freq_df <- repeat_cv_lasso_cox(train_expr = train_expr_filtered,
                                        train_clinical,
                                        significant_gene_vec = significant_gene2,
                                        repeats = 5,
                                        nfolds = 10,
                                        alpha = 1)
    message(sprintf("Bootstrap %d: LASSO Cox 完成，gene_freq_df 行数 = %d", 
                    b, nrow(gene_freq_df)))
    
    # 动态筛选变量，满足 EPV 要求
    gene_freq_df_best <- gene_freq_df %>%
      filter(freq >= 0.8)
    max_vars_allowed <- floor(events_train / min_epv)
    message(sprintf("Bootstrap %d: 动态筛选变量完成，gene_freq_df_best 行数 = %d, max_vars_allowed = %d", 
                    b, nrow(gene_freq_df_best), max_vars_allowed))
    
    if (max_vars_allowed == 0 || nrow(gene_freq_df) == 0) {
      message(sprintf("Bootstrap %d skipped: Too few events or no selected genes", b))
      next
    }
    
    selected_gene_df <- gene_freq_df[1:min(nrow(gene_freq_df_best), max_vars_allowed), ] %>%
      mutate(coef = mean_coef)
    message(sprintf("Bootstrap %d: 选择基因完成，selected_gene_df 行数 = %d, EPV = %.2f", 
                    b, nrow(selected_gene_df), events_train / nrow(selected_gene_df)))
    
    # 计算风险分数
    clinical_cleaned_risk_train <- compute_risk_score(
      gene_mat_scaled = train_expr_filtered,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = train_clinical_scaled,
      n_group = 3
    )
    message("训练集风险分数计算完成：clinical_cleaned_risk_train 维度 = ", 
            toString(dim(clinical_cleaned_risk_train)))
    
    predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
                     "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")
    predictors <- c(predictors0, colnames(clinical_cleaned_risk_train)[15:ncol(clinical_cleaned_risk_train)])
    message(sprintf("Bootstrap %d: 预测变量选择完成，predictors 长度 = %d", 
                    b, length(predictors)))
    
    # 相关性过滤
    predictor_data <- clinical_cleaned_risk_train[, predictors]
    numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
    filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.9)
    filtered_data <- cbind(
      predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
      filtered_numeric_vars
    )
    predictors_filtered <- colnames(filtered_data)
    df <- cbind(filtered_data, clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs")])
    message("预测变量相关性过滤完成：filtered_data 维度 = ", toString(dim(filtered_data)))
    
    # 拟合 Cox 模型
    results_train <- fit_cox_model(predictors_filtered, df)
    message(sprintf("Bootstrap %d: Cox 模型拟合完成，model 是否存在 = %s", 
                    b, !is.null(results_train$model)))
    
    if (is.null(results_train$model)) {
      message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
      next
    }
    
    # 发散检测
    if (any(abs(coef(results_train$model)) > coef_max)) {
      message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
      next
    }
    message("Cox 模型发散检测完成")
    
    # 检查完全分离
    concordance_val <- tryCatch({
      suppressWarnings(summary(results_train$model)$concordance[1])
    }, error = function(e) NA)
    message(sprintf("Bootstrap %d: 完全分离检查完成，concordance_val = %.4f", 
                    b, if (is.na(concordance_val)) NA else concordance_val))
    
    if (!is.na(concordance_val) && concordance_val >= min_concord) {
      message(sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
      next
    }
    
    # 测试集风险分数
    clinical_cleaned_risk_test <- compute_risk_score(
      gene_mat_scaled = test_expr_scaled,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = test_clinical_scaled,
      n_group = 3
    )
    message("测试集风险分数计算完成：clinical_cleaned_risk_test 维度 = ", 
            toString(dim(clinical_cleaned_risk_test)))
    
    # Cox 模型评估
    result_valid <- calculate_time_auc_cindex(
      "Cox", 
      fitted_model = results_train$model, 
      df = clinical_cleaned_risk_test
    )
    message(sprintf("Bootstrap %d: Cox 模型评估完成，iAUC = %.4f, c_index = %.4f", 
                    b, result_valid$iAUC, result_valid$c_index))
    
    # Random Survival Forest
    clinical_rsf <- df
    result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000)
    rsf_fit_best <- result_rsf_train$best_model
    message("RSF 模型拟合完成")
    
    result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
    message(sprintf("Bootstrap %d: RSF 模型评估完成，iAUC = %.4f, c_index = %.4f", 
                    b, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    
    # 保存性能和基因
    perf_list[[b]] <- list(
      cox = result_valid,
      rsf = result_rsf_valid
    )
    gene_list[[b]] <- selected_gene_df$gene
    rsf_predictor_list[[b]] <- predictors_filtered
    message("性能和基因保存完成")
    
    # 计算 Cox 和 RSF 综合评分
    cox_score <- if (!is.na(result_valid$iAUC) && !is.na(result_valid$c_index)) {
      result_valid$iAUC + result_valid$c_index
    } else {
      -Inf
    }
    rsf_score <- if (!is.na(result_rsf_valid$iAUC) && !is.na(result_rsf_valid$c_index)) {
      result_rsf_valid$iAUC + result_rsf_valid$c_index
    } else {
      -Inf
    }
    cat(sprintf("Cox_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
                b, cox_score, result_valid$iAUC, result_valid$c_index))
    cat(sprintf("rsf_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
                b, rsf_score, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    
    # 判断是否更优（Cox）
    if (cox_score > best_perf) {
      best_perf    <- cox_score
      best_model   <- results_train$model
      best_genes   <- selected_gene_df$gene
      best_iter    <- b
      best_seed    <- current_seed
      best_method  <- "Cox"
    }
    
    # 判断是否更优（RSF）
    if (rsf_score > best_perf) {
      best_perf    <- rsf_score
      best_model   <- rsf_fit_best
      best_genes   <- predictors_filtered
      best_iter    <- b
      best_seed    <- current_seed
      best_method  <- "RSF"
    }
    
    cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
    cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
    message(sprintf("Bootstrap %d: 迭代完成", b))
  }
  
  message("所有 Bootstrap 迭代完成")
  
  # 汇总
  perf_list_nz <- Filter(Negate(is.null), perf_list)
  
  cox_iAUC <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$iAUC)) x$cox$iAUC else NA_real_)
  cox_cidx <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$c_index)) x$cox$c_index else NA_real_)
  
  rsf_iAUC <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$iAUC)) x$rsf$iAUC else NA_real_)
  rsf_cidx <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$c_index)) x$rsf$c_index else NA_real_)
  
  mean_cox_iAUC <- mean(cox_iAUC, na.rm = TRUE)
  mean_cox_cidx <- mean(cox_cidx, na.rm = TRUE)
  mean_rsf_iAUC <- mean(rsf_iAUC, na.rm = TRUE)
  mean_rsf_cidx <- mean(rsf_cidx, na.rm = TRUE)
  
  message("性能指标汇总完成：mean_cox_iAUC = ", mean_cox_iAUC, 
          ", mean_cox_cidx = ", mean_cox_cidx,
          ", mean_rsf_iAUC = ", mean_rsf_iAUC,
          ", mean_rsf_cidx = ", mean_rsf_cidx)
  
  all_genes <- unlist(gene_list)
  gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  message("基因频率计算完成")
  
  return(list(
    mean_cox_iAUC = mean_cox_iAUC,
    mean_cox_cidx = mean_cox_cidx,
    mean_rsf_cidx = mean_rsf_cidx,
    mean_rsf_iAUC = mean_rsf_iAUC,
    gene_frequency = gene_freq,
    all_results = perf_list,
    best_model = best_model,
    best_perf = best_perf,
    best_predictors = best_genes,
    best_iter = best_iter,
    best_seed = best_seed,
    best_method = best_method,
    train_indices_list = train_indices_list
  ))
}


res39582_h <- run_bootstrap_validation_safe(expr_mat, clinical_df, 
                                            B = 20, 
                                            seed = 100000, 
                                            min_epv = 4, 
                                            coef_max = 10, 
                                            min_concord = 0.97)


saveRDS(res39582_h, file = "train_result_39582_h_best1_4702.rds")


res39582$best_model   # 最优 Cox 模型对象
res39582$best_perf    # 最优 iAUC
res39582$best_genes   # 对应的基因列表
res39582$mean_iAUC

res39582$mean_cindex
res39582$gene_frequency

res39582$best_predictors
res39582$best_iter
res39582_c$

  
  
  
  
  








##################################################################################
clinical_run_bootstrap_validation <- function(clinical_df, 
                                          B = 10, 
                                          seed = 90000,
                                          #min_epv = 2.5, 
                                          coef_max = 10,
                                          min_concord = 0.98) {
  
  n <- nrow(clinical_df)
  perf_list <- list()
  #gene_list <- list()
  #rsf_predictor_list <- list()
  
  # 保存最优模型信息
  best_model <- NULL
  best_perf  <- -Inf
  best_genes <- NULL
  best_iter <- NULL
  best_method <- NULL
  train_indices_list <- list()
  
  # message("初始化完成：expr_mat 维度 = ", toString(dim(expr_mat)), 
  #         ", clinical_df 维度 = ", toString(dim(clinical_df)))
  # 
  for (b in 1:B) {
    # Bootstrap 采样
    current_seed <- seed + b 
    set.seed(current_seed) 
    message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
    
    train_idx <- sample(seq_len(n), size = n, replace = TRUE)
    test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
    message(sprintf("Bootstrap %d: train_idx length = %d, test_idx length = %d", 
                    b, length(train_idx), length(test_idx)))
    
    if (length(test_idx) == 0) {
      message(sprintf("Bootstrap %d skipped: No OOB samples", b))
      next
    }
    
    
    train_clinical <- clinical_df[train_idx, , drop = FALSE]
    test_clinical  <- clinical_df[test_idx, , drop = FALSE]
    
    message("Bootstrap 采样完成：", #train_expr 维度 = ", toString(dim(train_expr)), 
            #", test_expr 维度 = ", toString(dim(test_expr)),
            ", train_clinical 维度 = ", toString(dim(train_clinical)),
            ", test_clinical 维度 = ", toString(dim(test_clinical)))
    
    
    
    # 事件数
    events_train <- sum(train_clinical$e_dmfs)
    message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
    message("事件数计算完成")
    
     
    predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
                     "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")
    
    message(sprintf("Bootstrap %d: 预测变量选择完成，predictors 长度 = %d", 
                    b, length(predictors0)))
    
    # 相关性过滤
    # predictor_data <- clinical_cleaned_risk_train[, predictors]
    # numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
    # filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.9)
    # filtered_data <- cbind(
    #   predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
    #   filtered_numeric_vars
    # )
    # predictors_filtered <- colnames(filtered_data)
    #df <- cbind(filtered_data, clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs")])
    df <- train_clinical %>%
      select(-geo_accession )
   # message("预测变量相关性过滤完成：filtered_data 维度 = ", toString(dim(filtered_data)))
    
    # 拟合 Cox 模型
    results_train <- fit_cox_model(predictors0, df)
    #vif_result <- vif(results_train$model)
    
    message(sprintf("Bootstrap %d: Cox 模型拟合完成，model 是否存在 = %s", 
                    b, !is.null(results_train$model)))
    
    if (is.null(results_train$model)) {
      message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
      next
    }
    
    # 发散检测
    if (any(abs(coef(results_train$model)) > coef_max)) {
      message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
      next
    }
    message("Cox 模型发散检测完成")
    
    # 检查完全分离
    concordance_val <- tryCatch({
      suppressWarnings(summary(results_train$model)$concordance[1])
    }, error = function(e) NA)
    message(sprintf("Bootstrap %d: 完全分离检查完成，concordance_val = %.4f", 
                    b, if (is.na(concordance_val)) NA else concordance_val))
    
    if (!is.na(concordance_val) && concordance_val >= min_concord) {
      message(sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
      next
    }
    
    
    # Cox 模型评估
    result_valid <- calculate_time_auc_cindex(
      "Cox", 
      fitted_model = results_train$model, 
      df = test_clinical
    )
    message(sprintf("Bootstrap %d: Cox 模型评估完成，iAUC = %.4f, c_index = %.4f", 
                    b, result_valid$iAUC, result_valid$c_index))
    
    # Random Survival Forest
    clinical_rsf <- df
    result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000)
    rsf_fit_best <- result_rsf_train$best_model
    message("RSF 模型拟合完成")
    
    result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = test_clinical)
    message(sprintf("Bootstrap %d: RSF 模型评估完成，iAUC = %.4f, c_index = %.4f", 
                    b, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    
    # 保存性能和基因
    perf_list[[b]] <- list(
      cox = result_valid,
      rsf = result_rsf_valid
    )
    #gene_list[[b]] <- selected_gene_df$gene
    #rsf_predictor_list[[b]] <- predictors_filtered
    #message("性能和基因保存完成")
    
    # 计算 Cox 和 RSF 综合评分
    cox_score <- if (!is.na(result_valid$iAUC) && !is.na(result_valid$c_index)) {
      result_valid$iAUC + result_valid$c_index
    } else {
      -Inf
    }
    rsf_score <- if (!is.na(result_rsf_valid$iAUC) && !is.na(result_rsf_valid$c_index)) {
      result_rsf_valid$iAUC + result_rsf_valid$c_index
    } else {
      -Inf
    }
    cat(sprintf("Cox_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
                b, cox_score, result_valid$iAUC, result_valid$c_index))
    cat(sprintf("rsf_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
                b, rsf_score, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    
    # 判断是否更优（Cox）
    if (cox_score > best_perf) {
      best_perf    <- cox_score
      best_model   <- results_train$model
      #best_genes   <- selected_gene_df$gene
      best_iter    <- b
      best_seed    <- current_seed
      best_method  <- "Cox"
    }
    
    # 判断是否更优（RSF）
    if (rsf_score > best_perf) {
      best_perf    <- rsf_score
      best_model   <- rsf_fit_best
      #best_genes   <- predictors_filtered
      best_iter    <- b
      best_seed    <- current_seed
      best_method  <- "RSF"
    }
    
    cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
    cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
    message(sprintf("Bootstrap %d: 迭代完成", b))
  }
  
  message("所有 Bootstrap 迭代完成")
  
  # 汇总
  perf_list_nz <- Filter(Negate(is.null), perf_list)
  
  cox_iAUC <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$iAUC)) x$cox$iAUC else NA_real_)
  cox_cidx <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$c_index)) x$cox$c_index else NA_real_)
  
  rsf_iAUC <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$iAUC)) x$rsf$iAUC else NA_real_)
  rsf_cidx <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$c_index)) x$rsf$c_index else NA_real_)
  
  mean_cox_iAUC <- mean(cox_iAUC, na.rm = TRUE)
  mean_cox_cidx <- mean(cox_cidx, na.rm = TRUE)
  mean_rsf_iAUC <- mean(rsf_iAUC, na.rm = TRUE)
  mean_rsf_cidx <- mean(rsf_cidx, na.rm = TRUE)
  
  message("性能指标汇总完成：mean_cox_iAUC = ", mean_cox_iAUC, 
          ", mean_cox_cidx = ", mean_cox_cidx,
          ", mean_rsf_iAUC = ", mean_rsf_iAUC,
          ", mean_rsf_cidx = ", mean_rsf_cidx)
  
  #all_genes <- unlist(gene_list)
  #gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  #message("基因频率计算完成")
  
  return(list(
    mean_cox_iAUC = mean_cox_iAUC,
    mean_cox_cidx = mean_cox_cidx,
    mean_rsf_cidx = mean_rsf_cidx,
    mean_rsf_iAUC = mean_rsf_iAUC,
    #gene_frequency = gene_freq,
    all_results = perf_list,
    best_model = best_model,
    best_perf = best_perf,
    #best_predictors = best_genes,
    best_iter = best_iter,
    best_seed = best_seed,
    best_method = best_method,
    train_indices_list = train_indices_list
  ))
}



res39582_clinical <- clinical_run_bootstrap_validation(clinical_df, 
                                            B = 20, 
                                            seed = 2000, 
                                           # min_epv = 4, 
                                            coef_max = 10, 
                                            min_concord = 0.97)

saveRDS(res39582_clinical, file = "train_result_39582_clinical_best1_5739.rds")
#################################################################################

gen_run_bootstrap_validation_safe <- function(expr_mat, clinical_df, 
                                          B = 10, 
                                          seed = 90000,
                                          min_epv = 2.5, 
                                          coef_max = 10,
                                          min_concord = 0.98) {
  
  n <- ncol(expr_mat)
  perf_list <- list()
  gene_list <- list()
  rsf_predictor_list <- list()
  
  # 保存最优模型信息
  best_model <- NULL
  best_perf  <- -Inf
  best_genes <- NULL
  best_iter <- NULL
  best_method <- NULL
  train_indices_list <- list()
  
  clinical_df <- clinical_df %>%
    select( "e_dmfs", "t_dmfs", "geo_accession")
  
  message("初始化完成：expr_mat 维度 = ", toString(dim(expr_mat)), 
          ", clinical_df 维度 = ", toString(dim(clinical_df)))
  
  for (b in 1:B) {
    # Bootstrap 采样
    current_seed <- seed + b 
    set.seed(current_seed) 
    message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
    
    train_idx <- sample(seq_len(n), size = n, replace = TRUE)
    test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
    message(sprintf("Bootstrap %d: train_idx length = %d, test_idx length = %d", 
                    b, length(train_idx), length(test_idx)))
    
    if (length(test_idx) == 0) {
      message(sprintf("Bootstrap %d skipped: No OOB samples", b))
      next
    }
    
    train_expr     <- expr_mat[, train_idx, drop = FALSE]
    test_expr      <- expr_mat[, test_idx, drop = FALSE]
    train_clinical <- clinical_df[train_idx, , drop = FALSE]
    test_clinical  <- clinical_df[test_idx, , drop = FALSE]
    
    message("Bootstrap 采样完成：train_expr 维度 = ", toString(dim(train_expr)), 
            ", test_expr 维度 = ", toString(dim(test_expr)),
            ", train_clinical 维度 = ", toString(dim(train_clinical)),
            ", test_clinical 维度 = ", toString(dim(test_clinical)))
    
    if (ncol(train_expr) == 0 || nrow(train_expr) == 0) {
      message(sprintf("Bootstrap %d skipped: Invalid train_expr dimensions", b))
      next
    }
    
    # 事件数
    events_train <- sum(train_clinical$e_dmfs)
    message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
    message("事件数计算完成")
    
    # MAD 过滤
    mad_train_expr <- apply(train_expr, 1, mad)
    cutoff <- quantile(mad_train_expr, 0.25)
    keep_mad <- mad_train_expr > cutoff
    message(sprintf("Bootstrap %d: MAD 过滤完成，保留基因数 = %d, cutoff = %.4f", 
                    b, sum(keep_mad), cutoff))
    
    if (sum(keep_mad) == 0) {
      message(sprintf("Bootstrap %d skipped: No genes passed MAD filter", b))
      next
    }
    
    train_expr2 <- train_expr[keep_mad, , drop = FALSE]
    test_expr2 <- test_expr[keep_mad, , drop = FALSE]
    message("MAD 过滤后矩阵生成完成：train_expr2 维度 = ", toString(dim(train_expr2)), 
            ", test_expr2 维度 = ", toString(dim(test_expr2)))
    
    # 单变量 Cox 回归
    sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = 0.01)
    message(sprintf("Bootstrap %d: 单变量 Cox 回归完成，显著基因数 = %d", 
                    b, if (is.null(sig_gene_df)) 0 else nrow(sig_gene_df)))
    
    if (is.null(sig_gene_df) || nrow(sig_gene_df) == 0) {
      message(sprintf("Bootstrap %d skipped: No significant genes from univariate Cox", b))
      next
    }
    
    significant_gene <- sig_gene_df$gene
    message(sprintf("Bootstrap %d: 提取显著基因完成，significant_gene 长度 = %d", 
                    b, length(significant_gene)))
    
    
    # 标准化训练集基因数据
    train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                                gene_mat_valid = train_expr2,
                                                significant_gene = significant_gene)
    message("标准化训练集基因数据完成：train_expr_scaled 维度 = ", 
            toString(dim(train_expr_scaled)))
    
    if (is.null(train_expr_scaled) || nrow(train_expr_scaled) == 0) {
      message(sprintf("Bootstrap %d skipped: Invalid standardized training data", b))
      next
    }
    
    
    # 相关性过滤
    train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = 0.90)
    message("相关性过滤完成：train_expr_filtered 维度 = ", 
            toString(dim(train_expr_filtered)))
    
    
    significant_gene2 <- rownames(train_expr_filtered)
    message(sprintf("Bootstrap %d: 提取过滤后基因完成，significant_gene2 长度 = %d", 
                    b, length(significant_gene2)))
    
    # 标准化测试集基因数据
    test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                               gene_mat_valid = test_expr2,
                                               significant_gene = significant_gene2)
    message("标准化测试集基因数据完成：test_expr_scaled 维度 = ", 
            toString(dim(test_expr_scaled)))
    
    # # 标准化训练集临床数据
    # train_clinical_scaled <- standardize_with_train_clinical(train_clinical,
    #                                                          train_clinical,
    #                                                          scale_cols = c("Age_at_Diagnosis"))
    # message("标准化训练集临床数据完成：train_clinical_scaled 维度 = ", 
    #         toString(dim(train_clinical_scaled)))
    # 
    # # 标准化测试集临床数据
    # test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
    #                                                         test_clinical,
    #                                                         scale_cols = c("Age_at_Diagnosis"))
    # message("标准化测试集临床数据完成：test_clinical_scaled 维度 = ", 
    #         toString(dim(test_clinical_scaled)))
    # 
    # LASSO Cox
    gene_freq_df <- repeat_cv_lasso_cox(train_expr = train_expr_filtered,
                                        train_clinical,
                                        significant_gene_vec = significant_gene2,
                                        repeats = 5,
                                        nfolds = 10,
                                        alpha = 1)
    message(sprintf("Bootstrap %d: LASSO Cox 完成，gene_freq_df 行数 = %d", 
                    b, nrow(gene_freq_df)))
    
    # 动态筛选变量，满足 EPV 要求
    gene_freq_df_best <- gene_freq_df %>%
      filter(freq >= 0.8)
    max_vars_allowed <- floor(events_train / min_epv)
    message(sprintf("Bootstrap %d: 动态筛选变量完成，gene_freq_df_best 行数 = %d, max_vars_allowed = %d", 
                    b, nrow(gene_freq_df_best), max_vars_allowed))
    
    if (max_vars_allowed == 0 || nrow(gene_freq_df) == 0) {
      message(sprintf("Bootstrap %d skipped: Too few events or no selected genes", b))
      next
    }
    
    selected_gene_df <- gene_freq_df[1:min(nrow(gene_freq_df_best), max_vars_allowed), ] %>%
      mutate(coef = mean_coef)
    message(sprintf("Bootstrap %d: 选择基因完成，selected_gene_df 行数 = %d, EPV = %.2f", 
                    b, nrow(selected_gene_df), events_train / nrow(selected_gene_df)))
    
    # 计算风险分数
    clinical_cleaned_risk_train <- compute_risk_score(
      gene_mat_scaled = train_expr_filtered,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = train_clinical,
      n_group = 3
    )
    # message("训练集风险分数计算完成：clinical_cleaned_risk_train 维度 = ", 
    #         toString(dim(clinical_cleaned_risk_train)))
    # 
    # predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
    #                  "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")
    
    #predictors <- c(predictors0, colnames(clinical_cleaned_risk_train)[15:ncol(clinical_cleaned_risk_train)])
    predictors <- colnames(clinical_cleaned_risk_train)[6:ncol(clinical_cleaned_risk_train)]
    message(sprintf("Bootstrap %d: 预测变量选择完成，predictors 长度 = %d", 
                    b, length(predictors)))
    
    
    # 相关性过滤
    #predictor_data <- clinical_cleaned_risk_train[, predictors]
    #numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
    # filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.9)
    # filtered_data <- cbind(
    #   predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
    #   filtered_numeric_vars
    # )
    # predictors_filtered <- colnames(filtered_data)
    
    
    
    # message("预测变量相关性过滤完成：filtered_data 维度 = ", toString(dim(filtered_data)))
    
    # 拟合 Cox 模型
    results_train <- fit_cox_model(predictors, clinical_cleaned_risk_train)
    message(sprintf("Bootstrap %d: Cox 模型拟合完成，model 是否存在 = %s", 
                    b, !is.null(results_train$model)))
    
    if (is.null(results_train$model)) {
      message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
      next
    }
    
    # 发散检测
    if (any(abs(coef(results_train$model)) > coef_max)) {
      message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
      next
    }
    message("Cox 模型发散检测完成")
    
    # 检查完全分离
    concordance_val <- tryCatch({
      suppressWarnings(summary(results_train$model)$concordance[1])
    }, error = function(e) NA)
    message(sprintf("Bootstrap %d: 完全分离检查完成，concordance_val = %.4f", 
                    b, if (is.na(concordance_val)) NA else concordance_val))
    
    if (!is.na(concordance_val) && concordance_val >= min_concord) {
      message(sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
      next
    }
    
    # 测试集风险分数
    clinical_cleaned_risk_test <- compute_risk_score(
      gene_mat_scaled = test_expr_scaled,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = test_clinical,
      n_group = 3
    )
    message("测试集风险分数计算完成：clinical_cleaned_risk_test 维度 = ", 
            toString(dim(clinical_cleaned_risk_test)))
    
    # Cox 模型评估
    result_valid <- calculate_time_auc_cindex(
      "Cox", 
      fitted_model = results_train$model, 
      df = clinical_cleaned_risk_test
    )
    message(sprintf("Bootstrap %d: Cox 模型评估完成，iAUC = %.4f, c_index = %.4f", 
                    b, result_valid$iAUC, result_valid$c_index))
    
    # Random Survival Forest
    
    clinical_rsf <- clinical_cleaned_risk_train %>%
      select(-c("geo_accession", "risk_score","risk_group" ))
    
    result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000)
    rsf_fit_best <- result_rsf_train$best_model
    message("RSF 模型拟合完成")
    
    result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
    message(sprintf("Bootstrap %d: RSF 模型评估完成，iAUC = %.4f, c_index = %.4f", 
                    b, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    
    # 保存性能和基因
    perf_list[[b]] <- list(
      cox = result_valid,
      rsf = result_rsf_valid
    )
    gene_list[[b]] <- selected_gene_df$gene
    rsf_predictor_list[[b]] <- predictors_filtered
    message("性能和基因保存完成")
    
    # 计算 Cox 和 RSF 综合评分
    cox_score <- if (!is.na(result_valid$iAUC) && !is.na(result_valid$c_index)) {
      result_valid$iAUC + result_valid$c_index
    } else {
      -Inf
    }
    rsf_score <- if (!is.na(result_rsf_valid$iAUC) && !is.na(result_rsf_valid$c_index)) {
      result_rsf_valid$iAUC + result_rsf_valid$c_index
    } else {
      -Inf
    }
    cat(sprintf("Cox_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
                b, cox_score, result_valid$iAUC, result_valid$c_index))
    cat(sprintf("rsf_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
                b, rsf_score, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    
    # 判断是否更优（Cox）
    if (cox_score > best_perf) {
      best_perf    <- cox_score
      best_model   <- results_train$model
      best_genes   <- selected_gene_df$gene
      best_iter    <- b
      best_seed    <- current_seed
      best_method  <- "Cox"
    }
    
    # 判断是否更优（RSF）
    if (rsf_score > best_perf) {
      best_perf    <- rsf_score
      best_model   <- rsf_fit_best
      best_genes   <- predictors_filtered
      best_iter    <- b
      best_seed    <- current_seed
      best_method  <- "RSF"
    }
    
    cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
    cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
    message(sprintf("Bootstrap %d: 迭代完成", b))
  }
  
  message("所有 Bootstrap 迭代完成")
  
  # 汇总
  perf_list_nz <- Filter(Negate(is.null), perf_list)
  
  cox_iAUC <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$iAUC)) x$cox$iAUC else NA_real_)
  cox_cidx <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$c_index)) x$cox$c_index else NA_real_)
  
  rsf_iAUC <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$iAUC)) x$rsf$iAUC else NA_real_)
  rsf_cidx <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$c_index)) x$rsf$c_index else NA_real_)
  
  mean_cox_iAUC <- mean(cox_iAUC, na.rm = TRUE)
  mean_cox_cidx <- mean(cox_cidx, na.rm = TRUE)
  mean_rsf_iAUC <- mean(rsf_iAUC, na.rm = TRUE)
  mean_rsf_cidx <- mean(rsf_cidx, na.rm = TRUE)
  
  message("性能指标汇总完成：mean_cox_iAUC = ", mean_cox_iAUC, 
          ", mean_cox_cidx = ", mean_cox_cidx,
          ", mean_rsf_iAUC = ", mean_rsf_iAUC,
          ", mean_rsf_cidx = ", mean_rsf_cidx)
  
  all_genes <- unlist(gene_list)
  gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  message("基因频率计算完成")
  
  return(list(
    mean_cox_iAUC = mean_cox_iAUC,
    mean_cox_cidx = mean_cox_cidx,
    mean_rsf_cidx = mean_rsf_cidx,
    mean_rsf_iAUC = mean_rsf_iAUC,
    gene_frequency = gene_freq,
    all_results = perf_list,
    best_model = best_model,
    best_perf = best_perf,
    best_predictors = best_genes,
    best_iter = best_iter,
    best_seed = best_seed,
    best_method = best_method,
    train_indices_list = train_indices_list
  ))
}



res39582_pur_gen <- gen_run_bootstrap_validation_safe(expr_mat, clinical_df, 
                                            B = 10, 
                                            seed = 100, 
                                            min_epv = 4, 
                                            coef_max = 10, 
                                            min_concord = 0.97)

