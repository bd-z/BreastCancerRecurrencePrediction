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

# English comments
col_mean <- mean(clinical_cleaned_39582_imputed$Age_at_Diagnosis, na.rm = TRUE)
clinical_cleaned_39582_imputed$Age_at_Diagnosis[is.na(clinical_cleaned_39582_imputed$Age_at_Diagnosis)] <- col_mean

md.pattern(df_miss)


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



run_bootstrap_validation_safe <- function(expr_mat, clinical_df, 
                                          B = 10, 
                                          #max_vars = 30, 
                                          seed = 10000,
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
  
  for (b in 1:B) {
    
    # Bootstrap 抽样
    current_seed <- seed + b 
    set.seed(current_seed) 
    message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
    
    train_idx <- sample(seq_len(n), size = n, replace = TRUE)
    test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
    
    train_indices_list[[b]] <- train_idx
    
    if (length(test_idx) == 0) next
    
    train_expr     <- expr_mat[, train_idx]
    test_expr      <- expr_mat[, test_idx]
    train_clinical <- clinical_df[train_idx, ]
    test_clinical  <- clinical_df[test_idx, ]
    
    # 事件数
    events_train <- sum(train_clinical$e_dmfs)
    message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
    
    # Step 1: 单变量 Cox
    sig_gene_df <- batch_univariate_cox_regression(train_expr, train_clinical)
    # sig_gene_df <- sig_gene_df[sig_gene_df$p.value < 0.05, ]
    # if (nrow(sig_gene_df) == 0) next
    
    #sig_gene_df_39582 <- sig_gene_df
    
    # 动态设置 p-value 阈值
    p_thresh <- if (events_train < 40) {
      0.01
    } else {
      0.05
    }
    
    sig_gene_df <- sig_gene_df[sig_gene_df$p.value < p_thresh, ]
    if (nrow(sig_gene_df) == 0) {
      message(sprintf("Bootstrap %d skipped: No genes passed p < %.2f (Events: %d)", b, p_thresh, events_train))
      next
    }
    
    
    
    
    significant_gene = sig_gene_df$gene
    
    # scale train gene data
    train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr,
                                                gene_mat_valid= train_expr,
                                                significant_gene = significant_gene)
    
    # scale test data
    test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr,
                                               gene_mat_valid= test_expr,
                                               significant_gene = significant_gene)
    
    
    # scale train clinical data
    train_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                             train_clinical,
                                                             scale_cols = c("Age_at_Diagnosis"))
    
    
    # scale test data
    test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                            test_clinical,
                                                            scale_cols = c("Age_at_Diagnosis"))
    
    # Step 2: LASSO Cox
    selected_gene_df <- lasso_cox_cv(train_expr_scaled, train_clinical, sig_gene_df)
    
    # 限制变量数
    # if (nrow(selected_gene_df) > max_vars) {
    #   selected_gene_df <- selected_gene_df[order(abs(selected_gene_df$coef), decreasing = TRUE), ]
    #   selected_gene_df <- selected_gene_df[1:max_vars, ]
    # }
    #############################
    # # 动态限制变量数：最多不超过 events_train / min_epv
    # max_vars_allowed <- floor(events_train / min_epv)
    # if (nrow(selected_gene_df) > max_vars_allowed) {
    #   selected_gene_df <- selected_gene_df[order(abs(selected_gene_df$coef), decreasing = TRUE), ]
    #   selected_gene_df <- selected_gene_df[1:max_vars_allowed, ]
    # }
    # 
    # 
    # # EPV 检查
    # message(sprintf("Bootstrap %d: Events = %d, Vars = %d, EPV = %.2f", 
    #                 b, events_train, nrow(selected_gene_df), events_train / nrow(selected_gene_df)))
    # 
    # if (events_train / nrow(selected_gene_df) < min_epv) {
    #   message(sprintf("Bootstrap %d skipped: EPV too low (%0.2f)", b, events_train / nrow(selected_gene_df)))
    #   next
    # }
    #############################
    # Step: 动态筛选变量，满足 EPV 要求
    max_vars_allowed <- floor(events_train / min_epv)
    
    if (max_vars_allowed == 0 || nrow(selected_gene_df) == 0) {
      message(sprintf("Bootstrap %d skipped: Too few events or no selected genes", b))
      next
    }
    
    # 限制变量数量（最多 max_vars_allowed 个）
    selected_gene_df <- selected_gene_df[order(abs(selected_gene_df$coef), decreasing = TRUE), ]
    selected_gene_df <- selected_gene_df[1:min(nrow(selected_gene_df), max_vars_allowed), ]
    
    # 打印信息
    message(sprintf("Bootstrap %d: Events = %d, Vars = %d, EPV = %.2f", 
                    b, events_train, nrow(selected_gene_df), events_train / nrow(selected_gene_df)))
    
    
    
    # Step 3: 风险分数
    clinical_cleaned_risk_train <- compute_risk_score(
      gene_mat_scaled = train_expr_scaled,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = train_clinical_scaled,
      n_group = 3
    )
    
    predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
                     "Chemotherapy_Adjuvant", "MMR_Status",
                     "KRAS_Mutation")
    #    "grade", "er", "age", "size"
    
    predictors  <- c(predictors0, colnames(clinical_cleaned_risk_train)[15:ncol(clinical_cleaned_risk_train)])
    
    # 相关性过滤
    predictor_data <- clinical_cleaned_risk_train[, predictors]
    numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
    filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.85)
    filtered_data <- cbind(
      predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
      filtered_numeric_vars
    )
    
    predictors_filtered <- colnames(filtered_data)
    df <- cbind(filtered_data, clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs")])
    
    # 拟合模型
    results_train <- fit_cox_model(predictors_filtered, df)
    
    if (is.null(results_train$model)) {
      message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
      next
    }
    
    
    # 发散检测
    if (any(abs(coef(results_train$model)) > coef_max)) {
      message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
      next
    }
    
    # 检查完全分离
    concordance_val <- tryCatch({
      suppressWarnings(summary(results_train$model)$concordance[1])
    }, error = function(e) NA)
    if (!is.na(concordance_val) && concordance_val >= min_concord) {
      message(
        sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
      next
    }
    
    # Step 4: 测试集评估
    clinical_cleaned_risk_test <- compute_risk_score(
      gene_mat_scaled = test_expr_scaled,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = test_clinical_scaled,
      n_group = 3
    )
    
    result_valid <- calculate_time_auc_cindex(
      "Cox", 
      fitted_model = results_train$model, 
      df = clinical_cleaned_risk_test
    )
    
    # perf_list[[b]] <- result_valid
    # gene_list[[b]] <- selected_gene_df$gene
    # 
    # # 保存最优模型
    # if (!is.na(result_valid$iAUC) && result_valid$iAUC > best_perf) {
    #   best_perf  <- result_valid$iAUC
    #   best_model <- results_train$model
    #   best_genes <- selected_gene_df$gene
    # }
    
    #
    # random forest
    clinical_rsf <- df
    
    # build RSF model
    result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5, ntree = 1000, seed = current_seed)
    # Best model
    rsf_fit_best <- result_rsf_train$best_model
    
    #Remove variables with negative/less importance from the predictors vector based on RSF model output
    imp <- result_rsf_train$importance
    vars_to_remove <- names(imp)[imp < 0.01]
    predictors_filtered_rsf <- setdiff(predictors_filtered, vars_to_remove)
    
    clinical_rsf <- df[, c("t_dmfs", "e_dmfs", predictors_filtered_rsf)]
    # build RSF model with predictors_filtered_rsf
    result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5, ntree = 1000, seed = current_seed)
    # Best model
    rsf_fit_best <- result_rsf_train$best_model
    result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
    
    
    # 把 Cox + RSF 一起保存
    perf_list[[b]] <- list(
      cox  = result_valid,
      rsf  = result_rsf_valid
    )
    
    # 保存基因
    gene_list[[b]] <- selected_gene_df$gene
    rsf_predictor_list[[b]] <- predictors_filtered_rsf
    
    
    # 当前 Cox 和 RSF 的综合评分（iAUC + C-index）
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
    cat(sprintf("Cox_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f",
                b,
                cox_score,
                result_valid$iAUC,
                result_valid$c_index
                ))
    cat(sprintf("rsf_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f",
                b,
                rsf_score,
                result_rsf_valid$iAUC,
                result_rsf_valid$c_index))
    
    
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
      best_genes   <- predictors_filtered_rsf
      best_iter    <- b
      best_seed    <- current_seed
      best_method  <- "RSF"
    }
    
    
    cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
    cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
  }
  
  # 汇总
 # iAUCs   <- sapply(perf_list, function(x) x$iAUC)
  #cindexs <- sapply(perf_list, function(x) x$c_index)
  
  
  
  ######################
  # 过滤掉为 NULL 的轮次
  perf_list_nz <- Filter(Negate(is.null), perf_list)
  
  # 取 Cox 指标
  cox_iAUC  <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$iAUC)) x$cox$iAUC else NA_real_)
  cox_cidx  <- sapply(perf_list_nz, function(x)
    if (!is.null(x$cox) && !is.null(x$cox$c_index)) x$cox$c_index else NA_real_)
  
  # 取 RSF 指标
  rsf_iAUC  <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$iAUC)) x$rsf$iAUC else NA_real_)
  rsf_cidx  <- sapply(perf_list_nz, function(x)
    if (!is.null(x$rsf) && !is.null(x$rsf$c_index)) x$rsf$c_index else NA_real_)
  
  # 汇总（去掉 NA）
  mean_cox_iAUC <- mean(cox_iAUC, na.rm = TRUE)
  mean_cox_cidx <- mean(cox_cidx, na.rm = TRUE)
  mean_rsf_iAUC <- mean(rsf_iAUC, na.rm = TRUE)
  mean_rsf_cidx <- mean(rsf_cidx, na.rm = TRUE)
  
  
  ######################################################
  
  all_genes <- unlist(gene_list)
  gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  
  return (list(
    #mean_iAUC   = mean(iAUCs, na.rm = TRUE),
    #mean_cindex = mean(cindexs, na.rm = TRUE),
    mean_cox_iAUC = mean_cox_iAUC,
    mean_cox_cidx = mean_cox_cidx,
    mean_rsf_cidx = mean_rsf_cidx,
    mean_rsf_iAUC = mean_rsf_iAUC,
    gene_frequency = gene_freq,
    all_results = perf_list,
    best_model  = best_model,   # ⬅ 最优模型对象
    best_perf   = best_perf,    # ⬅ 最优 iAUC
    best_predictors = best_genes,   # ⬅ 最优模型的基因列表或者RSF预测变量
    best_iter = best_iter,
    best_seed = best_seed,
    best_method = best_method,
    train_indices_list = train_indices_list
  ))
}


res39582_f <- run_bootstrap_validation_safe(expr_mat,
                                      clinical_df, 
                                      B = 30,
                                      seed = 80000,
                                      # max_vars = 20,
                                      min_epv = 4,
                                      coef_max = 10,
                                      min_concord = 0.98) 

saveRDS(res39582_c, file = "train_result_39582_d.rds")

res39582$best_model   # 最优 Cox 模型对象
res39582$best_perf    # 最优 iAUC
res39582$best_genes   # 对应的基因列表
res39582$mean_iAUC

res39582$mean_cindex
res39582$gene_frequency

res39582$best_predictors
res39582$best_iter


res39582_c$










# 加载必要的包
library(dplyr)
library(tidyr)
library(stringr)

# 假设 'clinical_selected' 是原始数据框
clinical_clean <- clinical_selected %>%
  # 选择以 'characteristics_ch' 开头的列
  select(starts_with("characteristics_ch")) %>%
  # 使用 str_replace_all 来去掉字段中的无关文字
  mutate(across(everything(), ~ str_replace_all(., ": |\\(", ""))) %>%
  # 分割每一列成两列，前部分为特征，后部分为值
  separate(col = characteristics_ch1.2, into = c("Sex", "Sex_value"), sep = "(?<=Sex)", extra = "merge") %>%
  separate(col = characteristics_ch1.3, into = c("Age_at_diagnosis", "Age_value"), sep = "(?<=diagnosis)", extra = "merge") %>%
  separate(col = characteristics_ch1.4, into = c("TNM_stage", "Stage_value"), sep = "(?<=stage)", extra = "merge") %>%
  separate(col = characteristics_ch1.5, into = c("TNM_t", "T_value"), sep = "(?<=t)", extra = "merge") %>%
  separate(col = characteristics_ch1.6, into = c("TNM_n", "N_value"), sep = "(?<=n)", extra = "merge") %>%
  separate(col = characteristics_ch1.7, into = c("TNM_m", "M_value"), sep = "(?<=m)", extra = "merge") %>%
  separate(col = characteristics_ch1.8, into = c("Tumor_location", "Location_value"), sep = "(?<=location)", extra = "merge") %>%
  separate(col = characteristics_ch1.9, into = c("Chemotherapy_adjuvant", "Adjuvant_value"), sep = "(?<=adjuvant)", extra = "merge") %>%
  separate(col = characteristics_ch1.13, into = c("OS_event", "OS_event_value"), sep = "(?<=event)", extra = "merge") %>%
  separate(col = characteristics_ch1.14, into = c("OS_delay", "OS_delay_value"), sep = "(?<=delay)", extra = "merge") %>%
  separate(col = characteristics_ch1.15, into = c("MMR_status", "MMR_status_value"), sep = "(?<=status)", extra = "merge") %>%
  separate(col = characteristics_ch1.16, into = c("CIMP_status", "CIMP_status_value"), sep = "(?<=status)", extra = "merge") %>%
  separate(col = characteristics_ch1.17, into = c("CIN_status", "CIN_status_value"), sep = "(?<=status)", extra = "merge") %>%
  separate(col = characteristics_ch1.18, into = c("TP53_mutation", "TP53_mutation_value"), sep = "(?<=mutation)", extra = "merge") %>%
  separate(col = characteristics_ch1.22, into = c("KRAS_mutation", "KRAS_mutation_value"), sep = "(?<=mutation)", extra = "merge") %>%
  separate(col = characteristics_ch1.26, into = c("BRAF_mutation", "BRAF_mutation_value"), sep = "(?<=mutation)", extra = "merge") %>%
  separate(col = characteristics_ch1.30, into = c("CIT_molecular_subtype", "CIT_subtype_value"), sep = "(?<=subtype)", extra = "merge")

# 进一步将数值转换为数值型（如适用）
clinical_clean <- clinical_clean %>%
  mutate(
    Age_value = as.numeric(Age_value),
    RFS_delay_value = as.numeric(RFS_delay_value),
    OS_delay_value = as.numeric(OS_delay_value)
  )

# 查看清理后的数据
head(clinical_clean)








# 加载必需的包
library(dplyr)
library(tidyr)

# 假设 `clinical_selected` 是原始数据框
clinical_clean <- clinical_selected %>%
  # 选取特定的列
  select(starts_with("characteristics_ch")) %>%
  # 将每列拆分为两个部分：特征和数值
  mutate(across(everything(), ~ gsub(": ", "", .))) %>%
  separate(col = characteristics_ch1.2, into = c("Sex", "Sex_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.3, into = c("Age_at_diagnosis", "Age_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.4, into = c("TNM_stage", "Stage_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.5, into = c("TNM_t", "T_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.6, into = c("TNM_n", "N_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.7, into = c("TNM_m", "M_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.8, into = c("Tumor_location", "Location_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.9, into = c("Chemotherapy_adjuvant", "Adjuvant_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.11, into = c("RFS_event", "RFS_event_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.12, into = c("RFS_delay", "RFS_delay_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.13, into = c("OS_event", "OS_event_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.14, into = c("OS_delay", "OS_delay_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.15, into = c("MMR_status", "MMR_status_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.16, into = c("CIMP_status", "CIMP_status_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.17, into = c("CIN_status", "CIN_status_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.18, into = c("TP53_mutation", "TP53_mutation_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.22, into = c("KRAS_mutation", "KRAS_mutation_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.26, into = c("BRAF_mutation", "BRAF_mutation_value"), sep = ": ", extra = "merge") %>%
  separate(col = characteristics_ch1.30, into = c("CIT_molecular_subtype", "CIT_subtype_value"), sep = ": ", extra = "merge")

# 查看清理后的数据
head(clinical_clean)
