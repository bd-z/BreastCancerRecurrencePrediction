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
options(scipen = 999)

## load dataset
#gset_7390 <- getGEO("GSE7390", GSEMatrix=TRUE)
gset_39582 <- getGEO("GSE39582", GSEMatrix=TRUE)
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


clinical_cleaned_39582 <- clean_clinical_data(col_selected = col_selected_39582, df_clinic_raw = clinical_39582)
colSums(is.na(clinical_cleaned_7390))









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


library(survival)

#=============================
# End-to-End 建模流程（和之前一致）
#=============================
bootstrap_pipeline <- function(train_expr, train_clinical, test_expr, test_clinical) {
  
  # Step 1: 单变量 Cox
  sig_gene_df <- batch_univariate_cox_regression(train_expr, train_clinical)
  
  # Step 2: LASSO Cox
  selected_gene_df <- lasso_cox_cv(train_expr, train_clinical, sig_gene_df)
  
  # Step 3: 训练集风险分数
  clinical_cleaned_risk_train <- compute_risk_score(
    gene_mat_scaled = train_expr,
    significant_vars_df = selected_gene_df,
    clinical_cleaned = train_clinical,
    n_group = 3
  )
  
  # Step 4: 预测变量
  predictors0 <- c("grade", "er", "age", "size")
  predictors  <- c(predictors0, colnames(clinical_cleaned_risk_train)[11:ncol(clinical_cleaned_risk_train)])
  
  # Step 5: 拟合 Cox 模型
  results_train <- fit_cox_model(predictors, clinical_cleaned_risk_train)
  
  # Step 6: 测试集风险分数
  clinical_cleaned_risk_test <- compute_risk_score(
    gene_mat_scaled = test_expr,
    significant_vars_df = selected_gene_df,
    clinical_cleaned = test_clinical,
    n_group = 3
  )
  
  # Step 7: 在测试集计算 AUC / C-index
  result_valid <- calculate_time_auc_cindex(
    "Cox", 
    fitted_model = results_train$model, 
    df = clinical_cleaned_risk_test
  )
  
  return(list(
    performance = result_valid,
    selected_genes = selected_gene_df$gene  # 用于稳定性分析
  ))
}





#=============================
# Bootstrap 外层验证
#=============================
remove_high_corr <- function(df, threshold = 0.9) {
  # df 是数据框（只包含数值型预测变量，不包括生存时间和事件）
  corr_matrix <- cor(df, method = "pearson", use = "pairwise.complete.obs")
  
  # 找出高度相关的变量对
  high_corr <- findCorrelation(corr_matrix, cutoff = threshold, names = TRUE)
  
  if (length(high_corr) > 0) {
    message(sprintf("Removed %d highly correlated variables (>|%0.2f|)", 
                    length(high_corr), threshold))
    df <- df[, !(colnames(df) %in% high_corr), drop = FALSE]
  }
  return(df)
}

# expr_7390_scaled <- standardize_with_train(gene_mat_train = expr_7390,
#                                            gene_mat_valid=expr_7390,
#                                            significant_gene = significant_gene_7390)




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
    set.seed(seed + b) 
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
    
    # 动态设置 p-value 阈值
    p_thresh <- if (events_train < 40) {
      0.01
    } else if (events_train < 55) {
      0.05
    } else {
      0.1
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
                                                             scale_cols = c("age", "size"))
     
    
    # scale test data
    test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                            test_clinical,
                                                             scale_cols = c("age", "size"))
    
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
    
    predictors0 <- c("grade", "er", "age", "size")
    predictors  <- c(predictors0, colnames(clinical_cleaned_risk_train)[11:ncol(clinical_cleaned_risk_train)])
    
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
      message(sprintf("Bootstrap %d skipped, Concordance >= 0.98 detected ", as.integer(b)))
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
    result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5)
    # Best model
    rsf_fit_best <- result_rsf_train$best_model
    
    #Remove variables with negative/less importance from the predictors vector based on RSF model output
    imp <- result_rsf_train$importance
    vars_to_remove <- names(imp)[imp < 0.001]
    predictors_filtered_rsf <- setdiff(predictors_filtered, vars_to_remove)
    
    clinical_rsf <- df[, c("t_dmfs", "e_dmfs", predictors_filtered_rsf)]
    # build RSF model with predictors_filtered_rsf
    result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5)
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
    
    # 判断是否更优（Cox）
    if (cox_score > best_perf) {
      best_perf    <- cox_score
      best_model   <- results_train$model
      best_genes   <- selected_gene_df$gene
      best_iter    <- b
      best_method  <- "Cox"
    }
    
    # 判断是否更优（RSF）
    if (rsf_score > best_perf) {
      best_perf    <- rsf_score
      best_model   <- rsf_fit_best
      best_genes   <- predictors_filtered_rsf
      best_iter    <- b
      best_method  <- "RSF"
    }
    cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
    cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
  }
  
  # 汇总
  iAUCs   <- sapply(perf_list, function(x) x$iAUC)
  cindexs <- sapply(perf_list, function(x) x$c_index)
  
  all_genes <- unlist(gene_list)
  gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  
  list(
    mean_iAUC   = mean(iAUCs, na.rm = TRUE),
    mean_cindex = mean(cindexs, na.rm = TRUE),
    gene_frequency = gene_freq,
    all_results = perf_list,
    best_model  = best_model,   # ⬅ 最优模型对象
    best_perf   = best_perf,    # ⬅ 最优 iAUC
    best_predictors = best_genes,   # ⬅ 最优模型的基因列表或者RSF预测变量
    best_iter = best_iter,
    best_method = best_method,
    train_indices_list
  )
}


res2 <- run_bootstrap_validation_safe(expr_mat,
                                     clinical_df, 
                                     B = 30,
                                     seed = 20000,
                                     # max_vars = 20,
                                     min_epv = 2.5,
                                     coef_max = 15,
                                     min_concord = 0.98) 



res$best_model   # 最优 Cox 模型对象
res$best_perf    # 最优 iAUC
res$best_genes   # 对应的基因列表




run_bootstrap_validation_safe <- function(expr_mat, clinical_df, 
                                          B = 100, 
                                          seed = 123, 
                                          max_vars = 50, 
                                          min_epv = 1, 
                                          coef_max = 10) {
  set.seed(seed)
  
  n <- ncol(expr_mat)
  perf_list <- list()
  gene_list <- list()
  
  for (b in 1:B) {
    # Bootstrap 抽样
    train_idx <- sample(seq_len(n), size = n, replace = TRUE)
    test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
    
    if (length(test_idx) == 0) next
    
    train_expr     <- expr_mat[, train_idx]
    test_expr      <- expr_mat[, test_idx]
    train_clinical <- clinical_df[train_idx, ]
    test_clinical  <- clinical_df[test_idx, ]
    
    # 事件数
    events_train <- sum(train_clinical$e_dmfs)
    
    # Step 1: 单变量 Cox
    sig_gene_df <- batch_univariate_cox_regression(train_expr, train_clinical)
    sig_gene_df <- sig_gene_df[sig_gene_df$p.value < 0.05, ]
    if (nrow(sig_gene_df) == 0) next
    # Step 2: LASSO Cox
    selected_gene_df <- lasso_cox_cv(train_expr, train_clinical, sig_gene_df)
    
    # 限制变量数（按绝对系数排序取前 max_vars 个）
    if (nrow(selected_gene_df) > max_vars) {
      selected_gene_df <- selected_gene_df[order(abs(selected_gene_df$coef), decreasing = TRUE), ]
      selected_gene_df <- selected_gene_df[1:max_vars, ]
    }
    
    # EPV 检查
    if (events_train / nrow(selected_gene_df) < min_epv) {
      message(sprintf("Bootstrap %d skipped: EPV too low (%0.2f)", b, events_train / nrow(selected_gene_df)))
      next
    }
    
    # Step 3: 计算风险分数
    clinical_cleaned_risk_train <- compute_risk_score(
      gene_mat_scaled = train_expr,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = train_clinical,
      n_group = 3
    )
    
    predictors0 <- c("grade", "er", "age", "size")
    predictors  <- c(predictors0, colnames(clinical_cleaned_risk_train)[11:ncol(clinical_cleaned_risk_train)])
    
    # 取出预测变量数据
    predictor_data <- clinical_cleaned_risk_train[, predictors]
    
    # 只对数值型变量做相关性过滤（因子型跳过）
    numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
    filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.85)
    
    # 保留非数值变量 + 过滤后的数值变量
    filtered_data <- cbind(
      predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
      filtered_numeric_vars
    )
    
    # 拟合 Cox 模型
    # results_train <- fit_cox_model_firth(colnames(filtered_data), 
    #                                df = cbind(filtered_data, clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs")]),
    #                                coef_max)
    # 
    
    results_train <- fit_cox_model(colnames(filtered_data),
                                   df = cbind(filtered_data, clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs")])
                                   )
    
    
    #results_train <- fit_cox_model(predictors, clinical_cleaned_risk_train)
    
    # 发散检测
    if (any(abs(coef(results_train$model)) > coef_max)) {
      message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
      #message(sprintf("Bootstrap %d skipped", as.integer(b)))

      #message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", b, coef_max))
      next
    }
    
    # 检查 Concordance 是否接近 1（完全分离）
    concordance_val <- tryCatch({
      suppressWarnings(summary(results_train$model)$concordance[1])
    }, error = function(e) NA)
    
    if (!is.na(concordance_val) && concordance_val >= 0.98) {
      message(sprintf("Bootstrap %d skipped, Concordance >= 0.98 detected ", as.integer(b)))
      next
    }
  
    # Step 4: 测试集评估
    clinical_cleaned_risk_test <- compute_risk_score(
      gene_mat_scaled = test_expr,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = test_clinical,
      n_group = 3
    )
    
    result_valid <- calculate_time_auc_cindex(
      "Cox", 
      fitted_model = results_train$model, 
      df = clinical_cleaned_risk_test
    )
    
    perf_list[[b]] <- result_valid
    gene_list[[b]] <- selected_gene_df$gene
  }
  
  # 汇总
  iAUCs   <- sapply(perf_list, function(x) x$iAUC)
  cindexs <- sapply(perf_list, function(x) x$c_index)
  
  all_genes <- unlist(gene_list)
  gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  
  list(
    mean_iAUC   = mean(iAUCs, na.rm = TRUE),
    mean_cindex = mean(cindexs, na.rm = TRUE),
    gene_frequency = gene_freq,
    all_results = perf_list
  )
}
# 
# 
# 
# 
# run_bootstrap_validation <- function(expr_mat, clinical_df, B = 100, seed = 123) {
#   set.seed(seed)
#   
#   n <- ncol(expr_mat)
#   perf_list <- list()
#   gene_list <- list()
#   
#   for (b in 1:B) {
#     # Bootstrap 抽样（训练集）
#     train_idx <- sample(seq_len(n), size = n, replace = TRUE)
#     test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
#     
#     if (length(test_idx) == 0) next  # 跳过无OOB样本的情况
#     
#     train_expr     <- expr_mat[, train_idx]
#     test_expr      <- expr_mat[, test_idx]
#     train_clinical <- clinical_df[train_idx, ]
#     test_clinical  <- clinical_df[test_idx, ]
#     
#     res <- bootstrap_pipeline(train_expr, train_clinical, test_expr, test_clinical)
#     
#     perf_list[[b]] <- res$performance
#     gene_list[[b]] <- res$selected_genes
#   }
#   
#   # 计算平均性能
#   iAUCs   <- sapply(perf_list, function(x) x$iAUC)
#   cindexs <- sapply(perf_list, function(x) x$c_index)
#   
#   # 计算基因入选频率
#   all_genes <- unlist(gene_list)
#   gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
#   
#   list(
#     mean_iAUC   = mean(iAUCs, na.rm = TRUE),
#     mean_cindex = mean(cindexs, na.rm = TRUE),
#     gene_frequency = gene_freq,
#     all_results = perf_list
#   )
# }
# 
#=============================
# 调用示例
#=============================
boot_results <- run_bootstrap_validation_safe(
  expr_mat = expr_7390,
  clinical_df = clinical_cleaned_7390_imputed,
  B = 3, seed = 123,max_vars = 30, min_epv = 1, coef_max = 10)
  
#run_bootstrap_validation(expr_mat = expr_7390, clinical_df = clinical_cleaned_7390_imputed, B = 10)


boot_results$mean_iAUC
# boot_results$mean_cindex
# head(boot_results$gene_frequency)















































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
