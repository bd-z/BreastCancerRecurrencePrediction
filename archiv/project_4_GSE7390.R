
library(GEOquery)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)

options(scipen = 999)

# library(GEOquery)
# gset <- getGEO("GSE7390", GSEMatrix=TRUE)
# expr <- exprs(gset[[1]])
# clinical <- pData(gset[[1]])

#saveRDS(expr_7390, file = "expr_7390.RDS")
#saveRDS(clinical_7390, file = "clinical_7390.rds")

expr <- readRDS("expr_7390.RDS")
clinical <- readRDS("clinical_7390.rds")

# 查看列名确定生存变量位置
colnames(clinical)

clinical_selected <- clinical[, c(
  "geo_accession",
  "t.dmfs:ch1", "e.dmfs:ch1",
  #"t.os:ch1", "e.os:ch1",
  #"t.rfs:ch1", "e.rfs:ch1",
  #"t.tdm:ch1", "e.tdm:ch1",
  "age:ch1",
  "er:ch1",
  "grade:ch1",
  "size:ch1",
  "node:ch1",
  "hospital:ch1", "samplename:ch1", "filename:ch1"
)]

original_names <- colnames(clinical_selected)
clean_names <- gsub(":ch1", "", original_names)
clean_names <- gsub("\\.", "_", clean_names)
colnames(clinical_selected) <- clean_names

clinical_cleaned <- as.data.frame(lapply(clinical_selected, function(x) {
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(x_num))) x else x_num
}))

for (i in 4:ncol(clinical_cleaned)) {
  if (length(unique(clinical_cleaned[, i])) <= 4) {
    clinical_cleaned[, i] <- factor(clinical_cleaned[, i])
  }
}

clinical_cleaned_7390 <- clinical_cleaned %>%
  select(geo_accession, t_dmfs, e_dmfs, age, er, grade, size, hospital)
hospital_counts <- table(clinical_cleaned_7390$hospital)
# 
# # 提取DMFS时间与事件
# clinical$time <- as.numeric(clinical$`t.dmfs:ch1`)
# clinical$status <- as.numeric(clinical$`e.dmfs:ch1`)
# 
# # 提取生存时间、事件状态及重要协变量
# clinical_selected <- clinical[, c(
#   "geo_accession",
#   "t.dmfs:ch1",
#   "e.dmfs:ch1",
#   
#   # DMFS 时间与事件
#   "t.os:ch1", "e.os:ch1",         # OS 时间与事件
#   "t.rfs:ch1", "e.rfs:ch1",       # RFS 时间与事件
#   "t.tdm:ch1", "e.tdm:ch1",       # 首次远处转移时间与事件
#   "age:ch1",                      # 年龄
#   "er:ch1",                       # 雌激素受体状态
#   "grade:ch1",                    # 肿瘤分级
#   "size:ch1",                     # 肿瘤大小
#   "node:ch1",                     # 淋巴结状态
#   #"Histtype:ch1",                 # 组织学类型
#   #"Angioinv:ch1",                 # 血管侵袭
#   #"Lymp_infil:ch1",               # 淋巴浸润
#   #"NPI:ch1",                      # 诺丁汉预后指数
#   #"veridex_risk:ch1",             # Veridex 风险评分
#   #"risknpi:ch1",                   # NPI 风险组\
#   "hospital:ch1", "samplename:ch1", "filename:ch1"
# )]
# 
# head(clinical_selected)
# 
# original_names <- colnames(clinical_selected)
# clean_names <- gsub(":ch1", "", original_names)
# clean_names <- gsub("\\.", "_", clean_names)
# colnames(clinical_selected) <- clean_names
# 
# 
# #write.csv(clinical_selected, "GSE7390_clinical_selected.csv", row.names = FALSE)
# 
# #clinical_selected <- read.csv("GSE7390_clinical_selected.csv")
# 
# 
# clinical_cleaned <- as.data.frame(lapply(clinical_selected, function(x) {
#   x_num <- suppressWarnings(as.numeric(as.character(x)))
#   if (all(is.na(x_num))) x else x_num
# })) #%>%
#   #select(-Angioinv, -Lymp_infil, -Histtype, - risknpi, -NPI)# %>%
#   #na.omit()
# colSums(is.na(clinical_cleaned))
# 
# 
# 
# 
# 
# 
# #table(clinical_cleaned$Histtype, useNA = "ifany")
# # Detect and define factor variables appropriately:
# for (i in 10:ncol(clinical_cleaned)) {
#   if (length(unique(clinical_cleaned[, i])) <= 4) {
#     clinical_cleaned[, i] <- factor(clinical_cleaned[, i])
#   }
# }
# 
# clinical_cleaned_7390 <- clinical_cleaned %>%
#   select(geo_accession, t_os, e_os, 
#          age, er, grade, size) # , hospital
# 

########################################
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
#saveRDS(gset_39582, file = "gset_39582.rds")
#gset_39582 <- readRDS(file = "gset_39582.rds")

expr_39582 <- exprs(gset_39582[[1]])
clinical_39582 <- pData(gset_39582[[1]])
##### 
# col_selected_39582 <- c(
#   "characteristics_ch1.2",
#   "characteristics_ch1.3",
#   "characteristics_ch1.4",
#   "characteristics_ch1.5",
#   "characteristics_ch1.6",
#   "characteristics_ch1.7",
#   "characteristics_ch1.8",
#   "characteristics_ch1.9",
#   "characteristics_ch1.11",
#   "characteristics_ch1.12",
#   "characteristics_ch1.13",
#   "characteristics_ch1.14",
#   "characteristics_ch1.15",
#   "characteristics_ch1.16",
#   "characteristics_ch1.17",
#   "characteristics_ch1.18",
#   "characteristics_ch1.22",
#   "characteristics_ch1.26",
#   "characteristics_ch1.30",
#   "source_name_ch1"
# )
# 
# 
# 
# clinical_selected_39582 <- clinical_39582[col_selected_39582] %>%
#   filter(
#     source_name_ch1 != "Frozen tissue of non tumoral colorectal mucosa"
#   )
# 
# 
# 
# 
# 
# #创建一个函数来清理前缀
# 
# clean_prefix <- function(x, prefix) { str_replace(x, paste0("^", prefix, ": "), "") }
# 
# #数据清理
# 
# cleaned_data_39582 <- clinical_selected_39582 %>% mutate( 
#   
#   Sex = as.factor(clean_prefix(characteristics_ch1.2, "Sex")),
#   # 清理诊断时年龄并转换为数值型
#   Age_at_Diagnosis = as.numeric(clean_prefix(characteristics_ch1.3, "age.at.diagnosis \\(year\\)")),
#   
#   # 清理TNM分期并转换为有序因子
#   TNM_Stage = factor(clean_prefix(characteristics_ch1.4, "tnm.stage"), 
#                      levels = c("1", "2", "3", "4"), ordered = TRUE),
#   
#   # 清理T分期并转换为有序因子
#   TNM_T = factor(clean_prefix(characteristics_ch1.5, "tnm.t"), 
#                  levels = c("T1", "T2", "T3", "T4"), ordered = TRUE),
#   
#   # 清理N分期并转换为有序因子
#   TNM_N = factor(clean_prefix(characteristics_ch1.6, "tnm.n"), 
#                  levels = c("N0", "N1", "N2"), ordered = TRUE),
#   
#   # 清理M分期并转换为因子
#   TNM_M = factor(clean_prefix(characteristics_ch1.7, "tnm.m"), 
#                  levels = c("M0", "M1")),
#   
#   # 清理肿瘤位置并转换为因子
#   Tumor_Location = factor(clean_prefix(characteristics_ch1.8, "tumor.location"), 
#                           levels = c("proximal", "distal")),
#   
#   # 清理辅助化疗状态并转换为因子
#   Chemotherapy_Adjuvant = factor(clean_prefix(characteristics_ch1.9, "chemotherapy.adjuvant"), 
#                                  levels = c("N", "Y")),
#   
#   # # 清理无复发生存事件并转换为因子
#   # RFS_Event = factor(clean_prefix(characteristics_ch1.11, "rfs.event"), 
#   #                    levels = c("0", "1")),
#   # 
#   # # 清理无复发生存时间并转换为数值型
#   # RFS_Delay = as.numeric(clean_prefix(characteristics_ch1.12, "rfs.delay")),
#   
#   # 清理总生存事件并转换为因子
#   OS_Event = as.numeric(clean_prefix(characteristics_ch1.13, "os.event")), 
#   # factor(clean_prefix(characteristics_ch1.13, "os.event"), 
#   #                 levels = c("0", "1")),
#   # 
#   # 清理总生存时间并转换为数值型
#   OS_Delay = as.numeric(clean_prefix(characteristics_ch1.14, "os.delay \\(months\\)")) * 30.42,
#   
#   # 清理MMR状态并转换为因子
#   MMR_Status = factor(clean_prefix(characteristics_ch1.15, "mmr.status"), 
#                       levels = c("pMMR", "dMMR")),
#   
#   # 清理CIMP状态并转换为因子
#   CIMP_Status = factor(clean_prefix(characteristics_ch1.16, "cimp.status"), 
#                        levels = c("-", "+")),
#   
#   # 清理CIN状态并转换为因子
#   CIN_Status = factor(clean_prefix(characteristics_ch1.17, "cin.status"), 
#                       levels = c("-", "+")),
#   
#   # 清理TP53突变状态并转换为因子
#   TP53_Mutation = factor(clean_prefix(characteristics_ch1.18, "tp53.mutation"), 
#                          levels = c("WT", "M")),
#   
#   # 清理KRAS突变状态并转换为因子
#   KRAS_Mutation = factor(clean_prefix(characteristics_ch1.22, "kras.mutation"), 
#                          levels = c("WT", "M")),
#   
#   # 清理BRAF突变状态并转换为因子
#   BRAF_Mutation = factor(clean_prefix(characteristics_ch1.26, "braf.mutation"), 
#                          levels = c("WT", "M")),
#   
#   # 清理分子亚型并转换为因子
#   Molecular_Subtype = factor(clean_prefix(characteristics_ch1.30, "cit.molecularsubtype"), 
#                              levels = c("C1", "C2", "C3", "C4", "C5", "C6"))) %>%
#   # 选择清理后的列 #RFS_Event, RFS_Delay,
#   select(Sex, Age_at_Diagnosis, TNM_T, TNM_N, TNM_M, Tumor_Location,
#          Chemotherapy_Adjuvant,
#          MMR_Status, KRAS_Mutation,
#          OS_Event, OS_Delay) # Molecular_Subtype, TNM_Stage, # TP53_Mutation,214个NA BRAF_Mutation, 71N, CIN_Status, 102NA
# #  CIMP_Status, NA70
# 
# # # %>%
# # #   # 处理缺失值：数值型变量用中位数填充，分类变量保持不变
# # #   mutate_if(is.numeric, ~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
# # #   # 确保分类变量的因子水平一致
# # #   mutate_if(is.factor, ~factor(., levels = levels(.)))
# 
# 
# colSums(is.na(cleaned_data_39582))
# 
# clinical_cleaned_39582 <- cleaned_data_39582 %>%
#   filter(!is.na(OS_Delay),
#          !is.na(OS_Event),
#          rowSums(is.na(.)) < 3,
#          OS_Delay > 0) 
# 
# colSums(is.na(clinical_cleaned_39582))
# sum(is.na(clinical_cleaned_39582))
# 
##### 

generate_box_plots(data=clinical_cleaned_7390, continuous_variables = c("age", "t_os", "size")) 

# 查找包含 NA 的列名
na_columns <- colnames(clinical_cleaned_7390)[colSums(is.na(clinical_cleaned_7390)) > 0]
na_columns
selected_col <- colnames(clinical_cleaned_7390)[1:(ncol(clinical_cleaned_7390) - 2)]

#clinical_cleaned_39582_imputed <- impute_missing_value(clinical_cleaned = clinical_cleaned_39582, selected_col, missed_col=na_columns)


#col_mean <- mean(clinical_cleaned_39582_imputed$Age_at_Diagnosis, na.rm = TRUE)
#clinical_cleaned_39582_imputed$Age_at_Diagnosis[is.na(clinical_cleaned_39582_imputed$Age_at_Diagnosis)] <- col_mean

#md.pattern(df_miss)


geo_accession_7390 <- clinical_cleaned_7390$geo_accession
expr_7390_aligned <- expr_7390[, geo_accession_7390, drop = FALSE]



expr_mat_7390=expr_7390_aligned

# clinical_df <- clinical_cleaned_39582_imputed %>%
#   mutate(
#     geo_accession = rownames(.)
#   )


# 
# clinical_df_0 <- clinical_cleaned_39582 %>%
#   rename(
#     e_dmfs = OS_Event,
#     t_dmfs = OS_Delay
#   )

clinical_df_7390 <- clinical_cleaned_7390 %>%
  mutate(
    e_dmfs = e_os,
    t_dmfs = t_os
    ) %>%
  select(-c(e_os, t_os))




# updated with inside imputation

run_bootstrap_validation_3_model <- function(expr_mat, clinical_df, 
                                          B = 10, 
                                          seed = 90000,
                                          min_epv = 2.5, 
                                          coef_max = 10,
                                          min_concord = 0.98) {
  
 
  km_p_list <- list() #save KM log rank test p value for test_data
  train_indices_list <- list()
  gene_list <- list()
  
  # model multimodal
  perf_list_mm <- list()
  rsf_predictor_list_mm <- list()
  best_perf_mm <- -Inf
  best_model_mm <- NULL
  best_iter_mm <- 0
  best_seed_mm <- NULL
  best_method_mm <- "None"
  
  # model gene
  perf_list_mg <- list()
  rsf_predictor_list_mg <- list()
  best_perf_mg <- -Inf
  best_model_mg <- NULL
  best_iter_mg <- 0
  best_seed_mg <- NULL
  best_method_mg <- "None"
  
  # model clinic
  perf_list_mc <- list()
  rsf_predictor_list_mc <- list()
  best_perf_mc <- -Inf
  best_model_mc <- NULL
  best_iter_mc <- 0
  best_seed_mc <- NULL
  best_method_mc <- "None"
  # # model multimodal
  # perf_list_mm <- list()
  # gene_list_mm <- list()
  # rsf_predictor_list_mm <- list()
  # best_perf_mm <- -Inf
  # best_model_mm <- NULL
  # best_genes_mm <- NULL
  # best_iter_mm <- 0
  # best_seed_mm <- NULL
  # best_method_mm <- "None"
  # 
  # # model gene
  # perf_list_mg <- list()
  # gene_list_mg <- list()
  # rsf_predictor_list_mg <- list()
  # best_perf_mg <- -Inf
  # best_model_mg <- NULL
  # best_genes_mg <- NULL
  # best_iter_mg <- 0
  # best_seed_mg <- NULL
  # best_method_mg <- "None"
  # 
  # # model clinic
  # perf_list_mc <- list()
  # gene_list_mc <- list()
  # rsf_predictor_list_mc <- list()
  # best_perf_mc <- -Inf
  # best_model_mc <- NULL
  # best_genes_mc <- NULL
  # best_iter_mc <- 0
  # best_seed_mc <- NULL
  # best_method_mc <- "None"
  # 
  # 
  # perf_list <- list()
  # gene_list <- list()
  # rsf_predictor_list <- list()
  # 
  # 
  # # 保存最优模型信息
  # best_model <- NULL
  # best_perf  <- -Inf
  # best_genes <- NULL
  # best_iter <- NULL
  # best_method <- NULL
  
  
  
  message("初始化完成：expr_mat 维度 = ", toString(dim(expr_mat)), 
          ", clinical_df 维度 = ", toString(dim(clinical_df)))
  
  for (b in 1:B) {
    # Bootstrap 采样
    current_seed <- seed + b 
    set.seed(current_seed) 
    message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
    
    # n <- ncol(expr_mat)
    # train_idx <- sample(seq_len(n), size = n, replace = TRUE)
    # train_indices_list[[b]] <- train_idx
    # test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
    # message(sprintf("Bootstrap %d: train_idx length = %d, test_idx length = %d", 
    #                 b, length(train_idx), length(test_idx)))
    # 
    # if (length(test_idx) == 0) {
    #   message(sprintf("Bootstrap %d skipped: No OOB samples", b))
    #   next
    # }
    # 
    # train_expr     <- expr_mat[, train_idx, drop = FALSE]
    # test_expr      <- expr_mat[, test_idx, drop = FALSE]
    # train_clinical <- clinical_df[train_idx, , drop = FALSE]
    # test_clinical  <- clinical_df[test_idx, , drop = FALSE]
    # 
    # columns_to_remove <- c("e_dmfs", "t_dmfs", "geo_accession", "hospital")
    # selected_col <- setdiff(colnames(clinical_df), columns_to_remove)
    # missed_col <- colnames(clinical_df)[colSums(is.na(clinical_df)) > 0]
    # 
    # res <- impute_fit_apply(
    #   train_df     = train_clinical,
    #   test_df      = test_clinical,
    #   selected_col = selected_col,
    #   missed_col   = missed_col
    # )
    # 
    # train_clinical <- res$train
    # test_clinical  <- res$test
    # 
    # 
    # message("Bootstrap 采样完成：train_expr 维度 = ", toString(dim(train_expr)), 
    #         ", test_expr 维度 = ", toString(dim(test_expr)),
    #         ", train_clinical 维度 = ", toString(dim(train_clinical)),
    #         ", test_clinical 维度 = ", toString(dim(test_clinical)))
    # 
    # if (ncol(train_expr) == 0 || nrow(train_expr) == 0) {
    #   message(sprintf("Bootstrap %d skipped: Invalid train_expr dimensions", b))
    #   next
    # }
    # 
    # # 事件数
    # events_train <- sum(train_clinical$e_dmfs)
    data_prepare_result <- bootstrap_sampl_split_miss_vale_imputation(expr_mat, clinical_df)
    
    if (is.null(data_prepare_result)) {
      message(sprintf("Bootstrap %d skipped: Invalid train_expr dimensions", b))
      next
    }
    train_expr <- data_prepare_result$train_expr
    test_expr <- data_prepare_result$test_expr
    train_clinical <- data_prepare_result$train_clinical
    test_clinical <- data_prepare_result$test_clinical
    events_train <- data_prepare_result$events_train
    
    message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
    message("事件数计算完成")
    
    # MAD 过滤
    # mad_train_expr <- apply(train_expr, 1, mad)
    # cutoff <- quantile(mad_train_expr, 0.1)
    # keep_mad <- mad_train_expr > cutoff
    # message(sprintf("Bootstrap %d: MAD 过滤完成，保留基因数 = %d, cutoff = %.4f", 
    #                 b, sum(keep_mad), cutoff))    # 
    # train_expr2 <- train_expr[keep_mad, , drop = FALSE]
    # test_expr2 <- test_expr[keep_mad, , drop = FALSE]
    filtered_data <- filter_genes_by_mad(train_expr, test_expr, quantile_cutoff = 0.25)
    train_expr2 <- filtered_data$train_expr
    test_expr3 <- filtered_data$test_expr  
    
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
    numeric_columns <- train_clinical %>%
      select(where(is.numeric)) %>%
      colnames
    
    numeric_columns <- setdiff(numeric_columns, columns_to_remove)
    
    train_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                             train_clinical,
                                                             scale_cols = numeric_columns)
    message("标准化训练集临床数据完成：train_clinical_scaled 维度 = ", 
            toString(dim(train_clinical_scaled)))
    
    # 标准化测试集临床数据
    test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                            test_clinical,
                                                            scale_cols = numeric_columns)
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
    
    gene_list[[b]] <- selected_gene_df$gene
    
    result <- compute_risk_score_train(
      gene_mat_scaled  = train_expr_filtered,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = train_clinical_scaled
      )
    # 提取各个组件
    q33_train <- result$q33
    q67_train <- result$q67
    clinical_cleaned_risk_train <- result$merged_df
    
    message("训练集风险分数计算完成：clinical_cleaned_risk_train 维度 = ", 
            toString(dim(clinical_cleaned_risk_train)))
    
    # 测试集风险分数
    clinical_cleaned_risk_test <- compute_risk_score_test(
      gene_mat_scaled = test_expr_scaled,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = test_clinical_scaled,
      q33_train,
      q67_train)
    
    message("测试集风险分数计算完成：clinical_cleaned_risk_test 维度 = ", 
            toString(dim(clinical_cleaned_risk_test)))
    
    km_p <- km_by_group(df = clinical_cleaned_risk_test, group_var = "risk_group")
    km_p_list[[b]] <- km_p
    
    # predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
    #                  "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")
    # 
    # predictors <- c(predictors0, colnames(clinical_cleaned_risk_train)[15:ncol(clinical_cleaned_risk_train)])
    
    columns_to_remove2 <- c("e_dmfs", "t_dmfs", "risk_score", "risk_group")
    character_columns <- clinical_cleaned_risk_train %>%
      select(where(is.character)) %>%
      colnames()
    predictors <- colnames(clinical_cleaned_risk_train) %>% setdiff(c(columns_to_remove2, character_columns))
    
    message(sprintf("Bootstrap %d: 预测变量选择完成，predictors 长度 = %d", 
                    b, length(predictors)))
    df <- clinical_cleaned_risk_train[, c(predictors, "t_dmfs", "e_dmfs")]
    
    # 相关性过滤
    # predictor_data <- clinical_cleaned_risk_train[, predictors]
    # numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
    # filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.9)
    # filtered_data <- cbind(
    #   predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
    #   filtered_numeric_vars
    # )
    # predictors_filtered <- colnames(filtered_data)
    # df <- cbind(filtered_data, clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs")])
    # message("预测变量相关性过滤完成：filtered_data 维度 = ", toString(dim(filtered_data)))
    
    #model multimodal
    model_m <- run_bootstrap_iteration(
      predictors_filtered = predictors,
      df = df,
      b = b,
      coef_max = coef_max,
      min_concord = min_concord,
      clinical_cleaned_risk_test = clinical_cleaned_risk_test,
      selected_gene_df = selected_gene_df,
      perf_list = perf_list_mm,
      rsf_predictor_list = rsf_predictor_list_mm,
      best_perf = best_perf_mm,
      best_model = best_model_mm,
      best_iter = best_iter_mm,
      best_seed = best_seed_mm,
      best_method = best_method_mm,
      current_seed = current_seed
    )
    
    # Update all variables
    perf_list_mm <- model_m$perf_list
    rsf_predictor_list_mm <- model_m$rsf_predictor_list
    best_perf_mm <- model_m$best_perf
    best_model_mm <- model_m$best_model
    best_iter_mm <- model_m$best_iter
    best_seed_mm <- model_m$best_seed
    best_method_mm <- model_m$best_method
    
    
    #model gene
    predictors_g <- predictors[5:length(predictors_filtered)]
    model_g <- run_bootstrap_iteration(
      predictors_filtered = predictors_g,
      df = df,
      b = b,
      coef_max = coef_max,
      min_concord = min_concord,
      clinical_cleaned_risk_test = clinical_cleaned_risk_test,
      selected_gene_df = selected_gene_df,
      perf_list = perf_list_mg,
      rsf_predictor_list = rsf_predictor_list_mg,
      best_perf = best_perf_mg,
      best_model = best_model_mg,
      best_iter = best_iter_mg,
      best_seed = best_seed_mg,
      best_method = best_method_mg,
      current_seed = current_seed
    )
    
    # Update all variables
    perf_list_mg <- model_g$perf_list
    rsf_predictor_list_mg <- model_g$rsf_predictor_list
    best_perf_mg <- model_g$best_perf
    best_model_mg <- model_g$best_model
    best_iter_mg <- model_g$best_iter
    best_seed_mg <- model_g$best_seed
    best_method_mg <- model_g$best_method
    
    
    #model clinic
    predictors_c <- predictors_filtered[1:4]
    model_c <- run_bootstrap_iteration(
      predictors_filtered = predictors_c,  # Fixed: should be predictors_filtered_c not predictors_filtered_g
      df = df,
      b = b,
      coef_max = coef_max,
      min_concord = min_concord,
      clinical_cleaned_risk_test = clinical_cleaned_risk_test,
      selected_gene_df = selected_gene_df,
      perf_list = perf_list_mc,
      rsf_predictor_list = rsf_predictor_list_mc,
      best_perf = best_perf_mc,
      best_model = best_model_mc,
      best_iter = best_iter_mc,
      best_seed = best_seed_mc,
      best_method = best_method_mc,
      current_seed = current_seed
    )
    
    # Update all variables
    perf_list_mc <- model_c$perf_list
    rsf_predictor_list_mc <- model_c$rsf_predictor_list
    best_perf_mc <- model_c$best_perf
    best_model_mc <- model_c$best_model
    best_iter_mc <- model_c$best_iter
    best_seed_mc <- model_c$best_seed
    best_method_mc <- model_c$best_method
#####    
    # #model multimodal
    # model_m <- run_bootstrap_iteration(
    #   predictors_filtered = predictors_filtered,
    #   df = df,
    #   b = b,
    #   coef_max = coef_max,
    #   min_concord = min_concord,
    #   clinical_cleaned_risk_test = clinical_cleaned_risk_test,
    #   selected_gene_df = selected_gene_df,
    #   perf_list = perf_list_mm,
    #   gene_list = gene_list_mm,
    #   rsf_predictor_list = rsf_predictor_list_mm,
    #   best_perf = best_perf_mm,
    #   best_model = best_model_mm,
    #   best_genes = best_genes_mm,
    #   best_iter = best_iter_mm,
    #   best_seed = best_seed_mm,
    #   best_method = best_method_mm,
    #   current_seed = current_seed
    # )
    # 
    # # 更新所有变量
    # perf_list_mm <- model_m$perf_list
    # gene_list_mm <- model_m$gene_list
    # rsf_predictor_list_mm <- model_m$rsf_predictor_list
    # best_perf_mm <- model_m$best_perf
    # best_model_mm <- model_m$best_model
    # best_genes_mm <- model_m$best_genes
    # best_iter_mm <- model_m$best_iter
    # best_seed_mm <- model_m$best_seed
    # best_method_mm <- model_m$best_method
    # 
    # 
    # #model gene
    # predictors_filtered_g <- predictors_filtered[5:length(predictors_filtered)]
    # model_g <- run_bootstrap_iteration(
    #   predictors_filtered = predictors_filtered_g,
    #   df = df,
    #   b = b,
    #   coef_max = coef_max,
    #   min_concord = min_concord,
    #   clinical_cleaned_risk_test = clinical_cleaned_risk_test,
    #   selected_gene_df = selected_gene_df,
    #   perf_list = perf_list_mg,
    #   gene_list = gene_list_mg,
    #   rsf_predictor_list = rsf_predictor_list_mg,
    #   best_perf = best_perf_mg,
    #   best_model = best_model_mg,
    #   best_genes = best_genes_mg,
    #   best_iter = best_iter_mg,
    #   best_seed = best_seed_mg,
    #   best_method = best_method_mg,
    #   current_seed = current_seed
    # )
    # 
    # # 更新所有变量
    # perf_list_mg <- model_g$perf_list
    # gene_list_mg <- model_g$gene_list
    # rsf_predictor_list_mg <- model_g$rsf_predictor_list
    # best_perf_mg <- model_g$best_perf
    # best_model_mg <- model_g$best_model
    # best_genes_mg <- model_g$best_genes
    # best_iter_mg <- model_g$best_iter
    # best_seed_mg <- model_g$best_seed
    # best_method_mg <- model_g$best_method
    # 
    # 
    # #model clinic
    # predictors_filtered_c <- predictors_filtered[1:4]
    # model_c <- run_bootstrap_iteration(
    #   predictors_filtered = predictors_filtered_g,
    #   df = df,
    #   b = b,
    #   coef_max = coef_max,
    #   min_concord = min_concord,
    #   clinical_cleaned_risk_test = clinical_cleaned_risk_test,
    #   selected_gene_df = selected_gene_df,
    #   perf_list = perf_list_mc,
    #   gene_list = gene_list_mc,
    #   rsf_predictor_list = rsf_predictor_list_mc,
    #   best_perf = best_perf_mc,
    #   best_model = best_model_mc,
    #   best_genes = best_genes_mc,
    #   best_iter = best_iter_mc,
    #   best_seed = best_seed_mc,
    #   best_method = best_method_mc,
    #   current_seed = current_seed
    # )
    # 
    # # 更新所有变量
    # perf_list_mc <- model_c$perf_list
    # gene_list_mc <- model_c$gene_list
    # rsf_predictor_list_mc <- model_c$rsf_predictor_list
    # best_perf_mc <- model_c$best_perf
    # best_model_mc <- model_c$best_model
    # best_genes_mc <- model_c$best_genes
    # best_iter_mc <- model_c$best_iter
    # best_seed_mc <- model_c$best_seed
    # best_method_mc <- model_c$best_method
    # 
    
    
    
  #####  
    # # 拟合 Cox 模型
    # results_train <- fit_cox_model(predictors_filtered, df)
    # message(sprintf("Bootstrap %d: Cox 模型拟合完成，model 是否存在 = %s",
    #                 b, !is.null(results_train$model)))
    # 
    # if (is.null(results_train$model)) {
    #   message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
    #   next
    # }
    # 
    # # 发散检测
    # if (any(abs(coef(results_train$model)) > coef_max)) {
    #   message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
    #   next
    # }
    # message("Cox 模型发散检测完成")
    # 
    # # 检查完全分离
    # concordance_val <- tryCatch({
    #   suppressWarnings(summary(results_train$model)$concordance[1])
    # }, error = function(e) NA)
    # message(sprintf("Bootstrap %d: 完全分离检查完成，concordance_val = %.4f",
    #                 b, if (is.na(concordance_val)) NA else concordance_val))
    # 
    # if (!is.na(concordance_val) && concordance_val >= min_concord) {
    #   message(sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
    #   next
    # }
    # 
    # 
    # 
    # # Cox 模型评估
    # result_valid <- calculate_time_auc_cindex(
    #   "Cox",
    #   fitted_model = results_train$model,
    #   df = clinical_cleaned_risk_test
    # )
    # message(sprintf("Bootstrap %d: Cox 模型评估完成，iAUC = %.4f, c_index = %.4f",
    #                 b, result_valid$iAUC, result_valid$c_index))
    # 
    # # Random Survival Forest
    # clinical_rsf <- df
    # result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000)
    # rsf_fit_best <- result_rsf_train$best_model
    # message("RSF 模型拟合完成")
    # 
    # result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
    # message(sprintf("Bootstrap %d: RSF 模型评估完成，iAUC = %.4f, c_index = %.4f",
    #                 b, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    # 
    # # 保存性能和基因
    # perf_list[[b]] <- list(
    #   cox = result_valid,
    #   rsf = result_rsf_valid
    # )
    # gene_list[[b]] <- selected_gene_df$gene
    # rsf_predictor_list[[b]] <- predictors_filtered
    # message("性能和基因保存完成")
    # 
    # # 计算 Cox 和 RSF 综合评分
    # cox_score <- if (!is.na(result_valid$iAUC) && !is.na(result_valid$c_index)) {
    #   result_valid$iAUC + result_valid$c_index
    # } else {
    #   -Inf
    # }
    # rsf_score <- if (!is.na(result_rsf_valid$iAUC) && !is.na(result_rsf_valid$c_index)) {
    #   result_rsf_valid$iAUC + result_rsf_valid$c_index
    # } else {
    #   -Inf
    # }
    # cat(sprintf("Cox_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
    #             b, cox_score, result_valid$iAUC, result_valid$c_index))
    # cat(sprintf("rsf_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
    #             b, rsf_score, result_rsf_valid$iAUC, result_rsf_valid$c_index))
    # 
    # # 判断是否更优（Cox）
    # if (cox_score > best_perf) {
    #   best_perf    <- cox_score
    #   best_model   <- results_train$model
    #   best_genes   <- selected_gene_df$gene
    #   best_iter    <- b
    #   best_seed    <- current_seed
    #   best_method  <- "Cox"
    # }
    # 
    # # 判断是否更优（RSF）
    # if (rsf_score > best_perf) {
    #   best_perf    <- rsf_score
    #   best_model   <- rsf_fit_best
    #   best_genes   <- predictors_filtered
    #   best_iter    <- b
    #   best_seed    <- current_seed
    #   best_method  <- "RSF"
    # }
    # 
    # cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
    # cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
    # message(sprintf("Bootstrap %d: 迭代完成", b))
  }
  
  message("所有 Bootstrap 迭代完成")
  
  ##### # 汇总
#####
    # perf_list_nz <- Filter(Negate(is.null), perf_list) # fill out null
  # 
  # cox_iAUC <- sapply(perf_list_nz, function(x)
  #   if (!is.null(x$cox) && !is.null(x$cox$iAUC)) x$cox$iAUC else NA_real_)
  # cox_cidx <- sapply(perf_list_nz, function(x)
  #   if (!is.null(x$cox) && !is.null(x$cox$c_index)) x$cox$c_index else NA_real_)
  # 
  # rsf_iAUC <- sapply(perf_list_nz, function(x)
  #   if (!is.null(x$rsf) && !is.null(x$rsf$iAUC)) x$rsf$iAUC else NA_real_)
  # rsf_cidx <- sapply(perf_list_nz, function(x)
  #   if (!is.null(x$rsf) && !is.null(x$rsf$c_index)) x$rsf$c_index else NA_real_)
  # 
  # mean_cox_iAUC <- mean(cox_iAUC, na.rm = TRUE)
  # mean_cox_cidx <- mean(cox_cidx, na.rm = TRUE)
  # mean_rsf_iAUC <- mean(rsf_iAUC, na.rm = TRUE)
  # mean_rsf_cidx <- mean(rsf_cidx, na.rm = TRUE)
  # 
  # message("性能指标汇总完成：mean_cox_iAUC = ", mean_cox_iAUC, 
  #         ", mean_cox_cidx = ", mean_cox_cidx,
  #         ", mean_rsf_iAUC = ", mean_rsf_iAUC,
  #         ", mean_rsf_cidx = ", mean_rsf_cidx
  #         )
  # 
  # all_genes <- unlist(gene_list)
  # gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  # message("基因频率计算完成")
  # 
#####  
  # summary_results_mm <- summarize_performance(
  #   perf_list = perf_list_mm,
  #   gene_list = gene_list_mm,
  #   B = B
  # )
  # summary_results_mg <- summarize_performance(
  #   perf_list = perf_list_mg,
  #   gene_list = gene_list_mg,
  #   B = B
  # )
  # 
  # summary_results_mc <- summarize_performance(
  #   perf_list = perf_list_mc,
  #   gene_list = gene_list_mc,
  #   B = B
  # )
#####
  
  # 创建数据框
  summary_df <- data.frame(
    model = c("mc", "mg", "mm"),
    rsf_predictors = I(list(
      rsf_predictor_list_mc[[1]],
      rsf_predictor_list_mg[[1]],
      rsf_predictor_list_mm[[1]]
    )),
    best_performance = c(best_perf_mc, best_perf_mg, best_perf_mm),
    best_iteration = c(best_iter_mc, best_iter_mg, best_iter_mm),
    best_seed = c(best_seed_mc, best_seed_mg, best_seed_mm),
    best_method = c(best_method_mc, best_method_mg, best_method_mm),
    stringsAsFactors = FALSE
  )
  
  summary_results_mm <- summarize_performance(
    perf_list = perf_list_mm,
    B = B
  )
  summary_results_mg <- summarize_performance(
    perf_list = perf_list_mg,
    B = B
  )
  
  summary_results_mc <- summarize_performance(
    perf_list = perf_list_mc,
    B = B
  )
  
  # 计算基因频率
  all_genes <- unlist(gene_list)
  gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  message("基因频率计算完成")
  
  # 计算km_p的平均值
  km_p_list_nz <- Filter(Negate(is.null), km_p_list) 
  km_p_values <- unlist(km_p_list_nz)
  mean_km_p <- mean(km_p_values, na.rm = TRUE)
  message("性能指标汇总完成: mean_km_p = ", mean_km_p)
  
  return(list(
    summary_df = summary_df,
    best_model_mc = best_model_mc,
    best_model_mg = best_model_mg,
    best_model_mm = best_model_mm,
    train_indices_list = train_indices_list,
    mean_km_p = mean_km_p,  
    km_p_values = km_p_list_nz, 
    gene_frequency = gene_freq,
    summary_results_mm = summary_results_mm,
    summary_results_mg = summary_results_mg,
    summary_results_mc = summary_results_mc
    
    
    # mean_cox_iAUC = mean_cox_iAUC,
    # mean_cox_cidx = mean_cox_cidx,
    # mean_rsf_cidx = mean_rsf_cidx,
    # mean_rsf_iAUC = mean_rsf_iAUC,
    # 
    # all_results = perf_list,
    # best_model = best_model,
    # best_perf = best_perf,
    # best_predictors = best_genes,
    # best_iter = best_iter,
    # best_seed = best_seed,
    # best_method = best_method,
    
  ))
}

# run_bootstrap_validation_3_model(expr_mat, clinical_df, 
#                                              B = 10, 
#                                              seed = 90000,
#                                              min_epv = 2.5, 
#                                              coef_max = 10,
#                                              min_concord = 0.98)




res7390_3_model3_2 <- run_bootstrap_validation_3_model(expr_mat = expr_mat_7390,
                                           clinical_df = clinical_df_7390, 
                                           B = 100, 
                                           seed = 82000000, 
                                           min_epv = 4, 
                                           coef_max = 10, 
                                           min_concord = 0.97)


saveRDS(res7390_3_model3, file = "res7390_3_model3.rds")


res39582$best_model   # 最优 Cox 模型对象
res39582$best_perf    # 最优 iAUC
res39582$best_genes   # 对应的基因列表
res39582$mean_iAUC

res39582$mean_cindex
res39582$gene_frequency

res39582$best_predictors
res39582$best_iter
res39582_c$
  
  
  # risk factor model
  
  run_bootstrap_validation_safe_riskscore <- function(expr_mat, clinical_df, 
                                                      B = 20, 
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
      
      columns_to_remove <- c("e_dmfs", "t_dmfs", "geo_accession", "hospital")
      selected_col <- setdiff(colnames(clinical_df), columns_to_remove)
      missed_col <- colnames(clinical_df)[colSums(is.na(clinical_df)) > 0]
      
      res <- impute_fit_apply(
        train_df     = train_clinical,
        test_df      = test_clinical,
        selected_col = selected_col,
        missed_col   = missed_col
      )
      
      train_clinical <- res$train
      test_clinical  <- res$test
      
      
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
                                                               scale_cols = numeric_columns)
      message("标准化训练集临床数据完成：train_clinical_scaled 维度 = ", 
              toString(dim(train_clinical_scaled)))
      
      # 标准化测试集临床数据
      test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                              test_clinical,
                                                              scale_cols = numeric_columns)
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
      
      # predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
      #                  "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")
      # predictors <- c(predictors0, "risk_score")
      
      columns_to_remove2 <- c("e_dmfs", "t_dmfs", "risk_group")
      
      character_columns <- clinical_cleaned_risk_train %>%
        select(where(is.character)) %>%
        colnames()
      
      columns_with_digits <- clinical_cleaned_risk_train %>%
        select(matches("\\d")) %>%
        names()
      
      predictors <- colnames(clinical_cleaned_risk_train) %>%
        setdiff(c(columns_to_remove2, character_columns, columns_with_digits))
      
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


res7390_riskscore <- run_bootstrap_validation_safe_riskscore(expr_mat_7390, clinical_df_7390, 
                                                              B = 20, 
                                                              seed = 9000010, 
                                                              min_epv = 4, 
                                                              coef_max = 10, 
                                                              min_concord = 0.97)

saveRDS(res7390_riskscore, "res7390_riskscore_1_35.rds") 

#res39582_riskscore <- readRDS("res39582_riskscore_1_562593.rds")

# updated with inside imputation  
####2025-8-24
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
    
    columns_to_remove <- c("e_dmfs", "t_dmfs", "geo_accession", "hospital")
    selected_col <- setdiff(colnames(clinical_df), columns_to_remove)
    missed_col <- colnames(clinical_df)[colSums(is.na(clinical_df)) > 0]
    
    res <- impute_fit_apply(
      train_df     = train_clinical,
      test_df      = test_clinical,
      selected_col = selected_col,
      missed_col   = missed_col
    )
    
    train_clinical <- res$train
    test_clinical  <- res$test
    
    
    
    message("Bootstrap 采样完成：", #train_expr 维度 = ", toString(dim(train_expr)), 
            #", test_expr 维度 = ", toString(dim(test_expr)),
            ", train_clinical 维度 = ", toString(dim(train_clinical)),
            ", test_clinical 维度 = ", toString(dim(test_clinical)))
    
    
    
    # 事件数
    events_train <- sum(train_clinical$e_dmfs)
    message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
    message("事件数计算完成")
    
    # 
    # predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
    #                  "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")
    
    columns_to_remove2 <- c("e_dmfs", "t_dmfs", "risk_group", "risk_score")
    
    character_columns <- clinical_cleaned_risk_train %>%
      select(where(is.character)) %>%
      colnames()
    
    columns_with_digits <- clinical_cleaned_risk_train %>%
      select(matches("\\d")) %>%
      names()
    
    predictors <- colnames(clinical_cleaned_risk_train) %>%
      setdiff(c(columns_to_remove2, character_columns, columns_with_digits))
    
    message(sprintf("Bootstrap %d: 预测变量选择完成，predictors 长度 = %d", 
                    b, length(predictors)))
    
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
    df <- train_clinical #%>%
    #select(-geo_accession )
    # message("预测变量相关性过滤完成：filtered_data 维度 = ", toString(dim(filtered_data)))
    
    # 拟合 Cox 模型
    results_train <- fit_cox_model(predictors, df)
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
    clinical_rsf <- df %>%
      select(-any_of("geo_accession"))
    
    result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000)
    rsf_fit_best <- result_rsf_train$best_model
    message("RSF 模型拟合完成")
    
    result_rsf_valid <- calculate_time_auc_cindex(
      "RSF",
      fitted_model = rsf_fit_best,
      df = test_clinical)
    
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



res7390_clinical_2025_9_03 <- clinical_run_bootstrap_validation(clinical_df_7390, 
                                                                 B = 20, 
                                                                 seed = 900000, 
                                                                 # min_epv = 4, 
                                                                 coef_max = 10, 
                                                                 min_concord = 0.98)

saveRDS(res7390_clinical_2025_9_03, file = "res7390_clinical_2025_9_03_1.3703.rds")
res39582_clinical_2025_8_30 <- readRDS("res39582_clinical_2025_8_30_1_5999.rds")
#################################################################################
# 
# gen_run_bootstrap_validation_safe <- function(expr_mat, clinical_df, 
#                                               B = 20, 
#                                               seed = 90000,
#                                               min_epv = 4, 
#                                               coef_max = 10,
#                                               min_concord = 0.98) {
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
#   clinical_df <- clinical_df %>%
#     mutate(geo_accession = rownames(clinical_df)) %>%
#     select( "e_dmfs", "t_dmfs", "geo_accession")
#   
#   message("初始化完成：expr_mat 维度 = ", toString(dim(expr_mat)), 
#           ", clinical_df 维度 = ", toString(dim(clinical_df)))
#   
#   for (b in 1:B) {
#     # Bootstrap 采样
#     current_seed <- seed + b 
#     set.seed(current_seed) 
#     message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
#     
#     train_idx <- sample(seq_len(n), size = n, replace = TRUE)
#     test_idx  <- setdiff(seq_len(n), unique(train_idx))  # OOB 样本
#     message(sprintf("Bootstrap %d: train_idx length = %d, test_idx length = %d", 
#                     b, length(train_idx), length(test_idx)))
#     
#     if (length(test_idx) == 0) {
#       message(sprintf("Bootstrap %d skipped: No OOB samples", b))
#       next
#     }
#     
#     train_expr     <- expr_mat[, train_idx, drop = FALSE]
#     test_expr      <- expr_mat[, test_idx, drop = FALSE]
#     train_clinical <- clinical_df[train_idx, , drop = FALSE]
#     test_clinical  <- clinical_df[test_idx, , drop = FALSE]
#     
#     message("Bootstrap 采样完成：train_expr 维度 = ", toString(dim(train_expr)), 
#             ", test_expr 维度 = ", toString(dim(test_expr)),
#             ", train_clinical 维度 = ", toString(dim(train_clinical)),
#             ", test_clinical 维度 = ", toString(dim(test_clinical)))
#     
#     if (ncol(train_expr) == 0 || nrow(train_expr) == 0) {
#       message(sprintf("Bootstrap %d skipped: Invalid train_expr dimensions", b))
#       next
#     }
#     
#     # 事件数
#     events_train <- sum(train_clinical$e_dmfs)
#     message(sprintf("Bootstrap %d: Number of events in training set = %d", b, events_train))
#     message("事件数计算完成")
#     
#     # MAD 过滤
#     mad_train_expr <- apply(train_expr, 1, mad)
#     cutoff <- quantile(mad_train_expr, 0.25)
#     keep_mad <- mad_train_expr > cutoff
#     message(sprintf("Bootstrap %d: MAD 过滤完成，保留基因数 = %d, cutoff = %.4f", 
#                     b, sum(keep_mad), cutoff))
#     
#     if (sum(keep_mad) == 0) {
#       message(sprintf("Bootstrap %d skipped: No genes passed MAD filter", b))
#       next
#     }
#     
#     train_expr2 <- train_expr[keep_mad, , drop = FALSE]
#     test_expr2 <- test_expr[keep_mad, , drop = FALSE]
#     message("MAD 过滤后矩阵生成完成：train_expr2 维度 = ", toString(dim(train_expr2)), 
#             ", test_expr2 维度 = ", toString(dim(test_expr2)))
#     
#     # 单变量 Cox 回归
#     sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = 0.01)
#     message(sprintf("Bootstrap %d: 单变量 Cox 回归完成，显著基因数 = %d", 
#                     b, if (is.null(sig_gene_df)) 0 else nrow(sig_gene_df)))
#     
#     if (is.null(sig_gene_df) || nrow(sig_gene_df) == 0) {
#       message(sprintf("Bootstrap %d skipped: No significant genes from univariate Cox", b))
#       next
#     }
#     
#     significant_gene <- sig_gene_df$gene
#     message(sprintf("Bootstrap %d: 提取显著基因完成，significant_gene 长度 = %d", 
#                     b, length(significant_gene)))
#     
#     
#     # 标准化训练集基因数据
#     train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
#                                                 gene_mat_valid = train_expr2,
#                                                 significant_gene = significant_gene)
#     message("标准化训练集基因数据完成：train_expr_scaled 维度 = ", 
#             toString(dim(train_expr_scaled)))
#     
#     if (is.null(train_expr_scaled) || nrow(train_expr_scaled) == 0) {
#       message(sprintf("Bootstrap %d skipped: Invalid standardized training data", b))
#       next
#     }
#     
#     
#     # 相关性过滤
#     train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = 0.90)
#     message("相关性过滤完成：train_expr_filtered 维度 = ", 
#             toString(dim(train_expr_filtered)))
#     
#     
#     significant_gene2 <- rownames(train_expr_filtered)
#     message(sprintf("Bootstrap %d: 提取过滤后基因完成，significant_gene2 长度 = %d", 
#                     b, length(significant_gene2)))
#     
#     # 标准化测试集基因数据
#     test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
#                                                gene_mat_valid = test_expr2,
#                                                significant_gene = significant_gene2)
#     message("标准化测试集基因数据完成：test_expr_scaled 维度 = ", 
#             toString(dim(test_expr_scaled)))
#     
#     # # 标准化训练集临床数据
#     # train_clinical_scaled <- standardize_with_train_clinical(train_clinical,
#     #                                                          train_clinical,
#     #                                                          scale_cols = c("Age_at_Diagnosis"))
#     # message("标准化训练集临床数据完成：train_clinical_scaled 维度 = ", 
#     #         toString(dim(train_clinical_scaled)))
#     # 
#     # # 标准化测试集临床数据
#     # test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
#     #                                                         test_clinical,
#     #                                                         scale_cols = c("Age_at_Diagnosis"))
#     # message("标准化测试集临床数据完成：test_clinical_scaled 维度 = ", 
#     #         toString(dim(test_clinical_scaled)))
#     # 
#     # LASSO Cox
#     gene_freq_df <- repeat_cv_lasso_cox(train_expr = train_expr_filtered,
#                                         train_clinical,
#                                         significant_gene_vec = significant_gene2,
#                                         repeats = 5,
#                                         nfolds = 10,
#                                         alpha = 1)
#     message(sprintf("Bootstrap %d: LASSO Cox 完成，gene_freq_df 行数 = %d", 
#                     b, nrow(gene_freq_df)))
#     
#     # 动态筛选变量，满足 EPV 要求
#     gene_freq_df_best <- gene_freq_df %>%
#       filter(freq >= 0.8)
#     max_vars_allowed <- floor(events_train / min_epv)
#     message(sprintf("Bootstrap %d: 动态筛选变量完成，gene_freq_df_best 行数 = %d, max_vars_allowed = %d", 
#                     b, nrow(gene_freq_df_best), max_vars_allowed))
#     
#     if (max_vars_allowed == 0 || nrow(gene_freq_df) == 0) {
#       message(sprintf("Bootstrap %d skipped: Too few events or no selected genes", b))
#       next
#     }
#     
#     selected_gene_df <- gene_freq_df[1:min(nrow(gene_freq_df_best), max_vars_allowed), ] %>%
#       mutate(coef = mean_coef)
#     message(sprintf("Bootstrap %d: 选择基因完成，selected_gene_df 行数 = %d, EPV = %.2f", 
#                     b, nrow(selected_gene_df), events_train / nrow(selected_gene_df)))
#     
#     # 计算风险分数
#     clinical_cleaned_risk_train <- compute_risk_score(
#       gene_mat_scaled = train_expr_filtered,
#       significant_vars_df = selected_gene_df,
#       clinical_cleaned = train_clinical,
#       n_group = 3
#     )
#     # message("训练集风险分数计算完成：clinical_cleaned_risk_train 维度 = ", 
#     #         toString(dim(clinical_cleaned_risk_train)))
#     # 
#     # predictors0 <- c("Sex", "Age_at_Diagnosis", "TNM_T", "TNM_N", "TNM_M", "Tumor_Location",
#     #                  "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")
#     
#     #predictors <- c(predictors0, colnames(clinical_cleaned_risk_train)[15:ncol(clinical_cleaned_risk_train)])
#     predictors <- colnames(clinical_cleaned_risk_train)[6:ncol(clinical_cleaned_risk_train)]
#     message(sprintf("Bootstrap %d: 预测变量选择完成，predictors 长度 = %d", 
#                     b, length(predictors)))
#     
#     
#     # 相关性过滤
#     #predictor_data <- clinical_cleaned_risk_train[, predictors]
#     #numeric_vars <- predictor_data[, sapply(predictor_data, is.numeric), drop = FALSE]
#     # filtered_numeric_vars <- remove_high_corr(numeric_vars, threshold = 0.9)
#     # filtered_data <- cbind(
#     #   predictor_data[, !sapply(predictor_data, is.numeric), drop = FALSE],
#     #   filtered_numeric_vars
#     # )
#     # predictors_filtered <- colnames(filtered_data)
#     
#     
#     
#     # message("预测变量相关性过滤完成：filtered_data 维度 = ", toString(dim(filtered_data)))
#     
#     # 拟合 Cox 模型
#     results_train <- fit_cox_model(predictors, clinical_cleaned_risk_train)
#     message(sprintf("Bootstrap %d: Cox 模型拟合完成，model 是否存在 = %s", 
#                     b, !is.null(results_train$model)))
#     
#     if (is.null(results_train$model)) {
#       message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
#       next
#     }
#     
#     # 发散检测
#     if (any(abs(coef(results_train$model)) > coef_max)) {
#       message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
#       next
#     }
#     message("Cox 模型发散检测完成")
#     
#     # 检查完全分离
#     concordance_val <- tryCatch({
#       suppressWarnings(summary(results_train$model)$concordance[1])
#     }, error = function(e) NA)
#     message(sprintf("Bootstrap %d: 完全分离检查完成，concordance_val = %.4f", 
#                     b, if (is.na(concordance_val)) NA else concordance_val))
#     
#     if (!is.na(concordance_val) && concordance_val >= min_concord) {
#       message(sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
#       next
#     }
#     
#     # 测试集风险分数
#     clinical_cleaned_risk_test <- compute_risk_score(
#       gene_mat_scaled = test_expr_scaled,
#       significant_vars_df = selected_gene_df,
#       clinical_cleaned = test_clinical,
#       n_group = 3
#     )
#     message("测试集风险分数计算完成：clinical_cleaned_risk_test 维度 = ", 
#             toString(dim(clinical_cleaned_risk_test)))
#     
#     # Cox 模型评估
#     result_valid <- calculate_time_auc_cindex(
#       "Cox", 
#       fitted_model = results_train$model, 
#       df = clinical_cleaned_risk_test
#     )
#     message(sprintf("Bootstrap %d: Cox 模型评估完成，iAUC = %.4f, c_index = %.4f", 
#                     b, result_valid$iAUC, result_valid$c_index))
#     
#     # Random Survival Forest
#     
#     clinical_rsf <- clinical_cleaned_risk_train %>%
#       select(-c("geo_accession", "risk_score","risk_group" ))
#     
#     result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000)
#     rsf_fit_best <- result_rsf_train$best_model
#     message("RSF 模型拟合完成")
#     
#     result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
#     message(sprintf("Bootstrap %d: RSF 模型评估完成，iAUC = %.4f, c_index = %.4f", 
#                     b, result_rsf_valid$iAUC, result_rsf_valid$c_index))
#     
#     # 保存性能和基因
#     perf_list[[b]] <- list(
#       cox = result_valid,
#       rsf = result_rsf_valid
#     )
#     gene_list[[b]] <- selected_gene_df$gene
#     rsf_predictor_list[[b]] <- predictors_filtered
#     message("性能和基因保存完成")
#     
#     # 计算 Cox 和 RSF 综合评分
#     cox_score <- if (!is.na(result_valid$iAUC) && !is.na(result_valid$c_index)) {
#       result_valid$iAUC + result_valid$c_index
#     } else {
#       -Inf
#     }
#     rsf_score <- if (!is.na(result_rsf_valid$iAUC) && !is.na(result_rsf_valid$c_index)) {
#       result_rsf_valid$iAUC + result_rsf_valid$c_index
#     } else {
#       -Inf
#     }
#     cat(sprintf("Cox_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
#                 b, cox_score, result_valid$iAUC, result_valid$c_index))
#     cat(sprintf("rsf_score in present iteration %d is %f, in which valid_iAUC is %f and valid_c_index is %f\n",
#                 b, rsf_score, result_rsf_valid$iAUC, result_rsf_valid$c_index))
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
#       best_genes   <- predictors_filtered
#       best_iter    <- b
#       best_seed    <- current_seed
#       best_method  <- "RSF"
#     }
#     
#     cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
#     cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
#     message(sprintf("Bootstrap %d: 迭代完成", b))
#   }
#   
#   message("所有 Bootstrap 迭代完成")
#   
#   # 汇总
#   perf_list_nz <- Filter(Negate(is.null), perf_list)
#   
#   cox_iAUC <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$cox) && !is.null(x$cox$iAUC)) x$cox$iAUC else NA_real_)
#   cox_cidx <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$cox) && !is.null(x$cox$c_index)) x$cox$c_index else NA_real_)
#   
#   rsf_iAUC <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$rsf) && !is.null(x$rsf$iAUC)) x$rsf$iAUC else NA_real_)
#   rsf_cidx <- sapply(perf_list_nz, function(x)
#     if (!is.null(x$rsf) && !is.null(x$rsf$c_index)) x$rsf$c_index else NA_real_)
#   
#   mean_cox_iAUC <- mean(cox_iAUC, na.rm = TRUE)
#   mean_cox_cidx <- mean(cox_cidx, na.rm = TRUE)
#   mean_rsf_iAUC <- mean(rsf_iAUC, na.rm = TRUE)
#   mean_rsf_cidx <- mean(rsf_cidx, na.rm = TRUE)
#   
#   message("性能指标汇总完成：mean_cox_iAUC = ", mean_cox_iAUC, 
#           ", mean_cox_cidx = ", mean_cox_cidx,
#           ", mean_rsf_iAUC = ", mean_rsf_iAUC,
#           ", mean_rsf_cidx = ", mean_rsf_cidx)
#   
#   all_genes <- unlist(gene_list)
#   gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
#   message("基因频率计算完成")
#   
#   return(list(
#     mean_cox_iAUC = mean_cox_iAUC,
#     mean_cox_cidx = mean_cox_cidx,
#     mean_rsf_cidx = mean_rsf_cidx,
#     mean_rsf_iAUC = mean_rsf_iAUC,
#     gene_frequency = gene_freq,
#     all_results = perf_list,
#     best_model = best_model,
#     best_perf = best_perf,
#     best_predictors = best_genes,
#     best_iter = best_iter,
#     best_seed = best_seed,
#     best_method = best_method,
#     train_indices_list = train_indices_list
#   ))
# }
# 
# 
# 
# res39582_pur_gen_2025_8_30 <- gen_run_bootstrap_validation_safe(expr_mat, clinical_df, 
#                                                                 B = 20, 
#                                                                 seed = 1000001, 
#                                                                 min_epv = 4, 
#                                                                 coef_max = 10, 
#                                                                 min_concord = 0.98)
# 
# saveRDS(res39582_pur_gen_2025_8_30, file = "res39582_pur_gen_2025_8_30_1_4964.rds")
