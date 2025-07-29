library(GEOquery)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(randomForestSRC)
library(pec)
library(riskRegression)
library(sva)
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


#### gene

################################################


identical(colnames(expr_7390), clinical_cleaned_7390_imputed$geo_accession)

# # 表达矩阵的样本名（列名）
# sample_ids_expr <- colnames(expr)
# 
# # 临床数据的 GEO ID（假设是 geo_accession 列）
# sample_ids_clin <- clinical_cleaned$geo_accession
# 
# # 找到两者共有的样本
# matched_ids <- intersect(sample_ids_expr, sample_ids_clin)
# 
# # 子集表达数据和临床数据
# expr_mat <- expr[, matched_ids]
# clin_matched <- clinical_cleaned[match(matched_ids, clinical_cleaned$geo_accession), ]
# 
# # 验证维度
# stopifnot(ncol(expr_mat) == nrow(clin_matched))
# 



corrected <- combat_three_batches(
  expr_A = expr_7390, clin_A = clinical_cleaned_7390_imputed,
  expr_B = expr_2990, clin_B = clinical_cleaned_2990_imputed,
  expr_C = NULL, clin_C = NULL,
  predictors <- c("grade", "size", "age", "er")
  #predictors = c("age", "grade", "er")
)

expr_7390_corrected <- corrected$A
expr_2990_corrected   <- corrected$B
expr_test_corrected  <- corrected$C










#2. 构建批量单变量 Cox 回归函数
#####
library(survival)

#keep <- rowSums(expr_mat != 0) >= 10  # 至少10个样本表达不为0
#expr_mat_filtered <- expr_mat[keep, ]


# 构建 Cox 分析函数
cox_results <- apply(expr_7390_corrected, 1, function(gene_expr) {
  df <- data.frame(
    expr = scale(as.numeric(gene_expr)),  # <--- 标准化这一步
    #expr = as.numeric(gene_expr),
    time = clinical_cleaned_7390_imputed$t_dmfs,
    status = clinical_cleaned_7390_imputed$e_dmfs
  )
  
  fit <- tryCatch(
    coxph(Surv(time, status) ~ expr, data = df),
    error = function(e) return(NULL)
  )
  
  if (is.null(fit)) return(c(NA, NA, NA, NA))
  
  s <- summary(fit)
  coef <- s$coefficients[1, "coef"]
  hr <- s$coefficients[1, "exp(coef)"]
  p <- s$coefficients[1, "Pr(>|z|)"]
  se <- s$coefficients[1, "se(coef)"]
  c(coef = coef, HR = hr, SE = se, p.value = p)
})

"218727_at"

#3. 整理结果表格并筛选显著基因
cox_df <- as.data.frame(t(cox_results))
cox_df$gene <- rownames(cox_df)

# 按 p 值排序
cox_df <- cox_df[order(cox_df$p.value), ]

# 筛选显著基因（比如 p < 0.05）
sig_genes_df <- subset(cox_df, p.value < 0.11)
sig_genes <-sig_genes_df$gene

# Step 2: 提取表达数据

expr_top_7390 <- t(expr_7390_corrected[sig_genes, ])  # 转置：样本 × 基因

expr_top_7390_scale <- scale(expr_top_7390) # scale() 函数在 R 中默认是对列进行操作的


# 
# # 将矩阵转为 data.frame，同时添加新的列
# expr_top_7390_with_survival <- data.frame(
#   time = clinical_cleaned_7390_imputed$t_dmfs,
#   status = clinical_cleaned_7390_imputed$e_dmfs,
#   expr_top_7390_scale
# )




# 1. 提取表达矩阵 (X) 和生存数据 (y)
X <- expr_top_7390_scale
y <- Surv(clinical_cleaned_7390_imputed$t_dmfs, clinical_cleaned_7390_imputed$e_dmfs)

# 2. Lasso Cox 回归 + 10折交叉验证
set.seed(123)
cvfit <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)

# 3. 查看最佳 lambda 值
best_lambda <- cvfit$lambda.min
cat("Best lambda:", best_lambda, "\n")
best_lambda1se <- cvfit$lambda.1se
plot(cvfit)



# 4. 提取非零系数的基因（Lasso选中的）
coef_min <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(coef_min)[which(coef_min != 0)]

cat("Selected genes:\n")
print(selected_genes)


# 可视化每个 lambda 的表现
plot(cvfit)
abline(v = log(cvfit$lambda.min), col = "red", lty = 2)   # 最小误差对应的 lambda
abline(v = log(cvfit$lambda.1se), col = "blue", lty = 2)  # 1-SE rule 对应的 lambda

# 查看所有 lambda 值和交叉验证误差
# cv_table <- data.frame(
#   lambda = cvfit$lambda,
#   cvm = cvfit$cvm,              # mean cross-validated deviance
#   cvsd = cvfit$cvsd,            # standard deviation
#   nzero = cvfit$nzero           # 每个lambda下的非零系数数
# )

print(head(cv_table))

cvfit$lambda.min
#[1] 0.09631261
#[1] 0.105703
#0.002
# 最佳 lambda 对应的非零系数
best_coef <- coef(cvfit, s = "lambda.min")


coef_min <- coef(cvfit, s = "lambda.min")

# 把稀疏矩阵转成普通 matrix，然后筛选非零
selected_genes <- as.matrix(coef_min)

# 只保留非零系数的基因
selected_df0 <- data.frame(
  gene = rownames(selected_genes),
  coef = selected_genes[, 1]
)
selected_df0 <- selected_df0[selected_df0$coef != 0, ]

## 7390data
significant_gene_7390 <- sub("^ge_", "", rownames(selected_df0))

expr_7390_scaled <- standardize_with_train(gene_mat_train = expr_7390_corrected,
                                           gene_mat_valid=expr_7390_corrected,
                                           significant_gene = significant_gene_7390)

clinical_cleaned_risk_7390 <- compute_risk_score(gene_mat_scaled = expr_7390_scaled,
                                                 significant_vars_df = selected_df0,
                                                 clinical_cleaned= clinical_cleaned_7390_imputed,
                                                 n_group = 4)

predictors <- c("grade", "er", "age", "size", "risk_score")

results_7390 <- fit_cox_model(predictors, clinical_cleaned_risk_7390)

predictors <- c("size", "risk_score")
results_7390_0 <- fit_cox_model(predictors, clinical_cleaned_risk_7390)

predictors <- colnames(clinical_cleaned_risk_7390)[11: length(colnames(clinical_cleaned_risk_7390))]
results_7390 <- fit_cox_model(predictors, clinical_cleaned_risk_7390)
fitted_model_7390 = results_7390$model
calculate_time_auc_cindex(model_type = "Cox",
                          fitted_model_7390,
                          df = clinical_cleaned_risk_7390) # not possible because NA in data


calculate_time_auc_cindex(model_type = "Cox",
                          fitted_model_7390,
                          df = clinical_cleaned_risk_2990) # not possible because NA in data





## 2990
significant_gene_7390 <- sub("^ge_", "", rownames(selected_df0))

expr_2990_scaled <- standardize_with_train(gene_mat_train = expr_7390_corrected,
                                           gene_mat_valid=expr_2990_corrected,
                                           significant_gene = significant_gene_7390)

clinical_cleaned_risk_2990 <- compute_risk_score(gene_mat_scaled = expr_2990_scaled,
                                                 significant_vars = selected_df0,
                                                 clinical_cleaned= clinical_cleaned_2990_imputed,
                                                 n_group = 4)

#predictors <- c("grade", "er", "age", "size", "risk_score")
predictors <- colnames(clinical_cleaned_risk_2990)[11: length(colnames(clinical_cleaned_risk_2990))]
results_2990 <- fit_cox_model(predictors, clinical_cleaned_risk_2990)

plot_km_by_group(clinical_cleaned_risk_2990, "risk_group")
fitted_model_2990 = results_2990$model
calculate_time_auc_cindex(model_type = "Cox",
                          fitted_model_2990,
                          df = clinical_cleaned_risk_2990) # not possible because NA in data




















expr_top_7390_scale_reselect <- expr_top_7390_scale[,selected_df0$gene]
colnames(expr_top_7390_scale_reselect) <- paste0("ge_", colnames(expr_top_7390_scale_reselect))


expr_top_7390_scale_reselect_with_survival <- data.frame(
  t_dmfs = clinical_cleaned_7390_imputed$t_dmfs,
  e_dmfs = clinical_cleaned_7390_imputed$e_dmfs,
   expr_top_7390_scale_reselect
 )


predictors <- colnames(expr_top_7390_scale_reselect)
results <- fit_cox_model(predictors, expr_top_7390_scale_reselect_with_survival)


# 提取系数表
cox_summary <- summary(results$model)$coefficients

# 转为 data.frame，方便处理
cox_df <- as.data.frame(cox_summary)

# 过滤 p < 0.05 的行
significant_vars <- subset(cox_df, `Pr(>|z|)` < 0.05)

# 查看结果
significant_vars[, c("coef", "Pr(>|z|)")]



# 
# significant_gene <- rownames(significant_vars)
# significant_gene <- sub("^ge_", "", significant_gene)
# gene_mat_train_s <- gene_mat_train[significant_gene, ]
# gene_mat_train_s_scaled <- scale(t(gene_mat_train_s))
# 
# risk_score <- as.numeric(as.matrix(expr_top_7390_scale_reselect_15) %*% significant_vars$coef)
# risk_score <-  gene_mat_train_s_scaled %*% significant_vars$coef
# clinical_cleaned_7390_imputed$risk_score <- risk_score


# expr_7390_scaled_test <- standardize_with_train(gene_mat_train = expr_7390_corrected,
#                                            gene_mat_valid = expr_7390_corrected,
#                                            significant_gene
# )
# 
# clinical_cleaned_risk <- compute_risk_score(gene_mat_scaled = expr_7390_scaled_test, significant_vars, clinical_cleaned= clinical_cleaned_7390_imputed, n_group = 4)


#expr_top_7390_scale_reselect_15 <- expr_top_7390_scale_reselect[, rownames(significant_vars)]

risk_score <- as.numeric(as.matrix(expr_top_7390_scale_reselect_15) %*% significant_vars$coef)

identical(rownames(expr_top_7390_scale_reselect_15), clinical_cleaned_7390_imputed$geo_accession)

clinical_cleaned_7390_imputed$risk_score_old <- risk_score


predictors <- colnames(clinical_cleaned_7390_imputed)[c(4,5,6,7,9)]
results <- fit_cox_model(predictors, clinical_cleaned_7390_imputed)



library(dplyr)

clinical_cleaned_7390_imputed <- clinical_cleaned_7390_imputed %>%
  mutate(risk_group = ntile(risk_score, 3)) %>%
  mutate(risk_group = factor(risk_group, labels = c("low", "medium", "high")))



plot_km_by_group(clinical_cleaned_7390_imputed, "risk_group")


significant_gene <- rownames(significant_vars)




expr_2990_scaled <- standardize_with_train(gene_mat_train = expr_7390_corrected,
                                           gene_mat_valid=expr_2990_corrected,
                                           significant_gene
                                           )


clinical_cleaned_risk_2990 <- compute_risk_score(gene_mat_scaled = expr_2990_scaled,
                                                 significant_vars,
                                                 clinical_cleaned= clinical_cleaned_2990,
                                                 n_group = 4)



plot_km_by_group(clinical_cleaned_risk_2990, "risk_group")
predictors <- colnames(clinical_cleaned_risk_2990)[c(4,5,6,7,8,9,10)]
df=clinical_cleaned_risk_2990
fitted_result = fit_cox_model(predictors, df=clinical_cleaned_risk_2990)
fitted_model = fitted_result$model
calculate_time_auc_cindex(model_type = "Cox", fitted_model, df) # not possible because NA in data



df0=clinical_cleaned_risk_2990 %>%
  na.omit()
fitted_result0 = fit_cox_model(predictors, df=df0)
fitted_model0 = fitted_result0$model
calculate_time_auc_cindex(model_type = "Cox", fitted_model0, df0)

df1=clinical_cleaned_risk_2990 %>%
  na.omit()
predictors <- colnames(clinical_cleaned_risk_2990)[c(4,5,6,7,8,9)]
fitted_result1 = fit_cox_model(predictors, df=df0)
fitted_model1 = fitted_result1$model
calculate_time_auc_cindex(model_type = "Cox", fitted_model1, df0)


# 
# # Step 3: LASSO-Cox 建模
# # library(glmnet)
# # x <- as.matrix(expr_top)
# # y <- Surv(clin_matched$t_dmfs, clin_matched$e_dmfs)
# # cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1)
# 
# expr_top_df <- as.data.frame(expr_top) 
# expr_top_df$geo_accession <- rownames(expr_top_df)
# 
# interested_col <- c("geo_accession", "grade", "size", "age", "er") 
# 
# 
# clinic_gene_7390 <- clin_matched %>%
#   select(all_of(interested_col)) %>%
#   na.omit() %>%
#   left_join(
#     expr_top_df, by = "geo_accession"
#   )


