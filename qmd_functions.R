generate_box_plots <- function(data, continuous_variables) {
  # Generate box plots using base R graphics
  for (var in continuous_variables) {
    
    var_data <- na.omit(data[[var]])
    par(mar = c(3, 4, 3, 2))  # Set margins
    boxplot(var_data,
            col = "lightblue",
            main = paste("Boxplot of", var),
            ylab = var,
            cex.main = 0.8,
            cex.lab = 0.8,
            outcol = "red",
            outpch = 19,
            cex = 0.6)
    
    grid(nx = NA, ny = NULL, lwd = 0.2)
  }
}

impute_missing_value <- function(clinical_cleaned, selected_col, missed_col) {
  # Perform multiple imputation using mice package
  # 
  # @param clinical_cleaned Cleaned clinical data frame
  # @param selected_col Columns used in imputation model
  # @param missed_col Target column(s) to impute
  # @return Data frame with imputed missing values
  
  # Extract relevant variables
  df_miss <- clinical_cleaned[, selected_col]
  
  # Set imputation methods
  methods <- make.method(df_miss)
  for (col in missed_col) {
    if (col %in% c("grade", "TNM_T", "TNM_N")) {
      methods[col] <- "polr"
    } else if (col %in% c("er", "TNM_M", "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")) {
      methods[col] <- "logreg"
    }
  }
  
  # Perform imputation
  imp <- mice(df_miss, m = 20, method = methods, seed = 123)
  
  # Majority vote function
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Replace missing values
  for (col in missed_col) {
    tbl <- imp$imp[[col]]
    if (is.null(tbl) || nrow(tbl) == 0) next
    
    filled_values <- apply(tbl, 1, Mode)
    ids <- names(filled_values)
    idx <- match(ids, rownames(clinical_cleaned))
    keep <- !is.na(idx) & !is.na(filled_values)
    
    if (any(keep)) {
      vals <- filled_values[keep]
      pos <- idx[keep]
      clinical_cleaned[[col]][pos] <- vals
    }
  }
  
  return(clinical_cleaned)
}

impute_fit_apply <- function(train_df, test_df, selected_col, missed_col, m = 20) {
  # subset columns to impute
  train_sel <- train_df[, selected_col, drop = FALSE]
  test_sel  <- test_df[,  selected_col, drop = FALSE]
  
  # ensure rownames exist for binding/matching
  #if (is.null(rownames(train_sel))) rownames(train_sel) <- paste0("train_", seq_len(nrow(train_sel)))
  if (is.null(rownames(train_sel)) || !any(grepl("train", rownames(train_sel)))) {
    rownames(train_sel) <- paste0("train_", seq_len(nrow(train_sel)))
  }
  #if (is.null(rownames(test_sel)))  rownames(test_sel)  <- paste0("test_" , seq_len(nrow(test_sel)))
  if (is.null(rownames(test_sel)) || !any(grepl("test", rownames(test_sel)))) {
    rownames(test_sel) <- paste0("test_", seq_len(nrow(test_sel)))
  }
  df_all  <- rbind(train_sel, test_sel)
  is_test <- grepl("^test_", rownames(df_all))
  
  # build mice methods (only columns in missed_col get imputed)
  methods <- make.method(df_all)
  methods[] <- ""  # default: no imputation
  # methods <- methods[, !names(methods) %in% missed_col, drop = FALSE]
  # 
  # methods[intersect(missed_col, names(methods))] <- c()  # initialize
  # # common_cols <- intersect(missed_col, names(methods))
  # # if (length(common_cols) > 0) {
  # #   methods[common_cols] <- ""
  # # }
  
  
  for (col in missed_col) {
    if (!col %in% names(methods)) next
    methods[col] <- if (col %in% c("grade","TNM_T","TNM_N")) "polr"
    else if (col %in% c("er","TNM_M")) "logreg"
    else if (col %in% c("age","size")) "pmm" # predictive mean matching for numeric
    else ""  # skip others
  }
  
  # "TNM_M","Chemotherapy_Adjuvant","MMR_Status","KRAS_Mutation" exist in GSE39582 dataset
  
  # fit on train rows, apply to test rows
  imp <- mice(df_all, m = m, method = methods, ignore = is_test, printFlag = FALSE, seed = 123)
  
  # majority vote across m imputations
  Mode <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA)
    ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]
  }
  
  # write back single imputed value per missing cell (no new factor levels)
  df_all_imp <- df_all
  for (col in intersect(missed_col, names(imp$imp))) {
    tbl <- imp$imp[[col]]
    if (is.null(tbl) || nrow(tbl) == 0) next
    vals <- apply(tbl, 1, Mode)                       # named by rownames
    idx  <- match(names(vals), rownames(df_all_imp))   # map IDs to positions
    keep <- !is.na(idx) & !is.na(vals)
    if (!any(keep)) next
    
    pos  <- idx[keep]
    v    <- vals[keep]
    
    if (is.factor(df_all_imp[[col]])) {
      # strictly refuse unseen levels
      ok <- v %in% levels(df_all_imp[[col]])
      df_all_imp[[col]][pos[ok]] <- v[ok]
      # leave non-ok as NA
    } else if (is.character(df_all_imp[[col]])) {
      df_all_imp[[col]][pos] <- as.character(v)
    } else if (is.numeric(df_all_imp[[col]])) {
      suppressWarnings(num <- as.numeric(v))
      ok <- !is.na(num)
      df_all_imp[[col]][pos[ok]] <- num[ok]
    } else {
      tmp <- as.character(df_all_imp[[col]]); tmp[pos] <- as.character(v)
      df_all_imp[[col]] <- tmp
    }
  }
  
  # split back
  train_imp <- df_all_imp[!is_test, , drop = FALSE]
  test_imp  <- df_all_imp[ is_test, , drop = FALSE]
  
  # train-driven fallback: numeric -> median; non-numeric -> mode
  for (col in colnames(train_imp)) {
    x_tr <- train_imp[[col]]
    if (is.numeric(x_tr)) {
      med <- median(x_tr, na.rm = TRUE)
      if (anyNA(train_imp[[col]])) train_imp[[col]][is.na(train_imp[[col]])] <- med
      if (anyNA(test_imp[[col]]))  test_imp[[col]][is.na(test_imp[[col]])]   <- med
    } else {
      md <- Mode(x_tr)
      if (anyNA(train_imp[[col]])) train_imp[[col]][is.na(train_imp[[col]])] <- md
      if (anyNA(test_imp[[col]]))  test_imp[[col]][is.na(test_imp[[col]])]   <- md
    }
  }
  
  # write back into original data.frames
  out_train <- train_df; out_test <- test_df
  out_train[, selected_col] <- train_imp
  out_test[,  selected_col] <- test_imp
  
  list(train = out_train, test = out_test)
}

bootstrap_sampl_split_miss_vale_imputation <- function(expr_mat, clinical_df) {
  # Bootstrap Sampling and Data Preparation#
  # Performs bootstrap sampling, data splitting, and missing value imputation
  # for expression and clinical data.
  #
  # @param expr_mat Expression matrix (genes x samples)
  # @param clinical_df Clinical data frame (samples x variables)
  # @return List containing processed training and test data, or NULL if skipped
  
  
  n <- ncol(expr_mat)
  train_idx <- sample(seq_len(n), size = n, replace = TRUE)
  test_idx <- setdiff(seq_len(n), unique(train_idx))  # OOB samples
  
  message(sprintf("train_idx length = %d, test_idx length = %d", 
                  length(train_idx), length(test_idx)))
  
  # Skip if no OOB samples
  if (length(test_idx) == 0) {
    message("skipped: No OOB samples")
    return(NULL)
  }
  
  # Split data
  train_expr <- expr_mat[, train_idx, drop = FALSE]
  test_expr <- expr_mat[, test_idx, drop = FALSE]
  train_clinical <- clinical_df[train_idx, , drop = FALSE]
  test_clinical <- clinical_df[test_idx, , drop = FALSE]
  
  # Prepare columns for imputation
  # remove not relevant columns for imputation
  columns_to_remove <- c("e_dmfs", "t_dmfs", "geo_accession", "hospital")
  selected_col <- setdiff(colnames(clinical_df), columns_to_remove)
  missed_col <- colnames(clinical_df)[colSums(is.na(clinical_df)) > 0]
  
  # Apply imputation
  res <- impute_fit_apply(
    train_df = train_clinical,
    test_df = test_clinical,
    selected_col = selected_col,
    missed_col = missed_col
  )
  
  train_clinical <- res$train
  test_clinical <- res$test
  
  # Validate dimensions
  message("Sampling completed: train_expr dim = ", toString(dim(train_expr)), 
          ", test_expr dim = ", toString(dim(test_expr)),
          ", train_clinical dim = ", toString(dim(train_clinical)),
          ", test_clinical dim = ", toString(dim(test_clinical)))
  
  # Skip if invalid dimensions
  if (ncol(train_expr) == 0 || nrow(train_expr) == 0) {
    message("skipped: Invalid train_expr dimensions")
    return(NULL)
  }
  
  # Count events in training set
  events_train <- sum(train_clinical$e_dmfs)
  
  return(list(
    train_expr = train_expr,
    test_expr = test_expr,
    train_clinical = train_clinical,
    test_clinical = test_clinical,
    events_train = events_train
  ))
}

filter_genes_by_mad <- function(train_expr, test_expr, quantile_cutoff = 0.1) {
  # Filter Genes by MAD (Median Absolute Deviation)  #
  # Filters genes based on MAD threshold, keeping top 90% most variable genes.  #
  # @param train_expr Training expression matrix (genes x samples)
  # @param test_expr Test expression matrix (genes x samples)
  # @param quantile_cutoff Quantile cutoff for MAD filtering (default: 0.1)
  # @return List containing filtered train and test expression matrices, or NULL if no genes pass filter
  
  # Calculate MAD for each gene
  mad_train_expr <- apply(train_expr, 1, mad)
  
  # Determine cutoff value
  cutoff <- quantile(mad_train_expr, quantile_cutoff)
  
  # Identify genes to keep (top 90% most variable)
  keep_mad <- mad_train_expr > cutoff
  
  message(sprintf("MAD filtering completed: %d genes retained, cutoff = %.4f", 
                  sum(keep_mad), cutoff))
  
  # Filter both train and test matrices
  train_expr2 <- train_expr[keep_mad, , drop = FALSE]
  test_expr2 <- test_expr[keep_mad, , drop = FALSE]
  
  return(list(
    train_expr = train_expr2,
    test_expr = test_expr2
  ))
}

batch_univariate_cox_regression <- function(train_expr, train_clinical, p_value) {
  # Description: Performs univariate Cox regression for each gene in the expression matrix
  # to evaluate associations with survival outcomes (DMFS time and status).
  # Input: train_expr (gene expression matrix), train_clinical (clinical data with t_dmfs and e_dmfs)
  # Output: sig_genes_df (data frame of significant genes with p < 0.05)
  
  
  # Perform Cox regression for each gene
  cox_results <- apply(train_expr, 1, function(gene_expr) {
    df <- data.frame(
      expr = gene_expr,                    # Gene expression
      time = train_clinical$t_dmfs,        # Survival time
      status = train_clinical$e_dmfs       # Event status (1=event, 0=censored)
    )
    
    fit <- tryCatch(
      coxph(Surv(time, status) ~ expr, data = df),  # Cox model
      error = function(e) return(NULL)              # Return NULL if fails
    )
    
    if (is.null(fit)) return(c(NA, NA, NA, NA))
    
    s <- summary(fit)
    c(coef = s$coefficients[1, "coef"],          # Log hazard ratio
      HR = s$coefficients[1, "exp(coef)"],       # Hazard ratio
      SE = s$coefficients[1, "se(coef)"],        # Standard error
      p.value = s$coefficients[1, "Pr(>|z|)"])   # P-value
  })
  
  # Organize results into data frame
  cox_df <- as.data.frame(t(cox_results))
  cox_df$gene <- rownames(cox_df)          # Add gene names
  cox_df <- cox_df[order(cox_df$p.value), ]  # Sort by p-value
  
  # Filter significant genes (p < 0.05)
  sig_genes_df <- subset(cox_df, p.value < p_value)
  
  return(sig_genes_df)
}

standardize_with_train <- function(gene_mat_train, gene_mat_valid, significant_gene) {
  # Standardizes gene_mat_valid by using the column means and standard
  # deviations of gene_mat_trainn.
  # gene_mat_train is the original training data. 
  # Note: In the data, columns represent samples and rows represent genes.
  # Standardization should be applied row-wise.
  
  significant_gene <- sub("^ge_", "", significant_gene)
  
  selected_gene_mat_train <- gene_mat_train[rownames(gene_mat_train) %in% significant_gene, ]
  
  selected_gene_mat_valid <- gene_mat_valid[rownames(gene_mat_valid) %in% significant_gene, ]
  
  means <- apply(selected_gene_mat_train, 1, mean) #RGIN = 2 means "operate over columns."
  sds <- apply(selected_gene_mat_train, 1, sd)
  
  
  selected_gene_mat_valid_scaled <- sweep(selected_gene_mat_valid, 1, means, "-")
  selected_gene_mat_valid_scaled <- sweep(selected_gene_mat_valid_scaled, 1, sds, "/")
  
  return(selected_gene_mat_valid_scaled)
}

standardize_with_train_clinical <- function(train_clinical, test_clinical, scale_cols) {
  # Description: Standardize selected columns in test_clinical using the mean 
  # and sd from train_clinical
  
  # Calculate column-wise means and standard deviations from the training data
  means <- sapply(train_clinical[, scale_cols, drop = FALSE], mean, na.rm = TRUE)
  sds   <- sapply(train_clinical[, scale_cols, drop = FALSE], sd, na.rm = TRUE)
  
  # Copy the test data
  test_scaled <- test_clinical
  
  # Subtract means
  test_scaled[, scale_cols] <- sweep(test_clinical[, scale_cols, drop = FALSE], 2, means, "-")
  # Divide by standard deviations
  test_scaled[, scale_cols] <- sweep(test_scaled[, scale_cols, drop = FALSE], 2, sds, "/")
  
  # Return standardized test data
  return(test_scaled)
}

remove_high_corr_genes <- function(expr_mat, cutoff = 0.90) {
  # expr_mat: genes x samples
  if (is.null(expr_mat) || !is.matrix(expr_mat)) return(expr_mat)
  if (nrow(expr_mat) < 2) return(expr_mat)  # when less than 2 gene, return
  cm <- suppressWarnings(cor(t(expr_mat), use = "pairwise.complete.obs"))
  if (!is.matrix(cm) || ncol(cm) < 2) return(expr_mat)
  
  idx <- caret::findCorrelation(cm, cutoff = cutoff)
  if (length(idx) == 0) return(expr_mat)
  
  expr_mat[-idx, , drop = FALSE]
}

fit_cox_model <- function(predictors, df) {
  tryCatch({
    formula <- as.formula(paste(
      "Surv(t_dmfs, e_dmfs) ~",
      paste(predictors, collapse = " + ")
    ))
    
    # Fit Cox model
    cox_model <- coxph(formula, data = df, x = TRUE, y = TRUE)
    
    # Model summary
    cat("Cox model summary:\n")
    model_summary <- summary(cox_model)
    print(model_summary)
    
    # Coefficient table
    coef_df <- data.frame(
      variable = rownames(model_summary$coefficients),
      coef = model_summary$coefficients[, "coef"],
      zscore = model_summary$coefficients[, "z"],
      pvalue = model_summary$coefficients[, "Pr(>|z|)"]
    )
    
    # Test PH assumption
    cat("Proportional hazards (PH) assumption test:\n")
    test_ph <- cox.zph(cox_model)
    print(test_ph)
    
    # Return a fixed structure
    list(
      model = cox_model,
      ph_test = test_ph,
      coef_table = coef_df
    )
  }, error = function(e) {
    message("Error in fit_cox_model: ", e$message)
    # Return a fixed structure with NULL values
    list(
      model = NULL,
      ph_test = NULL,
      coef_table = NULL
    )
  })
}

repeat_cv_lasso_cox <- function(train_expr = train_expr_filtered,
                                train_clinical,
                                significant_gene_vec= significant_gene2,
                                repeats = 20,
                                nfolds = 10,
                                alpha = 1) {
  # Performs repeated cross-validated LASSO Cox regression to identify 
  # consistently selected genes across multiple iterations. Returns genes ranked
  # by selection frequency and average coefficient magnitude.
  
  # Transpose expression matrix for glmnet compatibility
  X <- t(train_expr[significant_gene_vec, ])
  
  y <- Surv(train_clinical$t_dmfs, train_clinical$e_dmfs)
  
  # Initialize gene selection counters
  gene_names <- colnames(X)
  gene_counts <- setNames(rep(0, length(gene_names)), gene_names)
  gene_coef_sum <- setNames(rep(0, length(gene_names)), gene_names)
  
  #  Repeated cross-validation loop
  for (i in 1:repeats) {
    cvfit <- cv.glmnet(X, y, family = "cox", alpha = alpha, nfolds = nfolds)
    coef_nonzero <- as.numeric(coef(cvfit, s = "lambda.min"))
    names(coef_nonzero) <- rownames(coef(cvfit))
    
    # Update gene selection frequency and coefficient sums
    for (gene in gene_names) {
      if (coef_nonzero[gene] != 0) {
        gene_counts[gene] <- gene_counts[gene] + 1
        gene_coef_sum[gene] <- gene_coef_sum[gene] + coef_nonzero[gene]
      }
    }
  }
  
  # Calculate selection frequency and average coefficients
  gene_freq_df <- data.frame(
    gene = gene_names,
    freq = gene_counts / repeats,
    mean_coef = ifelse(gene_counts > 0, gene_coef_sum / gene_counts, 0)
  )
  # Sort by selection frequency and coefficient magnitude
  gene_freq_df <- gene_freq_df[order(-gene_freq_df$freq, -abs(gene_freq_df$mean_coef)), ]
  return(gene_freq_df)
}

compute_risk_score_train_2 <- function(gene_mat_scaled, significant_vars_df, clinical_cleaned) {
  # Computes risk scores based on scaled gene expression and model coefficients,
  # and assigns each sample to a risk group (low vs high) based on the median risk score.
  #
  # Arguments:
  # - gene_mat_scaled: matrix of scaled gene expression values (genes x samples)
  # - significant_vars_df: data frame with selected gene names (rownames) and coefficients
  # - clinical_cleaned: data frame with clinical data (samples as rows)
  #
  # Returns:
  # - A list containing:
  #     median_cut: the median cutoff used for risk grouping
  #     merged_df: clinical data with added columns: 'risk_score', 'risk_group' and gene data
  
  genes <- sub("^ge_", "", rownames(significant_vars_df))
  
  if(!"geo_accession" %in% colnames(clinical_cleaned)) {
    clinical_cleaned$geo_accession <- rownames(clinical_cleaned)
  }
  
  # Match sample order and subset by selected genes
  gene_mat_scaled <- gene_mat_scaled[, match(clinical_cleaned$geo_accession, colnames(gene_mat_scaled))]
  gene_mat_scaled <- gene_mat_scaled[genes, , drop = FALSE]
  
  if (!identical(rownames(gene_mat_scaled), rownames(significant_vars_df))){
    stop("Row names are not identical or not in the same order.")}
  
  # Compute risk score
  clinical_cleaned$risk_score <- as.vector(t(gene_mat_scaled) %*% significant_vars_df$coef)
  
  # median cutoff for binary grouping
  median_cut <- median(clinical_cleaned$risk_score, na.rm = TRUE)
  
  # Two groups: <= median = low risk, > median = high risk
  clinical_cleaned <- clinical_cleaned %>%
    dplyr::mutate(
      risk_group = dplyr::case_when(
        risk_score <= median_cut ~ 1,
        risk_score >  median_cut ~ 2,
        TRUE ~ NA_real_
      ),
      risk_group = factor(risk_group, levels = c(1, 2), labels = c("low", "high"), ordered = TRUE)
    )
  
  # Add gene expression data back to sample-level data
  gene_mat_scaled_t <- t(gene_mat_scaled)
  gene_df <- as.data.frame(gene_mat_scaled_t)
  
  # Add prefix "ge_" to column names starting with a digit
  colnames(gene_df)[grepl("^[0-9]", colnames(gene_df))] <-
    paste0("ge_", colnames(gene_df)[grepl("^[0-9]", colnames(gene_df))])
  
  gene_df$geo_accession <- rownames(gene_df)
  merged_df <- merge(clinical_cleaned, gene_df, by = "geo_accession", sort = FALSE)
  
  return(list(
    median_cut  = median_cut,
    merged_df = merged_df
  ))
}

compute_risk_score_test_2 <- function(gene_mat_scaled, significant_vars_df, clinical_cleaned, median_cut_train) {
  # Computes risk scores for the test set based on scaled gene expression and model coefficients,
  # and assigns each sample to a risk group (low vs high) using the training median cutoff.
  #
  # Arguments:
  # - gene_mat_scaled: matrix of scaled gene expression values (genes x samples)
  # - significant_vars_df: data frame with selected gene names (rownames) and coefficients
  # - clinical_cleaned: data frame with clinical data (samples as rows)
  # - median_cut_train: numeric, median cutoff value from the training set
  #
  # Returns:
  # - merged_df: clinical data with added columns: 'risk_score', 'risk_group' and gene data
  
  genes <- sub("^ge_", "", rownames(significant_vars_df))
  
  if(!"geo_accession" %in% colnames(clinical_cleaned)) {
    clinical_cleaned$geo_accession <- rownames(clinical_cleaned)
  }
  
  gene_mat_scaled <- gene_mat_scaled[, match(clinical_cleaned$geo_accession, colnames(gene_mat_scaled))]
  gene_mat_scaled <- gene_mat_scaled[genes, , drop = FALSE]
  
  if (!identical(rownames(gene_mat_scaled), rownames(significant_vars_df))){
    stop("Row names are not identical or not in the same order.")
  }
    
  
  # Compute risk score
  clinical_cleaned$risk_score <- as.vector(t(gene_mat_scaled) %*% significant_vars_df$coef)
  
  # Two groups: <= median_cut_train = low, > median_cut_train = high
  clinical_cleaned <- clinical_cleaned %>%
    dplyr::mutate(
      risk_group = dplyr::case_when(
        risk_score <= median_cut_train ~ 1,
        risk_score >  median_cut_train ~ 2,
        TRUE ~ NA_real_
      ),
      risk_group = factor(risk_group, levels = c(1, 2), labels = c("low", "high"), ordered = TRUE)
    )
  
  # Add gene data
  gene_mat_scaled_t <- t(gene_mat_scaled)
  gene_df <- as.data.frame(gene_mat_scaled_t)
  
  # Add prefix "ge_" to column names starting with a digit
  colnames(gene_df)[grepl("^[0-9]", colnames(gene_df))] <-
    paste0("ge_", colnames(gene_df)[grepl("^[0-9]", colnames(gene_df))])
  
  gene_df$geo_accession <- rownames(gene_df)
  merged_df <- merge(clinical_cleaned, gene_df, by = "geo_accession", sort = FALSE)
  
  return(merged_df)
}

km_by_group <- function(df, group_var) {
  # Convert group_var to factor
  df$group_var <- as.factor(df[[group_var]])
  
  # Create survival object
  surv_obj <- Surv(time = df$t_dmfs, event = df$e_dmfs)
  
  # Survival difference test
  surv_diff <- survdiff(surv_obj ~ group_var, data = df)
  
  # Calculate p-value
  p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  cat("Log-rank Kaplan-Meier test p-value:", p_value, "\n")
  
  fit <- survfit(surv_obj ~ group_var, data = df)
  
  return(list(
    p_value = p_value,
    fit = fit
    ))
}

rsf_kfold_cv_best <- function(data, K = 5, ntree = 1000) {
  # Performs K-fold cross-validation for a Random Survival Forest model.
  # Retrain the final model on the entire dataset to extract variable importance.
  # Returns the model with the highest C-index on the validation fold.
  # Output:
  #   Best model
  #   C-index per fold
  #   Best fold number
  #   Variable importance
  
  folds <- sample(rep(1:K, length.out = nrow(data)))
  cindex_vec <- numeric(K)
  model_list <- vector("list", K)
  
  for (k in 1:K) {
    cat("Fold", k, "\n")
    
    # Split data
    test_idx <- which(folds == k)
    train_data <- data[-test_idx, ]
    test_data  <- data[test_idx, ]
    
    # Train RSF
    model <- rfsrc(Surv(t_dmfs, e_dmfs) ~ ., data = train_data, ntree = ntree)
    model_list[[k]] <- model
    
    # Evaluate on test set
    c_idx <- cindex(object = model,
                    formula = Surv(t_dmfs, e_dmfs) ~ 1,
                    data = test_data)
    cindex_vec[k] <- round(as.numeric(c_idx$AppCindex), 4)
  }
  
  # Find best fold
  best_fold <- which.max(cindex_vec)
  best_model <- model_list[[best_fold]]
  
  cat("Best Fold:", best_fold, "\n")
  cat("Best C-index:", cindex_vec[best_fold], "\n")
  
  # retrain the final model on the entire dataset to extract variable importance.
  final_model <- rfsrc(Surv(t_dmfs, e_dmfs) ~ ., data = data,
                       ntree = 1000, importance = TRUE)
  
  var_importance <- final_model$importance
  print("Variable Importance(retrained on full data):")
  print(var_importance )
  
  return(list(
    cindex_per_fold = cindex_vec,
    best_model = best_model,
    best_fold = best_fold,
    importance = var_importance
  ))
}

calculate_time_auc_cindex <- function(model_type = c("Cox", "RSF"), fitted_model, df) {
  # This function evaluates the discrimination performance of a fitted survival 
  # model (either Cox or Random Survival Forest) on a given dataset. It performs
  # the following:
  # Computes time-dependent AUC using the Score() function.
  # Calculates the time-weighted average AUC (iAUC).
  # Calculates the concordance index (C-index).
  
  # Inputs:
  #   model_type: "Cox" or "RSF".
  # fitted_model: the trained model object.
  # df: dataset for evaluation.
  # Outputs:
  # AUC vs time plot.
  # Printed iAUC and C-index values.
  # Returns a list containing the AUC data frame, iAUC, and C-index.
  
  model <- match.arg(model_type)
  
  # extract unique event times
  time_points <- sort(unique(df$t_dmfs[df$e_dmfs == 1]))
  if (length(time_points) > 100) {
    time_points <- seq(min(time_points), max(time_points), length.out = 100)
  }
  
  # compute time-dependent 
  try_result <- try(
    Score(
      object = list(model = fitted_model),
      formula = Surv(t_dmfs, e_dmfs) ~ 1,
      data = df,
      metrics = "AUC",
      summary = "AUC",
      times = time_points,
      conf.int = FALSE,
      plots = FALSE,
      split.method = "none"
    ),
    silent = TRUE
  )
  
  if (inherits(try_result, "try-error")) {
    df = df %>% filter(t_dmfs < time_points[length(time_points)])
    time_points <- time_points[-length(time_points)]
    auc_result <- Score(
      object = list(model = fitted_model),
      formula = Surv(t_dmfs, e_dmfs) ~ 1,
      data = df,
      metrics = "AUC",
      summary = "AUC",
      times = time_points,
      conf.int = FALSE,
      plots = FALSE,
      split.method = "none"
    )
  } else {
    auc_result <- try_result
  }
  
  # extract AUC values
  auc_values <- auc_result$AUC$score
  auc_df <- data.frame(Time = auc_values$times, AUC = auc_values$AUC)
  
  # Bin event counts over small time intervals
  bin_width <- diff(range(time_points)) / length(time_points)
  event_counts <- sapply(time_points, function(t) {
    sum(df$e_dmfs == 1 & df$t_dmfs >= (t - bin_width/2) & df$t_dmfs < (t + bin_width/2))
  })
  
  time_diffs <- diff(c(auc_df$Time[1] - (auc_df$Time[2] - auc_df$Time[1]), auc_df$Time))
  weights <- time_diffs / sum(time_diffs)
  iAUC <- sum(auc_df$AUC * weights, na.rm = TRUE)
  
  # Compute C-index
  c_index_result <- cindex(
    object = fitted_model,
    formula = Surv(t_dmfs, e_dmfs) ~ 1,
    data = df
  )
  c_index <- round(as.numeric(c_index_result$AppCindex), 4)
  
  # Output summary
  cat("Time-weighted average AUC (iAUC):", round(iAUC, 4), "\n")
  cat("Concordance Index (C-index):", c_index, "\n")
  
  invisible(list(auc_df = auc_df, iAUC = iAUC, c_index = c_index))
}


run_bootstrap_iteration_1_model <- function(predictors_filtered,
                                    df,
                                    b,
                                    coef_max,
                                    min_concord, 
                                    clinical_cleaned_risk_test,
                                    selected_gene_df,
                                    perf_list,
                                    rsf_predictor_list,
                                    best_perf,
                                    best_model,
                                    best_iter,
                                    best_seed,
                                    best_method,
                                    current_seed) {
  
  # This function performs one iteration of a bootstrap sampling and model comparison process. 
  # Its main purpose is to:   
  # Train two types of survival models (a Cox Proportional Hazards model and a Random Survival Forest) on a sample of the data. 
  # Evaluate how well these models perform on a separate test dataset.
  # Check for problems like model divergence or overfitting.   
  # Keep track of the best-performing model across all bootstrap iterations, saving it if the current iteration's model is better.
  
  # Inputs:
  #   
  # predictors_filtered: Selected features for current iteration 
  # df: Training dataset with survival data  
  # b: Current bootstrap iteration number  # 
  # coef_max: Maximum allowed coefficient value for stability check  
  # min_concord: Minimum concordance threshold for overfitting detection 
  # clinical_cleaned_risk_test: Test dataset for model evaluation 
  # perf_list: Accumulated performance results from previous iterations  
  # rsf_predictor_list: Accumulated predictor lists from previous RSF models   
  # best_* parameters: Current best performance metrics and model details 
  # current_seed: Random seed for current iteration
  
  # Outputs: Returns a list containing:   
  # Updated perf_list and rsf_predictor_list with current results  
  # Updated best_perf, best_model, best_iter, best_seed, best_method if current models outperform 
  # cox_model and rsf_model objects from previous iteration  
  # iteration number for tracking
  
   
  
  # Fit Cox model
  results_train <- fit_cox_model(predictors_filtered, df)
  message(sprintf("Bootstrap %d: Cox model fitting completed, model exists = %s",
                  b, !is.null(results_train$model)))
  
  if (is.null(results_train$model)) {
    message(sprintf("Bootstrap %d skipped, Cox model is null ", as.integer(b)))
    return(list(perf_list = perf_list, rsf_predictor_list = rsf_predictor_list,
                best_perf = best_perf, best_model = best_model,
                best_iter = best_iter, best_seed = best_seed, best_method = best_method))
  }
  
  # Divergence detection
  if (any(abs(coef(results_train$model)) > coef_max)) {
    message(sprintf("Bootstrap %d skipped: coef > %0.1f detected", as.integer(b), coef_max))
    return(list(perf_list = perf_list, rsf_predictor_list = rsf_predictor_list,
                best_perf = best_perf, best_model = best_model,
                best_iter = best_iter, best_seed = best_seed, best_method = best_method))
  }
  message("Cox model divergence detection completed")
  
  # Check complete separation
  concordance_val <- tryCatch({
    suppressWarnings(summary(results_train$model)$concordance[1])
  }, error = function(e) NA)
  message(sprintf("Bootstrap %d: Complete separation check completed, concordance_val = %.4f",
                  b, if (is.na(concordance_val)) NA else concordance_val))
  
  if (!is.na(concordance_val) && concordance_val >= min_concord) {
    message(sprintf("Bootstrap %d skipped, Concordance >= %f detected ", as.integer(b), min_concord))
    return(list(perf_list = perf_list, rsf_predictor_list = rsf_predictor_list,
                best_perf = best_perf, best_model = best_model,
                best_iter = best_iter, best_seed = best_seed, best_method = best_method))
  }
  
  # Cox model evaluation
  result_valid <- calculate_time_auc_cindex(
    "Cox",
    fitted_model = results_train$model,
    df = clinical_cleaned_risk_test
  )
  message(sprintf("Bootstrap %d: Cox model evaluation completed, iAUC = %.4f, c_index = %.4f",
                  b, result_valid$iAUC, result_valid$c_index))
  
  # Random Survival Forest
  clinical_rsf <- df
  result_rsf_train <- rsf_kfold_cv_best(data = clinical_rsf, K = 5, ntree = 1000)
  rsf_fit_best <- result_rsf_train$best_model
  message("RSF model fitting completed")
  
  result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
  message(sprintf("Bootstrap %d: RSF model evaluation completed, iAUC = %.4f, c_index = %.4f",
                  b, result_rsf_valid$iAUC, result_rsf_valid$c_index))
  
  # Save performance
  perf_list[[b]] <- list(
    cox = result_valid,
    rsf = result_rsf_valid
  )
  rsf_predictor_list[[b]] <- predictors_filtered
  message("Performance saving completed")
  
  # Calculate Cox and RSF composite scores
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
  
  # Check if better (Cox)
  if (cox_score > best_perf) {
    best_perf    <- cox_score
    best_model   <- results_train$model
    best_iter    <- b
    best_seed    <- current_seed
    best_method  <- "Cox"
  }
  
  # Check if better (RSF)
  if (rsf_score > best_perf) {
    best_perf    <- rsf_score
    best_model   <- rsf_fit_best
    best_iter    <- b
    best_seed    <- current_seed
    best_method  <- "RSF"
  }
  
  cat(sprintf("Best model found in iteration %d using %s model\n", best_iter, best_method))
  cat(sprintf("Best iAUC + C-index = %.4f\n", best_perf))
  message(sprintf("Bootstrap %d: Iteration completed", b))
  
  # Return all updated variables
  return(list(
    perf_list = perf_list,
    rsf_predictor_list = rsf_predictor_list,
    best_perf = best_perf,
    best_model = best_model,
    best_iter = best_iter,
    best_seed = best_seed,
    best_method = best_method,
    cox_model = results_train$model,
    rsf_model = rsf_fit_best,
    iteration = b
  ))
}

run_bootstrap_validation_3_model <- function(expr_mat, clinical_df, 
                                             B = 10, 
                                             seed = 90000,
                                             min_epv = 2.5, 
                                             coef_max = 10,
                                             min_concord = 0.98) {
  # Description:
  # This function performs bootstrap validation of three survival models 
  # (clinical-only, gene-only, and multimodal) using expression and clinical data. 
  # The process includes:
  #   - Bootstrap sampling of training and test sets
  #   - Missing value imputation
  #   - Gene filtering (MAD, univariate Cox regression, correlation filtering)
  #   - Standardization of gene and clinical features
  #   - Variable selection via repeated cross-validated LASSO Cox
  #   - Risk score computation for training and test sets
  #   - Kaplan-Meier log-rank test on stratified risk groups
  #   - Model fitting and performance evaluation for multimodal, gene-only, 
  #     and clinic-only predictors
  #   - Aggregation of performance metrics, gene frequencies, and selected predictors
  #
  # Args:
  #   expr_mat: Gene expression matrix (genes x samples).
  #   clinical_df: Clinical data frame aligned with expression samples.
  #   B: Number of bootstrap iterations (default = 10).
  #   seed: Base random seed for reproducibility (default = 90000).
  #   min_epv: Minimum events per variable (default = 2.5).
  #   coef_max: Maximum allowed coefficient magnitude (default = 10).
  #   min_concord: Minimum concordance threshold for model selection (default = 0.98).
  #
  # Returns:
  #   A list containing:
  #     - summary_df: Summary of best models for each modality
  #     - best_model_mc: Best clinical-only model
  #     - best_model_mg: Best gene-only model
  #     - best_model_mm: Best multimodal model
  #     - train_indices_list: List of bootstrap training indices
  #     - mean_km_p: Mean p-value from Kaplan-Meier log-rank tests
  #     - km_p_values: All collected KM p-values
  #     - gene_frequency: Frequency of selected genes across bootstraps
  #     - summary_results_mm/mg/mc: Performance summaries for each model type
  
  km_p_list <- list() # save KM log-rank test p value for test_data
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
  
  
  
  message("Initialization completed: expr_mat dimensions = ", toString(dim(expr_mat)), 
          ", clinical_df dimensions = ", toString(dim(clinical_df)))
  
  for (b in 1:B) {
    # Bootstrap sampling
    current_seed <- seed + b 
    set.seed(current_seed) 
    message(sprintf("Bootstrap %d: seed = %d", b, seed + b))
    
    
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
    message("Event count completed")
    
    # MAD filtering
    filtered_data <- filter_genes_by_mad(train_expr, test_expr, quantile_cutoff = 0.25)
    train_expr2 <- filtered_data$train_expr
    test_expr2 <- filtered_data$test_expr  
    
    message("MAD filtering completed: train_expr2 dimensions = ", toString(dim(train_expr2)), 
            ", test_expr2 dimensions = ", toString(dim(test_expr2)))
    
    # Univariate Cox regression
    sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = 0.01)
    message(sprintf("Bootstrap %d: Univariate Cox regression completed, number of significant genes = %d", 
                    b, if (is.null(sig_gene_df)) 0 else nrow(sig_gene_df)))
    
    if (is.null(sig_gene_df) || nrow(sig_gene_df) == 0) {
      message(sprintf("Bootstrap %d skipped: No significant genes from univariate Cox", b))
      next
    }
    
    significant_gene <- sig_gene_df$gene
    message(sprintf("Bootstrap %d: Extracting significant genes completed, significant_gene length = %d", 
                    b, length(significant_gene)))
    
    if (length(significant_gene) < 2) {
      message(sprintf("Bootstrap %d skipped: significant genes < 2", b))
      next
    }
    
    # Standardize training gene data
    train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                                gene_mat_valid = train_expr2,
                                                significant_gene = significant_gene)
    message("Standardization of training gene data completed: train_expr_scaled dimensions = ", 
            toString(dim(train_expr_scaled)))
    
    if (is.null(train_expr_scaled) || nrow(train_expr_scaled) == 0) {
      message(sprintf("Bootstrap %d skipped: Invalid standardized training data", b))
      next
    }
    
    if (nrow(train_expr_scaled) < 2) {
      message(sprintf("Bootstrap %d skipped: nrow(train_expr_scaled) < 2", b))
      next
    }
    
    # Correlation filtering
    train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = 0.90)
    message("Correlation filtering completed: train_expr_filtered dimensions = ", 
            toString(dim(train_expr_filtered)))
    
    if (is.null(train_expr_filtered) || nrow(train_expr_filtered) < 2) {
      message(sprintf("Bootstrap %d skipped: No genes after correlation filtering", b))
      next
    }
    
    significant_gene2 <- rownames(train_expr_filtered)
    message(sprintf("Bootstrap %d: Extracting post-filtering genes completed, significant_gene2 length = %d", 
                    b, length(significant_gene2)))
    
    # Standardize testing gene data
    test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                               gene_mat_valid = test_expr2,
                                               significant_gene = significant_gene2)
    message("Standardization of testing gene data completed: test_expr_scaled dimensions = ", 
            toString(dim(test_expr_scaled)))
    
    # Standardize training clinical data
    numeric_columns <- train_clinical %>%
      select(where(is.numeric)) %>%
      colnames
    
    numeric_columns <- setdiff(numeric_columns, columns_to_remove)
    
    train_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                             train_clinical,
                                                             scale_cols = numeric_columns)
    message("Standardization of training clinical data completed: train_clinical_scaled dimensions = ", 
            toString(dim(train_clinical_scaled)))
    
    # Standardize testing clinical data
    test_clinical_scaled <- standardize_with_train_clinical(train_clinical,
                                                            test_clinical,
                                                            scale_cols = numeric_columns)
    message("Standardization of testing clinical data completed: test_clinical_scaled dimensions = ", 
            toString(dim(test_clinical_scaled)))
    
    # LASSO Cox
    gene_freq_df <- repeat_cv_lasso_cox(train_expr = train_expr_filtered,
                                        train_clinical,
                                        significant_gene_vec = significant_gene2,
                                        repeats = 5,
                                        nfolds = 10,
                                        alpha = 1)
    message(sprintf("Bootstrap %d: LASSO Cox completed, gene_freq_df rows = %d", 
                    b, nrow(gene_freq_df)))
    
    # Dynamic variable selection to satisfy EPV requirement
    gene_freq_df_best <- gene_freq_df %>%
      filter(freq >= 0.8)
    max_vars_allowed <- floor(events_train / min_epv)
    message(sprintf("Bootstrap %d: Dynamic variable selection completed, gene_freq_df_best rows = %d, max_vars_allowed = %d", 
                    b, nrow(gene_freq_df_best), max_vars_allowed))
    
    if (max_vars_allowed == 0 || nrow(gene_freq_df) == 0) {
      message(sprintf("Bootstrap %d skipped: Too few events or no selected genes", b))
      next
    }
    
    selected_gene_df <- gene_freq_df[1:min(nrow(gene_freq_df_best), max_vars_allowed), ] %>%
      mutate(coef = mean_coef)
    message(sprintf("Bootstrap %d: Gene selection completed, selected_gene_df rows = %d, EPV = %.2f", 
                    b, nrow(selected_gene_df), events_train / nrow(selected_gene_df)))
    
    gene_list[[b]] <- selected_gene_df$gene
    
    result <- compute_risk_score_train_2(
      gene_mat_scaled  = train_expr_filtered,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = train_clinical_scaled
    )
    # Extract components
    mean_cut_train <- result$mean_cut
    clinical_cleaned_risk_train <- result$merged_df
    
    message("Training risk score calculation completed: clinical_cleaned_risk_train dimensions = ", 
            toString(dim(clinical_cleaned_risk_train)))
    
    # Testing risk score
    clinical_cleaned_risk_test <- compute_risk_score_test_2(
      gene_mat_scaled = test_expr_scaled,
      significant_vars_df = selected_gene_df,
      clinical_cleaned = test_clinical_scaled,
      mean_cut_train
      )
    
    message("Testing risk score calculation completed: clinical_cleaned_risk_test dimensions = ", 
            toString(dim(clinical_cleaned_risk_test)))
    
    km_p <- km_by_group(df = clinical_cleaned_risk_test, group_var = "risk_group")
    km_p_list[[b]] <- km_p
    
    columns_to_remove2 <- c("e_dmfs", "t_dmfs", "risk_score", "risk_group")
    character_columns <- clinical_cleaned_risk_train %>%
      select(where(is.character)) %>%
      colnames()
    predictors <- colnames(clinical_cleaned_risk_train) %>% setdiff(c(columns_to_remove2, character_columns))
    
    message(sprintf("Bootstrap %d: Predictor selection completed, predictors length = %d", 
                    b, length(predictors)))
    df <- clinical_cleaned_risk_train[, c(predictors, "t_dmfs", "e_dmfs")]
    
    
    # model multimodal
    model_m <- run_bootstrap_iteration_1_model(
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
    
    
    # model gene
    predictors_g <- predictors[5:length(predictors_filtered)]
    model_g <- run_bootstrap_iteration_1_model(
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
    
    
    # model clinic
    predictors_c <- predictors_filtered[1:4]
    model_c <- run_bootstrap_iteration_1_model(
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
    
  }
  
  message("All Bootstrap iterations completed")
  
  ##### # Summary
  
  
  # Create summary data frame
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
  
  # Compute gene frequency
  all_genes <- unlist(gene_list)
  gene_freq <- sort(table(all_genes) / B, decreasing = TRUE)
  message("Gene frequency calculation completed")
  
  # Compute mean km_p
  km_p_list_nz <- Filter(Negate(is.null), km_p_list) 
  km_p_values <- unlist(km_p_list_nz)
  mean_km_p <- mean(km_p_values, na.rm = TRUE)
  message("Performance summary completed: mean_km_p = ", mean_km_p)
  
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
  ))
}

#geo_accession_7390 <- clinical_cleaned_7390$geo_accession
#expr_mat_7390 <- expr[, geo_accession_7390, drop = FALSE]
# res7390_3_model4 <- run_bootstrap_validation_3_model(expr_mat = expr_mat_7390,
#                                                      clinical_df = clinical_cleaned_7390, 
#                                                      B = 100, 
#                                                      seed = 92000000, 
#                                                      min_epv = 5, 
#                                                      coef_max = 10, 
#                                                      min_concord = 0.97)

#saveRDS(res7390_3_model4, "res7390_3_model4.rds")

##### PCA COX
train_cox_with_pca <- function(data = merged_df, 
                               time_col,
                               status_col,
                               clinical_cleaned_7390,
                               train_ratio = 0.7, 
                               variance_threshold = 0.7
                               ) {
  
  # Prepare feature matrix and survival data
  features <- as.matrix(data[, !(names(data) %in% c(time_col, status_col))])
  surv_data <- data[, c(time_col, status_col)]
  #colnames(surv_data) <- c("time", "status")
  
  # Split data into training and testing sets
  train_index <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
  all_indices <- 1:nrow(data)
  test_index <- setdiff(all_indices, unique(train_index))
  
  
  train_feat <- features[train_index, ]
  test_feat <- features[test_index, ]
  train_surv <- surv_data[train_index,]
  test_surv <- surv_data[test_index, ]
  
  # Perform PCA on training set
  pca_train <- prcomp(train_feat, center = TRUE, scale. = TRUE)
  pca_summary <- summary(pca_train)
  
  events_train <- sum(train_surv$e_dmfs)
  max_vars_allowed <- floor(events_train / 3)
  # Enforced an events-per-variable (EPV) rule of 3 in order to explain at least
  # 50% of the variance in the original data.
  
  cumulative_variance <- pca_summary$importance[3, max_vars_allowed]  # Third row is cumulative variance
  cat("The first", max_vars_allowed, "principal components explain", 
      round(cumulative_variance * 100, 2), "% of the variance\n")
  
  # Get PC scores for training set
  train_pcs <- pca_train$x[, 1:max_vars_allowed] #max_vars_allowed
  clinical_t_e_2 <- clinical_cleaned_7390 %>% select(geo_accession, hospital)
  
  train_pcs2 <- as.data.frame(train_pcs) %>%
    rownames_to_column("geo_accession") %>%
    left_join(clinical_t_e_2, 
              by = "geo_accession") %>%
    column_to_rownames("geo_accession")
  
  # Transform test set using training PCA parameters
  test_scaled <- scale(test_feat, center = pca_train$center, scale = pca_train$scale)
  test_pcs <- test_scaled %*% pca_train$rotation
  test_pcs <- test_pcs[, 1:max_vars_allowed] # n_pcs
  test_pcs2 <- as.data.frame(test_pcs) %>%
    rownames_to_column("geo_accession") %>%
    left_join(clinical_t_e_2, 
              by = "geo_accession") %>%
    column_to_rownames("geo_accession")
  
  
  # Create dataframes for modeling
  train_df <- data.frame(train_pcs2, t_dmfs = train_surv$t_dmfs, e_dmfs = train_surv$e_dmfs)
  test_df <- data.frame(test_pcs2, t_dmfs = test_surv$t_dmfs, e_dmfs = test_surv$e_dmfs)
  
  # Build COX model
  ## COX ME model
  library(coxme)
  cox_me_formula <- as.formula(paste("Surv(t_dmfs, e_dmfs) ~", 
                                     paste(colnames(train_pcs), collapse = " + "),"+ (1 | hospital)"))
  cox_me_model <- coxme(cox_me_formula, data = train_df, x = TRUE)
  cox_me_model_summ <- summary(cox_me_model)
  
  # this p value is to measue the random effect of hospital, smaller value indicates significient effect
  cox_me_p_value <- as.data.frame(cox_me_model_summ$chi) %>% select("p")
  coxme_cindex <- evaluate_coxme_cindex(cox_me_model, test_df, time_var = "t_dmfs", status_var = "e_dmfs")
  
  
  ## COX model
  cox_formula <- as.formula(paste("Surv(t_dmfs, e_dmfs) ~", 
                                  paste(colnames(train_pcs), collapse = " + ")))
  cox_model <- coxph(cox_formula, data = train_df, x = TRUE,  y=TRUE, model=TRUE)
  ph_test <- cox.zph(cox_model)
  global_ph_test <- ph_test$table[nrow(ph_test$table),3]
  # Make predictions on test set
  # test_risk <- predict(cox_model, newdata = test_df, type = "risk")
  
  test_result <- try(
    calculate_time_auc_cindex(model_type = "Cox", fitted_model= cox_model, df = test_df),
    silent = TRUE
  )
  
  failed <- FALSE
  if (inherits(test_result, "try-error")) {
    failed <- TRUE
  } else {
    # Safely extract results
    ci   <- test_result$c_index
    iauc <- test_result$iAUC
    # If either metric is non-finite, mark as failed
    if (!is.finite(ci) || !is.finite(iauc)) failed <- TRUE
  }
  
  if (failed) {
    message("Iteration: scoring failed or returned non-finite values. Skipping.")
    stop()  # <--- Key: skip to the next iteration
  }
  
  # Return results
  return(list(
    model = cox_model,
    train_pca = pca_train,
    max_vars_allowed = max_vars_allowed,
    cumulative_variance_explained = cumulative_variance,
    c_index = ci, #test_result$c_index,
    iAUC = iauc,#test_result$iAUC,
    total_score = ci + iauc, #test_result$c_index + test_result$iAUC,
    global_ph_test = global_ph_test,
    cox_me_model_summary = summary(cox_me_model),
    cox_me_p_value = cox_me_p_value,
    cox_me_cindex = coxme_cindex,
    cox_model_summary = summary(cox_model)
  ))
}

# Function to evaluate Cox model with PCA over multiple runs
bootstrap_train_cox_with_pca <- function(data, clinical_df, n_runs, seed) {
  
  # Initialize vectors to store performance metrics
  c_index_results <- numeric(n_runs)
  cox_me_cindex_results <- numeric(n_runs)
  iAUC_results <- numeric(n_runs)
  total_score_results <- numeric(n_runs)
  
  # Initialize list to store all results
  all_results_list <- list()
  
  # Initialize variables to track best results
  best_total_score <- -Inf
  best_run_id <- 0
  best_results <- NULL
  
  # Run iterations
  for (i in 1:n_runs) {
    # Set random seed for reproducibility
    current_seed <- seed + i
    set.seed(current_seed)
    
    cat(sprintf("Running iteration %d/%d...\n", i, n_runs))
    
    results <- tryCatch(
        {
          train_cox_with_pca(
            data = data,
            time_col = "t_dmfs",
            status_col = "e_dmfs",
            clinical_df,
            train_ratio = 0.7,
            variance_threshold = 0.7
          )
        },
        error = function(e) {
          message(paste("Error at iteration", i, ":", e$message))
          return(NULL)   
        }
      )
      
      if (is.null(results)) {
        next   
      }
    
    # Store performance metrics
    c_index_results[i] <- results$c_index
    cox_me_cindex_results[i] <- results$cox_me_cindex
    iAUC_results[i] <- results$iAUC
    total_score_results[i] <- results$total_score
    
    # Store complete results
    all_results_list[[i]] <- results
    
    # Update best results if current total score is higher
    if (results$total_score > best_total_score) {
      best_total_score <- results$total_score
      best_run_id <- i
      best_results <- results
      cat(sprintf("  New best result found! Run ID: %d, Total Score: %.4f\n", i, best_total_score))
    }
    
    # Print progress
    cat(sprintf("  Completed: C-index = %.3f, iAUC = %.3f, Total = %.3f\n",
                results$c_index, results$iAUC, results$total_score))
  }
  
  # Calculate average results
  avg_results <- data.frame(
    Metric = c("C-index", "CoxMe C-index", "iAUC", "Total Score"),
    Mean = c(
      mean(c_index_results),
      mean(cox_me_cindex_results),
      mean(iAUC_results),
      mean(total_score_results)
    ),
    SD = c(
      sd(c_index_results),
      sd(cox_me_cindex_results),
      sd(iAUC_results),
      sd(total_score_results)
    )
  )
  
  # Create best results dataframe
  best_results_df <- data.frame(
    Item = c("Best Run ID", "Best Seed", "C-index", "iAUC", "Total Score",
             "Number of PCs", "Cumulative Variance (%)", "Global PH Test p-value"),
    Value = c(
      best_run_id,
      seed + best_run_id,
      best_results$c_index,
      best_results$iAUC,
      best_results$total_score,
      best_results$max_vars_allowed,
      best_results$cumulative_variance_explained * 100,
      best_results$global_ph_test
    )
  )
  
  # Return results
  return(list(
    avg_results = avg_results,
    best_results_df = best_results_df,
    best_results = best_results
  ))
}

# pca_result_9_18_2 <- bootstrap_train_cox_with_pca(data = merged_df,
#                                            clinical_df = clinical_cleaned_7390,
#                                            n_runs = 100,
#                                            seed = 12395 )

#saveRDS(pca_result_9_18_2, "PCA_COX_result_model-inclusive_9_18_2.rds")

# EPV=3 model = FALSE in trainning
# pca_cox_result <- readRDS("PCA_COX_result_model-inclusive.rds")

#EPV=3 model = TRUE in trainning
#pca_cox_result <- readRDS("PCA_COX_result_model-inclusive_9_18_2.rds") 

#EPV=5
#pca_cox_result <- readRDS("PCA_COX_result_model-inclusive_9_18.rds") 

calculate_coxme_predictions <- function(cox_me_model, test_df) {
  # Extract fixed effects coefficients
  fixed_coef <- fixef(cox_me_model)
  
  # Create design matrix with model variables
  design_matrix <- as.matrix(test_df[, names(fixed_coef), drop = FALSE])
  
  # Calculate linear predictor (X * beta)
  linear_predictor <- design_matrix %*% fixed_coef
  
  return(as.numeric(linear_predictor))
}

evaluate_coxme_cindex <- function(cox_me_model, test_df, time_var, status_var) {
  
  predictions <- calculate_coxme_predictions(cox_me_model, test_df)
  
  cindex <- concordance(Surv(test_df[[time_var]], test_df[[status_var]]) ~ predictions)$concordance
  return(cindex)
}

#coxme_cindex <- evaluate_coxme_cindex(cox_me_model, test_df, time_var = "t_dmfs", status_var = "e_dmfs")

predict_pca_cox <- function(new_pat, pca_cox_result) {
  # @description 
  # This function takes new patient data and a trained PCA-COX model to predict
  # survival outcomes.
  #
  # @param new_pat Dataframe containing new patient data. include:
  #   - Gene expression features (variables used in original PCA)
  #   - t_dmfs: Time to event or last follow-up
  #   - e_dmfs: Event status (0=censored, 1=event)
  #   
  # @param pca_cox_result List containing trained PCA-COX model object with:
  #   - best_results$train_pca: PCA parameters (center, scale, rotation)
  #   - best_results$max_vars_allowed: Number of principal components to use
  #   - best_results$model: Trained COX regression model
  #
  # @return Dataframe with predictions containing:
  #   - risk_score: Linear predictor score from COX model
  #   - hazard_ratio: Relative hazard ratio compared to baseline
  #   - survival_P_at_val_time: Survival probability at patient's actual time point
  #   - val_time_dmfs: Actual observation time (for validation)
  #   - val_event_dmfs: Actual event status (for validation)
  
  
  # Extract gene features (exclude survival variables)
  test_feat <- as.matrix(new_pat[, !(names(new_pat) %in% c("t_dmfs", "e_dmfs"))])
  test_surv <- as.data.frame(new_pat[,(names(new_pat) %in% c("t_dmfs", "e_dmfs"))])
  predict_time_point <- test_surv$t_dmfs
  
  # Get PCA parameters from trained model
  pca_center <- pca_cox_result$best_results$train_pca$center
  pca_scale <- pca_cox_result$best_results$train_pca$scale
  pca_rotation <- pca_cox_result$best_results$train_pca$rotation
  pca_max_vars <- pca_cox_result$best_results$max_vars_allowed
  
  # Standardize new data using training parameters
  test_scaled <- scale(test_feat, center = pca_center, scale = pca_scale)
  
  # Transform to PCA space
  test_pcs <- test_scaled %*% pca_rotation
  test_pcs <- test_pcs[, 1:pca_max_vars, drop = FALSE]
  
  # Prepare data for COX prediction
  test_df <- as.data.frame(test_pcs)
  
  # Get COX model and predict
  cox_model <- pca_cox_result$best_results$model
  
  sfit <- survfit(cox_model, newdata = test_df)  
  out <- summary(sfit, times = predict_time_point)
  
  #survival_P_at_val_time = predict(cox_model, newdata = test_df, type = "survival", times = predict_time_point)
  
  predictions <- data.frame(
    risk_score = predict(cox_model, newdata = test_df, type = "lp"),
    hazard_ratio = predict(cox_model, newdata = test_df, type = "risk"),
    survival_P_at_val_time = out$surv,
    val_time_dmfs = new_pat$t_dmfs,
    val_event_dmfs = new_pat$e_dmfs
  )
  return(predictions)
}
# results <- predict_pca_cox(new_patient_data, pca_cox_model)

sanitize_latex <- function(x) {
  gsub("([%$&#_{}~^\\\\])", "\\\\\\1", x, perl = TRUE)
}

kable_table_func <- function(df, title, digit = 4){
  # convert a table into a kable function
  title <- sanitize_latex(title)
  kable(df, digits = digit, format = "latex", caption = title) %>%
    kable_styling(
      latex_options = c("striped", "hold_position", "scale_down"),
      full_width = FALSE,
      position = "center"
    ) %>%
    row_spec(0, bold = TRUE, background = "#D3D3D3") %>%
    column_spec(2, width = "8cm") 
}
