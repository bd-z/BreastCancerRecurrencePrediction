clean_clinical_data <- function(col_selected, df_clinic_raw) {
  # selects and cleans clinical columns, standardizes column names, converts
  # values to numeric where possible, and transforms low-cardinality columns to factors.
  # # Convert months to days if all t_dmfs values are within 0 to 240
  # Remove rows where t.dmfs or e.dmfs is NA, or where the number of NAs exceeds 3
  
  # Select specified columns from the raw clinical data
  clinical_selected <- df_clinic_raw[col_selected]
  
  # Clean column names: remove ":ch1" and replace "." with "_"
  original_names <- colnames(clinical_selected)
  clean_names <- gsub(":ch1", "", original_names)
  clean_names <- gsub("time", "t", clean_names)
  clean_names <- gsub("event", "e", clean_names)
  clean_names <- gsub("\\.", "_", clean_names)
  
  colnames(clinical_selected) <- clean_names
  
  # Convert character columns to numeric where possible
  clinical_cleaned <- as.data.frame(lapply(clinical_selected, function(x) {
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    if (all(is.na(x_num))) x else x_num
  }))
  
  # Convert columns with 4 or fewer unique values to factors (starting from column 3)
  for (i in 3:ncol(clinical_cleaned)) {
    if (length(unique(clinical_cleaned[, i])) <= 4) {
      clinical_cleaned[, i] <- factor(clinical_cleaned[, i])
    }
  }
  colnames(clinical_cleaned)[1] = "geo_accession"
  
  # Remove rows from clinical_cleaned where time_dmfs or event_dmfs is NA
  clinical_cleaned <- clinical_cleaned[!is.na(clinical_cleaned$t_dmfs) 
                                       & !is.na(clinical_cleaned$e_dmfs), ]
  
  clinical_cleaned$e_dmfs <- as.integer(as.character(clinical_cleaned$e_dmfs))
  
  # Check if all t_dmfs values are within 0 to 240
  if (all(clinical_cleaned$t_dmfs >= 0 & clinical_cleaned$t_dmfs <= 240)) {
    # Convert months to days using average month length
    clinical_cleaned$t_dmfs <- as.integer(clinical_cleaned$t_dmfs * 30.4375)
    message("t_dmfs assumed to be in months and converted to days.")
  } else {
    message("t_dmfs values exceed 240; assumed not to be in months. No conversion performed.")
  }
  
  if ("size_in_cm" %in% colnames(clinical_cleaned)) {
    clinical_cleaned <- clinical_cleaned %>%
      rename(size = size_in_cm)
  }
  
  # Remove rows where t.dmfs or e.dmfs is NA, or where the number of NAs exceeds 2
  clinical_cleaned_filtered <- clinical_cleaned %>%
    filter(!is.na(t_dmfs), !is.na(e_dmfs), rowSums(is.na(.)) < 2)
  
  return(clinical_cleaned_filtered)
}


# Function to generate box plots for continuous variables
generate_box_plots <- function(data, continuous_variables) {
  # Load required libraries
  library(ggplot2)
  library(patchwork)
  
  # Create a list to hold the plots
  plot_list <- list()
  
  # Loop through each continuous variable
  for (var in continuous_variables) {
    # Remove rows with NA values for the current variable
    df <- data[!is.na(data[[var]]), ]
    
    # Create the box plot
    p <- ggplot(df, aes(y = .data[[var]])) +
      geom_boxplot(fill = "lightblue", 
                   outlier.size = 5, 
                   outlier.colour = "red", 
                   outlier.shape = 1) +
      labs(title = paste("Boxplot of", var), y = var) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 8, face = "bold"),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black")
      )
    
    # Add the plot to the list
    plot_list[[var]] <- p
  }
  
  # Combine all plots and display
  combined_plot <- wrap_plots(plot_list, ncol = 2)  # you can adjust ncol
  print(combined_plot)
}

# # Example usage
# continuous_variables <- c("size", "age", "time.rfs")
# # Assuming 'df' is your data frame
# generate_box_plots(df_clean, continuous_variables)



library(mice)
# Define the imputation function
impute_missing_value <- function(clinical_cleaned, selected_col, missed_col) {
  # Description:
  # Performs multiple imputation on selected clinical variables using the mice package.
  # Supports imputing one or more target columns with missing values.
  #
  # Input:
  # - clinical_cleaned: A cleaned clinical data frame
  # - selected_col: A character vector of columns used in the imputation model
  # - missed_col: One or more target columns (character vector) to impute
  #
  # Output:
  # - A data frame with the missing values in missed_col imputed using majority vote across multiple imputations.
  #
  # Features:
  # - Automatically assigns imputation methods:
  #     - "polr" for ordered categorical variables (e.g., grade)
  #     - "logreg" for binary variables (e.g., er)
  #     - Defaults to mice's standard method otherwise
  # - Replaces missing values in the original dataset
  # - Returns the updated dataset
  
  
  # Extract relevant variables for imputation
  df_miss <- clinical_cleaned[, selected_col]
  
  # Show missing data pattern
  md.pattern(df_miss)
  
  # Initialize imputation method for each column
  methods <- make.method(df_miss)
  
  # Assign specific imputation methods for each missed_col
  for (col in missed_col) {
    if (col %in% c("grade", "TNM_T", "TNM_N")) {
      methods[col] <- "polr"     # Ordered categorical
    } else if (col %in% c("er", "TNM_M", "Chemotherapy_Adjuvant", "MMR_Status", "KRAS_Mutation")) {
      methods[col] <- "logreg"   # Binary categorical
    } else {
      methods[col] <- ""         # Let mice choose default method
    }
  }
  
  # Perform multiple imputation
  imp <- mice(df_miss, m = 20, method = methods, seed = 123)
  
  # # Define custom mode function for majority voting
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
   
  # # # Replace missing values for each target column
  # for (col in missed_col) {
  #   if (!is.null(imp$imp[[col]]) && nrow(imp$imp[[col]]) > 0) {
  #     filled_values <- apply(imp$imp[[col]], 1, Mode)
  #     print(paste("Filling column:", col))
  #     print(filled_values)
  #     clinical_cleaned[[col]][names(filled_values)] <- filled_values
  #     #clinical_cleaned[[col]][as.numeric(names(filled_values))] <- filled_values
  #   }
  # # }
  # # 
  
  
  for (col in missed_col) {
    tbl <- imp$imp[[col]]
    # skip if no imputations for this column
    if (is.null(tbl) || nrow(tbl) == 0) next
    
    # row-wise mode across m imputations -> one value per missing row
    filled_values <- apply(tbl, 1, Mode)  # named vector; names are IDs (GSM...)
    ids  <- names(filled_values)
    
    # map IDs to row positions
    idx <- match(ids, rownames(clinical_cleaned))
    keep <- !is.na(idx) & !is.na(filled_values)
    if (!any(keep)) next
    
    vals <- filled_values[keep]
    pos  <- idx[keep]
    
    # assign safely depending on column type
    if (is.factor(clinical_cleaned[[col]])) {
      # expand levels to include new values
      new_levels <- union(levels(clinical_cleaned[[col]]), unique(vals))
      clinical_cleaned[[col]] <- factor(clinical_cleaned[[col]], levels = new_levels)
      clinical_cleaned[[col]][pos] <- vals
    } else if (is.character(clinical_cleaned[[col]])) {
      clinical_cleaned[[col]][pos] <- as.character(vals)
    } else if (is.numeric(clinical_cleaned[[col]])) {
      suppressWarnings({
        num_vals <- as.numeric(vals)
      })
      ok <- !is.na(num_vals)
      clinical_cleaned[[col]][pos[ok]] <- num_vals[ok]
    } else {
      # fallback: coerce to character, assign, then factor back if needed
      tmp <- as.character(clinical_cleaned[[col]])
      tmp[pos] <- as.character(vals)
      clinical_cleaned[[col]] <- tmp
    }
  }
  return(clinical_cleaned)
  }
  



align_expr_clin <- function(expr, clin, sample_col = 1) {
  # Description: This function aligns the columns of a gene expression matrix
  # and the rows of a clinical data frame based on shared sample names. 
  # It ensures that the sample order is consistent between the two datasets.
  s <- as.character(clin[[sample_col]])
  s <- intersect(s, colnames(expr))
  list(
    expr = expr[, s, drop = FALSE],
    clin = clin[match(s, clin[[sample_col]]), , drop = FALSE]
  )
}


combat_three_batches <- function(expr_A, clin_A,
                                 expr_B, clin_B,
                                 expr_C = NULL, clin_C = NULL,
                                 predictors = NULL) {
  # Description: this function performs batch effect correction using ComBat on up to three gene expression datasets
  # (e.g., training set A, validation set B, and an optional prediction set C), while preserving
  # important clinical covariates. It ensures that expression and clinical data are aligned by sample
  # names before applying ComBat.
  
  # Inputs:
  #   - expr_A: expression matrix for batch A (genes x samples)
  #   - clin_A: clinical data for batch A (must include sample ID in the first column)
  #   - expr_B: expression matrix for batch B
  #   - clin_B: clinical data for batch B
  #   - expr_C (optional): expression matrix for batch C (e.g., new test data)
  #   - clin_C (optional): clinical data for batch C
  #   - predictors (optional): character vector of clinical covariate names to preserve during ComBat (e.g., c("age", "grade"))
  #
  # Output:
  #   A list containing corrected expression matrices:
  #     - $A : ComBat-corrected expression matrix for batch A
  #     - $B : ComBat-corrected expression matrix for batch B
  #     - $C : ComBat-corrected expression matrix for batch C (if provided)
  #
  # Note:
  #   - The function ensures expression and clinical data are aligned by sample name before correction.
  #   - Input expression matrices must have sample names as column names.
  #   - Clinical data must contain matching sample names in the first column.
  
  
  
  
  # Align expression and clinical data for batch A
  aligned_A <- align_expr_clin(expr_A, clin_A)
  expr_A_aligned <- aligned_A$expr
  clin_A_aligned <- aligned_A$clin
  
  # Align expression and clinical data for batch B
  aligned_B <- align_expr_clin(expr_B, clin_B)
  expr_B_aligned <- aligned_B$expr
  clin_B_aligned <- aligned_B$clin
  
  # Align expression and clinical data for batch C if provided
  if (!is.null(expr_C)) {
    aligned_C <- align_expr_clin(expr_C, clin_C)
    expr_C_aligned <- aligned_C$expr
    clin_C_aligned <- aligned_C$clin
  }
  
  # Sanity checks: number of samples must match between expr and clin
  stopifnot(ncol(expr_A_aligned) == nrow(clin_A_aligned))
  stopifnot(ncol(expr_B_aligned) == nrow(clin_B_aligned))
  if (!is.null(expr_C)) stopifnot(ncol(expr_C_aligned) == nrow(clin_C_aligned))
  
  # Create batch labels
  n_A <- ncol(expr_A_aligned)
  n_B <- ncol(expr_B_aligned)
  batch_A <- rep("A", n_A)
  batch_B <- rep("B", n_B)
  
  if (!is.null(expr_C)) {
    n_C <- ncol(expr_C_aligned)
    batch_C <- rep("C", n_C)
  }
  
  # Combine all expression and clinical data
  expr_all <- cbind(expr_A_aligned, expr_B_aligned, if (!is.null(expr_C)) expr_C_aligned)
  clin_all <- rbind(clin_A_aligned, clin_B_aligned, if (!is.null(expr_C)) clin_C_aligned)
  batch_all <- c(batch_A, batch_B, if (!is.null(expr_C)) batch_C)
  
  # Ensure sample order matches
  stopifnot(identical(colnames(expr_all), clin_all$geo_accession))
  
  # Build model matrix for covariates to preserve (if provided)
  if (!is.null(predictors)) {
    mod <- model.matrix(as.formula(
      paste("~", paste(predictors, collapse = "+"))
    ), data = clin_all)
  } else {
    mod <- NULL
  }
  
  # Apply ComBat batch correction (no transpose needed if input is genes x samples)
  expr_corrected <- ComBat(dat = expr_all, batch = batch_all, mod = mod)
  
  # Split back to original datasets
  expr_A_corrected <- expr_corrected[, 1:n_A ]
  expr_B_corrected <- expr_corrected[,(n_A + 1):(n_A + n_B)]
  if (!is.null(expr_C)) {
    expr_C_corrected <- expr_corrected[,(n_A + n_B + 1):(n_A + n_B + n_C)]
  } else {
    expr_C_corrected <- NULL
  }
  
  return(list(
    A = expr_A_corrected,
    B = expr_B_corrected,
    C = expr_C_corrected
  ))
}



# 
# 
# combat_three_batches <- function(expr_A, clin_A,
#                                  expr_B, clin_B,
#                                  expr_C = NULL, clin_C = NULL,
#                                  predictors = NULL) {
#   
#   
#   aligned_A <- align_expr_clin(expr_A, clin_A)
#   expr_A_aligned <- aligned_A$expr
#   clin_A_aligned <- aligned_A$clin
#   
#   aligned_B <- align_expr_clin(expr_B, clin_B)
#   expr_B_aligned <- aligned_B$expr
#   clin_B_aligned <- aligned_B$clin
#   
#   if (!is.null(expr_C)){
#     aligned_C <- align_expr_clin(expr_C, clin_C)
#     expr_C_aligned <- aligned_C$expr
#     clin_C_aligned <- aligned_C$clin
#   }
#   
#   
#   
#   stopifnot(ncol(expr_A_aligned) == nrow(clin_A_aligned))
#   stopifnot(ncol(expr_B_aligned) == nrow(clin_B_aligned))
#   if (!is.null(expr_C)) stopifnot(ncol(expr_C_aligned) == nrow(clin_C_aligned))
#   
#   # 添加 batch 标签
#   n_A <- ncol(expr_A_aligned)
#   n_B <- ncol(expr_B_aligned)
#   batch_A <- rep("A", n_A)
#   batch_B <- rep("B", n_B)
#   
#   if (!is.null(expr_C)) {
#     n_C <- ncol(expr_C_aligned)
#     batch_C <- rep("C", n_C)
#   }
#   
#   # 合并表达数据和临床数据
#   expr_all <- cbind(expr_A_aligned, expr_B_aligned, if (!is.null(expr_C)) expr_C_aligned)
#   clin_all <- rbind(clin_A_aligned, clin_B_aligned, if (!is.null(expr_C)) clin_C_aligned)
#   batch_all <- c(batch_A, batch_B, if (!is.null(expr_C)) batch_C)
#   stopifnot(identical(colnames(expr_all), clin_all$geo_accession))
#   
#   
#   # 提取表达部分（基因列）
#   # gene_cols <- grep("_at$", colnames(expr_all), value = TRUE)
#   # expr_mat <- expr_all[, gene_cols]
#   # rownames(expr_mat) <- rownames(expr_all)
#   # 
#   # 构造 mod 矩阵（如有）
#   if (!is.null(predictors)) {
#     mod <- model.matrix(as.formula(
#       paste("~", paste(predictors, collapse = "+"))
#     ), data = clin_all)
#   } else {
#     mod <- NULL
#   }
#   
#   # 转置表达矩阵（基因 × 样本） → ComBat 要求
#   #expr_mat_t <- t(as.matrix(expr_mat))
#   
#   # 执行 ComBat
#   expr_corrected <- ComBat(dat = expr_all, batch = batch_all, mod = mod)
#   #expr_corrected <- t(expr_corrected_t)  # 转回样本 × 基因
#   
#   # 拆分回 A, B, C
#   expr_A_corrected <- expr_corrected[, 1:n_A ]
#   expr_B_corrected <- expr_corrected[,(n_A + 1):(n_A + n_B)]
#   if (!is.null(expr_C)) {
#     expr_C_corrected <- expr_corrected[,(n_A + n_B + 1):(n_A + n_B + n_C)]
#   } else {
#     expr_C_corrected <- NULL
#   }
#   
#   return(list(
#     A = expr_A_corrected,
#     B = expr_B_corrected,
#     C = expr_C_corrected
#   ))
# }
# 







# Function to build and test a Cox proportional hazards model
fit_cox_model <- function(predictors, df) {
  tryCatch({
    # 构建公式
    formula <- as.formula(paste(
      "Surv(t_dmfs, e_dmfs) ~",
      paste(predictors, collapse = " + ")
    ))
    
    # 拟合 Cox 模型
    cox_model <- coxph(formula, data = df, x = TRUE, y = TRUE)
    
    # 模型摘要
    cat("Cox model summary:\n")
    model_summary <- summary(cox_model)
    print(model_summary)
    
    # 系数表
    coef_df <- data.frame(
      variable = rownames(model_summary$coefficients),
      coef = model_summary$coefficients[, "coef"],
      zscore = model_summary$coefficients[, "z"],
      pvalue = model_summary$coefficients[, "Pr(>|z|)"]
    )
    
    # PH 假设检验
    test_ph <- cox.zph(cox_model)
    print("Test proportional hazards assumption:")
    print(test_ph)
    
    # 返回固定结构
    list(
      model = cox_model,
      ph_test = test_ph,
      coef_table = coef_df
    )
  }, error = function(e) {
    message("Error in fit_cox_model: ", e$message)
    # 返回固定结构，值为 NULL
    list(
      model = NULL,
      ph_test = NULL,
      coef_table = NULL
    )
  })
}





# fit_cox_model <- function(predictors, df) {
#   # Build formula: Surv(...) ~ var1 + var2 + ...
#   formula <- as.formula(paste(
#     "Surv(t_dmfs, e_dmfs) ~",
#     paste(predictors, collapse = " + ")
#   ))
#   
#   # Fit Cox proportional hazards model
#   cox_model <- coxph(formula, data = df, x = TRUE, y = TRUE)
#   
#   # Print model summary
#   cat("Cox model summary:\n")
#   model_summary <- summary(cox_model)
#   print(model_summary)
#   
#   # Extract coefficients, z-score, and p-value
#   coef_df <- data.frame(
#     variable = rownames(model_summary$coefficients),
#     coef = model_summary$coefficients[, "coef"],
#     zscore = model_summary$coefficients[, "z"],
#     pvalue = model_summary$coefficients[, "Pr(>|z|)"]
#   )
#   
#   # Proportional hazards assumption test
#   test_ph <- cox.zph(cox_model)
#   print("Test proportional hazards assumption:")
#   print(test_ph)
#   
#   # Plot Schoenfeld residuals test
#   #print(ggcoxzph(test_ph))
#   
#   # Return everything
#   return(list(
#     model = cox_model,
#     ph_test = test_ph,
#     coef_table = coef_df
#   ))
# }




library(survival)
library(coxphf)

fit_cox_model_firth <- function(predictors, df, coef_max = 15) {
  # 构建公式
  formula <- as.formula(paste(
    "Surv(t_dmfs, e_dmfs) ~",
    paste(predictors, collapse = " + ")
  ))
  
  # 先尝试普通 Cox 拟合
  cox_model <- tryCatch({
    coxph(formula, data = df, x = TRUE, y = TRUE)
  }, error = function(e) NULL)
  
  use_firth <- FALSE
  
  # 检查是否需要切换到 Firth Cox
  if (is.null(cox_model)) {
    message("Cox model failed, switching to Firth Cox.")
    use_firth <- TRUE
  } else {
    # 检查系数是否发散
    if (any(abs(coef(cox_model)) > coef_max, na.rm = TRUE)) {
      message("Coefficient > ", coef_max, " detected, switching to Firth Cox.")
      use_firth <- TRUE
    }
    # 检查 Concordance 是否为 1（完全分离）
    concordance_val <- tryCatch({
      suppressWarnings(summary(cox_model)$concordance[1])
    }, error = function(e) NA)
    if (!is.na(concordance_val) && concordance_val >= 0.98) {
      message("Concordance = 1 detected, switching to Firth Cox.")
      use_firth <- TRUE
    }
  }
  
  # 如果触发 Firth Cox
  if (use_firth) {
    cox_model <- tryCatch({
      coxphf(formula, data = df)
    }, error = function(e) NULL)
    if (is.null(cox_model)) {
      warning("Firth Cox model also failed. Returning NULL.")
      return(NULL)
    }
    model_summary <- summary(cox_model)
  } else {
    model_summary <- summary(cox_model)
  }
  
  # 输出模型信息
  cat("Cox model summary:\n")
  print(model_summary)
  
  # 提取系数信息（Firth 和普通 Cox 格式略不同，这里做兼容）
  if ("coefficients" %in% names(model_summary)) {
    coef_df <- data.frame(
      variable = rownames(model_summary$coefficients),
      coef = model_summary$coefficients[, "coef"],
      zscore = model_summary$coefficients[, "z"],
      pvalue = model_summary$coefficients[, "Pr(>|z|)"]
    )
  } else {
    coef_df <- data.frame(
      variable = rownames(model_summary$coeff),
      coef = model_summary$coeff[, "coef"],
      pvalue = model_summary$coeff[, "p"],
      zscore = NA
    )
  }
  
  # PH 假设检验（只有普通 Cox 才能做）
  test_ph <- NULL
  if (!use_firth) {
    test_ph <- cox.zph(cox_model)
    print("Test proportional hazards assumption:")
    print(test_ph)
  }
  
  return(list(
    model = cox_model,
    ph_test = test_ph,
    coef_table = coef_df
  ))
}

# fit_cox_model <- function(predictors, df) {
#   
#   # Construct formula dynamically: Surv(...) ~ predictor1 + predictor2 + ...
#   formula <- as.formula(paste(
#     "Surv(t_dmfs, e_dmfs) ~",
#     paste(predictors, collapse = " + ")
#   ))
#   
#   # Fit the Cox proportional hazards model
#   cox_model <- coxph(formula, data = df, x = TRUE, y = TRUE)
#   
#   # Print model summary
#   print("cox model summary")
#   print(summary(cox_model))
#   
#   # Test proportional hazards assumption
#   test_ph <- cox.zph(cox_model)
#   print(ggcoxzph(test_ph))
#   
#   print("Test proportional hazards assumption:")
#   print(test_ph)
#   
#   # Return model and PH test result
#   return(list(model = cox_model,
#               ph_test = test_ph,
#               model_summary = summary(cox_model)))
# }

# Function to plot Kaplan-Meier survival curve by a target grouping variable
plot_km_by_group <- function(df, group_var) {
  # Convert group_var to factor
  df$group_var <- as.factor(df[[group_var]])
  
  # Fit survival model
  fit <- survfit(Surv(t_dmfs, e_dmfs) ~ group_var, data = df)
  
  # Create survival plot
  ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    xlab = "Days",
    ylab = "Survival Probability",
    palette = "Set1"
  )
}




calculate_time_auc_cindex <- function(model_type = c("Cox", "RSF"), fitted_model, df) {
  # This function evaluates the discrimination performance of a fitted survival 
  # model (either Cox or Random Survival Forest) on a given dataset. It performs
  # the following:
  # Computes time-dependent AUC using the Score() function.
  # Plot time-dependent AUC with smoothed event counts as a secondary axis using
  # time windows.
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
  
  # compute time-dependent AUC
  auc_result <- Score(
    object = list(model = fitted_model),
    formula = Surv(t_dmfs, e_dmfs) ~ 1,
    data = df,
    metrics = "AUC",
    summary = "AUC",
    times = time_points,
    conf.int = TRUE,
    plots = TRUE,
    split.method = "none"  # use external data, no cross-validation
  )
  
  # extract AUC values
  # 提取 AUC 数据框
  auc_values <- auc_result$AUC$score
  auc_df <- data.frame(Time = auc_values$times, AUC = auc_values$AUC)
  
  # Bin event counts over small time intervals
  bin_width <- diff(range(time_points)) / length(time_points)
  event_counts <- sapply(time_points, function(t) {
    sum(df$e_dmfs == 1 & df$t_dmfs >= (t - bin_width/2) & df$t_dmfs < (t + bin_width/2))
  })
  
  # # Plot AUC with event count as secondary axis
  # par(mar = c(5, 4, 4, 4) + 0.3)  # Leave room for right axis
  # plot(auc_df$Time, auc_df$AUC, type = "l", col = "green", lwd = 2,
  #      xlab = "Time", ylab = "AUC",
  #      main = paste("Time-Dependent AUC and Event Count (", model_type, ")", sep = ""))
  # grid()
  # 
  # par(new = TRUE)
  # plot(auc_df$Time, event_counts, type = "l", col = "red", lwd = 2, axes = FALSE,
  #      xlab = "", ylab = "", lty = 2)
  # axis(side = 4, col.axis = "red", col = "red")
  # mtext("Number of Events", side = 4, line = 3, col = "red")
  # legend("bottomleft", legend = c("AUC", "Event Count"),
  #        col = c("green", "red"), lty = c(1, 2), lwd = 2, bty = "n")
   
  # compute time-weighted average AUC (iAUC)
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

library(randomForestSRC)
library(survival)
library(riskRegression)

rsf_kfold_cv_best <- function(data, K = 5, ntree = 1000, seed = NULL) {
  # Performs K-fold cross-validation for a Random Survival Forest model.
  # Retrain the final model on the entire dataset to extract variable importance.
  # Returns the model with the highest C-index on the validation fold.
  # Output:
  #   Best model
  #   C-index per fold
  #   Best fold number
  #   Variable importance
   
  if (!is.null(seed)) set.seed
  
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
  final_model <- rfsrc(Surv(t_dmfs, e_dmfs) ~ ., data = clinical_rsf,
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





compute_risk_score <- function(gene_mat_scaled, significant_vars_df, clinical_cleaned, n_group = 3) {
  # Computes risk scores based on scaled gene expression and model coefficients,
  # and assigns each sample to a risk group (e.g., low, medium, high) based on quantiles.
  # 
  # Arguments:
  # - gene_mat_scaled: matrix of scaled gene expression values (genes x samples)
  # - significant_vars_df: data frame with selected gene names (rownames) and coefficients
  # - clinical_cleaned: data frame with clinical data (samples as rows)
  # - n_group: number of risk groups to divide samples into (default = 3)
  #
  # Returns:
  # - clinical_cleaned with added columns: 'risk_score', 'risk_group' and gene data
  
  # if risk_score/risk_group column already exists, rename it.
  names(clinical_cleaned)[names(clinical_cleaned) == "risk_score"] <- "lasso_risk_score"
  names(clinical_cleaned)[names(clinical_cleaned) == "risk_group"] <- "lasso_risk_group"
  
  genes <- sub("^ge_", "", rownames(significant_vars_df))
  # Keep only the samples in expr_2990_scaled that exist in clinical_cleaned$geo_accession
  gene_mat_scaled <- gene_mat_scaled[, match(clinical_cleaned$geo_accession, colnames(gene_mat_scaled))]
  gene_mat_scaled <- gene_mat_scaled[genes, ]
  
  if (!identical(rownames(gene_mat_scaled), rownames(significant_vars_df)))
    stop("Row names are not identical or not in the same order.")
  
  clinical_cleaned$risk_score <- as.vector(t(gene_mat_scaled) %*% significant_vars_df$coef)
  
  labels <- paste0("Risk", seq_len(n_group))
  
  clinical_cleaned <- clinical_cleaned %>%
    mutate(risk_group = ntile(risk_score, n_group),
           risk_group = factor(risk_group, labels = labels))
  
  # add gene data
  gene_mat_scaled_t <- t(gene_mat_scaled)
  gene_df <- as.data.frame(gene_mat_scaled_t)
  
  # Add prefix "ge_" to column names starting with a digit
  colnames(gene_df)[grepl("^[0-9]", colnames(gene_df))] <- 
    paste0("ge_", colnames(gene_df)[grepl("^[0-9]", colnames(gene_df))])
  
  gene_df$geo_accession <- rownames(gene_df)
  merged_df <- merge(clinical_cleaned, gene_df, by = "geo_accession", sort = FALSE)
  merged_df
}


split_expr_clinical <- function(expr_mat, clinical_df, 
                                stratify_col = "e_dmfs", train_frac = 0.7, seed = 123) {
  # Desprition: 该方程会自动将基因样本与医院样本对齐
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Stratified sampling based on the specified column (e.g., event indicator)
  train_index <- createDataPartition(clinical_df[[stratify_col]], p = train_frac, list = FALSE)
  
  # Split clinical data into training and testing sets
  train_clinical <- clinical_df[train_index, ]
  test_clinical  <- clinical_df[-train_index, ]
  
  # Get sample IDs (assumes geo_accession is the sample identifier)
  train_samples <- train_clinical$geo_accession
  test_samples  <- test_clinical$geo_accession
  
  # Subset expression matrix by sample IDs
  train_expr <- expr_mat[, train_samples]
  test_expr  <- expr_mat[, test_samples]
  
  # Return a named list containing all splits
  return(list(
    train_expr = train_expr,
    test_expr = test_expr,
    train_clinical = train_clinical,
    test_clinical = test_clinical
  ))
}




batch_univariate_cox_regression <- function(train_expr, train_clinical) {
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
  sig_genes_df <- subset(cox_df, p.value < 0.05)
  
  return(sig_genes_df)
}

# Example usage:
# sig_genes_df <- batch_univariate_cox_regression(train_expr, train_clinical)

lasso_cox_cv <- function(train_expr, train_clinical, sig_gene_df) {
  # Lasso Cox Regression with Cross-Validation
  # Description: Runs Lasso Cox with 10-fold CV to select significant genes.
  # Input: train_expr, train_clinical, sig_gene_df
  # Output: selected_gene_df
  
  sig_genes <- sig_gene_df$gene
  expr_transposed <- t(train_expr[sig_genes, ])
  X <- expr_transposed
  y <- Surv(train_clinical$t_dmfs, train_clinical$e_dmfs)
  
  set.seed(123)
  cvfit <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
  
  # plot(cvfit)
  # abline(v = log(cvfit$lambda.min), col = "red", lty = 2)   
  # abline(v = log(cvfit$lambda.1se), col = "blue", lty = 2)  
  
  coef_opt <- coef(cvfit, s = "lambda.min")
  selected_genes <- as.matrix(coef_opt)
  selected_gene_df <- data.frame(gene = rownames(selected_genes), coef = selected_genes[, 1])
  selected_gene_df <- selected_gene_df[selected_gene_df$coef != 0, ]
  
  return(selected_gene_df)
}

# Example: selected_gene_df <- lasso_cox_cv(train_expr, train_clinical, sig_gene_df)



cox_rsf_workflow <- function(gene_expr,
                         clinical_data_imputed,
                         train_frac,
                         seed,
                         clin_pred = TRUE,
                         riskscore_pred = FALSE,
                         gene_pred = TRUE,
                         clin_predictors = c("grade", "er", "age", "size")
                         ) {
  # Cox Workflow Function
  # Description: Executes a full Cox model workflow including data splitting, gene selection, risk scoring, and model evaluation.
  # Input: gene_expr (expression matrix), clinical_data_imputed (imputed clinical data), train_frac (training fraction), seed (random seed)
  # Output: results_train (Cox model results), result_valid (validation metrics)
  
  
  # Split data into train and test sets
  result <- split_expr_clinical(gene_expr, clinical_data_imputed, stratify_col = "e_dmfs", train_frac = train_frac, seed = seed)
  train_expr <- result$train_expr
  train_clinical <- result$train_clinical
  test_expr <- result$test_expr
  test_clinical <- result$test_clinical
  
  # Perform batch univariate Cox regression
  sig_gene_df <- batch_univariate_cox_regression(train_expr, train_clinical)
  
  # Apply Lasso Cox with cross-validation
  selected_gene_df <- lasso_cox_cv(train_expr, train_clinical, sig_gene_df)
  
  if (nrow(selected_gene_df) == 0) {
    message("No genes selected, skipping this run.")
    return(NULL)
  } else {
  
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
  
  # Fit Cox model on training data
  results_train <- fit_cox_model(predictors, clinical_cleaned_risk_train)
  
  # Compute risk scores for test data
  clinical_cleaned_risk_test <- compute_risk_score(gene_mat_scaled = test_expr,
                                                   significant_vars_df = selected_gene_df,
                                                   clinical_cleaned = test_clinical,
                                                   n_group = 3)
  
  # Calculate validation metrics
  result_valid <- calculate_time_auc_cindex("Cox", fitted_model = results_train$model, df = clinical_cleaned_risk_test)
  
  # random forest
  clinical_rsf <- clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs", predictors)]
  
  # build RSF model
  result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5)
  # Best model
  rsf_fit_best <- result_rsf_train$best_model
  
  #Remove variables with negative/less importance from the predictors vector based on RSF model output
  imp <- result_rsf_train$importance
  vars_to_remove <- names(imp)[imp < 0.01]
  predictors_filtered_rsf <- setdiff(predictors, vars_to_remove)
  
  clinical_rsf <- clinical_cleaned_risk_train[, c("t_dmfs", "e_dmfs", predictors_filtered_rsf)]
  # build RSF model with predictors_filtered_rsf
  result_rsf_train <- rsf_kfold_cv_best(clinical_rsf, K = 5)
  # Best model
  rsf_fit_best <- result_rsf_train$best_model
  
  result_rsf_valid <- calculate_time_auc_cindex("RSF", fitted_model = rsf_fit_best, df = clinical_cleaned_risk_test)
  
  return(list(results_train = results_train,
              result_valid = result_valid,
              result_rsf_train = result_rsf_train,
              result_rsf_valid = result_rsf_valid,
              predictors_filtered_rsf
              ))
  }
}

# Example usage:
# result <- cox_workflow(gene_expr, clinical_data_imputed, train_frac = 0.7, seed = 345)
