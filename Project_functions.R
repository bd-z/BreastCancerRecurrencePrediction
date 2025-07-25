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
  
  # Extract relevant variables for imputation
  df_miss <- clinical_cleaned[, selected_col]
  
  # Show missing data pattern (optional, informative only)
  md.pattern(df_miss)
  
  # Specify imputation methods
  methods <- make.method(df_miss)
  
  # Assign appropriate method based on the variable type
  if (missed_col == "grade") {
    methods["grade"] <- "polr"  # Ordered categorical: proportional odds logistic regression
  } else if (missed_col == "er") {
    methods["er"] <- "logreg"  # Binary categorical: logistic regression
  }
  
  # Perform multiple imputation (default m = 5)
  imp <- mice(df_miss, m = 5, method = methods, seed = 123)
  
  # Define custom mode function
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Use majority vote (mode) across imputations for the target variable
  filled_values <- apply(imp$imp[[missed_col]], 1, Mode)
  print(filled_values)
  
  # Replace missing values in the original dataset
  clinical_cleaned[[missed_col]][as.numeric(names(filled_values))] <- filled_values
  # Return updated dataset
  return(clinical_cleaned)
}

# Function to build and test a Cox proportional hazards model
fit_cox_model <- function(predictors, df) {
  
  # Construct formula dynamically: Surv(...) ~ predictor1 + predictor2 + ...
  formula <- as.formula(paste(
    "Surv(t_dmfs, e_dmfs) ~",
    paste(predictors, collapse = " + ")
  ))
  
  # Fit the Cox proportional hazards model
  cox_model <- coxph(formula, data = df, x = TRUE, y = TRUE)
  
  # Print model summary
  print("cox model summary")
  print(summary(cox_model))
  
  # Test proportional hazards assumption
  test_ph <- cox.zph(cox_model)
  print(ggcoxzph(test_ph))
  
  print("Test proportional hazards assumption:")
  print(test_ph)
  
  # Return model and PH test result
  return(list(model = cox_model, ph_test = test_ph))
}

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



library(survival)
library(randomForestSRC)
library(riskRegression)

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
  
  # Plot AUC with event count as secondary axis
  par(mar = c(5, 4, 4, 4) + 0.3)  # Leave room for right axis
  plot(auc_df$Time, auc_df$AUC, type = "l", col = "green", lwd = 2,
       xlab = "Time", ylab = "AUC",
       main = paste("Time-Dependent AUC and Event Count (", model_type, ")", sep = ""))
  grid()
  
  par(new = TRUE)
  plot(auc_df$Time, event_counts, type = "l", col = "red", lwd = 2, axes = FALSE,
       xlab = "", ylab = "", lty = 2)
  axis(side = 4, col.axis = "red", col = "red")
  mtext("Number of Events", side = 4, line = 3, col = "red")
  legend("bottomleft", legend = c("AUC", "Event Count"),
         col = c("green", "red"), lty = c(1, 2), lwd = 2, bty = "n")
   
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

rsf_kfold_cv_best <- function(data, K = 5, ntree = 1000, seed = 123) {
  # Performs K-fold cross-validation for a Random Survival Forest model.
  # Retrain the final model on the entire dataset to extract variable importance.
  # Returns the model with the highest C-index on the validation fold.
  # Output:
  #   Best model
  #   C-index per fold
  #   Best fold number
  #   Variable importance
   
  set.seed(seed)
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



compute_risk_score <- function(gene_mat_scaled, significant_vars, clinical_cleaned, n_group = 3) {
  # Computes risk scores based on scaled gene expression and model coefficients,
  # and assigns each sample to a risk group (e.g., low, medium, high) based on quantiles.
  # 
  # Arguments:
  # - gene_mat_scaled: matrix of scaled gene expression values (genes x samples)
  # - significant_vars: data frame with selected gene names (rownames) and coefficients
  # - clinical_cleaned: data frame with clinical data (samples as rows)
  # - n_group: number of risk groups to divide samples into (default = 3)
  #
  # Returns:
  # - clinical_cleaned with added columns: 'risk_score' and 'risk_group'
  
  genes <- sub("^ge_", "", rownames(significant_vars))
  # Keep only the samples in expr_2990_scaled that exist in clinical_cleaned$geo_accession
  gene_mat_scaled <- gene_mat_scaled[, colnames(gene_mat_scaled) %in% clinical_cleaned$geo_accession]
  clinical_cleaned$risk_score <- as.vector(t(gene_mat_scaled[genes, ]) %*% significant_vars$coef)
  
  labels <- paste0("Risk", seq_len(n_group))
  
  clinical_cleaned <- clinical_cleaned %>%
    mutate(risk_group = ntile(risk_score, n_group),
           risk_group = factor(risk_group, labels = labels))
  
  clinical_cleaned
}

