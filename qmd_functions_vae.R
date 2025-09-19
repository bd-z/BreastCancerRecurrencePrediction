#### VAE COX
vae_module <- nn_module(
  "VAE",
  initialize = function(input_dim = input_dim, latent_dim = latent_dim, hidden_dim = hidden_dim) {
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$fc_mean <- nn_linear(hidden_dim, latent_dim)
    self$fc_logvar <- nn_linear(hidden_dim, latent_dim)
    self$fc2 <- nn_linear(latent_dim, hidden_dim)
    self$output <- nn_linear(hidden_dim, input_dim)
  },
  
  reparameterize = function(mu, logvar) {
    std <- torch_exp(0.5 * logvar)
    eps <- torch_randn_like(std)
    mu + eps * std
  },
  
  forward = function(x) {
    h <- torch_relu(self$fc1(x))
    mu <- self$fc_mean(h)
    logvar <- self$fc_logvar(h)
    z <- self$reparameterize(mu, logvar)
    h_dec <- torch_relu(self$fc2(z))
    x_recon <- self$output(h_dec)
    list(x_recon = x_recon, mu = mu, logvar = logvar)
  }
)

loss_fn <- function(recon_x, x, mu, logvar) {
  recon_loss <- nnf_mse_loss(recon_x, x, reduction = "sum")
  kl_loss <- -0.5 * torch_sum(1 + logvar - mu^2 - torch_exp(logvar))
  recon_loss + kl_loss
}


find_best_latent_dim <- function(
    latent_dims,
    x_data,
    x_val_data,
    hidden_dim,
    lr = 0.001,
    max_epochs = 800,
    patience_train = 100,
    min_delta = 1e-6
) {
  
  input_dim <- ncol(x_data)   
  overall_best_val_loss <- Inf
  overall_best_latent_dim <- NA
  results <- list()
  
  for (latent_dim in latent_dims) {
    cat(sprintf("\n--- latent_dim = %d ---\n", latent_dim))
    
    model <- vae_module(input_dim, latent_dim, hidden_dim)
    optimizer <- optim_adam(model$parameters, lr = lr)
    
    best_train_loss <- Inf
    no_improve_count <- 0
    final_epoch <- max_epochs
    
    # training loop (early stop on train_loss)
    for (epoch in 1:max_epochs) {
      model$train(); optimizer$zero_grad()
      out <- model(x_data)
      loss <- loss_fn(out$x_recon, x_data, out$mu, out$logvar)
      loss$backward(); optimizer$step()
      
      train_loss <- loss$item()
      
      if (train_loss < best_train_loss - min_delta) {
        best_train_loss <- train_loss
        no_improve_count <- 0
      } else {
        no_improve_count <- no_improve_count + 1
      }
      if (no_improve_count >= patience_train) {
        final_epoch <- epoch
        break
      }
    }
    
    # validation after training
    model$eval()
    with_no_grad({
      val_out <- model(x_val_data)
      val_loss <- loss_fn(val_out$x_recon, x_val_data, val_out$mu, val_out$logvar)$item()
    })
    
    cat(sprintf("Finished latent_dim=%d | val_loss=%.6f | epoch=%d\n",
                latent_dim, val_loss, final_epoch))
    
    # store result
    results[[length(results) + 1]] <- list(
      latent_dim = latent_dim,
      best_val_loss = val_loss
    )
    
    # update best model
    if (val_loss < overall_best_val_loss) {
      overall_best_val_loss <- val_loss
      overall_best_latent_dim <- latent_dim
      torch_save(model$state_dict(), "best_model.pt")
    }
  }
  
  return(list(
    results = results,
    best_latent_dim = overall_best_latent_dim,
    best_overall_loss = overall_best_val_loss
  ))
}

# 
# best_2 <- find_best_latent_dim(
#   latent_dims = c(51, 61, 71, 81, 91, 101, 111, 121, 150, 200),
#   x_data = x_data,
#   x_val_data = x_val_data,
#   hidden_dim = 128
# )
# 
# results <- best_2$results
# best_latent_dim <- best_2$best_latent_dim
# best_overall_loss <- best_2$best_overall_loss
# 
# saveRDS(best_2, "best_latent_dim150.rds")

# best_3 <- find_best_latent_dim(
#   latent_dims = c(150, 160, 170, 200, 300, 400, 500, 100),
#   x_data = x_data,
#   x_val_data = x_val_data,
#   hidden_dim = 128
# )
# 
# results <- best_3$results
# 
# # Plot
# if (length(results) > 0) {
#   latent_dims <- sapply(results, function(x) x$latent_dim)
#   val_losses <- sapply(results, function(x) x$best_val_loss)
#   
#   plot(latent_dims, val_losses, type = "b", pch = 19, col = "blue",
#        xlab = "Latent Dimension", ylab = "Best Validation Loss",
#        main = "VAE Performance vs Latent Dimension")
#   points(best_latent_dim, best_overall_loss, pch = 19, col = "red", cex = 1.5)
#   legend("topright", legend = "Best", pch = 19, col = "red")
# }

extract_latent_variables <- function(model=model1, data=x_data, device) {
  # Extract latent variables mu from a VAE model (Cox model uses only mu)  
  # Arguments:
  # model: Trained VAE model
  # data: Input data matrix/data frame, should be tensor
  
  # Returns:
  # Latent variable mu matrix 
  
  # Convert to tensor
  #data_tensor <- torch_tensor(as.matrix(data), dtype = torch_float32())$to(device = device)
  
  model$eval()
  model$to(device = device)
  
  with_no_grad({
    output <- model(data)
    train_latent_mu <- as.array(output$mu$cpu())
  })
  
  return(as.data.frame(train_latent_mu))
}

run_vae_cox_iterations <- function(
    expr,                    # Gene expression matrix (genes x samples)
    clinical_df,             # Clinical data (must include t_dmfs, e_dmfs)
    n_iterations   = 50,
    seed           = 1234,
    latent_dim     = 11,
    hidden_dim     = 256,
    epochs         = 800,
    patience       = 50,
    lr             = 0.001,
    mad_quantile   = 0.25,
    uni_cox_p      = 0.01,
    corr_cutoff    = 0.90,
    save_dir       = NULL     # e.g. "best_models"; NULL means do not save
) {
  
  
  # --- Device setup ---
  device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
  
  # --- Initialization ---
  best_score        <- -Inf
  best_vae_model    <- NULL
  best_multiv_model <- NULL
  best_latent_dim   <- latent_dim
  best_iteration    <- 0
  iteration_results <- list()
  
  expr_mat     <- expr
  clinical_tbl <- clinical_df
  
  for (iteration in seq_len(n_iterations)) {
    set.seed(seed + iteration)
    message(sprintf("\n=== Iteration %d/%d ===\n", iteration, n_iterations))
    
    # 1) Bootstrap sampling + train/test split + missing value imputation
    data_prepare_result <- bootstrap_sampl_split_miss_vale_imputation(expr_mat, clinical_tbl)
    train_expr     <- data_prepare_result$train_expr
    test_expr      <- data_prepare_result$test_expr
    train_clinical <- data_prepare_result$train_clinical
    test_clinical  <- data_prepare_result$test_clinical
    
    # 2) MAD filtering
    filtered_data <- filter_genes_by_mad(train_expr, test_expr, quantile_cutoff = mad_quantile)
    train_expr2   <- filtered_data$train_expr
    test_expr2    <- filtered_data$test_expr
    
    # 3) Univariate Cox regression for gene selection
    sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = uni_cox_p)
    significant_gene <- sig_gene_df$gene
    message(sprintf("Number of significant genes = %d", length(significant_gene)))
    
    # 4) Standardization based on training data + removal of highly correlated genes
    train_expr_scaled <- standardize_with_train(
      gene_mat_train = train_expr2,
      gene_mat_valid = train_expr2,
      significant_gene = significant_gene
    )
    train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = corr_cutoff)
    significant_gene2 <- rownames(train_expr_filtered)
    
    test_expr_scaled <- standardize_with_train(
      gene_mat_train = train_expr2,
      gene_mat_valid = test_expr2,
      significant_gene = significant_gene2
    )
    
    # 5) Transpose to (samples x features)
    train_expr_df <- t(train_expr_filtered)
    test_expr_df  <- t(test_expr_scaled)
    input_dim     <- ncol(train_expr_df)
    
    # 6) Convert to torch tensors
    x_data     <- torch_tensor(train_expr_df, dtype = torch_float())$to(device = device)
    x_val_data <- torch_tensor(test_expr_df,  dtype = torch_float())$to(device = device)
    
    # 7) Define and train VAE
    model <- vae_module(
      input_dim   = input_dim,
      latent_dim  = latent_dim,
      hidden_dim  = hidden_dim
    )$to(device = device)
    
    optimizer <- optim_adam(model$parameters, lr = lr)
    
    best_train_loss  <- Inf
    best_epoch       <- 0
    patience_counter <- 0
    
    for (epoch in seq_len(epochs)) {
      model$train()
      optimizer$zero_grad()
      
      output <- model(x_data)
      loss   <- loss_fn(output$x_recon, x_data, output$mu, output$logvar)
      loss$backward()
      optimizer$step()
      
      train_loss <- loss$item()
      
      if (train_loss < best_train_loss) {
        best_train_loss  <- train_loss
        best_epoch       <- epoch
        patience_counter <- 0
      } else {
        patience_counter <- patience_counter + 1
      }
      
      if (patience_counter >= patience) {
        message(sprintf("Early stopping at epoch %d (best loss %.6f @ epoch %d)", 
                        epoch, best_train_loss, best_epoch))
        break
      }
    }
    
    # 8) Extract latent variables
    train_latent <- extract_latent_variables(model, x_data, device = device)
    val_latent   <- extract_latent_variables(model, x_val_data, device = device)
    
    # 9) Add survival outcomes
    train_latent$t_dmfs <- train_clinical$t_dmfs
    train_latent$e_dmfs <- train_clinical$e_dmfs
    val_latent$t_dmfs   <- test_clinical$t_dmfs
    val_latent$e_dmfs   <- test_clinical$e_dmfs
    
    # 10) Cox multivariable modeling with latent variables as predictors
    n_val <- ncol(train_latent) - 2
    multiv_formula <- as.formula(
      paste("Surv(t_dmfs, e_dmfs) ~", paste(paste0("V", seq_len(n_val)), collapse = " + "))
    )
    multiv_model <- survival::coxph(multiv_formula, data = train_latent, x = TRUE)
    
    # --- 11) Validation set evaluation (skip iteration on failure) ---
    
    result_valid <- try(
      calculate_time_auc_cindex(
        "Cox",
        fitted_model = multiv_model,
        df = val_latent
      ),
      silent = TRUE
    )
    
    failed <- FALSE
    if (inherits(result_valid, "try-error")) {
      failed <- TRUE
    } else {
      # Safely extract results
      ci   <- result_valid$c_index
      iauc <- result_valid$iAUC
      # If either metric is non-finite, mark as failed
      if (!is.finite(ci) || !is.finite(iauc)) failed <- TRUE
    }
    
    if (failed) {
      message(sprintf("Iteration %d: scoring failed or returned non-finite values. Skipping.", iteration))
      # Optional: record this iteration as NA 
      iteration_results[[iteration]] <- list(
        iteration = iteration,
        c_index   = NA_real_,
        iAUC      = NA_real_,
        score     = NA_real_,
        seed      = seed + iteration
      )
      next  # <--- Key: skip to the next iteration
    }
    
    # If this point is reached, valid c-index and iAUC were obtained
    ci   <- result_valid$c_index
    iauc <- result_valid$iAUC
    current_score <- ci + iauc
    
    message(sprintf("Iteration %d: C-index = %.4f, iAUC = %.4f, Score = %.4f",
                    iteration, ci, iauc, current_score))
    
    # Store results for this iteration
    iteration_results[[iteration]] <- list(
      iteration = iteration,
      c_index   = ci,
      iAUC      = iauc,
      score     = current_score,
      seed      = seed + iteration
    )
    
    # Update best model only if score is finite and better
    if (is.finite(current_score) && current_score > best_score) {
      best_score        <- current_score
      best_vae_model    <- model
      best_multiv_model <- multiv_model
      best_latent_dim   <- latent_dim
      best_iteration    <- iteration
      
      message(sprintf("New best model! Iteration %d, Score = %.4f", iteration, current_score))
      
      if (!is.null(save_dir)) {
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
        # Save VAE model
        torch_save(best_vae_model$state_dict(), file.path(save_dir, "best_vae_model.pt"))
        # Save Cox model
        saveRDS(best_multiv_model, file.path(save_dir, "best_cox_model.rds"))
        # Save model info
        model_info <- list(
          iteration  = best_iteration,
          score      = best_score,
          c_index    = ci,
          iAUC       = iauc,
          latent_dim = best_latent_dim,
          input_dim  = input_dim,
          hidden_dim = hidden_dim,
          seed       = seed + iteration,
          epochs     = epochs,
          patience   = patience,
          lr         = lr
        )
        saveRDS(model_info, file.path(save_dir, "model_info.rds"))
      }
    }
    
  }
  
  # --- Return objects for summary ---
  return(list(
    best_iteration    = best_iteration,
    best_score        = best_score,
    best_latent_dim   = best_latent_dim,
    iteration_results = iteration_results,
    best_multiv_model = best_multiv_model
  ))
}

#expr_mat <- expr
#clinical_df <- clinical_cleaned_7390

# res <- run_vae_cox_iterations(
#   expr          = expr,                   
#   clinical_df   = clinical_cleaned_7390,  
#   n_iterations  = 50,
#   seed          = 2025091520,
#   latent_dim    = 11,
#   hidden_dim    = 128,
#   epochs        = 800,
#   patience      = 50,
#   lr            = 0.001,
#   mad_quantile  = 0.25,
#   uni_cox_p     = 0.01,
#   corr_cutoff   = 0.90,
#   save_dir      = "best_models_dim11"           # NULL not save
# )
#saveRDS(res, "run_vae_cox_iterations_seed_2025091520.rds")
