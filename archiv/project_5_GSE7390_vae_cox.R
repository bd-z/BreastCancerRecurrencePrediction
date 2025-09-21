library(GEOquery)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(randomForestSRC)
library(pec)
library(mice)
options(scipen = 999)

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



expr_mat <- expr
clinical_df <- clinical_cleaned_7390
data_prepare_result <- bootstrap_sampl_split_miss_vale_imputation(expr_mat, clinical_df)
train_expr <- data_prepare_result$train_expr
test_expr <- data_prepare_result$test_expr
train_clinical <- data_prepare_result$train_clinical
test_clinical <- data_prepare_result$test_clinical
events_train <- data_prepare_result$events_train

# MAD filter
filtered_data <- filter_genes_by_mad(train_expr, test_expr, quantile_cutoff = 0.25)
train_expr2 <- filtered_data$train_expr
test_expr2 <- filtered_data$test_expr

# 单变量 Cox 回归
sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = 0.01)

significant_gene <- sig_gene_df$gene
message(sprintf(" extracted significant_gene length = %d",length(significant_gene)))

train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                            gene_mat_valid = train_expr2,
                                            significant_gene = significant_gene)
message("standardised train expr：train_expr_scaled dim = ", 
        toString(dim(train_expr_scaled)))

train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = 0.90)

significant_gene2 <- rownames(train_expr_filtered)

test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                           gene_mat_valid = test_expr2,
                                           significant_gene = significant_gene2)
message("standardised test expr：test_expr_scaled dim = ", 
        toString(dim(test_expr_scaled)))

library(torch)

train_expr_df <- t(train_expr_filtered)
test_expr_df <- t(test_expr_scaled)
input_dim <- ncol(train_expr_df)
#latent_dim <-50
hidden_dim <- 256 
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")



vae_module <- nn_module(
  "VAE",
  initialize = function(input_dim = input_dim, latent_dim = latent_dim, hidden_dim = hidden_dim) {
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$fc_mean <- nn_linear(hidden_dim, latent_dim)
    self$fc_logvar <- nn_linear(hidden_dim, latent_dim)
    self$fc2 <- nn_linear(latent_dim, hidden_dim)
    self$output <- nn_linear(hidden_dim, input_dim)
  },
  
  #mu, logvar 都是向量， torch_randn_like(std) 的结果也是向量
  
  reparameterize = function(mu, logvar) {
    std <- torch_exp(0.5 * logvar)
    eps <- torch_randn_like(std)
    #生成一个与 tensor 形状相同的张量，但这个新张量中的每个元素是从标准正态分布N(0,1) 中随机采样的
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


# 准备数据
x_data <- torch_tensor(train_expr_df, dtype = torch_float()) # expression_data
x_val_data <- torch_tensor(test_expr_df, dtype = torch_float())

#####
# model <- vae_module(
#   input_dim = input_dim,
#   latent_dim = latent_dim,
#   hidden_dim = hidden_dim
# )

#optimizer <- optim_adam(model$parameters, lr = 0.001)

#model$parameters 不是一个属性而是一个函数 即类似python中 model.parameters()，而这个函数的作用是：
# 返回模型中所有可学习的参数（张量）组成的列表
#####
# 定义损失函数
loss_fn <- function(recon_x, x, mu, logvar) {
  recon_loss <- nnf_mse_loss(recon_x, x, reduction = "sum")
  kl_loss <- -0.5 * torch_sum(1 + logvar - mu^2 - torch_exp(logvar))
  recon_loss + kl_loss
}
#####
# # 训练
# epochs <- 1000
# for (epoch in 1:epochs) {
#   model$train()
#   optimizer$zero_grad()
#   
#   output <- model(x_data)
#   loss <- loss_fn(output$x_recon, x_data, output$mu, output$logvar)
#   loss$backward()
#   optimizer$step()
#   
#   cat(sprintf("Epoch %d, Loss: %.2f\n", epoch, loss$item()))
# }


# #################################
# 
# library(torch)
# 
# # 假设您有训练数据 x_data 和验证数据 x_val_data
# #model <- VAE$new(input_dim = ncol(x_data), latent_dim = 50)
# #optimizer <- optim_adam(model$parameters, lr = 0.001)
# #loss_fn <- vae_loss_function
# 
# epochs <- 1000
# best_val_loss <- Inf
# best_epoch <- 0
# patience <- 50  # 早停法的耐心值
# patience_counter <- 0
# 
# # 记录损失历史
# train_loss_history <- numeric(epochs)
# val_loss_history <- numeric(epochs)
# 
# if (!dir.exists("models")) dir.create("models")
# 
# for (epoch in 1:epochs) {
#   # --- 训练阶段 ---
#   model$train()
#   optimizer$zero_grad()
#   
#   output <- model(x_data)
#   loss <- loss_fn(output$x_recon, x_data, output$mu, output$logvar)
#   loss$backward()
#   optimizer$step()
#   
#   train_loss <- loss$item()
#   train_loss_history[epoch] <- train_loss
#   
#   # --- 验证阶段 ---
#   model$eval()
#   with_no_grad({
#     val_output <- model(x_val_data)  # 使用验证数据！
#     val_loss <- loss_fn(val_output$x_recon, x_val_data, 
#                         val_output$mu, val_output$logvar)$item()
#   })
#   val_loss_history[epoch] <- val_loss
#   
#   cat(sprintf("Epoch %d/%d, Train Loss: %.4f, Val Loss: %.4f\n", 
#               epoch, epochs, train_loss, val_loss))
#   
#   # --- 基于验证损失保存最佳模型 ---
#   if (val_loss < best_val_loss) {
#     best_val_loss <- val_loss
#     best_epoch <- epoch
#     torch_save(model$state_dict(), "models/vae_model_best.pt")
#     patience_counter <- 0  # 重置耐心计数器
#     cat(sprintf("  -> New best model saved (val loss: %.4f)\n", best_val_loss))
#   } else {
#     patience_counter <- patience_counter + 1
#   }
#   
#   # 定期保存检查点
#   if (epoch %% 100 == 0) {
#     checkpoint <- list(
#       epoch = epoch,
#       model_state_dict = model$state_dict(),
#       optimizer_state_dict = optimizer$state_dict(),
#       train_loss = train_loss,
#       val_loss = val_loss
#     )
#     torch_save(checkpoint, sprintf("models/checkpoint_epoch_%d.pt", epoch))
#   }
#   
#   # --- 早停法 ---
#   if (patience_counter >= patience) {
#     cat(sprintf("Early stopping at epoch %d\n", epoch))
#     break
#   }
# }
# 
# # 保存最终模型和损失历史
# torch_save(model$state_dict(), "models/vae_model_final.pt")
# save(train_loss_history, val_loss_history, file = "models/loss_history.RData")
# 
# cat(sprintf("\nTraining completed! Best model from epoch %d with val loss: %.4f\n", 
#             best_epoch, best_val_loss))


##########################



# 定义要搜索的潜在维度范围
latent_dims <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
latent_dims <- c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50)
results <- list()

# 创建目录保存模型
if (!dir.exists("models")) dir.create("models")
if (!dir.exists("latent_dim_results")) dir.create("latent_dim_results")

for (latent_dim in latent_dims) {
  cat(sprintf("\n=== Training VAE with latent_dim = %d ===\n", latent_dim))
  
  # 初始化模型
  #model <- VAE$new(input_dim = ncol(x_data), latent_dim = latent_dim)
  model <- vae_module(
    input_dim = input_dim,
    latent_dim = latent_dim,
    hidden_dim = hidden_dim
  )
  
  
  optimizer <- optim_adam(model$parameters, lr = 0.001)
  #loss_fn <- vae_loss_function
  loss_fn <- function(recon_x, x, mu, logvar) {
    recon_loss <- nnf_mse_loss(recon_x, x, reduction = "sum")
    kl_loss <- -0.5 * torch_sum(1 + logvar - mu^2 - torch_exp(logvar))
    recon_loss + kl_loss
  }
  
  epochs <- 800
  best_val_loss <- Inf
  best_epoch <- 0
  patience <- 100
  patience_counter <- 0
  
  train_loss_history <- numeric(epochs)
  val_loss_history <- numeric(epochs)
  
  for (epoch in 1:epochs) {
    # --- 训练阶段 ---
    model$train()
    optimizer$zero_grad()
    
    output <- model(x_data)
    loss <- loss_fn(output$x_recon, x_data, output$mu, output$logvar)
    loss$backward()
    optimizer$step()
    
    train_loss <- loss$item()
    train_loss_history[epoch] <- train_loss
    
    # --- 验证阶段 ---
    model$eval()
    with_no_grad({
      val_output <- model(x_val_data)
      val_loss <- loss_fn(val_output$x_recon, x_val_data, 
                          val_output$mu, val_output$logvar)$item()
    })
    val_loss_history[epoch] <- val_loss
    
    if (epoch %% 100 == 0) {
      cat(sprintf("Epoch %d/%d, Train Loss: %.4f, Val Loss: %.4f\n", 
                  epoch, epochs, train_loss, val_loss))
    }
    
    # --- 基于验证损失保存最佳模型 ---
    if (val_loss < best_val_loss) {
      best_val_loss <- val_loss
      best_epoch <- epoch
      patience_counter <- 0
    } else {
      patience_counter <- patience_counter + 1
    }
    
    # --- 早停法 ---
    if (patience_counter >= patience) {
      cat(sprintf("Early stopping at epoch %d\n", epoch))
      break
    }
  }
  
  # 创建当前维度的结果对象
  current_result <- list(
    latent_dim = latent_dim,
    best_val_loss = best_val_loss,
    best_epoch = best_epoch,
    final_epoch = epoch,
    train_loss_history = train_loss_history[1:epoch],
    val_loss_history = val_loss_history[1:epoch]
  )
  
  # 保存到results列表
  results[[as.character(latent_dim)]] <- current_result
  
  # 保存模型和结果
  model_dir <- sprintf("models/latent_dim_%d", latent_dim)
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
  
  torch_save(model$state_dict(), file.path(model_dir, "vae_model_final.pt"))
  
  # 单独保存当前维度的结果
  save(current_result, file = file.path(model_dir, "training_results.RData"))
  
  cat(sprintf("Completed latent_dim %d: Best val loss = %.4f at epoch %d\n",
              latent_dim, best_val_loss, best_epoch))
}

# 保存所有结果
save(results, file = "latent_dim_results/all_results.RData")

# 分析结果并选择最佳latent_dim
best_latent_dim <- NULL
best_overall_loss <- Inf

cat("\n=== Results Summary ===\n")
for (latent_dim in latent_dims) {
  result <- results[[as.character(latent_dim)]]
  cat(sprintf("latent_dim %d: Best val loss = %.4f (epoch %d)\n",
              latent_dim, result$best_val_loss, result$best_epoch))
  
  if (result$best_val_loss < best_overall_loss) {
    best_overall_loss <- result$best_val_loss
    best_latent_dim <- latent_dim
  }
}

cat(sprintf("\nBest latent_dim: %d with validation loss: %.4f\n",
            best_latent_dim, best_overall_loss))

# 绘制结果比较图
if (length(results) > 0) {
  latent_dims <- sapply(results, function(x) x$latent_dim)
  val_losses <- sapply(results, function(x) x$best_val_loss)
  
  plot(latent_dims, val_losses, type = "b", pch = 19, col = "blue",
       xlab = "Latent Dimension", ylab = "Best Validation Loss",
       main = "VAE Performance vs Latent Dimension")
  points(best_latent_dim, best_overall_loss, pch = 19, col = "red", cex = 1.5)
  legend("topright", legend = "Best", pch = 19, col = "red")
}


######

library(torch)

# 加载模型
model_path <- "G:/clinical_data_science/final_project_github/BreastCancerRecurrencePrediction/models/latent_dim_30/vae_model_final.pt"

# 重新实例化模型（参数必须与训练时相同）
model1 <- vae_module(
  input_dim = input_dim,
  latent_dim = 30,
  hidden_dim = hidden_dim
)

#####
#model <- VAE$new(input_dim = ncol(x_data), latent_dim = 20)
#model1$load_state_dict(torch_load(model_path))

# # 提取训练数据潜在变量
# data_tensor <- torch_tensor(as.matrix(x_data), dtype = torch_float32())$to(device = device)
# with_no_grad({
#   output <- model1(data_tensor)
#   train_latent_mu <- as.array(output$mu$cpu())
#   train_latent_logvar <- as.array(output$logvar$cpu())
# })


# 保存结果
#save(train_latent_mu, train_latent_logvar, file = "train_latent_variables.RData")
#write.csv(train_latent_mu, "train_latent_mu.csv", row.names = FALSE)

# 查看前几个样本
# cat("前5个样本的前5个潜在特征:\n")
# print(train_latent_mu[1:5, 1:5])

#####








# 封装成函数
extract_latent_variables <- function(model=model1, data=x_data, device) {
  # """
  # 从VAE模型中提取潜在变量
  # 
  # 参数:
  # model: 训练好的VAE模型
  # data: 输入数据矩阵/数据框
  # use_mu: 是否使用均值mu（TRUE）还是重参数化z（FALSE）
  # 
  # 返回:
  # 潜在变量矩阵
  # """
  
  # 转换为tensor
  data_tensor <- torch_tensor(as.matrix(data), dtype = torch_float32())$to(device = device)
  
  model$eval()
  model$to(device = device)
  
  with_no_grad({
    output <- model(data_tensor)
    train_latent_mu <- as.array(output$mu$cpu())
  })
  
  return(as.data.frame(train_latent_mu))
}

# 使用函数提取训练数据潜在变量
train_latent <- extract_latent_variables(model1, x_data, device = device)
cat("训练数据潜在变量提取完成！\n")
cat(sprintf("维度: %d x %d\n", nrow(train_latent), ncol(train_latent)))

# 同样可以提取验证数据潜在变量
val_latent <- extract_latent_variables(model1, x_val_data, device = device)

cat("验证数据潜在变量提取完成！\n")
cat(sprintf("维度: %d x %d\n", nrow(val_latent), ncol(val_latent)))

n_val = ncol(train_latent)

#train_clinical$geo_accession

# 添加生存数据
train_latent$t_dmfs <- train_clinical$t_dmfs
train_latent$e_dmfs <- train_clinical$e_dmfs
val_latent$t_dmfs <- test_clinical$t_dmfs
val_latent$e_dmfs <- test_clinical$e_dmfs

multiv_formula <- as.formula(paste("Surv(t_dmfs, e_dmfs) ~", 
                                   paste(paste0("V", 1:n_val), collapse = " + ")))
multiv_model <- coxph(multiv_formula, data = train_latent, x=TRUE)

# 查看模型摘要
summary(multiv_model)
# 
# ph_test <- cox.zph(multiv_model)
# print(ph_test)


result_valid <- calculate_time_auc_cindex(
  "Cox", 
  fitted_model = multiv_model, 
  df = val_latent
)

message(sprintf("Bootstrap %d: Cox 模型评估完成，iAUC = %.4f, c_index = %.4f", 
                b, result_valid$iAUC, result_valid$c_index))



###############################################

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



# 设置迭代参数
n_iterations <- 10
best_score <- -Inf
best_vae_model <- NULL
best_multiv_model <- NULL
best_latent_dim <- NULL
best_iteration <- 0

# 创建结果存储
iteration_results <- list()
expr_mat <- expr
clinical_df <- clinical_cleaned_7390


for (iteration in 1:n_iterations) {
  # 设置随机种子
  set.seed(1234 + iteration)
  message(sprintf("\n=== 开始迭代 %d/%d ===\n", iteration, n_iterations))
  
  data_prepare_result <- bootstrap_sampl_split_miss_vale_imputation(expr_mat, clinical_df)
  train_expr <- data_prepare_result$train_expr
  test_expr <- data_prepare_result$test_expr
  train_clinical <- data_prepare_result$train_clinical
  test_clinical <- data_prepare_result$test_clinical
  events_train <- data_prepare_result$events_train
  
  # MAD filter
  filtered_data <- filter_genes_by_mad(train_expr, test_expr, quantile_cutoff = 0.25)
  train_expr2 <- filtered_data$train_expr
  test_expr2 <- filtered_data$test_expr
  
  # 单变量 Cox 回归
  sig_gene_df <- batch_univariate_cox_regression(train_expr2, train_clinical, p_value = 0.01)
  significant_gene <- sig_gene_df$gene
  message(sprintf("提取显著基因数量 = %d", length(significant_gene)))
  
  train_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                              gene_mat_valid = train_expr2,
                                              significant_gene = significant_gene)
  
  train_expr_filtered <- remove_high_corr_genes(train_expr_scaled, cutoff = 0.90)
  significant_gene2 <- rownames(train_expr_filtered)
  
  test_expr_scaled <- standardize_with_train(gene_mat_train = train_expr2,
                                             gene_mat_valid = test_expr2,
                                             significant_gene = significant_gene2)
  
  # 准备数据
  train_expr_df <- t(train_expr_filtered)
  test_expr_df <- t(test_expr_scaled)
  input_dim <- ncol(train_expr_df)
  latent_dim <- 16
  hidden_dim <- 256 
  device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
  
  # VAE 模型定义
    
  # 准备数据
  x_data <- torch_tensor(train_expr_df, dtype = torch_float())
  x_val_data <- torch_tensor(test_expr_df, dtype = torch_float())
  
  model <- vae_module(
    input_dim = input_dim,
    latent_dim = latent_dim,
    hidden_dim = hidden_dim
  )
  
  optimizer <- optim_adam(model$parameters, lr = 0.001)
    
  epochs <- 800
  best_train_loss <- Inf
  best_epoch <- 0
  patience <- 50
  patience_counter <- 0
  
  # 训练 VAE
  for (epoch in 1:epochs) {
    model$train()
    optimizer$zero_grad()
    
    output <- model(x_data)
    loss <- loss_fn(output$x_recon, x_data, output$mu, output$logvar)
    loss$backward()
    optimizer$step()
    
    train_loss <- loss$item()
    
    if (train_loss < best_train_loss) {
      best_train_loss <- train_loss
      best_epoch <- epoch
      patience_counter <- 0
    } else {
      patience_counter <- patience_counter + 1
    }
    
    if (patience_counter >= patience) {
      message(sprintf("早停于第 %d 轮", epoch))
      break
    }
  }
  
  # 提取潜在变量
  extract_latent_variables <- function(model, data, device = "cpu") {
    model$eval()
    with_no_grad({
      output <- model(data)
      latent_vars <- as.matrix(output$mu$cpu())
      as.data.frame(latent_vars)
    })
  }
  
  train_latent <- extract_latent_variables(model, x_data, device = device)
  val_latent <- extract_latent_variables(model, x_val_data, device = device)
  
  # 添加生存数据
  train_latent$t_dmfs <- train_clinical$t_dmfs
  train_latent$e_dmfs <- train_clinical$e_dmfs
  val_latent$t_dmfs <- test_clinical$t_dmfs
  val_latent$e_dmfs <- test_clinical$e_dmfs
  
  # 训练 Cox 模型
  n_val <- ncol(train_latent) - 2  # 减去两个生存变量
  multiv_formula <- as.formula(paste("Surv(t_dmfs, e_dmfs) ~", 
                                     paste(paste0("V", 1:n_val), collapse = " + ")))
  multiv_model <- coxph(multiv_formula, data = train_latent, x = TRUE)
  
  # 评估模型
  result_valid <- calculate_time_auc_cindex(
    "Cox", 
    fitted_model = multiv_model, 
    df = val_latent
  )
  
  # 计算评分
  current_score <- result_valid$c_index + result_valid$iAUC
  message(sprintf("迭代 %d: C-index = %.4f, iAUC = %.4f, 总分 = %.4f", 
                  iteration, result_valid$c_index, result_valid$iAUC, current_score))
  
  # 保存迭代结果
  iteration_results[[iteration]] <- list(
    iteration = iteration,
    c_index = result_valid$c_index,
    iAUC = result_valid$iAUC,
    score = current_score,
    seed = 1234 + iteration
  )
  
  # 更新最佳模型
  if (current_score > best_score) {
    best_score <- current_score
    best_vae_model <- model
    best_multiv_model <- multiv_model
    best_latent_dim <- latent_dim
    best_iteration <- iteration
    
    message(sprintf("新的最佳模型! 迭代 %d, 总分 = %.4f", iteration, current_score))
    
    # 保存最佳模型
    model_dir <- "best_models"
    if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
    
    # 保存 VAE 模型
    torch_save(best_vae_model$state_dict(), file.path(model_dir, "best_vae_model.pt"))
    
    # 保存 Cox 模型
    saveRDS(best_multiv_model, file.path(model_dir, "best_cox_model.rds"))
    
    # 保存模型信息
    model_info <- list(
      iteration = best_iteration,
      score = best_score,
      c_index = result_valid$c_index,
      iAUC = result_valid$iAUC,
      latent_dim = best_latent_dim,
      input_dim = input_dim,
      hidden_dim = hidden_dim,
      seed = 1234 + iteration
    )
    saveRDS(model_info, file.path(model_dir, "model_info.rds"))
  }
}

# 输出最终结果
message(sprintf("\n=== 完成所有迭代 ==="))
message(sprintf("最佳模型来自迭代 %d", best_iteration))
message(sprintf("最佳分数: %.4f (C-index: %.4f + iAUC: %.4f)", 
                best_score, 
                iteration_results[[best_iteration]]$c_index,
                iteration_results[[best_iteration]]$iAUC))

# 保存所有迭代结果
saveRDS(iteration_results, "iteration_results.rds")



# === 迭代結果彙總 ===
iter_df <- do.call(rbind, lapply(iteration_results, as.data.frame))
# 防止有 NA 或中途失敗的情況
mean_c_index <- mean(iter_df$c_index, na.rm = TRUE)
mean_iAUC   <- mean(iter_df$iAUC,   na.rm = TRUE)
sd_c_index  <- sd(iter_df$c_index,  na.rm = TRUE)
sd_iAUC     <- sd(iter_df$iAUC,     na.rm = TRUE)

message(sprintf("平均 C-index: %.4f (SD = %.4f)", mean_c_index, sd_c_index))
message(sprintf("平均 iAUC:   %.4f (SD = %.4f)", mean_iAUC, sd_iAUC))

# 也可以順手存成表格
iter_df$iteration <- sapply(iteration_results, function(x) x$iteration)
iter_df$seed      <- sapply(iteration_results, function(x) x$seed)
write.csv(iter_df, "iteration_results.csv", row.names = FALSE)

# 保存平均指標
avg_metrics <- list(mean_c_index = mean_c_index,
                    sd_c_index   = sd_c_index,
                    mean_iAUC    = mean_iAUC,
                    sd_iAUC      = sd_iAUC,
                    n_iterations = nrow(iter_df))
saveRDS(avg_metrics, "avg_metrics.rds")




























































##### 
# 
# #在 R 中可视化 VAE 隐藏变量与 grade 之间关系的代码
# library(ggplot2)
# 
# # 假设你的数据框已经存在并命名为 df_latent
# # 如果 grade 为字符型，建议先转换为因子型
# df_latent$grade <- as.factor(df_latent$grade)
# 
# # 绘制散点图
# ggplot(df_latent, aes(x = V1, y = V2, color = treatment)) +
#   geom_point(size = 3, alpha = 0.7) +
#   labs(
#     title = "Visualization of VAE Latent Variables by Grade",
#     x = "Latent Variable V1",
#     y = "Latent Variable V2",
#     color = "Grade"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 14),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12)
#   )
# 
# #对数据框中的两个隐藏变量 V1,V2 进行聚类分析
# #####
# # 加载需要的包
# library(cluster)
# 
# # 提取V1和V2列用于聚类
# latent_vars <- df_latent[, c("V1", "V2")]
# 
# # 使用k-means聚类方法 (这里假设聚成3类，可根据实际情况修改)
# set.seed(123)  # 设定随机种子，保证可复现性
# kmeans_result <- kmeans(latent_vars, centers = 3)
# 
# # 将聚类结果添加到原数据框中
# df_latent$cluster <- as.factor(kmeans_result$cluster)
# 
# # 可视化聚类结果
# ggplot(df_latent, aes(x = V1, y = V2, color = cluster)) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(
#     title = "K-means Clustering of VAE Latent Variables",
#     x = "Latent Variable V1",
#     y = "Latent Variable V2",
#     color = "Cluster"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 14),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12)
#   )
# 
# 
# # 合并数据框（按行合并，假设两个数据框的行数相同且顺序一致）
# final_df <- cbind(merged_df_filtered[, 1:13], df_latent[, c(1, 2, 6)])
# 
# ggplot(final_df, aes(x = V1, y = V2, color = er)) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(
#     title = "K-means Clustering of VAE Latent Variables",
#     x = "Latent Variable V1",
#     y = "Latent Variable V2",
#     color = ""
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 14),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12)
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(dplyr)
# 
# # 将cluster转换为数值型（便于相关性计算）
# final_df$cluster_numeric <- as.numeric(as.character(final_df$cluster))
# 
# # 仅选择数值型变量进行相关性计算
# numeric_cols <- final_df %>% 
#   select(where(is.numeric))
# 
# # 计算与cluster_numeric的相关性
# correlations <- cor(numeric_cols, use = "complete.obs")
# 
# # 提取cluster_numeric这一列的相关性并排序
# cluster_corr <- sort(abs(correlations[,"cluster_numeric"]), decreasing = TRUE)
# 
# # 显示相关性最高的列
# cluster_corr
# 
# 
# 
# 


#####
#cox
######################
# 构建多变量Cox模型
multiv_formula <- as.formula(paste("Surv(t_dmfs, e_dmfs) ~", 
                                   paste(paste0("V", 1:50), collapse = " + ")))
multiv_model <- coxph(multiv_formula, data = df_latent)

# 查看模型摘要
summary(multiv_model)













# 完整的训练数据转换流程
library(torch)

# 设置设备
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")

# 加载数据
#load("your_training_data.RData")  # 加载您的训练数据
# x_data 应该已经预处理过

# 加载模型
model_path <- "G:/clinical_data_science/final_project_github/BreastCancerRecurrencePrediction/models/latent_dim_20/vae_model_final.pt"
model1 <- vae_module(
  input_dim = input_dim,
  latent_dim = 20,
  hidden_dim = hidden_dim
)

#model <- VAE$new(input_dim = ncol(x_data), latent_dim = 20)
model1$load_state_dict(torch_load(model_path))
model1$to(device = device)
model1$eval()

with_no_grad({
  # 使用编码器获取潜在表示
  encoded_output <- model1$encode(x_data_tensor)
  
  # 通常使用均值mu作为潜在变量（更稳定）
  # 如果您在训练Cox模型时使用的是mu，这里也用mu
  latent_mu <- encoded_output[[1]]  # 均值
  latent_logvar <- encoded_output[[2]]  # 对数方差
  
  # 或者如果您需要重参数化后的z（通常用于训练，评估时用mu）
  # std <- torch_exp(0.5 * latent_logvar)
  # eps <- torch_randn_like(std)
  # latent_z <- latent_mu + eps * std
})


# 提取训练数据潜在变量
cat("正在提取训练数据潜在变量...\n")
train_latent <- extract_latent_variables(model1, x_data, use_mu = TRUE)

# 提取验证数据潜在变量（如果有）
if (exists("x_val_data")) {
  cat("正在提取验证数据潜在变量...\n")
  val_latent <- extract_latent_variables(model, x_val_data, use_mu = TRUE)
}

# 为Cox模型准备数据
# 假设您有生存数据：time_train, event_train, time_val, event_val

# 训练集
train_df <- data.frame(
  train_latent,
  time = time_train,
  event = event_train
)

# 验证集
val_df <- data.frame(
  val_latent,
  time = time_val, 
  event = event_val
)

# 设置列名
colnames(train_df)[1:ncol(train_latent)] <- paste0("z", 1:ncol(train_latent))
colnames(val_df)[1:ncol(val_latent)] <- paste0("z", 1:ncol(val_latent))

# 保存所有数据
save(train_df, val_df, file = "cox_ready_data.RData")

cat("数据准备完成！\n")
cat(sprintf("训练集: %d 样本, %d 特征\n", nrow(train_df), ncol(train_df) - 2))
cat(sprintf("验证集: %d 样本, %d 特征\n", nrow(val_df), ncol(val_df) - 2))
















######




















#   
#   # 保存该latent_dim的最佳结果
#   results[[as.character(latent_dim)]] <- list(
#     latent_dim = latent_dim,
#     best_val_loss = best_val_loss,
#     best_epoch = best_epoch,
#     final_epoch = epoch,
#     train_loss_history = train_loss_history[1:epoch],
#     val_loss_history = val_loss_history[1:epoch]
#   )
#   
#   # 保存模型
#   model_dir <- sprintf("models/latent_dim_%d", latent_dim)
#   if (!dir.exists(model_dir)) dir.create(model_dir)
#   
#   torch_save(model$state_dict(), file.path(model_dir, "vae_model_final.pt"))
#   save(results[[as.character(latent_dim)]], 
#        file = file.path(model_dir, "training_results.RData"))
#   
#   cat(sprintf("Completed latent_dim %d: Best val loss = %.4f at epoch %d\n",
#               latent_dim, best_val_loss, best_epoch))
# }

# 分析结果并选择最佳latent_dim
# best_latent_dim <- NULL
# best_overall_loss <- Inf
# 
# cat("\n=== Results Summary ===\n")
# for (latent_dim in latent_dims) {
#   result <- results[[as.character(latent_dim)]]
#   cat(sprintf("latent_dim %d: Best val loss = %.4f (epoch %d)\n",
#               latent_dim, result$best_val_loss, result$best_epoch))
#   
#   if (result$best_val_loss < best_overall_loss) {
#     best_overall_loss <- result$best_val_loss
#     best_latent_dim <- latent_dim
#   }
# }
# 
# cat(sprintf("\nBest latent_dim: %d with validation loss: %.4f\n",
#             best_latent_dim, best_overall_loss))
# 
# # 保存所有结果
# save(results, file = "latent_dim_results/all_results.RData")
# 
# # 绘制结果比较图
# plot_latent_dim_results <- function(results) {
#   latent_dims <- sapply(results, function(x) x$latent_dim)
#   val_losses <- sapply(results, function(x) x$best_val_loss)
#   
#   plot(latent_dims, val_losses, type = "b", pch = 19, col = "blue",
#        xlab = "Latent Dimension", ylab = "Best Validation Loss",
#        main = "VAE Performance vs Latent Dimension")
#   points(best_latent_dim, best_overall_loss, pch = 19, col = "red", cex = 1.5)
#   legend("topright", legend = "Best", pch = 19, col = "red")
# }
# 
# # 绘制比较图
# plot_latent_dim_results(results)
# 



#latent <- model(x_data)$mu$to(device = "cpu")$detach()$numpy()

# 1. 获取模型输出
output <- model(x_data)

# 2. 取出 mu（编码器输出的潜在变量）
mu <- output$mu

# 3. 移到 CPU，取消梯度，转为 R matrix
mu_cpu <- mu$to(device = "cpu")
mu_detached <- mu_cpu$detach()
latent <- as.matrix(mu_detached)
df_latent0 <- as.data.frame(latent)
df_latent <- as.data.frame(latent)
sum(rownames(train_expr_df) != train_clinical$geo_accession)



train_clinical$geo_accession

# 添加生存数据
df_latent$t_dmfs <- train_clinical$t_dmfs
df_latent$e_dmfs <- train_clinical$e_dmfs


##### 
# 
# #在 R 中可视化 VAE 隐藏变量与 grade 之间关系的代码
# library(ggplot2)
# 
# # 假设你的数据框已经存在并命名为 df_latent
# # 如果 grade 为字符型，建议先转换为因子型
# df_latent$grade <- as.factor(df_latent$grade)
# 
# # 绘制散点图
# ggplot(df_latent, aes(x = V1, y = V2, color = treatment)) +
#   geom_point(size = 3, alpha = 0.7) +
#   labs(
#     title = "Visualization of VAE Latent Variables by Grade",
#     x = "Latent Variable V1",
#     y = "Latent Variable V2",
#     color = "Grade"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 14),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12)
#   )
# 
# #对数据框中的两个隐藏变量 V1,V2 进行聚类分析
# #####
# # 加载需要的包
# library(cluster)
# 
# # 提取V1和V2列用于聚类
# latent_vars <- df_latent[, c("V1", "V2")]
# 
# # 使用k-means聚类方法 (这里假设聚成3类，可根据实际情况修改)
# set.seed(123)  # 设定随机种子，保证可复现性
# kmeans_result <- kmeans(latent_vars, centers = 3)
# 
# # 将聚类结果添加到原数据框中
# df_latent$cluster <- as.factor(kmeans_result$cluster)
# 
# # 可视化聚类结果
# ggplot(df_latent, aes(x = V1, y = V2, color = cluster)) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(
#     title = "K-means Clustering of VAE Latent Variables",
#     x = "Latent Variable V1",
#     y = "Latent Variable V2",
#     color = "Cluster"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 14),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12)
#   )
# 
# 
# # 合并数据框（按行合并，假设两个数据框的行数相同且顺序一致）
# final_df <- cbind(merged_df_filtered[, 1:13], df_latent[, c(1, 2, 6)])
# 
# ggplot(final_df, aes(x = V1, y = V2, color = er)) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(
#     title = "K-means Clustering of VAE Latent Variables",
#     x = "Latent Variable V1",
#     y = "Latent Variable V2",
#     color = ""
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 14),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12)
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(dplyr)
# 
# # 将cluster转换为数值型（便于相关性计算）
# final_df$cluster_numeric <- as.numeric(as.character(final_df$cluster))
# 
# # 仅选择数值型变量进行相关性计算
# numeric_cols <- final_df %>% 
#   select(where(is.numeric))
# 
# # 计算与cluster_numeric的相关性
# correlations <- cor(numeric_cols, use = "complete.obs")
# 
# # 提取cluster_numeric这一列的相关性并排序
# cluster_corr <- sort(abs(correlations[,"cluster_numeric"]), decreasing = TRUE)
# 
# # 显示相关性最高的列
# cluster_corr
# 
# 
# 
# 


#####
#cox
######################
# 构建多变量Cox模型
multiv_formula <- as.formula(paste("Surv(t_dmfs, e_dmfs) ~", 
                                   paste(paste0("V", 1:50), collapse = " + ")))
multiv_model <- coxph(multiv_formula, data = df_latent)

# 查看模型摘要
summary(multiv_model)









# 
# df_latent_t = t(df_latent0)
# a = df_latent_t[colnames(df_latent0),]
# 
# gene_freq_df1 <- repeat_cv_lasso_cox(train_expr = df_latent_t,
#                                     train_clinical,
#                                     significant_gene_vec = colnames(df_latent0),
#                                     repeats = 5,
#                                     nfolds = 10,
#                                     alpha = 1)
# 



predictors <- c(paste0("V", 1:5))


formula <- as.formula(paste(
  "Surv(time, status) ~",
  paste(predictors, collapse = " + ")
))

# ?M?? Cox ģ??
cox_model <- coxph(formula, data = df_latent)

# ?鿴ģ??ժҪ
summary(cox_model)


# 计算C-index

concordance(cox_model)



df_latent$time_day = df_latent$time
df_latent$status = as.integer(df_latent$status)

library(timeROC)
roc <- timeROC(T = df_latent$time_day, delta = df_latent$status,
               marker = predict(cox_model),
               cause = 1, times = c(365, 730), iid = TRUE)
#plot(roc)

plot(roc, time = 365, col = "blue", title = "1-year ROC")
plot(roc, time = 730, col = "red", title = "2-year ROC")















############################################



ggplot(df_latent, aes(x = dim1, y = dim2, color = as.factor(event))) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "VAE Latent Space by RFS Event", color = "Event (RFS)")




library(survival)
library(survminer)

#install.packages("survminer")

# 将 latent dim1 分为高低表达组
df_latent$group <- ifelse(df_latent$dim2 > median(df_latent$dim2), "High", "Low")

# 创建生存对象
surv_obj <- Surv(time = df_latent$time, event = df_latent$event)

# KM 曲线
fit <- survfit(surv_obj ~ group, data = df_latent)

ggsurvplot(fit, data = df_latent, pval = TRUE,
           title = "Kaplan-Meier by VAE latent dim1",
           legend.title = "Latent group")
























# 加上标签（比如治疗组或分级）
df_latent$group <- merged_df_filtered$treatment  # 替换为你关心的列，例如 grade, node 等

# 可视化
library(ggplot2)

ggplot(df_latent, aes(x = dim1, y = dim2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "VAE Latent Space", color = "Group")

########################
#install.packages("Rtsne")  # 若尚未安装
library(Rtsne)
set.seed(123)  # 保证可重复

pca_out <- prcomp(expression_data, center = TRUE, scale. = TRUE)
pca50 <- pca_out$x[, 1:50]  

tsne_out <- Rtsne(expression_data, dims = 2, perplexity = 30, verbose = TRUE)

df_tsne <- as.data.frame(tsne_out$Y)
colnames(df_tsne) <- c("TSNE1", "TSNE2")
df_tsne$event <- as.factor(merged_df_filtered$event.rfs)

library(ggplot2)
ggplot(df_tsne, aes(x = TSNE1, y = TSNE2, color = event)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "t-SNE on Expression Data", color = "Event RFS")









install.packages("uwot")
library(uwot)

set.seed(123)
umap_out <- umap(expression_data, n_neighbors = 15)

df_umap <- as.data.frame(umap_out)
colnames(df_umap) <- c("UMAP1", "UMAP2")



# 加入感兴趣的分组变量
df_umap$event <- as.factor(merged_df_filtered$event.rfs)
df_umap$grade <- as.factor(merged_df_filtered$grade)
df_umap$treatment <- as.factor(merged_df_filtered$treatment)

ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = grade)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP on Expression Data (colored by RFS event)", color = "treatment")




##### 
################################## CV


library(survival)
library(timeROC)  # or Hmisc for C-index

latent_dims <- c(2, 4, 5, 8, 16)
k_folds <- 5
input_dim <- ncol(expression_data)

results <- data.frame(latent_dim = latent_dims, mean_cindex = NA)

for (i in seq_along(latent_dims)) {
  dim <- latent_dims[i]
  cindexes <- c()
  
  folds <- sample(rep(1:k_folds, length.out = nrow(normalized_df)))
  
  for (k in 1:k_folds) {
    # 分出训练和验证集
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)
    
    x_train <- torch_tensor(as.matrix(normalized_df[train_idx, 14:147]), dtype = torch_float())
    x_test <- torch_tensor(as.matrix(normalized_df[test_idx, 14:147]), dtype = torch_float())
    
    # 训练 VAE 模型（关键：动态指定 latent_dim）
    model <- vae_module(input_dim, latent_dim = dim)
    optimizer <- optim_adam(model$parameters, lr = 0.001)
    
    # early stop
    best_loss <- Inf
    patience <- 5
    no_improve <- 0
    
    for (epoch in 1:200) {
      model$train()
      optimizer$zero_grad()
      output <- model(x_train)
      loss <- loss_fn(output$x_recon, x_train, output$mu, output$logvar)
      loss$backward()
      optimizer$step()
      cat(sprintf("Epoch %d, Loss: %.4f\n", epoch, loss$item()))
      
      if (loss$item() < best_loss - 1e-3) {
        best_loss <- loss$item()
        no_improve <- 0
      } else {
        no_improve <- no_improve + 1
        if (no_improve >= patience) {
          cat("Early stopping at epoch", epoch, "\n")
          break
        }
      }
    }
    
    # 提取训练 + 测试的 latent 表示（mu）
    mu_train <- as.matrix(model(x_train)$mu$detach()$cpu())
    mu_test <- as.matrix(model(x_test)$mu$detach()$cpu())
    
    # 生存数据
    surv_train <- normalized_df[train_idx, c("time","status")]
    surv_test <- normalized_df[test_idx, c("time","status")]
    
    # 用训练集拟合 Cox
    cox_model <- coxph(Surv(time, status) ~ ., data = data.frame(mu_train, surv_train))
    
    # 在验证集上预测
    risk_score <- predict(cox_model, newdata = data.frame(mu_test))
    
    # 评估 C-index
    cidx <- Hmisc::rcorr.cens(risk_score, Surv(surv_test$time, surv_test$status))["C Index"]
    cindexes <- c(cindexes, cidx)
  }
  # 平均 C-index
  results$mean_cindex[i] <- mean(cindexes)
}

best_dim <- results$latent_dim[which.max(results$mean_cindex)]
print(results)
cat("最佳潜变量维度是:", best_dim, "\n")

