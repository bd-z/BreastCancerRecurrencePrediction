[Predicting_Recurrence_Risk_in_Breast_Cancer.pdf](https://github.com/user-attachments/files/22452604/Predicting_Recurrence_Risk_in_Breast_Cancer.pdf)


This project evaluated multiple approaches to predicting distant recurrence-free survival (DRFS) in breast
cancer using gene expression and clinical data. Starting with more than 20,000 genes, univariate Cox
regression identified several thousand candidates, which were further reduced to a few dozen genes using
LASSO-Cox. To benchmark genomic against traditional clinical information for prognostic analysis, we
compared model performance across three input sets: (i) clinical variables only; (ii) clinical variables plus
gene-expression data; and (iii) gene-expression data only, using both Cox proportional hazards (PH) models
and Random Survival Forests (RSF). We also tested dimensionality-reduction methods—principal component
analysis (PCA) and variational autoencoder (VAE)—applied to the gene-expression data; the resulting
principal components and VAE latent variables were subsequently modeled with Cox regression.
Extensive bootstrap cross-validation was used to evaluate models by mean time-dependent AUC (iAUC)
and concordance index (C-index) across iterations. Data preprocessing included Multivariate Imputation by
Chained Equations (MICE) for missing clinical values, z-score normalization of gene-expression data, and
strict separation of training and validation sets to prevent data leakage. Across 100 bootstrap resamples
with replacement, none of the models achieved strong performance on the test set. By comparison, the
PCA-based Cox model modestly outperformed the others (mean iAUC = 0.57; mean C-index = 0.56). The
single best PCA–Cox model—identified on one train/test split by the highest score-was retained as the final
model for future prediction.
