# BreastCancerRecurrencePrediction
To build a multivariable survival model that integrates clinical features andselected gene expression profiles to predict long-term recurrence risk in breast cancer patients. 
Three models will be developed and compared: A clinical-only model; Ageneexpression-only model; A combined multimodal model. Performance will beevaluated to explore the potential improvement when combining clinical andgenomic features.

## Data Source
Dataset: GSE7390 from NCBI GEO
Cohort: 198 lymph node-negative breast cancer patients
Clinical variables: age, tumor grade, ER (estrogen receptor) status, tumor sizeetc. Gene expression: microarray data with thousands of probes
Outcome: recurrence status and time to recurrence (or overall survival time)

## Data Preprocessing and Exploratory Data Analysis (EDA)
Clean clinical data: handle missing values and outliers
Normalize gene expression data: e.g., log2 transformation, or z-score
standardization. Perform unsupervised clustering (e.g., K-means or hierarchical clustering) on gene expression data to identify potential molecular subgroups. Use t-SNE or
UMAP for 2D visualization of clustering patterns to aid interpretation of
subgroup structure. Feature Selection (if the proportional hazards (PH) assumption holds)
Method 1: Univariate Cox analysis to identify survival-associated genes
Method 2: LASSO-Cox regression to select key predictive genes

## Model Development
Three Cox proportional hazards models will be constructed:
Clinical model: age + tumor grade + ER status + tumor size + etc. Genomic model: selected gene expression features (e.g., top 30–100 genes)
Combined model: clinical variables + selected gene features
The risk score for patient will be calculated as:

<img width="477" height="53" alt="image" src="https://github.com/user-attachments/assets/a570af02-8982-4159-a379-de1e7f9ac2c5" />

Patients will be stratified into high- and low-risk groups (e.g., using median split). 

## Model Evaluation
C-index: concordance between predicted risk and actual survival time ranking; 
Time-dependent ROC/AUC: model sensitivity and specificity at various timepoints; 
Kaplan-Meier survival curves: visualize survival differences between risk groups; 
Use log-rank test to assess significance; 
External validation may be performed using an independent dataset with similar structure; 
Cross-validation will be performed.

## Challenges
High-dimensional genomic data: Requires careful feature selection to avoid overfitting, noise, and computational burden
PH assumption violation: May limit the use of Cox models;
## Further Exploration (options)
If the Proportional Hazards Assumption is Violated consider alternative models such as Random Survival Forests (RSF)
