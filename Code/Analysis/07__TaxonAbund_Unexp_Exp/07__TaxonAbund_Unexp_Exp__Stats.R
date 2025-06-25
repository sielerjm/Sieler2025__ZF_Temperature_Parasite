# 07__TaxonAbund_Unexp_Exp__Stats --------------------------------------------------

# Add start message
start_time <- Sys.time()
cat("Starting 07__TaxonAbund_Unexp_Exp__Stats.R at", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

tmp.psOBJ <- ps.list[["All"]]
tmp.resSubSection <- "All"

# Add Cluster Columns
tmp.psOBJ <- 
  tmp.psOBJ %>%
  
  # Group samples by Alpha Score
  ps_mutate(Cluster = if_else(
    Treatment == "Exposed" & Total.Worm.Count > 0,
    case_when(
      Simpson__Genus_norm <= 0.5 ~ "Low",
      Simpson__Genus_norm > 0.5 ~ "High",
      TRUE ~ "Other"
    ),
    "Other"
  ), .after = Treatment) %>%
  ps_mutate(Cluster = fct_relevel(factor(Cluster, levels = c("Other", "Low", "High")))) 




## TEMP_DPE_TREAT_PATH_WORM_CLUSTER ----------------------------------------


### MaAsLin2 ----------------------------------------------------------------

#### Run Maaslin2

# Methods derived from:
#   Elementary methods provide more replicable results in microbial differential abundance analysis
#   Juho Pelto, Kari Auranen, Janne Kujala, Leo Lahti
#   https://arxiv.org/abs/2404.02691

tmp.output.path.diffAbund <- paste0("All__TEMP_DPE_TREAT_PATH_WORM_CLUSTER")
tmp.file.path.diffAbund <- file.path(path.results, "Tables/MaAsLin2", tmp.output.path.diffAbund, "significant_results.tsv")

run_maaslin2(tmp.PS = tmp.psOBJ,
             tmp.fixed = c("Temperature", "DPE", "Treatment", "Pathology.Results", "Total.Worm.Count", "Cluster"),
             tmp.reference = c("Temperature,28", "Treatment,Control", "Pathology.Results,negative", "Cluster,Other"),
             tmp.output = tmp.output.path.diffAbund
)

#### Table

diffAbnd.stats[["All"]][["All__TEMP_DPE_TREAT_PATH_WORM_CLUSTER"]][["Maaslin2"]][["output"]] <- 
  readr::read_delim(tmp.file.path.diffAbund, 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>%
  dplyr::left_join(tmp.psOBJ %>% 
                     microViz::ps_melt() %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, " ", ".")) %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, "-", ".")) %>%
                     dplyr::select(Kingdom:Genus) %>%
                     dplyr::distinct(Genus, .keep_all = T), by = c("feature" = "Genus")) %>%
  dplyr::rename(Taxon = "feature")  %>%
  dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, "\\.", " ")) %>%
  dplyr::arrange(qval) %>%
  dplyr::ungroup()



## Exp__WORM ----------------------------------------


### MaAsLin2 ----------------------------------------------------------------

#### Run Maaslin2

# Methods derived from:
#   Elementary methods provide more replicable results in microbial differential abundance analysis
#   Juho Pelto, Kari Auranen, Janne Kujala, Leo Lahti
#   https://arxiv.org/abs/2404.02691

tmp.output.path.diffAbund <- paste0("Exp__WORM")
tmp.file.path.diffAbund <- file.path(path.results, "Tables/MaAsLin2", tmp.output.path.diffAbund, "significant_results.tsv")

run_maaslin2(tmp.PS = tmp.psOBJ %>% 
               microViz::ps_filter(Treatment == "Exposed"),
             tmp.fixed = c("Total.Worm.Count"),
             tmp.reference = c(),
             tmp.output = tmp.output.path.diffAbund
)

#### Table

diffAbnd.stats[["All"]][["Exp__WORM"]][["Maaslin2"]][["output"]] <- 
  readr::read_delim(tmp.file.path.diffAbund, 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>%
  dplyr::left_join(tmp.psOBJ %>% 
                     microViz::ps_melt() %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, " ", ".")) %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, "-", ".")) %>%
                     dplyr::select(Kingdom:Genus) %>%
                     dplyr::distinct(Genus, .keep_all = T), by = c("feature" = "Genus")) %>%
  dplyr::rename(Taxon = "feature")  %>%
  dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, "\\.", " ")) %>%
  dplyr::arrange(qval) %>%
  dplyr::ungroup()



## 28__TREAT_PATH_WORM ----------------------------------------


### MaAsLin2 ----------------------------------------------------------------

#### Run Maaslin2

# Methods derived from:
#   Elementary methods provide more replicable results in microbial differential abundance analysis
#   Juho Pelto, Kari Auranen, Janne Kujala, Leo Lahti
#   https://arxiv.org/abs/2404.02691

tmp.output.path.diffAbund <- paste0("28__TREAT_PATH_WORM")
tmp.file.path.diffAbund <- file.path(path.results, "Tables/MaAsLin2", tmp.output.path.diffAbund, "significant_results.tsv")

run_maaslin2(tmp.PS = tmp.psOBJ %>%
               microViz::ps_filter(Temperature == 28),
             tmp.fixed = c("Treatment", "Pathology.Results", "Total.Worm.Count"),
             tmp.reference = c("Treatment,Control", "Pathology.Results,negative"),
             tmp.output = tmp.output.path.diffAbund
)

#### Table

diffAbnd.stats[["All"]][["28__TREAT_PATH_WORM"]][["Maaslin2"]][["output"]] <- 
  readr::read_delim(tmp.file.path.diffAbund, 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>%
  dplyr::left_join(tmp.psOBJ %>% 
                     microViz::ps_melt() %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, " ", ".")) %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, "-", ".")) %>%
                     dplyr::select(Kingdom:Genus) %>%
                     dplyr::distinct(Genus, .keep_all = T), by = c("feature" = "Genus")) %>%
  dplyr::rename(Taxon = "feature")  %>%
  dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, "\\.", " ")) %>%
  dplyr::arrange(qval) %>%
  dplyr::ungroup()


## 32__TREAT_PATH_WORM ----------------------------------------


### MaAsLin2 ----------------------------------------------------------------

#### Run Maaslin2

# Methods derived from:
#   Elementary methods provide more replicable results in microbial differential abundance analysis
#   Juho Pelto, Kari Auranen, Janne Kujala, Leo Lahti
#   https://arxiv.org/abs/2404.02691

tmp.output.path.diffAbund <- paste0("32__TREAT_PATH_WORM")
tmp.file.path.diffAbund <- file.path(path.results, "Tables/MaAsLin2", tmp.output.path.diffAbund, "significant_results.tsv")

run_maaslin2(tmp.PS = tmp.psOBJ %>%
               microViz::ps_filter(Temperature == 32),
             tmp.fixed = c("Treatment", "Pathology.Results", "Total.Worm.Count"),
             tmp.reference = c("Treatment,Control", "Pathology.Results,negative"),
             tmp.output = tmp.output.path.diffAbund
)

#### Table

diffAbnd.stats[["All"]][["32__TREAT_PATH_WORM"]][["Maaslin2"]][["output"]] <- 
  readr::read_delim(tmp.file.path.diffAbund, 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>%
  dplyr::left_join(tmp.psOBJ %>% 
                     microViz::ps_melt() %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, " ", ".")) %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, "-", ".")) %>%
                     dplyr::select(Kingdom:Genus) %>%
                     dplyr::distinct(Genus, .keep_all = T), by = c("feature" = "Genus")) %>%
  dplyr::rename(Taxon = "feature")  %>%
  dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, "\\.", " ")) %>%
  dplyr::arrange(qval) %>%
  dplyr::ungroup()


## 35__TREAT_PATH_WORM ----------------------------------------


### MaAsLin2 ----------------------------------------------------------------

#### Run Maaslin2

# Methods derived from:
#   Elementary methods provide more replicable results in microbial differential abundance analysis
#   Juho Pelto, Kari Auranen, Janne Kujala, Leo Lahti
#   https://arxiv.org/abs/2404.02691

tmp.output.path.diffAbund <- paste0("35__TREAT_PATH_WORM")
tmp.file.path.diffAbund <- file.path(path.results, "Tables/MaAsLin2", tmp.output.path.diffAbund, "significant_results.tsv")

run_maaslin2(tmp.PS = tmp.psOBJ %>%
               microViz::ps_filter(Temperature == 35),
             tmp.fixed = c("Treatment", "Pathology.Results", "Total.Worm.Count"),
             tmp.reference = c("Treatment,Control", "Pathology.Results,negative"),
             tmp.output = tmp.output.path.diffAbund
)

#### Table

diffAbnd.stats[["All"]][["35__TREAT_PATH_WORM"]][["Maaslin2"]][["output"]] <- 
  readr::read_delim(tmp.file.path.diffAbund, 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>%
  dplyr::left_join(tmp.psOBJ %>% 
                     microViz::ps_melt() %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, " ", ".")) %>%
                     dplyr::mutate(Genus = stringr::str_replace_all(Genus, "-", ".")) %>%
                     dplyr::select(Kingdom:Genus) %>%
                     dplyr::distinct(Genus, .keep_all = T), by = c("feature" = "Genus")) %>%
  dplyr::rename(Taxon = "feature")  %>%
  dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, "\\.", " ")) %>%
  dplyr::arrange(qval) %>%
  dplyr::ungroup()



## 28 vs 32 Comparison -----------------------------------------------------

# ── file paths ───────────────────────────────────────────────────────────────
pth_28 <- file.path(
  path.results,
  "Tables/MaAsLin2/28__TREAT_PATH_WORM/28__significant_results.tsv"
)

pth_32 <- file.path(
  path.results,
  "Tables/MaAsLin2/32__TREAT_PATH_WORM/32__significant_results.tsv"
)

# ── helper 1: read & clean one MaAsLin2 sheet ────────────────────────────────
read_maaslin <- function(path, temp_label) {
  readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::filter(qval < 0.05) |>                  # keep FDR-significant hits
    dplyr::mutate(
      temp_group = temp_label,
      direction  = ifelse(coef > 0, "↑", "↓")
    ) |>
    dplyr::distinct(feature, metadata, .keep_all = TRUE) |>
    dplyr::select(temp_group, metadata, feature, direction, coef, qval)
}

df_28 <- read_maaslin(pth_28, "28 °C")
df_32 <- read_maaslin(pth_32, "32 °C")

# ── helper 2: comparison counts ──────────────────────────────────────────────
compare_temps <- function(df28, df32, metadata_var) {
  d28 <- df28 |> dplyr::filter(metadata == metadata_var)
  d32 <- df32 |> dplyr::filter(metadata == metadata_var)
  
  total_28 <- dplyr::n_distinct(d28$feature)
  total_32 <- dplyr::n_distinct(d32$feature)
  
  shared <- dplyr::inner_join(
    d28, d32, by = "feature", suffix = c(".28", ".32")
  )
  
  dplyr::tibble(
    category = c(
      "Total significant (28 °C)",
      "Total significant (32 °C)",
      "Shared, same direction",
      "Shared, opposite direction",
      "Unique to 28 °C",
      "Unique to 32 °C"
    ),
    n_taxa = c(
      total_28,
      total_32,
      sum(shared$direction.28 == shared$direction.32),
      sum(shared$direction.28 != shared$direction.32),
      total_28 - nrow(shared),
      total_32 - nrow(shared)
    )
  )
}

# ── comparison gt tables ────────────────────────────────────────────────────
# tbl_cmp_treatment <- compare_temps(df_28, df_32, "Treatment") |>
#   gt::gt() |>
#   gt::tab_header(
#     title    = md("Burden-Associated Genera: 28 °C vs 32 °C"),
#     subtitle = md("**Metadata:** Treatment  |  *FDR < 0.05*")
#   )

diffAbnd.stats[["All"]][["28_vs_32__TREAT_PATH_WORM"]][["Maaslin2"]][["Comparison__28_v_32"]] <- 
  compare_temps(df_28, df_32, "Total.Worm.Count") |>
  gt::gt() |>
  gt::tab_header(
    title    = md("Burden-Associated Genera: 28 °C vs 32 °C"),
    subtitle = md("**Metadata:** Total Worm Count  |  *FDR < 0.05*")
  )

# ── helper 3: top-N taxa with q-value display ───────────────────────────────
topN_taxa <- function(df28, df32, metadata_var, n_top = 10) {
  
  pooled <- dplyr::bind_rows(df28, df32) |>
    dplyr::filter(metadata == metadata_var)
  
  pooled_unique <- pooled |>
    dplyr::group_by(feature) |>
    dplyr::slice_min(qval, with_ties = FALSE) |>
    dplyr::ungroup()
  
  pooled_unique |>
    dplyr::arrange(qval) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::transmute(
      Rank        = dplyr::row_number(),
      Genus       = feature,
      Temperature = temp_group,
      Direction   = direction,
      Coef        = coef,
      qval        = qval
    )
}

fmt_qval <- function(gt_tbl) {
  gt_tbl |>
    # First: numeric formatting to 3 decimals
    gt::fmt_number(columns = qval, decimals = 3) |>
    # Second: override rows where qval < 0.001 with "<0.001"
    gt::fmt(
      columns = qval,
      rows    = qval < 0.001,
      fns     = function(x) rep("<0.001", length(x))
    )
}

# ── top-10 gt tables with custom q-value formatting ─────────────────────────
# tbl_top10_treatment <- topN_taxa(df_28, df_32, "Treatment") |>
#   gt::gt() |>
#   fmt_qval() |>
#   gt::fmt_number(columns = Coef, decimals = 3) |>
#   gt::tab_header(
#     title    = md("Top 10 Genera — Treatment"),
#     subtitle = md("*Smallest *q*; 28 °C + 32 °C*")
#   )

diffAbnd.stats[["All"]][["28_vs_32__TREAT_PATH_WORM"]][["Maaslin2"]][["Top10__28_v_32"]] <- 
  topN_taxa(df_28, df_32, "Total.Worm.Count") |>
  gt::gt() |>
  fmt_qval() |>
  gt::fmt_number(columns = Coef, decimals = 3) |>
  gt::tab_header(
    title    = md("Top 10 Genera — Total Worm Count"),
    subtitle = md("*Smallest *q*; 28 °C + 32 °C*")
  )



## Random Forest -----------------------------------------------------------



### Prep data ---------------------------------------------------------------



tmp.psOBJ <- ps.list[["All"]] %>%
  microViz::ps_filter(Treatment == "Exposed")

# Set seed for reproducibility
set.seed(42)

tmp.output.path.diffAbund <- paste0("Exp__WORM")
tmp.file.path.diffAbund <- read_tsv(file.path(path.results, "Tables/MaAsLin2", tmp.output.path.diffAbund, "significant_results.tsv"))

# Get significant taxa names
significant_taxa <- tmp.file.path.diffAbund %>%
  dplyr::filter(qval < 0.1) %>%  # Filter for significant taxa (q < 0.1)
  dplyr::pull(feature) %>%
  unique()

# Prepare microbiome data
microbiome_data <- tmp.psOBJ %>%
  microViz::ps_filter(Treatment == "Exposed") %>%
  microViz::tax_agg(rank = "Genus") %>%
  # Apply compositional normalization (total sum scaling) and log2 transformation
  microViz::tax_transform("compositional", rank = "Genus") %>%
  microViz::tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>%
  microViz::ps_melt() %>%
  dplyr::filter(Genus %in% significant_taxa) %>%
  dplyr::select(Sample, Genus, Abundance) %>%
  tidyr::pivot_wider(
    names_from = Genus,
    values_from = Abundance,
    values_fill = 0
  )

# Get metadata
metadata <- tmp.psOBJ %>%
  microViz::ps_filter(Treatment == "Exposed") %>%
  microViz::ps_melt() %>%
  dplyr::select(Sample, Total.Worm.Count, Temperature, DPE) %>%
  dplyr::distinct()

# Combine microbiome data with metadata
rf_data <- metadata %>%
  dplyr::left_join(microbiome_data, by = "Sample") %>%
  tidyr::drop_na()  # Remove any rows with missing values



### Test and Train ----------------------------------------------------------


# Split data into training and testing sets
set.seed(42)
train_index <- caret::createDataPartition(rf_data$Total.Worm.Count, p = 0.8, list = FALSE)
train_data <- rf_data[train_index, ]
test_data <- rf_data[-train_index, ]

# Train random forest model
set.seed(42)
rf_model <- ranger::ranger(
  formula = Total.Worm.Count ~ .,
  data = train_data %>% dplyr::select(-Sample),
  num.trees = 1000,
  importance = "permutation",
  num.threads = parallel::detectCores() - 1  # Use all but one core
)



### Evaluate Performance ----------------------------------------------------



# Evaluate model performance
predictions <- predict(rf_model, data = test_data %>% dplyr::select(-Sample))$predictions
rmse <- sqrt(mean((test_data$Total.Worm.Count - predictions)^2))
r2 <- 1 - sum((test_data$Total.Worm.Count - predictions)^2) / sum((test_data$Total.Worm.Count - mean(test_data$Total.Worm.Count))^2)

# Get variable importance
importance_data <- as.data.frame(rf_model$variable.importance) %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::rename(importance = `rf_model$variable.importance`) %>%
  dplyr::arrange(dplyr::desc(importance))

# Create performance metrics table
performance_tbl <- tibble::tibble(
  Metric = c("RMSE", "R-squared"),
  Value = c(round(rmse, 3), round(r2, 3))
)


#### Table

randomForest.stats[["Exposed"]][["RandomForest"]][["Performance"]] <- 

  # Create summary table
  gt::gt(performance_tbl) %>%
    gt::tab_header(
      title = "Random Forest Model Performance",
      subtitle = "Predicting Worm Burden from Microbiome Data"
    ) %>%
    gt::fmt_number(
      columns = Value,
      decimals = 3
  )

# # Create top variables table
# top_vars_tbl <- importance_data %>%
#   head(10) %>%
#   dplyr::mutate(
#     Rank = row_number(),
#     Importance = round(importance, 3)
#   ) %>%
#   dplyr::select(Rank, Variable = variable, Importance)
# 
# gt::gt(top_vars_tbl) %>%
#   gt::tab_header(
#     title = "Top 10 Most Important Variables",
#     subtitle = "Based on % Increase in MSE"
#   ) %>%
#   gt::fmt_number(
#     columns = Importance,
#     decimals = 3
#   )

# Create top 25 variables table
top25_vars_tbl <- importance_data %>%
  head(25) %>%
  dplyr::mutate(
    Rank = row_number(),
    Importance = round(importance, 3)
  ) %>%
  dplyr::select(Rank, Variable = variable, Importance)

#### Table

randomForest.stats[["Exposed"]][["RandomForest"]][["Top25_Features"]] <- 
  gt::gt(top25_vars_tbl) %>%
    gt::tab_header(
      title = "Top 25 Most Important Features",
      subtitle = "Based on % Increase in MSE"
    ) %>%
    gt::fmt_number(
      columns = Importance,
      decimals = 3
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "lightgray"),
      locations = gt::cells_body(
        rows = seq(1, 25, 2)  # Alternate row shading
      )
    )


### Stability Analysis ------------------------------------------------------

# Number of iterations
n_iter <- 100
# Store top 10 genera for each iteration
set.seed(42) # For reproducibility of the seeds themselves
seeds <- sample(1:10000, n_iter, replace = FALSE)

top10_list <- vector("list", n_iter)

# Progress bar
pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

for (i in seq_len(n_iter)) {
  set.seed(seeds[i])
  # Split data
  train_index <- caret::createDataPartition(rf_data$Total.Worm.Count, p = 0.8, list = FALSE)
  train_data <- rf_data[train_index, ]
  test_data <- rf_data[-train_index, ]
  # Train model
  rf_model <- ranger::ranger(
    formula = Total.Worm.Count ~ .,
    data = train_data %>% dplyr::select(-Sample),
    num.trees = 1000,
    importance = "permutation",
    num.threads = parallel::detectCores() - 1
  )
  # Get top 10 genera by importance
  importance_data <- as.data.frame(rf_model$variable.importance) %>%
    tibble::rownames_to_column("variable") %>%
    dplyr::rename(importance = `rf_model$variable.importance`) %>%
    dplyr::arrange(dplyr::desc(importance))
  top10 <- importance_data$variable[1:10]
  top10_list[[i]] <- top10
  utils::setTxtProgressBar(pb, i)
}
close(pb)

# Summarize how often each genus appears in the top 10
all_top10 <- unlist(top10_list)
top10_freq <- tibble::tibble(
  Genus = all_top10
) %>%
  dplyr::count(Genus, name = "Times_in_Top10") %>%
  dplyr::arrange(dplyr::desc(Times_in_Top10))

#### Table

randomForest.stats[["Exposed"]][["RandomForest"]][["StabilityAnalysis_Top10"]] <- 
  # Present as gt table
  gt::gt(top10_freq) %>%
    gt::tab_header(
      title = "Stability Analysis: Top Feature Frequency",
      subtitle = paste0("Number of times each feature appeared in the top 10 across ", n_iter, " runs.")
    )


### Cross-Validation Analysis -----------------------------------------------

# Set up cross-validation
set.seed(42)
n_folds <- 10
cv_folds <- caret::createFolds(rf_data$Total.Worm.Count, k = n_folds, list = TRUE)

# Initialize storage for metrics
cv_metrics <- tibble::tibble(
  Fold = integer(),
  RMSE = numeric(),
  R2 = numeric(),
  Null_RMSE = numeric(),
  Null_R2 = numeric()
)

# Progress bar
pb <- utils::txtProgressBar(min = 0, max = n_folds, style = 3)

# Run cross-validation
for (i in seq_along(cv_folds)) {
  # Split data
  test_indices <- cv_folds[[i]]
  train_cv <- rf_data[-test_indices, ]
  test_cv <- rf_data[test_indices, ]
  
  # Train RF model
  set.seed(42 + i)  # Different seed for each fold
  rf_cv <- ranger::ranger(
    formula = Total.Worm.Count ~ .,
    data = train_cv %>% dplyr::select(-Sample),
    num.trees = 1000,
    importance = "permutation",
    num.threads = parallel::detectCores() - 1
  )
  
  # Get predictions
  rf_pred <- predict(rf_cv, data = test_cv %>% dplyr::select(-Sample))$predictions
  
  # Calculate null model predictions (mean of training set)
  null_pred <- mean(train_cv$Total.Worm.Count)
  
  # Calculate metrics for RF model
  rf_rmse <- sqrt(mean((test_cv$Total.Worm.Count - rf_pred)^2))
  rf_r2 <- 1 - sum((test_cv$Total.Worm.Count - rf_pred)^2) / 
    sum((test_cv$Total.Worm.Count - mean(test_cv$Total.Worm.Count))^2)
  
  # Calculate metrics for null model
  null_rmse <- sqrt(mean((test_cv$Total.Worm.Count - null_pred)^2))
  null_r2 <- 1 - sum((test_cv$Total.Worm.Count - null_pred)^2) / 
    sum((test_cv$Total.Worm.Count - mean(test_cv$Total.Worm.Count))^2)
  
  # Store metrics
  cv_metrics <- cv_metrics %>%
    dplyr::add_row(
      Fold = i,
      RMSE = rf_rmse,
      R2 = rf_r2,
      Null_RMSE = null_rmse,
      Null_R2 = null_r2
    )
  
  utils::setTxtProgressBar(pb, i)
}
close(pb)

# Calculate summary statistics
cv_summary <- cv_metrics %>%
  dplyr::summarise(
    Mean_RMSE = mean(RMSE),
    SD_RMSE = sd(RMSE),
    Mean_R2 = mean(R2),
    SD_R2 = sd(R2),
    Mean_Null_RMSE = mean(Null_RMSE),
    SD_Null_RMSE = sd(Null_RMSE),
    Mean_Null_R2 = mean(Null_R2),
    SD_Null_R2 = sd(Null_R2)
  )

# Create summary table
cv_summary_gt <- tibble::tibble(
  Model = c("Random Forest", "Null Model"),
  RMSE = c(cv_summary$Mean_RMSE, cv_summary$Mean_Null_RMSE),
  R2 = c(cv_summary$Mean_R2, cv_summary$Mean_Null_R2)
) %>%
  # Create gt table
  gt::gt() %>%
  gt::tab_header(
    title = "Cross-Validation Summary Statistics",
    subtitle = paste0("Based on ", n_folds, "-fold cross-validation")
  ) %>%
  gt::fmt_number(
    columns = c(RMSE, R2),
    decimals = 3
  )

# Display tables
cv_metrics_gt <- cv_metrics %>%
  dplyr::mutate(
    RMSE = round(RMSE, 3),
    R2 = round(R2, 3),
    Null_RMSE = round(Null_RMSE, 3),
    Null_R2 = round(Null_R2, 3)
  ) %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Cross-Validation Metrics by Fold",
    subtitle = "Comparing Random Forest vs. Null Model (Mean Prediction)"
  ) %>%
  gt::fmt_number(
    columns = c(RMSE, R2, Null_RMSE, Null_R2),
    decimals = 3
  ) %>%
  gt::tab_spanner(
    label = "Random Forest",
    columns = c(RMSE, R2)
  ) %>%
  gt::tab_spanner(
    label = "Null Model",
    columns = c(Null_RMSE, Null_R2)
  )


#### Table

randomForest.stats[["Exposed"]][["RandomForest"]][["CrossValidation_Summary_10Fold"]] <- 
  cv_summary_gt

# Calculate improvement over null model
improvement <- cv_summary %>%
  dplyr::mutate(
    RMSE_Improvement = (Mean_Null_RMSE - Mean_RMSE) / Mean_Null_RMSE * 100,
    R2_Improvement = (Mean_R2 - Mean_Null_R2) * 100
  )

# Create improvement table
improvement_tbl <- tibble::tibble(
  Metric = c("RMSE", "R²"),
  Improvement = c(
    paste0(round(improvement$RMSE_Improvement, 1), "% reduction"),
    paste0(round(improvement$R2_Improvement, 1), "% increase")
  ),
  Description = c(
    "Reduction in Root Mean Square Error",
    "Increase in Coefficient of Determination"
  )
)

randomForest.stats[["Exposed"]][["RandomForest"]][["CrossValidation_Improvement"]] <- 
# Display as gt table
  gt::gt(improvement_tbl) %>%
    gt::tab_header(
      title = "Random Forest Model Improvement Over Null Model",
      subtitle = "Based on 10-fold cross-validation results"
    ) %>%
    gt::cols_label(
      Metric = "Performance Metric",
      Improvement = "Improvement",
      Description = "Description"
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "lightgreen"),
      locations = gt::cells_body(
        rows = 1:2
      )
    )



# Add at the end
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "secs")
cat("Completed 07__TaxonAbund_Unexp_Exp__Stats.R at", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time:", round(duration, 2), "seconds\n")