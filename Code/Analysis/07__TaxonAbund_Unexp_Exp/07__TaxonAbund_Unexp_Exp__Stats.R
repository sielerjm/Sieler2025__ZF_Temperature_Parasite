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


# ── paths ────────────────────────────────────────────────────────────────────
pth_28 <- "/Users/michaelsieler/Dropbox/Mac (2)/Documents/Sharpton_Lab/Projects_Repository/Rules_of_Life/Sieler2024__ZF_Temperature_Parasite/Results/Tables/MaAsLin2/28__TREAT_PATH_WORM/28__significant_results.tsv"
pth_32 <- "/Users/michaelsieler/Dropbox/Mac (2)/Documents/Sharpton_Lab/Projects_Repository/Rules_of_Life/Sieler2024__ZF_Temperature_Parasite/Results/Tables/MaAsLin2/32__TREAT_PATH_WORM/32__significant_results.tsv"

# ── helper to read & label one sheet ─────────────────────────────────────────
read_maaslin <- function(path, temp_label){
  readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::filter(qval < 0.05) |>              # keep only significant taxa
    dplyr::select(feature, coef) |>            # 'feature' = taxon name
    dplyr::mutate(
      direction = ifelse(coef > 0, "↑", "↓"),
      temp_group = temp_label
    )
}

df_28 <- read_maaslin(pth_28, "28 °C")
df_32 <- read_maaslin(pth_32, "32 °C")

# ── compare taxa between temperatures ────────────────────────────────────────
comparison_tbl <- dplyr::full_join(
  df_28 |> dplyr::select(feature, dir28 = direction),
  df_32 |> dplyr::select(feature, dir32 = direction),
  by = "feature"
) |>
  dplyr::mutate(
    category = dplyr::case_when(
      !is.na(dir28) & !is.na(dir32) & dir28 == dir32 ~ "Shared, same direction",
      !is.na(dir28) & !is.na(dir32) & dir28 != dir32 ~ "Shared, opposite direction",
      !is.na(dir28) &  is.na(dir32)                 ~ "Unique to 28 °C",
      is.na(dir28)  & !is.na(dir32)                 ~ "Unique to 32 °C"
    )
  )

# ── summary counts for quick inspection ─────────────────────────────────────
summary_tbl <- comparison_tbl |>
  dplyr::count(category, name = "n_taxa") |>
  dplyr::arrange(match(category,
                       c("Shared, same direction",
                         "Shared, opposite direction",
                         "Unique to 28 °C",
                         "Unique to 32 °C")))

# ── pretty table ─────────────────────────────────────────────────────────────
gt::gt(summary_tbl) |>
  gt::tab_header(
    title = "Overlap of FDR-Significant Genera (28 °C vs 32 °C)",
    subtitle = "Direction refers to sign of MaAsLin2 coefficient"
  )

# Add end message
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "secs")
cat("Completed 07__TaxonAbund_Unexp_Exp__Stats.R at", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time:", round(duration, 2), "seconds\n")



