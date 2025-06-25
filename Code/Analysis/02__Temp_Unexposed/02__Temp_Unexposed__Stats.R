# 02__Temp_Unexposed__Stats --------------------------------------------------

# Add start message with timestamp
start_time <- Sys.time()
cat("Starting 02__Temp_Unexposed__Stats.R at", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

tmp.psOBJ <- ps.list[["Unexposed"]]
tmp.resSubSection <- "Unexposed"

## TEMP -------------------------------------------------------------

### Alpha -------------------------------------------------------------------

#### GLM ---------------------------------------------------------------------

alpha.stats[[tmp.resSubSection]][["TEMP"]][["GLM"]] <- 
  tmp.psOBJ %>%
  
  # Convert phyloseq object into a dataframe and pivot longer by Alpha Metric and Score
  psObjToDfLong(div.score = "Alpha.Score", div.metric = "Alpha.Metric") %>%
  
  # Clean up cell value names by removing any strings including and after "__"
  cutCellNames(col = "Alpha.Metric", sep = "__") %>%
  
  # Run Tukey Test
  run_glm_models(formula_str = "Alpha.Score ~ Temperature")

### Table

alpha.stats[[tmp.resSubSection]][["TEMP"]][["GLM.Table"]] <-
  alpha.stats[[tmp.resSubSection]][["TEMP"]][["GLM"]] %>%
  
  # Create a column called "Alpha.Metric" for each metric's GLM model
  purrr::imap(., ~tidy(.x) %>% mutate(Alpha.Metric = .y)) %>%
  
  # Combine the metric dataframes
  dplyr::bind_rows() %>%
  
  # Add significance indicators in a new column
  SigStars() %>%
  
  # Create GT Table
  set_GT(var = "p.value", group.by = "Alpha.Metric") %>%
  
  # Title/caption
  gt::tab_header(
    title = "GLM Results",
    subtitle = "glm(Alpha.Score ~ Temperature); Unexposed fish"
  )


#### ANOVA -------------------------------------------------------------------

alpha.stats[[tmp.resSubSection]][["TEMP"]][["ANOVA"]] <-
  run_glm_anova(alpha.stats[[tmp.resSubSection]][["TEMP"]][["GLM"]])

### Table

alpha.stats[[tmp.resSubSection]][["TEMP"]][["ANOVA.Table"]] <-
  alpha.stats[[tmp.resSubSection]][["TEMP"]][["ANOVA"]] %>%
  
  # Create GT Table
  set_GT(var = "p.value", group.by = "Alpha.Metric") %>%
  
  # Title/caption
  gt::tab_header(
    title = "ANOVA of GLM",
    subtitle = "ANOVA(GLM(Alpha.Score ~ Temperature), type = 2); Unexposed fish"
  )



#### TUKEY -------------------------------------------------------------------

alpha.stats[[tmp.resSubSection]][["TEMP"]][["Tukey"]] <-
  tmp.psOBJ %>%
  
  # Convert phyloseq object into a dataframe and pivot longer by Alpha Metric and Score
  psObjToDfLong(div.score = "Alpha.Score", div.metric = "Alpha.Metric") %>%
  
  # Clean up cell value names by removing any strings including and after "__"
  cutCellNames(col = "Alpha.Metric", sep = "__") %>%
  
  # Run Tukey test on GLM models
  run_tukey_glm(., "Alpha.Score", "Alpha.Metric", c("Temperature"), 
                group_by_var = NULL) 

### Table

alpha.stats[[tmp.resSubSection]][["TEMP"]][["Tukey.Table"]] <-
  alpha.stats[[tmp.resSubSection]][["TEMP"]][["Tukey"]] %>%
  
  # Combine the different alpha diversity metrics into one dataframe
  dplyr::bind_rows() %>%
  
  # Create the table
  dplyr::group_by(Alpha.Metric, term) %>%
  set_GT(var = "adj.p.value", group.by = "Alpha.Metric") %>%
  
  # Title/caption
  gt::tab_header(
    title = "Pairwise Tukey's HSD, p.adj: Dunnett",
    subtitle = "Tukey(Alpha.Score ~ Temperature); Unexposed fish"
  )


### Beta --------------------------------------------------------------------



#### Capscale ----------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP"]][["CAP.mod"]] <-
  run_capscale(tmp.psOBJ, 
               dist.matrix = beta.dist.mat[[tmp.resSubSection]], 
               formula_str = "dist ~ Temperature")



##### ADONIS ------------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP"]][["CAP.ADONIS"]] <-
  run_cap_adonis(tmp.psOBJ,
                 dist.matrix = beta.dist.mat[[tmp.resSubSection]], 
                 formula_str = "dist ~ Temperature",
                 by.method = "terms") 

### Table

beta.stats[[tmp.resSubSection]][["TEMP"]][["CAP.ADONIS.Table"]] <-
  beta.stats[[tmp.resSubSection]][["TEMP"]][["CAP.ADONIS"]] %>%
    set_GT(var = "p.value", group.by = "Beta.Metric") %>%
  
  # Title/caption
  gt::tab_header(
    title = "ADONIS2",
    subtitle = "adonis2(Beta Distance ~ Temperature); Unexposed fish"
  )


#### Dispersion --------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.model"]] <-
  run_BetaDispersion(dist.matrix = beta.dist.mat[[tmp.resSubSection]], 
                     var = c("Temperature"))

##### ANOVA ------------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.ANOVA"]] <-
  get_HoD_anova(betaDisper = beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.model"]],
                var = c("Temperature"))

### Table

beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.ANOVA.Table"]] <-
  beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.ANOVA"]] %>%
  
    # Create the table
    dplyr::group_by(Beta.Metric) %>%
    set_GT(var = "p.value", group.by = "Beta.Metric")  %>%
  
    # Title/caption
    gt::tab_header(
      title = "ANOVA: Homogeneity of Dispersion",
      subtitle = "ANOVA(Beta Disperson ~ Temperature); Unexposed fish"
    )

##### Tukey ------------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.Tukey"]] <-
  get_HoD_tukey(betaDisper = beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.model"]],
                var = c("Temperature")) 

### Table

beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.Tukey.Table"]] <-
  beta.stats[[tmp.resSubSection]][["TEMP"]][["HoD.Tukey"]] %>%
  
  # Create the table
  dplyr::group_by(Beta.Metric) %>%
  set_GT(var = "adj.p.value", group.by = "Beta.Metric")  %>%
  
  # Title/caption
  gt::tab_header(
    title = "Tukey: Homogeneity of Dispersion",
    subtitle = "Tukey(Beta Disperson ~ Temperature); Unexposed fish"
  )


## TEMP:DPE ---------------------------------------------------------------

### Alpha -------------------------------------------------------------------

#### GLM ---------------------------------------------------------------------

alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]] <- 
  tmp.psOBJ %>%
  
  # Convert phyloseq object into a dataframe and pivot longer by Alpha Metric and Score
  psObjToDfLong(div.score = "Alpha.Score", div.metric = "Alpha.Metric") %>%
  
  # Clean up cell value names by removing any strings including and after "__"
  cutCellNames(col = "Alpha.Metric", sep = "__") %>%
  
  # Run Tukey Test
  run_glm_models(formula_str = "Alpha.Score ~ Temperature*DPE")

### Table

alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM.Table"]] <-
  alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]] %>%
  
  # Create a column called "Alpha.Metric" for each metric's GLM model
  purrr::imap(., ~tidy(.x) %>% mutate(Alpha.Metric = .y)) %>%
  
  # Combine the metric dataframes
  dplyr::bind_rows() %>%
  
  # Add significance indicators in a new column
  SigStars() %>%
  
  # Create GT Table
  set_GT(var = "p.value", group.by = "Alpha.Metric")%>%
  
  # Title/caption
  gt::tab_header(
    title = "GLM Results",
    subtitle = "glm(Alpha.Score ~ Temperature*DPE); Unexposed fish"
  )


#### ANOVA -------------------------------------------------------------------

alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["ANOVA"]] <-
  run_glm_anova(alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]])

### Table

alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["ANOVA.Table"]] <-
  alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["ANOVA"]] %>%
  
  # Create GT Table
  set_GT(var = "p.value", group.by = "Alpha.Metric") %>%
  
  # Title/caption
  gt::tab_header(
    title = "ANOVA of GLM",
    subtitle = "ANOVA(GLM(Alpha.Score ~ Temperature*DPE), type = 2); Unexposed fish"
  )


#### TUKEY -------------------------------------------------------------------

alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["Tukey"]] <-
  tmp.psOBJ %>%
  
  # Convert phyloseq object into a dataframe and pivot longer by Alpha Metric and Score
  psObjToDfLong(div.score = "Alpha.Score", div.metric = "Alpha.Metric") %>%
  
  # Clean up cell value names by removing any strings including and after "__"
  cutCellNames(col = "Alpha.Metric", sep = "__") %>%
  
  # Run Tukey test on GLM models
  run_tukey_glm(., "Alpha.Score", "Alpha.Metric", c("Temperature", "DPE"), 
                group_by_var = "Temperature") 

### Table

alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["Tukey.Table"]] <-
  alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["Tukey"]] %>%
  
  # Combine the different alpha diversity metrics into one dataframe
  dplyr::bind_rows() %>%
  
  # Create the table
  dplyr::group_by(Alpha.Metric, term) %>%
  set_GT(var = "adj.p.value", group.by = "Alpha.Metric") %>%
  
  # Title/caption
  gt::tab_header(
    title = "Pairwise Tukey's HSD, p.adj: Dunnett",
    subtitle = "Tukey(Alpha.Score ~ Temperature*DPE); Unexposed fish"
  )

### Beta --------------------------------------------------------------------



#### Capscale ----------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["CAP.mod"]] <-
  run_capscale(tmp.psOBJ, 
               dist.matrix = beta.dist.mat[[tmp.resSubSection]], 
               formula_str = "dist ~ Temperature*DPE")



##### ADONIS ------------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["CAP.ADONIS"]] <-
  run_cap_adonis(tmp.psOBJ,
                 dist.matrix = beta.dist.mat[[tmp.resSubSection]], 
                 formula_str = "dist ~ Temperature*DPE",
                 by.method = "terms") 

### Table

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["CAP.ADONIS.Table"]] <-
  beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["CAP.ADONIS"]] %>%
  set_GT(var = "p.value", group.by = "Beta.Metric") %>%
  
  # Title/caption
  gt::tab_header(
    title = "ADONIS2",
    subtitle = "adonis2(Beta Distance ~ Temperature*DPE); Unexposed fish"
  )


#### Dispersion --------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.model"]] <-
  run_BetaDispersion(dist.matrix = beta.dist.mat[[tmp.resSubSection]], 
                     var = c("Temp.DPE"))

##### ANOVA ------------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.ANOVA"]] <-
  get_HoD_anova(betaDisper = beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.model"]],
                var = c("Temp.DPE"))

### Table

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.ANOVA.Table"]] <-
  beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.ANOVA"]] %>%
  
  # Create the table
  dplyr::group_by(Beta.Metric) %>%
  set_GT(var = "p.value", group.by = "Beta.Metric")  %>%
  
  # Title/caption
  gt::tab_header(
    title = "ANOVA: Homogeneity of Dispersion",
    subtitle = "ANOVA(Beta Disperson ~ Temperature*DPE); Unexposed fish"
  )

##### Tukey ------------------------------------------------------------------

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.Tukey"]] <-
  get_HoD_tukey(betaDisper = beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.model"]],
                var = c("Temp.DPE")) 

### Table

beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.Tukey.Table"]] <-
  beta.stats[[tmp.resSubSection]][["TEMP:DPE"]][["HoD.Tukey"]] %>%
  
  # Create the table
  dplyr::group_by(Beta.Metric) %>%
  set_GT(var = "adj.p.value", group.by = "Beta.Metric")  %>%
  
  # Title/caption
  gt::tab_header(
    title = "Tukey: Homogeneity of Dispersion",
    subtitle = "Tukey(Beta Disperson ~ Temperature*DPE); Unexposed fish"
  )





# Add end message with timestamp and duration
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "secs")
cat("Completed 02__Temp_Unexposed__Stats.R at", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time:", round(duration, 2), "seconds\n")












