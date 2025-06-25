# Libraries ---------------------------------------------------------------
#   Organized by category and alphabet


# Style Guide -------------------------------------------------------------

#   library(<package_name>)  # Brief description


# General -----------------------------------------------------------------

library(forcats)  # suite of useful tools that solve common problems with factors
library(knitr)  # For knitting documents to HTML or PDF formats
library(parallel)  # Parallel processing functions
library(rcompanion)  # Allows us to run transformTukey()
library(scales)  # Normalizing values across a given scale
library(purrr)
library(furrr)


# Plotting ----------------------------------------------------------------

library(ggplot2)  # For plotting pretty graphs
library(ggbeeswarm)  # organizes the jitters better on boxplots
library(ggExtra)  # Additional ggplot settings
library(ggrepel)  # Adds label buffering to labels so they don't overlap, geom_label_repel()
library(ggpubr) # Statistical plotting 
library(gt)  # Making pretty tables
library(RColorBrewer)  # Making and using pretty color palettes
# library(rstatix) 
library(ComplexHeatmap)  # microViz rec
library(ggnewscale)  # Allows you to add an additional color/fill scale "new_scale_color()" https://github.com/eliocamp/ggnewscale


# Microbiome --------------------------------------------------------------

library(Maaslin2)  # Differential abundance analysis
library(phyloseq)  # Microbiome analysis package
library(phyloseqCompanion)  # Helper functions for phyloseq
library(vegan)  # Ecological analysis
library(microbiome) # microbiome analysis
library(picante) # alpha diversity analysis
library(microViz) # Microbiome analysis 


# Statistics --------------------------------------------------------------

library(broom)  # Convert Statistical Objects into Tidy Tibbles
library(car)  # Allows us to run Anova()
# library(easystats)  # Compilation of statistical packages
library(lme4)  # Linear mixed models
library(emmeans)  # Post-hoc model comparisons
library(multcomp)  # Multiple comparison tests
library(caret) # random forest
library(vip)  # random forest
library(ranger)  # random forest


# KEEP AS LAST LIBRARY TO LOAD --------------------------------------------

# Load last because of potential conflicting function names in previous libs
library(tidyverse)  # Making your code look pretty and tidy

cacheing=TRUE


# Install packages (placeholder)

# BiocManager::install(c("tidyverse"))


