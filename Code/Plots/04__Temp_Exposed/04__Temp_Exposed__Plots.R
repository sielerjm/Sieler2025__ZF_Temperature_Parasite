
# 04__Temp_Exposed__Plots -----------------------------------------------------

cat("04__Temp_Exposed__Plots \n")

tmp.psOBJ <- ps.list[["Exposed"]]
tmp.resSubSection <- "Exposed"


## Alpha ---------------------------------------------------------------


### TEMP --------------------------------------------------------------------

#### Sig Labels --------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP"]][["Tukey"]][["Sig.Labels"]] <- 
  alpha.stats[[tmp.resSubSection]][["TEMP"]][["Tukey"]] %>%
  SigStars(pval.var = "adj.p.value") %>%
  dplyr::filter(p.value < 0.05 ) %>%
  dplyr::filter(Alpha.Metric == "Simpson")


#### Plot --------------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP"]][["TUKEY"]][["Plot"]] <-
  tmp.psOBJ %>%
  microViz::samdat_tbl() %>%
  tidyr::pivot_longer(cols = contains("_norm"), names_to = "Alpha.Metric", values_to = "Alpha.Score") %>%
  cutCellNames(col = "Alpha.Metric") %>%
  filter(Alpha.Metric == "Simpson") %>%
  # filter(Treatment == "Control") %>%
  # dplyr::mutate(DPE = as.numeric(levels(DPE)[DPE])) %>%
  ggplot(aes(x = Temperature, y = Alpha.Score)) +
  
  
  geom_violin(aes(#color = Temperature,
    fill = Temperature,
    group = interaction(Temperature)),
    # position = position_dodge(width = 5),
    # width = 4,
    color = "white",
    draw_quantiles = c(.5),
    size = 3,
    alpha = .05,
    scale = "width",
    show.legend = F,
    na.rm = T) +
  
  geom_violin(aes(color = Temperature,
                  fill = Temperature,
                  group = interaction(Temperature)),
              # position = position_dodge(width = 5),
              # width = 4,
              size = 1.25,
              alpha = 0,
              scale = "width",
              show.legend = F,
              na.rm = T) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  
  geom_violin(aes(color = Temperature,
                  fill = Temperature,
                  group = interaction(Temperature)),
              # position = position_dodge(width = 5),
              # width = 4,
              draw_quantiles = c(.5),
              size = 1.25,
              alpha = 0,
              scale = "width",
              show.legend = F,
              na.rm = T) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  ggbeeswarm::geom_quasirandom(aes(
    color = Temperature),
    size = 3,
    fill = "white",
    shape = 23,
    stroke = 2,
    show.legend = T,
    # position = position_dodge2(width = 5),
    na.rm = T) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
  
  scale_x_discrete(
    labels = function(x) paste0(x, "°C"), # Append '°C' to x-axis labels
    breaks = c('28', '32', '35')
  ) +
  scale_y_continuous(limits = c(-0.05, 1.15), breaks = seq(0, 1, by = .25)) +
  labs(title = "Alpha Diversity Scores across Temperature", #caption = "Alpha; All",
       x = "Temperature",
       y = "Alpha Score (Normalized)"
  ) +
  guides(
    color = guide_legend(order = 1, title = "Temp (°C)",
                         override.aes = list(shape = 16, size = 5, linetype = NA))
  ) +
  theme(legend.position = "right",
        legend.justification = "center",
        legend.direction = "vertical",
        panel.grid.major.x = element_blank()) +
  
  # Statistical Notation
  ggpubr::stat_pvalue_manual(alpha.plots[[tmp.resSubSection]][["TEMP"]][["Tukey"]][["Sig.Labels"]],
                             label = "p.adj.sig",
                             y.position = c(1.075, 1.15),
                             size = 6,
                             label.size = 10,
                             bracket.size = 2,
                             hide.ns = T)


### TEMP:DPE ----------------------------------------------------------------

#### Sig Labels --------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP:DPE"]][["Tukey"]][["Sig.Labels"]] <- 
  alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["Tukey"]] %>%
  SigStars(pval.var = "adj.p.value") %>%
  dplyr::filter(p.value < 0.05 ) %>%
  dplyr::filter(Alpha.Metric == "Simpson") 


#### Plot --------------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP:DPE"]][["TUKEY"]][["Plot"]] <-
  tmp.psOBJ %>%
  microViz::samdat_tbl() %>%
  tidyr::pivot_longer(cols = contains("_norm"), names_to = "Alpha.Metric", values_to = "Alpha.Score") %>%
  cutCellNames(col = "Alpha.Metric") %>%
  filter(Alpha.Metric == "Simpson") %>%
  # filter(Treatment == "Control") %>%
  # dplyr::mutate(DPE = as.numeric(levels(DPE)[DPE])) %>%
  ggplot(aes(x = DPE, y = Alpha.Score, group = interaction(Temperature, DPE))) +
  
  geom_ribbon(data = modelbased::estimate_expectation(alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]][["Simpson"]], data = "grid"),
              aes(x = DPE, y = Predicted, color = Temperature, fill = Temperature, 
                  ymin = CI_low, ymax = CI_high), 
              alpha = 0.2,
              inherit.aes = F) + 
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_violin(aes(color = Temperature,
    fill = Temperature),
    # color = "white",
    draw_quantiles = c(.5),
    size = 1.5,
    alpha = .75,
    scale = "width",
    show.legend = F,
    na.rm = T) +
  
  scale_color_manual(values = c("white","white","white")) +
  scale_fill_manual(values = c("white","white","white")) +
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_violin(aes(color = Temperature,
                  fill = Temperature),
              size = .75,
              alpha = 0.05,
              scale = "width",
              show.legend = F,
              na.rm = T) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_violin(aes(color = Temperature,
                  fill = Temperature),
              draw_quantiles = c(.5),
              size = .75,
              alpha = 0,
              scale = "width",
              show.legend = F,
              na.rm = T) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  ggbeeswarm::geom_quasirandom(aes(#color = Temperature,
    color = Temperature),
    size = 1.5,
    fill = "white",
    shape = 23,
    stroke = .75,
    show.legend = T,
    dodge.width = 6,
    varwidth = T,
    na.rm = T) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_smooth(data = modelbased::estimate_expectation(alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]][["Simpson"]], data = "grid"),
              aes(x = DPE, y = Predicted, color = Temperature, fill = Temperature),
              size = 3,
              method = "glm", se = F,
              inherit.aes = F) +
  
  scale_color_manual(values = c("white","white","white")) +
  scale_fill_manual(values = c("white","white","white")) +
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_smooth(data = modelbased::estimate_expectation(alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]][["Simpson"]], data = "grid"),
              aes(x = DPE, y = Predicted, color = Temperature),
              linetype = "dashed",
              size = 1,
              method = "glm", se = F,
              inherit.aes = F) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
  
  # Scaling
  scale_x_continuous(limits = c(-5, 47), breaks = seq(0, 42, by = 7)) +
  scale_y_continuous(limits = c(-0.05, 1.15), breaks = seq(0, 1, by = .25)) + # Commented, since no sig comparisons, so no sig bars
  
  # Labels
  labs(#title = "Alpha Diversity Scores across Temperature and Time", #caption = "Alpha; All",
    x = "Days post exposure (dpe)",
    y = "Alpha Score (Normalized)"
  ) +
  
  # Facet
  facet_grid(.~ Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
  
  # Theme Settings
  theme(legend.position = "none",
        legend.justification = "center",
        legend.direction = "horizontal",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 12)
  ) #+
  # 
  # # Statistical Notation
  # - No Significant results, so no notation present
  # ggpubr::stat_pvalue_manual(alpha.plots[[tmp.resSubSection]][["TEMP:DPE"]][["Tukey"]][["Sig.Labels"]] %>%
  #                              dplyr::mutate(group1 = as.numeric(group1),
  #                                            group2 = as.numeric(group2)),
  #                            label = "p.adj.sig",
  #                            y.position = NULL, #c(1, 1.075, 1.15, 1, 1, 1.15, 1.075),
  #                            step.group.by = "term",
  #                            size = 6,
  #                            label.size = 7,
  #                            bracket.size = 1,
  #                            hide.ns = T)



## BETA --------------------------------------------------------------------


### TEMP --------------------------------------------------------------------

#### Plot --------------------------------------------------------------------

beta.plots[[tmp.resSubSection]][["TEMP"]][["CAP"]][["Plot"]] <-
  tmp.psOBJ %>%
  # ps_filter(Temperature == temp) %>%
  # ps_filter(DPE == 0) %>%
  ps_mutate(Cluster = if_else(
    Treatment == "Exposed" & Total.Worm.Count > 0,
    case_when(
      Simpson__Genus_norm %in% head(sort(Simpson__Genus_norm), 16) ~ "Low",
      Simpson__Genus_norm %in% tail(sort(Simpson__Genus_norm), 16) ~ "High",
      TRUE ~ "Other"
    ),
    "Other"
  )) %>%
  # ungroup()  %>%
  ps_mutate(Cluster = fct_relevel(factor(Cluster, levels = c("Other", "Low", "High")))) %>%
  # samdat_tbl()
  tax_agg("Genus") %>%
  dist_calc("bray") %>%
  dist_permanova(
    seed = 1,
    variables = c("Temperature"),
    n_processes = 8,
    n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
  ) %>% 
  ord_calc(method = "CAP", constraints = "Temperature.") %>%
  # ord_plot(color = "Temperature", #shape = "Treatment",
  #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
  #          size = 3) +
  microViz::ord_plot(
    color = "Temperature",
    fill = "white",
    size = 3,
    stroke = 1,
    shape = 23,
    constraint_vec_style = vec_constraint(colour = "black",
                                          size = 2,
                                          alpha = 1,
                                          arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
    constraint_lab_style = constraint_lab_style(
      type = "label",
      justify = "side",
      colour = "black",
      # max_angle = 90, 
      # perpendicular = TRUE, 
      size = 3.5,
      check_overlap = TRUE
    ), 
    auto_caption = NA
  ) +
  
  stat_ellipse(aes(color = Temperature, 
                   fill = Temperature),
               alpha = .05,
               linetype = "dashed",
               geom = "polygon"
  ) +
  
  # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
  
  scale_shape_manual(values = c(16, 23)) +
  scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
  scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C"))  +
  
  scale_x_continuous(limits = c(-2,2.5)) +
  
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        strip.text = element_text(size = 14))

# Reshuffle layers so ellipse is in background
beta.plots[[tmp.resSubSection]][["TEMP"]][["CAP"]][["Plot"]] <- 
  rearrange_layers(beta.plots[[tmp.resSubSection]][["TEMP"]][["CAP"]][["Plot"]])


### TEMP:DPE --------------------------------------------------------------------

#### Plot --------------------------------------------------------------------

beta.plots[[tmp.resSubSection]][["TEMP:DPE"]][["CAP"]][["Plot"]] <-
  tmp.psOBJ %>%
  # ps_filter(Temperature == temp) %>%
  # ps_filter(DPE == 0) %>%
  ps_mutate(Cluster = if_else(
    Treatment == "Exposed" & Total.Worm.Count > 0,
    case_when(
      Simpson__Genus_norm %in% head(sort(Simpson__Genus_norm), 16) ~ "Low",
      Simpson__Genus_norm %in% tail(sort(Simpson__Genus_norm), 16) ~ "High",
      TRUE ~ "Other"
    ),
    "Other"
  )) %>%
  # ungroup()  %>%
  ps_mutate(Cluster = fct_relevel(factor(Cluster, levels = c("Other", "Low", "High")))) %>%
  # samdat_tbl()
  tax_agg("Genus") %>%
  dist_calc("canberra") %>%
  dist_permanova(
    seed = 1,
    variables = c("Temperature"),
    n_processes = 8,
    n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
  ) %>% 
  ord_calc(constraints = c("Temperature.", "DPE", "Temp.DPE."),
           method = "CAP") %>%
  # ord_plot(color = "Temperature", #shape = "Treatment",
  #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
  #          size = 3) +
  microViz::ord_plot(
    axes = c(1, 2),
    color = "Temperature",
    fill = "white",
    size = 3,
    stroke = 1,
    shape = 23,
    constraint_vec_style = vec_constraint(colour = "black",
                                          size = 2,
                                          alpha = 1,
                                          arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
    constraint_lab_style = constraint_lab_style(
      type = "label",
      justify = "side",
      colour = "black",
      # max_angle = 90, 
      # perpendicular = TRUE, 
      size = 3.5,
      check_overlap = TRUE
    ), 
    auto_caption = NA
  ) +
  
  stat_ellipse(aes(color = Temperature, 
                   fill = Temperature),
               alpha = .05,
               linetype = "dashed",
               geom = "polygon"
  ) +
  
  # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
  
  scale_shape_manual(values = c(16, 23)) +
  scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
  scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C"))  +
  
  # scale_x_continuous(limits = c(-3.5,2.5)) +
  
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        strip.text = element_text(size = 14))


# Reshuffle layers so ellipse is in background
beta.plots[[tmp.resSubSection]][["TEMP:DPE"]][["CAP"]][["Plot"]] <- 
  rearrange_layers(beta.plots[[tmp.resSubSection]][["TEMP:DPE"]][["CAP"]][["Plot"]])


## SUPP --------------------------------------------------------------------


#### S4A: TEMP --------------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP"]][["TUKEY"]][["Plot_SUPP"]] <-
  {
    
    
    
    tmp.sig.bars <-
      alpha.stats[[tmp.resSubSection]][["TEMP"]][["Tukey"]] %>%
      SigStars(pval.var = "adj.p.value") %>%
      dplyr::filter(p.value < 0.05 ) 
    
    tmp.psOBJ %>%
      microViz::samdat_tbl() %>%
      tidyr::pivot_longer(cols = contains("_norm"), names_to = "Alpha.Metric", values_to = "Alpha.Score") %>%
      cutCellNames(col = "Alpha.Metric") %>%
      # filter(Alpha.Metric == "Simpson") %>%
      # filter(Treatment == "Control") %>%
      # dplyr::mutate(DPE = as.numeric(levels(DPE)[DPE])) %>%
      ggplot(aes(x = Temperature, y = Alpha.Score)) +
      
      
      geom_violin(aes(#color = Temperature,
        fill = Temperature,
        group = interaction(Temperature)),
        # position = position_dodge(width = 5),
        # width = 4,
        color = "white",
        draw_quantiles = c(.5),
        size = 3,
        alpha = .05,
        scale = "width",
        show.legend = F,
        na.rm = T) +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature,
                      group = interaction(Temperature)),
                  # position = position_dodge(width = 5),
                  # width = 4,
                  size = 1.25,
                  alpha = 0,
                  scale = "width",
                  show.legend = F,
                  na.rm = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature,
                      group = interaction(Temperature)),
                  # position = position_dodge(width = 5),
                  # width = 4,
                  draw_quantiles = c(.5),
                  size = 1.25,
                  alpha = 0,
                  scale = "width",
                  show.legend = F,
                  na.rm = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      ggbeeswarm::geom_quasirandom(aes(
        color = Temperature),
        size = 3,
        fill = "white",
        shape = 23,
        stroke = 2,
        show.legend = T,
        # position = position_dodge2(width = 5),
        na.rm = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
      
      scale_x_discrete(
        labels = function(x) paste0(x, "°C"), # Append '°C' to x-axis labels
        breaks = c('28', '32', '35')
      ) +
      scale_y_continuous(limits = c(-0.05, 1.2), breaks = seq(0, 1, by = .25)) +
      labs(title = "Alpha Diversity Scores across Temperature", #caption = "Alpha; All",
           x = "Temperature",
           y = "Alpha Score (Normalized)"
      ) +
      guides(
        color = guide_legend(order = 1, title = "Temp (°C)",
                             override.aes = list(shape = 16, size = 5, linetype = NA))
      ) +
      
      facet_grid(Alpha.Metric~.) +
      
      theme(legend.position = "none",
            legend.justification = "center",
            legend.direction = "vertical",
            panel.grid.major.x = element_blank(),
            strip.text.y = element_text(size = 12)) +  # Adjust the size for the facet labels on the y-axis (if applicable)) +
      
      # Statistical Notation
      ggpubr::stat_pvalue_manual(tmp.sig.bars,
                                 label = "p.adj.sig",
                                 y.position = c(1, 1.1, 1, 1.1, 1),
                                 size = 3,
                                 label.size = 6,
                                 bracket.size = 1,
                                 hide.ns = T)   
  }


### S4B: TEMP --------------------------------------------------------------------

beta.plots[[tmp.resSubSection]][["TEMP"]][["CAP"]][["Plot_SUPP"]] <-
  {
    # Canberra
    tmp.plot.canberra <-
      tmp.psOBJ %>%
      tax_agg("Genus") %>%
      dist_calc("canberra") %>%
      dist_permanova(
        seed = 1,
        variables = c("Temperature"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(method = "CAP", constraints = "Temperature.") %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_plot(
        color = "Temperature",
        fill = "white",
        size = 3,
        stroke = 1,
        shape = 23,
        constraint_vec_style = vec_constraint(colour = "black",
                                              size = 2,
                                              alpha = 1,
                                              arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
        constraint_lab_style = constraint_lab_style(
          type = "label",
          justify = "side",
          colour = "black",
          # max_angle = 90, 
          # perpendicular = TRUE, 
          size = 3.5,
          check_overlap = TRUE
        ), 
        auto_caption = NA
      ) +
      
      stat_ellipse(aes(color = Temperature, 
                       fill = Temperature),
                   alpha = .05,
                   linetype = "dashed",
                   geom = "polygon"
      ) +
      
      # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      scale_shape_manual(values = c(16, 23)) +
      scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
      scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C")) +
      
      # scale_x_continuous(limits = c(-2,1.5)) +
      labs(caption = "canberra") +
      
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.canberra <- 
      rearrange_layers(tmp.plot.canberra)
    
    # Unifrac
    tmp.plot.gunifrac <-
      tmp.psOBJ %>%
      # tax_agg("Genus") %>%
      dist_calc("gunifrac") %>%
      dist_permanova(
        seed = 1,
        variables = c("Temperature"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(method = "CAP", constraints = "Temperature.") %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_plot(
        color = "Temperature",
        fill = "white",
        size = 3,
        stroke = 1,
        shape = 23,
        constraint_vec_style = vec_constraint(colour = "black",
                                              size = 2,
                                              alpha = 1,
                                              arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
        constraint_lab_style = constraint_lab_style(
          type = "label",
          justify = "side",
          colour = "black",
          # max_angle = 90, 
          # perpendicular = TRUE, 
          size = 3.5,
          check_overlap = TRUE
        ), 
        auto_caption = NA
      ) +
      
      stat_ellipse(aes(color = Temperature, 
                       fill = Temperature),
                   alpha = .05,
                   linetype = "dashed",
                   geom = "polygon"
      ) +
      
      # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      scale_shape_manual(values = c(16, 23)) +
      scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
      scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C"))  +
      
      # scale_x_continuous(limits = c(-6,8)) +
      labs(caption = "gunifrac") +
      
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.gunifrac <- 
      rearrange_layers(tmp.plot.gunifrac)
    
    # Combine plots
    cowplot::plot_grid(tmp.plot.canberra, tmp.plot.gunifrac, ncol = 2)
  }


#### S4C: TEMP:DPE --------------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP:DPE"]][["TUKEY"]][["Plot_SUPP"]] <-
  {
    
    tmp.sig.bars <-
      alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["Tukey"]] %>%
      SigStars(pval.var = "adj.p.value") %>%
      dplyr::filter(p.value < 0.05 ) 
    
    tmp.psOBJ %>%
      microViz::samdat_tbl() %>%
      tidyr::pivot_longer(cols = contains("_norm"), names_to = "Alpha.Metric", values_to = "Alpha.Score") %>%
      cutCellNames(col = "Alpha.Metric") %>%
      # filter(Alpha.Metric == "Simpson") %>%
      # filter(Treatment == "Control") %>%
      # dplyr::mutate(DPE = as.numeric(levels(DPE)[DPE])) %>%
      ggplot(aes(x = DPE, y = Alpha.Score, group = interaction(Temperature, DPE))) +
      
      # geom_ribbon(data = modelbased::estimate_expectation(alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]][["Simpson"]], data = "grid"),
      #             aes(x = DPE, y = Predicted, color = Temperature, fill = Temperature, 
      #                 ymin = CI_low, ymax = CI_high), 
      #             alpha = 0.2,
      #             inherit.aes = F) + 
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature),
                  # color = "white",
                  draw_quantiles = c(.5),
                  size = 1.5,
                  alpha = .75,
                  scale = "width",
                  show.legend = F,
                  na.rm = T) +
      
      scale_color_manual(values = c("white","white","white")) +
      scale_fill_manual(values = c("white","white","white")) +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature),
                  size = .75,
                  alpha = 0.05,
                  scale = "width",
                  show.legend = F,
                  na.rm = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature),
                  draw_quantiles = c(.5),
                  size = .75,
                  alpha = 0,
                  scale = "width",
                  show.legend = F,
                  na.rm = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      ggbeeswarm::geom_quasirandom(aes(#color = Temperature,
        color = Temperature),
        size = 1.5,
        fill = "white",
        shape = 23,
        stroke = .75,
        show.legend = T,
        dodge.width = 6,
        varwidth = T,
        na.rm = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_smooth(#data = modelbased::estimate_expectation(alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]][["Simpson"]], data = "grid"),
        aes(x = DPE, y = Alpha.Score, color = Temperature, fill = Temperature, group = Temperature),
        size = 3,
        method = "glm", se = F,
        inherit.aes = T) +
      
      scale_color_manual(values = c("white","white","white")) +
      scale_fill_manual(values = c("white","white","white")) +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_smooth(#data = modelbased::estimate_expectation(alpha.stats[[tmp.resSubSection]][["TEMP:DPE"]][["GLM"]][["Simpson"]], data = "grid"),
        aes(x = DPE, y = Alpha.Score, color = Temperature, group = Temperature),
        size = 1.5,
        method = "glm", se = F,
        inherit.aes = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none") +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none") +
      
      # Scaling
      scale_x_continuous(limits = c(-5, 47), breaks = seq(0, 42, by = 7)) +
      scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, by = .25)) + # Commented, since no sig comparisons, so no sig bars
      
      # Labels
      labs(#title = "Alpha Diversity Scores across Temperature and Time", #caption = "Alpha; All",
        x = "Days post exposure (dpe)",
        y = "Alpha Score (Normalized)"
      ) +
      
      # Facet
      facet_grid(Alpha.Metric ~ Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      # Theme Settings
      theme(legend.position = "none",
            legend.justification = "center",
            legend.direction = "horizontal",
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(size = 12),
            strip.text.y = element_text(size = 12)
      ) #+
    
    # Statistical Notation
    # ggpubr::stat_pvalue_manual(tmp.sig.bars %>%
    #                              dplyr::mutate(group1 = as.numeric(group1),
    #                                            group2 = as.numeric(group2)),
    #                            label = "p.adj.sig",
    #                            y.position = c(1,1.1),
    #                            step.group.by = "term",
    #                            size = 3,
    #                            label.size = 5,
    #                            bracket.size = .75,
    #                            hide.ns = T)   
  }


### S4D: TEMP:DPE --------------------------------------------------------------------

beta.plots[[tmp.resSubSection]][["TEMP:DPE"]][["CAP"]][["Plot_SUPP"]] <-
  {
    # Canberra
    tmp.plot.bray <-
      tmp.psOBJ %>%
      tax_agg("Genus") %>%
      dist_calc("bray") %>%
      dist_permanova(
        seed = 1,
        variables = c("Temperature", "DPE", "Temp.DPE"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(method = "CAP", constraints = c("Temperature.", "DPE", "Temp.DPE.")) %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      ord_calc(constraints = c("Temperature.", "DPE", "Temp.DPE."),
               method = "CAP") %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_plot(
        axes = c(1, 2),
        color = "Temperature",
        fill = "white",
        size = 3,
        stroke = 1,
        shape = 23,
        constraint_vec_style = vec_constraint(colour = "black",
                                              size = 2,
                                              alpha = 1,
                                              arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
        constraint_lab_style = constraint_lab_style(
          type = "label",
          justify = "side",
          colour = "black",
          # max_angle = 90, 
          # perpendicular = TRUE, 
          size = 3.5,
          check_overlap = TRUE
        ), 
        auto_caption = NA
      ) +
      
      stat_ellipse(aes(color = Temperature, 
                       fill = Temperature),
                   alpha = .05,
                   linetype = "dashed",
                   geom = "polygon"
      ) +
      
      # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      scale_shape_manual(values = c(16, 23)) +
      scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
      scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C"))   +
      
      # scale_x_continuous(limits = c(-3.5,3)) +
      labs(caption = "bray-curtis") +
      
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.bray <- 
      rearrange_layers(tmp.plot.bray)
    
    # Unifrac
    tmp.plot.unifrac <-
      tmp.psOBJ %>%
      # tax_agg("Genus") %>%
      dist_calc("gunifrac") %>%
      dist_permanova(
        seed = 1,
        variables = c("Temperature", "DPE", "Temp.DPE"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(method = "CAP", constraints = c("Temperature.", "DPE", "Temp.DPE.")) %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_plot(
        axes = c(1, 2),
        color = "Temperature",
        fill = "white",
        size = 3,
        stroke = 1,
        shape = 23,
        constraint_vec_style = vec_constraint(colour = "black",
                                              size = 2,
                                              alpha = 1,
                                              arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
        constraint_lab_style = constraint_lab_style(
          type = "label",
          justify = "side",
          colour = "black",
          # max_angle = 90, 
          # perpendicular = TRUE, 
          size = 3.5,
          check_overlap = TRUE
        ), 
        auto_caption = NA
      ) +
      
      stat_ellipse(aes(color = Temperature, 
                       fill = Temperature),
                   alpha = .05,
                   linetype = "dashed",
                   geom = "polygon"
      ) +
      
      # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      scale_shape_manual(values = c(16, 23)) +
      scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
      scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C")) +
      
      scale_x_continuous(limits = c(-3,3)) +
      labs(caption = "gunifrac") +
      
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.unifrac <- 
      rearrange_layers(tmp.plot.unifrac)
    
    # Combine plots
    cowplot::plot_grid(tmp.plot.bray, tmp.plot.unifrac, ncol = 2)
  }


### S4D.2: TEMP:DPE (TIME) --------------------------------------------------------------------

beta.plots[[tmp.resSubSection]][["TEMP:DPE__TIME"]][["CAP"]][["Plot_SUPP"]] <-
  {
    
    
    # Bray
    tmp.plot.bray <-
      tmp.psOBJ %>%
      microViz::ps_mutate(
        DPE_labels = factor(
          DPE,                          # or use your existing DPE_labels column
          levels = c(14, 21, 28, 42),   # order you want in the legend
          labels = c("14", "21", "28", "42")
        )) %>%      
      microViz::tax_agg("Genus") %>%
      microViz::dist_calc("bray") %>%
      microViz::dist_permanova(
        seed = 1,
        variables = c("Temperature", "DPE", "Temp.DPE"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      microViz::ord_calc(method = "CAP", constraints = c("Temperature.", "DPE", "Temp.DPE.")) %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_calc(constraints = c("Temperature.", "DPE", "Temp.DPE."),
                         method = "CAP") %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_plot(
        axes = c(1, 2),
        fill = "DPE_labels",
        color = "white",
        size = 3,
        stroke = 1,
        shape = 23,
        constraint_vec_style = vec_constraint(colour = "black",
                                              size = 2,
                                              alpha = 1,
                                              arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
        constraint_lab_style = constraint_lab_style(
          type = "label",
          justify = "side",
          colour = "black",
          # max_angle = 90, 
          # perpendicular = TRUE, 
          size = 3.5,
          check_overlap = TRUE
        ), 
        auto_caption = NA
      ) +
      
      scale_fill_viridis_d(name = "Days Post Exposure (DPE)", 
                           option = "viridis", 
                           labels = c("14", "21", "28", "42")) + 
      
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale_color() +
      
      stat_ellipse(aes(color = Temperature, 
                       fill = Temperature),
                   alpha = .05,
                   linetype = "dashed",
                   geom = "polygon"
      ) +
      
      # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      scale_shape_manual(values = c(16, 23)) +
      scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
      scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C"))   +
      
      # scale_x_continuous(limits = c(-3.5,3)) +
      labs(caption = "bray-curtis") +
      
      theme(legend.position = "bottom",
            legend.direction = "vertical",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.bray <- 
      rearrange_layers(tmp.plot.bray)
    
    # Canberra
    
    tmp.plot.canberra <- 
      
      tmp.psOBJ %>%
      # ps_filter(Temperature == temp) %>%
      # ps_filter(DPE == 0) %>%
      microViz::ps_mutate(Cluster = if_else(
        Treatment == "Exposed" & Total.Worm.Count > 0,
        case_when(
          Simpson__Genus_norm %in% head(sort(Simpson__Genus_norm), 16) ~ "Low",
          Simpson__Genus_norm %in% tail(sort(Simpson__Genus_norm), 16) ~ "High",
          TRUE ~ "Other"
        ),
        "Other"
      )) %>%
      # ungroup()  %>%
      microViz::ps_mutate(Cluster = fct_relevel(factor(Cluster, levels = c("Other", "Low", "High")))) %>%
      # samdat_tbl()
      microViz::ps_mutate(
        DPE_labels = factor(
          DPE,                          # or use your existing DPE_labels column
          levels = c(14, 21, 28, 42),   # order you want in the legend
          labels = c("14", "21", "28", "42")
        )) %>%  
      microViz::tax_agg("Genus") %>%
      microViz::dist_calc("canberra") %>%
      microViz::dist_permanova(
        seed = 1,
        variables = c("Temperature", "DPE", "Temp.DPE"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      microViz::ord_calc(constraints = c("Temperature.", "DPE", "Temp.DPE."),
                         method = "CAP") %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_plot(
        axes = c(1, 2),
        fill = "DPE_labels",
        color = "white",
        size = 3,
        stroke = 1,
        shape = 23,
        constraint_vec_style = vec_constraint(colour = "black",
                                              size = 2,
                                              alpha = 1,
                                              arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
        constraint_lab_style = constraint_lab_style(
          type = "label",
          justify = "side",
          colour = "black",
          # max_angle = 90, 
          # perpendicular = TRUE, 
          size = 3.5,
          check_overlap = TRUE
        ), 
        auto_caption = NA
      ) +
      
      scale_fill_viridis_d(name = "Days Post Exposure (DPE)", 
                           option = "viridis", 
                           labels = c("14", "21", "28", "42")) + 
      
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale_color() +
      
      stat_ellipse(aes(color = Temperature, 
                       fill = Temperature),
                   alpha = .05,
                   linetype = "dashed",
                   geom = "polygon"
      ) +
      
      # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      scale_shape_manual(values = c(16, 23)) +
      scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
      scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C"))  +
      
      # scale_x_continuous(limits = c(-3.5,2.5)) +
      
      labs(caption = "canberra") +
      
      theme(legend.position = "bottom",
            legend.direction = "vertical",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.canberra <- 
      rearrange_layers(tmp.plot.canberra)
    
    # Unifrac
    tmp.plot.unifrac <-
      tmp.psOBJ %>%
      # tax_agg("Genus") %>%
      microViz::ps_mutate(
        DPE_labels = factor(
          DPE,                          # or use your existing DPE_labels column
          levels = c(14, 21, 28, 42),   # order you want in the legend
          labels = c("14", "21", "28", "42")
        )) %>%  
      microViz::dist_calc("gunifrac") %>%
      microViz::dist_permanova(
        seed = 1,
        variables = c("Temperature", "DPE", "Temp.DPE"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      microViz::ord_calc(method = "CAP", constraints = c("Temperature.", "DPE", "Temp.DPE.")) %>%
      # ord_plot(color = "Temperature", #shape = "Treatment",
      #          fill = "Temperature", #ifelse("Treatment" == "Control", "Temperature", "white"),
      #          size = 3) +
      microViz::ord_plot(
        axes = c(1, 2),
        fill = "DPE_labels",
        color = "white",
        size = 3,
        stroke = 1,
        shape = 23,
        constraint_vec_style = vec_constraint(colour = "black",
                                              size = 2,
                                              alpha = 1,
                                              arrow = grid::arrow(length = grid::unit(0.05, units = "npc"))),
        constraint_lab_style = constraint_lab_style(
          type = "label",
          justify = "side",
          colour = "black",
          # max_angle = 90, 
          # perpendicular = TRUE, 
          size = 3.5,
          check_overlap = TRUE
        ), 
        auto_caption = NA
      ) +
      
      scale_fill_viridis_d(name = "Days Post Exposure (DPE)", 
                           option = "viridis", 
                           labels = c("14", "21", "28", "42")) + 
      
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale_color() +
      
      stat_ellipse(aes(color = Temperature, 
                       fill = Temperature),
                   alpha = .05,
                   linetype = "dashed",
                   geom = "polygon"
      ) +
      
      # facet_grid(.~Temperature, labeller = labeller(Temperature = c("28" = "28°C", "32" = "32°C", "35" = "35°C"))) +
      
      scale_shape_manual(values = c(16, 23)) +
      scale_color_manual(values = col.Temp, labels = c("28°C", "32°C", "35°C")) +
      scale_fill_manual(values = c(col.Temp, "white"), labels = c("28°C", "32°C", "35°C")) +
      
      scale_x_continuous(limits = c(-3,3)) +
      labs(caption = "gunifrac") +
      
      theme(legend.position = "bottom",
            legend.direction = "vertical",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.unifrac <- 
      rearrange_layers(tmp.plot.unifrac)
    
    # Combine plots
    cowplot::plot_grid(tmp.plot.bray, tmp.plot.canberra,  tmp.plot.unifrac, ncol = 3)
  }


