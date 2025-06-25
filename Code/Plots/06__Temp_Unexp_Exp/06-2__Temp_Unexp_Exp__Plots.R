# 06-2__Temp_Unexp_Exp__Plots -----------------------------------------------------

cat("6-2__Temp_Unexp_Exp__Plots: Post-exposed \n")

## Post-Exposed

tmp.psOBJ <- ps.list[["PostExposed"]]
tmp.resSubSection <- "PostExposed"


## TEMP:TREAT (Post-Exposed) ----------------------------------------------------------------

### Alpha ---------------------------------------------------------------

#### Sig Labels --------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["Tukey"]][["Sig.Labels"]] <- 
  alpha.stats[[tmp.resSubSection]][["TEMP:TREAT"]][["Tukey"]] %>%
  SigStars(pval.var = "adj.p.value") %>%
  dplyr::filter(p.value < 0.05 ) %>%
  dplyr::filter(Alpha.Metric == "Simpson") 


#### Plot --------------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["TUKEY"]][["Plot"]] <-
  tmp.psOBJ %>%
  microViz::samdat_tbl() %>%
  tidyr::pivot_longer(cols = contains("_norm"), names_to = "Alpha.Metric", values_to = "Alpha.Score") %>%
  cutCellNames(col = "Alpha.Metric") %>%
  filter(Alpha.Metric == "Simpson") %>%
  # filter(Treatment == "Control") %>%
  # dplyr::mutate(DPE = as.numeric(levels(DPE)[DPE])) %>%
  ggplot(aes(x = Treatment, y = Alpha.Score, group = interaction(Temperature, Treatment))) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_violin(aes(color = Temperature,
                  fill = Temperature,
                  shape = Treatment,
                  alpha = ifelse(Treatment == "Control", Treatment, "white")),
              size = 1.25,
              # alpha = .05,
              scale = "width",
              show.legend = F,
              na.rm = T) +
  
  scale_alpha_manual(values = c(.5, .05)) +
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none", label = c("Control", "Post-Exposure")) +
  scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none", label = c("Control", "Post-Exposure")) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_violin(aes(color = Temperature,
                  fill = Temperature),
              # color = "white",
              draw_quantiles = c(.5),
              size = 3,
              alpha = 0,
              scale = "width",
              show.legend = F,
              na.rm = T) +
  
  scale_color_manual(values = c("white","white","white")) +
  scale_fill_manual(values = c("white","white","white")) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_violin(aes(color = Temperature,
                  fill = Temperature),
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
  
  ggbeeswarm::geom_quasirandom(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                                   fill = ifelse(Treatment == "Control", Temperature, "white"),
                                   shape = Treatment),
                               size = 2,
                               # fill = "white",
                               # shape = 23,
                               stroke = 1.25,
                               show.legend = T,
                               # dodge.width = 6,
                               # varwidth = T,
                               na.rm = T) +
  
  scale_fill_manual(values = c(col.Temp, "white"), name = "Treatment", guide = "none", label = c("Control", "Post-Exposure")) +
  scale_color_manual(values = c(col.Temp, "white"), name = "Temperature", guide = "none", label = c("Control", "Post-Exposure")) +
  scale_shape_manual(values = c(21, 23)) +
  
  # Scaling
  scale_x_discrete(labels = c("Control" = "Control", "Exposed" = "Exposed")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = .25)) +
  
  # Labels
  labs(#title = "Alpha Diversity Scores across Temperature and Time", #caption = "Alpha; All",
    x = "Treatment",
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
  )  +

# Statistical Notation
  ggpubr::stat_pvalue_manual(alpha.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["Tukey"]][["Sig.Labels"]],
                             label = "p.adj.sig",
                             y.position = c(1.075), #c(1, 1.075, 1.15, 1, 1, 1.15, 1.075),
                             step.group.by = "term",
                             size = 6,
                             label.size = 7,
                             bracket.size = 1,
                             hide.ns = T)




### Beta --------------------------------------------------------------------

#### Plot --------------------------------------------------------------------


beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot"]] <-
  tmp.psOBJ %>%
  tax_agg("Genus") %>%
  dist_calc("bray") %>%
  dist_permanova(
    seed = 1,
    variables = c("Temperature", "Treatment", "Temp.Treat"),
    n_processes = 8,
    n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
  ) %>% 
  ord_calc(
    constraints = c("Temperature.", "Treatment.", "Temp.Treat."),
    method = "CAP"
  ) %>%
  microViz::ord_plot(
    axes = c(1, 2),
    color =  NA,
    shape =  NA,
    fill = NA,
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
      check_overlap = F
    ), 
    auto_caption = NA,
    show.legend = F
  ) +
  
  geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                 fill = ifelse(Treatment == "Control", Temperature, "white"),
                 shape = Treatment),
             size = 3,
             stroke = 1,
             show.legend = T) +
  
  scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Exposed")) +
  scale_fill_manual(values = c(col.Temp, "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
  scale_color_manual(values = c(col.Temp, "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale("alpha") +
  
  stat_ellipse(aes(color = Temperature,
                   fill = Temperature,
                   linetype = Treatment,
                   alpha = ifelse(Treatment == "Control", "Temperature", "white")),
               # alpha = .05,
               # linetype = "dashed",
               size = .5,
               geom = "polygon",
               show.legend = F
  ) +
  
  scale_fill_manual(values = c(col.Temp, "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
  scale_color_manual(values = c(col.Temp, "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
  scale_alpha_manual(values = c(.5, .15), guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  
  scale_x_continuous(limits = c(-1.1, 3.1)) +
  
  guides(
    shape = guide_legend(order = , override.aes = list(shape = c(16, 23), size = 3))
  ) +
  
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top",
        strip.text = element_text(size = 14)) 

# Reshuffle layers so ellipse is in background and labels are on top
beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot"]] <- 
  rearrange_layers(beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot"]])

beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot"]] <- 
  move_label_layer_to_top(beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot"]])


### SUPP --------------------------------------------------------------------

#### S6C: TEMP:TREAT (Post-Exp) ---------------------------------------------------------------------

alpha.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["TUKEY"]][["Plot_SUPP"]] <-
  {
    
    tmp.sig.bars <-
      alpha.stats[[tmp.resSubSection]][["TEMP:TREAT"]][["Tukey"]] %>%
      SigStars(pval.var = "adj.p.value") %>%
      dplyr::filter(p.value < 0.05 ) 
    
    tmp.psOBJ %>%
      microViz::samdat_tbl() %>%
      tidyr::pivot_longer(cols = contains("_norm"), names_to = "Alpha.Metric", values_to = "Alpha.Score") %>%
      cutCellNames(col = "Alpha.Metric") %>%
      # filter(Alpha.Metric == "Simpson") %>%
      # filter(Treatment == "Control") %>%
      # dplyr::mutate(DPE = as.numeric(levels(DPE)[DPE])) %>%
      ggplot(aes(x = Treatment, y = Alpha.Score, group = interaction(Temperature, Treatment))) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature,
                      shape = Treatment,
                      alpha = ifelse(Treatment == "Control", Treatment, "white")),
                  size = 1.25,
                  # alpha = .05,
                  scale = "width",
                  show.legend = F,
                  na.rm = T) +
      
      scale_alpha_manual(values = c(.5, .05)) +
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temp (°C)", guide = "none", label = c("Control", "Post-Exposure")) +
      scale_color_manual(values = col.Temp, name = "Temp (°C)", guide = "none", label = c("Control", "Post-Exposure")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature),
                  # color = "white",
                  draw_quantiles = c(.5),
                  size = 3,
                  alpha = 0,
                  scale = "width",
                  show.legend = F,
                  na.rm = T) +
      
      scale_color_manual(values = c("white","white","white")) +
      scale_fill_manual(values = c("white","white","white")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      
      geom_violin(aes(color = Temperature,
                      fill = Temperature),
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
      
      ggbeeswarm::geom_quasirandom(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                                       fill = ifelse(Treatment == "Control", Temperature, "white"),
                                       shape = Treatment),
                                   size = 2,
                                   # fill = "white",
                                   # shape = 23,
                                   stroke = 1.25,
                                   show.legend = T,
                                   # dodge.width = 6,
                                   # varwidth = T,
                                   na.rm = T) +
      
      scale_fill_manual(values = c(col.Temp, "white"), name = "Treatment", guide = "none", label = c("Control", "Post-Exposure")) +
      scale_color_manual(values = c(col.Temp, "white"), name = "Temperature", guide = "none", label = c("Control", "Post-Exposure")) +
      scale_shape_manual(values = c(21, 23)) +
      
      # Scaling
      # scale_x_discrete(labels = c("Control" = "Control", "Exposed" = "Exposed")) +
      scale_y_continuous(limits = c(-0.05, 1.2), breaks = seq(0, 1, by = .25)) +
      
      # Labels
      labs(#title = "Alpha Diversity Scores across Temperature and Time", #caption = "Alpha; All",
        x = "Treatment",
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
      ) +
      
      # Statistical Notation
      ggpubr::stat_pvalue_manual(tmp.sig.bars,
                                 label = "p.adj.sig",
                                 y.position = c(1,1.1,1,1.1,1,1.1,1),
                                 step.group.by = "term",
                                 size = 3,
                                 label.size = 5,
                                 bracket.size = .75,
                                 hide.ns = T)
  }

#### S6D: TEMP:TREAT (Post-Exp) ---------------------------------------------------------------------

beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot_SUPP"]] <-
  {
    # Canberra
    tmp.plot.canberra <-
      tmp.psOBJ %>%
      tax_agg("Genus") %>%
      dist_calc("canberra") %>%
      dist_permanova(
        seed = 1,
        variables = c("Temperature", "Treatment", "Temp.Treat"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(method = "CAP", constraints = c("Temperature.", "Treatment.", "Temp.Treat.")) %>%
      
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp, "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp, "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp, "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      # scale_x_continuous(limits = c(-3.5,3)) +
      labs(caption = "canberra") +
      
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.canberra <- 
      rearrange_layers(tmp.plot.canberra)
    
    tmp.plot.canberra <- 
      move_label_layer_to_top(tmp.plot.canberra)
    
    # Unifrac
    tmp.plot.unifrac <-
      tmp.psOBJ %>%
      # tax_agg("Genus") %>%
      dist_calc("gunifrac") %>%
      dist_permanova(
        seed = 1,
        variables = c("Temperature", "Treatment", "Temp.Treat"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(method = "CAP", constraints = c("Temperature.", "Treatment.", "Temp.Treat.")) %>%
      
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp, "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp, "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp, "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp, "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      # scale_x_continuous(limits = c(-3,3)) +
      labs(caption = "gunifrac") +
      
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            strip.text = element_text(size = 14))
    
    # Reshuffle layers so ellipse is in background
    tmp.plot.unifrac <- 
      rearrange_layers(tmp.plot.unifrac)
    
    tmp.plot.unifrac <- 
      move_label_layer_to_top(tmp.plot.unifrac)
    
    # Combine plots
    cowplot::plot_grid(tmp.plot.canberra, tmp.plot.unifrac, ncol = 2)
  }

##### S6D.1 TREAT by TEMP -------------------------------------------------------------------


###### Bray -------------------------------------------------------------------


beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot_SUPP_6D.1_bray"]] <-
  {
    
    p1 <- tmp.psOBJ %>%
      ps_filter(Temperature == 28) %>%
      tax_agg("Genus") %>%
      dist_calc("bray") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[1], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[1], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[1], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[1], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "28°C (Bray-Curtis)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    p2 <- tmp.psOBJ %>%
      ps_filter(Temperature == 32) %>%
      tax_agg("Genus") %>%
      dist_calc("bray") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[2], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[2], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[2], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[2], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "32°C (Bray-Curtis)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    p3 <- tmp.psOBJ %>%
      ps_filter(Temperature == 35) %>%
      tax_agg("Genus") %>%
      dist_calc("bray") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[3], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[3], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[3], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[3], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "35°C (Bray-Curtis)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    cowplot::plot_grid(p1, p2, p3, ncol = 3)
  }

###### Canberra -------------------------------------------------------------------


beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot_SUPP_6D.1_canberra"]] <-
  {
    
    p1 <- tmp.psOBJ %>%
      ps_filter(Temperature == 28) %>%
      tax_agg("Genus") %>%
      dist_calc("canberra") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[1], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[1], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[1], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[1], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "28°C (Canberra)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    p2 <- tmp.psOBJ %>%
      ps_filter(Temperature == 32) %>%
      tax_agg("Genus") %>%
      dist_calc("canberra") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[2], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[2], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[2], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[2], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "32°C (Canberra)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    p3 <- tmp.psOBJ %>%
      ps_filter(Temperature == 35) %>%
      tax_agg("Genus") %>%
      dist_calc("canberra") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[3], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[3], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[3], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[3], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "35°C (Canberra)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    cowplot::plot_grid(p1, p2, p3, ncol = 3)
  }

###### Gunifrac -------------------------------------------------------------------


beta.plots[[tmp.resSubSection]][["TEMP:TREAT"]][["CAP"]][["Plot_SUPP_6D.1_gunifrac"]] <-
  {
    
    p1 <- tmp.psOBJ %>%
      ps_filter(Temperature == 28) %>%
      # tax_agg("Genus") %>%
      dist_calc("gunifrac") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[1], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[1], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[1], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[1], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "28°C (Gunifrac)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    p2 <- tmp.psOBJ %>%
      ps_filter(Temperature == 32) %>%
      # tax_agg("Genus") %>%
      dist_calc("gunifrac") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[2], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[2], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[2], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[2], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "32°C (Gunifrac)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    p3 <- tmp.psOBJ %>%
      ps_filter(Temperature == 35) %>%
      # tax_agg("Genus") %>%
      dist_calc("gunifrac") %>%
      dist_permanova(
        seed = 1,
        variables = c("Treatment"),
        n_processes = 8,
        n_perms = 999 # only 99 perms used in examples for speed (use 9999+!)
      ) %>% 
      ord_calc(
        constraints = c("Treatment."),
        method = "CAP"
      ) %>%
      microViz::ord_plot(
        axes = c(1, 2),
        color =  NA,
        shape =  NA,
        fill = NA,
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
          check_overlap = F
        ), 
        auto_caption = NA,
        show.legend = F
      ) +
      
      geom_point(aes(color = ifelse(Treatment != "Control", Temperature, "white"),
                     fill = ifelse(Treatment == "Control", Temperature, "white"),
                     shape = Treatment),
                 size = 3,
                 stroke = 1,
                 show.legend = T) +
      
      scale_shape_manual(values = c(21, 23), name = "Treatment", labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[3], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[3], "white"), name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale("alpha") +
      
      stat_ellipse(aes(color = Temperature,
                       fill = Temperature,
                       linetype = Treatment,
                       alpha = ifelse(Treatment == "Control", "Temperature", "white")),
                   # alpha = .05,
                   # linetype = "dashed",
                   size = .5,
                   geom = "polygon",
                   show.legend = F
      ) +
      
      scale_alpha_manual(values = c(.5, .15), guide = "none", labels = c("Control", "Pre-Exposed")) +
      scale_linetype_manual(values = c("solid", "dashed"), labels = c("Control", "Pre-Exposed")) +
      scale_fill_manual(values = c(col.Temp[3], "white"), guide = "none" ) + #name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      scale_color_manual(values = c(col.Temp[3], "white"), guide = "none" ) + #, name = "Temperature", labels = c("28°C", "32°C", "35°C", "")) +
      # scale_x_continuous(limits = c(-1.1, 3.1)) +
      
      guides(
        shape = guide_legend(override.aes = list(shape = c(16, 23), size = 3))
      ) +
      
      labs(title = "35°C (Gunifrac)") +
      
      theme(legend.position = "none",
            legend.direction = "horizontal",
            legend.title.position = "top",
            strip.text = element_text(size = 14)) 
    
    cowplot::plot_grid(p1, p2, p3, ncol = 3)
  }




