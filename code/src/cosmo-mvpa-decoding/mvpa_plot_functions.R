#plotting functions for MVPAs 

unimodal_mvpa_plot_ReadSpeech <- function(plot_data, summary_data){
    
    #Extract mean and se accuracy from summary stats
    mean_accu <- round(summary_data$mean_accuracy, digits = 2)
    se_accu <- round(summary_data$se_accuracy, digits = 2)
    
    plot_data %>% 
        ggplot(aes(x = conditions, y = accuracy, color = roiArea)) +
        geom_point(aes(alpha = I(0.7)), 
                   position = position_jitterdodge(jitter.width=0.2,
                                                   dodge.width = 0.3),
                   size = 4.5) +
        stat_summary(aes(x =conditions, y = accuracy, color = roiArea),
                     fun = "mean", 
                     geom = "crossbar", 
                     position = position_dodge(width = 0.5),
                     width = .55,
                     colour  = 'black',
                     inherit.aes = FALSE,
                     show.legend = FALSE) + 
        stat_summary(aes(x = conditions, y = accuracy, color = roiArea),
                     fun.max = function(x) mean(x) + (std.error(x)), 
                     fun.min = function(x) mean(x) - (std.error(x)), 
                     geom = "errorbar",
                     position = position_dodge(width = 0.5),
                     width = .15,
                     size = 0.8,
                     colour  = 'black',
                     inherit.aes = FALSE,
                     show.legend = FALSE) + 
        geom_hline(yintercept=c(50), 
                   linetype="dotted", 
                   colour="black", 
                   linewidth=.5) +
        #scale_color_manual(values = c(rep(dark_colors[1], 3))) +
        scale_y_continuous(name = "Clasisfier Accuracy [%]", 
                           limits = c(0,100)) +
        scale_x_discrete(name = 'Decoding conditions') +
        annotate("text", 
                 x = classifier, 
                 y = 5, 
                 label = mean_accu)+
        annotate("text", 
                 x = classifier, 
                 y = 0, 
                 label = se_accu) +
        theme_cowplot(font_size = 16, font_family = "Arial") +
        theme(axis.line = element_line(colour = 'black', size = 1),
              axis.ticks = element_line(colour = 'black', size = 1),
              axis.text = element_text(face="bold"),
              legend.position = 'none')
    
}

crossmodal_mvpa_plot_ReadSpeech <- function(plot_data, summary_data){
    
    #Extract mean and se accuracy from summary stats
    mean_accu <- round(summary_data$mean_accuracy, digits = 2)
    se_accu <- round(summary_data$se_accuracy, digits = 2)
    
    plot_data %>% 
        ggplot(aes(x = conditions, 
                   y = accuracy, 
                   color = TrainTest, 
                   alpha = TrainTest,
                   size = TrainTest,
                   group = interaction(conditions, TrainTest))) +
        geom_point(position = position_jitterdodge(jitter.width=0.2,
                                                   dodge.width = 0.9)) +
        
        #scale_alpha_manual(values = c(0.2,0.8, 0.2), 
        scale_alpha_manual(values = c(0.95,0.8, 0.95), 
                           name = "Decoding partitions") +
        # scale_color_manual(values = c(rep(dark_colors[1], 3)),
        #scale_color_manual(values = blind_triple_colors,
        #                   name = "Decoding partitions")+
        scale_size_manual(values = c(1.75,4,1.75),
                          name = "Decoding partitions")+
        scale_y_continuous(name = "Clasisfier Accuracy [%]", 
                           limits = c(0,100)) +
        scale_x_discrete(name = 'Decoding conditions') +
        ## ADD STAT SUMMARIES 
        stat_summary(fun = "mean",
                     geom = "crossbar",
                     position = position_dodge(width = 0.9),
                     width = .75,
                     size = 0.75,
                     #color = 'black',
                     alpha = 1,
                     show.legend = FALSE) +
        stat_summary(fun.max = function(x) mean(x) + (sd(x)/sqrt(20)),
                     fun.min = function(x) mean(x) - (sd(x)/sqrt(20)),
                     geom = "errorbar",
                     position = position_dodge(width = 0.9),
                     width = .15,
                     size = 0.8,
                     alpha = 1,
                     color = "black",
                     show.legend = FALSE) +
        #ADD CHANCE LINE
        geom_hline(yintercept=c(50), 
                   linetype="dotted", 
                   colour="black", 
                   linewidth=.5) +
        # ## TRY TO ANNOTATE THE PLOTS WITH MEAN AND SEM VALUES - only average!? 
        annotate("text",
                 x = 1:3,
                 y = 5,
                 label = mean_accu) +
        annotate("text",
                 x = 1:3,
                 y = 0,
                 label = se_accu) +
        theme_cowplot(font_size = 16, font_family = "Arial") +
        theme(axis.line = element_line(colour = 'black', size = 1),
              axis.ticks = element_line(colour = 'black', size = 1),
              axis.text = element_text(face="bold"),
              legend.position="bottom",
              legend.direction = "vertical",
              legend.justification = 0.60)
    
}


#plotting functions for univariate ROIs

univariate_fMRI_ROI_points <- function(data, roi, sub_group) {
    
    #Prepare triple colour schemes for groups
    blind_triple_colors <- c(brewer.pal(3, "Pastel2")[1], 
                             brewer.pal(3, "Set2")[1],   
                             brewer.pal(3, "Dark2")[1],  
                             brewer.pal(3, "Pastel2")[1], 
                             brewer.pal(3, "Set2")[1],   
                             brewer.pal(3, "Dark2")[1])
    
    sighted_triple_colors <- c(brewer.pal(3, "Pastel2")[2], 
                               brewer.pal(3, "Set2")[2],   
                               brewer.pal(3, "Dark2")[2],  
                               brewer.pal(3, "Pastel2")[2], 
                               brewer.pal(3, "Set2")[2],   
                               brewer.pal(3, "Dark2")[2])
    
    # Y Axis label for BOLD 
    bold_label <- "BOLD contrast estimate (a.u.)"
    
    # Prepare legend labels? 
    legend_labels <- if(sub_group == "Blind"){
        c("Control (Blind)", "Pseudowords (Blind)", "Words (Blind)")
    } else if (sub_group == "Sighted"){
        c("Control (Sighted)", "Pseudowords (Sighted)", "Words (Sighted)")
    }
    
    #Main plot
    data %>% 
        filter(Group == sub_group) %>% 
        filter(ROI == roi) %>% 
        mutate(Condition = factor(Condition, levels=c("Control", "Pseudowords", "Words"))) %>%
        ggplot(aes(x = Modality, y = Contrast_Estimate, group = interaction(Modality, Condition), 
                   colour = Condition)) +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_point(size = 3,
                   alpha = 0.8,
                   position = position_jitterdodge(jitter.width=0.2,
                                                   dodge.width = 1)) +
        stat_summary(fun = "mean",
                     geom = "crossbar",
                     position = position_dodge(width = 1),
                     width = .75,
                     size = 1,
                     show.legend = FALSE) +
        stat_summary(fun.max = function(x) mean(x) + (sd(x)/sqrt(20)),
                     fun.min = function(x) mean(x) - (sd(x)/sqrt(20)),
                     geom = "errorbar",
                     position = position_dodge(width = 1),
                     width = .15,
                     size = 0.8,
                     color = "black",
                     show.legend = FALSE) +
        scale_color_manual(values = if(sub_group == "Blind"){
            blind_triple_colors
        } else if (sub_group == "Sighted"){
            sighted_triple_colors
        },
        labels = legend_labels) +
        scale_y_continuous(limits = c(-15, 20), 
                           name = bold_label) +
        scale_x_discrete(name = "Modality") +
        theme_cowplot(font_size = 20, font_family = "Arial") +
        if(sub_group == "Blind") {
            theme(axis.line = element_line(colour = 'black', size = 1),
                  axis.ticks = element_line(colour = 'black', size = 1),
                  axis.text = element_text(face="bold"),
                  #    legend.position = "none",
                  legend.title = element_blank())}
    else if(sub_group == "Sighted") {
        theme(axis.line = element_line(colour = 'black', size = 1),
              axis.ticks = element_line(colour = 'black', size = 1),
              axis.text = element_text(face="bold"),
              axis.line.y = element_line(colour = 'black', size = 0),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title.y = element_blank(),
              # legend.position = "none",
              legend.title = element_blank())
    }
    
}

univariate_fMRI_ROI_bars <- function(data, roi, sub_group) {
    
    #Prepare triple colour schemes for groups
    blind_triple_colors <- c(brewer.pal(3, "Pastel2")[1], 
                             brewer.pal(3, "Set2")[1],   
                             brewer.pal(3, "Dark2")[1],  
                             brewer.pal(3, "Pastel2")[1], 
                             brewer.pal(3, "Set2")[1],   
                             brewer.pal(3, "Dark2")[1])
    
    sighted_triple_colors <- c(brewer.pal(3, "Pastel2")[2], 
                               brewer.pal(3, "Set2")[2],   
                               brewer.pal(3, "Dark2")[2],  
                               brewer.pal(3, "Pastel2")[2], 
                               brewer.pal(3, "Set2")[2],   
                               brewer.pal(3, "Dark2")[2])
    
    # Y Axis label for BOLD 
    bold_label <- "BOLD contrast estimate (a.u.)"
    
    # Prepare legend labels? 
    legend_labels <- if(sub_group == "Blind"){
        c("Control (Blind)", "Pseudowords (Blind)", "Words (Blind)")
    } else if (sub_group == "Sighted"){
        c("Control (Sighted)", "Pseudowords (Sighted)", "Words (Sighted)")
    }
    
    #Main plot
    data %>% 
        ## SUBSET THE DATA
        filter(Group == sub_group) %>% 
        filter(ROI_label == roi) %>% 
        mutate(Condition = factor(Condition, levels=c("Control", "Pseudowords", "Words"))) %>%
        ## ADD SUMMARIES - not needed for barplot?
        group_by(Modality, Condition) %>% 
        summarise(mean = mean(ContrastEstimate), sem = sd(ContrastEstimate)/sqrt(n())) %>% 
        ungroup() %>% 
        ## PLOT THAT PLOT
        ggplot(aes(x = Modality, y = mean, group = interaction(Modality, Condition), 
                   color = Condition)) +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_bar(stat = "summary",
                 fun = "mean",
                 position= "dodge",
                 color = "black",
                 size = 1,
                 aes(fill = Condition)) +
        geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),
                      position=position_dodge(0.9), 
                      width = 0.2, 
                      size = 1,
                      colour = 'black') +
        scale_fill_manual(values = if(sub_group == "Blind"){
            blind_triple_colors
        } else if (sub_group == "Sighted"){
            sighted_triple_colors
        },
        labels = legend_labels) +
        scale_y_continuous(name = bold_label
                           #limits = c(-15, 20), 
        ) +
        scale_x_discrete(name = "Modality") +
        #theme_cowplot(font_size = 20, font_family = "Arial") +
    #     if(sub_group == "Blind") {
    #         theme(axis.line = element_line(colour = 'black', size = 1),
    #               axis.ticks = element_line(colour = 'black', size = 1),
    #               axis.text = element_text(face="bold"),
    #               #    legend.position = "none",
    #               legend.title = element_blank(),
    #               legend.position = 'bottom')}
    # else if(sub_group == "Sighted") {
    #     theme(axis.line = element_line(colour = 'black', size = 1),
    #           axis.ticks = element_line(colour = 'black', size = 1),
    #           axis.text = element_text(face="bold"),
    #           axis.line.y = element_line(colour = 'black', size = 0),
    #           axis.text.y=element_blank(),
    #           axis.ticks.y=element_blank(),
    #           axis.title.y = element_blank(),
    #           legend.position = 'bottom',
    #           legend.title = element_blank())
    # }
        theme_cowplot(font_family = "Arial") +
        theme(axis.line = element_line(colour = 'black', size = 1),
                            axis.ticks = element_line(colour = 'black', size = 1),
                            axis.text = element_text(face="bold"),
                            #    legend.position = "none",
                            legend.title = element_blank(),
                            legend.position = 'bottom')
}


#NOW THESE ARE PRETTY COOL! THEY MAKE SUBPLOTS FOR BOTH GROUP, PASTE THEM 
#TOGETHER AND EQUALIZE THE Y AXES OF BOTH SUB-PLOTS TO BE IDENTICAL! 
# YOU CAN CHOOSE BARPLOTS OR POINTS + MEAN/SD! 
univariate_fMRI_ROI_bars_both_gr <- function(data, roi) {
    
    #Prepare triple colour schemes for groups
    blind_triple_colors <- c(rep(c(brewer.pal(3, "Pastel2")[1], 
                                   brewer.pal(3, "Set2")[1],   
                                   brewer.pal(3, "Dark2")[1]),2))
    
    sighted_triple_colors <- c(rep(c(brewer.pal(3, "Pastel2")[2], 
                                     brewer.pal(3, "Set2")[2],   
                                     brewer.pal(3, "Dark2")[2]),2))
    
    # Y Axis label for BOLD 
    bold_label <- "BOLD contrast estimate (a.u.)"
    legend_labels <- c('Control', 'Pseudowords', 'Words')
    
    groups <- unique(data$Group)
    
    for(iGr in seq_along(groups)){
        assign(paste0('p',iGr), 
               #Main plot
               data %>% 
                   ## SUBSET THE DATA
                   filter(Group == groups[iGr]) %>% 
                   filter(ROI_label == roi) %>% 
                   mutate(Condition = factor(Condition, levels=c("Control", "Pseudowords", "Words"))) %>%
                   ## ADD SUMMARIES - not needed for barplot?
                   group_by(Modality, Condition) %>% 
                   summarise(mean = mean(ContrastEstimate), 
                             sem = sd(ContrastEstimate)/sqrt(n()),
                             .groups = 'keep') %>% 
                   ungroup() %>% 
                   ## PLOT THAT PLOT
                   ggplot(aes(x = Modality, y = mean, group = interaction(Modality, Condition), 
                              color = Condition)) +
                   geom_hline(yintercept = 0, linetype = "dotted") +
                   geom_bar(stat = "summary",
                            fun = "mean",
                            position= "dodge",
                            color = "black",
                            size = 1,
                            aes(fill = Condition)) +
                   geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),
                                 position=position_dodge(0.9), 
                                 width = 0.2, 
                                 size = 1,
                                 colour = 'black') +
                   scale_fill_manual(values = if(iGr == 1){
                       blind_triple_colors
                   }else{
                       sighted_triple_colors
                   },labels = legend_labels) +
                   ggtitle(label = groups[iGr]) +
                   scale_x_discrete(name = "Modality") +
                   theme_cowplot(font_family = "Arial") +
                   theme(axis.line = element_line(colour = 'black', size = 1),
                         axis.ticks = element_line(colour = 'black', size = 1),
                         axis.text = element_text(face="bold"),
                         #    legend.position = "none",
                         legend.title = element_blank(),
                         legend.position = 'bottom',
                         #Center the plot title?
                         plot.title = element_text(hjust = 0.5)))
    }
        
        #Hard-coded patchworking???
        p_combined <- p1+p2

        #Get y axis limits
        p_ranges_y <- c(ggplot_build(p_combined[[1]])$layout$panel_scales_y[[1]]$range$range,
                        ggplot_build(p_combined[[2]])$layout$panel_scales_y[[1]]$range$range)

        #Update Y axes
        p_updated <- p_combined & scale_y_continuous(name = bold_label,
                                       limits = c(min(p_ranges_y), max(p_ranges_y)))
            
        return(p_updated)
    
}
univariate_fMRI_ROI_points_both_gr <- function(data, roi) {
    
    #Prepare triple colour schemes for groups
    blind_triple_colors <- c(rep(c(brewer.pal(3, "Pastel2")[1], 
                                   brewer.pal(3, "Set2")[1],   
                                   brewer.pal(3, "Dark2")[1]),2))
    
    sighted_triple_colors <- c(rep(c(brewer.pal(3, "Pastel2")[2], 
                                     brewer.pal(3, "Set2")[2],   
                                     brewer.pal(3, "Dark2")[2]),2))
    
    # Axis and legend labels
    bold_label <- "BOLD contrast estimate (a.u.)"
    legend_labels <- c('Control', 'Pseudowords', 'Words')
    
    groups <- unique(data$Group)
    
    for(iGr in seq_along(groups)){
        assign(paste0('p',iGr), 
               #Main plot
               data %>% 
                   ## SUBSET THE DATA
                   filter(Group == groups[iGr]) %>% 
                   filter(ROI_label == roi) %>% 
                   mutate(Condition = factor(Condition, levels=c("Control", "Pseudowords", "Words"))) %>%
                   ## PLOT THAT PLOT (summaries not needed, they are added as stat_summary)
                   ggplot(aes(x = Modality, y = ContrastEstimate, group = interaction(Modality, Condition), 
                              color = Condition)) +
                   geom_hline(yintercept = 0, linetype = "dotted") +
                   geom_point(size = 2,
                              alpha = 0.8,
                              position = position_jitterdodge(jitter.width=0.2,
                                                              dodge.width = 1)) +
                   stat_summary(fun = "mean",
                                geom = "crossbar",
                                position = position_dodge(width = 1),
                                width = 0.95,
                                size = 0.8,
                                show.legend = FALSE) +
                   stat_summary(fun.max = function(x) mean(x) + (sd(x)/sqrt(20)),
                                fun.min = function(x) mean(x) - (sd(x)/sqrt(20)),
                                geom = "errorbar",
                                position = position_dodge(width = 1),
                                width = .15,
                                size = 0.8,
                                color = "black",
                                show.legend = FALSE) +
                   scale_color_manual(values = if(iGr == 1){
                       blind_triple_colors}else{
                           sighted_triple_colors},
                       labels = legend_labels) +
                   ggtitle(label = groups[iGr]) +
                   scale_x_discrete(name = "Modality") +
                   theme_cowplot(font_family = "Arial") +
                   theme(axis.line = element_line(colour = 'black', size = 1),
                         axis.ticks = element_line(colour = 'black', size = 1),
                         axis.text = element_text(face="bold"),
                         #    legend.position = "none",
                         legend.title = element_blank(),
                         legend.position = 'bottom',
                         #Center the plot title?
                         plot.title = element_text(hjust = 0.5)))
    }
    
    #Hard-coded patchworking???
    p_combined <- p1+p2
    
    #Get y axis limits
    p_ranges_y <- c(ggplot_build(p_combined[[1]])$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(p_combined[[2]])$layout$panel_scales_y[[1]]$range$range)
    
    #Update Y axes
    p_updated <- p_combined & scale_y_continuous(name = bold_label,
                                                 limits = c(min(p_ranges_y), max(p_ranges_y)))
    
    return(p_updated)
    
}