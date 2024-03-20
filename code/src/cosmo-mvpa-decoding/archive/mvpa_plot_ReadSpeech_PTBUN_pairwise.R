#This stuff needs to be optimized in the following ways: 
# load and differentiate between groups
# loop across ROIs 
# keep plot as a function!

# perfect the plots (colours, labels, titels)

library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)




roi <- 'perVWFA' # lexVWFA, perVWFA

ffx <- 'beta' # beta tmap

featRatio <- '0.8'

space <-'IXI549Space'

FWHM <- 2

group <- 'blind' # blind or sighted

namePattern <- paste('sub-',group,'.*_label-VWFAwithinmodal_.*\\.csv$', sep = '')

# set te path to the task folder
cosmo_results <- '/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa'

#task_audVis <-  paste('task-audVisMotSpatialFreq_space-MNI152NLin2009cAsym_FWHM-', as.character(FWHM), '_node-mvpaBlockAverage', sep ='')
task_ReadSpeech <- paste0('task-MultimodalReadSpeech_space-IXI549Space_FWHM-',as.character(FWHM),'_node-mvpa6betas')

cosmo_results_path_ReadSpeech <- paste(cosmo_results, task_ReadSpeech, 'accuracy', sep ='/')

# read all the files 
allfiles_ReadSpeech <- list.files(cosmo_results_path_ReadSpeech, pattern = namePattern)

# allfiles_tact <- list.files(cosmo_results_path_tact, pattern = ".csv")


ReadSpeech <- lapply(paste(cosmo_results_path_ReadSpeech, allfiles_ReadSpeech, sep ='/'), read.csv, sep =',') %>% 
  rbindlist


mvpa <- ReadSpeech

# check that the headers are ok
head(mvpa)

# mvpa <- mvpa[,1:8]

mvpa$modality <- ifelse(mvpa$modality == 'reading', 'reading_WordsVPseudoVControl', 
                              ifelse(mvpa$modality == 'speech', 'speech_WordsVPseudoVControl', 'error'))


mvpa$modality_order <- ifelse(mvpa$modality == 'reading_WordsVPseudoVControl', 1, 
                        ifelse(mvpa$modality == 'speech_WordsVPseudoVControl', 2, 3))

#List all decoding types
#mvpa_conditions <- unique(mvpa$conditions)
#mvpa_rois <- unique(mvpa$roiArea)

mvpa_accuracy <-  mvpa %>%
  group_by(roiArea,
           ffxResults,
           modality,
          modality_order) %>%
          #conditions) %>%
  summarize(mean_accuracy = mean(accuracy),
            sd_accurarcy = sd(accuracy),
            se_accuracy = sd(accuracy)/sqrt(20),
            .groups = 'keep') 


## TRY THE MAIN LOOP ACROSS  ROIS and CONDITIONS
## ADD GROUP IN THE FUTURE

#for (r in seq_along(mvpa_rois)) {
#    for(c in seq_along(mvpa_conditions)) {
        
     ggplot() +
            geom_errorbar(data = subset(mvpa_accuracy, 
                                        ffxResults == ffx & roiArea == roi & modality_order !=4),
                          aes(x = reorder(modality, modality_order),
                              y = mean_accuracy,
                              ymin = mean_accuracy - se_accuracy, 
                              ymax = mean_accuracy + se_accuracy,
                              color = modality), 
                          width = .15,
                          position = position_dodge(1),
                          size = 1,
                          alpha = .8) +
            geom_point(data = subset(mvpa_accuracy, 
                                     ffxResults == ffx & roiArea == roi & modality_order !=4),
                       aes(x = reorder(modality, modality_order),
                           y = mean_accuracy,
                           color = modality),
                       fill = "white",
                       shape = 23,
                       size = 3,
                       stroke = 1.5) +
            geom_point(data = subset(mvpa, 
                                     ffxResults == ffx & roiArea == roi & modality_order !=4),
                       aes(x = reorder(modality, modality_order),
                           y = accuracy,
                           color = modality),
                       position = position_jitter(w = 0.3, h = 0.01),
                       alpha = 0.2,
                       size = 3.5,
                       stroke = 0,
                       show.legend = F) +
            geom_hline(yintercept=c(.33), 
                       linetype="dotted", 
                       colour="red", 
                       linewidth=.5) +
            theme_classic() +
            scale_y_continuous(limits = c(0, 1.01),
                               breaks = seq(0, 1, by = 0.2)) +
            scale_color_manual(values = c('#e35656',"#0091d4"), 
                               breaks = c('reading_WordsVPseudoVControl', 'speech_WordsVPseudoVControl'),
                               labels = c('reading_WordsVPseudoVControl', 'speech_WordsVPseudoVControl')) +
            xlab(roi) +
            ylab('classification acc.') +
            ggtitle(paste(group, ffx)) +
            theme(
                axis.text.x = element_blank(),
                axis.ticks = element_blank(),
                text=element_text(size=20),
                legend.position = 'right',
                legend.text=element_text(size=14),
                axis.line = element_line(size = 0.6),
                axis.text.y = element_text(size=14, colour='black'))
        
        #rename the plot 
#        assign(paste('plot',group,mvpa_rois[r],mvpa_conditions[c], sep = '_'), plot)
#        rm(plot)
        
        
 #   }
#}


    
    
    
    
    
#filename <- paste(cosmo_results, '/', 
#                  'unimod_space-', space, 
#                  '_FWHM-', FWHM, 
#                  '_ffx-', ffx, 
#                  '_featRation-', featRatio, 
#                  '_roi-', roi, '_', 
#                  format(Sys.time(), "%Y%m%d%H%M"), '.png', sep = '')

#ggsave(filename, device="png", units="in", width=6.54, height=4.54, dpi=300)   

  