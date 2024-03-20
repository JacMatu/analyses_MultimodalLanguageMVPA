

library('tidyverse')
library('data.table')

roi <- 'leftMTa' # leftMTp rightMTp

ffx <- 'beta' # beta tmap

featRatio <- '200'

space <-'MNI'

FWHM <- 0

namePattern <- paste(featRatio, '.*\\.csv$', sep = '')


# set te path to the task folder
cosmo_results_path <- '/Volumes/MICK/analisys_high-re_multiSensSpatFreq/outputs/derivatives/cosmo-mvpa/crossmodal'



# read all the files 
allfiles <- list.files(cosmo_results_path, pattern = namePattern)


mvpa <- lapply(paste(cosmo_results_path, allfiles, sep ='/'), read.csv, sep =',') %>% 
  rbindlist

# check that the headers are ok
head(mvpa)

# mvpa <- mvpa[,1:8]

# mvpa$modality <- ifelse(mvpa$modality == 'visual', 'visual_highVSlow', 
#                               ifelse(mvpa$modality == 'auditory', 'auditory_highVSlow', 
#                                      ifelse(mvpa$modality == 'fullMotionAudVis', 'fullMotion_AudVSVis', 
#                                             ifelse(mvpa$modality == 'tactile', 'tactile_highVSlow', 'error'))))
# 
mvpa$conditions_order <- ifelse(mvpa$conditions == 'trainVis_testAud', 1,
                                ifelse(mvpa$conditions == 'trainAud_testVis', 2,
                                       ifelse(mvpa$conditions == 'both', 3, 4)))

mvpa_accuracy <-  mvpa %>%
  group_by(roiArea,
           ffxResults,
           conditions,
           conditions_order) %>%
  summarize(mean_accuracy = mean(accuracy),
            sd_accurarcy = sd(accuracy),
            se_accuracy = sd(accuracy)/sqrt(15),
            .groups = 'keep') 

ggplot() +
  geom_point(data = subset(mvpa_accuracy, 
                           ffxResults == ffx & roiArea == roi),
             aes(x = reorder(conditions, conditions_order),
                 y = mean_accuracy,
                 color = conditions),
             shape = 5,
             size = 3,
             stroke = 2) +
  geom_errorbar(data = subset(mvpa_accuracy, 
                              ffxResults == ffx & roiArea == roi),
                aes(x = conditions,
                    y = mean_accuracy,
                    ymin = mean_accuracy - se_accuracy, 
                    ymax = mean_accuracy + se_accuracy,
                    color = conditions), 
                width = .15,
                position = position_dodge(1),
                size = 1,
                alpha = .8) +
  geom_point(data = subset(mvpa, 
                            ffxResults == ffx & roiArea == roi),
             aes(x = reorder(conditions, conditions_order),
                 y = accuracy,
                 color = conditions),
             position = position_jitter(w = 0.3, h = 0.01),
             alpha = 0.2,
             size = 3.5,
             stroke = 0,
             show.legend = F) +
  geom_hline(yintercept=c(.50), 
             linetype="dotted", 
             colour="red", 
             linewidth=.5) +
  theme_classic() +
  scale_y_continuous(limits = c(0.24, 0.86),
                     breaks = seq(0.3, 0.8, by = .1)) +
  scale_color_manual(values = c('#af5d1e',"#359c84", "#feb300"),
                     breaks = c('trainVis_testAud', 'trainAud_testVis', 'both'),
                     labels = c('trainVis_testAud', 'trainAud_testVis', 'both')) +
  xlab(roi) +
  ylab('classification acc.') +
  ggtitle(ffx) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    text=element_text(size=20),
    legend.position = 'right',
    legend.text=element_text(size=14),
    axis.line = element_line(size = 0.6),
    axis.text.y = element_text(size=14, colour='black'))

filename <- paste(cosmo_results_path, '/', 
                  'crossmod_space-', space, 
                  '_FWHM-', FWHM, 
                  '_ffx-', ffx, 
                  '_featRation-', featRatio, 
                  '_roi-', roi, '_', 
                  format(Sys.time(), "%Y%m%d%H%M"), '.png', sep = '')


ggsave(filename, device="png", units="in", width=6.54, height=4.54, dpi=300)   

  