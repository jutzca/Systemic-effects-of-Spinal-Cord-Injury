# ------------------------------------------------------------------------------------------------------------------
# Plots - blood markers in SCI
#
# November 2020
# L. Bourguignon
# ------------------------------------------------------------------------------------------------------------------

# Load packages ----
require(ggplot2)
require(data.table)
require(plyr)
require(dplyr)

# Load data ----

setwd('/Volumes/borgwardt/Projects/SCI_Blood_Biomarker/figures_labrotation_Lucie/code/App-1/')

# Age dataframe
df_age_combined <- read.csv('data/df_age_combined.csv')

# Bloodvalue dataframes - Sygen
CBC_sygen <- read.csv('data/df_final_CBC_sygen.csv')
kidney_sygen <- read.csv('data/df_final_kidney_sygen.csv')
liver_sygen <- read.csv('data/df_final_liver_sygen.csv')
muscle_sygen <- read.csv('data/df_final_muscle_sygen.csv')
pancreas_sygen <- read.csv('data/df_final_pancreas_sygen.csv')
rest_sygen <- read.csv('data/df_final_rest_sygen.csv')

# Bloodvalue dataframes - Murnau
CBC_murnau_long <- read.csv('data/df_final_CBC_long.csv')
kidney_murnau_long <- read.csv('data/df_final_kidney_long.csv')
liver_murnau_long <- read.csv('data/df_final_liver_long.csv')
muscle_murnau_long <- read.csv('data/df_final_muscle_long.csv')
pancreas_murnau_long <- read.csv('data/df_final_pancreas_long.csv')
rest_murnau_long <- read.csv('data/df_final_rest_long.csv')

# Bloodvalue dataframes - Combined
hb <- read.csv('data/df_combined/df_combined_hb.csv')
mcv <- read.csv('data/df_combined/df_combined_mcv.csv')

# Functions ----

make_plot_murnau <- function(mydata, list_norm_up, list_norm_down, times, grades, markers, levels_T6, levels_3cat, labels) {
  # Filter based on user input
  data_transformed <- mydata[mydata$bloodmarker %in% unlist(markers, use.names=FALSE), ] # filter blood markers
  data_transformed <- data_transformed[data_transformed$time %in% unlist(times, use.names=FALSE), ] # filter time points
  data_transformed <- data_transformed[data_transformed$AIS %in% unlist(grades, use.names=FALSE), ] # filter AIS grades
  data_transformed <- data_transformed[data_transformed$level_T6 %in% unlist(levels_T6, use.names=FALSE), ] # filter above/below T6
  data_transformed <- data_transformed[data_transformed$level_3cat %in% unlist(levels_3cat, use.names=FALSE), ] # filter cervical/thoracic/lumbar injuries

  # Make time column as factor variable
  data_transformed$time_f = factor(data_transformed$time, levels=unlist(times, use.names=FALSE)) # make time a factor
  data_transformed$AIS_f = factor(data_transformed$AIS, levels=unlist(grades, use.names=FALSE)) # make AIS grades factors
  data_transformed$bloodmarker_f = factor(data_transformed$bloodmarker, levels=unlist(markers, use.names=FALSE)) # make blood markers factors
  
  # Remove NAs = if no bloodvalue at a certain time point for a certain patient
  data_transformed = data_transformed[complete.cases(data_transformed),] 
  
  # Create grid and select the normal levels for the markers selected by user
  norm_lvl <- expand.grid(bloodmarker = unique(data_transformed$bloodmarker_f), time = unique(data_transformed$time_f))
  norm_lvl$up <- list_norm_up
  norm_lvl$down <- list_norm_down

  p1 <- data_transformed %>% 
    dplyr::group_by(time_f, bloodmarker) %>%  # do the same calcs for each box
    dplyr::mutate(value2 = filter_lims(bloodvalue)) %>%  # new variable (value2) so as not to displace first one)
    ggplot(aes(x=AIS, y=value2)) + 
    geom_boxplot(na.rm = TRUE, coef = 5) +  # remove NAs, and set the whisker length to all included points
    geom_point(fill=alpha("black", 0.01), alpha = 0.2, shape=21, colour="grey20", position=position_jitter(width=0.2)) +
    facet_grid(bloodmarker ~ time_f, scales="free", labeller = labeller(bloodmarker=labels)) + 
    theme_bw() +
    geom_hline(data = norm_lvl, aes(yintercept = up), linetype="dashed") + # horizontal line for upper limit of the norm values
    geom_hline(data = norm_lvl, aes(yintercept = down), linetype="dashed") + # horizontal line for lower limit of the norm values
    geom_rect(data = norm_lvl, aes(xmin=0, xmax=Inf, ymin=down, ymax=up, fill="Normal range"), inherit.aes = FALSE, alpha = 0.01) +
    scale_fill_manual('', values = 'blue', guide = guide_legend(override.aes = list(alpha = 0.07))) +
    theme(legend.key = element_rect(color="black", linetype="dashed")) +
    scale_x_discrete(labels=c("AIS A" = "A", "AIS B" = "B", "AIS C" = "C", 'AIS D' = 'D')) + 
    ylab("Bloodvalue")
  
  # Output plot
  return(p1)
}

make_plot_murnau_line <- function(mydata, list_norm_up, list_norm_down, times, grades, markers, levels_T6, levels_3cat, labels) {
  # Filter based on user input
  data_transformed <- mydata[mydata$bloodmarker %in% unlist(markers, use.names=FALSE), ] # filter blood markers
  data_transformed <- data_transformed[data_transformed$time %in% unlist(times, use.names=FALSE), ] # filter time points
  data_transformed <- data_transformed[data_transformed$AIS %in% unlist(grades, use.names=FALSE), ] # filter AIS grades
  data_transformed <- data_transformed[data_transformed$level_T6 %in% unlist(levels_T6, use.names=FALSE), ] # filter above/below T6
  data_transformed <- data_transformed[data_transformed$level_3cat %in% unlist(levels_3cat, use.names=FALSE), ] # filter cervical/thoracic/lumbar injuries
  
  # Make time column as factor variable
  data_transformed$time_f = factor(data_transformed$time, levels=unlist(times, use.names=FALSE)) # make time a factor
  data_transformed$AIS_f = factor(data_transformed$AIS, levels=unlist(grades, use.names=FALSE)) # make AIS grades factors
  data_transformed$bloodmarker_f = factor(data_transformed$bloodmarker, levels=unlist(markers, use.names=FALSE)) # make blood markers factors
  
  # Remove NAs = if no bloodvalue at a certain time point for a certain patient
  data_transformed = data_transformed[complete.cases(data_transformed),] 
  
  # Create grid and select the normal levels for the markers selected by user
  norm_lvl <- expand.grid(bloodmarker = unique(data_transformed$bloodmarker_f), time = unique(data_transformed$time_f))
  norm_lvl$up <- list_norm_up
  norm_lvl$down <- list_norm_down
  
  p1 <- #data_transformed %>% 
    #dplyr::group_by(time_f, bloodmarker) %>%  # do the same calcs for each box
    #dplyr::mutate(value2 = filter_lims(bloodvalue)) %>% 
    ggplot(data=data_transformed, aes(x=time_f, 
                                      y=bloodvalue, 
                                      color=AIS, 
                                      fill=AIS, 
                                      group=AIS)) +
    stat_summary(fun.data = "mean_cl_boot", 
                 geom="smooth", 
                 se = TRUE,  
                 size=0.5, 
                 linetype=1, 
                 alpha=0.2) +
    facet_grid(bloodmarker ~ . , 
               scales = 'free',
               labeller = labeller(bloodmarker=labels)) + 
    scale_fill_manual(name = "AIS grade",
                      labels = c('A', 'B', 'C', 'D'),
                      values = c('#218317',"#457fe1", "#b30099", "#ffba00")) + 
    scale_color_manual(name = "AIS grade",
                       labels = c('A', 'B', 'C', 'D'),
                       values = c('#218317',"#457fe1", "#b30099", "#ffba00")) +
    theme_bw() + 
    ylab('Bloodvalue') + 
    xlab("Time after injury") +
    geom_hline(data = norm_lvl, aes(yintercept = up), linetype="dashed") + # horizontal line for upper limit of the norm values
    geom_hline(data = norm_lvl, aes(yintercept = down), linetype="dashed") + # horizontal line for lower limit of the norm values
    geom_rect(data = norm_lvl, aes(xmin=0, xmax=Inf, ymin=down, ymax=up), inherit.aes = FALSE, alpha = 0.01) #+ # fill in between the lines for normal ranges
  
  p2 = p1 + geom_point(data = data_transformed, aes(size="Normal ranges", shape = NA), colour = "lightgray", alpha = 0.3)
  fig <- p2 + guides(size=guide_legend("", override.aes=list(shape=15, size = 10)))
  
  # Output plot
  return(fig)
}

#----------------------------------------------------------------------------------------------

make_plot_combined_test <- function(mydata, list_norm_up, list_norm_down, times, grades, levels_T6, levels_anat, outliers) {
  # Filter based on user input
  data_transformed <- mydata[mydata$time2 %in% unlist(times, use.names=FALSE), ] # filter time points
  data_transformed <- data_transformed[data_transformed$AIS %in% unlist(grades, use.names=FALSE), ] # filter AIS grades
  data_transformed <- data_transformed[data_transformed$level_anat %in% unlist(levels_anat, use.names=FALSE), ] # filter above/below T6
  data_transformed <- data_transformed[data_transformed$level_T6 %in% unlist(levels_T6, use.names=FALSE), ] # filter cervical/thoracic/lumbar injuries
  
  # Make time column as factor variable
  data_transformed$time_f = factor(data_transformed$time2, levels=unlist(times, use.names=FALSE)) # make time a factor
  data_transformed$AIS_f = factor(data_transformed$AIS, levels=unlist(grades, use.names=FALSE)) # make AIS grades a factor
  
  # Remove NAs = if no bloodvalue at a certain time point for a certain patient
  data_transformed = data_transformed[complete.cases(data_transformed),]

  # Create grid and select the normal levels for the markers selected by user
  norm_lvl <- expand.grid(AIS = unique(data_transformed$AIS), time = unique(data_transformed$time2))
  norm_lvl$up <- list_norm_up
  norm_lvl$down <- list_norm_down
  
  # Create plot
  # Make different plots depending on if the user chose to display outliers or not
  if (outliers == 'False'){
    p1 <- ggplot(aes(x=database, y=bloodvalue), data=data_transformed) + # x-axis = databases; y-axis = blood values; fill boxplot according to number of patients per AIS grade per time point per blood marker
      geom_boxplot() + # create the box plots displaying outliers
      geom_point(fill=alpha("black", 0.01), alpha = 0.2, shape=21, colour="grey20", position=position_jitter(width=0.2, height=0.01)) +
      facet_grid(AIS~time2,scales="free") + # create facet grids : x-axis = time points; y-axis = AIS grades
      theme_bw() +
      geom_hline(data = norm_lvl, aes(yintercept = up), linetype="dashed") + # horizontal line for upper limit of the norm values
      geom_hline(data = norm_lvl, aes(yintercept = down), linetype="dashed") + # horizontal line for lower limit of the norm values
      geom_rect(data = norm_lvl, aes(xmin=0, xmax=Inf, ymin=down, ymax=up, fill="Normal range"), inherit.aes = FALSE, alpha = 0.01) + # fill in between the lines for normal ranges
      scale_fill_manual('', values = 'blue', guide = guide_legend(override.aes = list(alpha = 0.07))) +
      theme(legend.key = element_rect(color="black", linetype="dashed"))
    
  } else {
    
    ylim1 = boxplot.stats(data_transformed$bloodvalue)$stats[c(1, 5)] # determine y-axis limits
    p1 <- ggplot(aes(x=database, y=bloodvalue), data=data_transformed) + # x-axis = databases; y-axis = blood values; fill boxplot according to number of patients per AIS grade per time point per blood marker
      geom_boxplot(outlier.color = NA) + # create the box plots not displaying outliers
      geom_point(fill=alpha("black", 0.01), alpha = 0.2, shape=21, colour="grey20", position=position_jitter(width=0.2, height=0.5)) +
      facet_grid(AIS~time2,scales="free", labeller = labeller(AIS=c(A = 'AIS A', B = 'AIS B', C = 'AIS C', D = 'AIS D'))) + # create facet grids : x-axis = time points; y-axis = AIS grades
      theme_bw() +
      geom_hline(data = norm_lvl, aes(yintercept = up), linetype="dashed") + # horizontal line for upper limit of the norm values
      geom_hline(data = norm_lvl, aes(yintercept = down), linetype="dashed") + # horizontal line for lower limit of the norm values
      geom_rect(data = norm_lvl, aes(xmin=0, xmax=Inf, ymin=down, ymax=up, fill="Normal range"), inherit.aes = FALSE, alpha = 0.01) + # fill in between the lines for normal ranges
      scale_fill_manual('',values = 'blue', guide = guide_legend(override.aes = list(alpha = 0.07))) +
      theme(legend.key = element_rect(color="black", linetype="dashed")) + 
      coord_cartesian(ylim = ylim1*1.25) # rescale y-axis of each facets
  }
  
  return(p1)
}


completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

#----------------------------------------------------------------------------------------------
grades_Sygen <- c("AIS A","AIS B","AIS C","AIS D")
time_Sygen <- c("Week 0", "Week 1", "Week 2", "Week 4", "Week 8", "Week 52")
grades_Murnau <- c("A","B","C","D")
time_Murnau <- c("Week 0", "Week 1", "Week 2", "Week 4", "Week 8", "Month 3", "Month 4", "Month 5", "Month 6")
levels_T6 <- c('above_T6', 'below_T6')
levels_anat <- c('cervical', 'thoracic', 'lumbar')
list_norm_up_normalized <- c(100.0, 66.01941747572816, 9.01639344262294, 33.333333333333336, 21.64948453608249, 100.0, 17.98561151079137,
                             18.000000000000007, 100.0, 100.0, 66.00000000000001, 20.0, 19.54022988505747, 100, 18.604651162790695,
                             43.04635761589404, 5.882352941176471, 11.11111111111111, 4.225352112676056, 51.724137931034484)

list_norm_down_normalized <- c(-100.0, -66.01941747572816, -9.016393442622958, -33.333333333333336, -21.649484536082472, -100.0, -17.985611510791358,
                               -18.000000000000007, -100.0, -100.0, -66.0, -20.0, -19.54022988505747, -100, -18.604651162790695,
                               -43.04635761589405, -5.882352941176471, -11.11111111111111, -4.225352112676056, -51.724137931034484)
#----------------------------------------------------------------------------------------------
list_norm_up = c(36, 60, 40, 8, 3, 1, 52, 18, 34, 100, 440, 5.9, 8, 10.8)
list_norm_down = c(32, 45, 20, 4, 1, 0, 35, 12, 26, 80, 140, 3.8, 6, 4.3)
list_all_markers = c("MCHC",
                     "Neutrophils",
                     "Lymphocytes",
                     "Monocytes",
                     "Eosinophils",
                     "Basophils",
                     "Hematocrit",
                     "Hemoglobin",
                     "MCH",
                     "MCV",
                     "Thrombocytes",
                     "Erythrocytes",
                     "Total serum",
                     "Leucocytes")
labeller = c(MCHC = "MCHC \n g/dL",
             Neutrophils = "Neutrophils \n %",
             Lymphocytes = "Lymphocytes \n %",
             Monocytes = "Monocytes \n %",
             Eosinophils = "Eosinophils \n %",
             Basophils = "Basophils \n %",
             Hematocrit = "Hematocrit \n %",
             Hemoglobin = "Hemoglobin \n g/dL",
             MCH = "MCH \n pg",
             MCV = "MCV\n fL",
             Thrombocytes = "Thrombocytes \n thou/mcL",
             Erythrocytes = "Erythrocytes \n _/pL",
             `Total serum` = "Total serum concentration \n g/dL",
             Leucocytes = "Leucocytes \n _/nL")

subset(CBC_sygen, bloodmarker== "Thrombocytes")
for (i in c(1:dim(CBC_sygen)[1])){
  if (CBC_sygen[i, "bloodmarker"] == "Thrombocytes"){
    CBC_sygen[i, "bloodvalue"] <- CBC_sygen[i, "bloodvalue"]/1000
  }
}
subset(CBC_sygen, bloodmarker== "Thrombocytes")

fig1 <- make_plot_murnau(CBC_sygen, 
                         list_norm_up[c(7, 11, 12, 14)], 
                         list_norm_down[c(7, 11, 12, 14)], 
                         time_Sygen, 
                         grades_Sygen, 
                         list_all_markers[c(7, 11, 12, 14)], 
                         levels_T6, 
                         levels_anat, 
                         labeller[c(7, 11, 12, 14)])
fig1

#----------------------------------------------------------------------------------------------

list_all_markers = c("Erythrocytes",
                     "Hemoglobin",
                     "Hb per RBC",
                     "Hematocrit",
                     "Leucocytes",
                     "MCHC",
                     "MCV",
                     "Thrombocytes")
list_norm_up = c(5.9, 18, 34, 52, 10.8, 36, 100, 440)
list_norm_down = c(3.80, 12, 27, 35, 4.3, 32, 80, 140)
labeller = c(Erythrocytes = "Erythrocytes \n _/pL", 
             Hemoglobin = "Hemoglobin \n g/dL", 
             `Hb per RBC` = "Hb per RBC \n pg", 
             Hematocrit = "Hematocrit \n %",
             Leucocytes = "Leucocytes \n _/nL",
             MCHC = "MCHC \n g/dL",
             MCV = "MCV \n fL", 
             Thrombocytes = "Thrombocytes \n thou/mcL")

fig2 <- make_plot_murnau(CBC_murnau_long, 
                         list_norm_up[c(1,4,5,8)],
                         list_norm_down[c(1,4,5,8)],
                         time_Murnau,
                         grades_Murnau,
                         list_all_markers[c(1,4,5,8)],
                         levels_T6, 
                         levels_anat, 
                         labeller[c(1,4,5,8)])
fig2

#----------------------------------------------------------------------------------------------

list_all_times <- c("Week 0", "Week 1", "Week 2", "Week 4", "Week 8")

hb$database <- as.factor(hb$database)
hb$database <- factor(hb$database, levels=rev(levels(hb$database)))
fig3A <- make_plot_combined_test(mydata = hb,
                                list_norm_up = list_norm_up_normalized[12],
                                list_norm_down = list_norm_down_normalized[12],
                                times = list_all_times,
                                grades = grades_Murnau,
                                levels_T6,levels_anat,
                                outliers = 'True')
fig3A

# data_transformed <- mydata[mydata$time2 %in% unlist(times, use.names=FALSE), ] # filter time points
# data_transformed <- data_transformed[data_transformed$AIS %in% unlist(grades, use.names=FALSE), ] # filter AIS grades
# data_transformed <- data_transformed[data_transformed$level_anat %in% unlist(levels_anat, use.names=FALSE), ] # filter above/below T6
# data_transformed <- data_transformed[data_transformed$level_T6 %in% unlist(levels_T6, use.names=FALSE), ] # filter cervical/thoracic/lumbar injuries
# 
# # Make time column as factor variable
# data_transformed$time_n = as.numeric(data_transformed$time_f) # make time a factor
# 
# data_transformed$time_f = factor(data_transformed$time2, levels=unlist(times, use.names=FALSE)) # make time a factor
# data_transformed$AIS_f = factor(data_transformed$AIS, levels=unlist(grades, use.names=FALSE)) # make AIS grades a factor
# 
# # Remove NAs = if no bloodvalue at a certain time point for a certain patient
# data_transformed = data_transformed[complete.cases(data_transformed),]
# 
# # Create grid and select the normal levels for the markers selected by user
# norm_lvl <- expand.grid(AIS = unique(data_transformed$AIS), time = unique(data_transformed$time2))
# norm_lvl$up <- list_norm_up
# norm_lvl$down <- list_norm_down
# 
# p1 <- ggplot(data=data_transformed, aes(x=time_f, 
#                                         y=bloodvalue, 
#                                         color=database, 
#                                         fill=database, 
#                                         group=database)) +
#   stat_summary(fun.data = "mean_cl_boot", 
#                 geom="smooth", 
#                 se = TRUE,  
#                 size=0.5, 
#                 linetype=1, 
#                 alpha=0.2) +
#   facet_grid(AIS ~ . , 
#              scales = 'free',
#              labeller = labeller(AIS=c(A = 'AIS A',B = 'AIS B',C = 'AIS C',D = 'AIS D'))) + 
#   scale_fill_manual(name = "Database",
#                     labels = c('Murnau', 'Sygen'),
#                     values = c('#CC0000', "#457fe1")) + 
#   scale_color_manual(name = "Database",
#                      labels = c('Murnau', 'Sygen'),
#                      values = c('#CC0000',"#457fe1")) +
#   theme_bw() + 
#   ylab('Bloodvalue') + 
#   xlab("Time after injury") +
#   coord_cartesian(ylim = c(-40,-5)) +
#   geom_hline(data = norm_lvl, aes(yintercept = up), linetype="dashed") + # horizontal line for upper limit of the norm values
#   geom_hline(data = norm_lvl, aes(yintercept = down), linetype="dashed") + # horizontal line for lower limit of the norm values
#   geom_rect(data = norm_lvl, aes(xmin=0, xmax=Inf, ymin=down, ymax=up), inherit.aes = FALSE, alpha = 0.01) #+ # fill in between the lines for normal ranges
#   #scale_fill_manual('',values = 'blue', guide = guide_legend(override.aes = list(alpha = 0.07))) +
#   #theme(legend.key = element_rect(color="black", linetype="dashed")) 
# 
# p2 = p1 + geom_point(data = data_transformed, aes(size="Normal range", shape = NA), colour = "lightgray", alpha = 0.3)
# fig3A <- p2 + guides(size=guide_legend("", override.aes=list(shape=15, size = 10)))
# fig3A
#----------------------------------------------------------------------------------------------

mcv$database <- as.factor(mcv$database)
mcv$database <- factor(mcv$database, levels=rev(levels(mcv$database)))

fig3B <- make_plot_combined_test(mydata = mcv,
                                 list_norm_up = list_norm_up_normalized[18],
                                 list_norm_down = list_norm_down_normalized[18],
                                 times = list_all_times,
                                 grades_Murnau,levels_T6,levels_anat,
                                 outliers = 'True')
fig3B

#----------------------------------------------------------------------------------------------

list_norm_up = c(36, 60, 40, 8, 3, 1, 52, 18, 34, 100, 440, 5.9, 8, 10.8)
list_norm_down = c(32, 45, 20, 4, 1, 0, 35, 12, 26, 80, 140, 3.8, 6, 4.3)
list_all_markers = c("MCHC",
                     "Neutrophils",
                     "Lymphocytes",
                     "Monocytes",
                     "Eosinophils",
                     "Basophils",
                     "Hematocrit",
                     "Hemoglobin",
                     "MCH",
                     "MCV",
                     "Thrombocytes",
                     "Erythrocytes",
                     "Total serum",
                     "Leucocytes")
labeller = c(MCHC = "MCHC \n g/dL",
             Neutrophils = "Neutrophils \n %",
             Lymphocytes = "Lymphocytes \n %",
             Monocytes = "Monocytes \n %",
             Eosinophils = "Eosinophils \n %",
             Basophils = "Basophils \n %",
             Hematocrit = "Hematocrit \n %",
             Hemoglobin = "Hemoglobin \n g/dL",
             MCH = "MCH \n pg",
             MCV = "MCV\n fL",
             Thrombocytes = "Thrombocytes \n thou/mcL",
             Erythrocytes = "Erythrocytes \n _/pL",
             `Total serum` = "Total serum concentration \n g/dL",
             Leucocytes = "Leucocytes \n _/nL")
suppfig2A <- make_plot_murnau(CBC_sygen, 
                              list_norm_up[c(1,2,3,4,5,6,8,9,10,13)], 
                              list_norm_down[c(1,2,3,4,5,6,8,9,10,13)], 
                              time_Sygen, 
                              grades_Sygen, 
                              list_all_markers[c(1,2,3,4,5,6,8,9,10,13)], 
                              levels_T6, 
                              levels_anat, 
                              labeller[c(1,2,3,4,5,6,8,9,10,13)])
suppfig2A

#----------------------------------------------------------------------------------------------

list_all_markers = c("Erythrocytes",
                     "Hemoglobin",
                     "Hb per RBC",
                     "Hematocrit",
                     "Leucocytes",
                     "MCHC",
                     "MCV",
                     "Thrombocytes")
list_norm_up = c(5.9, 18, 34, 52, 10.8, 36, 100, 440)
list_norm_down = c(3.80, 12, 27, 35, 4.3, 32, 80, 140)
labeller = c(Erythrocytes = "Erythrocytes \n _/pL", 
             Hemoglobin = "Hemoglobin \n g/dL", 
             `Hb per RBC` = "Hb per RBC \n pg", 
             Hematocrit = "Hematocrit \n %",
             Leucocytes = "Leucocytes \n _/nL",
             MCHC = "MCHC \n g/dL",
             MCV = "MCV \n fL", 
             Thrombocytes = "Thrombocytes \n thou/mcL")
suppfig2B <- make_plot_murnau(CBC_murnau_long, 
                              list_norm_up[c(2,3,6,7)],
                              list_norm_down[c(2,3,6,7)],
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers[c(2,3,6,7)],
                              levels_T6, 
                              levels_anat, 
                              labeller[c(2,3,6,7)])
suppfig2B

library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
plottest <- ggarrange(suppfig2A, suppfig2B, ncol=1, nrow=2, labels = c("(A)", "(B)",family="TT Times New Roman"))
plottest
#----------------------------------------------------------------------------------------------

list_all_markers = c("Alkaline phosphatase",
                     "Total bilirubin",
                     "Chloride",
                     "ASAT",
                     "ALAT")
list_norm_up = c(171, 1, 106, 35, 45)
list_norm_down = c(35, 0.3, 96, 0, 0)
labeller = c(`Alkaline phosphatase` = "Alkaline phosphatase \n U/L",
             `Total bilirubin` = "Total bilirubin \n mg/dL",
             Chloride = "Chloride \n meq/L",
             ASAT = "ASAT \n U/L",
             ALAT = "ALAT \n U/L")
suppfig3A <- make_plot_murnau(liver_sygen, 
                         list_norm_up,
                         list_norm_down,
                         time_Sygen,
                         grades_Sygen,
                         list_all_markers,
                         levels_T6, 
                         levels_anat, 
                         labeller)
suppfig3A

#----------------------------------------------------------------------------------------------

list_all_markers = c("Alkaline phosphatase",
                     "Gamma-GT",
                     "Total bilirubin",
                     "ASAT",
                     "ALAT",
                     "Lactate dehydrogenase")
list_norm_up = c(171, 65, 1.1, 35, 45, 248)
list_norm_down = c(35, 0, 0, 0, 0, 0)
labeller = c(`Alkaline phosphatase` = "Alkaline phosphatase \n U/L", 
             `Gamma-GT` = "Gamma-GT \n U/L", 
             `Total bilirubin` = "Total bilirubin \n mg/dL", 
             ASAT = "ASAT \n U/L",
             ALAT = "ALAT \n U/L",
             `Lactate dehydrogenase` = "Lactate dehydrogenase \n U/L")
suppfig3B <- make_plot_murnau(liver_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig3B

#----------------------------------------------------------------------------------------------

list_all_markers = c("Albumin",
                     "Creatinin",
                     "Uric acid",
                     "Blood urea nitrogen",
                     "Calcium")
list_norm_up = c(5.4, 1, 7, 8.3, 2.66)
list_norm_down = c(3.4, 0.5, 2.3, 1.7, 2.22)
labeller = c(Albumin = "Albumin \n g/dL",
             Creatinin = "Creatinine \n mg/dL",
             `Uric acid` = "Uric acid \n mg/dL",
             `Blood urea nitrogen` = "Blood urea nitrogen \n mmol/L",
             Calcium = "Calcium \n mmol/L")

kidney_sygen <- read.csv('data/df_final_kidney_sygen.csv')
subset(kidney_sygen, bloodmarker== "Blood urea nitrogen")
for (i in c(1:dim(kidney_sygen)[1])){
  if (kidney_sygen[i, "bloodmarker"] == "Blood urea nitrogen"){
    kidney_sygen[i, "bloodvalue"] <- kidney_sygen[i, "bloodvalue"]*0.3571
  }
  if (kidney_sygen[i, "bloodmarker"] == "Calcium"){
    kidney_sygen[i, "bloodvalue"] <- kidney_sygen[i, "bloodvalue"]*0.2495
  }
}
subset(kidney_sygen, bloodmarker== "Blood urea nitrogen")

suppfig4A <- make_plot_murnau(kidney_sygen, 
                              list_norm_up,
                              list_norm_down,
                              time_Sygen,
                              grades_Sygen,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig4A

#----------------------------------------------------------------------------------------------

list_all_markers = c("Calcium",
                     "Creatinin",
                     "Total proteins",
                     "Blood urea nitrogen")
list_norm_up = c(2.66, 1, 8.2, 8.3)
list_norm_down = c(2.22, 0.5, 5.7, 1.7)
labeller = c(Calcium = "Calcium \n mmol/L", 
             Creatinin = "Creatinine \n mg/dL", 
             `Total proteins` = "Total proteins \n g/dL",
             `Blood urea nitrogen` = "Blood urea nitrogen \n mmol/L")
suppfig4B <- make_plot_murnau(kidney_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig4B

plottest <- ggarrange(suppfig4A, suppfig4B, ncol=1, nrow=2, labels = c("(A)", "(B)",family="TT Times New Roman"))
plottest

#----------------------------------------------------------------------------------------------

list_all_markers = c("Amylase")
list_norm_up = c(115)
list_norm_down = c(0)
labeller = c(Amylase = "Amylase \n U/dL")
suppfig6A <- make_plot_murnau(pancreas_sygen, 
                              list_norm_up,
                              list_norm_down,
                              time_Sygen,
                              grades_Sygen,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig6A

#----------------------------------------------------------------------------------------------

list_all_markers = c("Potassium",
                     "Sodium")
list_norm_up = c(5.1, 148)
list_norm_down = c(3.5, 136)
labeller = c(Potassium = "Potassium \n mmol/L",
             Sodium = "Sodium \n mmol/L")
suppfig5A <- make_plot_murnau(muscle_sygen, 
                              list_norm_up,
                              list_norm_down,
                              time_Sygen,
                              grades_Sygen,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig5A

#----------------------------------------------------------------------------------------------

list_all_markers = c("Prothrombin time",
                     "Cholesterol",
                     "Carbon dioxide",
                     "Glucose",
                     "Triglycerides")
list_norm_up = c(26, 200, 29, 5.9, 200)
list_norm_down = c(37, 0, 23, 4.1, 0)
labeller = c(`Prothrombin time` = "Prothrombin time \n s",
             Cholesterol = "Cholesterol \n mg/dL",
             `Carbon dioxide` = "Carbon dioxide \n meq/L",
             Glucose = "Glucose \n mg/dL",
             Triglycerides = "Triglycerides \n mg/dL")

subset(rest_sygen, bloodmarker== "Glucose")
for (i in c(1:dim(rest_sygen)[1])){
  if (rest_sygen[i, "bloodmarker"] == "Glucose"){
    rest_sygen[i, "bloodvalue"] <- rest_sygen[i, "bloodvalue"]*0.0555
  }
}
subset(rest_sygen, bloodmarker== "Glucose")

suppfig7A <- make_plot_murnau(rest_sygen, 
                              list_norm_up,
                              list_norm_down,
                              time_Sygen,
                              grades_Sygen,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig7A

#----------------------------------------------------------------------------------------------

list_all_markers = c("Amylase",
                     "Lipase")
list_norm_up = c(115, 80)
list_norm_down = c(0, 0)
labeller = c(Amylase = "Amylase \n U/L", 
             Lipase = "Lipase \n U/L")
suppfig6B <- make_plot_murnau(pancreas_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig6B

plottest <- ggarrange(suppfig6A, suppfig6B, ncol=1, nrow=2, labels = c("(A)", "(B)",family="TT Times New Roman"))
plottest

#----------------------------------------------------------------------------------------------

list_all_markers = c("Cholinesterase",
                     "Potassium",
                     "Sodium")
list_norm_up = c(12920, 5.1, 148)
list_norm_down = c(5320, 3.5, 136)
labeller = c(Cholinesterase = "Cholinesterase \n U/L", 
             Potassium = "Potassium \n mmol/L", 
             Sodium = "Sodium \n mmol/L")
suppfig5B <- make_plot_murnau(muscle_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig5B

plottest <- ggarrange(suppfig5A, suppfig5B, ncol=1, nrow=2, labels = c("(A)", "(B)",family="TT Times New Roman"))
plottest

#----------------------------------------------------------------------------------------------

list_all_markers = c("Protein C reactive",
                     "Glucose",
                     "Prothrombin time",
                     "Partial thromboplastin time",
                     "Quick test")
list_norm_up = c(0.5, 5.9, 1.1, 40, 127)
list_norm_down = c(0, 4.1, 0, 26, 80)
labeller = c(`Protein C reactive` = "CRP \n mg/dL", 
             Glucose = "Glucose \n mmol/L", 
             `Prothrombin time` = "INR \n without unit",
             `Partial thromboplastin time` = "Partial thromboplastin time \n s",
             `Quick test` = "Quick test \n %")
suppfig7B <- make_plot_murnau(rest_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig7B

plottest <- ggarrange(suppfig7A, suppfig7B, ncol=1, nrow=2, labels = c("(A)", "(B)",family="TT Times New Roman"))
plottest

#----------------------------------------------------------------------------------------------

suppfig1 <- ggplot(df_age_combined, aes(x=age)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.1, color="#FF6666") +
  facet_grid(. ~ dataset) +
  theme_bw()+
  theme(axis.text.y.left =element_text(colour="red"),
        axis.title.y.left =element_text(colour="red"))
suppfig1 <- suppfig1 + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "proportion (%)"))
suppfig1
