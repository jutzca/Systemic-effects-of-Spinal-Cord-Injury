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
  
  # Create a error message to deliver when filters applied end up removing all data when taken together
  if (dim(data_transformed)[1] == 0){
    rects <- data.frame(x = 1:1,
                        colors = c("red"),
                        text = "No data left")
    p <- ggplot(rects, aes(x, y = 0, fill = colors, label = text)) +
      geom_tile(width = .25, height = .1) + # make square tiles
      geom_text(color = "white") + # add white text in the middle
      scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
      coord_fixed() + # make sure tiles are square
      theme_void() # remove any axis markings
    return(p)
  }
  
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
  
  # Create a error message to deliver when filters applied end up removing all data when taken together
  if (dim(data_transformed)[1] == 0){
    rects <- data.frame(x = 1:1,
                        colors = c("red"),
                        text = "No data left")
    p <- ggplot(rects, aes(x, y = 0, fill = colors, label = text)) +
      geom_tile(width = .25, height = .1) + # make square tiles
      geom_text(color = "white") + # add white text in the middle
      scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
      coord_fixed() + # make sure tiles are square
      theme_void() # remove any axis markings
    return(p)
  }
  
  # Create grid and select the normal levels for the markers selected by user
  norm_lvl <- expand.grid(AIS = unique(data_transformed$AIS), time = unique(data_transformed$time2))
  norm_lvl$up <- list_norm_up
  norm_lvl$down <- list_norm_down
  
  # Create plot
  # Make different plots depending on if the user chose to display outliers or not
  if (outliers == 'False'){
    # p1 <- ggplot(aes(x=database, y=bloodvalue, fill=nb), data=data_transformed) + # x-axis = databases; y-axis = blood values; fill boxplot according to number of patients per AIS grade per time point per blood marker
    #   geom_boxplot() + # create the box plots displaying outliers
    #   facet_grid(AIS~time2,scales="free") + # create facet grids : x-axis = time points; y-axis = AIS grades
    #   scale_fill_gradientn(colors = colors) + # fill according to a gradient of colours
    #   theme_bw() +
    #   labs(fill = "Number of patients") + # set legend of filling for box plots
    #   geom_hline(data = norm_lvl, aes(yintercept = up), linetype="dashed") + # horizontal line for upper limit of the norm values
    #   geom_hline(data = norm_lvl, aes(yintercept = down), linetype="dashed") + # horizontal line for lower limit of the norm values
    #   geom_rect(data = norm_lvl, aes(xmin=0, xmax=Inf, ymin=down, ymax=up), inherit.aes = FALSE, fill="blue", alpha = 0.01) # fill in between the lines for normal ranges
    # 
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

#----------------------------------------------------------------------------------------------
grades_Sygen <- c("AIS A","AIS B","AIS C","AIS D")
time_Sygen <- c("Week 0", "Week 1", "Week 2", "Week 4", "Week 8", "Week 52")
grades_Murnau <- c("A","B","C","D")
time_Murnau <- c("Week 0", "Week 1", "Week 2", "Week 4", "Week 8", "Month 3", "Month 4", "Month 5", "Month 6")
levels_T6 <- c('above_T6', 'below_T6')
levels_anat <- c('cervical', 'thoracic', 'lumbar')
#----------------------------------------------------------------------------------------------

fig1 <- make_plot_murnau(CBC_sygen, 
                         list_norm_up[c(1, 7, 8, 10, 11, 12, 14)], 
                         list_norm_down[c(1, 7, 8, 10, 11, 12, 14)], 
                         time_Sygen, 
                         grades_Sygen, 
                         list_all_markers[c(1, 7, 8, 10, 11, 12, 14)], 
                         levels_T6, 
                         levels_anat, 
                         labeller[c(1, 7, 8, 10, 11, 12, 14)])
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
             Thrombocytes = "Thrombocytes \n thou/μL")

fig2 <- make_plot_murnau(CBC_murnau_long, 
                         list_norm_up,
                         list_norm_down,
                         time_Murnau,
                         grades_Murnau,
                         list_all_markers,
                         levels_T6, 
                         levels_anat, 
                         labeller)
fig2

#----------------------------------------------------------------------------------------------

list_all_times <- c("Week 0", "Week 1", "Week 2", "Week 4", "Week 8")

hb$database <- as.factor(hb$database)
hb$database <- factor(hb$database, levels=rev(levels(hb$database)))
fig3A <- make_plot_combined_test(mydata = hb,
                                list_norm_up = list_norm_up_normalized[12],
                                list_norm_down = list_norm_down_normalized[12],
                                times = list_all_times,
                                grades_Murnau,levels_T6,levels_anat,
                                outliers = 'True')
fig3A

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
suppfig2A <- make_plot_murnau(liver_sygen, 
                         list_norm_up,
                         list_norm_down,
                         time_Sygen,
                         grades_Sygen,
                         list_all_markers,
                         levels_T6, 
                         levels_anat, 
                         labeller)
suppfig2A

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
suppfig2B <- make_plot_murnau(liver_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig2B

#----------------------------------------------------------------------------------------------

list_all_markers = c("Albumin",
                     "Creatinin",
                     "Uric acid",
                     "Blood urea nitrogen",
                     "Calcium")
list_norm_up = c(5.4, 1, 7, 23.24, 10.66)
list_norm_down = c(3.4, 0.5, 2.3, 4.76, 8.9)
labeller = c(Albumin = "Albumin \n g/dL",
             Creatinin = "Creatinine \n mg/dL",
             `Uric acid` = "Uric acid \n mg/dL",
             `Blood urea nitrogen` = "Blood urea nitrogen \n mg/dL",
             Calcium = "Calcium \n mg/dL")
suppfig3A <- make_plot_murnau(kidney_sygen, 
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
suppfig3B <- make_plot_murnau(kidney_murnau_long, 
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

list_all_markers = c("Amylase")
list_norm_up = c(115)
list_norm_down = c(0)
labeller = c(Amylase = "Amylase \n U/dL")
suppfig5A <- make_plot_murnau(pancreas_sygen, 
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

list_all_markers = c("Potassium",
                     "Sodium")
list_norm_up = c(5.1, 148)
list_norm_down = c(3.5, 136)
labeller = c(Potassium = "Potassium \n meq/L",
             Sodium = "Sodium \n meq/L")
suppfig4A <- make_plot_murnau(muscle_sygen, 
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

list_all_markers = c("Prothrombin time",
                     "Cholesterol",
                     "Carbon dioxide",
                     "Glucose",
                     "Triglycerides")
list_norm_up = c(26, 200, 29, 106, 200)
list_norm_down = c(37, 0, 23, 74, 0)
labeller = c(`Prothrombin time` = "Prothrombin time \n s",
             Cholesterol = "Cholesterol \n mg/dL",
             `Carbon dioxide` = "Carbon dioxide \n meq/L",
             Glucose = "Glucose \n mg/dL",
             Triglycerides = "Triglycerides \n mg/dL")
suppfig6A <- make_plot_murnau(rest_sygen, 
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

list_all_markers = c("Amylase",
                     "Lipase")
list_norm_up = c(115, 80)
list_norm_down = c(0, 0)
labeller = c(Amylase = "Amylase \n U/L", 
             Lipase = "Lipase \n U/L")
suppfig5B <- make_plot_murnau(pancreas_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig5B

#----------------------------------------------------------------------------------------------

list_all_markers = c("Cholinesterase",
                     "Potassium",
                     "Sodium")
list_norm_up = c(12920, 5.1, 148)
list_norm_down = c(5320, 3.5, 136)
labeller = c(Cholinesterase = "Cholinesterase \n U/L", 
             Potassium = "Potassium \n mmol/L", 
             Sodium = "Sodium \n mmol/L")
suppfig4B <- make_plot_murnau(muscle_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig4B

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
suppfig6B <- make_plot_murnau(rest_murnau_long, 
                              list_norm_up,
                              list_norm_down,
                              time_Murnau,
                              grades_Murnau,
                              list_all_markers,
                              levels_T6, 
                              levels_anat, 
                              labeller)
suppfig6B

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
