# ------------------------------------------------------------------------------------------------------------------
# Statistical anaylsis - blood markers in SCI
#
# July 2020
# L. Bourguignon
# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------

# Load packages ----
library(lme4)
library(sjPlot)
library(stats)
library(rlist)
library(multcomp)
library(xtable)
library(stringr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(lsmeans)
library(car)

# Set working directory (where to store results) ----


# ------------------------------------------------------------------------------------------------------------------

# Define functions ----

# # Create tab_model .html files summarizing results from linear mixed effects models
stat_test <- function(df, name)
{
  df = df[complete.cases(df),] #remove rows (= 1 time point for 1 individual for 1 blood marker) with missing values
  remove_AIS = c('E', '', 'AIS E', 'ND')
  df = df[!df$AIS %in% remove_AIS,] #filter out patients with AIS E or missing grades at baseline (Sygen) or acute 1 stage (Murnau)
  
  test <- lmer(bloodvalue ~ time + sex + AIS + age + level +
                 (1|random_effect), data = df, na.action = na.omit)
  
  return (tab_model(test, title = name, file = paste("results_stats/summary_table_", name, ".html", sep = "")))
}

# # ANOVA test to assess the overall contribution of each variable in explaining blood values for inidividual dataset
anova_test <- function(df)
{
  df = df[complete.cases(df),] #remove rows (= 1 time point for 1 individual for 1 blood marker) with missing values
  remove_AIS = c('E', '', 'AIS E', 'ND')
  df = df[!df$AIS %in% remove_AIS,] #filter out patients with AIS E or missing grades at baseline (Sygen) or acute 1 stage (Murnau)
  
  #test <- aov(bloodvalue ~ age + level + sex + AIS, data = df, na.action = na.omit) # fit ANOVA test 
  test <- Anova(lm1 <- lmer(bloodvalue ~ time + AIS + age + level + sex + multiple_trauma + AIS*time + (1|random_effect), data = df, na.action = na.omit), type="III")
  
  return (test)
}

slope <- function(df)
{
  df = df[complete.cases(df),] #remove rows (= 1 time point for 1 individual for 1 blood marker) with missing values
  remove_AIS = c('E', '', 'AIS E', 'ND')
  df = df[!df$AIS %in% remove_AIS,] #filter out patients with AIS E or missing grades at baseline (Sygen) or acute 1 stage (Murnau)
  
  df$time = as.numeric(as.character(df$time))
  
  test <- lmer(bloodvalue ~ time*AIS + age + level + sex+(1|random_effect), data = df, na.action = na.omit) # fit ANOVA test
  m.lst <- lstrends(test, "AIS", var = "time")
  
  return(summary(pairs(m.lst)))
}

# # ANOVA test to assess the overall contribution of each variable in explaining blood values when comparing the 2 datasets
anova_combined <- function(df)
{
  df = df[complete.cases(df),] #remove rows (= 1 time point for 1 individual for 1 blood marker) with missing values
  remove_AIS = c('E', '', 'AIS E', 'ND')
  df = df[!df$AIS %in% remove_AIS,] #filter out patients with AIS E or missing grades at baseline (Sygen) or acute 1 stage (Murnau)
  
  #test <- aov(bloodvalue ~ time*database + AIS + level_T6 + sex + age + (1|random_effect), data = df, na.action = na.omit)
  test <- Anova(lm1 <- lmer(bloodvalue ~ time*database + AIS + level + sex + age + multiple_trauma + (1|random_effect), data = df, na.action = na.omit), type="III")
  return (test)
}

# # Adjusted lmer to compare pairs of severity grades in explaining blood values
pairwise_test <- function(df, name)
{
  df = df[complete.cases(df),] #remove rows (= 1 time point for 1 individual for 1 blood marker) with missing values
  remove_AIS = c('E', '', 'AIS E', 'ND')
  df = df[!df$AIS %in% remove_AIS,] #filter out patients with AIS E or missing grades at baseline (Sygen) or acute 1 stage (Murnau)
  
  test <- lmer(bloodvalue ~ time + sex + AIS + age + level + multiple_trauma + (1|random_effect), data = df, na.action = na.omit)
  test_comp <- glht(test, linfct = mcp(AIS = "Tukey"))
  
  return(summary(test_comp))
}

anova_test_1time <- function(df)
{
  df = df[complete.cases(df),] #remove rows (= 1 time point for 1 individual for 1 blood marker) with missing values
  remove_AIS = c('E', '', 'AIS E', 'ND')
  df = df[!df$AIS %in% remove_AIS,] #filter out patients with AIS E or missing grades at baseline (Sygen) or acute 1 stage (Murnau)
  
  test <- aov(bloodvalue ~ AIS + age + level + sex, data = df, na.action = na.omit) # fit ANOVA test 
  result <- TukeyHSD(test, "AIS")
  
  return (result)
}

add_age_murnau <- function(df){
  list_age <- c(71, 33, 56, 28, 85, 61, 19, 51, 75, 79, 51, 40)
  count = 0
  for (i in c(60389, 60423, 61276, 60459, 60748, 60523, 60495, 60504, 60731, 60444, 60445, 60447)){
    count = count +1
    index <- df$random_effect == i
    df$age[index] <- list_age[count]
  }
  return (df)
}

# ------------------------------------------------------------------------------------------------------------------

# Load data ----

# # Murnau
df_amylase = read.csv("stats/df_stat_analysis_amylase_review.csv", header = TRUE, sep = ",")
df_ap = read.csv("stats/df_stat_analysis_ap_review.csv", header = TRUE, sep = ",")
df_calcium = read.csv("stats/df_stat_analysis_calcium_review.csv", header = TRUE, sep = ",")
df_che = read.csv("stats/df_stat_analysis_che_review.csv", header = TRUE, sep = ",")
df_creatinin = read.csv("stats/df_stat_analysis_creatinin_review.csv", header = TRUE, sep = ",")
df_crp = read.csv("stats/df_stat_analysis_crp_review.csv", header = TRUE, sep = ",")
df_ery = read.csv("stats/df_stat_analysis_ery_review.csv", header = TRUE, sep = ",")
df_gamma_gt = read.csv("stats/df_stat_analysis_gamma_gt_review.csv", header = TRUE, sep = ",")
df_bilirubin = read.csv("stats/df_stat_analysis_gesamt_bilirubin_review.csv", header = TRUE, sep = ",")
df_gesamteiweiss = read.csv("stats/df_stat_analysis_gesamteiweiss_review.csv", header = TRUE, sep = ",")
df_glucose = read.csv("stats/df_stat_analysis_glucose_review.csv", header = TRUE, sep = ",")
df_harnstoff = read.csv("stats/df_stat_analysis_harnstoff_review.csv", header = TRUE, sep = ",")
df_hb = read.csv("stats/df_stat_analysis_hb_review.csv", header = TRUE, sep = ",")
df_hbe = read.csv("stats/df_stat_analysis_hbe_review.csv", header = TRUE, sep = ",")
df_hk = read.csv("stats/df_stat_analysis_hk_review.csv", header = TRUE, sep = ",")
df_inr = read.csv("stats/df_stat_analysis_inr_review.csv", header = TRUE, sep = ",")
df_kalium = read.csv("stats/df_stat_analysis_kalium_review.csv", header = TRUE, sep = ",")
df_ldh = read.csv("stats/df_stat_analysis_ldh_review.csv", header = TRUE, sep = ",")
df_leuco_nl = read.csv("stats/df_stat_analysis_leuco_nl_review.csv", header = TRUE, sep = ",")
df_lipase = read.csv("stats/df_stat_analysis_lipase_review.csv", header = TRUE, sep = ",")
df_mchc = read.csv("stats/df_stat_analysis_mchc_review.csv", header = TRUE, sep = ",")
df_mcv_m = read.csv("stats/df_stat_analysis_mcv_review.csv", header = TRUE, sep = ",")
df_natrium = read.csv("stats/df_stat_analysis_natrium_review.csv", header = TRUE, sep = ",")
df_ptt = read.csv("stats/df_stat_analysis_ptt_review.csv", header = TRUE, sep = ",")
df_quick = read.csv("stats/df_stat_analysis_quick_review.csv", header = TRUE, sep = ",")
df_thrombo = read.csv("stats/df_stat_analysis_thrombo_review.csv", header = TRUE, sep = ",")
df_got = read.csv("stats/df_stat_analysis_got_review.csv", header = TRUE, sep = ",")
df_gpt = read.csv("stats/df_stat_analysis_gpt_review.csv", header = TRUE, sep = ",")

list_murnau = list(df_amylase, df_ap, df_calcium, df_che, df_creatinin, df_crp, df_ery, df_gamma_gt, df_bilirubin, df_gesamteiweiss,
                   df_glucose, df_harnstoff, df_hb, df_hbe, df_hk, df_inr, df_kalium, df_ldh, df_leuco_nl, df_lipase,
                   df_mchc, df_mcv_m, df_natrium, df_ptt, df_quick, df_thrombo, df_got, df_gpt)

list_names_murnau = c("Amylase", "Alkaline phosphatase", "Calcium", "Cholinesterase", "Creatinine",
                      "CRP", "Erythrocytes", "Gamma GT", "Bilirubin", "Proteins", "Glucose", "Blood urea",
                      "Hemoglobin", "Hemoglobin per erythrocyte", "Hematocrit", "INR", "Potassium",
                      "Lactate dehydrogenase", "Leucocytes", "Lipase","MCHC", "MCV", "Sodium", "Prothrombin time",
                      "Quick test", "Thrombocytes", "ASAT", "ALAT")

# # Sygen
df_alb = read.csv("stats/df_stat_analysis_sygen_alb_review.csv", header = TRUE, sep = ",")
df_alk = read.csv("stats/df_stat_analysis_sygen_alk_review.csv", header = TRUE, sep = ",")
df_amy = read.csv("stats/df_stat_analysis_sygen_amy_review.csv", header = TRUE, sep = ",")
df_ap1 = read.csv("stats/df_stat_analysis_sygen_ap1_review.csv", header = TRUE, sep = ",")
df_bc9 = read.csv("stats/df_stat_analysis_sygen_bc9_review.csv", header = TRUE, sep = ",")
df_bt = read.csv("stats/df_stat_analysis_sygen_bt_review.csv", header = TRUE, sep = ",")
df_bua = read.csv("stats/df_stat_analysis_sygen_bua_review.csv", header = TRUE, sep = ",")
df_bun = read.csv("stats/df_stat_analysis_sygen_bun_review.csv", header = TRUE, sep = ",")
df_cab = read.csv("stats/df_stat_analysis_sygen_cab_review.csv", header = TRUE, sep = ",")
df_chc = read.csv("stats/df_stat_analysis_sygen_chc_review.csv", header = TRUE, sep = ",")
df_cho = read.csv("stats/df_stat_analysis_sygen_chi_review.csv", header = TRUE, sep = ",")
df_ck0 = read.csv("stats/df_stat_analysis_sygen_ck0_review.csv", header = TRUE, sep = ",")
df_clb = read.csv("stats/df_stat_analysis_sygen_clb_review.csv", header = TRUE, sep = ",")
df_co2 = read.csv("stats/df_stat_analysis_sygen_co2_review.csv", header = TRUE, sep = ",")
df_dbn = read.csv("stats/df_stat_analysis_sygen_dbn_review.csv", header = TRUE, sep = ",")
df_ddn = read.csv("stats/df_stat_analysis_sygen_ddn_review.csv", header = TRUE, sep = ",")
df_djn = read.csv("stats/df_stat_analysis_sygen_djn_review.csv", header = TRUE, sep = ",")
df_dln = read.csv("stats/df_stat_analysis_sygen_dln_review.csv", header = TRUE, sep = ",")
df_dnn = read.csv("stats/df_stat_analysis_sygen_dnn_review.csv", header = TRUE, sep = ",")
df_glu = read.csv("stats/df_stat_analysis_sygen_glu_review.csv", header = TRUE, sep = ",")
df_hct = read.csv("stats/df_stat_analysis_sygen_hct_review.csv", header = TRUE, sep = ",")
df_hgb = read.csv("stats/df_stat_analysis_sygen_hgb_review.csv", header = TRUE, sep = ",")
df_kb = read.csv("stats/df_stat_analysis_sygen_kb_review.csv", header = TRUE, sep = ",")
df_mch = read.csv("stats/df_stat_analysis_sygen_mch_review.csv", header = TRUE, sep = ",")
df_mcv = read.csv("stats/df_stat_analysis_sygen_mcv_review.csv", header = TRUE, sep = ",")
df_nab = read.csv("stats/df_stat_analysis_sygen_nab_review.csv", header = TRUE, sep = ",")
df_plt = read.csv("stats/df_stat_analysis_sygen_plt_review.csv", header = TRUE, sep = ",")
df_pt1 = read.csv("stats/df_stat_analysis_sygen_pt1_review.csv", header = TRUE, sep = ",")
df_rbc = read.csv("stats/df_stat_analysis_sygen_rbc_review.csv", header = TRUE, sep = ",")
df_sgo = read.csv("stats/df_stat_analysis_sygen_sgo_review.csv", header = TRUE, sep = ",")
df_sgp = read.csv("stats/df_stat_analysis_sygen_sgp_review.csv", header = TRUE, sep = ",")
df_stp = read.csv("stats/df_stat_analysis_sygen_stp_review.csv", header = TRUE, sep = ",")
df_tri = read.csv("stats/df_stat_analysis_sygen_tri_review.csv", header = TRUE, sep = ",")
df_wbc = read.csv("stats/df_stat_analysis_sygen_wbc_review.csv", header = TRUE, sep = ",")

for (i in c(1:dim(df_plt)[1])){
  df_plt[i, "bloodvalue"] <- df_plt[i, "bloodvalue"]/1000
  df_bun[i, "bloodvalue"] <- df_bun[i, "bloodvalue"]*0.3571
  df_cab[i, "bloodvalue"] <- df_cab[i, "bloodvalue"]*0.2495
  df_glu[i, "bloodvalue"] <- df_glu[i, "bloodvalue"]*0.0555
}

list_sygen = list(df_alb, df_alk, df_amy, df_ap1, df_bc9, df_bt, df_bua, df_bun, df_cab, df_chc, df_cho,
                  df_ck0, df_clb, df_co2, df_dbn, df_ddn, df_djn, df_dln, df_dnn, df_glu, df_hct, df_hgb,
                  df_kb, df_mch, df_mcv, df_nab, df_plt, df_rbc, df_sgo, df_sgp, df_stp, df_tri, df_wbc)

list_names_sygen = c("Albumin", "Alkaline phosphatase", "Amylase", "Prothrombin time", "Creatinine",
                     "Total bilirubin", "Uric acid", "Blood urea nitrogen", "Calcium", "MCHC",
                     "Cholesterol", "Creatin phosphokinase", "Chloride", "Carbon dioxide", "Neutrophils",
                     "Lymphocytes", "Monocytes", "Eosinophils", "Basophils", "Glucose", "Hematocrit", 
                     "Hemoglobin", "Potassium", "MCH", "MCV", "Sodium", "Thrombocytes", "Erythrocytes",
                     "ASAT", "ALAT", "Total serum", "Triglycerides", "Leucocytes")

# ------------------------------------------------------------------------------------------------------------------

# Statistical analysis 

# ------------------------------------------------------------------------------------------------------------------
# # Murnau
# ------------------------------------------------------------------------------------------------------------------

# # # ANOVA tests
test <- anova_test(df_amylase)

# Create table to store results
df_results_ANOVA_Murnau <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Bloodmarker", "Variable", "Chisq", "Df", "Pr(>Chisq)")
colnames(df_results_ANOVA_Murnau) <- x

# List of variables for which one wants to report the ANOVA results
list_variable = c('','time', 'AIS', 'age', 'level', 'sex', 'multiple_trauma', 'time:AIS')


for (i in 1:length(list_names_murnau)){ # go through each haematological marker
  df = list_murnau[[i]] # select corresponding data frame
  df = add_age_murnau(df) # add the age information for 12 patients which data have been updated
  colnames(df)[which(names(df) == "level_T6")] <- "level"
  row_name <- c(list_names_murnau[i]) # select corresponding haematological marker name
  #df$time = as.factor(df$time) # can choose to take time as a continuous variable or a factor variable (for factor, uncomment this line)
  temp <- anova_test(df) # run ANOVA test
  for (j in c(2, 3, 4, 5, 6, 7, 8)) { # for each variable of interest, report results
    row_value <- c(formatC(temp[["Chisq"]][j]), # extract degree of freedom
                   formatC(temp[["Df"]][j], digits = 3), # extract F-value with 3 digits
                   formatC(temp[["Pr(>Chisq)"]][j], format = "e", digits = 2)) # extract p-value with format X.XXE-X
    df_results_ANOVA_Murnau[nrow(df_results_ANOVA_Murnau) + 1,"Bloodmarker"] <- row_name # create a new row and give it the name of the haematological marker
    df_results_ANOVA_Murnau[nrow(df_results_ANOVA_Murnau),"Variable"] <- list_variable[j] # report name of variable j
    df_results_ANOVA_Murnau[nrow(df_results_ANOVA_Murnau), 3:5] <- row_value # report results for haematological marker i, variable j
  }
}

write.csv(df_results_ANOVA_Murnau,'results_stats/df_results_ANOVA_Murnau_interactions_type3_RE_agemodified_review_bintrauma.csv') # save ANOVA results in .csv format

# ------------------------------------------------------------------------------------------------------------------

# # # Pairwise comparison of severuty grades in explaining blood values
# # Can do the same analysis for time by replacing all mentions of 'AIS' by 'time', including in the pairwise_test function
# pairwise_test(df_amylase, "Amylase")

# Create table to store results
df_results_glht_Murnau_AIS <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Bloodmarker", "Comparison", "Effect", "p-value", "adjusted p-value")
colnames(df_results_glht_Murnau_AIS) <- x

for (i in 1:length(list_names_murnau)){ # go through each haematological marker
  df = list_murnau[[i]] # select corresponding data frame
  df = add_age_murnau(df) # add the age information for 12 patients which data have been updated
  colnames(df)[which(names(df) == "level_T6")] <- "level"
  row_name <- c(list_names_murnau[i]) # select corresponding haematological marker name
  df$AIS <- as.factor(df$AIS) # make sure AIS grade/time column is a factor
  #df$time <- as.factor(df$time) # make sure AIS grade/time column is a factor
  temp <- pairwise_test(df) # run lmer and glht functions to get the pairwise differences between different AIS grade/time point
  for (j in 1:length(names(temp$test$coefficients))){ # for each pairs compared, report results
    row_value <- c(names(temp$test$coefficients)[j], # extract name of pair
                   round(temp$test$coefficients[[j]], digits = 3), # extract coefficient
                   formatC(temp$test$pfunction()[[j]], format = "e", digits = 2), # extract p function
                   formatC(temp$test$pvalues[[j]], format = "e", digits = 2)) # extract p values
    df_results_glht_Murnau_AIS[nrow(df_results_glht_Murnau_AIS) + 1, "Bloodmarker"] <- row_name # create a new row and give it the name of the haematological marker
    df_results_glht_Murnau_AIS[nrow(df_results_glht_Murnau_AIS), 2:5] <- row_value # report results for haematological marker i, pair j
  }
}

write.csv(df_results_glht_Murnau_AIS,'results_stats/df_results_glht_Murnau_AIS_agecorrected_review_bintrauma.csv') # save lmer results in .csv format

# ------------------------------------------------------------------------------------------------------------------

# # # Slope analysis - comparing bloodvalues for different AIS grades over time
# slope(df_amylase)

# Create table to store results
df_results_slope_Murnau <- data.frame(matrix(ncol = 7, nrow = 0))
x <- c("Bloodmarker", "Contrast", "Estimate", "SE", "df","t.ratio", "p- value")
colnames(df_results_slope_Murnau) <- x

for (i in 1:length(list_names_murnau)){ # go through each haematological marker
  df = list_murnau[[i]] # select corresponding data frame
  row_name <- c(list_names_murnau[i]) # select corresponding haematological marker name
  #df$AIS <- as.factor(df$AIS) # make sure AIS grade/time is a factor
  temp <- slope(df) # run slope analysis
  for (j in 1:length(levels(temp$contrast))){ # for each pair compared
    row_value <- c(levels(temp$contrast)[j], # report name of the pair
                   round(temp$estimate[[j]], digits = 3), # report the estimate of effect
                   round(temp$SE[[j]], digits = 3), # 
                   round(temp$df[[j]], digits = 3), # report degree of freedom 
                   round(temp$t.ratio[[j]], digits = 3), # report t value
                   formatC(temp$p.value[[j]], format = "e", digits = 2)) # report p- value
    df_results_slope_Murnau[nrow(df_results_slope_Murnau) + 1, "Bloodmarker"] <- row_name # create a new row and give it the name of the haematological marker
    df_results_slope_Murnau[nrow(df_results_slope_Murnau), 2:7] <- row_value # report results for haematological marker i, pair j
  }
}

write.csv(df_results_slope_Murnau,'results_stats/df_results_slope_Murnau.csv') # save slope analysis results in .csv format

# ------------------------------------------------------------------------------------------------------------------

# # # Pairwise comparison - comparing bloodvalues for different AIS grades at end point
# slope(df_amylase)

# Create table to store results
df_results_glht_Murnau_AIS_endpoint <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Bloodmarker", "Comparison", "Effect", "p-value", "adjusted p-value")
colnames(df_results_glht_Murnau_AIS_endpoint) <- x

for (i in 1:length(list_names_murnau)){ # go through each haematological marker
  df = list_murnau[[i]] # select corresponding data frame
  df_subset = subset(df, time == 7)
  df = add_age_murnau(df) # add the age information for 12 patients which data have been updated
  row_name <- c(list_names_murnau[i]) # select corresponding haematological marker name
  df$AIS <- as.factor(df$AIS) # make sure AIS grade/time column is a factor
  #df$time <- as.factor(df$time) # make sure AIS grade/time column is a factor
  temp <- pairwise_test(df) # run lmer and glht functions to get the pairwise differences between different AIS grade/time point
  for (j in 1:length(names(temp$test$coefficients))){ # for each pairs compared, report results
    row_value <- c(names(temp$test$coefficients)[j], # extract name of pair
                   round(temp$test$coefficients[[j]], digits = 3), # extract coefficient
                   formatC(temp$test$pfunction()[[j]], format = "e", digits = 2), # extract p function
                   formatC(temp$test$pvalues[[j]], format = "e", digits = 2)) # extract p values
    df_results_glht_Murnau_AIS_endpoint[nrow(df_results_glht_Murnau_AIS_endpoint) + 1, "Bloodmarker"] <- row_name # create a new row and give it the name of the haematological marker
    df_results_glht_Murnau_AIS_endpoint[nrow(df_results_glht_Murnau_AIS_endpoint), 2:5] <- row_value # report results for haematological marker i, pair j
  }
}

write.csv(df_results_glht_Murnau_AIS_endpoint,'results_stats/df_results_glht_Murnau_AIS_endpoint.csv') # save lmer results in .csv format

# ------------------------------------------------------------------------------------------------------------------

df_results_glht_Murnau_AIS_endpoint <- data.frame(matrix(ncol = 6, nrow = 0))
x <- c("Bloodmarker", "Variable", "diff", "lwr", "upr", "p.adj")
colnames(df_results_glht_Murnau_AIS_endpoint) <- x

for (i in 1:length(list_names_murnau)){ # go through each haematological marker
  df = list_murnau[[i]] # select corresponding data frame
  df_subset = subset(df, time == 7) # add the age information for 12 patients which data have been updated
  df_subset$AIS <- as.factor(df_subset$AIS)
  temp_fct <- anova_test_1time(df_subset)
  temp <- data.frame(temp_fct$AIS)
  row_name <- c(list_names_murnau[i])
  for (j in c(1, 2, 3, 4, 5, 6)){
    row_value <- c(formatC(temp["diff"][[1]][j], digits = 3), 
                   formatC(temp["lwr"][[1]][j], digits = 3),
                   formatC(temp["upr"][[1]][j], digits = 3),
                   formatC(temp["p.adj"][[1]][j], format = "e", digits = 2))
    df_results_glht_Murnau_AIS_endpoint[nrow(df_results_glht_Murnau_AIS_endpoint) + 1,"Bloodmarker"] <- row_name
    df_results_glht_Murnau_AIS_endpoint[nrow(df_results_glht_Murnau_AIS_endpoint),"Variable"] <- row.names(result)[j]
    df_results_glht_Murnau_AIS_endpoint[nrow(df_results_glht_Murnau_AIS_endpoint), 3:6] <- row_value
  }
}

write.csv(df_results_glht_Murnau_AIS_endpoint,'results_stats/df_results_glht_Murnau_AIS_endpoint.csv')

# ------------------------------------------------------------------------------------------------------------------
# # Sygen
# ------------------------------------------------------------------------------------------------------------------

# # # ANOVA tests
# df_test = df_glu
# colnames(df_test)[which(names(df_test) == "level_T6")] <- "level"
# df_test$multiple_trauma <- ifelse(df_test$nb_inj>=2, 1, 0)
# anova_test(df_test)
# df_test$bloodvalue <- df_test$bloodvalue * 0.0555
# anova_test(df_test)

anova_test(df_alb)

# Create table to store results
df_results_ANOVA_Sygen <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Bloodmarker", "Variable", 'Chisq', "df", "pvalue")
colnames(df_results_ANOVA_Sygen) <- x

list_variable = c('', 'time', 'AIS', 'age', 'level', 'sex', 'multiple_trauma', 'time:AIS')

for (i in 1:length(list_names_sygen)){ # go through each haematological marker
  df = list_sygen[[i]] # select corresponding data frame
  colnames(df)[which(names(df) == "level_T6")] <- "level"
  df$multiple_trauma <- ifelse(df$nb_inj>=2, 1, 0)
  #df$time = as.factor(df$time)
  temp <- anova_test(df)
  row_name <- c(list_names_sygen[i])
  for (j in c(2, 3, 4, 5, 6, 7, 8)){
    # row_value <- c(formatC(temp[["Df"]][j]), 
    #                formatC(temp[["F value"]][j], digits = 3),
    #                formatC(temp[["Pr(>F)"]][j], format = "e", digits = 2))
    row_value <- c(formatC(temp[["Chisq"]][j]), # extract degree of freedom
                   formatC(temp[["Df"]][j], digits = 3), # extract F-value with 3 digits
                   formatC(temp[["Pr(>Chisq)"]][j], format = "e", digits = 2)) # extract p-value with format X.XXE-X
    df_results_ANOVA_Sygen[nrow(df_results_ANOVA_Sygen) + 1,"Bloodmarker"] <- row_name
    df_results_ANOVA_Sygen[nrow(df_results_ANOVA_Sygen),"Variable"] <- list_variable[j]
    df_results_ANOVA_Sygen[nrow(df_results_ANOVA_Sygen), 3:5] <- row_value
  }
}

write.csv(df_results_ANOVA_Sygen,'results_stats/df_results_ANOVA_Sygen_interactions_type3_RE_review_bintrauma.csv')

# ------------------------------------------------------------------------------------------------------------------

# # # Pairwise comparison of severity grades in explaining blood values
pairwise_test(df_alb, "Albumin")
colnames(df_alb)[which(names(df_alb) == "level_T6")] <- "level"

# Create table to store results
df_results_glht_Sygen_AIS <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Bloodmarker", "Comparison", "Effect", "p-value", "adjusted p-value")
colnames(df_results_glht_Sygen_AIS) <- x

for (i in 1:length(list_names_sygen)){ # go through each haematological marker
  df = list_sygen[[i]] # select corresponding data frame
  colnames(df)[which(names(df) == "level_T6")] <- "level"
  df$AIS <- as.factor(df$AIS)
  df$multiple_trauma <- ifelse(df$nb_inj>=2, 1, 0)
  temp <- pairwise_test(df, list_names_sygen[i])
  row_name <- c(list_names_sygen[i])
  for (j in 1:length(names(temp$test$coefficients))){
    row_value <- c(names(temp$test$coefficients)[j],
                   round(temp$test$coefficients[[j]], digits = 3), 
                   formatC(temp$test$pfunction()[[j]], format = "e", digits = 2),
                   formatC(temp$test$pvalues[[j]], format = "e", digits = 2))
    #df_results_glht_Sygen_AIS[i+j,"Bloodmarker"] <- row_name
    #df_results_glht_Sygen_AIS[i+j,2:5] <- row_value
    #df_results_glht_Sygen_AIS <- rbind(df_results_glht_Sygen_AIS, )
    df_results_glht_Sygen_AIS[nrow(df_results_glht_Sygen_AIS) + 1, "Bloodmarker"] <- row_name
    df_results_glht_Sygen_AIS[nrow(df_results_glht_Sygen_AIS), 2:5] <- row_value
  }
}

write.csv(df_results_glht_Sygen_AIS,'results_stats/df_results_glht_Sygen_AIS_review_bintrauma.csv')

# ------------------------------------------------------------------------------------------------------------------

# # # Slope analysis - comparing bloodvalues for different AIS grades over time
# slope(df_amylase)

# Create table to store results
df_results_slope_sygen <- data.frame(matrix(ncol = 7, nrow = 0))
x <- c("Bloodmarker", "Contrast", "Estimate", "SE", "df","t.ratio", "p- value")
colnames(df_results_slope_sygen) <- x

for (i in 1:length(list_names_sygen)){ # go through each haematological marker
  df = list_sygen[[i]] # select corresponding data frame
  #df$AIS <- as.factor(df$AIS)
  temp <- slope(df)
  row_name <- c(list_names_sygen[i])
  for (j in 1:length(levels(temp$contrast))){
    row_value <- c(levels(temp$contrast)[j],
                   round(temp$estimate[[j]], digits = 3), 
                   round(temp$SE[[j]], digits = 3),
                   round(temp$df[[j]], digits = 3),
                   round(temp[5][[1]][j], digits = 3),
                   formatC(temp$p.value[[j]], format = "e", digits = 2))
    df_results_slope_sygen[nrow(df_results_slope_sygen) + 1, "Bloodmarker"] <- row_name
    df_results_slope_sygen[nrow(df_results_slope_sygen), 2:7] <- row_value
  }
}

write.csv(df_results_slope_sygen,'results_stats/df_results_slope_sygen.csv')

# ------------------------------------------------------------------------------------------------------------------

# # # Pairwise comparison of severity grades in explaining blood values
try <- anova_test_1time(df_subset)

# Create table to store results
df_results_glht_Sygen_AIS_baseline <- data.frame(matrix(ncol = 6, nrow = 0))
x <- c("Bloodmarker", "Variable", "diff", "lwr", "upr", "p.adj")
colnames(df_results_glht_Sygen_AIS_baseline) <- x

for (i in 1:length(list_names_sygen)){ # go through each haematological marker
  df = list_sygen[[i]] # select corresponding data frame
  df_subset = subset(df, time == 0)
  df_subset$AIS <- as.factor(df_subset$AIS)
  temp_fct <- anova_test_1time(df_subset)
  temp <- data.frame(temp_fct$AIS)
  row_name <- c(list_names_sygen[i])
  for (j in c(1, 2, 3, 4, 5, 6)){
    row_value <- c(formatC(temp["diff"][[1]][j], digits = 3), 
                   formatC(temp["lwr"][[1]][j], digits = 3),
                   formatC(temp["upr"][[1]][j], digits = 3),
                   formatC(temp["p.adj"][[1]][j], format = "e", digits = 2))
    df_results_glht_Sygen_AIS_baseline[nrow(df_results_glht_Sygen_AIS_baseline) + 1,"Bloodmarker"] <- row_name
    df_results_glht_Sygen_AIS_baseline[nrow(df_results_glht_Sygen_AIS_baseline),"Variable"] <- row.names(result)[j]
    df_results_glht_Sygen_AIS_baseline[nrow(df_results_glht_Sygen_AIS_baseline), 3:6] <- row_value
  }
}

write.csv(df_results_glht_Sygen_AIS_baseline,'results_stats/df_results_glht_Sygen_AIS_baseline.csv')


# ------------------------------------------------------------------------------------------------------------------
# # Combined
# ------------------------------------------------------------------------------------------------------------------

list_murnau_combined = list(df_amylase, df_ap, df_calcium, df_creatinin, df_ery, df_bilirubin, df_gesamteiweiss,
                            df_glucose, df_harnstoff, df_hb, df_hk, df_kalium, df_leuco_nl,
                            df_mchc, df_mcv_m, df_natrium, df_ptt, df_thrombo, df_got, df_gpt)
list_sygen_combined = list(df_amy, df_alk, df_cab, df_bc9, df_rbc, df_bua, df_alb, df_glu, df_bun, df_hgb, df_hct,
                           df_kb, df_wbc, df_chc, df_mcv, df_nab, df_ap1, df_plt, df_sgo, df_sgp)
list_names_combined = c("Amylase", "Alkaline phosphatase", "Calcium", "Creatinin", "Erythrocytes", "Bilirubin",
                        "Proteins/Albumin", "Glucose", "Blood urea","Hemoglobin", "Hematocrit", "Potassium",
                        "Leucocytes", "MCHC", "MCV", "Sodium", "Prothrombin time", "Thrombocytes", "ASAT", "ALAT")

df_results_ANOVA_combined <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("Bloodmarker", "df", "Chisq", "pvalue")
colnames(df_results_ANOVA_combined) <- x

for (i in 1:length(list_names_combined)){
  df_murnau = list_murnau_combined[[i]] # select df of interest for Murnau data
  df_murnau = add_age_murnau(df_murnau)
  colnames(df_murnau)[which(names(df_murnau) == "level_T6")] <- "level"
  df_sygen = list_sygen_combined[[i]] # slect df of interest for Sygen data
  df_sygen$multiple_trauma <- ifelse(df_sygen$nb_inj>=2, 1, 0)
  df_sygen$nb_inj <- NULL
  df_murnau$time <- as.factor(df_murnau$time) # make time a factor
  df_sygen$time <- as.factor(df_sygen$time)
  df_sygen$AIS <- as.character(df_sygen$AIS)
  colnames(df_sygen)[which(names(df_sygen) == "level_T6")] <- "level"
  #df_sygen$AIS <- as.factor(df_sygen$AIS) # make AIS a factor
  
  # add a column identifying the database where the row come from
  for (j in 1:dim(df_murnau)[1]){
    df_murnau[j, "database"] <- 'Murnau'
  }
  for (j in 1:dim(df_sygen)[1]){
    df_sygen[j, "database"] <- 'Sygen'
  }
  
  # Unify the values for AIS grades
  df_sygen$AIS[df_sygen$AIS=="AIS A"]<-"A"
  df_sygen$AIS[df_sygen$AIS=="AIS B"]<-"B"
  df_sygen$AIS[df_sygen$AIS=="AIS C"]<-"C"
  df_sygen$AIS[df_sygen$AIS=="AIS D"]<-"D"
  # Unify the values for sex
  df_sygen$sex[df_sygen$sex==2]<-"m"
  df_sygen$sex[df_sygen$sex==1]<-"w"
  
  df_sygen$AIS <- as.factor(df_sygen$AIS)
  df_sygen$time <- as.factor(df_sygen$time)
  df_sygen$level <- as.factor(df_sygen$level)
  df_sygen$sex <- as.factor(df_sygen$sex)
  df_murnau$AIS <- as.factor(df_murnau$AIS)
  df_murnau$time <- as.factor(df_murnau$time)
  df_murnau$level <- as.factor(df_murnau$level)
  df_murnau$sex <- as.factor(df_murnau$sex)
  
  # combine both dfs
  df_combined = rbind(df_murnau, df_sygen)
  
  temp <- anova_combined(df_combined) # perform ANOVA analysis
  
  row_name <- c(list_names_combined[i]) # select haematological marker name
  
  row_value <- c(formatC(temp[["Df"]][3]), # report degree of freedom
                 formatC(temp[["Chisq"]][3], digits = 3), # report f value
                 formatC(temp[["Pr(>Chisq)"]][3], format = "e", digits = 2)) # report p- value
  
  df_results_ANOVA_combined[i,"Bloodmarker"] <- row_name
  df_results_ANOVA_combined[i,2:4] <- row_value
}

write.csv(df_results_ANOVA_combined,'results_stats/df_results_ANOVA_combined_type3_review_bintrauma.csv')


# ------------------------------------------------------------------------------------------------------------------

# Create heatmaps summarising the pvalue from lmer/glht tests 

data_test <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("Blood markers", "Time", "pvalue")
colnames(data_test) <- x

for (i in 1:length(list_names_murnau)){
  df = list_murnau[[i]]
  df$time <- as.factor(df$time)
  temp <- pairwise_test(df, list_names_murnau[i])
  row_name <- c(list_names_murnau[i])
  for (j in 1:length(names(temp$test$coefficients))){
    row_value <- c(names(temp$test$coefficients)[j],
                   formatC(temp$test$pvalues[[j]]))
    data_test[nrow(data_test) + 1, "Blood markers"] <- row_name
    data_test[nrow(data_test), 2:3] <- row_value
  }
}
#data_test

# Make p-value numerics
data_test$pvalue <- as.numeric(as.character(data_test$pvalue))

# Select which p-values to write in the heatmaps (only significant ones)
for (i in 1:dim(data_test)[1]){
  if (data_test[i,"pvalue"] < 0.05){
    #data_test[i,"sign"] = formatC(data_test[i,"pvalue"], format = "e", digits = 2)
    data_test[i,"sign"] = '*'
  } else {
    data_test[i,"sign"] = ''
  }
}

# Reshape names of the variables to display it in a more readable fashion 
for (i in 1:dim(data_test)[1]){
  x <- data_test[i, 'Time']
  x <- strsplit(x, '[ - ]')
  x <- x[[1]][c(3,2,1)]
  try <- paste0(x[1], x[2], x[3])
  data_test[i, 'Times compared'] <- try
}

# Create heatmap
data_test$`Times compared` <- f(data_test$`Times compared`)
colors <- c("#D7191C", "#FDAE61", "#ABD9E9", "#2C7BB6")
ggplot(data = data_test, aes(x = `Times compared`, y = `Blood markers`)) + # x axis : to adapt (time or AIS grades)
  geom_tile(aes(fill = pvalue)) + # fill cells of heatmap according to p-value
  geom_text(aes(label=sign), color="white") + # add value in cells when p-value is significant
  scale_fill_gradientn(colors = colors, breaks=c(0.001, 0.25, 0.5, 0.75, 1)) + # using a continuous gradient for filling
  ggtitle("Statistical results from liner mixed model \n Time in murnau study") + # title : to adapt
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.text.x = element_text(angle = 90))

