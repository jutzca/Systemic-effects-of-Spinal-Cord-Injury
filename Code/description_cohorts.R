# ------------------------------------------------------------------------------------------------------------------
# Description of included and excluded cohorts
#
# January 2021
# L. Bourguignon
# ------------------------------------------------------------------------------------------------------------------

# Load packages ----
library("readxl")
library(data.table)
library(dplyr)
library(plyr)
library(tidyr)

# Load data ----
df_sygen_raw <- read.csv('~/data/clinical_trial/JohnKramersProject_DATA_2019-10-07_0111.csv')
df_murnau_raw <- read.csv('~/data/HematologicalBiomark_DATA_2021-01-07_0729.csv')
df_murnau_demog <- read_excel('~/data/qry_CATHERINE_Demogr_60000.xlsx')
df_murnau_demog <- as.data.frame(df_murnau_demog)
df_murnau_age <- read.csv('~/df_age_updated_murnau.csv')
df_amylase <- read.csvl('~/data/df_amylase_murnau.csv')

# Sygen cohorts ----

df_sygen_subcol <- df_sygen_raw[c('ptid', 'age', 'sexcd', 'ais1', 'splvl', 'lower01', 'lower52')]
# Number of NA for age
sum(is.na(df_sygen_subcol$age))
# Number of NA for sex
sum(is.na(df_sygen_subcol$sexcd))
# Number of NA for AIS grade at baseline
sum(is.na(df_sygen_subcol$ais1))
# Number of NA for level of injury
sum(is.na(df_sygen_subcol$splvl))
# Number of NA for lems at week 1
sum(is.na(df_sygen_subcol$lower01))

# Remove patients with missing AIS grade at baseline
df_sygen_subcol <- subset(df_sygen_subcol, ais1 != '')
# Remove patients with NA for age, sex, AIS grade at baseline and/or level of injury
df_sygen_subcol_nona <- df_sygen_subcol %>% drop_na(age, sexcd, ais1, splvl)

# Count patients per sex in included cohort
table(df_sygen_subcol_nona$sexcd)
# Count patients per AIS grade in included cohort
table(df_sygen_subcol_nona$ais1)
# Mean age in included cohort
mean(df_sygen_subcol_nona$age)
# SD age in included cohort
sd(df_sygen_subcol_nona$age)
# Mean lems at baseline in included cohort
mean(df_sygen_subcol_nona$lower01, na.rm = TRUE)
# SD lems at baseline in included cohort
sd(df_sygen_subcol_nona$lower01, na.rm = TRUE)
# Mean lems at 1year post-trauma in included cohort
mean(df_sygen_subcol_nona$lower52, na.rm = TRUE)
# SD lems at 1year post-trauma in included cohort
sd(df_sygen_subcol_nona$lower52, na.rm = TRUE)
# Number of NA for lems at 1year post-trauma in included cohort
sum(is.na(df_sygen_subcol_nona$lower52))


# Murnau cohorts ----

df_murnau_demog <- setnames(df_murnau_demog, "Patientennummer", "patientennummer")
df_murnau_merge <- df_murnau_raw %>% inner_join(df_murnau_demog, by=c("patientennummer"))
df_murnau_merge2 <- df_murnau_merge %>% inner_join(df_murnau_age, by=c("patientennummer"))

df_amylase2 <- df_amylase[c('random_effect', 'AIS', 'level')]
df_amylase2 <- setnames(df_amylase2, "random_effect", "patientennummer")
df_amylase3 <- df_amylase2[!duplicated(df_amylase2),]
df_murnau_merge3 <- df_murnau_merge2 %>% inner_join(df_amylase3, by=c("patientennummer"))

df_murnau_subcol <- df_murnau_merge3[c('patientennummer', 'Sex', 'va_lems', 'c_lems', 'AIS', 'level', 'age')]
# Number of NA for age
sum(is.na(df_murnau_subcol$age))
# Number of NA for sex
sum(is.na(df_murnau_subcol$Sex))
# Count patients per AIS grade
table(df_murnau_subcol$AIS)
# Number of NA for level of injury
sum(is.na(df_murnau_subcol$level))
# Number of NA for lems at baseline
sum(is.na(df_murnau_subcol$va_lems))
# Number of NA for lems at 1year post-injury
sum(is.na(df_murnau_subcol$c_lems))

# Remove patients with missing or normal AIS grade at baseline
df_murnau_subcol <- subset(df_murnau_subcol, AIS != '')
df_murnau_subcol <- subset(df_murnau_subcol, AIS != 'E')
df_murnau_subcol <- subset(df_murnau_subcol, AIS != 'ND')
df_murnau_subcol <- subset(df_murnau_subcol, age != 'NA')
# Remove patients with NA for age, sex, AIS grade at baseline and/or level of injury
df_murnau_subcol_nona <- df_murnau_subcol %>% drop_na(Sex, va_lems, AIS, level)
df_murnau_subcol_nona$Sex[df_murnau_subcol_nona$Sex=='w'] <- 'f'
# Number of NA for lems at 1year post-injury in included cohort
sum(is.na(df_murnau_subcol_nona$c_lems))
# Count patients per sex in included cohort
table(df_murnau_subcol_nona$Sex)
# Count patients per AIS grade in included cohort
table(df_murnau_subcol_nona$AIS)
# Mean age in included cohort
mean(df_murnau_subcol_nona$age)
# SD age in included cohort
sd(df_murnau_subcol_nona$age)
# Mean lems at baseline in included cohort
mean(df_murnau_subcol_nona$va_lems)
# SD lems at baseline in included cohort
sd(df_murnau_subcol_nona$va_lems)

df_murnau_temp <- subset(df_murnau_subcol_nona, va_lems != '')
df_murnau_temp <- subset(df_murnau_temp, va_lems != 'ND')
df_murnau_temp <- df_murnau_temp %>% drop_na(va_lems)
# Mean lems at baseline in included cohort
mean(as.numeric(df_murnau_temp$va_lems), na.rm = TRUE)
# SD lems at baseline in included cohort
sd(as.numeric(df_murnau_temp$va_lems), na.rm = TRUE)
# Mean lems at 1year post-injury in included cohort
mean(as.numeric(df_murnau_temp$c_lems), na.rm = TRUE)
# SD lems at 1year post-injury in included cohort
sd(as.numeric(df_murnau_temp$c_lems), na.rm = TRUE)

patientennummer <- df_murnau_merge3$patientennummer [!df_murnau_merge3$patientennummer %in% df_murnau_subcol_nona$patientennummer ]
# Create dataframe for excluded cohort
excluded_murnau <- as.data.frame(patientennummer)
excluded_murnau <- excluded_murnau %>% inner_join(df_murnau_merge3, by=c("patientennummer"))
excluded_murnau <- excluded_murnau[c('patientennummer', 'Sex', 'va_lems', 'c_lems', 'AIS', 'level', 'age')]
# Number of NA for lems at 1year post-injury in excluded cohort
sum(is.na(excluded_murnau$c_lems))
# Number of NA for age in excluded cohort
sum(is.na(excluded_murnau$age))
# Count patients per sex in excluded cohort
table(excluded_murnau$Sex)
# Count patients per AIS grade in excluded cohort
table(excluded_murnau$AIS)
# Mean age in excluded cohort
mean(excluded_murnau$age, na.rm = T)
# SD age in excluded cohort
sd(excluded_murnau$age, na.rm = T)
# Mean lems at baseline in excluded cohort
mean(as.numeric(excluded_murnau$va_lems), na.rm = TRUE)
# SD lems at baseline in excluded cohort
sd(as.numeric(excluded_murnau$va_lems), na.rm = TRUE)


# Statistical tests ----

# Compare Sygen and Murnau in terms of sex proportions
M <- as.table(rbind(c(46, 143), c(193, 560)))
print(chisq.test(M))

# Compare Sygen and Murnau in terms of age distributions
t.test(as.numeric(df_murnau_subcol_nona$age), as.numeric(df_sygen_subcol_nona$age), alternative = c("two.sided"))

# Compare Sygen and Murnau in terms of AIS grades proportions
M <- as.table(rbind(c(81, 446), c(22, 77), c(26,149), c(110,31)))
print(chisq.test(M))

# Compare Sygen and Murnau in terms of lems at baseline distributions
t.test(as.numeric(df_murnau_subcol_nona$va_lems), as.numeric(df_sygen_subcol_nona$lower01), alternative = c("two.sided"), na.rm = TRUE)
# Compare Sygen and Murnau in terms of lems at 1year post-injury distributions
t.test(as.numeric(df_murnau_subcol_nona$c_lems), as.numeric(df_sygen_subcol_nona$lower52), alternative = c("two.sided"), na.rm = TRUE)

# Compare included and excluded cohorts from Murnau in terms of age distributions
t.test(as.numeric(df_murnau_subcol_nona$age), as.numeric(excluded_murnau$age), alternative = c("two.sided"))

# Compare included and excluded cohorts from Murnau in terms of sex proportions
M <- as.table(rbind(c(46, 42), c(193, 82)))
print(chisq.test(M))

# Compare included and excluded cohorts from Murnau in terms of AIS grades proportions
M <- as.table(rbind(c(94, 72), c(0, 21), c(0,5), c(0,26)))
print(chisq.test(M))

