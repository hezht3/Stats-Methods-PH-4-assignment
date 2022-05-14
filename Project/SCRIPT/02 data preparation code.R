######################################### CLHLS biostats 624 project ########################################
# Aim: Survival data setup

require(tidyverse)
require(haven)
require(lubridate)
require(gtsummary)
require(survival)

setwd("D:/OneDrive - Johns Hopkins/Course/140.624.71 - Statistical Methods in Public Health IV/assignments/Stats-Methods-PH-4-assignment/Project")


# --------------------------------------------------------------------------------------------------------- #
##################################### import cleaned covariates dataset #####################################
# --------------------------------------------------------------------------------------------------------- #

data <- read_rds("./DATA/CLEAN/hearloss_ana18_repeated.rds")


# --------------------------------------------------------------------------------------------------------- #
############### dealing with missing values: time-fixed variables, time-varying variables LOCF ##############
# --------------------------------------------------------------------------------------------------------- #
# Summary statistics on missing values in covariates, stratified by complete case exposure and eligibility
# Time-fixed variables: missing values filled with available measurements from cross-sectional surveys
# Time-varying variables: missing values filled with last observation carried forward

########################### time-fixed variables filling with repeated measurements #########################

Mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
}

# gender
data_fill <- data %>%
    mutate(gender = if_else(is.na(gender) == TRUE, as.character(gender_r), as.character(gender))) %>% 
    mutate(gender = recode_factor(gender, "1" = "male", "2" = "female"))

# ethnicity
data_fill <- data_fill %>%
    mutate(ethnicity = if_else(is.na(ethnicity) == TRUE, as.character(ethnicity_r), as.character(ethnicity))) %>% 
    group_by(id) %>% 
    mutate(ethnicity = Mode(ethnicity)) %>% 
    mutate(ethnicity = recode_factor(ethnicity, "1" = "han", "0" = "others"))

# education
data_fill <- data_fill %>%
    mutate(edug = if_else(is.na(edug) == TRUE, as.character(edug_r), as.character(edug))) %>% 
    group_by(id) %>% 
    mutate(edug = Mode(edug)) %>% 
    mutate(edug = recode_factor(edug, "1" = "none",
                                "2" = "primary school",
                                "3" = "middle school or higher"))

# occupataion
data_fill <- data_fill %>%
    mutate(occupation = if_else(is.na(occupation) == TRUE, as.character(occupation_r), as.character(occupation))) %>% 
    group_by(id) %>% 
    mutate(occupation = Mode(occupation)) %>% 
    mutate(occupation = recode_factor(occupation, "1" = "mannual", "0" = "non-mannual"))

#################### time-varying variables filling with last observation carried forward ###################

data_fill <- data_fill %>% 
    group_by(id) %>% 
    fill(hearloss,
         smkl_bi, dril_bi, pa_bi, diet_bi,
         residence, coresidence, marital,
         hypertension, diabetes, heartdisea, strokecvd, copd, tb, cancer, ulcer, parkinson, bedsore,
         bathing, dressing, toileting, transferring, continence, feeding, ci)


# --------------------------------------------------------------------------------------------------------- #
############### survival data setup: time origin, enrollment; time metric: length of follow-up ##############
# --------------------------------------------------------------------------------------------------------- #
# Event: death
# Time origin: study enrollment (first survey date after age >= 80 with complete measurements in exposure and covariates)
# Time metric: time since enrollment, duration of follow-up
# Study entry: same as time origin
# Study exit: death
# Censoring:
#     1. Lost to F/U defined as 1 day after last survey date
#     2. Administrative censoring defined as 1 day after last survey date (wave 18 survey date)

################################## for main adjustment model: sex, age, SES #################################

# Step 1: identify time origin
#   Step 1.1: drop event status problematic participants (`censor_18` == 1 & `dthdate` <= `survey_date`)
problematic_id <- data_fill %>% 
    filter(censor_18 == 1 & dthdate <= survey_date) %>% 
    distinct(id) %>% 
    pull(id)   # 183 participants, 389 person-periods

data_main <- data_fill %>% filter(!(id %in% problematic_id))

#   Step 1.2: drop participants with missing endpoint status
#      Case 1: censor_18 == NA & lost_18 == NA & interview2018 == NA
#      Case 2: censor_18 == 0 & lost_18 == NA & interview2018 == NA
data_main <- data_main %>% 
    filter(!(censor_18 == 0 & is.na(lost_18) == TRUE & is.na(interview2018) == TRUE))   # 62 participants, 86 person-periods

#   Step 1.3: exclude observations before study entry - less than 80 years old (based on `trueage`, not at risk)
data_main <- data_main %>% filter(trueage >= 80)   # 20175 person-periods

#   Step 1.4: exclude observations before study entry - first complete case in exposure and covariates
#   After this step, the first row of each participant is the identified time origin/study entry
data_main <- data_main %>%    # missing values indicator variable, 1 = no missing, 0 = missing
    mutate(na.indi = case_when(if_all(c(hearloss,
                                        smkl_bi, dril_bi, pa_bi, diet_bi,
                                        gender, trueage, ethnicity, edug, occupation, residence, coresidence, marital,
                                        hypertension, diabetes, heartdisea, strokecvd, copd, tb, cancer, ulcer, parkinson, bedsore,
                                        bathing, dressing, toileting, transferring, continence, feeding, ci),
                                      ~ is.na(.x) == FALSE) ~ 1,
                               TRUE ~ 0))

data_main <- data_main %>%    # cumulative sum of `na.cum`, observations at and after first complete case are positive numbers
    group_by(id) %>% 
    mutate(na.cum = cumsum(na.indi)) %>% 
    ungroup() %>% 
    filter(na.cum > 0)   # 1451 participants, 2722 person-periods
# 35836 participants, 62122 person-periods remained in the dataset

data_main <- data_main %>% 
    mutate(across(c("hearloss"),
                  ~ forcats::fct_explicit_na(.x)))   # code missing values in exposure as a separate level

# Step 2: generate entry time `sur_time1` for each person-period (based on `survey_date`)
data_main <- data_main %>% mutate(sur_time1 = survey_date)   # no NA values in `sur_time1`

# Step 3: generate exit time `sur_time2` for each person-period

#   Step 3.1: Examine missing study visits
problematic_idwave <- data_main %>% 
    mutate(wave = as.numeric(wave)) %>% 
    group_by(id) %>% 
    mutate(wavediff = wave - lag(wave)) %>% 
    ungroup() %>% 
    filter(wavediff > 1) %>% 
    mutate(idwave = paste(id, wave, sep = "_")) %>% 
    pull(idwave)   # 1 participant: `50001898` missed wave 2

#   Step 3.2: for person-periods prior to the last row for each participant
#      1. If the participant does not miss any study visit, exit time = entry time for next person-period
#      2. If the participant misses any study visit, exit time = mid-point survey date for next survey wave
#         (Assuming the participant continues to contribute person-time to the mid-point of next survey wave,
#         and does not contribute person-time till return survey date ("gaps"))
#      3. No missing values in covariates
data_main <- data_main %>% 
    arrange(id, sur_time1) %>% 
    mutate(sur_time1_r = sur_time1) %>%    # backup `sur_time1`
    mutate(idwave = paste(id, wave, sep = "_")) %>% 
    mutate(sur_time1_r = if_else(idwave %in% problematic_idwave,   # participants misses visits
                                 as.character("9999-12-31"),       # code entry time for first return visit as "9999-12-31"
                                 as.character(sur_time1_r))) %>% 
    group_by(id) %>% 
    mutate(sur_time2 = lead(sur_time1_r, 1)) %>%    # exit time = entry time for next person-period
    ungroup() %>% 
    mutate(sur_time2 = case_when(sur_time2 == "9999-12-31" & wave == "1" ~ "2000-07-16",
                                 sur_time2 == "9999-12-31" & wave == "2" ~ "2002-05-17",
                                 sur_time2 == "9999-12-31" & wave == "3" ~ "2005-05-10",
                                 sur_time2 == "9999-12-31" & wave == "4" ~ "2008-08-02",
                                 sur_time2 == "9999-12-31" & wave == "5" ~ "2011-08-20",
                                 sur_time2 == "9999-12-31" & wave == "6" ~ "2014-05-27",
                                 sur_time2 == "9999-12-31" & wave == "7" ~ "2018-08-01",
                                 TRUE ~ sur_time2)) %>%    # exit time = mid-point survey date for next survey wave
    mutate(sur_time2 = ymd(sur_time2)) %>% 
    mutate(sur_time1_r = NULL)

#   Step 3.3: for last person-periods (last row) for each participant
#      1. If the participant has event, exit time = death date `dthdate`
#      2. If the participant lost to F/U, exit time = last survey date + 1 day
#      3. Administrative censoring at last survey date + 1 day (wave 18 survey date + 1 day)
data_main <- data_main %>% 
    group_by(id) %>% 
    # person-periods prior to last row
    mutate(censor_18 = ifelse(is.na(sur_time2) == FALSE, 0, censor_18)) %>%
    mutate(sur_time2 = case_when(
        # event (failure)
        is.na(sur_time2) == TRUE & censor_18 == 1 ~ as.character(dthdate),
        # lost to F/U
        is.na(sur_time2) == TRUE & censor_18 == 0 & lost_18 == 1 ~ as.character(sur_time1 + 1),
        # administrative censoring
        is.na(sur_time2) == TRUE & censor_18 == 0 & is.na(lost_18) == TRUE ~ as.character(interview2018),
        # person-periods prior to last row
        is.na(sur_time2) == FALSE ~ as.character(sur_time2),
        TRUE ~ as.character(NA))) %>% 
    mutate(sur_time2 = ymd(sur_time2)) %>% 
    mutate(sur_time2 = as.numeric(sur_time2 - min(sur_time1))/365.25) %>% 
    mutate(sur_time1 = as.numeric(sur_time1 - min(sur_time1))/365.25) %>% 
    mutate(id_num = 1:n())

# Step 4: examine person-periods with missing values ("gaps"), correct measurment errors, examine gaps
#    Observations excluded are middle waves with missing values in exposure stats
data_main %>% filter(lifesty_bi == "(Missing)") %>% nrow()   # 2125 person-periods

#    Correct measurement errors in `trueage`
#    Correct observations with `trueage` in later wave lower than earlier wave
data_main <- data_main %>% 
    mutate(calage = floor((survey_date - bthdate)/365.25)) %>% 
    mutate(calage = as.numeric(calage))

data_main %>%
    ggplot(aes(x = trueage, y = calage)) +
    geom_point()

problematic_id <- data_main %>%
    group_by(id) %>%
    mutate(agediff = trueage - lag(trueage)) %>%
    ungroup() %>%
    filter(agediff <= 0) %>%
    pull(id)

data_main <- data_main %>% 
    mutate(trueage = ifelse(id %in% problematic_id, calage, trueage))

#   Examine "gaps"
data_main %>% 
    mutate(wave = as.numeric(wave)) %>% 
    group_by(id) %>% 
    mutate(wavediff = wave - lag(wave)) %>% 
    ungroup() %>% 
    summarise(`wave` = table(wavediff))   # including missing exposures: 756 participants with 1 gap, 46 participants with 2 gaps, 6 participants with 3 gaps

# Step 5: create risk sets
data_main <- data_main %>% 
    mutate(SurvObj = Surv(time = sur_time1, time2 = sur_time2, event = censor_18))

write_rds(data_main, file = "./DATA/CLEAN/hearloss_ana18_repeated_main (enroll).rds")
write_dta(data_main %>% select(- na.indi, - na.cum), "./DATA/CLEAN/hearloss_ana18_repeated_main (enroll).dta")


survfit(SurvObj ~ hearloss, data = data_main)
survminer::ggsurvplot(survfit(SurvObj ~ hearloss, data = data_main)[c(1:4)])

coxph(SurvObj ~ hearloss, data = data_main)


data_main <- data_main %>% 
    mutate(lifesty_3 = case_when(lifesty_bi == "0" | lifesty_bi == "1" ~ "0-1",
                                 lifesty_bi == "2" ~ "2",
                                 lifesty_bi == "3" | lifesty_bi == "4" ~ "3-4")) %>% 
    mutate(lifesty_3 = factor(lifesty_3, levels = c("0-1", "2", "3-4")))

coxph(SurvObj ~ hearloss * lifesty_3, data = data_main)

data_main <- data_main %>% 
    mutate(smkl_c = case_when(
        smkl_cat == "heavy current smoker" | smkl_cat == "light current smoker" ~ "current",
        smkl_cat == "former smoker" | smkl_cat == "never smoke" ~ "not current")) %>% 
    mutate(smkl_c = factor(smkl_c, levels = c("not current", "current")))
