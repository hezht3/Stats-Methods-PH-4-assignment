######################################### CLHLS biostats 624 project ########################################
# Aim: Table 1

require(tidyverse)
require(survival)
require(survminer)
require(gtsummary)
require(flextable)

setwd("D:/OneDrive - Johns Hopkins/Course/140.624.71 - Statistical Methods in Public Health IV/assignments/Stats-Methods-PH-4-assignment/Project")


# --------------------------------------------------------------------------------------------------------- #
##################################### import cleaned covariates dataset #####################################
# --------------------------------------------------------------------------------------------------------- #

data_main <- read_rds("./DATA/CLEAN/hearloss_ana18_repeated_main (enroll).rds")


# --------------------------------------------------------------------------------------------------------- #
########################################### kaplan-meier estimates ##########################################
# --------------------------------------------------------------------------------------------------------- #

data_main %>% 
    group_by(id) %>% 
    filter(id_num == n()) %>% 
    ungroup() %>% 
    summarise(median = median(sur_time2),
              mean = mean(sur_time2),
              min = min(sur_time2),
              max = max(sur_time2),
              p25 = quantile(sur_time2, 0.25),
              p75 = quantile(sur_time2, 0.75))

data_main %>% 
    group_by(id) %>% 
    filter(id_num == n()) %>% 
    ungroup() %>% 
    summarise(event = sum(censor_18)/n(),
              lost = sum(lost_18, na.rm = TRUE)/n())

pyears(SurvObj ~ hearloss, data = data_main, scale = 1) %>% summary()

survfit(SurvObj ~ hearloss, data = data_main)

tiff("./OUTPUT/Figure. survival curve.tiff",
     width = 8000, height = 7230, pointsize = 12, res = 600)
ggsurvplot(survfit(SurvObj ~ hearloss, data = data_main), data = data_main,
           palette = "jama", linetype = 1, censor.size = 0.8, # change the width of the curve by this command
           xlab = "Duration of follow-up (years)", ylab = "Probability of survival",
           break.x.by = 4, xlim = c(0,20.5), conf.int = TRUE,
           legend.title = "Hearing status",
           legend.labs = c("Able to hear, without hearing aid", "Able to hear, with hearing aid",
                           "Partly able to hear, despite hearing aid", "Not able to hear"),
           risk.table = TRUE, cumevents = TRUE, fontsize = 6,
           risk.table.title = "No. at Risk", cumevents.title = "Cumulative No. of Deaths",
           risk.table.height = 0.15, cumevents.height = 0.15,
           risk.table.x.text = FALSE, cumevents.x.text = FALSE, risk.table.y.text = FALSE, cumevents.y.text = FALSE,
           font.main = c(15, "bold", "black"), font.x = c(15, "bold", "black"), font.y = c(15, "bold", "black"),
           font.tickslab = c(15, "bold", "black"), font.legend = c(15, "bold", "black"),
           tables.theme = theme_cleantable()) +
    guides(colour = guide_legend(nrow = 2))
dev.off()


# --------------------------------------------------------------------------------------------------------- #
####################################### cox proportional hazard model #######################################
# --------------------------------------------------------------------------------------------------------- #

crude.model <- coxph(SurvObj ~ hearloss, data = data_main)
adj1.model <- coxph(SurvObj ~ hearloss + gender + factor(trueage), data = data_main)
adj2.model <- coxph(SurvObj ~ hearloss + gender + factor(trueage)
                    + ethnicity + edug + occupation + residence + coresidence + marital, data = data_main)
adj3.model <- coxph(SurvObj ~ hearloss + gender + factor(trueage)
                    + ethnicity + edug + occupation + residence + coresidence + marital
                    + smkl_bi + dril_bi + pa_bi + diet_bi, data = data_main)
adj4.model <- coxph(SurvObj ~ hearloss + gender + factor(trueage)
                    + ethnicity + edug + occupation + residence + coresidence + marital
                    + smkl_bi + dril_bi + pa_bi + diet_bi
                    + hypertension + diabetes + heartdisea + strokecvd + copd + tb + cancer + ulcer + parkinson
                    + bedsore + bathing + dressing + toileting + transferring + continence + feeding
                    + factor(ci), data = data_main)

tbl_merge(list(crude.model, adj1.model, adj2.model, adj3.model, adj4.model) %>%
              map(~ .x %>%
                      tbl_regression(exponentiate = TRUE, include = c("hearloss")) %>%
                      bold_labels),
          tab_spanner = c("**Unadjusted**", "**Model 1**", "**Model 2**", "**Model 3**", "**Model 4**")) %>%
    as_flex_table() %>%
    save_as_docx(path = "./OUTPUT/Table 2. cox model.docx")

coxph(SurvObj ~ as.numeric(hearloss) + gender + factor(trueage)
      + ethnicity + edug + occupation + residence + coresidence + marital, data = data_main)   # trend test
