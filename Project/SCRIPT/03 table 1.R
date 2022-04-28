######################################### CLHLS biostats 624 project ########################################
# Aim: Table 1

require(tidyverse)
require(gtsummary)
require(flextable)

setwd("D:/OneDrive - Johns Hopkins/Course/140.624.71 - Statistical Methods in Public Health IV/assignments/Stats-Methods-PH-4-assignment/Project")


# --------------------------------------------------------------------------------------------------------- #
##################################### import cleaned covariates dataset #####################################
# --------------------------------------------------------------------------------------------------------- #

data_main <- read_rds("./DATA/CLEAN/hearloss_ana18_repeated_main (enroll).rds")


# --------------------------------------------------------------------------------------------------------- #
######################################### baseline characteristics ##########################################
# --------------------------------------------------------------------------------------------------------- #

data_main %>%
    filter(id_num == 1) %>%
    select(hearloss, gender, trueage, edug, residence, coresidence) %>%
    tbl_summary(by = "hearloss",
                label = list(gender ~ "Sex",
                             trueage ~ "Age at baseline",
                             edug ~ "Education background",
                             residence ~ "Residential status",
                             coresidence ~ "Living arrangement"),
                statistic = list(all_continuous() ~ "{mean} ({sd}) {median} ({p25}, {p75})"),
                missing = "no") %>%
    add_overall() %>%
    bold_labels() %>%
    as_flex_table() %>%
    save_as_docx(path = "./OUTPUT/Table 1. baseline summary statistics.docx")
