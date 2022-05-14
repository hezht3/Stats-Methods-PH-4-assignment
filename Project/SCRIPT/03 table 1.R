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
    as.data.frame() %>% 
    select(hearloss, gender, trueage, ethnicity, edug, occupation, residence, coresidence, marital) %>%
    tbl_summary(by = "hearloss",
                label = list(gender ~ "Sex",
                             trueage ~ "Age at baseline",
                             ethnicity ~ "Ethnicity",
                             edug ~ "Education background",
                             occupation ~ "Main occupation before age 60",
                             residence ~ "Residential status",
                             coresidence ~ "Living arrangement",
                             marital ~ "Current marital status"),
                statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
                missing = "no") %>%
    add_overall() %>%
    bold_labels() %>%
    as_flex_table() %>%
    save_as_docx(path = "./OUTPUT/Table 1. baseline summary statistics.docx")
