######################################### CLHLS biostats 624 project ########################################

# Aim: exploratory data cleaning and analysis on hearing loss and mortality


require(tidyverse)
require(data.table)
require(haven)
require(sjlabelled)
require(lubridate)

setwd("D:/OneDrive - Johns Hopkins/Course/140.624.71 - Statistical Methods in Public Health IV/assignments/Stats-Methods-PH-4-assignment/Project")


# --------------------------------------------------------------------------------------------------------- #
################################# import raw crs and long data across waves #################################
# --------------------------------------------------------------------------------------------------------- #

############################################ cross-sectional data ###########################################

filenames <- list.files("./DATA/CRS")
file <- substr(filenames, 11, 12)

crs <- list()
for(i in 1:length(filenames)) {
    crs[[file[i]]] <- read_dta(paste0("./DATA/CRS/", filenames[i]))
}

# wave indicator variable
for(wave in file) {
    crs[[wave]] <- crs[[wave]] %>%
        mutate(wave = wave) %>%
        mutate(wave = recode(wave, "98" = "1", "00" = "2", "02" = "3", "05" = "4",
                             "08" = "5", "11" = "6", "14" = "7"))
}

############################################# longitudinal data #############################################

filenames <- list.files("./DATA/LONG")
file <- substr(filenames, 4, 5)

long <- list()
for(i in 1:length(filenames)) {
    long[[file[i]]] <- read_dta(paste0("./DATA/LONG/", filenames[i]))
}

# extract new participants in current wave
for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    if(wave != "11" & wave != "08" & wave != "02") {
        long[[wave]] <- long[[wave]] %>% filter(id %% 100 == as.numeric(wave))
    } else if(wave == "11") {
        long[[wave]] <- long[[wave]] %>% filter(id %% 100 == 12)
    } else if(wave == "02") {
        long[[wave]] <- long[[wave]] %>% filter(id %% 100 == as.numeric(wave) | id == 21009300)   # missing in `00` wave
    } else {
        long[[wave]] <- long[[wave]] %>% filter(id %% 100 == 8 | id %% 100 == 9)
    }
}

# wave indicator variable
for(wave in file) {
    long[[wave]] <- long[[wave]] %>%
        mutate(wave = wave) %>%
        mutate(wave = recode(wave, "98" = "1", "00" = "2", "02" = "3", "05" = "4",
                             "08" = "5", "11" = "6", "14" = "7"))
}


# --------------------------------------------------------------------------------------------------------- #
############################## time-fixed variables, study exit date and status #############################
# --------------------------------------------------------------------------------------------------------- #

############### demographics: sex, ethnicity, education, occupation before age 60 at baseline ###############

for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    long[[wave]] <- long[[wave]] %>% 
        # gender
        mutate(gender = ifelse(a1 != 1 & a1 != 2, NA, a1)) %>% 
        mutate(gender = recode_factor(gender, "1" = "male", "2" = "female")) %>% 
        # ethnicity
        mutate(ethnicity = case_when(a2 == 1 ~ "1",
                                     a2 == 9 ~ as.character(NA),
                                     is.na(a2) == TRUE ~ as.character(NA),
                                     TRUE ~ "0")) %>%
        mutate(ethnicity = recode_factor(ethnicity, "1" = "han", "0" = "others")) %>% 
        # education
        mutate(edu = ifelse(f1 == 88 | f1 == 99, NA, f1)) %>%
        mutate(edug = case_when(edu == 0 ~ "1",
                                edu > 0 & edu <= 5 ~ "2",
                                edu > 5 ~ "3")) %>%
        mutate(edug = recode_factor(edug, "1" = "none",
                                    "2" = "primary school",
                                    "3" = "middle school or higher")) %>% 
        mutate_at(vars(edug), funs(setattr(., "label", "education category"))) %>%
        # occupation
        mutate(occupation = case_when(f2 == 9 ~ as.character(NA),
                                      f2 == 0 | f2 == 1 ~ "0",
                                      is.na(f2) == TRUE ~ as.character(NA),
                                      TRUE ~ "1")) %>%
        mutate(occupation = recode_factor(occupation, "1" = "mannual", "0" = "non-mannual"))
}

############################################# birth year, month #############################################
# change month > 12 to 12 (wave 98, 1 changed)
# assume day = 15 for all participants

for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    long[[wave]] <- long[[wave]] %>% 
        mutate(v_bthmon = ifelse(v_bthmon > 12, 12, v_bthmon)) %>% 
        # birthdate
        mutate(bthdate = paste(v_bthyr, v_bthmon, 15, sep = "-")) %>% 
        mutate(bthdate = ymd(bthdate))
}   # 2 missing in wave 98, 127 missing in wave 11

########################################## study exit date, status ##########################################

# [event date]: `dthdate`
for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    long[[wave]] <- long[[wave]] %>% 
        mutate(dthdate = stataXml::fromStataTime(dthdate, '%td'))
}

# [LOFU date]: `lostdate`, assuming lost to F/U at mid-point between adjacent two interview waves
for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    long[[wave]] <- long[[wave]] %>% 
        mutate(lostdate = stataXml::fromStataTime(lostdate, '%td'))
}

# [administrative censoring]: `interview2018`, interview date in wave 2018
for(wave in c("98", "00", "02", "05", "08", "11")) {
    long[[wave]] <- long[[wave]] %>% 
        mutate(interview2018 = stataXml::fromStataTime(interview2018, '%td'))
}
long[["14"]] <- long[["14"]] %>% rename("interview2018" = "in18")

# [event status]: `censor_18`

# [LOFU]: `lost_18`


# --------------------------------------------------------------------------------------------------------- #
############################# time-varying variables, study exit date and status ############################
# --------------------------------------------------------------------------------------------------------- #

#################################### height, weight, for constructing BMI ###################################

# wave 98
crs[["98"]] <- crs[["98"]] %>%
    # g12: weight
    mutate(weight = g12) %>%
    mutate(weight = ifelse(weight == 999, NA, weight)) %>%
    # g81: acromion process styloideus
    mutate(armlength = g81) %>%
    mutate(armlength = ifelse(armlength == 999, NA, armlength)) %>%
    # g82: right knee to the floor
    mutate(kneelength = g82) %>%
    mutate(kneelength = ifelse(kneelength == 999, NA, kneelength))

# wave 00
crs[["00"]] <- crs[["00"]] %>%
    # g10: weight
    mutate(weight = g10) %>%
    mutate(weight = ifelse(weight == 999, NA, weight))

# wave 02
crs[["02"]] <- crs[["02"]] %>%
    # g101: weight
    mutate(weight = g101) %>%
    mutate(weight = ifelse(weight == 999, NA, weight)) %>%
    # g102a: length from wrist to shoulder
    mutate(armlength = g102a) %>%
    mutate(armlength = ifelse(armlength == 999, NA, armlength)) %>%
    # g102b: length from kneel to floor
    mutate(kneelength = g102b) %>%
    mutate(kneelength = ifelse(kneelength == 999, NA, kneelength))

# wave 05
crs[["05"]] <- crs[["05"]] %>%
    # g101: weight
    mutate(weight = g101) %>%
    mutate(weight = ifelse(weight == 999, NA, weight)) %>%
    # g102: stature when the elder was young
    mutate(youngheight = g102) %>%
    mutate(youngheight = ifelse(youngheight == 999, NA, youngheight)) %>%
    mutate_at(vars(youngheight), funs(setattr(., "label", "stature when the elder was young")))

# wave 08
crs[["08"]] <- crs[["08"]] %>%
    # g101: weight
    mutate(weight = g101) %>%
    mutate(weight = ifelse(weight == 888 | weight == 999, NA, weight)) %>%
    # g1021: directly measured height of the interviewee
    mutate(meaheight = g1021) %>%
    mutate(meaheight = ifelse(meaheight == 888 | meaheight == 999, NA, meaheight)) %>%
    mutate_at(vars(meaheight), funs(setattr(., "label", "directly measured height of the interviewee"))) %>%
    # g122: height measured from the top of the right arm bone to the top of the right shoulder
    mutate(armlength = g122) %>%
    mutate(armlength = ifelse(armlength == 888 | armlength == 999, NA, armlength)) %>%
    # g123: height measured from the right knee joint to the ground
    mutate(kneelength = g123) %>%
    mutate(kneelength = ifelse(kneelength == 888 | kneelength == 999, NA, kneelength))

# wave 11, 14
for(wave in c("11", "14")) {
    crs[[wave]] <- crs[[wave]] %>%
        # g101: weight
        mutate(weight = g101) %>%
        mutate(weight = ifelse(weight == 888 | weight == 999, NA, weight)) %>%
        # g1021: directly measured height of the interviewee
        mutate(meaheight = g1021) %>%
        mutate(meaheight = ifelse(meaheight == 888 | meaheight == 999, NA, meaheight)) %>%
        mutate_at(vars(meaheight), funs(setattr(., "label", "directly measured height of the interviewee"))) %>%
        # g122: height measured from the top of the right arm bone to the top of the right shoulder
        mutate(armlength = g122) %>%
        mutate(armlength = ifelse(armlength == 88 | armlength == 99, NA, armlength)) %>%
        # g123: height measured from the right knee joint to the ground
        mutate(kneelength = g123) %>%
        mutate(kneelength = ifelse(kneelength == 88 | kneelength == 99, NA, kneelength))
}

############################### demographics: residence, coresidence, marital ###############################

# wave 98
crs[["98"]] <- crs[["98"]] %>%
    # residence
    mutate(residence = recode_factor(residenc, "2" = "rural", "1" = "urban (city or town)")) %>% 
    # co-residence
    mutate(coresidence = ifelse(a51 == 9, NA, a51)) %>%
    mutate(coresidence = recode_factor(coresidence, "1" = "with household member(s)",
                                       "2" = "alone",
                                       "3" = "in an institution")) %>% 
    # marital
    mutate(marital = case_when(f41 == 1 ~ "1",
                               f41 == 9 ~ as.character(NA),
                               is.na(f41) == TRUE ~ as.character(NA),
                               TRUE ~ "2")) %>%
    mutate(marital = recode_factor(marital, "2" = "others", "1" = "married"))

# wave 00, 02, 05, 08
for(wave in c("00", "02", "05", "08")) {
    crs[[wave]] <- crs[[wave]] %>%
        # residence
        mutate(residence = case_when(residenc == 3 ~ 2,
                                     residenc == 1 | residenc == 2 ~ 1)) %>% 
        mutate(residence = recode_factor(residence, "2" = "rural", "1" = "urban (city or town)")) %>% 
        # co-residence
        mutate(coresidence = ifelse(a51 == 9, NA, a51)) %>%
        mutate(coresidence = recode_factor(coresidence, "1" = "with household member(s)",
                                           "2" = "alone",
                                           "3" = "in an institution")) %>% 
        # marital
        mutate(marital = case_when(f41 == 1 ~ "1",
                                   f41 == 9 ~ as.character(NA),
                                   is.na(f41) == TRUE ~ as.character(NA),
                                   TRUE ~ "2")) %>%
        mutate(marital = recode_factor(marital, "2" = "others", "1" = "married"))
}

# wave 11, 14
for(wave in c("11", "14")) {
    crs[[wave]] <- crs[[wave]] %>%
        # residence
        mutate(residence = case_when(residenc == 3 ~ 2,
                                     residenc == 1 | residenc == 2 ~ 1)) %>% 
        mutate(residence = recode_factor(residence, "2" = "rural", "1" = "urban (city or town)")) %>% 
        # co-residence
        mutate(coresidence = ifelse(a51 == 8 | a51 == 9, NA, a51)) %>%
        mutate(coresidence = recode_factor(coresidence, "1" = "with household member(s)",
                                           "2" = "alone",
                                           "3" = "in an institution")) %>% 
        # marital
        mutate(marital = case_when(f41 == 1 ~ "1",
                                   f41 == 8 | f41 == 9 ~ as.character(NA),
                                   is.na(f41) ~ as.character(NA),
                                   TRUE ~ "2")) %>%
        mutate(marital = recode_factor(marital, "2" = "others", "1" = "married"))
}

############### demographics: sex, ethnicity, education, occupation before age 60 at baseline ###############
# time-varying measurements of time-fixed demographics [for missing values]

# wave 98
crs[["98"]] <- crs[["98"]] %>%
    # gender
    mutate(gender_r = ifelse(a1 != 1 & a1 != 2, NA, a1)) %>% 
    mutate(gender_r = recode_factor(gender_r, "1" = "male", "2" = "female")) %>% 
    # ethnicity
    mutate(ethnicity_r = case_when(a2 == 1 ~ "1",
                                   a2 == 9 ~ as.character(NA),
                                   is.na(a2) == TRUE ~ as.character(NA),
                                   TRUE ~ "0")) %>%
    mutate(ethnicity_r = recode_factor(ethnicity_r, "1" = "han", "0" = "others")) %>% 
    # education
    # 00 wave: f1 == -2 means "already asked in 98 survey"
    # However, no "-2" value was found in table(crs[["00"]]$f1)
    mutate(edu_r = ifelse(f1 == 88 | f1 == 99, NA, f1)) %>%
    mutate(edug_r = case_when(edu_r == 0 ~ "1",
                              edu_r > 0 & edu_r <= 5 ~ "2",
                              edu_r > 5 ~ "3")) %>%
    mutate(edug_r = recode_factor(edug_r, "1" = "none",
                                  "2" = "primary school",
                                  "3" = "middle school or higher")) %>% 
    mutate_at(vars(edug_r), funs(setattr(., "label", "education category"))) %>%
    # occupation before 60
    # Codebook of 08/11/14 waves shows f2 (occupation before 60) == -1 means "not applicable"
    # However, no "-1" value was found in table(crs[["crs##"]]$f2)
    mutate(occupation_r = case_when(f2 == 9 ~ as.character(NA),
                                    f2 == 0 | f2 == 1 ~ "0",
                                    is.na(f2) == TRUE ~ as.character(NA),
                                    TRUE ~ "1")) %>%
    mutate(occupation_r = recode_factor(occupation_r, "1" = "mannual", "0" = "non-mannual"))

# wave 00, 02, 05, 08
for(wave in c("00", "02", "05", "08")) {
    crs[[wave]] <- crs[[wave]] %>%
        # gender
        mutate(gender_r = ifelse(a1 != 1 & a1 != 2, NA, a1)) %>% 
        mutate(gender_r = recode_factor(gender_r, "1" = "male", "2" = "female")) %>% 
        # ethnicity
        mutate(ethnicity_r = case_when(a2 == 1 ~ "1",
                                       a2 == 9 ~ as.character(NA),
                                       is.na(a2) == TRUE ~ as.character(NA),
                                       TRUE ~ "0")) %>%
        mutate(ethnicity_r = recode_factor(ethnicity_r, "1" = "han", "0" = "others")) %>% 
        # education
        # 00 wave: f1 == -2 means "already asked in 98 survey"
        # However, no "-2" value was found in table(crs[["00"]]$f1)
        mutate(edu_r = ifelse(f1 == 88 | f1 == 99, NA, f1)) %>%
        mutate(edug_r = case_when(edu_r == 0 ~ "1",
                                  edu_r > 0 & edu_r <= 5 ~ "2",
                                  edu_r > 5 ~ "3")) %>%
        mutate(edug_r = recode_factor(edug_r, "1" = "none",
                                      "2" = "primary school",
                                      "3" = "middle school or higher")) %>% 
        mutate_at(vars(edug_r), funs(setattr(., "label", "education category"))) %>%
        # occupation before 60
        # Codebook of 08/11/14 waves shows f2 (occupation before 60) == -1 means "not applicable"
        # However, no "-1" value was found in table(crs[["crs##"]]$f2)
        mutate(occupation_r = case_when(f2 == 9 ~ as.character(NA),
                                        f2 == 0 | f2 == 1 ~ "0",
                                        is.na(f2) == TRUE ~ as.character(NA),
                                        TRUE ~ "1")) %>%
        mutate(occupation_r = recode_factor(occupation_r, "1" = "mannual", "0" = "non-mannual"))
}

# resolve `a2` and `f2` missing issue using longitudinal dataset
crs[["14"]] <- crs[["14"]] %>% 
    mutate(a2 = NULL,
           f2 = NULL) %>% 
    left_join(long[["14"]] %>% select(id, a2, f2), by = "id")

# wave 11, 14
for(wave in c("11", "14")) {
    crs[[wave]] <- crs[[wave]] %>%
        # gender
        mutate(gender_r = ifelse(a1 != 1 & a1 != 2, NA, a1)) %>% 
        mutate(gender_r = recode_factor(gender_r, "1" = "male", "2" = "female")) %>%
        # ethnicity
        mutate(ethnicity_r = case_when(a2 == 1 ~ "1",
                                       a2 == 9 ~ as.character(NA),
                                       is.na(a2) ~ as.character(NA),
                                       TRUE ~ "0")) %>%
        mutate(ethnicity_r = recode_factor(ethnicity_r, "1" = "han", "0" = "others")) %>% 
        # education
        # 00 wave: f1 == -2 means "already asked in 98 survey"
        # However, no "-2" value was found in table(crs[["00"]]$f1)
        mutate(edu_r = ifelse(f1 == 88 | f1 == 99, NA, f1)) %>%
        mutate(edug_r = case_when(edu_r == 0 ~ "1",
                                  edu_r > 0 & edu_r <= 5 ~ "2",
                                  edu_r > 5 ~ "3")) %>%
        mutate(edug_r = recode_factor(edug_r, "1" = "none",
                                      "2" = "primary school",
                                      "3" = "middle school or higher")) %>% 
        mutate_at(vars(edug_r), funs(setattr(., "label", "education category"))) %>%
        # occupation before 60
        # Codebook of 08/11/14 waves shows f2 (occupation before 60) == -1 means "not applicable"
        # However, no "-1" value was found in table(crs[["crs##"]]$f2)
        mutate(occupation_r = case_when(f2 == 9 ~ as.character(NA),
                                        f2 == 0 | f2 == 1 ~ "0",
                                        is.na(f2) ~ as.character(NA),
                                        TRUE ~ "1")) %>%
        mutate(occupation_r = recode_factor(occupation_r, "1" = "mannual", "0" = "non-mannual"))
}

################################################## smoking ##################################################

# wave 98
crs[["98"]] <- crs[["98"]] %>% 
    mutate(d63 = ifelse(d63 == -1 | d63 == 888 | d63 == 999, NA, d63),
           d64 = ifelse(d64 == -1 | d64 == 888 | d64 == 999, NA, d64),
           d65 = ifelse(d65 == -1 | d65 == 88 | d65 == 99, NA, d65)) %>% 
    # smoke years
    mutate(smk_year = d64 - d63) %>% 
    mutate_at(vars(smk_year), funs(setattr(., "label", "smoking year"))) %>%
    # smoke category
    mutate(smkl = case_when(d61 == 2 & d62 == 2 ~ "1",             # never smoking
                            d61 != 1 & d62 == 1 ~ "2",             # former smoking
                            d61 == 1 & d65 < 20 ~ "3",             # current smoking & times/day < 20
                            d61 == 1 & d65 >= 20 ~ "4")) %>%       # current smoking & times/day >= 20
    mutate(smkl = recode_factor(smkl, "1" = "never smoke", "2" = "former smoker",
                                "3" = "light current smoker", "4" = "heavy current smoker")) %>% 
    mutate_at(vars(smkl), funs(setattr(., "label", "smoker category")))

# wave 00, 02, 05, 08, 11, 14
for(wave in c("00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>% 
        mutate(d73 = ifelse(d73 == -1 | d73 == 888 | d73 == 999, NA, d73),
               d74 = ifelse(d74 == -1 | d74 == 888 | d74 == 999, NA, d74),
               d75 = ifelse(d75 == -1 | d75 == 88 | d75 == 99, NA, d75)) %>%
        # smoke years
        mutate(smk_year = d74 - d73) %>%
        mutate_at(vars(smk_year), funs(setattr(., "label", "smoking year"))) %>%
        # smoke category
        mutate(smkl = case_when(d71 == 2 & d72 == 2 ~ "1",               # never smoking
                                d71 != 1 & d72 == 1 ~ "2",               # former smoking
                                d71 == 1 & d75 < 20 ~ "3",               # current smoking & d65_wt<20
                                d71 == 1 & d75 >= 20 ~ "4")) %>%         # current smoking & d65_wt>=20
        mutate(smkl = recode_factor(smkl, "1" = "never smoke", "2" = "former smoker",
                                    "3" = "light current smoker", "4" = "heavy current smoker")) %>% 
        mutate_at(vars(smkl), funs(setattr(., "label", "smoker category")))
}

################################################## drinking #################################################

# The unit of wine is liang (Chinese unit) 1 liang = 50 grams
# wave 98: liquor-45.5%; wine-12%; rice wine-15%
# else wave: heavy spirit liquor-53%; light spirit liquor-38%; wine-12%; rice wine-15%; beer-4%; others-24.4%

# wave 98
crs[["98"]] <- crs[["98"]] %>%   # d75 = drinkcatg, d76 = drinkvol
    mutate(d76 = ifelse(d76 == -1 | d76 == 88 | d76 == 99, NA, d76)) %>%
    # alcohol
    mutate(alcohol = case_when(d75 == 1 ~ d76 * 50 * 0.455,            # Very strong liquor
                               d75 == 2 ~ d76 * 50 * 0.12,             # Wine
                               d75 == 3 ~ d76 * 50 * 0.15)) %>%        # Rice wine
    # dril
    mutate(dril = case_when(d71 == 2 & d72 == 2 ~ "1",                                # never drinker
                            d71 != 1 & d72 == 1 ~ "2",                                # former drinker
                            gender_r == "male" & d71 == 1 & alcohol <= 25 ~ "3",        # light current drinker
                            gender_r == "female" & d71 == 1 & alcohol <= 15 ~ "3",      # light current drinker
                            gender_r == "male" & d71 == 1 & alcohol > 25 ~ "4",         # heavy current drinker
                            gender_r == "female" & d71 == 1 & alcohol > 15 ~ "4")) %>%  # heavy current drinker
    mutate(dril = recode_factor(dril, "1" = "never drinker", "2" = "former drinker",
                                "3" = "current & light drinker", "4" = "current & heavy drinker")) %>% 
    mutate_at(vars(dril), funs(setattr(., "label", "drinker category")))

# wave 00, 02, 05, 08, 11, 14
for(wave in c("00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>%   # d85 = drinkcatg, d86 = drinkvol
        mutate(d86 = ifelse(d86 == -1 | d86 == 88 | d86 == 99, NA, d86)) %>%
        # alcohol
        mutate(alcohol = case_when(d85 == 1 ~ d86 * 50 * 0.53,             # Very strong liquor
                                   d85 == 2 ~ d86 * 50 * 0.38,             # Not very strong liquor
                                   d85 == 3 ~ d86 * 50 * 0.12,             # Wine
                                   d85 == 4 ~ d86 * 50 * 0.15,             # Rice wine
                                   d85 == 5 ~ d86 * 50 * 0.04,             # Beer
                                   d85 == 6 ~ d86 * 50 * 0.244)) %>%       # Others
        # dril
        mutate(dril = case_when(d81 == 2 & d82 == 2 ~ "1",                                # never drinker
                                d81 != 1 & d82 == 1 ~ "2",                                # former drinker
                                gender_r == "male" & d81 == 1 & alcohol <= 25 ~ "3",        # light current drinker
                                gender_r == "female" & d81 == 1 & alcohol <= 15 ~ "3",      # light current drinker
                                gender_r == "male" & d81 == 1 & alcohol > 25 ~ "4",         # heavy current drinker
                                gender_r == "female" & d81 == 1 & alcohol > 15 ~ "4")) %>%  # heavy current drinker
        mutate(dril = recode_factor(dril, "1" = "never drinker", "2" = "former drinker",
                                    "3" = "current & light drinker", "4" = "current & heavy drinker")) %>% 
        mutate_at(vars(dril), funs(setattr(., "label", "drinker category")))
}

######################################### regular physical activity #########################################

# wave 98
crs[["98"]] <- crs[["98"]] %>% 
    mutate(d83 = ifelse(d83 == -1 | d83 == 888 | d83 == 999, NA, d83)) %>% 
    mutate(pa = case_when(d81 == 1 & d83 < 50 ~ "1",
                          d81 == 1 & d83 >= 50 ~ "2",
                          d81 != 1 & d82 == 1 ~ "3",
                          d81 == 2 & d82 == 2 ~ "4")) %>%
    mutate(pa = recode_factor(pa, "1" = "current & start age < 50", "2" = "current & start age >= 50",
                              "3" = "former", "4" = "never")) %>% 
    mutate_at(vars(pa), funs(setattr(., "label", "physical activity")))

# wave 00, 02, 05, 08, 11, 14
for(wave in c("00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>%
        mutate(d93 = ifelse(d93 == -1 | d93 == 888 | d93 == 999, NA, d93)) %>%
        mutate(pa = case_when(d91 == 1 & d93 < 50 ~ "1",
                              d91 == 1 & d93 >= 50 ~ "2",
                              d91 != 1 & d92 == 1 ~ "3",
                              d91 == 2 & d92 == 2 ~ "4")) %>%
        mutate(pa = recode_factor(pa, "1" = "current & start age < 50", "2" = "current & start age >= 50",
                                  "3" = "former", "4" = "never")) %>% 
        mutate_at(vars(pa), funs(setattr(., "label", "physical activity")))
}

#################################################### diet ###################################################

# wave 98, 00, 02, 05
for(wave in c("98", "00", "02", "05")) {
    crs[[wave]] <- crs[[wave]] %>% 
        # fruit, vegetable
        mutate(across(c("d31", "d32"),
                      ~ case_when(.x == 9 ~ as.numeric(NA),
                                  .x == 1 | .x == 2 ~ 1,
                                  .x == 3 ~ 2,
                                  .x == 4 ~ 3),
                      .names = "{.col}_new")) %>% 
        # fish, egg, bean, salt vegetable, sugar, tea, garlic
        mutate(across(c("d4fish2", "d4egg2", "d4bean2", "d4veg2", "d4suga2", "d4tea2", "d4garl2"),
                      ~ .x,
                      .names = "{.col}_new")) %>% 
        # add value label for all diet variables
        mutate(across(ends_with("_new"),
                      ~ recode_factor(.x,
                                      "1" = "everyday or excepte winter",
                                      "2" = "occasionally",
                                      "3" = "rarely or never"))) %>% 
        # rename variables
        rename("fruit" = "d31_new", "veg" = "d32_new", "fish" = "d4fish2_new", "egg" = "d4egg2_new",
               "bean" = "d4bean2_new", "saltveg" = "d4veg2_new", "sugar" = "d4suga2_new",
               "tea" = "d4tea2_new", "garlic" = "d4garl2_new")
}

 
####################################### self-reported disease history #######################################

# wave 98
crs[["98"]] <- crs[["98"]] %>% 
    mutate(across(c("g17a1", "g17b1", "g17c1", "g17d1", "g17e1", "g17f1", "g17i1", "g17k1", "g17l1", "g17m1"),
                  ~ case_when(.x == 1 & get(glue::glue(paste(substr({cur_column()}, 1, 4), "2", sep = ""))) == 1 ~ 1,
                              .x == 1 & get(glue::glue(paste(substr({cur_column()}, 1, 4), "2", sep = ""))) != 1 ~ 2,
                              .x == 2 ~ 3,
                              TRUE ~ as.numeric(NA)),
                  .names = "{.col}_new")) %>% 
    mutate(across(ends_with("_new"),
                  ~ recode_factor(.x, "3" = "no disease", "2" = "disease", "1" = "disease & serious diasbility"))) %>% 
    rename("hypertension" = "g17a1_new", "diabetes" = "g17b1_new", "heartdisea" = "g17c1_new", "strokecvd" = "g17d1_new",
           "copd" = "g17e1_new", "tb" = "g17f1_new", "cancer" = "g17i1_new", "ulcer" = "g17k1_new", "parkinson" = "g17l1_new",
           "bedsore" = "g17m1_new")

# wave 00, 02, 05, 08, 11, 14
for(wave in c("00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>% 
        mutate(across(c("g15a1", "g15b1", "g15c1", "g15d1", "g15e1", "g15f1", "g15i1", "g15k1", "g15l1", "g15m1"),
                      ~ case_when(.x == 1 & get(glue::glue(paste(substr({cur_column()}, 1, 4), "3", sep = ""))) == 1 ~ 1,
                                  .x == 1 & get(glue::glue(paste(substr({cur_column()}, 1, 4), "3", sep = ""))) != 1 ~ 2,
                                  .x == 2 ~ 3,
                                  TRUE ~ as.numeric(NA)),
                      .names = "{.col}_new")) %>% 
        mutate(across(ends_with("_new"),
                      ~ recode_factor(.x, "3" = "no disease", "2" = "disease", "1" = "disease & serious diasbility"))) %>% 
        rename("hypertension" = "g15a1_new", "diabetes" = "g15b1_new", "heartdisea" = "g15c1_new", "strokecvd" = "g15d1_new",
               "copd" = "g15e1_new", "tb" = "g15f1_new", "cancer" = "g15i1_new", "ulcer" = "g15k1_new", "parkinson" = "g15l1_new",
               "bedsore" = "g15m1_new")
}

########################################## activity of daily living #########################################

for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>%
        mutate(across(c("e1", "e2", "e3", "e4", "e5", "e6"),
                      ~ case_when(.x == 1 ~ 0,
                                  .x == 2 | .x == 3 ~ 1,
                                  TRUE ~ as.numeric(NA)),
                      .names = "{.col}_new")) %>% 
        # mutate(across(ends_with("new"),
        #               ~ recode_factor(.x, "0" = "do not need help", "1" = "need help"))) %>% 
        rename("bathing" = "e1_new", "dressing" = "e2_new", "toileting" = "e3_new",
               "transferring" = "e4_new", "continence" = "e5_new", "feeding" = "e6_new")
}

############################################ self-reported health ###########################################

for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>% 
        mutate(srhealth = ifelse(b12 == 8 | b12 == 9, NA, b12)) %>% 
        mutate(srhealth = recode_factor(srhealth, "1" = "very good", "2" = "good", "3" = "so so", "4" = "bad", "5" = "very bad")) %>% 
        mutate_at(vars(srhealth), funs(setattr(., "label", "self-reported health")))
}

################################################# cognition #################################################

# orientation section
for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>% 
        mutate(across(c("c11", "c12", "c13", "c14", "c15"),
                      ~ case_when(.x == 8 | .x == 9 ~ as.numeric(NA),
                                  .x == 1 ~ 1,
                                  is.na(.x) ~ as.numeric(NA),
                                  TRUE ~ 0),
                      .names = "{.col}_new")) %>% 
        # mutate(across(ends_with("new"),
        #               ~ recode_factor(.x, "1" = "correct", "0" = "unable to answer or wrong"))) %>% 
        rename("time_orientation1" = "c11_new", "time_orientation2" = "c12_new", "time_orientation3" = "c13_new",
               "time_orientation4" = "c14_new", "place_orientation" = "c15_new")
}

# naming foods
for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>% 
        mutate(namefo = case_when((c16 >= 7 & c16 < 88) | (c16 > 88 & c16 < 99) | c16 > 99 ~ 7, # c16: # of kinds of food named in 1 minute
                                  c16 == 88 ~ 0,                                                # unable to answer
                                  c16 == 99 ~ as.numeric(NA),                                   # missing
                                  is.na(c16) == TRUE ~ as.numeric(NA),                          # missing
                                  TRUE ~ as.numeric(c16))) %>% 
        # mutate(namefo = recode_factor(namefo, "0" = "unable to answer", "1" = "1", "2" = "2", "3" = "3", "4" = "4",
        #                                       "5" = "5", "6" = "6", "7" = ">=7")) %>% 
        mutate_at(vars(namefo), funs(setattr(., "label", "# of kinds of food named in 1 minute")))
}

# registration calculation etc.
for(wave in c("98", "00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>% 
        mutate(across(c("c21a", "c21b", "c21c", "c31a", "c31b", "c31c", "c31d", "c31e", "c32", "c41a", "c41b", "c41c", "c51a", "c51b",
                        "c52", "c53a", "c53b", "c53c"),
                      ~ case_when(.x == 8 | .x == 9 ~ 0,
                                  is.na(.x) == TRUE ~ 0,
                                  TRUE ~ as.numeric(.x)),
                      .names = "{.col}_new")) %>% 
        # mutate(across(ends_with("new"),
        #               ~ recode_factor(.x, "0" = "wrong", "1" = "correct"))) %>% 
        rename("registration1" = "c21a_new", "registration2" = "c21b_new", "registration3" = "c21c_new",
               "calculation1" = "c31a_new", "calculation2" = "c31b_new", "calculation3" = "c31c_new",
               "calculation4" = "c31d_new", "calculation5" = "c31e_new",
               "copyf" = "c32_new", "delayed_recall1" = "c41a_new", "delayed_recall2" = "c41b_new", "delayed_recall3" = "c41c_new",
               "naming_objects1" = "c51a_new", "naming_objects2" = "c51b_new", "repeating_sentence" = "c52_new",
               "listening_obeying1" = "c53a_new", "listening_obeying2" = "c53b_new", "listening_obeying3" = "c53c_new")
}

################################################ hearing loss ###############################################

# based on interviewer's answer of "Was the interviewee able to hear?"
# wave 98
crs[["98"]] <- crs[["98"]] %>% 
    mutate(hearloss = case_when(h1a == 1 ~ "yes, without hearing aid",
                                h1a == 2 ~ "yes, but needs hearing aid",
                                h1a == 3 ~ "partly, despite hearing aid",
                                h1a == 4 ~ "no",
                                h1a == 9 ~ as.character(NA))) %>% 
    mutate(hearloss = factor(hearloss,
                             levels = c("yes, without hearing aid",
                                        "yes, but needs hearing aid",
                                        "partly, despite hearing aid",
                                        "no")))

# wave 00, 02, 05, 08, 11, 14
for(wave in c("00", "02", "05", "08", "11", "14")) {
    crs[[wave]] <- crs[[wave]] %>% 
        mutate(hearloss = case_when(h1 == 1 ~ "yes, without hearing aid",
                                    h1 == 2 ~ "yes, but needs hearing aid",
                                    h1 == 3 ~ "partly, despite hearing aid",
                                    h1 == 4 ~ "no",
                                    h1 == 8 | h1 == 9 ~ as.character(NA))) %>% 
        mutate(hearloss = factor(hearloss,
                                 levels = c("yes, without hearing aid",
                                            "yes, but needs hearing aid",
                                            "partly, despite hearing aid",
                                            "no")))
}

################################################ survey date ################################################

# crs98: year9899, month98, day98
# crs00: 2000, month00, day00
# crs02: 2002, month02, day02
# crs05: 2005, monthin, dayin
# crs08: 2008, monthin, dayin
# crs11: yearin, monthin, dayin
# crs14: yearin, monthin, dayin

# Replacement of missing values rule:
# a. If only the interview day is missing, then the day is assumed to be 15th
# b. If both month and day are missing and the year isn’t missing, or only the month is missing, the month/day is assumed to be that of the mid-point between the earliest interview date and the latest interview date of that year.
# c. No interview year is missing.

# Modify input mistakes:
# a. Change day 29/max of Feb to 28 for years 99, 01, 02, 03, 05, 06, 07, 09, 10, 11, 13, 14 (non-leap year)
# b. Change day 30/max of Feb to 29 for years 00, 04, 08, 12 (leap year)
# c. Change day 31 to 30 for months 4, 6, 9, 11


# wave 98

# crs[["98"]] %>%
#     filter(month98 != 99) %>%
#     filter(day98 != 99) %>%
#     mutate(survey_date = mdy(paste(month98, day98, year9899, sep = "-"))) %>%
#     summarise(median_survey_date = median(survey_date, na.rm = TRUE)) # 1998-05-04

crs[["98"]] <- crs[["98"]] %>%
    mutate(day98 = ifelse(day98 == 99, NA, day98)) %>% 
    mutate(month98 = ifelse(month98 == 99, NA, month98)) %>% 
    mutate(day98 = ifelse(is.na(day98) == TRUE & is.na(month98) == FALSE, 15, day98),
           day98 = ifelse(is.na(day98) == TRUE & is.na(month98) == TRUE, 4, day98),
           month98 = ifelse(is.na(month98) == TRUE, 5, month98)) %>% # no missing month for year9899 == 1999
    mutate(day98 = ifelse(day98 >= 29 & month98 == 2, 28, day98),
           day98 = ifelse(day98 >= 31 & (month98 == 4 | month98 == 6 | month98 == 9 | month98 == 11), 30, day98)) %>%
    mutate(survey_date = mdy(paste(month98, day98, year9899, sep = "-")))


# wave 00

# crs[["00"]] %>%
#     filter(month00 != 99) %>%
#     filter(day00 != 99) %>%
#     mutate(survey_date = mdy(paste(month00, day00, 2000, sep = "-"))) %>%
#     summarise(median_survey_date = median(survey_date, na.rm = TRUE)) # 2000-07-16

# no missing value for month or day
crs[["00"]] <- crs[["00"]] %>%
    mutate(day00 = ifelse(day00 >= 29 & month00 == 2, 28, day00),
           day00 = ifelse(day00 >= 31 & (month00 == 4 | month00 == 6 | month00 == 9 | month00 == 11), 30, day00)) %>%
    mutate(survey_date = mdy(paste(month00, day00, 2000, sep = "-")))


# wave 02

# crs[["02"]] %>%
#     filter(month02 != 99) %>%
#     filter(day02 != 99) %>%
#     mutate(survey_date = mdy(paste(month02, day02, 2002, sep = "-"))) %>%
#     summarise(median_survey_date = median(survey_date, na.rm = TRUE)) # 2002-05-17

# no missing value for month
crs[["02"]] <- crs[["02"]] %>%
    mutate(day02 = ifelse(day02 == 99, NA, day02)) %>% 
    mutate(day02 = ifelse(is.na(day02) == TRUE, 15, day02)) %>%
    mutate(day02 = ifelse(day02 >= 29 & month02 == 2, 28, day02),
           day02 = ifelse(day02 >= 31 & (month02 == 4 | month02 == 6 | month02 == 9 | month02 == 11), 30, day02)) %>%
    mutate(survey_date = mdy(paste(month02, day02, 2002, sep = "-")))


# wave 05

# crs[["05"]] %>%
#     filter(monthin != 99) %>%
#     filter(dayin != 99) %>%
#     mutate(survey_date = mdy(paste(monthin, dayin, 2005, sep = "-"))) %>%
#     summarise(median_survey_date = median(survey_date, na.rm = TRUE)) # 2005-05-10

# no missing value for month or day
crs[["05"]] <- crs[["05"]] %>% 
    mutate(dayin = ifelse(dayin >= 29 & monthin == 2, 28, dayin),
           dayin = ifelse(dayin >= 31 & (monthin == 4 | monthin == 6 | monthin == 9 | monthin == 11), 30, dayin)) %>%
    mutate(survey_date = mdy(paste(monthin, dayin, 2005, sep = "-")))


# wave 08

# crs[["08"]] %>%
#     filter(monthin != 99) %>%
#     filter(dayin != 99) %>%
#     mutate(survey_date = mdy(paste(monthin, dayin, 2008, sep = "-"))) %>%
#     summarise(median_survey_date = median(survey_date, na.rm = TRUE)) # 2008-08-02

# no missing value for month or day
crs[["08"]] <- crs[["08"]] %>% 
    mutate(dayin = ifelse(dayin >= 29 & monthin == 2, 28, dayin),
           dayin = ifelse(dayin >= 31 & (monthin == 4 | monthin == 6 | monthin == 9 | monthin == 11), 30, dayin)) %>%
    mutate(survey_date = mdy(paste(monthin, dayin, 2008, sep = "-")))


# wave 11

# crs[["11"]] %>%
#     filter(monthin != 99) %>%
#     filter(dayin != 99) %>%
#     mutate(survey_date = mdy(paste(monthin, dayin, yearin, sep = "-"))) %>%
#     summarise(median_survey_date = median(survey_date, na.rm = TRUE)) # 2011-08-20

# no missing value for year, month or day
crs[["11"]] <- crs[["11"]] %>% 
    mutate(dayin = ifelse(dayin >= 29 & monthin == 2, 28, dayin),
           dayin = ifelse(dayin >= 31 & (monthin == 4 | monthin == 6 | monthin == 9 | monthin == 11), 30, dayin)) %>%
    mutate(survey_date = mdy(paste(monthin, dayin, yearin, sep = "-")))


# wave 14

# crs[["14"]] %>%
#     filter(monthin != 99) %>%
#     filter(dayin != 99) %>%
#     mutate(survey_date = mdy(paste(monthin, dayin, yearin, sep = "-"))) %>%
#     summarise(median_survey_date = median(survey_date, na.rm = TRUE)) # 2014-05-27

# no missing value for year, month or day
crs[["14"]] <- crs[["14"]] %>% 
    mutate(dayin = ifelse(dayin >= 29 & monthin == 2, 28, dayin),
           dayin = ifelse(dayin >= 31 & (monthin == 4 | monthin == 6 | monthin == 9 | monthin == 11), 30, dayin)) %>%
    mutate(survey_date = mdy(paste(monthin, dayin, yearin, sep = "-")))


# --------------------------------------------------------------------------------------------------------- #
################################### pool all waves into analytic data set ###################################
# --------------------------------------------------------------------------------------------------------- #

################################## longitudinal data: time-fixed variables ##################################

keep.var <- c("id", "wave", "gender", "ethnicity", "edu", "edug", "occupation", "bthdate", "dthdate",
              "lostdate", "interview2018", "censor_18", "lost_18")

data.long <- long[["98"]][, names(long[["98"]]) %in% keep.var]
data.long <- copy_labels(data.long, long[["98"]])
for(wave in c("00", "02", "05", "08", "11", "14")) {
    data.long <- bind_rows(data.long, long[[wave]][, names(long[[wave]]) %in% keep.var])
    data.long <- copy_labels(data.long, long[[wave]])
}

################################ cross-sectional data: time-varying variables ###############################

keep.var <- c("id", "wave", "survey_date", "bthdate_r",
              "smkl", "smk_year", "dril", "alcohol", "pa",
              "fruit", "veg", "fish", "egg", "bean", "saltveg", "sugar", "tea", "garlic",
              "weight", "armlength", "kneelength", "youngheight", "meaheight",
              "gender_r", "ethnicity_r", "edu_r", "edug_r", "occupation_r",
              "trueage", "residence", "coresidence", "marital",
              "hypertension", "diabetes", "heartdisea", "strokecvd", "copd", "tb", "cancer", "ulcer", "parkinson", "bedsore",
              "bathing", "dressing", "toileting", "transferring", "continence", "feeding",
              "srhealth",
              "time_orientation1", "time_orientation2", "time_orientation3", "time_orientation4", "place_orientation",
              "namefo", "registration1", "registration2", "registration3", "calculation1", "calculation2", "calculation3",
              "calculation4", "calculation5", "copyf", "delayed_recall1", "delayed_recall2", "delayed_recall3",
              "naming_objects1", "naming_objects2", "repeating_sentence", "listening_obeying1", "listening_obeying2",
              "listening_obeying3", "hearloss")

data.crs <- crs[["98"]][, names(crs[["98"]]) %in% keep.var]
data.crs <- copy_labels(data.crs, crs[["98"]])
for(wave in c("00", "02", "05", "08", "11", "14")) {
    data.crs <- bind_rows(data.crs, crs[[wave]][, names(crs[[wave]]) %in% keep.var])
    data.crs <- copy_labels(data.crs, crs[[wave]])
}

################################ merge time-fixed with time-varying variables ###############################

data <- data.crs %>% 
    left_join(data.long %>% select(- wave), by = "id")

############################################### missing number ##############################################

data$diet_miss <- apply(data[, c('fruit', 'veg', 'fish', 'egg', 'bean', 'saltveg', 'sugar', 'tea', 'garlic')],
                        1, function(x) length(which(is.na(x) == TRUE)))
attr(data$diet_miss, "label") <- "missing num for diet variables"

data$adl_miss <- apply(data[, c('bathing', 'dressing', 'toileting', 'transferring', 'continence', 'feeding')],
                       1, function(x) length(which(is.na(x) == TRUE)))
attr(data$adl_miss, "label") <- "missing num for activity of daily living variables"

data$ci_miss <- apply(data[, c('time_orientation1', 'time_orientation2', 'time_orientation3', 'time_orientation4',
                               'place_orientation', 'namefo', 'registration1', 'registration2', 'registration3',
                               'calculation1', 'calculation2', 'calculation3', 'calculation4', 'calculation5',
                               'copyf', 'delayed_recall1', 'delayed_recall2', 'delayed_recall3', 'naming_objects1',
                               'naming_objects2', 'repeating_sentence', 'listening_obeying1', 'listening_obeying2',
                               'listening_obeying3')],
                      1, function(x) length(which(is.na(x) == TRUE)))
attr(data$ci_miss, "label") <- "missing num for cognitive impairment variables"


# --------------------------------------------------------------------------------------------------------- #
######################################## generate composite variable ########################################
# --------------------------------------------------------------------------------------------------------- #

#################################################### diet ###################################################

data <- data %>% 
    mutate(across(c("fruit", "veg", "garlic", "fish", "egg", "bean", "tea"),
                  ~ case_when(.x == "everyday or excepte winter" ~ 2,
                              .x == "occasionally" ~ 1,
                              .x == "rarely or never" ~ 0),
                  .names = "{.col}_new")) %>% 
    mutate(across(c("sugar", "saltveg"),
                  ~ case_when(.x == "everyday or excepte winter" ~ 0,
                              .x == "occasionally" ~ 1,
                              .x == "rarely or never" ~ 2),
                  .names = "{.col}_new"))

data <- data %>% 
    mutate(she = rowSums(data %>% select(ends_with("_new")), na.rm = TRUE)) %>%
    mutate(she = ifelse(diet_miss > 2, NA, she)) %>% 
    mutate_at(vars(she), funs(setattr(., "label", "simplified healthy eating index"))) # based on final Diet paper, needs to be verified

#################################################### BMI ####################################################

data <- data %>% 
    mutate(height = case_when(wave == "5" | wave == "6" | wave == "7" ~ meaheight/100,
                              TRUE ~ youngheight/100)) %>% 
    mutate(bmi = weight/(height^2)) %>% 
    mutate(bmi = ifelse(bmi < 12 | bmi >= 40, NA, bmi)) %>%
    mutate(bmi_cat = case_when(bmi < 18.5 ~ 1,
                               bmi >= 18.5 & bmi < 24 ~ 2,
                               bmi >= 24 ~ 3)) %>% 
    mutate(bmi_cat = recode_factor(bmi_cat, "1" = "underweight", "2" = "normal weight", "3" = "overweight & obesity")) %>% 
    mutate_at(vars(bmi_cat), funs(setattr(., "label", "bmi category")))

########################################## activity of daily living #########################################

data <- data %>% 
    mutate(adl_sum = rowSums(data %>% select("bathing", "dressing", "toileting", "transferring", "continence", "feeding"), na.rm = TRUE)) %>% 
    mutate(adl_sum = ifelse(adl_miss > 1, NA, adl_sum)) %>%
    mutate_at(vars(adl_sum), funs(setattr(., "label", "activity of daily living_sum"))) %>%
    mutate(adl = case_when(adl_sum > 0 ~ 1,
                           adl_sum == 0 ~ 0)) %>% 
    mutate(adl = recode_factor(adl, "0" = "do not need help", "1" = "need help")) %>% 
    mutate_at(vars(adl), funs(setattr(., "label", "activity of daily living_sum")))

############################################ cognitive impairment ###########################################

data <- data %>% 
    mutate(ci = rowSums(data %>% 
                            select('time_orientation1', 'time_orientation2', 'time_orientation3', 'time_orientation4', 'place_orientation', 'namefo',
                                   'registration1', 'registration2', 'registration3', 'calculation1', 'calculation2', 'calculation3', 'calculation4',
                                   'calculation5', 'copyf', 'delayed_recall1', 'delayed_recall2', 'delayed_recall3', 'naming_objects1',
                                   'naming_objects2', 'repeating_sentence', 'listening_obeying1', 'listening_obeying2', 'listening_obeying3'),
                        na.rm = TRUE)) %>% 
    mutate_at(vars(ci), funs(setattr(., "label", "cognition impairment_sum")))

data <- data %>% 
    mutate(ci = ifelse(ci_miss > 3, NA, ci)) %>% 
    mutate(mmse = case_when(ci >= 25 ~ 1,
                            ci >= 18 & ci < 25 ~ 2,
                            ci >= 10 & ci < 18 ~ 3,
                            ci >= 0 & ci < 10 ~ 4)) %>% 
    mutate(mmse = recode_factor(mmse, "1" = "no impairment", "2" = "mild", "3" = "moderate", "4" = "severe")) %>% 
    mutate_at(vars(mmse), funs(setattr(., "label", "cognition impairment category")))

############################################# Composite exposure ############################################

# dichotomized healthy lifestyle

data <- data %>%
    # dichotomized smoking
    mutate(smkl_bi = case_when(smkl == "never smoke" ~ 1,
                               smkl %in% c("former smoker", "light current smoker", "heavy current smoker") ~ 0)) %>%
    mutate_at(vars(smkl_bi), funs(setattr(., "label", "dichotomized main measure: smoke"))) %>%
    # dichotomized drinking
    mutate(dril_bi = case_when(dril == "never drinker" ~ 1,
                               dril %in% c("former drinker", "current & light drinker", "current & heavy drinker") ~ 0)) %>%
    mutate_at(vars(dril_bi), funs(setattr(., "label", "dichotomized main measure: drink"))) %>%
    # dichotomized physical activity
    mutate(pa_bi = case_when(pa == "current & start age < 50" | pa == "current & start age >= 50" ~ 1,
                             pa == "former" | pa == "never" ~ 0)) %>%
    mutate_at(vars(pa_bi), funs(setattr(., "label", "dichotomized main measure: physical activity"))) %>%
    # dichotomized diet
    mutate(diet_bi = case_when(she %in% c(11:18) ~ 1,
                               she %in% c(0:10) ~ 0)) %>%
    mutate_at(vars(diet_bi), funs(setattr(., "label", "dichotomized main measure: diet"))) %>% 
    # dichotomized bmi
    mutate(bmi_bi = case_when(bmi_cat == "underweight" | bmi_cat == "overweight & obesity" ~ 0,
                              bmi_cat == "normal weight" ~ 1)) %>%
    mutate_at(vars(bmi_bi), funs(setattr(., "label", "dichotomized main measure: bmi")))

data <- data %>% 
    mutate(lifesty_bi = rowSums(data[, c("smkl_bi", "dril_bi", "pa_bi", "diet_bi")])) %>%
    mutate(lifesty_bi = factor(lifesty_bi, levels = c("0", "1", "2", "3", "4"))) %>%
    mutate_at(vars(lifesty_bi), funs(setattr(., "label", "healthy lifestyle profile")))

data <- data %>%
    mutate(lifesty5_bi = rowSums(data[, c("smkl_bi", "dril_bi", "pa_bi", "diet_bi", "bmi_bi")])) %>%
    mutate(lifesty5_bi = factor(lifesty5_bi, levels = c("0", "1", "2", "3", "4", "5"))) %>%
    mutate_at(vars(lifesty5_bi), funs(setattr(., "label", "healthy lifestyle profile, with bmi")))

data <- data %>% 
    mutate(smkl_bi = recode_factor(smkl_bi, "0" = "unhealthy", "1" = "healthy")) %>% 
    mutate(dril_bi = recode_factor(dril_bi, "0" = "unhealthy", "1" = "healthy")) %>% 
    mutate(pa_bi = recode_factor(pa_bi, "0" = "unhealthy", "1" = "healthy")) %>% 
    mutate(diet_bi = recode_factor(diet_bi, "0" = "unhealthy", "1" = "healthy")) %>% 
    mutate(bmi_bi = recode_factor(bmi_bi, "0" = "unhealthy", "1" = "healthy"))

# quartile healthy lifestyle

data <- data %>%
    mutate(smkl_cat = case_when(smkl == "never smoke" ~ 3,
                                smkl == "former smoker" ~ 2,
                                smkl == "light current smoker" ~ 1,
                                smkl == "heavy current smoker" ~ 0)) %>%
    mutate_at(vars(smkl_cat), funs(setattr(., "label", "categorized main measure: smoke"))) %>%
    mutate(dril_cat = case_when(dril == "never drinker" ~ 3,
                                dril == "former drinker" ~ 2,
                                dril == "current & light drinker" ~ 1,
                                dril == "current & heavy drinker" ~ 0)) %>%
    mutate_at(vars(dril_cat), funs(setattr(., "label", "categorized main measure: drink"))) %>%
    mutate(pa_cat = case_when(pa == "current & start age < 50" ~ 3,
                              pa == "current & start age >= 50" ~ 2,
                              pa == "former" ~ 1,
                              pa == "never" ~ 0)) %>%
    mutate_at(vars(pa_cat), funs(setattr(., "label", "categorized main measure: physical activity"))) %>%
    mutate(diet_cat = case_when(she %in% c(13:18) ~ 3,
                                she %in% c(11:12) ~ 2,
                                she %in% c(9:10) ~ 1,
                                she %in% c(0:8) ~ 0)) %>% # Min: 0, 1st Q: 8, median: 10, 3rd Q: 12, Max: 18
    mutate_at(vars(diet_cat), funs(setattr(., "label", "categorized main measure: diet")))

data <- data %>%
    mutate(lifesty_cat = rowSums(data[, c("smkl_cat", "dril_cat", "pa_cat", "diet_cat")])) %>% 
    mutate(lifesty_cat = case_when(lifesty_cat %in% c(9:12) ~ "Q4",
                                   lifesty_cat == 8 ~ "Q3",
                                   lifesty_cat == 7 ~ "Q2",
                                   lifesty_cat %in% c(0:6) ~ "Q1")) %>%
    mutate(lifesty_cat = factor(lifesty_cat, levels = c("Q1", "Q2", "Q3", "Q4"))) %>%
    mutate_at(vars(lifesty_cat), funs(setattr(., "label", "quartile healthy lifestyle profile")))

data <- data %>% 
    mutate(smkl_cat = recode_factor(smkl_cat, "0" = "heavy current smoker", "1" = "light current smoker",
                                    "2" = "former smoker", "3" = "never smoke")) %>% 
    mutate(dril_cat = recode_factor(dril_cat, "0" = "current & heavy drinker", "1" = "current & light drinker",
                                    "2" = "former drinker", "3" = "never drinker")) %>% 
    mutate(pa_cat = recode_factor(pa_cat, "0" = "never", "1" = "former",
                                  "2" = "current & start age >= 50", "3" = "current & start age < 50")) %>% 
    mutate(diet_cat = recode_factor(diet_cat, "0" = "0-8", "1" = "9-10", "2" = "11-12", "3" = "13-18"))


# --------------------------------------------------------------------------------------------------------- #
############################### export data set for survival data preparation ###############################
# --------------------------------------------------------------------------------------------------------- #

data <- data %>% arrange(id, wave, survey_date)

write_dta(data, "./DATA/CLEAN/hearloss_ana18_repeated.dta")
saveRDS(data, "./DATA/CLEAN/hearloss_ana18_repeated.rds")