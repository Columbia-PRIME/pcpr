# File: queens.R
# Date: 07/15/2022
# Author: Lawrence Chillrud <lgc2139@cumc.columbia.edu>
# Description: Prepares queens speciated PM2.5 data for pcpr using raw data from
# EPA's AQS datamart.

#--------------------------#
####      CONTENTS      ####
#--------------------------#
# N. Notes
# 0. Package imports
# 1. Download data
# 2. Prepare data
# 3. Save data

#--------------------------#
####      N. NOTES      ####
#--------------------------#
# This script prepares the queens speciated PM2.5 data for the pcpr package
# using the raw daily concentrations from the EPA's AQS data mart (The Queens
# monitor's AQS Site ID is: 36-081-0124). This script relies on the the URLs
# defined in the "urls" vector below as raw inputs, generating data/queens.rda
# as output. Documentation for the final queens dataset can be found in R/data.R

#--------------------------#
#### 0. PACKAGE IMPORTS ####
#--------------------------#
library(dplyr)
library(stringr)
library(tidyr)

#--------------------------#
####  1. DOWNLOAD DATA  ####
#--------------------------#
# 1a. Define vector of URLs (years 2001-2021) for downloading data:
urls <- paste0("https://www3.epa.gov/cgi-bin/broker?_service=data&_program=dataprog.Daily.sas&check=site&debug=0&year=", 2001:2021, "&site=36-081-0124")

# 1b. Read the URLs in (this should take ~15-20min depending on OS!):
raw_data <- purrr::map_dfr(urls, readr::read_csv)

#--------------------------#
####  2. PREPARE DATA   ####
#--------------------------#
# 2a. Define speciated PM2.5 exposures of interest:
chems <- c(
  "aluminum", "ammonium ion", "arsenic", "barium", "bromine",
  "cadmium", "calcium", "chlorine", "chromium", "copper",
  "elemental carbon", "iron", "lead", "magnesium", "manganese",
  "nickel", "organic carbon", "potassium ion", "selenium", "silicon",
  "sodium", "sulfur", "titanium", "total nitrate", "vanadium", "zinc"
)

# 2b. Define abbreviations for the column names of the final cleaned data:
abbreviations <- c(
  "Al", "NH4", "As", "Ba", "Br", "Cd", "Ca", "Cl", "Cr", "Cu",
  "EC", "Fe", "Pb", "Mg", "Mn", "Ni", "OC", "K", "Se", "Si",
  "Na", "S", "Ti", "NO3", "V", "Zn"
)

# 2c. Extract daily (24h) PM2.5 species of interest, pivot to wider format:
clean_data <- raw_data %>%
  mutate(`Parameter Name` = str_to_lower(`Parameter Name`)) %>%
  filter(
    str_detect(`Parameter Name`, "pm2.5"),
    `Units of Measure` == "Micrograms/cubic meter (LC)",
    `Duration Description` == "24 HOUR"
  ) %>%
  select(`Parameter Name`, `Date (Local)`, `Arithmetic Mean`) %>%
  mutate(
    `Parameter Name` = ifelse(`Parameter Name` == "ec csn_rev unadjusted pm2.5 lc tot", "elemental carbon", `Parameter Name`),
    `Parameter Name` = ifelse(`Parameter Name` == "oc csn_rev unadjusted pm2.5 lc tot", "organic carbon", `Parameter Name`),
    `Parameter Name` = str_remove(`Parameter Name`, " pm2.5(.*)")
  ) %>%
  filter(`Parameter Name` %in% chems) %>%
  pivot_wider(names_from = `Parameter Name`, values_from = `Arithmetic Mean`, values_fill = NA) %>%
  arrange(`Date (Local)`) %>%
  select(`Date (Local)`, all_of(chems))

# 2d. Rename the columns of the data to the abbreviated chemical names:
queens <- setNames(clean_data, c("Date", abbreviations))

#--------------------------#
####    3. SAVE DATA    ####
#--------------------------#
# 3a. Save as a .csv file:
readr::write_csv(queens, here::here("data-raw/queens.csv"))

# 3b. Save as a .rda file:
save(queens, file = here::here("data/queens.rda"))

# 3c. For administrative use (i.e., for the pcpr package):
# usethis::use_data(queens, overwrite = TRUE)
