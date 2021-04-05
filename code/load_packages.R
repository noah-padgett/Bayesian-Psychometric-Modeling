# ============================================= #
# script: load_packages.R
# Project: Code Guide BPM
# Author(s): R.N. Padgett
# ============================================= #
# Data Created: 2020-09-01
# Date Modified: 2020-09-04
# By: R. Noah Padgett
# ============================================= #
# Stems from Padgett's Independent Study
# ============================================= #
# Purpose:
# This R script is for loading all necessary
#   R packages
#
# No output - just loading packages into the
#   environment
# ============================================= #
# Set up directory and libraries
# rm(list=ls())
# list of packages
packages <- c("patchwork", "tidyr", "dplyr",
              "dtplyr", "data.table", "ggplot2",
              "R2jags","R2WinBUGS", "blavaan",
              "rstan", "brms", "bayesplot", "ggmcmc",
              "kableExtra", "extraDistr",
              "dagitty", "ggdag", "ggraph", "cowplot",
              "pdftools")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load packages
lapply(packages, library, character.only = TRUE)

w.d <- getwd()

