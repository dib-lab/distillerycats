library(ranger)
library(dplyr)
library(readr)
library(feather)

remotes::install_github("silkeszy/Pomona")
library(Pomona)
file.create(snakemake@output[['pomona']])
