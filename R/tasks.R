#!/usr/bin/env Rscript
# Prints the model x phenotype x sample-set task list (one per line) for the stage 3 and 4 job arrays.
source("R/config.R")
an <- load_config()$pipeline$analysis
tk <- CJ(model = an$models, phenotype = an$phenotypes, relatedness = an$relatedness, sorted = FALSE)
writeLines(paste(tk$model, tk$phenotype, tk$relatedness))
