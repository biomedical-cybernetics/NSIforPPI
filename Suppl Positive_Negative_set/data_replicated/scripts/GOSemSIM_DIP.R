# R Script Documentation
# Author: Ilyes Abdelhamid, 2022
# Description: Compute Wang Gene Ontology (GO) semantic similarity for BP and CC
# for all protein pairs (existing and missing) in the Yeast DIP network.

# Set working directory where the R script is

library(GOSemSim)
library(tidyverse)
library(readxl)
library(data.table)


name = "Yeast"
organism = "yeast"

df <- fread(file = paste0('../Uniprot_to_EntrezID_', name, '_DIP_net.txt'))
genes <- df$To


# Start the timer
start_time <- Sys.time()
d <- godata('org.Sc.sgd.db', ont="BP", computeIC=FALSE)
BP = mgeneSim(genes, semData = d, measure="Wang", verbose = T)
write.table(BP, file = paste('../', name,"_DIP_BP.txt", sep=""), quote = FALSE, sep = "\t")
# End the timer
end_time <- Sys.time()
# Calculate the elapsed time
elapsed_time <- end_time - start_time

cat("Elapsed Time for BP: ", elapsed_time, " minutes\n")

# Start the timer
start_time <- Sys.time()
d <- godata('org.Sc.sgd.db', ont="CC", computeIC=FALSE)
CC = mgeneSim(genes, semData = d, measure="Wang", verbose = T)
write.table(CC, file = paste('../',name,"_DIP_CC.txt", sep=""), quote = FALSE, sep = "\t")
# End the timer
end_time <- Sys.time()
# Calculate the elapsed time
elapsed_time <- end_time - start_time

cat("Elapsed Time for CC: ", elapsed_time, " minutes\n")


