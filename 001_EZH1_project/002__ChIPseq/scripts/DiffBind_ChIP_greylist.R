#! /usr/bin/env Rscript

library("DiffBind") 

# Load the raw counts
load("output/DiffBind/sample_count.RData")

# Apply greylist

sample_dba_greylist = dba.blacklist(sample_count, blacklist=FALSE, greylist=TRUE, cores=1)

# Save as R object
save(sample_dba_greylist, file = "output/DiffBind/sample_dba_greylist.RData")
