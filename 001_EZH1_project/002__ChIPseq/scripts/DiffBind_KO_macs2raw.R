#! /usr/bin/env Rscript

library("DiffBind") 

# Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba_KO = dba(sampleSheet=read.table("output/DiffBind/meta_sample_KO_macs2raw.txt", header = TRUE, sep = "\t"))

# Count
sample_count_KO = dba.count(sample_dba_KO, bParallel=TRUE) # bParallel=TRUE is to count using multiple processors

# Save the sample count as R object
save(sample_count_KO, file = "output/DiffBind/sample_count_KO_macs2raw.RData")


# Apply greylist


sample_count_KO_blackgreylist = dba.blacklist(sample_count_KO, blacklist=TRUE, greylist=TRUE, cores=1)


# Save as R object
save(sample_count_KO_blackgreylist, file = "output/DiffBind/sample_count_KO_macs2raw_greylist.RData")



