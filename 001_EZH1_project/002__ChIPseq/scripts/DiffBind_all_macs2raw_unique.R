#! /usr/bin/env Rscript

library("DiffBind") 

# Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_all_macs2raw_unique.txt", header = TRUE, sep = "\t"))

# Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba, bParallel=TRUE) # bParallel=TRUE is to count using multiple processors

# Save the sample count as R object
save(sample_count, file = "output/DiffBind/count_all_macs2raw_unique.RData")



# Apply greylist
sample_count_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE, cores=1)


# Save as R object
save(sample_count_blackgreylist, file = "output/DiffBind/count_all_macs2raw_unique_blackgreylist.RData")


