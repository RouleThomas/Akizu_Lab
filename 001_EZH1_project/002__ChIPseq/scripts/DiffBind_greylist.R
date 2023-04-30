#! /usr/bin/env Rscript

library("DiffBind") 

# Load the raw counts
load("output/DiffBind/sample_count_2dN.RData")
sample_count_2dN = sample_count
rm(sample_count)

load("output/DiffBind/sample_count_all.RData")
sample_count_all = sample_count
rm(sample_count)

load("output/DiffBind/sample_count_ESC.RData")
sample_count_ESC = sample_count
rm(sample_count)

load("output/DiffBind/sample_count_HET.RData")
sample_count_HET = sample_count
rm(sample_count)

load("output/DiffBind/sample_count_KO.RData")
sample_count_KO = sample_count
rm(sample_count)

load("output/DiffBind/sample_count_NPC.RData")
sample_count_NPC = sample_count
rm(sample_count)

load("output/DiffBind/sample_count_WT.RData")
sample_count_WT = sample_count
rm(sample_count)

# Apply greylist

sample_count_2dN_greylist = dba.blacklist(sample_count_2dN, blacklist=FALSE, greylist=TRUE, cores=1)
sample_count_all_greylist = dba.blacklist(sample_count_all, blacklist=FALSE, greylist=TRUE, cores=1)
sample_count_ESC_greylist = dba.blacklist(sample_count_ESC, blacklist=FALSE, greylist=TRUE, cores=1)
sample_count_HET_greylist = dba.blacklist(sample_count_HET, blacklist=FALSE, greylist=TRUE, cores=1)
sample_count_KO_greylist = dba.blacklist(sample_count_KO, blacklist=FALSE, greylist=TRUE, cores=1)
sample_count_NPC_greylist = dba.blacklist(sample_count_NPC, blacklist=FALSE, greylist=TRUE, cores=1)
sample_count_WT_greylist = dba.blacklist(sample_count_WT, blacklist=FALSE, greylist=TRUE, cores=1)



# Save as R object
save(sample_count_2dN_greylist, file = "output/DiffBind/sample_count_2dN_greylist.RData")
save(sample_count_all_greylist, file = "output/DiffBind/sample_count_all_greylist.RData")
save(sample_count_ESC_greylist, file = "output/DiffBind/sample_count_ESC_greylist.RData")
save(sample_count_HET_greylist, file = "output/DiffBind/sample_count_HET_greylist.RData")
save(sample_count_KO_greylist, file = "output/DiffBind/sample_count_KO_greylist.RData")
save(sample_count_NPC_greylist, file = "output/DiffBind/sample_count_NPC_greylist.RData")
save(sample_count_WT_greylist, file = "output/DiffBind/sample_count_WT_greylist.RData")

