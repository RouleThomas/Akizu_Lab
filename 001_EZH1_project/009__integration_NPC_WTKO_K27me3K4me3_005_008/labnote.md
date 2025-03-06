# Project goal

Dataset integration of 2 biol rep:
- `005__CutRun` NPC WT vs KO H3K27me3 and H3K4me3 (native)
- `008__CutRun` NPC WT vs KO H3K27me3 and H3K4me3 (FA)

# Pipeline integration

- copy the `*.unique.dupmark.sorted.bam` and `*.bai` files
- re -do masc2 peak calling 
- Calculate the MG1655 SF
- DiffBind_TMM SF calculation
- THOR to scale bigwig


# Copy and rename files


--> Copy `*.unique.dupmark.sorted.bam` and `*.bai` files and manually added `_008` or `_005` respectively

*NOTE: For 005 there was sample inversion for `NPC_KO_H3K4me3` (eg. NPC_KO_H3K27me1 is NPC_KO_H3K4me3); I made sure to cp the good file by checking the library size*


# macs2 peak calling

Let's just re-do this as there was file renaiming, to avoid fail with copy/paste..



**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# AB per AB
sbatch scripts/macs2_broad_H3K27me3.sh # 13878669 ok
sbatch scripts/macs2_broad_H3K4me3.sh # 13878889 ok

# pool
sbatch scripts/macs2_broad_H3K27me3_pool.sh # 15518679 ok
sbatch scripts/macs2_broad_H3K4me3_pool.sh # 15518691 ok
```



```bash
conda activate bowtie2 # for bedtools
sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive
sbatch scripts/macs2_raw_peak_signif_pool.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive


# quick command to print median size of peak within a bed
awk '{print $3-$2}' your_bed_file.bed | sort -n | awk 'BEGIN {c=0; sum=0;} {a[c++]=$1; sum+=$1;} END {if (c%2) print a[int(c/2)]; else print (a[c/2-1]+a[c/2])/2;}'
```

Then keep only the significant peaks (re-run the script to test different qvalue cutoff) and remove peaks overlapping with blacklist regions. MACS2 column9 output is -log10(qvalue) format so if we want 0.05; 
- q0.05: `q value = -log10(0.05) = 1.30103`
- q0.01 = 2
- q0.005 = 2.30103
- q0.001 = 3
- q0.0001 = 4
- q0.00001 = 5


**Optimal qvalue** according to IGV:
- NPC_H3K27me3: 1.30103 (**2.3 more true peaks**)
- NPC_H3K4me3: 1.30103 (**2.3 more true peaks**)

--> Let's go with 2.3




# Ecoli scaling factor
## Mapping E coli

--> Already done in `005__CutRun` and `008__CutRun`
----> Create new xls file with all these values `SpikeIn_MG1655_005_008.xlsx`

*NOTE: For 005 there was sample inversion for `NPC_KO_H3K4me3` (eg. NPC_KO_H3K27me1 is NPC_KO_H3K4me3); I made sure to use the correct E coli count*

Now calculate SF in R, as for histone SF:


```R
# package
library("tidyverse")
library("readxl")
library("ggpubr")

# SF H3K27me3
spikein <- read_excel("output/spikein/SpikeIn_MG1655_005_008.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "H3K27me3")
# Total reads per IP
spikein_H3K27me3_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_H3K27me3_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_H3K27me3_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)


# SF H3K4me3
spikein <- read_excel("output/spikein/SpikeIn_MG1655_005_008.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "H3K4me3")
# Total reads per IP
spikein_H3K4me3_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_H3K4me3_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_H3K4me3_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)


```

--> OK



# Spike in scaling
## With MG1655 spike in for PTM CutRun

--> Let's use MG1655 as default method for spike in normalization! Seems more accurate and can be used for all AB!

**Using our scaling factor, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/)) 


### Adjust library size with MG1655 scaling factor and apply normalization
Total number of reads is our library size (used `samtools flagstat` to double check) :

`samtools flagstat output/bowtie2/*.dupmark.sorted.bam` used to obtain library size (first value=library size)
--> Values save in GoogleDrive `009__*/samples_009.xlsx`. Histone-norm-library-size = library-size * SF. Using the non-reciprocal scaling factor, we increase the library-size; the more histone enriched, the more library size is increased, thus the more signal will decrease.

Now let's use these new histone-scaled library size and normalize with library-size,TMM or RLE. Let's use the **unique bam files** together with the **unique bam MACS2 raw files (xlsx, not the bed with pre-filtered qvalue)**

***Key points:***
- **Let's do 1 DiffBind per AB (H3K27me3, H3K4me3,...) and tissue (PSC, NPC); otherwise the TMM normalization may take all, unrelated, samples into account!** --> Files are `meta_sample_macs2raw_unique*.txt`
- **For the non-histone CutRun, I will use the library size non histone scaled in DiffBind to collect TMM normalized SF**; I tested with and without specifying library size; and it does not change a lot the SF... Let's better use the one RiP method w/o providing the library size! Should provide BETTER correction



```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# ONE PER ONE
## NPC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(12342151,11660741,9417384,21433114), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_H3K27me3.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_H3K27me3_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_H3K27me3_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()






# NPC_H3K4me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K4me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K4me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K4me3.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_H3K4me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_H3K4me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(7793480, 11433838, 9269965, 8910719), normalize = DBA_NORM_TMM) # 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_H3K4me3.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_H3K4me_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_H3K4me_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()



# ALL TOGETHER FOR PCA/HEATMAP PLOT
## NPC
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(12342151,11660741,9417384,21433114,7793480, 11433838, 9269965, 8910719), normalize = DBA_NORM_TMM) # 
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()


```

--> experimental batch effect.. Samples cluster more by experiment (005; 008) than per genotypes. However, applying Greylist/blacklist/LibSpikeIn scaled increase the genotype clustering.
----> GOOD


Let's try to **reduce batch effect using `dba.contrast block**`, discuss [here](https://support.bioconductor.org/p/96441/)




```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# ONE PER ONE
## NPC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
## apply blocking replicate factor 
sample_count_blackgreylist_block =  dba.contrast(sample_count_blackgreylist, block=DBA_REPLICATE, minMembers = 2)

sample_count_blackgreylist_block <- dba.contrast(sample_count_blackgreylist,
                           design="~Replicate + Treatment", minMembers =2, 
                           contrast=c("Treatment", "WT", "KO"))


### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist_block, library = c(12342151,11660741,9417384,21433114), normalize = DBA_NORM_TMM) 


sample_count_blackgreylist_LibHistoneScaled_TMM_analyze <- dba.analyze(sample_count_blackgreylist_LibHistoneScaled_TMM)

sample_count_blackgreylist_LibHistoneScaled_TMM_analyze.DB <- dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_analyze, contrast=1, method=DBA_DESEQ2_BLOCK, bCounts=TRUE)

dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM_analyze, bRetrieve=TRUE)






### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_H3K27me3_blackgreylist_block_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_H3K27me3_blackgreylist_block_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
#### Here is to retrieve the scaling factor value
#sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
#console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
#writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_H3K27me3.txt")



```

--> Playing with block factor do not change anything on the SF values, nor PCA/heatmap plots...


#### E coli bam spike in scaling method


Lets try to use the **E coli bam spike in scaling method** (add MG1655 bam files in the meta file):
--> For NPC samples only
--> Vignete 7.6 from [DiffBind doc](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) to check how apply spikein

**spikein_LIB normalization**

```bash
conda activate DiffBind
```

```R
library("DiffBind") 


# ONE PER ONE
## NPC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K27me3_MG1655bam.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3_MG1655bam.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3_MG1655bam.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### Spike in normalization 
sample_count_blackgreylist_LIB_spikein = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_LIB, spikein=TRUE) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LIB_spikein_SF = dba.normalize(sample_count_blackgreylist_LIB_spikein, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LIB_spikein_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LIB_spikein_SF_NPC_H3K27me3.txt")

### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_H3K27me3_blackgreylist_LIB_spikein_SF.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LIB_spikein)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_H3K27me3_blackgreylist_LIB_spikein_SF.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LIB_spikein,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
```


*NOTE: if apply the TMM, it seems to 'erase' the spike in bam effect; because in samples_005.xlsx, the DiffBind TMM (classic) SF = TMM_spikein_SF. So maybe TMM_spikein_SF is the way to go? NEED RNASEQ!!!*





# Generate Spike in scaled bigwig

--> Reciprocal from DiffBind_TMM is to be used when converting bam to bigwig!


### MG1655/E coli scaled bigwig


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM.sh # 14461837 ok
```


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
Check some known target regulated in 2months neurons:
--> NEUROG2 seems less in KO which is good.
--> EFNA5 tiny decrease in KO (only in normalized data!)
--> GRIK3 tiny increase in KO

--> Something is WEIRD... When samples have MORE spike in, their signal should be reduced, as they overall have more DNA; but if I used the reciprocal from DiffBind_TMM; this is not respected (ie. sample with more spike in, we increased their signal...!)... That is true for both histone/MG1655-spike in DiffBind TMM norm...
----> What should be the BEST to use, is then the NON-reciprocal_DiffBind_TMM !!!
------> Let's try and compare with gene expression...! Maybe it is still good as we take into account the library size with the DiffBind_TMM method??
--------> YESSS in the end we correct the library size with the SF!!!!! So we 're good!!! reciprocal DiffBind_TMM IS TO BE USED!!
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# THOR

Let's use THOR, notably to have IGG scaled bigwig...!

Comparison to do; NPC WT vs KO:
- H3K27me3
- H3K4me3


--> SF to use in THOR are the **reciprocal of MG1655_DiffBind_TMM**
--> Configs file created manually as `output/THOR/NPC_H3K27me3.config`


*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

*Follow these parameters: `WTvsHET_unique_Keepdup` (perform best in previous CutRun)*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

# classic DiffBind TMM method
sbatch scripts/THOR_NPC_H3K27me3.sh # 14462899 ok
sbatch scripts/THOR_NPC_H3K4me3.sh # 14462901 ok


# TMM from THOR soleley (THOR default)
sbatch scripts/THOR_NPC_H3K27me3_TMM.sh # 19142097 ok


# LIB_spikein_SF method (bam spike in)
sbatch scripts/THOR_NPC_H3K27me3_LIB_spikein.sh # 16224494 ok



```

Comparing **raw vs DiffBind_TMM vs THOR** on IGV
- H3K4me3: THOR is the cleaner, replicates looks more similar
- H3K27me3: THOR is the cleaner, replicates looks more similar 


Generate median tracks:
```bash
conda activate BedToBigwig
# classic DiffBind TMM method
sbatch scripts/bigwigmerge_THOR.sh # 15516016 ok
sbatch scripts/bigwigmerge_THOR_H3K4me3.sh # 15517815 ok

# LIB_spikein_SF method (bam spike in)
sbatch --dependency=afterany:16224494 scripts/bigwigmerge_THOR_H3K27me3_LIB_spikein.sh # 16224515 ok

```

--> TMM and LIB spike in method look very similar, at least for the few genes I checked 



## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R

# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# H3K27me3 TMM method (classic)
diffpeaks <- read_tsv("output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_WTvsKO_H3K27me3/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_NPC_WTvsKO_H3K27me3/log2FC_qval40.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 40) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO_qval40") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 10) %>%
  write_tsv("output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval10.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())




# H3K4me3
diffpeaks <- read_tsv("output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_WTvsKO_H3K4me3/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_NPC_WTvsKO_H3K4me3/log2FC_qval40.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 40) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO_qval40") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval50.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 20) %>%
  group_by(X6) %>%
  summarise(n = n())


# H3K27me3 LIB method
diffpeaks <- read_tsv("output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/NPCWTvsKOH3K27me3LIBspikein-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO_qval50") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 40) %>%
  write_tsv("output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval40.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())


```

- *NOTE: FC positive = less in KO; negative = more in KO*

**Optimal qvalue:**
--> *H3K27me3*; qval 30 for TMM and LIB (qval40 ok also,more stringeant)
--> *H3K4me3*; qval 20 (qval30 ok also,more stringeant)





# deepTools plot

On all genes, compare raw, DiffBind_TMM, THOR bigwigs


```bash
conda activate deeptools

# All genes
## H3K27me3
sbatch scripts/matrix_TSS_10kb_H3K27me3_raw_allGenes.sh # 14466946 fail; 14477113 ok
sbatch --dependency=afterany:14461837 scripts/matrix_TSS_10kb_H3K27me3_DiffBindTMM_allGenes.sh # 14470202 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_allGenes.sh # 14471783 fail; 14677072 ok

## H3K4me3
sbatch scripts/matrix_TSS_10kb_H3K4me3_raw_allGenes.sh # 14467901 fail; 14477116 fail (missabotated sample); 14516566 ok
sbatch --dependency=afterany:14461837 scripts/matrix_TSS_10kb_H3K4me3_DiffBindTMM_allGenes.sh # 14470728 ok
sbatch --dependency=afterany:14462899:14462901 scripts/matrix_TSS_10kb_H3K4me3_THOR_allGenes.sh # 14472822 ok


# median tracks
## all genes
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_allGenes.sh # 15516179 ok
sbatch scripts/matrix_TSS_5kb_H3K27me3_median_THOR_allGenes.sh # 15516394 ok
sbatch scripts/matrix_TSS_5kb_H3K4me3_median_THOR_allGenes.sh # 15516489 ok
sbatch scripts/matrix_TSS_2kb_H3K4me3_median_THOR_allGenes.sh # 15518489 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THORLIBspikein_allGenes.sh # 17135510 ok

## THOR diff. regions
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene.sh # 17490761 ok


## only genes with peak in WT and or KO qval 2.3
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks.sh # 15523775 disapear so relaunche as 15649755 ok
sbatch scripts/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks.sh # 15523776 ok
sbatch scripts/matrix_TSS_5kb_H3K4me3_median_THOR_genePeaks.sh # 15523800 ok
sbatch scripts/matrix_TSS_2kb_H3K4me3_median_THOR_genePeaks.sh # 15523832 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THORLIBspikein_genePeaks.sh # 17135552 ok



## only genes with peak in WT and or KO qval 3
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q3.sh # 15529967 ok
sbatch scripts/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks_macs2q3.sh # 15529966 ok
#XX sbatch scripts/matrix_TSS_5kb_H3K4me3_median_THOR_genePeaks_macs2q3.sh #  x
#XX sbatch scripts/matrix_TSS_2kb_H3K4me3_median_THOR_genePeaks_macs2q3.sh #  x

## only genes with peak in WT and or KO qval 4
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q4.sh # 15529977 ok
sbatch scripts/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks_macs2q4.sh # 15529988 ok
#XX sbatch scripts/matrix_TSS_5kb_H3K4me3_median_THOR_genePeaks_macs2q4.sh #  x
#XX sbatch scripts/matrix_TSS_2kb_H3K4me3_median_THOR_genePeaks_macs2q4.sh #  x

## only genes with peak in WT and or KO qval 5
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q5.sh # 15635588 ok

## only genes with peak in WT and or KO qval 1.3
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q1.30103.sh # 15641090 ok

## only genes with peak in WT and or KO qval 8
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q8.sh # 15641207/15641968 fail not sure why; 15642663 ok

## only genes with peak in WT and or KO qval 10
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q10.sh # 15641571 ok

# pearson corr plots
sbatch scripts/multiBigwigSummary_H3K27me3_raw.sh # 14480543 ok fail erase; 14554827 ok
sbatch scripts/multiBigwigSummary_H3K27me3_DiffBindTMM.sh # 14481446 ok fail erase; 14554967 ok
sbatch scripts/multiBigwigSummary_H3K27me3_THOR.sh # 14482864 ok fail erase; 14554971 ok
sbatch scripts/multiBigwigSummary_H3K27me3_THORLIBspikein.sh # 17135572 ok

sbatch scripts/multiBigwigSummary_H3K4me3_raw.sh # 14481079 fail (missabotated sample); 14516643 ok fail erase; 14554975 ok
sbatch scripts/multiBigwigSummary_H3K4me3_DiffBindTMM.sh # 14482286 ok fail erase; 14554976 ok
sbatch scripts/multiBigwigSummary_H3K4me3_THOR.sh # 14522864 ok fail erase; 14554977 ok


## H3K27me3 signal in EZH2/SUZ12 overlapping H3K27me3 for WT and KO (EZH2 peak from 001005)
### WT KO separately
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_EZH2Peaks001005OverlapH3K27me3_WT.sh # 19610146 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_EZH2Peaks001005OverlapH3K27me3_KO.sh # 19610179 ok

sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT.sh # 19610226 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_KO.sh # 19610248 ok

### WT KO together within WT peaks
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_EZH2Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks.sh # 19615119 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks.sh # 19615129 ok


### WT KO together within their respective genotype peaks FAIL
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_EZH2Peaks001005OverlapH3K27me3_WTKO.sh # 19611003 ok but not great for comparison


```


--> experimental batch effect.. Cluster more per experiment than per genotype for raw, DiffBindTMM and THOR...
----> Maybe batch effect because of overall aspecific signal like in intergenic region? Let's do pearson corr plot in gene body region only, and then in gene body and promoters (add 2kb upstream TSS)
-----> **Batch effet removed with THOR TMM, still present with LIBspikein**

--> H3K27me3 signal spreading around PRC2 sites (EZH2/SUZ12) stronger in WT as compared to KO


```bash
conda activate deeptools


# generate bed file from the gtf
cat /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$7}}' | tr -d '";' > /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_gene.bed

cat /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' > /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol.bed

# pearson corr plots - gene body
sbatch scripts/multiBigwigSummary_H3K27me3_raw_BEDgene.sh # 14551329 fail; 14555735 ok
sbatch scripts/multiBigwigSummary_H3K27me3_DiffBindTMM_BEDgene.sh # 14551977 fail; 14555790 ok
sbatch scripts/multiBigwigSummary_H3K27me3_THOR_BEDgene.sh # 14552069 fail; 14555796 ok

sbatch scripts/multiBigwigSummary_H3K4me3_raw_BEDgene.sh # 14551893 fail; 14555902 ok
sbatch scripts/multiBigwigSummary_H3K4me3_DiffBindTMM_BEDgene.sh # 14552266 fail; 14555905 ok
sbatch scripts/multiBigwigSummary_H3K4me3_THOR_BEDgene.sh # 14552400 fail; 14556059 ok
```


--> using gene only or whole genome do NOT change anything..






# deepTools plot on THOR diff peak genes


- isolate peak gain / lost --> deepTool plot
- assign diff. peak to genes (done in ChIPseeker)
- generate gtf from diff gene list --> deepTool plot



```bash
# PEAK
# THOR diff peaks
## isolate gain / lost peaks
awk -F'\t' '$16 > 1' output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval40.bed > output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval40_gain.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval40.bed > output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval40_lost.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval20.bed > output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval20_gain.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval20.bed > output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval20_lost.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval30.bed > output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval30_gain.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval30.bed > output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval30_lost.bed


## deeptools plot
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_q30_peak.sh # 14898621 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_q40_peak.sh # 14898622 ok

sbatch scripts/matrix_TSS_5kb_H3K4me3_THOR_q20_peak.sh # 14898678 ok
sbatch scripts/matrix_TSS_5kb_H3K4me3_THOR_q30_peak.sh # 14898706 ok

### with median not separating up and down
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_q30_peak.sh # 15530182 ok
sbatch scripts/matrix_TSS_5kb_H3K4me3_median_THOR_q20_peak.sh # 15530189 ok



# GENE
### create gtf from gene list
#### create gtf of diff bound genes 
#### Modify the .txt file that list all genes so that it match gtf structure

## Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt > output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836_as_gtf_geneSymbol.txt

## Filter the gtf
grep -Ff output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_H3K27me3_q30_pos_Promoter_5.gtf


## files
output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram194.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram355.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol_Venndiagram462.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol_Venndiagram636.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram220.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt

meta/ENCFF159KBI_H3K27me3_q30_pos_Promoter_5.gtf
meta/ENCFF159KBI_H3K27me3_q30_neg_Promoter_5.gtf
meta/ENCFF159KBI_H3K27me3_q40_pos_Promoter_5.gtf
meta/ENCFF159KBI_H3K27me3_q40_neg_Promoter_5.gtf
meta/ENCFF159KBI_H3K4me3_q20_pos_Promoter_5.gtf
meta/ENCFF159KBI_H3K4me3_q20_neg_Promoter_5.gtf
meta/ENCFF159KBI_H3K4me3_q30_pos_Promoter_5.gtf
meta/ENCFF159KBI_H3K4me3_q30_neg_Promoter_5.gtf
meta/ENCFF159KBI_H3K4me3_q30_pos_Promoter_5.gtf
meta/ENCFF159KBI_H3K4me3_q30_neg_Promoter_5.gtf

## deeptools plot
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_q30_gene.sh # 14900049 fail gtf; 15121974 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_q40_gene.sh # 14900082 fail gtf; 15121985 ok
sbatch scripts/matrix_TSS_5kb_H3K4me3_THOR_q20_gene.sh # 14900149 fail gtf; 15122006 ok
sbatch scripts/matrix_TSS_5kb_H3K4me3_THOR_q30_gene.sh # 14900199 fail gtf; 15122007 ok

### with median not separating up and down
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_q30_gene.sh # 15647072 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_q40_gene.sh # 15647203 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THOR_q50_gene.sh # 15647306 ok

```

--> Mybe qvalue treshold could be increase to make gian lost more clear



# deepTools plot to compare bigwig Ferguson/Local maxima vs THOR

Lets's compare THOR and Ferguson method (ie. the identified peaks, and how the bigwig behave). We notably found with THOR more region that gain than regions that lost H3K27me3, but the opposite is found with Ferguson DESEQ2/EDGER method...

**Isolate gain lost peaks**
```bash
# PEAK
# EDGER diff peaks (GlmfitLikelihoodRatioRawCounts) 121 / 40 = 
output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.txt
## Extract all diff bound regions, in bed format
awk 'NR>1 {print $7 "\t" $8 "\t" $9 "\t" $6 "\t" $2}' output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.txt > output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.bed

awk '$4 < 0.05' output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.bed > output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.FDR05.bed
awk '$5 > 0' output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.FDR05.bed > output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.FDR05FCpos.bed
awk '$5 < 0' output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.FDR05.bed > output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.FDR05FCneg.bed


# DESEQ2 diff peaks (lfcShrinkNORMAL) 132 / 45
output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.txt
## Extract all diff bound regions, in bed format
awk 'NR>1 {print $8 "\t" $9 "\t" $10 "\t" $7 "\t" $3}' output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.txt > output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.bed

awk '$4 < 0.05' output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.bed > output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.FDR05.bed
awk '$5 > 0' output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.FDR05.bed > output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.FDR05FCpos.bed
awk '$5 < 0' output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.FDR05.bed > output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.FDR05FCneg.bed


#--> All files of interest below:
# Regions
## THOR gain lost
output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30_gain.bed
output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30_lost.bed
## EDGER gain lost - GlmfitLikelihoodRatioRawCounts
output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.FDR05FCpos.bed
output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.FDR05FCneg.bed
## DESEQ2 gain lost - lfcShrinkNORMAL
output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.FDR05FCpos.bed
output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.FDR05FCneg.bed

#Bigwigs
## THOR
output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw # WT
output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw # KO
## Ferguson/Local Maxima
output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median.bw # WT
output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median.bw # KO

```


Lets do **deepTools plot**:
- gain lost (THOR and DESEQ2)
- all genes
- consensus H3K27me3 peaks (q2.3)


```bash
conda activate deeptools
# deeptools plot

## gain lost (THOR and DESEQ2)
sbatch scripts/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-THORq30PosNegPeaks.sh # 38161517 ok
sbatch scripts/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-EDGERGlmfitLikelihoodRatioRawCountsFDR05FCposNegPeaks.sh # 38161537 ok
sbatch scripts/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-DESEQ2lfcShrinkNORMALFDR05FCposNegPeaks.sh # 38161566 ok

## all genes
sbatch scripts/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-allGenes.sh # 38615466 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3-LocalMaxima_THOR-allGenes.sh # 38616206 ok

## consensus H3K27me3 peaks
sbatch scripts/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-consensusPeaks.sh # 38616352 ok 
sbatch scripts/matrix_TSS_10kb_H3K27me3-LocalMaxima_THOR-consensusPeaks.sh # 38616376 ok
```


--> THOR diff regions also look different in Ferguson/LocalMaxima bigwigs! Good news, we see many regions gaining H3K27me3 as before!!

--> DESEQ2/EDGER diff regions are not detected looking at THOR bigwig!
  --> THOR may miss these ones, they are usually smaller in size! But real!

--> So sliding window method look cool for bigger regions, but DESEQ2 cool too, for smaller regions... Let's try another method (csaw?) to have the best of both worlds

--> For all genes and consensus peaks, WT and KO present same H3K27me3 levels, for both LocalMaxima and  THOR bigwigs (so it is only using the THOR diff sites that I see an increase of H3K27me3 level...)




# ChIPseeker peak gene assignment

## From optimal qval bed files peaks
Let's assign **peak to genes from MACS2 peak**:

**Optimal qvalue** according to IGV:
- NPC_H3K27me3 WT and KO; 2.30103
- NPC_H3K4me3 WT and KO; 2.30103


--> Assign peak to genes for NPC:

```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")


# Import macs2 peaks
## H3K27me3
H3K27me3_WT_005 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_005_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3_WT_008 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_008_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3_KO_005 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_005_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3_KO_008 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_008_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

H3K27me3_WT_pool = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval10/NPC_WT_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3_KO_pool = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval10/NPC_KO_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

## H3K4me3
H3K4me3_WT_005 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K4me3_005_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K4me3_WT_008 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K4me3_008_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K4me3_KO_005 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K4me3_005_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K4me3_KO_008 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K4me3_008_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

H3K4me3_WT_pool = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval4/NPC_WT_H3K4me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K4me3_KO_pool = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval4/NPC_KO_H3K4me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

# Tidy peaks 
## H3K27me3
H3K27me3_WT_005_gr = makeGRangesFromDataFrame(H3K27me3_WT_005,keep.extra.columns=TRUE)
H3K27me3_WT_008_gr = makeGRangesFromDataFrame(H3K27me3_WT_008,keep.extra.columns=TRUE)
H3K27me3_KO_005_gr = makeGRangesFromDataFrame(H3K27me3_KO_005,keep.extra.columns=TRUE)
H3K27me3_KO_008_gr = makeGRangesFromDataFrame(H3K27me3_KO_008,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3_WT_005=H3K27me3_WT_005_gr, H3K27me3_WT_008=H3K27me3_WT_008_gr,  H3K27me3_KO_005=H3K27me3_KO_005_gr, H3K27me3_KO_008=H3K27me3_KO_008_gr)

H3K27me3_WT_pool_gr = makeGRangesFromDataFrame(H3K27me3_WT_pool,keep.extra.columns=TRUE)
H3K27me3_KO_pool_gr = makeGRangesFromDataFrame(H3K27me3_KO_pool,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3_WT_pool=H3K27me3_WT_pool_gr, H3K27me3_KO_pool=H3K27me3_KO_pool_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_H3K27me3_pool.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_H3K27me3_pool.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
H3K27me3_WT_005_annot <- as.data.frame(peakAnnoList[["H3K27me3_WT_005"]]@anno)
H3K27me3_WT_008_annot <- as.data.frame(peakAnnoList[["H3K27me3_WT_008"]]@anno)
H3K27me3_KO_005_annot <- as.data.frame(peakAnnoList[["H3K27me3_KO_005"]]@anno)
H3K27me3_KO_008_annot <- as.data.frame(peakAnnoList[["H3K27me3_KO_008"]]@anno)

H3K27me3_WT_pool_annot <- as.data.frame(peakAnnoList[["H3K27me3_WT_pool"]]@anno)
H3K27me3_KO_pool_annot <- as.data.frame(peakAnnoList[["H3K27me3_KO_pool"]]@anno)

## Convert entrez gene IDs to gene symbols
H3K27me3_WT_005_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_WT_005_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_WT_005_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_WT_005_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_WT_008_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_WT_008_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_WT_008_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_WT_008_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K27me3_WT_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_WT_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_WT_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_WT_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


H3K27me3_KO_005_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_KO_005_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_KO_005_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_KO_005_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_KO_008_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_KO_008_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_KO_008_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_KO_008_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K27me3_KO_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_KO_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_KO_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_KO_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(H3K27me3_WT_005_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_WT_005_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_WT_008_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_WT_008_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_KO_005_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_KO_005_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_KO_008_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_KO_008_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_WT_pool_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_WT_pool_qval10.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_KO_pool_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_KO_pool_qval10.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K27me3_WT_005_annot_promoterAnd5 = tibble(H3K27me3_WT_005_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_WT_008_annot_promoterAnd5 = tibble(H3K27me3_WT_008_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_KO_005_annot_promoterAnd5 = tibble(H3K27me3_KO_005_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_KO_008_annot_promoterAnd5 = tibble(H3K27me3_KO_008_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

H3K27me3_WT_pool_annot_promoterAnd5 = tibble(H3K27me3_WT_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_KO_pool_annot_promoterAnd5 = tibble(H3K27me3_KO_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
H3K27me3_WT_005_annot_promoterAnd5_geneSymbol = H3K27me3_WT_005_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_WT_008_annot_promoterAnd5_geneSymbol = H3K27me3_WT_008_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_KO_005_annot_promoterAnd5_geneSymbol = H3K27me3_KO_005_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_KO_008_annot_promoterAnd5_geneSymbol = H3K27me3_KO_008_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

H3K27me3_WT_pool_annot_promoterAnd5_geneSymbol = H3K27me3_WT_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_KO_pool_annot_promoterAnd5_geneSymbol = H3K27me3_KO_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(H3K27me3_WT_005_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_WT_005_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_WT_008_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_WT_008_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_KO_005_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_KO_005_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_KO_008_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_KO_008_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

write.table(H3K27me3_WT_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_WT_pool_qval10_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_KO_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_KO_pool_qval10_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)




## H3K4me3
H3K4me3_WT_005_gr = makeGRangesFromDataFrame(H3K4me3_WT_005,keep.extra.columns=TRUE)
H3K4me3_WT_008_gr = makeGRangesFromDataFrame(H3K4me3_WT_008,keep.extra.columns=TRUE)
H3K4me3_KO_005_gr = makeGRangesFromDataFrame(H3K4me3_KO_005,keep.extra.columns=TRUE)
H3K4me3_KO_008_gr = makeGRangesFromDataFrame(H3K4me3_KO_008,keep.extra.columns=TRUE)
gr_list <- list(H3K4me3_WT_005=H3K4me3_WT_005_gr, H3K4me3_WT_008=H3K4me3_WT_008_gr,  H3K4me3_KO_005=H3K4me3_KO_005_gr,  H3K4me3_KO_008=H3K4me3_KO_008_gr)

H3K4me3_WT_pool_gr = makeGRangesFromDataFrame(H3K4me3_WT_pool,keep.extra.columns=TRUE)
H3K4me3_KO_pool_gr = makeGRangesFromDataFrame(H3K4me3_KO_pool,keep.extra.columns=TRUE)
gr_list <- list(H3K4me3_WT_pool=H3K4me3_WT_pool_gr, H3K4me3_KO_pool=H3K4me3_KO_pool_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_H3K4me3_pool.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_H3K4me3_pool.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
H3K4me3_WT_005_annot <- as.data.frame(peakAnnoList[["H3K4me3_WT_005"]]@anno)
H3K4me3_WT_008_annot <- as.data.frame(peakAnnoList[["H3K4me3_WT_008"]]@anno)
H3K4me3_KO_005_annot <- as.data.frame(peakAnnoList[["H3K4me3_KO_005"]]@anno)
H3K4me3_KO_008_annot <- as.data.frame(peakAnnoList[["H3K4me3_KO_008"]]@anno)

H3K4me3_WT_pool_annot <- as.data.frame(peakAnnoList[["H3K4me3_WT_pool"]]@anno)
H3K4me3_KO_pool_annot <- as.data.frame(peakAnnoList[["H3K4me3_KO_pool"]]@anno)
## Convert entrez gene IDs to gene symbols
H3K4me3_WT_005_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_WT_005_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_WT_005_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_WT_005_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_WT_008_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_WT_008_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_WT_008_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_WT_008_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K4me3_WT_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_WT_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_WT_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_WT_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



H3K4me3_KO_005_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_KO_005_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_KO_005_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_KO_005_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_KO_008_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_KO_008_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_KO_008_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_KO_008_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K4me3_KO_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_KO_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_KO_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_KO_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(H3K4me3_WT_005_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_WT_005_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_WT_008_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_WT_008_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_KO_005_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_KO_005_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_KO_008_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_KO_008_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_WT_pool_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_WT_pool_qval4.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_KO_pool_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_KO_pool_qval4.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K4me3_WT_005_annot_promoterAnd5 = tibble(H3K4me3_WT_005_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_WT_008_annot_promoterAnd5 = tibble(H3K4me3_WT_008_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_KO_005_annot_promoterAnd5 = tibble(H3K4me3_KO_005_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_KO_008_annot_promoterAnd5 = tibble(H3K4me3_KO_008_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

H3K4me3_WT_pool_annot_promoterAnd5 = tibble(H3K4me3_WT_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_KO_pool_annot_promoterAnd5 = tibble(H3K4me3_KO_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
H3K4me3_WT_005_annot_promoterAnd5_geneSymbol = H3K4me3_WT_005_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_WT_008_annot_promoterAnd5_geneSymbol = H3K4me3_WT_008_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_KO_005_annot_promoterAnd5_geneSymbol = H3K4me3_KO_005_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_KO_008_annot_promoterAnd5_geneSymbol = H3K4me3_KO_008_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

H3K4me3_WT_pool_annot_promoterAnd5_geneSymbol = H3K4me3_WT_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_KO_pool_annot_promoterAnd5_geneSymbol = H3K4me3_KO_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(H3K4me3_WT_005_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_WT_005_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_WT_008_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_WT_008_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_KO_005_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_KO_005_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_KO_008_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_KO_008_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

write.table(H3K4me3_WT_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_WT_pool_qval4_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_KO_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_KO_pool_qval4_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

```


**NPC**:

--> Export gene list (`output/ChIPseeker/annotation_macs2_H3K*me3_*_00*_qval2.30103_promoterAnd5_geneSymbol.txt`) to Online Venn diagram to check replicate homonegeity between `CutRun__005` and `CutRun__008`




## From consensus peak

### Consensus peak H3K27me3 no extension/extension, qvalue 2.3
- consensus peak H3K27me3, WT KO, no extension: `output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge.bed`
- consensus peak H3K27me3, WT KO, 100bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge100bp.bed`
- consensus peak H3K27me3, WT KO, 500bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge500bp.bed`



```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")


# Import consensus peak
NPC_WTKO_H3K27me3_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
NPC_WTKO_H3K27me3_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
NPC_WTKO_H3K27me3_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K27me3
NPC_WTKO_H3K27me3_pool_peaks_merge_gr = makeGRangesFromDataFrame(NPC_WTKO_H3K27me3_pool_peaks_merge,keep.extra.columns=TRUE)
NPC_WTKO_H3K27me3_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(NPC_WTKO_H3K27me3_pool_peaks_merge100bp,keep.extra.columns=TRUE)
NPC_WTKO_H3K27me3_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(NPC_WTKO_H3K27me3_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(NPC_WTKO_H3K27me3_pool_peaks_merge=NPC_WTKO_H3K27me3_pool_peaks_merge_gr, NPC_WTKO_H3K27me3_pool_peaks_merge100bp=NPC_WTKO_H3K27me3_pool_peaks_merge100bp_gr, NPC_WTKO_H3K27me3_pool_peaks_merge500bp=NPC_WTKO_H3K27me3_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_NPC_WTKO_H3K27me3_pool_peaks_merge.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_NPC_WTKO_H3K27me3_pool_peaks_merge.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
NPC_WTKO_H3K27me3_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["NPC_WTKO_H3K27me3_pool_peaks_merge"]]@anno)
NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["NPC_WTKO_H3K27me3_pool_peaks_merge100bp"]]@anno)
NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["NPC_WTKO_H3K27me3_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
NPC_WTKO_H3K27me3_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K27me3_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
NPC_WTKO_H3K27me3_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K27me3_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(NPC_WTKO_H3K27me3_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
NPC_WTKO_H3K27me3_pool_peaks_merge_annot_promoterAnd5 = tibble(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
NPC_WTKO_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol = NPC_WTKO_H3K27me3_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(NPC_WTKO_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```







### Consensus peak H3K4me3 no extension/extension, qvalue 2.3
- consensus peak H3K4me3, WT KO KOEF1aEZH1, no extension: `output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge.bed`
- consensus peak H3K4me3, WT KO KOEF1aEZH1, 100bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge100bp.bed`
- consensus peak H3K4me3, WT KO KOEF1aEZH1, 500bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge500bp.bed`



```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")


# Import consensus peak
NPC_WTKO_H3K4me3_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
NPC_WTKO_H3K4me3_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
NPC_WTKO_H3K4me3_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K4me3
NPC_WTKO_H3K4me3_pool_peaks_merge_gr = makeGRangesFromDataFrame(NPC_WTKO_H3K4me3_pool_peaks_merge,keep.extra.columns=TRUE)
NPC_WTKO_H3K4me3_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(NPC_WTKO_H3K4me3_pool_peaks_merge100bp,keep.extra.columns=TRUE)
NPC_WTKO_H3K4me3_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(NPC_WTKO_H3K4me3_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(NPC_WTKO_H3K4me3_pool_peaks_merge=NPC_WTKO_H3K4me3_pool_peaks_merge_gr, NPC_WTKO_H3K4me3_pool_peaks_merge100bp=NPC_WTKO_H3K4me3_pool_peaks_merge100bp_gr, NPC_WTKO_H3K4me3_pool_peaks_merge500bp=NPC_WTKO_H3K4me3_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_NPC_WTKO_H3K4me3_pool_peaks_merge.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_NPC_WTKO_H3K4me3_pool_peaks_merge.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
NPC_WTKO_H3K4me3_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["NPC_WTKO_H3K4me3_pool_peaks_merge"]]@anno)
NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["NPC_WTKO_H3K4me3_pool_peaks_merge100bp"]]@anno)
NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["NPC_WTKO_H3K4me3_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
NPC_WTKO_H3K4me3_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K4me3_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
NPC_WTKO_H3K4me3_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K4me3_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(NPC_WTKO_H3K4me3_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_NPC_WTKO_H3K4me3_pool_peaks_merge_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
NPC_WTKO_H3K4me3_pool_peaks_merge_annot_promoterAnd5 = tibble(NPC_WTKO_H3K4me3_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
NPC_WTKO_H3K4me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol = NPC_WTKO_H3K4me3_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(NPC_WTKO_H3K4me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_NPC_WTKO_H3K4me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_NPC_WTKO_H3K4me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_NPC_WTKO_H3K4me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```




# create gtf of gene with peak in WT and/or KO

```bash
# isolate all the genes bound with H3K27me3/H3K4me3 in WT and or KO
cat output/ChIPseeker/annotation_macs2_H3K4me3_WT_pool_qval4_promoterAnd5_geneSymbol.txt output/ChIPseeker/annotation_macs2_H3K4me3_KO_pool_qval4_promoterAnd5_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_macs2_H3K4me3_WTKO_pool_qval4_promoterAnd5_geneSymbol.txt

cat output/ChIPseeker/annotation_macs2_H3K27me3_WT_pool_qval8_promoterAnd5_geneSymbol.txt output/ChIPseeker/annotation_macs2_H3K27me3_KO_pool_qval8_promoterAnd5_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_macs2_H3K27me3_WTKO_pool_qval8_promoterAnd5_geneSymbol.txt


### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure

## Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_macs2_H3K4me3_WTKO_pool_qval4_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_macs2_H3K4me3_WTKO_pool_qval4_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_macs2_H3K27me3_WTKO_pool_qval8_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_macs2_H3K27me3_WTKO_pool_qval8_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_THOR_H3K27me3_q30_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_THOR_H3K27me3_q30_promoterAnd5_as_gtf_geneSymbol.txt

## Filter the gtf
grep -Ff output/ChIPseeker/annotation_macs2_H3K4me3_WTKO_pool_qval4_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_macs2_H3K4me3_WTKO_pool_qval4_Promoter_5.gtf

grep -Ff output/ChIPseeker/annotation_macs2_H3K27me3_WTKO_pool_qval8_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_macs2_H3K27me3_WTKO_pool_qval8_Promoter_5.gtf

grep -Ff output/ChIPseeker/annotation_THOR_H3K27me3_q30_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_THOR_H3K27me3_q30_Promoter_5.gtf


```





## From THOR diff bound peaks
Let's assign **peak to genes from THORs peak**:

**Optimal qvalue** according to IGV:
- NPC_H3K27me3 WT vs KO; qval 30
- NPC_H3K4me3 WT vs KO; qval 20


--> Assign peak to genes for NPC:

```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")


# Import THOR peaks
# TMM
## H3K27me3 _ q30
H3K27me3_q30_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K27me3_q30_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 
H3K27me3_q30 = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16)

## H3K27me3 _ q40
H3K27me3_q40_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval40.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K27me3_q40_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval40.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 
H3K27me3_q40 = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval40.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16)

## H3K27me3 _ q50
H3K27me3_q50_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval50.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K27me3_q50_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval50.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 
H3K27me3_q50 = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval50.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16)


# LIBspikein
## H3K27me3 _ q30
H3K27me3_q30_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval30.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K27me3_q30_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval30.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 

## H3K27me3 _ q40
H3K27me3_q40_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval40.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K27me3_q40_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval40.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 

## H3K27me3 _ q50
H3K27me3_q50_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval50.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K27me3_q50_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein/THOR_qval50.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 



# Tidy peaks 
## H3K27me3
H3K27me3_q30_pos_gr = makeGRangesFromDataFrame(H3K27me3_q30_pos,keep.extra.columns=TRUE)
H3K27me3_q30_neg_gr = makeGRangesFromDataFrame(H3K27me3_q30_neg,keep.extra.columns=TRUE)
H3K27me3_q40_pos_gr = makeGRangesFromDataFrame(H3K27me3_q40_pos,keep.extra.columns=TRUE)
H3K27me3_q40_neg_gr = makeGRangesFromDataFrame(H3K27me3_q40_neg,keep.extra.columns=TRUE)
H3K27me3_q50_pos_gr = makeGRangesFromDataFrame(H3K27me3_q50_pos,keep.extra.columns=TRUE)
H3K27me3_q50_neg_gr = makeGRangesFromDataFrame(H3K27me3_q50_neg,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3_q30_pos=H3K27me3_q30_pos_gr, H3K27me3_q30_neg=H3K27me3_q30_neg_gr,  H3K27me3_q40_pos=H3K27me3_q40_pos_gr, H3K27me3_q40_neg=H3K27me3_q40_neg_gr,
H3K27me3_q50_pos=H3K27me3_q50_pos_gr, H3K27me3_q50_neg=H3K27me3_q50_neg_gr)

H3K27me3_q30_gr = makeGRangesFromDataFrame(H3K27me3_q30,keep.extra.columns=TRUE)
H3K27me3_q40_gr = makeGRangesFromDataFrame(H3K27me3_q40,keep.extra.columns=TRUE)
H3K27me3_q50_gr = makeGRangesFromDataFrame(H3K27me3_q50,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3_q30=H3K27me3_q30_gr, H3K27me3_q40=H3K27me3_q40_gr,  H3K27me3_q50=H3K27me3_q50_gr)


# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_THOR_H3K27me3.pdf", width = 8, height = 3)
pdf("output/ChIPseeker/plotAnnoBar_THORLIBspikein_H3K27me3.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_THOR_H3K27me3.pdf", width = 8, height = 3)
pdf("output/ChIPseeker/plotDistToTSS_THORLIBspikein_H3K27me3.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
H3K27me3_q30_pos_annot <- as.data.frame(peakAnnoList[["H3K27me3_q30_pos"]]@anno)
H3K27me3_q30_neg_annot <- as.data.frame(peakAnnoList[["H3K27me3_q30_neg"]]@anno)
H3K27me3_q40_pos_annot <- as.data.frame(peakAnnoList[["H3K27me3_q40_pos"]]@anno)
H3K27me3_q40_neg_annot <- as.data.frame(peakAnnoList[["H3K27me3_q40_neg"]]@anno)
H3K27me3_q50_pos_annot <- as.data.frame(peakAnnoList[["H3K27me3_q50_pos"]]@anno)
H3K27me3_q50_neg_annot <- as.data.frame(peakAnnoList[["H3K27me3_q50_neg"]]@anno)

H3K27me3_q30_annot <- as.data.frame(peakAnnoList[["H3K27me3_q30"]]@anno)
H3K27me3_q40_annot <- as.data.frame(peakAnnoList[["H3K27me3_q40"]]@anno)
H3K27me3_q50_annot <- as.data.frame(peakAnnoList[["H3K27me3_q50"]]@anno)


## Convert entrez gene IDs to gene symbols
H3K27me3_q30_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q30_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q30_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q30_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_q30_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q30_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q30_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q30_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K27me3_q40_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q40_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q40_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q40_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_q40_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q40_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q40_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q40_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K27me3_q50_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q50_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_q50_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q50_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K27me3_q30_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q30_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q30_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q30_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_q40_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q40_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q40_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q40_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_q50_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q50_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(H3K27me3_q30_pos_annot, file="output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q30_pos.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q30_neg_annot, file="output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q30_neg.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q40_pos_annot, file="output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q40_pos.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q40_neg_annot, file="output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q40_neg.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q50_pos_annot, file="output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q50_pos.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q50_neg_annot, file="output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q50_neg.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

write.table(H3K27me3_q30_annot, file="output/ChIPseeker/annotation_THOR_H3K27me3_q30.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q40_annot, file="output/ChIPseeker/annotation_THOR_H3K27me3_q40.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q50_annot, file="output/ChIPseeker/annotation_THOR_H3K27me3_q50.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K27me3_q30_pos_annot_promoterAnd5 = tibble(H3K27me3_q30_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q30_neg_annot_promoterAnd5 = tibble(H3K27me3_q30_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q40_pos_annot_promoterAnd5 = tibble(H3K27me3_q40_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q40_neg_annot_promoterAnd5 = tibble(H3K27me3_q40_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q50_pos_annot_promoterAnd5 = tibble(H3K27me3_q50_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q50_neg_annot_promoterAnd5 = tibble(H3K27me3_q50_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

H3K27me3_q30_annot_promoterAnd5 = tibble(H3K27me3_q30_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q40_annot_promoterAnd5 = tibble(H3K27me3_q40_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q50_annot_promoterAnd5 = tibble(H3K27me3_q50_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))



### Save output gene lists
H3K27me3_q30_pos_annot_promoterAnd5_geneSymbol = H3K27me3_q30_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q30_neg_annot_promoterAnd5_geneSymbol = H3K27me3_q30_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q40_pos_annot_promoterAnd5_geneSymbol = H3K27me3_q40_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q40_neg_annot_promoterAnd5_geneSymbol = H3K27me3_q40_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q50_pos_annot_promoterAnd5_geneSymbol = H3K27me3_q50_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q50_neg_annot_promoterAnd5_geneSymbol = H3K27me3_q50_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

H3K27me3_q30_annot_promoterAnd5_geneSymbol = H3K27me3_q30_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q40_annot_promoterAnd5_geneSymbol = H3K27me3_q40_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q50_annot_promoterAnd5_geneSymbol = H3K27me3_q50_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(H3K27me3_q30_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q30_pos_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q30_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q30_neg_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q40_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q40_pos_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q40_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q40_neg_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q50_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q50_pos_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q50_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORLIBspikein_H3K27me3_q50_neg_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

write.table(H3K27me3_q30_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K27me3_q30_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q40_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K27me3_q40_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q50_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K27me3_q50_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



## H3K4me3
## H3K4me3 _ q20
H3K4me3_q20_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval20.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K4me3_q20_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval20.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 
## H3K4me3 _ q30
H3K4me3_q30_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval30.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K4me3_q30_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval30.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 
## H3K4me3 _ q40
H3K4me3_q40_pos = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval40.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K4me3_q40_neg = as_tibble(read.table('output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval40.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 



# Tidy peaks 
## H3K4me3
H3K4me3_q20_pos_gr = makeGRangesFromDataFrame(H3K4me3_q20_pos,keep.extra.columns=TRUE)
H3K4me3_q20_neg_gr = makeGRangesFromDataFrame(H3K4me3_q20_neg,keep.extra.columns=TRUE)
H3K4me3_q30_pos_gr = makeGRangesFromDataFrame(H3K4me3_q30_pos,keep.extra.columns=TRUE)
H3K4me3_q30_neg_gr = makeGRangesFromDataFrame(H3K4me3_q30_neg,keep.extra.columns=TRUE)
H3K4me3_q40_pos_gr = makeGRangesFromDataFrame(H3K4me3_q40_pos,keep.extra.columns=TRUE)
H3K4me3_q40_neg_gr = makeGRangesFromDataFrame(H3K4me3_q40_neg,keep.extra.columns=TRUE)

gr_list <- list(H3K4me3_q20_pos=H3K4me3_q20_pos_gr, H3K4me3_q20_neg=H3K4me3_q20_neg_gr, H3K4me3_q30_pos=H3K4me3_q30_pos_gr, H3K4me3_q30_neg=H3K4me3_q30_neg_gr,  H3K4me3_q40_pos=H3K4me3_q40_pos_gr, H3K4me3_q40_neg=H3K4me3_q40_neg_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_THOR_H3K4me3.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_THOR_H3K4me3.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
H3K4me3_q20_pos_annot <- as.data.frame(peakAnnoList[["H3K4me3_q20_pos"]]@anno)
H3K4me3_q20_neg_annot <- as.data.frame(peakAnnoList[["H3K4me3_q20_neg"]]@anno)
H3K4me3_q30_pos_annot <- as.data.frame(peakAnnoList[["H3K4me3_q30_pos"]]@anno)
H3K4me3_q30_neg_annot <- as.data.frame(peakAnnoList[["H3K4me3_q30_neg"]]@anno)
H3K4me3_q40_pos_annot <- as.data.frame(peakAnnoList[["H3K4me3_q40_pos"]]@anno)
H3K4me3_q40_neg_annot <- as.data.frame(peakAnnoList[["H3K4me3_q40_neg"]]@anno)


## Convert entrez gene IDs to gene symbols
H3K4me3_q20_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_q20_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_q20_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_q20_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_q20_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_q20_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_q20_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_q20_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K4me3_q30_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_q30_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_q30_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_q30_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_q30_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_q30_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_q30_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_q30_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

H3K4me3_q40_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_q40_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_q40_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_q40_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_q40_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_q40_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_q40_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_q40_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(H3K4me3_q20_pos_annot, file="output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_q20_neg_annot, file="output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_q30_pos_annot, file="output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_q30_neg_annot, file="output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_q40_pos_annot, file="output/ChIPseeker/annotation_THOR_H3K4me3_q40_pos.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_q40_neg_annot, file="output/ChIPseeker/annotation_THOR_H3K4me3_q40_neg.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K4me3_q20_pos_annot_promoterAnd5 = tibble(H3K4me3_q20_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_q20_neg_annot_promoterAnd5 = tibble(H3K4me3_q20_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_q30_pos_annot_promoterAnd5 = tibble(H3K4me3_q30_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_q30_neg_annot_promoterAnd5 = tibble(H3K4me3_q30_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_q40_pos_annot_promoterAnd5 = tibble(H3K4me3_q40_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_q40_neg_annot_promoterAnd5 = tibble(H3K4me3_q40_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))




### Save output gene lists
H3K4me3_q20_pos_annot_promoterAnd5_geneSymbol = H3K4me3_q20_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_q20_neg_annot_promoterAnd5_geneSymbol = H3K4me3_q20_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_q30_pos_annot_promoterAnd5_geneSymbol = H3K4me3_q30_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_q30_neg_annot_promoterAnd5_geneSymbol = H3K4me3_q30_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_q40_pos_annot_promoterAnd5_geneSymbol = H3K4me3_q40_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_q40_neg_annot_promoterAnd5_geneSymbol = H3K4me3_q40_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()



write.table(H3K4me3_q20_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_q20_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_q30_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_q30_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_q40_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K4me3_q40_pos_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_q40_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_H3K4me3_q40_neg_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


```


# Functional analysis with enrichR

Functional analysis **enrichR with THOR diff bound genes**; using webtool Venn diagram I isolated the specific genes that gain / lost H3K*me3 `output/ChIPseeker/annotation_THOR_H3K*me3_q*_pos_promoterAnd5_geneSymbol_Venndiagram*.txt`


**IMPOPRTANT NOTE: Run the reading and processing ONE BY ONE !!! Otherwise, lead to bug!!!!**

```R
# library
library("tidyverse")
library("enrichR")
library("ggrepel")

# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 

### GeneSymbol list of signif gain/lost H3K27me3 in each genotypes
output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram194.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram355.txt

output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol_Venndiagram462.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol_Venndiagram636.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram220.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt


### GeneSymbol list of signif gain/lost H3K27me3 in each genotypes and DEGs (Carolina RNAseq)
output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt

output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"
# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string
# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K27me3_q30.pdf", width=8, height=5)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K4me3_q20.pdf", width=8, height=6)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K4me3_q30.pdf", width=8, height=6)

pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K27me3_q40.pdf", width=8, height=6)

pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K27me3_q30_Lost_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=6)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=6)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="DEGs and H3K27me3",   # H3K27me3  H3K4me3
                    labels = c("down-reg and Gain", "up-reg and Lost"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "GO_Biological_Process_2023") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Molecular_Function_2023") # 
### GeneSymbol list of signif up/down genes in each genotypes
output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram194.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram355.txt

output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol_Venndiagram462.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol_Venndiagram636.txt


output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram220.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt


### GeneSymbol list of signif gain/lost H3K27me3 in each genotypes and DEGs (Carolina RNAseq)
output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt

output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)
# Extracting KEGG data and assigning types
up <- eup$GO_Molecular_Function_2023
down <- edown$GO_Molecular_Function_2023
up$type <- "up"
down$type <- "down"
# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)
# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string
# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K4me3_q30.pdf", width=8, height=3)
pdf("output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K4me3_q20.pdf", width=8, height=5)
pdf("output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K27me3_q30.pdf", width=8, height=6)

pdf("output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K27me3_q40.pdf", width=8, height=5)
pdf("output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K27me3_q30_Lost_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)
pdf("output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K4me3",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "GO_Molecular_Function_2023") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()

## save output
write.table(gos, "output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Cellular_Component_2023") # 

### GeneSymbol list of signif up/down genes in each genotypes
output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram194.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram355.txt

output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol_Venndiagram462.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol_Venndiagram636.txt


output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram220.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt

### GeneSymbol list of signif gain/lost H3K27me3 in each genotypes and DEGs (Carolina RNAseq)
output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt

output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Cellular_Component_2023
down <- edown$GO_Cellular_Component_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)



# Plotting with enhanced aesthetics


pdf("output/GO/enrichR_GO_Cellular_Component_2023_THOR_H3K27me3_q30.pdf", width=8, height=4)
pdf("output/GO/enrichR_GO_Cellular_Component_2023_THOR_H3K4me3_q20.pdf", width=8, height=5)
pdf("output/GO/enrichR_GO_Cellular_Component_2023_THOR_H3K4me3_q30.pdf", width=8, height=3)

pdf("output/GO/enrichR_GO_Cellular_Component_2023_THOR_H3K27me3_q40.pdf", width=8, height=4)

pdf("output/GO/enrichR_GO_Cellular_Component_2023_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=4)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K4me3",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "GO_Cellular_Component_2023") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Cellular_Component_2023_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("KEGG_2021_Human") # 

### GeneSymbol list of signif up/down genes in each genotypes
output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram194.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram355.txt

output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol_Venndiagram462.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol_Venndiagram636.txt


output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram220.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt

### GeneSymbol list of signif gain/lost H3K27me3 in each genotypes and DEGs (Carolina RNAseq)
output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt

output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
down <- edown$KEGG_2021_Human
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)
# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string
# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics

pdf("output/GO/enrichR_KEGG_2021_Human_THOR_H3K4me3_q30.pdf", width=8, height=4)
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_H3K4me3_q20.pdf", width=8, height=5)
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_H3K27me3_q30.pdf", width=8, height=3)

pdf("output/GO/enrichR_KEGG_2021_Human_THOR_H3K27me3_q40.pdf", width=8, height=2)

pdf("output/GO/enrichR_KEGG_2021_Human_THOR_H3K27me3_q30_Lost_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=2)
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=3)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K4me3",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "KEGG_2021_Human") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X") # 

### GeneSymbol list of signif up/down genes in each genotypes
output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram194.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram355.txt

output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol_Venndiagram462.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol_Venndiagram636.txt


output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram220.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt

### GeneSymbol list of signif gain/lost H3K27me3 in each genotypes and DEGs (Carolina RNAseq)
output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt

output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt

### negative control 202 random genes
meta/ENCFF159KBI_geneSymbol_202random_1.bed
meta/ENCFF159KBI_geneSymbol_202random_2.bed
meta/ENCFF159KBI_geneSymbol_202random_3.bed
meta/ENCFF159KBI_geneSymbol_202random_4.bed
meta/ENCFF159KBI_geneSymbol_202random_5.bed
meta/ENCFF159KBI_geneSymbol_202random_6.bed
meta/ENCFF159KBI_geneSymbol_202random_7.bed
meta/ENCFF159KBI_geneSymbol_202random_8.bed
meta/ENCFF159KBI_geneSymbol_202random_9.bed
meta/ENCFF159KBI_geneSymbol_202random_10.bed

meta/ENCFF159KBI_geneSymbol_500random_1.bed
meta/ENCFF159KBI_geneSymbol_500random_2.bed
meta/ENCFF159KBI_geneSymbol_500random_3.bed
meta/ENCFF159KBI_geneSymbol_500random_4.bed
meta/ENCFF159KBI_geneSymbol_500random_5.bed
meta/ENCFF159KBI_geneSymbol_500random_6.bed
meta/ENCFF159KBI_geneSymbol_500random_7.bed
meta/ENCFF159KBI_geneSymbol_500random_8.bed
meta/ENCFF159KBI_geneSymbol_500random_9.bed
meta/ENCFF159KBI_geneSymbol_500random_10.bed


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X`
down <- edown$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X`
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)
# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string
# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics


pdf("output/GO/enrichR_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.pdf", width=8, height=3)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K4me3",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)


### TF plot like JC paper
#### import all TF genes
TF = read.csv("output/GO/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_TFonly.txt") %>%
                               as_tibble()

#### re do gos without filtering for pvalue
up <- eup$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X`
down <- edown$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X`
up$type <- "up"
down$type <- "down"
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)
gos <- rbind(down, up)

### combine TF with gos
gos_TF = TF %>%
  left_join(as_tibble(gos) ) %>%
  filter(type == "down") 

gos_TF_tidy = gos_TF %>%
  separate(Term, into = c("TF", "db"))


### plot


pdf("output/GO/VolcanoPlotTF_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT_down.pdf", width=5, height=4)
ggplot(gos_TF_tidy, aes(x = Odds.Ratio, y = -logAdjP, color = db)) +
  geom_point(aes(color = ifelse(-logAdjP < 1.3, "not signif.", db))) +
  scale_color_manual(values = c("lightgreen", "darkgreen", "grey")) + # Replace with your actual colors
  theme_bw() +
  labs(x = "Odds Ratio", y = "-log10(adjusted p-value)") +
  geom_text_repel(data = subset(gos_TF_tidy, logAdjP < 1.3),
                  aes(label = TF),
                  nudge_x = 0.2,  # Adjust this value to nudge labels to the right
                  size = 3,
                  max.overlaps = 25)  +
  guides(color = guide_legend(override.aes = list(label = ""))) # just to remove the "a" added in fig legend
dev.off()


pdf("output/GO/VolcanoPlotTF_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_THOR_H3K27me3_q30_Gain_q05FC05_NPC_KO_vs_NPC_WT_down_posterCHOP1.pdf", width=5, height=4)
ggplot(gos_TF_tidy, aes(x = Odds.Ratio, y = -logAdjP, color = db)) +
  geom_point(aes(color = ifelse(-logAdjP < 1.3, "not signif.", db))) +
  scale_color_manual(values = c("blue", "lightblue", "grey")) + # Replace with your actual colors
  theme_bw() +
  labs(x = "Odds Ratio", y = "-log10(adjusted p-value)") +
  geom_text_repel(data = subset(gos_TF_tidy, logAdjP < 1.3),
                  aes(label = TF),
                  nudge_x = 0.2,  # Adjust this value to nudge labels to the right
                  size = 6,
                  max.overlaps = 30)  +
  guides(color = guide_legend(override.aes = list(label = ""))) # just to remove the "a" added in fig legend
dev.off()

```

--> `ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X` db identified *EZH2* and *SUZ12* as genes upreg in KO
----> plot `x = odd.ratio` and `y = logadjPval` like JC done after.


Let's add a negative control for ENCODE_ChEA analysis; by selecting 202 random genes and see if we find EZH2 SUZ12 too:

```bash

awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_1.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_2.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_3.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_4.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_5.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_6.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_7.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_8.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_9.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 202 > meta/ENCFF159KBI_geneSymbol_202random_10.bed

awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_1.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_2.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_3.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_4.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_5.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_6.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_7.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_8.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_9.bed
awk '!seen[$4]++ {print $4}' meta/ENCFF159KBI_geneSymbol.bed | shuf -n 500 > meta/ENCFF159KBI_geneSymbol_500random_10.bed

```

# Functional analysis with enrichGO (single list of genes dotplot)


We will use clusterProfile package. Tutorial [here](https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html).

Let's do a test of the pipeline with genes from cluster4 amd cluster14 from the rlog counts. Our background list will be all genes tested for differential expression.

**IMPORTANT NOTE: When doing GO, do NOT set a universe (background list of genes) it perform better!**

```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")

## Read GTF file
gtf_file <- "../../Master/meta/ENCFF159KBI.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


# Genes that gain H3K27me3 in NPC (009)
## Files
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt

gain_H3K27me3_KO = read_csv("output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt", col_names = "gene_name")
  

ego <- enrichGO(gene = as.character(gain_H3K27me3_KO$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_BP_annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836_top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP_annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836_top10.pdf", width=6, height=5)
dotplot(ego, showCategory=10, font.size = 15)
dev.off()

pdf("output/GO/dotplot_BP_annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836_top5.pdf", width=7, height=3)
dotplot(ego, showCategory=5) 
dev.off()

pdf("output/GO/dotplot_BP_annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836_top5v2.pdf", width=10, height=10)
dotplot(ego, showCategory=5) 
dev.off()



# Genes that gain H3K27me3 in NPC (009) and that are downregulated RNAseq (001)
## Files
output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt # 202 genes

gain_H3K27me3_downreg_KO = read_csv("output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", col_names = "gene_name")
  

ego <- enrichGO(gene = as.character(gain_H3K27me3_downreg_KO$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_BP_downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP_downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_top10.pdf", width=6, height=5)
dotplot(ego, showCategory=10, font.size = 15)
dev.off()

pdf("output/GO/dotplot_BP_downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_top5.pdf", width=7, height=3)
dotplot(ego, showCategory=5) 
dev.off()

pdf("output/GO/dotplot_BP_downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_top5v2.pdf", width=6, height=3)
dotplot(ego, showCategory=5, font.size = 15) 
dev.off()


```






# Does the peak that gain H3K27me3 in KO is because of EZH2 binding in NPC?

## Check EZH2 binding profile in NPC (from 005__CutRun)
Check bigwig files of H3K27me3 (`009__CutRun`; THOR q30) and EZH2 (`005__CutRun`; THOR qval10) in WT and KO for the peaks that gain and lost H3K27me3 (2 separate bed). If EZH2 binding increase in the peak that gain H3K27me3 then it strenghten EZH2 compensation in EZH1 KO.

```bash
conda activate deeptools

# PEAKS
## Using classic Recirpocl DiffBind TMM method
sbatch scripts/matrix_TSS_10kb_H3K27me3_EZH2_median_THOR_gainLost.sh # 16222103 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_EZH2_SUZ12_median_THOR_gainLost.sh # 16222113 ok

## Using DiffBind MG1655 bam scaling factor (see 005); use bam MG1655 in DiffBind and Library seq normalization
sbatch scripts/matrix_TSS_10kb_H3K27me3_THORDiffBindTMM_EZH2SUZ12_THORLIBspikein_gainLost.sh # 16223845 ok FAIL because H3K27me3 data with TMM (classic) and PRC2 subunits with LIB norm.
sbatch scripts/matrix_TSS_10kb_H3K27me3_THORDiffBindTMM_EZH2SUZ12_THORLIBspikein_gainLost_corr.sh # 17139835 ok

## Using Default THOR TMM norm for EZH2 and SUZ12, but DiffBindTMM reciprocal THOR for H3K27me3
sbatch scripts/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.sh # 19142544 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.sh # 19142545 ok

sbatch scripts/matrix_TSS_25kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.sh # 19383876 ok
sbatch scripts/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.sh # 19383961 ok




# GENES
## Using classic Recirpocl DiffBind TMM method
sbatch scripts/matrix_TSS_10kb_H3K27me3_EZH2_SUZ12_median_THOR_gainLost_gene.sh # 17197231 ok


```

--> No changes of EZH2 and SUZ12 in peak that gain H3K27me3 in KO... But only 1 replicate... Also looking at all genes strong decrease of EZH2 and SUZ12 binding in KO, here, same level... Another replicate may show changes!!!

--> Also all peaks that gain H3K27me3 in KO, seems already bound by EZH2 and SUZ12 in the WT! So maybe there is an increase EZH2 activity upon EZH1 KO, rather than an increase binding

--> Showing gene or peaks lead to same result (EZH2, SUZ12, barely follow H3K27me3 WTvsKO pattern)

--> `*THORLIBspikein*` norm is weird as show always gain of EZH2/SUZ12, even in regions that lose H3K27me3...

--> Seems KO lead to more H3K27me3 sharp peak; maybe less H3K27me3 spreading in KO, or just less noise/bakground --> increase window to 25/50kb to see if signal disapear: still signal; so likely **more noise in WT vs KO**




## Collect promoter sequence and run MEME

Collect fasta of promoter regions 1kb upstream and 250bp downstream TSS of genes that gain H3K27me3 and are downregulated in KO:
- convert all genes gtf to bed (already done here: `/scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol.bed`)
- change size of each genes to -1kb to +250bp around TSS
- convert bed to fasta
- go in MEME-suite to do MEME motif discovery and TOMTOM (check if motif ressemble known TF)

```bash
# collect promoter region of each genes (-1kb to +250bp)
awk 'BEGIN{OFS="\t"} {
  if($5 == "+") {
    start = ($2 - 1000 < 0) ? 0 : $2 - 1000; # Prevent negative start for forward strand
    end = $2 + 250;
  } else {
    start = $3 - 250;
    end = $3 + 1000;
  }
  print $1, start, end, $4, $5;
}' /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol.bed > /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol_prom1kb250bp.bed
# collect promoter region of each genes (-2.5kb to +2.5kb)
awk 'BEGIN{OFS="\t"} {
  if($5 == "+") {
    start = ($2 - 2500 < 0) ? 0 : $2 - 2500; # Prevent negative start for forward strand
    end = $2 + 2500;
  } else {
    start = $3 - 2500;
    end = $3 + 2500;
  }
  print $1, start, end, $4, $5;
}' /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol.bed > /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol_prom2500bp2500bp.bed


# filter the ENCFF159KBI_geneSymbol_prom1kb250bp.bed to only keep the genes that gain H3K27me3 in KO and are downregulated
grep -Fwf output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol_prom1kb250bp.bed > output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_prom1kb250bp.bed

# convert bed to fasta
conda activate BedToBigwig

bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>

bedtools getfasta -fi /scr1/users/roulet/Akizu_Lab/Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_prom1kb250bp.bed -fo output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_prom1kb250bp.fasta 
```


--> `downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30_prom1kb250bp.fasta` run into [MEME](https://meme-suite.org/meme/tools/meme); in `Classic/ANR`






# THOR with ChIPseeker

Generate xls file with gene name and THOR peak metrics (Naiara Slack task 20240314):
- import THOR bed output (correct qvalue) `output/THOR/THOR*/THOR_qval*.bed` (q30)
- import ChIPseeker output `output/ChIPseeker/annotation_THOR*.txt`
- combine THOR and ChIPseeker with Peak name
- output all metrics in `THOR` folder! Then filter them in xls to keep only the relevant ones


```R
library("tidyverse")

# import files
THOR = read_tsv("output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character())) %>%
       rename("name" = "X4", "qval" = "X15", "FC" = "X16", "WT_count_1" = "X11", "WT_count_2" = "X12", "KO_count_1" = "X13", "KO_count_2" = "X14") %>%
       mutate(WT_count = (WT_count_1)+(WT_count_2) / 2,
              KO_count = (KO_count_1)+(KO_count_2) / 2) %>%
       dplyr::select(name, WT_count, KO_count, FC, qval)

chipseeker_gain = read_tsv("output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos.txt",
                      col_names = TRUE, trim_ws = TRUE) %>%
       dplyr::select(seqnames, start, end, width, name, annotation, geneChr, geneStart, geneEnd, geneLength, geneStrand, geneId, transcriptId, distanceToTSS, geneSymbol, gene) %>%
       add_column(peak = "gain")
chipseeker_lost = read_tsv("output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg.txt",
                      col_names = TRUE, trim_ws = TRUE) %>%
       dplyr::select(seqnames, start, end, width, name, annotation, geneChr, geneStart, geneEnd, geneLength, geneStrand, geneId, transcriptId, distanceToTSS, geneSymbol, gene) %>%
       add_column(peak = "lost")
chipseeker = chipseeker_gain %>%
    bind_rows(chipseeker_lost)
# combine and filter

THOR_chipseeker = THOR %>% left_join(chipseeker)


write.table(THOR_chipseeker, file="output/THOR/THOR_NPC_WTvsKO_H3K27me3/THORq30_chipseeker-NPC_WTvsKO.txt", sep="\t", quote=FALSE, row.names=FALSE)



```





# Pilot grant 20240204

Patient1 has EZH1 GOF mutation + KMT2A (H3K4 methyltransferase) KO
Hypothesis: Gene that gain H3K27me3 in patient1 are because loss of H3K4m3 so more room for H3K27me3 expansion

Isolate bivalent genes in WT: 
- Macs2 peak calling in the 2 pool replicate (qval 2.3) 
- ChIPseeker gene peak assignment
- Filter –in promoter/TSS peaks `annotation_macs2_H3K*me3_WT_pool_qval2.30103_promoterAnd5_geneSymbol`
- Venn diagram of peak enriched genes


**Metrics qval 2.3**:
- WT_H3K27me3:
    - 10,479 peaks
    - 4,877 genes
- WT_H3K4me3
    - 12,248 peaks
    - 11,013 genes
--> 2,138 bivalent genes

**Metrics qval 4**:
- WT_H3K27me3:
    - 5,016 peaks
    - 3,433 genes
- WT_H3K4me3
    - 11,017 peaks
    - 10,177 genes

--> 1,016 bivalent genes


## Generate deeptool plots for the H3K4me3, H3K27me3 and bivalent genes


```bash
# Generate gtf file from gene list:

### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3only.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3only_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3only.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3only_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K4me3only.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K4me3only_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3only.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3only_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3andH3K4me3.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3andH3K4me3_as_gtf_geneSymbol.txt

## Filter the gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3only_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3only.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3only_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3only.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3.gtf

grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K4me3only_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K4me3only.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3only_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3only.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3andH3K4me3_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3andH3K4me3.gtf



# deeptool plots
## qval 2.3
sbatch scripts/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3.sh # 17478801 ok
sbatch scripts/matrix_TSS_10kb_Venn_overlap_WT_H3K27me3H3K4me3.sh # 17478808 ok
## qval 4
sbatch scripts/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4.sh # 17490336 ok


```


--> The Venn overlap gene filtering worked great; H3K4me3/H3K27me3/bivalent genes are clearly identified at q2.3. But some genes show no signal, even when decreasing the zMax scale; let's try more stringeant qvalues
----> XXX With more stringeat qvalues XXX

--> Optimal qvalue is XXX




## enrichR functional analysis

On the **bivalent** - H3K4me3 and H3K27me3 (NPC) at qval 2.3 and 4 with the g**enes that gain H3K27me3 in GOF/HET** (`003__CutRun`)




--> Below code modified to show only 1 set of genes (no up and down, only 1 set)


```R
# packages
library("tidyverse")
library("enrichR")


# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 

### GeneSymbol list of DEGs per tissue
output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt
output/ChIPseeker/Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.txt


# IF starting with geneSymbol

## Read and preprocess data for DEGs genes
gene_names_up <- read.csv("output/ChIPseeker/Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
up$type <- "up"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 50) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)

# Combine the two dataframes
gos <- up
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
new_order <- up_pathways
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.pdf", width=20, height=12)
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.pdf", width=20, height=10)
ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Purple")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.txt", sep="\t", row.names=FALSE, quote=FALSE)



# Define databases for enrichment
dbs <- c("KEGG_2021_Human") # 

### GeneSymbol list of DEGs per tissue
output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt
output/ChIPseeker/Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.txt


# IF starting with geneSymbol

## Read and preprocess data for DEGs genes
gene_names_up <- read.csv("output/ChIPseeker/Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
up$type <- "up"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 50) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)

# Combine the two dataframes
gos <- up
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
new_order <- up_pathways
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_KEGG_2021_Human_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.pdf", width=20, height=12)
pdf("output/GO/enrichR_KEGG_2021_Human_Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.pdf", width=20, height=10)
ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Purple")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__bivalentOnly.txt", sep="\t", row.names=FALSE, quote=FALSE)



```




## Expression of bivalent genes



Collect all bivalent genes (qval2.3 and qval4) and check their log2FC expression level in WT vs GOF in NPC (`001*/001__RNAseq`)+ Highlight in red the bivalent and that gain H3K27me3 in HET



```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

# import RNAseq data WT vs HET NPC (001) #############################################################

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_HET_vs_WT", type="apeglm")


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols




# import bivalent genes ############################################
## gain H3K27me3 in HET 8wN
output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15_gain_Promoter5_geneSymbol.txt
## bivalent only
output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3.txt
output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3andH3K4me3.txt
## bivalent and that gain H3K27me3
output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__WTbivalentand003GainHETTHORq15.txt
output/ChIPseeker/Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__WTbivalentmacs2qval4and003GainHETTHORq15.txt



bivalent_qval2 = read.table("output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()
bivalent_qval4 = read.table("output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4__H3K27me3andH3K4me3.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()

bivalent_qval2_gainH3K27me3HET = read.table("output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__WTbivalentand003GainHETTHORq15.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()
bivalent_qval4_gainH3K27me3HET = read.table("output/ChIPseeker/Venn_overlap_WTbivalentmacs2qval4_003GainHETTHORq15__WTbivalentmacs2qval4and003GainHETTHORq15.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()


## isolate from res (all DEGs genes) the bivalent genes ###############################################
bivalent_qval2
bivalent_qval4
bivalent_qval2_gainH3K27me3HET
bivalent_qval4_gainH3K27me3HET
#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene")

res_bivalent_qval2 = bivalent_qval2 %>% 
  inner_join(res_tibble)
res_bivalent_qval4 = bivalent_qval4 %>% 
  inner_join(res_tibble)
res_bivalent_qval2_gainH3K27me3HET = bivalent_qval2_gainH3K27me3HET %>% 
  inner_join(res_tibble)
res_bivalent_qval4_gainH3K27me3HET = bivalent_qval4_gainH3K27me3HET %>% 
  inner_join(res_tibble)

### highlight genes that gain H3K27me3 in HET ###############################################
gain_H3K27me3_HET = read.table("../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15_gain_Promoter5_geneSymbol.txt", 
                                           header = TRUE) %>%
                               as_tibble() %>%
    dplyr::select(GeneSymbol) %>%
    unique()

highlight_genes <- gain_H3K27me3_HET$GeneSymbol # gain H3K27me3 in HET 8wN (003__CutRun)

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_bivalent_qval2_gainH3K27me3HET$log2FoldChange < -0.5 & res_bivalent_qval2_gainH3K27me3HET$padj < 5e-2, 'Sky Blue',     ######### CAHNGE NAME HERE
    ifelse(res_bivalent_qval2_gainH3K27me3HET$log2FoldChange > 0.5 & res_bivalent_qval2_gainH3K27me3HET$padj < 5e-2, 'Orange',   ######### CAHNGE NAME HERE
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'






pdf("../001__RNAseq/output/deseq2_hg38/plotVolcano_res_bivalent_qval2.pdf", width=8, height=8)  
pdf("../001__RNAseq/output/deseq2_hg38/plotVolcano_res_bivalent_qval4.pdf", width=8, height=8)  
pdf("../001__RNAseq/output/deseq2_hg38/plotVolcano_res_bivalent_qval2_gainH3K27me3HET.pdf", width=8, height=8)  

EnhancedVolcano(res_bivalent_qval2_gainH3K27me3HET,            ######### CAHNGE NAME HERE
  lab = NA,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, NPC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 3,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()




```

--> The code seems bugy, too many genes are shown! They should be around 1k bivalent genes and the plot mention there is 10k genes. The gene name converions may lead to the issue!





# RNAseq integration 

**IMPORTANT NOTE: Here it is advisable to REMOVE all genes from chromosome X and Y BEFORE doing the DEGs analysis (X chromosome re-activation occurs in some samples, notably these with more cell passage; in our case, the HET and KO)**
--> It is good to do this on the count matrix see [here](https://support.bioconductor.org/p/119932/)
### 'one-by-one' comparison
Comparison WT vs mutant (KO or HET) for each time-points:
- NPC KO vs WT

--> DEGs is redo as in `001__RNAseq`


### NPC KO vs WT
Take ressource
```bash
srun --mem=100g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols



## import gene list gain / lost H3K27me3

# H3K27me3
### GAIN

H3K27me3_qval30_Gain = read.table("../009__integration_NPC_WTKO_K27me3K4me3_005_008/output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = H3K27me3_qval30_Gain %>% 
  left_join(res_tibble) 


### LOST
H3K27me3_qval30_Lost = read.table("../009__integration_NPC_WTKO_K27me3K4me3_005_008/output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()



#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = H3K27me3_qval30_Lost %>% 
  left_join(res_tibble)

## PLOT
### GAIN
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

pdf("output/deseq2/plotVolcano_res_001009_Gain_H3K27me3_qval30_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 202
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 73

# Save as gene list for GO analysis:

upregulated <- res_Gain[res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, ]
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, ]
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Gain_H3K27me3_qval30.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### LOST
highlight_genes <- c("") # NA

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.5)'

pdf("output/deseq2/plotVolcano_res_001009_Lost_H3K27me3_qval30_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Lost[res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, ]
upregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Lost[res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, ]
downregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05FC05_NPC_KO_vs_NPC_WT_001009_Lost_H3K27me3_qval30.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


```


--> Functional analysis of Gain/Lost and DEGs done at `# Functional analysis with enrichR`


# Assess H3K27me3 spreading around PRC2 peak (SUZ12, EZH2)

- Collect EZH2 or SUZ12 peaks in WT (`005__*/macs2/*1.3`) overlapping with H3K27me3 (`009__*/macs2/*2.3`)
- Generate bed with regions of 1kb-1.2kb-1.5-2kb length up and downstream of EZH2/SUZ12 peaks --> `bedtools flank`
- Quantify H3K27me3 read density in the peak (convert bigwig to bedgraph) --> `bedtools`
- Quantify H3K27me3 read density 1.2kb up the peak; down the peak --> Try [bin.bw](https://rdrr.io/github/jmonlong/PopSV/man/bin.bw.html) 
- Represent data in R `ggplot`


--> Option1: Do same in KO and compare

--> Option2: Use peak in WT and check signal in WT vs KO





```bash
conda activate BedToBigwig

# Bigwig to bedGraph

## convert bigwig to bedGraph (File already there; with median file generation); bigWigToBedGraph in.bigwig out.bedGraph
output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.sorted.bedGraph # WT median
output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.sorted.bedGraph # KO median

# macs2 H3K27me3 peaks
output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_peaks.broadPeak # WT EZH2 H3K27me3 peak
output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_peaks.broadPeak # KO EZH2 H3K27me3 peak

# macs2 EZH2 peaks overlapping with H3K27me3
../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks.broadPeak # WT peak
../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks.broadPeak # KO peak
## keep only EZH2 peak overlapping
bedtools intersect -wa -a ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks.broadPeak -b output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_pool_peaks.broadPeak > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak
bedtools intersect -wa -a ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks.broadPeak -b output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_pool_peaks.broadPeak > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak




## collect flanking regions (bedtools flank [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> [-b or (-l and -r)])
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kb.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1200 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1.2kb.broadPeak

bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank2kb.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1500 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank1.5kb.broadPeak



## collect flanking regions - separating upstream / downstream (see NOTE at the end) (-l upstream -r downstream)
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 2000 -r 0 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbupstream.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 0 -r 1000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbdownstream.broadPeak

bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 2000 -r 0 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank2kbupstream.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 0 -r 2000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank2kbdownstream.broadPeak






# macs2 SUZ12 peaks

output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_peaks.broadPeak # WT EZH2 H3K27me3 peak
output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_peaks.broadPeak # KO EZH2 H3K27me3 peak

# macs2 SUZ12 peaks overlapping with H3K27me3
../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks.broadPeak # WT peak
../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks.broadPeak # KO peak
## keep only SUZ12 peak overlapping
bedtools intersect -wa -a ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks.broadPeak -b output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_pool_peaks.broadPeak > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak
bedtools intersect -wa -a ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks.broadPeak -b output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_pool_peaks.broadPeak > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak




## collect flanking regions (bedtools flank [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> [-b or (-l and -r)])
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kb.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1500 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1.5kb.broadPeak

bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1kb.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1200 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1.2kb.broadPeak

## collect flanking regions - separating upstream / downstream (see NOTE at the end) (-l upstream -r downstream)
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 2000 -r 0 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbupstream.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 0 -r 2000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbdownstream.broadPeak

bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 1000 -r 0 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1kbupstream.broadPeak
bedtools flank -i ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak -g ../../Master/meta/GRCh38_chrom_sizes.tab -l 0 -r 1000 > ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1kbdownstream.broadPeak



```


## Set up bin.bw
Create conda env to use bin.bw (copy one with devtool installed; `ChIPseqSpikeInFree` )

```bash
conda create --name binBw --clone ChIPseqSpikeInFree
conda activate binBw
# --> FAIL, try with scRNAseq env

conda create --name binBw_v2 --clone scRNAseq
conda activate binBw_v2
```

Go in R and install bin.bw

```R
BiocManager::install("DNAcopy")
# Install PopSV from GitHub
devtools::install_github("jmonlong/PopSV")
```


--> WORK!!

## Use bin.bw


**Option1: Do same in KO and compare**



```bash
conda activate binBw_v2
```

```R
library("PopSV")
library("tidyverse")
library("Rsamtools")
library("ggpubr")
library("data.table")

# EZH2 #####################################################

## WT EZH2_in Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end")
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2")
# WT EZH2_ upstream / downtream - 1kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2upstream1kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2downstream1kb") #


# WT EZH2_ upstream / downtream - 2kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2upstream2kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2downstream2kb") #




## KO EZH2_in Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end") %>% unique()
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2")
## KO EZH2_ upstream / downtream - 1kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2upstream1kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2downstream1kb") # 


## KO EZH2_ upstream / downtream - 2kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2upstream2kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2downstream2kb") # 



## Boxplot - downstream peak upstream
# WT 1 kb EZH2
### WT_EZH2 peak
WT_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_EZH2.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  WT_EZH2 upstream / downstream 1kb
WT_EZH2upstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2upstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_EZH2downstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2downstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_EZH2_tidy_1kb = WT_EZH2 %>%
  bind_rows(WT_EZH2upstream1kb) %>%
  bind_rows(WT_EZH2downstream1kb) %>%
  mutate(bc_norm = bc / length) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/WT_EZH2flank1kb.pdf", width=4, height=4)  
ggplot(WT_EZH2_tidy_1kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  WT_EZH2 upstream / downstream 2kb
WT_EZH2upstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2upstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_EZH2downstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2downstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_EZH2_tidy_2kb = WT_EZH2 %>%
  bind_rows(WT_EZH2upstream2kb) %>%
  bind_rows(WT_EZH2downstream2kb) %>%
  mutate(bc_norm = bc / length) 

pdf("output/binBw/WT_EZH2flank2kb.pdf", width=4, height=4)  
ggplot(WT_EZH2_tidy_2kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream"))) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))
dev.off()


### KO_EZH2
############### slight test
# KO_EZH2_with_direction = as_tibble(fread(cmd = "gunzip -c output/binBw/KO_EZH2.bgz") ) %>%
#  add_column(direction = "peak") 
# write.table(KO_EZH2_with_direction, file = "output/binBw/KO_EZH2_with_direction.txt", sep = "\t", row.names = FALSE, quote = FALSE)
################

# KO 1 kb EZH2
### KO_EZH2 peak
KO_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/KO_EZH2.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  KO_EZH2 upstream / downstream 1kb
KO_EZH2upstream1kb <- fread(cmd = "gunzip -c output/binBw/KO_EZH2upstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_EZH2downstream1kb <- fread(cmd = "gunzip -c output/binBw/KO_EZH2downstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_EZH2_tidy_1kb = KO_EZH2 %>%
  bind_rows(KO_EZH2upstream1kb) %>%
  bind_rows(KO_EZH2downstream1kb) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_EZH2flank1kb.pdf", width=4, height=4)  
ggplot(KO_EZH2_tidy_1kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  KO_EZH2 upstream / downstream 2kb
KO_EZH2upstream2kb <- fread(cmd = "gunzip -c output/binBw/KO_EZH2upstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_EZH2downstream2kb <- fread(cmd = "gunzip -c output/binBw/KO_EZH2downstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_EZH2_tidy_2kb = KO_EZH2 %>%
  bind_rows(KO_EZH2upstream2kb) %>%
  bind_rows(KO_EZH2downstream2kb) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_EZH2flank2kb.pdf", width=4, height=4)  
ggplot(KO_EZH2_tidy_2kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()



## Boxplot - ratio Inside (peak) vs Outside (5 or 3')
## 1kb WT vs KO
WT_EZH2_tidy_1kb_ratio = WT_EZH2_tidy_1kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_EZH2_tidy_1kb_ratio = KO_EZH2_tidy_1kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

EZH2_tidy_1kb_ratio = WT_EZH2_tidy_1kb_ratio %>%
  bind_rows(KO_EZH2_tidy_1kb_ratio) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_EZH2upstream1kb_ratio.pdf", width=3, height=3)  
ggplot(EZH2_tidy_1kb_ratio, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()


pdf("output/binBw/WTvsKO_EZH2downstream1kb_ratio.pdf", width=3, height=3)  
ggplot(EZH2_tidy_1kb_ratio, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()





## 2kb WT vs KO
WT_EZH2_tidy_2kb_ratio = WT_EZH2_tidy_2kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_EZH2_tidy_2kb_ratio = KO_EZH2_tidy_2kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

EZH2_tidy_2kb_ratio = WT_EZH2_tidy_2kb_ratio %>%
  bind_rows(KO_EZH2_tidy_2kb_ratio) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_EZH2upstream2kb_ratio.pdf", width=3, height=3)  
ggplot(EZH2_tidy_2kb_ratio, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()



pdf("output/binBw/WTvsKO_EZH2downstream2kb_ratio.pdf", width=3, height=3)  
ggplot(EZH2_tidy_2kb_ratio, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()




# SUZ12 #####################################################

## WT SUZ12_in Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end")
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12")
# WT SUZ12_ upstream / downtream - 1kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12upstream1kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12downstream1kb") #


# WT SUZ12_ upstream / downtream - 2kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12upstream2kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12downstream2kb") #




## KO SUZ12_in Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end") %>% unique()
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12")
## KO SUZ12_ upstream / downtream - 1kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12upstream1kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12downstream1kb") # 


## KO SUZ12_ upstream / downtream - 2kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12upstream2kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12downstream2kb") # 



## Boxplot - downstream peak upstream
# WT 1 kb SUZ12
### WT_SUZ12 peak
WT_SUZ12 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_SUZ12.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  WT_SUZ12 upstream / downstream 1kb
WT_SUZ12upstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12upstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_SUZ12downstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12downstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_SUZ12_tidy_1kb = WT_SUZ12 %>%
  bind_rows(WT_SUZ12upstream1kb) %>%
  bind_rows(WT_SUZ12downstream1kb) %>%
  mutate(bc_norm = bc / length) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/WT_SUZ12flank1kb.pdf", width=4, height=4)  
ggplot(WT_SUZ12_tidy_1kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  WT_SUZ12 upstream / downstream 2kb
WT_SUZ12upstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12upstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_SUZ12downstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12downstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_SUZ12_tidy_2kb = WT_SUZ12 %>%
  bind_rows(WT_SUZ12upstream2kb) %>%
  bind_rows(WT_SUZ12downstream2kb) %>%
  mutate(bc_norm = bc / length) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/WT_SUZ12flank2kb.pdf", width=4, height=4)  
ggplot(WT_SUZ12_tidy_2kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()


### KO_SUZ12
############### slight test
# KO_SUZ12_with_direction = as_tibble(fread(cmd = "gunzip -c output/binBw/KO_SUZ12.bgz") ) %>%
#  add_column(direction = "peak") 
# write.table(KO_SUZ12_with_direction, file = "output/binBw/KO_SUZ12_with_direction.txt", sep = "\t", row.names = FALSE, quote = FALSE)
################

# KO 1 kb SUZ12
### KO_SUZ12 peak
KO_SUZ12 <- as_tibble(fread(cmd = "gunzip -c output/binBw/KO_SUZ12.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_overlap_001009_NPC_KO_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  KO_SUZ12 upstream / downstream 1kb
KO_SUZ12upstream1kb <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12upstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_SUZ12downstream1kb <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12downstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_SUZ12_tidy_1kb = KO_SUZ12 %>%
  bind_rows(KO_SUZ12upstream1kb) %>%
  bind_rows(KO_SUZ12downstream1kb) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_SUZ12flank1kb.pdf", width=4, height=4)  
ggplot(KO_SUZ12_tidy_1kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  KO_SUZ12 upstream / downstream 2kb
KO_SUZ12upstream2kb <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12upstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_SUZ12downstream2kb <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12downstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_KO_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_SUZ12_tidy_2kb = KO_SUZ12 %>%
  bind_rows(KO_SUZ12upstream2kb) %>%
  bind_rows(KO_SUZ12downstream2kb) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_SUZ12flank2kb.pdf", width=4, height=4)  
ggplot(KO_SUZ12_tidy_2kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()



## Boxplot - ratio Inside (peak) vs Outside (5 or 3')
## 1kb WT vs KO
WT_SUZ12_tidy_1kb_ratio = WT_SUZ12_tidy_1kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_SUZ12_tidy_1kb_ratio = KO_SUZ12_tidy_1kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

SUZ12_tidy_1kb_ratio = WT_SUZ12_tidy_1kb_ratio %>%
  bind_rows(KO_SUZ12_tidy_1kb_ratio) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_SUZ12upstream1kb_ratio.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_1kb_ratio, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "red3")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()


pdf("output/binBw/WTvsKO_SUZ12downstream1kb_ratio.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_1kb_ratio, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "red3")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()





## 2kb WT vs KO
WT_SUZ12_tidy_2kb_ratio = WT_SUZ12_tidy_2kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_SUZ12_tidy_2kb_ratio = KO_SUZ12_tidy_2kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

SUZ12_tidy_2kb_ratio = WT_SUZ12_tidy_2kb_ratio %>%
  bind_rows(KO_SUZ12_tidy_2kb_ratio) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_SUZ12upstream2kb_ratio.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_2kb_ratio, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "red3")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()



pdf("output/binBw/WTvsKO_SUZ12downstream2kb_ratio.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_2kb_ratio, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "red3")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()

```

- *NOTE: for `output/binBw/WT_EZH2flank1kb.bgz` I initially indicated whether upstream and downstream the peak with `direction`; the method I used was to say row 1 = upstream and row 2 = downstream, etc... BUT that lead to issue when peaks are close together!!!! Let's better generate two separate file; one for downstream; one for upstream*





**Option2: Use peak in WT and check signal in WT vs KO**



```bash
conda activate binBw_v2
```

```R
library("PopSV")
library("tidyverse")
library("Rsamtools")
library("ggpubr")
library("data.table")

# EZH2 #####################################################

## WT EZH2_in Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end")
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2")
# WT EZH2_ upstream / downtream - 1kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2upstream1kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2downstream1kb") #


# WT EZH2_ upstream / downtream - 2kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2upstream2kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2downstream2kb") #




## KO EZH2_in WT Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end")
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2_WTpeaks")
# KO EZH2_ upstream / downtream - 1kb  in WT Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2upstream1kb_WTpeaks") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2downstream1kb_WTpeaks") #


# KO EZH2_ upstream / downtream - 2kb  in WT Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2upstream2kb_WTpeaks") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_EZH2downstream2kb_WTpeaks") #






## Boxplot - downstream peak upstream
# WT 1 kb EZH2
### WT_EZH2 peak
WT_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_EZH2.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  WT_EZH2 upstream / downstream 1kb
WT_EZH2upstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2upstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_EZH2downstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2downstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_EZH2_tidy_1kb = WT_EZH2 %>%
  bind_rows(WT_EZH2upstream1kb) %>%
  bind_rows(WT_EZH2downstream1kb) %>%
  mutate(bc_norm = bc / length) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))


pdf("output/binBw/WT_EZH2flank1kb.pdf", width=4, height=4)  
ggplot(WT_EZH2_tidy_1kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  WT_EZH2 upstream / downstream 2kb
WT_EZH2upstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2upstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_EZH2downstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_EZH2downstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_EZH2_tidy_2kb = WT_EZH2 %>%
  bind_rows(WT_EZH2upstream2kb) %>%
  bind_rows(WT_EZH2downstream2kb) %>%
  mutate(bc_norm = bc / length) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/WT_EZH2flank2kb.pdf", width=4, height=4)  
ggplot(WT_EZH2_tidy_2kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()


### KO_EZH2
############### slight test
# KO_EZH2_with_direction = as_tibble(fread(cmd = "gunzip -c output/binBw/KO_EZH2.bgz") ) %>%
#  add_column(direction = "peak") 
# write.table(KO_EZH2_with_direction, file = "output/binBw/KO_EZH2_with_direction.txt", sep = "\t", row.names = FALSE, quote = FALSE)
################

# KO 1 kb EZH2
### KO_EZH2 peak
KO_EZH2_WTpeaks <- as_tibble(fread(cmd = "gunzip -c output/binBw/KO_EZH2_WTpeaks.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  KO_EZH2 upstream / downstream 1kb
KO_EZH2upstream1kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_EZH2upstream1kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_EZH2downstream1kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_EZH2downstream1kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_EZH2_tidy_1kb_WTpeaks = KO_EZH2_WTpeaks %>%
  bind_rows(KO_EZH2upstream1kb_WTpeaks) %>%
  bind_rows(KO_EZH2downstream1kb_WTpeaks) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_EZH2flank1kb_WTpeaks.pdf", width=4, height=4)  
ggplot(KO_EZH2_tidy_1kb_WTpeaks, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  KO_EZH2 upstream / downstream 2kb
KO_EZH2upstream2kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_EZH2upstream2kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_EZH2downstream2kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_EZH2downstream2kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_EZH2_tidy_2kb_WTpeaks = KO_EZH2_WTpeaks %>%
  bind_rows(KO_EZH2upstream2kb_WTpeaks) %>%
  bind_rows(KO_EZH2downstream2kb_WTpeaks) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_EZH2flank2kb_WTpeaks.pdf", width=4, height=4)  
ggplot(KO_EZH2_tidy_2kb_WTpeaks, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()



## Boxplot - ratio Inside (peak) vs Outside (5 or 3')
## 1kb WT vs KO
WT_EZH2_tidy_1kb_ratio = WT_EZH2_tidy_1kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_EZH2_tidy_1kb_ratio_WTpeaks = KO_EZH2_tidy_1kb_WTpeaks %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

EZH2_tidy_1kb_ratio_WTpeaks = WT_EZH2_tidy_1kb_ratio %>%
  bind_rows(KO_EZH2_tidy_1kb_ratio_WTpeaks) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_EZH2upstream1kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(EZH2_tidy_1kb_ratio_WTpeaks, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()


pdf("output/binBw/WTvsKO_EZH2downstream1kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(EZH2_tidy_1kb_ratio_WTpeaks, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()





## 2kb WT vs KO
WT_EZH2_tidy_2kb_ratio = WT_EZH2_tidy_2kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_EZH2_tidy_2kb_ratio_WTpeaks = KO_EZH2_tidy_2kb_WTpeaks %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

EZH2_tidy_2kb_ratio_WTpeaks = WT_EZH2_tidy_2kb_ratio %>%
  bind_rows(KO_EZH2_tidy_2kb_ratio_WTpeaks) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_EZH2upstream2kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(EZH2_tidy_2kb_ratio_WTpeaks, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 1.9 ) +
  ylim(0,2.5)
dev.off()



pdf("output/binBw/WTvsKO_EZH2downstream2kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(EZH2_tidy_2kb_ratio_WTpeaks, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 1.9 ) +
  ylim(0,2.5)
dev.off()



## pearson corr ratio
### 1kb
#### upstream
EZH2_tidy_1kb_ratio_WTpeaks

EZH2_tidy_1kb_ratio_WTpeaks_upstream <- EZH2_tidy_1kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_upstream) %>%
  spread(key = genotype, value = ratio_upstream)

# Calculate Pearson correlation
correlation <- cor(EZH2_tidy_1kb_ratio_WTpeaks_upstream$WT, EZH2_tidy_1kb_ratio_WTpeaks_upstream$KO)

# Plot
pdf("output/binBw/scatter_plot_EZH2_1kb_upstream.pdf", width=5, height=5)
ggplot(EZH2_tidy_1kb_ratio_WTpeaks_upstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 4, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 4, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()

#### downstream

EZH2_tidy_1kb_ratio_WTpeaks_downstream <- EZH2_tidy_1kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_downstream) %>%
  spread(key = genotype, value = ratio_downstream)

# Calculate Pearson correlation
correlation <- cor(EZH2_tidy_1kb_ratio_WTpeaks_downstream$WT, EZH2_tidy_1kb_ratio_WTpeaks_downstream$KO)

# Plot
pdf("output/binBw/scatter_plot_EZH2_1kb_downstream.pdf", width=5, height=5)
ggplot(EZH2_tidy_1kb_ratio_WTpeaks_downstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 4, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 4, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()



### 2kb
#### upstream
EZH2_tidy_2kb_ratio_WTpeaks

EZH2_tidy_2kb_ratio_WTpeaks_upstream <- EZH2_tidy_2kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_upstream) %>%
  spread(key = genotype, value = ratio_upstream)

# Calculate Pearson correlation
correlation <- cor(EZH2_tidy_2kb_ratio_WTpeaks_upstream$WT, EZH2_tidy_2kb_ratio_WTpeaks_upstream$KO)

# Plot
pdf("output/binBw/scatter_plot_EZH2_2kb_upstream.pdf", width=5, height=5)
ggplot(EZH2_tidy_2kb_ratio_WTpeaks_upstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 4, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 4, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()

#### downstream

EZH2_tidy_2kb_ratio_WTpeaks_downstream <- EZH2_tidy_2kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_downstream) %>%
  spread(key = genotype, value = ratio_downstream)

# Calculate Pearson correlation
correlation <- cor(EZH2_tidy_2kb_ratio_WTpeaks_downstream$WT, EZH2_tidy_2kb_ratio_WTpeaks_downstream$KO)

# Plot
pdf("output/binBw/scatter_plot_EZH2_2kb_downstream.pdf", width=5, height=5)
ggplot(EZH2_tidy_2kb_ratio_WTpeaks_downstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 4, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 4, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()







# SUZ12 #####################################################

## WT SUZ12_in Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end")
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12")
# WT SUZ12_ upstream / downtream - 1kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12upstream1kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12downstream1kb") #


# WT SUZ12_ upstream / downtream - 2kb
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12upstream2kb") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_SUZ12downstream2kb") #




## KO SUZ12_in WT Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%
  dplyr::select("chr", "start", "end")
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12_WTpeaks")
# KO SUZ12_ upstream / downtream - 1kb  in WT Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12upstream1kb_WTpeaks") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12downstream1kb_WTpeaks") #


# KO SUZ12_ upstream / downtream - 2kb  in WT Peaks
bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() #
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12upstream2kb_WTpeaks") # 

bwFile <- "output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw"
regions <- read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")) %>%   dplyr::select("chr", "start", "end") %>% unique() # 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/KO_SUZ12downstream2kb_WTpeaks") #






## Boxplot - downstream peak upstream
# WT 1 kb SUZ12
### WT_SUZ12 peak
WT_SUZ12 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_SUZ12.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  WT_SUZ12 upstream / downstream 1kb
WT_SUZ12upstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12upstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_SUZ12downstream1kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12downstream1kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_SUZ12_tidy_1kb = WT_SUZ12 %>%
  bind_rows(WT_SUZ12upstream1kb) %>%
  bind_rows(WT_SUZ12downstream1kb) %>%
  mutate(bc_norm = bc / length) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))


pdf("output/binBw/WT_SUZ12flank1kb.pdf", width=4, height=4)  
ggplot(WT_SUZ12_tidy_1kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  WT_SUZ12 upstream / downstream 2kb
WT_SUZ12upstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12upstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
WT_SUZ12downstream2kb <- fread(cmd = "gunzip -c output/binBw/WT_SUZ12downstream2kb.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


WT_SUZ12_tidy_2kb = WT_SUZ12 %>%
  bind_rows(WT_SUZ12upstream2kb) %>%
  bind_rows(WT_SUZ12downstream2kb) %>%
  mutate(bc_norm = bc / length) %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/WT_SUZ12flank2kb.pdf", width=4, height=4)  
ggplot(WT_SUZ12_tidy_2kb, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "lightgrey")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()


### KO_SUZ12
############### slight test
# KO_SUZ12_with_direction = as_tibble(fread(cmd = "gunzip -c output/binBw/KO_SUZ12.bgz") ) %>%
#  add_column(direction = "peak") 
# write.table(KO_SUZ12_with_direction, file = "output/binBw/KO_SUZ12_with_direction.txt", sep = "\t", row.names = FALSE, quote = FALSE)
################

# KO 1 kb SUZ12
### KO_SUZ12 peak
KO_SUZ12_WTpeaks <- as_tibble(fread(cmd = "gunzip -c output/binBw/KO_SUZ12_WTpeaks.bgz") ) %>%
 add_column(direction = "peak") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2"))) %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 

##  KO_SUZ12 upstream / downstream 1kb
KO_SUZ12upstream1kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12upstream1kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_SUZ12downstream1kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12downstream1kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank1kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_SUZ12_tidy_1kb_WTpeaks = KO_SUZ12_WTpeaks %>%
  bind_rows(KO_SUZ12upstream1kb_WTpeaks) %>%
  bind_rows(KO_SUZ12downstream1kb_WTpeaks) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_SUZ12flank1kb_WTpeaks.pdf", width=4, height=4)  
ggplot(KO_SUZ12_tidy_1kb_WTpeaks, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()

##  KO_SUZ12 upstream / downstream 2kb
KO_SUZ12upstream2kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12upstream2kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "upstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbupstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 
KO_SUZ12downstream2kb_WTpeaks <- fread(cmd = "gunzip -c output/binBw/KO_SUZ12downstream2kb_WTpeaks.bgz") %>% as_tibble()  %>%
 add_column(direction = "downstream") %>%
 left_join(read.table("../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_flank2kbdownstream.broadPeak", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand", "trash", "trash1", "trash2")))  %>%
 mutate(length = end - start ) %>%
 dplyr::select("name", "length", "direction", "bc") 


KO_SUZ12_tidy_2kb_WTpeaks = KO_SUZ12_WTpeaks %>%
  bind_rows(KO_SUZ12upstream2kb_WTpeaks) %>%
  bind_rows(KO_SUZ12downstream2kb_WTpeaks) %>%
  mutate(bc_norm = bc / length)  %>%
  unique() %>%
  mutate(direction = factor(direction, levels = c("upstream", "peak", "downstream")))

pdf("output/binBw/KO_SUZ12flank2kb_WTpeaks.pdf", width=4, height=4)  
ggplot(KO_SUZ12_tidy_2kb_WTpeaks, aes(x = direction, y = bc_norm, fill = direction)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "red4") +
  theme_bw() +
  labs(x = NULL, y = "H3K27me3 density") +
  scale_fill_manual(values = c("red2", "red3", "red2")) +
  ggpubr::stat_compare_means(aes(group = direction), label = "p.format",
                                    comparisons = list(c("peak", "upstream"), c("peak", "downstream")))
dev.off()



## Boxplot - ratio Inside (peak) vs Outside (5 or 3')
## 1kb WT vs KO
WT_SUZ12_tidy_1kb_ratio = WT_SUZ12_tidy_1kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_SUZ12_tidy_1kb_ratio_WTpeaks = KO_SUZ12_tidy_1kb_WTpeaks %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

SUZ12_tidy_1kb_ratio_WTpeaks = WT_SUZ12_tidy_1kb_ratio %>%
  bind_rows(KO_SUZ12_tidy_1kb_ratio_WTpeaks) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_SUZ12upstream1kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_1kb_ratio_WTpeaks, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()


pdf("output/binBw/WTvsKO_SUZ12downstream1kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_1kb_ratio_WTpeaks, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 2 ) +
  ylim(0,2.5)
dev.off()





## 2kb WT vs KO
WT_SUZ12_tidy_2kb_ratio = WT_SUZ12_tidy_2kb %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "WT")

KO_SUZ12_tidy_2kb_ratio_WTpeaks = KO_SUZ12_tidy_2kb_WTpeaks %>%
  dplyr::select(name,bc_norm,direction) %>%
  pivot_wider(names_from = direction, values_from = bc_norm, names_prefix = "bc_norm_") %>%
  mutate(ratio_upstream = bc_norm_upstream / bc_norm_peak,
         ratio_downstream = bc_norm_downstream / bc_norm_peak) %>%
  add_column(genotype = "KO")

SUZ12_tidy_2kb_ratio_WTpeaks = WT_SUZ12_tidy_2kb_ratio %>%
  bind_rows(KO_SUZ12_tidy_2kb_ratio_WTpeaks) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")))



pdf("output/binBw/WTvsKO_SUZ12upstream2kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_2kb_ratio_WTpeaks, aes(x = genotype, y = ratio_upstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 1.9 ) +
  ylim(0,2.5)
dev.off()



pdf("output/binBw/WTvsKO_SUZ12downstream2kb_ratio_WTpeaks.pdf", width=3, height=3)  
ggplot(SUZ12_tidy_2kb_ratio_WTpeaks, aes(x = genotype, y = ratio_downstream, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
 # geom_jitter(width = 0.15, size = 0.3, alpha = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "Out/In H3K27me3 ratio") +
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  ggpubr::stat_compare_means(aes(group = genotype), label = "p.format",
                                    comparisons = list(c("WT", "KO")),
                                    label.y = 1.9 ) +
  ylim(0,2.5)
dev.off()



## pearson corr ratio
### 1kb
#### upstream
SUZ12_tidy_1kb_ratio_WTpeaks

SUZ12_tidy_1kb_ratio_WTpeaks_upstream <- SUZ12_tidy_1kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_upstream) %>%
  spread(key = genotype, value = ratio_upstream) %>%
  filter(!is.na(WT), !is.na(KO), is.finite(WT), is.finite(KO))  # Ensure there are no NA or infinite values


# Calculate Pearson correlation
correlation <- cor(SUZ12_tidy_1kb_ratio_WTpeaks_upstream$WT, SUZ12_tidy_1kb_ratio_WTpeaks_upstream$KO)

# Plot
pdf("output/binBw/scatter_plot_SUZ12_1kb_upstream.pdf", width=5, height=5)
ggplot(SUZ12_tidy_1kb_ratio_WTpeaks_upstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 2, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 2, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()

#### downstream

SUZ12_tidy_1kb_ratio_WTpeaks_downstream <- SUZ12_tidy_1kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_downstream) %>%
  spread(key = genotype, value = ratio_downstream)%>%
  filter(!is.na(WT), !is.na(KO), is.finite(WT), is.finite(KO))  # Ensure there are no NA or infinite values

# Calculate Pearson correlation
correlation <- cor(SUZ12_tidy_1kb_ratio_WTpeaks_downstream$WT, SUZ12_tidy_1kb_ratio_WTpeaks_downstream$KO)

# Plot
pdf("output/binBw/scatter_plot_SUZ12_1kb_downstream.pdf", width=5, height=5)
ggplot(SUZ12_tidy_1kb_ratio_WTpeaks_downstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 2, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 2, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()



### 2kb
#### upstream
SUZ12_tidy_2kb_ratio_WTpeaks

SUZ12_tidy_2kb_ratio_WTpeaks_upstream <- SUZ12_tidy_2kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_upstream) %>%
  spread(key = genotype, value = ratio_upstream)%>%
  filter(!is.na(WT), !is.na(KO), is.finite(WT), is.finite(KO))  # Ensure there are no NA or infinite values

# Calculate Pearson correlation
correlation <- cor(SUZ12_tidy_2kb_ratio_WTpeaks_upstream$WT, SUZ12_tidy_2kb_ratio_WTpeaks_upstream$KO)

# Plot
pdf("output/binBw/scatter_plot_SUZ12_2kb_upstream.pdf", width=5, height=5)
ggplot(SUZ12_tidy_2kb_ratio_WTpeaks_upstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 3, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 3, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()

#### downstream

SUZ12_tidy_2kb_ratio_WTpeaks_downstream <- SUZ12_tidy_2kb_ratio_WTpeaks %>%
  filter(genotype %in% c("WT", "KO")) %>%
  select(name, genotype, ratio_downstream) %>%
  spread(key = genotype, value = ratio_downstream)%>%
  filter(!is.na(WT), !is.na(KO), is.finite(WT), is.finite(KO))  # Ensure there are no NA or infinite values

# Calculate Pearson correlation
correlation <- cor(SUZ12_tidy_2kb_ratio_WTpeaks_downstream$WT, SUZ12_tidy_2kb_ratio_WTpeaks_downstream$KO)

# Plot
pdf("output/binBw/scatter_plot_SUZ12_2kb_downstream.pdf", width=5, height=5)
ggplot(SUZ12_tidy_2kb_ratio_WTpeaks_downstream, aes(x = WT, y = KO)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "Wild type IN/OUT H3K27me3 ratio",
       y = "Ezh1 KO IN/OUT H3K27me3 ratio") +
  theme_bw() +
  annotate("text", x = 3, y = 6, label = paste("r =", round(correlation, 2)), size = 5, hjust = 0) +
  annotate("text", x = 3, y = 5, label = paste("r² =", round(correlation^2, 2)), size = 5, hjust = 0)
dev.off()






```




--> H3K27me3 spreading is affected upon EZH1 KO. Not as much as EZH2 KO as shown in mice, but still a significant decrease of spreading around PRC2 sites (EZH2, SUZ12; in a 1-2kb window)




# Normalization method from Ferguson et al

Summary pipeline, improved after troubleshooting `001*/016*`:
- mapping Bowtie2 only keep uniquely aligned reads
- Convert bam to bigwig with 1bp bins
- Identify local maxima (norm99) and generate blacklist region from the bigiwg
  - let's instead use ENCODE blacklist regions
- Apply SF to bam to bigwig
- smooth the bigiwg to 50bp bins for better vizualization




```bash
conda activate deeptools

# Convert bam to bigwig 
#--> Already done

# Convert bigwig to bedgraph
conda activate BedToBigwig
## Unique bigwig (1bp resolution)
sbatch scripts/BedToBigwig_Ferguson_unique.sh # 38030781 ok


# Remove blacklist regions
conda activate BedToBigwig
## Unique bigwig (1bp resolution)
sbatch scripts/BedintersectBlacklist_Ferguson_unique.sh # 38030863 ok
```

Use Python to identify local maxima, quantify the height for the 99th percentile peak

```bash
srun --mem=250g --pty bash -l

# Identify local maxima
## Unique bigwig (1bp resolution)
python scripts/LocalMaxima_Ferguson_unique.py

#  calculate the 99th percentile of the signal heights (score) in the local maxima files.
## Unique bigwig (1bp resolution)
python scripts/Percentile99_Ferguson_unique.py

# normalize AB per AB (using WT sample 1st replicate as reference)
## Unique bigwig (1bp resolution)
### 99th percentile
python scripts/norm_H3K27me3_Ferguson_Perc99_unique.py
python scripts/norm_H3K4me3_Ferguson_Perc99_unique.py



```
--> Works!

- *NOTE: **Local Maxima** = value is higher than its neighboring points. In the context of your CUT&RUN data (or other genomic data), local maxima refer to genomic positions where the signal intensity (e.g., read depth or coverage in the bedGraph file) is greater than the signal in the surrounding regions.*
- *NOTE: **Percentile 99** = signal level that is greater than 99% of all other signal values in the dataset.*


Convert normalized bedGraph back to bigwig

```bash
conda activate BedToBigwig

# Unique bigwig (1bp resolution)
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique.sh # 38033297 ok


```
--> XXX





## Bigwig Ferguson

### Median tracks

Let's generate median tracks for Ferguson and THOR_Ferguson bigwigs (Norm99 unique)

**Run wiggletools:**
```bash
conda activate BedToBigwig

# Calculate median
## bigwig_Ferguson
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique-H3K27me3.sh # 38035621 ok
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique-H3K4me3.sh # 38035824 ok


## smooth bigwig
### calculate bin signal with multiBigwigSummary
conda activate deeptools

sbatch scripts/bigwigsmooth_Norm99_Ferguson_unique.sh # 38037234 ok
### Re-convert to bigwig
conda activate BedToBigwig

sbatch --dependency=afterany:38037234 scripts/bigwigsmooth_Norm99_Ferguson_unique_part2.sh # 38040571 xxx
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*

-->  Smoothing bigwig using `multiBigwigSummary` work great! 







# Ferguson method Diff binding


Ferguson method:
- Call peaks with SEACR or MACS2
  --> Done with macs2 (broad, default, no qvalue filtering)
- Calculate length-normalize signal for each locus (computeMatrix –scale-regions (gene or peak) )
  --> Not sure how to deal with the peak; so let's instead do per gene promoter (1kb up 250bp down); save in `output/edgeR`
- Diff. log-normalize Counts using edgeR and limma on consensus peak or genes (consensus peak?)


## Diff binding on consensus peak

--> `001*/016*` been copied; Here start directly with **optimal qvalue for macs2 consensus peak** --> 2.3 for both H3K27me3 and H3K4me3




NPC_WT_H3K27me3_005
NPC_WT_H3K27me3_008
NPC_KO_H3K27me3_005
NPC_KO_H3K27me3_008

```bash
conda activate BedToBigwig

# concatenate and sort bed files
## qvalue 2.3 ##############
cat output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval2.30103/NPC_WT_H3K4me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/NPC_KO_H3K4me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.broadPeak




# merge = consensus peak identification
## qvalue 2.3 ##############
### no merge extension
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge.bed

### with 100bp peak merging
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge100bp.bed

### with 500bp peak merging
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge500bp.bed
```

--> All good; consensus peak files are: `output/macs2/broad/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge[SIZE].bed`


Now calculate **signal in consensus peak**:


```bash
conda activate deeptools

### H3K27me3

# sample per sample (replicate per replicate)
## no merge extension
## qvalue 2.3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.sh # 38043666 ok
sbatch scripts/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.sh # 38043672 ok

#### KO
sbatch scripts/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.sh # 38043680 ok
sbatch scripts/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.sh # 38043689 ok





### H3K4me3

## qvalue 2.3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKO_H3K4me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K4me3_005-FergusonUniqueNorm99.sh # 38043792 ok
sbatch scripts/LengthNormSignal_WTKO_H3K4me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K4me3_008-FergusonUniqueNorm99.sh # 38043797 ok

#### KO
sbatch scripts/LengthNormSignal_WTKO_H3K4me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K4me3_005-FergusonUniqueNorm99.sh # 38043913 ok
sbatch scripts/LengthNormSignal_WTKO_H3K4me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K4me3_008-FergusonUniqueNorm99.sh # 38043933 ok


```

--> I set here `--binSize 100 --regionBodyLength 100`; seems it give 1 value per row/peak. Look good.



### H3K27me3, no extension, qval 2.3 - R DESEQ2


```R
library("tidyverse")
library("DESeq2")
#library("edgeR")
library("EnhancedVolcano")



set.seed(42)

# import bed reference to collect gene name
NPC_WTKO_H3K27me3_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 
SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())





# import BED position from matrix
BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())


# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")  

SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")  



# Tidy into a single tibble
SCORE_BED_WTKO_H3K27me3_pool_peaks = SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008) %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005) %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008)


######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO = SCORE_BED_WTKO_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  


# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -peakID), pull(countData_WTvsKO, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKO_raw <- SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKO_raw, -sample), pull(colData_WTvsKO_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5 # below 2000 look like noise on IGV
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.1 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.1 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, NPC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0.1,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(res_tibble$log2FoldChange > 0.1 & res_tibble$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_tibble$log2FoldChange < -0.1 & res_tibble$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.1 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.1 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, pvalue) %>%
  filter(pvalue < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, pvalue) %>%
  filter(pvalue < 0.05, log2FoldChange < -0.1)

```

--> Seems that there are fewer Diff. bound regions as compare to when using THOR; here; 40 lost 35 gain 

Let's compare with Venn diagram the gene that gain/lost H3K27me3 with Ferguson and THOR method; files:
- *Ferguson, Gain*: `output/edgeR/upregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3.txt` # 35 genes
- *Ferguson, Lost*: `output/edgeR/downregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3.txt` # 40 genes
- *THOR, Gain*: `output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol.txt` # 840 genes
- *THOR, Lost*: `output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol.txt` # 224 genes


Let's try edgeR for diff. binding
--> After tested, it seems edgeR identify more sites, but mostly noise, very small peaks. So I would rather use DESEQ2 with relax parameters. 





### H3K27me3, no extension, qval 2.3 - R DESEQ2 RELAX PARAMETERS

Things to adapt:
- Extract results with relaxed threshold (FDR < 0.1 instead of 0.05):
res <- results(dds, alpha=0.1)  # Adjusted p-value threshold qval 0.1 instead of 0.05
- Reduce the log2 fold-change threshold (e.g., detect smaller effects):
res <- results(dds, alpha=0.1, lfcThreshold=0.5)  # Default is 1.0
- Use different shrinkage methods (apeglm = conservative, ashr = relaxed)
res_shrunk <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="ashr")

```bash
conda activate deseq2
```


```R
library("tidyverse")
library("DESeq2")
library("ashr")
library("EnhancedVolcano")



set.seed(42)

# import bed reference to collect gene name
NPC_WTKO_H3K27me3_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 
SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())





# import BED position from matrix
BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())


# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")  

SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")  



# Tidy into a single tibble
SCORE_BED_WTKO_H3K27me3_pool_peaks = SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008) %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005) %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008)


######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO = SCORE_BED_WTKO_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  


# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -peakID), pull(countData_WTvsKO, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKO_raw <- SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKO_raw, -sample), pull(colData_WTvsKO_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5 # below 2000 look like noise on IGV
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

## HERE PICK METHOD !!!
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="ashr") # HERE CHANGE normal ashr apeglm = `lfcShrinkASHR`
res <- results(dds, alpha=0.1, pAdjustMethod = "fdr") # = `resultsFDR`



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.1 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.1 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-resultsFDR.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-resultsFDR.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, NPC, H3K27me3',
  pCutoff = 5e-2,         # 5e-2
  FCcutoff = 0.1,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(res_tibble$log2FoldChange > 0.1 & res_tibble$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_tibble$log2FoldChange < -0.1 & res_tibble$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.1 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.1 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-resultsFDR.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-resultsFDR.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.1, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.1, log2FoldChange < -0.1)


res_tibble %>%
  mutate(FDR = p.adjust(pvalue, method = "fdr")) %>% 
  dplyr::select(peakID, geneSymbol, log2FoldChange, FDR) %>%
  filter(FDR < 0.05, log2FoldChange < -0.1)
res_tibble %>%
  mutate(FDR = p.adjust(pvalue, method = "fdr")) %>% 
  dplyr::select(peakID, geneSymbol, log2FoldChange, FDR) %>%
  filter(FDR < 0.05, log2FoldChange > 0.1)



```

--> ashr, apeglm, or using results() instead of lfcshrinkage() works great, gave more diff bound regions. Surprinsigly, still more region that lose H3K27me3 than one gaining H3K27me3... Opposite result as when using THOR. 
  --> Let's check deepTool profile on these genes and check both THOR and Ferguson bigwigs





### H3K27me3, no extension, qval 2.3 - R EDGER

Let's use the [edgeR guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

```bash
conda activate monocle3
```


```R
# packages
library('tidyverse')
library('edgeR')
library("EnhancedVolcano")



set.seed(42)

# import bed reference to collect gene name
NPC_WTKO_H3K27me3_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_NPC_WTKO_H3K27me3_pool_peaks_merge_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 
SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())





# import BED position from matrix
BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())


# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")  

SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 = SCORE_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 %>%
  left_join(BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008 ) %>%
  left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")  



# Tidy into a single tibble
SCORE_BED_WTKO_H3K27me3_pool_peaks = SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_005 %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_WT_H3K27me3_008) %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_005) %>%
  bind_rows(SCORE_BED_WTKO_H3K27me3_pool_peaks__NPC_KO_H3K27me3_008)


######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO = SCORE_BED_WTKO_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKO_H3K27me3_pool_peaks_WTvsKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -peakID), pull(countData_WTvsKO, peakID)) 

### Test with log2 normalize counts
log2_counts <- log2(counts_all_matrix + 1)
log2_counts_rounded <- round(log2_counts)
head(log2_counts_rounded)
############


# CREATE DGEList object ####
## Define sample groups (adjust based on your design)
group <- factor(c("WT", "WT", "KO", "KO"))
## Create DGEList
y <- DGEList(counts=counts_all_matrix, group=group) # CHANGE TO LOG COUNT OR NOT use: log2_counts_rounded OR counts_all_matrix
## Check DGEList
y # 13005 genes
## fIlter out lowly express genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE] # 12984 genes

# y <- normLibSizes(y) # TO RUN OR NOT !!!!!!!!!!!!!!!!!!!!!

y <- estimateDisp(y)

########### Diff binding with the Exact Test #########################
et <- exactTest(y)
topTags(et)
######################################################################

########### Diff binding with Glmfit limma --> PREFERRED OPTION !!!!!! #########################
design <- model.matrix(~group)
## quasi-likelihood (QL) F-test
fit <- glmQLFit(y, design)
qlf.2vs1 <- glmQLFTest(fit, coef=2) # KO vs WT
topTags(qlf.2vs1)

## likelihood ratio test ######--> PREFERRED OPTION !!!!!! 
fit <- glmFit(y, design)
lrt.2vs1 <- glmLRT(fit, coef=2)
topTags(lrt.2vs1)
######################################################################




##############  TESTING    ###########################################################





########### Diff binding with the limma followed by the eBayes t-statistics #########################
## Define sample groups
group <- factor(c("WT", "WT", "KO", "KO"))
## Create design matrix for limma
design <- model.matrix(~ 0 + group)  # Avoids intercept
colnames(design) <- levels(group)

## Compute mean signal per peak (log2 transformed counts)
mean_signal <- rowMeans(log2_counts)

## Set a threshold to remove low-signal peaks (adjust as needed)
threshold <- 2  # You can try 0.5, 1, or 2 depending on the noise level
filtered_log2_counts <- log2_counts[mean_signal > threshold, ]


## Print the design matrix
design

## Fit linear model
fit <- lmFit(filtered_log2_counts, design)      # CHANGE TO LOG COUNT OR NOT use: filtered_log2_counts OR log2_counts_rounded OR counts_all_matrix
## Define the contrast (KO vs WT)
contrast_matrix <- makeContrasts(KO - WT, levels=design)
## Apply the contrast
fit2 <- contrasts.fit(fit, contrast_matrix)
## Apply empirical Bayes moderation
fit2 <- eBayes(fit2)
## Extract results
results <- topTable(fit2, adjust="fdr", number=Inf) # BH BY holm
## Check the first few results
head(results)

results = as_tibble(results, rownames = "peakID")
results %>% filter(adj.P.Val < 0.05)

results_all <- results %>% left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot)



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  results_all$logFC < -0 & results_all$adj.P.Val < 5e-2, 'Sky Blue',
    ifelse(results_all$logFC > 0 & results_all$adj.P.Val < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'


# Export result
write.table(results_all, file="output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-limmaEbayesLog2Count.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_q05fc0-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-limmaEbayesLog2Count.pdf", width=3, height=4)    
EnhancedVolcano(results_all,
  lab = results_all$geneSymbol,
  x = 'logFC',
  y = 'adj.P.Val',
  title = 'KO vs WT, NPC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(results_all$logFC > 0 & results_all$adj.P.Val < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(results_all$logFC < -0 & results_all$adj.P.Val < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- results_all[!is.na(results_all$logFC) & !is.na(results_all$adj.P.Val) & results_all$logFC > 0.1 & results_all$adj.P.Val < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- results_all[!is.na(results_all$logFC) & !is.na(results_all$adj.P.Val) & results_all$logFC < -0.1 & results_all$adj.P.Val < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-limmaEbayesLog2Count.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-limmaEbayesLog2Count.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


results_all %>% dplyr::select(peakID, geneSymbol, logFC, adj.P.Val) %>%
  filter(adj.P.Val < 0.05, logFC > 0)
results_all %>% dplyr::select(peakID, geneSymbol, logFC, adj.P.Val) %>%
  filter(adj.P.Val < 0.05, logFC < -0)



topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, PValue) %>%
  filter(PValue < 0.05, logFC > 0)
topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, PValue) %>%
  filter(PValue < 0.05, logFC < -0)




######################################################################






########################################
## likelihood ratio test without TMM normalization, with raw count (NOT log2 norm) ########
####################################
topTags_all = as_tibble(topTags(lrt.2vs1, n = Inf)$table, rownames = "peakID") %>%
  mutate(logFC = -logFC)


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, 'Sky Blue',
    ifelse(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'


topTags_all <- topTags_all %>% left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(topTags_all, file="output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_q05fc0-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.pdf", width=3, height=4)    
EnhancedVolcano(topTags_all,
  lab = topTags_all$geneSymbol,
  x = 'logFC',
  y = 'FDR',
  title = 'KO vs WT, NPC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC > 0.1 & topTags_all$FDR < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC < -0.1 & topTags_all$FDR < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCounts.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC > 0)
topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC < -0)



topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, PValue) %>%
  filter(PValue < 0.05, logFC > 0)
topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, PValue) %>%
  filter(PValue < 0.05, logFC < -0)






########################################
## Exact test, without TMM normalization, with raw counts (NOT log2 norm) ########
####################################
topTags_all = as_tibble(topTags(et, n = Inf)$table, rownames = "peakID") %>%
  mutate(logFC = -logFC)


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, 'Sky Blue',
    ifelse(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'


topTags_all <- topTags_all %>% left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(topTags_all, file="output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCounts.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_q05fc0-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCounts.pdf", width=3, height=4)    
EnhancedVolcano(topTags_all,
  lab = topTags_all$geneSymbol,
  x = 'logFC',
  y = 'FDR',
  title = 'KO vs WT, NPC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC > 0.1 & topTags_all$FDR < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC < -0.1 & topTags_all$FDR < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCounts.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCounts.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC > 0)
topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC < -0)






########################################
## Exact test, WITH TMM normalization, with raw counts (NOT log2 norm) ########
####################################
topTags_all = as_tibble(topTags(et, n = Inf)$table, rownames = "peakID") %>%
  mutate(logFC = -logFC)


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, 'Sky Blue',
    ifelse(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'


topTags_all <- topTags_all %>% left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(topTags_all, file="output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCountsTMM.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_q05fc0-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCountsTMM.pdf", width=3, height=4)    
EnhancedVolcano(topTags_all,
  lab = topTags_all$geneSymbol,
  x = 'logFC',
  y = 'FDR',
  title = 'KO vs WT, NPC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC > 0.1 & topTags_all$FDR < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC < -0.1 & topTags_all$FDR < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCountsTMM.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-ExactTestRawCountsTMM.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC > 0)
topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC < -0)




########################################
## likelihood ratio test, WITH TMM normalization, with raw counts (NOT log2 norm) ########
####################################
topTags_all = as_tibble(topTags(lrt.2vs1, n = Inf)$table, rownames = "peakID") %>%
  mutate(logFC = -logFC)


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, 'Sky Blue',
    ifelse(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'


topTags_all <- topTags_all %>% left_join(NPC_WTKO_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(topTags_all, file="output/edgeR/EDGER-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCountsRawCountsTMM.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_q05fc0-WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCountsRawCountsTMM.pdf", width=3, height=4)    
EnhancedVolcano(topTags_all,
  lab = topTags_all$geneSymbol,
  x = 'logFC',
  y = 'FDR',
  title = 'KO vs WT, NPC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(topTags_all$logFC > 0 & topTags_all$FDR < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(topTags_all$logFC < -0 & topTags_all$FDR < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC > 0.1 & topTags_all$FDR < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- topTags_all[!is.na(topTags_all$logFC) & !is.na(topTags_all$FDR) & topTags_all$logFC < -0.1 & topTags_all$FDR < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCountsRawCountsTMM.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc0_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-GlmfitLikelihoodRatioRawCountsRawCountsTMM.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC > 0)
topTags_all %>% dplyr::select(peakID, geneSymbol, logFC, FDR) %>%
  filter(FDR < 0.05, logFC < -0)


```


Use log2 norm count or raw count?
--> Not sure, test both; log2 norm count gave very few DEG, so **prefer raw count**

Which diff binding test: **Exact test, or Glmfit limma**
--> Ferguson used limma, so **focus on limma**

With limma use QL test or likelihood test
--> **Likelihood test** give more DEG, always, whatever test done


*Testing*:
- *likelihood ratio test, without TMM normalization, with raw counts (NOT log2 norm)*: 121 down / 40 up; looks good on IGV --> Look good
- *Exact test, without TMM normalization, with raw counts (NOT log2 norm)*: 112 down / 31 up; looks good on IGV --> Look good, like likelihood ratio test
- *likelihood ratio test, without TMM normalization, with log2 norm counts*: 1 up / 0 down --> Look bad
- *Exact test, without TMM normalization, with log2 norm counts*: no DEGs

- *Exact test, WITH TMM normalization, with raw counts (NOT log2 norm)*: 116 down / 33 up; looks good on IGV --> Look good
- *likelihood ratio test, WITH TMM normalization, with raw counts (NOT log2 norm)*: 123 down / 42 up; looks good on IGV --> Look good
- *Diff binding with the limma followed by the eBayes t-statistics* = Paper method: 50 down / 26 up 


--> Prefer **likelihood ratio test, without TMM normalization, with raw counts (NOT log2 norm)** = `GlmfitLikelihoodRatioRawCounts`; as I am not sure what the TMM is doing and that do not change much results. Let's export genes and check overlap with THOR:
  --> *For Lost*: Very few overlap. Interestingly, method-specific genes seems true!! Seems THOR misses complete loss in mutant (like signal to complete absence); and EDGER misses slight differences but still visible. For the ones detected in THOR, but not in edgeR, the pvalue is low but padj much higher...
  --> *For gain*: same



--> It may be worst changing pval trehsold or method for pvalue adjustment? 
  --> It did not help that much, detect more but still miss many from THOR...



# CSAW sliding window - Ferguson/Local Maxima


Let's use the [csaw/sliding window method](https://bioconductor.org/books/release/csawBook/counting-reads-into-windows.html) to check for diff binding. CSAW uses BAM as input, so we will need instead to use conmputeMatrix from 
deeptools to generate matrix.

Let's check which window size to set:
- CSAW_H3K27me3: csaw guide recommend bin of 2000bp every 500 bp (overlapping window)
- THOR_H3K27me3: THOR used per default bin of 100bp every 50 bp (overlapping window)

From [this paper](https://www.nature.com/articles/s41467-018-03538-9): Window widths were set to reflect broad or more narrow distribution of the investigated histone modifications: width 1000 bp and a spacing interval of 100 bp were used for H3K27me3 and H3K9me3, whereas width 150 bp and spacing 50 bp were used for H3K4me3


Lets try both:
- 1000bp every 100bp = `bin1000space100`
- 150bp every 50bp = `bin150space50`


## Generate bed file that cover the whole genome = window bed

```bash
conda activate BedToBigwig
# To calculate signal
bedtools makewindows -g ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -w 1000 -s 100 > ../../Master/meta/GRCh38_bin1000space100.bed
bedtools makewindows -g ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -w 150 -s 50 > ../../Master/meta/GRCh38_bin150space50.bed
bedtools makewindows -g ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -w 150 -s 150 > ../../Master/meta/GRCh38_bin150space150.bed
# To filter out low abundance window
bedtools makewindows -g ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -w 10000 -s 10000 > ../../Master/meta/GRCh38_bin10000space10000.bed
```
--> only chr 1-21 X,Y,M are included in `../../Master/meta/GRCh38_chrom_sizes_MAIN.tab`

--> Looks good!


## Calculate signal in the whole genome

- 1000bp every 100bp = `bin1000space100`
- 150bp every 50bp = `bin150space50`
- 150bp every 150bp = `bin150space150`


This correspond to that part of the code in *csaw*: `win.data <- windowCounts(h3k27me3data$Path, param=param, width=2000, spacing=500, ext=200)`


```bash
conda activate deeptools

#H3K27me3
## sample per sample (replicate per replicate)
#### WT
sbatch scripts/LengthNormSignal-bin1000space100-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.sh # 38173952 different rows from the other; 38316733
sbatch scripts/LengthNormSignal-bin1000space100-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.sh # 38173981 fail misname sample; 38239043 ok

sbatch scripts/LengthNormSignal-bin150space50-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.sh # 38174112 ok
sbatch scripts/LengthNormSignal-bin150space50-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.sh # 38174128 ok

sbatch scripts/LengthNormSignal-bin150space150-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.sh # 38544052 xxx
sbatch scripts/LengthNormSignal-bin150space150-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.sh # 38544120 xxx

#### KO
sbatch scripts/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.sh # 38174041 ok
sbatch scripts/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.sh # 38174071 ok

sbatch scripts/LengthNormSignal-bin150space50-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.sh # 38174151 ok
sbatch scripts/LengthNormSignal-bin150space50-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.sh # 38174227 ok

sbatch scripts/LengthNormSignal-bin150space150-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.sh # 38544188 xxx
sbatch scripts/LengthNormSignal-bin150space150-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.sh # 38544203 xxx
```
--> All good


Let's calculate *broad* signal; **to filter out low abundance window**. This corresp

This correspond to that part of the code in *csaw*: `bins <- windowCounts(h3k27me3data$Path, bin=TRUE, width=10000, param=param)`. And it is basically:
- 10000bp every 10000bp = `bin10000space10000` = NO OVERLAP OF THE WINDOW HERE
  --> The ../../Master/meta/GRCh38_bin1000space100.bed as been updated to `../../Master/meta/GRCh38_bin10000space10000.bed` below


```bash
conda activate deeptools

#H3K27me3
## sample per sample (replicate per replicate)
#### WT
sbatch scripts/LengthNormSignal-bin10000space10000-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.sh # 38250761 ok
sbatch scripts/LengthNormSignal-bin10000space10000-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.sh # 38250811 ok

#### KO
sbatch scripts/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.sh # 38251323 ok
sbatch scripts/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.sh # 38251502 ok
```
--> All good


## Run csaw - bin1000space100

Let's follow the [csaw guide](https://bioconductor.org/books/release/csawBook/counting-reads-into-windows.html)

```bash
conda activate DiffBind
```

```R
# packages
library("tidyverse")
library("csaw")
library("edgeR")
library("statmod")

set.seed(42)


##################################################################
# import samples (bin1000space100) ######################
####################################################################


# import SCORE 
SCORE_NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, scoer per row, coordinate and row
SCORE_BED_NPC_WT_H3K27me3_005 = SCORE_NPC_WT_H3K27me3_005 %>%
  left_join(BED_NPC_WT_H3K27me3_005 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_NPC_WT_H3K27me3_008 = SCORE_NPC_WT_H3K27me3_008 %>%
  left_join(BED_NPC_WT_H3K27me3_008 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")

SCORE_BED_NPC_KO_H3K27me3_005 = SCORE_NPC_KO_H3K27me3_005 %>%
  left_join(BED_NPC_KO_H3K27me3_005 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_NPC_KO_H3K27me3_008 = SCORE_NPC_KO_H3K27me3_008 %>%
  left_join(BED_NPC_KO_H3K27me3_008 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")



# Convert to RangedSummarizedExperiment Object
## Convert to GRange
samples <- list(
  "WT_005" = SCORE_BED_NPC_WT_H3K27me3_005,
  "WT_008" = SCORE_BED_NPC_WT_H3K27me3_008,
  "KO_005" = SCORE_BED_NPC_KO_H3K27me3_005,
  "KO_008" = SCORE_BED_NPC_KO_H3K27me3_008
)
## Ensure all samples have the same genomic windows (assuming identical structure)
gr <- GRanges(
  seqnames = samples$WT_005$chr, 
  ranges = IRanges(start = samples$WT_005$start, end = samples$WT_005$end)
)
## Create a matrix of counts where each column corresponds to a sample
counts_matrix <- do.call(cbind, lapply(samples, function(df) df$score))
## Define colData (metadata for samples)
col_data <- data.frame(
  sample_name = names(samples),
  genotype = c("WT", "WT", "KO", "KO"),
  replicate = c("R1", "R2", "R1", "R2"),
  totals = colSums(counts_matrix)  # Library size per sample
)
## Create SummarizedExperiment object
RSE_SCORE_BED_NPC_H3K27me3 <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  rowRanges = gr,
  colData = col_data
)
## Check object
RSE_SCORE_BED_NPC_H3K27me3


##################################################################
# import background sample (bin10000space10000) ######################
####################################################################

# import SCORE 
SCORE_NPC_WT_H3K27me3_005_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_WT_H3K27me3_008_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_005_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_008_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_NPC_WT_H3K27me3_005_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_WT_H3K27me3_008_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_005_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_008_background <- read.delim("output/edgeR/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, scoer per row, coordinate and row
SCORE_BED_NPC_WT_H3K27me3_005_background = SCORE_NPC_WT_H3K27me3_005_background %>%
  left_join(BED_NPC_WT_H3K27me3_005_background ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")

SCORE_BED_NPC_WT_H3K27me3_008_background = SCORE_NPC_WT_H3K27me3_008_background %>%
  left_join(BED_NPC_WT_H3K27me3_008_background ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")

SCORE_BED_NPC_KO_H3K27me3_005_background = SCORE_NPC_KO_H3K27me3_005_background %>%
  left_join(BED_NPC_KO_H3K27me3_005_background ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_NPC_KO_H3K27me3_008_background = SCORE_NPC_KO_H3K27me3_008_background %>%
  left_join(BED_NPC_KO_H3K27me3_008_background ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")



# Convert to RangedSummarizedExperiment Object
## Convert to GRange
samples <- list(
  "WT_005" = SCORE_BED_NPC_WT_H3K27me3_005_background,
  "WT_008" = SCORE_BED_NPC_WT_H3K27me3_008_background,
  "KO_005" = SCORE_BED_NPC_KO_H3K27me3_005_background,
  "KO_008" = SCORE_BED_NPC_KO_H3K27me3_008_background
)
## Ensure all samples have the same genomic windows (assuming identical structure)
gr <- GRanges(
  seqnames = samples$WT_005$chr, 
  ranges = IRanges(start = samples$WT_005$start, end = samples$WT_005$end)
)
## Create a matrix of counts where each column corresponds to a sample
counts_matrix <- do.call(cbind, lapply(samples, function(df) df$score))
## Define colData (metadata for samples)
col_data <- data.frame(
  sample_name = names(samples),
  genotype = c("WT", "WT", "KO", "KO"),
  replicate = c("R1", "R2", "R1", "R2"),
  totals = colSums(counts_matrix)  # Library size per sample
)
## Create SummarizedExperiment object
RSE_SCORE_BED_NPC_H3K27me3_background <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  rowRanges = gr,
  colData = col_data
)
## Check object
RSE_SCORE_BED_NPC_H3K27me3_background


#--> Lets forget about the backgruond part for now

# Filter out low quality windows by count size
abundances <- aveLogCPM(asDGEList(RSE_SCORE_BED_NPC_H3K27me3))
summary(abundances)


## Convert to data frame for ggplot2
abundance_df <- data.frame(logCPM = abundances)
pdf("output/csaw/hist_aveLogCPM-bin1000space100.pdf", width=6, height=5)
ggplot(abundance_df, aes(x = logCPM)) +
  geom_histogram(bins = 50, fill = "gray", color = "black") +
  ggtitle("Distribution of Window Abundance (aveLogCPM)") +
  xlab("Average Log CPM") +
  ylab("Frequency")+
  geom_vline(xintercept = -11, color = "red", linetype = "dashed", size = 1) +  
  theme_bw()
dev.off()

keep <- abundances > -11 # Here treshold of 10 apply. CAN BE CHANGED!!!
summary(keep)

## Apply filtering
filtered.data <- RSE_SCORE_BED_NPC_H3K27me3[keep,]

# Normalization
# --> LETS TRY WITHOUT NORMALIZATION, as theorcially already done...

# Statistical modelling
y <- asDGEList(filtered.data)
str(y)

genotype <- RSE_SCORE_BED_NPC_H3K27me3$genotype
## Convert to factor
genotype <- factor(genotype, levels = c("WT", "KO"))  # Ensures WT is first
## Create design matrix
design <- model.matrix(~0 + genotype)
colnames(design) <- levels(genotype)
design

# estimate the negative binomial (NB) and quasi-likelihood (QL) dispersions for each window
y <- estimateDisp(y, design) # Model biological variability across replicates./Reduce technical noise and improve variance structure.
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE) # Shrink extreme dispersions, ensuring more accurate differential windows.
summary(fit$df.prior)

pdf("output/csaw/plotMDS_glmQLFit-bin1000space100.pdf", width=6, height=5)
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("black", "red")[as.integer(genotype)])
dev.off()


# Save - LOAD ####################################
save.image(file = "output/csaw/bin1000space100_v1.RData") # keep <- abundances > -10
load("output/csaw/bin1000space100_v1.RData")
##################################################

# test for DB between conditions in each window using the QL F-test
contrast <- makeContrasts(KO-WT, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

# Grouping windows into regions
merged <- mergeResultsList(list(filtered.data), 
    tab.list=list(res$table),
    equiweight=TRUE, tol=100, merge.args=list(max.width=30000))
merged$regions


tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
table(tabcom$direction[is.sig])
tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

out.ranges <- merged$regions
mcols(out.ranges) <- data.frame(tabcom, best.logFC=tabbest$rep.logFC)

# Add annotation
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")

anno <- detailRanges(out.ranges, orgdb=org.Hs.eg.db,
    txdb=TxDb.Hsapiens.UCSC.hg38.knownGene)
mcols(out.ranges) <- cbind(mcols(out.ranges), anno)

# SAVE OUTPUT
out_tibble =  as_tibble(as.data.frame(out.ranges)) %>%
  dplyr::rename(seqnames = seqnames, start = start, end = end, strand = strand)
write_tsv(out_tibble , "output/csaw/results_bin1000space100.tsv")

# Check some genes
out_tibble %>%
  dplyr::select(seqnames, start, end, rep.logFC, best.logFC, FDR, direction, overlap) %>%
  filter(FDR < 0.05, best.logFC > 0)


```

- NOTE: `SummarizedExperiment()` need same library size, otheriwse error. However my totals are different between background and my counts because in background the counts are counted without overlapping bins. Versus in my count samples, bins are overlapping so some reads are counted mutliple times, so we end up with a much higher library size! So I **manually normalize library size**: `win.data$totals`=`RSE_SCORE_BED_NPC_H3K27me3` to match `bins$totals`=`RSE_SCORE_BED_NPC_H3K27me3_background`
  - Then error with `SummarizedExperiment()`:  `Error in if (prop.seen > 1 { : missing value where TRUE/FALSE needed`. Not sure why, there is no NA, and I try removing all the windowns with 0 but still same error. So let's instead **manually filtered low-abundance windows**
--> Lets forget about the background part for now and let's use another moethdo to filter out low quality windows. Instead, lets simply **Filter out low quality windows by count size**

- NOTE: # `tol` at `mergeWindows()` or `mergeResultsList()` sets the minimum distance to merge binding sites: large values (500-1000 bp) reduce redundancy, while small values (<200 bp) resolve individual sites. max.width limits cluster size, splitting larger ones into equal subclusters.

--> My bigiwgs are already normalized, so let's try to NOT apply any normalization. 
  --> Few sites obtain, and most down... Let's try to apply TMM normalization. But for **TMM normalization I need background and sample with same total.. So for that I will count into non-overlapping bins!**





## Run csaw - bin1000space100_genePromoters - signal in gene and promoters only - Without TMM normalization

Let's follow the [csaw guide](https://bioconductor.org/books/release/csawBook/counting-reads-into-windows.html)

```bash
conda activate DiffBind
```

```R
# packages
library("tidyverse")
library("csaw")
library("edgeR")
library("statmod")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")

set.seed(42)

# Load bin1000space100
load("output/csaw/bin1000space100_v1.RData")

# Lets only keep windows in gene and promoters

broads <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
broads <- resize(broads, width(broads)+5000, fix="end")
head(broads)
suppressWarnings(keep <- overlapsAny(rowRanges(RSE_SCORE_BED_NPC_H3K27me3), broads))
sum(keep)


## Apply filtering
filtered.data <- RSE_SCORE_BED_NPC_H3K27me3[keep,]

# Normalization
# --> LETS TRY WITHOUT NORMALIZATION, as theorcially already done...

# Statistical modelling
y <- asDGEList(filtered.data)
str(y)

genotype <- RSE_SCORE_BED_NPC_H3K27me3$genotype
## Convert to factor
genotype <- factor(genotype, levels = c("WT", "KO"))  # Ensures WT is first
## Create design matrix
design <- model.matrix(~0 + genotype)
colnames(design) <- levels(genotype)
design

# estimate the negative binomial (NB) and quasi-likelihood (QL) dispersions for each window
y <- estimateDisp(y, design) # Model biological variability across replicates./Reduce technical noise and improve variance structure.
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE) # Shrink extreme dispersions, ensuring more accurate differential windows.
summary(fit$df.prior)

pdf("output/csaw/plotMDS_glmQLFit-bin1000space100_genePromoters.pdf", width=6, height=5)
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("black", "red")[as.integer(genotype)])
dev.off()


# Save - LOAD ####################################
save.image(file = "output/csaw/bin1000space100_genePromoters_v1.RData") # keep windows in gene promoter 5kb upstream
load("output/csaw/bin1000space100_genePromoters_v1.RData")
##################################################

# test for DB between conditions in each window using the QL F-test
contrast <- makeContrasts(KO-WT, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

# Grouping windows into regions
merged <- mergeResultsList(list(filtered.data), 
    tab.list=list(res$table),
    equiweight=TRUE, tol=100, merge.args=list(max.width=30000))
merged$regions


tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
table(tabcom$direction[is.sig])
tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

out.ranges <- merged$regions
mcols(out.ranges) <- data.frame(tabcom, best.logFC=tabbest$rep.logFC)

# Add annotation
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")

anno <- detailRanges(out.ranges, orgdb=org.Hs.eg.db,
    txdb=TxDb.Hsapiens.UCSC.hg38.knownGene)
mcols(out.ranges) <- cbind(mcols(out.ranges), anno)

# SAVE OUTPUT
out_tibble =  as_tibble(as.data.frame(out.ranges)) %>%
  dplyr::rename(seqnames = seqnames, start = start, end = end, strand = strand)
write_tsv(out_tibble , "output/csaw/results_bin1000space100_genePromoters.tsv")

# Check some genes
out_tibble %>%
  dplyr::select(seqnames, start, end, rep.logFC, best.logFC, FDR, direction, overlap) %>%
  filter(FDR < 0.05, best.logFC > 0)



```







## Run csaw - bin150space150 - Proportion background + TMM norm

Let's do like in the H3K27me3 tutorial; use large bins to remove low quality windows, and uses them to perform TMM normalization.

BUT, use non-overlapping bins for both sample and background to have the same library size!



```bash
conda activate DiffBind
```

```R
# packages
library("tidyverse")
library("csaw")
library("edgeR")
library("statmod")

set.seed(42)


##################################################################
# import background (bin10000space10000) ######################
####################################################################

# Load background = bin10000space10000
load("output/csaw/bin1000space100_v1.RData")
RSE_SCORE_BED_NPC_H3K27me3_background # = bin10000space10000


##################################################################
# import samples (bin150space150) ######################
####################################################################


# import SCORE 
SCORE_NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space150-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, scoer per row, coordinate and row
SCORE_BED_NPC_WT_H3K27me3_005 = SCORE_NPC_WT_H3K27me3_005 %>%
  left_join(BED_NPC_WT_H3K27me3_005 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_NPC_WT_H3K27me3_008 = SCORE_NPC_WT_H3K27me3_008 %>%
  left_join(BED_NPC_WT_H3K27me3_008 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")

SCORE_BED_NPC_KO_H3K27me3_005 = SCORE_NPC_KO_H3K27me3_005 %>%
  left_join(BED_NPC_KO_H3K27me3_005 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_NPC_KO_H3K27me3_008 = SCORE_NPC_KO_H3K27me3_008 %>%
  left_join(BED_NPC_KO_H3K27me3_008 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")



# Convert to RangedSummarizedExperiment Object
## Convert to GRange
samples <- list(
  "WT_005" = SCORE_BED_NPC_WT_H3K27me3_005,
  "WT_008" = SCORE_BED_NPC_WT_H3K27me3_008,
  "KO_005" = SCORE_BED_NPC_KO_H3K27me3_005,
  "KO_008" = SCORE_BED_NPC_KO_H3K27me3_008
)
## Ensure all samples have the same genomic windows (assuming identical structure)
gr <- GRanges(
  seqnames = samples$WT_005$chr, 
  ranges = IRanges(start = samples$WT_005$start, end = samples$WT_005$end)
)
## Create a matrix of counts where each column corresponds to a sample
counts_matrix <- do.call(cbind, lapply(samples, function(df) df$score))
## Define colData (metadata for samples)
col_data <- data.frame(
  sample_name = names(samples),
  genotype = c("WT", "WT", "KO", "KO"),
  replicate = c("R1", "R2", "R1", "R2"),
  totals = colSums(counts_matrix)  # Library size per sample
)
## Create SummarizedExperiment object
RSE_SCORE_BED_NPC_H3K27me3 <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  rowRanges = gr,
  colData = col_data
)
## Check object
RSE_SCORE_BED_NPC_H3K27me3


# Save - LOAD ####################################
save.image(file = "output/csaw/bin150space150_PropBackground_TMM_v1.RData") # Just sample bin150space150 and backgruond are loaded
load("output/csaw/bin150space150_PropBackground_TMM_v1.RData")
##################################################



# Filtering of low-abundance windows
## Scale/normalize our background to have same totals as in sample
scale_factor <- sum(RSE_SCORE_BED_NPC_H3K27me3$totals) / sum(RSE_SCORE_BED_NPC_H3K27me3_background$totals)
assay(RSE_SCORE_BED_NPC_H3K27me3_background) <- assay(RSE_SCORE_BED_NPC_H3K27me3_background) * scale_factor
RSE_SCORE_BED_NPC_H3K27me3_background$totals <- RSE_SCORE_BED_NPC_H3K27me3$totals # force them to be identical, rounding issue lead to slight differences

# Add chr size to our RSE objects
chrom_sizes <- read.delim("../../Master/meta/GRCh38_chrom_sizes_MAIN.tab", header=FALSE, sep="\t", col.names = c("chromosome", "size"))

chrom_size_vector <- setNames(chrom_sizes$size, chrom_sizes$chromosome)

## Assign chromosome sizes to rowRanges of RSE objects
seqlevelsStyle(rowRanges(RSE_SCORE_BED_NPC_H3K27me3)) <- "UCSC"  # Ensure UCSC style if needed
seqlengths(rowRanges(RSE_SCORE_BED_NPC_H3K27me3)) <- chrom_size_vector

seqlevelsStyle(rowRanges(RSE_SCORE_BED_NPC_H3K27me3_background)) <- "UCSC"  # Ensure UCSC style if needed
seqlengths(rowRanges(RSE_SCORE_BED_NPC_H3K27me3_background)) <- chrom_size_vector

# Remove windows without any counts

log2_cpm <- cpm(assay(RSE_SCORE_BED_NPC_H3K27me3), log=TRUE, prior.count=1) # Compute log2-CPM for `RSE_SCORE_BED_NPC_H3K27me3`
log2_cpm_means <- rowMeans(log2_cpm)

## Generate the histogram plot
pdf("output/csaw/hist_log2_cpm_means-RSE_SCORE_BED_NPC_H3K27me3-bin150space150_PropBackground_TMM.pdf", width=6, height=5)
hist(log2_cpm_means, main="Log2 CPM Distribution", breaks=50, col="grey",
    xlab="Log2 CPM", ylab="Frequency", border="black")
abline(v=-8.5, col="red", lwd=2, lty=2)  # Threshold line in red
dev.off()

keep <- log2_cpm_means > -8.5
summary(keep)

## Apply the filter to the RSE object
RSE_SCORE_BED_NPC_H3K27me3_filtered <- RSE_SCORE_BED_NPC_H3K27me3[keep, ]




## FILTER
#filter.stat <- filterWindowsGlobal(RSE_SCORE_BED_NPC_H3K27me3, RSE_SCORE_BED_NPC_H3K27me3_background) # All windows kept
filter.stat <- filterWindowsGlobal(RSE_SCORE_BED_NPC_H3K27me3_filtered, RSE_SCORE_BED_NPC_H3K27me3_background) # window with very few reads removed

min.fc <- 1

#pdf("output/csaw/hist_filter.stat-bin150space150_PropBackground_TMM.pdf", width=6, height=5)
pdf("output/csaw/hist_filter.stat_log2CPMmeans8.5-bin150space150_PropBackground_TMM.pdf", width=6, height=5)
hist(filter.stat$filter, main="", breaks=50,
    xlab="Background abundance (log2-CPM)")
abline(v=log2(min.fc), col="red")
dev.off()

keep2 <- filter.stat$filter > log2(min.fc)
summary(keep2)
## Aply filtering
filtered.data <- RSE_SCORE_BED_NPC_H3K27me3_filtered[keep2,]


# Normalization for composition biases
RSE_SCORE_BED_NPC_H3K27me3 <- normFactors(RSE_SCORE_BED_NPC_H3K27me3_background, se.out=RSE_SCORE_BED_NPC_H3K27me3)
(normfacs <- RSE_SCORE_BED_NPC_H3K27me3$norm.factors)


# Statistical modelling
y <- asDGEList(filtered.data)
str(y)


genotype <- RSE_SCORE_BED_NPC_H3K27me3$genotype
## Convert to factor
genotype <- factor(genotype, levels = c("WT", "KO"))  # Ensures WT is first
## Create design matrix
design <- model.matrix(~0 + genotype)
colnames(design) <- levels(genotype)
design
 


y <- estimateDisp(y, design)
summary(y$trended.dispersion)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)


pdf("output/csaw/plotMDS_glmQLFit-bin150space150_PropBackground_TMM.pdf", width=6, height=5)
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("black", "red")[as.integer(genotype)])
dev.off()


contrast <- makeContrasts(KO-WT, levels=design)

res <- glmQLFTest(fit, contrast=contrast)


# Consolidating results from multiple window sizes
#--> not done



# Grouping windows into regions
merged <- mergeResultsList(list(filtered.data), 
    tab.list=list(res$table),
    equiweight=TRUE, tol=100, merge.args=list(max.width=30000))
merged$regions


tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
table(tabcom$direction[is.sig])
tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

out.ranges <- merged$regions
mcols(out.ranges) <- data.frame(tabcom, best.logFC=tabbest$rep.logFC)

# Add annotation
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")

anno <- detailRanges(out.ranges, orgdb=org.Hs.eg.db,
    txdb=TxDb.Hsapiens.UCSC.hg38.knownGene)
mcols(out.ranges) <- cbind(mcols(out.ranges), anno)

# SAVE OUTPUT
out_tibble =  as_tibble(as.data.frame(out.ranges)) %>%
  dplyr::rename(seqnames = seqnames, start = start, end = end, strand = strand)
write_tsv(out_tibble , "output/csaw/results_bin150space150_PropBackground_TMM.tsv")

# Check some genes
out_tibble %>%
  dplyr::select(seqnames, start, end, rep.logFC, best.logFC, FDR, direction, overlap) %>%
  filter(FDR < 0.05, best.logFC > 0)












```


- NOTE: With this new method, totals is almost identical, but not exactly, lets manually scale the totals from background to be same as totals from sample. 
- NOTE: As I generated the RSE object manually, it misses chr lenght, I thus have to add them manually at `## Assign chromosome sizes to rowRanges of RSE objects`





## Run csaw - bin150space50

Let's follow the [csaw guide](https://bioconductor.org/books/release/csawBook/counting-reads-into-windows.html)

```bash
conda activate DiffBind
```

```R
# packages
library("tidyverse")
library("csaw")
library("edgeR")
library("statmod")

set.seed(42)


##################################################################
# import samples (bin150space50) ######################
####################################################################


# import SCORE 
SCORE_NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_NPC_WT_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_WT_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_WT_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_WT_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_005 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_NPC_KO_H3K27me3_008 <- read.delim("output/edgeR/LengthNormSignal-bin150space50-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, scoer per row, coordinate and row
SCORE_BED_NPC_WT_H3K27me3_005 = SCORE_NPC_WT_H3K27me3_005 %>%
  left_join(BED_NPC_WT_H3K27me3_005 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_NPC_WT_H3K27me3_008 = SCORE_NPC_WT_H3K27me3_008 %>%
  left_join(BED_NPC_WT_H3K27me3_008 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")

SCORE_BED_NPC_KO_H3K27me3_005 = SCORE_NPC_KO_H3K27me3_005 %>%
  left_join(BED_NPC_KO_H3K27me3_005 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_NPC_KO_H3K27me3_008 = SCORE_NPC_KO_H3K27me3_008 %>%
  left_join(BED_NPC_KO_H3K27me3_008 ) %>%
  dplyr::select(chr, start, end, score) %>%
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")



# Convert to RangedSummarizedExperiment Object
## Convert to GRange
samples <- list(
  "WT_005" = SCORE_BED_NPC_WT_H3K27me3_005,
  "WT_008" = SCORE_BED_NPC_WT_H3K27me3_008,
  "KO_005" = SCORE_BED_NPC_KO_H3K27me3_005,
  "KO_008" = SCORE_BED_NPC_KO_H3K27me3_008
)
## Ensure all samples have the same genomic windows (assuming identical structure)
gr <- GRanges(
  seqnames = samples$WT_005$chr, 
  ranges = IRanges(start = samples$WT_005$start, end = samples$WT_005$end)
)
## Create a matrix of counts where each column corresponds to a sample
counts_matrix <- do.call(cbind, lapply(samples, function(df) df$score))
## Define colData (metadata for samples)
col_data <- data.frame(
  sample_name = names(samples),
  genotype = c("WT", "WT", "KO", "KO"),
  replicate = c("R1", "R2", "R1", "R2"),
  totals = colSums(counts_matrix)  # Library size per sample
)
## Create SummarizedExperiment object
RSE_SCORE_BED_NPC_H3K27me3 <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  rowRanges = gr,
  colData = col_data
)
## Check object
RSE_SCORE_BED_NPC_H3K27me3



# Filter out low quality windows by count size
abundances <- aveLogCPM(asDGEList(RSE_SCORE_BED_NPC_H3K27me3))
summary(abundances)


## Convert to data frame for ggplot2
abundance_df <- data.frame(logCPM = abundances)
pdf("output/csaw/hist_aveLogCPM-bin150space50.pdf", width=6, height=5)
ggplot(abundance_df, aes(x = logCPM)) +
  geom_histogram(bins = 50, fill = "gray", color = "black") +
  ggtitle("Distribution of Window Abundance (aveLogCPM)") +
  xlab("Average Log CPM") +
  ylab("Frequency")+
  geom_vline(xintercept = -9, color = "red", linetype = "dashed", size = 1) +  
  theme_bw()
dev.off()



keep <- abundances > -9 # Here treshold of 10 apply. CAN BE CHANGED!!!
summary(keep)

## Apply filtering
filtered.data <- RSE_SCORE_BED_NPC_H3K27me3[keep,]

# Normalization
# --> LETS TRY WITHOUT NORMALIZATION, as theorcially already done...

# Statistical modelling
y <- asDGEList(filtered.data)
str(y)

genotype <- RSE_SCORE_BED_NPC_H3K27me3$genotype
## Convert to factor
genotype <- factor(genotype, levels = c("WT", "KO"))  # Ensures WT is first
## Create design matrix
design <- model.matrix(~0 + genotype)
colnames(design) <- levels(genotype)
design

pdf("output/csaw/plotMDS-bin150space50.pdf", width=6, height=5)
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("black", "red")[as.integer(genotype)])
dev.off()

# estimate the negative binomial (NB) and quasi-likelihood (QL) dispersions for each window
y <- estimateDisp(y, design) # Model biological variability across replicates./Reduce technical noise and improve variance structure.
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE) # Shrink extreme dispersions, ensuring more accurate differential windows.
summary(fit$df.prior)

pdf("output/csaw/plotMDS_glmQLFit-bin150space50.pdf", width=6, height=5)
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("black", "red")[as.integer(genotype)])
dev.off()


# Save - LOAD ####################################
save.image(file = "output/csaw/bin150space50_v1.RData")
load("output/csaw/bin150space50_v1.RData")
##################################################


# test for DB between conditions in each window using the QL F-test
contrast <- makeContrasts(KO-WT, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

# Grouping windows into regions
merged <- mergeResultsList(list(filtered.data), 
    tab.list=list(res$table),
    equiweight=TRUE, tol=100, merge.args=list(max.width=30000))
merged$regions


tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
table(tabcom$direction[is.sig])
tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

# Save rds
out.ranges <- merged$regions
mcols(out.ranges) <- data.frame(tabcom, best.logFC=tabbest$rep.logFC)
saveRDS(file="output/csaw/results_bin150space50.rds", out.ranges)

# Add annotation
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")

anno <- detailRanges(out.ranges, orgdb=org.Hs.eg.db,
    txdb=TxDb.Hsapiens.UCSC.hg38.knownGene)
mcols(out.ranges) <- cbind(mcols(out.ranges), anno)

# SAVE OUTPUT
out_tibble =  as_tibble(as.data.frame(out.ranges)) %>%
  dplyr::rename(seqnames = seqnames, start = start, end = end, strand = strand)
write_tsv(out_tibble , "output/csaw/results_bin150space50.tsv")

# Check some genes
out_tibble %>%
  dplyr::select(seqnames, start, end, rep.logFC, best.logFC, FDR, direction, overlap) %>%
  filter(FDR < 0.05, best.logFC > 0)


```

--> Very few Diff. bound regions here

--> Maybe the filtering of low quality window is to harsh? Or not enough?


## Run csaw - bin150space50 + bin1000space100



```bash
conda activate DiffBind
```

```R
# packages
library("tidyverse")
library("csaw")
library("edgeR")
library("statmod")

set.seed(42)

# Load `y` and `fit` from bin1000space100 and bin150space50
load("output/csaw/bin1000space100_v1.RData")
fit_bin1000space100 <- glmQLFit(y, design, robust=TRUE) # Shrink extreme dispersions, ensuring more accurate differential windows.
summary(fit$df.prior)
y_bin1000space100 = y

load("output/csaw/bin150space50_v1.RData")
fit_bin150space50 = fit
y_bin150space5 = y


```






# diffreps - differential binding analysis

[diffreps](https://github.com/shenlab-sinai/diffreps) uses sliding window method on bed file. 

## Install diffreps


```bash
cd ../../Master/software

# download repo
git clone https://github.com/shenlab-sinai/diffreps.git
cd diffreps

# use useless conda env to install diffreps
conda activate ChIPseqSpikeInFree
conda install bioconda::perl-cpan-shell # to allow Perl CPAN installation

# install the perl dependencies --> Will be installed in ChIPseqSpikeInFree conda env
perl -MCPAN -e 'install Statistics::TTest'
perl -MCPAN -e 'install Math::CDF'
perl -MCPAN -e 'install Parallel::ForkManager'

# install
perl Makefile.PL (Optional, PREFIX=your_perl_directory)
make
make test
make install
```
--> All good




## Run diffreps


Input files needed:
- output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.sorted.bedGraph
- output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.sorted.bedGraph
- output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.sorted.bedGraph
- output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.sorted.bedGraph


Troubleshoot [discussion](https://groups.google.com/forum/?fromgroups\#!forum/diffreps-discuss)
--> diffreps required BED6 file (= A BED file where each feature is described by chrom, start, end, name, score, and strand.)



```bash
conda activate ChIPseqSpikeInFree

# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.sorted.bedGraph > output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.sorted.bedGraph > output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.sorted.bedGraph > output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.sorted.bedGraph > output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed



########### 1000bp every 100bp #################################

# 1000bp every 100bp (Default histone) - Negative binomial
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb-diff.nb.txt --window 1000 --step 100
# --> 111 diff sites

# 1000bp every 100bp (Default histone) - Negative binomial Pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb_pval05-diff.nb.txt --window 1000 --step 100 --pval 0.05
# --> 111 diff sites (same as just pval cutoff is changed...)

# 1000bp every 100bp (Default histone) - Negative binomial Pval 0.001
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb_pval001-diff.nb.txt --window 1000 --step 100 --pval 0.001
# --> xxx diff sites (same as just pval cutoff is changed...)


# 1000bp every 100bp (Default histone) - Negative binomial with nsd 10
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb_nsd10-diff.nb.txt --window 1000 --step 100 --nsd 10
# --> 89 diff sites



# 1000bp every 100bp (Default histone) - t test
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_tt-diff.nb.txt --window 1000 --step 100 --meth tt
# --> 17 diff sites

# 1000bp every 100bp (Default histone) - G test
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_gt-diff.nb.txt --window 1000 --step 100 --meth gt
# --> 202 diff sites



# 1000bp every 100bp (Default histone) - G test with pval 0.001
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval001-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.001
# --> xxx diff sites 

# 1000bp every 100bp (Default histone) - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05
# --> xxx diff sites 



########### 500bp every 100bp #################################

# 500bp every 100bp (Default histone) - Negative binomial
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_nb-diff.nb.txt --window 500 --step 100
# --> 77 diff sites

# 500bp every 100bp (Default histone) - Negative binomial with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_nb_pval05-diff.nb.txt --window 500 --step 100 --pval 0.05
# --> xxx diff sites

# 500bp every 100bp (Default histone) - Negative binomial with pval 0.001
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_nb_pval001-diff.nb.txt --window 500 --step 100 --pval 0.001
# --> xxx diff sites

# 500bp every 100bp (Default histone) - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05
# --> xxx diff sites

# 500bp every 100bp (Default histone) - G test with pval 0.001
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_gt_pval001-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.001
# --> xxx diff sites




########### 250bp every 50bp #################################

# 250bp every 50bp (Default histone) - Negative binomial
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_nb-diff.nb.txt --window 250 --step 50
# --> 10 diff sites

# 250bp every 50bp (Default histone) - Negative binomial with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_nb_pval05-diff.nb.txt --window 250 --step 50 --pval 0.05

# 250bp every 50bp (Default histone) - Negative binomial with pval 0.001
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_nb_pval001-diff.nb.txt --window 250 --step 50 --pval 0.001




# 250bp every 50bp (Default histone) - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05

# 250bp every 50bp (Default histone) - G test with pval 0.001
diffReps.pl -tr output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bed -co output/bigwig_Ferguson/NPC_WT_H3K27me3_005_unique_norm99.bed output/bigwig_Ferguson/NPC_WT_H3K27me3_008_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_gt_pval001-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.001



```


--> *1000bp every 100bp looks optimal* (more diff site and recommended by the tool)

--> Changing *Pval treshold to 0.001* looks good; then investigate in R more in detail; but it collect more diff sites. That look true for the most. = This correspond to the nb of windows kept!!

--> *NB test looks ok* (and recommended to use). G-test looks more 'relax', could use it too!



## Explore diffreps results in R

These versions looks good:

- 1000bp every 100bp (Default histone) - Negative binomial Pval 0.05: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb_pval05-diff.nb.txt`
- 1000bp every 100bp (Default histone) - Negative binomial Pval 0.001: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb_pval001-diff.nb.txt`
- 1000bp every 100bp (Default histone) - G test with pval 0.05: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval05-diff.nb.txt`
- 1000bp every 100bp (Default histone) - G test with pval 0.001: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval001-diff.nb.txt`

- 500bp every 100bp (Default histone) - Negative binomial Pval 0.05: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_nb_pval05-diff.nb.txt`
- 500bp every 100bp (Default histone) - Negative binomial Pval 0.001: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_nb_pval001-diff.nb.txt`
- 500bp every 100bp (Default histone) - G test with pval 0.05: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_gt_pval05-diff.nb.txt`
- 500bp every 100bp (Default histone) - G test with pval 0.001: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_gt_pval001-diff.nb.txt`



- 250bp every 50bp (Default histone) - Negative binomial Pval 0.05: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_nb_pval05-diff.nb.txt`
- 250bp every 50bp (Default histone) - Negative binomial Pval 0.001: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_nb_pval001-diff.nb.txt`
- 250bp every 50bp (Default histone) - G test with pval 0.05: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_gt_pval05-diff.nb.txt`
- 250bp every 50bp (Default histone) - G test with pval 0.001: `output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_gt_pval001-diff.nb.txt`



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")


# import files
bin1000space100_nb_pval05 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_nb_pval001 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_nb_pval001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval001 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 

bin500space100_nb_pval05 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_nb_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_nb_pval001 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_nb_pval001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval001 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin500space100_gt_pval001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 


bin250space50_nb_pval05 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_nb_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_nb_pval001 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_nb_pval001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval001 <- read.delim("output/diffreps/NPC_H3K27me3_unique_norm99.bed-bin250space50_gt_pval001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 

# Replace Inf by min/max values

bin1000space100_nb_pval05$log2FC[bin1000space100_nb_pval05$log2FC == Inf] <- max(bin1000space100_nb_pval05$log2FC[is.finite(bin1000space100_nb_pval05$log2FC)], na.rm = TRUE)
bin1000space100_nb_pval05$log2FC[bin1000space100_nb_pval05$log2FC == -Inf] <- min(bin1000space100_nb_pval05$log2FC[is.finite(bin1000space100_nb_pval05$log2FC)], na.rm = TRUE)

bin1000space100_nb_pval001$log2FC[bin1000space100_nb_pval001$log2FC == Inf] <- max(bin1000space100_nb_pval001$log2FC[is.finite(bin1000space100_nb_pval001$log2FC)], na.rm = TRUE)
bin1000space100_nb_pval001$log2FC[bin1000space100_nb_pval001$log2FC == -Inf] <- min(bin1000space100_nb_pval001$log2FC[is.finite(bin1000space100_nb_pval001$log2FC)], na.rm = TRUE)

bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == Inf] <- max(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == -Inf] <- min(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin1000space100_gt_pval001$log2FC[bin1000space100_gt_pval001$log2FC == Inf] <- max(bin1000space100_gt_pval001$log2FC[is.finite(bin1000space100_gt_pval001$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval001$log2FC[bin1000space100_gt_pval001$log2FC == -Inf] <- min(bin1000space100_gt_pval001$log2FC[is.finite(bin1000space100_gt_pval001$log2FC)], na.rm = TRUE)


bin500space100_nb_pval05$log2FC[bin500space100_nb_pval05$log2FC == Inf] <- max(bin500space100_nb_pval05$log2FC[is.finite(bin500space100_nb_pval05$log2FC)], na.rm = TRUE)
bin500space100_nb_pval05$log2FC[bin500space100_nb_pval05$log2FC == -Inf] <- min(bin500space100_nb_pval05$log2FC[is.finite(bin500space100_nb_pval05$log2FC)], na.rm = TRUE)

bin500space100_nb_pval001$log2FC[bin500space100_nb_pval001$log2FC == Inf] <- max(bin500space100_nb_pval001$log2FC[is.finite(bin500space100_nb_pval001$log2FC)], na.rm = TRUE)
bin500space100_nb_pval001$log2FC[bin500space100_nb_pval001$log2FC == -Inf] <- min(bin500space100_nb_pval001$log2FC[is.finite(bin500space100_nb_pval001$log2FC)], na.rm = TRUE)

bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == Inf] <- max(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)
bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == -Inf] <- min(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)

bin500space100_gt_pval001$log2FC[bin500space100_gt_pval001$log2FC == Inf] <- max(bin500space100_gt_pval001$log2FC[is.finite(bin500space100_gt_pval001$log2FC)], na.rm = TRUE)
bin500space100_gt_pval001$log2FC[bin500space100_gt_pval001$log2FC == -Inf] <- min(bin500space100_gt_pval001$log2FC[is.finite(bin500space100_gt_pval001$log2FC)], na.rm = TRUE)


bin250space50_nb_pval05$log2FC[bin250space50_nb_pval05$log2FC == Inf] <- max(bin250space50_nb_pval05$log2FC[is.finite(bin250space50_nb_pval05$log2FC)], na.rm = TRUE)
bin250space50_nb_pval05$log2FC[bin250space50_nb_pval05$log2FC == -Inf] <- min(bin250space50_nb_pval05$log2FC[is.finite(bin250space50_nb_pval05$log2FC)], na.rm = TRUE)

bin250space50_nb_pval001$log2FC[bin250space50_nb_pval001$log2FC == Inf] <- max(bin250space50_nb_pval001$log2FC[is.finite(bin250space50_nb_pval001$log2FC)], na.rm = TRUE)
bin250space50_nb_pval001$log2FC[bin250space50_nb_pval001$log2FC == -Inf] <- min(bin250space50_nb_pval001$log2FC[is.finite(bin250space50_nb_pval001$log2FC)], na.rm = TRUE)

bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == Inf] <- max(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)
bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == -Inf] <- min(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)

bin250space50_gt_pval001$log2FC[bin250space50_gt_pval001$log2FC == Inf] <- max(bin250space50_gt_pval001$log2FC[is.finite(bin250space50_gt_pval001$log2FC)], na.rm = TRUE)
bin250space50_gt_pval001$log2FC[bin250space50_gt_pval001$log2FC == -Inf] <- min(bin250space50_gt_pval001$log2FC[is.finite(bin250space50_gt_pval001$log2FC)], na.rm = TRUE)


# List of dataset names
file_names <- c("bin1000space100_nb_pval05", "bin1000space100_nb_pval001", 
                "bin1000space100_gt_pval05", "bin1000space100_gt_pval001",
                "bin500space100_nb_pval05", "bin500space100_nb_pval001",
                "bin500space100_gt_pval05", "bin500space100_gt_pval001",
                "bin250space50_nb_pval05", "bin250space50_nb_pval001",
                "bin250space50_gt_pval05", "bin250space50_gt_pval001")

## Function to read and format each file
read_and_process <- function(file) {
  df <- get(file)  # Load dataset from environment
  df$dataset <- file  # Add dataset identifier
  return(df)
}

## Combine all datasets into one
combined_data <- bind_rows(lapply(file_names, read_and_process)) 

combined_data_counts <- combined_data %>% 
  filter(padj<0.05) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")

## plot

pdf("output/diffreps/hist-log2FC_distribution-padj05.pdf", width=7, height=4)
combined_data %>% 
  filter(padj<0.05) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 3) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 7, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()


# Focus bin1000space100_gt_pval05 and bin500space100_gt_pval05
combined_data_counts_filt = combined_data_counts %>%  filter(dataset %in% c("bin1000space100_gt_pval05" , "bin500space100_gt_pval05"))


pdf("output/diffreps/hist-log2FC_distribution-padj05_gt_pval05.pdf", width=7, height=4)
combined_data %>% 
  filter(padj<0.05) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!! 
  filter(dataset %in% c("bin1000space100_gt_pval05" , "bin500space100_gt_pval05")) %>% 
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 7, face = "bold")) +
  geom_text(data = combined_data_counts_filt , 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts_filt$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()



# Investigate some unique genes


combined_data %>% 
  filter(padj<0.05, log2FC > 0) %>%   
  filter(dataset %in% c("bin1000space100_gt_pval05")) %>%
  dplyr::select(-Length, -dataset) 


combined_data %>% 
  filter(padj<0.05, log2FC > 0) %>%   
  filter(dataset %in% c("bin500space100_gt_pval05")) %>%
  dplyr::select(-Length, -dataset) 





```


--> changing bin size gave different results:
  - bin 1000 seems to give more diff bound regions, and overall lost
  - bin 500 seems to give less diff bound regions, and overall gain

--> The pvalur for window gave more site at 0.05 than 0.001 --> Prefer pval 0.05

--> gt gave more diff. regions than nb 

--> To maximize diff bound regions: pval05 gt = `bin1000space100_gt_pval05` and `bin500space100_gt_pval05`
  --> These two seems good on IGV too. At padj 0.05. Lets try to find a way to integrate both into one file; to avoid overlap/duplicate
    --> Check what csaw did?

--> 







