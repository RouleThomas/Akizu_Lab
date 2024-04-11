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

```


--> experimental batch effect.. Cluster more per experiment than per genotype for raw, DiffBindTMM and THOR...
----> Maybe batch effect because of overall aspecific signal like in intergenic region? Let's do pearson corr plot in gene body region only, and then in gene body and promoters (add 2kb upstream TSS)
-----> **Batch effet removed with THOR TMM, still present with LIBspikein**


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
gene_names_down <- read.csv("meta/ENCFF159KBI_geneSymbol_500random_9.bed", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("meta/ENCFF159KBI_geneSymbol_500random_10.bed", header=FALSE, stringsAsFactors=FALSE)
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

# GENES
## Using classic Recirpocl DiffBind TMM method
sbatch scripts/matrix_TSS_10kb_H3K27me3_EZH2_SUZ12_median_THOR_gainLost_gene.sh # 17197231 ok


```

--> No changes of EZH2 and SUZ12 in peak that gain H3K27me3 in KO... But only 1 replicate... Also looking at all genes strong decrease of EZH2 and SUZ12 binding in KO, here, same level... Another replicate may show changes!!!

--> Also all peaks that gain H3K27me3 in KO, seems already bound by EZH2 and SUZ12 in the WT! So maybe there is an increase EZH2 activity upon EZH1 KO, rather than an increase binding

--> Showing gene or peaks lead to same result (EZH2, SUZ12, barely follow H3K27me3 WTvsKO pattern)

--> `*THORLIBspikein*` norm is weird as show always gain of EZH2/SUZ12, even in regions that lose H3K27me3...


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


