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
sbatch scripts/matrix_TSS_10kb_H3K27me3_median_THORLIBspikein_allGenes.sh # 17135510 xxx


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

# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 

### GeneSymbol list of signif up/down genes in each genotypes
output/ChIPseeker/annotation_THOR_H3K4me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram194.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram355.txt

output/ChIPseeker/annotation_THOR_H3K4me3_q20_neg_promoterAnd5_geneSymbol_Venndiagram462.txt
output/ChIPseeker/annotation_THOR_H3K4me3_q20_pos_promoterAnd5_geneSymbol_Venndiagram636.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt

output/ChIPseeker/annotation_THOR_H3K27me3_q30_neg_promoterAnd5_geneSymbol_Venndiagram220.txt
output/ChIPseeker/annotation_THOR_H3K27me3_q30_pos_promoterAnd5_geneSymbol_Venndiagram836.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt", header=FALSE, stringsAsFactors=FALSE)
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

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K4me3",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
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
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_THOR_H3K27me3_q40.txt", sep="\t", row.names=FALSE, quote=FALSE)





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

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt", header=FALSE, stringsAsFactors=FALSE)
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
write.table(gos, "output/GO/enrichR_GO_Molecular_Function_2023_THOR_H3K27me3_q40.txt", sep="\t", row.names=FALSE, quote=FALSE)





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

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt", header=FALSE, stringsAsFactors=FALSE)
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
write.table(gos, "output/GO/enrichR_GO_Cellular_Component_2023_THOR_H3K27me3_q40.txt", sep="\t", row.names=FALSE, quote=FALSE)







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

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_neg_promoterAnd5_geneSymbol_Venndiagram145.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_THOR_H3K27me3_q40_pos_promoterAnd5_geneSymbol_Venndiagram554.txt", header=FALSE, stringsAsFactors=FALSE)
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
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_THOR_H3K27me3_q40.txt", sep="\t", row.names=FALSE, quote=FALSE)


```



# Does the peak that gain H3K27me3 in KO is because of EZH2 binding in NPC?


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
sbatch scripts/matrix_TSS_10kb_H3K27me3_EZH2_SUZ12_median_THOR_gainLost_gene.sh # 17197231 xxx


```

--> No changes of EZH2 and SUZ12 in peak that gain H3K27me3 in KO... But only 1 replicate... Also looking at all genes strong decrease of EZH2 and SUZ12 binding in KO, here, same level... Another replicate may show changes!!!

--> Also all peaks that gain H3K27me3 in KO, seems already bound by EZH2 and SUZ12 in the WT! So maybe there is an increase EZH2 activity upon EZH1 KO, rather than an increase binding

--> Showing gene or peaks lead to same result (EZH2, SUZ12, barely follow H3K27me3 WTvsKO pattern)

--> `*THORLIBspikein*` norm is weird as show always gain of EZH2/SUZ12, even in regions that lose H3K27me3...




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


