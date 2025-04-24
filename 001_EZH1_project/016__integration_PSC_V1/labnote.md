# Project goal

Data integration PSC WT, KO and KOEF1aEZH1:
    - EZH1
    - EZH2
    - SUZ12
    - H3K27me3

Detail of sample used:
- *WT_EZH1*: 1 Bio Rep (`006__CutRun`)
- *WT_EZH2*: 3 Bio Rep (`006__CutRun`, `010__CutRun`, `014__CutRun`)
- *WT_SUZ12*: 3 Bio Rep (`006__CutRun`, `013__CutRun` --> 2 Rep, use the best one!, `014__CutRun`)
- *WT_H3K27me3*: 3 Bio Rep (`006__CutRun`, `010__CutRun`, `013__CutRun` --> 2 Rep, use the best one!)
- *KO_EZH1*: 3 Bio Rep (`006__CutRun`, `013__CutRun` --> 2 Rep, use the best one!, `014__CutRun`)
- *KO_EZH2*: 3 Bio Rep (`013__CutRun` --> 2 Rep, use the best one!, `014__CutRun` *2 Rep)
- *KO_SUZ12*: 3 Bio Rep (`013__CutRun` --> 2 Rep, use the best one!, `014__CutRun` *2 Rep --> 3 Rep, use the best two ones!)
- *KO_H3K27me3*: 3 Bio Rep (`006__CutRun`, `013__CutRun` --> 2 Rep, use the best one!, `014__CutRun`)
- *KOEF1aEZH1_EZH1*: 3 Bio Rep (`005__CutRun`, `006__CutRun`, `013__CutRun` --> 2 Rep, use the best one!)
- *KOEF1aEZH1_EZH2*: 3 Bio Rep (`006__CutRun`, `013__CutRun` --> 2 Rep, use the best one!, `014__CutRun`)
- *KOEF1aEZH1_SUZ12*: 3 Bio Rep (`005__CutRun`, `006__CutRun`,`013__CutRun` --> 2 Rep, use the best one!)
- *KOEF1aEZH1_H3K27me3*: 3 Bio Rep (`005__CutRun`, `006__CutRun`,`013__CutRun` --> 2 Rep, use the best one!)

--> *NOTE: for `013__CutRun` R1 is better and has been used; for `014__CutRun` KO_SUZ12, I used R1 and R2: they are more similar, R3 show a lot more signal, outlier as compare to R1 R2.*

# Pipeline integration

- copy the `*.unique.dupmark.sorted.bam` and `*.bai` files
- re -do masc2 peak calling 
- Calculate the MG1655 SF
- DiffBind_TMM SF calculation
- THOR to scale bigwig


# Copy and rename files


--> Copy `*.unique.dupmark.sorted.bam` and `*.bai` files and manually rename them. File per file, not to make mistake, with their associated IGG! Use following nomenclature: [SAMPLENAME]_[005R1005R2008R3] for replicates.

--> I keep replicate number from the experiment; if no number after `R`, means experiment has only 1 rep.


# THOR diff binding
## Run THOR



```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge


# Default THOR TMM normalization (no E coli spike in norm)
sbatch scripts/THOR_PSC_WTvsKO_EZH1_TMM.sh # 29755426 fail; 29757697 fail
sbatch scripts/THOR_PSC_WTvsKO_EZH2_TMM.sh # 29755604 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_TMM.sh # 29755723 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_TMM.sh # 29755816 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_TMM.sh # 29755881 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_TMM.sh # 29755963 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_TMM.sh # 29756081 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_TMM.sh # 29756122 ok

# Default THOR genes normalization (no E coli spike in norm) - HK genes from THOR tutorial
sbatch scripts/THOR_PSC_WTvsKO_EZH1_housekeep.sh # 29758500 ok
sbatch scripts/THOR_PSC_WTvsKO_EZH2_housekeep.sh # 29758509 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_housekeep.sh # 29758553 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_housekeep.sh # 29758588 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeep.sh # 29758609 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeep.sh # 29758622 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeep.sh # 29758674 ok 
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeep.sh # 29758780 ok  


# Default THOR genes normalization (no E coli spike in norm) - HOX genes used
sbatch scripts/THOR_PSC_WTvsKO_EZH1_housekeepHOX.sh # 29760504 ok
sbatch scripts/THOR_PSC_WTvsKO_EZH2_housekeepHOX.sh # 29760650 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_housekeepHOX.sh # 29760682 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX.sh # 29760701 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOX.sh # 29760823 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOX.sh # 29760852 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOX.sh # 29761039 ok 
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX.sh # 29761092 ok  

# Default THOR genes normalization (no E coli spike in norm) - HOX genes used without input
sbatch scripts/THOR_PSC_WTvsKO_EZH1_housekeepHOXnoInput.sh # 29829147 ok
sbatch scripts/THOR_PSC_WTvsKO_EZH2_housekeepHOXoInput.sh # 29829152 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_housekeepHOXoInput.sh # 29829163 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXoInput.sh # 29829179 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOXoInput.sh # 29829181 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOXoInput.sh # 29829193 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOXoInput.sh # 29829206 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOXoInput.sh # 29829211 ok 

# Default THOR genes normalization (no E coli spike in norm) - HOX + HK genes used
sbatch scripts/THOR_PSC_WTvsKO_EZH1_housekeepHOXHK.sh # 29823147 ok
sbatch scripts/THOR_PSC_WTvsKO_EZH2_housekeepHOXHK.sh # 29823148 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_housekeepHOXHK.sh # 29823150 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXHK.sh # 29823152 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOXHK.sh # 29823155 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOXHK.sh # 29823157 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOXHK.sh # 29823160 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOXHK.sh # 29823164 ok

# Default THOR genes normalization (no E coli spike in norm) - HOX + HK genes used without input
sbatch scripts/THOR_PSC_WTvsKO_EZH1_housekeepHOXHKnoInput.sh # 29829033 ok
sbatch scripts/THOR_PSC_WTvsKO_EZH2_housekeepHOXHKnoInput.sh # 29829038 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_housekeepHOXHKnoInput.sh # 29829040 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXHKnoInput.sh # 29829049 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOXHKnoInput.sh # 29829051 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOXHKnoInput.sh # 29829060 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOXHKnoInput.sh # 29829066 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOXHKnoInput.sh # 29829072 ok

# THOR scaling factor manually calculated following EpiCypher (E coli spike in norm) 
sbatch scripts/THOR_PSC_WTvsKO_EZH1_EpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKO_EZH2_EpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_EpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_EpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_EpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_EpiCypher.sh # 29869761 ok; NO NEED to do the other as this one is so bad!
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_EpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_EpiCypher.sh #  xxx

# THOR scaling factor DiffBind_TMM_EpiCypher (E coli spike in norm) --> SF to use in THOR are the **reciprocal of MG1655_DiffBind_TMM**
sbatch scripts/THOR_PSC_WTvsKO_EZH1_DiffBindTMMEpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKO_EZH2_DiffBindTMMEpiCypher.sh # 32559674 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_DiffBindTMMEpiCypher.sh # 32560104 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMMEpiCypher.sh # 29871946 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_DiffBindTMMEpiCypher.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_DiffBindTMMEpiCypher.sh # 32559686 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_DiffBindTMMEpiCypher.sh # 32560204 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher.sh # 29872019 ok



# THOR scaling factor DiffBind_TMM (not Epicypher, not E coli spike in corrected) --> SF to use in THOR are the **reciprocal of DiffBind_TMM**
sbatch scripts/THOR_PSC_WTvsKO_EZH1_DiffBindTMM.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKO_EZH2_DiffBindTMM.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_DiffBindTMM.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMM.sh # 33908796 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1_DiffBindTMM.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_DiffBindTMM.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_DiffBindTMM.sh #  xxx
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMM.sh #  xxx


# THOR scaling factor from Ferguson unique norm99 - without IGG (SF directly collected from python code and copied to xlsx)
sbatch scripts/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput.sh # 35227285 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput.sh # 35227286 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput.sh # 35227288 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput.sh # 35227289 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput.sh # 35227291 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput.sh # 35227292 ok


# THOR scaling factor from Ferguson unique norm99 - with IGG (SF directly collected from python code and copied to xlsx)
sbatch scripts/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99.sh # 35227426 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99.sh # 35227427 ok
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99.sh # 35227428 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99.sh # 35227430 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99.sh # 35227431 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99.sh # 35227432 ok




```

**Conclusion for replicate similarity**:
- *Default THOR TMM normalization*: very bad, replicate very different (potential batch effect remaining, different signal noise ratio)
- *Housekeeping gene normalization with HK genes*: very bad, replicate very different (fail likely because HK genes lowly H3K27me3)
- *Housekeeping gene normalization with HOX genes*: very good, replicate cluster very well together
- *Housekeeping gene normalization with HK + HOX genes*: bad, replicate very different (fail likely because HK genes lowly H3K27me3 again, add noises)
- *Housekeeping gene normalization with HOX genes without input*: very good, replicate cluster very well together
- *Housekeeping gene normalization with HK + HOX genes without input*: bad, replicate very different
- *SpikeIn EpiCypher normalization*: very bad, replicate very different (potential over correction, SF are huge values)
- *SpikeIn EpiCypher DiffBind TMM normalization*: perform OK, less homogeneous than HOX normalization...
- *DiffBindTMM*: perform OK bad, less homogeneous than HOX normalization...
- *FergusonUniqueNorm99*: very good, replicate cluster very well together (*FergusonUniqueNorm99_noInput* perform similarly), let's prefer not using IGG, to be more in agreement with Ferguson method that do not use it.


--> *Error* `IndexError: cannot do a non-empty take from an empty axes.` on `scripts/THOR_PSC_WTvsKO_EZH1_TMM.sh`. Probably because WT vs KO, and KO EZH1 as no signal at all... Weird comparison! That is a control...

--> *Housekeeping genes* has been generated in `001*/002*` at `#### THOR with housekeeping genes normalization`. Collected from the [rgt-THOR tutorial](https://reg-gen.readthedocs.io/en/latest/thor/tool_usage.html)
    --> Let's also try **housekeeping gene normalization using the HOX genes**. Generate in `meta/`: Works great!!





```bash
# Generate gtf file for the HOX genes
nano meta/HOX_genes.txt
HOXA1
HOXA2
HOXA3
HOXA-AS3
HOXA4
HOXA5
HOXA6
HOXA7
HOXA9
HOXA10
HOXA10-AS
HOXA11
HOXA11-AS
HOXA13


### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' meta/HOX_genes.txt > meta/HOX_genes_as_gtf_geneSymbol.txt

## Filter the gtf
grep -Ff meta/HOX_genes_as_gtf_geneSymbol.txt ../015__RNAseq_PSC/meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_HOX_genes.gtf

## convert gtf to bed
awk 'BEGIN {OFS="\t"} $3 == "gene" {print $1, $4-1, $5, $9, ".", $7}' meta/ENCFF159KBI_HOX_genes.gtf > meta/ENCFF159KBI_HOX_genes.bed


```

- *NOTE: HOX genes manually selected, identified from IGV*




### Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R
# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# H3K27me3 WTvsKO housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/PSCWTvsKOH3K27me3housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 20) %>%
  group_by(X6) %>%
  summarise(n = n())



# H3K27me3 WTvsKOEF1aEZH1 housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX/PSCWTvsKOEF1aEZH1H3K27me3housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 20) %>%
  group_by(X6) %>%
  summarise(n = n())





# SUZ12 WTvsKO housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_SUZ12_housekeepHOX/PSCWTvsKOSUZ12housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
#--> NONE!



# SUZ12 WTvsKOEF1aEZH1 housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOX/PSCWTvsKOEF1aEZH1SUZ12housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
#--> NONE!



# EZH2 WTvsKO housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_EZH2_housekeepHOX/PSCWTvsKOEZH2housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
#--> NONE!


# EZH2 WTvsKOEF1aEZH1 housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOX/PSCWTvsKOEF1aEZH1EZH2housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
#--> NONE!




# EZH1 WTvsKO housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_EZH1_housekeepHOX/PSCWTvsKOEZH1housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
#--> NONE!





# EZH1 WTvsKOEF1aEZH1 housekeepHOX
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOX/PSCWTvsKOEF1aEZH1EZH1housekeepHOX-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOX/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOX/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOX/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 20) %>%
  group_by(X6) %>%
  summarise(n = n())



# H3K27me3 WTvsKO housekeepHOX without input
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXnoInput/PSCWTvsKOH3K27me3housekeepHOXnoInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXnoInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXnoInput/log2FC_qval10.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 10) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval10") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 15) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXnoInput/THOR_qval15.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 15) %>%
  group_by(X6) %>%
  summarise(n = n())




# SUZ12 WTvsKO housekeepHOX without input
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_SUZ12_housekeepHOXnoInput/PSCWTvsKOSUZ12housekeepHOXnoInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_SUZ12_housekeepHOXnoInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_SUZ12_housekeepHOXnoInput/log2FC_qval15.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 15) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval15") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 15) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_SUZ12_housekeepHOXnoInput/THOR_qval15.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 15) %>%
  group_by(X6) %>%
  summarise(n = n())






# EZH2 WTvsKO housekeepHOX without input
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_EZH2_housekeepHOXnoInput/PSCWTvsKOEZH2housekeepHOXnoInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_EZH2_housekeepHOXnoInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_EZH2_housekeepHOXnoInput/log2FC_qval15.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 15) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval15") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_EZH2_housekeepHOXnoInput/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 15) %>%
  group_by(X6) %>%
  summarise(n = n())



# H3K27me3 WTvsKOEF1aEZH1 housekeepHOX without input
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOXnoInput/PSCWTvsKOEF1aEZH1H3K27me3housekeepHOXnoInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOXnoInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOXnoInput/log2FC_qval15.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 15) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval15") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOXnoInput/THOR_qval20.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 15) %>%
  group_by(X6) %>%
  summarise(n = n())


# SUZ12 WTvsKOEF1aEZH1 housekeepHOX without input
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOXnoInput/PSCWTvsKOEF1aEZH1SUZ12housekeepHOXnoInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOXnoInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOXnoInput/log2FC_qval15.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 15) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval15") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_housekeepHOXnoInput/THOR_qval20.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 15) %>%
  group_by(X6) %>%
  summarise(n = n())




# EZH2 WTvsKOEF1aEZH1 housekeepHOX without input
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOXnoInput/PSCWTvsKOEF1aEZH1EZH2housekeepHOXnoInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOXnoInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOXnoInput/log2FC_qval15.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 15) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval15") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeepHOXnoInput/THOR_qval20.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 15) %>%
  group_by(X6) %>%
  summarise(n = n())





# H3K27me3 WTvsKO DiffBindTMMEpiCypher

diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMMEpiCypher/PSCWTvsKOH3K27me3DiffBindTMMEpiCypher-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMMEpiCypher/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMMEpiCypher/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMMEpiCypher/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 20) %>%
  group_by(X6) %>%
  summarise(n = n())




# H3K27me3 WTvsKOEF1aEZH1 DiffBindTMMEpiCypher

diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher/PSCWTvsKOEF1aEZH1H3K27me3DiffBindTMMEpiCypher-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())





# SUZ12 WTvsKO DiffBindTMMEpiCypher
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_SUZ12_DiffBindTMMEpiCypher/PSCWTvsKOSUZ12DiffBindTMMEpiCypher-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_SUZ12_DiffBindTMMEpiCypher/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_SUZ12_DiffBindTMMEpiCypher/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_SUZ12_DiffBindTMMEpiCypher/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())





# H3K27me3 WTvsKO Ferguson Unique Norm99 (no Input) THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/log2FC_qval30.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval30") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())




# H3K27me3 WTvsKOEF1aEZH1 Ferguson Unique Norm99 (no Input) THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/log2FC_qval30.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval30") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())




# SUZ12 WTvsKOEF1aEZH1 Ferguson Unique Norm99 (no Input) THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1SUZ12FergusonUniqueNorm99noInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())




# SUZ12 WTvsKO Ferguson Unique Norm99 (no Input) THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput/PSCWTvsKOSUZ12FergusonUniqueNorm99noInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput/THOR_qval20.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 20) %>%
  group_by(X6) %>%
  summarise(n = n())





# EZH2 WTvsKOEF1aEZH1 Ferguson Unique Norm99 (no Input) THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1EZH2FergusonUniqueNorm99noInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2", "count_KOEF1aEZH1_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2+count_KOEF1aEZH1_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput/log2FC_qval20.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval20") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 20) %>%
  group_by(X6) %>%
  summarise(n = n())




# EZH2 WTvsKO Ferguson Unique Norm99 (no Input) THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEZH2FergusonUniqueNorm99noInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/log2FC_qval10.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 10) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval10") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())





XXXY PURSUE FOR EZH2 and SUZ12 XXX

```



**Optimal qvalue:**
HousekeepHOXInput
- WTvsKO_H3K27me3: qval20
- WTvsKOEF1aEZH1_H3K27me3: qval20
HousekeepHOXnoInput
- WTvsKO_H3K27me3: qval10
- WTvsKO_SUZ12: qval15 (SUZ12 R1 s1-rep0, very noisy)
- WTvsKO_EZH2: qval15 
- WTvsKOEF1aEZH1_H3K27me3: qval15
- WTvsKOEF1aEZH1_SUZ12: qval15
DiffBindTMMEpiCypher
- WTvsKO_H3K27me3: qval20
- WTvsKOEF1aEZH1_H3K27me3: qval20
FergusonUniqueNorm99_noInput
- WTvsKO_H3K27me3: qval30
- WTvsKOEF1aEZH1_H3K27me3: qval30
- WTvsKO_SUZ12: qval20
- WTvsKOEF1aEZH1_SUZ12: qval20
- WTvsKOEF1aEZH1_EZH2: qval10



# DiffBind diff binding

Follow recommendation [here](https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf); let's try:
- without using spike in normalization
- using spike in normalization, same method as usual

I modified the meta file and put the replicate in the Tissue column to be used with BLOCK after (eg. Replicate column is ignore in previous meta file); seems need to use [design instead](https://support.bioconductor.org/p/9138520/)

--> Seems the order is: blacklist > count > normalize > contrast > analyze.

--> For `dba.count()` I could put summit= 250 for H3K27me3 broad peaks; and keep default for the other mark. Used in this [paper](https://www.cell.com/developmental-cell/pdf/S1534-5807(20)30551-7.pdf) and discuss in DiffBind Bioconductor


```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# ONE PER ONE
## H3K27me3
### Load dba
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_H3K27me3_RepTissue.txt", header = TRUE, sep = "\t"))
### Blacklist Greylist
sample_dba_blackgreylist = dba.blacklist(sample_dba, blacklist=TRUE, greylist=TRUE)
### Count

sample_blackgreylist_count = dba.count(sample_dba_blackgreylist) 
## This take time, here is checkpoint command to save/load:
save(sample_blackgreylist_count, file = "output/DiffBind/sample_blackgreylist_count_macs2raw_unique_H3K27me3_RepTissue.RData")
load("output/DiffBind/sample_blackgreylist_count_macs2raw_unique_H3K27me3_RepTissue.RData")


### Data normalization


#### TMM ################################################
sample_blackgreylist_count_TMM = dba.normalize(sample_blackgreylist_count, normalize = DBA_NORM_TMM)
##### Here is to retrieve the scaling factor value
sample_blackgreylist_count_TMM_SF = dba.normalize(sample_blackgreylist_count_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_blackgreylist_count_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_blackgreylist_count_TMM_SF_SF_H3K27me3.txt")
##### Contrast our Replicate (=Tissue)

### contrast
sample_blackgreylist_count_TMM_contrast = dba.contrast(sample_blackgreylist_count_TMM, categories=DBA_TREATMENT, design="~Tissue + Treatment")

sample_blackgreylist_count_TMM_contrast = dba.contrast(sample_blackgreylist_count_TMM, categories=DBA_TREATMENT, design="~Treatment")

# Diff bind. all method
sample_blackgreylist_count_TMM_contrast_analyze <- dba.analyze(sample_blackgreylist_count_TMM_contrast, method= DBA_ALL_METHODS)
# Diff bind.
sample_blackgreylist_count_TMM_contrast_analyze <- dba.analyze(sample_blackgreylist_count_TMM_contrast, method= DBA_DESEQ2)

pdf("output/DiffBind/plotPCA_DESEQ2_H3K27me3_contrast1.pdf", width=5, height=5)
dba.plotPCA(sample_blackgreylist_count_TMM_contrast_analyze, contrast=1)
dev.off()
pdf("output/DiffBind/plotPCA_DESEQ2_H3K27me3_contrast2.pdf", width=5, height=5)
dba.plotPCA(sample_blackgreylist_count_TMM_contrast_analyze, contrast=2)
dev.off()
## Retrieving the differentially bound sites
sample_blackgreylist_count_TMM_contrast_analyze_DESEQ2= dba.report(sample_blackgreylist_count_TMM_contrast_analyze, contrast=2)
sum(sample_blackgreylist_count_TMM_contrast_analyze_DESEQ2$Fold>0)

### Export
df= as.data.frame(sample_blackgreylist_count_TMM_contrast_analyze_DESEQ2)
df$start <- df$start - 1
write.table(df, file = "output/DiffBind/TMM_DESEQ2_H3K27me3_contrast2.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##### Here is to retrieve the scaling factor value

sample_blackgreylist_count_TMM_contrast_analyze_SF = dba.normalize(sample_blackgreylist_count_TMM_contrast_analyze, bRetrieve=TRUE)
console_output <- capture.output(print(sample_dba_blackgreylist_RLE_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_DBA_NORM_RLE_unique_SF_H3K27me3.txt")

```

--> The SF has been applied to bigwig but resutl in very different replicates.. Let's try [blocking factor to reduce batch effect between replicates](https://support.bioconductor.org/p/96441/)



- NOTE: Using blacklist before or after count change output. Prefer using before count
- NOTE: Very few DB when using design= ~Tissue + Treatment in dba.contrast(). Prefer using design= ~Treatment instead
  - Fore H3K27me3, DESEQ2 show more DB than EDGER
- NOTE: Big issue is that I can only collect SF after normalization; so it seems the Batch effect cannot be corrected with DiffBind, as the contrast is not taken into account in the SF collected! So cannot generate bigiwg...






# Ferguson method Diff binding


Ferguson method:
- Call peaks with SEACR or MACS2
  --> Done with macs2 (broad, default, no qvalue filtering)
- Calculate length-normalize signal for each locus (computeMatrix –scale-regions (gene or peak) )
  --> Not sure how to deal with the peak; so let's instead do per gene promoter (1kb up 250bp down); save in `output/edgeR`
- Diff. log-normalize Counts using edgeR and limma on consensus peak or genes (consensus peak?)


## Diff binding on consensus peak

- Option1: Identify peak in WT and in KO, separately using MACS2, then merge overlapping peak = consensus peak. Then calculate signal in these regions
- Option2: Identify peak in WT and KO sample together in MACS2. Like as input in MACS2 WT and KO samples. Automatically peak region should be adjusted
  -->Issue is that if peak only present in one condition; region might be missed…

Let's do option1 first, that look better overall; so let's simply bedtool intersect and keep the overlapping entire regions = consensus peak list (without pek extension, merging peak region if 100bp apart)

```bash
conda activate BedToBigwig

# concatenate and sort bed files
## Raw - non qvalue filtered ##############
cat output/macs2/broad/PSC_WT_H3K27me3_pool_peaks.broadPeak output/macs2/broad/PSC_KO_H3K27me3_pool_peaks.broadPeak output/macs2/broad/PSC_KOEF1aEZH1_H3K27me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak
cat output/macs2/broad/PSC_WT_EZH2_pool_peaks.broadPeak output/macs2/broad/PSC_KO_EZH2_pool_peaks.broadPeak output/macs2/broad/PSC_KOEF1aEZH1_EZH2_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak
cat output/macs2/broad/PSC_WT_SUZ12_pool_peaks.broadPeak output/macs2/broad/PSC_KO_SUZ12_pool_peaks.broadPeak output/macs2/broad/PSC_KOEF1aEZH1_SUZ12_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak
## qvalue 2.3 ##############
### WT KO KOEF
cat output/macs2/broad/broad_blacklist_qval2.30103/PSC_WT_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/PSC_KO_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/PSC_KOEF1aEZH1_H3K27me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval2.30103/PSC_WT_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/PSC_KO_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/PSC_KOEF1aEZH1_EZH2_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval2.30103/PSC_WT_SUZ12_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/PSC_KO_SUZ12_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/PSC_KOEF1aEZH1_SUZ12_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak
### WT KO
cat output/macs2/broad/broad_blacklist_qval2.30103/PSC_WT_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/PSC_KO_EZH2_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKO_EZH2_pool_peaks.sorted.broadPeak

## qvalue 3 ##############
cat output/macs2/broad/broad_blacklist_qval3/PSC_WT_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/PSC_KO_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/PSC_KOEF1aEZH1_H3K27me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval3/PSC_WT_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/PSC_KO_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/PSC_KOEF1aEZH1_EZH2_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval3/PSC_WT_SUZ12_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/PSC_KO_SUZ12_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/PSC_KOEF1aEZH1_SUZ12_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak




# merge = consensus peak identification
## Raw - non qvalue filtered ##############
### no merge extension
bedtools merge -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge.bed
### with 100bp peak merging
bedtools merge -d 100 -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp.bed
### with 500bp peak merging
bedtools merge -d 500 -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge500bp.bed

## qvalue 2.3 ##############
### no merge extension
#### WT KO KOEF
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge.bed
#### WT KO
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKO_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKO_EZH2_pool_peaks.sorted.merge.bed
### with 100bp peak merging
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp.bed
### with 500bp peak merging
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge500bp.bed

## qvalue 3 ##############
### no merge extension
bedtools merge -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge.bed
### with 100bp peak merging
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp.bed
### with 500bp peak merging
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge500bp.bed




```

--> All good; consensus peak files are: `output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge[SIZE].bed`


Now calculate **signal in consensus peak**:


```bash
conda activate deeptools

# condition per condition (not optimal in the end, as I need replicate for stats)
## no merge extension
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-FergusonUniqueNorm99smooth50bp.sh # 35839533 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-FergusonUniqueNorm99smooth50bp.sh # 35839677 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-FergusonUniqueNorm99smooth50bp.sh # 35839774 ok
## with 100bp peak merging
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp-FergusonUniqueNorm99smooth50bp.sh # 35839846 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp-FergusonUniqueNorm99smooth50bp.sh # 35839893 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp-FergusonUniqueNorm99smooth50bp.sh # 35839944 ok
## with 500bp peak merging
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp-FergusonUniqueNorm99smooth50bp.sh # 35839996 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge500bp-FergusonUniqueNorm99smooth50bp.sh # 35840048 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp-FergusonUniqueNorm99smooth50bp.sh # 35840149 ok



# sample per sample (replicate per replicate)
## no merge extension
## Raw - non qvalue filtered ##############
### H3K27me3
#### WT
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.sh # 35902923 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.sh # 35902926 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.sh # 35902928 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.sh # 35902932 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.sh # 35902933 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.sh # 35902936 ok
#### KOEF1aEZH1
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.sh # 35902941 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.sh # 35902943 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.sh # 35902945 ok

## qvalue 2.3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.sh # 38014865 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.sh # 38014873 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.sh # 38014875 xxx
#### KO
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.sh # 38014886 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.sh # 38014897 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.sh # 38014907 xxx
#### KOEF1aEZH1
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.sh # 38014939 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.sh # 38014967 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.sh # 38014991 xxx

## qvalue 3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.sh # 38015087 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.sh # 38015107 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.sh # 38015128 xxx
#### KO
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.sh # 38015162 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.sh # 38015212 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.sh # 38015235 xxx
#### KOEF1aEZH1
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.sh # 38015308 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.sh # 38015343 xxx
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.sh # 38015351 xxx





### EZH2
## Raw - non qvalue filtered ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_006R-FergusonUniqueNorm99.sh # 36718018 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_010R-FergusonUniqueNorm99.sh # 36718023 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.sh # 36718029 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.sh # 36718032 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.sh # 36718036 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.sh # 36718041 ok
#### KOEF1aEZH1
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.sh # 36718053 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.sh # 36718078 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.sh # 36718139 ok

## qvalue 2.3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_006R-FergusonUniqueNorm99.sh # 38015714 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_010R-FergusonUniqueNorm99.sh # 38015834 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.sh # 38015978 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.sh # 38015991 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.sh # 38015994 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.sh # 38015995 ok
#### KOEF1aEZH1
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.sh # 38016016 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.sh # 38016020 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.sh # 38016024 ok

## qvalue 3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_006R-FergusonUniqueNorm99.sh # 38016116 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_010R-FergusonUniqueNorm99.sh # 38016125 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.sh # 38016129 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.sh # 38016131 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.sh # 38016132 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.sh # 38016138 ok
#### KOEF1aEZH1
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.sh # 38016143 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.sh # 38016145 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.sh # 38016150 ok






### SUZ12
#### WT
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_WT_SUZ12_006R-FergusonUniqueNorm99.sh # 36717819 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_WT_SUZ12_013R1-FergusonUniqueNorm99.sh # 36717847 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_WT_SUZ12_014R1-FergusonUniqueNorm99.sh # 36717854 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_KO_SUZ12_013R1-FergusonUniqueNorm99.sh # 36717896 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_KO_SUZ12_014R1-FergusonUniqueNorm99.sh # 36717925 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_KO_SUZ12_014R2-FergusonUniqueNorm99.sh # 36717942 ok
#### KOEF1aEZH1
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_SUZ12_005R-FergusonUniqueNorm99.sh # 36717958 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_SUZ12_006R-FergusonUniqueNorm99.sh # 36717971 ok
sbatch scripts/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_SUZ12_013R1-FergusonUniqueNorm99.sh # 36717989 ok




```

--> I set here `--binSize 100 --regionBodyLength 100`; seems it give 1 value per row/peak. Look good.



### H3K27me3, no extension, Raw - non qvalue filtered - R DESEQ2


```R
library("tidyverse")
library("DESeq2")
library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
  
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())


# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R ) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")


SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1)



######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO %>%
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
colData_WTvsKO_raw <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO %>%
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
keep <- rowSums(counts(dds)) >= 100 # below 2000 look like noise on IGV
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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)





######################################################
### WT vs KOEF1aEZH1 ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KOEF1aEZH1"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKOEF1aEZH1 <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKOEF1aEZH1, -peakID), pull(countData_WTvsKOEF1aEZH1, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKOEF1aEZH1_raw <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKOEF1aEZH1_raw, -sample), pull(colData_WTvsKOEF1aEZH1_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 100 # below 2000 look like noise on IGV
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")



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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot)
#export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)



pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)
```


### H3K27me3, no extension,  qvalue 2.3 - R DESEQ2

```bash
conda activate deseq2
```

```R
library("tidyverse")
library("DESeq2")
# library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot-qval2.30103.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
  
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())


# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R ) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")


SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1)



######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO %>%
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
colData_WTvsKO_raw <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO %>%
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
keep <- rowSums(counts(dds)) >= 5 # before was 100
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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)




######################################################
### WT vs KOEF1aEZH1 ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KOEF1aEZH1"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKOEF1aEZH1 <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKOEF1aEZH1, -peakID), pull(countData_WTvsKOEF1aEZH1, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKOEF1aEZH1_raw <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKOEF1aEZH1_raw, -sample), pull(colData_WTvsKOEF1aEZH1_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5 # before was 100
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")



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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot)
#export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)



pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)
```




### H3K27me3, no extension,  qvalue 3 - R DESEQ2

```bash
conda activate deseq2
```

```R
library("tidyverse")
library("DESeq2")
# library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot-qval3.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
  
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())


# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R ) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1 %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")


SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1 = SCORE_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_006R %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_010R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_WT_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KO_H3K27me3_014R2) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks__PSC_KOEF1aEZH1_H3K27me3_013R1)



######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO %>%
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
colData_WTvsKO_raw <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKO %>%
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
keep <- rowSums(counts(dds)) >= 5 # before was 100
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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KO_vs_PSC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)




######################################################
### WT vs KOEF1aEZH1 ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KOEF1aEZH1"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKOEF1aEZH1 <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKOEF1aEZH1, -peakID), pull(countData_WTvsKOEF1aEZH1, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKOEF1aEZH1_raw <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKOEF1aEZH1_raw, -sample), pull(colData_WTvsKOEF1aEZH1_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5 # before was 100
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")



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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot)
#export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)



pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)
```







### EZH2, no extension,  Raw - non qvalue filtered- R DESEQ2




```R
library("tidyverse")
library("DESeq2")
#library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_010R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_010R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R ) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")


SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1)




######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO %>%
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
colData_WTvsKO_raw <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO %>%
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
keep <- rowSums(counts(dds)) >= 1500 # below 2000 look like noise on IGV
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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_EZH2_pool_peaks-PSC_KO_vs_PSC_WT-EZH2.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_EZH2_pool_peaks-PSC_KO_vs_PSC_WT-EZH2.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, EZH2',
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


XXXXXXXXX BELOW NOT MOPDIFIED XXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXX


upregulated_genes <- sum(res_tibble$log2FoldChange > 0.1 & res_tibble$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_tibble$log2FoldChange < -0.1 & res_tibble$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.1 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.1 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)





######################################################
### WT vs KOEF1aEZH1 ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 = SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks %>%
  filter(genotype %in% c("WT", "KOEF1aEZH1"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKOEF1aEZH1 <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKOEF1aEZH1, -peakID), pull(countData_WTvsKOEF1aEZH1, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKOEF1aEZH1_raw <- SCORE_BED_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_WTvsKOEF1aEZH1 %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKOEF1aEZH1_raw, -sample), pull(colData_WTvsKOEF1aEZH1_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 100 # below 2000 look like noise on IGV
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")



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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot)
#export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)



pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)



```

--> For EZH2, it seems the peak regions are sometime noise; add *noise to the DGB analysis, perform badly; we compare noise to noise*. Let's be more stringeant on peak identification. **Increase qvalue at MACS2**





### EZH2, no extension,  qvalue 2.3 - R DESEQ2




```R
library("tidyverse")
library("DESeq2")
#library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot-qval2.30103.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_010R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_010R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval2.30103-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R ) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")


SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1)




######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO %>%
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
colData_WTvsKO_raw <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO %>%
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
keep <- rowSums(counts(dds)) >= 5 # befor 5
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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-EZH2.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-EZH2.pdf", width=3, height=4)    
pdf("output/edgeR/plotVolcano_res_q05fc01-test.pdf", width=3, height=4)    

EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, EZH2',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KO_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)





######################################################
### WT vs KOEF1aEZH1 ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKOEF1aEZH1 = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks %>%
  filter(genotype %in% c("WT", "KOEF1aEZH1"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKOEF1aEZH1 <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKOEF1aEZH1 %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKOEF1aEZH1, -peakID), pull(countData_WTvsKOEF1aEZH1, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKOEF1aEZH1_raw <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKOEF1aEZH1 %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKOEF1aEZH1_raw, -sample), pull(colData_WTvsKOEF1aEZH1_raw, sample))

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

res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")



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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot)
#export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.txt", sep="\t", row.names=FALSE, quote=FALSE)



pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs WT, PSC, EZH2',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval2.30103-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)



```







### EZH2, no extension,  qvalue 3 - R DESEQ2




```R
library("tidyverse")
library("DESeq2")
#library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot <- read.delim("output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot-qval3.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)


# import SCORE 

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_010R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_010R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_WT_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KO_EZH2_014R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge-qval3-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, scoer per row, coordinate and row


SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R ) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1 %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")


SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1 = SCORE_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1  %>%
  left_join(BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1) %>%
  left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_006R %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_010R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_WT_EZH2_014R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KO_EZH2_014R2) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_006R) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_013R1) %>%
  bind_rows(SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks__PSC_KOEF1aEZH1_EZH2_014R1)




######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO %>%
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
colData_WTvsKO_raw <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKO %>%
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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot)
# Export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KO_vs_PSC_WT-EZH2.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KO_vs_PSC_WT-EZH2.pdf", width=3, height=4)    

EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, EZH2',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KO_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KO_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)





######################################################
### WT vs KOEF1aEZH1 ####################################
######################################################

SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKOEF1aEZH1 = SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks %>%
  filter(genotype %in% c("WT", "KOEF1aEZH1"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKOEF1aEZH1 <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKOEF1aEZH1 %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKOEF1aEZH1, -peakID), pull(countData_WTvsKOEF1aEZH1, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKOEF1aEZH1_raw <- SCORE_BED_WTKOKOEF1aEZH1_EZH2_pool_peaks_WTvsKOEF1aEZH1 %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKOEF1aEZH1_raw, -sample), pull(colData_WTvsKOEF1aEZH1_raw, sample))

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

res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")



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


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot)
#export result
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.txt", sep="\t", row.names=FALSE, quote=FALSE)



pdf("output/edgeR/plotVolcano_res_q05fc01-WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs WT, PSC, EZH2',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_WTKOKOEF1aEZH1_EZH2_pool_peaks-qval3-PSC_KOEF1aEZH1_vs_PSC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.1)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.1)



```






## Diff binding on all genes


```bash
# Calculate length-normalize signal for each gene promoters (1kb up 250bp down); here no lenght normalization as all region are same length
conda activate deeptools
sbatch scripts/LengthNormSignal_prom1kb250bp-FergusonUniqueNorm99smooth50bp.sh # 35838217 ok

# condition per condition (not optimal in the end, as I need replicate for stats)
sbatch scripts/LengthNormSignal_prom1kb250bp-WT_H3K27me3-FergusonUniqueNorm99smooth50bp.sh # 35849866 ok; fail need --missingDataAsZero 35891003; fail need --averageTypeBins sum  and --binSize 1250 (size promoter region!) 35893049 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-KO_H3K27me3-FergusonUniqueNorm99smooth50bp.sh # 35849895 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-KOEF1aEZH1_H3K27me3-FergusonUniqueNorm99smooth50bp.sh # 35849902 ok

sbatch scripts/LengthNormSignal_prom1kb250bp-WT_SUZ12-FergusonUniqueNorm99smooth50bp.sh # 35849929 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-KO_SUZ12-FergusonUniqueNorm99smooth50bp.sh # 35849939 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-KOEF1aEZH1_SUZ12-FergusonUniqueNorm99smooth50bp.sh # 35849944 ok

sbatch scripts/LengthNormSignal_prom1kb250bp-WT_EZH2-FergusonUniqueNorm99smooth50bp.sh # 35849953 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-KO_EZH2-FergusonUniqueNorm99smooth50bp.sh # 35849959 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-KOEF1aEZH1_EZH2-FergusonUniqueNorm99smooth50bp.sh # 35849971 ok

# sample per sample (replicate per replicate)
## H3K27me3
### WT
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.sh # 35894716 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.sh # 35894772 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.sh # 35894791 ok
### KO
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.sh # 35894808 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.sh # 35894813 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99smooth50bp.sh # 35894816 ok
### KOEF1aEZH1
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.sh # 35894822 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.sh # 35894834 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.sh # 35894855 ok


## SUZ12
### WT
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_SUZ12_006R-FergusonUniqueNorm99smooth50bp.sh # 35895792 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_SUZ12_013R1-FergusonUniqueNorm99smooth50bp.sh # 35895818 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_SUZ12_014R1-FergusonUniqueNorm99smooth50bp.sh # 35895823 ok
### KO
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_SUZ12_013R1-FergusonUniqueNorm99smooth50bp.sh # 35895847 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_SUZ12_014R1-FergusonUniqueNorm99smooth50bp.sh # 35895855 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_SUZ12_014R2-FergusonUniqueNorm99smooth50bp.sh # 35895870 ok
### KOEF1aEZH1
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_SUZ12_005R-FergusonUniqueNorm99smooth50bp.sh # 35895881 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_SUZ12_006R-FergusonUniqueNorm99smooth50bp.sh # 35895893 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_SUZ12_013R1-FergusonUniqueNorm99smooth50bp.sh # 35895913 ok


## EZH2
### WT
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_EZH2_006R-FergusonUniqueNorm99smooth50bp.sh # 35896656 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_EZH2_010R-FergusonUniqueNorm99smooth50bp.sh # 35896683 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_WT_EZH2_014R1-FergusonUniqueNorm99smooth50bp.sh # 35896690 ok
### KO
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_EZH2_013R1-FergusonUniqueNorm99smooth50bp.sh # 35896710 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_EZH2_014R1-FergusonUniqueNorm99smooth50bp.sh # 35896727 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KO_EZH2_014R2-FergusonUniqueNorm99smooth50bp.sh # 35896739 ok
### KOEF1aEZH1
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_EZH2_006R-FergusonUniqueNorm99smooth50bp.sh # 35896749 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_EZH2_013R1-FergusonUniqueNorm99smooth50bp.sh # 35896760 ok
sbatch scripts/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_EZH2_014R1-FergusonUniqueNorm99smooth50bp.sh # 35896775 ok



# Calculate length-normalize signal for each gene promoters (2.5kb up 2.5kb down); here no lenght normalization as all region are same length = 5kb
## H3K27me3
### WT
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.sh # 35902595 ok
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.sh # 35902601 ok
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.sh # 35902605 ok
### KO
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.sh # 35902609 ok
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.sh # 35902612 ok
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99smooth50bp.sh # 35902621 ok
### KOEF1aEZH1
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.sh # 35902624 ok
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.sh # 35902631 ok
sbatch scripts/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.sh # 35902638 ok

## SUZ12
XXX Maybe another region size is more adapted
## EZH2
XXX Maybe another region size is more adapted



```

--> So I had to add the following parameters to make it work. Here all regions from my bed are around promoter and are 1250bp length:
  -  `--missingDataAsZero` = replace NA per 0 when no data
  -  `--averageTypeBins sum` = here the bigwig are normalized, so I want the sum per bin of all my signal (does not matter because I do have 1 bin)
  -  `--binSize 1250` = Sets the bin size to exactly the length of each region (1250bp). This ensures that each gene is represented by only one bin, rather than being split into multiple bins.
  -  `--regionBodyLength 1250` = Defines the length to which each region should be scaled. Since all regions in my BED file are already 1250bp, I use the exact same value here to prevent any artificial resizing.

    --> Then, the `--outFileNameMatrix` contain one value per row = signal score per region. And `--outFileSortedRegions` contains the bed file region --> Usefull to group gene coordinate and signal 



--> According to deepTool plot, **optiomal size to calculate signal around gene TSS**:
  - *H3K27me3*: Most signal around -2.5 and +2.5kb; new bed file generated in `001*/009*`: `/scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_geneSymbol_prom2500bp2500bp.bed`
  - *SUZ12*: XXX
  - *EZH2*: XXX

Let's load the count matrix for all genes into R and perform Diff Bind analysis with edgeR



### H3K27me3

```bash
conda activate monocle3
```

#### H3K27me3 - 1kb upstream 250bp downstream TSS - R DESEQ2

```R
library("tidyverse")
library("DESeq2")
library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
ENCFF159KBI_geneSymbol_prom1kb250bp <- read.delim("../../Master/meta/ENCFF159KBI_geneSymbol_prom1kb250bp.bed", header=FALSE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = V1, start =V2 , end = V3, geneSymbol= V4)


# import SCORE 
SCORE_prom1kb250bp_PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())


# import BED position from matrix
BED_prom1kb250bp_PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom1kb250bp-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())




# Put together, gene name, scoer per row, coordinate and row

SCORE_BED_geneSymbol__prom1kb250bp_WT_H3K27me3_006R = SCORE_prom1kb250bp_PSC_WT_H3K27me3_006R %>%
  left_join(BED_prom1kb250bp_PSC_WT_H3K27me3_006R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_geneSymbol__prom1kb250bp_WT_H3K27me3_010R = SCORE_prom1kb250bp_PSC_WT_H3K27me3_010R %>%
  left_join(BED_prom1kb250bp_PSC_WT_H3K27me3_010R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_geneSymbol__prom1kb250bp_WT_H3K27me3_013R1 = SCORE_prom1kb250bp_PSC_WT_H3K27me3_013R1 %>%
  left_join(BED_prom1kb250bp_PSC_WT_H3K27me3_013R1) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")

SCORE_BED_geneSymbol__prom1kb250bp_KO_H3K27me3_006R = SCORE_prom1kb250bp_PSC_KO_H3K27me3_006R  %>%
  left_join(BED_prom1kb250bp_PSC_KO_H3K27me3_006R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_geneSymbol__prom1kb250bp_KO_H3K27me3_013R1 = SCORE_prom1kb250bp_PSC_KO_H3K27me3_013R1  %>%
  left_join(BED_prom1kb250bp_PSC_KO_H3K27me3_013R1) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_geneSymbol__prom1kb250bp_KO_H3K27me3_014R2 = SCORE_prom1kb250bp_PSC_KO_H3K27me3_014R2  %>%
  left_join(BED_prom1kb250bp_PSC_KO_H3K27me3_014R2) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_geneSymbol__prom1kb250bp_KOEF1aEZH1_H3K27me3_005R = SCORE_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_005R  %>%
  left_join(BED_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_geneSymbol__prom1kb250bp_KOEF1aEZH1_H3K27me3_006R = SCORE_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_006R  %>%
  left_join(BED_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_geneSymbol__prom1kb250bp_KOEF1aEZH1_H3K27me3_013R1 = SCORE_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_013R1  %>%
  left_join(BED_prom1kb250bp_PSC_KOEF1aEZH1_H3K27me3_013R1) %>%
  left_join(ENCFF159KBI_geneSymbol_prom1kb250bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_geneSymbol__prom1kb250bp_H3K27me3 = SCORE_BED_geneSymbol__prom1kb250bp_WT_H3K27me3_006R %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_WT_H3K27me3_010R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_WT_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_KO_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_KO_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_KO_H3K27me3_014R2) %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_KOEF1aEZH1_H3K27me3_005R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_KOEF1aEZH1_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom1kb250bp_KOEF1aEZH1_H3K27me3_013R1)



######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_geneSymbol__prom1kb250bp_H3K27me3_WTvsKO = SCORE_BED_geneSymbol__prom1kb250bp_H3K27me3 %>%
  filter(genotype %in% c("WT", "KO")) %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_geneSymbol__prom1kb250bp_H3K27me3_WTvsKO %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -geneSymbol), pull(countData_WTvsKO, geneSymbol)) 


## Create colData file that describe all our samples
colData_WTvsKO_raw <- SCORE_BED_geneSymbol__prom1kb250bp_H3K27me3_WTvsKO %>%
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
keep <- rowSums(counts(dds)) >= 100 # below 2000 look like noise on IGV
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")

## Export result as 'raw_NPC_KO_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/edgeR/raw_PSC_KO_vs_PSC_WT_H3K27me3.txt")
### If need to import: res <- read_csv("output/edgeR/raw_PSC_KO_vs_PSC_WT_H3K27me3.txt") #To import



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


res_tibble <- as_tibble(res, rownames = "geneSymbol")


pdf("output/edgeR/plotVolcano_res_q05fc01_PSC_KO_vs_PSC_WT_H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc01_PSC_KO_vs_PSC_WT_H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc01_PSC_KO_vs_PSC_WT_H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





```

--> Some geneSymbol were duplicated (likely different isoforms), for these I took median signal at the step  `# Put together, gene name, scoer per row, coordinate and row`






#### H3K27me3 - 2.5kb upstream 2.5kb downstream TSS - R DESEQ2

```R
library("tidyverse")
library("DESeq2")
library("edgeR")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
ENCFF159KBI_geneSymbol_prom2500bp2500bp <- read.delim("../../Master/meta/ENCFF159KBI_geneSymbol_prom2500bp2500bp.bed", header=FALSE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = V1, start =V2 , end = V3, geneSymbol= V4)


# import SCORE 
SCORE_prom2500bp2500bp_PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())


# import BED position from matrix
BED_prom2500bp2500bp_PSC_WT_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_WT_H3K27me3_010R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_WT_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_KO_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_KO_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_KO_H3K27me3_014R2 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KO_H3K27me3_014R2-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_005R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_006R <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_013R1 <- read.delim("output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_013R1-FergusonUniqueNorm99smooth50bp.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())




# Put together, gene name, scoer per row, coordinate and row

SCORE_BED_geneSymbol__prom2500bp2500bp_WT_H3K27me3_006R = SCORE_prom2500bp2500bp_PSC_WT_H3K27me3_006R %>%
  left_join(BED_prom2500bp2500bp_PSC_WT_H3K27me3_006R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_geneSymbol__prom2500bp2500bp_WT_H3K27me3_010R = SCORE_prom2500bp2500bp_PSC_WT_H3K27me3_010R %>%
  left_join(BED_prom2500bp2500bp_PSC_WT_H3K27me3_010R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_geneSymbol__prom2500bp2500bp_WT_H3K27me3_013R1 = SCORE_prom2500bp2500bp_PSC_WT_H3K27me3_013R1 %>%
  left_join(BED_prom2500bp2500bp_PSC_WT_H3K27me3_013R1) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")

SCORE_BED_geneSymbol__prom2500bp2500bp_KO_H3K27me3_006R = SCORE_prom2500bp2500bp_PSC_KO_H3K27me3_006R  %>%
  left_join(BED_prom2500bp2500bp_PSC_KO_H3K27me3_006R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_geneSymbol__prom2500bp2500bp_KO_H3K27me3_013R1 = SCORE_prom2500bp2500bp_PSC_KO_H3K27me3_013R1  %>%
  left_join(BED_prom2500bp2500bp_PSC_KO_H3K27me3_013R1) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_geneSymbol__prom2500bp2500bp_KO_H3K27me3_014R2 = SCORE_prom2500bp2500bp_PSC_KO_H3K27me3_014R2  %>%
  left_join(BED_prom2500bp2500bp_PSC_KO_H3K27me3_014R2) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_geneSymbol__prom2500bp2500bp_KOEF1aEZH1_H3K27me3_005R = SCORE_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_005R  %>%
  left_join(BED_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_005R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R1")
SCORE_BED_geneSymbol__prom2500bp2500bp_KOEF1aEZH1_H3K27me3_006R = SCORE_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_006R  %>%
  left_join(BED_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_006R) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R2")
SCORE_BED_geneSymbol__prom2500bp2500bp_KOEF1aEZH1_H3K27me3_013R1 = SCORE_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_013R1  %>%
  left_join(BED_prom2500bp2500bp_PSC_KOEF1aEZH1_H3K27me3_013R1) %>%
  left_join(ENCFF159KBI_geneSymbol_prom2500bp2500bp) %>%
  dplyr::select(geneSymbol, score) %>%
  group_by(geneSymbol) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KOEF1aEZH1", replicate = "R3")



# Tidy into a single tibble
SCORE_BED_geneSymbol__prom2500bp2500bp_H3K27me3 = SCORE_BED_geneSymbol__prom2500bp2500bp_WT_H3K27me3_006R %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_WT_H3K27me3_010R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_WT_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_KO_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_KO_H3K27me3_013R1) %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_KO_H3K27me3_014R2) %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_KOEF1aEZH1_H3K27me3_005R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_KOEF1aEZH1_H3K27me3_006R) %>%
  bind_rows(SCORE_BED_geneSymbol__prom2500bp2500bp_KOEF1aEZH1_H3K27me3_013R1)



######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_geneSymbol__prom2500bp2500bp_H3K27me3_WTvsKO = SCORE_BED_geneSymbol__prom2500bp2500bp_H3K27me3 %>%
  filter(genotype %in% c("WT", "KO")) %>%
  mutate(median_score = round(median_score))


# Convert to wide format
countData_WTvsKO <- SCORE_BED_geneSymbol__prom2500bp2500bp_H3K27me3_WTvsKO %>%
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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -geneSymbol), pull(countData_WTvsKO, geneSymbol)) 


## Create colData file that describe all our samples
colData_WTvsKO_raw <- SCORE_BED_geneSymbol__prom2500bp2500bp_H3K27me3_WTvsKO %>%
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
keep <- rowSums(counts(dds)) >= 100 # below 2000 look like noise on IGV
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")

## Export result as 'raw_NPC_KO_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/edgeR/raw_prom2500bp2500bp_PSC_KO_vs_PSC_WT_H3K27me3.txt")
### If need to import: res <- read_csv("output/edgeR/raw_PSC_KO_vs_PSC_WT_H3K27me3.txt") #To import



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


res_tibble <- as_tibble(res, rownames = "geneSymbol")


pdf("output/edgeR/plotVolcano_res_prom2500bp2500bp_q05fc01_PSC_KO_vs_PSC_WT_H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC, H3K27me3',
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
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_prom2500bp2500bp_q05fc01_PSC_KO_vs_PSC_WT_H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_prom2500bp2500bp_q05fc01_PSC_KO_vs_PSC_WT_H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





```

--> Some geneSymbol were duplicated (likely different isoforms), for these I took median signal at the step  `# Put together, gene name, scoer per row, coordinate and row`







# deepTool plots

## THOR TMM with DEGs

Let's generate deepTool plot with THOR Default TMM method (no spike in) on DEGs from `015__RNAseq`


```bash
# gtf file generated in 015:
meta/ENCFF159KBI_upregulated_q05fc05_PSC_KO_vs_PSC_WT.gtf
meta/ENCFF159KBI_downregulated_q05fc05_PSC_KO_vs_PSC_WT.gtf

meta/ENCFF159KBI_upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.gtf
meta/ENCFF159KBI_downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.gtf

# WT vs KO ####################
## Keeping rep individuals, 1 plot per IP
sbatch scripts/matrix_TSS_5kb_PSC_EZH1_WTKO_DEGWTvsKOq05fc05_TMM.sh # 29756558 fail THOR fail before;
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKO_DEGWTvsKOq05fc05_TMM.sh # 29756559 xxx
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKO_DEGWTvsKOq05fc05_TMM.sh # 29756560 xxx
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3_WTKO_DEGWTvsKOq05fc05_TMM.sh # 29756561 xxx

sbatch scripts/matrix_TSS_5kb_PSC_EZH1_WTKO_DEGWTvsKOq05fc05_housekeepHOX.sh # 29817709 ok
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKO_DEGWTvsKOq05fc05_housekeepHOX.sh # 29817862 ok
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKO_DEGWTvsKOq05fc05_housekeepHOX.sh # 29817949 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3_WTKO_DEGWTvsKOq05fc05_housekeepHOX.sh # 29817749 ok
## Pulling rep, all IP together
sbatch scripts/matrix_TSS_5kb_PSC_EZH1EZH2SUZ12H3K27me3_WTKOmerge_DEGWTvsKOq05fc05_TMM.sh


# WT vs KOEF1aEZH1 ###############
## Keeping rep individuals, 1 plot per IP
sbatch scripts/matrix_TSS_5kb_PSC_EZH1_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_TMM.sh # 29759887 xxx
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_TMM.sh # 29759890 xxx
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_TMM.sh # 29759917 xxx
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_TMM.sh # 29759919 ok

sbatch scripts/matrix_TSS_5kb_PSC_EZH1_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_housekeepHOX.sh # 29818027 ok
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_housekeepHOX.sh # 29818108 ok
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_housekeepHOX.sh # 29818113 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3_WTKOEF1aEZH1_DEGWTvsKOEF1aEZH1q05fc05_housekeepHOX.sh # 29817853 ok
## Pulling rep, all IP together
sbatch scripts/matrix_TSS_5kb_PSC_EZH1EZH2SUZ12H3K27me3_WTKOEF1aEZH1merge_DEGWTvsKOEF1aEZH1q05fc05_TMM.sh
```

--> *TMM Default normalization* leads to very different Replicate! It does not work well at all...



## Bigwig Ferguson

### Median tracks

Let's generate median tracks for Ferguson and THOR_Ferguson bigwigs (Norm99 unique)


**Not recommended method here, see IMPORTANT NOTE**


**Run wiggletools:**
```bash
conda activate BedToBigwig

# Calculate median
## bigwig_Ferguson
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique-H3K27me3.sh # 35440488 ok
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique-EZH2.sh # 35440565 ok
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique-SUZ12.sh # 35440572 ok
## bigwig_THOR_Ferguson
sbatch scripts/bigwigmerge_THOR_FergusonUniqueNorm99-H3K27me3.sh # 35440624 ok
sbatch scripts/bigwigmerge_THOR_FergusonUniqueNorm99-EZH2.sh # 35440625 ok
sbatch scripts/bigwigmerge_THOR_FergusonUniqueNorm99-SUZ12.sh # 35440635 ok



## smooth bigwig - smooth50bp
### calculate bin signal with multiBigwigSummary
conda activate deeptools

sbatch scripts/bigwigsmooth_Norm99_Ferguson_unique.sh # 35640884 fail; miss some samples; rerun 35652328 ok
### Re-convert to bigwig
conda activate BedToBigwig

sbatch --dependency=afterany:35652328 scripts/bigwigsmooth_Norm99_Ferguson_unique_part2.sh # 35652519 ok

## smooth bigwig - smooth250bp
### calculate bin signal with multiBigwigSummary
conda activate deeptools

sbatch scripts/bigwigsmooth250bp_Norm99_Ferguson_unique.sh # 39632034 ok
### Re-convert to bigwig
conda activate BedToBigwig

sbatch --dependency=afterany:39632034 scripts/bigwigsmooth250bp_Norm99_Ferguson_unique_part2.sh # 39632090 ok
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*

-->  Smoothing bigwig using `multiBigwigSummary` work great! 

--> ***IMPORTANT NOTE**: with this method, sharp CutRun like EZH2 SUZ12, signal is not smooth enough, too much sharp. That is because I used the 1bp resolution and did median on them, and then perform the smoothing. I should do instead: perform the smoothing from 1bp to 50bp (or 250bp), and then generate the median. --> File name as `-IndividualSampleSmoothing`*
  --> Not better in the end; still very sharp...

Let's do **smoothing per sample and then median**:



```bash

## smooth bigwig - smooth50bp SAMPLE PER SAMPLE
### calculate bin signal with multiBigwigSummary
conda activate deeptools

sbatch scripts/bigwigsmooth50bp_Norm99_Ferguson_unique-IndividualSampleSmoothing-H3K27me3.sh # 39634388 ok
sbatch scripts/bigwigsmooth50bp_Norm99_Ferguson_unique-IndividualSampleSmoothing-EZH2.sh # 39634451 ok
sbatch scripts/bigwigsmooth50bp_Norm99_Ferguson_unique-IndividualSampleSmoothing-SUZ12.sh # 39634541 ok



### Re-convert to bigwig
conda activate BedToBigwig

sbatch --dependency=afterany:39634388 scripts/bigwigsmooth50bp_Norm99_Ferguson_unique-IndividualSampleSmoothing-H3K27me3_part2.sh # 39634917 ok
sbatch --dependency=afterany:39634451 scripts/bigwigsmooth50bp_Norm99_Ferguson_unique-IndividualSampleSmoothing-EZH2_part2.sh # 39635023 ok
sbatch --dependency=afterany:39634541 scripts/bigwigsmooth50bp_Norm99_Ferguson_unique-IndividualSampleSmoothing-SUZ12_part2.sh # 39635126 ok
```


--> Its not better.. Still very sharp... 
  --> I found the issue, I was applying SF to local maxima bed... And not the initial bigwig. Correct one below:


```bash
conda activate BedToBigwig


# Calculate median
## bigwig initialBigwig
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique_initialBigwig-H3K27me3.sh # 39916201 ok
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique_initialBigwig-EZH2.sh # 39916207 ok
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique_initialBigwig-SUZ12.sh # 39916209 ok
```


--> Looks good!




### deepTool plots

Let's generate deepTool plot for all genes; peak/gene with changes of H3K27me3 levels.


```bash
conda activate deeptools
# All genes
## H3K27me3
sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-FergusonUniqueNorm99.sh # 35441130 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3_WTKOKOEF1aEZH1-FergusonUniqueNorm99.sh # 354411422 ok

sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.sh # 354411514 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.sh # 354411517 ok

sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp.sh # 35683527 ok


## SUZ12
sbatch scripts/matrix_TSS_10kb_PSC_SUZ12_WTKOKOEF1aEZH1-FergusonUniqueNorm99.sh # 354411439 ok
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-FergusonUniqueNorm99.sh # 354411452 ok

sbatch scripts/matrix_TSS_10kb_PSC_SUZ12_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.sh # 35441521 ok
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.sh # 35441525 ok

sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp.sh # 35683918 ok


## EZH2
sbatch scripts/matrix_TSS_10kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99.sh # 354411499 ok
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99.sh # 354411506 ok

sbatch scripts/matrix_TSS_10kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.sh # 35441551 ok
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.sh # 35441554 ok

sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp.sh # 35684066 ok


# Peak with THOR H3K27me3 changes
## separate peak gain lost
### WT vs KO
awk -F'\t' '$18 > 1' output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30.bed > output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_positive.bed
awk -F'\t' '$18 < 1' output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30.bed > output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_negative.bed
### WT vs KOEF1aEZH1
awk -F'\t' '$18 > 1' output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30.bed > output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_positive.bed
awk -F'\t' '$18 < 1' output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30.bed > output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_negative.bed

## Check signal in region/peak with H3K27me3 changes btwn WT vs KO
sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-THOR_FergusonUniqueNorm99.sh # 35441983 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-THOR_FergusonUniqueNorm99.sh # 35441990 ok

sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-FergusonUniqueNorm99.sh # 35442396 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-FergusonUniqueNorm99.sh # 35442397 ok

sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-FergusonUniqueNorm99smooth50bp.sh # 35685990 ok



## Check signal in region/peak with H3K27me3 changes btwn WT vs KOEF1aEZH1
sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-THOR_FergusonUniqueNorm99.sh # 35442072 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-THOR_FergusonUniqueNorm99.sh # 35442074 ok

sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99.sh # 35442504 ok
sbatch scripts/matrix_TSS_5kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99.sh # 35442538 ok

sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp.sh # 35685802 ok



# All macs2 WT H3K27me3 peaks (no trehsold)
sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-macs2_WT_H3K27me3_pool-THOR_FergusonUniqueNorm99smooth50bp.sh # 35696938 ok

# All consensus WTKOKOEF1aEZH1 peaks
## H3K27me3 consensus peak
sbatch scripts/matrix_TSS_10kb_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-FergusonUniqueNorm99smooth50bp_H3K27me3.sh # 36189289 ok
sbatch scripts/matrix_TSS_10kb_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.sh # 36391767 ok
sbatch scripts/matrix_TSS_10kb_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-FergusonUniqueNorm99smooth50bp_EZH2.sh # 36391796 ok
sbatch scripts/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.sh # 36392530 ok
sbatch scripts/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks-FergusonUniqueNorm99smooth50bp_EZH2.sh # 36392535 ok

## SUZ12 consensus peak
sbatch scripts/matrix_TSS_10kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_H3K27me3.sh # 36635203 ok
sbatch scripts/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.sh # 36636059 ok
sbatch scripts/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_EZH2.sh # 36636660 ok

## EZH2 consensus peak
sbatch scripts/matrix_TSS_10kb_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks-FergusonUniqueNorm99smooth50bp_H3K27me3.sh # 36637271 ok
sbatch scripts/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.sh # 36637940 ok
sbatch scripts/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks-FergusonUniqueNorm99smooth50bp_EZH2.sh # 36638223 ok






# Peak with DESEQ2 H3K27me3 changes
## separate peak gain lost
### WT vs KO
awk -F'\t' '$3 > 0.1 && $6 <= 0.05 {print $7, $8, $9, $1, $3, $6}' OFS='\t' output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt > output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3_positivepadj05FC01.bed # Filter log2FoldChange > 0.1 and padj <0.05 and print genomic coordinates and save as bed
awk -F'\t' '$3 < -0.1 && $6 <= 0.05 {print $7, $8, $9, $1, $3, $6}' OFS='\t' output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt > output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3_negativepadj05FC01.bed # Filter log2FoldChange > 0.1 and padj <0.05 and print genomic coordinates and save as bed

### WT vs KOEF1aEZH1
awk -F'\t' '$3 > 0.1 && $6 <= 0.05 {print $7, $8, $9, $1, $3, $6}' OFS='\t' output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt > output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3_positivepadj05FC01.bed # Filter log2FoldChange > 0.1 and padj <0.05 and print genomic coordinates and save as bed
awk -F'\t' '$3 < -0.1 && $6 <= 0.05 {print $7, $8, $9, $1, $3, $6}' OFS='\t' output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3.txt > output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3_negativepadj05FC01.bed # Filter log2FoldChange > 0.1 and padj <0.05 and print genomic coordinates and save as bed


## Check signal in region/peak with H3K27me3 changes btwn WT vs KO
sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.sh # 36688729 ok
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.sh # 36688748 ok
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.sh # 36688866 ok
## Check signal in region/peak with H3K27me3 changes btwn WT vs KOEF1aEZH1
sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.sh # 36689206 ok
sbatch scripts/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.sh # 36689213 ok
sbatch scripts/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.sh # 36689220 ok




## Check EZH1 and EZH2 signal in WT and OE in EZH1 peak (identified from KOEF1aEZH1)
sbatch scripts/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH1_EZH2-peakKOEF1aEZH1_EZH1macs2q3.sh # 42096114 ok
sbatch scripts/matrix_TSS_5kb_PSC_WTKO-initialBigwig_EZH1_EZH2-peakKOEF1aEZH1_EZH1macs2q3.sh # 42097105 xxx




```


--> Bigwig_Ferguson vs THOR_bigwig_Ferguson: Look mostly similar; but:
  - THOR plot is more smoothed, cleaner to show than bigwig_Ferguson
  - H3K27me3 and SUZ12 same
  - EZH2 is different!! THOR show less EZH2 in mutants; Ferguson show more...

  --> Lets prefer Ferguson method; conclusion more in agreement with our expectations. And we have more control on the normalization than THOR (not clear what THOR do). So generate smoothed Ferguson (bin 50 as in the paper)

--> I did not check THOR plot q30 carefully as NOT sure I can 'trust' THOR normalization method.


**Conclusion** (*Ferguson unique norm99*):
- More H3K27me3 in KO and KOEF1aEZH1
  - Changes in H3K27me3 seems to be occuring in both KO and KOEF1aEZH1 at same place (weird?)
- A bit more EZH2 in KO
  - notably when looking consensus EZH2 peaks



## QC PCA plots

Let's do PCA plots to compare the bigwig, how they perform to put same samples together.




```bash
conda activate deeptools

# pearson corr plots - bigwig regions GTF as BED
### All 3 genotypes
sbatch scripts/multiBigwigSummary_H3K27me3_THOR_TMM_bed.sh # 32405260 ok
sbatch scripts/multiBigwigSummary_H3K27me3_THOR_housekeepHOX_bed.sh # 32405464 ok
sbatch scripts/multiBigwigSummary_H3K27me3_THOR_DiffBindTMMEpiCypher_bed.sh # 32405785 ok

sbatch scripts/multiBigwigSummary_SUZ12_THOR_TMM_bed.sh # 32409383 ok
sbatch scripts/multiBigwigSummary_SUZ12_THOR_housekeepHOX_bed.sh # 32409439 ok
#sbatch scripts/multiBigwigSummary_SUZ12_THOR_DiffBindTMMEpiCypher_bed.sh # SUZ12_THOR_DiffBindTMMEpiCypher not generated

sbatch scripts/multiBigwigSummary_EZH2_THOR_TMM_bed.sh # 32409934 ok
sbatch scripts/multiBigwigSummary_EZH2_THOR_housekeepHOX_bed.sh # 32410183 ok
 #sbatch scripts/multiBigwigSummary_EZH2_THOR_DiffBindTMMEpiCypher_bed.sh # EZH2_THOR_DiffBindTMMEpiCypher not generated


# pearson corr plots - bigwig bins
### All 3 genotypes
sbatch scripts/multiBigwigSummary_H3K27me3_THOR_TMM.sh # 32396020 ok
sbatch scripts/multiBigwigSummary_H3K27me3_THOR_housekeepHOX.sh # 32396689 ok
sbatch scripts/multiBigwigSummary_H3K27me3_THOR_DiffBindTMMEpiCypher.sh # 32397002 ok

### 2 genotypes only
sbatch scripts/multiBigwigSummary_H3K27me3_WTvsKO_THOR_TMM.sh # 323400765 ok
sbatch scripts/multiBigwigSummary_H3K27me3_WTvsKO_THOR_housekeepHOX.sh # 323400857 ok
sbatch scripts/multiBigwigSummary_H3K27me3_WTvsKO_THOR_DiffBindTMMEpiCypher.sh # 323400961 ok

sbatch scripts/multiBigwigSummary_H3K27me3_WTvsKOEF1aEZH1_THOR_TMM.sh # 323401491 ok
sbatch scripts/multiBigwigSummary_H3K27me3_WTvsKOEF1aEZH1_THOR_housekeepHOX.sh # 323401667 ok
sbatch scripts/multiBigwigSummary_H3K27me3_WTvsKOEF1aEZH1_THOR_DiffBindTMMEpiCypher.sh # 323401695 ok
```

--> None normalization show good clustering... (WT and KOEF1aEZH1 OK; but KO_R2 weird outlier). Same when treating only two genotype at a time.
  --> Same result when using the bed-file (eg. clustering on gene only)...



## DIFFREPS

### peaks - DIFFREPS
Let's check gain / lost H3K27me3 regions in WT vs KO; and check H3K27me3 and EZH2 signal in:
- H3K27me3 binding changes peak
- H3K27me3 binding changes peak, overlapping with EZH2 consensus peaks (WT/KO)



```bash
conda activate deeptools


# PEAK
# DIFFREPS diff peaks
## isolate gain / lost peaks
### Parameters FAIL - bigwig local maxima...
awk -F'\t' '$4 == "Gain" ' output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO.txt > output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Gain.txt
awk -F'\t' '$4 == "Lost" ' output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO.txt > output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Lost.txt
### Good `*initialBigwig` - nb combined window
awk -F'\t' '$4 == "Gain" ' output/diffreps/merged_intervals-padj05_nb_pval0001-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt > output/diffreps/merged_intervals-padj05_nb_pval0001-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Gain.txt
awk -F'\t' '$4 == "Lost" ' output/diffreps/merged_intervals-padj05_nb_pval0001-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt > output/diffreps/merged_intervals-padj05_nb_pval0001-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Lost.txt

awk -F'\t' '$4 == "Gain" ' output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt > output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Gain.txt
awk -F'\t' '$4 == "Lost" ' output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt > output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Lost.txt





# DEEPTOOL PLOTS
## H3K27me3 binding changes
### Parameters FAIL - bigwig local maxima...
sbatch scripts/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.sh # interactive 
sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.sh # interactive 

### GOOD `*initialBigwig` - nb combined window
sbatch scripts/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak.sh # 39926889 ok 
sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak.sh # 39926902 ok 

sbatch scripts/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1-initialBigwig_H3K27me3_EZH2-peak.sh # 40025589 ok

sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1-initialBigwig_H3K27me3_EZH2-peak.sh # 40025596 ok


### Good `*initialBigwig` - gt pval05 padj 001
sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3_EZH2-peak.sh # 40810626 ok
sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3-bin1000space100_gt_pval05_padj001_FC1-initialBigwig_H3K27me3_EZH2-peak.sh # 40810627 ok



########################################################################################
## H3K27me3 binding changes overlapping with EZH2 consensus peak (WT/KO) ######################
########################################################################################

### Filter H3K27me3 binding changes to keep the one overlapping with EZH2
### Parameters FAIL - bigwig local maxima...
bedtools intersect -wa -a output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Gain.txt -b output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKO_EZH2_pool_peaks.sorted.merge.bed > output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Gain-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed
bedtools intersect -wa -a output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Lost.txt -b output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKO_EZH2_pool_peaks.sorted.merge.bed > output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Lost-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed
#### bigwig Ferguson smooth50bp
sbatch scripts/matrix_TSS_5kb-DIFFREPS-WTKO_H3K27me3_5kb2kb1kb500bp250bp__padj05_gt_pval05-macs2qval2.3_WTKO_EZH2-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.sh 

### GOOD `*initialBigwig`
conda activate BedToBigwig

bedtools intersect -wa -a output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Gain.txt -b output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKO_EZH2_pool_peaks.sorted.merge.bed > output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1_initialBigwig-5kb2kb1kb500bp250bp-WTvsKO-Gain-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed
bedtools intersect -wa -a output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Lost.txt -b output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKO_EZH2_pool_peaks.sorted.merge.bed > output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1_initialBigwig-5kb2kb1kb500bp250bp-WTvsKO-Lost-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed

sbatch scripts/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak.sh # 40044647 ok



```

--> H3K27me3 show clear changes, EZH2 follow the pattern only for Lost region. For gain, EZH2 not clean enough, or no EZH2 peak!
  --> Let's *keep Only H3K27me3 binding changes region that overlap with EZH2 consensus (WT and/or KO) peak*
  

--> Seems that my smooth50bp is not smooth enough, many sharp peaks; may affect vizualization; lets smooth to 250bp instead


### genes - DIFFREPS



**Generate GTF file from these geneSymbol list** of genes:


- DIFFREPS genes - NB default
- DIFFREPS gene overlapping with EZH2 consensus peak WT/KO
- DIFFREPS genes - GT pval05 padj 001 with without FC tresh 1



```bash
# Generate gtf file from gene list:

### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt

## Filter the gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5.gtf

grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5.gtf


grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5.gtf

grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5.gtf



# deeptool plots
sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1-initialBigwig_H3K27me3_EZH2-gene.sh # 40044612 ok

sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3-padj001_nb_pval0001_log2FC1-macs2qval2EZH2_WTKO-initialBigwig_H3K27me3_EZH2-gene.sh # 40045139 ok

sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3_EZH2-gene.sh # 40812453 ok
sbatch scripts/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3-bin1000space100_gt_pval05_padj001_FC1-initialBigwig_H3K27me3_EZH2-gene.sh # 40813239 ok




```






# Bigwig
## Generate raw bigwig from unique.dupmark.sorted


Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch scripts/bamtobigwig_unique_1.sh # 29756540 ok
sbatch scripts/bamtobigwig_unique_2.sh # 29756546 ok
sbatch scripts/bamtobigwig_unique_3.sh # 29756550 ok
```

## Generate norm bigwig from DiffBind SF

Let's use directly the SF from DiffBind (not the reciprocal).



```bash
conda activate deeptools

# directly SF from DiffBind
sbatch scripts/bamtobigwig-DBA_NORM_TMM_unique_SF_H3K27me3.sh # 32649115 ok
sbatch scripts/bamtobigwig-DBA_NORM_RLE_unique_SF_H3K27me3.sh # 32649206 ok
sbatch scripts/bamtobigwig-DBA_NORM_LIB_unique_SF_H3K27me3.sh # 32649234 ok

# Reciprocal SF (=1/SF) from DiffBind
sbatch scripts/bamtobigwig-DBA_NORM_TMM_unique_reciprocalSF_H3K27me3.sh # 32650322 ok

```

--> using SF directly from DiffBind is bad, very heterogeneous replicates, for each norm method... (Reciprocal also bad)

## Pearson correlation heatmap on bigwig signals

```bash
conda activate deeptools


# Generate compile bigwig (.npz) files
## H3K27me3
sbatch scripts/multiBigwigSummary_H3K27me3_unique_raw.sh # 35225590 ok
sbatch scripts/multiBigwigSummary_H3K27me3_unique_norm99.sh # 35225604 ok
## EZH2
sbatch scripts/multiBigwigSummary_EZH2_unique_raw.sh # 35225618 ok
sbatch scripts/multiBigwigSummary_EZH2_unique_norm99.sh # 35225652 ok
## SUZ12
sbatch scripts/multiBigwigSummary_SUZ12_unique_raw.sh # 35225673 ok
sbatch scripts/multiBigwigSummary_SUZ12_unique_norm99.sh # 35225733 ok




# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_H3K27me3_unique_raw.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_H3K27me3_005R PSC_KOEF1aEZH1_H3K27me3_006R PSC_KOEF1aEZH1_H3K27me3_013R1 PSC_KO_H3K27me3_006R PSC_KO_H3K27me3_013R1 PSC_KO_H3K27me3_014R2 PSC_WT_H3K27me3_006R PSC_WT_H3K27me3_010R PSC_WT_H3K27me3_013R1 \
    --markers o o o o o o o o o \
    --colors blue blue blue red red red black black black \
    -o output/bigwig/multiBigwigSummary_H3K27me3_unique_raw_plotPCA.pdf

plotPCA -in output/bigwig_Ferguson/multiBigwigSummary_H3K27me3_unique_norm99.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_H3K27me3_005R PSC_KOEF1aEZH1_H3K27me3_006R PSC_KOEF1aEZH1_H3K27me3_013R1 PSC_KO_H3K27me3_006R PSC_KO_H3K27me3_013R1 PSC_KO_H3K27me3_014R2 PSC_WT_H3K27me3_006R PSC_WT_H3K27me3_010R PSC_WT_H3K27me3_013R1 \
    --markers o o o o o o o o o \
    --colors blue blue blue red red red black black black \
    -o output/bigwig_Ferguson/multiBigwigSummary_H3K27me3_unique_norm99_plotPCA.pdf
```

--> Unique norm 99 whose work great show quite bad correlation in PCA plot... Maybe that's not a big deal...




## Merge bigiwg files, replicate

XXX


# peak calling
## MACS2 peak calling on bam unique


--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad` and `narrow` **


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_1.sh # 29773632 ok
sbatch scripts/macs2_broad_2.sh # 29773633 ok
sbatch scripts/macs2_broad_3.sh # 29773661 ok

# pool replicate
sbatch scripts/macs2_broad_pool_1.sh # 35693862 ok
sbatch scripts/macs2_broad_pool_2.sh # 35693955 ok
sbatch scripts/macs2_broad_pool_3.sh # 35694072 ok
```

--> all good



## MACS2 peak qvalue filtering

For **consensus peak** counting (ie Ferguson / local maxima method); I used the **pool peak** to generate the consensus peak file for the three genotype comparison

```bash
conda activate bowtie2 # for bedtools

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
- WT_H3K27me3: 2.3 or 3
- KO_H3K27me3: 2.3 or 3
- KOEF1aEZH1_H3K27me3: 3
--> Let's setup *optimal qvalue to 3 (look more true) for H3K27me3*
- WT_EZH2: 2.3 or 3
- KO_EZH2: 2.3 or 3
- KOEF1aEZH1_EZH2: 2.3 or 3
--> Let's setup *optimal qvalue to 3 (look more true) for EZH2*
- KOEF1aEZH1_EZH1: 2.3 or 3
--> Let's setup *optimal qvalue to 3 (look more true) for EZH1*












# ChIPseeker peak gene assignment

## From THOR diff peaks

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


# Import THOR peaks - housekeepHOX and DiffBindTMMEpiCypher
## H3K27me3
WTvsKO_H3K27me3_housekeepHOX = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/THOR_qval20.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
WTvsKOEF1aEZH1_H3K27me3_housekeepHOX = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX/THOR_qval20.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

WTvsKO_H3K27me3_DiffBindTMMEpiCypher = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMMEpiCypher/THOR_qval20.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher/THOR_qval20.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 




# Tidy peaks 
## H3K27me3
WTvsKO_H3K27me3_housekeepHOX_gr = makeGRangesFromDataFrame(WTvsKO_H3K27me3_housekeepHOX,keep.extra.columns=TRUE)
WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_gr = makeGRangesFromDataFrame(WTvsKOEF1aEZH1_H3K27me3_housekeepHOX,keep.extra.columns=TRUE)
WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gr = makeGRangesFromDataFrame(WTvsKO_H3K27me3_DiffBindTMMEpiCypher,keep.extra.columns=TRUE)
WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_gr = makeGRangesFromDataFrame(WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher,keep.extra.columns=TRUE)
gr_list <- list(WTvsKO_H3K27me3_housekeepHOX=WTvsKO_H3K27me3_housekeepHOX_gr, WTvsKOEF1aEZH1_H3K27me3_housekeepHOX=WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_gr,  WTvsKO_H3K27me3_DiffBindTMMEpiCypher=WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gr, WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher=WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_gr)



# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_H3K27me3_housekeepHOX_DiffBindTMMEpiCypher.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_H3K27me3_housekeepHOX_DiffBindTMMEpiCypher.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
WTvsKO_H3K27me3_housekeepHOX_annot <- as.data.frame(peakAnnoList[["WTvsKO_H3K27me3_housekeepHOX"]]@anno)
WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot <- as.data.frame(peakAnnoList[["WTvsKOEF1aEZH1_H3K27me3_housekeepHOX"]]@anno)
WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot <- as.data.frame(peakAnnoList[["WTvsKO_H3K27me3_DiffBindTMMEpiCypher"]]@anno)
WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot <- as.data.frame(peakAnnoList[["WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher"]]@anno)

## Convert entrez gene IDs to gene symbols
WTvsKO_H3K27me3_housekeepHOX_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_housekeepHOX_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKO_H3K27me3_housekeepHOX_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_housekeepHOX_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot, file="output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
WTvsKO_H3K27me3_housekeepHOX_annot_promoterAnd5 = tibble(WTvsKO_H3K27me3_housekeepHOX_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot_promoterAnd5 = tibble(WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5 = tibble(WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5 = tibble(WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
WTvsKO_H3K27me3_housekeepHOX_annot_promoterAnd5_geneSymbol = WTvsKO_H3K27me3_housekeepHOX_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot_promoterAnd5_geneSymbol = WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5_geneSymbol = WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5_geneSymbol = WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)





# Import THOR peaks - Ferguson Unique Norm 99 (no input)
## H3K27me3
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30 = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30 = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 


WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50 = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval50.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50 = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval50.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 


# Tidy peaks 
## H3K27me3
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_gr = makeGRangesFromDataFrame(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30,keep.extra.columns=TRUE)
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_gr = makeGRangesFromDataFrame(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30,keep.extra.columns=TRUE)

WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_gr = makeGRangesFromDataFrame(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50,keep.extra.columns=TRUE)
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_gr = makeGRangesFromDataFrame(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50,keep.extra.columns=TRUE)

gr_list <- list(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30=WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_gr, WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30=WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_gr, WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50=WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_gr, WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50=WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_H3K27me3_FergusonUniqueNorm99_noInput.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_H3K27me3_FergusonUniqueNorm99_noInput.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot <- as.data.frame(peakAnnoList[["WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30"]]@anno)
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot <- as.data.frame(peakAnnoList[["WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30"]]@anno)
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot <- as.data.frame(peakAnnoList[["WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50"]]@anno)
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot <- as.data.frame(peakAnnoList[["WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50"]]@anno)

## Convert entrez gene IDs to gene symbols
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$gene <- mapIds(org.Hs.eg.db, keys = WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot, file="output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot, file="output/ChIPseeker/annotation_THORq30_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", sep="\t", quote=F, row.names=F) 

write.table(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot, file="output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot, file="output/ChIPseeker/annotation_THORq50_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", sep="\t", quote=F, row.names=F) 

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5 = tibble(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5 = tibble(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5 = tibble(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5 = tibble(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5_geneSymbol = WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5_geneSymbol = WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5_geneSymbol = WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5_geneSymbol = WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq30_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORq30_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

write.table(WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_THORq50_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORq50_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

```





## From consensus peak

### Consensus peak H3K27me3 no extension/extension, Raw - non qvalue filtered
- consensus peak H3K27me3, WT KO KOEF1aEZH1, no extension: `output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed`
- consensus peak H3K27me3, WT KO KOEF1aEZH1, 100bp merge: `output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed`
- consensus peak H3K27me3, WT KO KOEF1aEZH1, 500bp merge: `output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed`


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
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K27me3
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_gr, PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_gr, PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge"]]@anno)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp"]]@anno)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```


### Consensus peak H3K27me3 no extension/extension, qvalue 2.3
- consensus peak H3K27me3, WT KO KOEF1aEZH1, no extension: `output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed`
- consensus peak H3K27me3, WT KO KOEF1aEZH1, 100bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed`
- consensus peak H3K27me3, WT KO KOEF1aEZH1, 500bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed`


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
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K27me3
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_gr, PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_gr, PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge-qval2.30103.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge-qval2.30103.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge"]]@anno)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp"]]@anno)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot-qval2.30103.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot-qval2.30103.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot-qval2.30103.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot-qval2.30103-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot-qval2.30103-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot-qval2.30103-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```


### Consensus peak H3K27me3 no extension/extension, qvalue 3
- consensus peak H3K27me3, WT KO KOEF1aEZH1, no extension: `output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed`
- consensus peak H3K27me3, WT KO KOEF1aEZH1, 100bp merge: `output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed`
- consensus peak H3K27me3, WT KO KOEF1aEZH1, 500bp merge: `output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed`


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
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K27me3
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_gr, PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_gr, PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp=PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge-qval3.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge-qval3.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge"]]@anno)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp"]]@anno)
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot-qval3.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot-qval3.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot-qval3.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge_annot-qval3-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge100bp_annot-qval3-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks_merge500bp_annot-qval3-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```





### Consensus peak EZH2 no extension/extension, Raw - non qvalue filtered
- consensus peak EZH2, WT KO KOEF1aEZH1, no extension: `output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed`
- consensus peak EZH2, WT KO KOEF1aEZH1, 100bp merge: `output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed`
- consensus peak EZH2, WT KO KOEF1aEZH1, 500bp merge: `output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed`


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
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K27me3
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_gr, PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_gr, PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge"]]@anno)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp"]]@anno)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```



### Consensus peak EZH2 no extension/extension, qvalue 2.3
- consensus peak EZH2, WT KO KOEF1aEZH1, no extension: `output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed`
- consensus peak EZH2, WT KO KOEF1aEZH1, 100bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed`
- consensus peak EZH2, WT KO KOEF1aEZH1, 500bp merge: `output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed`


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
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K27me3
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_gr, PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_gr, PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge-qval2.30103.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge-qval2.30103.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge"]]@anno)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp"]]@anno)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot-qval2.30103.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot-qval2.30103.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot-qval2.30103.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot-qval2.30103-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot-qval2.30103-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot-qval2.30103-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```




### Consensus peak EZH2 no extension/extension, qvalue 3
- consensus peak EZH2, WT KO KOEF1aEZH1, no extension: `output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed`
- consensus peak EZH2, WT KO KOEF1aEZH1, 100bp merge: `output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed`
- consensus peak EZH2, WT KO KOEF1aEZH1, 500bp merge: `output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed`


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
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge100bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed')) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3) 

# Tidy peaks 
## H3K27me3
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp,keep.extra.columns=TRUE)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_gr = makeGRangesFromDataFrame(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_gr, PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_gr, PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp=PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge-qval3.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge-qval3.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge"]]@anno)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp"]]@anno)
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot <- as.data.frame(peakAnnoList[["PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp"]]@anno)


## Convert entrez gene IDs to gene symbols
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot-qval3.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot-qval3.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot, file="output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot-qval3.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5 = tibble(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol = PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge_annot-qval3-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge100bp_annot-qval3-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks_merge500bp_annot-qval3-promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```


## From DIFFREPS


### On combine windows 5kb-250bp - initialBigwig

Let's assign peak to genes on the two best windowns/parameters:
- merge interval from **Negative binomial (nb) pval 0.0001 for window and padj 0.001**; of 5kb, 2kb, 1kb, 500bp, 250bp: `output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt`



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


# Import diff peaks
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1 <- read.delim("output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt", sep = "\t", header = TRUE) %>%
  as_tibble()
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1 %>%
  filter(direction == "Gain") %>% mutate(log2FC = as.numeric(log2FC)) %>% filter(log2FC > 1)
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1 %>%
  filter(direction == "Lost") %>% mutate(log2FC = as.numeric(log2FC)) %>% filter(log2FC < -1)
### SAVE Gain and Lost peaks
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain, file="output/diffreps/merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost, file="output/diffreps/merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost.txt", sep="\t", quote=F, row.names=F) 
########

# Tidy peaks 
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_gr = makeGRangesFromDataFrame(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain,keep.extra.columns=TRUE)
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_gr = makeGRangesFromDataFrame(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost,keep.extra.columns=TRUE)


gr_list <- list(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain=merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_gr,merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost=merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_initialBigwig.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_initialBigwig.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot <- as.data.frame(peakAnnoList[["merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain"]]@anno)
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot <- as.data.frame(peakAnnoList[["merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost"]]@anno)


## Convert entrez gene IDs to gene symbols
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot, file="output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_initialBigwig_log2FC1_Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot, file="output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_initialBigwig_log2FC1_Lost_annot.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5 = tibble(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5 = tibble(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5_geneSymbol = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5_geneSymbol = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



```




### On bin1000space100_gt_pval05_padj001 - WT vs KO H3K27me3

Let's assign peak to genes on the two best windowns/parameters as in `001*/009*`:
- Bin 1000bp space 100bp, G test, pval 0.05 and padj 0.001: done without FC and with FC 1 treshold



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


# Import diff peaks
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>%
  filter(padj < 0.001)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC>0)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC<0)

PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>%
  filter(padj < 0.001)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC>1)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC< (-1) )


### SAVE Gain and Lost peaks
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost.txt", sep="\t", quote=F, row.names=F) 

write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost.txt", sep="\t", quote=F, row.names=F) 
########

# Tidy peaks 
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain,keep.extra.columns=TRUE)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost,keep.extra.columns=TRUE)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain,keep.extra.columns=TRUE)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost,keep.extra.columns=TRUE)

gr_list <- list(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_gr,PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_gr, PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_gr,PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain"]]@anno)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost"]]@anno)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain"]]@anno)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost"]]@anno)

## Convert entrez gene IDs to gene symbols
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot.txt", sep="\t", quote=F, row.names=F)  


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



```






### On bin1000space100_gt_pval05_padj001 - WT vs KOEF1aEZH1 H3K27me3

Let's assign peak to genes on the two best windowns/parameters as in `001*/009*`:
- Bin 1000bp space 100bp, G test, pval 0.05 and padj 0.001: done without FC and with FC 1 treshold


XXXY


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


# Import diff peaks
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>%
  filter(padj < 0.001)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC>0)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC<0)

PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>%
  filter(padj < 0.001)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC>1)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001 %>%
  filter(log2FC< (-1) )


### SAVE Gain and Lost peaks
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost.txt", sep="\t", quote=F, row.names=F) 

write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost, file="output/diffreps/PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost.txt", sep="\t", quote=F, row.names=F) 
########

# Tidy peaks 
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain,keep.extra.columns=TRUE)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost,keep.extra.columns=TRUE)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain,keep.extra.columns=TRUE)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_gr = makeGRangesFromDataFrame(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost,keep.extra.columns=TRUE)

gr_list <- list(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_gr,PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_gr, PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_gr,PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost=PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain"]]@anno)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost"]]@anno)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain"]]@anno)
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost"]]@anno)

## Convert entrez gene IDs to gene symbols
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot, file="output/ChIPseeker/annotation_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot.txt", sep="\t", quote=F, row.names=F)  


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5 = tibble(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_geneSymbol = PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001_FC1__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



```







## DIFFREPS overlapping with EZH2 consensus peaks




output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1_initialBigwig-5kb2kb1kb500bp250bp-WTvsKO-Lost-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed




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


# Import diff peaks
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain <- read.delim("output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1_initialBigwig-5kb2kb1kb500bp250bp-WTvsKO-Gain-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename(Chr = V1, Start = V2, End = V3)
  merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost <- read.delim("output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1_initialBigwig-5kb2kb1kb500bp250bp-WTvsKO-Lost-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename(Chr = V1, Start = V2, End = V3)


# Tidy peaks 
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_gr = makeGRangesFromDataFrame(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain,keep.extra.columns=TRUE)
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_gr = makeGRangesFromDataFrame(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost,keep.extra.columns=TRUE)


gr_list <- list(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain=merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_gr,merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost=merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_initialBigwig__macs2qval2EZH2_WTKO.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_initialBigwig__macs2qval2EZH2_WTKO.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot <- as.data.frame(peakAnnoList[["merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain"]]@anno)
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot <- as.data.frame(peakAnnoList[["merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost"]]@anno)


## Convert entrez gene IDs to gene symbols
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot, file="output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_initialBigwig_Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot, file="output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_initialBigwig_Lost_annot.txt", sep="\t", quote=F, row.names=F)  



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5 = tibble(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5 = tibble(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5_geneSymbol = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5_geneSymbol = merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1__macs2qval2EZH2_WTKO_Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



```



# RNAseq integration 

## Using 001015 RNAseq (Jasmine)

--> DEGs is redo as in `015__RNAseq`


### PSC KO vs WT


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

## collect all samples ID
samples <- c("PSC_WT_R1", "PSC_WT_R2" ,"PSC_WT_R3" ,"PSC_KO_R1" ,"PSC_KO_R2", "PSC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../015__RNAseq_PSC/output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
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
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
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
res$geneSymbol <- gene_symbols



## import gene list gain / lost H3K27me3

################################################################################################
##### H3K27me3 ##################################################################################
################################################################################################

### housekeepHOX WT vs KO ###############################################################
THORq20_WTvsKO_H3K27me3_housekeepHOX_gain = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()
THORq20_WTvsKO_H3K27me3_housekeepHOX_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique() 

### DESEQ2 WT vs KO ###############################################################
DESEQ2_WTvsKO_H3K27me3_gain = read.table("output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(log2FoldChange > 0.1,
                                      padj < 0.05,
                                      annotation != "Distal intergenic") %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = DESEQ2_WTvsKO_H3K27me3_gain %>% 
  left_join(res_tibble) 

### LOST
### housekeepHOX WT vs KO ###############################################################
THORq20_WTvsKO_H3K27me3_housekeepHOX_lost = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKO_H3K27me3_housekeepHOX_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

### DESEQ2 WT vs KO ###############################################################
DESEQ2_WTvsKO_H3K27me3_lost = read.table("output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(log2FoldChange < -0.1,
                                      padj < 0.05,
                                      annotation != "Distal intergenic") %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = DESEQ2_WTvsKO_H3K27me3_lost %>% 
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


#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8) 
pdf("output/deseq2/plotVolcano_DESEQ2_WTvsKO_H3K27me3_gain_noIntergenic__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  

EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 10
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 13


## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)
pdf("output/deseq2/plotVolcano_DESEQ2_WTvsKO_H3K27me3_lost_noIntergenic__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  

EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 22
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 0


### spikein WT vs KO ###############################################################
THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5 %>% 
  left_join(res_tibble) 


### LOST

THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5 %>% 
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

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 2
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 7






## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 22
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 0



### Ferguson unique WT vs KO H3K27me3 (no input) ###############################################################
THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()
THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique() 

THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()
THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique() 

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain_promoterAnd5 %>% 
  left_join(res_tibble) 

### LOST
THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost_promoterAnd5 %>% 
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


#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  

EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 10
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 13


## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  

EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 22
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 0





```







### PSC KOEF1aEZH1 vs WT


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

## collect all samples ID
samples <- c("PSC_WT_R1", "PSC_WT_R2" ,"PSC_WT_R3" ,"PSC_KOEF1aEZH1_R1" ,"PSC_KOEF1aEZH1_R2", "PSC_KOEF1aEZH1_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../015__RNAseq_PSC/output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
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
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
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
res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$geneSymbol <- gene_symbols



## import gene list gain / lost H3K27me3

# H3K27me3
### housekeepHOX WT vs KOEF1aEZH1 ###############################################################
THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_gain = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_gain_promoterAnd5 %>% 
  left_join(res_tibble) 


### LOST

THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_lost = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_lost_promoterAnd5 %>% 
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

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_gain__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_gain_promoterAnd5__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KOEF1aEZH1 vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 10
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 13






## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_lost__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_housekeepHOX_lost_promoterAnd5__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KOEF1aEZH1 vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 22
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 0










# H3K27me3
### spikein WT vs KOEF1aEZH1 ###############################################################
THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_gain = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()




#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5 %>% 
  left_join(res_tibble) 


### LOST

THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_lost = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()



#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5 %>% 
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


#pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_gain__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KOEF1aEZH1 vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 2
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 7






## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_lost__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5__PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KOEF1aEZH1 vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 127
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 15


```




## Using 001001 RNAseq (Carolina)



### PSC KO vs WT


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

## collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1" ,"ESC_KO_R2" ,"ESC_KO_R3")

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
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
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
res$geneSymbol <- gene_symbols



## import gene list gain / lost H3K27me3

################################################################################################
##### H3K27me3 ##################################################################################
################################################################################################


### Ferguson unique WT vs KO H3K27me3 (no input) ###############################################################
THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()
THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique() 

THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()
THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique() 

### DESEQ2 WT vs KO ###############################################################
DESEQ2_WTvsKO_H3K27me3_gain = read.table("output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(log2FoldChange > 0.1,
                                      padj < 0.05,
                                      annotation != "Distal intergenic") %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = DESEQ2_WTvsKO_H3K27me3_gain %>% 
  left_join(res_tibble) 

### LOST
THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq30_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

### DESEQ2 WT vs KO ###############################################################
DESEQ2_WTvsKO_H3K27me3_lost = read.table("output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KO_vs_PSC_WT-H3K27me3.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(log2FoldChange < -0.1,
                                      padj < 0.05,
                                      annotation != "Distal intergenic") %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = DESEQ2_WTvsKO_H3K27me3_lost %>% 
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


#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_gain_promoterAnd5__ESC_KO_vs_ESC_WT_001001.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_DESEQ2_WTvsKO_H3K27me3__gain_noIntergenic__ESC_KO_vs_ESC_WT_001001.pdf", width=8, height=8)  

EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 10
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 13


## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
#pdf("output/deseq2/plotVolcano_THORq50_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput_lost_promoterAnd5__ESC_KO_vs_ESC_WT_001001.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_DESEQ2_WTvsKO_H3K27me3__lost_noIntergenic__ESC_KO_vs_ESC_WT_001001.pdf", width=8, height=8)  

EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 22
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 0






XXXY BELOW NOT MOD !!!!!!!!!!!!!!! XXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXX



### housekeepHOX WT vs KO ###############################################################
THORq20_WTvsKO_H3K27me3_housekeepHOX_gain = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()
THORq20_WTvsKO_H3K27me3_housekeepHOX_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique() 

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = THORq20_WTvsKO_H3K27me3_housekeepHOX_gain_promoterAnd5 %>% 
  left_join(res_tibble) 

### LOST
THORq20_WTvsKO_H3K27me3_housekeepHOX_lost = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKO_H3K27me3_housekeepHOX_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_housekeepHOX_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = THORq20_WTvsKO_H3K27me3_housekeepHOX_lost_promoterAnd5 %>% 
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


#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_gain_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 10
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 13


## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_housekeepHOX_lost_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 22
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 0


### spikein WT vs KO ###############################################################
THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 > 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5 %>% 
  left_join(res_tibble) 


### LOST

THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5 = read.table("output/ChIPseeker/annotation_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_annot.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
                               as_tibble() %>%
                               filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
                               filter(V18 < 1) %>% # FILTER FC positive here!!
                               dplyr::select(geneSymbol) %>%
                               unique()

#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5 %>% 
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

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_gain_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 2
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE) # 7






## PLOT
### LOST
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

#pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_THORq20_WTvsKO_H3K27me3_DiffBindTMMEpiCypher_lost_promoterAnd5__PSC_KO_vs_PSC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
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


upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 22
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE) # 0






```













# Ecoli scaling factor

From [EpiCypher](https://support.epicypher.com/docs/normalizing-to-e-coli-spike-in-dna), it seems uniquely aligned reads should be used, from both human, and E coli!

## Count the nb of uniquely aligned reads

Let's count the nb of uniquely aligned reads in our sample (human)

```bash
conda activate bowtie2

sbatch scripts/samtools_unique_count.sh # 29775711 ok
```

--> Output of read counts for all samples in `samples_001016.xlsx`

Let's check the nb and proportion of E coli vs human reads for each sample. Created file with output metrics in `output/QC/samples_001016_readCount.txt`

```R
# packages
library("tidyverse")
library("ggpubr")


# import
readCount = read_tsv("output/QC/samples_001016_readCount.txt") %>%
  filter(uniqMappedReads >0 ) %>%
  tidyr::pivot_longer(cols = c(uniqMappedReads, uniqMappedMG1655Reads), 
                      names_to = "Read_Type", values_to = "Read_Count") %>%
  group_by(new_sample_ID) %>%
  mutate(Ecoli_prop = round((Read_Count[Read_Type == "uniqMappedMG1655Reads"] / 
                             Read_Count[Read_Type == "uniqMappedReads"]) * 100 , 3)) %>%
  tidyr::pivot_longer(cols = c(Read_Count), 
                      names_to = "Read_Measure", values_to = "Count_Value") %>%
  mutate(new_sample_ID = factor(new_sample_ID, levels = unique(new_sample_ID[order(exp)])))
 

sample_order <- c(
"PSC_KOEF1aEZH1_EZH1_005R",
"PSC_KOEF1aEZH1_SUZ12_005R",
"PSC_KOEF1aEZH1_H3K27me3_005R",

"PSC_WT_EZH1_006R",
"PSC_WT_EZH2_006R",
"PSC_WT_SUZ12_006R",
"PSC_WT_H3K27me3_006R",
"PSC_KO_EZH1_006R",
"PSC_KO_H3K27me3_006R",
"PSC_KOEF1aEZH1_EZH1_006R",
"PSC_KOEF1aEZH1_EZH2_006R",
"PSC_KOEF1aEZH1_SUZ12_006R",
"PSC_KOEF1aEZH1_H3K27me3_006R",

"PSC_WT_EZH2_010R",
"PSC_WT_H3K27me3_010R",
"PSC_WT_SUZ12_013R1",
"PSC_WT_H3K27me3_013R1",
"PSC_KO_EZH1_013R1",
"PSC_KO_EZH2_013R1",
"PSC_KO_SUZ12_013R1",
"PSC_KO_H3K27me3_013R1",
"PSC_KOEF1aEZH1_EZH1_013R1",
"PSC_KOEF1aEZH1_EZH2_013R1",
"PSC_KOEF1aEZH1_SUZ12_013R1",
"PSC_KOEF1aEZH1_H3K27me3_013R1",

"PSC_WT_EZH2_014R1",
"PSC_WT_SUZ12_014R1",
"PSC_KO_EZH1_014R2",
"PSC_KO_EZH2_014R1",
"PSC_KO_EZH2_014R2",
"PSC_KO_SUZ12_014R1",
"PSC_KO_SUZ12_014R2",
"PSC_KO_H3K27me3_014R2",
"PSC_KOEF1aEZH1_EZH2_014R1"
)

# Set new_sample_ID as a factor with the specified order
readCount <- readCount %>%
  mutate(new_sample_ID = factor(new_sample_ID, levels = sample_order))


# plot
pdf("output/QC/histogram_ReadCount.pdf", width=12, height=6)  # Increase the width for better spacing
ggplot(readCount, aes(x = new_sample_ID, y = Count_Value, fill = exp)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Read_Type, scales = "free_y") +  # Place facet_wrap before theme adjustments
  geom_text(data = readCount, aes(x = new_sample_ID, y = 0, label = Ecoli_prop), 
            angle = 90, vjust = 0.5, hjust = -0.5, size = 3, color = "blue") +  # Adjust hjust for alignment
  labs(x = "Sample ID (Ecoli Proportion in Blue)", y = "Read Count", fill = "Experiment") +  # Label axes
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 10)  # Adjust facet label size
  )
dev.off()

```



## Calculate E Coli SF (scaling factor)

--> I did it in `sample_001016.xlsx`. To adjust library size, I did:
  - spike in proportion = (uniqMappedMG1655Reads / uniqMappedReads) *100
  - adjustedLibSize = spike in proportion * uniqMappedReads
  


**Using our spike in proportion, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/)) 



```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# ONE PER ONE
## H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_H3K27me3.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(9441400,4470000,1140000,4401600,1244000,2761200,5406400,9551400,1025200), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_H3K27me3.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_H3K27me3_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_H3K27me3_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()




# ONE PER ONE
## EZH2

### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_EZH2.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_EZH2.RData")
load("output/DiffBind/sample_count_macs2raw_unique_EZH2.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_EZH2.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_EZH2.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(15963200,5584400,4618000,1898400,4033200,3976200,18113200,1707800,4936300), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_EZH2.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_EZH2_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_EZH2_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()


## SUZ12

### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_SUZ12.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_SUZ12.RData")
load("output/DiffBind/sample_count_macs2raw_unique_SUZ12.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_SUZ12.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_SUZ12.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(19901200,1859800,4436600,1674200,4206000,4051800,11078200,27481600,1739600), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_SUZ12.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_SUZ12_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_SUZ12_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
```



# Normalization method from Ferguson et al 

Summary pipeline:
- mapping Bowtie2 and sort bam removing unmapped fragments with SAMtools (`Phred>33 --dovetail`)
  - Collect the `output/bowtie2/*filter.bam` reads
- Convert bam to bigwig (keeping duplicate!) with 50bp bins
- Identify local maxima and generate blacklist region from the bigiwg
  - let's instead use ENCODE blacklist regions





```bash
conda activate deeptools

# Convert bam to bigwig (keeping duplicate!)
sbatch scripts/bamtobigwig_Ferguson_1.sh # 33969577 ok
sbatch scripts/bamtobigwig_Ferguson_2.sh # 33969578 ok
sbatch scripts/bamtobigwig_Ferguson_3.sh # 33969581 ok


# Convert bam to bigwig (unique! 50bp resolution)
sbatch scripts/bamtobigwig_Ferguson50bp_1.sh # 35637227 ok
sbatch scripts/bamtobigwig_Ferguson50bp_2.sh # 35637229 ok
sbatch scripts/bamtobigwig_Ferguson50bp_3.sh # 35637230 ok 

# Convert bigwig to bedgraph
conda activate BedToBigwig
## Default bigwig bin50
sbatch scripts/BedToBigwig_Ferguson.sh # 33973549 ok
## IGG subtracted bigwig bin50
sbatch scripts/BedToBigwig_Ferguson_subtractIGG.sh # 34123124 ok
## Unique bigwig (1bp resolution)
sbatch scripts/BedToBigwig_Ferguson_unique.sh # 35195181 ok
## Unique bigwig (1bp resolution) - IGG subtracted
sbatch scripts/BedToBigwig_Ferguson_subtractIGG_unique.sh # 35204218 ok
## Unique bigwig (50bp resolution)
sbatch scripts/BedToBigwig_Ferguson50bp_unique.sh # 35637290 ok




# Remove blacklist regions
conda activate BedToBigwig
## Default bigwig bin50
sbatch scripts/BedintersectBlacklist_Ferguson.sh # 33974981 ok
## IGG subtracted bigwig bin50
sbatch scripts/BedintersectBlacklist_Ferguson_subtractIGG.sh # 34124084 ok
## Unique bigwig (1bp resolution)
sbatch scripts/BedintersectBlacklist_Ferguson_unique.sh # 35195293 ok
## Unique bigwig (1bp resolution) - IGG subtracted
sbatch scripts/BedintersectBlacklist_Ferguson_subtractIGG_unique.sh # 35204304 ok
## Unique bigwig (50bp resolution)
sbatch scripts/BedintersectBlacklist_Ferguson50bp_unique.sh # 35637868 ok


```

Use Python to identify local maxima, quantify the height for the 75-99th percentile peak

```bash
srun --mem=250g --pty bash -l

# Identify local maxima
## Default bigwig bin50
python scripts/LocalMaxima_Ferguson.py
## IGG subtracted bigwig bin50
python scripts/LocalMaxima_Ferguson_subtractIGG.py
## Unique bigwig (1bp resolution)
python scripts/LocalMaxima_Ferguson_unique.py
## Unique bigwig (1bp resolution) - IGG subtracted
python scripts/LocalMaxima_Ferguson_subtractIGG_unique.py
## Unique bigwig (50bp resolution)
python scripts/LocalMaxima_Ferguson50bp_unique.py




#  calculate the 99th percentile of the signal heights (score) in the local maxima files.
## Default bigwig bin50
python scripts/Percentile99_Ferguson.py
python scripts/Percentile75_Ferguson.py
python scripts/Percentile90_Ferguson.py
python scripts/Percentile95_Ferguson.py
python scripts/Percentile98_Ferguson.py
## IGG subtracted bigwig bin50
python scripts/Percentile99_Ferguson_subtractIGG.py
python scripts/Percentile95_Ferguson_subtractIGG.py
python scripts/Percentile90_Ferguson_subtractIGG.py
## Unique bigwig (1bp resolution)
python scripts/Percentile99_Ferguson_unique.py
python scripts/Percentile95_Ferguson_unique.py
python scripts/Percentile90_Ferguson_unique.py
## Unique bigwig (1bp resolution) - IGG subtracted
python scripts/Percentile99_Ferguson_subtractIGG_unique.py
python scripts/Percentile95_Ferguson_subtractIGG_unique.py
python scripts/Percentile90_Ferguson_subtractIGG_unique.py
## Unique bigwig (50bp resolution)
python scripts/Percentile99_Ferguson50bp_unique.py




# normalize AB per AB (using WT sample 1st replicate as reference)
## Default bigwig bin50
### default 99th percentile
python scripts/norm_H3K27me3_Ferguson.py
python scripts/norm_SUZ12_Ferguson.py
python scripts/norm_EZH2_Ferguson.py
python scripts/norm_IGG_Ferguson.py
### 90th percentile
python scripts/norm_H3K27me3_Ferguson_Perc90.py
python scripts/norm_SUZ12_Ferguson_Perc90.py
python scripts/norm_EZH2_Ferguson_Perc90.py
python scripts/norm_IGG_Ferguson_Perc90.py
### 95th percentile
python scripts/norm_H3K27me3_Ferguson_Perc95.py
python scripts/norm_SUZ12_Ferguson_Perc95.py
python scripts/norm_EZH2_Ferguson_Perc95.py
python scripts/norm_IGG_Ferguson_Perc95.py
### 98th percentile
python scripts/norm_H3K27me3_Ferguson_Perc98.py
python scripts/norm_SUZ12_Ferguson_Perc98.py
python scripts/norm_EZH2_Ferguson_Perc98.py
python scripts/norm_IGG_Ferguson_Perc98.py

## IGG subtracted bigwig bin50
### 90th percentile
python scripts/norm_H3K27me3_Ferguson_Perc90_subtractIGG.py
python scripts/norm_SUZ12_Ferguson_Perc90_subtractIGG.py
python scripts/norm_EZH2_Ferguson_Perc90_subtractIGG.py
### 95th percentile
python scripts/norm_H3K27me3_Ferguson_Perc95_subtractIGG.py
python scripts/norm_SUZ12_Ferguson_Perc95_subtractIGG.py
python scripts/norm_EZH2_Ferguson_Perc95_subtractIGG.py
### 99th percentile
python scripts/norm_H3K27me3_Ferguson_Perc99_subtractIGG.py
python scripts/norm_SUZ12_Ferguson_Perc99_subtractIGG.py
python scripts/norm_EZH2_Ferguson_Perc99_subtractIGG.py

## Unique bigwig (1bp resolution)
### 90th percentile
python scripts/norm_H3K27me3_Ferguson_Perc90_unique.py
python scripts/norm_SUZ12_Ferguson_Perc90_unique.py
python scripts/norm_EZH2_Ferguson_Perc90_unique.py
### 95th percentile
python scripts/norm_H3K27me3_Ferguson_Perc95_unique.py
python scripts/norm_SUZ12_Ferguson_Perc95_unique.py
python scripts/norm_EZH2_Ferguson_Perc95_unique.py
### 99th percentile
python scripts/norm_H3K27me3_Ferguson_Perc99_unique.py
python scripts/norm_SUZ12_Ferguson_Perc99_unique.py
python scripts/norm_EZH2_Ferguson_Perc99_unique.py

## Unique bigwig (1bp resolution) - IGG subtracted
### 90th percentile
python scripts/norm_H3K27me3_Ferguson_Perc90_subtractIGG_unique.py
python scripts/norm_SUZ12_Ferguson_Perc90_subtractIGG_unique.py
python scripts/norm_EZH2_Ferguson_Perc90_subtractIGG_unique.py
### 95th percentile
python scripts/norm_H3K27me3_Ferguson_Perc95_subtractIGG_unique.py
python scripts/norm_SUZ12_Ferguson_Perc95_subtractIGG_unique.py
python scripts/norm_EZH2_Ferguson_Perc95_subtractIGG_unique.py
### 99th percentile
python scripts/norm_H3K27me3_Ferguson_Perc99_subtractIGG_unique.py
python scripts/norm_SUZ12_Ferguson_Perc99_subtractIGG_unique.py
python scripts/norm_EZH2_Ferguson_Perc99_subtractIGG_unique.py


## Unique bigwig (50bp resolution)
### 99th percentile
python scripts/norm_H3K27me3_Ferguson50bp_Perc99_unique.py
python scripts/norm_SUZ12_Ferguson50bp_Perc99_unique.py
python scripts/norm_EZH2_Ferguson50bp_Perc99_unique.py


#python scripts/norm_EZH1_Ferguson.py

# v2= test using scaling_factor = reference_value / percentile_value
python scripts/norm_H3K27me3_Ferguson_v2.py
python scripts/norm_SUZ12_Ferguson_v2.py
python scripts/norm_EZH2_Ferguson_v2.py
python scripts/norm_IGG_Ferguson_v2.py



## Unique bigwig (1bp resolution) APPLYING SF TO INITIAL BIGWIG
### 99th percentile
python scripts/norm_H3K27me3_Ferguson_Perc99_unique_initialBigwig.py
python scripts/norm_SUZ12_Ferguson_Perc99_unique_initialBigwig.py
python scripts/norm_EZH2_Ferguson_Perc99_unique_initialBigwig.py


```
--> Works!

- *NOTE: **Local Maxima** = value is higher than its neighboring points. In the context of your CUT&RUN data (or other genomic data), local maxima refer to genomic positions where the signal intensity (e.g., read depth or coverage in the bedGraph file) is greater than the signal in the surrounding regions.*
- *NOTE: **Percentile 99** = signal level that is greater than 99% of all other signal values in the dataset.*
- *NOTE: in `norm_*_v2.py` I tested `scaling_factor = percentile_value / reference_value` instead of `scaling_factor = reference_value / percentile_value` and it is not good!! **V2 is NOT GOOD***
- *NOTE: SF obtained from unique 1bp and unique 50bp are very different!*


--> **IMPORTANT NOTE: I noticed that I applied the SF to the local maxima bigwig, BUT I should have applied it to the initial bigwig files!!!**
  --> This has been corrected `APPLYING SF TO INITIAL BIGWIG`; `*_initialBigwig`




Convert normalized bedGraph back to bigwig

```bash
conda activate BedToBigwig

sbatch scripts/BedToBigwig_Norm_Ferguson.sh # 33980900 ok
sbatch scripts/BedToBigwig_Norm_IGG_Ferguson.sh # 33984013 ok

# default
sbatch scripts/BedToBigwig_Normv2_Ferguson.sh # 34088401 ok
sbatch scripts/BedToBigwig_Norm90_Ferguson.sh # 34115498 ok
sbatch scripts/BedToBigwig_Norm95_Ferguson.sh # 34091046 ok
sbatch scripts/BedToBigwig_Norm98_Ferguson.sh # 34115552 ok
# subtract IGG signal in raw file
sbatch scripts/BedToBigwig_Norm90_Ferguson_subtractIGG.sh # 34142015 ok
sbatch scripts/BedToBigwig_Norm95_Ferguson_subtractIGG.sh # 34142064 ok
sbatch scripts/BedToBigwig_Norm99_Ferguson_subtractIGG.sh # 34142137 ok
# Unique bigwig (1bp resolution)
sbatch scripts/BedToBigwig_Norm90_Ferguson_unique.sh # 35196105 ok
sbatch scripts/BedToBigwig_Norm95_Ferguson_unique.sh # 35196140 ok
sbatch scripts/BedToBigwig_Norm99_Ferguson_uniqu.sh # 35196143 ok
## Unique bigwig (1bp resolution) - IGG subtracted
sbatch scripts/BedToBigwig_Norm90_Ferguson_subtractIGG_unique.sh # 35213761 ok
sbatch scripts/BedToBigwig_Norm95_Ferguson_subtractIGG_unique.sh # 35213763 ok
sbatch scripts/BedToBigwig_Norm99_Ferguson_subtractIGG_unique.sh # 35213766 ok
# Unique bigwig (50bp resolution)
sbatch scripts/BedToBigwig_Norm99_Ferguson50bp_unique.sh # 35639537 ok


# Unique bigwig (1bp resolution) APPLYING SF TO INITIAL BIGWIG
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_1.sh # 39902700 xxx
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_2.sh # 39902701 xxx
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_3.sh # 39902702 xxx




# Subtract Igg signal (after normalization - likely not recommended)
conda activate deeptools

sbatch scripts/bigwigCompare_Norm_Ferguson_subtractIGG.sh # 33995282 ok
sbatch scripts/bigwigCompare_Norm_Ferguson_subtractIGG_unique.sh # 35197178 ok

```
--> Replicates are very heterogeneous... Subtracting processed Igg same... Test with subtracting IGG from raw files after.

--> Using 75 percentile give same SF; tried 90, 95, 98
  --> 90 perform best! Almost identical replicate!!! Let's try subtracting IGG on raw before applying normalization see if improvement

--> Unique reads at 99percentile works GREAT! Replicate homogeneous (subtract IGG perform badly)
  --> 50bp perform badly! Let's try instead to smooth the 1bp version
    --> NO need smoothing, instead use the `_initialBigwig` one --> These are all good! I mistakenly applied SF to local maxima previously...



Let's try to use sample-specific blacklist regions, for that I will use [Greenscreen](https://github.com/sklasfeld/GreenscreenProject) to generate a blacklist for all our samples:
- Call peaks in all Igg files



```bash
conda activate macs2

sbatch scripts/macs2_broad_greenscreen.sh # 34088154 ok

#sbatch scripts/generate_greenscreenBed.sh [QVAL] [MERGE_DISTANCE] [DISTINCT_NINPUTS]
sbatch scripts/generate_greenscreenBed.sh 10 5000 5 # 34088405 ok
sbatch scripts/generate_greenscreenBed.sh 10 5000 3 # 34088536 ok


```
--> We have 12 IGG, so if region call in >= 5/3 IGG it is accounted as a Greenscreen region.
  --> Greenscreen region file: `output/macs2/broad/qval10/gs_merge5000bp_call5_12inputs.bed`
    --> Only 13/17 regions removed when using >= 5/3 IGG



Let's try to remove IGG signal BEFORE normalization:


```bash
conda activate deeptools


# Subtract Igg signal from raw bedgraph file
sbatch scripts/bigwigCompare_raw_subtractIGG.sh # 34088708 ok
```

--> Overall method work great. However, the 1bp resolution I started with make the EZH2/SUZ12 bigwig weird looking; very sharp... Let's instead generate 50bp resolution bigwig (from the bam)














# DIFFREPS


Lets use DIFFREPS to identify gain lost regions

From `001*/009*` best parameters to use were:
- *gt test, pval 0.05 for window, padj 0.05 FDR*
- in window of *5kb, 2kb, 1kb, 500bp, 250bp*;
--> and then combine into a single file. More detail at `## Run diffreps` 



## WT vs KO - DIFFREPS - FAIL

--> FAIL here as I applied SF to local maxima and not initial bigwigs...


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bedGraph 
output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bedGraph
output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bedGraph 

output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bedGraph 
output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bedGraph
output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bedGraph 

# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05


# 2000bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05



# 1000bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05


# 250bp every 50bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05

```


### Explore diffreps results in R




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin5000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin2000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 







# Replace Inf by min/max values
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == Inf] <- max(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == -Inf] <- min(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == Inf] <- max(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == -Inf] <- min(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == Inf] <- max(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == -Inf] <- min(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == Inf] <- max(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)
bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == -Inf] <- min(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)

bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == Inf] <- max(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)
bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == -Inf] <- min(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)



# List of dataset names
file_names <- c("bin5000space100_gt_pval05", "bin2000space100_gt_pval05", "bin1000space100_gt_pval05", "bin500space100_gt_pval05", "bin250space50_gt_pval05")

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

pdf("output/diffreps/hist-log2FC_distribution-padj05.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.05) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 4, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()



# Combine windows - pval05
combined_data_select = combined_data %>% 
  filter(padj<0.05)

## Convert to GRanges
## Convert combined_data_filt to GRanges
gr_combined <- GRanges(
  seqnames = combined_data_select$Chrom,
  ranges = IRanges(start = combined_data_select$Start, end = combined_data_select$End),
  log2FC = combined_data_select$log2FC,
  padj = combined_data_select$padj,
  dataset = combined_data_select$dataset
)
## Merge overlapping windows across all datasets
merged_gr <- reduce(gr_combined, ignore.strand = TRUE)
## Find overlaps with original intervals
ov <- findOverlaps(merged_gr, gr_combined)
## Summarize merged regions and assign labels
merged_df <- as.data.frame(merged_gr) %>%
  mutate(
    log2FC_list = lapply(seq_along(merged_gr), function(i) gr_combined$log2FC[subjectHits(ov)[queryHits(ov) == i]]),
    dataset_list = lapply(seq_along(merged_gr), function(i) gr_combined$dataset[subjectHits(ov)[queryHits(ov) == i]]),
    direction = sapply(log2FC_list, function(fc) {
      if (all(fc > 0)) return("Gain")
      if (all(fc < 0)) return("Lost")
      return("Mixed")
    }),
    Largest_window = sapply(dataset_list, function(ds) {
      if ("bin5000space100_gt_pval05" %in% ds) return("5kb")
      if ("bin2000space100_gt_pval05" %in% ds) return("2kb")
      if ("bin1000space100_gt_pval05" %in% ds) return("1kb")
      if ("bin500space100_gt_pval05" %in% ds) return("1kb")
      return("250bp")
    }),
    log2FC = sapply(seq_along(log2FC_list), function(i) {
      ds <- dataset_list[[i]]
      fc <- log2FC_list[[i]]
      
      # Mixed: both negative and positive log2FC
      if(any(fc > 0) && any(fc < 0)) {
        return(paste(min(fc), max(fc), sep = "_"))
      }
      
      # Non-mixed: safely find log2FC of Largest_window
      idx <- which(ds == Largest_window[i])
      if(length(idx) > 0) return(fc[idx[1]])
      
      # Fallback if for some reason largest window is missing (rare)
      return(round(mean(fc), 2))
    })
  ) %>%
  select(seqnames, start, end, direction, Largest_window, log2FC) %>%
  as_tibble()
#--> 395 Gain, 331 Lost, 3 Mixed

# PLOT combine windows
merged_df_counts <- merged_df %>%  
  filter(direction != "Mixed") %>%
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(direction) %>%
  summarise(count = n(), .groups = "drop")

pdf("output/diffreps/hist-log2FC_distribution-padj05_gt_pval05-WindowCombine_5kb2kb1kb500bp250bp.pdf", width=3, height=3)
merged_df %>%  
  filter(direction != "Mixed") %>%
  mutate(log2FC = as.numeric(log2FC)) %>%
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  labs(title = "Log2FC Distribution",
       x = "Log2 Fold Change",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 7, face = "bold")) +
  geom_text(data = merged_df_counts , 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(merged_df_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()

## Save output
write.table(merged_df, "output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO.txt", sep = "\t", quote = FALSE, row.names = FALSE)


```

--> Looks great!! Nice combination of gain/lost, with a bit more gain.








## WT vs KO - DIFFREPS - initialBigwig - G test pval 05

--> GOOD here as I applied SF to initial bigwigs


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bedGraph 

# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05


# 2000bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05



# 1000bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05


# 250bp every 50bp - G test with pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05

```




### Explore diffreps results in R




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 







# Replace Inf by min/max values
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == Inf] <- max(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == -Inf] <- min(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == Inf] <- max(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == -Inf] <- min(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == Inf] <- max(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == -Inf] <- min(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == Inf] <- max(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)
bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == -Inf] <- min(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)

bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == Inf] <- max(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)
bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == -Inf] <- min(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)



# List of dataset names
file_names <- c("bin5000space100_gt_pval05", "bin2000space100_gt_pval05", "bin1000space100_gt_pval05", "bin500space100_gt_pval05", "bin250space50_gt_pval05")

## Function to read and format each file
read_and_process <- function(file) {
  df <- get(file)  # Load dataset from environment
  df$dataset <- file  # Add dataset identifier
  return(df)
}

## Combine all datasets into one
combined_data <- bind_rows(lapply(file_names, read_and_process)) 

combined_data_counts <- combined_data %>% 
  filter(padj<0.01) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")

## plot

pdf("output/diffreps/hist-log2FC_distribution-padj01_initialBigwig.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.01) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 4, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()



# Combine windows - pval05
combined_data_select = combined_data %>% 
  filter(padj<0.001)

## Convert to GRanges
## Convert combined_data_filt to GRanges
gr_combined <- GRanges(
  seqnames = combined_data_select$Chrom,
  ranges = IRanges(start = combined_data_select$Start, end = combined_data_select$End),
  log2FC = combined_data_select$log2FC,
  padj = combined_data_select$padj,
  dataset = combined_data_select$dataset
)
## Merge overlapping windows across all datasets
merged_gr <- reduce(gr_combined, ignore.strand = TRUE)
## Find overlaps with original intervals
ov <- findOverlaps(merged_gr, gr_combined)
## Summarize merged regions and assign labels
merged_df <- as.data.frame(merged_gr) %>%
  mutate(
    log2FC_list = lapply(seq_along(merged_gr), function(i) gr_combined$log2FC[subjectHits(ov)[queryHits(ov) == i]]),
    dataset_list = lapply(seq_along(merged_gr), function(i) gr_combined$dataset[subjectHits(ov)[queryHits(ov) == i]]),
    direction = sapply(log2FC_list, function(fc) {
      if (all(fc > 0)) return("Gain")
      if (all(fc < 0)) return("Lost")
      return("Mixed")
    }),
    Largest_window = sapply(dataset_list, function(ds) {
      if ("bin5000space100_gt_pval05" %in% ds) return("5kb")
      if ("bin2000space100_gt_pval05" %in% ds) return("2kb")
      if ("bin1000space100_gt_pval05" %in% ds) return("1kb")
      if ("bin500space100_gt_pval05" %in% ds) return("1kb")
      return("250bp")
    }),
    log2FC = sapply(seq_along(log2FC_list), function(i) {
      ds <- dataset_list[[i]]
      fc <- log2FC_list[[i]]
      
      # Mixed: both negative and positive log2FC
      if(any(fc > 0) && any(fc < 0)) {
        return(paste(min(fc), max(fc), sep = "_"))
      }
      
      # Non-mixed: safely find log2FC of Largest_window
      idx <- which(ds == Largest_window[i])
      if(length(idx) > 0) return(fc[idx[1]])
      
      # Fallback if for some reason largest window is missing (rare)
      return(round(mean(fc), 2))
    })
  ) %>%
  select(seqnames, start, end, direction, Largest_window, log2FC) %>%
  as_tibble()


# PLOT combine windows
## FC pos/neg 0 treshold
merged_df_counts <- merged_df %>%  
  mutate(log2FC = as.numeric(log2FC)) %>%
  filter(direction != "Mixed") %>%
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(direction) %>%
  summarise(count = n(), .groups = "drop")
## FC pos/neg 1 treshold
merged_df_counts <- merged_df %>%  
  mutate(log2FC = as.numeric(log2FC)) %>%
  filter(log2FC > 1 | log2FC < -1) %>%
  mutate(direction = ifelse(log2FC < -1, "Negative", "Positive")) %>%
  group_by(direction) %>%
  summarise(count = n(), .groups = "drop")

pdf("output/diffreps/hist-log2FC_distribution-padj001_gt_pval05-WindowCombine_5kb2kb1kb500bp250bp_initialBigwig.pdf", width=3, height=3)
merged_df %>%  
  filter(direction != "Mixed") %>%
  mutate(log2FC = as.numeric(log2FC)) %>%
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  labs(title = "Log2FC Distribution",
       x = "Log2 Fold Change",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 7, face = "bold")) +
  geom_text(data = merged_df_counts , 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(merged_df_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()

## Save output
write.table(merged_df, "output/diffreps/merged_intervals-padj001_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)


```

--> Looks great!! Nice combination of gain/lost, with a bit more gain. AND much more Diff. bound site here!!








## WT vs KO - DIFFREPS - initialBigwig - NB test Default

--> GOOD here as I applied SF to initial bigwigs


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bedGraph 

# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp -  Negative binomial pval 0.0001
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_nb_pval0001-diff.nb.txt --window 5000 --step 100 --meth nb --pval 0.0001


# 2000bp every 100bp -  Negative binomial pval 0.0001
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_nb_pval0001-diff.nb.txt --window 2000 --step 100 --meth nb --pval 0.0001



# 1000bp every 100bp -  Negative binomial pval 0.0001
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_nb_pval0001-diff.nb.txt --window 1000 --step 100 --meth nb --pval 0.0001

# 500bp every 100bp -  Negative binomial pval 0.0001
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_nb_pval0001-diff.nb.txt --window 500 --step 100 --meth nb --pval 0.0001


# 250bp every 50bp -  Negative binomial pval 0.0001
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_nb_pval0001-diff.nb.txt --window 250 --step 50 --meth nb --pval 0.0001

```




### Explore diffreps results in R




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_nb_pval0001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_nb_pval0001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_nb_pval0001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_nb_pval0001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_nb_pval0001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_nb_pval0001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_nb_pval0001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_nb_pval0001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_nb_pval0001 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_nb_pval0001-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 



# Replace Inf by min/max values
bin5000space100_nb_pval0001$log2FC[bin5000space100_nb_pval0001$log2FC == Inf] <- max(bin5000space100_nb_pval0001$log2FC[is.finite(bin5000space100_nb_pval0001$log2FC)], na.rm = TRUE)
bin5000space100_nb_pval0001$log2FC[bin5000space100_nb_pval0001$log2FC == -Inf] <- min(bin5000space100_nb_pval0001$log2FC[is.finite(bin5000space100_nb_pval0001$log2FC)], na.rm = TRUE)

bin2000space100_nb_pval0001$log2FC[bin2000space100_nb_pval0001$log2FC == Inf] <- max(bin2000space100_nb_pval0001$log2FC[is.finite(bin2000space100_nb_pval0001$log2FC)], na.rm = TRUE)
bin2000space100_nb_pval0001$log2FC[bin2000space100_nb_pval0001$log2FC == -Inf] <- min(bin2000space100_nb_pval0001$log2FC[is.finite(bin2000space100_nb_pval0001$log2FC)], na.rm = TRUE)

bin1000space100_nb_pval0001$log2FC[bin1000space100_nb_pval0001$log2FC == Inf] <- max(bin1000space100_nb_pval0001$log2FC[is.finite(bin1000space100_nb_pval0001$log2FC)], na.rm = TRUE)
bin1000space100_nb_pval0001$log2FC[bin1000space100_nb_pval0001$log2FC == -Inf] <- min(bin1000space100_nb_pval0001$log2FC[is.finite(bin1000space100_nb_pval0001$log2FC)], na.rm = TRUE)

bin500space100_nb_pval0001$log2FC[bin500space100_nb_pval0001$log2FC == Inf] <- max(bin500space100_nb_pval0001$log2FC[is.finite(bin500space100_nb_pval0001$log2FC)], na.rm = TRUE)
bin500space100_nb_pval0001$log2FC[bin500space100_nb_pval0001$log2FC == -Inf] <- min(bin500space100_nb_pval0001$log2FC[is.finite(bin500space100_nb_pval0001$log2FC)], na.rm = TRUE)

bin250space50_nb_pval0001$log2FC[bin250space50_nb_pval0001$log2FC == Inf] <- max(bin250space50_nb_pval0001$log2FC[is.finite(bin250space50_nb_pval0001$log2FC)], na.rm = TRUE)
bin250space50_nb_pval0001$log2FC[bin250space50_nb_pval0001$log2FC == -Inf] <- min(bin250space50_nb_pval0001$log2FC[is.finite(bin250space50_nb_pval0001$log2FC)], na.rm = TRUE)



# List of dataset names
file_names <- c("bin5000space100_nb_pval0001", "bin2000space100_nb_pval0001", "bin1000space100_nb_pval0001", "bin500space100_nb_pval0001", "bin250space50_nb_pval0001")

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

pdf("output/diffreps/hist-log2FC_distribution-padj05_nb_pval0001_initialBigwig.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.05) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 4, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()



# Combine windows 
combined_data_select = combined_data %>% 
  filter(padj<0.001)

## Convert to GRanges
## Convert combined_data_filt to GRanges
gr_combined <- GRanges(
  seqnames = combined_data_select$Chrom,
  ranges = IRanges(start = combined_data_select$Start, end = combined_data_select$End),
  log2FC = combined_data_select$log2FC,
  padj = combined_data_select$padj,
  dataset = combined_data_select$dataset
)
## Merge overlapping windows across all datasets
merged_gr <- reduce(gr_combined, ignore.strand = TRUE)
## Find overlaps with original intervals
ov <- findOverlaps(merged_gr, gr_combined)
## Summarize merged regions and assign labels
merged_df <- as.data.frame(merged_gr) %>%
  mutate(
    log2FC_list = lapply(seq_along(merged_gr), function(i) gr_combined$log2FC[subjectHits(ov)[queryHits(ov) == i]]),
    dataset_list = lapply(seq_along(merged_gr), function(i) gr_combined$dataset[subjectHits(ov)[queryHits(ov) == i]]),
    direction = sapply(log2FC_list, function(fc) {
      if (all(fc > 0)) return("Gain")
      if (all(fc < 0)) return("Lost")
      return("Mixed")
    }),
    Largest_window = sapply(dataset_list, function(ds) {
      if ("bin5000space100_nb_pval0001" %in% ds) return("5kb")
      if ("bin2000space100_nb_pval0001" %in% ds) return("2kb")
      if ("bin1000space100_nb_pval0001" %in% ds) return("1kb")
      if ("bin500space100_nb_pval0001" %in% ds) return("1kb")
      return("250bp")
    }),
    log2FC = sapply(seq_along(log2FC_list), function(i) {
      ds <- dataset_list[[i]]
      fc <- log2FC_list[[i]]
      
      # Mixed: both negative and positive log2FC
      if(any(fc > 0) && any(fc < 0)) {
        return(paste(min(fc), max(fc), sep = "_"))
      }
      
      # Non-mixed: safely find log2FC of Largest_window
      idx <- which(ds == Largest_window[i])
      if(length(idx) > 0) return(fc[idx[1]])
      
      # Fallback if for some reason largest window is missing (rare)
      return(round(mean(fc), 2))
    })
  ) %>%
  select(seqnames, start, end, direction, Largest_window, log2FC) %>%
  as_tibble()
#--> 395 Gain, 331 Lost, 3 Mixed

# PLOT combine windows
## FC pos/neg 0 treshold
merged_df_counts <- merged_df %>%  
  mutate(log2FC = as.numeric(log2FC)) %>%
  filter(direction != "Mixed") %>%
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(direction) %>%
  summarise(count = n(), .groups = "drop")
## FC pos/neg 1 treshold
merged_df_counts <- merged_df %>%  
  mutate(log2FC = as.numeric(log2FC)) %>%
  filter(log2FC > 1 | log2FC < -1) %>%
  mutate(direction = ifelse(log2FC < -1, "Negative", "Positive")) %>%
  group_by(direction) %>%
  summarise(count = n(), .groups = "drop")

pdf("output/diffreps/hist-log2FC_distribution-padj001_nb_pval0001_log2FC1-WindowCombine_5kb2kb1kb500bp250bp_initialBigwig.pdf", width=3, height=3)
merged_df %>%  
  filter(direction != "Mixed") %>%
  mutate(log2FC = as.numeric(log2FC)) %>%
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
  labs(title = "Log2FC Distribution",
       x = "Log2 Fold Change",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 7, face = "bold")) +
  geom_text(data = merged_df_counts , 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(merged_df_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()

## Save output
write.table(merged_df, "output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)


```

--> Looks great!! Nice combination of gain/lost, with a bit more gain. AND much more Diff. bound site here!!
  --> `padj001_nb_pval0001_log2FC1` may be optimal








## WT vs KO - DIFFREPS - initialBigwig - G test pval 0.05

--> Focus on `pval 0.05 padj 0.001 1000bp 100bp`: parameter that works best in `001*/009*`.


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bedGraph 

# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05


# 2000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05



# 1000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05


# 250bp every 50bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05

```





### Explore diffreps results in R




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/PSC_WT_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 



# Replace Inf by min/max values
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == Inf] <- max(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == -Inf] <- min(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == Inf] <- max(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == -Inf] <- min(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == Inf] <- max(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == -Inf] <- min(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == Inf] <- max(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)
bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == -Inf] <- min(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)

bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == Inf] <- max(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)
bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == -Inf] <- min(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)



# List of dataset names
file_names <- c("bin5000space100_gt_pval05", "bin2000space100_gt_pval05", "bin1000space100_gt_pval05", "bin500space100_gt_pval05", "bin250space50_gt_pval05")

## Function to read and format each file
read_and_process <- function(file) {
  df <- get(file)  # Load dataset from environment
  df$dataset <- file  # Add dataset identifier
  return(df)
}

## Combine all datasets into one
combined_data <- bind_rows(lapply(file_names, read_and_process)) 

combined_data_counts <- combined_data %>% 
  filter(padj<0.001) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")

## plot

pdf("output/diffreps/hist-log2FC_distribution-padj001_gt_pval05_initialBigwig.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.05) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 4, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()



# Combine windows 
combined_data_select = combined_data %>% 
  filter(padj<0.001)

## Convert to GRanges
## Convert combined_data_filt to GRanges
gr_combined <- GRanges(
  seqnames = combined_data_select$Chrom,
  ranges = IRanges(start = combined_data_select$Start, end = combined_data_select$End),
  log2FC = combined_data_select$log2FC,
  padj = combined_data_select$padj,
  dataset = combined_data_select$dataset
)
## Merge overlapping windows across all datasets
merged_gr <- reduce(gr_combined, ignore.strand = TRUE)
## Find overlaps with original intervals
ov <- findOverlaps(merged_gr, gr_combined)
## Summarize merged regions and assign labels
merged_df <- as.data.frame(merged_gr) %>%
  mutate(
    log2FC_list = lapply(seq_along(merged_gr), function(i) gr_combined$log2FC[subjectHits(ov)[queryHits(ov) == i]]),
    dataset_list = lapply(seq_along(merged_gr), function(i) gr_combined$dataset[subjectHits(ov)[queryHits(ov) == i]]),
    direction = sapply(log2FC_list, function(fc) {
      if (all(fc > 0)) return("Gain")
      if (all(fc < 0)) return("Lost")
      return("Mixed")
    }),
    Largest_window = sapply(dataset_list, function(ds) {
      if ("bin5000space100_gt_pval05" %in% ds) return("5kb")
      if ("bin2000space100_gt_pval05" %in% ds) return("2kb")
      if ("bin1000space100_gt_pval05" %in% ds) return("1kb")
      if ("bin500space100_gt_pval05" %in% ds) return("1kb")
      return("250bp")
    }),
    log2FC = sapply(seq_along(log2FC_list), function(i) {
      ds <- dataset_list[[i]]
      fc <- log2FC_list[[i]]
      
      # Mixed: both negative and positive log2FC
      if(any(fc > 0) && any(fc < 0)) {
        return(paste(min(fc), max(fc), sep = "_"))
      }
      
      # Non-mixed: safely find log2FC of Largest_window
      idx <- which(ds == Largest_window[i])
      if(length(idx) > 0) return(fc[idx[1]])
      
      # Fallback if for some reason largest window is missing (rare)
      return(round(mean(fc), 2))
    })
  ) %>%
  select(seqnames, start, end, direction, Largest_window, log2FC) %>%
  as_tibble()
#--> 395 Gain, 331 Lost, 3 Mixed

# PLOT combine windows
## FC pos/neg 0 treshold
merged_df_counts <- merged_df %>%  
  mutate(log2FC = as.numeric(log2FC)) %>%
  filter(direction != "Mixed") %>%
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(direction) %>%
  summarise(count = n(), .groups = "drop")
## FC pos/neg 1 treshold
merged_df_counts <- merged_df %>%  
  mutate(log2FC = as.numeric(log2FC)) %>%
  filter(log2FC > 1 | log2FC < -1) %>%
  mutate(direction = ifelse(log2FC < -1, "Negative", "Positive")) %>%
  group_by(direction) %>%
  summarise(count = n(), .groups = "drop")

pdf("output/diffreps/hist-log2FC_distribution-padj001_gt_pval05_log2FC1-WindowCombine_5kb2kb1kb500bp250bp_initialBigwig.pdf", width=3, height=3)
merged_df %>%  
  filter(direction != "Mixed") %>%
  mutate(log2FC = as.numeric(log2FC)) %>%
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
  labs(title = "Log2FC Distribution",
       x = "Log2 Fold Change",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 7, face = "bold")) +
  geom_text(data = merged_df_counts , 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(merged_df_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()

## Save output
write.table(merged_df, "output/diffreps/merged_intervals-padj001_gt_pval05_log2FC1-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)


```

XXX 
--> Looks great!! Nice combination of gain/lost, with a bit more gain. AND much more Diff. bound site here!!
  --> `padj001_gt_pval05_log2FC1` may be optimal
XXX







## WT vs KOEF1aEZH1 - DIFFREPS - initialBigwig - G test pval 0.05

--> Focus on `pval 0.05 padj 0.001 1000bp 100bp`: parameter that works best in `001*/009*`.


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph
output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph 

# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 1000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/PSC_WTKOEF1aEZH1_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05



```





### Explore diffreps results in R




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin1000space100_gt_pval05 <- read.delim("output/diffreps/PSC_WTKOEF1aEZH1_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 


# Replace Inf by min/max values
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == Inf] <- max(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == -Inf] <- min(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)


# List of dataset names
file_names <- c("bin1000space100_gt_pval05")

## Function to read and format each file
read_and_process <- function(file) {
  df <- get(file)  # Load dataset from environment
  df$dataset <- file  # Add dataset identifier
  return(df)
}

## Combine all datasets into one
combined_data <- bind_rows(lapply(file_names, read_and_process)) 

combined_data_counts <- combined_data %>% 
  filter(padj<0.001) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
  mutate(direction = ifelse(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")

## plot

pdf("output/diffreps/hist-PSC_WTvsKOEF1aEZH1_H3K27me3-log2FC_distribution-padj001_gt_pval05_initialBigwig.pdf", width=4, height=4)
combined_data %>% 
  filter(padj<0.001) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 4, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()

```


--> KOEF1aEZH1 show more H3K27me3 lost. It is great, as we have OE EZH1 less EZH2, and less H3K27me3.














# SICER2



## WT vs KO - SICER2 - FAIL

--> FAIL, as I applied SF to local maxima bed and not bigwig


--> SICER2 work with simplicate. BUT I do not have always corresponding WT/KO from the same experiment!!!
  --> Here, I will do **006R WT with 006R KO, 013R1 WT with 013R1 KO**, I will *use only two bio reps*



```bash
conda activate sicer2

##################################
# DATA PREP #################
##################################

output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed
output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed
output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed

output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed
output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed
output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed

# Scale up and Round score
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bed > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bed > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bed > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bed > output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bed > output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bed > output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_scaleUpRounded.bed


##################################
# Run SICER2 #################
##################################

## FDR 0.05 Rep 005 window 200 gap 600 e-value 50000
sicer_df -t output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_scaleUpRounded.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_scaleUpRounded.bed -s hg38 --window_size 200 -fdr_df 0.05 --gap_size 600 --e_value 50000 -o output/sicer2/window200gap600fdr05evalue50000
sicer_df -t output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_scaleUpRounded.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_scaleUpRounded.bed -s hg38 --window_size 200 -fdr_df 0.05 --gap_size 600 --e_value 50000 -o output/sicer2/window200gap600fdr05evalue50000

## FDR 0.05 Rep 005 window 200 gap 600 e-value 1000
sicer_df -t output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_scaleUpRounded.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_scaleUpRounded.bed -s hg38 --window_size 200 -fdr_df 0.05 --gap_size 600 --e_value 1000 -o output/sicer2/window200gap600fdr05evalue1000
sicer_df -t output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_scaleUpRounded.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_scaleUpRounded.bed -s hg38 --window_size 200 -fdr_df 0.05 --gap_size 600 --e_value 1000 -o output/sicer2/window200gap600fdr05evalue1000


```

- `window200gap600fdr05evalue50000`: 006R show almost only Lost H3K27me3, but 013R1 is more homogeneous (even tho more Lost)
- `window200gap600fdr05evalue1000`: 006R show almost only Lost H3K27me3, but 013R1 is more homogeneous (even tho more Lost)

--> I would recommend using DIFFREPS, and not SICER2, for this experiment; likely because replicate are treated individually, I think DIFFREPS decrease batch effects by treating replicates together.











## WT vs KO - SICER2 - initialBigwig



--> SICER2 work with simplicate. BUT I do not have always corresponding WT/KO from the same experiment!!!
  --> Here, I will do **006R WT with 006R KO, 013R1 WT with 013R1 KO**, I will *use only two bio reps*



```bash
conda activate sicer2

##################################
# DATA PREP #################
##################################

output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed

output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed

# Scale up and Round score
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig.bed > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig.bed > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_initialBigwig_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig.bed > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig.bed > output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig.bed > output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig_scaleUpRounded.bed
awk 'OFS="\t" {print $1, $2, $3, $4, int($5*1000 + 0.5), "+"}' output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig.bed > output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_initialBigwig_scaleUpRounded.bed


##################################
# Run SICER2 #################
##################################

## Default parameters Rep 005
sicer_df -t output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig_scaleUpRounded.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_initialBigwig_scaleUpRounded.bed -s hg38 --window_size 200 -fdr_df 0.01 --gap_size 600 --e_value 1000 -o output/sicer2/window200gap600fdr01evalue1000_initialBigwig
#--> 2132 increase, 14394 decreased
sicer_df -t output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig_scaleUpRounded.bed output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_initialBigwig_scaleUpRounded.bed -s hg38 --window_size 200 -fdr_df 0.01 --gap_size 600 --e_value 1000 -o output/sicer2/window200gap600fdr01evalue1000_initialBigwig
#--> 5138 increase, 6315 decreased
```

- `window200gap600fdr01evalue1000_initialBigwig`: 006R show almost only Lost H3K27me3, but 013R1 is more homogeneous (even tho more Lost)




### Explore result in R

Lets follow this workflow to integrate bio rep:
- Load each replicate into R.
- Keep only significant regions (rows) from each replicate --> Use `*increased*`/`*decreased*` files
- Convert significant peaks from each replicate into GenomicRanges objects.
- Intersect the peaks across replicates
- Export final intersect peaks as a BED file




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")



###############################################################
####### FDR 0.01 window 200 gap 600 E-value 1000 ######
###############################################################
# import file
PSC_KO_H3K27me3_006R__decreased <- read.table("output/sicer2/window200gap600fdr01evalue1000_initialBigwig/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig_scaleUpRounded-W200-G600-decreased-islands-summary-FDR0.01", header=FALSE, sep="\t") %>%
  dplyr::rename(chr= V1, start= V2, end= V3, Readcount_KO= V4, Normalized_Readcount_KO= V5, Readcount_WT= V6, Normalized_Readcount_WT= V7) %>%
  as_tibble()
PSC_KO_H3K27me3_006R__increased <- read.table("output/sicer2/window200gap600fdr01evalue1000_initialBigwig/PSC_KO_H3K27me3_006R_unique_norm99_initialBigwig_scaleUpRounded-W200-G600-increased-islands-summary-FDR0.01", header=FALSE, sep="\t") %>%
  dplyr::rename(chr= V1, start= V2, end= V3, Readcount_KO= V4, Normalized_Readcount_KO= V5, Readcount_WT= V6, Normalized_Readcount_WT= V7) %>%
  as_tibble()

PSC_KO_H3K27me3_013R1__decreased <- read.table("output/sicer2/window200gap600fdr01evalue1000_initialBigwig/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig_scaleUpRounded-W200-G600-decreased-islands-summary-FDR0.01", header=FALSE, sep="\t") %>%
  dplyr::rename(chr= V1, start= V2, end= V3, Readcount_KO= V4, Normalized_Readcount_KO= V5, Readcount_WT= V6, Normalized_Readcount_WT= V7) %>%
  as_tibble()
PSC_KO_H3K27me3_013R1__increased <- read.table("output/sicer2/window200gap600fdr01evalue1000_initialBigwig/PSC_KO_H3K27me3_013R1_unique_norm99_initialBigwig_scaleUpRounded-W200-G600-increased-islands-summary-FDR0.01", header=FALSE, sep="\t") %>%
  dplyr::rename(chr= V1, start= V2, end= V3, Readcount_KO= V4, Normalized_Readcount_KO= V5, Readcount_WT= V6, Normalized_Readcount_WT= V7) %>%
  as_tibble()

# DECREASED / LOST #####################
PSC_KO_H3K27me3_006R__decreased_GR <- GRanges(seqnames=PSC_KO_H3K27me3_006R__decreased$"chr",
                   ranges=IRanges(start=PSC_KO_H3K27me3_006R__decreased$start, end=PSC_KO_H3K27me3_006R__decreased$end))
PSC_KO_H3K27me3_013R1__decreased_GR <- GRanges(seqnames=PSC_KO_H3K27me3_013R1__decreased$"chr",
                   ranges=IRanges(start=PSC_KO_H3K27me3_013R1__decreased$start, end=PSC_KO_H3K27me3_013R1__decreased$end))
# find overlap
overlap <- findOverlaps(PSC_KO_H3K27me3_006R__decreased_GR, PSC_KO_H3K27me3_013R1__decreased_GR)
rep1_overlap <- PSC_KO_H3K27me3_006R__decreased_GR[queryHits(overlap)]
rep2_overlap <- PSC_KO_H3K27me3_013R1__decreased_GR[subjectHits(overlap)]
# Merge overlapping regions into consensus peaks
consensus_peaks <- reduce(c(rep1_overlap, rep2_overlap))
# Export bed file
rtracklayer::export.bed(consensus_peaks, "output/sicer2/window200gap600fdr01evalue1000_initialBigwig/PSC_KO_H3K27me3_006R013R1-decreased.bed")

# INCREASED / GAIN #####################
PSC_KO_H3K27me3_006R__increased_GR <- GRanges(seqnames=PSC_KO_H3K27me3_006R__increased$"chr",
                   ranges=IRanges(start=PSC_KO_H3K27me3_006R__increased$start, end=PSC_KO_H3K27me3_006R__increased$end))
PSC_KO_H3K27me3_013R1__increased_GR <- GRanges(seqnames=PSC_KO_H3K27me3_013R1__increased$"chr",
                   ranges=IRanges(start=PSC_KO_H3K27me3_013R1__increased$start, end=PSC_KO_H3K27me3_013R1__increased$end))
# find overlap
overlap <- findOverlaps(PSC_KO_H3K27me3_006R__increased_GR, PSC_KO_H3K27me3_013R1__increased_GR)
rep1_overlap <- PSC_KO_H3K27me3_006R__increased_GR[queryHits(overlap)]
rep2_overlap <- PSC_KO_H3K27me3_013R1__increased_GR[subjectHits(overlap)]
# Merge overlapping regions into consensus peaks
consensus_peaks <- reduce(c(rep1_overlap, rep2_overlap))
# Export bed file
rtracklayer::export.bed(consensus_peaks, "output/sicer2/window200gap600fdr01evalue1000_initialBigwig/PSC_KO_H3K27me3_006R013R1-increased.bed")



```

- *Window 200, Gap 600, E-value 1000* --> This one seems the best + recommended by SICER2
  - Lost= 2689
  - Gain= 130
--> WEIRD!!
  




# Functional analysis with enrichR

Functional analysis **enrichR with DIFFREPS diff bound genes** --> Code copied from `001*/009*` `# Functional analysis with enrichR`


**IMPOPRTANT NOTE: Run the reading and processing ONE BY ONE !!! Otherwise, lead to bug!!!!**

```R
# library
library("tidyverse")
library("enrichR")
library("ggrepel")

# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 

### GeneSymbol list of signif gain/lost H3K27me3 from DIFFREPS intitialBigwig bin1000space100 gt pval05 padj001
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt




# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
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
pdf("output/GO/enrichR_GO_Biological_Process_2023_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.pdf", width=8, height=3)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K27me3",   # H3K27me3  H3K4me3 DEGs and H3K27me3
                    labels = c("Lost", "Gain"), # down-reg and Gain up-reg and Lost
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
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Molecular_Function_2023") # 
### GeneSymbol list of signif gain/lost H3K27me3 from DIFFREPS intitialBigwig bin1000space100 gt pval05 padj001
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt



# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
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


pdf("output/GO/enrichR_GO_Molecular_Function_2023_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.pdf", width=8, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K27me3",   # H3K27me3  H3K4me3
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
write.table(gos, "output/GO/enrichR_GO_Molecular_Function_2023_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Cellular_Component_2023") # 

### GeneSymbol list of signif gain/lost H3K27me3 from DIFFREPS intitialBigwig bin1000space100 gt pval05 padj001
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
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


pdf("output/GO/enrichR_GO_Cellular_Component_2023_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.pdf", width=8, height=4)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="H3K27me3",   # H3K27me3  H3K4me3
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
write.table(gos, "output/GO/enrichR_GO_Cellular_Component_2023_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("KEGG_2021_Human") # 

### GeneSymbol list of signif gain/lost H3K27me3 from DIFFREPS intitialBigwig bin1000space100 gt pval05 padj001
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt



# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
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


pdf("output/GO/enrichR_KEGG_2021_Human_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.pdf", width=8, height=4)

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
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X") # 

### GeneSymbol list of signif gain/lost H3K27me3 from DIFFREPS intitialBigwig bin1000space100 gt pval05 padj001
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt

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
gene_names_down <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Lost_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/annotation_PSC_WTKO_H3K27me3_PSC_WT_H3K27me3_bin1000space100_gt_pval05_padj001__Gain_annot_promoterAnd5_geneSymbol.txt", header=FALSE, stringsAsFactors=FALSE)
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

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)



# Plotting with enhanced aesthetics


pdf("output/GO/enrichR_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.pdf", width=8, height=3)

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
write.table(gos, "output/GO/enrichR_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_DIFFREPS_bin1000space100_gt_pval05_padj001_GainLost_promoterAnd5.txt", sep="\t", row.names=FALSE, quote=FALSE)


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

pdf("output/GO/VolcanoPlotTF_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_DIFFREPS_bin1000space100_gt_pval05_padj001_Gain_promoterAnd5.pdf", width=5, height=4)
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


pdf("output/GO/VolcanoPlotTF_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_DIFFREPS_bin1000space100_gt_pval05_padj001_Gain_promoterAnd5_posterCHOP1.pdf", width=5, height=4)
ggplot(gos_TF_tidy, aes(x = Odds.Ratio, y = -logAdjP, color = db)) +
  geom_point(aes(color = ifelse(-logAdjP < 1.3, "not signif.", db))) +
  scale_color_manual(values = c("blue", "lightblue", "grey")) + # Replace with your actual colors
  theme_bw() +
  labs(x = "Odds Ratio", y = "-log10(adjusted p-value)") +
  geom_text_repel(data = subset(gos_TF_tidy, logAdjP < 1.3),
                  aes(label = TF),
                  nudge_x = 0.2,  # Adjust this value to nudge labels to the right
                  size = 6,
                  max.overlaps = 50)  + # 30
  guides(color = guide_legend(override.aes = list(label = ""))) # just to remove the "a" added in fig legend
dev.off()



pdf("output/GO/VolcanoPlotTF_ENCODE_and_ChEA_Consensus_TFs_from_ChIP_DIFFREPS_bin1000space100_gt_pval05_padj001_Gain_promoterAnd5_posterCHOP1.pdf", width=5, height=4)
ggplot(gos_TF_tidy, aes(x = Odds.Ratio, y = -logAdjP, color = db)) +
  geom_point(aes(color = ifelse(-logAdjP < 1.3, "not signif.", db))) +
  scale_color_manual(values = c("blue", "lightblue", "grey")) + # Replace with your actual colors
  theme_bw() +
  labs(x = "Odds Ratio", y = "-log10(adjusted p-value)") +
  geom_text_repel(data = subset(gos_TF_tidy, logAdjP < 1.3),
                  aes(label = TF),
                  nudge_x = 0.2,  # Adjust this value to nudge labels to the right
                  size = 6,
                  max.overlaps = 50)  + # 30
  guides(color = guide_legend(override.aes = list(label = ""))) # just to remove the "a" added in fig legend
dev.off()


```

--> `ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X` db identified *EZH2* and *SUZ12* as genes upreg in KO
----> plot `x = odd.ratio` and `y = logadjPval` like JC done after.


Let's add a negative control for ENCODE_ChEA analysis; by selecting 202 random genes and see if we find EZH2 SUZ12 too:
XXX Below not ran; should make for the same number as I have for Gain H3K27me3 here = 1,195 genes

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




