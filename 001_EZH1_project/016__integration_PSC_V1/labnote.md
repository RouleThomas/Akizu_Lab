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
sbatch scripts/bamtobigwig-DBA_NORM_TMM_unique_SF_H3K27me3.sh # 32649115 xxx
sbatch scripts/bamtobigwig-DBA_NORM_RLE_unique_SF_H3K27me3.sh # 32649206 xxx
sbatch scripts/bamtobigwig-DBA_NORM_LIB_unique_SF_H3K27me3.sh # 32649234 xxx

# Reciprocal SF (=1/SF) from DiffBind
sbatch scripts/bamtobigwig-DBA_NORM_TMM_unique_reciprocalSF_H3K27me3.sh # 32650322 xxx

```

--> using SF directly from DiffBind is bad, very heterogeneous replicates, for each norm method...


XXX If awful test reciprocal jut in case



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
```

--> all good




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


# Import THOR peaks
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


```



# RNAseq integration 



--> DEGs is redo as in `015__RNAseq`


## PSC KO vs WT


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

# H3K27me3
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






















# H3K27me3
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







## PSC KOEF1aEZH1 vs WT


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
# Convert bam to bigwig (keeping duplicate!)
conda activate deeptools

sbatch scripts/bamtobigwig_Ferguson_1.sh # 33969577 ok
sbatch scripts/bamtobigwig_Ferguson_2.sh # 33969578 ok
sbatch scripts/bamtobigwig_Ferguson_3.sh # 33969581 ok

# Convert bigwig to bedgraph
conda activate BedToBigwig

sbatch scripts/BedToBigwig_Ferguson.sh # 33973549 ok

# Remove blacklist regions
sbatch scripts/BedintersectBlacklist_Ferguson.sh # 33974981 xxx
```

Use Python to identify local maxima, quantify the height fir the 99th percentile peak

```bash
srun --mem=250g --pty bash -l

# Identify local maxima
python scripts/LocalMaxima_Ferguson.py
#  calculate the 99th percentile of the signal heights (score) in the local maxima files.
python scripts/Percentile99_Ferguson.py

# normalize AB per AB (using WT sample 1st replicate as reference)
python scripts/norm_H3K27me3_Ferguson.py
python scripts/norm_SUZ12_Ferguson.py
python scripts/norm_EZH2_Ferguson.py
python scripts/norm_IGG_Ferguson.py

#python scripts/norm_EZH1_Ferguson.py
```
--> Works!

- *NOTE: **Local Maxima** = value is higher than its neighboring points. In the context of your CUT&RUN data (or other genomic data), local maxima refer to genomic positions where the signal intensity (e.g., read depth or coverage in the bedGraph file) is greater than the signal in the surrounding regions.*
- *NOTE: **Percentile 99** = signal level that is greater than 99% of all other signal values in the dataset.*


Convert normalized bedGraph back to bigwig

```bash
conda activate BedToBigwig

sbatch scripts/BedToBigwig_Norm_Ferguson.sh # 33980900 ok
sbatch scripts/BedToBigwig_Norm_IGG_Ferguson.sh # 33984013 ok

# Subtract Igg signal
conda activate deeptools

sbatch scripts/bigwigCompare_Norm_Ferguson_subtractIGG.sh # 33995282 xxx

xxxy check subtract
```
--> Replicates are very heterogeneous...


Let's try to remove IGG signal.












