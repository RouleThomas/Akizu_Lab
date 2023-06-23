# Identify EZH1 target

## Objective and data finding
Objective is to identify EZH1 target; without EHZ1 ChIP.

Ideal would be to focus on genes H3K27me3-bound, non bound with EZH2, but bound with SUZ12. --> As 6/16/2023; I did not find NPC/neurons with such ChiP performed; however I found:
- ChIP H9 5 days NPC: H3K27me3, EZH2 from [ENCODE](https://www.encodeproject.org/biosamples/ENCBS018TPT/ 
) 

### ENCODE ChIP H3K27me3 and EZH2

#### Complete re-analysis

YYY Not priority, let's try with the already processed files first.


#### Using already processed files

Let's:
- collect the bed of the peaks from the ChIPs
- assign peak to genes with ChIPseeker
- Filter genes bound with H3K27me3 but NOT with EZH2 = putative EZH1 target 


Here is GEO for [H3K27me3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123199) and [EZH2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95944). There are many files available, let's pick the following and transfer to `input/` folder.


Look at the various bigwig/bed on IGV and decide which bed to use.

--> The optimal files to use seems to be:
- H3K27me3: `input/ENCODE/H3K27me3/GSE95944_ENCFF778VZK_signal_p-value_GRCh38.bigWig` and `input/ENCODE/H3K27me3/GSE95944_ENCFF053EMF_replicated_peaks_GRCh38.bed`
- EZH2: `input/ENCODE/EZH2/GSE123199_ENCFF253ZHN_signal_p-value_GRCh38.bigWig` and `input/ENCODE/EZH2/GSE123199_ENCFF575KEH_conservative_idr_thresholded_peaks_GRCh38.bed` (HERE NOT SURE for the bed, conservative or optimal; I pick conservatei at seems there is more peak... Could be more stringeant)


*NOTE: Not clear optimal bigwig file to use, but found clear ID from the selected file [here](https://github.com/kundajelab/encode-ui-sandbox/blob/master/data/ENCODE2018Core.tsv) so I pick them. For the bed, the replicated peak look more broad; less messy*


#### Assign peak to genes with ChIPseeker



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
library(VennDiagram)


# Import peaks
H3K27me3 = read.table('input/ENCODE/H3K27me3/GSE95944_ENCFF053EMF_replicated_peaks_GRCh38.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, V5=V5, V6=V6, V7=V7, V8=V8, FC=V9, V10=V10) 
EZH2 = read.table('input/ENCODE/EZH2/GSE123199_ENCFF575KEH_conservative_idr_thresholded_peaks_GRCh38.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, V4=V4, V5=V5, V6=V6, V7=V7, V8=V8, FC=V9, V10=V10) 


# Tidy peaks
H3K27me3_gr = makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=TRUE)
EZH2_gr = makeGRangesFromDataFrame(EZH2,keep.extra.columns=TRUE)


gr_list <- list(H3K27me3=H3K27me3_gr, EZH2=EZH2_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_H3K27me3_EZH2.pdf", width=7, height=7)
vennplot(genes)
dev.off()

## Remove the distal intergenic
# Filter out peaks with 'Distal Intergenic' annotation
peakAnnoList_noIntergenic <- lapply(peakAnnoList, function(x) {
  x <- as.data.frame(x) # convert GRanges to dataframe
  x <- x[x$annotation != "Distal Intergenic", ] # remove rows with 'Distal Intergenic'
  return(x)
})

genes = lapply(peakAnnoList_noIntergenic, function(i) unique(i$geneId))

pdf("output/ChIPseeker/overlap_genes_H3K27me3_EZH2_noIntergenic.pdf", width=7, height=7) # 
vennplot(genes)
dev.off()



# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
H3K27me3_annot <- as.data.frame(peakAnnoList[["H3K27me3"]]@anno)
EZH2_annot <- as.data.frame(peakAnnoList[["EZH2"]]@anno)

## Peak distance to TSS
pdf("output/ChIPseeker/DistToTSS.pdf", width=14, height=5)
plotDistToTSS(peakAnnoList)
dev.off()

## Convert entrez gene IDs to gene symbols
H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")

H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(H3K27me3_annot, file="output/ChIPseeker/annotation_H3K27me3.txt", sep="\t", quote=F, row.names=F)
write.table(EZH2_annot, file="output/ChIPseeker/annotation_EZH2.txt", sep="\t", quote=F, row.names=F)
```




#### Identify H3K27me3 non-EZH2 bound genes with deepTools

Convert the assign gene files to ChIPseeker and use bedtools to identify the non-overlapping genes


```bash
conda activate BedToBigwig
# File where the diffbound sites has been assigned to genes
output/ChIPseeker/annotation_H3K27me3.txt
output/ChIPseeker/annotation_EZH2.txt


# filter out the intergenic and convert to bed (remove header manually)
grep -v "Intergenic" output/ChIPseeker/annotation_H3K27me3.txt > output/ChIPseeker/annotation_H3K27me3_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_EZH2.txt > output/ChIPseeker/annotation_EZH2_noIntergenic.bed

# Collect gene ID
awk -F'\t' '(NR==1 || FNR>1) {print $22}' output/ChIPseeker/annotation_H3K27me3_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_H3K27me3_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $22}' output/ChIPseeker/annotation_EZH2_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_EZH2_noIntergenic_geneSymbol.txt


# Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_H3K27me3_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_H3K27me3_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_EZH2_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_EZH2_noIntergenic_as_gtf_geneSymbol.txt


# Filter the gtf
grep -Ff output/ChIPseeker/annotation_H3K27me3_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_H3K27me3_noIntergenic_ENCODE.gtf
grep -Ff output/ChIPseeker/annotation_EZH2_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_EZH2_noIntergenic_ENCODE.gtf


# Isolate genes
## H3K27me3+EZH2 bound
bedtools intersect -wa -a meta/ENCFF159KBI_EZH2_noIntergenic_ENCODE.gtf -b meta/ENCFF159KBI_H3K27me3_noIntergenic_ENCODE.gtf > meta/ENCFF159KBI_EZH2_H3K27me3_noIntergenic_ENCODE.gtf
### Sort and remove dupplicates
sort meta/ENCFF159KBI_EZH2_H3K27me3_noIntergenic_ENCODE.gtf | uniq > meta/ENCFF159KBI_EZH2_H3K27me3_noIntergenic_ENCODE_sort.gtf
## H3K27me3 NON-EZH2 bound
bedtools intersect -v -a meta/ENCFF159KBI_H3K27me3_noIntergenic_ENCODE.gtf -b meta/ENCFF159KBI_EZH2_noIntergenic_ENCODE.gtf > meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE.gtf
### Sort and remove dupplicates
sort meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE.gtf | uniq > meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_sort.gtf


# Count how many genes
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
```
**Nb of genes:**
- Bound witgh H3K27me3: 8,763
- Bound with EZH2: 5,860
- Bound with EZH2 and H3K27me3: 4,702
- Bound with H3K27me3 but NOT EZH2: 4,104


#### fine-tune the putative EZH1-target

To fine-tune the list of putative EZH1-target genes, let's look at the GO of the H3K27me3-bound/non-EZH2-bound and keep the one related to neurons = *stringeant putative-EZH1 target*


I collected from [here](http://geneontology.org/) the list of genes related to neurons. I just typed "neurons" and select "Associations between GO terms and genes or gene products." and output all the list of genes for human.; Also did it with "brain"

```bash
conda activate deseq2

# select the gene list name only
awk -F'\t' '{print $3}' meta/neurons_related_GO_geneList.txt | sort | uniq > meta/neurons_related_GO_geneList_geneSymbol.txt
awk -F'\t' '{print $3}' meta/brain_related_GO_geneList.txt | sort | uniq > meta/brain_related_GO_geneList_geneSymbol.txt
# combine both file
cat meta/neurons_related_GO_geneList_geneSymbol.txt meta/brain_related_GO_geneList_geneSymbol.txt | sort | uniq > meta/neurons_brain_related_GO_geneList_geneSymbol.txt

# Filter our GTF putative EZH1 (H3K27me3-nonEZH2) with only these neurons/brain genes
## Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' meta/neurons_brain_related_GO_geneList_geneSymbol.txt > meta/neurons_brain_related_GO_geneList_as_gtf_geneSymbol.txt
## Filter the gtf
grep -Ff meta/neurons_brain_related_GO_geneList_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_neurons_brain_related_GO.gtf


# Isolate genes  H3K27me3 NON-EZH2 bound AND GO neurons/brain
conda activate BedToBigwig

bedtools intersect -wa -a meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_sort.gtf -b meta/ENCFF159KBI_neurons_brain_related_GO.gtf > meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_neurons_brain_related_GO.gtf
### Sort and remove dupplicates
sort meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_neurons_brain_related_GO.gtf | uniq >  meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_neurons_brain_related_GO_sort.gtf
### Count how many genes
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_neurons_brain_related_GO_sort.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
```
--> `neurons_related_GO_geneList.txt` = 1,212 unique genes with GO neurons

--> `brain_related_GO_geneList.txt` = 1,085 unique genes with GO brain

--> `neurons_brain_related_GO_geneList_geneSymbol.txt` = 1,887 unique genes with GO related to neurons and/or brain

--> `ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_neurons_brain_related_GO_sort.gtf` = 391 unique genes with GO related to neurons and/or brain AND H3K27me3-bound Not-EZH2 bound = **stringeant putative-EZH1 target** 


#### deepTools profile of these genes in our CutRun

Let's check deepTools profile in our CutRun data; so at 8weeks neurons

```bash
# deepTools plot
conda activate deeptools
## Bound with EZH2 and H3K27me3
sbatch scripts/matrix_gene_1kb_THOR_EZH2_H3K27me3_ENCODE.sh # 1368994
## Bound with H3K27me3 but NOT EZH2
sbatch scripts/matrix_gene_1kb_THOR_H3K27me3_NotEZH2_ENCODE.sh # 1369033
## Bound with H3K27me3 but NOT EZH2 AND GO neurons/brain
sbatch scripts/matrix_gene_1kb_THOR_H3K27me3_NotEZH2_ENCODE_GO_neurons_brain.sh # 1488762
```

--> Bound with EZH2 and H3K27me3 show alsmost similar profile WT vs mutants

--> Bound with H3K27me3 but NOT EZH2; Both mutants are much more H3K27me3 than WT... 


#### Proportion of genes EZH-1-2 target ~ EZH1 mutations

Let's try to strenghten the following hypothesis:
- Genes that gain H3K27me3 in KO are EZH2 target (EZH1 goes out so EZH2 is coming and methylate)
- Genes that gain H3K27me3 in HET are EZH1 target 

For that, look at the proportion of each genes Gain / Lost in HET / KO; the following:
- Unknown
- EZH2-target
- EZH1-target
- EZH1-stringeant-target


By doing:
- In Bash (awk) output the list of genes for each of the **H3K27me3 categories** (ie. Gain in HET, Lost in HET,...)
- In Bash (awk) output the list of genes for each of the **EZH categories** (ie. EZH1 target, EZH2 target,...)
- In R, import these tables and check for the **proportion of EZH categories with each of the H3K27me3 categories**

--> Analysis in `output/proportion`
```bash
conda activate deseq2
# cp gene list H3K27me3 from 003__CutRun
cp ../003__CutRun/output/ChIPseeker/annotation_WTvsHET_qval10_DownHET_noIntergenic_geneSymbol.txt output/proportion
cp ../003__CutRun/output/ChIPseeker/annotation_WTvsHET_qval10_UpHET_noIntergenic_geneSymbol.txt output/proportion
cp ../003__CutRun/output/ChIPseeker/annotation_WTvsKO_qval10_UpKO_noIntergenic_geneSymbol.txt output/proportion
cp ../003__CutRun/output/ChIPseeker/annotation_WTvsKO_qval10_DownKO_noIntergenic_geneSymbol.txt output/proportion
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("rtracklayer")

# Extract gene list EZH1 (raw and stringenat) categories
## Raw EZH1 target
#### Read GTF file
gtf_file <- "meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_sort.gtf"
gtf_data <- import(gtf_file)
#### Extract gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_name <- elementMetadata(gene_data)$gene_name
#### Save
write(gene_name, file = "output/proportion/H3K27me3_NotEZH2_ENCODE_geneSymbol.txt")


## Stringeant EZH1 target
#### Read GTF file
gtf_file <- "meta/ENCFF159KBI_H3K27me3_NotEZH2_noIntergenic_ENCODE_neurons_brain_related_GO_sort.gtf"
gtf_data <- import(gtf_file)
#### Extract gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_name <- elementMetadata(gene_data)$gene_name
#### Save
write(gene_name, file = "output/proportion/H3K27me3_NotEZH2_ENCODE_neurons_brain_related_GO_geneSymbol.txt")
```
*NOTE: extract gene list from gtf*

**H3K27me3 categories**:
- `output/proportion/annotation_WTvsHET_qval10_DownHET_noIntergenic_geneSymbol.txt` = LostHET
- `output/proportion/annotation_WTvsHET_qval10_UpHET_noIntergenic_geneSymbol.txt` = GainHET
- `output/proportion/annotation_WTvsKO_qval10_DownKO_noIntergenic_geneSymbol.txt` = LostKO
- `output/proportion/annotation_WTvsKO_qval10_UpKO_noIntergenic_geneSymbol.txt` = GainKO

**EZH categories**:
- `output/ChIPseeker/annotation_EZH2_noIntergenic_geneSymbol.txt` = EZH2_Target
- `output/proportion/H3K27me3_NotEZH2_ENCODE_geneSymbol.txt` = EZH1_rawTarget
- `output/proportion/H3K27me3_NotEZH2_ENCODE_neurons_brain_related_GO_geneSymbol.txt` = EZH1_stringeantTarget

Now in R generate the **proportion of EZH categories with each of the H3K27me3 categories** plot:

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("rtracklayer")

# Import all gene categories
LostHET = read_csv("output/proportion/annotation_WTvsHET_qval10_DownHET_noIntergenic_geneSymbol.txt", col_names = c("LostHet"))
GainHET = read_csv("output/proportion/annotation_WTvsHET_qval10_UpHET_noIntergenic_geneSymbol.txt", col_names = c("GainHET"))
LostKO = read_csv("output/proportion/annotation_WTvsKO_qval10_DownKO_noIntergenic_geneSymbol.txt", col_names = c("LostKO"))
GainKO = read_csv("output/proportion/annotation_WTvsKO_qval10_UpKO_noIntergenic_geneSymbol.txt", col_names = c("GainKO"))

EZH2_Target = read_csv("output/ChIPseeker/annotation_EZH2_noIntergenic_geneSymbol.txt", col_names = c("EZH2_Target"))
EZH1_Target = read_csv("output/proportion/H3K27me3_NotEZH2_ENCODE_geneSymbol.txt", col_names = c("EZH1_rawTarget"))
EZH1_stringeantTarget = read_csv("output/proportion/H3K27me3_NotEZH2_ENCODE_neurons_brain_related_GO_geneSymbol.txt", col_names = c("EZH1_stringeantTarget"))

# Create a data frame for each gene category and annotate the category
LostHET <- LostHET %>% rename(Gene = LostHet) %>% mutate(Category = "LostHET")
GainHET <- GainHET %>% rename(Gene = GainHET) %>% mutate(Category = "GainHET")
LostKO <- LostKO %>% rename(Gene = LostKO) %>% mutate(Category = "LostKO")
GainKO <- GainKO %>% rename(Gene = GainKO) %>% mutate(Category = "GainKO")

# Combine all gene categories into one data frame
all_genes <- bind_rows(LostHET, GainHET, LostKO, GainKO)

# Annotate each gene with the EZH category it falls into
all_genes <- all_genes %>%
  mutate(
    EZH_Category = case_when(
      .$Gene %in% EZH2_Target$EZH2_Target ~ "EZH2_Target",
      .$Gene %in% EZH1_Target$EZH1_rawTarget ~ "EZH1_Target",
      .$Gene %in% EZH1_stringeantTarget$EZH1_stringeantTarget ~ "EZH1_stringeantTarget",
      TRUE ~ "Unknown"
    )
  )

# Compute proportions
proportions <- all_genes %>%
  group_by(Category, EZH_Category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Category) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

proportions$Category <- factor(proportions$Category, 
                                   levels = c("LostHET", "GainHET", "LostKO", "GainKO"))



# Plot the proportions
pdf("output/proportion/H3K27me3_EZH_Category.pdf", width=7, height=6)
ggplot(proportions, aes(x = Category, y = Proportion, fill = EZH_Category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(Proportion)), 
            position = position_fill(vjust = 0.5),
            size = 5) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Proportion", x = "Category", fill = "EZH Category") +
  theme_bw()
dev.off()

# Correct the code to also display EZH1-stringeant
# Import all gene categories
LostHET = read_csv("output/proportion/annotation_WTvsHET_qval10_DownHET_noIntergenic_geneSymbol.txt", col_names = c("LostHet"))
GainHET = read_csv("output/proportion/annotation_WTvsHET_qval10_UpHET_noIntergenic_geneSymbol.txt", col_names = c("GainHET"))
LostKO = read_csv("output/proportion/annotation_WTvsKO_qval10_DownKO_noIntergenic_geneSymbol.txt", col_names = c("LostKO"))
GainKO = read_csv("output/proportion/annotation_WTvsKO_qval10_UpKO_noIntergenic_geneSymbol.txt", col_names = c("GainKO"))

EZH2_Target = read_csv("output/ChIPseeker/annotation_EZH2_noIntergenic_geneSymbol.txt", col_names = c("EZH2_Target"))
EZH1_Target = read_csv("output/proportion/H3K27me3_NotEZH2_ENCODE_geneSymbol.txt", col_names = c("EZH1_rawTarget"))
EZH1_stringeantTarget = read_csv("output/proportion/H3K27me3_NotEZH2_ENCODE_neurons_brain_related_GO_geneSymbol.txt", col_names = c("EZH1_stringeantTarget"))

# Remove EZH1_stringeantTarget genes from EZH1_Target
EZH1_Target <- EZH1_Target[!(EZH1_Target$EZH1_rawTarget %in% EZH1_stringeantTarget$EZH1_stringeantTarget), ]

# Create a data frame for each gene category and annotate the category
LostHET <- LostHET %>% rename(Gene = LostHet) %>% mutate(Category = "LostHET")
GainHET <- GainHET %>% rename(Gene = GainHET) %>% mutate(Category = "GainHET")
LostKO <- LostKO %>% rename(Gene = LostKO) %>% mutate(Category = "LostKO")
GainKO <- GainKO %>% rename(Gene = GainKO) %>% mutate(Category = "GainKO")


# Combine all gene categories into one data frame
all_genes <- bind_rows(LostHET, GainHET, LostKO, GainKO)

# Annotate each gene with the EZH category it falls into
all_genes <- all_genes %>%
  mutate(
    EZH_Category = case_when(
      .$Gene %in% EZH2_Target$EZH2_Target ~ "EZH2_Target",
      .$Gene %in% EZH1_Target$EZH1_rawTarget ~ "EZH1_Target",
      .$Gene %in% EZH1_stringeantTarget$EZH1_stringeantTarget ~ "EZH1_stringeantTarget",
      TRUE ~ "Unknown"
    )
  )

# Compute proportions
proportions <- all_genes %>%
  group_by(Category, EZH_Category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Category) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

proportions$Category <- factor(proportions$Category, 
                                   levels = c("LostHET", "GainHET", "LostKO", "GainKO"))



# Plot the proportions
pdf("output/proportion/H3K27me3_EZH_Category_all.pdf", width=7, height=6)
ggplot(proportions, aes(x = Category, y = Proportion, fill = EZH_Category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(Proportion)), 
            position = position_fill(vjust = 0.5),
            size = 5) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Proportion", x = "Category", fill = "EZH Category") +
  theme_bw()
dev.off()
```
*NOTE: It fail to show EZH1_stringeantTarget as it falls within EZH1_Target; so I corrected the script afterward to remove the EZH1_stringeant that are within EZH1_Target*




