# Project and goals 

- Identify mutations in single or multiple CR genes leading to NDD:
    - Collect all CR genes
        - Identify the brain express ones (`GTEx`, `BrainSpan Atlas`)
        - Collect their allele frequency score (`gnomAD`)
    - Filter in the one known as invovled in NDD, cehck db:
        - SFARI
            lists 1,231 genes implicated in autism, with annotations and links to published papers.
        - ARCCUS at CHOP
        - PennMedecineBioBank
        - UK Biobank
        - Decipher
        - GeneMatcher
        - Matchmaker Exchange
    - Find individuals with CR mutations and NDD
        - Enter genes and check if NDD phenotype (`PennMedecineBioBank`, `Decipher`, `GeneMatcher`, `Matchmaker Exchange`)
            - For individual with phenotype, check whether other CR genes is mutated
        - Check if disease associated genes; already known (`OMIM`, `SFARI`)
    - Classify individuals for syndrome severity and look for correlation with single or multiple mutations (hope to find relationship here)
    - Identify group of CR mutated genes leading to NDD
        - Generate summary tables with patient's symptom severity and # of genes mutated

*Additional databases*:
- `STRING`, `GeneMANIA` for gene regulatory networks (see if link between NDD genes and CR genes)  To add in summary tables
- `eQTL` analysis to check whether mutation in NDD is associated with gene expression changes in the brain (GTEx)
- *machine learning models* to predict NDD risk based on mutation patterns in CR genes




# gene list

I combined both manually curated gene list; guided with ChatGPT; **conv archived in my account as: *Human PRC2 Core Subunits*.** with gene list from databasesl; msigdb; C5, C2

--> Xlsx table genereted as `007_NDD_CRmut/CRgenes_*.xlsx`

Generate a clean complete table of unique gene list in R:



```R
# packages
library("tidyverse")
library("readxl")
library("biomaRt")

# import file
clean_manual <- read_excel("meta/CRgenes_compil.xlsx", sheet = 1) %>%
    filter(list %in% c("manual")) %>%
    dplyr::select(geneSymbol, list)
clean_manual_msigdb <- read_excel("meta/CRgenes_compil.xlsx", sheet = 1) %>%
    filter(list %in% c("manual_msigdb")) %>%
    dplyr::select(geneSymbol, list)
clean_msigdb <- read_excel("meta/CRgenes_compil.xlsx", sheet = 1) %>%
    filter(list %in% c("msigdb")) %>%
    dplyr::select(geneSymbol, list)

raw <- read_excel("meta/CRgenes_compil.xlsx", sheet = 2)


# compile information from raw
clean_manual_raw = clean_manual %>%
    left_join(raw %>% dplyr::select(geneSymbol, Class, Complex_Family, Type, ProteinFunction)) %>%
    add_column(Term = "NA")

clean_msigdb_raw = clean_msigdb %>%
    left_join(raw %>% dplyr::select(geneSymbol,Term))
# --> Here same genes in different Term, so need combine Term ID into the same value!
clean_msigdb_combined_raw <- clean_msigdb_raw %>%
  group_by(geneSymbol) %>%
  summarize(Term = paste(unique(Term), collapse = "-"), .groups = 'drop')

clean_manual_msigdb_raw = clean_manual_msigdb %>%
    left_join(raw %>% dplyr::select(geneSymbol,Term)) %>%
    filter(Term != "NA") %>%
    group_by(geneSymbol) %>%
    summarize(Term = paste(unique(Term), collapse = "-"), .groups = 'drop') %>%
    left_join(raw %>% dplyr::select(geneSymbol, Class, Complex_Family, Type, ProteinFunction)) %>%
    filter(Class != "NA")

geneList_CR_V1 = clean_manual_msigdb_raw %>%
    bind_rows(clean_msigdb_combined_raw) %>%
    bind_rows(clean_manual_raw)

### save output: write.table(geneList_CR_V1, "meta/geneList_CR_V1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
### load: geneList_CR_V1 = read_tsv("meta/geneList_CR_V1.txt")

# Add ensembl gene ID for all genes
geneList_CR_V1

## convert gene Ensembl to symbol 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

### Convert Ensembl gene IDs to gene symbols
Ensembl <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                     filters = "external_gene_name",
                     values = geneList_CR_V1$geneSymbol,
                     mart = ensembl)  %>%
    dplyr::rename("geneSymbol" = "external_gene_name") %>%
    as_tibble()
### Merge gene symbols to your dataframe
geneList_CR_Ensembl_V1 <- geneList_CR_V1 %>%
    left_join(Ensembl)

### save output: write.table(geneList_CR_Ensembl_V1, "meta/geneList_CR_Ensembl_V1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
### load: geneList_CR_Ensembl_V1 = read_tsv("meta/geneList_CR_Ensembl_V1.txt")
```

- *NOTE: 842 unique `geneSymbol` corresponding to 938 unique `ensembl_gene_id`*

# Expression in the brain

Let's filter our gene list `meta/geneList_CR_V1.txt` to keep only genes expressed in the brain.
- [GTEx](https://www.gtexportal.org/home): has TPM -expression in many human tissue, including brain. [Here](https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression) median expression gene tissue --> Let's filter out and keep only the genes with tpm>1 = express in the brain: `input/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`



```R
# packages
library("tidyverse")


# import files
GTEx = read_delim("input/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
geneList_CR_Ensembl_V1 = read_tsv("meta/geneList_CR_Ensembl_V1.txt")


# filter brain express genes
## keep all column of brain gene expression (12 tissues)
GTEx_brain = GTEx %>%
    dplyr::select("Name", "Brain - Anterior cingulate cortex (BA24)","Brain - Caudate (basal ganglia)","Brain - Cerebellar Hemisphere","Brain - Cerebellum",	"Brain - Cortex","Brain - Frontal Cortex (BA9)","Brain - Hippocampus","Brain - Hypothalamus","Brain - Nucleus accumbens (basal ganglia)","Brain - Putamen (basal ganglia)","Brain - Spinal cord (cervical c-1)","Brain - Substantia nigra") %>%
    dplyr::rename("ensembl_gene_id" = "Name") %>%
    separate(ensembl_gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>%
    dplyr::select(-version)



## keep the gene if at least one of the tissue show a tpm value > 1
GTEx_all = GTEx_brain %>%
    dplyr::select("ensembl_gene_id") %>% unique()

GTEx_brain_tpm1 = GTEx_brain %>%
  filter(if_any(c(-ensembl_gene_id), ~.x >= 1)) %>%
    dplyr::select("ensembl_gene_id") %>% unique()

GTEx_brain_tpm5 = GTEx_brain %>%
  filter(if_any(c(-ensembl_gene_id), ~.x >= 5)) %>%
    dplyr::select("ensembl_gene_id") %>% unique()

# combine with our CR gene list
geneList_CR_Ensembl_V1 %>% inner_join(GTEx_all)

geneList_CR_Ensembl_V1_tpm1 = geneList_CR_Ensembl_V1 %>% inner_join(GTEx_brain_tpm1)
geneList_CR_Ensembl_V1 %>% inner_join(GTEx_brain_tpm5)



### save output: write.table(geneList_CR_Ensembl_V1_tpm1, "output/GTEx/geneList_CR_Ensembl_V1_tpm1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
### load: geneList_CR_Ensembl_V1 = read_tsv("output/GTEx/geneList_CR_Ensembl_V1_tpm1.txt")


```

- *NOTE: I join gene lsit with GTEx expression data using ensembl gene ID to recover most of the genes; but some are still missing in the GTEx db*


--> 842 unique `geneSymbol` corresponding to 938 unique `ensembl_gene_id`
----> 823 unique `ensembl_gene_id` have GTEx expression data
------> 694 unique `ensembl_gene_id` show tpm >1 in at least one brain tissue
------> 630 unique `ensembl_gene_id` show tpm >5 in at least one brain tissue



# described as NDD related

Let's filter my list of candidate CR genes by keeping only the one already defined as NDD related (autism, or other neurolodevlopmental issue):
- OMIM: Provide disease associated with the genes.
- [SFARI](https://gene.sfari.org/database/human-gene/): provide a lsit of genes implicated in autism, with annotations and links to published papers


## Select OMIM genes related to NDD





```R
library("tidyverse")
library("readxl")

# import OMIM genes + disease
genemap2 = read_excel("input/OMIM/genemap2.xlsx", sheet = 1) %>%
    dplyr::select("Approved Gene Symbol","MIM Number",  "Comments", "Phenotypes", "Mouse Gene Symbol/ID") %>%
    rename("geneSymbol" = "Approved Gene Symbol")



# join with gene CR gene list express in brain
geneList_CR_Ensembl_V1 = read_tsv("output/GTEx/geneList_CR_Ensembl_V1_tpm1.txt")

geneList_CR_Ensembl_V1_genemap2 = geneList_CR_Ensembl_V1 %>%
    inner_join(genemap2) %>%
    filter(Phenotypes != "NA")


### save output: write.table(geneList_CR_Ensembl_V1_genemap2, "output/OMIM/geneList_CR_Ensembl_V1_genemap2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
### load: geneList_CR_Ensembl_V1 = read_tsv("output/OMIM/geneList_CR_Ensembl_V1_genemap2.txt")


```


--> Among the 842 CR genes, 218 has OMIM phenotypes; among them 103 are linked to NDD (detail in GoogleDrive `output/OMIM/CR_OMIM.xlsx`)


Putting together CR_brain express `output/GTEx/geneList_CR_Ensembl_V1_tpm1.txt` with OMIM_NDD `output/OMIM/CR_OMIM.xlsx` and SFARI genes; we end up with 160 unique CR genes realted to NDD




# PennMedecineBioBank
- electronic heatlh record with WGS!

Let's check whether our 160 CR-brain-NDD genes are present in the PennMedecineBioBank.







