# SNX14 Vanessa project

Quick RNAseq analysis using already processed output
- GO with different qval / FC treshold (test different db; obj lipid related Terms)
- GSEA on all terms and output the significant one related to lipid; just GSEA plot


## Filter data

From Shared Folders Google drive; `RNAseq of CB&CX_1mon&1yr.xlsx` --> Filtered log2fc and qval and output gene list

## GO on gene lists


```bash
conda activate deseq2
```

```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")
library("enrichR")
library("biomaRt")




# import DEG gene list
output/GO/geneSymbol_1month_CB_qval05FCless0.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.txt

output/GO/geneSymbol_1month_CX_qval05FCless0.txt
output/GO/geneSymbol_1month_CX_qval05FCmore0.txt
output/GO/geneSymbol_1year_CX_qval05FCless0.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.txt

output/GO/geneSymbol_1month_CB_qval05FCless1.txt
output/GO/geneSymbol_1month_CB_qval05FCmore1.txt
output/GO/geneSymbol_1year_CB_qval05FCless1.txt
output/GO/geneSymbol_1year_CB_qval05FCmore1.txt

output/GO/geneSymbol_1year_CX_qval05FCless1.txt
output/GO/geneSymbol_1year_CX_qval05FCmore1.txt


output/GO/geneSymbol_1month_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.5.txt

output/GO/geneSymbol_1year_CX_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.5.txt

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt
# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Adjusted.P.value (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Adjusted.P.value, decreasing = TRUE), ], 50)  #---> this not run for `*_complete.txt`
down <- head(down[order(down$Adjusted.P.value, decreasing = TRUE), ], 50) #---> this not run for `*_complete.txt`
 
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos <- gos %>% filter(P.value <= 0.05)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC0.pdf", width=12, height=7)
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_GO_BP_1month_CX_qval05FC0.pdf", width=12, height=6)
pdf("output/GO/enrichR_GO_BP_1year_CX_qval05FC0.pdf", width=12, height=6)

pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC1.pdf", width=12, height=9)
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC1.pdf", width=12, height=4)

pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC0.5.pdf", width=12, height=11)
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC0.5.pdf", width=12, height=8)

pdf("output/GO/enrichR_GO_BP_1year_CX_qval05FC0.5.pdf", width=12, height=4)

pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC0.5_filtUp.pdf", width=12, height=11)
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC0.5_filtUp.pdf", width=12, height=8)

pdf("output/GO/enrichR_GO_BP_1month_CB_corr.pdf", width=12, height=5)
pdf("output/GO/enrichR_GO_BP_1month_CX_corr.pdf", width=12, height=5)
pdf("output/GO/enrichR_GO_BP_1year_CX_corr.pdf", width=12, height=3)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="#165CAA", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_BP_1month_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1month_CX_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CX_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_BP_1month_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_BP_1month_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CX_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_BP_1month_CB_qval05FC0.5_complete.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CB_qval05FC0.5_complete.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CX_qval05FC0.5_complete.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_BP_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)



# Define databases for enrichment
dbs <- c("GO_Molecular_Function_2023") 

# import DEG gene list

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt


# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Molecular_Function_2023
down <- edown$GO_Molecular_Function_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

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

pdf("output/GO/enrichR_GO_MF_1month_CB_corr.pdf", width=12, height=5)

pdf("output/GO/enrichR_GO_MF_1year_CB_corr.pdf", width=8, height=2)

pdf("output/GO/enrichR_GO_MF_1month_CX_corr.pdf", width=12, height=3)

pdf("output/GO/enrichR_GO_MF_1year_CX_corr.pdf", width=8, height=2)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO MF pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output

write.table(gos, "output/GO/enrichR_GO_MF_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("GO_Cellular_Component_2023") 

# import DEG gene list
output/GO/geneSymbol_1month_CB_qval05FCless0.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.txt

output/GO/geneSymbol_1month_CX_qval05FCless0.txt
output/GO/geneSymbol_1month_CX_qval05FCmore0.txt
output/GO/geneSymbol_1year_CX_qval05FCless0.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.txt

output/GO/geneSymbol_1month_CB_qval05FCless1.txt
output/GO/geneSymbol_1month_CB_qval05FCmore1.txt
output/GO/geneSymbol_1year_CB_qval05FCless1.txt
output/GO/geneSymbol_1year_CB_qval05FCmore1.txt

output/GO/geneSymbol_1month_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.5.txt

output/GO/geneSymbol_1year_CX_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.5.txt
## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt


# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Cellular_Component_2023
down <- edown$GO_Cellular_Component_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

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
pdf("output/GO/enrichR_GO_CC_1month_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_GO_CC_1year_CB_qval05FC0.pdf", width=12, height=6)
pdf("output/GO/enrichR_GO_BP_1month_CX_qval05FC0.pdf", width=12, height=6)
pdf("output/GO/enrichR_GO_BP_1year_CX_qval05FC0.pdf", width=12, height=6)

pdf("output/GO/enrichR_GO_CC_1month_CB_qval05FC1.pdf", width=12, height=4)
pdf("output/GO/enrichR_GO_CC_1year_CB_qval05FC1.pdf", width=12, height=3)

pdf("output/GO/enrichR_GO_CC_1month_CB_qval05FC0.5.pdf", width=12, height=4)
pdf("output/GO/enrichR_GO_CC_1year_CB_qval05FC0.5.pdf", width=12, height=4)

pdf("output/GO/enrichR_GO_CC_1month_CB_corr.pdf", width=12, height=4)
pdf("output/GO/enrichR_GO_CC_1year_CB_corr.pdf", width=12, height=4)

pdf("output/GO/enrichR_GO_CC_1year_CX_corr.pdf", width=12, height=2)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_CC_1month_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_CC_1year_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_CC_1month_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_CC_1year_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_CC_1month_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_CC_1year_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)


write.table(gos, "output/GO/enrichR_GO_CC_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("Reactome_2022") 

# import DEG gene list
## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt


# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Reactome_2022
down <- edown$Reactome_2022
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("R-HSA-[0-9]+$", "", gos$Term)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_Reactome_2022_1month_CB_corr.pdf", width=10, height=4)

pdf("output/GO/enrichR_Reactome_2022_1year_CB_corr.pdf", width=10, height=5)

pdf("output/GO/enrichR_Reactome_2022_1month_CX_corr.pdf", width=10, height=4)

pdf("output/GO/enrichR_Reactome_2022_1year_CX_corr.pdf", width=10, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output

write.table(gos, "output/GO/enrichR_Reactome_2022_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("KEGG_2019_Mouse") 

# import DEG gene list
output/GO/geneSymbol_1month_CB_qval05FCless0.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.txt

output/GO/geneSymbol_1month_CX_qval05FCless0.txt
output/GO/geneSymbol_1month_CX_qval05FCmore0.txt
output/GO/geneSymbol_1year_CX_qval05FCless0.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.txt

output/GO/geneSymbol_1month_CB_qval05FCless1.txt
output/GO/geneSymbol_1month_CB_qval05FCmore1.txt
output/GO/geneSymbol_1year_CB_qval05FCless1.txt
output/GO/geneSymbol_1year_CB_qval05FCmore1.txt

output/GO/geneSymbol_1month_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.5.txt

output/GO/geneSymbol_1year_CX_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.5.txt

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt


# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_KEGG_2019_Mouse_1month_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_KEGG_2019_Mouse_1year_CB_qval05FC0.pdf", width=12, height=6)
pdf("output/GO/enrichR_GO_BP_1month_CX_qval05FC0.pdf", width=12, height=6)
pdf("output/GO/enrichR_GO_BP_1year_CX_qval05FC0.pdf", width=12, height=6)

pdf("output/GO/enrichR_KEGG_2019_Mouse_1month_CB_qval05FC1.pdf", width=12, height=4)
pdf("output/GO/enrichR_KEGG_2019_Mouse_1year_CB_qval05FC1.pdf", width=12, height=4)

pdf("output/GO/enrichR_KEGG_2019_Mouse_1month_CB_qval05FC0.5.pdf", width=12, height=4)
pdf("output/GO/enrichR_KEGG_2019_Mouse_1year_CB_qval05FC0.5.pdf", width=12, height=4)
pdf("output/GO/enrichR_KEGG_2019_Mouse_1year_CX_qval05FC0.5.pdf", width=12, height=2)

pdf("output/GO/enrichR_KEGG_2019_Mouse_1month_CB_corr.pdf", width=10, height=4)
pdf("output/GO/enrichR_KEGG_2019_Mouse_1year_CB_corr.pdf", width=10, height=4)

pdf("output/GO/enrichR_KEGG_2019_Mouse_1month_CX_corr.pdf", width=10, height=2)
pdf("output/GO/enrichR_KEGG_2019_Mouse_1year_CX_corr.pdf", width=10, height=3)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1month_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1year_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1year_CX_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1month_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1year_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1month_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1year_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1year_CX_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("WikiPathways_2019_Mouse") 

# import DEG gene list
output/GO/geneSymbol_1month_CB_qval05FCless0.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.txt

output/GO/geneSymbol_1month_CX_qval05FCless0.txt
output/GO/geneSymbol_1month_CX_qval05FCmore0.txt
output/GO/geneSymbol_1year_CX_qval05FCless0.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.txt

output/GO/geneSymbol_1month_CB_qval05FCless1.txt
output/GO/geneSymbol_1month_CB_qval05FCmore1.txt
output/GO/geneSymbol_1year_CB_qval05FCless1.txt
output/GO/geneSymbol_1year_CB_qval05FCmore1.txt


output/GO/geneSymbol_1month_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.5.txt

output/GO/geneSymbol_1year_CX_qval05FCless0.5.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.5.txt

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$WikiPathways_2019_Mouse
down <- edown$WikiPathways_2019_Mouse
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("WP[0-9]+$", "", gos$Term)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.pdf", width=12, height=4)

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC1.pdf", width=12, height=4)

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.5.pdf", width=12, height=4)
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.5.pdf", width=12, height=4)

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_corr.pdf", width=12, height=4)

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_corr.pdf", width=12, height=4)

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CX_corr.pdf", width=12, height=4)

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1year_CX_corr.pdf", width=12, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("HDSigDB_Mouse_2021") 

# import DEG gene list

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$HDSigDB_Mouse_2021
down <- edown$HDSigDB_Mouse_2021
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics

pdf("output/GO/enrichR_HDSigDB_Mouse_2021_1month_CB_corr.pdf", width=12, height=4)

pdf("output/GO/enrichR_HDSigDB_Mouse_2021_1year_CB_corr.pdf", width=12, height=7)

pdf("output/GO/enrichR_HDSigDB_Mouse_2021_1month_CX_corr.pdf", width=12, height=2)

pdf("output/GO/enrichR_HDSigDB_Mouse_2021_1year_CX_corr.pdf", width=12, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output

write.table(gos, "output/GO/enrichR_HDSigDB_Mouse_2021_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("Descartes_Cell_Types_and_Tissue_2021") 

# import DEG gene list

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Descartes_Cell_Types_and_Tissue_2021
down <- edown$Descartes_Cell_Types_and_Tissue_2021
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics



pdf("output/GO/enrichR_Descartes_Cell_Types_and_Tissue_2021_1month_CB_corr.pdf", width=12, height=7)
pdf("output/GO/enrichR_Descartes_Cell_Types_and_Tissue_2021_1year_CB_corr.pdf", width=12, height=7)

pdf("output/GO/enrichR_Descartes_Cell_Types_and_Tissue_2021_1month_CX_corr.pdf", width=12, height=7)

pdf("output/GO/enrichR_Descartes_Cell_Types_and_Tissue_2021_1year_CX_corr.pdf", width=12, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output

write.table(gos, "output/GO/enrichR_Descartes_Cell_Types_and_Tissue_2021_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)








# Define databases for enrichment
dbs <- c("Tabula_Muris") 

# import DEG gene list

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Tabula_Muris
down <- edown$Tabula_Muris
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("CL:[0-9]+$", "", gos$Term)

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

pdf("output/GO/enrichR_Tabula_Muris_1month_CB_corr.pdf", width=12, height=7)

pdf("output/GO/enrichR_Tabula_Muris_1year_CB_corr.pdf", width=12, height=9)

pdf("output/GO/enrichR_Tabula_Muris_1month_CX_corr.pdf", width=12, height=7)

pdf("output/GO/enrichR_Tabula_Muris_1year_CX_corr.pdf", width=12, height=9)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output

write.table(gos, "output/GO/enrichR_Tabula_Muris_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("DisGeNET") 

# import DEG gene list

## re analysis (`output/deseq2_corr`)
output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt

output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt
output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$DisGeNET
down <- edown$DisGeNET
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)


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

pdf("output/GO/enrichR_DisGeNET_1month_CB_corr.pdf", width=12, height=12)
pdf("output/GO/enrichR_DisGeNET_1year_CB_corr.pdf", width=12, height=12)

pdf("output/GO/enrichR_DisGeNET_1month_CX_corr.pdf", width=12, height=12)
pdf("output/GO/enrichR_DisGeNET_1year_CX_corr.pdf", width=12, height=12)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output

write.table(gos, "output/GO/enrichR_DisGeNET_1year_CX_corr.txt", sep="\t", row.names=FALSE, quote=FALSE)








# Define databases for enrichment
dbs <- c("KOMP2_Mouse_Phenotypes_2022") # nothing
dbs <- c("KOMP2_Molecular_Function_2022") # nothing
dbs <- c("Elsevier_Pathway_Collection") # nothing lipd
dbs <- c("KEGG_2021_Human")


# import DEG gene list
output/GO/geneSymbol_1month_CB_qval05FCless0.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.txt

output/GO/geneSymbol_1month_CX_qval05FCless0.txt
output/GO/geneSymbol_1month_CX_qval05FCmore0.txt
output/GO/geneSymbol_1year_CX_qval05FCless0.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCless0.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCmore0.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
down <- edown$KEGG_2021_Human
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_GO_BP_1month_CX_qval05FC0.pdf", width=12, height=6)
pdf("output/GO/enrichR_GO_BP_1year_CX_qval05FC0.pdf", width=12, height=6)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="dodgerblue2", "up"="firebrick2")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18),
    legend.position = "none"
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_BP_1month_CX_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_1year_CX_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)




```

--> qval 05 FC 0.5 is the best parameter; the only one displaying lipid-related terms
----> filter to keep when overlap >1 gene



Let's try **gene network** vizualization (may be better in our case); **emaplot**:

```R

# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("org.Mm.eg.db")



# Import genes_cluster list and background list
output/GO/geneSymbol_1month_CB_qval05FCless0.txt
output/GO/geneSymbol_1month_CB_qval05FCmore0.txt
output/GO/geneSymbol_1year_CB_qval05FCless0.txt
output/GO/geneSymbol_1year_CB_qval05FCmore0.txt

output/GO/geneSymbol_1month_CX_qval05FCless0.txt
output/GO/geneSymbol_1month_CX_qval05FCmore0.txt
output/GO/geneSymbol_1year_CX_qval05FCless0.txt
output/GO/geneSymbol_1year_CX_qval05FCmore0.txt

geneSymbol_1month_CB_qval05FCmore05 <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCmore0.5.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1month_CB_qval05FCless05 <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCless0.5.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1month_CB_qval05FC05 = geneSymbol_1month_CB_qval05FCmore05 %>%
  bind_rows(geneSymbol_1month_CB_qval05FCless05)



# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(geneSymbol_1month_CB_qval05FC05$V1), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Mm.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Save GO analyses
GO_summary <- data.frame(ego)

write.csv(GO_summary, "output/GO/geneSymbol_1month_CB_qval05FC05.csv")



# Vizualization

pdf("output/GO/dotplot_geneSymbol_1month_CB_qval05FC05.pdf", width=8, height=15)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_geneSymbol_1month_CB_qval05FC05.pdf", width=8, height=11)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()


```

--> Using the enrichGO method we have nothing related to lipid (maybe the GO database is outdated?)



# GSEA



Let's do GSEA and look for signficant hit related to lipid:
- C2; curated gene list (KEGG, wikipathway, reactome, …)
- C5; ontology
- C8; cell type signature

Then do gsea plots for *term that are significant in 1mo CB but not in 1mo CTX* from the GSEA_GO_BP_V1 analysis (log2fc+pvalue ranking):
*POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION*
*UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS*
*UNSATURATED_FATTY_ACID_METABOLIC_PROCESS*
*LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS*
FATTY_ACID_METABOLIC_PROCESS
CELLULAR_RESPONSE_TO_LIPID
*RESPONSE_TO_LIPID* (this may be redundant to cellular response to lipid, if so just choose one)
*CELLULAR_RESPONSE_TO_OXYGEN_LEVELS*
CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND
RESPONSE_TO_OXYGEN_LEVELS (this may be redundant to cellular response to oxygen containing compound, if so just choose one)



Use `conda activate deseq2` and R:



```R
# Packages
library("tidyverse")
library("clusterProfiler")
library("msigdbr") # BiocManager::install("msigdbr")
library("org.Mm.eg.db")
library("enrichplot") # for gseaplot2()
library("pheatmap")
library("readxl")
library("ggrepel")
library("forcats") 

# import DEGs


CB_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 1) %>%
  dplyr::select(gene_name, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue))
CX_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 3) %>%
  dplyr::select(gene_name, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue))
CB_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 2) %>%
  dplyr::select(Gene.name, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue)) %>%
  dplyr::rename(gene_name = Gene.name)
CX_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 4) %>%
  dplyr::select(Gene.name, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue)) %>%
  dplyr::rename(gene_name = Gene.name)

# import DEGs - re-alanalysis (_corr)
CB_1month <- read.table("output/deseq2_corr/filtered_CB_1month__KO_vs_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()  %>%
  filter(!is.na(GeneSymbol)) %>% # filter to keep only the geneSymbol gene
  dplyr::select(GeneSymbol, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue)) %>%
  unique()

CB_1year <- read.table("output/deseq2_corr/filtered_CB_1year__KO_vs_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()  %>%
  filter(!is.na(GeneSymbol)) %>% # filter to keep only the geneSymbol gene
  dplyr::select(GeneSymbol, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue)) %>%
  unique()

CX_1month <- read.table("output/deseq2_corr/filtered_CX_1month__KO_vs_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()  %>%
  filter(!is.na(GeneSymbol)) %>% # filter to keep only the geneSymbol gene
  dplyr::select(GeneSymbol, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue)) %>%
  unique()

CX_1year <- read.table("output/deseq2_corr/filtered_CX_1year__KO_vs_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()  %>%
  filter(!is.na(GeneSymbol)) %>% # filter to keep only the geneSymbol gene
  dplyr::select(GeneSymbol, log2FoldChange, pvalue) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue)) %>%
  unique()


# import msigdbr cell marker db 
hs_hallmark_sets <- msigdbr(
  species = "Mus musculus", # Replace with species name relevant to your data
  category = "C2"   # From C2 pathways
)
hs_hallmark_sets <- msigdbr(
  species = "Mus musculus", # Replace with species name relevant to your data
  category = "C5"   # From C5 ontology
)


# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- CX_1year$combined_score           ################ CHANGE !!!!!!!!!!!!!!!! ######################
names(lfc_vector) <- CX_1year$GeneSymbol           ################ CHANGE !!!!!!!!!!!!!!!! ######################
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA

## without pvalue cutoff
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    hs_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

gsea_result_df <- data.frame(gsea_results@result)


# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_CX_1year_C2_complete.tsv"        ################ CHANGE !!!!!!!!!!!!!!!! ######################
  )
)

readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_CB_1year_C5_complete.tsv"        ################ CHANGE !!!!!!!!!!!!!!!! ######################
  )
)


readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_CX_1year_C2_complete_corr.tsv"        ################ CHANGE !!!!!!!!!!!!!!!! ######################
  )
)



# lipid-containg term
REACTOME_METABOLISM_OF_LIPIDS
GOBP_RESPONSE_TO_LIPID


# plots


pdf("output/gsea/gsea_CB_1month_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1month_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)
pdf("output/gsea/gsea_CB_1year_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1year_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)

pdf("output/gsea/gsea_CB_1month_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1month_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1year_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)
pdf("output/gsea/gsea_CB_1year_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)

enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GOBP_RESPONSE_TO_LIPID",
  title = "GOBP_RESPONSE_TO_LIPID",
  color.line = "#0d76ff"
)
dev.off()


# plot for term that are significant in 1mo CB but not in 1mo CTX (all condition in same plot; 1 plot per term)
## 20240124 ppt
## Find corresponding gene ID number
gsea_results$ID # TO CHECK WHICH NB (=gene lists genes) TO USE
which(gsea_results$ID == "GOBP_RESPONSE_TO_OXYGEN_LEVELS")

### CX_1month IDs; 250, 1145, 1505, 538, 8267, 1708, 7148, 1108
gsea_results_CX_1month = gsea_results

250 FATTY_ACID_METABOLIC_PROCESS
1145 UNSATURATED_FATTY_ACID_METABOLIC_PROCESS
1505 LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS
538 UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS
8267 POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION
1708 RESPONSE_TO_LIPID (this may be redundant to cellular response to lipid, if so just choose one)
7148 RESPONSE_TO_OXYGEN_LEVELS (this may be redundant to cellular response to oxygen containing compound, if so just choose one)
1108 CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND

4716 CELLULAR_RESPONSE_TO_OXYGEN_LEVELS
10714 CELLULAR_RESPONSE_TO_LIPID

### CX_1year IDs; 4657, 5573, 9494, 688, 4932, 7407, 12234, 6616
gsea_results_CX_1year = gsea_results


4657 FATTY_ACID_METABOLIC_PROCESS
5573 UNSATURATED_FATTY_ACID_METABOLIC_PROCESS
9494 LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS
688 UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS
4932 POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION
7407 RESPONSE_TO_LIPID (this may be redundant to cellular response to lipid, if so just choose one)
12234 RESPONSE_TO_OXYGEN_LEVELS (this may be redundant to cellular response to oxygen containing compound, if so just choose one)
6616 CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND

8738 CELLULAR_RESPONSE_TO_OXYGEN_LEVELS
3563 CELLULAR_RESPONSE_TO_LIPID



### CB_1month IDs; 391, 607, 684, 1070, 851, 37, 675, 42
gsea_results_CB_1month = gsea_results


391 FATTY_ACID_METABOLIC_PROCESS
607 UNSATURATED_FATTY_ACID_METABOLIC_PROCESS
684 LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS
1070 UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS
851 POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION
37 RESPONSE_TO_LIPID (this may be redundant to cellular response to lipid, if so just choose one)
675 RESPONSE_TO_OXYGEN_LEVELS (this may be redundant to cellular response to oxygen containing compound, if so just choose one)
42 CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND

571 CELLULAR_RESPONSE_TO_OXYGEN_LEVELS
82 CELLULAR_RESPONSE_TO_LIPID


### CB_1year IDs; 2223, 10335, 11523, 4074, 9118, 22, 3692, 574
gsea_results_CB_1year = gsea_results


2223 FATTY_ACID_METABOLIC_PROCESS
10335 UNSATURATED_FATTY_ACID_METABOLIC_PROCESS
11523 LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS
4074 UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS
9118 POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION
22 RESPONSE_TO_LIPID (this may be redundant to cellular response to lipid, if so just choose one)
3692 RESPONSE_TO_OXYGEN_LEVELS (this may be redundant to cellular response to oxygen containing compound, if so just choose one)
574 CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND

1865 CELLULAR_RESPONSE_TO_OXYGEN_LEVELS
205 CELLULAR_RESPONSE_TO_LIPID



pdf("output/gsea/gsea_CX_1month-curatedTerms_V1.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1year-curatedTerms_V1.pdf", width=14, height=8)
pdf("output/gsea/gsea_CB_1month-curatedTerms_V1.pdf", width=14, height=8)
pdf("output/gsea/gsea_CB_1year-curatedTerms_V1.pdf", width=14, height=8)
### CX_1month IDs; 250, 1145, 1505, 538, 8267, 1708, 7148, 1108
### CX_1year IDs; 4657, 5573, 9494, 688, 4932, 7407, 12234, 6616
### CB_1month IDs; 391, 607, 684, 1070, 851, 37, 675, 42
### CB_1year IDs; 2223, 10335, 11523, 4074, 9118, 22, 3692, 574
pdf("output/gsea/gsea_CB_1year-curatedTerms_V2.pdf", width=10, height=8)
gseaplot2(gsea_results_CB_1year, geneSetID = c(2223, 10335, 11523, 4074, 9118, 22, 3692, 574))
dev.off()

### CX_1month; 7148, 1708, 8267, 1505
### CX_1year; 12234, 7407, 4932, 9494
### CB_1month; 675, 37, 851, 684
### CB_1year; 3692, 22, 9118, 11523
pdf("output/gsea/gsea_CB_1year-curatedTerms_V3.pdf", width=7, height=7)
gseaplot2(gsea_results_CB_1year, geneSetID = c(3692, 22, 9118, 11523), base_size = 30) 
dev.off()



## heatmap Norm Enrichment Score
### import gsea results

CX_1month_C5 = readr::read_tsv("output/gsea/gsea_results_CX_1month_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1month")
CX_1year_C5 = readr::read_tsv("output/gsea/gsea_results_CX_1year_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1year")
CB_1month_C5 = readr::read_tsv("output/gsea/gsea_results_CB_1month_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1month")
CB_1year_C5 = readr::read_tsv("output/gsea/gsea_results_CB_1year_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1year")

### import gsea results - re-analysis C5
CX_1month_C5 = readr::read_tsv("output/gsea/gsea_results_CX_1month_C5_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1month")
CX_1year_C5 = readr::read_tsv("output/gsea/gsea_results_CX_1year_C5_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1year")
CB_1month_C5 = readr::read_tsv("output/gsea/gsea_results_CB_1month_C5_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1month")
CB_1year_C5 = readr::read_tsv("output/gsea/gsea_results_CB_1year_C5_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1year")


gsea_result_df_tidy = CX_1month_C5 %>%
  bind_rows(CX_1year_C5) %>%
  bind_rows(CB_1month_C5) %>%
  bind_rows(CB_1year_C5)
gsea_result_df_tidy_C5 = gsea_result_df_tidy
### import gsea results - re-analysis C2
CX_1month_C2 = readr::read_tsv("output/gsea/gsea_results_CX_1month_C2_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1month")
CX_1year_C2 = readr::read_tsv("output/gsea/gsea_results_CX_1year_C2_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1year")
CB_1month_C2 = readr::read_tsv("output/gsea/gsea_results_CB_1month_C2_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1month")
CB_1year_C2 = readr::read_tsv("output/gsea/gsea_results_CB_1year_C2_complete_corr.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1year")


gsea_result_df_tidy = CX_1month_C2 %>%
  bind_rows(CX_1year_C2) %>%
  bind_rows(CB_1month_C2) %>%
  bind_rows(CB_1year_C2)
gsea_result_df_tidy_C2 = gsea_result_df_tidy
### Set up the heatmap
desired_ids <- c(
"GOBP_RESPONSE_TO_LIPID",
"GOBP_CELLULAR_RESPONSE_TO_LIPID",
"GOBP_LIPID_METABOLIC_PROCESS"
)

desired_ids <- c(
"GOBP_FATTY_ACID_METABOLIC_PROCESS",
"GOBP_UNSATURATED_FATTY_ACID_METABOLIC_PROCESS",
"GOBP_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS",
"GOBP_UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS",
"GOBP_POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION",
"GOBP_RESPONSE_TO_LIPID",
"GOBP_RESPONSE_TO_OXYGEN_LEVELS",
"GOBP_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND"
)


desired_ids <- c(
"GOBP_RESPONSE_TO_OXYGEN_LEVELS",
"GOBP_RESPONSE_TO_LIPID",
"GOBP_POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION",
"GOBP_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS"
)

desired_ids <- c(
"GOBP_REGULATION_OF_SYNAPSE_MATURATION",
"GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION",
"GOCC_NEURON_TO_NEURON_SYNAPSE",
"GOCC_NEURON_PROJECTION",
"GOBP_GABAERGIC_NEURON_DIFFERENTIATION",
"HP_NEUROMUSCULAR_DYSPHAGIA",
"GOBP_LIPID_LOCALIZATION",
"GOBP_LIPID_EXPORT_FROM_CELL",
"GOBP_REGULATION_OF_LIPID_TRANSPORT",
"GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS",
"GOBP_FATTY_ACID_HOMEOSTASIS",
"GOBP_POSITIVE_REGULATION_OF_UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS",
"GOBP_FATTY_ACID_TRANSPORT",
"GOMF_OXYGEN_CARRIER_ACTIVITY",
"GOBP_OXYGEN_TRANSPORT",
"GOMF_OXYGEN_BINDING",
"GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
"GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES"
)


desired_ids <- c(
"WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER",
"LEIN_CEREBELLUM_MARKERS",
"KONDO_HYPOXIA",
"BURTON_ADIPOGENESIS_11",
"REACTOME_FATTY_ACID_METABOLISM",
"KIM_BIPOLAR_DISORDER_OLIGODENDROCYTE_DENSITY_CORR_DN",
"LEIN_CEREBELLUM_MARKERS",
"LEIN_ASTROCYTE_MARKERS",
"WP_COMPLEMENT_SYSTEM_IN_NEURONAL_DEVELOPMENT_AND_PLASTICITY",
"WP_LIPID_PARTICLES_COMPOSITION"
)


# Filter the data for desired IDs
filtered_data <- gsea_result_df_tidy %>%
  filter(ID %in% desired_ids)

filtered_data$genotype <-
  factor(filtered_data$genotype,
         c("CX_1month", "CX_1year", "CB_1month", "CB_1year"))


pdf("output/gsea/heatmap_GOBP_LIPID.pdf", width=3, height=4)

pdf("output/gsea/heatmap_GOBP_LIPID-curatedTerms_V2.pdf", width=8, height=4)
pdf("output/gsea/heatmap_GOBP_LIPID-curatedTerms_V3.pdf", width=8, height=4)

pdf("output/gsea/heatmap_GOBP_LIPID-curatedTerms_corr.pdf", width=8, height=4)

pdf("output/gsea/heatmap_GOBP_LIPID_corr.pdf", width=12, height=6)

pdf("output/gsea/heatmap_C5_examplesV1.pdf", width=18, height=10)
pdf("output/gsea/heatmap_C2_examplesV1.pdf", width=18, height=10)

ggplot(filtered_data, aes(x=genotype, y=ID, fill=NES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, vjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="NES") +
  geom_text(aes(label=sprintf("%.2f", NES)), 
            color = ifelse(filtered_data$pvalue <= 0.05, "black", "grey50"), 
            size=4) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



## re-analysis select relevant terms
# V1
### --> output all terms significant in CB but not significant in CX
#### Step 1: Filter for significant CB_1month or CB_1year
significant_CB <- gsea_result_df_tidy %>%
  filter((genotype == "CB_1month" | genotype == "CB_1year") & pvalue <= 0.05)
#### Step 2: Identify IDs significant in both CX_1month and CX_1year
significant_CX_both <- gsea_result_df_tidy %>%
  filter(genotype == "CX_1month" & pvalue <= 0.05) %>%
  select(ID) %>%
  intersect(gsea_result_df_tidy %>%
              filter(genotype == "CX_1year" & pvalue <= 0.05) %>%
              select(ID))

# Final step: Exclude IDs that are significant in both CX_1month and CX_1year from the CB results
final_filtered_results <- significant_CB %>%
  filter(!ID %in% significant_CX_both$ID)

# V2 corrected
#### Step 1: Filter for significant CB_1month or CB_1year
significant_CB <- gsea_result_df_tidy %>%
  filter((genotype == "CB_1month" | genotype == "CB_1year") & pvalue <= 0.05)
# Identify IDs significant in either CX_1month or CX_1year
significant_CX_either <- gsea_result_df_tidy %>%
  filter((genotype == "CX_1month" | genotype == "CX_1year") & pvalue <= 0.05) %>%
  dplyr::select(ID) %>%
  distinct()
# Final step: Exclude IDs that are significant in either CX_1month or CX_1year from the CB results
final_filtered_results <- significant_CB %>%
  filter(!ID %in% significant_CX_either$ID)



readr::write_tsv(
  final_filtered_results,
  file.path("output/gsea/final_filtered__results_C2_corr.tsv"        ################ CHANGE !!!!!!!!!!!!!!!! ######################
  )
)



# dotplot with term highlighted
## Term to highlight C2 db
gsea_result_df_tidy_C2

desired_ids <- c(
"REACTOME_FATTY_ACID_METABOLISM",
"WP_LIPID_PARTICLES_COMPOSITION",
"REACTOME_TRANSPORT_OF_FATTY_ACIDS"
)

gsea_result_df_tidy_C2_sample = gsea_result_df_tidy_C2 %>% 
  filter(genotype == "CB_1year",
         NES >0,
         pvalue <0.15) %>%
  dplyr::select(-genotype)


# Creating a new column for color
gsea_result_df_tidy_C2_sample$color <- ifelse(gsea_result_df_tidy_C2_sample$ID %in% desired_ids, "red", "grey")

pdf("output/gsea/dotplot_C2_CB_1year_pos.pdf", width=10, height=8)
ggplot(gsea_result_df_tidy_C2_sample, aes(x = pvalue, y = NES, label = ID)) +
  # Add grey points for all IDs
  geom_point(aes(color = color), size = 2, alpha = 0.5) +
  # Add red points for desired IDs on top
  geom_point(data = subset(gsea_result_df_tidy_C2_sample, ID %in% desired_ids),
             aes(color = color), size = 4, alpha = 0.8) +
  scale_color_identity() + # Use the colors as is
  # Add text labels with repelling
  geom_text_repel(
    data = subset(gsea_result_df_tidy_C2_sample, ID %in% desired_ids),
    aes(label = ID), 
    box.padding = unit(0.35, "lines"), 
    point.padding = unit(0.5, "lines"),
    size = 5,
    color = "black"
  ) +
  # Add theme and vertical line
  theme_bw() +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  labs(title = "", x = "p.value", y = "Normalized Enrichment Score (NES)") +
  theme(
    legend.position = "none", # Remove legend
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )+
  ylim(1,1.7) # value pos: 1,1.7; value neg: -1.4,-1.20
dev.off()

## Term to highlight C5 db
gsea_result_df_tidy_C5

desired_ids <- c(
"GOBP_OXYGEN_TRANSPORT",
"GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
"GOBP_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
"GOBP_LIPID_LOCALIZATION",
"GOBP_REGULATION_OF_LIPID_LOCALIZATION",
"GOBP_LIPID_EXPORT_FROM_CELL",
"GOBP_REGULATION_OF_LIPID_TRANSPORT",
"GOBP_NEGATIVE_REGULATION_OF_LIPID_METABOLIC_PROCESS",
"GOBP_NEGATIVE_REGULATION_OF_LIPID_BIOSYNTHETIC_PROCESS",
"GOBP_POSITIVE_REGULATION_OF_LIPID_LOCALIZATION",
"GOBP_POSITIVE_REGULATION_OF_LIPID_TRANSPORT",
"GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS",
"GOMF_IRON_ION_BINDING",
"GOBP_POSITIVE_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
"GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
"GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
"GOBP_RESPONSE_TO_LIPID",
"GOBP_SEQUESTERING_OF_IRON_ION",
"GOMF_LIPID_KINASE_ACTIVITY",
"GOBP_REGULATION_OF_IRON_ION_TRANSMEMBRANE_TRANSPORT",
"GOBP_CELLULAR_RESPONSE_TO_LIPID",
"GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
"GOBP_POSITIVE_REGULATION_OF_UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS",
"GOBP_FATTY_ACID_TRANSPORT",
"GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT",
"GOBP_NEGATIVE_REGULATION_OF_FATTY_ACID_BETA_OXIDATION",
"GOBP_FATTY_ACID_HOMEOSTASIS",
"GOBP_NEGATIVE_REGULATION_OF_FATTY_ACID_METABOLIC_PROCESS"
)

gsea_result_df_tidy_C5_sample = gsea_result_df_tidy_C5 %>% 
  filter(genotype == "CX_1month",
         NES <0,
         pvalue <0.15) %>%
  dplyr::select(-genotype)


# Creating a new column for color
gsea_result_df_tidy_C5_sample$color <- ifelse(gsea_result_df_tidy_C5_sample$ID %in% desired_ids, "red", "grey")

pdf("output/gsea/dotplot_C5_CX_1month_neg.pdf", width=10, height=8)
ggplot(gsea_result_df_tidy_C5_sample, aes(x = pvalue, y = NES, label = ID)) +
  # Add grey points for all IDs
  geom_point(aes(color = color), size = 2, alpha = 0.5) +
  # Add red points for desired IDs on top
  geom_point(data = subset(gsea_result_df_tidy_C5_sample, ID %in% desired_ids),
             aes(color = color), size = 4, alpha = 0.8) +
  scale_color_identity() + # Use the colors as is
  # Add text labels with repelling
  geom_text_repel(
    data = subset(gsea_result_df_tidy_C5_sample, ID %in% desired_ids),
    aes(label = ID), 
    box.padding = unit(0.35, "lines"), 
    point.padding = unit(0.5, "lines"),
    size = 3,
    color = "black",
    max.overlaps = 50
  ) +
  # Add theme and vertical line
  theme_bw() +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  labs(title = "", x = "p.value", y = "Normalized Enrichment Score (NES)") +
  theme(
    legend.position = "none", # Remove legend
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  ylim(-1.4,-1.20) # value pos: 1,1.7 ; value neg: -1.4,-1.20
dev.off()




# Remove some term for C5 and add class (option 1) )

# Oxygen-related processes
Oxygen = c(
  "GOBP_OXYGEN_TRANSPORT",
  "GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
  "GOBP_POSITIVE_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
  "GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS"
)

# Lipid-related processes
Lipid = c(
  "GOBP_LIPID_LOCALIZATION",
  "GOBP_LIPID_EXPORT_FROM_CELL",
  "GOBP_REGULATION_OF_LIPID_TRANSPORT",
  "GOBP_RESPONSE_TO_LIPID",
  "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS",
  "GOBP_FATTY_ACID_TRANSPORT",
  "GOMF_LIPID_KINASE_ACTIVITY",
  "GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT",
  "GOBP_FATTY_ACID_HOMEOSTASIS"
)

# Iron-related processes
Iron = c(
  "GOMF_IRON_ION_BINDING",
  "GOBP_SEQUESTERING_OF_IRON_ION",
  "GOBP_REGULATION_OF_IRON_ION_TRANSMEMBRANE_TRANSPORT"
)



gsea_result_df_tidy_C5_sample = gsea_result_df_tidy_C5 %>% 
  filter(genotype == "CB_1year",
         NES >0,
         pvalue <0.15) %>%
  dplyr::select(-genotype)

special_ids <- c(Oxygen, Lipid, Iron)

# Creating a new column for color based on category
gsea_result_df_tidy_C5_sample$color <- case_when(
  gsea_result_df_tidy_C5_sample$ID %in% Lipid ~ "red",
  gsea_result_df_tidy_C5_sample$ID %in% Oxygen ~ "blue",
  gsea_result_df_tidy_C5_sample$ID %in% Iron ~ "black",
  TRUE ~ "grey" # Default color
)

# Creating a new column for size based on category membership
gsea_result_df_tidy_C5_sample$size <- ifelse(gsea_result_df_tidy_C5_sample$ID %in% special_ids, 4, 2)

#pdf("output/gsea/dotplot_C5_CB_1year_neg_filtV1.pdf", width=10, height=8)
pdf("output/gsea/dotplot_C5_CB_1year_pos_filtV1.pdf", width=10, height=8)
ggplot(gsea_result_df_tidy_C5_sample, aes(x = pvalue, y = NES, label = ID)) +
  geom_point(aes(color = color, size = size), alpha = 0.5) +
  scale_color_identity() +
  scale_size_identity() +
  geom_text_repel(
    data = subset(gsea_result_df_tidy_C5_sample, ID %in% special_ids),
    aes(label = ID), 
    box.padding = unit(0.35, "lines"), 
    point.padding = unit(0.5, "lines"),
    size = 4,
    color = "black",
    max.overlaps = 50
  ) +
  theme_bw() +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  labs(title = "", x = "p.value", y = "Normalized Enrichment Score (NES)") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  ylim(1,1.7) # value pos: 1,1.7 ; value neg: -1.4,-1.20
dev.off()



# Remove some term for C5 and add class (option 1) and only display the top 150 signficant term in waterfall plot )

# Oxygen-related processes
Oxygen = c(
  "GOBP_OXYGEN_TRANSPORT",
  "GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
  "GOBP_POSITIVE_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
  "GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS"
)

# Lipid-related processes
Lipid = c(
  "GOBP_LIPID_LOCALIZATION",
  "GOBP_LIPID_EXPORT_FROM_CELL",
  "GOBP_REGULATION_OF_LIPID_TRANSPORT",
  "GOBP_RESPONSE_TO_LIPID",
  "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS",
  "GOBP_FATTY_ACID_TRANSPORT",
  "GOMF_LIPID_KINASE_ACTIVITY",
  "GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT",
  "GOBP_FATTY_ACID_HOMEOSTASIS"
)

# Iron-related processes
Iron = c(
  "GOMF_IRON_ION_BINDING",
  "GOBP_SEQUESTERING_OF_IRON_ION",
  "GOBP_REGULATION_OF_IRON_ION_TRANSMEMBRANE_TRANSPORT"
)
special_ids <- c(Oxygen, Lipid, Iron)

### Filter and prepare the dataset
gsea_result_df_tidy_C5_sample <- gsea_result_df_tidy_C5 %>%
  filter(genotype == "CX_1year", NES > 0, pvalue < 0.05) %>%
  mutate(
    Category = case_when(
      ID %in% Lipid ~ "Lipid",
      ID %in% Oxygen ~ "Oxygen",
      ID %in% Iron ~ "Iron",
      TRUE ~ "Other"
    ),
    color = case_when(
      ID %in% Lipid ~ "red",
      ID %in% Oxygen ~ "blue",
      ID %in% Iron ~ "black",
      TRUE ~ "grey" # Default color
    )
  ) %>%
  mutate(Label = ifelse(ID %in% special_ids, as.character(ID), NA_character_)) %>%
  arrange(desc(NES))

# Plotting
pdf("output/gsea/waterfall_C5_CX_1year_pos.pdf", width=8, height=5)
ggplot(gsea_result_df_tidy_C5_sample, aes(x = reorder(ID, -NES), y = NES)) +
  geom_point(aes(color = color), size = 2, alpha = 0.6) +
  scale_color_identity() +
  geom_text_repel(
    aes(label = Label, color = color), 
    box.padding = unit(0.35, "lines"), 
    point.padding = unit(0.5, "lines"), 
    size = 3, 
    na.rm = TRUE,
    max.overlaps = 50
  ) +
  theme_classic() +
  labs(title = "", x = "", y = "NES") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(), # Remove the x-axis text
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) 
dev.off()



















# improve readability with geom_line (option 2); not the chosen one)
Oxygen = c(
  "GOBP_OXYGEN_TRANSPORT",
  "GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
  "GOBP_POSITIVE_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
  "GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS"
)

# Lipid-related processes
Lipid = c(
  "GOBP_LIPID_LOCALIZATION",
  "GOBP_LIPID_EXPORT_FROM_CELL",
  "GOBP_REGULATION_OF_LIPID_TRANSPORT",
  "GOBP_RESPONSE_TO_LIPID",
  "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS",
  "GOBP_FATTY_ACID_TRANSPORT",
  "GOMF_LIPID_KINASE_ACTIVITY",
  "GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT",
  "GOBP_FATTY_ACID_HOMEOSTASIS"
)

# Iron-related processes
Iron = c(
  "GOMF_IRON_ION_BINDING",
  "GOBP_SEQUESTERING_OF_IRON_ION",
  "GOBP_REGULATION_OF_IRON_ION_TRANSMEMBRANE_TRANSPORT"
)
special_ids <- c(Oxygen, Lipid, Iron)

## Adjusting the data frame
gsea_result_df_tidy_C5_sample <- gsea_result_df_tidy_C5 %>% 
  filter(genotype == "CB_1year", NES < 0, pvalue < 0.15) %>%
  dplyr::select(-genotype) %>%
  mutate(
    color = case_when(
      ID %in% Lipid ~ "red",
      ID %in% Oxygen ~ "blue",
      ID %in% Iron ~ "black",
      TRUE ~ "grey" # Default color
    ),
    size = ifelse(ID %in% special_ids, 4, 1) # Assign size based on ID being in special_ids
  )

  
#pdf("output/gsea/dotplot_C5_CX_1month_pos_filtV2.pdf", width=10, height=8)
pdf("output/gsea/dotplot_C5_CB_1year_neg_filtV2.pdf", width=10, height=8)

ggplot(gsea_result_df_tidy_C5_sample, aes(x = pvalue, y = NES)) +
  geom_point(aes(color = color, size = size), alpha = 0.2) + # All points
  scale_color_identity() +
  scale_size_identity() +
  geom_smooth(aes(group = 1), method = "loess", color = "darkgrey", size = 1) + # Smooth line
  geom_point(data = subset(gsea_result_df_tidy_C5_sample, ID %in% special_ids),
             aes(color = color), size = 4, alpha = 0.8) + # Special points with correct colors
  geom_text_repel(
    data = subset(gsea_result_df_tidy_C5_sample, ID %in% special_ids),
    aes(label = ID, color = color), 
    box.padding = unit(0.35, "lines"), 
    point.padding = unit(0.5, "lines"),
    size = 4,
    max.overlaps = 50
  ) +
  guides(color = guide_legend(title = "Category")) + # Add a legend for color
  theme_bw() +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  labs(title = "", x = "p.value", y = "Normalized Enrichment Score (NES)") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  ylim(-1.4,-1.20) # value pos: 1,1.7 ; value neg: -1.4,-1.20

dev.off()





```

--> The GSEA dotplot where very dense, notably for C5, so I generated a more condensed version `*_filt.pdf` and add color per category




Now GSEA with the log2fc ranking (save under `GoogleDrive/*/gsea/gsea_GO_BP_V2_log2FoldChangeRanking.xlsx`)


```R
# Packages
library("tidyverse")
library("clusterProfiler")
library("msigdbr") # BiocManager::install("msigdbr")
library("org.Mm.eg.db")
library("enrichplot") # for gseaplot2()
library("pheatmap")
library("readxl")

# import DEGs


CB_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 1) %>%
  dplyr::select(gene_name, log2FoldChange)
CX_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 3) %>%
  dplyr::select(gene_name, log2FoldChange)
CB_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 2) %>%
  dplyr::select(Gene.name, log2FoldChange)%>%
  dplyr::rename(gene_name = Gene.name)
CX_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 4) %>%
  dplyr::select(Gene.name, log2FoldChange)%>%
  dplyr::rename(gene_name = Gene.name)


# import msigdbr cell marker db 
hs_hallmark_sets <- msigdbr(
  species = "Mus musculus", # Replace with species name relevant to your data
  category = "C2"   # From C2 pathways
)
hs_hallmark_sets <- msigdbr(
  species = "Mus musculus", # Replace with species name relevant to your data
  category = "C5"   # From C5 ontology
)


# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- CX_1year$log2FoldChange           ################ CHANGE !!!!!!!!!!!!!!!! ######################
names(lfc_vector) <- CX_1year$gene_name           ################ CHANGE !!!!!!!!!!!!!!!! ######################
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA

## without pvalue cutoff
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    hs_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

gsea_result_df <- data.frame(gsea_results@result)



readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_CX_1year_C5_complete-log2FoldChangeRanking.tsv"        ################ CHANGE !!!!!!!!!!!!!!!! ######################
  )
)




#### BELOW CODE TO MODIFY!!!!!!!

# lipid-containg term
REACTOME_METABOLISM_OF_LIPIDS
GOBP_RESPONSE_TO_LIPID


# plots


pdf("output/gsea/gsea_CB_1month_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1month_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)
pdf("output/gsea/gsea_CB_1year_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1year_C2-REACTOME_METABOLISM_OF_LIPIDS.pdf", width=14, height=8)

pdf("output/gsea/gsea_CB_1month_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1month_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)
pdf("output/gsea/gsea_CX_1year_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)
pdf("output/gsea/gsea_CB_1year_C5-GOPB_RESPONSE_TO_LIPID.pdf", width=14, height=8)

enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GOBP_RESPONSE_TO_LIPID",
  title = "GOBP_RESPONSE_TO_LIPID",
  color.line = "#0d76ff"
)
dev.off()







## heatmap Norm Enrichment Score
### import gsea results

CX_1month_C5 = readr::read_tsv("output/gsea/gsea_results_CX_1month_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1month")
CX_1year_C5 = readr::read_tsv("output/gsea/gsea_results_CX_1year_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CX_1year")
CB_1month_C5 = readr::read_tsv("output/gsea/gsea_results_CB_1month_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1month")
CB_1year_C5 = readr::read_tsv("output/gsea/gsea_results_CB_1year_C5_complete.tsv") %>%
  dplyr::select(ID, NES, pvalue) %>%
  add_column(genotype = "CB_1year")





gsea_result_df_tidy = CX_1month_C5 %>%
  bind_rows(CX_1year_C5) %>%
  bind_rows(CB_1month_C5) %>%
  bind_rows(CB_1year_C5)

### Set up the heatmap
desired_ids <- c(
"GOBP_RESPONSE_TO_LIPID",
"GOBP_CELLULAR_RESPONSE_TO_LIPID",
"GOBP_LIPID_METABOLIC_PROCESS"
)





# Filter the data for desired IDs
filtered_data <- gsea_result_df_tidy %>%
  filter(ID %in% desired_ids)

filtered_data$genotype <-
  factor(filtered_data$genotype,
         c("CX_1month", "CX_1year", "CB_1month", "CB_1year"))


pdf("output/gsea/heatmap_GOBP_LIPID.pdf", width=3, height=4)



ggplot(filtered_data, aes(x=genotype, y=ID, fill=NES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="NES") +
  geom_text(aes(label=sprintf("%.2f", NES)), 
            color = ifelse(filtered_data$pvalue <= 0.01, "black", "grey50"), 
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()


```

--> **gsea log2fc + pvalue ranking is better**








## GSEA using [webtool](https://www.webgestalt.org/)


Let's try to reproduce Ying analysis; she may have used only the signficiant DEGs as input for GSEA (*this method is NOT recommended*)

Let's create a `*.rnk` file with 2 columns; geneSymbol and rank score --> **let's filter the DEGs only and rank with log2fc**

--> Let's create 2 input files; using pvalue treshold 0.05 and qvalue 0.05

```R
# packages
library("tidyverse")
library("readxl")


# pvalue filtering

# import DEGs
CB_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 1) %>%
  dplyr::select(gene_name, log2FoldChange, pvalue) %>%
  filter(pvalue <=0.05) %>%
  dplyr::select(-pvalue)
CX_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 3) %>%
  dplyr::select(gene_name, log2FoldChange, pvalue) %>%
  filter(pvalue <=0.05) %>%
  dplyr::select(-pvalue)
CB_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 2) %>%
  dplyr::select(Gene.name, log2FoldChange, pvalue) %>%
  dplyr::rename(gene_name = Gene.name)%>%
  filter(pvalue <=0.05) %>%
  dplyr::select(-pvalue)
CX_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 4) %>%
  dplyr::select(Gene.name, log2FoldChange, pvalue) %>%
  dplyr::rename(gene_name = Gene.name)%>%
  filter(pvalue <=0.05) %>%
  dplyr::select(-pvalue)


# Create .rnk file
## Order our DEG
CB_1month_rank = CB_1month %>%
  arrange(desc(log2FoldChange))
CB_1year_rank = CB_1year %>%
  arrange(desc(log2FoldChange))
CX_1month_rank = CX_1month %>%
  arrange(desc(log2FoldChange))
CX_1year_rank = CX_1year %>%
  arrange(desc(log2FoldChange))
### save output

write.table(CB_1month_rank, "output/gsea/CB_1month_rank_pvalue.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(CB_1year_rank, "output/gsea/CB_1year_rank_pvalue.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(CX_1month_rank, "output/gsea/CX_1month_rank_pvalue.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(CX_1year_rank, "output/gsea/CX_1year_rank_pvalue.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# qvalue filtering

# import DEGs
CB_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 1) %>%
  dplyr::select(gene_name, log2FoldChange, padj) %>%
  filter(padj <=0.05) %>%
  dplyr::select(-padj)
CX_1month = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 3) %>%
  dplyr::select(gene_name, log2FoldChange, padj) %>%
  filter(padj <=0.05) %>%
  dplyr::select(-padj)
CB_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 2) %>%
  dplyr::select(Gene.name, log2FoldChange, padj) %>%
  dplyr::rename(gene_name = Gene.name)%>%
  filter(padj <=0.05) %>%
  dplyr::select(-padj)
CX_1year = read_excel("output/gsea/RNAseq of CB&CX_1mon&1yr.xlsx", sheet = 4) %>%
  dplyr::select(Gene.name, log2FoldChange, padj) %>%
  dplyr::rename(gene_name = Gene.name)%>%
  filter(padj <=0.05) %>%
  dplyr::select(-padj)


# Create .rnk file
## Order our DEG
CB_1month_rank = CB_1month %>%
  arrange(desc(log2FoldChange))
CB_1year_rank = CB_1year %>%
  arrange(desc(log2FoldChange))
CX_1month_rank = CX_1month %>%
  arrange(desc(log2FoldChange))
CX_1year_rank = CX_1year %>%
  arrange(desc(log2FoldChange))
### save output

write.table(CB_1month_rank, "output/gsea/CB_1month_rank_padj.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(CB_1year_rank, "output/gsea/CB_1year_rank_padj.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(CX_1month_rank, "output/gsea/CX_1month_rank_padj.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(CX_1year_rank, "output/gsea/CX_1year_rank_padj.rnk", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



```













# Let's generate TPM from the bam

- Transfer bam files to cluster (one done with Novogene and the other with something else; hopefully comparable...)
- Do featureCounts from the bam
- Do tpm from featureCounts


Sample name for 1year:
- CTRL= 171HetCB	175HetCB	474WTCB	(or CX)
- KO= 174MTCB	 177MTCB	(or CX)


```bash
conda activate featurecounts

sbatch scripts/featurecounts_1month_CB.sh # 7865267 ok
sbatch scripts/featurecounts_1month_CX.sh # 7865273 ok
sbatch scripts/featurecounts_1year_CB.sh # 7865290 ok
sbatch scripts/featurecounts_1year_CX.sh # 7865295 ok
```

--> All > 75% alignment. Good!



## Calculate TPM and RPKM

Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
conda activate deseq2
# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM.sh #
# mv all output to output/tpm or rpkm folder
mv output/featurecounts/*tpm* output/tpm/
mv output/featurecounts/*rpkm* output/rpkm/
```

All good. 

If needed to **display gene with TPM**:

Let's generate a heatmap to check expression level of [RESPONSE_TO_LIPID GOBP genes](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/GOBP_RESPONSE_TO_LIPID.html)



```bash
nano output/tpm/GOBP_RESPONSE_TO_LIPID_raw.txt # here copy paste gene list

```


```R
# library
library("biomaRt")
library("readr")
library("stringr")
library("tidyverse")
library("dplyr")
library("ComplexHeatmap")
library("circlize")
library("viridis")
library("reshape2")

# import gene list BP lipid
GOBP_RESPONSE_TO_LIPID_raw <- read_lines("output/tpm/GOBP_RESPONSE_TO_LIPID_raw.txt") 
GOBP_RESPONSE_TO_LIPID =  tibble(geneSymbol = unlist(str_split(GOBP_RESPONSE_TO_LIPID_raw, ",")))

# create df for file renaming
ID	genotype	time	tissue	replicate
S_CB_KO1	KO	1month	CB	R1
S_CB_KO2	KO	1month	CB	R2
S_CB_KO3	KO	1month	CB	R3
S_CB_WT1	WT	1month	CB	R1
S_CB_WT2	WT	1month	CB	R2
S_CB_WT3	WT	1month	CB	R3
S_CX_KO1	KO	1month	CX	R1
S_CX_KO2	KO	1month	CX	R2
S_CX_KO3	KO	1month	CX	R3
S_CX_WT1	WT	1month	CX	R1
S_CX_WT2	WT	1month	CX	R2
S_CX_WT3	WT	1month	CX	R3
171HetCB	WT	1year	CB	R1
174MTCB	KO	1year	CB	R1
175HetCB	WT	1year	CB	R2
177MTCB	KO	1year	CB	R2
474WTCB	WT	1year	CB	R3
171HetCX	WT	1year	CX	R1
174MTCX	KO	1year	CX	R1
175HetCX	WT	1year	CX	R2
177MTCX	KO	1year	CX	R2
474WTCX	WT	1year	CX	R3



fileName <- tibble(
  ID = c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", 
         "S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3", 
         "171HetCB", "174MTCB", "175HetCB", "177MTCB", "474WTCB", 
         "171HetCX", "174MTCX", "175HetCX", "177MTCX", "474WTCX"),
  genotype = c("KO", "KO", "KO", "WT", "WT", "WT", 
               "KO", "KO", "KO", "WT", "WT", "WT", 
               "WT", "KO", "WT", "KO", "WT", 
               "WT", "KO", "WT", "KO", "WT"),
  time = c(rep("1month", 12), rep("1year", 10)),
  tissue = c(rep("CB", 6), rep("CX", 6), rep("CB", 5), rep("CX", 5)),
  replicate = c(rep("R1", 2), rep("R2", 2), rep("R3", 2), rep("R1", 2), rep("R2", 2), rep("R3", 2), "R1", "R1", "R2", "R2", "R3", "R1", "R1", "R2", "R2", "R3"), 
  new_ID = c("1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
             "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
) 


## import tpm
#### Generate TPM for ALL samples
#### collect all samples ID
samples_1month <- c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", "S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3")

## Make a loop for importing all tpm data and keep only ID and count column
sample_data_1month <- list()
for (sample in samples_1month) {
  sample_data_1month[[sample]] <- read_delim(paste0("../001__RNAseq/output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.bam.1month_rnaseqyz072420.")) %>%
    rename(!!sample := starts_with("output.bam.1month_rnaseqyz072420."))
}

samples_1year <- c("171HetCB", "174MTCB", "175HetCB", "177MTCB", "474WTCB", "171HetCX", "174MTCX", "175HetCX", "177MTCX", "474WTCX")

sample_data_1year <- list()
for (sample in samples_1year) {
  sample_data_1year[[sample]] <- read_delim(paste0("../001__RNAseq/output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.bam.1year_AL1804271_R2_new_analysis.")) %>%
    rename(!!sample := starts_with("output.bam.1year_AL1804271_R2_new_analysis."))
}



## Merge all dataframe into a single one
tpm_all_sample_1month <- purrr::reduce(sample_data_1month, full_join, by = "Geneid")
tpm_all_sample_1year <- purrr::reduce(sample_data_1year, full_join, by = "Geneid")

tpm_all_sample = tpm_all_sample_1month %>%
  left_join(tpm_all_sample_1year)



# write.csv(tpm_all_sample, file="../001__RNAseq/output/tpm/tpm_all_sample.txt")
### If need to import: tpm_all_sample <- read_csv("../001__RNAseq/output/tpm/tpm_all_sample.txt") %>% dplyr::select(-("...1"))#To import

# plot some genes
tpm_all_sample_tidy <- tpm_all_sample %>%
  gather(key = 'ID', value = 'tpm', -Geneid) %>%
  left_join(fileName)
  
  
tpm_all_sample_tidy$Geneid <- gsub("\\..*", "", tpm_all_sample_tidy$Geneid)


tpm_all_sample_tidy_geneOnly = tpm_all_sample_tidy %>%
  dplyr::select(Geneid)

## convert gene Ensembl to symbol 
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
############ IF ensembl fail; try using mirror or previosu version
ensembl <-  useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl", 
                   mirror = "useast") # uswest useast asia

ensembl <-  useEnsembl(biomart = 'ensembl', 
                       dataset = 'mmusculus_gene_ensembl',
                       version = 110)
############ (more info: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html)


# Convert Ensembl gene IDs to gene symbols
genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = tpm_all_sample_tidy_geneOnly$Geneid,
                     mart = ensembl)

# Merge gene symbols to your dataframe
tpm_all_sample_tidy <- left_join(tpm_all_sample_tidy, genesymbols, 
                                 by = c("Geneid" = "ensembl_gene_id")) %>%
                       unique() 
## save output:  write.table(tpm_all_sample_tidy, sep = "\t", quote = FALSE, row.names=FALSE, file="output/tpm/tpm_all_sample_tidy.txt")


## HERE GOOD!!!!!!!! read: tpm_all_sample_tidy <- read.csv("output/tpm/tpm_all_sample_tidy.txt", sep = "\t", header = TRUE)



tpm_all_sample_tidy_GOBP_RESPONSE_TO_LIPID = GOBP_RESPONSE_TO_LIPID %>%
  rename("external_gene_name" = "geneSymbol") %>%
  left_join(tpm_all_sample_tidy) %>%
  mutate(TPM = log2(tpm +1)) %>%
  filter(!is.na(external_gene_name) & external_gene_name != "") %>% # remove geneanme NA
  unique()  %>%
  group_by(external_gene_name, new_ID) %>%
  summarise(TPM = mean(TPM, na.rm = TRUE), .groups = 'drop') # clean  data


# heatmap
## CB_1month
CB_1month <- tpm_all_sample_tidy_GOBP_RESPONSE_TO_LIPID %>%
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3")) 

CB_1month$new_ID_ordered <- factor(CB_1month$new_ID, levels = c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                                                  "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3"))
### Order the rows based on the new ID order and TPM values
CB_1month <- CB_1month %>%
  arrange(new_ID_ordered, desc(TPM))
### Pivot the data to a long format suitable for ggplot
long_df <- CB_1month %>%
  pivot_longer(cols = TPM, names_to = "Condition", values_to = "Expression")
### Raw uncluseted heatmap
pdf("output/tpm/heatmap_CB_1month-GOBP_RESPONSE_TO_LIPID.pdf", width=4, height=8)
ggplot(long_df, aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()
### Raw uncluseted heatmap with only DEGs
#### import DEGs

geneSymbol_1month_CB_qval05FCmore0 <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCmore0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1month_CB_qval05FCless0 <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCless0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1month_CB_qval05FC0 = geneSymbol_1month_CB_qval05FCmore0 %>%
  bind_rows(geneSymbol_1month_CB_qval05FCless0) %>%
  rename("external_gene_name" = "V1")

pdf("output/tpm/heatmap_CB_1month-GOBP_RESPONSE_TO_LIPID-DEGs.pdf", width=4, height=8)
ggplot(long_df %>% inner_join(geneSymbol_1month_CB_qval05FC0), aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()



## CB_1year
CB_1year <- tpm_all_sample_tidy_GOBP_RESPONSE_TO_LIPID %>%
  filter(new_ID %in% c("1year_CB_WT_R1", "1year_CB_WT_R2", "1year_CB_WT_R3", 
                       "1year_CB_KO_R1", "1year_CB_KO_R2", "1year_CB_KO_R3")) 

CB_1year$new_ID_ordered <- factor(CB_1year$new_ID, levels = c("1year_CB_WT_R1", "1year_CB_WT_R2", "1year_CB_WT_R3", 
                                                  "1year_CB_KO_R1", "1year_CB_KO_R2", "1year_CB_KO_R3"))
### Order the rows based on the new ID order and TPM values
CB_1year <- CB_1year %>%
  arrange(new_ID_ordered, desc(TPM))
### Pivot the data to a long format suitable for ggplot
long_df <- CB_1year %>%
  pivot_longer(cols = TPM, names_to = "Condition", values_to = "Expression")
### Raw uncluseted heatmap
pdf("output/tpm/heatmap_CB_1year-GOBP_RESPONSE_TO_LIPID.pdf", width=4, height=8)
ggplot(long_df, aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()
### Raw uncluseted heatmap with only DEGs
#### import DEGs

geneSymbol_1year_CB_qval05FCmore0 <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCmore0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1year_CB_qval05FCless0 <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCless0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1year_CB_qval05FC0 = geneSymbol_1year_CB_qval05FCmore0 %>%
  bind_rows(geneSymbol_1year_CB_qval05FCless0) %>%
  rename("external_gene_name" = "V1")

pdf("output/tpm/heatmap_CB_1year-GOBP_RESPONSE_TO_LIPID-DEGs.pdf", width=4, height=8)
ggplot(long_df %>% inner_join(geneSymbol_1year_CB_qval05FC0), aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()



## CX_1month
CX_1month <- tpm_all_sample_tidy_GOBP_RESPONSE_TO_LIPID %>%
  filter(new_ID %in% c("1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
                       "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3")) 

CX_1month$new_ID_ordered <- factor(CX_1month$new_ID, levels = c("1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
                                                  "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3"))
### Order the rows based on the new ID order and TPM values
CX_1month <- CX_1month %>%
  arrange(new_ID_ordered, desc(TPM))
### Pivot the data to a long format suitable for ggplot
long_df <- CX_1month %>%
  pivot_longer(cols = TPM, names_to = "Condition", values_to = "Expression")
### Raw uncluseted heatmap
pdf("output/tpm/heatmap_CX_1month-GOBP_RESPONSE_TO_LIPID.pdf", width=4, height=8)
ggplot(long_df, aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()
### Raw uncluseted heatmap with only DEGs
#### import DEGs

geneSymbol_1month_CX_qval05FCmore0 <- read.csv("output/GO/geneSymbol_1month_CX_qval05FCmore0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1month_CX_qval05FCless0 <- read.csv("output/GO/geneSymbol_1month_CX_qval05FCless0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1month_CX_qval05FC0 = geneSymbol_1month_CX_qval05FCmore0 %>%
  bind_rows(geneSymbol_1month_CX_qval05FCless0) %>%
  rename("external_gene_name" = "V1")

pdf("output/tpm/heatmap_CX_1month-GOBP_RESPONSE_TO_LIPID-DEGs.pdf", width=4, height=2)
ggplot(long_df %>% inner_join(geneSymbol_1month_CX_qval05FC0), aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()


## CX_1year
CX_1year <- tpm_all_sample_tidy_GOBP_RESPONSE_TO_LIPID %>%
  filter(new_ID %in% c("1year_CX_WT_R1", "1year_CX_WT_R2", "1year_CX_WT_R3", 
                       "1year_CX_KO_R1", "1year_CX_KO_R2", "1year_CX_KO_R3")) 

CX_1year$new_ID_ordered <- factor(CX_1year$new_ID, levels = c("1year_CX_WT_R1", "1year_CX_WT_R2", "1year_CX_WT_R3", 
                                                  "1year_CX_KO_R1", "1year_CX_KO_R2", "1year_CX_KO_R3"))
### Order the rows based on the new ID order and TPM values
CX_1year <- CX_1year %>%
  arrange(new_ID_ordered, desc(TPM))
### Pivot the data to a long format suitable for ggplot
long_df <- CX_1year %>%
  pivot_longer(cols = TPM, names_to = "Condition", values_to = "Expression")
### Raw uncluseted heatmap
pdf("output/tpm/heatmap_CX_1year-GOBP_RESPONSE_TO_LIPID.pdf", width=4, height=8)
ggplot(long_df, aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()
### Raw uncluseted heatmap with only DEGs
#### import DEGs

geneSymbol_1year_CX_qval05FCmore0 <- read.csv("output/GO/geneSymbol_1year_CX_qval05FCmore0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1year_CX_qval05FCless0 <- read.csv("output/GO/geneSymbol_1year_CX_qval05FCless0.txt", header=FALSE, stringsAsFactors=FALSE)
geneSymbol_1year_CX_qval05FC0 = geneSymbol_1year_CX_qval05FCmore0 %>%
  bind_rows(geneSymbol_1year_CX_qval05FCless0) %>%
  rename("external_gene_name" = "V1")

pdf("output/tpm/heatmap_CX_1year-GOBP_RESPONSE_TO_LIPID-DEGs.pdf", width=4, height=2)
ggplot(long_df %>% inner_join(geneSymbol_1year_CX_qval05FC0), aes(x = new_ID_ordered, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = -1))
dev.off()







## clustered heatmap _ v1 with all replicate

desired_samples <- tpm_all_sample_tidy_GOBP_RESPONSE_TO_LIPID %>% 
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3")) 

### Create a matrix for the heatmap
heatmap_data <- dcast(desired_samples, external_gene_name ~ new_ID, value.var = "TPM")

### Convert to matrix and transpose
heatmap_matrix <- as.matrix(heatmap_data[, -1]) # remove the gene names column
rownames(heatmap_matrix) <- heatmap_data$external_gene_name

### Remove all rows which have 0 or 1 across all samples
heatmap_matrix <- heatmap_matrix[!(rowSums(heatmap_matrix != 0) <= 1), ]

### plot


pdf("output/tpm/heatmap_CB_1month-GOBP_RESPONSE_TO_LIPID-cluster.pdf", width=5, height=8)
pheatmap(heatmap_matrix, 
         scale = "none", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         color = viridis::viridis(100), 
         border_color = NA)
dev.off()


## clustered heatmap _ v2 with median replicates



desired_samples <- tpm_all_sample_tidy_GOBP_RESPONSE_TO_LIPID %>% 
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3")) %>%
  mutate(new_ID_grouped = sub("_R[0-9]+$", "", new_ID)) %>%
  group_by(new_ID_grouped, external_gene_name) %>%
  mutate(tpm_median = median(TPM)) %>%
  dplyr::select(-TPM, -new_ID) %>%
  unique()


### Create a matrix for the heatmap
heatmap_data <- dcast(desired_samples, external_gene_name ~ new_ID_grouped, value.var = "tpm_median")

### Convert to matrix and transpose
heatmap_matrix <- as.matrix(heatmap_data[, -1]) # remove the gene names column
rownames(heatmap_matrix) <- heatmap_data$external_gene_name

### Remove all rows which have 0 or 1 across all samples
heatmap_matrix <- heatmap_matrix[!(rowSums(heatmap_matrix != 0) <= 1), ]

### plot


pdf("output/tpm/heatmap_CB_1month-GOBP_RESPONSE_TO_LIPID-cluster_median.pdf", width=4, height=8)
pheatmap(heatmap_matrix, 
         scale = "none", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         color = viridis::viridis(100), 
         border_color = NA)
dev.off()





## heatmap _ v2 with median replicates - geneLists 20240124
tpm_all_sample_tidy <- read.csv("output/tpm/tpm_all_sample_tidy.txt", sep = "\t", header = TRUE)
### import gene list
geneLists = read_csv("output/tpm/geneListsCurated_heatmap_20240124_geneSymbol.txt") # all leading edge genes
geneLists = read_csv("output/tpm/geneListsCurated_heatmap_20240124_geneSymbol_CB1monthGenes.txt") # CB 1 month leading genes only

### combine with expression
tpm_all_sample_tidy_geneLists <- geneLists %>%
  left_join(tpm_all_sample_tidy, by = "external_gene_name")

desired_samples <- tpm_all_sample_tidy_geneLists %>%
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3",
                       "1year_CB_WT_R1", "1year_CB_WT_R2", "1year_CB_WT_R3", 
                       "1year_CB_KO_R1", "1year_CB_KO_R2", "1year_CB_KO_R3",
                       "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
                       "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3",
                       "1year_CX_WT_R1", "1year_CX_WT_R2", "1year_CX_WT_R3", 
                       "1year_CX_KO_R1", "1year_CX_KO_R2", "1year_CX_KO_R3")) %>%
  mutate(new_ID_grouped = sub("_R[0-9]+$", "", new_ID)) %>%
  group_by(new_ID_grouped, external_gene_name) %>%
  summarise(tpm_median = median(log2(tpm + 1))) %>%
  ungroup()

# Pivot the data to a long format suitable for ggplot
long_df <- desired_samples %>%
  pivot_longer(cols = tpm_median, names_to = "Condition", values_to = "Expression")


long_df$new_ID_grouped <-
  factor(long_df$new_ID_grouped,
         c("1month_CX_WT", "1month_CX_KO", "1year_CX_WT", "1year_CX_KO",
           "1month_CB_WT", "1month_CB_KO", "1year_CB_WT", "1year_CB_KO"))

pdf("output/tpm/heatmap-geneListsCurated-median.pdf", width=5, height=5)
pdf("output/tpm/heatmap-geneListsCurated_CB1monthGenes-median.pdf", width=3, height=8)

ggplot(long_df, aes(x = new_ID_grouped, y = reorder(external_gene_name, Expression), fill = Expression) )+
  geom_tile() +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=3, name="log2(tpm+1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        legend.position = "right")
dev.off()





## Z - score plot ######################
# Pivot the data to a long format suitable for ggplot

desired_samples <- tpm_all_sample_tidy_geneLists %>%
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3",
                       "1year_CB_WT_R1", "1year_CB_WT_R2", "1year_CB_WT_R3", 
                       "1year_CB_KO_R1", "1year_CB_KO_R2", "1year_CB_KO_R3",
                       "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
                       "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3",
                       "1year_CX_WT_R1", "1year_CX_WT_R2", "1year_CX_WT_R3", 
                       "1year_CX_KO_R1", "1year_CX_KO_R2", "1year_CX_KO_R3")) %>%
  mutate(new_ID_grouped = sub("_R[0-9]+$", "", new_ID)) %>%
  group_by(new_ID_grouped, external_gene_name) %>%
  summarise(tpm_median = median(log2(tpm + 1))) %>%
  ungroup() %>%
  mutate(z_score = (tpm_median - mean(tpm_median)) / sd(tpm_median))


long_df <- desired_samples %>%
  pivot_longer(cols = z_score, names_to = "Condition", values_to = "Expression")


long_df$new_ID_grouped <-
  factor(long_df$new_ID_grouped,
         c("1month_CX_WT", "1month_CX_KO", "1year_CX_WT", "1year_CX_KO",
           "1month_CB_WT", "1month_CB_KO", "1year_CB_WT", "1year_CB_KO"))

pdf("output/tpm/heatmap-geneListsCurated_CB1monthGenes-median_Zscore.pdf", width=4, height=8)

ggplot(long_df, aes(x = new_ID_grouped, y = reorder(external_gene_name, Expression), fill = Expression) )+
  geom_tile() +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="red", midpoint=0, name="z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        legend.position = "right")
dev.off()

######################################## STAT boxplot

  

library("ggpubr")

pdf("output/tpm/boxplot-geneListsCurated-median.pdf", width=6, height=4)
ggboxplot(long_df, x = "new_ID_grouped", y = "Expression",
  add.params = list(size = 1, alpha = 0.5),
      fill = "new_ID_grouped", palette = c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"),  add = "jitter") + theme_classic() +
  stat_compare_means(comparisons = list( c("1month_CX_WT", "1month_CX_KO"), c("1year_CX_WT", "1year_CX_KO"), c("1month_CB_WT", "1month_CB_KO"), c("1year_CB_WT", "1year_CB_KO") )) + # Add pairwise comparisons p-value
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5) )
dev.off()


pdf("output/tpm/boxplot_nojitter-geneListsCurated-median.pdf", width=6, height=4)
ggboxplot(long_df, x = "new_ID_grouped", y = "Expression",
  add.params = list(size = 1, alpha = 0.5),
      fill = "new_ID_grouped", palette = c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")) + theme_classic() +
  stat_compare_means(comparisons = list( c("1month_CX_WT", "1month_CX_KO"), c("1year_CX_WT", "1year_CX_KO"), c("1month_CB_WT", "1month_CB_KO"), c("1year_CB_WT", "1year_CB_KO") )) + # Add pairwise comparisons p-value
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5) )
dev.off()

########################################


## heatmap _ from top5 GSEA terms

tpm_all_sample_tidy <- read.csv("output/tpm/tpm_all_sample_tidy.txt", sep = "\t", header = TRUE)
### import gene list
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1year_corr.txt") # top 5
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1month_Top5_corr.txt") # top 5
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1month_Top10_corr.txt") # top 5
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1month_Top20_corr.txt") # top 5

## from top5 GSEA terms TPM

### combine with expression
tpm_all_sample_tidy_geneLists <- geneLists %>%
  left_join(tpm_all_sample_tidy, by = "external_gene_name")

desired_samples <- tpm_all_sample_tidy_geneLists %>%
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3",
                       "1year_CB_WT_R1", "1year_CB_WT_R2", "1year_CB_WT_R3", 
                       "1year_CB_KO_R1", "1year_CB_KO_R2", "1year_CB_KO_R3",
                       "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
                       "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3",
                       "1year_CX_WT_R1", "1year_CX_WT_R2", "1year_CX_WT_R3", 
                       "1year_CX_KO_R1", "1year_CX_KO_R2", "1year_CX_KO_R3")) %>%
  mutate(new_ID_grouped = sub("_R[0-9]+$", "", new_ID)) %>%
  group_by(new_ID_grouped, external_gene_name) %>%
  summarise(tpm_median = median(log2(tpm + 1))) %>%
  ungroup() 


# Pivot the data to a long format suitable for ggplot
long_df <- desired_samples %>%
  pivot_longer(cols = tpm_median, names_to = "Condition", values_to = "Expression")


long_df$new_ID_grouped <-
  factor(long_df$new_ID_grouped,
         c("1month_CX_WT", "1month_CX_KO", "1year_CX_WT", "1year_CX_KO",
           "1month_CB_WT", "1month_CB_KO", "1year_CB_WT", "1year_CB_KO"))

pdf("output/tpm/heatmap-leadingEdgeGenes_CB_1month_Top20_corr-median.pdf", width=5, height=5)

ggplot(long_df, aes(x = new_ID_grouped, y = reorder(external_gene_name, Expression), fill = Expression) )+
  geom_tile() +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=2.5, name="log2(tpm+1)") +  # mid point4 CB1year / CB 1 month Top5/10 2.5 / Top10   / Top20
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        legend.position = "right")
dev.off()

########################## boxplot stat

library("ggpubr")

pdf("output/tpm/boxplot-leadingEdgeGenes_CB_1year_corr-median.pdf", width=6, height=4)
ggboxplot(long_df, x = "new_ID_grouped", y = "Expression",
  add.params = list(size = 1, alpha = 0.5),
      fill = "new_ID_grouped", palette = c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"),  add = "jitter") + theme_classic() +
  stat_compare_means(comparisons = list( c("1month_CX_WT", "1month_CX_KO"), c("1year_CX_WT", "1year_CX_KO"), c("1month_CB_WT", "1month_CB_KO"), c("1year_CB_WT", "1year_CB_KO") )) + # Add pairwise comparisons p-value
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5) )
dev.off()


##########################

## Zscore 

## from top5 GSEA terms TPM
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1year_corr.txt") # top 5


geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1month_Top5_corr.txt") # top 5
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1month_Top10_corr.txt") # top 10
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1month_Top20_corr.txt") # top 20


### combine with expression
tpm_all_sample_tidy_geneLists <- geneLists %>%
  left_join(tpm_all_sample_tidy, by = "external_gene_name")

desired_samples <- tpm_all_sample_tidy_geneLists %>%
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3",
                       "1year_CB_WT_R1", "1year_CB_WT_R2", "1year_CB_WT_R3", 
                       "1year_CB_KO_R1", "1year_CB_KO_R2", "1year_CB_KO_R3",
                       "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
                       "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3",
                       "1year_CX_WT_R1", "1year_CX_WT_R2", "1year_CX_WT_R3", 
                       "1year_CX_KO_R1", "1year_CX_KO_R2", "1year_CX_KO_R3")) %>%
  mutate(new_ID_grouped = sub("_R[0-9]+$", "", new_ID)) %>%
  group_by(new_ID_grouped, external_gene_name) %>%
  summarise(tpm_median = median(log2(tpm + 1))) %>%
  ungroup()%>%
  mutate(z_score = (tpm_median - mean(tpm_median)) / sd(tpm_median))




# Pivot the data to a long format suitable for ggplot
long_df <- desired_samples %>%
  pivot_longer(cols = z_score, names_to = "Condition", values_to = "Expression")


long_df$new_ID_grouped <-
  factor(long_df$new_ID_grouped,
         c("1month_CX_WT", "1month_CX_KO", "1year_CX_WT", "1year_CX_KO",
           "1month_CB_WT", "1month_CB_KO", "1year_CB_WT", "1year_CB_KO"))

# pdf("output/tpm/heatmap-leadingEdgeGenes_CB_1year_corr-Zscore.pdf", width=5, height=5)
pdf("output/tpm/heatmap-leadingEdgeGenes_CB_1month_Top20_corr-Zscore.pdf", width=5, height=5)

ggplot(long_df, aes(x = new_ID_grouped, y = reorder(external_gene_name, Expression), fill = Expression) )+
  geom_tile() +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        legend.position = "right")
dev.off()






## zscore on tpm 

## from top5 GSEA terms TPM

### combine with expression
tpm_all_sample_tidy_geneLists <- geneLists %>%
  left_join(tpm_all_sample_tidy, by = "external_gene_name")

desired_samples <- tpm_all_sample_tidy_geneLists %>%
  filter(new_ID %in% c("1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
                       "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3",
                       "1year_CB_WT_R1", "1year_CB_WT_R2", "1year_CB_WT_R3", 
                       "1year_CB_KO_R1", "1year_CB_KO_R2", "1year_CB_KO_R3",
                       "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
                       "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3",
                       "1year_CX_WT_R1", "1year_CX_WT_R2", "1year_CX_WT_R3", 
                       "1year_CX_KO_R1", "1year_CX_KO_R2", "1year_CX_KO_R3")) %>%
  mutate(new_ID_grouped = sub("_R[0-9]+$", "", new_ID)) %>%
  group_by(new_ID_grouped, external_gene_name) %>%
  summarise(tpm_median = median(tpm + 1)) %>%
  ungroup()%>%
  mutate(z_score = (tpm_median - mean(tpm_median)) / sd(tpm_median))




# Pivot the data to a long format suitable for ggplot
long_df <- desired_samples %>%
  pivot_longer(cols = z_score, names_to = "Condition", values_to = "Expression")


long_df$new_ID_grouped <-
  factor(long_df$new_ID_grouped,
         c("1month_CX_WT", "1month_CX_KO", "1year_CX_WT", "1year_CX_KO",
           "1month_CB_WT", "1month_CB_KO", "1year_CB_WT", "1year_CB_KO"))

pdf("output/tpm/heatmap-leadingEdgeGenes_CB_1year_corr-tpmZscore.pdf", width=5, height=5)

ggplot(long_df, aes(x = new_ID_grouped, y = reorder(external_gene_name, Expression), fill = Expression) )+
  geom_tile() +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=2, name="Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        legend.position = "right")
dev.off()





## heatmap _ from top5 GSEA terms and using log2fc instead of tpm!

CX_1month <- as_tibble(read.csv("output/deseq2_corr/filtered_CX_1month__KO_vs_WT.txt", sep = "\t", header = TRUE) ) %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj) %>%
  add_column(group = "CX_1month") %>%
  drop_na()
CX_1year <- as_tibble(read.csv("output/deseq2_corr/filtered_CX_1year__KO_vs_WT.txt", sep = "\t", header = TRUE) ) %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj) %>%
  add_column(group = "CX_1year")  %>%
  drop_na()
CB_1month <- as_tibble(read.csv("output/deseq2_corr/filtered_CB_1month__KO_vs_WT.txt", sep = "\t", header = TRUE) ) %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj) %>%
  add_column(group = "CB_1month")  %>%
  drop_na()
CB_1year <- as_tibble(read.csv("output/deseq2_corr/filtered_CB_1year__KO_vs_WT.txt", sep = "\t", header = TRUE) ) %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj) %>%
  add_column(group = "CB_1year")  %>%
  drop_na()

DEGs_tidy = CX_1month %>%
  bind_rows(CX_1year) %>%
  bind_rows(CB_1month) %>%
  bind_rows(CB_1year) %>%
  dplyr::rename("external_gene_name" = "GeneSymbol") %>%
  mutate(log2FoldChange = ifelse(padj > 0.05, 0, log2FoldChange)) # Set log2FC to 0 for non-significant genes


# Import gene list
geneLists <- read_tsv("output/tpm/GSEA_leadingEdgeGenes_CB_1year_corr.txt")

# Combine with expression
DEGs_tidy_geneLists <- geneLists %>%
  left_join(DEGs_tidy, by = "external_gene_name")

# Pivot the data for heatmap plotting
long_heatmap_data <- DEGs_tidy_geneLists %>%
  pivot_longer(cols = starts_with("log2FoldChange"), names_to = "Condition", values_to = "Expression") %>%
  mutate(group = factor(group, levels = c("CX_1month", "CX_1year", "CB_1month", "CB_1year")))
  


# Plot the heatmap
pdf("output/deseq2_corr/heatmap-leadingEdgeGenes_CB_1year_corr-DEGslog2FC.pdf", width = 10, height = 8)
ggplot(long_heatmap_data, aes(x = group, y = reorder(external_gene_name, Expression), fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, name = "log2(FoldChange)", 
                       limit = c(-2, 2), space = "Lab", na.value = "white") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_blank(), # Remove y-axis labels
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    legend.position = "right"
  ) +
  labs(x = "", y = "Gene")
dev.off()



```


--> `*top5` correspond to the selection of the top 5 leading edge genes from each Term of the signficant GSEA at 1year CB; excell file for selection in `Fig_V2/GSEA_leadingEdgeGenes_CB_1year_corr.xlsx` and output gene list as `Fig_V2/GSEA_leadingEdgeGenes_CB_1year_corr.txt` cp into `tpm/GSEA_leadingEdgeGenes_CB_1year_corr.txt` in HPC

--> To strenghten difference I tried tpm solely or zscore on tpm solely but  it is bad, not visible because value are too different







## bigwig coverage files

Generate bigwig coverage files from the bam.


Let's generate **TPM coverage**:

```bash
conda activate deeptools
# run time-per-time:
sbatch scripts/TPM_bw_1.sh # 11761814 xxx
sbatch scripts/TPM_bw_2.sh # 11761885 xxx
```



# DEGs with deseq2


**IMPORTANT NOTE: Here it is advisable to REMOVE all genes from chromosome X and Y BEFORE doing the DEGs analysis (X chromosome re-activation occurs in some samples, notably these with more cell passage; in our case, the HET and KO)**


## PCA and clustering with deseq2 mm10

In R; Import all counts and combine into one matrix, deseq2 dataframe (dds)
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
### 1month data
samples <- c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", "S_CX_KO1" ,"S_CX_KO2" ,"S_CX_KO3" ,"S_CX_WT1", "S_CX_WT2" ,"S_CX_WT3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/bam/1month_rnaseqyz072420/")) %>%
    rename(!!sample := starts_with("output/bam/1month_rnaseqyz072420/"))
}

# Merge all dataframe into a single one
counts_all_1 <- reduce(sample_data, full_join, by = "Geneid")

### 1year data
samples <- c("171HetCB", "174MTCB", "175HetCB" ,"177MTCB" ,"474WTCB", "171HetCX", "174MTCX" ,"175HetCX", "177MTCX", "474WTCX")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/bam/1year_AL1804271_R2_new_analysis/")) %>%
    rename(!!sample := starts_with("output/bam/1year_AL1804271_R2_new_analysis/"))
}

counts_all_2 <- reduce(sample_data, full_join, by = "Geneid")

### combine 1month and 1year
counts_all = counts_all_1 %>%
  left_join(counts_all_2)

### rename sample
fileName <- tibble(
  ID = c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", 
         "S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3", 
         "171HetCB", "174MTCB", "175HetCB", "177MTCB", "474WTCB", 
         "171HetCX", "174MTCX", "175HetCX", "177MTCX", "474WTCX"),
  genotype = c("KO", "KO", "KO", "WT", "WT", "WT", 
               "KO", "KO", "KO", "WT", "WT", "WT", 
               "WT", "KO", "WT", "KO", "WT", 
               "WT", "KO", "WT", "KO", "WT"),
  time = c(rep("1month", 12), rep("1year", 10)),
  tissue = c(rep("CB", 6), rep("CX", 6), rep("CB", 5), rep("CX", 5)),
  replicate = c(rep("R1", 2), rep("R2", 2), rep("R3", 2), rep("R1", 2), rep("R2", 2), rep("R3", 2), "R1", "R1", "R2", "R2", "R3", "R1", "R1", "R2", "R2", "R3"), 
  new_ID = c("1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
             "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
) 

id_to_new_id <- fileName[, c("ID", "new_ID")]
new_column_names <- setNames(id_to_new_id$new_ID, id_to_new_id$ID)
counts_all_rename <- counts_all
names(counts_all_rename)[-1] <- new_column_names[match(names(counts_all)[-1], names(new_column_names))]



# ALL sample together
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
counts_all_matrix = make_matrix(select(counts_all_rename, -Geneid), pull(counts_all_rename, Geneid)) 

## Create colData file that describe all our samples
### Including replicate


coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID)

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -new_ID), pull(coldata_raw, new_ID))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
### Genotype only
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype )
#### Tissue and genotype and time
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype + tissue + time)



# Data normalization
vsd <- vst(dds, blind=TRUE)


# Vizualization for quality metrics
## Heatmap of the sample-to-sample distances
### vsd 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$time, vsd$tissue, vsd$genotype, vsd$replicate, sep="_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("output/deseq2_corr/heatmap_cluster_vsd.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

## PCA
### vsd 
pdf("output/deseq2_corr/PCA_vsd.pdf", width=7, height=5)

pcaData <- plotPCA(vsd, intgroup=c("time", "tissue", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=tissue, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()
dev.off()

pdf("output/deseq2_corr/PCA_vsd_1.pdf", width=7, height=5)

pcaData <- plotPCA(vsd, intgroup=c("time", "tissue", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()
dev.off()



# Separated by time
counts_all_rename_1month = counts_all_rename %>%
  dplyr::select("Geneid", "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3")

counts_all_rename_1year = counts_all_rename %>%
  dplyr::select("Geneid", "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
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
counts_all_matrix = make_matrix(select(counts_all_rename_1year, -Geneid), pull(counts_all_rename_1year, Geneid)) 

## Create colData file that describe all our samples
### Including replicate


coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(time == "1month")

coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(time == "1year")

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -new_ID), pull(coldata_raw, new_ID))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
#### Tissue and genotype and time
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype + tissue)



# Data normalization
vsd <- vst(dds, blind=TRUE)


# Vizualization for quality metrics
## Heatmap of the sample-to-sample distances
### vsd 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$tissue, vsd$genotype, vsd$replicate, sep="_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("output/deseq2_corr/heatmap_cluster_vsd_1month.pdf", width=5, height=6)
pdf("output/deseq2_corr/heatmap_cluster_vsd_1year.pdf", width=5, height=6)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

## PCA
### vsd 
pdf("output/deseq2_corr/PCA_vsd_1month.pdf", width=7, height=5)
pdf("output/deseq2_corr/PCA_vsd_1year.pdf", width=7, height=5)

pcaData <- plotPCA(vsd, intgroup=c("tissue", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=tissue, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()
dev.off()



# Separated by tissue
counts_all_rename_CB = counts_all_rename %>%
  dplyr::select("Geneid", "1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3")

counts_all_rename_CX = counts_all_rename %>%
  dplyr::select("Geneid", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
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
counts_all_matrix = make_matrix(select(counts_all_rename_CX, -Geneid), pull(counts_all_rename_CX, Geneid)) 

## Create colData file that describe all our samples
### Including replicate


coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(tissue == "CB")

coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(tissue == "CX")

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -new_ID), pull(coldata_raw, new_ID))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
#### Tissue and genotype and time
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype + time)



# Data normalization
vsd <- vst(dds, blind=TRUE)


# Vizualization for quality metrics
## Heatmap of the sample-to-sample distances
### vsd 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$tissue, vsd$genotype, vsd$replicate, sep="_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("output/deseq2_corr/heatmap_cluster_vsd_CB.pdf", width=5, height=6)
pdf("output/deseq2_corr/heatmap_cluster_vsd_CX.pdf", width=5, height=6)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

## PCA
### vsd 

pdf("output/deseq2_corr/PCA_vsd_CB.pdf", width=7, height=5)
pdf("output/deseq2_corr/PCA_vsd_CX.pdf", width=7, height=5)

pcaData <- plotPCA(vsd, intgroup=c("time", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()
dev.off()
```



## 'one-by-one' comparison
Comparison tisse/time per tisse/time:
- CB_1month _ WT vs KO
- CB_1year _ WT vs KO
- CX_1month _ WT vs KO
- CX_1year _ WT vs KO

*NOTE: For the relax1 featureCounts version, number are with decimal as I used the `fraction` option, so I need to round them to make it work with deseq2 at the dds step.*

###  CB_1month _ WT vs KO

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
### 1month data
samples <- c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3")


## Make a loop for importing all featurecounts data and keep only ID and count column

sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/bam/1month_rnaseqyz072420/")) %>%
    dplyr::rename(!!sample := starts_with("output/bam/1month_rnaseqyz072420/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")


fileName <- tibble(
  ID = c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", 
         "S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3", 
         "171HetCB", "174MTCB", "175HetCB", "177MTCB", "474WTCB", 
         "171HetCX", "174MTCX", "175HetCX", "177MTCX", "474WTCX"),
  genotype = c("KO", "KO", "KO", "WT", "WT", "WT", 
               "KO", "KO", "KO", "WT", "WT", "WT", 
               "WT", "KO", "WT", "KO", "WT", 
               "WT", "KO", "WT", "KO", "WT"),
  time = c(rep("1month", 12), rep("1year", 10)),
  tissue = c(rep("CB", 6), rep("CX", 6), rep("CB", 5), rep("CX", 5)),
  replicate = c(rep("R1", 2), rep("R2", 2), rep("R3", 2), rep("R1", 2), rep("R2", 2), rep("R3", 2), "R1", "R1", "R2", "R2", "R3", "R1", "R1", "R2", "R2", "R3"), 
  new_ID = c("1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
             "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
) 

id_to_new_id <- fileName[, c("ID", "new_ID")]
new_column_names <- setNames(id_to_new_id$new_ID, id_to_new_id$ID)
counts_all_rename <- counts_all
names(counts_all_rename)[-1] <- new_column_names[match(names(counts_all)[-1], names(new_column_names))]

counts_all = counts_all_rename

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(tissue == "CB",
         time == "1month")


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -new_ID), pull(coldata_raw, new_ID))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

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
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

###### save output
res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2_corr/filtered_CB_1month__KO_vs_WT.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
######

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################

keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, '#165CAA',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Red',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Red'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == '#165CAA'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

### DEG_CB_1month
highlight_genes <- c("Fabp5", "Fabp1", "Acsl5", "Gstm2", "Dcn", "Hmox1", "Nol3", "Snx14") # 1month_CB
highlight_genes <- c("Fabp5", "Dcn", "Snx14") # 1month_CB_signifOnly

pdf("output/deseq2_corr/plotVolcano_DEG_CB_1month.pdf", width=8, height=8)  
pdf("output/deseq2_corr/plotVolcano_DEG_CB_1month_signifOnly.pdf", width=8, height=8)  


EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = ' ',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 6,
  labSize = 11,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  ylim(0,110)

dev.off()


# count genes


upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)



# Save as gene list for GO analysis:
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_corr/upregulated_q05FC05_DEG_CB_1month.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_corr/downregulated_q05FC05_DEG_CB_1month.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```




###  CB_1year _ WT vs KO

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
### 1year data
samples <- c("171HetCB", "174MTCB", "175HetCB" ,"177MTCB" ,"474WTCB")

## Make a loop for importing all featurecounts data and keep only ID and count column

sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/bam/1year_AL1804271_R2_new_analysis/")) %>%
    dplyr::rename(!!sample := starts_with("output/bam/1year_AL1804271_R2_new_analysis/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")


fileName <- tibble(
  ID = c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", 
         "S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3", 
         "171HetCB", "174MTCB", "175HetCB", "177MTCB", "474WTCB", 
         "171HetCX", "174MTCX", "175HetCX", "177MTCX", "474WTCX"),
  genotype = c("KO", "KO", "KO", "WT", "WT", "WT", 
               "KO", "KO", "KO", "WT", "WT", "WT", 
               "WT", "KO", "WT", "KO", "WT", 
               "WT", "KO", "WT", "KO", "WT"),
  time = c(rep("1month", 12), rep("1year", 10)),
  tissue = c(rep("CB", 6), rep("CX", 6), rep("CB", 5), rep("CX", 5)),
  replicate = c(rep("R1", 2), rep("R2", 2), rep("R3", 2), rep("R1", 2), rep("R2", 2), rep("R3", 2), "R1", "R1", "R2", "R2", "R3", "R1", "R1", "R2", "R2", "R3"), 
  new_ID = c("1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
             "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
) 

id_to_new_id <- fileName[, c("ID", "new_ID")]
new_column_names <- setNames(id_to_new_id$new_ID, id_to_new_id$ID)
counts_all_rename <- counts_all
names(counts_all_rename)[-1] <- new_column_names[match(names(counts_all)[-1], names(new_column_names))]

counts_all = counts_all_rename

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(tissue == "CB",
         time == "1year")


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -new_ID), pull(coldata_raw, new_ID))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

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
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

###### save output
res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2_corr/filtered_CB_1year__KO_vs_WT.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
######

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################

keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, '#165CAA',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Red',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Red'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == '#165CAA'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

### DEG_CB_1year
highlight_genes <- c("Acacb", "Acsm5", "Cp","Lyz2", "C4b", "Cd68", "Trem2", "Apoe","Gfap","Casp3","Lcn2","Snx14","Rgs8" ,"Pcp2","Car8","Calb1") # 1year_CB 



pdf("output/deseq2_corr/plotVolcano_DEG_CB_1year.pdf", width=8, height=8)  


EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = ' ',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 6,
  labSize = 11,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  ylim(0,110)

dev.off()


# count genes


upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)



# Save as gene list for GO analysis:
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_corr/upregulated_q05FC05_DEG_CB_1year.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_corr/downregulated_q05FC05_DEG_CB_1year.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```


###  CX_1month _ WT vs KO

I had to update conda env deseq2 due to lfcShrink error (discuss [here](https://www.biostars.org/p/9559740/)): 
```bash
Error in optimHess(par = init, fn = nbinomFn, gr = nbinomGr, x = x, y = y,  : 
  non-finite value supplied by optim; 
```

- clone conda env with `conda create --name deseq2V2 --clone deseq2`
- install devtools in R: `install.packages("devtools")`
- install last version of apeglm: `devtools::install_github("azhu513/apeglm")`



Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
### 1month data
samples <- c("S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3")


## Make a loop for importing all featurecounts data and keep only ID and count column

sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/bam/1month_rnaseqyz072420/")) %>%
    dplyr::rename(!!sample := starts_with("output/bam/1month_rnaseqyz072420/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")


fileName <- tibble(
  ID = c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", 
         "S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3", 
         "171HetCB", "174MTCB", "175HetCB", "177MTCB", "474WTCB", 
         "171HetCX", "174MTCX", "175HetCX", "177MTCX", "474WTCX"),
  genotype = c("KO", "KO", "KO", "WT", "WT", "WT", 
               "KO", "KO", "KO", "WT", "WT", "WT", 
               "WT", "KO", "WT", "KO", "WT", 
               "WT", "KO", "WT", "KO", "WT"),
  time = c(rep("1month", 12), rep("1year", 10)),
  tissue = c(rep("CB", 6), rep("CX", 6), rep("CB", 5), rep("CX", 5)),
  replicate = c(rep("R1", 2), rep("R2", 2), rep("R3", 2), rep("R1", 2), rep("R2", 2), rep("R3", 2), "R1", "R1", "R2", "R2", "R3", "R1", "R1", "R2", "R2", "R3"), 
  new_ID = c("1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
             "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
) 

id_to_new_id <- fileName[, c("ID", "new_ID")]
new_column_names <- setNames(id_to_new_id$new_ID, id_to_new_id$ID)
counts_all_rename <- counts_all
names(counts_all_rename)[-1] <- new_column_names[match(names(counts_all)[-1], names(new_column_names))]

counts_all = counts_all_rename

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(tissue == "CX",
         time == "1month")


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -new_ID), pull(coldata_raw, new_ID))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

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
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

###### save output
res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2_corr/filtered_CX_1month__KO_vs_WT.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
######

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################

keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, '#165CAA',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Red',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Red'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == '#165CAA'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

### DEG_CX_1month
highlight_genes <- c("Snx14" ) # 1month_CX 

pdf("output/deseq2_corr/plotVolcano_DEG_CX_1month.pdf", width=8, height=8)  


EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = ' ',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 6,
  labSize = 11,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  ylim(0,110)

dev.off()


# count genes


upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)



# Save as gene list for GO analysis:
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_corr/upregulated_q05FC05_DEG_CX_1month.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_corr/downregulated_q05FC05_DEG_CX_1month.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```





###  CX_1year _ WT vs KO

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
### 1year data
samples <- c("171HetCX", "174MTCX", "175HetCX" ,"177MTCX" ,"474WTCX")

## Make a loop for importing all featurecounts data and keep only ID and count column

sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/bam/1year_AL1804271_R2_new_analysis/")) %>%
    dplyr::rename(!!sample := starts_with("output/bam/1year_AL1804271_R2_new_analysis/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")


fileName <- tibble(
  ID = c("S_CB_KO1", "S_CB_KO2", "S_CB_KO3", "S_CB_WT1", "S_CB_WT2", "S_CB_WT3", 
         "S_CX_KO1", "S_CX_KO2", "S_CX_KO3", "S_CX_WT1", "S_CX_WT2", "S_CX_WT3", 
         "171HetCB", "174MTCB", "175HetCB", "177MTCB", "474WTCB", 
         "171HetCX", "174MTCX", "175HetCX", "177MTCX", "474WTCX"),
  genotype = c("KO", "KO", "KO", "WT", "WT", "WT", 
               "KO", "KO", "KO", "WT", "WT", "WT", 
               "WT", "KO", "WT", "KO", "WT", 
               "WT", "KO", "WT", "KO", "WT"),
  time = c(rep("1month", 12), rep("1year", 10)),
  tissue = c(rep("CB", 6), rep("CX", 6), rep("CB", 5), rep("CX", 5)),
  replicate = c(rep("R1", 2), rep("R2", 2), rep("R3", 2), rep("R1", 2), rep("R2", 2), rep("R3", 2), "R1", "R1", "R2", "R2", "R3", "R1", "R1", "R2", "R2", "R3"), 
  new_ID = c("1month_CB_KO_R1", "1month_CB_KO_R2", "1month_CB_KO_R3", 
             "1month_CB_WT_R1", "1month_CB_WT_R2", "1month_CB_WT_R3", 
             "1month_CX_KO_R1", "1month_CX_KO_R2", "1month_CX_KO_R3", 
             "1month_CX_WT_R1", "1month_CX_WT_R2", "1month_CX_WT_R3", 
             "1year_CB_WT_R1", "1year_CB_KO_R1", "1year_CB_WT_R2", 
             "1year_CB_KO_R2", "1year_CB_WT_R3", "1year_CX_WT_R1", 
             "1year_CX_KO_R1", "1year_CX_WT_R2", "1year_CX_KO_R2", 
             "1year_CX_WT_R3")
) 

id_to_new_id <- fileName[, c("ID", "new_ID")]
new_column_names <- setNames(id_to_new_id$new_ID, id_to_new_id$ID)
counts_all_rename <- counts_all
names(counts_all_rename)[-1] <- new_column_names[match(names(counts_all)[-1], names(new_column_names))]

counts_all = counts_all_rename

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(fileName) %>%
  dplyr::select(-ID) %>%
  filter(tissue == "CX",
         time == "1year")


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -new_ID), pull(coldata_raw, new_ID))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

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
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

###### save output
res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2_corr/filtered_CX_1year__KO_vs_WT.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
######

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################

keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, '#165CAA',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Red',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Red'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == '#165CAA'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

### DEG_CX_1year
highlight_genes <- c("Sv2c", "Doc2b","Snx14") # 1year_CX


pdf("output/deseq2_corr/plotVolcano_DEG_CX_1year.pdf", width=8, height=8)  


EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = ' ',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 6,
  labSize = 11,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  ylim(0,110)

dev.off()


# count genes


upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)



# Save as gene list for GO analysis:
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_corr/upregulated_q05FC05_DEG_CX_1year.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_corr/downregulated_q05FC05_DEG_CX_1year.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```













# volcano plot labelling genes _ V1 with output from Novogene & Genewiser

Let's label some genes in the volcano plots
- import RNAseq Novogen data
- label genes with enhancedvolcano


Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("biomaRt")
library("readxl")

# import DEGs output
DEG_CB_1month = read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 1) %>%
  dplyr::select(gene_name, log2FoldChange,padj) %>%
  unique()
DEG_CB_1year= read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 2) %>%
  dplyr::select(Gene.name, log2FoldChange,padj) %>%
  unique() %>%
  dplyr::rename("gene_name" = "Gene.name")
DEG_CX_1month = read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 3) %>%
  dplyr::select(gene_name, log2FoldChange,padj) %>%
  unique()
DEG_CX_1year= read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 4)%>%
  dplyr::select(Gene.name, log2FoldChange,padj) %>%
  unique() %>%
  dplyr::rename("gene_name" = "Gene.name")


# re-calculate pvalue adj
DEG_CB_1month = read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 1) %>%
  dplyr::select(gene_name, log2FoldChange,pvalue) %>%
  unique() %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method
DEG_CB_1year= read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 2) %>%
  dplyr::select(Gene.name, log2FoldChange,pvalue) %>%
  unique() %>%
  dplyr::rename("gene_name" = "Gene.name")%>%
  mutate(padj = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method
DEG_CX_1month = read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 3) %>%
  dplyr::select(gene_name, log2FoldChange,pvalue) %>%
  unique()%>%
  mutate(padj = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method
DEG_CX_1year= read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 4)%>%
  dplyr::select(Gene.name, log2FoldChange,pvalue) %>%
  unique() %>%
  dplyr::rename("gene_name" = "Gene.name")%>%
  mutate(padj = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method





# re-calculate pvalue adj and compare with previous   ########################
DEG_CB_1month = read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 1) %>%
  dplyr::select(gene_name, log2FoldChange,pvalue,padj) %>%
  unique() %>%
  mutate(padj_new = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method
DEG_CB_1year= read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 2) %>%
  dplyr::select(Gene.name, log2FoldChange,pvalue,padj) %>%
  unique() %>%
  dplyr::rename("gene_name" = "Gene.name")%>%
  mutate(padj_new = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method
DEG_CX_1month = read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 3) %>%
  dplyr::select(gene_name, log2FoldChange,pvalue,padj) %>%
  unique()%>%
  mutate(padj_new = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method
DEG_CX_1year= read_excel('output/deseq2/RNAseq of CB&CX_1mon&1yr.xlsx', sheet = 4)%>%
  dplyr::select(Gene.name, log2FoldChange,pvalue,padj) %>%
  unique() %>%
  dplyr::rename("gene_name" = "Gene.name")%>%
  mutate(padj_new = p.adjust(pvalue, method = "BH")) # Benjamini & Hochberg method


write.table(DEG_CB_1month, sep = "\t", quote = FALSE, row.names=FALSE, file="output/deseq2/DEG_CB_1month.txt")
write.table(DEG_CB_1year, sep = "\t", quote = FALSE, row.names=FALSE, file="output/deseq2/DEG_CB_1year.txt")
write.table(DEG_CX_1month, sep = "\t", quote = FALSE, row.names=FALSE, file="output/deseq2/DEG_CX_1month.txt")
write.table(DEG_CX_1year, sep = "\t", quote = FALSE, row.names=FALSE, file="output/deseq2/DEG_CX_1year.txt")

########################################################################


### DEG_CB_1month
highlight_genes <- c("Fabp5", "Fabp1", "Acsl5", "Gstm2", "Dcn", "Hmox1", "Nol3", "Snx14") # 1month_CB
highlight_genes <- c("Fabp5", "Dcn", "Snx14") # 1month_CB_signifOnly

highlight_genes <- c("Acacb", "Acsm5", "Cp","Lyz2", "C4b", "Cd68", "Trem2", "Apoe","Gfap","Casp3","Lcn2","Snx14" ) # 1year_CB 
highlight_genes <- c("Snx14" ) # 1month_CX 
highlight_genes <- c("Sv2c", "Doc2b","Snx14") # 1year_CX


# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  DEG_CB_1month$log2FoldChange < -0.5 & DEG_CB_1month$padj < 5e-2, 'Sky Blue',
    ifelse(DEG_CB_1month$log2FoldChange > 0.5 & DEG_CB_1month$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Red'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

pdf("output/deseq2/plotVolcano_DEG_CB_1month.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_DEG_CB_1month_signifOnly.pdf", width=8, height=8)  

pdf("output/deseq2/plotVolcano_DEG_CX_1month.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_DEG_CX_1year.pdf", width=8, height=8)  

pdf("output/deseq2/plotVolcano_DEG_CB_1year.pdf", width=8, height=8)  


pdf("output/deseq2/plotVolcano_DEG_CB_1month_pvalueadj_signifOnly.pdf", width=8, height=8)  


pdf("output/deseq2/plotVolcano_DEG_CB_1year_pvalueadj.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_DEG_CX_1month_pvalueadj.pdf", width=8, height=8)  
pdf("output/deseq2/plotVolcano_DEG_CX_1year_pvalueadj.pdf", width=8, height=8)  

EnhancedVolcano(DEG_CX_1year,
  lab = DEG_CX_1year$gene_name,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 1yearCB',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 6,
  labSize = 11,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  ylim(0,100)

dev.off()


# count genes


upregulated_genes <- sum(DEG_CX_1year$log2FoldChange > 0.5 & DEG_CX_1year$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(DEG_CX_1year$log2FoldChange < -0.5 & DEG_CX_1year$padj < 5e-2, na.rm = TRUE)



# Save as gene list for GO analysis:
upregulated <- DEG_CX_1year[!is.na(DEG_CX_1year$log2FoldChange) & !is.na(DEG_CX_1year$padj) & DEG_CX_1year$log2FoldChange > 0.5 & DEG_CX_1year$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- DEG_CX_1year[!is.na(DEG_CX_1year$log2FoldChange) & !is.na(DEG_CX_1year$padj) & DEG_CX_1year$log2FoldChange < -0.5 & DEG_CX_1year$padj < 5e-2, ]
#### Save
write.table(upregulated$gene_name, file = "output/deseq2/upregulated_q05FC05_DEG_CX_1year.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$gene_name, file = "output/deseq2/downregulated_q05FC05_DEG_DEG_CX_1year.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

```




--> Needed to re-calculate the padj as some genes were having an NA value... Including Lcn2...
----> Plot corrected as `*pvalueadj.pdf`








# Upload files to GEO


Let's upload additional files to the current FTP server (see email 2/2/2024)

- upload the processed file to my personal space (`uploads/thomasroule@orcid_A787EGG4`) and transfer files
- send email to GEO to confirm (here no need to send another metadata file!)

*NOTE: For count files, I had to copy all of them into a new folder `001__RNAseq/count/` as I have to upload the entire content to geo and do not want to add the summary files.*


```bash
module load lftp


# connect to ftp
lftp -u geoftp,inAlwokhodAbnib5 ftp-private.ncbi.nlm.nih.gov # geoftp = username; inAlwokhodAbnib5 = pwd
cd uploads/thomasroule@orcid_A787EGG4

mirror -R ../001__RNAseq/output/bigwig/
mirror -R ../001__RNAseq/count/
mirror -R ../001__RNAseq/meta/ # contain the meta_additionalProcessedFiles.xlsx file
```






