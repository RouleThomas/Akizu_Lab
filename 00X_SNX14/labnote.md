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

# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCless0.5.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCmore0.5.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 50)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 50)

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
pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC0.pdf", width=12, height=7)
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_GO_BP_1month_CX_qval05FC0.pdf", width=12, height=6)
pdf("output/GO/enrichR_GO_BP_1year_CX_qval05FC0.pdf", width=12, height=6)

pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC1.pdf", width=12, height=9)
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC1.pdf", width=12, height=4)

pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC0.5.pdf", width=12, height=11)

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
write.table(gos, "output/GO/enrichR_GO_BP_1month_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CB_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1month_CX_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CX_qval05FC0.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_BP_1month_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CB_qval05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(gos, "output/GO/enrichR_GO_BP_1month_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)



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


# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCless1.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1month_CB_qval05FCmore1.txt", header=FALSE, stringsAsFactors=FALSE)
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

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCless1.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCmore1.txt", header=FALSE, stringsAsFactors=FALSE)
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


# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCless1.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCmore1.txt", header=FALSE, stringsAsFactors=FALSE)
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

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.pdf", width=12, height=4)
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.pdf", width=12, height=4)

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC1.pdf", width=12, height=4)



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


