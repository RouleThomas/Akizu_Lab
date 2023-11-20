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


# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCless0.5.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCmore0.5.txt", header=FALSE, stringsAsFactors=FALSE)
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
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC0.5.pdf", width=12, height=8)

pdf("output/GO/enrichR_GO_BP_1year_CX_qval05FC0.5.pdf", width=12, height=4)

pdf("output/GO/enrichR_GO_BP_1month_CB_qval05FC0.5_filtUp.pdf", width=12, height=11)
pdf("output/GO/enrichR_GO_BP_1year_CB_qval05FC0.5_filtUp.pdf", width=12, height=8)


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
write.table(gos, "output/GO/enrichR_GO_BP_1year_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_BP_1year_CX_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)



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



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1year_CX_qval05FCless0.5.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1year_CX_qval05FCmore0.5.txt", header=FALSE, stringsAsFactors=FALSE)
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



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1year_CX_qval05FCless0.5.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1year_CX_qval05FCmore0.5.txt", header=FALSE, stringsAsFactors=FALSE)
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

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCless0.5.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/GO/geneSymbol_1year_CB_qval05FCmore0.5.txt", header=FALSE, stringsAsFactors=FALSE)
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

pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1month_CB_qval05FC0.5.pdf", width=12, height=4)
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.5.pdf", width=12, height=4)


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

write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_1year_CB_qval05FC0.5.txt", sep="\t", row.names=FALSE, quote=FALSE)

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
lfc_vector <- CB_1year$combined_score           ################ CHANGE !!!!!!!!!!!!!!!! ######################
names(lfc_vector) <- CB_1year$gene_name           ################ CHANGE !!!!!!!!!!!!!!!! ######################
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


