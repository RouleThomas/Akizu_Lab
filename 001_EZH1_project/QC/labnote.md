# QC metrics on current CutRuns

- File in Google: `001__EZH1_project/metrics`
- Importing into HPC: `input/QC_Metrics_CutRuns.xlsx`


```bash
srun --mem=5g --pty bash -l
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("readxl")
library("ggpubr")

# import QC file
QC_file <- read_excel("input/QC_Metrics_CutRuns.xlsx", sheet = 1)%>%
    filter(Note == "NA") 

QC_file$QC = factor(QC_file$QC, c("HIGH", "MEDIUM", "LOW"))


# All IP except IGG
QC_file_filt = QC_file %>%
    filter(IP != "IGG",
    QC != "MEDIUM")

## insert size
pdf("output/ggboxplot_all.pdf", width=7, height=7)
ggboxplot(QC_file_filt, x = "QC", y = "Fastp_InsertSizePeak",
  add.params = list(size = 1, alpha = 0.5),
      fill = "QC", palette = c("darkgrey","darkgrey","darkgrey"),  add = "jitter") + theme_classic() +
  stat_compare_means(comparisons = list( c("LOW", "HIGH"), c("MEDIUM", "HIGH") ))  # Add pairwise comparisons p-value
dev.off()


## Fastp_sizeDistributionPattern

pdf("output/SizeDistributionPattern_all.pdf", width=7, height=7)
QC_file_filt %>%
  ggplot(., aes(x = QC, fill = Fastp_sizeDistributionPattern)) +
  geom_bar(position = "fill") + # 'fill' stacks the bars and scales them to 1
  scale_y_continuous(labels = scales::percent_format()) + # Convert y-axis to percentage
  labs(x = "QC Category", y = "Percentage (%)", fill = "Size Distribution Pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("output/SizeDistributionPattern_all_ABsep.pdf", width=7, height=7)
QC_file_filt %>%
  ggplot(., aes(x = QC, fill = Fastp_sizeDistributionPattern)) +
  geom_bar(position = "fill") + # 'fill' stacks the bars and scales them to 1
  scale_y_continuous(labels = scales::percent_format()) + # Convert y-axis to percentage
  labs(x = "QC Category", y = "Percentage (%)", fill = "Size Distribution Pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~IP)
dev.off()


# only H3K27me3 IP
# All IP except IGG
QC_file_filt = QC_file %>%
    filter(IP == "H3K27me3")


pdf("output/ggboxplot_H3K27me3.pdf", width=7, height=7)
ggboxplot(QC_file_filt, x = "QC", y = "Fastp_InsertSizePeak",
  add.params = list(size = 1, alpha = 0.5),
      fill = "QC", palette = c("darkgrey","darkgrey","darkgrey"),  add = "jitter") + theme_classic() +
  stat_compare_means(comparisons = list( c("LOW", "HIGH"), c("MEDIUM", "HIGH") ))  # Add pairwise comparisons p-value
dev.off()




```