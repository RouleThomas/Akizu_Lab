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


# All histone marks
QC_file_filt = QC_file %>%
    filter(IP %in% c("H3K27me3", "H3K4me3", "H3K27ac", "H3K27me1"),
    QC != "MEDIUM",
    Exp_ID != "003_CutRun")


## Novogene_BioanalyzerPattern

pdf("output/BioanalyzerPattern_histone.pdf", width=7, height=7)
QC_file_filt %>%
  ggplot(., aes(x = QC, fill = Novogene_BioanalyzerPattern)) +
  geom_bar(position = "fill") + # 'fill' stacks the bars and scales them to 1
  scale_y_continuous(labels = scales::percent_format()) + # Convert y-axis to percentage
  labs(x = "QC Category", y = "Percentage (%)", fill = "Size Distribution Pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## Fastp_sizeDistributionPattern

pdf("output/Fastp_sizeDistributionPattern_histone.pdf", width=7, height=7)
QC_file_filt %>%
  ggplot(., aes(x = QC, fill = Fastp_sizeDistributionPattern)) +
  geom_bar(position = "fill") + # 'fill' stacks the bars and scales them to 1
  scale_y_continuous(labels = scales::percent_format()) + # Convert y-axis to percentage
  labs(x = "QC Category", y = "Percentage (%)", fill = "Size Distribution Pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




# PRC2 subunits
QC_file_filt = QC_file %>%
    filter(IP %in% c("EZH1cs", "EZH1pt", "EZH2", "SUZ12", "HA"),
    QC != "MEDIUM",
    Exp_ID != "003_CutRun")


## Novogene_BioanalyzerPattern

pdf("output/BioanalyzerPattern_PRC2.pdf", width=7, height=7)
QC_file_filt %>%
  ggplot(., aes(x = QC, fill = Novogene_BioanalyzerPattern)) +
  geom_bar(position = "fill") + # 'fill' stacks the bars and scales them to 1
  scale_y_continuous(labels = scales::percent_format()) + # Convert y-axis to percentage
  labs(x = "QC Category", y = "Percentage (%)", fill = "Size Distribution Pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## Fastp_sizeDistributionPattern

pdf("output/Fastp_sizeDistributionPattern_PRC2.pdf", width=7, height=7)
QC_file_filt %>%
  ggplot(., aes(x = QC, fill = Fastp_sizeDistributionPattern)) +
  geom_bar(position = "fill") + # 'fill' stacks the bars and scales them to 1
  scale_y_continuous(labels = scales::percent_format()) + # Convert y-axis to percentage
  labs(x = "QC Category", y = "Percentage (%)", fill = "Size Distribution Pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()








```