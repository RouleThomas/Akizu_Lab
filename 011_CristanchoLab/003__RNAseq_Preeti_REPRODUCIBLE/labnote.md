--> This folder is the **clean + reproducible** version of: `001_CristanchoLab/002__RNAseq_Preeti` (original exploration)

# Data Overview

- **Data type**: Bulk RNA-seq
- **Library**: Directional mRNA (poly(A) enrichment)
- **Platform**: Illumina NovaSeq X Plus (paired-end, 150 bp)
- **Factors**:
  - **Cell type**: ReN, PSC
  - **Condition**: Normoxia (24h), Hypoxia (24h)
  - **Replicates**: 4 biological replicates per group (total n = 16)



# Data processing

## Sample Renaming

| Key | Original sample name | New sample name |
| --: | -------------------- | --------------- |
|  A1 | ReN_Nor_24h_Exp3     | ReN_Norm_Rep1   |
|  A2 | ReN_Nor_24h_Exp4     | ReN_Norm_Rep2   |
|  A3 | ReN_Nor_24h_Exp5     | ReN_Norm_Rep3   |
|  A4 | ReN_Nor_24h_Exp6     | ReN_Norm_Rep4   |
|  A5 | ReN_Hyp_24h_Exp3     | ReN_Hypo_Rep1   |
|  A6 | ReN_Hyp_24h_Exp4     | ReN_Hypo_Rep2   |
|  A7 | ReN_Hyp_24h_Exp5     | ReN_Hypo_Rep3   |
|  A8 | ReN_Hyp_24h_Exp6     | ReN_Hypo_Rep4   |
|  B1 | hiPSCs_Nor_24h_Exp3  | PSC_Norm_Rep1   |
|  B2 | hiPSCs_Nor_24h_Exp4  | PSC_Norm_Rep2   |
|  B3 | hiPSCs_Nor_24h_Exp5  | PSC_Norm_Rep3   |
|  B4 | hiPSCs_Nor_24h_Exp6  | PSC_Norm_Rep4   |
|  B5 | hiPSCs_Hyp_24h_Exp3  | PSC_Hypo_Rep1   |
|  B6 | hiPSCs_Hyp_24h_Exp4  | PSC_Hypo_Rep2   |
|  B7 | hiPSCs_Hyp_24h_Exp5  | PSC_Hypo_Rep3   |
|  B8 | hiPSCs_Hyp_24h_Exp6  | PSC_Hypo_Rep4   |



## Read pre-processing & Quality filtering

Reads were trimmed using **fastp (VERSION)** with default adapter detection and quality filtering parameters.

The following command was run for each (paired-end) sample:
```bash
fastp \
  -i <sample>_1.fq.gz \
  -I <sample>_2.fq.gz \
  -o <sample>_1.fq.gz \
  -O <sample>_2.fq.gz \
  -j <sample>.json \
  -h <sample>.html
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/fastp.sh`


## Read mapping

Reads were aligned to the human reference genome (**hg38**, *GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta*) using **STAR** (v2.7.3a).  
Resulting BAM files were sorted by coordinate and indexed using **SAMtools** (v1.16.1).

The following command was run for each (paired-end) sample:
```bash
STAR \
  --genomeDir <genome_index> \
  --readFilesIn <sample>_1.fq.gz <sample>_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate

samtools index <sample>_Aligned.sortedByCoord.out.bam
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/STAR_mapping_fastp.sh`


## Gene-level quantification

Gene-level read counts were generated using **featureCounts** (v2.0.1) with the **GENCODE v47** gene annotation.

The following command was run for each sample:
```bash
featureCounts -p -C -O -M --fraction -s 2
	-a <annotation.gtf> \
	-o <sample>.txt <sample>_Aligned.sortedByCoord.out.bam
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/featurecounts_multi.sh`



## Isoform-level quantification


Isoform-level read counts were generated using **kallisto** (v0.44.0) with the **GENCODE v47** gene annotation.


The following command was run for each sample:
```bash
kallisto quant \
  -i <transcripts.idx> \
  -o <sample>_quant \
  -b 100 \
  -g <annotation.gtf> \
  --rf-stranded \
  --genomebam \
  --chromosomes GRCh38_chrom_sizes.tab \
  <sample>_1.fq.gz <sample>_2.fq.gz
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/kallisto_count_gtf.sh`




## Transcript abundance estimation (TPM/RPKM)

Gene-level TPM and RPKM values were computed from the featureCounts output using a custom R script `RPKM_TPM_featurecounts.R`. For each input featureCounts table, the script converts raw counts to TPM/RPKM using the feature length column (Length) and outputs two tab-separated files: `*_tpm.txt` and `*_rpkm.txt`.

The following command was run for each sample:
```bash
Rscript scripts/RPKM_TPM_featurecounts.R <featureCounts_output>.txt <output_prefix>
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/scripts/featurecounts_TPM.sh`



## Coverage track & QC

Genome-wide coverage tracks (.bigWig) were generated from coordinate-sorted BAM files using **bamCoverage** (v3.5.1). Coverage was normalized using BPM (bins per million mapped reads) with single–base resolution.

The following command was run for each sample:
```bash
bamCoverage \
  --bam <sample>_Aligned.sortedByCoord.out.bam \
  --outFileName <sample>.bw \
  --outFileFormat bigwig \
  --normalizeUsing BPM \
  --binSize 1
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/bigwigmerge_STAR_TPM_bw.sh `

For each biological condition, replicate bigWig files were aggregated by computing the median coverage signal across using **wiggletools** (v1.0.0). Median coverage tracks were first generated in bedGraph format, subsequently sorted using **bedtools** (v2.30.0), and finally converted back to bigWig format using **bedGraphToBigWig** (v4).

```bash
# Compute median coverage across replicates
wiggletools write_bg <condition>_median.bedGraph \
  median <rep1>.bw <rep2>.bw <rep3>.bw <rep4>.bw
# Sort bedGraph
bedtools sort -i <condition>_median.bedGraph > <condition>_median.sorted.bedGraph
# Convert bedGraph to bigWig
bedGraphToBigWig \
  <condition>_median.sorted.bedGraph \
  GRCh38_chrom_sizes.tab \
  <condition>_median.bw
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/bigwigmerge_STAR_TPM_bw.sh`






# Data analysis



## Gene expression (TPM)




## Differential Gene Expression 



## Alternative splicing 



## Functional enrichment 





