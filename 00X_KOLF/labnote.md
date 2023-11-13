

# Resp to reviewers; check if KOLF CNVs are present in Pangenome.

Pangenome [paper](https://www.nature.com/articles/s41586-023-05896-x) and [data](https://github.com/human-pangenomics/hpp_pangenome_resources).

--> From the data, VCF decomposed seems to be the best to use.

- Download the data
- Transfer data to HPC
- in HPC add .gz at the end of both files (vcf and the vcf index) and gunzip
- load the vcf in IGV and IGV  create index


The data is analyzed on hg38 genome. However our variant are called on hg19. Coordinate as been corrected with this [webtool](https://genome.ucsc.edu/cgi-bin/hgLiftOver)



**Hg19 coordinates**
Chr3p14.2. Chr3:61242795-61528270 (chr3:61135039-61911274)
Chr18q22.1. Chr18:62115075-62264706
Chr3p13.  Chr3:72289657-72404542
Chr9q33. Chr9:119244942-119355584
Chr6p22. Chr6:15487649-15722102

**Hg38 coordinates**

Chr3p14.2. chr3:61257121-61542596
Chr18q22.1. chr18:64447840-64597471
Chr3p13. chr3:72240506-72355391
Chr9q33. chr9:116482663-116593305
Chr6p22. chr6:15487418-15721871





To check for **dupplication**:
- Extract sequence (of the duplicated CNV)
- In the pangenome VCF filter to keep all Dupplication event
- Retrieve all dupplicated sequence from the VCF
- BLAST extracted sequence to pangenome duplicated sequence


To check for **deletion**:
- In the pangenome VCF filter to keep all Deletion event
- Check on IGV region of interest for large deletion






# Check our variant with bcftools

## Deletion 
Check what variant are found within our genomic coordinates:

```bash
module load BCFtools/1.15.1-GCC-11.3.0
 
# compress vcf with bgzip to use bcftools
bgzip -c input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf > input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz

# generate index
bcftools index input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz

# Collect all event in each coordinates
## Chr3p14.2
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' -r chr3:61228469-61513944 input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz > output/bcftools/Chr3p14.2.txt
## Chr3p13
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' -r chr3:72338808-72453693 input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz > output/bcftools/Chr3p13.txt


```


To facilitate observation of Deletion events; lets generate a new VCF file containing only deletion

```bash
# Extract coordinates of the deletion event
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz | awk -F"\t" 'length($3) > length($4)' > output/bcftools/deletions.txt
# Filter the original vcf for these deletion coordinates
bcftools view -T output/bcftools/deletions.txt input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz -Oz -o output/bcftools/deletions.vcf.gz
# decompress the vcf for vizualization in IGV
cp output/bcftools/deletions.vcf.gz output/bcftools/IGV_deletions.vcf.gz
gunzip output/bcftools/IGV_deletions.vcf.gz

```

- *NOTE: To isolate deletion we compare the size of refernece versus alternate allele, if reference > alternate = Deletion event*

Conclusion:
- Chr3p13: no big deletion event detected
- Chr9q33: no big deletion event detected
- Chr6p22: no big deletion event detected





# Dupplication


- Extract FASTA sequence of these coordinates:
Chr3p14.2. chr3:61257121-61542596
Chr18q22.1. chr18:64447840-64597471
- Filter VCF to keep only dupplication
- Extract DNA sequence as FASTA format
- BLAST

Install Blast+ to perform blast:

```bash

conda create -n blast -c bioconda blast # command can be launch from anywhere (directory and node)

conda create -n blast blast==2.15.0


```

```bash
conda activate bowtie2 

# create bed of coordinate
nano output/fasta/Chr3p14.bed 
nano output/fasta/Chr18q22.bed 

# use getfasta to get fasta seq
bedtools getfasta -fi ../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed output/fasta/Chr3p14.bed > output/fasta/Chr3p14.fa
bedtools getfasta -fi ../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed output/fasta/Chr18q22.bed > output/fasta/Chr18q22.fa


# Extract all dupplicated sequence from the VCF
## Extract coordinates of the deletion event
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz | awk -F"\t" 'length($3) < length($4)' > output/bcftools/insertion.txt
## Filter the original vcf for these deletion coordinates
bcftools view -T output/bcftools/insertion.txt input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz -Oz -o output/bcftools/insertion.vcf.gz
## Generate a FASTA of all insertions sequence from the output output/bcftools/insertion.txt
awk 'BEGIN{OFS="\n"}{print ">"$1":"$2, $4}' output/bcftools/insertion.txt > output/bcftools/insertion.fasta
## clean the FASTA (for when there is A,T in a row)
awk '
BEGIN {FS="\n"; RS=">"; ORS=""}
NR > 1 {
    header = $1;
    gsub(/\r/, "", header);
    split($2, seqs, ",");
    for (i in seqs) {
        if (length(seqs[i]) > 0) {
            printf(">%s\n%s\n", header, seqs[i]);
        }
    }
}' output/bcftools/insertion.fasta | sed 's/>$//' > output/bcftools/insertion_clean.fasta
## Create a BLAST db
conda activate blast 
makeblastdb -in output/bcftools/insertion_clean.fasta -dbtype nucl -out output/blast/insertion_db
## BLAST

blastn -query output/fasta/Chr3p14.fa -db output/blast/insertion_db -out output/blast/Chr3p14_blast.txt
blastn -query output/fasta/Chr18q22.fa -db output/blast/insertion_db -out output/blast/Chr18q22_blast.txt

## Extract the Length infoprmation (lenght of signficant overlap)

grep "Length=" output/blast/Chr3p14_blast.txt | sed 's/.*Length=//' | sed 's/[^0-9]*//g' > output/blast/Chr3p14_blast_Length.txt
grep "Length=" output/blast/Chr18q22_blast.txt | sed 's/.*Length=//' | sed 's/[^0-9]*//g' > output/blast/Chr18q22_blast_Length.txt

```




















