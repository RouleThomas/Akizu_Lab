

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
Chr3p14.2. chr3:61228469-61513944 (chr3:61149366-61925600)
Chr18q22.1. chr18:59782308-59931939
Chr3p13. chr3:72338808-72453693
Chr9q33. chr9:122007220-122117862
Chr6p22. chr6:15487880-15722333



# Check our variant with bcftools


Check what variant are found within our genomic coordinates:

```bash
module load BCFtools/1.15.1-GCC-11.3.0
 
# compress vcf with bgzip to use bcftools
bgzip -c input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf > input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz

# Hg38 coordinates
bcftools view -r chr3:61228469-61513944 input/pangenome/minigraph-cactus/decomposed_test/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.vcf.gz



```



