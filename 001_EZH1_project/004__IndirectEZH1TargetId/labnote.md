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
- H3K27me3: XXX
- EZH2: XXX



