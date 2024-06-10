# Project

Study KCNC1 mutation (Kcnc1-Arg320His). In human HT KCNC1 missense mutation (dominant negative=mutated non functional prot compete with the functional prot.) associated with epilepsy and cerebellar ataxia (Cerebellum/CB not working well, but histology do not show anything --> scRNAseq may). Generated mice model with same mutation; also got phentoype (HM mice die very young, HT are used here). So KCNC1 is a potassium channel and if it is not working properly, the excitatory neurons will not be slow-down, too much active and lead to epilepsy. KCNC work in tetramer with other KCNC proteins.

A drug has been design to promote potassium channel activation and it rescue the phenotype in mice (currently, tested in human): paper [here](https://pubmed.ncbi.nlm.nih.gov/38266642/)

KCNC1 express day12 in mice neuron (important for fast-spike interneuron (inhibitory); express in Cortex/Ctx and CB. In CB expressed in unipolar brush cell (UBC), purkinje layer interneuron (PLI), Golgi, Molecular Layer Interneuron 2 (MLI2), Purkinje). Histology no marker for MLI2. --> [webtool](https://portal.nemoarchive.org/) for cell type in CB. This [paper](https://www.sciencedirect.com/science/article/pii/S0896627324002484) show that MLI2 where KCNC1 is express, inhibit MLI1 who inhibit Purkinje cell activity. Loss of Purkinje activtiy or Purkinje cells leads to ataxia. 

2 genotype (WT, KCNC1-mutant); 2 brain regions (Ctx, CB); 3 ages (p14,35,180 = 14 is when KCNC1 start being express in intenreunon; p35= 1st phenotype very mild; p180= degeneration)



# Objectives

- Check cell population changes in WT vs mutant; of Ctx and CB. Then DEGs analysis
- check where KCNC1 is express (compare with CB [paper](https://pubmed.ncbi.nlm.nih.gov/34616064/)). Check trajectory of expression too.
- consequence of interneuron not working properly; may be increased excitatory neurons activity (leading to epilepsy); check if more excitatory neurons, or more activtiy (check marker genes of activity maybe? Channel?)
- Check KCNC1 expression over time and within the different cell types. 
    - The weird thing is that the KCNC1 mutation disease is progressive. Not clear why. Maybe where KCNC1 expression pattern change, or the cell types appearition themself change (Like MLI2 appear later?)

--> Focus the analysis on the cell types expressing KCNC1



# Paper for cell type annotation

- CB: [paper](https://pubmed.ncbi.nlm.nih.gov/34616064/); cell types here should be found at p35-p180
- Ctx: XXX



# Data acess

To access data folow email
```
Service Request : Fileshare_Access_2801 
Requester : goldberge 
Fileshare : goldberg_lab_scb 
AccessType : ReadWrite 
Fileshare Path : \\ressmb05.research.chop.edu\goldberg_lab_scb 
```

Open files app/Connect to server; add Fileshare path and add my credentials

Then, right click on the server and open with terminal

--> Copy all files to `002*/005*`

```bash
cp -r * /scr1/users/roulet/Akizu_Lab/002_scRNAseq/005__Goldberg/input_raw/

```

--> ALL GOOD

File used in `/snRNAseq_Kcnc1_R320H/snRNAseq_Kcnc1_reorganized/*`

For the analysis I will follow the YAP1 scRNAseq in `002003`


## Counting with cellranger count

Within each folder in `/snRNAseq_Kcnc1_R320H/snRNAseq_Kcnc1_reorganized/*` I have two lanes L001 and L002 with I1/I2 and R1/R2 fastq. 








```bash 
conda activate scRNAseq
which cellranger

# Run count using mice genome
## p14 _ Kcnc1 _ CB and CX
sbatch scripts/cellranger_count_Kcnc1_p14_CB_Rep1.sh # 20177866 FAIL corrupted; 20178043 FAIL should have delete previous; 20203124 xxx
sbatch scripts/cellranger_count_Kcnc1_p14_CB_Rep2.sh # 20177913 FAIL corrupted; 20178032 FAIL should have delete previous; 20203125 xxx
sbatch scripts/cellranger_count_Kcnc1_p14_CB_Rep3.sh # 20177921 ok
sbatch scripts/cellranger_count_Kcnc1_p14_CX_Rep1.sh # 20177943 FAIL corrupted; 20178594 FAIL should have delete previous; 20203127 xxx
sbatch scripts/cellranger_count_Kcnc1_p14_CX_Rep2.sh # 20177952 ok
sbatch scripts/cellranger_count_Kcnc1_p14_CX_Rep3.sh # 20177959 ok


## p14 _ WT _ CB and CX
sbatch scripts/cellranger_count_WT_p14_CB_Rep1.sh # 20178099 ok
sbatch scripts/cellranger_count_WT_p14_CB_Rep2.sh # 20178116 ok
sbatch scripts/cellranger_count_WT_p14_CB_Rep3.sh # 20178127 ok
sbatch scripts/cellranger_count_WT_p14_CX_Rep1.sh # 20178154 ok 
sbatch scripts/cellranger_count_WT_p14_CX_Rep2.sh # 20178162 ok
sbatch scripts/cellranger_count_WT_p14_CX_Rep3.sh # 20178170 ok


## p35 _ Kcnc1 _ CB and CX
sbatch scripts/cellranger_count_Kcnc1_p35_CB_Rep1.sh # 20178221 fail fastq path; 20203342 xxx
sbatch scripts/cellranger_count_Kcnc1_p35_CB_Rep2.sh # 20178240 fail fastq path; 20203343 xxx
sbatch scripts/cellranger_count_Kcnc1_p35_CB_Rep3.sh # 20178250 fail fastq path; 20203344 xxx
sbatch scripts/cellranger_count_Kcnc1_p35_CX_Rep1.sh # 20178258 fail fastq path; 20203351 xxx 
sbatch scripts/cellranger_count_Kcnc1_p35_CX_Rep2.sh # 20178280 fail fastq path; 20203352 xxx
sbatch scripts/cellranger_count_Kcnc1_p35_CX_Rep3.sh # 20178344 fail fastq path; 20203353 xxx


## p35 _ WT _ CB and CX
sbatch scripts/cellranger_count_WT_p35_CB_Rep1.sh # 20178351 fail fastq path; 20203366 xxx
sbatch scripts/cellranger_count_WT_p35_CB_Rep2.sh # 20178354 fail fastq path; 20203367 xxx
sbatch scripts/cellranger_count_WT_p35_CB_Rep3.sh # 20178357 fail fastq path; 20203368 xxx
sbatch scripts/cellranger_count_WT_p35_CX_Rep1.sh # 20178372 ok 
sbatch scripts/cellranger_count_WT_p35_CX_Rep2.sh # 20178379 fail fastq path; 20203383 xxx
sbatch scripts/cellranger_count_WT_p35_CX_Rep3.sh # 20178386 fail fastq path; 20203384 xxx



## p180 _ Kcnc1 _ CB and CX
sbatch scripts/cellranger_count_Kcnc1_p180_CB_Rep1.sh # 20178463 fail fastq path; 20203418 xxx
sbatch scripts/cellranger_count_Kcnc1_p180_CB_Rep2.sh # 20178466 fail fastq path; 20203419 xxx
sbatch scripts/cellranger_count_Kcnc1_p180_CB_Rep3.sh # 20178471 fail fastq path; 20203420 xxx
sbatch scripts/cellranger_count_Kcnc1_p180_CX_Rep1.sh # 20178478 fail fastq path; 20203432 xxx
sbatch scripts/cellranger_count_Kcnc1_p180_CX_Rep2.sh # 20178481 fail fastq path; 20203437 xxx
sbatch scripts/cellranger_count_Kcnc1_p180_CX_Rep3.sh # 20178485 fail fastq path; 20203438 xxx


## p180 _ WT _ CB and CX
sbatch scripts/cellranger_count_WT_p180_CB_Rep1.sh # 20178494 fail fastq path; 20203461 xxx
sbatch scripts/cellranger_count_WT_p180_CB_Rep2.sh # 20178497 fail fastq path; 20203464 xxx
sbatch scripts/cellranger_count_WT_p180_CB_Rep3.sh # 20178501 fail fastq path; 20203465 xxx
sbatch scripts/cellranger_count_WT_p180_CX_Rep1.sh # 20178518 fail fastq path; 20203505 xxx
sbatch scripts/cellranger_count_WT_p180_CX_Rep2.sh # 20178525 fail fastq path; 20203508 xxx
sbatch scripts/cellranger_count_WT_p180_CX_Rep3.sh # 20178587 fail fastq path; 20203509 xxx
```


--> Bug for:
- `cellranger_count_Kcnc1_p14_CB_Rep1`: `Sequence and quality length mismatch: file: "/scr1/users/roulet/Akizu_Lab/002_scRNAseq/003__YAP1/input/24hgastruloidhumanUN/24hgastruloidhumanUN_S1_L001_R2_001.fastq.gz", line: 292599520`
    - Another file available in `8_31_23 core reanalyzed corrupted files`

-----> To avoid error, let's directly use the reanalyzed one for the concerned samples: 
- Kcnc1_p14_CB_Rep1
- Kcnc1_p14_CB_Rep2
- Kcnc1_p14_CX_Rep1
- WT_p35_CX_Rep1
























## RNA contamination and doublet detection
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object

```bash
srun --mem=500g --pty bash -l
conda deactivate # base environment needed
python3 scrublet.py [input_path] [output_path]
# Run doublet detection/scrublet sample per sample
python3 scripts/scrublet_doublets.py E7mousecontrolQCNOTfail/outs/filtered_feature_bc_matrix output/doublets/embryo_E7_control.tsv
python3 scripts/scrublet_doublets.py E7mousecYAPKO/outs/filtered_feature_bc_matrix output/doublets/embryo_E7_cYAPKO.tsv

python3 scripts/scrublet_doublets.py 24hgastruloidhumanUNQCNOTfail/outs/filtered_feature_bc_matrix output/doublets/humangastruloid_UNTREATED24hr.tsv
python3 scripts/scrublet_doublets.py 24hgastruloidhumanDASA/outs/filtered_feature_bc_matrix output/doublets/humangastruloid_DASATINIB24hr.tsv

```
Doublet detection score:
- humangastruloid_UNTREATED24hr: 42% (previous time: 0% doublet)
- humangastruloid_DASATINIB24hr: 0.2% (previous time: 34.2% doublet)
- embryo_E7_control: 7.2% (previous time: 0.1% doublet)
- embryo_E7_cYAPKO: 6.2% (previous time: 3.6%)
--> Successfully assigned doublet





# embryo E7 (second sample) analysis in Seurat



