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


# XXX






