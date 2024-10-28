# Project and goal

Part of gastrulation paper with Conchi Estaras.

--> See `002*/003*/gastrulation paper/Figure 4` for detail: We need to know the overlap of the 199 QSER1:YAP cobound peaks (identified in `008*/001*` and `008*/003*`) with H3K27ac, H3K4me1, H3K27me3 and H3K36me3 modification marks

Pipeline:
- Identify data from ENCODE
- Download ENCODE data (bigwig and peak file)
- Follow `Figure 4` ppt

# Data download from ENCODE

Prioritize most recent datasets, or all datasets from same lab:
- **Bing Ren**, UCSD_Project ROADMAP; [2013 2 Bio Rep; processed 2020](https://www.encodeproject.org/experiments/ENCSR928HYM/):
    - **H3K27me3**: *ENCFF395GVR* bigwig signal p-value; *ENCFF599KDF* bed
    - **H3K4me1**:  *ENCFF164XHJ* bigwig signal p-value; *ENCFF613QAB* bed
    - **H3K36me3**: *ENCFF483UMR* bigwig signal p-value; *ENCFF681CEO* bed
    - **H3K27ac**:  *ENCFF390JIZ* bigwig signal p-value; *ENCFF045CUG* bed
- **Bradley Bernstein**, Broad_Project ENCODE; 2013
    - **H3K27me3**: *ENCFF193PKI* bigwig signal p-value; *ENCFF305KNA* bed
    - **H3K4me1**:  *ENCFF706CHK* bigwig signal p-value; *ENCFF984DGO* bed
    - **H3K36me3**: *ENCFF985CVI* bigwig signal p-value; *ENCFF504KOV* bed
    - **H3K27ac**:  *ENCFF771GNB* bigwig signal p-value; *ENCFF317QGQ* bed


--> File transfer to HPC cluster at `output/ENCODE`

--> Let's check both.





