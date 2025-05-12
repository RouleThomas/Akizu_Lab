
# Projects

Collaboration with Joan from Fox Chase Cancer.

Simulate known mutation signature (SBS, DBS..) to the human reference genome/exome, generate mutated peptides and assess their binding affinity to all known MHC-I variants (ie. HLA variants)



# Backgrounds



MHC-I, can be produced by three genes HLA-A/B/C each individual has 6 alleles  Produce 6 different MHC-I complex per ind.
Each allele are prone to polymorphism! Notably HLA-A  Many different possibilities of sequence combination of HLAs. (>3k alleles variants)
Cancer mutations  abnormal proteins  broken down into peptides. MHC-I presents these peptides (neoantigens) on the cell surface for recognition and killing by CD8+ T cell 
MHC-I peptide binding recognition is variable; some people can present certain peptides/neoantigen that others cannot
Different mutagens (e.g. tobacco, UV, chemo), and cancer type, create distinct mutation signatures


--> Interesting docs:
- [Mutation signature](https://medium.com/@hylke.donker/mutational-signatures-explained-1dc435b2d7b7)

For SBS-96:
    6 mutation types
    4 possible 5’ flanking bases (A/C/G/T)
    4 possible 3’ flanking bases (A/C/G/T)
--> In total, there are 96 = 4 x 6 x 4 singlet classes when sorted by three letter — trinucleotide — motifs.



# Pipeline


- Collect ALL mutation signatures
    COSMIC Mutation Signature Database (include Alexandrov et al 2018)
    Kucab et al 2019; 41 environmental agents (6 DBS, 8 ID)
    Pich et al 2019 (therapies induced mutation)
- Simulate a 1000 random mutation signature profile for each known signature on human reference exome
- Confirm profile by generating mutation profile plot
- Collect mutated peptides for each signature 
- Assess binding affinity of MHC-I to these mutated peptide; for each HLA/MHC-I known variants (NetMHCPan3.0 ; PHBR score (Marty et al 2017))



# Tool


To simulate the mutation signature to genome/exome:
- [MSA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04450-8) --> More recent, test this one first!
- [SigProfilerSimulator](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03772-3)










