# 5. Ensamblaje de captura dirigida y recuperación de secuencias

[← Volver al README](../README.md)

## 5.1 Ensamblaje de loci con `HybPiper`

```bash
conda activate hybpiper
while read F; do hybpiper assemble --readfiles ${F}_R1_fastp.fq ${F}_R2_fastp.fq --targetfile_dna /media/sf_VM_phylogenomics_Workshop_2026/data/HybPiper/RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta --cpu 2 --cov_cutoff 8 --bwa --timeout_assemble 80 --timeout_extract_contigs 100 --hybpiper_output ~/Documents/Phylo_work_2026/HybPiper --force_overwrite --prefix $F; done < namelist.temp
```

## 5.2 Resumen de recuperación de loci

```bash
hybpiper stats --targetfile_dna RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta gene namelist.temp
```

## 5.3 Recuperación de secuencias ensambladas

```bash
hybpiper retrieve_sequences dna --targetfile_dna RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta --sample_names namelist.temp
```
