# B. Ensamblaje de captura dirigida y recuperación de secuencias

[← Volver al README](../README.md)

### Ensamblaje de loci con `HybPiper`

HybPiper es un pipeline bioinformático diseñado para recuperar genes a partir de lecturas de secuenciación obtenidas mediante target enrichment o Hyb-Seq. Su función principal es asignar las lecturas a genes de referencia, ensamblarlas por locus y extraer secuencias listas para análisis filogenómicos. Entre sus salidas más útiles están las secuencias codificantes recuperadas, estadísticas de éxito de captura por muestra y por gen, y archivos que luego pueden alinearse e incorporarse a inferencias filogenéticas. 


Activa primero el entorno conda de HybPiper:

```bash
conda activate hybpiper
```

Después ejecuta el bucle de ensamblaje:

```bash
while read F; do hybpiper assemble --readfiles ${F}_R1_fastp.fq ${F}_R2_fastp.fq --targetfile_dna /media/sf_VM_phylogenomics_Workshop_2026/data/HybPiper/RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta --cpu 2 --cov_cutoff 8 --bwa --timeout_assemble 80 --timeout_extract_contigs 100 --hybpiper_output ~/Documents/Phylo_work_2026/HybPiper --force_overwrite --prefix $F; done < namelist.temp
```

### 5.2 Resumen de recuperación de loci

Tras el ensamblaje, resume las estadísticas de recuperación con:

```bash
hybpiper stats --targetfile_dna RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta gene namelist.temp
```

Este paso ayuda a evaluar el éxito del enriquecimiento y del ensamblaje entre muestras y loci.

### 5.3 Recuperación de secuencias ensambladas

Recupera las secuencias de ADN para todas las muestras y loci con:

```bash
hybpiper retrieve_sequences dna --targetfile_dna RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta --sample_names namelist.temp
```

Esto produce archivos FASTA por gen para los alineamientos posteriores.
