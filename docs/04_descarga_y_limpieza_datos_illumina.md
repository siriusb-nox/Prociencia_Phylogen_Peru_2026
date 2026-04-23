# 4. Descarga y limpieza de datos Illumina

[← Volver al README](../README.md)

## 4.1 Descarga de lecturas crudas desde NCBI SRA

```bash
chmod +x download_sra_fastq.sh
./download_sra_fastq.sh accessions.txt ~/Documents/Phylo_work_2026/raw_reads 4
```

## 4.2 Renombrar archivos FASTQ con nombres de especies y códigos de acceso

```bash
#!/usr/bin/env bash
set -euo pipefail

MAPFILE="accession_to_species.tsv"
declare -A species_map

while IFS=$'\t' read -r acc species; do
    [[ -z "${acc}" ]] && continue
    species_map["$acc"]="$species"
done < "$MAPFILE"

for f in *_*.fastq.gz; do
    [[ -e "$f" ]] || continue
    acc="${f%%_*}"
    rest="${f#*_}"
    pair="${rest%%.*}"

    if [[ -n "${species_map[$acc]:-}" ]]; then
        new="${species_map[$acc]}_${acc}_${pair}.fastq.gz"
        mv -- "$f" "$new"
    fi
done
```

## 4.3 Recorte de lecturas con `fastp`

```bash
while read F; do fastp -i ${F}_1.fastq.gz -I ${F}_2.fastq.gz -o ../trimmed/${F}_R1_fastp.fq -O ../trimmed/${F}_R2_fastp.fq -q 30 -w 8 -R ${F}.report -h ../trimmed/${F}.webpage; done < namelist.temp
```
