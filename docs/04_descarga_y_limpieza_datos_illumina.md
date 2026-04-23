# A. Descarga y limpieza de datos Illumina

[← Volver al README](../README.md)

En este paso, vamos a acceder datos de secuenciacion de illumina a partir del repositorio [NCBI](https://www.ncbi.nlm.nih.gov/). La unica manera de acceder a estos datos es a traves de linea de comando y con la herramienta 'fasterq-dump'.

Prepara un archivo de texto plano con un acceso por línea, por ejemplo `accessions.txt`:

```text
ERR5006156
ERR7622624
ERR7622662
SRR26876932
SRR26876810
SRR26876866
SRR6127594
SRR6127593
ERR7599286
ERR4180009
```

El comando para descargar un solo set de datos a partir de un código de accesión es:

```bash
fasterq-dump ACCESSION --split-files --threads 2 
```
Para bajar múltiples archivos con un solo comando, podemos ejecutar el siguiente script (guarda el siguiente texto como `download_sra_fastq.sh`):

```bash
#!/usr/bin/env bash
set -euo pipefail

ACC_FILE="${1:-accessions.txt}"
OUTDIR="${2:-$PWD/sra_fastq}"
THREADS="${3:-4}"

TMPDIR="${OUTDIR}/tmp"
SRA_DIR="${OUTDIR}/sra_cache"
FASTQ_DIR="${OUTDIR}/fastq"

mkdir -p "${TMPDIR}" "${SRA_DIR}" "${FASTQ_DIR}"

while IFS= read -r ACC || [[ -n "$ACC" ]]; do
    [[ -z "$ACC" ]] && continue

    echo "Processing ${ACC} ..."

    prefetch "${ACC}" --output-directory "${SRA_DIR}"

    fasterq-dump "${SRA_DIR}/${ACC}" \
        --split-files \
        --threads "${THREADS}" \
        --temp "${TMPDIR}" \
        --outdir "${FASTQ_DIR}"

    gzip -f "${FASTQ_DIR}/${ACC}"*.fastq

    echo "Finished ${ACC}"
done < "${ACC_FILE}"

echo "All downloads completed."
echo "Compressed FASTQ files are in: ${FASTQ_DIR}"
```

Ejecútalo así:

```bash
chmod +x download_sra_fastq.sh
./download_sra_fastq.sh accessions.txt ~/Documents/Phylo_work_2026/raw_reads 4
```

### 4.2 Renombrar archivos FASTQ con nombres de especies y códigos de acceso

Supón que tienes un archivo tabulado de dos columnas llamado `accession_to_species.tsv`:

```text
ERR5006156	Apostasia_nuda
ERR7622624	Aa_maderoi
ERR7622662	Anacamptis_pyramidalis
SRR26876932	Cycnoches_amparoanum
SRR26876810	Coryanthes_hunteriana
SRR26876866	Vanilla_hartii
SRR6127594	Phragmipedium_lindleyanum
SRR6127593	Paphiopedilum_callosum
ERR7599286	Erythrorchis_cassythoides
ERR4180009	Dioscorea_caucasica
```

Puedes renombrar los FASTQ descargados con un script de Bash como este:

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

### Recorte de lecturas con `fastp`

El conjunto final de muestras usado en la práctica tras el recorte fue:

```text
Aa_maderoi_ERR7622624_R1_fastp.fq
Aa_maderoi_ERR7622624_R2_fastp.fq
Anacamptis_pyramidalis_ERR7622662_R1_fastp.fq
Anacamptis_pyramidalis_ERR7622662_R2_fastp.fq
Apostasia_nuda_ERR5006156_R1_fastp.fq
Apostasia_nuda_ERR5006156_R2_fastp.fq
Coryanthes_hunteriana_SRR26876810_R1_fastp.fq
Coryanthes_hunteriana_SRR26876810_R2_fastp.fq
Cycnoches_amparoanum_SRR26876932_R1_fastp.fq
Cycnoches_amparoanum_SRR26876932_R2_fastp.fq
Dioscorea_caucasica_ERR4180009_R1_fastp.fq
Dioscorea_caucasica_ERR4180009_R2_fastp.fq
Paphiopedilum_callosum_SRR6127593_R1_fastp.fq
Paphiopedilum_callosum_SRR6127593_R2_fastp.fq
Vanilla_hartii_SRR26876866_R1_fastp.fq
Vanilla_hartii_SRR26876866_R2_fastp.fq
```

El recorte se ejecutó con el siguiente bucle:

```bash
while read F; do fastp -i ${F}_1.fastq.gz -I ${F}_2.fastq.gz -o ../trimmed/${F}_R1_fastp.fq -O ../trimmed/${F}_R2_fastp.fq -q 30 -w 8 -R ${F}.report -h ../trimmed/${F}.webpage; done < namelist.temp
```

Este comando:

- lee los nombres de muestra desde `namelist.temp`
- recorta lecturas pareadas
- escribe las lecturas recortadas en `../trimmed/`
- mantiene lecturas con umbral de calidad `q=30`
- genera un informe HTML por muestra


Como output, el programa produce reportes de los resultados de filtrado y de calidad, antes y despues del filtrado, en formato html.
