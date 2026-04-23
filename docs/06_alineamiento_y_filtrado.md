# 6. Alineamiento por locus y filtrado

[← Volver al README](../README.md)

## 6.1 Alineamiento de loci con `MAFFT`

```bash
ls *.FNA | parallel 'mafft --maxiterate 1000 --globalpair --adjustdirectionaccurately --thread 2 {} > ../alingSeqs_3xCov/{.}.fasta'
```

## 6.2 Filtrado de sitios mal alineados con `TAPER`

```bash
mkdir -p Taper_filt
for F in *.fasta; do julia correction_multi.jl -m N -a N -c 1 $F > Taper_filt/$F; done
```

## 6.3 Eliminación de archivos vacíos

```bash
find . -type f -empty
find . -type f -empty -delete
```

## 6.4 Eliminación de loci con menos de 3 secuencias

```bash
./filter_bySeqnum_v0.04102026.sh
```
