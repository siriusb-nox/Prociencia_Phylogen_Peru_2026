# 7. Inferencia de árboles génicos y postprocesamiento

[← Volver al README](../README.md)

## 7.1 Inferencia de árboles génicos con `RAxML-NG`

```bash
mkdir -p ../raxml_gene_trees
ls *.fasta | parallel -j 2 'raxml-ng --all --msa {} --model GTR+G --bs-trees 100 --threads 2 --seed 4321 --prefix ../raxml_gene_trees/{.}'
```

## 7.2 Organización de las salidas de RAxML-NG

```bash
mkdir -p output_byproduct MLtrees_support
mv *raxml.* output_byproduct/
mv *.support MLtrees_support/
```

## 7.3 Renombrar y concatenar árboles con soporte

```bash
rename 's/\.support$/.support.tre/' *.support
cat *.tre > Orchidaceae_allSub_63g_MLtree.tre
```

## 7.4 Colapsar ramas con soporte débil

```bash
nw_ed Orchidaceae_allSub_63g_MLtree.tre 'i & b<=20' o > Orchidaceae_allSub_63g_MLtree_20LBScoll.tre
```
