# 7. Inferencia de árboles génicos y postprocesamiento

[← Volver al README](../README.md)

Existen diferentes maneras de inferir arboles filogeneticos: Los principales metodos son:

**Máxima parsimonia**: busca el árbol filogenético que explique los datos con el menor número posible de cambios evolutivos. En otras palabras, prefiere la hipótesis más simple entre las alternativas compatibles con el alineamiento.  

**Máxima verosimilitud**: busca el árbol que hace que los datos observados sean más probables bajo un modelo explícito de evolución molecular. Es decir, compara distintas topologías y elige aquella con la mayor verosimilitud dado el alineamiento y el modelo evolutivo.

**Inferencia bayesiana**: estima una distribución de probabilidad de árboles a partir de los datos, un modelo evolutivo y supuestos previos. Normalmente utiliza MCMC para explorar muchas topologías y calcular la probabilidad posterior de clados y relaciones filogenéticas

![Figure 1](https://github.com/siriusb-nox/Prociencia_Phylogen_Peru_2026/blob/main/IMG/sub_models.png?raw=true)
**Figura 1**: Algunos modelos markovianos de substitucion de nucleotidos


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
