# Taller de "Principios en Filogenómica" - 15-19 Abril (Madre de Dios) y 22-24 Abril (UNTRM) 
## Instructors: [Oscar Alejandro Perez-Escobar](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.70484) [(Royal Botanic Gardens, Kew)](https://scholar.google.com/citations?user=tSzyp6QAAAAJ&hl=es) & [Alexander Damian](https://cameron.botany.wisc.edu/staff/damian-l-alexander/) [(University of Wisconsin-Madison)](https://scholar.google.co.uk/citations?user=RwjJKcwAAAAJ&hl=th)
##### Este tayer esta apoyado por [Prociencia](https://prociencia.gob.pe/) & la [Swiss Orchid Foundation](https://novagenetix.ch/swiss-orchid-foundation/)
## Introduccion
Este repositorio contiene un tutorial introductorio para una práctica de filogenómica ejecutada en una máquina virtual de Ubuntu sobre portátiles con Windows. El tutorial cubre un flujo de trabajo mínimo para descargar lecturas Illumina, limpiar los datos, ensamblar loci in-silico, alinear genes, filtrar loci y sitios problemáticos, inferir árboles génicos de máxima verosimilitud y arboles de especies usando un metodo de "Multispecies Coalescent" (MSC).

En este taller, emplearemos una pequeña muestra de datos de la filogenia de orquideas de Perez-Escobar et al. (2024, [_New. Phyt._](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.19580)), que representa las cinco subfamilias mas un grupo externo (_Dioscorea_). El objetivo del taller es lograr ensamblar una filogenia a partir de sets de datos illumina los cuales primero deben ser obtenidos de la base de datos SRA del [NCBI](https://www.ncbi.nlm.nih.gov/) asi como tambien aprender a diseñar proyectos de investigacion en filogenomica. 

<p align="center">
  <img src="https://github.com/siriusb-nox/Prociencia_Phylogen_Peru_2026/blob/main/IMG/Orchid_ToL3.0_v5_OP_distproof_LR_96dpi_MQ_online_EN.jpg" alt="Orchid Tree of Life v2.0"/>
</p>


---

## 1. Entorno computacional del taller

La práctica fue diseñada para ejecutarse en una máquina virtual de Ubuntu en VirtualBox. Alternativamente, los ejercicios del taller se podran implementar en un ambiente Linux, asumiendo que todas las dependencias han sido instaladas.

### 1.1 Información de la máquina virtual

- **Sistema operativo invitado:** Ubuntu
- **Software de virtualización:** VirtualBox
- **Usuario:** `vboxuser`

> **La contraseña de cada sesion Ubuntu sera entregada el dia de la practica**.

### 1.2 Programas requeridos
Para ejecutar este taller, se requieren los siguientes programas:

- [VirtualBox](https://www.virtualbox.org/wiki/Downloads)
- Ubuntu (máquina virtual)
- Bash
- [SRA Toolkit: `prefetch` y `fasterq-dump`](https://github.com/ncbi/sra-tools)
- [fastp](https://github.com/opengene/fastp)
- [conda](https://docs.conda.io/)
- [HybPiper](https://github.com/mossmatters/HybPiper)
- [BWA](https://github.com/lh3/BWA)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- [Julia](https://julialang.org/)
- [TAPER](https://github.com/chaoszhang/TAPER)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- [Newick Utilities](https://github.com/tjunier/newick_utils)
- utilidades estándar de UNIX: `grep`, `find`, `rename`, `cat`, `mv`, `ls`, `gzip`


### 1.3 Estructura del pipeline
Este tutorial esta dividio en cuatro principales (Figura 1):

A. [Descarga y limpieza de datos Illumina](docs/04_descarga_y_limpieza_datos_illumina.md)

B. [Ensamblaje de captura dirigida y recuperación de secuencias](docs/05_ensamblaje_y_recuperacion_hybpiper.md)

C. [Alineamiento por locus y filtrado](docs/06_alineamiento_y_filtrado.md)

D. [Inferencia de árboles génicos y postprocesamiento](docs/07_arboles_genicos_raxml_ng.md)

E. [Inferencia del árbol de especies con ASTRAL](docs/08_arbol_de_especies_astral.md)

![Figure 1](https://github.com/siriusb-nox/Prociencia_Phylogen_Peru_2026/blob/main/IMG/image.png?raw=true)
**Figura 1**: Vista simplificada del tutorial/pipeline


---

## 2. Introducción muy breve a Bash

Esta sección está adaptada del tutorial de Bash usado en el repositorio del taller de ONT y simplificada para esta práctica de filogenómica.

### 2.1 ¿Qué es Bash?

Bash es la shell usada en muchos sistemas operativos de tipo UNIX, como Ubuntu y Linux. En bioinformática, la shell es especialmente útil porque muchas herramientas están diseñadas para ejecutarse desde la terminal y porque Bash facilita el trabajo con archivos de texto, secuencias y directorios.

Un comando típico tiene esta estructura:

```bash
programa opcion1 opcion2 entrada > salida
```

Por ejemplo:

```bash
ls
ls -lrt
```

El primer comando lista archivos. El segundo los muestra con más detalle y los ordena por fecha.

### 2.2 Comandos básicos

#### `pwd`

Muestra el directorio actual:

```bash
pwd
```

#### `ls`

Lista los archivos del directorio actual:

```bash
ls
ls -l
```

#### `cd`

Cambia de directorio:

```bash
cd
cd ..
cd /home/vboxuser/Documents
```

#### `mkdir`

Crea directorios:

```bash
mkdir nueva_carpeta
mkdir carpeta1 carpeta2 carpeta3
mkdir -p ~/Documents/Phylo_work_2026
```

#### `cat`

Muestra o concatena archivos de texto:

```bash
cat archivo.txt
cat archivo1.txt archivo2.txt > combinado.txt
```

#### `less`

Permite inspeccionar archivos de texto página por página:

```bash
less archivo.txt
```

#### `grep`

Busca patrones de texto:

```bash
grep '^>' archivo.fasta
grep -c '^>' archivo.fasta
```

#### `man`

Muestra la página de ayuda de un comando:

```bash
man ls
man grep
```

### 2.3 Ideas útiles de Bash para este taller

#### Bucles

Un bucle repite un comando para muchos archivos o muestras:

```bash
for F in *.fasta; do echo "$F"; done
```

#### Redirección de salida

El símbolo `>` envía la salida a un archivo:

```bash
cat entrada.txt > salida.txt
```

#### Tuberías

Las tuberías envían la salida de un programa a otro:

```bash
ls *.FNA | wc -l
```

#### Ejecución en paralelo

GNU Parallel permite correr muchos trabajos al mismo tiempo:

```bash
ls *.FNA | parallel 'echo {}'
```

---

## 3. Estructura sugerida de directorios de trabajo

Una estructura ordenada de carpetas hace que la práctica sea mucho más fácil de seguir. 

```bash
mkdir -p /media/foldername/Phylo_work_2026_Peru/{raw_reads,trimmed,HybPiper,alignSeqs_3xCov,Taper_filt,raxml_gene_trees}
```

También puedes crear carpetas separadas para scripts, logs y archivos de referencia.

---

## 3. Notas para estudiantes

- Comprueba siempre que estás en el directorio correcto con `pwd`.
- Usa `ls` con frecuencia para seguir la pista de tus archivos.
- Crea las carpetas de salida antes de redirigir archivos hacia ellas.
- Guarda los entornos y programas en el sistema de archivos local de Ubuntu y no dentro de carpetas compartidas de VirtualBox.
- Para trabajos largos, guarda los comandos en scripts de shell para poder relanzarlos fácilmente.
- Para resolver problemas, inspecciona la salida de `less`, `cat` y `grep` antes de borrar archivos.

---

## 4. Agradecimientos

La introducción breve a Bash usada en este taller fue adaptada del tutorial de Bash empleado en el repositorio del taller de ONT del mismo autor y reformulada aquí para una práctica de filogenómica.
