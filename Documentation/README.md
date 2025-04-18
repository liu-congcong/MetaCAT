# MetaCAT

Metagenome Clustering and Association Tool.

## Documentation

The [MetaCAT](https://github.com/liu-congcong/MetaCAT) algorithm.

Try `MetaCAT [command] -h|--help` for full help.

### Commands

- `coverage`            Generate COVERAGE file from bam files.
- `seed`                Generate SEED file from an assembly file.
- `cluster`             Cluster sequences based on ASSEMBLY, COVERAGE and SEED files.
- `cluster2mapping`     Generate mapping file from fasta formatted cluster files.
- `mapping2cluster`     Generate fasta formatted cluster files from assembly and mapping file.
- `benchmarkGT`         Benchmark with ground truth.
- `benchmarkRW`         Benchmark for real-world datasets.
- `indexBam`            Index bam files.
- `representative`      Select representatives from fasta formatted cluster files.
- `abundance`           Generate abundance table from metagenomic data.
- `abundanceTest`       Significance test for abundance.
- `variant`             Call SNPs from metagenomic data.
- `mwas`                MWAS for metagenomic data.
- `plotMWAS`            QQ and Manhattan plots for MWAS results.
- `checkm2`             Estimate the quality of clusters using CheckM2.
- `gtdbtk`              Classify the clusters using GTDB-Tk.

### `coverage`

```TEXT
usage: MetaCAT coverage [options] -b <BAMs> -o <COVERAGE>

options:
  -h, --help            show this help message and exit
  -b <str> [<str> ...], --bam <str> [<str> ...]
                        Path to the sorted bam files.
                        All bam files must have the same header.
  -o <str>, --output <str>
                        Path to the coverage file.
  -ti <int>, --threads-index <int>
                        Number of threads for indexing files.
                        Works only with multiple bam files.
                        The value should be a positive integer.
                        Default: all threads.
  -tc <int>, --threads-count <int>
                        Number of threads for counting depth.
                        The value should be a positive integer.
                        Default: all threads.
  -m <int>, --min-sequence-length <int>
                        Calculations applies only to the sequences with length being greater than or equal to this value.
                        The value should be greater than 2 * "trim" + 1.
                        Default: 100.
  -l, --long            Compute coverage for long-read bam files.
                        Default: False.
  --trim <int>          Ignore the regions at both ends of a sequence when do calculations.
                        The value should be a non-negative integer.
                        Default: 0.
  --mapq <int>          MAPQ threshold.
                        The value should be a non-negative integer.
                        Default: 0.
  --identity <float>    Identity threshold.
                        The value should be from 0 to 1.
                        Default: 0.90.
```

### `seed`

```text
usage: MetaCAT seed [options] -f <ASSEMBLY> -o <SEED>

options:
  -h, --help            show this help message and exit
  -f <str>, --fasta <str>
                        Path to the fasta formatted assembly file.
  -o <str>, --output <str>
                        Path to the seed file.
  -t <int>, --threads <int>
                        Number of threads.
                        The value should be a positive integer.
                        Default: all threads.
  --fraggenescan <str>  Path to the "FragGeneScan".
                        Default: environmental FragGeneScan.
  --hmmsearch <str>     Path to the "hmmsearch".
                        Default: environmental hmmsearch.
  --hmmpress <str>      Path to the "hmmpress".
                        Default: environmental hmmpress.
```

### `cluster`

```text
usage: MetaCAT cluster [options] -f <ASSEMBLY> -c <COVERAGE> -s <SEED> -o <METACAT>

options:
  -h, --help            show this help message and exit
  -f <str>, --fasta <str>
                        Path to the fasta formatted assembly file.
  -c <str>, --coverage <str>
                        Path to the coverage file.
  -s <str>, --seed <str>
                        Path to the seed file.
  -o <str>, --output <str>
                        Prefix to the output files.
                        "*.*.fasta" the fasta formatted cluster file.
                        "*.mapping" the mapping of sequences to clusters.
  -m <int>, --min-sequence-length <int>
                        Sequences with the length being greater than or equal to this value can be involved in the clustering algorithm.
                        The value should be a positive integer.
                        Default: 1500 (#samples = 1), 500 (#samples > 1).
  -t <int>, --threads <int>
                        Number of threads.
                        The value should be a positive integer.
                        Default: all threads.
  --no-fasta-clusters   Do not output fasta formatted cluster files.
                        Default: False.
  --min-pca-variance-ratio <float>
                        Minimum percentage of variance explained in PCA in SWDPGMM.
                        The value should be from 0 to 1.
                        Default: 0.90.
  --min-pca-components <int>
                        Minimum number of components retained in PCA in SWDPGMM.
                        The value should be a positive integer.
                        Default: 20.
  --swdpgmm-engine {auto,cpu,gpu}
                        Device used to run SWDPGMM.
                        Default: auto (SWDPGMM will automatically run on the GPU if it is available).
  --min-swdpgmm-clusters <int>
                        SWDPGMM is enabled if the estimated number of genomes is greater than or equal to this value.
                        The value should be an integer.
                        Default: 100.
  --min-swdpgmm-sequences <int>
                        SWDPGMM is enabled if the number of sequences that meet the length threshold is greater than or equal to this value.
                        The value should be an integer.
                        Default: 500000.
  --max-swdpgmm-iterations <int>
                        Maximum number of SWDPGMM iterations to perform.
                        The value should be a positive integer.
                        Default: 30.
  --kmer-frequence-weight <float> [<float> ...]
                        Weights of kmer frequency probabilistic model for constructing the affinity graph.
                        The values should be from 0 to 1.
                        Default: 0.10 0.30 0.50 0.70 0.90.
  --min-kmer-frequence-probability <float>
                        Probability of kmer frequence probabilistic model lower than this value will be set to 0.
                        The value should be from 0 to 1.
                        Default: 0.50.
  --coverage-weight <float> [<float> ...]
                        Weight of read coverage probabilistic model for constructing the affinity graph.
                        The values should be from 0 to 1.
                        Default: 0.90 0.70 0.50 0.30 0.10.
  --min-coverage-probability <float>
                        Probability of coverage probabilistic model lower than this value will be set to 0.
                        The value should be from 0 to 1.
                        Default: 0.50.
  --neighbors <int>     Number of neighbors in affinity model.
                        The values should be positive integer.
                        Default: 10.
  --min-cluster <int>   Minimum size of a cluster to output.
                        The value should be a positive integer.
                        Default: 100000 bp.
  --max-seeds <int>     Maximum number of seed sequences involved in the seed model.
                        The value should be a positive integer.
                        Default: 30000.
  --random-number <int>
                        Random number generator seeded by the given integer.
                        Default: 0.
```

### `cluster2mapping`

```text
usage: MetaCAT cluster2mapping -f <FASTAs> -o <MAPPING>

options:
  -h, --help            show this help message and exit
  -f <str> [<str> ...], --fasta <str> [<str> ...]
                        Path to the fasta formatted cluster files.
  -o <str>, --output <str>
                        Path to the output file.
```

### `mapping2cluster`

```text
usage: MetaCAT mapping2cluster [options] -f <ASSEMBLY> -m <MAPPING>

options:
  -h, --help            show this help message and exit
  -f <str>, --fasta <str>
                        Path to the fasta formatted assembly file.
  -m <str>, --mapping <str>
                        Path to the mapping file.
                        A header line "Sequence ID<tab>Cluster ID" should be present.
  -o <str>, --output <str>
                        Directory to output fasta formatted cluster files.
                        Default: /home/liucongcong/data/z1.
  --checkm2 <str>       Path to the file generated by MetaCAT's checkm2, only available if "--cluster" is disable.
                        A header line "Name<tab>Completeness<tab>Contamination ..." should be present.
  --contamination <float>
                        Contamination threshold, only available if "--checkm2" is enabled.
                        The value should be from 0 to 1.
                        Default: 0.10.
  --completeness <float>
                        Completeness threshold, only available if "--checkm2" is enabled.
                        The value should be from 0 to 1.
                        Default: 0.70.
  --min-cluster <int>   Minimum size of a cluster to output.
                        The value should be a positive integer.
                        Default: 1.
  --cluster <str> [<str> ...]
                        Output only the selected clusters.
                        Default: all clusters.
```

### `benchmarkGT`

```text
usage: MetaCAT benchmarkGT [options] -gt <GROUND-TRUTH> -m <MAPPINGs> -o <BENCHMARK>

options:
  -h, --help            show this help message and exit
  -gt <str>, --ground-truth <str>
                        Ground truth file is tab-delimited, contains 3 columns and a single header line as follows.
                        Sequence ID<tab>Genome ID<tab>Length.
  -m <str> [<str> ...], --mapping <str> [<str> ...]
                        Path to the mapping files.
                        A header line "Sequence ID<tab>Cluster ID" should be present.
  -o <str>, --output <str>
                        Prefix to the output files.
  --dataset <str>       Title of dataset.
                        Default: value of "--ground-truth".
  --label <str> [<str> ...]
                        Labels of the mapping files.
                        If "--label" is selected, the number of values should be equal to the number of mapping files.
  --precision <float> [<float> ...]
                        Precision thresholds to plot.
                        The value should be from 0 to 1.
                        Default: 0.95 0.90.
  --recall <float> [<float> ...]
                        Recall thresholds to plot.
                        The value should be from 0 to 1.
                        Default: 0.95 0.90 0.80 0.70.
  --color <str>         The color string.
                        Default: Reds_r.
  --color-min <float>   Minimum color.
                        The value should be from 0 to 1.
                        Default: 0.20.
  --color-max <float>   Maximum color.
                        The value should be from 0 to 1.
                        Default: 0.75.
```

### `benchmarkRW`

```text
usage: MetaCAT benchmarkRW [options] --combine -i <CHECKM2> -o <BENCHMARK>

options:
  -h, --help            show this help message and exit
  -i <str>, --input <str>
                        Path to the file generated by MetaCAT's checkm2.
                        A header line "Name<tab>Completeness<tab>Contamination ..." should be present.
                        The format of the 1st column must be "Dataset.Program.ClusterID".
                        "Dataset" and "ClusterID" fields must contain no dots.
  -o <str>, --output <str>
                        Prefix to the output files.
  -c, --combine         Combine all datasets.
                        Default: False.
  --combined-dataset <str>
                        Title of combined dataset.
                        Default: Total Datasets.
  --contamination <float> [<float> ...]
                        Contamination thresholds to plot.
                        The value should be from 0 to 1.
                        Default: 0.05 0.10.
  --completeness <float> [<float> ...]
                        Completeness thresholds to plot.
                        The value should be from 0 to 1.
                        Default: 0.95 0.90 0.80 0.70.
  --rows <int>          Number of rows of the subplot grid.
                        The value should be a positive integer.
                        Default: auto.
  --columns <int>       Number of columns of the subplot grid.
                        The value should be a positive integer.
                        Default: auto.
  --program <str> [<str> ...]
                        Order of programs in the figure.
  --label <str> [<str> ...]
                        Labels of the programs.
                        If "--label" is selected, the number of values should be equal to the number of programs.
  --color <str>         The color string.
                        Default: Reds_r.
  --color-min <float>   Minimum color.
                        The value should be from 0 to 1.
                        Default: 0.20.
  --color-max <float>   Maximum color.
                        The value should be from 0 to 1.
                        Default: 0.75.
```

### `indexBam`

```text
usage: MetaCAT indexBam [options] -i <BAMs>

options:
  -h, --help            show this help message and exit
  -i <str> [<str> ...], --input <str> [<str> ...]
                        Path to the sorted bam files.
                        All bam files must have the same header.
  -t <int>, --threads <int>
                        Number of threads.
                        Works only with multiple bam files.
                        The value should be a positive integer.
                        Default: all threads.
```

### `representative`

```text
usage: MetaCAT representative [options] -f <FASTAs> -c <CHECKM2> -g <GTDBTK> -o <REPRESENTATIVE>

options:
  -h, --help            show this help message and exit
  -f <str> [<str> ...], --fasta <str> [<str> ...]
                        Path to the fasta files.
  -c <str>, --checkm2 <str>
                        Path to the file generated by MetaCAT's checkm2.
                        A header line "Name<tab>Completeness<tab>Contamination<tab>...<tab>...<tab>...<tab>Contig_N50 ..." should be present.
  -g <str>, --gtdbtk <str>
                        Path to the file generated by MetaCAT's gtdbtk.
                        A header line "user_genome<tab>classification ..." should be present.
  -o <str>, --output <str>
                        Prefix to the output files.
                        "*.annotation": the classifications of all clusters.
                        "*.assembly": the combined representative assembly.
                        "*.mapping": the mapping of sequences to clusters.
  -t <int>, --threads <int>
                        Number of threads.
                        The value should be a positive integer.
                        Default: all threads.
  --mash <str>          Path to the "mash".
                        Default: environmental mash.
  --distance <float>    Maximum mash distance.
                        The value should be from 0 to 1.
                        Default: 0.05.
  --kmer-size <int>     K-mer size used in Mash.
                        The value should be a positive integer.
                        Default: 21.
  --sketch-size <int>   Number of non-redundant min-hashes used in Mash.
                        The value should be a positive integer.
                        Default: 1000.
  --contamination <float>
                        Contamination threshold.
                        The value should be from 0 to 1.
                        Default: 0.10.
  --completeness <float>
                        Completeness threshold.
                        The value should be from 0 to 1.
                        Default: 0.70.
  --rename              Sequences are stored in the assembly file from 1 to N.
                        Default: False.
```

### `abundance`

```text
usage: MetaCAT abundance [options] -a <ANNOTATION> -b <BAMs> -m <MAPPINGs> -o <ABUNDANCE>

options:
  -h, --help            show this help message and exit
  -a <str>, --annotation <str>
                        Path to the annotation file generated by MetaCAT's representative.
                        A header line "Cluster ID<tab>Classification" should be present.
  -b <str> [<str> ...], --bam <str> [<str> ...]
                        Path to the sorted bam files.
  -m <str> [<str> ...], --mapping <str> [<str> ...]
                        Path to the mapping files.
                        A header line "Sequence ID<tab>Cluster ID" should be present.
                        The number of values should be equal to the number of bam files, if the number of values is not one.
                        Otherwise all bam files must have the same header.
  -o <str>, --output <str>
                        Path to the output file.
  -ti <int>, --threads-index <int>
                        Number of threads for indexing files.
                        Works only with multiple bam files.
                        The value should be a positive integer.
                        Default: all threads.
  -tc <int>, --threads-count <int>
                        Number of threads for counting depth.
                        The value should be a positive integer.
                        Default: all threads.
  -l, --long            Compute abundance for long-read bam files.
                        Default: False.
  --trim <int>          Ignore the regions at both ends of a sequence when do calculations.
                        The value should be a non-negative integer.
                        Default: 0.
  --mapq <int>          MAPQ threshold.
                        The value should be a non-negative integer.
                        Default: 0.
  --identity <float>    Identity threshold.
                        The value should be from 0 to 1.
                        Default: 0.90.
  --min-abundance <float>
                        Absolute abundance less than this value will be set to 0.
                        The value should be from 0 to 1.
                        Default: 0.10.
```

### `abundanceTest`

```text
usage: MetaCAT abundanceTest [options] -a <ABUNDANCE> -g <GROUP> -o <TEST>

options:
  -h, --help            show this help message and exit
  -a <str>, --abundance <str>
                        Path to the abundance file.
                        A header line "Abundance<tab>ID1<tab>ID2 ..." should be present.
  -g <str>, --group <str>
                        Path to the group file.
                        A header line "ID<tab>Group" should be present.
  -o <str>, --output <str>
                        Prefix to the output files.
                        "*.G1-G2.*: significance test for group G1 and G2."
  --coverage <float>    Minimum coverage in each of groups.
                        The value should be from 0 to 1.
                        Default: 0.50.
  --rows <int>          Number of rows of the subplot grid.
                        The value should be a positive integer.
                        Default: auto.
  --columns <int>       Number of columns of the subplot grid.
                        The value should be a positive integer.
                        Default: auto.
  --classification <str> [<str> ...]
                        Classification levels used for comparison.
                        The values can be s, g, f, o, c, p, d.
                        Default: s.
  --comparison <str> <str>
                        Paired groups for comparison.
                        Default: all paired groups.
  --alpha <float>       The probability of making the wrong decision when the null hypothesis is true.
                        Default: 0.05.
  --multiple-test <str>
                        Method for multiple hypothesis testing.
                        The values can be bonferroni, benjamini-hochberg, none.
                        Default: bonferroni.
```

### `variant`

```text
usage: MetaCAT variant [options] -i <BAMs> -o <VARIANT>

options:
  -h, --help            show this help message and exit
  -i <str> [<str> ...], --input <str> [<str> ...]
                        Path to the sorted bam files.
                        All bam files must have the same header.
  -o <str>, --output <str>
                        Path to the output file.
  -tc <int>, --threads-call <int>
                        Number of threads for calling variants.
                        The value should be a positive integer.
                        Default: all threads.
  -ti <int>, --threads-index <int>
                        Number of threads for indexing files.
                        Works only with multiple bam files.
                        The value should be a positive integer.
                        Default: all threads.
  --coverage <float>    Minimum population coverage of a SNP.
                        Default: 0.50.
  --mapq <int>          MAPQ threshold.
                        The value should be a non-negative integer.
                        Default: 0.
  --identity <float>    Identity threshold.
                        The value should be from 0 to 1.
                        Default: 0.90.
  --depth <int>         Read depth for a sample less than this value will be set to 0.
                        The value should be a positive integer.
                        Default: 1.
```

### `mwas`

```text
usage: MetaCAT mwas [options] -p <PHENOTYPE> -m <MAPPING> -ad <ABUNDANCE> -at <ANNOTATION> -v <VARIANT> -o <MWAS>

options:
  -h, --help            show this help message and exit
  -p <str>, --phenotype <str>
                        Path to the phenotype file.
                        A header line "ID<tab>Phenotype1<tab>Phenotype2 ..." should be present.
                        Each phenotype is a numeric value.
  -m <str>, --mapping <str>
                        Path to the mapping file.
                        A header line "Sequence ID<tab>Cluster ID" should be present.
  -ad <str>, --abundance <str>
                        Path to the abundance file.
                        A header line "Abundance<tab>ID1<tab>ID2 ..." should be present.
  -at <str>, --annotation <str>
                        Path to the annotation file generated by MetaCAT's representative.
                        A header line "Cluster ID<tab>Classification" should be present.
  -v <str>, --variant <str>
                        Path to the variant file.
                        A header line "Chromosome<tab>Position<tab>Population<tab>ID1<tab>ID2 ..." should be present.
  -o <str>, --output <str>
                        Path to the output file.
  -c <str>, --covariate <str>
                        Path to the covariate file.
                        A header line "ID<tab>Covariate1<tab>Covariate2 ..." should be present.
  -t <int>, --threads <int>
                        Number of threads.
                        Default: 56.
  --phenotype-column <int>
                        Column of phenotype.
                        Default: 2.
  --mwas <str>          Path to the mwas file.
                        Major alleles will be determined in this file rather than in the variant file.
                        A header line "Classification<tab>Chromosome<tab>Position<tab>Allele<tab>Beta<tab>SE<tab>P<tab>Significance" should be present.
  --size <float>        Minimum proportion of samples for mwas.
                        Default: 0.50.
  --pca-components <int>
                        Number of components of GDM to keep.
                        Default: 5.
  --alpha <float>       The probability of making the wrong decision when the null hypothesis is true.
                        Default: 0.05.
```

### `plotMWAS`

```text
usage: MetaCAT plotMWAS [options] -m <MWAS> -o <MWAS>

options:
  -h, --help            show this help message and exit
  -m <str>, --mwas <str>
                        Path to the mwas file.
                        A header line "Classification<tab>Chromosome<tab>Position<tab>Allele<tab>Beta<tab>SE<tab>P<tab>Significance" should be present.
  -o <str>, --output <str>
                        Prefix to the output files.
                        "*.qq.pdf": the qq plot for the mwas.
                        "*.manhattan": the manhattan plot for the mwas.
  --variant-annotation <str>
                        Path to the variant annotation file.
                        A header line "Chromosome<tab>Position<tab>Annotation ..." should be present.
  --alpha <float>
                        The probability of making the wrong decision when the null hypothesis is true.
                        Default: 0.05.
  --max-negative-log-p <float>
                        Maximum value of -log10(P).
                        Default: 20.
```

### `checkm2`

```text
usage: MetaCAT checkm2 [options] -i <FASTAs> -o <CHECKM2>

options:
  -h, --help            show this help message and exit
  -i <str> [<str> ...], --input <str> [<str> ...]
                        Path to the fasta files.
  -o <str>, --output <str>
                        Path to the output file.
  -t <int>, --threads <int>
                        Number of threads.
                        The value should be a positive integer.
                        Default: all threads.
  --checkm2 <str>       Path to the "checkm2".
                        Default: environmental checkm2.
```

### `gtdbtk`

```text
usage: MetaCAT gtdbtk [options] -i <FASTAs> -o <GTDBTk>

options:
  -h, --help            show this help message and exit
  -i <str> [<str> ...], --input <str> [<str> ...]
                        Path to the fasta files.
  -o <str>, --output <str>
                        Path to the output file.
  --mash-db <str>       Path to the Mash reference sketch database used by GTDB-Tk.
                        If not specified, it will be automatically created as "GTDB-Tk.msh".
                        Default: GTDB-Tk.msh.
  -c <str>, --checkm2 <str>
                        Path to the file generated by MetaCAT's checkm2.
                        A header line "Name<tab>Completeness<tab>Contamination ..." should be present.
  --contamination <float>
                        Contamination threshold, only available if "--checkm2" is enabled.
                        The value should be from 0 to 1.
                        Default: 0.10.
  --completeness <float>
                        Completeness threshold, only available if "--checkm2" is enabled.
                        The value should be from 0 to 1.
                        Default: 0.70.
  -t <int>, --threads <int>
                        Number of threads.
                        The value should be a positive integer.
                        Default: all threads.
  --gtdbtk <str>        Path to the "gtdbtk".
                        Default: environmental gtdbtk.
```

## File Formats

### `abundance`

| `Abundance` | *`ID1`* | *`ID2`* | `...` |
| :---------: | :-----  | ------: | :---: |
|    `...`    |  `...`  |  `...`  | `...` |

### `annotation`

| `Cluster ID` | `Classification` |
| :----------: | :--------------: |
|    `...`     |      `...`       |

### `checkm2`

| `Name` | `Completeness` | `Contamination` | `*` | `*` | `*` | `Contig_N50` | `...` |
| :----: | :------------: | :-------------: | :-: | :-: | :-: | :----------: | :---: |
| `...`  |     `...`      |      `...`      |`...`|`...`|`...`|    `...`     | `...` |

### `covariate`

| `ID` | *`Covariate1`* | *`Covariate2`* | `...` |
| :--: | :------------: | :------------: | :---: |
|`...` |     `...`      |     `...`      | `...` |

### `ground-truth`

| `Sequence ID` | `Genome ID` | `Length` |
| :-----------: | :---------: | :------: |
|     `...`     |    `...`    |  `...`   |

### `group`

| `ID` | `Group` |
| :--: | :-----: |
|`...` |  `...`  |

### `gtdbtk`

| `user_genome` | `classification` | `...` |
| :-----------: | :--------------: | :---: |
|     `...`     |      `...`       | `...` |

### `mapping`

| `Sequence ID` | `Cluster ID` |
| :-----------: | :----------: |
|     `...`     |    `...`     |

### `mwas`

| `Classification` | `Chromosome` | `Position` | `Allele` | `Beta` | `SE` | `P` | `Significance` | `...` |
| :--------------: | :----------: | :--------: | :------: | :----: | :--: | :-: | :------------: | :---: |
|      `...`       |    `...`     |   `...`    |  `...`   | `...`  |`...` |`...`|     `...`      | `...` |

### `phenotype`

| `ID` | `Phenotype` |
| :--: | :---------: |
|`...` |    `...`    |

### `variant`

| `Chromosome` | `Position` | `Population` | *`ID1`* | *`ID2`* | `...` |
| :----------: | :--------: | :----------: | :-----  | ------: | :---: |
|    `...`     |   `...`    |    `...`     |  `...`  |  `...`  | `...` |
