# MetaCAT

Metagenome Clustering and Association Tool.

## Getting started

- [`Group metagenomic sequences into clusters`](./README.md#Group-metagenomic-sequences-into-clusters)
- [`Evaluate the clusters`](./README.md#Evaluate-the-clusters)
- [`Classify the clusters`](./README.md#Classify-the-clusters)
- [`Define a set of representative genomes`](./README.md#Define-a-set-of-representative-genomes)
- [`Compute the relative abundance table`](./README.md#Compute-the-relative-abundance-table)
- [`Determine the differential marker species`](./README.md#Determine-the-differential-marker-species)
- [`Call microbial SNPs`](./README.md#Call-microbial-SNPs)
- [`Perform MWAS for metagenomic data`](./README.md#Perform-MWAS-for-metagenomic-data)

### Group metagenomic sequences into clusters

#### Preparations

- A fasta formatted assembly file `dataset.assembly`.

- A set of sorted bam files with the same header `dataset.*.bam`.

#### Generate COVERAGE file from bam files

```bash
MetaCAT coverage --bam dataset.*.bam --output dataset.metacat.coverage
```

#### Generate SEED file from an assembly file

```bash
MetaCAT seed --fasta dataset.assembly --output dataset.metacat.seed
```

#### Cluster sequences based on ASSEMBLY, COVERAGE and SEED files

```bash
MetaCAT cluster --fasta dataset.assembly --coverage dataset.metacat.coverage --seed dataset.metacat.seed --output dataset.metacat
```

This step will generate a mapping file `dataset.metacat.mapping`,

and a set of fasta formatted cluster files `dataset.metacat.*.fasta`.

### Evaluate the clusters

#### Preparations

- A set of fasta formatted cluster files `dataset.metacat.*.fasta`.

- `checkm2` and its dataset are installed and available in the system PATH.

- Activate the `checkm2` environment first if it is installed!

#### Estimate the quality of clusters using CheckM2

```bash
MetaCAT checkm2 --fasta dataset.metacat.*.fasta --output dataset.metacat.checkm2
```

The output file is identical to the `quality_report.tsv` file produced by `checkm2 predict`.

#### Benchmark for real-world datasets

```bash
MetaCAT benchmarkRW --combine --checkm2 dataset.metacat.checkm2 --output dataset.metacat.benchmarkRW
```

The format of the 1st column in `checkm2` must be "Dataset.Program.ClusterID".

"Dataset" and "ClusterID" fields must contain no dots.

This step will generate a benchmark file `dataset.metacat.benchmarkRW.tsv`,

and a set of pdf formatted plot files `dataset.metacat.benchmarkRW.*.pdf`.

### Classify the clusters

#### Preparations

- A set of fasta formatted cluster files `dataset.metacat.*.fasta`.

- `gtdbtk` and its database are installed and available in the system PATH.

- Activate the `gtdbtk` environment first if it is installed!

#### Classify the clusters using GTDB-Tk

```bash
MetaCAT gtdbtk --fasta dataset.metacat.*.fasta --output dataset.metacat.gtdbtk
```

The output file is identical to the concatenated `*.summary.tsv` files produced by `gtdbtk classify_wf`.

If the file is generated manually, ensure that it includes only one header line.

### Define a set of representative genomes

#### Preparations

- A set of fasta formatted cluster files `dataset.metacat.*.fasta`.

- Quality assignment results `dataset.metacat.checkm2`.

- Classification results `dataset.metacat.gtdbtk`.

- `fastANI` is installed and available in the system PATH, FastANI is used to compute the similarity between paired genomes.

#### Select representatives from fasta formatted cluster files

```bash
MetaCAT representative --fasta dataset.metacat.*.fasta --checkm2 dataset.metacat.checkm2 --gtdbtk dataset.metacat.gtdbtk --output dataset.metacat.representative
```

This step will generate a mapping file `dataset.metacat.representative.mapping`,

a fasta formatted assembly file `dataset.metacat.representative.assembly`,

and an annotation file `dataset.metacat.representative.annotation`.

After this step, all bam files must be regenerated using the new assembly `dataset.metacat.representative.assembly`.

### Compute the relative abundance table

#### Preparations

- A representative annotation file `dataset.metacat.representative.annotation`.

- A representative mapping file `dataset.metacat.representative.mapping`.

- A set of sorted bam files with the same header `dataset.*.bam`.

#### Generate abundance table from metagenomic data

```bash
MetaCAT abundance --annotation dataset.metacat.representative.annotation --mapping dataset.metacat.representative.mapping \
--bam dataset.*.bam --output dataset.metacat.abundance
```

Sample IDs are defined as the basenames of the bam files, excluding file extensions.

If bam files are stored on an HDD instead of an NVMe drive,

it is recommended to reduce the number of threads (e.g., -tc 20 -ti 20).

### Determine the differential marker species

#### Preparations

- An abundance file `dataset.metacat.abundance`.

- A group file [`group`](./README.md#group).

#### Significance test for abundance

```bash
MetaCAT abundanceTest --abundance dataset.metacat.abundance --group group --output dataset.metacat.abundanceTest
```

This step will generate a set of pdf formatted plot files `dataset.metacat.abundanceTest.*.pdf`.

### Call microbial SNPs

#### Preparations

- A representative annotation file `dataset.metacat.representative.annotation`.

- A representative mapping file `dataset.metacat.representative.mapping`.

- A set of sorted bam files with the same header `dataset.*.bam`.

#### Call SNPs from metagenomic data

```bash
MetaCAT variant --annotation dataset.metacat.representative.annotation --mapping dataset.metacat.representative.mapping \
--bam dataset.*.bam --output dataset.metacat.variant
```

Sample IDs are defined as the basenames of the bam files, excluding file extensions.

If bam files are stored on an HDD instead of an NVMe drive,

it is recommended to reduce the number of threads (e.g., -tc 20 -ti 20).

### Perform MWAS for metagenomic data

#### Preparations

- An abundance file `dataset.metacat.abundance`.

- A phenotype file [`phenotype`](./README.md#phenotype).

- A covariate file [`covariate`](./README.md#covariate).

- A variant file `dataset.metacat.variant`.

#### MWAS for metagenomic data

```bash
MetaCAT mwas --abundance dataset.metacat.abundance --phenotype phenotype \
--covariate covariate --variant dataset.metacat.variant --output dataset.metacat.mwas
```

#### QQ and Manhattan plots for MWAS results

```bash
MetaCAT plotMWAS --mwas dataset.metacat.mwas --output dataset.metacat.plotMWAS
```

This step will generate two pdf formatted plot files `dataset.metacat.plotMWAS.*.pdf`.

---

## Commands

- [`coverage`](./README.md#coverage)           Generate COVERAGE file from bam files.
- [`seed`](./README.md#seed)               Generate SEED file from an assembly file.
- [`cluster`](./README.md#cluster)            Cluster sequences based on ASSEMBLY, COVERAGE and SEED files.
- [`cluster2mapping`](./README.md#cluster2mapping)    Generate mapping file from fasta formatted cluster files.
- [`mapping2cluster`](./README.md#mapping2cluster)    Generate fasta formatted cluster files from assembly and mapping file.
- [`benchmarkGT`](./README.md#benchmarkGT)        Benchmark with ground truth.
- [`benchmarkRW`](./README.md#benchmarkRW)        Benchmark for real-world datasets.
- [`indexBam`](./README.md#indexBam)           Index bam files.
- [`representative`](./README.md#representative)     Select representatives from fasta formatted cluster files.
- [`abundance`](./README.md#abundance)          Generate abundance table from metagenomic data.
- [`abundanceTest`](./README.md#abundanceTest)      Significance test for abundance.
- [`variant`](./README.md#variant)            Call SNPs from metagenomic data.
- [`mwas`](./README.md#mwas)               MWAS for metagenomic data.
- [`plotMWAS`](./README.md#plotMWAS)           QQ and Manhattan plots for MWAS results.
- [`checkm2`](./README.md#checkm2)            Estimate the quality of clusters using CheckM2.
- [`gtdbtk`](./README.md#gtdbtk)             Classify the clusters using GTDB-Tk.

Try `MetaCAT [command] -h|--help` for full help.

### coverage

```text
usage: MetaCAT coverage [options] -b <BAMs> -o <COVERAGE>

options:
  -h, --help            show this help message and exit
  -b <str> [<str> ...], --bam <str> [<str> ...]
                        Path to the sorted bam files or directories containing them.
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
  --bam-suffix <str>    If directories are provided with "-b|--bam", only files with the specified suffix will be selected.
                        Default: bam.
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

### seed

```text
usage: MetaCAT seed [options] -f <ASSEMBLY> -o <SEED>

optional arguments:
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
                        Default: internal FragGeneScan.
  --hmmsearch <str>     Path to the "hmmsearch".
                        Default: internal hmmsearch.
```

### cluster

```text
usage: MetaCAT cluster [options] -f <ASSEMBLY> -c <COVERAGE> -s <SEED> -o <METACAT>

optional arguments:
  -h, --help            show this help message and exit
  -f <str>, --fasta <str>
                        Path to the fasta formatted assembly file.
  -c <str>, --coverage <str>
                        Path to the coverage file.
  -s <str>, --seed <str>
                        Path to the seed file.
  -o <str>, --output <str>
                        Prefix of the output files.
                        "*.*.fasta" the fasta formatted cluster file.
                        "*.mapping" the mapping of sequences to clusters.
  -m <int>, --min-sequence-length <int>
                        Sequences with the length being greater than or equal to this value can be involved in the clustering algorithm.
                        The value should be a positive integer.
                        Default: 1000 (#samples = 1), 500 (#samples > 1).
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
                        Default: 200.
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
                        Default: 0.90 0.70 0.50 0.30 0.10.
  --min-kmer-frequence-probability <float>
                        Probability of kmer frequence probabilistic model lower than this value will be set to 0.
                        The value should be from 0 to 1.
                        Default: 0.50.
  --coverage-weight <float> [<float> ...]
                        Weight of read coverage probabilistic model for constructing the affinity graph.
                        The values should be from 0 to 1.
                        Default: 0.10 0.30 0.50 0.70 0.90.
  --min-coverage-probability <float>
                        Probability of coverage probabilistic model lower than this value will be set to 0.
                        The value should be from 0 to 1.
                        Default: 0.50.
  --seed-neighbors <int>
                        Number of neighbors in seed affinity model.
                        The values should be positive integer.
                        Default: 10.
  --sequence-neighbors <int>
                        Number of neighbors in sequence affinity model.
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

### cluster2mapping

```text
usage: MetaCAT cluster2mapping [options] -f <FASTAs> -o <MAPPING>

optional arguments:
  -h, --help            show this help message and exit
  -f <str> [<str> ...], --fasta <str> [<str> ...]
                        Path to the fasta formatted cluster files or directories containing them.
  -o <str>, --output <str>
                        Path to the output file.
  --fasta-suffix <str>  If directories are provided with "-f|--fasta", only files with the specified suffix will be selected.
                        Default: fasta.
  --min-cluster <int>   Minimum size of a cluster to output.
                        The value should be a positive integer.
                        Default: 1.
```

### mapping2cluster

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

### benchmarkGT

```text
usage: MetaCAT benchmarkGT [options] -gt <GROUND-TRUTH> -m <MAPPINGs> -o <BENCHMARK>

optional arguments:
  -h, --help            show this help message and exit
  -gt <str>, --ground-truth <str>
                        Ground truth file is tab-delimited, contains 3 columns and a single header line as follows.
                        Sequence ID<tab>Genome ID<tab>Length.
  -m <str> [<str> ...], --mapping <str> [<str> ...]
                        Path to the mapping files.
                        A header line "Sequence ID<tab>Cluster ID" should be present.
  -o <str>, --output <str>
                        Prefix of the output files.
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

### benchmarkRW

```text
usage: MetaCAT benchmarkRW [options] --combine -c <CHECKM2> -o <BENCHMARK>

optional arguments:
  -h, --help            show this help message and exit
  -c <str>, --checkm2 <str>
                        Path to the file generated by MetaCAT's checkm2.
                        A header line "Name<tab>Completeness<tab>Contamination ..." should be present.
                        The format of the 1st column must be "Dataset.Program.ClusterID".
                        "Dataset" and "ClusterID" fields must contain no dots.
  -o <str>, --output <str>
                        Prefix of the output files.
  --combine             Combine all datasets.
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

### indexBam

```text
usage: MetaCAT indexBam [options] -b <BAMs>

optional arguments:
  -h, --help            show this help message and exit
  -b <str> [<str> ...], --bam <str> [<str> ...]
                        Path to the sorted bam files or directories containing them.
                        All bam files must have the same header.
  -t <int>, --threads <int>
                        Number of threads.
                        Works only with multiple bam files.
                        The value should be a positive integer.
                        Default: all threads.
  --bam-suffix <str>    If directories are provided with "-b|--bam", only files with the specified suffix will be selected.
                        Default: bam.
```

### representative

```text
usage: MetaCAT representative [options] -f <FASTAs> -c <CHECKM2> -g <GTDBTK> -o <REPRESENTATIVE>

optional arguments:
  -h, --help            show this help message and exit
  -f <str> [<str> ...], --fasta <str> [<str> ...]
                        Path to the fasta formatted files or directories containing them.
  -c <str>, --checkm2 <str>
                        Path to the file generated by MetaCAT's checkm2.
                        A header line "Name<tab>Completeness<tab>Contamination ..." should be present.
  -g <str>, --gtdbtk <str>
                        Path to the file generated by MetaCAT's gtdbtk.
                        A header line "user_genome<tab>classification ..." should be present.
  -o <str>, --output <str>
                        Prefix of the output files.
                        "*.annotation": the classifications of all clusters.
                        "*.assembly": the combined representative assembly.
                        "*.mapping": the mapping of sequences to clusters.
  -t <int>, --threads <int>
                        Number of threads.
                        The value should be a positive integer.
                        Default: all threads.
  -e {fastANI,mash}, --engine {fastANI,mash}
                        Engine for computing the similarity between paired genomes.
                        Default: mash.
  --mash <str>          Path to the "mash".
                        Default: internal mash.
  --mash-distance <float>
                        Maximum mash distance.
                        The value should be from 0 to 1.
                        Default: 0.05.
  --mash-kmer-size <int>
                        K-mer size used in Mash.
                        The value should be a positive integer.
                        Default: 21.
  --mash-sketch-size <int>
                        Number of non-redundant min-hashes used in Mash.
                        The value should be a positive integer.
                        Default: 5000.
  --fastani <str>       Path to the "fastANI".
                        Default: environmental fastANI.
  --fastani-ani <float>
                        Minimum ANI to consider genomes as the same species.
                        The value should be from 0 to 100.
                        Default: 95.0.
  --fastani-kmer-size <int>
                        K-mer size used in fastANI.
                        The value should be a positive integer.
                        Default: 16.
  --fastani-fragment-length <int>
                        Fragment length used in Mash.
                        The value should be a positive integer.
                        Default: 3000.
  --fasta-suffix <str>  If directories are provided with "-f|--fasta", only files with the specified suffix will be selected.
                        Default: fasta.
  --contamination <float>
                        Contamination threshold.
                        The value should be from 0 to 1.
                        Default: 0.10.
  --completeness <float>
                        Completeness threshold.
                        The value should be from 0 to 1.
                        Default: 0.70.
```

### abundance

```text
usage: MetaCAT abundance [options] -a <ANNOTATION> -b <BAMs> -m <MAPPINGs> -o <ABUNDANCE>

optional arguments:
  -h, --help            show this help message and exit
  -a <str>, --annotation <str>
                        Path to the annotation file generated by MetaCAT's representative.
                        A header line "Cluster ID<tab>Classification" should be present.
  -b <str> [<str> ...], --bam <str> [<str> ...]
                        Path to the sorted bam files or directories containing them.
                        All bam files must have the same header.
  -m <str>, --mapping <str>
                        Path to the mapping file generated by MetaCAT's representative.
                        A header line "Sequence ID<tab>Cluster ID" should be present.
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
  --bam-suffix <str>    If directories are provided with "-b|--bam", only files with the specified suffix will be selected.
                        Default: bam.
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

### abundanceTest

```text
usage: MetaCAT abundanceTest [options] -a <ABUNDANCE> -g <GROUP> -o <TEST>

optional arguments:
  -h, --help            show this help message and exit
  -a <str>, --abundance <str>
                        Path to the abundance file.
                        A header line "Abundance<tab>ID1<tab>ID2 ..." should be present.
  -g <str>, --group <str>
                        Path to the group file.
                        A header line "ID<tab>Group" should be present.
  -o <str>, --output <str>
                        Prefix of the output files.
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

### variant

```text
usage: MetaCAT variant [options] -a <ANNOTATION> -b <BAMs> -m <MAPPING> -o <VARIANT>

optional arguments:
  -h, --help            show this help message and exit
  -a <str>, --annotation <str>
                        Path to the annotation file generated by MetaCAT's representative.
                        A header line "Cluster ID<tab>Classification" should be present.
  -b <str> [<str> ...], --bam <str> [<str> ...]
                        Path to the sorted bam files or directories containing them.
                        All bam files must have the same header.
  -m <str>, --mapping <str>
                        Path to the mapping file generated by MetaCAT's representative.
                        A header line "Sequence ID<tab>Cluster ID" should be present.
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
  --bam-suffix <str>    If directories are provided with "-b|--bam", only files with the specified suffix will be selected.
                        Default: bam.
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

### mwas

```text
usage: MetaCAT mwas [options] -p <PHENOTYPE> -a <ABUNDANCE> -v <VARIANT> -o <MWAS>

optional arguments:
  -h, --help            show this help message and exit
  -a <str>, --abundance <str>
                        Path to the abundance file generated by MetaCAT's abundance.
                        A header line "Abundance<tab>ID1<tab>ID2 ..." should be present.
  -p <str>, --phenotype <str>
                        Path to the phenotype file.
                        A header line "ID<tab>Phenotype1<tab>Phenotype2 ..." should be present.
                        Each phenotype is a numeric value.
  -v <str>, --variant <str>
                        Path to the variant file generated by MetaCAT's variant.
                        A header line "Classification<tab>Chromosome<tab>Position<tab>Population<tab>ID1<tab>ID2 ..." should be present.
  -o <str>, --output <str>
                        Path to the output file.
  -c <str>, --covariate <str>
                        Path to the covariate file.
                        A header line "ID<tab>Covariate1<tab>Covariate2 ..." should be present.
  -t <int>, --threads <int>
                        Number of threads.
                        Default: all threads.
  --mwas <str>          Path to the mwas file generated by MetaCAT's mwas.
                        Major alleles will be determined in this file rather than in the variant file.
                        A header line "Classification<tab>Chromosome<tab>Position<tab>Allele<tab>Beta<tab>SE<tab>P<tab>Significance" should be present.
  --alpha <float>       The probability of making the wrong decision when the null hypothesis is true.
                        Default: 0.05.
  --pca-components <int>
                        Number of components of GDM to keep.
                        Default: 5.
  --phenotype-column <int>
                        Column of phenotype.
                        Default: 2.
  --size <float>        Minimum proportion of samples for mwas.
                        Default: 0.50.
```

### plotMWAS

```text
usage: MetaCAT plotMWAS [options] -m <MWAS> -o <MWAS>

optional arguments:
  -h, --help            show this help message and exit
  -m <str>, --mwas <str>
                        Path to the mwas file generated by MetaCAT's mwas.
                        Major alleles will be determined in this file rather than in the variant file.
                        A header line "Classification<tab>Chromosome<tab>Position<tab>Allele<tab>Beta<tab>SE<tab>P<tab>Significance" should be present.
  -o <str>, --output <str>
                        Prefix of the output files.
                        "*.qq.pdf": the qq plot for the mwas.
                        "*.manhattan.pdf": the manhattan plot for the mwas.
  --variant-annotation <str>
                        Path to the variant annotation file.
                        A header line "Classification<tab>Chromosome<tab>Position<tab>Annotation ..." should be present.
  --alpha <float>       The probability of making the wrong decision when the null hypothesis is true.
                        Default: 0.05.
  --max-negative-log-p <float>
                        Maximum value of -log10(P).
                        Default: 20.
  --qq-width <float>    Width of the QQ plot in inches.
                        The value should be a positive float.
                        Default: 3.0
  --qq-height <float>   Height of the QQ plot in inches.
                        The value should be a positive float.
                        Default: 3.0
  --manhattan-width <float>
                        Width of the QQ plot in inches.
                        The value should be a positive float.
                        Default: 9.0
  --manhattan-height <float>
                        Height of the QQ plot in inches.
                        The value should be a positive float.
                        Default: 3.0
```

### checkm2

```text
usage: MetaCAT checkm2 [options] -f <FASTAs> -o <CHECKM2>

optional arguments:
  -h, --help            show this help message and exit
  -f <str> [<str> ...], --fasta <str> [<str> ...]
                        Path to the fasta files or directories containing them.
  -o <str>, --output <str>
                        Path to the output file.
  -t <int>, --threads <int>
                        Number of threads.
                        The value should be a positive integer.
                        Default: all threads.
  --checkm2 <str>       Path to the "checkm2".
                        Default: environmental checkm2.
  --fasta-suffix <str>  If directories are provided with "-f|--fasta", only files with the specified suffix will be selected.
                        Default: fasta.
```

### gtdbtk

```text
usage: MetaCAT gtdbtk [options] -f <FASTAs> -o <GTDBTk>

optional arguments:
  -h, --help            show this help message and exit
  -f <str> [<str> ...], --fasta <str> [<str> ...]
                        Path to the fasta files or directories containing them.
  -o <str>, --output <str>
                        Path to the output file.
  --fasta-suffix <str>  If directories are provided with "-f|--fasta", only files with the specified suffix will be selected.
                        Default: fasta.
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

### abundance

|Abundance|*ID1*|*ID2*|...|
|:-------:|:---:|:---:|:-:|
|   ...   | ... | ... |...|

Description:

The relative abundance file generated by the `abundance` command in `MetaCAT`.

*`ID1`*, *`ID2`*, etc. represent the basenames of the input BAM files, excluding file extensions.

### annotation

|Cluster ID|Classification|
|:--------:|:------------:|
|   ...    |     ...      |

Description:

The annotation file generated by the `representative` command in `MetaCAT`.

### checkm2

|Name|Completeness|Contamination|...|
|:--:|:----------:|:-----------:|:-:|
|... |    ...     |     ...     |...|

Description:

The genome quality assessments generated by the `checkm2` command in `MetaCAT`.

Entries in the `Name` column represent the basenames of the input FASTA files, excluding file extensions.

### covariate

|ID |*Covariate1*|*Covariate2*|...|
|:-:|:----------:|:----------:|:-:|
|...|    ...     |    ...     |...|

Description:

The covariate file used by the `mwas` command in `MetaCAT`.

Entries in the `ID` column represent sample IDs.

Missing values should be indicated using **NA** (case-sensitive).

### ground-truth

|Sequence ID|Genome ID|Length|
|:---------:|:-------:|:----:|
|    ...    |   ...   | ...  |

Description:

The ground-truth file used by the `benchmarkGT` command in `MetaCAT`.

### group

|ID |Group|
|:-:|:---:|
|...| ... |

Description:

The group file used by the `abundanceTest` command in `MetaCAT`.

Entries in the `ID` column represent sample IDs.

Missing values should be indicated using **NA** (case-sensitive).

### gtdbtk

|user_genome|classification|...|
|:---------:|:------------:|:-:|
|    ...    |     ...      |...|

Description:

The genome classification generated by the `gtdbtk` command in `MetaCAT`.

Entries in the `user_genome` column represent the basenames of the input FASTA files, excluding file extensions.

### mapping

|Sequence ID|Cluster ID|
|:---------:|:--------:|
|    ...    |   ...    |

Description:

The mapping between each sequence ID and cluster ID generated by the `cluster` and `representative` commands in `MetaCAT`.

### mwas

|Classification|Chromosome|Position|Allele|Beta|SE | P |Significance|...|
|:------------:|:--------:|:------:|:----:|:--:|:-:|:-:|:----------:|:-:|
|     ...      |   ...    |  ...   | ...  |... |...|...|    ...     |...|

Description:

The mwas results generated by the `mwas` command in `MetaCAT`.

Entries `+` and `−` in the `Significance` column indicate significant and non-significant association, respectively.

### phenotype

| ID |Phenotype|
|:--:|:-------:|
|... |   ...   |

Description:

The phenotype file used by the `mwas` command in `MetaCAT`.

Entries in the `ID` column represent sample IDs.

Binary traits should be encoded using only 0 and 1.

Missing values should be indicated using **NA** (case-sensitive).

### variant

|Classification|Chromosome|Position|Population|*ID1*|*ID2*|...|
|:------------:|:--------:|:------:|:--------:|:---:|:---:|:-:|
|     ...      |   ...    |  ...   |   ...    | ... | ... |...|

Description:

The variants generated by the `variant` command in `MetaCAT`.

*`ID1`*, *`ID2`*, etc. represent the basenames of the input BAM files, excluding file extensions.
