# MetaCAT

Metagenome Clustering and Association Tool.

## Installation

**MetaCAT uses the following dependencies and assumes they are on your system path.**

[`CheckM2`](https://github.com/chklovski/CheckM2)
[`FragGeneScan`](https://sourceforge.net/projects/fraggenescan/)
[`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk)
[`Hmmer`](http://hmmer.org)
[`Mash`](https://github.com/marbl/Mash)

`FragGeneScan` (Linux version), `Hmmer` (Linux version) and `Mash` (Linux version) are also available in MetaCAT.

Install `MetaCAT`.

### v1.0.0

```TEXT
pip3 install https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.0/metacat-1.0.0-py3-none-any.whl
```

### v1.0.2

```TEXT
pip3 install https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.2/metacat-1.0.2-py3-none-any.whl
```

### GPU version for MetaCAT

To enable the `GPU` for MetaCAT, you need to install [`CuPy`](https://cupy.dev).

```TEXT
pip3 install cupy-cudaXXX
```

## Benchmarks

[HERE](https://github.com/liu-congcong/MetaCAT/tree/main/Benchmarks).

## Documentation

[HERE](https://github.com/liu-congcong/MetaCAT/tree/main/Documentation).

## Examples

[HERE](https://github.com/liu-congcong/MetaCAT/tree/main/Examples).

## Datasets for Review

[HERE](https://github.com/liu-congcong/MetaCAT/tree/main/Datasets).

## Updates

### v1.0.0

The first release.

### v1.0.1

Add `--bam-suffix` in `coverage` command.

Add `--fasta-suffix` in `cluster2mapping` command.

Add `--input-suffix` in `indexBam` command.

Add `--fasta-suffix` in `representative` command.

Add `--bam-suffix` and `--mapping-suffix` in `abundance` command.

Add `--input-suffix` in `variant` command.

Add `--input-suffix` in `checkm2` command.

Add `--input-suffix` in `gtdbtk` command.

### v1.0.2

Fixed a bug in the `checkm2` command that occurred when handling a large number of input files.
