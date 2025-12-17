# MetaCAT

Metagenome Clustering and Association Tool.

## Dependencies

- [`CheckM2`](https://github.com/chklovski/CheckM2)

- [`CuPy`](https://cupy.dev)

- [`FragGeneScan`](https://sourceforge.net/projects/fraggenescan/)

- [`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk)

- [`Hmmer`](http://hmmer.org)

- [`Mash`](https://github.com/marbl/Mash)

- [`FastANI`](https://github.com/ParBLiSS/FastANI)

- [`skani`](https://github.com/bluenote-1577/skani)

### Dependencies used for clutering sequneces

`FragGeneScan`, `hmmpress` and `hmmsearch` (from the `Hmmer` tool) are used in the `seed` command of MetaCAT.

MetaCAT uses built-in versions of `FragGeneScan`, `hmmpress`, and `hmmsearch` by default, \
if you are running on a supported platform (i.e., x86_64 Linux or arm macOS).

Otherwise, you need to compile them manually and ensure they are callable by MetaCAT.

You can specify their paths using `--fraggenescan`, `--hmmpress` and `--hmmsearch` options.

### Dependencies used for estimating cluster quality

`CheckM2` is used in the `checkm2` command of MetaCAT.

By default, MetaCAT uses the `checkm2` from the environment.

You can specify its path using the `--checkm2` option.

`CheckM2` and its database (*.dmnd) must be installed manually, for examples:

```bash
pip3 install checkm2 # CheckM2 v1.0.1 has been tested #
checkm2 database --download --path /path/to/database
checkm2 database --setdblocation /path/to/*.dmnd
```

For more details, please visit the [`CheckM2`](https://github.com/chklovski/CheckM2) repository.

### Dependencies used for classifying clusters

`GTDB-Tk` is used in the `gtdbtk` command of MetaCAT.

By default, MetaCAT uses `gtdbtk` from the environment.

You can specify its path using the `--gtdbtk` option.

`GTDB-Tk` must be installed manually, for example:

```bash
pip3 install gtdbtk # GTDB-TK v2.3.2 has been tested #
```

Please note that you also need to install the `GTDB-Tk` database and other required dependencies.

```bash
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar xvf gtdbtk_r226_data.tar.gz
mv release226 ~/gtdb-226
echo "export GTDBTK_DATA_PATH="$HOME/gtdb-226"" >> ~/.bashrc
```

Use the following command to verify the installation:

```bash
gtdbtk --check_install
```

For more details, please visit the [`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk) repository.

### Dependencies used for selecting representatives

`Mash` or `FastANI` or `Skani` is used in the `representative` command of MetaCAT.

By default, MetaCAT uses the built-in `Mash` binary \
if you are running on a supported platform (i.e., x86_64 Linux or arm macOS).

Otherwise, you need to compile `Mash` manually and ensure it is callable by MetaCAT.

You can specify its path using the `--mash` option.

Additionally, we provide `FastANI` and `Skani` as alternatives (`--engine fastANI|skani`) in `representative` command.

To use `FastANI` or `Skani`, you need to compile it manually.

MetaCAT uses `FastANI` or `Skani` from the environment by default.

You can specify its path using the `--fastani` or `--skani` option.

## Install MetaCAT

```TEXT
pip3 install https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.3/metacat-1.0.3-py3-none-any.whl
```

### GPU version for MetaCAT

To enable the `GPU` for MetaCAT, you need to install [`CuPy`](https://cupy.dev).

```TEXT
pip3 install cupy-cudaXXX
```

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

### v1.0.3

Fixed a bug that could potentially cause a segmentation fault.

Support `skani`.

Rewrite indexBam.
