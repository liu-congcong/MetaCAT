# MetaCAT

Metagenome Clustering and Association Tool.

## Dependencies

[`CheckM2`](https://github.com/chklovski/CheckM2)
[`CuPy`](https://cupy.dev)
[`FragGeneScan`](https://sourceforge.net/projects/fraggenescan/)
[`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk)
[`Hmmer`](http://hmmer.org)
[`Mash`](https://github.com/marbl/Mash)
[`FastANI`](https://github.com/ParBLiSS/FastANI)
[`skani`](https://github.com/bluenote-1577/skani)

### Dependencies used for clutering sequneces

`FragGeneScan`, `hmmpress` and `hmmsearch` (from the `Hmmer` tool) are used in the `seed` command of MetaCAT.

MetaCAT uses built-in versions of `FragGeneScan`, `hmmpress`, and `hmmsearch` by default, \
if you are running on a supported platform (i.e., x86_64 Linux or arm macOS).

Otherwise, you need to compile them manually and ensure they are callable by MetaCAT.

You can specify their paths using `--fraggenescan`, `--hmmpress` and `--hmmsearch` options.

Note that the built-in binaries are only accessible to the `seed` command and cannot be called by external programs.

### Dependencies used for estimating cluster quality

`CheckM2` is used in the `checkm2` command of MetaCAT.

By default, MetaCAT uses the `checkm2` from the environment.

You can specify its path using the `--checkm2` option.

`CheckM2` and its database must be installed manually.

For more details, please visit the [`CheckM2`](https://github.com/chklovski/CheckM2) repository.

### Dependencies used for classifying clusters

`GTDB-Tk` is used in the `gtdbtk` command of MetaCAT.

By default, MetaCAT uses `gtdbtk` from the environment.

You can specify its path using the `--gtdbtk` option.

`GTDB-Tk` and its database must be installed manually.

For more details, please visit the [`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk) repository.

### Dependencies used for selecting representatives

`Mash` or `FastANI` or `Skani` is used in the `representative` command of MetaCAT.

By default, MetaCAT uses the built-in `Mash` binary \
if you are running on a supported platform (i.e., x86_64 Linux or arm macOS).

Otherwise, you need to compile `Mash` manually and ensure it is callable by MetaCAT.

You can specify its path using the `--mash` option.

Note that the built-in binaries are only accessible to the `representative` command and cannot be called by external programs.

Additionally, we provide `FastANI` and `Skani` as alternatives (`--engine fastANI|skani`) in `representative` command.

To use `FastANI` or `Skani`, you need to compile it manually.

MetaCAT uses `FastANI` or `Skani` from the environment by default.

You can specify its path using the `--fastani` or `--skani` option.

## Install MetaCAT and dependencies

If only the core modules of MetaCAT (e.g., clustering) are required, installing MetaCAT alone is sufficient.

`CheckM2`, `GTDB-Tk`, and their dependencies may require specific Python versions.

The following installation includes MetaCAT and all required dependencies, and has been tested with Python 3.12.

- `MetaCAT` will be installed in `./metacat`
- `CheckM2` will be installed in `./checkm2`
- `GTDB-Tk` will be installed in `./gtdbtk`
- Binary dependencies (`CheckM2` & `GTDB-Tk`) will be installed in `./metacat/bin`

You may adjust the installation paths as needed.

### Install MetaCAT

```text
python3 -m venv metacat
source metacat/bin/activate
pip3 install https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.4/metacat-1.0.4-py3-none-any.whl
deactivate
```

**If CheckM2, GTDB-Tk, and their dependencies are already available in the system PATH, the following steps can be skipped.**

### Install dependencies used by CheckM2 and GTDB-Tk

All dependencies will be installed or linked into "metacat/bin/".

```text
for i in diamond FastTree FastTreeMP guppy hmmalign hmmsearch pplacer prodigal skani
do
wget "https://raw.githubusercontent.com/liu-congcong/MetaCAT/main/Dependencies/${i}"
mv ${i} metacat/bin/
chmod 755 metacat/bin/${i}
done
```

### Install CheckM2

```text
python3 -m venv checkm2
source checkm2/bin/activate
pip3 install setuptools wheel pandas numpy scipy requests tqdm lightgbm scikit-learn==1.6.1 tensorflow==2.17
git clone --recursive https://github.com/liu-congcong/CheckM2.git checkm2-src
cd checkm2-src
python3 setup.py install
cd ..
rm -rf checkm2-src
deactivate
```

### Download CheckM2's database from [Zenodo](https://zenodo.org/records/14897628)

```text
wget "https://zenodo.org/records/14897628/files/checkm2_database.tar.gz"
tar xvf checkm2_database.tar.gz
rm checkm2_database.tar.gz
mv CheckM2_database/uniref100.KO.1.dmnd checkm2/
rm -rf CheckM2_database
checkm2/bin/checkm2 database --setdblocation checkm2/uniref100.KO.1.dmnd
```

### Link CheckM2 to "metacat/bin/"

```text
ln -s $(pwd)/checkm2/bin/checkm2 metacat/bin/checkm2
```

### Install GTDB-Tk

Note that GTDB-Tk is frequently updated and may not be compatible with earlier database releases.

Therefore, each GTDB-Tk version should be used with its corresponding database to ensure compatibility.

```text
python3 -m venv gtdbtk
source gtdbtk/bin/activate
pip3 install gtdbtk==2.7.1
deactivate
```

### Download GTDB-Tk's database from [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing-gtdbtk-reference-data)

```text
wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
tar xvf gtdbtk_r232_data.tar.gz
rm gtdbtk_r232_data.tar.gz
mv release232 gtdbtk/
echo "export GTDBTK_DATA_PATH=\"$(pwd)/gtdbtk/release232\"" >> ~/.bashrc
source ~/.bashrc
```

### Link GTDB-Tk to "metacat/bin/"

```text
ln -s $(pwd)/gtdbtk/bin/gtdbtk metacat/bin/gtdbtk
```

### GPU version for MetaCAT

To enable the **GPU** for MetaCAT, you need to install [`CuPy`](https://cupy.dev).

```text
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

### v1.0.4

Support `GTDB-Tk v2.7.1`.

Minor bug fixes.
