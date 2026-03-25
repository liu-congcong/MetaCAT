# MetaCAT for Ubuntu

## Install a new Ubuntu on Windows and some dependencies

```text
wsl --install -d Ubuntu-24.04
sudo apt update
sudo apt install -y python3.12-venv libgomp1
```

## Install MetaCAT

```text
python3 -m venv metacat && echo "export PATH=\"$(pwd)/metacat/bin:\$PATH\"" >> ~/.bashrc
source metacat/bin/activate
pip3 install ...
deactivate && source ~/.bashrc
```

## Install CheckM2

```text
python3 -m venv checkm2 && echo "export PATH=\"$(pwd)/checkm2/bin:\$PATH\"" >> ~/.bashrc
source checkm2/bin/activate
pip3 install setuptools wheel pandas numpy scipy requests tqdm lightgbm scikit-learn==1.6.1 tensorflow==2.17
git clone --recursive https://github.com/liu-congcong/CheckM2.git checkm2-src && cd checkm2-src
python3 setup.py install && cd .. && rm -rf checkm2-src && deactivate && source ~/.bashrc
```

### Download database from [Zenodo](https://zenodo.org/records/14897628)

```text
wget "https://zenodo.org/records/14897628/files/checkm2_database.tar.gz"
tar xvf checkm2_database.tar.gz && rm checkm2_database.tar.gz
mv CheckM2_database/uniref100.KO.1.dmnd checkm2 && rm -rf CheckM2_database
checkm2 database --setdblocation checkm2/uniref100.KO.1.dmnd
```

### Install dependencies

```text
for i in diamond prodigal
do
wget "https://raw.githubusercontent.com/liu-congcong/MetaCAT/main/Dependencies/${i}"
mv ${i} gtdbtk/bin && chmod 755 gtdbtk/bin/${i}
done
```

### Run test

```text
checkm2 testrun
```

### Run checkm2 command in MetaCAT

```text
MetaCAT checkm2 ...
```

## Install GTDB-Tk

```text
python3 -m venv gtdbtk && echo "export PATH=\"$(pwd)/gtdbtk/bin:\$PATH\"" >> ~/.bashrc
source gtdbtk/bin/activate
pip3 install gtdbtk
deactivate && source ~/.bashrc
```

### Download database from [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing-gtdbtk-reference-data)

```text
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar xvf gtdbtk_r226_data.tar.gz && rm gtdbtk_r226_data.tar.gz
mv release226 gtdbtk
echo "export GTDBTK_DATA_PATH=\"$(pwd)/gtdbtk/release226\"" >> ~/.bashrc && source ~/.bashrc
```

### Install dependencies

```text
for i in FastTree FastTreeMP guppy hmmalign hmmsearch pplacer prodigal skani
do
wget "https://raw.githubusercontent.com/liu-congcong/MetaCAT/main/Dependencies/${i}"
mv ${i} gtdbtk/bin && chmod 755 gtdbtk/bin/${i}
done
```

### Run test

```text
gtdbtk test
```

### Run gtdbtk command in MetaCAT

```text
MetaCAT gtdbtk ...
```
