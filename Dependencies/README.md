# MetaCAT for Ubuntu

## Install a new Ubuntu on Windows and some dependencies

```text
wsl --install -d Ubuntu-24.04
sudo apt update
sudo apt install -y python3.12-venv libgomp1
```

## Install MetaCAT

```text
python3 -m venv metacat
# All dependencies will be installed or linked into "metacat/bin/". #
source metacat/bin/activate
# Install MetaCAT using pip #
pip3 install https://github.com/liu-congcong/MetaCAT/releases/download/vx.x.x/metacat-x.x.x-py3-none-any.whl
deactivate
```

### Install dependencies used by CheckM2 and GTDB-Tk

```text
for i in diamond FastTree FastTreeMP guppy hmmalign hmmsearch pplacer prodigal skani
do
wget "https://raw.githubusercontent.com/liu-congcong/MetaCAT/main/Dependencies/${i}" && \
mv ${i} metacat/bin/ && \
chmod 755 metacat/bin/${i}
done
```

## Install CheckM2

```text
python3 -m venv checkm2
source checkm2/bin/activate
pip3 install setuptools wheel pandas numpy scipy requests tqdm lightgbm scikit-learn==1.6.1 tensorflow==2.17
git clone --recursive https://github.com/liu-congcong/CheckM2.git checkm2-src && \
cd checkm2-src && \
python3 setup.py install && \
cd .. && \
rm -rf checkm2-src
deactivate
```

### Download database from [Zenodo](https://zenodo.org/records/14897628)

```text
wget "https://zenodo.org/records/14897628/files/checkm2_database.tar.gz" && \
tar xvf checkm2_database.tar.gz && \
rm checkm2_database.tar.gz && \
mv CheckM2_database/uniref100.KO.1.dmnd checkm2 && \
rm -rf CheckM2_database && \
checkm2/bin/checkm2 database --setdblocation checkm2/uniref100.KO.1.dmnd
```

### Link CheckM2 to "metacat/bin/"

```text
ln -s $(pwd)/checkm2/bin/checkm2 metacat/bin/checkm2
```

## Install GTDB-Tk

```text
python3 -m venv gtdbtk
source gtdbtk/bin/activate
pip3 install gtdbtk
deactivate
```

### Download database from [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing-gtdbtk-reference-data)

```text
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz && \
tar xvf gtdbtk_r226_data.tar.gz && \
rm gtdbtk_r226_data.tar.gz && \
mv release226 gtdbtk && \
echo "export GTDBTK_DATA_PATH=\"$(pwd)/gtdbtk/release226\"" >> ~/.bashrc
source ~/.bashrc
```

### Link GTDB-Tk to "metacat/bin/"

```text
ln -s $(pwd)/gtdbtk/bin/gtdbtk metacat/bin/gtdbtk
```

### Run MetaCAT

```text
source metacat/bin/activate
MetaCAT ...
deactivate
```
