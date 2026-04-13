# Install MetaCAT for Linux 

## 1. Install metacat using python in system environment

```text
python3 -m venv metacat && echo "export PATH=\"\$PATH:$(pwd)/metacat/bin/\"" >> ~/.bashrc
source metacat/bin/activate
pip3 install https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.3/metacat-1.0.3-py3-none-any.whl
deactivate && source ~/.bashrc
```

### Install dependencies

* prodigal, diamond, FastTree, FastTreeMP, guppy, hmmalign, hmmsearch, pplacer, skani
```text
for i in FastTree FastTreeMP guppy hmmalign hmmsearch pplacer skani prodigal diamond
do
wget "https://raw.githubusercontent.com/liu-congcong/MetaCAT/main/Dependencies/${i}"
mv ${i} gtdbtk/bin && chmod 755 gtdbtk/bin/${i}
done
```

* checkm2
```text
python3 -m venv checkm2 && echo "export PATH=\"\$PATH:$(pwd)/checkm2/bin/\"" >> ~/.bashrc
source checkm2/bin/activate
pip3 install setuptools wheel pandas numpy scipy requests tqdm lightgbm scikit-learn==1.6.1 tensorflow==2.17
git clone --recursive https://github.com/liu-congcong/CheckM2.git checkm2-src && cd checkm2-src
python3 setup.py install && cd .. && rm -rf checkm2-src && deactivate && source ~/.bashrc
checkm2 testrun
MetaCAT checkm2 -h
```

* GTDB-Tk
```text
python3 -m venv gtdbtk && echo "export PATH=\"\$PATH:$(pwd)/gtdbtk/bin/\"" >> ~/.bashrc
source gtdbtk/bin/activate
pip3 install gtdbtk
deactivate && source ~/.bashrc
MetaCAT gtdbtk -h
```

## 2. Install dependencies using Conda/Mamba package manager
### Create a virtual environment and activate:
```text
conda create -n MetaCAT python=3.12
conda activate MetaCAT
```
### Install metacat using pip command in a virtual environment:
```text
wget https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.3/metacat-1.0.3-py3-none-any.whl
pip install metacat-1.0.3-py3-none-any.whl
echo "export PATH=\"\$PATH:$(pwd)/metacat/bin/\"" >> ~/.bashrc
source ~/.bashrc && conda activate MetaCAT
```
### Install dependencies using package manager
* prodigal: `pip install prodigal` or `conda install -c conda-forge -c bioconda prodigal`
* diamond: `conda install -c conda-forge -c bioconda diamond`
* HMMER: `pip install HMMER` or `conda install -c bioconda HMMER`
* pplacer: `conda install -c bioconda pplacer`
* Skani: `conda install -c bioconda Skani`
* FastTree and FastTreeMP: `conda install -c conda-forge -c bioconda fasttree`
* Guppy: `conda install -c conda-forge -c bioconda guppy`
* gtdbtk: `pip install gtdbtk` or `conda install -c conda-forge -c bioconda gtdbtk`
* checkm2: `conda install -c conda-forge -c bioconda checkm2`

## 3. Download database from [Zenodo](https://zenodo.org/records/14897628)

```text
wget "https://zenodo.org/records/14897628/files/checkm2_database.tar.gz"
tar xvf checkm2_database.tar.gz && rm checkm2_database.tar.gz
mv CheckM2_database/uniref100.KO.1.dmnd checkm2 && rm -rf CheckM2_database
checkm2 database --setdblocation checkm2/uniref100.KO.1.dmnd
```

## 4. Download database from [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing-gtdbtk-reference-data)

```text
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar xvf gtdbtk_r226_data.tar.gz && rm gtdbtk_r226_data.tar.gz
mv release226 gtdbtk
echo "export GTDBTK_DATA_PATH=\"$(pwd)/gtdbtk/release226\"" >> ~/.bashrc && source ~/.bashrc
gtdbtk test
```
