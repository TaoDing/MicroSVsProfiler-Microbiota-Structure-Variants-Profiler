# Installation and software requirements

## 1. Tested computing environment
| Component | Original environment |
| --- | --- |
| OS | CentOS Linux 7.9.2009 (Core), 64-bit |
| CPU | Intel Xeon Gold 6130 @ 2.10 GHz |
| RAM | 512 GB |

These values describe the historical analysis environment, **not minimum system requirements**.

To record the environment on a new server:

```bash
cat /etc/redhat-release
uname -m
lscpu | grep -E 'Model name|^CPU\(s\):'
free -h
```
## 2. Command-line tools
The workflow uses the following software:
| Software | Version |
| --- | --- |
| KneadData | 0.7.5 |
| Trimmomatic | 0.39 |
| Bowtie2 | 2.2.3 |
| seqtk | 1.3 |
| Kraken2 | 2.1.2 |
| Bracken | 2.6.2 |
| MetaPhlAn | 3.0 |
| HUMAnN | 3.0 |
| SGV-Finder | 1.0 |
| Bracken | 2.6.2 |
## 3. Installation instructions

### 3.1 GEM Mapper
```bash
# Extract a bzip2-compressed tar archive
tar -xjf GEM-binaries.tbz2

# Add extracted executable directory to PATH
export PATH="/absolute/path/to/GEM-binaries:$PATH"
```
### 3.2 Python environment
Created a Python 3.8.1 environment:
```bash
conda create -n python3.8.1 python=3.8.1
conda activate python3.8.1
```
### 3.3 SGV-Finder and reference files
- Source repository: https://github.com/segalab/SGVFinder
- Data archive: https://zenodo.org/record/3237975/files/DataFiles.tar.gz
```bash
# Download the tools and databases, and install the tools
wget https://zenodo.org/record/3237975/files/DataFiles.tar.gz
git clone https://github.com/segalab/SGVFinder.git

# Inspect the archive before extracting
tar -tzf DataFiles.tar.gz | head

# Follow SGV-Finder's instructions for placing reference data.
cd SGVFinder
python setup.py build
```
### 3.4 KneadData
```bash
conda create -n kneaddata
conda activate kneaddata
conda install -c biobakery kneaddata
mkdir -p kneaddata_database
cd kneaddata_database
kneaddata_database --download human_genome bowtie2 ./
```
### 3.5 seqtk
```bash
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
```
Add the compiled executable to `PATH` if necessary, and record the Git revision used.
### 3.6 HUMAnN 3.0
```bash
conda create --name mpa -c bioconda python=3.7
conda activate mpa
conda install humann
humann_test
```
**Verify:** This does not pin HUMAnN 3.0 or the MetaPhlAn database version.
### 3.7 Kraken2
```bash
conda install -y kraken2
```
**Verify:** Pin Kraken2 version `2.1.2` if reproducing the reported analysis and document the Kraken2 reference database separately.
### 3.8 GraPhlAn
```bash
conda install graphlan
```
## 4. R dependencies
```r
library(data.table)
library(flashClust)
library(ppcor)
library(gplots)
library(cowplot)
library(ggplot2)
library(plyr)
```