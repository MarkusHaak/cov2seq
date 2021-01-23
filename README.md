# cov2seq

Automation and QC tool for Nanopore sequencing experiments of SARS-CoV2.

## Installation

Cov2seq needs to be installed on a UNIX based operating system.

It is highly recommended to install cov2seq and all its dependencies in an anaconda3 environment with python version 3.8:

```
conda create -n cov2seq python=3.8
conda activate cov2seq
conda config --add channels conda-forge bioconda
```

Then install the dependencies using conda and pip:

```
conda install mafft iqtree augur bedtools samtools longshot
```

At the point of writing, setuptools is complaining that Nextstrain augur requires bioconda version 1.76 when running the installation of cov2seq, but version 1.78 is installed by conda. Use pip to fix this issue if you encounter it:

```
pip uninstall bioconda
pip install bioconda==1.76
```

Cov2seq additionally requires some binaries and files from different git repositories. You can clone them to whatever location you want, but need to change the corresponding entry in the configuration file as explained later. In order to let cov2seq find its location by default, setup the default directory tree as follows:

```
cd
mkdir -p cov2seq/nanopore/samples
mkdir -p cov2seq/results/samples
mkdir -p cov2seq/illumina/samples
mkdir -p cov2seq/sanger/samples
```

Then change into the root of the directory tree and clone the nextstrain ncov github repository:

```
cd ~/cov2seq
git clone https://github.com/nextstrain/ncov.git
sed -i 's/>MN908947/>MN908947.3/' ncov/defaults/reference_seq.fasta
```

Keep this repository up to date by running `git pull` on a regular basis to receive the latest nextstrain clade information. The reference sequence information under ncov/defaults/reference_seq.* are used in the pipeline. Possibly, the fasta header in ncov/defaults/reference_seq.fasta needs to be changed from `>MN908947` to `>MN908947.3`, since the latter is used as the reference id in the ARTIC pipeline.

Cov2seq also needs the information about primer schemes used in the artic pipeline. Either skip the following step and change the corresponding path in the config file to the location of your primer_schemes directory, or clone the artic-ncov2019 repository to the cov2seq root directory:

```
git clone https://github.com/artic-network/artic-ncov2019.git
```

Additionally, the SnpEff binaries are required to annotate vcf files. Download and unzip the binaries and create a MN908947.3 reference:

```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
./scripts/buildDbNcbi.sh MN908947.3
cd ..
```

Now clone the cov2seq repository (to any location you like) and install cov2seq:

```
git clone https://github.com/MarkusHaak/cov2seq.git
cd cov2seq
python setup.py install
```

The last lines output to the console after a successful installation should indicate the location of the basic configuration file used by cov2seq. You can copy this file to your home directory and make changes that will be prioritized over the entries in the default configuration file. You can also specify the location of any configuration file when running cov2seq or alternatively change all options with their respective command line arguments.
