# cov2seq

Automation and QC tool for Nanopore sequencing experiments of SARS-CoV2.

## Installation

Cov2seq needs to be installed on a UNIX based operating system.

It is highly recommended to install cov2seq and all its dependencies in an anaconda3 environment. To setup such an environment and install most of the dependencies, ake use of the environment.yml file:

```
conda env create -f environment.yml
conda activate cov2seq
```

Now clone and install the ARTIC fieldbioinformatics software suit. In our lab, we currently use a slightly modified version of the v1.2.1 branch that is using a different approach to normalization suitable for sequencing SARS-CoV2 with rapid kits. If you want to use the official distribution you can do so, they are fully compatible.

```
git clone https://github.com/MarkusHaak/fieldbioinformatics.git
cd fieldbioinformatics
python setup.py install
```

Cov2seq additionally requires some binaries and files from different git repositories. You can clone them to whatever location you want, but need to change the corresponding entry in the configuration file as explained later. In order to let cov2seq find its location by default, setup the default directory tree as follows:

```
mkdir -p ~/cov2seq/nanopore/samples
mkdir -p ~/cov2seq/results/samples
mkdir -p ~/cov2seq/illumina/samples
mkdir -p ~/cov2seq/sanger/samples
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
