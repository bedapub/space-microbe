# 'Spatial Transcriptomics meets Microbiome' Snakemake Workflow

## Description
TBD

## Graph
TBD

## Citation
TBD

## Table of contents
[Quick Start](#quick-start)  
    * [Installation](#installation)  
    * [Test Installation](#test-installation)  
    * [Basic Usage](#basic-usage)  
    * [Main Output](#main-output)  
    * [Environment Creation with Mamba](#environment-creation-with-mamba)  



## Quick Start
### Installation

**Step 1 - Clone the repository**

First, clone the repository to the desired location on your system by using the following command:

```bash
git clone https://github.com/bedapub/xyz.git
```

**Step 2 - Install conda, if not already installed**

Next, ensure that conda is installed on your system, otherwise install using the instructions provided [here](https://developers.google.com/earth-engine/guides/python_install-conda/).

**Step 3 - Create and activate a new conda environment**

Then, create a new conda environment via the following command, replacing `/full/path/to/cloned/repo` with the appropriate path to your cloned repository:

```bash
conda env create -n st_microbiome_env --file /full/path/to/cloned/repo/environment.yaml
```

And activate the environment by executing:

```bash
conda activate st_microbiome_env
```

If you have trouble creating the environment using the above commands, you can alternatively follow the instructions [here](#environment-creation-with-mamba).

**Step 4 - Download Kraken database**

Building the Kraken standard database requires a lot of memory and disk space.
See manual https://ccb.jhu.edu/software/kraken/MANUAL.html#standard-kraken-database

The Snakemake workflow is using the "Standard plus protozoa & fungi" database (PlusPF, 53Gb archive size) provided by Ben Langmead et al. 
See https://benlangmead.github.io/aws-indexes/k2

```bash
mkdir -p k2_pluspf_20230314 && cd k2_pluspf_20230314
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230314.tar.gz
tar xvzf k2_pluspf_20230314.tar.gz
```

**Step 5 - Prepare input data**

All input data, ie `BAM` files, must be located in local directories, specified by an input file, the so called _"file list"_.
The _"file list"_ is a tab-delimited text file and contains two columns: 
the first column denotes a sample alias (no white space, slash, etc. allowed) and the second column contains the paths to the corresponding `BAM` files.
Then, Snakemake will read this "file list" and process the `BAM` files according to the workflow protocol.

Example file
```
sampleA    /path/to/sampleA/possorted_genome_bam.bam
sampleB    /path/to/sampleB/possorted_genome_bam.bam
sampleC    /path/to/sampleC/possorted_genome_bam.bam
```


### Test Installation

In order to test that the Snakemake workflow can be run on your data, navigate to your cloned repository using `cd` and download the `BAM` files required for testing via the following command:

```bash
wget <url/to/repo/with/test_dataset.tar.gz>
tar xvzf test_dataset.tar.gz
```

Then, run the wrapper script `run.py` which will create a new `config.yaml` and launch the Snakemake workflow either locally (use `--cores <int>` option) or submitted to the cluster (use `--profile <path>` option). In the example below, Snakemake will be run locally and by using 4 cores for data processing.

```bash
conda activate st_microbiome_env

python run.py --file-list <path to FILE_LIST> \
              --outdir <path to output folder> \
              --kraken-db <path to KRAKEN_DB> \
              --cores 4
```

### Basic Usage

For basic usage, first activate the conda environment and then run the wrapper script with the appropriate arguments, i.e. path to the input file list, output directory, path to kraken database, and use the optional parameters as required.

By the parameter `--profile <path>` the workflow will be submitted to the cluster with runnin up to `--jobs <int>` in parralle. However, if the option `--cores <int>` is used then the workflow will be executed locally.

```
conda activate st_microbiome_env

python run.py --file-list <path to FILE_LIST> \
              --outdir <path to output folder> \
              --kraken-db <path to KRAKEN_DB> \
              [--kraken-threads KRAKEN_THREADS] \
              [--cores CORES] \
              [--jobs JOBS] \
              [--profile <path to configuration file for cluster PROFILE>]
```

The wrapper script has the following parameters:

```
python run.py --help

optional arguments:
  -h, --help            show this help message and exit
  --file-list FILE_LIST, -f FILE_LIST
                        Path to file list (tabular wiht Sample-ID and BAM-filepath)
  --outdir OUTDIR, -o OUTDIR
                        Path to output directory, may not exist
  --kraken-db KRAKEN_DB, -d KRAKEN_DB
                        Path to the Kraken database
  --kraken-threads KRAKEN_THREADS, -k KRAKEN_THREADS
                        Number of cores to use for the Kraken classification step
  --cores CORES, -c CORES
                        Number of cores to use for local run
  --jobs JOBS, -j JOBS  Number of jobs for running on the cluster
  --profile PROFILE, -p PROFILE
                        Path to cluster profile, if omitted Snakemake will be run locally
```                        

### Main Output

TBD


### Environment Creation with Mamba

This method may be quicker than the one described above. Here we create the conda environment with `mamba` and without the provided `environment.yml` file. At the end, we use `conda-minify` to create a new `environment.yml` file.

```bash
DIR=/path/to/your/preferred/destination/folder/st_microbiome_env

conda create -y -f -p $DIR -c conda-forge python=3.9 mamba
conda activate $DIR
mamba install -c bioconda -c conda-forge -c jamespreed -y conda-minify snakemake samtools multiqc cutadapt umi_tools 10x_bamtofastq fastp kraken2

conda-minify --name $DIR -f environment.yml 
```

----------------------------------
# INTERNAL VERSION (may be deleted)
# Snakemake Version of the Pipeline 'Spatial Transcriptomics meets Microbiome'


## Prerequisites

### Activate conda environment (only for pRED)

To activate the conda environment for st_microbiome, please run the following code in the sHPC shell:

```
ml purge && ml Anaconda3 && conda activate /projects/site/pred/ngs/envs/st_microbiome
```

### Create new conda environment

It may be that you need to create a new conda environment from scratch. To do this please run one of the following code (`mamba` may be the fastest):

```bash
DIR=/path/to/your/preferred/destination/folder/st_microbiome
NAME=st_microbiome

# load module (only pRED)
ml purge && ml Anaconda3

# create conda env with mamba and without environment.yml file
conda create -y -f -p $DIR/$NAME -c conda-forge python=3.9 mamba
conda activate $DIR/$NAME
mamba install -c bioconda -c conda-forge -c jamespreed -y conda-minify snakemake samtools multiqc cutadapt umi_tools 10x_bamtofastq fastp kraken2

# then create the environment.yml file
conda-minify --name $DIR/$NAME -f environment.yml 

# create conda env with environment.yml file
conda env create -p $DIR/$NAME -f environment.yml 
conda activate $DIR/$NAME
```

### Kraken Database

Building the Kraken standard database requires a lot of memory and disk space.
See manual https://ccb.jhu.edu/software/kraken/MANUAL.html#standard-kraken-database


The pipeline is using the "Standard plus protozoa & fungi" database (PlusPF, 53Gb archive size) provided by Ben Langmead et al. 
See https://benlangmead.github.io/aws-indexes/k2

```bash
mkdir -p k2_pluspf_20230314 && cd k2_pluspf_20230314
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230314.tar.gz
tar xvfz k2_pluspf_20230314.tar.gz
```

### Input Data (BAM files)

All input data, ie `BAM` files, must be located in local directories, specified by an input file, the so called "file list".
The "file list" is a tab-delimited text file and contains two columns: 
the first column denotes a sample alias (no white space, slash, etc. allowed) and the second column contains the paths to the corresponding `BAM` files.
Then, Snakemake will read this "file list" and process the `BAM` files according to the workflow protocol.


### Workflow Configuration File

In order to run the Snakemake workflow, one has to specify several parameters by a configuration file, e.g. `config.yaml`.

The structure and format of the yaml file is as follows (white-space delimited).
```
file_list: 'my/input/file_list.txt'                                                  # Path to file list (tabular wiht Sample-ID and BAM-filepath)
results: 'output'                                                                    # Path to output directory, may not exist
kraken_threads: 12                                                                   # Number of cores to use for the Kraken classification step
kraken_db: '/projects/site/pred/ngs/pipelines/st_microbiome/kraken2_Standard_PlusPF' # Path to the Kraken database
```

The workflow will process all BAM files that are in the input directory `bam_dir` (must have extenstion `.bam`).

## How to run (pRED only)

There are two ways to run the workflow:
1. locally in the terminal
2. submit to a cluster

To run it locally, use the following command where the first input argument is the path to the workflow configuration file, `config.yaml` (see above) and the second, optional argument is to specify the maximal number of cores to use.

```bash
./run_locally.sh <config.yaml> [cores]
```

If the workflow should be executed on a cluster, then use the following command.

```bash
./run_hpc.sh <config.yaml>
```



### Development

```bash
LSF_PROFILE=/apps/rocs/etc/apps/snakemake/lsf/v1.4_memfix
SNAKE_FILE=Snakefile
CONFIG_FILE=config.yaml
JOBS=100

ml purge && ml Anaconda3 && conda activate /projects/site/pred/ngs/envs/st_microbiome

snakemake --snakefile $SNAKE_FILE --rulegraph all --configfile $CONFIG_FILE | dot -Tpdf > rulegraph.pdf

# Submit to cluster
snakemake --snakefile $SNAKE_FILE \
    --configfile $CONFIG_FILE \
    --jobs $JOBS \
    --profile $LSF_PROFILE \
    --notemp \
    --latency-wait 6 \
    --rerun-incomplete \
    --keep-going \
    --verbose
    
    
# Run it locally
snakemake --cores 4 

conda deactivate

```