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

The structure and format of the yaml file is as follows.
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