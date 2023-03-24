# Sample GitLab Project

This sample project shows how a project in GitLab looks for demonstration purposes. It contains issues, merge requests and Markdown files in many branches,
named and filled with lorem ipsum.

You can look around to get an idea how to structure your project and, when done, you can safely delete this project.

[Learn more about creating GitLab projects.](https://docs.gitlab.com/ee/gitlab-basics/create-project.html)


# Snakemake Version of the Microbiome Pipeline


## Prerequisites

### Activate conda environment (only for pRED)

To activate the conda environment for st_microbiome, please run the following code in the sHPC shell:

```

ml purge && ml Anaconda3 && conda activate /projects/site/pred/ngs/envs/st_microbiome

```

### Create new conda environment

It may be that you need to create a new conda environment from scratch. To do this please run the following code:

```bash
DIR=/path/to/your/preferred/destination/folder/st_microbiome
NAME=st_microbiome

# load module (only pRED)
ml purge && ml Anaconda3

# without environment.yml file
conda create -y -f -p $DIR/$NAME python=3.9 mamba
conda activate $DIR/$NAME
mamba install -y snakemake samtools multiqc cutadapt umi_tools 10x_bamtofastq fastp kraken2
conda install -y conda-minify -c jamespreed

# then create the environment.yml file
conda-minify --name $DIR/$NAME -f environment.yml 

# with environment.yml file
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

### Workflow Configuration File

In order to run the Snakemake workflow, one has to specify several parameters by a configuration file, e.g. `config.yaml`.

The structure and format of the yaml file is as follows
```
results: 'output'                                                                    # Path to output directory, may not exist
bam_dir: 'st_microbiome'                                                             # Path to input directory containing BAM files (must exist)
samples: ['OSCC_2_possorted_genome_bam', 'CRC_16_possorted_genome_bam']              # List of BAM file names without file extension that are present in the input directory
kraken_threads: 12                                                                   # Number of cores to use for the Kraken classification step
kraken_db: '/projects/site/pred/ngs/pipelines/st_microbiome/kraken2_Standard_PlusPF' # Path to the Kraken database
```

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