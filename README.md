# Sample GitLab Project

This sample project shows how a project in GitLab looks for demonstration purposes. It contains issues, merge requests and Markdown files in many branches,
named and filled with lorem ipsum.

You can look around to get an idea how to structure your project and, when done, you can safely delete this project.

[Learn more about creating GitLab projects.](https://docs.gitlab.com/ee/gitlab-basics/create-project.html)


# Snakemake Version of the Microbiome Pipeline


## Prerequisites

### Activate conda environment

To activate the conda environment for st_microbiome, please run the following code in the sHPC shell:

```

ml purge && ml Anaconda3 && conda activate /projects/site/pred/ngs/envs/st_microbiome

```

### Create new conda environment

It may be the case that you need to create a new conda environment. To do this please run the following code:

```

ml purge && ml Anaconda3
conda env create -p </path/to/your/preferred/destination/folder/st_microbiome> -f environment.yml 
conda activate </path/to/your/preferred/destination/folder/st_microbiome> 


# without environment.yml

conda create -y -f -p /projects/site/pred/ngs/envs/st_microbiome python=3.9 mamba
conda activate /projects/site/pred/ngs/envs/st_microbiome
mamba install snakemake samtools multiqc cutadapt umi_tools 10x_bamtofastq fastp kraken2

```

### Kraken Database

Building the Kraken standard database requires a lot of memory and disk space.
See manual https://ccb.jhu.edu/software/kraken/MANUAL.html#standard-kraken-database


The pipeline is using the "Standard plus protozoa & fungi" database (PlusPF, 53Gb archive size) provided by Ben Langmead et al. 
See https://benlangmead.github.io/aws-indexes/k2

```bash

wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230314.tar.gz

```

### How to run

```bash

ml purge && ml Anaconda3 && conda activate /projects/site/pred/ngs/envs/st_microbiome


snakemake

conda deactivate
```



### Development

```bash

ml purge && ml Anaconda3 && conda activate /projects/site/pred/ngs/envs/st_microbiome


snakemake --cores 4 output/{CRC_16,OSCC_2}_possorted_genome_bam_unm_srt.bam
snakemake --cores 4 output/{OSCC_2,CRC_16}_possorted_genome_bam_bamtofastq.done

conda deactivate


```