#!/bin/bash

N_CORES_DEFAULT=8

function usage()
{
  echo ""
  echo "1 input arguments required:"
  echo "1. input config.yaml file"
  echo ""
  echo "1 optional argument"
  echo "1. maximal numbers of cores to use (default $N_CORES_DEFAULT)"
  echo "Contact roland.schmucki@roche.com / tel 71330"
  echo ""
  echo ""
}

get_abs_filename() {
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

if [[ $# -lt 1 ]];then usage && exit 1; fi

# Define variable (pRED)
CONFIG_FILE=$(get_abs_filename $1)
SNAKE_FILE=./Snakefile
GRAPH_FILE=rulegraph.pdf
CONDA_ENV=/projects/site/pred/ngs/envs/st_microbiome

if [[ $# -gt 1 ]]; then
  N_CORES=$2
else
  N_CORES=$N_CORES_DEFAULT
fi

if [ ! -e $CONFIG_FILE ] ; then usage && exit 1; fi

# Activate environment (pRED)
ml purge && ml Anaconda3 && conda activate $CONDA_ENV

# Create workflow diagram
snakemake --snakefile $SNAKE_FILE --rulegraph all --configfile $CONFIG_FILE | dot -Tpdf > rulegraph.pdf

# Run the workflow on the cluster
snakemake --snakefile $SNAKE_FILE \
    --configfile $CONFIG_FILE \
    --cores $N_CORES \
    --notemp \
    --latency-wait 6 \
    --rerun-incomplete \
    --keep-going \
    --verbose
