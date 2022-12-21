#!/usr/bin/env bash

# This creates a job template for LSF on sHPC
# Program: bamtofastq
# Type: Array
# script by galvezbe @ROCHE microbiome-squad - modified by anzboecs

# Call: bash extract_unmapped_bam2fastq_template.sh FILELIST SETTINGS OUTDIR

# Input file list should be in Format with no header:
# Name possorted_genome_bam.bam

if [ "$#" -ne 3 ]; then
        echo "Incorrect number of input files."
        exit 1
else
        FILELIST=$1
        SETTINGS=$2
        OUTDIR=$3
        # extend OUTDIR path if relative
        if [[ ! "$OUTDIR" = /* ]]; then
                OUTDIR=`pwd -P`/$OUTDIR
        fi
        if [[ ! "$FILELIST" = /* ]]; then
                FILELIST=`pwd -P`/$FILELIST
        fi
        echo "Input file list: $FILELIST"
        echo "Input settings: $SETTINGS"
        echo "Output directory: $OUTDIR"
fi
#==================================== arguments ==========================#

FILELIST_BASE=`basename ${FILELIST/.txt/}`
OUTFILE=$OUTDIR/bsub/extract_unmapped.$FILELIST_BASE.run.bsub
NFILES=`cat $FILELIST | wc -l`

if [ `grep -c '$ARVDIR' $FILELIST` -ne 0 ]; then
ARVLOAD="export ARVDIR=/local/\$USER-arv/\${LSB_JOBID}_\${LSB_JOBINDEX}/ \n
function finish { \n
  ### Unmount Arvados on exit \n
  ml arvados
  arv-mount --unmount \$ARVDIR/ \n
  rmdir \$ARVDIR \n
  echo 'Exit from finish' \n
} \n
trap finish EXIT \n
\n
source /apps/rocs/init.sh \n
ml arvados \n
\n
mkdir -p \$ARVDIR \n
arv-mount --read-only \$ARVDIR/"
else
ARVLOAD=''
fi

mkdir -p $OUTDIR $OUTDIR/bsub/

# printing lsf job header 
echo "#!/bin/bash
#BSUB -J bash_pipeline.$FILELIST_BASE[1-$NFILES]       # Job array name and job indexes
#BSUB -n 6                 # number of tasks per job, Use 4, 8, 12, 24 
#BSUB -q preempt                  # Select queue
#BSUB -M 30000                    # RAM (i.e 10Gb per job)
#BSUB -R \"span[hosts=1]\"        # Allocate all tasks in 1 host
#BSUB -o `realpath $OUTDIR`/bsub/extract_unmapped.$FILELIST_BASE.%J_%I.bsub.stdout   # Output files
#BSUB -e `realpath $OUTDIR`/bsub/extract_unmapped.$FILELIST_BASE.%J_%I.bsub.stderr   # Error files

`echo -e $ARVLOAD`

## parsing file list 

THREADS=6
LINE=\`sed -n \${LSB_JOBINDEX}p `realpath $FILELIST`\`
NAME=\`echo \$LINE | awk '{print \$1}'\`
BAM=\`echo \$LINE | awk '{print \$2}' | sed \"s,[$]ARVDIR,\$ARVDIR,g\"\`
OUTDIR=`realpath $OUTDIR`

echo "[\$NAME] Processing \$BAM "

##  Loading tool enviroment and sHPC modules
" > $OUTFILE

cat $SETTINGS >> $OUTFILE


echo -e "\nRun: bsub < $OUTFILE\n"
