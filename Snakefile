"""
Snakemake workflow for ST microbiome data analysis pipeline
Spatial Transcriptomics meets Microbiome
"""
import sys
import pandas as pd
import subprocess # to get git hash
from os.path import dirname, join


# ----------------------------------------------------------
# Define the pipeline configuration file
configfile: 'config.yaml'


# ----------------------------------------------------------
# Define variables
OD = config['results']
INDIR = config['input_directory']

# ----------------------------------------------------------
"""
main rule
"""
rule all:
    input:
        fq = os.path.join(OD, '{sample}_trim.fq.gz'),
        html = os.path.join(OD, '{sample}_fastp.html'),
        json = os.path.join(OD, '{sample}_fastp.json'),
        txt = os.path.join(OD, '{sample}_kraken-output.txt'),
        report = os.path.join(OD, '{sample}_kraken-report.txt')

 
# ----------------------------------------------------------
"""
extracting unmapped reads from possorted_genome_bam.bam file
the rule processes all .bam files in input directory INDIR
"""
rule samtools_view:
    input:
        os.path.join(INDIR, '{sample}.bam'),
    output:
        os.path.join(OD, '{sample}_unm.bam')
    threads: 4
    resources:
        mem_mb=10000
    shell:
        """
        samtools view \
            --output-fmt BAM \
            --bam \
            --require-flags 4 \
            --threads {threads} \
            {input} > {output}
        """
# ----------------------------------------------------------
"""
sort extracted reads
"""
rule samtools_sort:
    input:
        os.path.join(OD, '{sample}_unm.bam')
    output:
        os.path.join(OD, '{sample}_unm_srt.bam')
    threads: 4
    resources:
        mem_mb=10000
    shell:
        """
        samtools sort -n \
            --output-fmt BAM \
            --threads {threads} \
            {input} > {output}
        """
        
# ----------------------------------------------------------
"""
converting bam file to R1 and R2 fastq files (in original format)
"""
rule bamtofastq:
    input:
        os.path.join(OD, '{sample}_unm_srt.bam')
    output:
        os.path.join(OD, '{sample}_bamtofastq.done')
    params:
        dir = os.path.join(OD, '{sample}')
    threads: 4
    resources:
        mem_mb=10000
    shell:
        """
        mkdir -p {params.dir}
        bamtofastq \
            --nthreads={threads} \
            {input} \
            {params.dir}/fastq \
        && touch {output}
        """

# ----------------------------------------------------------
"""
concatenating multiple fastq files into one R1
"""
rule concatenate1:
    input:
        os.path.join(OD, '{sample}_bamtofastq.done')
    output:
        os.path.join(OD, '{sample}_R1.fastq.gz')
    params:
        dir = os.path.join(OD, '{sample}', 'fastq')
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        rm -f {output}
        for name in {params.dir}/*/*_R1_*.fastq.gz; do 
           cat $name >> {output}
        done
        """
# ----------------------------------------------------------
"""
concatenating multiple fastq files into one R2
"""
rule concatenate2:
    input:
        os.path.join(OD, '{sample}_bamtofastq.done')
    output:
        os.path.join(OD, '{sample}_R2.fastq.gz')
    params:
        dir = os.path.join(OD, '{sample}', 'fastq')
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        rm -f {output}
        for name in {params.dir}/*/*_R2_*.fastq.gz; do 
           cat $name >> {output}
        done
        """        

# ----------------------------------------------------------
"""
extracting BC & UMI from R1 and placing on R2 read header with umi_tools
"""
rule umi:
    input:
        fq1 = os.path.join(OD, '{sample}_R1.fastq.gz'),
        fq2 = os.path.join(OD, '{sample}_R2.fastq.gz')
    output:
        fq1 = os.path.join(OD, '{sample}_R1_extracted.fastq.gz'),
        fq2 = os.path.join(OD, '{sample}_R2_extracted.fastq.gz')
    params:
        bc = 'CCCCCCCCCCCCCCCCNNNNNNNNNNNN'
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        umi_tools extract \
        --bc-pattern={params.bc} \
        --stdin {input.fq1} \
        --stdout {output.fq1} \
        --read2-in {input.fq2} \
        --read2-out={output.fq2}
        """


# ----------------------------------------------------------
"""
Quality controol on R2 fastq files

cutadapt
 -g trims adapters from 5'end (TSO adapter sequence) and polyA
 -m removing reads shorter than minimum read length of 31
 -n removing up to 2 adapters
"""
rule cutadapt:
    input:
        os.path.join(OD, '{sample}_R2_extracted.fastq.gz')
    output:
        os.path.join(OD, '{sample}_trimTSO.fq.gz')
    params:
        adaptor = 'AAGCAGTGGTATCAACGCAGAGTACATGGG',
        polyA = '^A\{30\}'
    threads: 6
    resources:
        mem_mb=10000
    shell:
        """
        cutadapt \
            --front {params.adaptor} \
            --front {params.polyA} \
            --times 2 \
            --minimum-length 31 \
            --cores {threads} \
            --output {output} \
            {input}
        """


# ----------------------------------------------------------
"""
fastp
 -l minimum read length = 31 (kraken kmer length)
 -A disables automatic adapter trimming
 --trim_poly_X trims polyX sequences from 3' end
 trimmomatic-style sliding window quality trimming
 
NOTE the tools outputs this:
WARNING: you specified the options for cutting by quality, but forogt to enable any of cut_front/cut_tail/cut_right. This will have no effect.
"""
rule fastp:
    input:
        os.path.join(OD, '{sample}_trimTSO.fq.gz')
    output:
        fq = os.path.join(OD, '{sample}_trim.fq.gz'),
        html = os.path.join(OD, '{sample}_fastp.html'),
        json = os.path.join(OD, '{sample}_fastp.json')
    params:
    threads: 6
    resources:
        mem_mb=10000
    shell:
        """
         fastp --in1 {input} \
             --cut_right_window_size 4 \
             --cut_right_mean_quality 20 \
             --disable_adapter_trimming \
             --length_required 31 \
             --thread {threads} \
             --dont_eval_duplication \
             --trim_poly_x \
             --out1 {output.fq} \
             --html {output.html} \
             --json {output.json}
        """

# ----------------------------------------------------------
"""
Kraken2 profiling

Define DB path 
K2DB=/projects/site/pred/microbiome/database/kraken2_Standard_PlusPF/

Step 1: Classify reads with kraken2 (with KrakenUniq option)

singularity exec --home /projects /projects/site/pred/microbiome/microbiome-toolbox/apps/containers/kraken2/latest \
  kraken2 --db $K2DB \
  --memory-mapping \
  --confidence 0.1 \
  --threads 24 \
  --use-names \
  --gzip-compressed \
  --report-minimizer-data \
  --output $OUTDIR/${NAME}/kraken/${NAME}_output.txt \
  --report $OUTDIR/${NAME}/kraken/${NAME}_report.txt \
  $OUTDIR/${NAME}/${NAME}_trim.fq.gz
"""
rule kraken_classify:
    input:
        os.path.join(OD, '{sample}_trim.fq.gz')
    output:
        txt = os.path.join(OD, '{sample}_kraken-output.txt'),
        report = os.path.join(OD, '{sample}_kraken-report.txt')
    params:
        db = config['kraken_db']
    threads: 24
    resources:
        mem_mb=10000
    shell:
        """
          kraken2 --db {params.db} \
              --memory-mapping \
              --confidence 0.1 \
              --threads {threads} \
              --use-names \
              --gzip-compressed \
              --report-minimizer-data \
              --output {output.txt} \
              --report {output.report} \
              {input}
        """


# ----------------------------------------------------------
"""
Kraken2 profiling

Define DB path 
K2DB=/projects/site/pred/microbiome/database/kraken2_Standard_PlusPF/

Step 2: modifying kraken output as input for R package
"""
### extract only classified reads & remve all human reads
#sed -e '/^U/d' -e '/sapiens/d' $OUTDIR/${NAME}/kraken/${NAME}_output.txt > $OUTDIR/${NAME}/${NAME}_filtered_output.txt
### extract spatial BC & UMIs into separate tabs
#sed -i 's/\_/\t/g' $OUTDIR/${NAME}/${NAME}_filtered_output.txt
### extract taxid into separate tab
#sed -i -e 's/ (taxid /\t/g' -e 's/)//g' $OUTDIR/${NAME}/${NAME}_filtered_output.txt
### extract only CB, UMI, taxid columns
#awk -F "\t" '{ print $3,"\t",$4,"\t",$6 }' $OUTDIR/${NAME}/${NAME}_filtered_output.txt > $OUTDIR/${NAME}/${NAME}_kraken_output.txt
