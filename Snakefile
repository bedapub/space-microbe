"""
Snakemake workflow for ST microbiome data analysis pipeline
Spatial Transcriptomics meets Microbiome
"""
import sys
from os.path import join

# ----------------------------------------------------------
# Define the pipeline configuration file
configfile: 'config.yaml'


# ----------------------------------------------------------
# Define variables
OD = config['results']


# ----------------------------------------------------------
# From input FileList create SAMPLES array and temp BAM directory
# and create symbolic links to all BAM files
file_list = config['file_list']
SAMPLES = []
INDIR = os.path.join(OD, 'bam')
BAMDIR = INDIR
if not os.path.exists(BAMDIR):
    os.makedirs(BAMDIR)

with open(config['file_list']) as f:
    for line in f:
       (sample, src) = line.split()
       SAMPLES.append(sample)
       #a[sample] = src
       dest = os.path.join(BAMDIR, sample+'.bam')
       if os.path.exists (src):
           if not os.path.islink(dest):
               os.symlink(src, dest)
       else:
            raise SystemExit('File does not exist:'+src)
#print(SAMPLES)


# ----------------------------------------------------------
# Declare local rules, ie not submitted to the cluster
localrules: all, profiling


# ----------------------------------------------------------
"""
main rule
"""
rule all:
    input:
        os.path.join(OD, 'multiqc_QC', 'multiqc_report.html'),
        expand(os.path.join(OD, '{sample}', '{sample}_profiling-output.txt'), sample=SAMPLES)

        

# ----------------------------------------------------------
"""
extracting unmapped reads from possorted_genome_bam.bam file
the rule processes all .bam files in input directory INDIR
"""
rule samtools_view:
    input:
        os.path.join(INDIR, '{sample}.bam'),
    output:
        temp(os.path.join(OD, '{sample}', '{sample}_unm.bam'))
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
        os.path.join(OD, '{sample}', '{sample}_unm.bam')
    output:
        temp(os.path.join(OD, '{sample}', '{sample}_unm_srt.bam'))
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

Mark directory as TEMP ?
"""
rule bamtofastq:
    input:
        os.path.join(OD, '{sample}', '{sample}_unm_srt.bam')
    output:
        temp(os.path.join(OD, '{sample}', '{sample}_bamtofastq.done')),
        temp(directory(os.path.join(OD, '{sample}', 'fastq')))
    params:
        dir = os.path.join(OD, '{sample}')
    threads: 4
    resources:
        mem_mb=10000
    shell:
        """
        rm -rf {params.dir}/fastq && \
        mkdir -p {params.dir} && \
        bamtofastq \
            --nthreads={threads} \
            {input} \
            {params.dir}/fastq \
        && touch {output[0]}
        """

# ----------------------------------------------------------
"""
concatenating multiple fastq files into one R1
"""
rule concatenate1:
    input:
        os.path.join(OD, '{sample}', '{sample}_bamtofastq.done'),
        os.path.join(OD, '{sample}', 'fastq')
    output:
        temp(os.path.join(OD, '{sample}', '{sample}_R1.fastq.gz'))
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
        os.path.join(OD, '{sample}', '{sample}_bamtofastq.done'),
        os.path.join(OD, '{sample}', 'fastq')
    output:
        temp(os.path.join(OD, '{sample}', '{sample}_R2.fastq.gz'))
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
        fq1 = os.path.join(OD, '{sample}', '{sample}_R1.fastq.gz'),
        fq2 = os.path.join(OD, '{sample}', '{sample}_R2.fastq.gz')
    output:
        fq1 = temp(os.path.join(OD, '{sample}', '{sample}_R1_extracted.fastq.gz')),
        fq2 = temp(os.path.join(OD, '{sample}', '{sample}_R2_extracted.fastq.gz'))
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
Quality control on R2 fastq files

cutadapt
 -g trims adapters from 5'end (TSO adapter sequence) and polyA
 -m removing reads shorter than minimum read length of 31
 -n removing up to 2 adapters
"""
rule cutadapt:
    input:
        os.path.join(OD, '{sample}', '{sample}_R2_extracted.fastq.gz')
    output:
        fq = temp(os.path.join(OD, '{sample}', '{sample}_trimTSO.fq.gz')),
        report = temp(os.path.join(OD, '{sample}', '{sample}_cutadapt.txt'))
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
            --output {output.fq} \
            {input} > {output.report}
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
        os.path.join(OD, '{sample}', '{sample}_trimTSO.fq.gz')
    output:
        html = os.path.join(OD, '{sample}', '{sample}_fastp.html'),
        fq = os.path.join(OD, '{sample}', '{sample}_trim.fq.gz'),
        json = temp(os.path.join(OD, '{sample}', '{sample}_fastp.json'))
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

How much memory is required ?

Step 1: Classify reads with kraken2 (with KrakenUniq option)
"""
rule classify:
    input:
        os.path.join(OD, '{sample}', '{sample}_trim.fq.gz')
    output:
        txt = protected(os.path.join(OD, '{sample}', '{sample}_kraken-output.txt')),
        report = protected(os.path.join(OD, '{sample}', '{sample}_kraken-report.txt'))
    params:
        db = config['kraken_db']
    threads: config['kraken_threads']
    resources:
        mem_mb=100000
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

Step 2: modifying kraken output as input for R package
"""
rule profiling:
    input:
        os.path.join(OD, '{sample}', '{sample}_kraken-output.txt')
    output:
        f = temp(os.path.join(OD, '{sample}', '{sample}_filtered-output.txt')),
        p = os.path.join(OD, '{sample}', '{sample}_profiling-output.txt')
    shell:
        """
        ### extract only classified reads & remove all human reads
        sed -e '/^U/d' -e '/sapiens/d' {input} > {output.f}
        ### extract spatial BC & UMIs into separate tabs
        sed -i 's/\_/\t/g' {output.f}
        ### extract taxid into separate tab
        sed -i -e 's/ (taxid /\t/g' -e 's/)//g' {output.f}
        ### extract only CB, UMI, taxid columns
        awk -F \"\t\" '{{ print $3,\"\t\",$4,\"\t\",$6 }}' {output.f} > {output.p}
        """


# ----------------------------------------------------------
"""
MultiQC for cutadapt and fastp
and remove BAM directory
"""
rule multiqc:
    input:
        expand(os.path.join(OD, '{sample}', '{sample}_fastp.html'), sample=SAMPLES),
        expand(os.path.join(OD, '{sample}', '{sample}_fastp.json'), sample=SAMPLES),
        expand(os.path.join(OD, '{sample}', '{sample}_cutadapt.txt'), sample=SAMPLES),
        expand(os.path.join(OD, '{sample}', '{sample}_kraken-report.txt'), sample=SAMPLES)
    output:
        os.path.join(OD, 'multiqc_QC', 'multiqc_report.html')
    params:
        indir = OD,
        outdir = os.path.join(OD, 'multiqc_QC')
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        multiqc \
            --force \
            --module fastp \
            --module kraken \
            --module cutadapt \
            --outdir {params.outdir} \
            {params.indir} \
        && rm -rf {BAMDIR}
        """
      
        
        