"""
Wrapper script for the "Spatials Transcriptomics meets Microbiome" Snakemake workflow
"""
import os
import sys
import logging
import argparse
import subprocess
from pathlib import Path

SNAKEFILE = 'Snakefile'
GRAPHFILE = 'rulegraph.pdf'
CONFIGFILE = 'config.yaml'

"""
Create the config file
"""
def config(args):
    configfile = os.path.join(args.outdir, CONFIGFILE)
    f = open(configfile, 'w')
    f.write('file_list: \''+args.file_list
            +'\'\nresults: \''+args.outdir
            +'\'\nkraken_db: \''+args.kraken_db
            +'\'\nkraken_threads: '+str(args.kraken_threads)+'\n')
    f.close()
    
    
"""
Create Snakemake Graph file
"""
def graph(args):
    snakefile = os.path.join('.', SNAKEFILE)
    outfile = os.path.join(args.outdir, GRAPHFILE)
    configfile = os.path.join(args.outdir, CONFIGFILE)
                        
    cmd = 'snakemake'\
          +' --rulegraph all'\
          +' --snakefile '+snakefile\
          +' --configfile '+configfile\
          +' | dot -Tpdf > '+outfile
    logging.info('COMMAND='+cmd)
    subprocess.run(cmd, shell=True, check=True)
    
    
"""
Run Snakemake workflow
"""
def workflow(args):
    snakefile = os.path.join('.', SNAKEFILE)
    configfile = os.path.join(args.outdir, CONFIGFILE)
                        
    cmd = 'snakemake'\
          +' --snakefile '+snakefile\
          +' --configfile '+configfile\
          +' --latency-wait 6'\
          +' --rerun-incomplete'\
          +' --keep-going'\
          +' --verbose'
    if args.profile:
        cmd = cmd+' --jobs '+str(args.jobs)+' --profile '+args.profile
        logging.info('Submit workflow to the cluster')
    else:
        cmd = cmd+' --cores '+str(args.cores)
        logging.info('Run workflow locally')
    logging.info('COMMAND='+cmd)
    subprocess.run(cmd, shell=True, check=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Wrapper script for the \"Spatial Transcriptomics meets Microbiome\" Snakemake workflow.'
    )
    parser.add_argument('--file-list', '-f', help='Path to file list (tabular with Sample-ID and BAM-filepath)', required=True)
    parser.add_argument('--outdir', '-o', help='Path to output directory, may not exist', required=True)
    parser.add_argument('--kraken-db', '-d', help='Path to the Kraken database', required=True)
    parser.add_argument('--kraken-threads', '-k', help='Number of cores to use for the Kraken classification step', default=12)
    parser.add_argument('--cores', '-c', help='Number of cores to use for local run', required=False, default=8)
    parser.add_argument('--jobs', '-j', help='Number of jobs for running on the cluster', required=False, default=100)
    parser.add_argument('--profile', '-p', help='Path to cluster profile, if omitted Snakemake will be run locally', 
                        required=False, default=None) #LSF_PROFILE=/apps/rocs/etc/apps/snakemake/lsf/v1.4_memfix
    args = parser.parse_args()

    # Create output directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Create config.yaml file
    config(args)
    
    # Create graph file
    graph(args)
    
    # Launch workflow
    workflow(args)