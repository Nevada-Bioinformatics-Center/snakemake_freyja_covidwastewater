# snakemake_freyja_covidwastewater

This snakemake pipeline is designed to automate periodic Sars-CoV-2 wastewater sequencing data. It uses PE or SE Illumina sequencing data, trims the reads with fastp, maps with minimap2, classifies the reads with Kraken2, then processes the BAM files with Freyja to create an aggregated summary file and figure. 

This pipeline was developed for Drs. Subhash Verma and Krishna Pagilla at the University of Nevada, Reno. 

### Configuration

1. Clone repo

2. Copy config.template.yaml to config.yaml

3. Edit config.yaml for all parameters

`directory` Input directory with sequencing data is organized in a specific manner. Sub-directories should be named with the wastewater sample collection date in YYYY-MM-DD format. Within this sub-directory, only one collected sample may be present with R1 and R2 reads labeled in the filename. If only an R1 filename is found, it will assume SE Illumina sequencing.

Examole of directory tree:
```
INPUT_DIR
├── 2022-10-03
│   ├── WW106_S13_L001_R1_001.fastq.gz -> /data/gpfs/assoc/raw_seq_10202022/WW106_L1_ds.06365230dcfa473a8fb6b8098ab760c4/WW106_S13_L001_R1_001.fastq.gz
│   └── WW106_S13_L001_R2_001.fastq.gz -> /data/gpfs/assoc/raw_seq_10202022/WW106_L1_ds.06365230dcfa473a8fb6b8098ab760c4/WW106_S13_L001_R2_001.fastq.gz
├── 2022-10-10
│   ├── WW113_S17_L001_R1_001.fastq.gz -> /data/gpfs/assoc/raw_seq_11072022/WW113_L1_ds.f0035b32a09a4ba5b32d9f51df134902/WW113_S17_L001_R1_001.fastq.gz
│   └── WW113_S17_L001_R2_001.fastq.gz -> /data/gpfs/assoc/raw_seq_11072022/WW113_L1_ds.f0035b32a09a4ba5b32d9f51df134902/WW113_S17_L001_R2_001.fastq.gz
└── 2022-10-24
    ├── WW117_S21_L001_R1_001.fastq.gz -> /data/gpfs/assoc/raw_seq_11072022/WW117_L1_ds.e1dbdd81b7244f7a9ab5db98e5ddcfdb/WW117_S21_L001_R1_001.fastq.gz
    └── WW117_S21_L001_R2_001.fastq.gz -> /data/gpfs/assoc/raw_seq_11072022/WW117_L1_ds.e1dbdd81b7244f7a9ab5db98e5ddcfdb/WW117_S21_L001_R2_001.fastq.gz
```


`krakendb` Location for the uncompressed kraken database; can be obtained here: https://benlangmead.github.io/aws-indexes/k2

`ref` Location of the Wuhan genome assembly; can be downloaded from Freyja's project page here: https://raw.githubusercontent.com/andersen-lab/Freyja/main/freyja/data/NC_045512_Hu-1.fasta

`projectname` The project name will be added to filenames allowing for multiple sampling sites in an individual report. 


4. Create a "freyja" conda environment

```
##Create environment name
conda create -n freyja
conda activate freyja

##Add the necessary channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

##Install freyja
conda install freyja
```

5. Run `freyja update` to pull the latest covid variant classifications.


### Periodic Updates

The pipeline is designed to create a new aggregated object and figure every time you run it with a date in the filenames. You may want to run `freyja update` within the freyja conda environment before running the pipeline to update the classifications.

### Run the pipeline

Run the pipeline using snakemake.  Installation instructions here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

The following commands assume the snakemake pipeline is ran on a local computer, thus be sure to configure the `--cores` variable with the proper amount available on the local system. Snakemake will automatically parallelize jobs that can be ran at once.

```
conda activate snakemake
snakemake --use-conda -prn --cores 16  ## This command tests and does a dry-run of the pipeline
snakemake --use-conda -pr --cores 16   ## This command actually runs the pipeline
```


### Acknowledgement

This work was supported by funds from the US Treasury through the Coronavirus Aid, Relief, and Economic Security (CARES) Act and grants from the National Institute of General Medical Sciences (GM103440 and GM104944).
