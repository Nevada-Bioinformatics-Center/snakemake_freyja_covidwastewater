# snakemake_freyja_covidwastewater

With collaboration with Dr. Subhash Verma and Dr. Krishna Pagilla's labs, we have developed the following pipeline to process covid wastewater samples.

This snakemake pipeline is designed to automate periodic Sars-CoV-2 wastewater sequencing data. It uses PE or SE illumina sequencing data, trims the reads with fastp, maps with minimap2, classifies the reads with Kraken2, then processes the BAM files with Freyja to create an aggregated summary file and figure. 

### Configuration

1. Clone this repo

2. Copy the config.template.yaml to config.yaml.

3. Edit the config.yaml for all parameters

`directory` Input Directory with sequencing data organized in a specific manner. Sub-Directories should be named with the date of wastewater sample collection in YYYY-MM-DD format. Within this sub-directory, only 1 sample may be present with R1 and R2 reads labeled in the filename. If only an R1 filename is found, it will assume SE illumina sequencing.

I prefer creating a "raw_seq" folder with the date I downloaded the data. Then symlink the sample data into the appropriate dated folder. The directory tree below may help visualize this structure.

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


`krakendb` The location for the uncompressed kraken database. These can be obtained here: https://benlangmead.github.io/aws-indexes/k2

`ref` Location of the Wuhan genome assembly which can be downloaded from Freyja's project page here: https://raw.githubusercontent.com/andersen-lab/Freyja/main/freyja/data/NC_045512_Hu-1.fasta

`projectname` Project Name will be added to filenames. This allows for multiple sampling sites to be used when generating reports. 


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

The pipeline is designed to create a new aggregate and figure every time you run it with a date in the filenames. You may want to run `freyja update` within the freyja conda environment before running the pipeline to update the classifications.

### Run the pipeline

Run the pipeline using snakemake.  Installation instructions here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html


The following commands assume you are running snakemake on a local computer. 
Be sure to configure the --cores variable with the proper amount on your local system. Snakemake will automatically parallelize the jobs it can run at once.

```
conda activate snakemake
snakemake --use-conda -prn --cores 16  ## This command tests and does a dry-run of the pipeline
snakemake --use-conda -pr --cores 16  ## This command actually runs the pipeline
```

