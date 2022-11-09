# snakemake_freyja_covidwastewater

This snakemake pipeline is designed to automate periodic Sars-CoV-2 wastewater sequencing data. It uses PE or SE illumina sequencing data, trims the reads with fastp, maps with minimap2, classifies the reads with Kraken2, then processes the BAM files with Freyja to create an aggregated summary file and figure. 

### Configuration

1. Clone this repo

2. Copy the config.template.yaml to config.yaml.

3. Edit the config.yaml for all parameters
`directory` Input Directory with sequencing data organized in a specific manner. Sub-Directories should be named with the date of wastewater sample collection in YYYY-MM-DD format. Within this sub-directory, only 1 sample may be present with R1 and R2 reads labeled in the filename. If only an R1 filename is found, it will assume SE illumina sequencing.`

`krakendb` The location for the uncompressed kraken database. These can be obtained here: https://benlangmead.github.io/aws-indexes/k2

`ref` Location of the Wuhan genome assembly which can be downloaded from Freyja's project page here: https://raw.githubusercontent.com/andersen-lab/Freyja/main/freyja/data/NC_045512_Hu-1.fasta

`projectname` Project Name will be added to filenames. This allows for multiple sampling sites to be used when generating reports. 


4. Create a "freyja" conda environment

5. Run `freyja update` to pull the latest covid variant classifications.


### Periodic Updates

The pipeline is designed to create a new aggregate and figure every time you run it with a date in the filenames. You may want to run `freyja update` within the freyja conda environment before running the pipeline to update the classifications.


