import glob
import datetime

configfile: "config.yaml"

inputdirectory=config["directory"]
ref=config["ref"]
projectname=config["projectname"]
wrappers_version="v1.3.0"
wildcard_constraints:
    date="\d\d\d\d-\d\d-\d\d",
    sample="[a-zA-Z]+.+",



##Get PE data
DATE, SAMPLES, = glob_wildcards(inputdirectory+"/{date}/{sample}_R2_001.fastq.gz", followlinks=True)
print("PE Sample List")
print(SAMPLES)
print(DATE)

PE_dict = {DATE[i]: SAMPLES[i] for i in range(len(DATE))}
print("PE:", str(PE_dict))

## Get SE data
ALL_DATE, ALL_SAMPLES, = glob_wildcards(inputdirectory+"/{date}/{sample}_R1_001.fastq.gz", followlinks=True)
ALL_dict = {ALL_DATE[i]: ALL_SAMPLES[i] for i in range(len(ALL_DATE))}
print("ALL:", str(ALL_dict))

##Create SE lists
SE_SAMPLES = list()
SE_DATE = list()
for date in ALL_dict:
    if date not in PE_dict:
        lsample = ALL_dict[date]
        SE_DATE.append(date)
        SE_SAMPLES.append(lsample)

print("SE Sample List")
print(SE_SAMPLES)
print(SE_DATE)


today = datetime.date.today().strftime("%Y-%m-%d")
print("Today's date: %s" % today)

freyja_updatedir="/data/gpfs/assoc/inbre/projects/subhash_verma/covid_wastewater/freyja_db_update"


def get_map_reads_input(wildcards):
    #print("Wildcards:", wildcards.sample)
    if wildcards.sample in SAMPLES and wildcards.date in DATE:
        return [
            "trimmed/fastp/{date}/{sample}.1.fastq.gz",
            "trimmed/fastp/{date}/{sample}.2.fastq.gz",
        ]
    return ["trimmed/fastp/{date}/{sample}.single.fastq.gz"]


##### target rules #####
rule all:
    input: 
       expand("trimmed/fastp/{date}/{sample}.2.fastq.gz", zip, date=DATE, sample=SAMPLES),
       expand("report/fastp/{date}_{sample}.fastp.json", zip, date=DATE, sample=SAMPLES),
       expand("trimmed/fastp/{date}/{sample}.single.fastq.gz", zip, date=SE_DATE, sample=SE_SAMPLES),
       expand("report/fastp/se/{date}_{sample}.fastp.json", zip, date=SE_DATE, sample=SE_SAMPLES),
       expand("minimap2/{date}/{date}_{sample}_aln.bam", zip, date=ALL_DATE, sample=ALL_SAMPLES),
       expand("freyja/{date}/{date}_{sample}.variants.tsv", zip, date=ALL_DATE, sample=ALL_SAMPLES),
       expand("freyja/demix/{date}_{sample}.output", zip, date=ALL_DATE, sample=ALL_SAMPLES),
       expand("freyja/dates/{date}_{sample}.csv", zip, date=DATE, sample=SAMPLES),
       "freyja/times_metadata.csv",
       #"%s/update_%s.flag" % (freyja_updatedir, today ),
       "freyja/demix/aggregate_%s_%s.tsv" % (projectname, today),
       "freyja/wastewater_timecourse_%s_%s.pdf" % (projectname, today),
       "freyja/wastewater_timecourse_lineages_%s_%s.pdf" % (projectname, today),
       "qc/multiqc_report_%s_%s.html" % (projectname, today),


rule fastp_se:
    input:
        sample=[inputdirectory+"/{date}/{sample}_R1_001.fastq.gz"]
    output:
        trimmed=["trimmed/fastp/{date}/{sample}.single.fastq.gz"],
        html="report/fastp/se/{date}_{sample}.html",
        json="report/fastp/se/{date}_{sample}.fastp.json"
        #html="report/fastp/se/{date}_{sample}.html",
        #json="report/fastp/se/{date}_{sample}.fastp.json"
    log:
        "logs/fastp/{date}_{sample}.log"
    params:
        adapters="",
        extra=""
    threads: 16
    wrapper:
        f"{wrappers_version}/bio/fastp"

rule fastp_pe:
    input:
        sample=[inputdirectory+"/{date}/{sample}_R1_001.fastq.gz", inputdirectory+"/{date}/{sample}_R2_001.fastq.gz"]
    output:
        trimmed=["trimmed/fastp/{date}/{sample}.1.fastq.gz", "trimmed/fastp/{date}/{sample}.2.fastq.gz"],
        html="report/fastp/{date}_{sample}.html",
        json="report/fastp/{date}_{sample}.fastp.json"
    log:
        "logs/fastp/{date}_{sample}.log"
    params:
        adapters="--detect_adapter_for_pe",
        extra=""
    threads: 16
    resources: time_min=480, mem_mb=40000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/fastp"

rule minimap2_bam:
    input:
        target=ref, # can be either genome index or genome fasta
        query=get_map_reads_input,
    output:
        "minimap2/{date}/{date}_{sample}_aln.bam"
    log:
        "logs/minimap2/{date}_{sample}.log"
    params:
        extra="-x sr",           # optional
        sorting="coordinate",           # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra=""                # optional: extra arguments for samtools/picard
    threads: 16
    wrapper:
        f"{wrappers_version}/bio/minimap2/aligner"

rule freyja_update:
    input:
        bam=expand("minimap2/{date}/{date}_{sample}_aln.bam", zip, date=ALL_DATE, sample=ALL_SAMPLES),
    output:
        "%s/update_%s.flag" % (freyja_updatedir, today ),
    log:
        "logs/freyja/update_%s.log" % today
    threads: 1
    resources: time_min=480, mem_mb=40000, cpus=1
    conda:
        "freyja"
    shell:
        "freyja update > {log} 2>&1 && touch {output}"

rule freyja_variant:
    input:
        bam="minimap2/{date}/{date}_{sample}_aln.bam",
        ref=ref, 
        #update="%s/update_%s.flag" % (freyja_updatedir, today ),
    output:
        variants="freyja/{date}/{date}_{sample}.variants.tsv",
        depths="freyja/{date}/{date}_{sample}.depths",
    log:
        "logs/freyja/{date}_{sample}_variant.log"
    params:
        "freyja/{date}"
    threads: 1
    resources: time_min=480, mem_mb=40000, cpus=1
    conda:
        "freyja"
    shell:
        "mkdir -p {params} && freyja variants {input.bam} --variants {output.variants} --depths {output.depths} --ref {input.ref} > {log} 2>&1"

rule freyja_demix:
    input:
        variants="freyja/{date}/{date}_{sample}.variants.tsv",
        depths="freyja/{date}/{date}_{sample}.depths",
    output:
        out="freyja/demix/{date}_{sample}.output",
        #odir=directory("freyja/demix/"),
    log:
        "logs/freyja/{date}_{sample}_demix.log"
    resources: time_min=480, mem_mb=40000, cpus=1
    threads: 1
    conda:
        "freyja"
    shell:
        "freyja demix {input.variants} {input.depths} --output {output.out} > {log} 2>&1"


rule freyja_timescsv_sample:
    input:
        "freyja/{date}/{date}_{sample}.variants.tsv",
    output:
        "freyja/dates/{date}_{sample}.csv",
    params:
        "{date}_{sample}.variants.tsv,{date}",
    threads: 1
    shell:
        """
        echo {params} > {output}
        """

rule freyja_timescsv_combine:
    input:
        csv=expand("freyja/dates/{date}_{sample}.csv", zip, date=ALL_DATE, sample=ALL_SAMPLES)
        #timestemp="freyja/times_metadata_temp.csv"
    output:
        "freyja/times_metadata.csv"
    threads: 1
    shell:
        """
        echo  "Sample,sample_collection_datetime" > {output}
        cat {input.csv} >> {output}
        """

rule freyja_aggregate_all:
    input:
        demix=expand("freyja/demix/{date}_{sample}.output", zip, date=ALL_DATE, sample=ALL_SAMPLES)
    output:
        "freyja/demix/aggregate_%s_%s.tsv" % (projectname, today),
    params:
        indir="freyja/demix/",
    log:
        "logs/freyja/aggregate_all.log"
    threads: 1
    resources: time_min=480, mem_mb=40000, cpus=1
    conda:
        "freyja"
    shell:
        "freyja aggregate {params.indir} --output {output} --ext output > {log} 2>&1"

rule freyja_plot:
    input:
        aggr="freyja/demix/aggregate_%s_%s.tsv" % (projectname, today),
        times="freyja/times_metadata.csv",
    output:
        "freyja/wastewater_timecourse_%s_%s.pdf" % (projectname, today),
    log:
        "logs/freyja/aggregate_plot.log"
    threads: 1
    conda:
        "freyja"
    shell:
        "freyja plot {input.aggr} --output {output} --times {input.times} --interval D --windowsize 14 2> {log}"

rule freyja_plot_lineages:
    input:
        aggr="freyja/demix/aggregate_%s_%s.tsv" % (projectname, today),
        times="freyja/times_metadata.csv",
    output:
        "freyja/wastewater_timecourse_lineages_%s_%s.pdf" % (projectname, today),
    log:
        "logs/freyja/aggregate_plot_lineages.log"
    threads: 1
    conda:
        "freyja"
    shell:
        #"freyja plot {input.aggr} --output {output} --times {input.times} --interval D --windowsize 14 --colors colors.csv --lineages 2> {log}"
        "freyja plot {input.aggr} --output {output} --times {input.times} --interval D --windowsize 14 --lineages 2> {log}"


rule samtools_index:
    input:
        "minimap2/{date}/{date}_{sample}_aln.bam",
    output:
        "minimap2/{date}/{date}_{sample}_aln.bam.bai",
    params:
        "" # optional params string
    resources: time_min=480, mem_mb=2000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/samtools/index"

##QC
rule qualimap:
    input:
        bam="minimap2/{date}/{date}_{sample}_aln.bam",
        bai="minimap2/{date}/{date}_{sample}_aln.bam.bai"
    output:
        outdir=directory("qc/qualimap/{date}_{sample}/"),
        genomestats="qc/qualimap/{date}_{sample}/genome_results.txt"
    log:
        "logs/qualimap/{date}_{sample}.log"
    params:
    resources: time_min=480, mem_mb=30000, cpus=16
    threads: 16  # Use at least two threads
    conda:
        "envs/qualimap.yaml"
    shell:
        "(qualimap bamqc -bam {input.bam} -nt {threads} -outdir {output.outdir}) 2> {log}"

rule kraken2_pe:
    input:
        query=get_map_reads_input,
        db=config["krakendb"]
    output:
        out="qc/kraken2/pe/{date}_{sample}.txt",
        report="qc/kraken2/pe/{date}_{sample}_report.txt",
    log:
        "logs/kraken2/{date}_{sample}.log"
    params:
    resources: time_min=980, mem_mb=90000, cpus=16
    threads: 16  # Use at least two threads
    conda:
        "envs/kraken2.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} --report {output.report} --output {output.out} --paired {input.query} ) 2> {log}"

rule kraken2_se:
    input:
        query=get_map_reads_input,
        db=config["krakendb"]
    output:
        out="qc/kraken2/se/{date}_{sample}.txt",
        report="qc/kraken2/se/{date}_{sample}_report.txt",
    log:
        "logs/kraken2/{date}_{sample}.log"
    params:
    resources: time_min=980, mem_mb=90000, cpus=16
    threads: 16  # Use at least two threads
    conda:
        "envs/kraken2.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} --report {output.report} --output {output.out} {input.query} ) 2> {log}"


rule multiqc:
    input:
        expand("minimap2/{date}/{date}_{sample}_aln.bam", zip, date=ALL_DATE, sample=ALL_SAMPLES),
        expand("qc/kraken2/pe/{date}_{sample}_report.txt", zip, date=DATE, sample=SAMPLES),
        expand("qc/kraken2/se/{date}_{sample}_report.txt", zip, date=SE_DATE, sample=SE_SAMPLES),
        expand("qc/qualimap/{date}_{sample}/", zip, date=ALL_DATE, sample=ALL_SAMPLES),
        expand("report/fastp/{date}_{sample}.fastp.json", zip, date=DATE, sample=SAMPLES),
        expand("report/fastp/se/{date}_{sample}.fastp.json", zip, date=SE_DATE, sample=SE_SAMPLES),
    output:
        "qc/multiqc_report_%s_%s.html" % (projectname, today),
    log:
        "logs/multiqc_fastp.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

