#Basecalling 
#Prior to runing this snakemake, the fast5 reads were demultiplexed at the fast5 level by deepbinner version 0.2.0. 
#With the following command: 
# "deepbinner realtime --in_dir {fast5} --out_dir {demultiplexed_fast5} --rapid"

#Fats5s are then basecalled by strain by Guppy, the multple fastq per strain are then merged into a
#single file and this is moved to a dir called ONT and then zipped before reads are filtered with filtlong version 0.2.0


configfile:
    "strains.yml"

rule all:
    input:
        expand("data/fast5/{strain}_basecalled", strain=config["strain"]),
        expand("data/reads/filtered_ONT/{strain}_zipped", strain=config["strain"])

rule guppy_basecaller:
    input:
        "data/fast5/{strain}_fast5_reads"
    output:
        touch("data/fast5/{strain}_basecalled")
    params:
        "data/reads/{strain}"
    shell:
         "guppy_basecaller -i {input}  -s {params}  --flowcell FLO-MIN106 --kit SQK-RBK004 -x \"cuda:0\" --fast5_out"
     
rule merge_fastq:
    input:
        marker="data/fast5/{strain}_basecalled"
    output:
        "data/reads/raw_ONT/{strain}.fastq"
    params:
        "data/reads/{strain}"
    shell:
        "cat  {params}/*.fastq  >  {output}" 

rule zip_filtered_fastq:
    input:
        "data/reads/filtered_ONT/{strain}.fastq"
    output:
        touch("data/reads/filtered_ONT/{strain}_zipped") 
    run:        
        shell("gzip {input}")


rule filtlong:
    input: 
        R1="/data/Georgia/SC12A_reads/illumina/{strain}_R1.fastq.gz",
        R2="/data/Georgia/SC12A_reads/illumina/{strain}_R2.fastq.gz",
        ONT="data/reads/raw_ONT/{strain}.fastq"
    output:
        "data/reads/filtered_ONT/{strain}.fastq"
    shell:
        "filtlong -1 {input.R1} -2 {input.R2} --min_length 1000 --mean_q_weight 10 --trim --split 500 --target_bases 500000000 {input.ONT} > {output}"





