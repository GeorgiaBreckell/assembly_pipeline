#Basecalling Tims Strains. 
#Fast5 reads will be in the dir data/fast5/{strain} after being depmultiplexed by deepbinner. 
#Fats5s are then basecalled by strain by Guppy, the multple fastq per strain are then merged into a
#single file and this is moved to a dir called ONT and then zipped for easier storage. 

#run_ID=config["run_ID"]

configfile:
    "strains.yml"

rule all:
    input:
        expand("data/fast5/{strain}_basecalled", strain=config["strain"]),
        #expand("{run_ID}/{strain}.fastq.gz", strain=config["strain"], run_ID=config["run_ID"]),
        #expand("/data/{shared_path}/{run_ID}.fastq.tar.gz", shared_path=config["fastq_shared_path"], run_ID=config["run_ID"]),
        #expand("/data/{personal_path}/{strain}.fastq.gz", personal_path=config["fastq_personal_path"], strain=config["strain"]),
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

#rule all_fastq_zipped:
#    input:
#        "{run_ID}/{strain}.fastq.gz"
#    output:
#        temp("all_fastq_zipped.txt"),
#    shell:
#        "if [ ls {run_ID} | wc -l == {config[num_samples]}] then touch {output}"


#rule tar_zip_fastq:
#    input:
#        "all_fastq_zipped.txt"
#    output: 
#        "{run_ID}.fastq.tar.gz"
#    params:
#        "{run_ID}/"
#    shell:
#        "tar -cvzf {output} {params}"

#rule move_shared_fastq: 
#    input:
#        "{run_ID}.fastq.tar.gz"
#    output: 
#        "/data/{shared_path}/{run_ID}.fastq.tar.gz"
#    shell:
#        "cp {input} {output}"

rule filtlong:
    input: 
        R1="/data/Georgia/SC12A_reads/illumina/{strain}_R1.fastq.gz",
        R2="/data/Georgia/SC12A_reads/illumina/{strain}_R2.fastq.gz",
        ONT="data/reads/raw_ONT/{strain}.fastq"
    output:
        "data/reads/filtered_ONT/{strain}.fastq"
    shell:
        "filtlong -1 {input.R1} -2 {input.R2} --min_length 1000 --mean_q_weight 10 --trim --split 500 --target_bases 500000000 {input.ONT} > {output}"

#rule zip_fastq:
#    input:
#        "filtered_ONT/{strain}.fastq"
#    output:
#        zipped="filtered_ONT/{strain}.fastq.gz", 
#    shell:        
#        "gzip {input}"

#rule move_personal_fastq: 
#    input:
#        "filtered_ONT/{strain}.fastq.gz"
#    output: 
#        "/data/{personal_path}/{strain}.fastq.gz"
#    shell:
#        "cp {input} {output}"




