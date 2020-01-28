##################################################################################
###     Assembly pipeline for hybrid assemblies of ONT and illumina data
###
###     Part 3 of 3 To be used with assembly.Snakefile and 
###     assembly_polishing.Snakefile  
###     
###     
###     This snakemake is the third step in the assembly pipeline, 
###     it takes as input a polished assembly generated from the second snakemake
###     for the strains and assemblers listed in the config file. 
###     This snakemake generates metrics describing the accuracy of each assembly. 
###     The snakemake should be run in the same dir as the 1st and 2nd Snakefiles were run.  
###     Each assembly will be assessed as follows: (results saved in a QC directory) 
###     Nanopore and Illumina reads are mapped to the finished assembly and pdf plots of read coverage
###     are produced. 
###     Illumina withheld reads are mapped to the finished assembly and the fraction of concordantly
###     mapping reads is calculated, and output to a single file per assembler. 
###     A Socru assesment of rRNA operon arrangment is performed to assess global assembly structure.
###     An ORF assesment is is performed, as suggested by Mick Watson.
### 
###         
###     Georgia Breckell 28.01.2020
###
###
###################################################################################

configfile:
    "SC11A_config.yml"

rule all:
    input:
        expand("{strain}/{assembler}/QC/all_nanopore_mapping_coverage.pdf", strain=config["ONT"], assembler=config["assembler"]),
        expand("{strain}/{assembler}/QC/all_illumina_mapping_coverage.pdf", strain=config["ONT"], assembler=config["assembler"]),
        expand("{strain}/{assembler}/QC/concordent_mapping.txt", strain=config["ONT"], assembler=config["assembler"]),
        expand("{strain}/{assembler}/QC/ORF.pdf", strain=config["ONT"], assembler=config["assembler"]),
        expand("{strain}/{assembler}/QC/socru.txt",strain=config["ONT"], assembler=config["assembler"]),



###Coverage and Concordent read assesment
#Index the polished genome for mapping metrics
rule BWA_index: 
    input:
        "{strain}/{assembler}/polished_genome.fasta"
    output:
        touch("{strain}/{assembler}/polished_genome.IDXDONE")
    shell:
        "bwa index {input}" 

#Mapping all illumina reads to the genome
rule BWA_mem_illumina_all: 
    input:
        R1="Illumina/{strain}_R1.fastq.gz", 
        R2="Illumina/{strain}_R2.fastq.gz",
        idxdone="{strain}/{assembler}/polished_genome.IDXDONE"
    output:
        temp("{strain}/{assembler}/QC/all_illumina_mapping.sam")
    params:
        indexed_fasta="{strain}/{assembler}/polished_genome.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.R1} {input.R2} > {output}"

#Mapping the 10,000 withheld Illumina reads to the genome
rule BWA_mem_illumina_WH: 
    input:
        R1="Illumina_subsampled/{strain}_R1_WH.fastq.gz", 
        R2="Illumina_subsampled/{strain}_R2_WH.fastq.gz",
        idxdone="{strain}/{assembler}/polished_genome.IDXDONE"
    output:
        temp("{strain}/{assembler}/QC/WH_illumina_mapping.sam")
    params:
        indexed_fasta="{strain}/{assembler}/polished_genome.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.R1} {input.R2} > {output}"

#Mapping nanopore reads to the genome
rule BWA_mem_nanopore: 
    input:
        reads="ONT/{strain}.fastq.gz",
        idxdone="{strain}/{assembler}/polished_genome.IDXDONE"
    output:
        temp("{strain}/{assembler}/QC/nanopore_mapping.sam")               
    params:              
        indexed_fasta="{strain}/{assembler}/polished_genome.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.reads}  > {output}"

#Sorting the alingment files to allow extraction of coverage and concordant mapping data
rule samtools_sort_all:
    input:
        "{strain}/{assembler}/QC/all_illumina_mapping.sam"
    output:
        "{strain}/{assembler}/QC/all_illumina_mapping_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_sort_WH:
    input:
        "{strain}/{assembler}/QC/WH_illumina_mapping.sam"
    output:
        "{strain}/{assembler}/QC/WH_illumina_mapping_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_sort_nanopore:
    input:
        "{strain}/{assembler}/QC/nanopore_mapping.sam" 
    output:
        "{strain}/{assembler}/QC/nanopore_mapping_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

#Extracting coverage stats for the alignments.
rule samtools_coverage_illumina:
    input:
        "{strain}/{assembler}/QC/all_illumina_mapping_sorted.bam"
    output:
        "{strain}/{assembler}/QC/all_illumina_mapping_coverage.txt"
    shell:
        "samtools depth {input} > {output}"

rule samtools_coverage_nanopore:
    input:
        "{strain}/{assembler}/QC/nanopore_mapping_sorted.bam"
    output:
        "{strain}/{assembler}/QC/all_nanopore_mapping_coverage.txt"
    shell:
        "samtools depth {input} > {output}"

#Calculating concordant reads for Illumina alignments        
rule illumina_WH_mapping_concordent_reads:
    input:  
        "{strain}/{assembler}/QC/WH_illumina_mapping_sorted.bam"
    output:
        "{strain}/{assembler}/QC/illumina_WH_mapping_concordent_reads.txt"
    shell:
        "samtools view -c -f 0x2 {input} > {output}"

rule illumina_all_mapping_concordent_reads:
    input:  
        "{strain}/{assembler}/QC/all_illumina_mapping_sorted.bam"
    output:
        "{strain}/{assembler}/QC/illumina_all_mapping_concordent_reads.txt"
    shell:
        "samtools view -c -f 0x2 {input} > {output}"

#Make one combined file of concordent reads for all strains 
rule Illumina_concoordent_reads_file:
    input:
        "{strain}/{assembler}/QC/illumina_WH_mapping_concordent_reads.txt"
    output:
        touch("{strain}/{assembler}/QC/concordent_mapping.txt")
    params:
        "{assembler}_concordent_mapping.txt"
    run:
        shell("echo {input} >> {params}"),
        shell("cat {input} >> {params}")

#Plotting coverage stats
rule plot_illumina_coverage:
    input:
        "{strain}/{assembler}/QC/all_illumina_mapping_coverage.txt"
    output:
        "{strain}/{assembler}/QC/all_illumina_mapping_coverage.pdf"
    shell:
        "Rscript ./scripts/coverage_plots_g.R {input}"

rule plot_nanopore_coverage:
    input:
        "{strain}/{assembler}/QC/all_nanopore_mapping_coverage.txt"
    output:
        "{strain}/{assembler}/QC/all_nanopore_mapping_coverage.pdf"
    shell:
        "Rscript ./scripts/coverage_plots_g.R {input}"

###Socru assesment
#Run socru assesment on assemblies 
rule socru:
    input:
        "{strain}/{assembler}/polished_genome.fasta"
    output:
        "{strain}/{assembler}/QC/socru.txt"
    shell:
        "socru Escherichia_coli {input} > {output}"

#Make one combined file of Socru output for all strains
rule Scoru_file:
    input:
        "{strain}/{assembler}/QC/socru.txt"
    output:
        touch("{strain}/{assembler}/QC/socru.txt")
    params:
        "{assembler}_socru_mapping.txt"
    run:
        shell("echo {input} >> {params}"),
        shell("cat {input} >> {params}")

###ORF assesment
#Run prodigal step of ORF assesment
rule ORF_prodigal:
        input: 
            "{strain}/{assembler}/polished_genome.fasta"
        output: 
            "{strain}/{assembler}/QC/ORF.faa"
        benchmark: 
            "{strain}/{assembler}/benchmarks/{strain}.prodigal.bench.txt"
        shell: 
            "prodigal -a {output} -q -i {input}"

#Run Diamond step of ORF assesment
rule ORF_diamond:
        input: 
            "{strain}/{assembler}/QC/ORF.faa"
        output: 
            "{strain}/{assembler}/QC/ORF.data"
        threads: 
            32
        params:
                db="/data/databases/diamond/uniprot.dmnd",
                of="6 qlen slen"
        benchmark: 
            "{strain}/{assembler}/benchmarks/{strain}.diamond.bench.txt"
        shell: 
            "diamond blastp --threads 32 --max-target-seqs 1 --db {params.db} --query {input} --outfmt {params.of} --out {output}"

#Run plotting step of ORF assesment
rule ORF_hist:
        input: 
            "{strain}/{assembler}/QC/ORF.data"
        output: 
            "{strain}/{assembler}/QC/ORF.pdf"
        shell: 
            "R --slave --no-restore --file=hist_orf_lengths.R --args {input} {output}"







       
