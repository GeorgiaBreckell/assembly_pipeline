##################################################################################
###    Assembly polishing pipeline for hybrid assemblies of ONT and illumina data
###
###     Part 2 of 3 To be used with assembly.Snakefile and 
###     assesment.Snakefile
###                     
###		This snakemake is the second step in the assembly pipeline, 
###		it takes as input an assembly generated from the first snakemake
###		for the strains and assemblers listed in the config file.
###     The snakemake should be run in the same dir as the assembly.Snakefile was run.  
###     Each assembly will be polished 4 times iteratively with Racon using ONT reads, 
###		followed by illumina read polishing using Pilon. Finally the assembly is polished 
###     using Racon and adapted Illumina reads.
###     The final output is a polished genome fasta file and a genome stats file for each 
###		assembly, as well as genome stats pooled by assembler. 
### 
###         
###     Georgia Breckell 19.11.2019
###
###
###
###################################################################################

configfile:
    "SC11A_config.yml"

rule all:
    input: 
        expand("{strain}/{assembler}_genome_stats.txt", strain=config["ONT"],assembler=config["assembler"])
        
###GENOME POLISHING

#First alignment for the first round of polishing. Assembly reads are aligned to the assembly using Minimap2
#This needs to be a seperate input for each assembler as they all have different 
#output formats, following this first round of polishing the different assemblers can be treated the same 

rule Minimap2_round1:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="{strain}/{assembler}/assembly.fasta",
    output:
        alignment=temp("{strain}/{assembler}/polishing/Minimap2_round1.sam"),
    conda:
        "python3_6.yml"
    shell:
        "minimap2 -a {input.assembly} {input.reads} > {output.alignment}"

#First round of Racon polishing using Minimap2 alignments, assembly reads and the assembly
rule Racon_round1:
    input:
        reads="ONT/{strain}.fastq.gz", 
        alignment="{strain}/{assembler}/polishing/Minimap2_round1.sam",
        assembly="{strain}/{assembler}/assembly.fasta",
    output:
        Racon_round1="{strain}/{assembler}/polishing/Racon_round1.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round1}"),

#Second alignment for polishing, the output from the first round of polishing is used as the input assembly. 
rule Minimap2_round2:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="{strain}/{assembler}/polishing/Racon_round1.fasta",

    output:
        alignment=temp("{strain}/{assembler}/polishing/Minimap2_round2.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Second round of Racon polishing. The output from the first round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to first racon output.
rule Racon_round2:
    input:
        reads="ONT/{strain}.fastq.gz",
        alignment="{strain}/{assembler}/polishing/Minimap2_round2.sam",
        assembly="{strain}/{assembler}/polishing/Racon_round1.fasta",
    output:
        Racon_round2="{strain}/{assembler}/polishing/Racon_round2.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round2}")

#Third alignment for polishing, the output from the second round of polishing is used as the input assembly.
rule Minimap2_round3:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="{strain}/{assembler}/polishing/Racon_round2.fasta",
    output:
        alignment=temp("{strain}/{assembler}/polishing/Minimap2_round3.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Third round of Racon polishing. The output from the second round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to second racon output.
rule Racon_round3:
    input:
        reads="ONT/{strain}.fastq.gz",
        alignment="{strain}/{assembler}/polishing/Minimap2_round3.sam",
        assembly="{strain}/{assembler}/polishing/Racon_round2.fasta",
    output:
        Racon_round3="{strain}/{assembler}/polishing/Racon_round3.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round3}"),

#Fouth alignment for polishing, the output from the third round of polishing is used as the input assembly.
rule Minimap2_round4:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="{strain}/{assembler}/polishing/Racon_round3.fasta",
    output:
        alignment=temp("{strain}/{assembler}/polishing/Minimap2_round4.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Fouth round of Racon polishing. The output from the third round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to third racon output.
rule Racon_round4:
    input:
        reads="ONT/{strain}.fastq.gz",
        alignment="{strain}/{assembler}/polishing/Minimap2_round4.sam",
        assembly="{strain}/{assembler}/polishing/Racon_round3.fasta",
    output:
        Racon_round4="{strain}/{assembler}/polishing/Racon_round4.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round4}"),

#BWA indexing of fouth Racon output.  
rule BWA_index_Racon_polished: 
    input:
        assembly="{strain}/{assembler}/polishing/Racon_round4.fasta",
    output:
        touch("{strain}/{assembler}/polishing/racon_round4_makeidx.done"),
    run:
        shell("bwa index {input.assembly}"),

#BWA mem aligns Illumina reads to the fouth racon output (partially polished genome)
rule BWA_mem_illumina: 
    input:
        R1="Illumina/{strain}_R1.fastq.gz", 
        R2="Illumina/{strain}_R2.fastq.gz", 
        idxdone="{strain}/{assembler}/polishing/racon_round4_makeidx.done",
    output:
        mapping=temp("{strain}/{assembler}/polishing/Illumina_mapping.sam"),
    params:
        indexed_fasta="{strain}/{assembler}/polishing/Racon_round4.fasta",
    run:
        shell("bwa mem {params.indexed_fasta} {input.R1} {input.R2} > {output.mapping}"),

#BWA mem alignment is sorted for input into Pilon
rule sorted_BWA_mem_alignment: 
    input:
        mapping="{strain}/{assembler}/polishing/Illumina_mapping.sam",   
    output:
        "{strain}/{assembler}/polishing/Illumina_mapping_sorted.bam",   
    run:
        shell("samtools sort {input.mapping} -o {output}"),

#BWA mem sorted alignment is indexed 
rule indexed_BWA_mem_alignment:
    input:
        "{strain}/{assembler}/polishing/Illumina_mapping_sorted.bam",
    output:
        touch("{strain}/{assembler}/polishing/BWA_alignment_makeidx.done"),
    run:
        shell("samtools index {input}"),

#Pilon polish the genome using Illumina short reads
rule Pilon_polish: 
    input:
        assembly="{strain}/{assembler}/polishing/Racon_round4.fasta",
        bam="{strain}/{assembler}/polishing/Illumina_mapping_sorted.bam",
        idxdone="{strain}/{assembler}/polishing/BWA_alignment_makeidx.done",
    output:
        pilon="{strain}/{assembler}/polishing/pilon/pilon.fasta",
    params:
        outdir="{strain}/{assembler}/polishing/pilon/",
    run:
        shell("pilon --genome {input.assembly} --bam {input.bam} --outdir {params.outdir}"),

#Index the pilon polished output
rule BWA_index_pilon: 
    input:
        "{strain}/{assembler}/polishing/pilon/pilon.fasta",
    output:
        touch("{strain}/{assembler}/polishing/pilon/makeidx.done"),
    run:
        shell("bwa index {input}"),

#BWA mem is used to align the Illumina reads to the indexed Pilon output
#For this step the Illumina reads need to be combined, this is NOT included in this snakefile
rule Illumina_on_pilon_alignment: 
    input:
        combined_illumina="Illumina/{strain}_illumina.fastq.gz",
        idxdone_pilon="{strain}/{assembler}/polishing/pilon/makeidx.done",
    output:
        temp("{strain}/{assembler}/polishing/pilon/all_illumina_mapping_{strain}.sam")
    params:
        indexed_fasta="{strain}/{assembler}/polishing/pilon/pilon.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.combined_illumina} > {output}"

#Polish the Pilon output with Racon using Illumina reads. Generates the Final polished genome, output to the main directory
rule Racon_illumina:
    input:
        combined_illumina="Illumina/{strain}_illumina.fastq.gz",
        alignment="{strain}/{assembler}/polishing/pilon/all_illumina_mapping_{strain}.sam",
        assembly="{strain}/{assembler}/polishing/pilon/pilon.fasta"
    output:
        Racon_illumina="{strain}/{assembler}/polished_genome.fasta"
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.combined_illumina} {input.alignment} {input.assembly} > {output.Racon_illumina}")

#Making Final Genome stats file
rule genome_stats: 
        input:
            "{strain}/{assembler}/polished_genome.fasta"
        output:
            "{strain}/{assembler}/genome_stats.txt"
        shell:
            "seqkit stats {input} {output}" 

rule pooled_stats:
        input:
            "{strain}/{assembler}/genome_stats.txt"
        output:
            "{strain}/{assembler}_genome_stats.txt"
        params:
            "{assembler}_genome_stats.txt"
        shell:
            "cat {input} >> {params}"
