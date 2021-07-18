##################################################################################
###     Assembly pipeline for hybrid assemblies of ONT and illumina data
###
###     Reads should be in a directoy named ONT and Illumina respectively
###     and output will be in an individual directory for each strain, 
###     organsied by assembler. A table of results including genome stats and QC measures is
###	the final output. 
###
###     Assemblers and Strains can be specifed in the Config.yml file 
###         
###     GB 18.07.2021
###
###
###
###################################################################################

configfile:
    "configs/assemblies.yml"

rule all:
    input:
        expand("results/{strain}/{assembler}/assembly.fasta", strain=config["strain"], assembler=config["assembler"]),
        expand("results/{strain}/{assembler}/output_dir_removed", strain=config["strain"],assembler=config["assembler"]),
        expand("results/{strain}/{assembler}/bandage_plot.png", strain=config["strain"],assembler=config["assembler_gfa"]),
        expand("results/{strain}/{assembler}/genome_stats.pooled", strain=config["strain"],assembler=config["assembler"]),
        expand("results/{strain}/{assembler}/QC/all_nanopore_mapping_coverage.pdf", strain=config["strain"], assembler=config["assembler"]),
        expand("results/{strain}/{assembler}/QC/all_illumina_mapping_coverage.pdf", strain=config["strain"], assembler=config["assembler"]),
        expand("results/{strain}/{assembler}/QC/ORF.pdf", strain=config["strain"], assembler=config["assembler"]),
        expand("results/final_table/{strain}_{assembler}_all_metrics_merged.marker", strain=config["strain"],assembler=config["assembler"])

#Randomly subsample 1000 Nanopore long reads. Generates a new fastq containing the subsampled reads.
rule subsample_ONT:
    input:
        original="ONT/{strain}.fastq.gz"
    output:
        "ONT_subsampled/{strain}_withheld.fastq.gz"
    shell:
        "seqtk sample -s 11 {input} 1000 > {output}"

#Identify which long reads were subsampled by comparing with main dataset. Generates a list of read IDs for the subsampled reads.
rule ID_subsampled_reads_ONT: 
    input:
        subsampled="ONT_subsampled/{strain}_withheld.fastq.gz",
        original="ONT/{strain}.fastq.gz"
    output:
        "ONT_subsampled/{strain}_common.list"
    shell:
        "seqkit common {input.subsampled} {input.original} | grep '@' | cut -c 2-37 > {output} "
        
#Removes subsampled reads from the main dataset, according to the read IDs on the list above.
rule remove_subsampled_reads_ONT:
    input: 
        common_list="ONT_subsampled/{strain}_common.list",
        original="ONT/{strain}.fastq.gz"
    output:
        "ONT_subsampled/{strain}_assembly.fastq.gz"
    shell:
        "seqkit grep -f {input.common_list} -v {input.original} -o {output}"

#Randomly subsample 5000 paired Illumina reads, Generates a new fastq containing the subsampled reads
rule Subsample_Illumina:
    input:
        R1="illumina/{strain}_R1.fastq.gz",
        R2="illumina/{strain}_R2.fastq.gz"
    output:
        subsampled_R1="illumina_subsampled/{strain}_R1_WH.fastq.gz",
        subsampled_R2="illumina_subsampled/{strain}_R2_WH.fastq.gz"
    run:
        shell("seqtk sample -s 11 {input.R1} 5000 > {output.subsampled_R1}")
        shell("seqtk sample -s 11 {input.R2} 5000 > {output.subsampled_R2}")

#Identify which long reads were subsampled by comparing with main dataset. Generates a list of read IDs for the subsampled reads.
rule ID_subsampled_reads_Illumina: 
    input:
        subsampled_R1="illumina_subsampled/{strain}_R1_WH.fastq.gz",
        R1="illumina/{strain}_R1.fastq.gz",
        subsampled_R2="illumina_subsampled/{strain}_R2_WH.fastq.gz",
        R2="illumina/{strain}_R2.fastq.gz"
    output:
        R1="illumina_subsampled/{strain}_R1_common.list",
        R2="illumina_subsampled/{strain}_R2_common.list"
    run:
        shell("seqkit common {input.subsampled_R1} {input.R1} | grep '@' | cut -c 2-37 > {output.R1}") 
        shell("seqkit common {input.subsampled_R2} {input.R2} | grep '@' | cut -c 2-37 > {output.R2}") 
 
##Removes subsampled reads from the main dataset, according to the read IDs on the list above.
#Generates an "Assembly" dataset from which all assemblies are created.    
rule remove_subsampled_reads_Illumina:
    input: 
        R1_common="illumina_subsampled/{strain}_R1_common.list",
        R1="illumina/{strain}_R1.fastq.gz",
        R2_common="illumina_subsampled/{strain}_R2_common.list",
        R2="illumina/{strain}_R2.fastq.gz"
    output:
        R1="illumina_subsampled/{strain}_R1_assembly.fastq.gz",
        R2="illumina_subsampled/{strain}_R2_assembly.fastq.gz"
    run:
        shell("seqkit grep -f {input.R1_common} -v {input.R1} -o {output.R1}")
        shell("seqkit grep -f {input.R2_common} -v {input.R2} -o {output.R2}")

#Unicycler Hybrid assembly using the Assembly dataset of reads, (subsampled reads have been removed). 
rule unicycler_assembly:
    input:
        ont="ONT_subsampled/{strain}_assembly.fastq.gz",
        r1="illumina_subsampled/{strain}_R1_assembly.fastq.gz",
        r2="illumina_subsampled/{strain}_R2_assembly.fastq.gz"
    output:
        fasta="results/{strain}/unicycler/unicycler_output/assembly.fasta",
        gfa="results/{strain}/unicycler/unicycler_output/assembly.gfa"
    params:
        out_prefix="results/{strain}/unicycler/unicycler_output/"
    log:
        "results/{strain}/logs/unicycler.log"
    benchmark:
        "results/{strain}/benchmarks/unicycler.assembly.benchmark.txt"
    threads: 8
    shell:
        "unicycler -1 {input.r1} -2 {input.r2} -l {input.ont} -o {params.out_prefix}" 
        " 1> {log}"

#Copies the unicycler assembly out of the unicycler output dir
rule copy_unicycler:
    input:
        assembly="results/{strain}/unicycler/unicycler_output/assembly.fasta",
        gfa="results/{strain}/unicycler/unicycler_output/assembly.gfa"
    output:
        assembly="results/{strain}/unicycler/assembly.fasta",
        gfa="results/{strain}/unicycler/assembly.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"), 
        shell("cp {input.gfa} {output.gfa}")   

#remove unicycler assembly directory   
rule remove_unicycler_output_dir:
    input:
        assembly="results/{strain}/unicycler/assembly.fasta",
        gfa="results/{strain}/unicycler/assembly.gfa"
    output:
        touch("results/{strain}/unicycler/output_dir_removed")
    params:
        output_dir="results/{strain}/unicycler/unicycler_output/"
    shell:
        "rm -r {params.output_dir}"

#Canu long read only assembly with the "Assembly dataset of reads"
rule canu_assembly:
    input:
        "ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        fasta="results/{strain}/canu/canu_output/{strain}.contigs.fasta",
        gfa="results/{strain}/canu/canu_output/{strain}.unitigs.gfa"
    log:
        "results/{strain}/logs/canu.log"
    benchmark:
        "results/{strain}/benchmarks/canu.assembly.benchmark.txt"
    params:
        directory="results/{strain}/canu/canu_output/",
        prefix="{strain}"
    shell:
        "canu -d {params.directory} -p {params.prefix} genomeSize=5m -nanopore-raw {input} -maxMemory=32g -maxThreads=10" 
        "1>{log}" 

#Copies and renames the canu output fasta to generic naming convention 
rule copy_canu:
    input:
        assembly="results/{strain}/canu/canu_output/{strain}.contigs.fasta",
        gfa="results/{strain}/canu/canu_output/{strain}.unitigs.gfa"
    output:
        assembly="results/{strain}/canu/{strain}.contigs.fasta",
        gfa="results/{strain}/canu/{strain}.unitigs.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"),
        shell("cp {input.gfa} {output.gfa}")

rule rename_canu:
    input:
        assembly="results/{strain}/canu/{strain}.contigs.fasta",
        gfa="results/{strain}/canu/{strain}.unitigs.gfa"
    output:
        assembly="results/{strain}/canu/assembly.fasta",
        gfa="results/{strain}/canu/assembly.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"),
        shell("cp {input.gfa} {output.gfa}")

#remove canu assembly directory   
rule remove_canu_output_dir:
    input:
        assembly="results/{strain}/canu/assembly.fasta",
        gfa="results/{strain}/canu/assembly.gfa"
    output:
        touch("results/{strain}/canu/output_dir_removed")
    params:
        output_dir="results/{strain}/canu/canu_output/"
    shell:
        "rm -r {params.output_dir}"

#Raven long read only assembly with the "Assembly dataset of reads", no renmaing needed as output is in desired convention. 
#Raven output as single fasta file so no movement out of the "output_dir" is needed either 
rule raven_assembly:
    input:
        "ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        assembly="results/{strain}/raven/assembly.fasta",
        dir_removed=touch("results/{strain}/raven/output_dir_removed")
    log:
        "results/{strain}/logs/raven.log"    
    benchmark:
        "results/{strain}/benchmarks/raven.assembly.benchmark.txt"
    shell:
        "raven {input} > {output.assembly}"
        "1> {log} "

#flye long read only assembly with the "Assembly dataset of reads", no renaming is needed because flye output is in desired convention 
rule flye_assembly:
    input:
        "ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        fasta="results/{strain}/flye/flye_output/assembly.fasta",
        gfa="results/{strain}/flye/flye_output/assembly_graph.gfa"
    params:
        out_prefix="results/{strain}/flye/flye_output/"
    log:
        "results/{strain}/logs/flye.log"
    benchmark:
        "results/{strain}/benchmarks/flye.assembly.benchmark.txt"
    conda:
        "environments/assemblies_2_7.yml"
    shell:
        "flye --nano-raw {input} --out-dir {params.out_prefix} --plasmids" 
        "1>{log}"

#Copies the flye assembly out of the flye output dir
rule copy_flye:
    input:
        assembly="results/{strain}/flye/flye_output/assembly.fasta",
        gfa="results/{strain}/flye/flye_output/assembly_graph.gfa"
    output:
        assembly="results/{strain}/flye/assembly.fasta",
        gfa="results/{strain}/flye/assembly_graph.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"),
        shell("cp {input.gfa} {output.gfa}")

#rename GFA file to standard convention
rule rename_flye:
    input:
        gfa="results/{strain}/flye/assembly_graph.gfa"
    output:
        "results/{strain}/flye/assembly.gfa"
    run:
        shell("mv {input.gfa} {output}")

#rule to remove flye assembly directory   
rule remove_flye_output_dir:
    input:
        assembly="results/{strain}/flye/assembly.fasta",
    output:
        touch("results/{strain}/flye/output_dir_removed")
    params:
        output_dir="results/{strain}/flye/flye_output/"
    shell:
        "rm -r {params.output_dir}"

#Wtdbg2(redbean) long read only assembly with assembly dataset
rule wtdbg2_assembly_1:
        input:
            ont="ONT_subsampled/{strain}_assembly.fastq.gz"
        output:
            "results/{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.lay.gz"
        params:
            out_prefix="results/{strain}/wtdbg2/wtdbg2_output/{strain}"
        log:
            "results/{strain}/logs/wtdgb2_1.log"
        benchmark:
            "results/{strain}/benchmarks/wtdbg2_1.assembly.benchmark.txt"
        shell:
            "wtdbg2 -x ont -g 4.8m -t 16 -i {input.ont} -o {params.out_prefix}" 
            "1>{log}"

#Second half of the redbean assembly
rule wtdbg2_assembly_2:
        input:
            wtdbg2_config="results/{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.lay.gz"
        output:
            "results/{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.fa" 
        params:
            out_prefix="results/{strain}/wtdbg2/wtdbg2_output/{strain}"
        log:
            "results/{strain}/logs/wtdgb2_2.log"
        benchmark:
            "results/{strain}/benchmarks/wtdbg2_2.assembly.benchmark.txt"
        shell:
            "wtpoa-cns -t 16 -i {input} -o {params.out_prefix}.ctg.fa" 
            "1>{log}"

#Copies and renames the wtdbg2 output fasta to generic naming convention 
rule copy_wtdbg2_1:
    input:
        "results/{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.fa"
    output:
        "results/{strain}/wtdbg2/{strain}.ctg.fa" 
    shell:
        "cp {input} {output}"

rule rename_wtdbg2_1:
    input:
        "results/{strain}/wtdbg2/{strain}.ctg.fa"
    output:
        "results/{strain}/wtdbg2/assembly.fasta"
    shell:
        "mv {input} {output}"

#rule to remove wtdbg2 assembly directory   
rule remove_wtdbg2_output_dir:
    input:
        assembly="results/{strain}/wtdbg2/assembly.fasta",
    output:
        touch("results/{strain}/wtdbg2/output_dir_removed")
    params:
        output_dir="results/{strain}/wtdbg2/wtdbg2_output/"
    shell:
        "rm -r {params.output_dir}"

#Bandage plots of the assembly graph for each assembler (except Ra and wtdbg2 - no .gfa file produced) 
rule Bandage_plot:
    input:
        gfa="results/{strain}/{assembler}/assembly.gfa"
    output:
        bandage_plot="results/{strain}/{assembler}/bandage_plot.png"
    run:
        shell("Bandage image {input.gfa} {output.bandage_plot}")

###GENOME POLISHING

#First alignment for the first round of polishing. Assembly reads are aligned to the assembly using Minimap2
#This needs to be a seperate input for each assembler as they all have different 
#output formats, following this first round of polishing the different assemblers can be treated the same 

rule Minimap2_round1:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/assembly.fasta",
    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round1.sam"),
    conda:
        "environments/base.yml"
    shell:
        "minimap2 -a {input.assembly} {input.reads} > {output.alignment}"

#First round of Racon polishing using Minimap2 alignments, assembly reads and the assembly
rule Racon_round1:
    input:
        reads="ONT/{strain}.fastq.gz", 
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round1.sam",
        assembly="results/{strain}/{assembler}/assembly.fasta",
    output:
        Racon_round1="results/{strain}/{assembler}/polishing/Racon_round1.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round1}"),

#Second alignment for polishing, the output from the first round of polishing is used as the input assembly. 
rule Minimap2_round2:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/polishing/Racon_round1.fasta",

    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round2.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Second round of Racon polishing. The output from the first round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to first racon output.
rule Racon_round2:
    input:
        reads="ONT/{strain}.fastq.gz",
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round2.sam",
        assembly="results/{strain}/{assembler}/polishing/Racon_round1.fasta",
    output:
        Racon_round2="results/{strain}/{assembler}/polishing/Racon_round2.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round2}")

#Third alignment for polishing, the output from the second round of polishing is used as the input assembly.
rule Minimap2_round3:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/polishing/Racon_round2.fasta",
    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round3.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Third round of Racon polishing. The output from the second round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to second racon output.
rule Racon_round3:
    input:
        reads="ONT/{strain}.fastq.gz",
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round3.sam",
        assembly="results/{strain}/{assembler}/polishing/Racon_round2.fasta",
    output:
        Racon_round3="results/{strain}/{assembler}/polishing/Racon_round3.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round3}"),

#Fouth alignment for polishing, the output from the third round of polishing is used as the input assembly.
rule Minimap2_round4:
    input:
        reads="ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/polishing/Racon_round3.fasta",
    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round4.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Fouth round of Racon polishing. The output from the third round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to third racon output.
rule Racon_round4:
    input:
        reads="ONT/{strain}.fastq.gz",
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round4.sam",
        assembly="results/{strain}/{assembler}/polishing/Racon_round3.fasta",
    output:
        Racon_round4="results/{strain}/{assembler}/polishing/Racon_round4.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round4}"),

#BWA indexing of fouth Racon output.  
rule BWA_index_Racon_polished: 
    input:
        assembly="results/{strain}/{assembler}/polishing/Racon_round4.fasta",
    output:
        touch("results/{strain}/{assembler}/polishing/racon_round4_makeidx.done"),
    run:
        shell("bwa index {input.assembly}"),

#BWA mem aligns Illumina reads to the fouth racon output (partially polished genome)
rule BWA_mem_illumina: 
    input:
        R1="illumina/{strain}_R1.fastq.gz", 
        R2="illumina/{strain}_R2.fastq.gz", 
        idxdone="results/{strain}/{assembler}/polishing/racon_round4_makeidx.done",
    output:
        mapping=temp("results/{strain}/{assembler}/polishing/Illumina_mapping.sam"),
    params:
        indexed_fasta="results/{strain}/{assembler}/polishing/Racon_round4.fasta",
    run:
        shell("bwa mem {params.indexed_fasta} {input.R1} {input.R2} > {output.mapping}"),

#BWA mem alignment is sorted for input into Pilon
rule sorted_BWA_mem_alignment: 
    input:
        mapping="results/{strain}/{assembler}/polishing/Illumina_mapping.sam",   
    output:
        "results/{strain}/{assembler}/polishing/Illumina_mapping_sorted.bam",   
    run:
        shell("samtools sort {input.mapping} -o {output}"),

#BWA mem sorted alignment is indexed 
rule indexed_BWA_mem_alignment:
    input:
        "results/{strain}/{assembler}/polishing/Illumina_mapping_sorted.bam",
    output:
        touch("results/{strain}/{assembler}/polishing/BWA_alignment_makeidx.done"),
    run:
        shell("samtools index {input}"),

#Pilon polish the genome using Illumina short reads
rule Pilon_polish: 
    input:
        assembly="results/{strain}/{assembler}/polishing/Racon_round4.fasta",
        bam="results/{strain}/{assembler}/polishing/Illumina_mapping_sorted.bam",
        idxdone="results/{strain}/{assembler}/polishing/BWA_alignment_makeidx.done",
    output:
        pilon="results/{strain}/{assembler}/polishing/pilon/pilon.fasta",
    params:
        outdir="results/{strain}/{assembler}/polishing/pilon/",
    run:
        shell("pilon --genome {input.assembly} --bam {input.bam} --outdir {params.outdir}"),

#Index the pilon polished output
rule BWA_index_pilon: 
    input:
        "results/{strain}/{assembler}/polishing/pilon/pilon.fasta",
    output:
        touch("results/{strain}/{assembler}/polishing/pilon/makeidx.done"),
    run:
        shell("bwa index {input}"),

#BWA mem is used to align the Illumina reads to the indexed Pilon output
#For this step the Illumina reads need to be combined, this is NOT included in this snakefile
rule Illumina_on_pilon_alignment: 
    input:
        combined_illumina="illumina/{strain}_illumina.fastq.gz",
        idxdone_pilon="results/{strain}/{assembler}/polishing/pilon/makeidx.done",
    output:
        temp("results/{strain}/{assembler}/polishing/pilon/all_illumina_mapping_{strain}.sam")
    params:
        indexed_fasta="results/{strain}/{assembler}/polishing/pilon/pilon.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.combined_illumina} > {output}"

#Polish the Pilon output with Racon using Illumina reads. Generates the Final polished genome, output to the main directory
rule Racon_illumina:
    input:
        combined_illumina="illumina/{strain}_illumina.fastq.gz",
        alignment="results/{strain}/{assembler}/polishing/pilon/all_illumina_mapping_{strain}.sam",
        assembly="results/{strain}/{assembler}/polishing/pilon/pilon.fasta"
    output:
        Racon_illumina="results/{strain}/{assembler}/polished_genome.fasta"
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.combined_illumina} {input.alignment} {input.assembly} > {output.Racon_illumina}")

#Making Final Genome stats file
rule genome_stats: 
        input:
            "results/{strain}/{assembler}/polished_genome.fasta"
        output:
            "results/{strain}/{assembler}/genome_stats.txt"
        shell:
            "seqkit stats -T {input} > {output}" 

rule pooled_stats:
        input:
            "results/{strain}/{assembler}/genome_stats.txt"
        output:
            touch("results/{strain}/{assembler}/genome_stats.pooled"),
        params:
            "results/{assembler}_genome_stats.txt"
        shell:
           "cat {input} >> {params}"

####Assembly QC 

###Coverage and Concordent read assesment
#Index the polished genome for mapping metrics
rule BWA_index: 
    input:
        "results/{strain}/{assembler}/polished_genome.fasta"
    output:
        touch("results/{strain}/{assembler}/polished_genome.IDXDONE")
    shell:
        "bwa index {input}" 

#Mapping all illumina reads to the genome
rule BWA_mem_illumina_all: 
    input:
        R1="illumina/{strain}_R1.fastq.gz", 
        R2="illumina/{strain}_R2.fastq.gz",
        idxdone="results/{strain}/{assembler}/polished_genome.IDXDONE"
    output:
        temp("results/{strain}/{assembler}/QC/all_illumina_mapping.sam")
    params:
        indexed_fasta="results/{strain}/{assembler}/polished_genome.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.R1} {input.R2} > {output}"

#Mapping the 10,000 withheld Illumina reads to the genome
rule BWA_mem_illumina_WH: 
    input:
        R1="illumina_subsampled/{strain}_R1_WH.fastq.gz", 
        R2="illumina_subsampled/{strain}_R2_WH.fastq.gz",
        idxdone="results/{strain}/{assembler}/polished_genome.IDXDONE"
    output:
        temp("results/{strain}/{assembler}/QC/WH_illumina_mapping.sam")
    params:
        indexed_fasta="results/{strain}/{assembler}/polished_genome.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.R1} {input.R2} > {output}"

#Mapping nanopore reads to the genome
rule BWA_mem_nanopore: 
    input:
        reads="ONT/{strain}.fastq.gz",
        idxdone="results/{strain}/{assembler}/polished_genome.IDXDONE"
    output:
        temp("results/{strain}/{assembler}/QC/nanopore_mapping.sam")               
    params:              
        indexed_fasta="results/{strain}/{assembler}/polished_genome.fasta"
    shell:
        "bwa mem {params.indexed_fasta} {input.reads}  > {output}"

#Sorting the alingment files to allow extraction of coverage and concordant mapping data
rule samtools_sort_all:
    input:
        "results/{strain}/{assembler}/QC/all_illumina_mapping.sam"
    output:
        "results/{strain}/{assembler}/QC/all_illumina_mapping_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_sort_WH:
    input:
        "results/{strain}/{assembler}/QC/WH_illumina_mapping.sam"
    output:
        "results/{strain}/{assembler}/QC/WH_illumina_mapping_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_sort_nanopore:
    input:
        "results/{strain}/{assembler}/QC/nanopore_mapping.sam" 
    output:
        "results/{strain}/{assembler}/QC/nanopore_mapping_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

#Extracting coverage stats for the alignments.
rule samtools_coverage_illumina:
    input:
        "results/{strain}/{assembler}/QC/all_illumina_mapping_sorted.bam"
    output:
        "results/{strain}/{assembler}/QC/all_illumina_mapping_coverage.txt"
    shell:
        "samtools depth {input} > {output}"

rule samtools_coverage_nanopore:
    input:
        "results/{strain}/{assembler}/QC/nanopore_mapping_sorted.bam"
    output:
        "results/{strain}/{assembler}/QC/all_nanopore_mapping_coverage.txt"
    shell:
        "samtools depth {input} > {output}"

#Calculating concordant reads for Illumina alignments        
rule illumina_WH_mapping_concordent_reads:
    input:  
        "results/{strain}/{assembler}/QC/WH_illumina_mapping_sorted.bam"
    output:
        "results/{strain}/{assembler}/QC/illumina_WH_mapping_concordent_reads.txt"
    shell:
        "samtools view -c -f 0x2 {input} > {output}"

rule illumina_all_mapping_concordent_reads:
    input:  
        "results/{strain}/{assembler}/QC/all_illumina_mapping_sorted.bam"
    output:
        "results/{strain}/{assembler}/QC/illumina_all_mapping_concordent_reads.txt"
    shell:
        "samtools view -c -f 0x2 {input} > {output}"

rule illumina_WH_non_mapping_concordent_reads:
    input:
        "results/{strain}/{assembler}/QC/WH_illumina_mapping_sorted.bam"
    output:
        "results/{strain}/{assembler}/QC/illumina_WH_non_mapping_concordent_reads.txt"
    shell:
        "samtools view -c -f 0x12 {input} > {output}"

#Make one combined file of concordent reads for all strains 
rule Illumina_concoordent_reads_file:
    input:
        "results/{strain}/{assembler}/QC/illumina_WH_mapping_concordent_reads.txt"
    output:
        touch("results/{strain}/{assembler}/QC/concordent_mapping.txt")
    params:
        "{assembler}_concordent_mapping.txt"
    run:
        shell("echo {input} >> {params}"),
        shell("cat {input} >> {params}")

rule Illumina_non_mapping_reads_file:
    input:
        "results/{strain}/{assembler}/QC/illumina_WH_non_mapping_concordent_reads.txt"
    output:
        touch("results/{strain}/{assembler}/QC/non_mapping.txt")
    params:
        "{assembler}_non_mapping.txt"
    run:
        shell("echo {input} >> {params}"),
        shell("cat {input} >> {params}")

#Plotting coverage stats
rule plot_illumina_coverage:
    input:
        "results/{strain}/{assembler}/QC/all_illumina_mapping_coverage.txt"
    output:
        "results/{strain}/{assembler}/QC/all_illumina_mapping_coverage.pdf"
    shell:
        "Rscript ./scripts/coverage_plots_g.R {input}"

rule plot_nanopore_coverage:
    input:
        "results/{strain}/{assembler}/QC/all_nanopore_mapping_coverage.txt"
    output:
        "results/{strain}/{assembler}/QC/all_nanopore_mapping_coverage.pdf"
    shell:
        "Rscript ./scripts/coverage_plots_g.R {input}"

###Socru assesment
#Run socru assesment on assemblies 
rule socru:
    input:
        "results/{strain}/{assembler}/polished_genome.fasta"
    output:
        "results/{strain}/{assembler}/QC/socru.txt"
    shell:
        "socru Escherichia_coli {input} > {output}"

#Make one combined file of Socru output for all strains
rule Scoru_file:
    input:
        "results/{strain}/{assembler}/QC/socru.txt"
    output:
        touch("results/{strain}/{assembler}/QC/socru.txt")
    params:
        "{assembler}_socru_mapping.txt"
    run:
        shell("echo {input} >> {params}"),
        shell("cat {input} >> {params}")

###ORF assesment
#Run prodigal step of ORF assesment
rule ORF_prodigal:
        input: 
            "results/{strain}/{assembler}/polished_genome.fasta"
        output: 
            "results/{strain}/{assembler}/QC/ORF.faa"
        benchmark: 
            "results/{strain}/{assembler}/benchmarks/{strain}.prodigal.bench.txt"
        shell: 
            "prodigal -a {output} -q -i {input}"

#Run Diamond step of ORF assesment
rule ORF_diamond:
        input: 
            "results/{strain}/{assembler}/QC/ORF.faa"
        output: 
            "results/{strain}/{assembler}/QC/ORF.data"
        threads: 
            32
        params:
                db="/data/databases/diamond/uniprot.dmnd",
                of="6 qlen slen"
        benchmark: 
            "results/{strain}/{assembler}/benchmarks/{strain}.diamond.bench.txt"
        shell: 
            "diamond blastp --threads 32 --max-target-seqs 1 --db {params.db} --query {input} --outfmt {params.of} --out {output}"

#Run plotting step of ORF assesment
rule ORF_hist:
        input: 
            "results/{strain}/{assembler}/QC/ORF.data"
        output: 
            "results/{strain}/{assembler}/QC/ORF.pdf"
        shell: 
            "R --slave --no-restore --file=./scripts/hist_orf_lengths.R --args {input} {output}"

#Run plasmid finder program 
rule plasmid_finder:
        input:
            "results/{strain}/{assembler}/polished_genome.fasta"
        output:
            touch("results/{strain}/{assembler}/QC/plasmidfinder.marker")
        params:
            "results/{strain}/{assembler}/QC/"
        shell: 
            "plasmidfinder.py -i {input} -o {params} -p /data/databases/plasmidfinder_db/"

#Make CSV file of results 
rule plasmid_finder_csv:
        input:
            "results/{strain}/{assembler}/QC/plasmidfinder.marker"
        params:
            "results/{strain}/{assembler}/QC/data.json"
        output:
            touch("results/{strain}/{assembler}/QC/plasmid_finder_data_csv")
        shell:
            "python ./scripts/plasmid_finder_json_parsing.py {params}"


rule mlplasmids:
        input: 
            "results/{strain}/{assembler}/polished_genome.fasta"
        output:
            "results/{strain}/{assembler}/QC/mlplasmids_tab"
        shell:
            "Rscript mlplasmids/scripts/run_mlplasmids.R {input} {output} 0.01 'Escherichia coli'"

rule mlplasmids_by_assembler:
        input:
            "results/{strain}/{assembler}/QC/mlplasmids_tab"
        output:
            touch("results/{strain}/{assembler}/QC/mlplasmids_tab_merged")
        shell: 
            "python scripts/merge_mlplasmids.py {input}"

rule make_tab_output:
        input: 
            ORF= "results/{strain}/{assembler}/QC/ORF.data",
            genome_stats= "results/{strain}/{assembler}/genome_stats.txt",
            mlplasmids= "results/{strain}/{assembler}/QC/mlplasmids_tab",
            plasmid_finder_marker="results/{strain}/{assembler}/QC/plasmid_finder_data_csv",
            concordant_reads="results/{strain}/{assembler}/QC/illumina_WH_mapping_concordent_reads.txt",
            non_mapping_reads="results/{strain}/{assembler}/QC/illumina_WH_non_mapping_concordent_reads.txt",
            socru="results/{strain}/{assembler}/QC/socru.txt"
        output: 
            touch("results/final_table/{strain}_{assembler}_all_metrics_merged.marker")
        params: 
            plasmid_finder="results/{strain}/{assembler}/QC/plasmid_finder.csv"

        shell:
            "python scripts/All_parsing.py {input.ORF} {input.genome_stats} {input.mlplasmids} {params.plasmid_finder} {input.concordant_reads} {input.non_mapping_reads} {input.socru}"


