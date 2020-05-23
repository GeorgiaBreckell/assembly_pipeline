##################################################################################
###     Assembly pipeline for hybrid assemblies of ONT and illumina data
###
###     Part 1 of 3 To be used with assembly_polishing.Snakefile and 
###     assesment.Snakefile. 
###     
###     Produces hybrid assemblies using Unicycler V0.4.4
###     long read only assemblies are produced using:
###     Canu, Flye, Ra, Wtdbg2 (redbean)                        
###     The first step of this pipeline,subsets both the ONT 
###     and Illumina data, the majority of reads are used in 
###     the assembly, while withheld reads are used for mapping
###     during QC in part 3 (assesment.Snakefile). 
###     
###
###     Reads should be in a directoy named ONT and Illumina respectively
###     and output will be in an individual directory for each strain, 
###     organsied by assembler. A fasta file, GFA file and bandage plot
###     are produced from running this snakemake. 
###
###     Assemblers and Strains can be specifed in the Config.yml file 
###         
###     Georgia Breckell  19.11.2019
###
###
###
###################################################################################

configfile:
    "assemblies_config.yml"

rule all:
    input:
        expand("{strain}/{assembler}/assembly.fasta", strain=config["strain"], assembler=config["assembler"]),
        expand("{strain}/{assembler}/output_dir_removed", strain=config["strain"],assembler=config["assembler"]),
        expand("{strain}/{assembler}/bandage_plot.png", strain=config["strain"],assembler=config["assembler_gfa"]), 

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
        fasta="{strain}/unicycler/unicycler_output/assembly.fasta",
        gfa="{strain}/unicycler/unicycler_output/assembly.gfa"
    params:
        out_prefix="{strain}/unicycler/unicycler_output/"
    log:
        "{strain}/logs/unicycler.log"
    benchmark:
        "{strain}/benchmarks/unicycler.assembly.benchmark.txt"
    threads: 8
    shell:
        "unicycler -1 {input.r1} -2 {input.r2} -l {input.ont} -o {params.out_prefix}" 
        " 1> {log}"

#Copies the unicycler assembly out of the unicycler output dir
rule copy_unicycler:
    input:
        assembly="{strain}/unicycler/unicycler_output/assembly.fasta",
        gfa="{strain}/unicycler/unicycler_output/assembly.gfa"
    output:
        assembly="{strain}/unicycler/assembly.fasta",
        gfa="{strain}/unicycler/assembly.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"), 
        shell("cp {input.gfa} {output.gfa}")   

#remove unicycler assembly directory   
rule remove_unicycler_output_dir:
    input:
        assembly="{strain}/unicycler/assembly.fasta",
        gfa="{strain}/unicycler/assembly.gfa"
    output:
        touch("{strain}/unicycler/output_dir_removed")
    params:
        output_dir="{strain}/unicycler/unicycler_output/"
    shell:
        "rm -r {params.output_dir}"

#Canu long read only assembly with the "Assembly dataset of reads"
rule canu_assembly:
    input:
        "ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        fasta="{strain}/canu/canu_output/{strain}.contigs.fasta",
        gfa="{strain}/canu/canu_output/{strain}.unitigs.gfa"
    log:
        "{strain}/logs/canu.log"
    benchmark:
        "{strain}/benchmarks/canu.assembly.benchmark.txt"
    params:
        directory="{strain}/canu/canu_output/",
        prefix="{strain}"
    shell:
        "canu -d {params.directory} -p {params.prefix} genomeSize=5m -nanopore-raw {input} -maxMemory=32g -maxThreads=10 1>{log}" 

#Copies and renames the canu output fasta to generic naming convention 
rule copy_canu:
    input:
        assembly="{strain}/canu/canu_output/{strain}.contigs.fasta",
        gfa="{strain}/canu/canu_output/{strain}.unitigs.gfa"
    output:
        assembly="{strain}/canu/{strain}.contigs.fasta",
        gfa="{strain}/canu/{strain}.unitigs.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"),
        shell("cp {input.gfa} {output.gfa}")

rule rename_canu:
    input:
        assembly="{strain}/canu/{strain}.contigs.fasta",
        gfa="{strain}/canu/{strain}.unitigs.gfa"
    output:
        assembly="{strain}/canu/assembly.fasta",
        gfa="{strain}/canu/assembly.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"),
        shell("cp {input.gfa} {output.gfa}")

#remove canu assembly directory   
rule remove_canu_output_dir:
    input:
        assembly="{strain}/canu/assembly.fasta",
        gfa="{strain}/canu/assembly.gfa"
    output:
        touch("{strain}/canu/output_dir_removed")
    params:
        output_dir="{strain}/canu/canu_output/"
    shell:
        "rm -r {params.output_dir}"

#ra long read only assembly with the "Assembly dataset of reads", no renmaing needed as output is in desired convention. 
#ra output as single fasta file so no movement out of the "output_dir" is needed either 
rule ra_assembly:
    input:
        "ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        assembly="{strain}/ra/assembly.fasta",
	dir_removed=touch("{strain}/ra/output_dir_removed")
    #conda:
    #    "python2_7.yml"
    log:
        "{strain}/logs/ra.log"
    benchmark:
        "{strain}/benchmarks/ra.assembly.benchmark.txt"
    run:
       shell( "./ra/build/bin/ra -x ont -t 4 {input} > {output.assembly} 1>{log}")

#flye long read only assembly with the "Assembly dataset of reads", no renaming is needed because flye output is in desired convention 
rule flye_assembly:
    input:
        "ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        fasta="{strain}/flye/flye_output/assembly.fasta",
        gfa="{strain}/flye/flye_output/assembly_graph.gfa"
    params:
        out_prefix="{strain}/flye/flye_output/"
    log:
        "{strain}/logs/flye.log"
    benchmark:
        "{strain}/benchmarks/flye.assembly.benchmark.txt"
    conda:
        "environments/assemblies_2_7.yml"
    shell:
        "flye --nano-raw {input} --genome-size 5g --out-dir {params.out_prefix} --plasmids 1>{log}"

#Copies the flye assembly out of the flye output dir
rule copy_flye:
    input:
        assembly="{strain}/flye/flye_output/assembly.fasta",
        gfa="{strain}/flye/flye_output/assembly_graph.gfa"
    output:
        assembly="{strain}/flye/assembly.fasta",
        gfa="{strain}/flye/assembly_graph.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"),
        shell("cp {input.gfa} {output.gfa}")

#rename GFA file to standard convention
rule rename_flye:
    input:
        gfa="{strain}/flye/assembly_graph.gfa"
    output:
        "{strain}/flye/assembly.gfa"
    run:
        shell("mv {input.gfa} {output}")

#rule to remove flye assembly directory   
rule remove_flye_output_dir:
    input:
        assembly="{strain}/flye/assembly.fasta",
    output:
        touch("{strain}/flye/output_dir_removed")
    params:
        output_dir="{strain}/flye/flye_output/"
    shell:
        "rm -r {params.output_dir}"

#Wtdbg2(redbean) long read only assembly with assembly dataset
rule wtdbg2_assembly_1:
        input:
            ont="ONT_subsampled/{strain}_assembly.fastq.gz"
        output:
            "{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.lay.gz"
        params:
            out_prefix="{strain}/wtdbg2/wtdbg2_output/{strain}"
        log:
            "{strain}/log/wtdgb2_1.log"
        benchmark:
            "{strain}/benchmarks/wtdbg2_1.assembly.benchmark.txt"
        shell:
            "wtdbg2 -x ont -g 4.8m -t 16 -i {input.ont} -o {params.out_prefix} 1>{log}"

#Second half of the redbean assembly
rule wtdbg2_assembly_2:
        input:
            wtdbg2_config="{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.lay.gz"
        output:
            "{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.fa" 
        params:
            out_prefix="{strain}/wtdbg2/wtdbg2_output/{strain}"
        log:
            "{strain}/log/wtdgb2_2.log"
        benchmark:
            "{strain}/benchmarks/wtdbg2_2.assembly.benchmark.txt"
        shell:
            "wtpoa-cns -t 16 -i {input} -o {params.out_prefix}.ctg.fa 1>{log}"

#Copies and renames the wtdbg2 output fasta to generic naming convention 
rule copy_wtdbg2_1:
    input:
        "{strain}/wtdbg2/wtdbg2_output/{strain}.ctg.fa"
    output:
        "{strain}/wtdbg2/{strain}.ctg.fa" 
    shell:
        "cp {input} {output}"

rule rename_wtdbg2_1:
    input:
        "{strain}/wtdbg2/{strain}.ctg.fa"
    output:
        "{strain}/wtdbg2/assembly.fasta"
    shell:
        "mv {input} {output}"

#rule to remove wtdbg2 assembly directory   
rule remove_wtdbg2_output_dir:
    input:
        assembly="{strain}/wtdbg2/assembly.fasta",
    output:
        touch("{strain}/wtdbg2/output_dir_removed")
    params:
        output_dir="{strain}/wtdbg2/wtdbg2_output/"
    shell:
        "rm -r {params.output_dir}"

#Bandage plots of the assembly graph for each assembler (except Ra and wtdbg2 - no .gfa file produced) 
rule Bandage_plot:
    input:
        gfa="{strain}/{assembler}/assembly.gfa"
    output:
        bandage_plot="{strain}/{assembler}/bandage_plot.png"
    run:
        shell("Bandage image {input.gfa} {output.bandage_plot}")






