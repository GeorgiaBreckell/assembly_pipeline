# ONT/Illumina hybrid assembly pipeline

This is a snakemake based assembly pipeline for assembling, polishing and assessing hybrid genomes from Illumina short reads and Nanopore long reads. Developed in an effort to establish which assembly option was best suited for producing a high quality de-novo genome assembly for natural isolates of E.coli. 

The pipeline takes as input paired Illumina short reads, along with Nanopore long reads for each isolate and produces up to 5 hybrid assemblies, using the assemblers Flye, Canu, Ra, Wtdbg2 and Unicycler. 

Prior to assembly a small subset of Nanopore and Illumina reads are removed from the dataset. These reads are withheld from the assemblies and are used during the assembly assessment. We believe using these *"novel"* reads when assessing read coverage and the fraction of concordantly mapping reads improves the validity of these metrics.  

Each assembly is polished iteratively by Racon and Pilon using both the long and short reads. We therefore consider all final assemblies to be hybrid assemblies, regardless of whether the original assembly was derived entirely from long reads. 

The quality of each polished assembly is then assessed across a range of metrics including: 
* The fraction of short open reading frames
* The fraction of concordantly mapping illumina reads
* A Socru assessment of rRNA operon arrangement, to assess global assembly structure
* A mlPlasmids assessment to predict short contigs of plasmid origin 
* Use of plasmidFinder to identify plasmids within isolates 
* Read coverage plots for Nanopore and Illumina withheld reads.

Tested and used on a Linux system  

## How to use this pipeline

1. Clone the repository where you will run the pipeline 
```
git clone https://github.com/GeorgiaBreckell/assembly_pipeline.git
cd assembly_pipeline
```
2. Set up the environment and install required tools
```
conda create --name assemblies --file environments/base.yaml
conda activate assemblies
``` 
   - Ra must be installed into the working dir seperately from [here](https://github.com/lbcb-sci/ra)
   - mlplasmids is an R based tool installed from [here](https://gitlab.com/sirarredondo/mlplasmids), tested in R v4.0.0

3. Write your config file
``` 
cp assemblies_config.yml.example assemblies_config.yml
```
Add appropriate details to the config file for your isolates,

ie for *sample_ID*.fastq add *sample_ID* to the strains section of the config. 

An output dir will then be produced for each *sample_ID* containing seperate results dir for each assembler run. 

Assemblers can be removed from the config file if a particular assembler is not required. 

4. Start the Snakefile
```
snakemake -s assemblies.Snakefile --use-conda --cores 
```
Specifiy the number of cores you have available for the pipeline after the --cores option. 

A dry run (initiated with the -np flag) can be useful to check the installation and config file is correct. 

