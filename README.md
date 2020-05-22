# ONT/Illumina hybrid assembly pipeline

This is a snakemake based assembly pipeline for assembling, polishing and assessing hybrid genomes from Illumina short reads and Nanopore long reads. Developed in an effort to establish which assembly option was best suited for producing a high quality de-novo genome assembly for natural isolates of E.coli. 

The pipeline takes as input paired Illumina short reads, along with Nanopore long reads for each isolate and produces up to 5 hybrid assemblies, using the assemblers Flye, Canu, Ra, Wtdbg2 and Unicycler. Each assembly is polished iteratively by Racon and Pilon using both the long and short reads. We therefore consider all final assemblies to be hybrid assemblies, regardless of whether the original assembly was derived entirely from long reads. 

The quality of each polished assembly is then assessed across a range of metrics including: 
* The fraction of short open reading frames
* The fraction of concordantly mapping illumina reads
* A Socru assessment of rRNA operon arrangement, to assess global assembly structure
* A mlPlasmids assessment to predict short contigs of plasmid origin 
* Use of plasmidFinder to identify plasmids within isolates 
* Read coverage plots across the genome

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
   - Ra and Wtdbg2 must be installed into the working dir seperately from here: [ra](https://github.com/lbcb-sci/ra), [wtdbg2](https://github.com/ruanjue/wtdbg2)  
   - mlplasmids is an R based tool installed from [here](https://gitlab.com/sirarredondo/mlplasmids), tested in R v4.0.0
   - If run on a linux system, Bandage plots will be produced for Flye, Unicycler and Canu assemblies. For OSX, the Bandage GUI can be installed from [here](http://rrwick.github.io/Bandage/) and is useful for visualising contig relationships. 

3. Write your config file
``` 
cp assemblies_config.yml.example assemblies_config.yml
```
Add appropriate details to config file for your isolates:
