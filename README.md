# Phylogenetic pipeline #
# Author:
Aleksandra Cupriak 

## Contents ##
  * [General Description](#general-description)
  * [Requirments](#requirments)
  * [General Information](#general-information)
  * [Usage](#usage)
  * [Examples](#examples)

## General Description ##

Phylogenetic pipeline leading to the calculation of a set of gene trees and a genome tree based on proteomes

## Requirments ##
  * [Python3](https://docs.python-guide.org/starting/install3/linux/)
    * [Biopython](https://biopython.org/)
    * [ete3](http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html)
    * [toytree](https://toytree.readthedocs.io/en/latest/)
  * [R](https://www.r-project.org/) 
    * [ape](https://cran.r-project.org/web/packages/ape/index.html)
    * [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html)
  * [snakemake](https://snakemake.readthedocs.io/en/stable/)
  * [MMseqs2](https://github.com/soedinglab/MMseqs2)
  * [fasturec](https://bitbucket.org/pgor17/fasturec)
  * [clustalo](http://www.clustal.org/omega/)

## General Information ##
  * Pipeline consists of the following steps: proteomes downloading, clustering, clusters selecting, multialignments calculation, neighbor joining trees building, bootsrapping (optionally), supertree building, consensus tree building
  * The correct path to the input file for running the entire pipeline: species/species.txt
  * Input file for the entire pipeline (species/species.txt) in each line must contain the UniProt Proteome ID and organism name (separated by tab)
  * Running the entire pipeline will lead to the calculation of the genome tree/s based on the downloaded proteomes for the species given in the input file species/species.txt
  * Specifying the output file for any step in the pipeline when run snakemake allows to run only a part of the pipeline if the corresponding input files are present
  * Parameters to pipeline can be passed through the configuration file
  * Pipeline should be run from the main directory of the project
  * Binary for fasturec should be located in usr/bin directory

## Configuration file ##
  * The correct template of the configuration file is presented in the file 'configuration.yaml'
  * Parameters specification:
      + names: output file name fragments
      + results: specification of methods for generating a genome trees
      + to_remove: list of organisms that will be excluded from the analysis (before clustering)
      + c: c parameter for MMseq2 clustering ()
      + min_cluster_size: int, minimum acceptable cluster size (cluster selection)
      + max_cluster_size: int, maximum acceptable cluster size (cluster selection); to ignore can be set as False
      + max_repeat: float, the maximum acceptable frequency of the most common organism in a given cluster (cluster selection)
      + orthological: bool, whether the cluster must be otological (cluster selection)
      + choose_random: bool, whether paralogous clusters will be used for the artificial construction of orthologous clusters (cluster selection)
      + excluded_species: list of organisms that will be excluded from the analysis (after clustering)
      + bootstrap: bool, whether bootsrap will be executed
      + bootstrap_thresh: float, the rejection threshold of NJ trees after bootstrap
      + p: float, parameter p to construct the majority consensus tree


## Usage ##
Usage:

To run the entire pipeline:

    % snakemake -c [arg] [snakemake parameters]

To run the pipline for specific output file:

    % snakemake -c [arg] [output] [snakemake parameters]

To know more about usage of snakemake workflow management system read snakemake documentation.
    
## Examples ##
    Running the entire pipeline:

    > snakemake -c8 --configfile configuration.yaml

    Running piepelin only to the clustering step:

    > snakemake -c8  "results/mmseq_res/mmseq_result_allorganisms_all_seqs.fasta" --configfile configuration.yaml