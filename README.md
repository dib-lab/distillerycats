# distillerycats
Distill sourmash with spacegraphcats.
For metaanalysis of metagenomic datasets.
Currently takes as input metagenome sequencing reads from multiple studies, the study from which each sample originated, and a classification variable of interest. 

This repository is under active development. 
The code in the master branch is not guarunteed to be completely documented or functional yet. 
We plan to remove this message soon when this is no longer true :)

## Assumptions made by the pipeline

Some or all of the assumptions may be removed in future iterations

+ Assumes paired-end shotgun metagenome reads
+ Assumes data are already downloaded and in `inputs/raw`
+ Currently focused on human microbiome samples
  + Removes human reads as contaminant host
  + Classifies sequences against [human microbiome MAG databases](https://osf.io/hza89/)
+ Requires samples from at least two studies.
+ Reads must be > 31 bp in length
+ Requires snakemake-minimal, sourmash, feather, and pandas.

## Running the pipeline
```
git clone https://github.com/dib-lab/distillerycats.git
cd distillerycats
conda env create -f environment.yml -n dcats
conda activate dcats
```
