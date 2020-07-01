# distillerycats
Distill sourmash with spacegraphcats!

Find disease associations across metagenomes with k-mers using sourmash, and then recover pangenomic accessory elements using spacegraphcats.

distillerycats currently takes as input metagenome sequencing reads from multiple studies, the study from which each sample originated, and a classification variable of interest.

**Warning:** This repository is under active development. 
The code in the master branch is not guaranteed to be completely functional yet, and our documentation is still being written.
We will update this message as we make progress :)

## Assumptions made by the pipeline

Some or all of the assumptions may be removed in future iterations

+ Assumes paired-end shotgun metagenome reads
+ Assumes data are already downloaded and in `inputs/raw`
+ Currently focused on human microbiome samples
  + Removes human reads as contaminant host
  + Classifies sequences against [human microbiome MAG databases](https://osf.io/hza89/)
+ Requires samples from at least two studies.
+ Reads must be > 31 bp in length.
+ Requires snakemake-minimal, [sourmash](http://sourmash.rtfd.io/), [feather](https://blog.rstudio.com/2016/03/29/feather/), and pandas.

## Getting started

This pipeline uses conda for software and environment management. To get started, install [miniconda](https://docs.conda.io/en/latest/miniconda.html). If you're new to miniconda, see [this tutorial](https://angus.readthedocs.io/en/2019/conda_tutorial.html). 

To get started, clone this repository. All computation will take place within the cloned github directory.:

```
git clone https://github.com/dib-lab/distillerycats.git
```

`cd` into the directory and create the environment:

```
cd distillerycats
conda env create -f environment.yml -n dcats
conda activate dcats
```

Make sure your input data is in the `inputs/raw` directory. This pipeline assumes that all input data paths follow this format:

```
inputs/raw/{sample}_R1.fq.gz
inputs/raw/{sample}_R2.fq.gz
```

Where all input samples are paired-end reads with and `R1` and `R2` file, reads are gzipped, and files end with `.fq.gz`.` 

The `{sample}` root of the input files should match `sample` column in the metadata file:
```
sample,study,var
PSM7J199,iHMP,CD
PSM7J1BJ,iHMP,CD
PSM7J177,iHMP,CD
```

The metadata file must be located in the `inputs` directory, and should be called `test_metadata.csv`. It needs to be in csv format, and should have the sample columns `sample`, `study`, and `var`. The names must match exactly in capitalization/spelling. 

If you would like to try distillerycats on a (large) test data set, we have uploaded one to OSF. 

To run the full pipeline including read preprocessing, including adapter trimming, human DNA removal, and k-mer abundance trimming, download this following dataset, untar it, and make sure the fastq files are in `inputs/raw` as below. Note that k-mer trimming and human DNA removal both take ~64GB of RAM. 

```
mkdir -p inputs/raw
wget -O inputs/test_data.tar.gz https:...
tar xf inputs/test_data.tar.gz -C inputs/raw
```

If you would like to skip read preprocessing but would still like to try the pipeline, you can either download the preprocessed reads, or you can download the signatures. These files should be placed in `ouputs/abundtrim` and `outputs/sigs`, respectively. The first step in the pipeline after preprocessing is to calculate signatures, and signatures are much lighter weight than preprocessed fastq files. 

preprocessed reads:

```
mkdir -p outputs
wget -O outputs/abundtrim.tar.gz ...
tar xf outputs/abundtrim.tar.gz -C outputs/
```

signatures:

```
mkdir -p outputs
wget -O outputs/sigs.tar.gz https://osf.io/scbk8/
tar xf outputs/sigs.tar.gz -C outputs/
```

---

@taylorreiter
@ctb
