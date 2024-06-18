# Microexon Code

## Installation

### Prerequisites

#### Hardware
- A x64 machine running Linux (tested on Ubuntu 18.04). Microexon Code does not run on macOS or Windows.
- An Nvidia GPU capable of CUDA 9.0+
- An internet connection as some required data is downloaded on demand
- If you want to run the code on large-scale data, you may need to run the code on a cluster or in the cloud. Snakemake is compatible with many schedulers and cloud APIs. See https://snakemake.readthedocs.io for details.

#### Software
Install Snakemake according to the instructions at https://snakemake.readthedocs.io/en/stable/getting_started/installation.html. This requires a Conda-based Python3 installation.

In short:

```shell
$ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
$ bash Miniforge3-$(uname)-$(uname -m).sh
```

Next, install Snakemake into a fresh environment and activate it:
```shell
$ mamba create -c conda-forge -c bioconda -n snakemake=8.14 snakemake
$ mamba activate snakemake
```

You will also need to add the `microexon_code/` directory to the local PYTHONPATH so the local packages can be found:
```shell
$ export PYTHONPATH=`pwd`  # or substitute the actual path to microexon_code/
```

## Overview
Microexon Code is a lightweight machine learning model using gradient boosted trees that is built on top of third-party feature detectors. To train the model or for inference, features need to be computed on the inputs. We use Snakemake to orchestrate this workflow. Some components require some inputs (e.g. model snapshots) which are supplied in the `resources/` folder.

The scripts supplied in this workflow can be used to train the model (even though we provide fully trained models that can be used for inference), replicate the analysis presented in the manuscript (such as _in silico_ saturation mutagenesis) and run inference on use supplied datasets.

## Usage
Microexon Code is run by executing Snakemake commands. Snakemake will then build a computational graph and run only those steps of the workflow required to obtain the desired output. The output can be specified by either supplying a named rule or supplying an _output_ file name.

### Inference for wild-type sequences or individual variants
In this mode, you may score sequences from a reference genome and optionally apply individual variants to each sequence. This is designed for a small number of variants which are applied to the wild-type sequence on the fly. Each input can have different variants applied. We use this mode to score SNVs such as our saturation _in silico_ mutagenesis experiments. 

To run inference on your own inputs you need to follow these steps exactly. The data format as well as the file names and paths must be exactly as specified as otherwise, Snakemake cannot deduce the correct workflow. For examples of correct inputs, see `resources/inputs/`

Create an input file in gzipped tab format with the following columns:
  - `event`: An event ID you can choose
  - `chrom`: Chromosome in UCSC format (e.g. "chr5", "chrX")
  - `strand`: `+` or `-`
  - `upIntStart`, `upIntEnd`, `dnIntStart`, `dnIntEnd`: Genomic coordinates of the upstream and downstream introns of the microexon, 1-based indexing
  - `geneName`: Gene name
  - `lengthDiff`: The length of the microexon
  - `label`: Not used for inference, can be set to any value
  - `variant` (optional column): Comma-delimited list of variants (e.g. `chr10:26770836:T:A`) to apply to the sequence. Variants must be non-overlapping. 

Store the input file under `resources/inputs/{task}/{prefix}.{species}.tab.gz`, where you can choose your own values for `{task}` and `{prefix}` and specify the desired genome build for `{species}`. Do not use any dots (`.`) in the prefix. An example of a valid path is `resources/inputs/wt/MicEvents.hg38.tab.gz`.
 
Run inference with the following command:
```shell
  $ snakemake --cores {N} --use-conda --conda-not-block-search-path-envvars results/predictions/{model_name}/{task}/predictions.{prefix}.{species}.csv.gz
``` 
Substitute `{N}` for the number of cores on your machine or the number of parallel jobs you'd like to run

Substitute the values for `{task}`, `{prefix}` and `{species}` from your input file. For `{model_name}`, choose which model you want to use for the predictions:
  - `known_microexons`: This model was trained on known microexons from seven mammalian species
  - `known_microexons_hg38`: This model was only trained on human known microexons
  - `novel_microexons`: This model was trained known microexons from seven mammalian species and the human novel microexons we detected in our study

As an example, you can run predictions for the known human wild-type microexons with the command
```shell
snakemake --use-conda --conda-not-block-search-path-envvars results/predictions/known_microexons/wt/predictions.MicEvents.hg38.csv.gz
```
