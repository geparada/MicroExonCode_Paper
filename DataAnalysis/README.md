# Data Analysis for Microexon Code Paper

This directory contains the scripts used for data analysis in the study "An expanded repertoire of brain microexons is directly impacted by autism-associated genetic variation."

## Contents

### Notebooks

All the original Jupyter notebooks used to conduct the data analysis are stored inside the `Notebooks/` directory. These notebooks contain the core workflow and code used in the analysis process.

### Setup

To run these notebooks, you can create a virtual environment that contains most of the necessary dependencies using the following command:

```sh
mamba create --name test numpy pandas pybedtools pybigwig pysam r-data.table r-dbplyr r-ggplot2 r-ggsignif r-reshape rpy2 scipy seaborn bedtools bioconductor-biobase bioconductor-annotationdbi bioconductor-iranges bioconductor-pcamethods r-cowplot
```

An exhaustive list of all the libraries included in the virtual environment we used to carry out these analyses and their corresponding versions is found in the `env/environment_clone.yml` file.

###  Input Tables
Some input tables necessary to run these notebooks are available in the `Input/` directory. Additional files can be shared upon request.


