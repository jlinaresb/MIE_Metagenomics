# Machine Learning algorithms reveals country-specific metagenomic taxa from American Gut Project data

This repository includes code based on article published in [MIE](https://efmi.org/2020/12/10/31st-medical-informatics-europe-conference-mie2021-athens-greece/) congress to be held on 29th-31st May, 2021.

DOI will be added when available.


## American Gut Project (AGP) Data
The data used in this work was downloaded from the AGP (ftp://ftp.microbio.me/AmericanGut/). Three files were downloaded: otu table, otu tree and clinical data. [Phyloseq](https://joey711.github.io/phyloseq/) package was used to manage this data and filtered by feces samples.

## Abstract
In recent years, microbiota has become an increasingly relevant factor for the understanding and potential treatment of diseases. A major problem has been the variability that exists between different countries of origin when designing a global treatment. For this reason, in this work, based on the data reported by the largest study of microbioma in the world, a classification model has been developed based on Machine Learning (ML) capable of predicting the country of origin (United Kingdom vs United States) according to metagenomic data. The data were used for the training of a glmnet algorithm and a Random Forest algorithm. Both algorithms obtained similar results (0.698 and 0.672 in AUC, respectively). Furthermore, thanks to the application of a multivariate feature selection algorithm, eleven metagenomic genres highly correlated with the country of origin were obtained. An in-depth study of the variables used in each model is shown in the present work.

# Getting Started
There are two main files in this repository: preprocessing analyis and machine learning code. In addition, utils functions are included in a utils.r file. To run this code, make sure you have the following:

## Prerequisites
The packages we´ve used:

``` r
install.packages(c('mlr', 'h2o', 'FCBF', 'parallelMap', 'BiocManager'))

BiocManager::install("phyloseq")
```

## Authors and affiliations
Jose Liñares Blanco U+00B9
Carlos Fernandez Lozano
Jose A. Seoane
Guillermo Lopez Campos

$$_^1$$ Department of Computer Science and Information Technologies, Faculty of Computer Science, University of A Coruña, CITIC, Campus Elviña s/n, A Coruña, 15071, Spain

## Questions?
If you have any questions, please feel fre to contact (j.linares@udc.es)