# bigmelon
Memory Efficient Methods for DNA Methylation Analysis

## Installation
Currently there is a bug in the release version of wateRmelon that has yet to be pushed to bioConductor. This bug prevents bigmelon from parsing in EPIC idat files. As a result we recommend installing wateRmelon (an important dependency straight from github.

```
library(devtools)
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19') # Optional
devtools::install_git('https://github.com/schalkwyk/wateRmelon')
BiocManager::install('bigmelon')
```


## Authors

Tyler J. Gorrie-Stone and Leonard C. Schalkwyk

School of Biological Sciences
University of Essex
Colchester, UK

## Abstract
DNA methylation analyses are getting ever bigger.With the release of the HumanMethylationEpic microarray by Illumina and datasets reaching into the thousands, analysis of these large datasets using popular R packages is becoming impractical due to memory requirements and even the time required to read the data from disk. As such there is an increasing need for computationally efficient methods to perform meaningful analysis on high dimension data.

The bigmelon R package provides a memory-efficient work-flow that enables users to perform complex, large scale analyses required in EWAS without huge RAM. Building on the CoreArray Genome Data Structure (.gds) file format and libraries packaged in 'gdsfmt', we provide a familiar wateRmelon-like work flow that facilitates reading-in, preprocessing, quality control and statistical analysis.

To demonstrate large-scale data analyses, we stored the entire contents of the marmal-aid database (>14,500 samples) in a .gds file and demonstrate quality measures and principal components analysis.

Overall, bigmelon provides a familiar environment for users to perform large-scale analyses where  convention methods would run out of memory. Bigmelon shows reasonable performance in speed compared to conventional methods.
