# bigmelon
Memory Efficient Methods for DNA Methylation Analysis

## Software status

| Resource:     | Bioconductor        | Travis CI     |
| ------------- | ------------------- | ------------- |
| _Platforms:_  | _Multiple_          | _Linux_       |
| R CMD check   | <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/bigmelon/"><img border="0" src="http://bioconductor.org/shields/build/release/bioc/bigmelon.svg" alt="Build status"></a> (release)</br><a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/bigmelon/"><img border="0" src="http://bioconductor.org/shields/build/devel/bioc/bigmelon.svg" alt="Build status"></a> (devel) | <a href="https://travis-ci.org/tjgorrie/bigmelon"><img src="https://travis-ci.org/tjgorrie/bigmelon.svg" alt="Build status"></a> |

## Installation
We are constantly updating both bigmelon and wateRmelon and will try to keep the github inline with the bioconductor mirror as much as possible. However if you want the latest features before they are committed to bioconductor you can install both wateRmelon and bigmelon using the code below. 

When we introduce functionality into bigmelon this will also be mirrored in wateRmelon and will often require the absolute latest version.

```
library(devtools)
devtools::install_git('https://github.com/schalkwyk/wateRmelon')
devtools::install_git('https://github.com/tjgorrie/bigmelon')
```

## Authors

Tyler J. Gorrie-Stone and Leonard C. Schalkwyk

School of Biological Sciences
University of Essex
Colchester, UK

## Citation
If you use bigmelon for your analyses, please use citation() within R or cite [our paper](10.1093/bioinformatics/bty713) as:

Tyler J Gorrie-Stone, Melissa C Smart, Ayden Saffari, Karim Malki, Eilis Hannon, Joe Burrage, Jonathan Mill, Meena Kumari, Leonard C Schalkwyk, Bigmelon: tools for analysing large DNA methylation datasets, Bioinformatics, Volume 35, Issue 6, 15 March 2019, Pages 981–986, https://doi.org/10.1093/bioinformatics/bty713 

## Abstract
DNA methylation analyses are getting ever bigger.With the release of the HumanMethylationEpic microarray by Illumina and datasets reaching into the thousands, analysis of these large datasets using popular R packages is becoming impractical due to memory requirements and even the time required to read the data from disk. As such there is an increasing need for computationally efficient methods to perform meaningful analysis on high dimension data.

The bigmelon R package provides a memory-efficient work-flow that enables users to perform complex, large scale analyses required in EWAS without huge RAM. Building on the CoreArray Genome Data Structure (.gds) file format and libraries packaged in 'gdsfmt', we provide a familiar wateRmelon-like work flow that facilitates reading-in, preprocessing, quality control and statistical analysis.

To demonstrate large-scale data analyses, we stored the entire contents of the marmal-aid database (>14,500 samples) in a .gds file and demonstrate quality measures and principal components analysis.

Overall, bigmelon provides a familiar environment for users to perform large-scale analyses where  convention methods would run out of memory. Bigmelon shows reasonable performance in speed compared to conventional methods.
