# bigmelon
Memory Efficient Methods for DNA Methylation Analysis

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

## Note
An early version of this was reviewed for BioC in 2014 (issue 1044).  The reviewer made the reasonable point that it was essentially just a stub at that point. Now it is much better developed so we have submitted it again.  Because of the time 
passed and the fact the submission process has changed this is a new submission.

