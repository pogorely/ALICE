# ALICE
Detecting TCR involved in immune responses from single RepSeq datasets.

## About
Here we provide an R implementation of ALICE approach, see [paper](https://www.biorxiv.org/content/early/2018/07/23/375162) for details.
## Software requirements
Any OS where R is available (Linux, OS X, Windows), however parallel computing is currently not available on Windows.  

## Installation

1. Install R distribution of choice (i.e. from [R Core team](https://cloud.r-project.org/), or [Microsoft R Open](https://mran.microsoft.com/open/) )
2. Install BioStrings package from bioconductor: open R console and execute following commands: 
```R
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```
3. Install data.table, stringdistr, igraph packages:
```R
install.packages("igraph")
install.packages("data.table")
install.packages("stringdistr")
```

## Quick start
Lets load some data first and organize it to the list. 

This is one VJ combination (TRBV9-TRBJ2-7) from S1 donor from link on day 0 and day 15 after yellow fever immunization.

```R
S1d15<-fread("sample/S1_d15_V9_J2_7.tsv")
S1d0<-fread("sample/S1_d0_V9_J2_7.tsv")
S1<-list(d0=S1d0,d15=S1d15)
```

Now lets run ALICE pipeline. Note, that algorithm will create a folder for files produced, if it does not exist:
```R
source("ALICE.R")
S1_alice<-ALICE_pipeline(S1,folder="S1_res",cores=1,iter=10,nrec=5e5) 
sapply(S1_alice,nrow)
```
For this VJ-combination we have no significant results in day 0 timepoint, and 34 significant hits in day 15 timepoint. 

Note, that for demo purposes we run it on 1 core with 10 iterations (0.5 mln sequences in each iteration) for generative probability estimation simulation. Total number of simulated TCR sequences (both inframe and out-of-frame) in this case is 5 million. 
In [paper](https://www.biorxiv.org/content/early/2018/07/23/375162) we used 200 mln simulated sequences for each VJ-class, and this takes a lot of time. 

## Input file format
Algorithm operates on R dataset with following mandatory columns: 

*CDR3.amino.acid.sequence*, *bestVGene*, *bestJGene*, *Read.count*. 

See sample datasets in *sample* folder.
