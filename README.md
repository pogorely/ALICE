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
3. Install *data.table*, *stringdist*, *igraph* packages:
```R
install.packages("igraph")
install.packages("data.table")
install.packages("stringdist")
```
4. Download this github repository.

## Quick start
Lets load some data first and organize it to the list. 

This is one VJ combination (TRBV9-TRBJ2-7) from S1 donor from [link](https://www.biorxiv.org/content/early/2018/04/13/300343) on day 0 and day 15 after yellow fever immunization.

```R
library(data.table)
S1d15<-fread("sample/S1_d15_V9_J2_7.tsv")
S1d0<-fread("sample/S1_d0_V9_J2_7.tsv")
S1<-list(d0=S1d0,d15=S1d15)
```

Now lets run ALICE pipeline. Note, that algorithm will create a folder for files produced, if it does not exist:
```R
source("ALICE.R")
S1_alice<-ALICE_pipeline(DTlist=S1,folder="S1_res",cores=1,iter=10,nrec=5e5) #this takes few minutes to run
sapply(S1_alice,nrow)
```
For this VJ-combination we have no significant results in day 0 timepoint, and 34 significant hits in day 15 timepoint. 

Note, that for demo purposes we run it on 1 core with 10 iterations (0.5 mln sequences in each iteration) for generative probability estimation simulation. Total number of simulated TCR sequences (both inframe and out-of-frame) in this case is 5 million. 

In [paper](https://www.biorxiv.org/content/early/2018/07/23/375162) we used 100 mln simulated sequences for each VJ-class, and this takes a lot of time. 
Datasets from the paper are available [here](https://github.com/pogorely/ALICE_sample_data).

###Additional parameters
*Read_count_filter*(default 0) and *Read_count_neighbour*(default 1) parameters are two conceptually different count threshold for clones considered by the algorithm. 

Algorithm discards all clones with count  *Read_count_filter* or less prior to analysis, and it does not consider as neighbours clone with count *Read_count_neighbour* or less (but such clones are not discarded, so if they have a lot of high count neighbours, it could be significant hit).

*qL*(default is FALSE) uses different selection factor for different lengths instead of universal selection factor. 

*cor_method* (default "BH") specifies multiple testing correction method supplied to *p.adjust* function from *stats* package. *P_thres* (default 0.001) determines p-value threshold after correction.

## Input file format
Algorithm operates on R *data.table* with following mandatory columns: 

*CDR3.amino.acid.sequence*, *bestVGene*, *bestJGene*, *Read.count*. 

See sample datasets in *sample* folder.

There is a pipeline to import data in immunoSEQ format. It imports all files in given *folder* with counts larger than *Read_thres*, and also identifies CDR3nt and collapses clones with unique CDR3nt+V+J. immunoSEQ data with short reads do not contain complete CDR3nt sequence. For such data specify *trim* parameter (for reads of length 60 typically *trim=3*), in this case algorithm will identified CDR3nt trimmed from J side by certain number of amino acids.  

```R
source("import.R")
dtlist<-import_immunoseq_pipeline(folder="path",trim=3,Read_thres=0)
results<-ALICE_pipeline(dtlist)
```

## (Experimental) Using pipeline with OLGA for Pgen estimation
Install OLGA first (see [OLGA github](https://github.com/zsethna/OLGA) for details).

ALICE will call OLGA in background. Multi-core usage is not available on Windows.  

```R
source("ALICE.R")

#import data

S1d15<-fread("sample/S1_d15_V9_J2_7.tsv")
S1d0<-fread("sample/S1_d0_V9_J2_7.tsv")
S1<-list(d0=S1d0,d15=S1d15)

S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,folder="S1_res",cores=1)
sapply(S1_alice,nrow)

```