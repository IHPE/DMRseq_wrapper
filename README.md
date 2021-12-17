# DMRseq Wrapper

## Description
DMRseq ([github](https://github.com/kdkorthauer/dmrseq) or [bioconductor](https://bioconductor.org/packages/release/bioc/html/dmrseq.html)) 
performs differential methylation analysis: it aims at detecting
differentially methylated regions. It is adapted to organisms that have a
genome assembly at the scaffold level.
This script wraps DMRseq, prepares the data and parallelises the analysis.

In order to reproduce the results from the article you need to follow the following steps:

- Create a working folder and download the script in this repository into this folder.
- Inside the working folder, create a folder named 'data' that will hold the input data.
- Download the associated data of the article into the 'data/' folder.
- Launch the script:

```bash
Rscript script_dmrseq_nosim_SC_01.07.19.R
```

Or open the folder in rstudio as the working directory, open the script and source it or step through the code.
