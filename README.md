# asp_lncrna
Investigating the Transcriptional Response to Itraconazole in Aspergillus Fumigatus

**Validation of lncRNAs using SRA data**

We validated the presence of lncRNAs in other datasets by utilising sequencing reads from public experiments. Three experiments were used to validate this: a time-series experiment on _Aspergillus fumigatus_ after exposure to itraconazole (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA482512), expression analysis of _ppzA_ mutant, and wild-type, after iron-starvation (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA381768) and transcriptional response of mitochondrial mutant compared to wild-type (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA627614). 

SRA read accessions were downloaded and named accordingly (e.g. time, fe, mito).

The SRAToolkit was then used to download these reads. Corresponding accessions and HPC script is found within `SRA/accessions`

**Testing strandedness**

Strandendess was tested for by using Salmon to map reads, from each of the SRA data sets, onto an indexed version of the newly generated GTF. 



