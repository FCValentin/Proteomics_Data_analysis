# __Nolwenn Proteomics Data__

Visual output of normalized data matrix

## Getting Started

This repository reports all scripts used for bio-informatics analysis.  

### Prerequisites and genome version
Genome version and annotation :
* Orthologs of Xenopus laevis 9.2 from [Xenbase](https://download.xenbase.org/xenbase/Genomics/JGI/Xenla9.2/)
* SampleAnnot file to describe samples
* SamplePerCond.tsv file that contains all samples to project for a comparison
  
### Tools used for analysis

* R libraries used are shown in each R files

### Scripts used for analysis

* [PCA](RScript/PCA.R) analysis
* [Heatmap](RScript/Heatmap.R) analysis
* All the analysis performed from samplesPerCond.tsv => [Analysis.R](RScript/Analysis.R) 
* Binarization of the datamatrix => [TableIniatialization.R](RScript/TableIniatialization.R)

## Authors

**Valentin FRANCOIS--CAMPION** 
