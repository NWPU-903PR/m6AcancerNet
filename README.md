## m6AcancerNet
A novel network-based approach for identifying m6A-mediated driver genes which are defined as genes mediated by m6A methylation, significantly mutated and functionally interacted in cancer.
## Software requirments
All the code are based on R language, and R version ≥ 3.6.3  is suggested.
## Data Acquisition
* RNA-seq and mutation data are downloaded from TCGA by using "TCGAbiolinks" R package.
* Gene annotation data are downloaded from GENCODE (https://www.gencodegenes.org/). We recommend using **GENCODE v22**, because it's consistent with the reference files used by GDC(https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files).
* m6A sites co-methylation data are gained from the article **m6Acomet**(https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2840-3).
* The reference gene functional network data are gained from **GeneMANIA**.
## Usage Example
### Take AML dataset as an example
* 
