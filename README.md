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
1. **Gene co-expression network.** Use *GeneCoexpressionNet.R* to generate the gene co-expression network;
2. **Gene functional interaction network.** Use *CorrectNet.R* to generate the gene functional interaction network;
3. **Gene-site heterogeneous network.** Use *GeneSiteNet.R* to obtain the gene-site heterogeneous network;
4. **Transition matrix of the gene-site heterogeneous network.** Use *Normalization.R* to generate the transition matrix of the gene-site heterogeneous network;
5. **Functional m6A-mediated genes.** Use *m6ARWRH.R* to perform RWRH for each of the m6A regulator and finally gain potential functional m6A-mediated gene.

Then, the mutation frequency of these functional m6A-mediated genes are input as heat vector to perform HotNet2 in three PPI networks. The code and network data used by HotNet2 can be downloaded from https://github.com/raphael-group/hotnet2.
