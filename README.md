# PhosMap
A Comprehensive R Package For Analyzing Quantitative Phosphoproteomics Data.
- PhosMap is a comprehensive R package for analyzing quantitative phosphoproteomics data. Modules in PhosMap were classified into two major categories: (1) data pre-processing and (2) data analysis and presentation. 
    - An intact data pre-processing procedure of phosphoproteomics data covered three main steps: merging input files after quality control, mapping phosphorylation sites (p-sites) to the corresponding protein sequence and data normalization. 
    - PhosMap incorporated four analysis modules, including clustering and differential expression analysis, time course analysis, kinase-substrate enrichment analysis to find activated/deactivated kinases and motif enrichment analysis.

* In PhosMap, "ksea" and "rmotif" packages were imported for kinase-substrate enrichment analysis to find activated/deactivated kinases and motif enrichment analysis, respectively. The two packages are not denpendent on Bioconductor or CRAN, so users need to install them by the following methods prior to installing PhosMap.
  - (1) devtools::install_github(https://github.com/omarwagih/rmotifx/') and devtools::install_github(https://github.com/evocellnet/ksea/')
  - (2) Downloading souce code package of 'ksea' and 'rmotif' for local installation via above-mentioned URLs
