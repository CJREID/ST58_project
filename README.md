# ST58 Manuscript
Scripts, data and outputs of a paper on the genomic epidemiology of *E. coli* ST58

## Overview
The following repository allows conscientious readers of the ST58 manuscript to reproduce all the data processing, statistics and figures presented in the paper.

It comprises two directories: __`scripts`__ and __`raw_data`__ (the contents of which should be self explanatory) and generates an __`outputs`__ folder with __`figures`__ and __`data`__ subdirectories.

## Installation
### Software requirements
These scripts are currently functional on mac OS Big Sur 11.5.2 using RStudio 1.4.1106 and R version 4.0.5. We cannot guarantee they will work on other distributions of R or RStudio. Your OS should not be an issue provided you use these versions of R and RStudio though.

### R Packages
You will need to install the following packages and versions to work with the scripts:
- data.table_1.14.0
- tidyverse_1.3.1
- magrittr_2.0.1
- RColorBrewer_1.1-2
- ggtree_3.1.0
- pheatmap_1.0.12
- reshape2_1.4.4
- ggpubr_0.4.0
- ggplot2_3.3.5
- tibble_3.1.4
- purrr_0.3.4
- readr_2.0.1
- stringr_1.4.0
- forcats_0.5.1
- tidyr_1.1.3
- dplyr_1.0.7

### Issues with ggtree
There have recently been some issues with ggtree in the way it interacts with dplyr. The solution is to install the latest version of ggtree directly from github instead of via BiocManager. You can do this in the console on RStudio with the __`remotes`__ package like so:
```
install.packages("remotes")
remotes::install_github("YuLab-SMU/ggtree")
```

## Usage
Clone this repository
```
git clone https://github.com/CJREID/ST58_project.git
cd ST58_project
pwd
```
Open the data_vis.R script in a text editor or RStudio and set the variable __`wrkdir`__ on line 18 to the output of __`pwd`__ above and save the script.

Run the data_vis.R script and watch your __`outputs`__ folder magically fill with goodies.

## Outputs
### Figures
1. Figure 1. Metadata summary
2. Figure 2. Core gene phylogenetic tree with metadata (Additional manual annotation in manuscript)
3. Figure 3. Source and F RST distributions by BAP cluster
4. Figure 4. pCERC4 heatmap (Exported in 3 parts; manually edited for publication)
5. Figure 5. Source distribution by ColV status in ST58 Collection
6. Figure 6. Tree-heatmap of presence/absence of BAP2 (ColV) positively- and negatively-associated genes (Additional manual annotation in manuscript)
7. Figure 7. Virulence and resistance gene carriage rates by BAP cluster
8. Figure 8. Source distributions of ColV in Enterobase genome collection
9. Figure 9. Heatmap of low SNP distances between epidemiologically unrelated sequences

### Supplementary
#### Tables
1. Supplementary Data 1. Processed metadata and accession numbers for all sequences
2. Supplementary Data 2. Gene presence/absence data
3. Supplementary Data 3. Genes postively/negatively-associated with BAP clusters as identified by Roary and Scoary
4. Supplementary Data 4. Enterobase metadata and ColV gene screening

#### Figures
1. Fig. S1. Source distribution across all collection years
2. Fig. S2. Serotype distributions by source and ColV+/-
3. Fig. S3. FimH distributions by source and ColV+/-
4. Fig. S4. F RST by source and ColV carriage
5. Fig. S5. Heatmap of Liu ColV criteria gene carriage for ColV+ sequences (Manually edited in manuscript)
6. Fig. S6. Absolute and relative ColV carriage by sources
7. Fig. S7. BAP6 Scoary Heatmap (Additional manual annotation in manuscript)
8. Fig. S8. Tree-heatmap of antimicrobial resistance genes
9. Fig. S9. Tree-heatmap of virulence-associated genes
10. Fig. S10. Tree-heatmap of plasmid replicon genes
11. Fig. S11. Tree-heatmap of pairwise SNP distance between all sequences, mapped to the core gene phylogeny

