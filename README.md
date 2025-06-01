# An Expanded Subventricular Zone Supports Postnatal Cortical Interneuron Migration in Gyrencephalic Brains
Repository hosting code used for transcriptomics analysis used in "An Expanded Subventricular Zone Supports Postnatal Cortical Interneuron Migration in Gyrencephalic Brains".
Kim, JaeYeon<sup>1</sup>; Poddar, Aunoy<sup>1</sup>; Sandoval, Kadellyn; Chu, Julia; Horton, Emma; Cui, Di...; Paredes, Mercedes et al. (2025). _An Expanded Subventricular Zone Supports Postnatal Cortical Interneuron Migration in Gyrencephalic Brains._

# Description
In large, gyrencephalic brains, the stem cell niche known as the subventricular zone (SVZ) is expanded and arranged into previously unknown structures presumably to facilitate the expanded migration of inhibitory interneurons. To identify the cell types and relevant molecular mechanisms relevant for the development of inhibitory interneurons in the expanded SVZ, also known as the "Arc", we conducted single-nucleus RNA sequencing on GW30 and GW39 microdissected Arc regions from human brain samples. The code used to process and analyze this data, which is used to support the findings presented in Kim et al., are in this repository.

# Abstract
Cortical GABAergic interneurons generated in the ventral developing brain travel long distances to their final destinations. While there are examples of interneuron migration in the neonatal human brain, the extent of postnatal migration across species and how it contributes to cortical interneuron composition remains unknown. Here, we demonstrate that neonatal gyrencephalic brains, including humans and piglets, harbor an elaborate subventricular zone (SVZ), termed the Arc due to its curved morphology and expanded neuroblast populations. The Arc is absent in lissencephalic marmoset and mouse brains. Transcriptomic and histological approaches revealed that Arc neurons are diverse interneurons from the medial and caudal ganglionic eminences that migrate into the frontal, cingulate, and temporal cortex. Arc cortical targets in human and piglet brains exhibit an increase in VIP+ neuronal density compared to other regions. Our findings reveal that the Arc is a developmental structure that supports the expansion of postnatal neuronal migration for cortical interneuron patterning in gyrencephalic brains.


# Single-Nucleus RNA Sequencing of Human Arc
![Figure_2B_Arc_UMAP_Representation](https://github.com/aunoyp/arc/blob/main/readme_figs/Fig_2B_UMAP.png?raw=true)

## Preprocessing Files
FASTQ Files uploaded to GEO Browser (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255968). 
Post-Cellbender corrected count matrices filtered to remove cells expressing high percentages of mitochondrial and ribosomal RNA and run through the DoubletFinder 2.0.3 R package to remove doublets.
Preprocessing.Rmd loads matrices and filters for valid barcodes that pass filter thresholds saved in <THIS DIRECTORY>).

## CellRank
[Figure 2G](https://github.com/aunoyp/arc/blob/main/snRNA_seq/fig2_cellrank.ipynb)
[Figure 5B](https://github.com/aunoyp/arc/blob/final_version_edits/snRNA_seq/fig5_cellrank.ipynb)

## Monocle
[Monocle Graph Script](https://github.com/aunoyp/arc/blob/final_version_edits/snRNA_seq/fig2_graphtest.R)

## WGCNA
[WGCNA Analysis](https://github.com/aunoyp/arc/blob/final_version_edits/snRNA_seq/wgcna.Rmd)

# Spatial Transcriptomic Analysis of Piglet Arc Migratory Streams

## Loading Image Data
[Loading Image Data](https://github.com/aunoyp/arc/blob/final_version_edits/hiplex/st_all.Rmd)

## Plotting Fxns
[Helper Functions](https://github.com/aunoyp/arc/blob/final_version_edits/hiplex/st_functions.R)

## Nearest Neighbor Analysis
[Adapated from the histoCAT neighborhood analysis interaction score.](https://github.com/aunoyp/arc/blob/final_version_edits/hiplex/st_nn.Rmd)

