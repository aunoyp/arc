# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 10)

#interneurons = readRDS("/wynton/group/paredes/Aunoy/integrations_w_IPC/reviewer_fig2.rds")
#interneurons = readRDS("/wynton/group/paredes/Aunoy/integrations_w_IPC/reviewer_fig2_wgcna_subtype.rds")
#interneurons = readRDS("/wynton/group/paredes/Aunoy/integrations_w_IPC/cgewgcna.rds")
interneurons = readRDS("/wynton/group/paredes/Aunoy/integrations_w_IPC/cgewgcna.rds")

#interneurons = interneurons[, interneurons$lineage == "CGE_Cortical"]


#interneurons@misc$active_wgcna = NULL
#interneurons@misc$inter_fraction_0.05 = NULL

## Setup
#print("setting up WGCNA")
#interneurons <- SetupForWGCNA(
#  interneurons,
#  gene_select = "fraction", # the gene selection approach
#  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
#  group.by = "dataset",
#  wgcna_name = "interneuron_expression_modules" # the name of the hdWGCNA experiment
#)
#interneurons$lineage = "cge_lineage"

# construct metacells:

#print("creating metacells")
#interneurons <- MetacellsByGroups(
#  seurat_obj = interneurons,
#  group.by = c("lineage", "dataset"), # specify the columns in seurat_obj@meta.data to group by
#  reduction = 'integrated.scvi', # select the dimensionality reduction to perform KNN on
#  k = 25, # nearest-neighbors parameter
#  max_shared = 10, # maximum number of shared cells between two metacells
#  ident.group = 'lineage' # set the Idents of the metacell seurat object
#)

#interneurons = NormalizeMetacells(interneurons)

# set up a gene expression matrix list for each study
#interneurons <- SetMultiExpr(
#  interneurons,
#  group_name = "cge_lineage",
#  group.by = "lineage",
#  multi.group.by = "dataset",
#  multi_groups = NULL
#)
#print("searching for soft-power threshold")
## Select soft-power threshold
#interneurons <- TestSoftPowersConsensus(
#  interneurons,
#    group.by = 'lineage',
#group_name = 'cge_lineage',
#setDatExpr = FALSE
#)

print("constructing network...")
interneurons<- ConstructNetwork(
  interneurons,
  soft_power=c(6,6,6,5),
  consensus=TRUE,
  overwrite_tom = TRUE,
  tom_name = "TOM_CGE",
  setDatExpr=FALSE,
  detectCutHeight=0.999999
)

# plot the dendrogram
fig_dir = "~/"
pdf(paste0(fig_dir, "CGE_dendro.pdf"),height=3, width=6)
PlotDendrogram(interneurons, main='CGE Consensus Dendrogram')
dev.off()
print("computing eigengenes...")
interneurons <- ModuleEigengenes(
  interneurons,
  group.by.vars="dataset", # snRNAseq batch
  verbose=FALSE
)

print("computing module connectivity...")
# compute module connectivity:
interneurons <- ModuleConnectivity(interneurons)

## run RenameModules
interneurons <- ResetModuleNames(
  interneurons,
  new_name = "CGE"
)
print(names(GetModules(interneurons)))


# save data:
print("saving object")
saveRDS(interneurons, "/wynton/group/paredes/Aunoy/integrations_w_IPC/cgewgcna.rds")
