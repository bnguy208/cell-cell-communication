rm(list=ls())

library(anndata)
library(tibble)
library(CellChat)
library(patchwork)
library(circlize)

#  Part I: Data input & processing and initialization of CellChat object
## Prepare required input data for CellChat analysis

# Read the data into R using anndata R package
ad <- read_h5ad("ad_merged_neuron_epithelial_data.h5ad")

# Access count data matrix
counts <- t(as.matrix(ad$X))

# Access meta data
meta <- ad$obs

# # Ensure barcodes match
meta <- meta[rownames(meta) %in% colnames(counts), ]
counts <- counts[, colnames(counts) %in% rownames(meta)]

# Save the data for CellChat analysis
save(counts, meta, file = "ad_merged_data.RData")

## Create a CellChat object
cellchat <- createCellChat(object = counts, meta = meta, group.by = "celltype")

## Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# Use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Set the used database in the object
cellchat@DB <- CellChatDB.use

## Pre-processing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 2 * 1024^3)  # Increase the maximum allowed size for globals (2 GB)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network
## Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(5,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Part III: Visualization of cell-cell communication network
## Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("BMP")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[3,] # show one ligand-receptor pair

# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

save.image('cell_chat_environment.RData')

# Tuft and NMUens subset
cellchat_subset <- subsetCellChat(cellchat, idents.use = c("Tuft", "ENC6"))

# Check communication probabilities between Tuft and ENC6
subsetCommunication(cellchat, sources.use = "Tuft", targets.use = "ENC6")
subsetCommunication(cellchat, sources.use = "ENC6", targets.use = "Tuft")

# Visualize interactions between Tuft and ENC6
netVisual_circle(cellchat_subset, 
                 sources.use = "Tuft", 
                 targets.use = "ENC6", 
                 title.name = "Interactions between Tuft and ENC6")
