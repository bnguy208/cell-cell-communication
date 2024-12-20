library(anndata)
library(tibble)
library(CellChat)
library(patchwork)
library(circlize)

#####################################
# PART 1: Prepare input data matrix 
#####################################
# Read the data into R using anndata R package
ad <- read_h5ad("ad_mn.h5ad")

# Access count data matrix
counts <- t(as.matrix(ad$X))

# Access meta data
meta <- ad$obs

# Ensure barcodes match
meta <- meta[rownames(meta) %in% colnames(counts), ]
counts <- counts[, colnames(counts) %in% rownames(meta)]

# Save the data for CellChat analysis
save(counts, meta, file = "ad_mn_data.RData")

#####################################
# PART 2: Create CellChat object
#####################################
# Create CellChat object                           
cellchat <- createCellChat(object = counts, meta = meta, group.by = "celltype")

#####################################
# PART 3: Set & update CellChatDB
#####################################
# Load the mouse database
CellChatDB <- CellChatDB.mouse

# Show ligand-receptor categories
showDatabaseCategory(CellChatDB.mouse)

# Remove problematic interactions
which(CellChatDB[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB[["interaction"]] <- CellChatDB[["interaction"]][-1887,]

which(CellChatDB[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB[["interaction"]] <- CellChatDB[["interaction"]][-1900,]

# Set all CellChatDB for cell-cell communication analysis 
cellchat@DB <- CellChatDB

#####################################
# PART 4: Preprocess Data
#####################################
cellchat <- subsetData(cellchat) # Subset the data based on the database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#####################################
# PART 5: Run CellChat Pipeline
#####################################
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10) # Adjust 'min.cells' as needed
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#####################################
# PART 6: Visualize the results
#####################################
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xp = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F, title.name = "Number of Interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
