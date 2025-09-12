suppressMessages(library(Seurat))
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))

### prepare urd root & tip
s_vRG <- subset(sob, anno2tp_fin == 'vRG')
s_vRG <- FindVariableFeatures(s_vRG, assay = 'SCT', nfeatures = 3000)
s_vRG <- RunPCA(s_vRG, assay = 'SCT', npcs = 50, features = VariableFeatures(s_vRG))
s_vRG <- RunHarmony(s_vRG, group.by.vars = 'orig.ident')
s_vRG <- RunUMAP(s_vRG, dims = 1:40, reduction = 'harmony', reduction.name = 'umap_harmony')
s_vRG <- FindNeighbors(s_vRG, dims = 1:40, reduction = 'harmony', graph.name = c('harmony_nn','harmony_snn'))

s_vRG <- FindClusters(s_vRG, graph.name = 'harmony_snn', resolution = 0.5)

sob@meta.data$Root_Tip <- NA
sob@meta.data$Root_Tip[(Cells(sob) %in% Cells(s_vRG)[s_vRG$seurat_clusters %in% c(4)]) & (sob@meta.data$timepoint == 'GW11')] <- "Root"
sob@meta.data$Root_Tip[which((sob@meta.data$anno2tp_fin %in% c("CFuPN","L6B","oRG")) & (sob@meta.data$timepoint == 'GW13'))] <- "Tips"

saveRDS(sob@meta.data, 'meta.urd_11-13.root_tip.rds')

### read in sob
sob <- readRDS('sob.mnn.11-13.rds')
sob@meta.data <- readRDS('sob_meta.anno_harmony_fin.0717.rds')

s <- subset(sob, anno2tp_fin %in% c('vRG','oRG','IPC_RG','IPC_Neuron','nascent ExN', 'Newborn CPN','CFuPN','L6B'))
s <- subset(s, timepoint %in% c('GW11','GW13'))

s@meta.data <- readRDS('meta.urd_11-13.root_tip.rds')[colnames(s),]
s <- subset(s, features = intersect(unique(unlist(readRDS('../hvglist_normdata_3000filtered.11-13.rds'))), rownames(s)))

### urd
axial <- createURD(count.data = s@assays$RNA@counts, meta = s@meta.data) #, min.cells=3, min.counts=3
axial <- calcDM(axial, knn = 226, sigma = 15, dcs.store = 50)

axial@group.ids$anno2tp_fin <- axial@meta$anno2tp_fin
axial@group.ids$Root_Tip <- as.character(meta[rownames(axial@group.ids),"Root_Tip"])
# pseudotime
root.cells <- cellsInCluster(axial, "Root_Tip", "Root")
axial.floods <- floodPseudotime(axial, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)
axial <- floodPseudotimeProcess(axial, axial.floods, floods.name="pseudotime")

axial@group.ids$Root_Tip <- as.character(meta[rownames(axial@group.ids),"Root_Tip"])
axial@group.ids$timepoint <- as.character(axial@meta[rownames(axial@group.ids),"timepoint"])
axial@group.ids$location <- as.character(axial@meta[rownames(axial@group.ids),"location"])
axial@group.ids$cellTypes <- as.character(axial@meta[rownames(axial@group.ids),"anno2tp_fin"])
root.cells <- cellsInCluster(axial, "Root_Tip", "Root")
#
s_rt <- readRDS('s_rt.rds')
tiplist <- c("1","2","3","4")
axial@group.ids[rownames(s_rt@meta.data)[which(s_rt@meta.data$tip.clusters %in% tiplist)], "tip.clusters"] <- s_rt@meta.data$tip.clusters[which(s_rt@meta.data$tip.clusters %in% tiplist)]

axial.ptlogistic <- pseudotimeDetermineLogistic(axial, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)
axial.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(axial, "pseudotime", logistic.params=axial.ptlogistic))

axial <- urdSubset(axial, cells.keep=rownames(axial.biased.tm))

axial.walks <- simulateRandomWalksFromTips(axial, tip.group.id="tip.clusters", root.cells=root.cells, 
                                           transition.matrix = axial.biased.tm, n.per.tip = 25000, 
                                           root.visits = 1, max.steps = 5000, verbose = F)
axial <- processRandomWalksFromTips(axial, axial.walks, verbose = F)

# annotation
tiplist <- c("1","2","3")
axial.tree <- loadTipCells(axial, "tip.clusters")
axial.tree <- buildTree(axial.tree, pseudotime = "pseudotime", tips.use=as.numeric(tiplist), divergence.method = "preference", cells.per.pseudotime.bin = 25, 
                        bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)
#axial.tree <- buildTree(object = axial.tree, pseudotime="pseudotime", divergence.method = "ks", tips.use=as.numeric(tiplist), weighted.fusion = T, use.only.original.tips = T, cells.per.pseudotime.bin=80, bins.per.pseudotime.window = 5, minimum.visits = 1, visit.threshold = 0.7, p.thresh = 0.025, save.breakpoint.plots = NULL, dendro.node.size = 100, min.cells.per.segment = 10, min.pseudotime.per.segment = .01, verbose = F)
axial.tree <- nameSegments(axial.tree, segments=tiplist, 
                           segment.names = c("oRG","CFuPN","L6B"), short.names = tiplist)
axial.tree@meta$cellTypes <- sob@meta.data[rownames(axial.tree@meta),'anno2tp_fin']


