library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 8)
set.seed(1234)
library(irlba)
library(reshape2)

#remove certain batches from ref data, set up hdWGCNA assay
combined.RG = subset(combined, Subclass == "Radial glia")
combined.RG = subset(combined.RG, Estimated_postconceptional_age_in_days != 54 & Estimated_postconceptional_age_in_days != 60  & Estimated_postconceptional_age_in_days != 147)
DefaultAssay(combined.RG) = "RNA"
#combined.RG = JoinLayers(combined.RG, assay = "RNA")
combined.RG = NormalizeData(combined.RG, assay = "RNA")
combined.RG = FindVariableFeatures(combined.RG, assay = "RNA", nfeatures = 3000)
combined.RG = SetupForWGCNA(
    combined.RG,
    features = VariableFeatures(combined.RG),
    wgcna_name = "wgcna"
)

combined.RG = ScaleData(combined.RG, assay = "RNA", features = VariableFeatures(combined.RG), vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
combined.RG = RunPCA(combined.RG, assay = "RNA", verbose = FALSE)
combined.RG = MetacellsByGroups(
    seurat_obj = combined.RG,
    group.by = "Estimated_postconceptional_age_in_days",
    reduction = "pca",
    k = 20,
    max_shared = 10,
    ident.group = "Estimated_postconceptional_age_in_days",
    assay = "RNA"
)

combined.RG = NormalizeMetacells(combined.RG)
combined.RG <- SetDatExpr(
  combined.RG,
  group_name = sort(as.numeric(unique(combined.RG@meta.data$metacell_grouping))),#c(54, 63, 77, 91, 98, 112, 126), #kept "Estimated_postconceptional_age_in_days"
  group.by='Estimated_postconceptional_age_in_days', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

#calculate softpower
combined.RG = TestSoftPowers(combined.RG, networkType = "signed")
wrap_plots(PlotSoftPowers(combined.RG), ncol = 2)

#construct gene network and calculate modules, minModuleSize = 5 & mergeCutHeight = 0 to make the module small enough.
combined.RG <- ConstructNetwork(
  combined.RG,
  tom_name = 'AGE_RG',
  minModuleSize = 5,
  mergeCutHeight = 0,
  overwrite_tom = TRUE
)
PlotDendrogram(combined.RG, main='RG Modules')

combined.RG = ModuleEigengenes(combined.RG, assay = "RNA")

#Select modules whose eigengene value highly correlated with age.
MEs = GetMEs(combined.RG, harmonized=FALSE)
MEcorAge = as.data.frame(cor(MEs, combined.RG$Estimated_postconceptional_age_in_days))
MEcorAge = subset(MEcorAge, V1 >= 0.6 | V1<= -0.6)
MEs = MEs[, rownames(MEcorAge)]

#perform PCA for the module eigengenes
pcaresults = prcomp(MEs)
PCA = as.data.frame(pcaresults$x)
PCA$Estimated_postconceptional_age_in_days = combined.RG$Estimated_postconceptional_age_in_days
PCA$Estimated_postconceptional_age_in_days = factor(PCA$Estimated_postconceptional_age_in_days)

#extract modules
Modules = GetModules(combined.RG)
Modulelist = split(Modules, f = Modules$color)
Modulenames = list()
for(i in names(Modulelist)){
    Modulenames[[i]] = rownames(Modulelist[[i]])
}
Modulenames = Modulenames[rownames(MEcorAge)]

#re-calculate module eigengenes (this time we can get feature loadings)
mepca.list = list()
combined.RG = ScaleData(combined.RG, assay = "RNA", features = VariableFeatures(combined.RG)) #scale the data once again to get unregressed matrix
for(i in names(Modulenames)){
    npcs = min(50, length(Modulenames[[i]]) - 1)
    mepca.list[[i]] = irlba(A = t(x = combined.RG@assays$RNA@scale.data[Modulenames[[i]],]), nv = npcs) #PCA
}
MEs_reference = data.frame(matrix(ncol = 0, nrow = length(colnames(combined.RG))))
for(i in names(mepca.list)){
    MEs_reference[[i]] = mepca.list[[i]]$u[, 1] * mepca.list[[i]]$d[1]
}
rownames(MEs_reference) = colnames(combined.RG)

hvg <- VariableFeatures(combined.RG)
Folding_dat$Estimated_postconceptional_age_in_days = 63
Folding_dat$Study = "Folding"
Folding_dat$batch = "spatial"
Folding_dat = subset(Folding_dat, layer_anno == "VZ")
#remove certain batches from ref data, set up hdWGCNA assay
combined.RG = subset(combined, Subclass == "Radial glia")
DefaultAssay(combined.RG) = "RNA"
DefaultAssay(Folding_dat) = "Spatial"
#merge folding data with ref data and scale them together
combined.RG.Folding = Folding_dat
DefaultAssay(combined.RG.Folding) = "Spatial"
combined.RG.Folding[["SCT"]] = NULL
combined.RG[["SCT"]] = NULL
combined.RG.Folding$Study = "folding"
combined.RG$Study = "ref"
combined.RG$batch = combined.RG$Estimated_postconceptional_age_in_days
combined.RG.Folding = merge(combined.RG.Folding, combined.RG)
combined.RG.Folding = JoinLayers(combined.RG.Folding)
combined.RG.Folding = NormalizeData(combined.RG.Folding)
combined.RG.Folding = ScaleData(combined.RG.Folding, features = hvg)

mat = GetAssayData(combined.RG.Folding,assay = 'RNA', slot = 'scale.data')
#project the merged data to the module feature loadings
MEs_folding = data.frame(matrix(ncol = 0, nrow = ncol(combined.RG.Folding)))
for(i in names(Modulenames)){
    genelist <- intersect(Modulenames[[i]], rownames(mat))
    feature_loading = mepca.list[[i]]$v[, 1][Modulenames[[i]] %in% rownames(mat)]
    MEs_folding[[i]] = t(as.matrix(t(feature_loading)) %*% mat[Modulenames[[i]][Modulenames[[i]] %in% rownames(mat)], ])
}
rownames(MEs_folding) = colnames(combined.RG.Folding)

pcaresults_ref = prcomp(MEs_reference)
pseudoage_ref = as.data.frame(pcaresults_ref$x)
pseudoage_ref$Estimated_postconceptional_age_in_days = combined$Estimated_postconceptional_age_in_days[rownames(MEs_reference)]
pseudoage_ref$Estimated_postconceptional_age_in_days = factor(pseudoage_ref$Estimated_postconceptional_age_in_days)

pseudoage_folding = data.frame(matrix(ncol = 0, nrow = nrow(combined.RG.Folding@meta.data)))
pseudoage_folding[["PC1"]] = t(t(as.matrix(pcaresults_ref$rotation[, 1])) %*% t(MEs_folding))[,1]
rownames(pseudoage_folding) = colnames(combined.RG.Folding)

