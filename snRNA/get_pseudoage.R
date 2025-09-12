library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(irlba)
library(reshape2)
theme_set(theme_cowplot())

enableWGCNAThreads(nThreads = 8)
set.seed(1234)

combined = readRDS("/Reference Data/integrated.v4.rds")
Folding_dat = readRDS("../sob.pseudoage.objlist.rds")

### WGCNA
combined.RG = subset(combined, Subclass == "Radial glia")
combined.RG = subset(combined.RG, Estimated_postconceptional_age_in_days != 54 & Estimated_postconceptional_age_in_days != 60  & Estimated_postconceptional_age_in_days != 147)
DefaultAssay(combined.RG) = "RNA"
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
  group_name = sort(as.numeric(unique(combined.RG@meta.data$metacell_grouping))),
  group.by='Estimated_postconceptional_age_in_days',
  assay = 'RNA',
  slot = 'data'
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
combined.RG = ScaleData(combined.RG, assay = "RNA", features = VariableFeatures(combined.RG))
for(i in names(Modulenames)){
    npcs = min(50, length(Modulenames[[i]]) - 1)
    mepca.list[[i]] = irlba(A = t(x = combined.RG@assays$RNA@scale.data[Modulenames[[i]],]), nv = npcs)
}
MEs_reference = data.frame(matrix(ncol = 0, nrow = length(colnames(combined.RG))))
for(i in names(mepca.list)){
    MEs_reference[[i]] = mepca.list[[i]]$u[, 1] * mepca.list[[i]]$d[1]
}
rownames(MEs_reference) = colnames(combined.RG)

### apply
hvg <- VariableFeatures(combined.RG)
for(i in 1:length(Folding_dat)){
    Folding_dat[[i]]@meta.data <- readRDS('../sob_meta.anno_harmony_fin.0722.rds')[colnames(Folding_dat[[i]]),]
    Folding_dat[[i]]$location_ori = NULL
    Folding_dat[[i]]$Type = Folding_dat[[i]]$anno2tp_fin 
    Folding_dat[[i]]$Estimated_postconceptional_age_in_days = ifelse(Folding_dat[[i]]$timepoint == "GW11", 63, 
                                                                     ifelse(Folding_dat[[i]]$timepoint == "GW13", 77, 105))
    Folding_dat[[i]]$timepoint = NULL
    Folding_dat[[i]]$Region = Folding_dat[[i]]$location
    Folding_dat[[i]]$location = NULL
    Folding_dat[[i]]$anno2tp_fin = NULL
    Folding_dat[[i]]$Study = "Folding"
    Folding_dat[[i]] = subset(Folding_dat[[i]], Type == "vRG")
}

#remove certain batches from ref data, set up hdWGCNA assay
combined.RG = subset(combined, Subclass == "Radial glia")
DefaultAssay(combined.RG) = "RNA"

#merge folding data with ref data and scale them together
combined.RG.Folding = merge(Folding_dat[[1]], Folding_dat[-1])
DefaultAssay(combined.RG.Folding) = "RNA"
combined.RG.Folding[["SCT"]] = NULL
combined.RG[["SCT"]] = NULL
combined.RG.Folding$Study = "folding"
combined.RG$Study = "ref"
combined.RG$batch = combined.RG$Estimated_postconceptional_age_in_days
combined.RG.Folding = merge(combined.RG.Folding, combined.RG)
combined.RG.Folding = NormalizeData(combined.RG.Folding)
combined.RG.Folding = ScaleData(combined.RG.Folding, features = hvg)

#project the merged data to the module feature loadings
MEs_folding = data.frame(matrix(ncol = 0, nrow = ncol(combined.RG.Folding)))
for(i in names(Modulenames)){
    feature_loading = mepca.list[[i]]$v[, 1]
    MEs_folding[[i]] = t(as.matrix(t(feature_loading)) %*% combined.RG.Folding@assays$RNA@scale.data[Modulenames[[i]], ])
}
rownames(MEs_folding) = colnames(combined.RG.Folding)

#perform PCA for the module eigengenes of the reference data
pcaresults_ref = prcomp(MEs_reference)
pseudoage_ref = as.data.frame(pcaresults_ref$x)
pseudoage_ref$Estimated_postconceptional_age_in_days = combined$Estimated_postconceptional_age_in_days[rownames(MEs_reference)]

pseudoage_ref$Estimated_postconceptional_age_in_days = factor(pseudoage_ref$Estimated_postconceptional_age_in_days)

#project module eigengenes of the folding data to the pseudoage axis
pseudoage_folding = data.frame(matrix(ncol = 0, nrow = ncol(combined.RG.Folding)))
pseudoage_folding[["PC1"]] = t(t(as.matrix(pcaresults_ref$rotation[, 1])) %*% t(MEs_folding))[,1]
rownames(pseudoage_folding) = colnames(combined.RG.Folding)

pseudoage_folding$batch = combined.RG.Folding$batch
pseudoage_folding$Study = combined.RG.Folding$Study
pseudoage_folding$PC1 = -pseudoage_folding$PC1 # pseudoage
pseudoage_folding$batch <- factor(pseudoage_folding$batch,
    levels = c(sort(unique(as.numeric(combined.RG$Estimated_postconceptional_age_in_days))), 'GW11_C', 'GW11_G', 'GW11_D', 'GW13_C', 'GW13_G', 'GW13_D'), ordered = T)

