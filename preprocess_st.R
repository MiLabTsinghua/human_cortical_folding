library(Seurat)

tag = 'GW11'
s <- Load10X_Spatial(paste0('/home/baiyq/projects/fold_dev/data/10X_Visium/', tag, '/outs'))

# normalize1
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
s <- subset(s, subset = (nFeature_Spatial > 800) & (percent.mt < 5))
s <- SCTransform(s, assay = 'Spatial', new.assay.name = 'SCT', vars.to.regress = c('percent.mt'), verbose = FALSE)
# normalize2
s <- CellCycleScoring(s, assay = 'SCT', s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
s <- SCTransform(s, assay = 'Spatial', new.assay.name = 'SCT', vars.to.regress = c('percent.mt','S.Score','G2M.Score'), verbose = FALSE)
