suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(batchelor))

### merge & SCT norm
t <- read.delim('../Filtered_cell_metadata_23rd_Jan_2024.csv',sep = ',')
t$location <- unlist(lapply(t[,1], function(item) strsplit(item, split = '[_]')[[1]][3]))
t$tag <- paste0(t$Week, '_',t$location)

taglist <- list.files('/data/snrna/')[1:6]
soblist <- lapply(taglist, function(tag) {
    print(tag)
    scrna <- Read10X(data.dir = paste0('/data/snrna/', tag, '/outs/filtered_feature_bc_matrix'))
    s <- CreateSeuratObject(count = scrna, project = tag, min.cells = 0, min.features = 0)
    # QC
    cell_list <- t[t$tag == tag, 'cell_ID']
    s <- subset(s, cells = cell_list)
    # rename
    s <- RenameCells(s, new.names = paste0(tag, '_', colnames(s)))
    s@meta.data$batch <- tag
    return(s)
})

genelist <- lapply(soblist, function(item) {
        item <- NormalizeData(item)
        item <- FindVariableFeatures(item, nfeatures = 3000)
        return(VariableFeatures(item))
})
saveRDS(genelist, 'hvglist_normdata_3000filtered.11-13.rds')

sob.sct <- soblist[[1]]
for (i in 2:length(soblist)) {
    sob.sct <- merge(sob.sct, y = soblist[[i]])
}

s <- sob.sct
# normalize1
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
s <- SCTransform(s, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt'), verbose = FALSE)
# normalize2
s <- CellCycleScoring(s, assay = 'SCT', s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
s <- SCTransform(s, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt','S.Score','G2M.Score'), verbose = FALSE)

saveRDS(s, 'sob.merge_sct.11-13.rds')
