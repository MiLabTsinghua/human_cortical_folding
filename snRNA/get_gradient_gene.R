suppressMessages(library(Seurat))

tag <- 'vRG'
ctlist <- list('vRG' = c('vRG'), 'IPC' = c('IPC1','IPC2'), 'EX' = c("Nascent ExN","CFuPN","L6B"))
s <- subset(sob, anno2tp_fin %in% ctlist[[tag]])
s <- subset(s, location %in% c('Sulcus','Distant'))
s <- subset(s, timepoint == 'GW11')
s <- SetIdent(s, value = s@meta.data$location)
s <- PrepSCTFindMarkers(s, assay = "SCT", verbose = TRUE)
ref.markers <- FindAllMarkers(object = s, only.pos = TRUE, min.pct = 0.25, assay = "SCT", logfc.threshold = log2(1.2))

s <- subset(sob, anno2tp_fin %in% ctlist[[tag]])  
s <- subset(s, genes = ref.markers$gene)
mat_mean <- as.data.frame(t(apply(s@assays$SCT@data, 1, function(entry) tapply(entry, as.character(s$location), mean))))
mat_mean$gene <- rownames(mat_mean)
rownames(mat_mean) <- NULL
head(mat_mean)

mat_snrna <- merge(ref.markers, mat_mean, by = 'gene')
mat_snrna <- mat_snrna[mat_snrna$p_val_adj < 0.05,]
flaglist <- apply(mat_snrna, 1, function(entry) {
    if (as.character(entry['cluster']) == 'Distant') {
        flag <- ifelse((entry['Distant'] > entry['Adjacent']) & (entry['Adjacent'] > entry['Sulcus']), 'Pass', 'Fail')
    }
    if (as.character(entry['cluster']) == 'Sulcus') {
        flag <- ifelse((entry['Distant'] < entry['Adjacent']) & (entry['Adjacent'] < entry['Sulcus']), 'Pass', 'Fail')
    }
    return(flag)
})
mat_snrna$flag <- unlist(flaglist)

