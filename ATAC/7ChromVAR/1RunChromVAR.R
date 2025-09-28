library(ArchR)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(RColorBrewer)  
library(TFBSTools)

setwd('~/workspace/Mida_collab/')
crfd_clean <- loadArchRProject('CorticalFolding_Clean/')
addArchRThreads(threads = 4)
pwm_path <- '/home/gaojie/workspace/Mida_collab/Motif_PWM' ###all PW matrix are downloaded from GimmeMotifs




all_motifs <- list.files(pwm_path)
lapply(all_motifs, function(motif_file_name){
    path <- paste0(pwm_path,'/',motif_file_name)
    line1 <- readLines(path)[1]
    motif_name <- substr(line1,2,nchar(line1))
    pfm_matrix <- read.csv(path,skip = 1,sep='\t',header=FALSE,col.names = c('A','C','G','T')) %>%
                      as.matrix() %>%
                      t()
    pwm_normalized = apply(pfm_matrix, 2, function(col) col / sum(col))
    print(colSums(pwm_normalized))
    pwm =  PWMatrix(
        ID = motif_name,
        name = motif_name,
        profileMatrix = pwm_normalized
    )
    return(pwm)
}) %>%
do.call(PWMatrixList,.) -> pwm_list
names(pwm_list) <- ID(pwm_list)
saveRDS(pwm_list,'ChromVar/pwm_list.rds')
pwm_list <- readRDS('ChromVar/pwm_list.rds')
length(pwm_list)


    
crfd <- subsetArchRProject(ArchRProj = crfd_clean, 
                                cells = rownames(data.frame(crfd@cellColData) %>% filter(celltype == 'vRG'),
                                outputDirectory = 'CorticalFolding_Clean_vRG/',
                                threads=4)
crfd <- addMotifAnnotations(
    ArchRProj = crfd,
    motifPWMs = pwm_list,
    name = 'GimmeMotif',
    force = TRUE
    )

crfd <- addBgdPeaks(crfd,force=TRUE)

crfd <- addDeviationsMatrix(
  ArchRProj = crfd, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(crfd,'CorticalFolding_Clean_vRG/')
                           
plotVarDev <- getVarDeviations(crfd, name = "MotifMatrix", plot = FALSE,n = 10)
colors <- colorRampPalette(rev(brewer.pal(name = 'RdYlBu', n = 5)))(250)
motif_matrix <- getMatrixFromProject(crfd,'MotifMatrix')
res_motif <- data.frame(row.names = rownames(motif_matrix@assays@data$deviations))


                           
for (sample in unique(crfd@cellColData$Sample)){
    cells <- crfd@cellColData %>% data.frame() %>% filter(Sample==sample) %>% rownames()
    tmp_mat <- motif_matrix@assays@data$deviations[,cells]
    res_motif[[sample]] <- rowSums(tmp_mat)/length(cells)
}
res_motif <- res_motif[,c('Sulcus','Adjacent','Distant')]
lapply(rownames(res_motif),function(motif){
    S_cells <- crfd@cellColData %>% data.frame() %>% filter(Sample=='Sulcus') %>% rownames()
    D_cells <- crfd@cellColData %>% data.frame() %>% filter(Sample=='Distant') %>% rownames()
    S_data <- motif_matrix@assays@data$deviations[motif,S_cells]
    D_data <- motif_matrix@assays@data$deviations[motif,D_cells]
    wilc <- wilcox.test(c(S_data,D_data),alternative = "two.sided") ###using two sided wilcox test to see if there is a significant difference
    return(wilc$p.value)
}) -> wilc_list
res_motif <- res_motif %>% mutate(p.value = wilc_list)

res_motif <- res_motif %>% mutate(p.adjust = p.adjust(res_motif$p.value,method = 'fdr'))

res_motif <- res_motif %>% filter(p.adjust <= 0.01) ###set the filtering threshold to 0.01
res_motif <- res_motif[((res_motif$Sulcus > res_motif$Adjacent) & (res_motif$Adjacent > res_motif$Distant) | (res_motif$Sulcus < res_motif$Adjacent) & (res_motif$Adjacent < res_motif$Distant)),] %>%  ###making sure they are gradient shifts 
                arrange(p.adjust)

res_motif %>% head()

head(res_motif[,c('Sulcus','Adjacent','Distant')])

GSEA_bi <- read.csv('CellOracle/output/GSEA_mat_bi_exp0.04pct0.05size20.csv',row.names = 1)
TF_dict <- read.csv('GimmeMotif_motif2tf.csv',row.names = 1)
selected_motifs <- c('GM.5.0.Sox.0013', ##LEF1
                     'GM.5.0.SMAD.0003', ##NFIC
                     'GM.5.0.SMAD.0004', ##NFIX
                     'GM.5.0.C2H2_ZF.0313', ##MAZ
                     'GM.5.0.Homeodomain.0145', ##MEIS1
                     'GM.5.0.Homeodomain_POU.0016' ##POU3F2
                    )

chrom_TFs <- unique(unlist(
      lapply(TF_dict[rownames(res_motif),'direct'], function(x) {
        cleaned <- gsub("\\[|\\]|\\s", "", x)
        no_quotes <- gsub("'", "", unlist(strsplit(cleaned, ",")))
        return(no_quotes)
      })
    ))

chrom_intersect <- chrom_TFs[chrom_TFs %in% rownames(GSEA_bi)]
res_motif$direct <- TF_dict[rownames(res_motif),'direct']
res_motif_vis <- res_motif[grepl(paste(paste0("'",c('MAZ','LEF1','NFIX','NFIC','MEIS1','POU3F2'),"'"),collapse = '|'),res_motif$direct),]
res_motif_vis %>% head()

TF_heat <- Heatmap(
  res_motif_vis[selected_motifs,c('Sulcus','Adjacent','Distant')] %>% t() %>% scale() %>% t(),
  name = "TFs_Zscore",
    col = colors,
        height = unit(6,'cm'),
    width = unit(3,'cm'),
    column_order = factor(c('Sulcus','Adjacent','Distant'),levels = c('Sulcus','Adjacent','Distant')),
    row_order = selected_motifs,
    cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8)
)

pdf("ChromVar/Selected_motifs_heatmap.pdf", width = 10, height = 8)
draw(TF_heat)
dev.off()