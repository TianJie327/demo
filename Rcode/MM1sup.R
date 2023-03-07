########################################################################################

                                   #1、celltype

########################################################################################


library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggraph)
library(limma)
library(magrittr)
library(tibble)

setwd("~/zhuangle_file/MM1sup")
MM_1data <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
MM_1 <- CreateSeuratObject(counts = MM_1data, min.cells = 3, min.features = 400)
save(MM_1, file = "1_MM1.Rda")


table(grepl("^MT-", row.names(MM_1)))
MM_1 <- PercentageFeatureSet(MM_1, pattern = "^MT-", col.name = "percent.mt")

MM_1 <- subset(MM_1, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20) 
MM_1 <- CellCycleScoring(object = MM_1, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
MM_1 <- SCTransform(MM_1, vars.to.regress = c("percent.mt"), verbose = FALSE)
MM_1 <- RunPCA(MM_1, verbose = FALSE)
MM_1 <- RunUMAP(MM_1, dims = 1:50, verbose = FALSE)
MM_1 <- FindNeighbors(MM_1, dims = 1:50, verbose = FALSE)
MM_1 <- FindClusters(MM_1, verbose = FALSE, resolution = 0.1)

P1 <- DimPlot(MM_1, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 10) + NoLegend()
pdf("P_umap.pdf", width = 10)
P1
dev.off()
save(MM_1, file = "2_MM_1_SCT.Rda")


DefaultAssay(object = MM_1) <- "SCT"

#Melanoma 
MM_1gene <- c("MLANA", "MITF", "TYR","PMEL")
pdf("P_Melanoma_FeaturePlot.pdf", width = 15, height = 8)
FeaturePlot(MM_1, reduction = "umap", features = MM_1gene, cols = c("gray", "red"), label = TRUE)
dev.off()
pdf("P_Melanoma_VlnPlot.pdf")
VlnPlot(MM_1, features = MM_1gene, log=F, ncol = 2)
dev.off()

#Epidermal appendage keratinocytes 
MM_1gene <- c("KRT18", "KRT8")
pdf("P_Epidermal appendage keratinocytes_FeaturePlot.pdf", width = 10, height = 4)
FeaturePlot(MM_1, reduction = "umap", features = MM_1gene, cols = c("gray", "red"), label = TRUE)
dev.off()
pdf("P_Epidermal appendage keratinocytes_VlnPlot.pdf")
VlnPlot(MM_1, features = MM_1gene, log=F)
dev.off()

#Epidermal keratinocytes
MM_1gene <- c("KRT5", "KRT14")
pdf("P_Epidermal keratinocytes.pdf", width = 10, height = 1)
FeaturePlot(MM_1, reduction = "umap", features = MM_1gene, cols = c("gray", "red"), label = TRUE)
dev.off()
pdf("P_Epidermal keratinocytes_VlnPlot.pdf")
VlnPlot(MM_1, features = MM_1gene, log=F)
dev.off()

#Endothelial
MM_1gene <- c("PECAM1", "VWF")
pdf("P_Endothelial.pdf", width = 10, height = 4)
FeaturePlot(MM_1, reduction = "umap", features = MM_1gene, cols = c("gray", "red"), label = TRUE)
dev.off()
pdf("P_Endothelial_VlnPlot.pdf")
VlnPlot(MM_1, features = MM_1gene, log=F)
dev.off()

#Cancer-associated fibroblasts
MM_1gene <- c("FAP", "MME")
pdf("P_Cancer-associated fibroblasts.pdf", width = 10, height = 4)
FeaturePlot(MM_1, reduction = "umap", features = MM_1gene, cols = c("gray", "red"), label = TRUE)
dev.off()
pdf("P_Cancer-associated fibroblasts_VlnPlot.pdf")
VlnPlot(MM_1, features = MM_1gene, log=F)
dev.off()

#Fibroblasts
MM_1gene <- c("ACTA2", "DCN", "LUM")
pdf("P_Fibroblasts.pdf", width = 10, height = 8)
FeaturePlot(MM_1, reduction = "umap", features = MM_1gene, cols = c("gray", "red"), label = TRUE)
dev.off()
pdf("P_Fibroblasts_VlnPlot.pdf")
VlnPlot(MM_1, features = MM_1gene, log=F)
dev.off()

#Immune cells
MM_1gene <- c("PTPRC")
pdf("P_Immune cells.pdf", width = 10, height = 8)
FeaturePlot(MM_1, reduction = "umap", features = MM_1gene, cols = c("gray", "red"), label = TRUE)
dev.off()
pdf("P_Immune cells_VlnPlot.pdf")
VlnPlot(MM_1, features = MM_1gene, log=F)
dev.off()



new.cluster.ids <- c("Melanoma_1",                        #0
                     "Melanoma_2",                        #1
                     "Endothelial" ,                      #2
                     "Fibroblasts",                       #3
                     "Melanoma_3",                        #4
                     "Epidermal_appendage_keratinocytes", #5
                     "Immune_cells",                      #6
                     "Epidermal_keratinocytes",           #7
                     "Cancer_associated_fibroblasts"      #8     
)                     
                     
)
names(new.cluster.ids) <- levels(MM_1)
MM_1 <- RenameIdents(MM_1, new.cluster.ids)

pdf(file = "Figure4.pdf", width = 11)
DimPlot(MM_1, reduction = "umap", label = TRUE, pt.size = 0.5 , label.size = 6, label.box = T, repel = T, label.color = "white") + NoLegend()
dev.off()

save(MM_1, file = "MM_1_celltype.Rda")

########################################################################################

                                    #2、monocle2

########################################################################################

library(monocle)


MM_1$cellType <- Idents(MM_1)
Cells.sub <- MM_1[,MM_1$seurat_clusters %in% c(0, 1 ,4)]
table(Cells.sub$seurat_clusters)
MM_1sub <- Cells.sub
MM_1sub$cellType <- factor(MM_1sub$cellType, levels = c("Melanoma_1", "Melanoma_2", "Melanoma_3"))
table(MM_1sub$cellType)
data <- as(as.matrix(MM_1sub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = MM_1sub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        negbinomial.size()
)

mycds <- detectGenes(mycds, min_expr = 0.1) 
print(head(fData(mycds)))
expressed_genes <- row.names(subset(fData(mycds), num_cells_expressed >= 10)) 
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds)
save(mycds, file = "1_mycds.Rda") 
diff_test_res <- differentialGeneTest(mycds[expressed_genes,],
                                      fullModelFormulaStr = "~cellType")
save(diff_test_res, file = "2_diff_test_res_fullModelFormulaStr_cellType.Rda")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
mycds <- setOrderingFilter(mycds, ordering_genes)
plot_ordering_genes(mycds)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
save(mycds, file = "3_af_order_mycds.Rda") 

load("~/3_af_order_mycds.Rda")
p1 <- plot_cell_trajectory(mycds, color_by = "State", cell_size = 0.1, cell_name_size = 1, state_number_size = 1)
p2 <- plot_cell_trajectory(mycds, color_by = "cellType", cell_size = 0.1, cell_name_size =1, state_number_size = 1)
p3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime", cell_size = 0.1, cell_name_size = 1, state_number_size = 5)
P1 <- p1|p2|p3

pdf("Fig6.pdf", width = 15, height = 5)
P1
dev.off()

BEAM_res <- BEAM(mycds, branch_point = 1)
save(BEAM_res, file = "BEAM_res.Rda")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
save(mycds, file = "P8_mycds.Rda")


load("~/BEAM_res.Rda")
load("~/P8_mycds.Rda")

MMg <- c("MITF", "MLANA", "TYR", "PMEL",
         "RAB4A", "RAB5A", "RAB7A", "RAB11A",
         "LAMP1", "LAMP2", "M6PR")

Pmmg <- plot_genes_branched_heatmap(mycds[MMg,],
                                    branch_point = 1,
                                    num_clusters = 2,
                                    branch_labels = c("Non-Pigmentation", "Pigmentation"),
                                    branch_colors = c(
                                      "#b80c09", #Melanoma_3
                                      "#d3d3d3", #Melanoma_2
                                      "#000000"  #Melanoma_1
                                    ),
                                    use_gene_short_name = T,
                                    show_rownames = T,
                                    return_heatmap = T,
                                    cores = 4
)
ggsave("fig7A.pdf", Pmmg$ph_res, width = 6.5, height = 10)


########################################################################################

                                    #3、GOBP

########################################################################################


library(AnnotationHub)
library(BiocGenerics)
library(parallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(readr)
library(enrichplot)
library(DO.db)
library(KEGGREST)
library(dplyr)
library(ggnewscale)
library(ggupset)
library(ggridges)
library(patchwork)
library(tidyverse)

MMCluster <- MM_1[,MM_1$seurat_clusters %in% c("0",
                                               "1",
                                               "4")]

cluster.markers <- FindAllMarkers(MMCluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene1<- filter(cluster.markers, cluster == "Melanoma_1")
fivenum(gene1$avg_log2FC)
gene <- filter(cluster.markers, cluster == "Melanoma_1" & avg_log2FC >= 0.35 & p_val_adj < 0.05)
SYMBOL <- rownames(gene)
ENID <- bitr(SYMBOL , fromType= "SYMBOL" , toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
go1_MM1 <-enrichGO(gene = ENID$ENTREZID,
                   OrgDb          = org.Hs.eg.db, 
                   keyType        = "ENTREZID",
                   ont            = "BP",         
                   pAdjustMethod  = "BH",
                   pvalueCutoff   = 0.1,         
                   qvalueCutoff   = 0.05,
                   minGSSize      = 10,           
                   maxGSSize      = 500,          
                   readable       = FALSE,
                   pool           = FALSE
)
go1_MM1_symbol <- setReadable(go1_MM1, 'org.Hs.eg.db',keyType = "ENTREZID")
go1_MM1top <- as.data.frame(go1_MM1_symbol)
BP_MM1 <- go1_MM1top$Description[1:20]

gene1<- filter(cluster.markers, cluster == "Melanoma_2")
fivenum(gene1$avg_log2FC)
gene <- filter(cluster.markers, cluster == "Melanoma_2" & avg_log2FC >= 0.37 & p_val_adj < 0.05)
SYMBOL <- rownames(gene)
ENID <- bitr(SYMBOL , fromType= "SYMBOL" , toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
go1_MM2 <-enrichGO(gene = ENID$ENTREZID,
                   OrgDb          = org.Hs.eg.db, 
                   keyType        = "ENTREZID",
                   ont            = "BP",        
                   pAdjustMethod  = "BH",
                   pvalueCutoff   = 0.1,         
                   qvalueCutoff   = 0.05,
                   minGSSize      = 10,          
                   maxGSSize      = 500,          
                   readable       = FALSE,
                   pool           = FALSE
)
go1_MM2_symbol <- setReadable(go1_MM2, 'org.Hs.eg.db',keyType = "ENTREZID")
go1_MM2top <- as.data.frame(go1_MM2_symbol)
BP_MM2 <- go1_MM2top$Description[1:20]

gene1<- filter(cluster.markers, cluster == "Melanoma_3")
fivenum(gene1$avg_log2FC)
gene <- filter(cluster.markers, cluster == "Melanoma_3" & avg_log2FC >= 0.7 & p_val_adj < 0.05)
SYMBOL <- rownames(gene)
ENID <- bitr(SYMBOL , fromType= "SYMBOL" , toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
go1_MM3 <-enrichGO(gene = ENID$ENTREZID,
                   OrgDb          = org.Hs.eg.db, 
                   keyType        = "ENTREZID",
                   ont            = "BP",         
                   pAdjustMethod  = "BH",
                   pvalueCutoff   = 0.1,         
                   qvalueCutoff   = 0.05,
                   minGSSize      = 10,          
                   maxGSSize      = 500,         
                   readable       = FALSE,
                   pool           = FALSE
)
go1_MM3_symbol <- setReadable(go1_MM3, 'org.Hs.eg.db',keyType = "ENTREZID")
go1_MM3top <- as.data.frame(go1_MM3_symbol)
BP_MM3 <- go1_MM3top$Description[1:20]

goMM1 <- go1_MM1top[1:20,]
goMM1$Phenotype <- "Melanoma_1"
goMM2 <- go1_MM2top[1:20,]
goMM2$Phenotype <- "Melanoma_2"
goMM3 <- go1_MM3top[1:20,]
goMM3$Phenotype <- "Melanoma_3"
go.df <- rbind(goMM1,goMM2,goMM3)
go.df$Description <- factor(go.df$Description, levels = rev(go.df$Description))

COLS <- c("#000000", 
          "#d3d3d3", 
          "#b80c09" 
)
go_bar <- ggplot(data = go.df, 
                 aes(x = Description, y = Count, fill = Phenotype))+ 
  geom_bar(stat = "identity", width = 0.9)+ 
  scale_fill_manual(values = COLS)+
  theme_bw()+ 
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ 
  labs(x = "GO terms", y = "GeneNumber", title = "Barplot of Enriched GO Terms")+ 
  theme(axis.title = element_text(size = 20), 
        axis.text = element_text(size = 11), 
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 11), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm") ,
        axis.text.x = element_text(face = "bold", color="gray50", angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "top",
        panel.grid = element_blank()
  ) 
ggsave(go_bar,filename = "fig5A.pdf", width = 20, height = 10)


########################################################################################

                                      #4、VlnPlot

########################################################################################


scRNA <- subset(MM_1, idents=c("Melanoma_1", "Melanoma_2", "Melanoma_3"))
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)
gene <- c("MLANA", "MITF", "TYR","PMEL")

pdf("fig5B.pdf")
VlnPlot(scRNA, features = gene, log=T)
dev.off()


gene <- c("RB1", "PCNA", "CDK2","CDK4", "CDK6", "CENPF")

pdf("fig7B.pdf")
VlnPlot(scRNA, features = gene, log=T)
dev.off()


