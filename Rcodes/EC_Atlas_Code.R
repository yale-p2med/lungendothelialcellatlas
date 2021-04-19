############################################
# Analysis and Visualizations of the EC project

library(dplyr)
library(Seurat)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyr)
library(parallel)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(future)
library(seriation)
  
# load necessary all objects
FullDataset_Object <- readRDS(file = "FullDataset_Object.rds")
EC_Subset_Object <- readRDS(file = "EC_Subset_Object.rds")

load("mouse.Robj")
load("mouse_EC.Robj")
load("CombinedSpeciesECs_obj.Robj")
load("color_vectors.Robj")

# reorder for plotting
FullDataset_Object$cell.type.ident <- factor(FullDataset_Object$cell.type.ident, 
	levels=c("AM", "M", "cMono", "ncMono", "cDC1", "cDC2", "DC_mature", "DC_Langerhans", "pDC", "Mast",
	"T_helper", "T_cytotox", "T_reg", "NK","ILC", "B", "Plasma", 
	"AT1", "AT2", "Basal", "Suprabasal", "Goblet", "Club", "Ciliated","Ionocyte","PNEC",
	"Fibroblast_adventitial", "Fibroblast_alveolar", "SMC","Pericyte","Mesothel",
	"EC_arterial",  "EC_capillary_A", "EC_capillary_B", "EC_venous","EC_bronchial","EC_lymph", "Cell_Cycle"))


#############################################
######### UMAPs
# UMAP of EC object by cell type, cohort, subject
pdf(file = "EC_Subset_Object.UMAP.cohort.ident.pdf", useDingbats=FALSE)
DimPlot(EC_Subset_Object, reduction = "umap", group.by="cohort.ident", cols=cohort_colors, cells = sample(colnames(EC_Subset_Object)))
dev.off()

pdf(file = "EC_Subset_Object.UMAP.cell.type.ident.pdf", useDingbats=FALSE)
DimPlot(EC_Subset_Object, reduction = "umap", group.by="cell.type.ident", cols=celltype_colors, cells = sample(colnames(EC_Subset_Object)))
dev.off()

pdf(file = "EC_Subset_Object.UMAP.subject.ident.pdf", useDingbats=FALSE)
DimPlot(EC_Subset_Object, reduction = "umap", group.by="subject.ident", cols=subject_colors, cells = sample(colnames(EC_Subset_Object))) + NoLegend()
dev.off()

pdf(file = "EC_Subset_Object.UMAP.sample.location.pdf", useDingbats=FALSE)
DimPlot(EC_Subset_Object, reduction = "umap", group.by="location.ident", cols=location_colors, cells = sample(colnames(EC_Subset_Object)))
dev.off()


#########
# UMAP of Full object by cell type, cohort, subject
# dont use pdfs as it is to big	   
p1 <- DimPlot(FullDataset_Object, reduction = "umap", group.by="cell.type.ident", label=FALSE, cols=celltype_colors, cells = sample(colnames(FullDataset_Object)))
png(file = "FullDataset_Object.newUMAP.cell.type.ident.png", width = 10, height = 10,units = 'in', res = 300)
p1 + NoLegend()
dev.off()
p2 <- DimPlot(FullDataset_Object, reduction = "umap", group.by="cohort.ident", cols=cohort_colors, label=FALSE, cells = sample(colnames(FullDataset_Object)))
png(file = "FullDataset_Object.newUMAP.cohort.ident.png", width = 10, height = 10,units = 'in', res = 300)
p2 + NoLegend()
dev.off()
p3 <- DimPlot(FullDataset_Object, reduction = "umap", group.by="subject.ident", cols=subject_colors, label=FALSE, cells = sample(colnames(FullDataset_Object)))
png(file = "FullDataset_Object.newUMAP.subject.ident.png", width = 10, height = 10,units = 'in', res = 300)
p3 + NoLegend()
dev.off()
png(file = "FullDataset_Object.newUMAP.cell.type.ident.withLabels.png", width = 10, height = 10,units = 'in', res = 300)
DimPlot(FullDataset_Object, reduction = "umap", group.by="cell.type.ident", label=TRUE, cols=celltype_colors, 
	cells = sample(colnames(FullDataset_Object))) + NoLegend()
dev.off()

pdf("FullDataset_Object.newUMAP.Legends.pdf", useDingbats=FALSE)
grid.draw(cowplot::get_legend(p1))
grid.newpage()
grid.draw(cowplot::get_legend(p2))
grid.newpage()
grid.draw(cowplot::get_legend(p3))
dev.off()


#########
# UMAP of full mouse object by cell type, cohort, subject
# dont use pdfs as it is to big	   
p1 <- DimPlot(mouse, reduction = "umap", group.by="lineage.ident", label=FALSE, cols=lineage_colors, cells = sample(colnames(mouse)))
png(file = "mouse.UMAP.lineage.ident.png", width = 10, height = 10,units = 'in', res = 300)
p1 + NoLegend()
dev.off()
p2 <- DimPlot(mouse, reduction = "umap", group.by="orig.ident", label=FALSE, cols=mouse_subject_colors, cells = sample(colnames(mouse)))
png(file = "mouse.UMAP.orig.ident.png", width = 10, height = 10,units = 'in', res = 300)
p2 + NoLegend()
dev.off()
p3 <- DimPlot(mouse, reduction = "umap", group.by="cohort.ident", label=FALSE, cols=mouse_cohort_colors, cells = sample(colnames(mouse)))
png(file = "mouse.UMAP.cohort.ident.png", width = 10, height = 10,units = 'in', res = 300)
p3 + NoLegend()
dev.off()
p4 <- DimPlot(mouse, reduction = "umap", group.by="cell.type.ident", label=FALSE, cols=c(celltype_colors), 
    cells = sample(colnames(mouse)))
png(file = "mouse.UMAP.celltype.ident.png", width = 10, height = 10,units = 'in', res = 300)
p4 + NoLegend()
dev.off()

pdf("mouse.UMAP.Legends.pdf", useDingbats=FALSE)
grid.draw(cowplot::get_legend(p1))
grid.newpage()
grid.draw(cowplot::get_legend(p2))
grid.newpage()
grid.draw(cowplot::get_legend(p3))
grid.newpage()
grid.draw(cowplot::get_legend(p4))
dev.off()


pdf(file = "mouse_EC.UMAP.cell.type.ident.pdf", useDingbats=FALSE)
DimPlot(mouse_EC, reduction = "umap", group.by="cell.type.ident", label=FALSE, cols=celltype_colors, cells = sample(colnames(mouse_EC)))
dev.off()
pdf(file = "mouse_EC.UMAP.orig.ident.pdf", useDingbats=FALSE)
DimPlot(mouse_EC, reduction = "umap", group.by="orig.ident", label=FALSE, cols=mouse_subject_colors, cells = sample(colnames(mouse_EC)))
dev.off()
pdf(file = "mouse_EC.UMAP.cohort.ident.pdf", useDingbats=FALSE)
DimPlot(mouse_EC, reduction = "umap", group.by="cohort.ident", label=FALSE, cols=mouse_cohort_colors, cells = sample(colnames(mouse_EC)))
dev.off()

png(file = "mouse_EC.UMAP.cell.type.ident.png", width = 10, height = 10,units = 'in', res = 300)
DimPlot(mouse_EC, reduction = "umap", group.by="cell.type.ident", label=FALSE, cols=celltype_colors, cells = sample(colnames(mouse_EC)), pt.size=1) + NoLegend()
dev.off()
png(file = "mouse_EC.UMAP.orig.ident.png", width = 10, height = 10,units = 'in', res = 300)
DimPlot(mouse_EC, reduction = "umap", group.by="orig.ident", label=FALSE, cols=mouse_subject_colors, cells = sample(colnames(mouse_EC)), pt.size=1) + NoLegend()
dev.off()
png(file = "mouse_EC.UMAP.cohort.ident.png", width = 10, height = 10,units = 'in', res = 300)
DimPlot(mouse_EC, reduction = "umap", group.by="cohort.ident", label=FALSE, cols=mouse_cohort_colors, cells = sample(colnames(mouse_EC)), pt.size=1) + NoLegend()
dev.off()


######
# UMAPs of combined species ECs
pdf("UMAP.CombinedSpeciesECs_obj.CellTypeIdent.pdf", useDingbats=FALSE)
DimPlot(CombinedSpeciesECs_obj, group.by="cell.type.ident", cells = sample(colnames(CombinedSpeciesECs_obj)), cols=celltype_colors)
dev.off()
pdf("UMAP.CombinedSpeciesECs_obj.CohortIdent.pdf", useDingbats=FALSE)
DimPlot(CombinedSpeciesECs_obj, group.by="cohort.ident", cells = sample(colnames(CombinedSpeciesECs_obj)), cols=c(cohort_colors, mouse_cohort_colors))
dev.off()
pdf("UMAP.CombinedSpeciesECs_obj.Species.pdf", useDingbats=FALSE)
DimPlot(CombinedSpeciesECs_obj, group.by="species", cells = sample(colnames(CombinedSpeciesECs_obj)), cols=species_colors)
dev.off()


####################
# get differentials including DOR
# group.name is a metadat column name, e.g. "celltype", logFC.min.filter ist the minimal log FC filter
myMultithreadFindMarkersDOR <- function(seurat.object, group.name, logFC.min.filter, ncores){
    library(parallel)
    seurat.object$cellBarcode <- colnames(seurat.object)
    num.cellsTotal <- nrow(seurat.object@meta.data)
    group.levels <- as.vector(unique(seurat.object@meta.data[[group.name]])) #maybe as vector necessary in case of character data
    # now as function to multithread, input is a vector of group.levels
    MultithreadFindMarkersDOR <- function(x){
        #### get the index of cell barcodes that belong to that cell type
            temp.cells <- seurat.object@meta.data %>% filter(seurat.object@meta.data[,group.name] == x) %>% pull(cellBarcode)
            temp.avgExpAllGenes.inGroup <- rowMeans(seurat.object@assays$RNA@data[,temp.cells])
            temp.avgExpAllGenes.outGroup <- rowMeans(seurat.object@assays$RNA@data[,!colnames(seurat.object@assays$RNA@data) %in% temp.cells])
            temp.all.logFC <- temp.avgExpAllGenes.inGroup - temp.avgExpAllGenes.outGroup
        ##### Filter out genes with a FC < logFC.min.filter
            temp.filt.logFC <- temp.all.logFC[temp.all.logFC >= logFC.min.filter]
            temp.numGenesFilt <-  length(temp.filt.logFC)
        ##### gather the rest of the information for each of these genes    
            cat("Computing statistics for ", temp.numGenesFilt, " genes for ", x,".\n", sep="")
        ##### get basic info on number of cells    
            num.cellsInGroup <- length(temp.cells)
            num.cellsOutGroup <- num.cellsTotal - length(temp.cells)
        #### get number of cells within & outside the group with these genes    
            num.TruePos <- rowSums(seurat.object@assays$RNA@data[names(temp.filt.logFC),temp.cells] > 0)
            num.FalsePos <- rowSums(seurat.object@assays$RNA@data[names(temp.filt.logFC), !colnames(seurat.object@assays$RNA@data) %in% temp.cells] > 0)
        #### get number of cells in and outside group without genes
            num.FalseNeg <- num.cellsInGroup - num.TruePos
            num.TrueNeg <- num.cellsOutGroup - num.FalsePos
        ###### use these values to calculate log(DOR) w/ a pseudocount of 0.5 to avoid +/- infinity values         
            temp.logDOR <- log((num.TruePos+0.5)/(num.FalsePos+0.5)/((num.FalseNeg+0.5)/(num.TrueNeg+0.5)))
            temp.data <-seurat.object@assays$RNA@data[names(temp.filt.logFC),]
            p_val <- sapply(X = 1:nrow(temp.data),
                FUN=function(x){
                    return(wilcox.test(x=temp.data[x,temp.cells], 
                        y=temp.data[x,seurat.object$cellBarcode[!seurat.object$cellBarcode %in% temp.cells]])$p.value)})
            names(p_val) <- names(temp.filt.logFC)
            p_val_adj <- p.adjust(p = p_val, method = "bonferroni", n = nrow(x = seurat.object))
        #### put everything I'm interested in into one big table
            temp.df <- data.frame(group.name = rep(x, temp.numGenesFilt),
                gene=names(temp.filt.logFC),
                logFC=temp.filt.logFC,
                logDOR=temp.logDOR,
                pct.1=num.TruePos/num.cellsInGroup,
                pct.2=num.FalsePos/num.cellsOutGroup,
                p_val=p_val,
                p_val_adj=p_val_adj,
                nCells.in=num.cellsInGroup,
                nCells.out=num.cellsOutGroup)
            colnames(temp.df)[colnames(temp.df)=="group.name"] <- group.name
        #### rearrange in order of logDOR; highest to lowest
            return(temp.df %>% arrange(desc(logDOR)))
    }
    outputList <- mclapply(group.levels, MultithreadFindMarkersDOR, mc.cores=8)
    return(Reduce(rbind, outputList))
}

### EC
EC.markers.DOR <- myMultithreadFindMarkersDOR(
    seurat.object=EC_Subset_Object, 
    group.name="cell.type.ident", 
    logFC.min.filter=0.25, ncores=8)
write.table(EC.markers.DOR, "EC.markers.DOR.txt", row.names = FALSE, sep="\t")

###  FullDataset_Object
FullDataset_Object.markers.DOR <- myMultithreadFindMarkersDOR(
    seurat.object=FullDataset_Object, 
    group.name="cell.type.ident", 
    logFC.min.filter=0.25, ncores=8)

EC_general_genes <- Reduce(intersect, list(
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_venous") %>% pull(gene),
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_A") %>% pull(gene),
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_B") %>% pull(gene),
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_arterial") %>% pull(gene),
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_lymph") %>% pull(gene)))

FullDataset_Object.markers.DOR$EC_general_marker <- "no"
FullDataset_Object.markers.DOR$EC_general_marker[ FullDataset_Object.markers.DOR$gene %in% EC_general_genes & 
	FullDataset_Object.markers.DOR$cell.type.ident %in% c("EC_venous", "EC_capillary_A", "EC_capillary_B", "EC_arterial", "EC_lymph")] <- "yes" 

EC_vasc_genes <- Reduce(intersect, list(
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_venous") %>% pull(gene),
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_A") %>% pull(gene),
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_B") %>% pull(gene),
	FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_arterial") %>% pull(gene)))
	
EC_vasc_genes <- EC_vasc_genes[! EC_vasc_genes %in%  (FullDataset_Object.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_lymph") %>% pull(gene))]

FullDataset_Object.markers.DOR$EC_vasc_marker <- "no"
FullDataset_Object.markers.DOR$EC_vasc_marker[ FullDataset_Object.markers.DOR$gene %in% EC_vasc_genes & 
	FullDataset_Object.markers.DOR$cell.type.ident %in% c("EC_venous", "EC_capillary_A", "EC_capillary_B", "EC_arterial")] <- "yes" 

write.table(FullDataset_Object.markers.DOR, "FullDataset_Object.markers.DOR.txt", row.names = FALSE, sep="\t")

### mouse_EC
mouse_EC.markers.DOR <- myMultithreadFindMarkersDOR(
    seurat.object=mouse_EC, 
    group.name="cell.type.ident", 
    logFC.min.filter=0.25, ncores=8)
write.table(mouse_EC.markers.DOR, "mouse_EC.markers.DOR.txt", row.names = FALSE, sep="\t")

###  mouse full
mouse.markers.DOR <- myMultithreadFindMarkersDOR(
    seurat.object=mouse, 
    group.name="cell.type.ident", 
    logFC.min.filter=0.25, ncores=8)

mouse_EC_general_genes_woBronchial <- Reduce(intersect, list(
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_venous") %>% pull(gene),
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_A") %>% pull(gene),
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_B") %>% pull(gene),
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_arterial") %>% pull(gene),
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_lymph") %>% pull(gene)))

mouse.markers.DOR$EC_general_marker <- "no"
mouse.markers.DOR$EC_general_marker[ mouse.markers.DOR$gene %in% mouse_EC_general_genes_woBronchial & 
	mouse.markers.DOR$cell.type.ident %in% c("EC_venous", "EC_capillary_A", "EC_capillary_B", "EC_arterial","EC_lymph")] <- "yes" 

mouse_EC_vasc_genes_woBronchial <- Reduce(intersect, list(
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_venous") %>% pull(gene),
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_A") %>% pull(gene),
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_capillary_B") %>% pull(gene),
	mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_arterial") %>% pull(gene)))

mouse_EC_vasc_genes_woBronchial <- mouse_EC_vasc_genes_woBronchial[! mouse_EC_vasc_genes_woBronchial %in%  (mouse.markers.DOR %>% filter(p_val_adj<0.05 & logFC>0.25 & cell.type.ident=="EC_lymph") %>% pull(gene))]

mouse.markers.DOR$EC_vasc_marker <- "no"
mouse.markers.DOR$EC_vasc_marker[ mouse.markers.DOR$gene %in% mouse_EC_vasc_genes_woBronchial & 
	mouse.markers.DOR$cell.type.ident %in% c("EC_venous", "EC_capillary_A", "EC_capillary_B", "EC_arterial")] <- "yes" 
	
write.table(mouse.markers.DOR, "mouse.markers.DOR.txt", row.names = FALSE, sep="\t")


#################################################
#### make heatmaps
#### a) FullDataset_Object marker genes per average per subject and per celltype
FullDataset_Object$lineage.ident <- rep(NA)
FullDataset_Object$lineage.ident[FullDataset_Object$cell.type.ident %in% c("AM","cDC1", "cDC2", "cMono", "DC_mature","M","Mast","ncMono", 
	"pDC","DC_Langerhans")] <- "myeloid"
FullDataset_Object$lineage.ident[FullDataset_Object$cell.type.ident %in% c("B","ILC","NK","Plasma", "T_cytotox", "T_helper", "T_reg","T_gammadelta")] <- "lymphoid"
FullDataset_Object$lineage.ident[FullDataset_Object$cell.type.ident %in% c("AT1", "AT2", "Basal", "Ciliated","Club" ,"Goblet", "Suprabasal", 
	"Ionocyte", "PNEC")] 	<- "epithelial"
FullDataset_Object$lineage.ident[FullDataset_Object$cell.type.ident %in% c("EC_arterial", "EC_bronchial", "EC_capillary_A", "EC_capillary_B","EC_lymph", 
	"EC_venous")] <- "endothelial"
FullDataset_Object$lineage.ident[FullDataset_Object$cell.type.ident %in% c("Fibroblast_adventitial", "Fibroblast_alveolar","Mesothel", "Pericyte",
	"SMC")] <- "stromal"

EC_vasc_genes_to_plot <- c("SLCO2A1", "CLEC14A", "PCDH17", "EPAS1", "LIMS2", "ADGRL4","AQP1", "IFI27", "FLT1", "BMPR2", 
	  "SHROOM4","ENG", "ESAM", "ITM2A")
	
EC_general_genes_to_plot <- c("CALCRL", "LDB2", "PCAT19", "CLDN5","RAMP2","ERG","PECAM1", "CDH5", "EGFL7", 
	 "TIE1", "GNG11", "HYAL2", "TM4SF1", "PODXL", "CAVIN2", "ARHGAP29", "ST6GALNAC3", "DOCK9", "VAMP5", "SASH1")
	
genes_for_heatmap <- c( "FABP4","C1QB",  # AM
	"TREM2", "SPP1",    # M     
	"S100A12", "VCAN", # cMono
	"LILRB1", "LILRB2", # ncMono
	"CLEC9A", "CLNK", # cDC1
	"TREM2", "CLEC10A",    # cDC2 
	"CCR7", "CCL22", # DC_mature
	"FCER1A", "CD1E", # DC_Langerhans
	"SCT", "IRF4", # pDC
	"TPSAB1", "TPSB2", # Mast
	"THEMIS", "IL7R",  # T_helper
	"CD8A", "CD8B",  # T_cytotox
	"FOXP3", "CTLA4",  # T_reg
	"NKG7", "KLRD1",  # NK
	"KLRC1", "XCL1",  # ILC
	"MS4A1","CD79A",  # B
	"IGHG1", "MZB1",  # Plasma
	"AGER","CLDN18", # AT1
	"SFTPC","SFTPA1", # AT2
	"TP63","MIR205HG", # Basal
	"SERPINB3","SERPINB4", # Suprabasal
	"BPIFB1","SCGB3A1", # Goblet	
	"SCGB3A2","CLIC6", # Club	
	"HYDIN","CFAP299", # Ciliated
	"FOXI1","PDE1C", # Ionocytes  
	"CHGA","GRP", # PNECs
	"MFAP5","SCARA5", # Fibro
	"ITGA8","SCN7A", # Myofibro
	"MYH11","TAGLN", # SMC
	"COX4I2","PDGFRB", # Pericyte	
	"MSLN","CALB2", # Mesothel
	"MMRN1", "CCL21", "PROX1", "PKHD1L1", "SEMA3D", "TFF3","TM4SF18","FLT4", "LYVE1","PDPN", "SNCG", "SCN3B", "TBX1", "LINC02147", "RELN","KLHL4", # EC_lymph	
	sample(EC_general_genes_to_plot),# EC_vasc+EC_lymph
	sample(EC_vasc_genes_to_plot))	# only EC_vasc

# also determines the order
celltypes_to_plot <- c("AM", "M", "cMono", "ncMono", "cDC1", "cDC2", "DC_mature", "DC_Langerhans", "pDC", "Mast",
	"T_helper", "T_cytotox", "T_reg", "NK","ILC", "B", "Plasma", 
	"AT1", "AT2", "Basal", "Suprabasal", "Goblet", "Club", "Ciliated","Ionocyte","PNEC",
	"Fibroblast_adventitial", "Fibroblast_alveolar", "SMC","Pericyte","Mesothel",
	"EC_lymph","EC_arterial",  "EC_capillary_A", "EC_capillary_B", "EC_venous","EC_bronchial")

### with multithreading, first make all possible combinations
FullDataset_Object$cellBarcode <- colnames(FullDataset_Object) 
cellTypes <- levels(as.factor(FullDataset_Object$cell.type.ident)) 
cellTypes <- cellTypes[cellTypes %in% celltypes_to_plot]
meta.data.sub <- FullDataset_Object@meta.data[,c("cohort.ident", "cell.type.ident", "subject.ident", "cellBarcode")]

get.CT.DS.subj.vector <- function(cellTypes){
	tmp.meta.data <- meta.data.sub %>% filter(cell.type.ident== cellTypes)
	cohorts <- unique(tmp.meta.data$cohort.ident)
	subjects <- unique(tmp.meta.data$subject.ident)
	tmp.CT.DS.subj.vector <- vector()
	for(j in 1:length(cohorts)){ 
		for(k in 1:length(subjects)){
		temp.cells <- tmp.meta.data %>% filter(cohort.ident==cohorts[j] & subject.ident==subjects[k]) %>% pull(cellBarcode)
		if ( length(temp.cells) >1 ) { # >=
			tmp.CT.DS.subj.vector <- c(tmp.CT.DS.subj.vector, paste(cellTypes, cohorts[j], subjects[k], sep="__"))
			}
		}
	}
	cat("Completed for ", cellTypes, ".\n", sep="")
	return(tmp.CT.DS.subj.vector)
}
celltype_cohort_subject.list <- parallel::mclapply(cellTypes, get.CT.DS.subj.vector, mc.cores=8)
celltype_cohort_subject <- unlist(celltype_cohort_subject.list)

DefaultAssay(FullDataset_Object) <- "RNA"
get.SubjectcohortCellTypeAvg <- function(celltype_cohort_subject){
	temp.cell.type <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][1]
	temp.cohort <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][2]
	temp.subject <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][3]
	temp.meta.data <- FullDataset_Object@meta.data[,c("cohort.ident", "cell.type.ident", "subject.ident", "cellBarcode")]
	temp.cells <- temp.meta.data %>% filter(cell.type.ident==temp.cell.type & cohort.ident==temp.cohort & 
		subject.ident==temp.subject) %>% pull(cellBarcode) 
	if (length(temp.cells) > 1) { 
		tmp.df <- as.data.frame(rowMeans(GetAssayData(FullDataset_Object)[,temp.cells]))	   
    } else { 
		tmp.df <- as.data.frame(GetAssayData(FullDataset_Object)[,temp.cells])
		cat("Subject",temp.subject,"only has 1",temp.cell.type,"cell, using singlet for",temp.cohort,"representation...\n",sep=" ") 
    } 
	colnames(tmp.df) <- paste(celltype_cohort_subject)
	return(tmp.df)
}

collapsed.mtx.list <- parallel::mclapply(celltype_cohort_subject, get.SubjectcohortCellTypeAvg, mc.cores=8)
collapsed.SubjectcohortCellTypeAvg.mtx <- Matrix(as.matrix(do.call(cbind, collapsed.mtx.list)), sparse = TRUE) 
dim(collapsed.SubjectcohortCellTypeAvg.mtx) 

heatmap_metadata <- as.data.frame(cbind(colnames(collapsed.SubjectcohortCellTypeAvg.mtx), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 1), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 2), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 3)))
colnames(heatmap_metadata) <- c("cell.ident","cell.type.ident", "cohort.ident", "subject.ident")

heatmap_metadata$cell.type.ident <- factor(heatmap_metadata$cell.type.ident, levels=celltypes_to_plot)
heatmap_metadata$cohort.ident <- factor(heatmap_metadata$cohort.ident, levels=c("Yale_BWH","Vanderbilt_TGen","WSI_Groningen", 
	"Northwestern","Leuven_VIB" ))

cell_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cell.ident)
cohort_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cohort.ident)
celltype_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cell.type.ident)
subject_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(subject.ident )

heatmap_df1 <-  as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[genes_for_heatmap,cell_order])

lineage_for_heatmap <- c( "LYZ",  "PTPRC", "CD69",  "EPCAM", "CDH1", "PDGFRB", "COL1A2", "PECAM1" ,"CLDN5")
heatmap_df2 <- as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[lineage_for_heatmap,cell_order])

heatmap_df <- rbind(heatmap_df2, heatmap_df1)

heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))

heatmap_cohort_annotation <- HeatmapAnnotation(cell_type=celltype_order, cohort=cohort_order, subject=subject_order, 
    col = list(cohort=cohort_colors, cell_type=celltype_colors, subject=subject_colors),
	show_legend = c("subject" = FALSE))
	
pdf("FullDataset_Object.Heatmap.Markers.AveragePerSubject.pdf", height=10, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=celltype_order, column_title = NULL, use_raster=FALSE)
dev.off()



####
# first EC only
ECgenes_for_heatmap <- c("MMRN1", "CCL21", "PROX1", "PKHD1L1", "SEMA3D", "TFF3","TM4SF18","FLT4", "LYVE1","PDPN", "SNCG", "SCN3B", "TBX1", "LINC02147", "RELN","KLHL4", # EC_lymph	
	sample(EC_general_genes_to_plot),# EC_vasc+EC_lymph
	sample(EC_vasc_genes_to_plot))
heatmap_df <-  as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[ECgenes_for_heatmap,cell_order])
heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))
heatmap_cohort_annotation <- HeatmapAnnotation(cell_type=celltype_order, cohort=cohort_order, subject=subject_order, 
    col = list(cohort=cohort_colors, cell_type=celltype_colors, subject=subject_colors),
	show_legend = c("subject" = FALSE))
pdf("FullDataset_Object.Heatmap.EConly.Markers.AveragePerSubject.8_14_20.pdf", height=5, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=celltype_order, column_title = NULL, use_raster=FALSE)
dev.off()

# now lineage heatmap
heatmap_metadata$lineage.ident <- rep(NA)
heatmap_metadata$lineage.ident[heatmap_metadata$cell.type.ident %in% c("AM","cDC1", "cDC2", "cMono", "DC_mature","M","Mast","ncMono", 
	"pDC","DC_Langerhans")] <- "myeloid"
heatmap_metadata$lineage.ident[heatmap_metadata$cell.type.ident %in% c("B","ILC","NK","Plasma", "T_cytotox", "T_helper", "T_reg","T_gammadelta")] <- "lymphoid"
heatmap_metadata$lineage.ident[heatmap_metadata$cell.type.ident %in% c("AT1", "AT2", "Basal", "Ciliated","Club" ,"Goblet", "Suprabasal", 
	"Ionocyte", "PNEC")] 	<- "epithelial"
heatmap_metadata$lineage.ident[heatmap_metadata$cell.type.ident %in% c("EC_arterial", "EC_bronchial", "EC_capillary_A", "EC_capillary_B","EC_lymph", 
	"EC_venous")] <- "endothelial"
heatmap_metadata$lineage.ident[heatmap_metadata$cell.type.ident %in% c("Fibroblast_adventitial", "Fibroblast_alveolar","Mesothel", "Pericyte",
	"SMC")] <- "stromal"
heatmap_metadata$lineage.ident <- factor(heatmap_metadata$lineage.ident, levels=c("myeloid", "lymphoid", "epithelial", "stromal", "endothelial"))

cell_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cell.ident)
cohort_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cohort.ident)
celltype_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cell.type.ident)
subject_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(subject.ident )
lineage_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(lineage.ident)

heatmap_df1 <-  as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[genes_for_heatmap,cell_order])
lineage_for_heatmap <- c( "LYZ",  "PTPRC", "CD69",  "EPCAM", "CDH1", "PDGFRB", "COL1A2", "PECAM1" ,"CLDN5")
heatmap_df2 <- as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[lineage_for_heatmap,cell_order])
heatmap_df <- rbind(heatmap_df2, heatmap_df1)

heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))

heatmap_cohort_annotation <- HeatmapAnnotation(Lineage=lineage_order, Cell_type=celltype_order, Cohort=cohort_order, Subject=subject_order, 
    col = list(Cohort=cohort_colors, Cell_type= celltype_colors,Lineage=lineage_colors, Subject=subject_colors),
	show_legend = c("Subject" = FALSE))

pdf("FullDataset_Object.Heatmap.Lineage.Markers.AveragePerSubject.8_14_20.pdf", height=10, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=lineage_order, column_title = NULL, use_raster=FALSE)
dev.off()
pdf("FullDataset_Object.Heatmap.CellType.Markers.AveragePerSubject.8_14_20.pdf", height=10, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=celltype_order, column_title = NULL, use_raster=FALSE)
dev.off()


#### b) EC_Subset_Object marker genes per average per subject and per celltype
celltypes_to_plot <- as.factor(c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))

genes_for_heatmap <- c("DKK2", "IGFBP3","FBLN5", "SERPINE2", "CLDN10", "GJA5", "CXCL12","BMX", "LTBP4", "HEY1","SOX5","SEMA3G",  # arterial
	"SOSTDC1", "EDNRB", "HPGD", "CYP3A5", "PRKG1", "TBX2", "RCSD1", "EDA", "B3GALNT1", "EXPH5", "NCALD","S100A4", # cap A
	"CA4", "AFF3", "ADGRL2","BTNL9","RGCC", "ADGRF5", "KIAA1217",#common capA and capB
	"FCN3", "IL7R", "CD36",  "NRXN3","SLC6A4", "GPIHBP1","ARHGAP18", "IL18R1", # cap B
	"CPE", "CLU", "C7","PTGS1","EFEMP1","MMRN1","PKHD1L1", "PDZRN4", "DKK3","PLAT", "CDH11","HDAC9", # venous
	"ACKR1","IGFBP7","MCTP1", "VWF", #common bronch and venous
	"COL15A1",  "ZNF385D", "EBF1", "TSHZ2","FLRT2","OLFM1", "CPXM2","PLVAP", "TPD52L1","PDE7B","VWA1", "SPRY1") # peribronchial
	
###  with multithreading, first make all possible combinations
EC_Subset_Object$cellBarcode <- colnames(EC_Subset_Object) 
cellTypes <- levels(as.factor(EC_Subset_Object$cell.type.ident)) 
cellTypes <- cellTypes[cellTypes %in% celltypes_to_plot]
meta.data.sub <- EC_Subset_Object@meta.data[,c("cohort.ident", "cell.type.ident", "subject.ident", "cellBarcode")]

get.CT.DS.subj.vector <- function(cellTypes){
	tmp.meta.data <- meta.data.sub %>% filter(cell.type.ident== cellTypes)
	cohorts <- unique(tmp.meta.data$cohort.ident)
	subjects <- unique(tmp.meta.data$subject.ident)
	tmp.CT.DS.subj.vector <- vector()
	for(j in 1:length(cohorts)){ 
		for(k in 1:length(subjects)){
		temp.cells <- tmp.meta.data %>% filter(cohort.ident==cohorts[j] & subject.ident==subjects[k]) %>% pull(cellBarcode)
		if ( length(temp.cells) >1 ) { # >=
			tmp.CT.DS.subj.vector <- c(tmp.CT.DS.subj.vector, paste(cellTypes, cohorts[j], subjects[k], sep="__"))
			}
		}
	}
	cat("Completed for ", cellTypes, ".\n", sep="")
	return(tmp.CT.DS.subj.vector)
}
celltype_cohort_subject.list <- parallel::mclapply(cellTypes, get.CT.DS.subj.vector, mc.cores=8)
celltype_cohort_subject <- unlist(celltype_cohort_subject.list)

DefaultAssay(EC_Subset_Object) <- "RNA"
get.SubjectcohortCellTypeAvg <- function(celltype_cohort_subject){
	temp.cell.type <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][1]
	temp.cohort <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][2]
	temp.subject <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][3]
	temp.meta.data <- EC_Subset_Object@meta.data[,c("cohort.ident", "cell.type.ident", "subject.ident", "cellBarcode")]
	temp.cells <- temp.meta.data %>% filter(cell.type.ident==temp.cell.type & cohort.ident==temp.cohort & 
		subject.ident==temp.subject) %>% pull(cellBarcode) 
	if (length(temp.cells) > 1) { 
		tmp.df <- as.data.frame(rowMeans(GetAssayData(EC_Subset_Object)[,temp.cells]))	   
     } else { 
		tmp.df <- as.data.frame(GetAssayData(EC_Subset_Object)[,temp.cells])
		cat("Subject",temp.subject,"only has 1",temp.cell.type,"cell, using singlet for",temp.cohort,"representation...\n",sep=" ") 
    } 
	colnames(tmp.df) <- paste(celltype_cohort_subject)
	return(tmp.df)
}

collapsed.mtx.list <- parallel::mclapply(celltype_cohort_subject, get.SubjectcohortCellTypeAvg, mc.cores=8)
collapsed.SubjectcohortCellTypeAvg.mtx <- Matrix(as.matrix(do.call(cbind, collapsed.mtx.list)), sparse = TRUE) 
dim(collapsed.SubjectcohortCellTypeAvg.mtx) 

heatmap_metadata <- as.data.frame(cbind(colnames(collapsed.SubjectcohortCellTypeAvg.mtx), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 1), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 2), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 3)))
colnames(heatmap_metadata) <- c("cell.ident","cell.type.ident", "cohort.ident", "subject.ident")

heatmap_metadata$cell.type.ident <- factor(heatmap_metadata$cell.type.ident, levels=celltypes_to_plot)
heatmap_metadata$cohort.ident <- factor(heatmap_metadata$cohort.ident, levels=c("Yale_BWH","Vanderbilt_TGen","WSI_Groningen", 
	"Northwestern","Leuven_VIB", "Leuven_LKI"))

cell_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cell.ident)
cohort_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cohort.ident)
celltype_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(cell.type.ident)
subject_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, subject.ident ) %>% pull(subject.ident )

heatmap_df <-  as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[genes_for_heatmap,cell_order])
heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))

heatmap_cohort_annotation <- HeatmapAnnotation(cell_type=celltype_order, cohort=cohort_order, subject=subject_order, 
    col = list(cohort=cohort_colors, cell_type=celltype_colors, subject=subject_colors),
	show_legend = c("subject" = FALSE))

pdf("EC_Subset_Object.Heatmap.Markers.AveragePerSubject.pdf", height=10, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=celltype_order, column_title = NULL, use_raster=FALSE)
dev.off()



########################
####  EC Mouse marker genes per average per animal and per celltype
celltypes_to_plot <- as.factor(c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))

genes_for_heatmap <- c("Dkk2", "Gja5", "Cxcl12","Bmx", "Efnb2", "Fbln2", "Gja4", "Dll4", "Hey1", "Htra1", # arterial
	"Car4", "Emp2", "Ednrb","Rcsd1", "Prx", "Tbx2", "Fibin", "Igfbp7", "Itga3", "Kdr", # cap A
	"Gpihbp1","Plcb1", "Glp1r", "Bmp6", "Cd74","Aplnr", "Adgrl3", "Sntb1", "Kit", "AC165951.1", # cap B
	"Cpe", "Slc6a2", "Fgl2","Ptgs1", "Vcam1", "Nr2f2", "Hdac9", "Prcp", "Ackr3", "Rgs5", # venous
    "Vwf", #common bronch and venous
	"Ackr1", "Col15a1", "Meox2", "Ebf3", "Tshz2","Gja1","C1qtnf9", "Aqp7","Timp4", "Vwa1") # peribronchial

### with multithreading, first make all possible combinations
mouse_EC$cellBarcode <- colnames(mouse_EC) 
cellTypes <- levels(as.factor(mouse_EC$cell.type.ident)) 
cellTypes <- cellTypes[cellTypes %in% celltypes_to_plot]
meta.data.sub <- mouse_EC@meta.data[,c("cohort.ident", "cell.type.ident", "orig.ident", "cellBarcode")]

get.CT.DS.subj.vector <- function(cellTypes){
	tmp.meta.data <- meta.data.sub %>% filter(cell.type.ident== cellTypes)
	cohorts <- unique(tmp.meta.data$cohort.ident)
	subjects <- unique(tmp.meta.data$orig.ident)
	tmp.CT.DS.subj.vector <- vector()
	for(j in 1:length(cohorts)){ 
		for(k in 1:length(subjects)){
		temp.cells <- tmp.meta.data %>% filter(cohort.ident==cohorts[j] & orig.ident==subjects[k]) %>% pull(cellBarcode)
		if ( length(temp.cells) >1 ) { # >=
			tmp.CT.DS.subj.vector <- c(tmp.CT.DS.subj.vector, paste(cellTypes, cohorts[j], subjects[k], sep="__"))
			}
		}
	}
	cat("Completed for ", cellTypes, ".\n", sep="")
	return(tmp.CT.DS.subj.vector)
}
celltype_cohort_subject.list <- parallel::mclapply(cellTypes, get.CT.DS.subj.vector, mc.cores=8)
celltype_cohort_subject <- unlist(celltype_cohort_subject.list)

DefaultAssay(mouse_EC) <- "RNA"
get.SubjectcohortCellTypeAvg <- function(celltype_cohort_subject){
	temp.cell.type <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][1]
	temp.cohort <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][2]
	temp.subject <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][3]
	temp.meta.data <- mouse_EC@meta.data[,c("cohort.ident", "cell.type.ident", "orig.ident", "cellBarcode")]
	temp.cells <- temp.meta.data %>% filter(cell.type.ident==temp.cell.type & cohort.ident==temp.cohort & 
		orig.ident==temp.subject) %>% pull(cellBarcode) 
	if (length(temp.cells) > 1) { 
		tmp.df <- as.data.frame(rowMeans(GetAssayData(mouse_EC)[,temp.cells]))	   
     } else { 
		tmp.df <- as.data.frame(GetAssayData(mouse_EC)[,temp.cells])
		cat("Subject",temp.subject,"only has 1",temp.cell.type,"cell, using singlet for",temp.cohort,"representation...\n",sep=" ") 
    } 
	colnames(tmp.df) <- paste(celltype_cohort_subject)
	return(tmp.df)
}

collapsed.mtx.list <- parallel::mclapply(celltype_cohort_subject, get.SubjectcohortCellTypeAvg, mc.cores=8)
collapsed.SubjectcohortCellTypeAvg.mtx <- Matrix(as.matrix(do.call(cbind, collapsed.mtx.list)), sparse = TRUE) 
dim(collapsed.SubjectcohortCellTypeAvg.mtx) 

heatmap_metadata <- as.data.frame(cbind(colnames(collapsed.SubjectcohortCellTypeAvg.mtx), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 1), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 2), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 3)))
colnames(heatmap_metadata) <- c("cell.ident","cell.type.ident", "cohort.ident", "orig.ident")

heatmap_metadata$cell.type.ident <- factor(heatmap_metadata$cell.type.ident, levels=celltypes_to_plot)

cell_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(cell.ident)
cohort_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(cohort.ident)
celltype_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(cell.type.ident)
subject_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(orig.ident )

heatmap_df <-  as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[genes_for_heatmap,cell_order])
heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))

heatmap_cohort_annotation <- HeatmapAnnotation(cell_type=celltype_order, cohort=cohort_order, animal=subject_order, 
    col = list(cohort=mouse_cohort_colors, cell_type=celltype_colors, animal=mouse_subject_colors)) #show_legend = c("animal" = FALSE)

pdf("mouse_EC.Heatmap.Markers.AveragePerAnimal.pdf", height=10, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=celltype_order, column_title = NULL, use_raster=FALSE)
dev.off()


####  Full Mouse Dataset lineage marker genes per average per subject and per celltype
EC_vasc_genes_to_plot <- c("Flt1", "Clec14a", "Slco2a1", "Bmpr2", "Acvrl1", "Esam","Tmem100", "Gata2",  "Epas1")
EC_general_genes_to_plot <- c("Calcrl", "Erg", "Pecam1", "Cldn5","Cdh5","Tie1","Cav1", "Cavin2", "Ramp2", "Sash1")

genes_for_heatmap <- c("Chil3","Ear1",  # AM
	"Ms4a7", "C1qc",    # M     
	"Ms4a6c", "Ly6a2", # cMono
	"Cd300e", "Adgre4", # ncMono
	"Clec9a", "Xcr1", # cDC1
	"Clec10a", "Cd209a",    # cDC2 
	"Ccr7", "Ccl22", # DC_mature
	"Siglech", "Klk1", # pDC
	"Hdc", "Mcpt8", # Mast
	"Cd4", "Il7r",  # T_helper
	"Cd8a", "Cd8b1",  # T_cytotox
	"Ctla4","Tnfrsf4",   # T_reg
	"Klrb1c", "Nkg7",  # NK
	"Il23r", "Cxcr6",  # ILC
	"Ms4a1","Cd79a",  # B
	"Iglv1", "Tnfrsf17",  # Plasma
	"Ager","Cldn18", # AT1
	"Sftpc","Sftpd", # AT2
	"Trp63","Krt17", # Basal
	"Scgb1a1","Muc5b", # Secretory	
	"Hydin","Cfap299", # Ciliated
	"Mfap5","Scara5", # Fibro
	"Itga8","Inmt", # Myofibro
	"Myh11","Tagln", # SMC
	"Cox4i2","Pdgfrb", # Pericyte
	"Msln","Upk3b", # Mesothel
	"Mmrn1", "Reln", "Ccl21a", "Prox1", "Tbx1", "Klhl4","Scn3a", "Sema3d", "Pdpn", "Flt4", # EC_lymph	
	sample(EC_general_genes_to_plot),# EC_vasc+EC_lymph
	sample(EC_vasc_genes_to_plot))	# only EC_vasc

celltypes_to_plot <- c("AM", "M", "cMono", "ncMono", "cDC1", "cDC2", "DC_mature",  "pDC", "Mast",
	"T_helper", "T_cytotox", "T_reg", "NK","ILC", "B", "Plasma", 
	"AT1", "AT2", "Basal",  "Secretory", "Ciliated",
	"Fibroblast_adventitial", "Fibroblast_alveolar", "SMC","Pericyte","Mesothel",
	"EC_lymph","EC_arterial",  "EC_capillary_A", "EC_capillary_B", "EC_venous","EC_bronchial")

### with multithreading, first make all possible combinations
mouse$cellBarcode <- colnames(mouse) 
cellTypes <- levels(as.factor(mouse$cell.type.ident)) 
cellTypes <- cellTypes[cellTypes %in% celltypes_to_plot]
meta.data.sub <- mouse@meta.data[,c("cohort.ident", "cell.type.ident", "orig.ident", "cellBarcode")]

get.CT.DS.subj.vector <- function(cellTypes){
	tmp.meta.data <- meta.data.sub %>% filter(cell.type.ident== cellTypes)
	cohorts <- unique(tmp.meta.data$cohort.ident)
	subjects <- unique(tmp.meta.data$orig.ident)
	tmp.CT.DS.subj.vector <- vector()
	for(j in 1:length(cohorts)){ 
		for(k in 1:length(subjects)){
		temp.cells <- tmp.meta.data %>% filter(cohort.ident==cohorts[j] & orig.ident==subjects[k]) %>% pull(cellBarcode)
		if ( length(temp.cells) >1 ) { # >=
			tmp.CT.DS.subj.vector <- c(tmp.CT.DS.subj.vector, paste(cellTypes, cohorts[j], subjects[k], sep="__"))
			}
		}
	}
	cat("Completed for ", cellTypes, ".\n", sep="")
	return(tmp.CT.DS.subj.vector)
}
celltype_cohort_subject.list <- parallel::mclapply(cellTypes, get.CT.DS.subj.vector, mc.cores=8)
celltype_cohort_subject <- unlist(celltype_cohort_subject.list)

DefaultAssay(mouse) <- "RNA"
get.SubjectcohortCellTypeAvg <- function(celltype_cohort_subject){
	temp.cell.type <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][1]
	temp.cohort <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][2]
	temp.subject <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][3]
	temp.meta.data <- mouse@meta.data[,c("cohort.ident", "cell.type.ident", "orig.ident", "cellBarcode")]
	temp.cells <- temp.meta.data %>% filter(cell.type.ident==temp.cell.type & cohort.ident==temp.cohort & 
		orig.ident==temp.subject) %>% pull(cellBarcode) 
	if (length(temp.cells) > 1) { 
		tmp.df <- as.data.frame(rowMeans(GetAssayData(mouse)[,temp.cells]))	   
    } else { 
		tmp.df <- as.data.frame(GetAssayData(mouse)[,temp.cells])
		cat("Subject",temp.subject,"only has 1",temp.cell.type,"cell, using singlet for",temp.cohort,"representation...\n",sep=" ") 
    } 
	colnames(tmp.df) <- paste(celltype_cohort_subject)
	return(tmp.df)
}

collapsed.mtx.list <- parallel::mclapply(celltype_cohort_subject, get.SubjectcohortCellTypeAvg, mc.cores=8)
collapsed.SubjectcohortCellTypeAvg.mtx <- Matrix(as.matrix(do.call(cbind, collapsed.mtx.list)), sparse = TRUE) 
dim(collapsed.SubjectcohortCellTypeAvg.mtx) 

heatmap_metadata <- as.data.frame(cbind(colnames(collapsed.SubjectcohortCellTypeAvg.mtx), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 1), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 2), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 3)))
colnames(heatmap_metadata) <- c("cell.ident","cell.type.ident", "cohort.ident", "orig.ident")

heatmap_metadata$cell.type.ident <- factor(heatmap_metadata$cell.type.ident, levels=celltypes_to_plot)

cell_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(cell.ident)
cohort_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(cohort.ident)
celltype_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(cell.type.ident)
subject_order <- heatmap_metadata %>% arrange(cell.type.ident, cohort.ident, orig.ident ) %>% pull(orig.ident )

heatmap_df <-  as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[genes_for_heatmap,cell_order])
heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))

heatmap_cohort_annotation <- HeatmapAnnotation(cell_type=celltype_order, cohort=cohort_order, animal=subject_order, 
    col = list(cohort=mouse_cohort_colors, cell_type=celltype_colors, animal=mouse_subject_colors),
	show_legend = c("animal" = FALSE))
	
pdf("mouse.Heatmap.Markers.AveragePerSubject.pdf", height=10, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=celltype_order, column_title = NULL, use_raster=FALSE)
dev.off()

####
### Violin Plots EC adventitial in mice
genes_to_plot <- c("Col15a1", "Vwa1", "Ebf1",
	"Meox2", "Tshz2",  "Sparcl1",  
	"Nrp2","Hspg2","Magi1",
	"Plvap", "Igfbp7", "Ackr1")
	
mouse_EC$cell.type.ident <- 	factor(mouse_EC$cell.type.ident, levels=c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))
	
# without points
figure_list <- list()
for (i in 1:length(genes_to_plot)){
figure_list[[i]] <- VlnPlot(mouse_EC, genes_to_plot[i], group.by="cell.type.ident", pt.size = 0, col=celltype_colors)
figure_list[[i]] <- figure_list[[i]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
		axis.text.x=element_blank(), axis.ticks.x=element_blank()) + NoLegend()
}
# get common legend
figure_list[[length(figure_list)+1]] <- cowplot::get_legend(VlnPlot(mouse_EC, genes_to_plot[i], pt.size = 0, group.by="cell.type.ident",  col=celltype_colors))
mega_vln <- cowplot::plot_grid(plotlist=figure_list, ncol =3 ) #labels = LETTERS[1:(length(figure_list)-1)]
pdf("mouse_EC.mega_vln_bronchialEC.pdf", compress = FALSE, useDingbats = FALSE, height = 5, width = 9)
mega_vln
dev.off()


###################
# Sam's Connectome
library(connectome)
load("connectome.Robj")

connectome.filtered <- FilterConnectome(connectome.NicheNet,  min.pct =0.20) 
# remove self mappings & integrins
connectome.filtered <- connectome.filtered[connectome.filtered$ligand != connectome.filtered$receptor,]

integrins <- grep( "^ITG", unique(connectome.filtered$receptor), value=TRUE)

pdf(paste0("CircosPlot.ECs.Senders.05_20_20.pdf"))	
	connectome_to_plot <- connectome.filtered %>% 
		filter(source %in% c("EC_arterial", "EC_capillary_A", "EC_capillary_B", "EC_venous", "EC_bronchial", "EC_lymph")) %>% 
		filter(! target %in% c("EC_arterial", "EC_capillary_A", "EC_capillary_B", "EC_venous", "EC_bronchial", "EC_lymph")) %>% 
		filter(recept.scale > 0 &  ligand.scale > 0) %>%
		filter(! receptor %in% integrins) %>%
		filter(p_val_adj.lig < 1e-5 & p_val_adj.rec < 1e-5)  %>% 
		top_n(n=75, wt=weight_sc) 
	CircosPlot(connectome_to_plot, cols.use = celltype_colors)
dev.off()

pdf(paste0("CircosPlot.ECs.Receivers.03_09_21.pdf"))	
	connectome_to_plot <- connectome.filtered %>% 
		filter(target %in% c("EC_arterial", "EC_capillary_A", "EC_capillary_B", "EC_venous", "EC_bronchial", "EC_lymph")) %>% 
		filter(! source %in% c("EC_arterial", "EC_capillary_A", "EC_capillary_B", "EC_venous", "EC_bronchial", "EC_lymph")) %>% 
		filter(recept.scale > 0 &  ligand.scale > 0) %>%
		filter(! receptor %in% integrins) %>%
		filter(p_val_adj.lig < 1e-5 & p_val_adj.rec < 1e-5)  %>% 
		top_n(n=75, wt=weight_sc) 
	CircosPlot(connectome_to_plot, cols.use = celltype_colors)
dev.off()


########################################################################################################
# mouse human comparison
###

# input is the metadata slot which is going to get sampled from, the metadata slot used for balancing, and the seurat object,
# numberToPull number of cell ids to extract, identAttribute is a value in toSampleFrom
myBalancedRandomSampleling <- function(numberToPull, identAttribute,toSampleFrom, toBalance, seurat.object){
filtered.cell.ids <- rownames(seurat.object@meta.data[seurat.object@meta.data[,toSampleFrom]==identAttribute,])
sorted.count.table <- seurat.object@meta.data[filtered.cell.ids,] %>% count_(toBalance) %>% arrange(n) # need to add "_" to count so it uses the standard evalution
output.vector <- vector()
number.of.cells.left.to.distribute <- numberToPull
for (i in 1:nrow(sorted.count.table)){
	# if number of cells is equal or smaller than average to be sampled from, pick all
	temp.Balance.Attribute <- as.character(sorted.count.table[i,toBalance])
	temp.filtered.cell.ids <- rownames(seurat.object@meta.data[seurat.object@meta.data[,toBalance]== temp.Balance.Attribute,])[
			rownames(seurat.object@meta.data[seurat.object@meta.data[,toBalance]== temp.Balance.Attribute,]) %in% filtered.cell.ids]
	if ( sorted.count.table$n[i]   <= number.of.cells.left.to.distribute/(nrow(sorted.count.table)-i+1) ) {
		temp.output <- temp.filtered.cell.ids
		number.of.cells.left.to.distribute <- number.of.cells.left.to.distribute - length(temp.output)
		}
	# else randomly select the number of cells left to distribute	
	else{
		set.seed(7)
		temp.output <- sample(x=temp.filtered.cell.ids, size=floor(number.of.cells.left.to.distribute/(nrow(sorted.count.table)-i+1)))
		number.of.cells.left.to.distribute <- number.of.cells.left.to.distribute - length(temp.output)
		}
	output.vector <- c(output.vector, temp.output)
	}
return(output.vector)
}

table(mouse_EC$cell.type.ident)

# downsample to smallest population (number to DownSample) in mouse and human
nDS <- 500

mCapA <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_capillary_A",toSampleFrom="cell.type.ident", toBalance="orig.ident", seurat.object=mouse_EC)
mCapB <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_capillary_B",toSampleFrom="cell.type.ident", toBalance="orig.ident", seurat.object=mouse_EC)
mArt <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_arterial",toSampleFrom="cell.type.ident", toBalance="orig.ident", seurat.object=mouse_EC)
mVen <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_venous",toSampleFrom="cell.type.ident", toBalance="orig.ident", seurat.object=mouse_EC)

hCapA <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_capillary_A",toSampleFrom="cell.type.ident", toBalance="subject.ident", seurat.object=EC_Subset_Object)
hCapB <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_capillary_B",toSampleFrom="cell.type.ident", toBalance="subject.ident", seurat.object=EC_Subset_Object)
hArt <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_arterial",toSampleFrom="cell.type.ident", toBalance="subject.ident", seurat.object=EC_Subset_Object)
hVen <- myBalancedRandomSampleling(numberToPull=nDS, identAttribute="EC_venous",toSampleFrom="cell.type.ident", toBalance="subject.ident", seurat.object=EC_Subset_Object)

mouseEC_sub <- subset(mouse_EC, cells=c(mCapA,mCapB,mVen,mArt)) 
humanEC_sub <- subset(EC_Subset_Object, cells=c(hCapA,hCapB,hVen,hArt))   

# function, see above
mouse_ECsub_markers <- myMultithreadFindMarkersDOR(seurat.object=mouseEC_sub, group.name="cell.type.ident", logFC.min.filter=0.2, ncores=8)
human_ECsub_markers <- myMultithreadFindMarkersDOR(seurat.object=humanEC_sub, group.name="cell.type.ident", logFC.min.filter=0.2, ncores=8)

###########
# translate human gene to mouse and vice versa
library(biomaRt)
# GRCh38.p12 (this is what i mapped against) is Ensemble97 (July 2019)
# listEnsemblArchives()
# http://jul2019.archive.ensembl.org
# listDatasets(useMart("ensembl", host="http://jul2019.archive.ensembl.org"))
human_mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="http://jul2019.archive.ensembl.org")
mouse_mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="http://jul2019.archive.ensembl.org")

translator <- getLDS(attributes = c("ensembl_gene_id", "external_gene_name", "mgi_symbol"), mart = mouse_mart, 
	attributesL = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"), martL = human_mart, uniqueRows=T)
colnames(translator)[colnames(translator) %in% c("Gene.name", "Gene.name.1")] <- c("Gene.name.mouse", "Gene.name.human")
translator <- translator[,c("Gene.name.mouse", "Gene.name.human")]
translator <- distinct(translator)
head(translator)

# only translate if there is exactly one unique gene in the species to be translated
duplicated_mouse_genes <- unique(translator$Gene.name.mouse[duplicated(translator$Gene.name.mouse)])
duplicated_human_genes <- unique(translator$Gene.name.human[duplicated(translator$Gene.name.human)])
translator <- translator %>% filter(! Gene.name.mouse %in% duplicated_mouse_genes & ! Gene.name.human %in% duplicated_human_genes) 

## translate!
colnames(mouse_ECsub_markers)[colnames(mouse_ECsub_markers)=="gene"] <- "Gene.name.mouse"
colnames(human_ECsub_markers)[colnames(human_ECsub_markers)=="gene"] <- "Gene.name.human"

mouse_ECsub_markers <- mouse_ECsub_markers %>% plyr::join(x=., y=translator, by="Gene.name.mouse") 
human_ECsub_markers <- human_ECsub_markers %>% plyr::join(x=., y=translator, by="Gene.name.human") 

mouse_ECsub_markers$origin <- "mouse"
human_ECsub_markers$origin <- "human"

markers_ECsub <- bind_rows(human_ECsub_markers, mouse_ECsub_markers)
markers_ECsub$Gene.name.merged <- paste(markers_ECsub$Gene.name.human , markers_ECsub$Gene.name.mouse, sep="_")

# select specific marker genes per cell type
markers_ECsub <- markers_ECsub %>% filter(p_val_adj < 0.05)

Art_genes_human <- markers_ECsub %>% filter(cell.type.ident == "EC_arterial" & origin=="human") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
Art_genes_mouse <- markers_ECsub %>% filter(cell.type.ident == "EC_arterial" & origin=="mouse") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
common_art_genes <- intersect(Art_genes_human,Art_genes_mouse)

capA_genes_human <- markers_ECsub %>% filter(cell.type.ident == "EC_capillary_A" & origin=="human") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
capA_genes_mouse <- markers_ECsub %>% filter(cell.type.ident == "EC_capillary_A" & origin=="mouse") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
common_capA_genes <- intersect(capA_genes_human,capA_genes_mouse)

capB_genes_human <- markers_ECsub %>% filter(cell.type.ident == "EC_capillary_B" & origin=="human") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
capB_genes_mouse <- markers_ECsub %>% filter(cell.type.ident == "EC_capillary_B" & origin=="mouse") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
common_capB_genes <- intersect(capB_genes_human,capB_genes_mouse)

Ven_genes_human <- markers_ECsub %>% filter(cell.type.ident == "EC_venous" & origin=="human") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
Ven_genes_mouse <- markers_ECsub %>% filter(cell.type.ident == "EC_venous" & origin=="mouse") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
common_ven_genes <- intersect(Ven_genes_human,Ven_genes_mouse)

# remove non-specific genes
common_art_genes <- common_art_genes[! common_art_genes %in% c(capA_genes_human, capA_genes_mouse,
	capB_genes_human, capB_genes_mouse, Ven_genes_human, Ven_genes_mouse)]
common_capA_genes <- common_capA_genes[! common_capA_genes %in% c(Art_genes_human, Art_genes_mouse, 
	capB_genes_human, capB_genes_mouse, Ven_genes_human, Ven_genes_mouse)]
common_capB_genes <- common_capB_genes[! common_capB_genes %in% c(Art_genes_human, Art_genes_mouse, capA_genes_human, capA_genes_mouse,
	 Ven_genes_human, Ven_genes_mouse)]
common_ven_genes <- common_ven_genes[! common_ven_genes %in% c(Art_genes_human, Art_genes_mouse, capA_genes_human, capA_genes_mouse,
	capB_genes_human, capB_genes_mouse)]

genes_to_plot <- c(common_art_genes, common_capA_genes, common_capB_genes, common_ven_genes) #common_bronch_genes

# extract average and percent expressed data frames
DefaultAssay(FullDataset_Object) <- "RNA"
avg.exp.human.df <- data.frame(
	"EC_arterial"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident=="EC_arterial")]),
	"EC_capillary_A"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident=="EC_capillary_A")]),
	"EC_capillary_B"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident=="EC_capillary_B")]),
	"EC_venous"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident=="EC_venous")]),
	"EC_bronchial"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident=="EC_bronchial")]),
	"EC_lymph"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident=="EC_lymph")]),
	"Epithelial"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident %in% 
		c("AT1","Basal", "AT2", "Ciliated", "Club" , "Goblet", "Ionocyte", "PNEC", "Suprabasal"))]),
	"Stromal"=rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident %in% 
		c("Fibroblast_adventitial", "Mesothel","Fibroblast_alveolar", "Pericyte", "SMC"))]),
	"Myeloid" = rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident %in% 
		c("AM", "cMono", "ncMono", "cDC2", "cDC1", "DC_Langerhans", "DC_mature", "M", "Mast","pDC"))]),
	"Lymphoid" =rowMeans(GetAssayData(FullDataset_Object, slot="data")[,WhichCells(FullDataset_Object,expression=cell.type.ident %in% 
		c("T_helper", "B", "ILC","NK", "Plasma", "T_cytotox", "T_reg"))])
)

avg.exp.mouse.df <- data.frame(
	"EC_arterial"=rowMeans(GetAssayData(mouse_EC, slot="data")[,WhichCells(mouse_EC,expression=cell.type.ident=="EC_arterial")]),
	"EC_capillary_A"=rowMeans(GetAssayData(mouse_EC, slot="data")[,WhichCells(mouse_EC,expression=cell.type.ident=="EC_capillary_A")]),
	"EC_capillary_B"=rowMeans(GetAssayData(mouse_EC, slot="data")[,WhichCells(mouse_EC,expression=cell.type.ident=="EC_capillary_B")]),
	"EC_venous"=rowMeans(GetAssayData(mouse_EC, slot="data")[,WhichCells(mouse_EC,expression=cell.type.ident=="EC_venous")]),
	"EC_bronchial"=rowMeans(GetAssayData(mouse_EC, slot="data")[,WhichCells(mouse_EC,expression=cell.type.ident=="EC_bronchial")]),
	"EC_lymph"=rowMeans(GetAssayData(mouse, slot="data")[,WhichCells(mouse,expression=cell.type.ident=="EC_lymph")]),
	"Epithelial"=rowMeans(GetAssayData(mouse, slot="data")[,WhichCells(mouse,expression=lineage.ident=="epithelial")]),
	"Stromal"=rowMeans(GetAssayData(mouse, slot="data")[,WhichCells(mouse,expression=lineage.ident=="stromal")]),
	"Myeloid" = rowMeans(GetAssayData(mouse, slot="data")[,WhichCells(mouse,expression=lineage.ident=="myeloid")]),
	"Lymphoid" =rowMeans(GetAssayData(mouse, slot="data")[,WhichCells(mouse,expression=lineage.ident=="lymphoid")])
)

# subset to genes in translator, move rownames
avg.exp.human.df <- avg.exp.human.df %>% .[translator$Gene.name.human,] %>% tibble::rownames_to_column("Gene.name.human") %>%
	 inner_join(x=., y=translator, by="Gene.name.human")  
avg.exp.mouse.df <- avg.exp.mouse.df %>% .[translator$Gene.name.mouse,] %>% tibble::rownames_to_column("Gene.name.mouse") %>%
	 inner_join(x=., y=translator, by="Gene.name.mouse")  
	 		 	 
avg.exp.human.df$Gene.name.merged <- paste(avg.exp.human.df$Gene.name.human , avg.exp.human.df$Gene.name.mouse, sep="_")	 
avg.exp.mouse.df$Gene.name.merged <- paste(avg.exp.mouse.df$Gene.name.human , avg.exp.mouse.df$Gene.name.mouse, sep="_")	 
	 
# subset stuff to plot, then scale 
avg.exp.mouse.df.EC <- avg.exp.mouse.df %>% tibble::column_to_rownames("Gene.name.merged") %>%	
	.[genes_to_plot, c("EC_arterial", "EC_capillary_A", "EC_capillary_B", "EC_venous")]  #, "EC_bronchial"
avg.exp.human.df.EC <- avg.exp.human.df %>% tibble::column_to_rownames("Gene.name.merged") %>%	
	.[genes_to_plot, c("EC_arterial", "EC_capillary_A", "EC_capillary_B", "EC_venous")]  #, "EC_bronchial"

#### heatmaps unity normalized:
avg.exp.human.df.EC.normalized <- apply(avg.exp.human.df.EC , MARGIN=1, FUN= myUnityNormalize) 
avg.exp.mouse.df.EC.normalized <- apply(avg.exp.mouse.df.EC , MARGIN=1, FUN= myUnityNormalize) 
rownames(avg.exp.human.df.EC.normalized) <- paste(rownames(avg.exp.human.df.EC.normalized), "human", sep="_")
rownames(avg.exp.mouse.df.EC.normalized) <- paste(rownames(avg.exp.mouse.df.EC.normalized), "mouse", sep="_")

normalized.merged.ECsub <- rbind(avg.exp.human.df.EC.normalized, avg.exp.mouse.df.EC.normalized)
# sort by average unity normalized scpre, groupbed by cell type
sort_scores_ECsub <- c(sort(apply(normalized.merged.ECsub[,common_art_genes], 2, function(x) median(x)), decreasing=FALSE),
	sort(apply(normalized.merged.ECsub[,common_capA_genes], 2, function(x) median(x)), decreasing=FALSE),
	sort(apply(normalized.merged.ECsub[,common_capB_genes], 2, function(x) median(x)), decreasing=FALSE),
	sort(apply(normalized.merged.ECsub[,common_ven_genes], 2, function(x) median(x)), decreasing=FALSE))


normalized.merged.ECsub <- normalized.merged.ECsub[c("EC_arterial_human", "EC_arterial_mouse","EC_capillary_A_human", "EC_capillary_A_mouse",
	"EC_capillary_B_human", "EC_capillary_B_mouse","EC_venous_human", "EC_venous_mouse")

# change orientation of heatmap
normalized.merged.ECsub <- t(normalized.merged.ECsub)

pdf("Species.comparison.merged.EC.black.UnityNomalized.08_06_20.pdf", width=3, height=6)
Heatmap(normalized.merged.ECsub, col = viridis(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
    use_raster=FALSE, column_title = "UnityNomalized", column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6))
dev.off()


genes_to_show <- c("DKK2_Dkk2", "GJA5_Gja5", "CXCL12_Cxcl12", "BMX_Bmx", "EFNB2_Efnb2", "FBLN2_Fbln2", "GJA4_Gja4", "DLL4_Dll4",  # arterial
	"EDNRB_Ednrb", "RCSD1_Rcsd1", "TBX2_Tbx2", "PRX_Prx", "CAVIN2_Cavin2", "ITGA3_Itga3", "KDR_Kdr","TBX3_Tbx3","CAV2_Cav2","PXDN_Pxdn", # capA
	"GPIHBP1_Gpihbp1", "EPAS1_Epas1",  "ZNF608_Zfp608",  "ADGRE5_Adgre5", # capB
	"CPE_Cpe",  "HDAC9_Hdac9", "PTGS1_Ptgs1",  "ACKR3_Ackr3", "NR2F2_Nr2f2", "RGS5_Rgs5", "PRCP_Prcp") #venous
# sort them
genes_to_show <- names(sort_scores_ECsub)[names(sort_scores_ECsub) %in% genes_to_show ]
	
row_annotation <- rowAnnotation(markers = anno_mark(at = which(names(sort_scores_ECsub)  %in% genes_to_show ), labels=genes_to_show,
	labels_gp = gpar(fontsize = 6)))

pdf("Species.comparison.merged.EC.black.UnityNomalized.some.labeled.08_06_20.pdf", width=3, height=6)
Heatmap(normalized.merged.ECsub, col = viridis(256), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    use_raster=FALSE, column_title = "UnityNomalized",right_annotation =row_annotation, column_names_gp = gpar(fontsize = 6))
dev.off()



###################################
## now for general EC and vascular EC marker

# first downsample to get similar representation
myBalancedRandomSamplelingBig <- function(numberToPullperLevel, toSampleFrom, toBalance, seurat.object){
identAttribute <- unique(seurat.object@meta.data[,toSampleFrom])
output.vector <- list()
for (j in 1:length(identAttribute)){
	filtered.cell.ids <- rownames(seurat.object@meta.data[seurat.object@meta.data[,toSampleFrom]==identAttribute[j],])
	sorted.count.table <- seurat.object@meta.data[filtered.cell.ids,] %>% count_(toBalance) %>% arrange(n) # need to add "_" to count so it uses the standard evalution
	output.vector[[j]] <- vector()
	number.of.cells.left.to.distribute <- numberToPullperLevel
	for (i in 1:nrow(sorted.count.table)){
		# if number of cells is equal or smaller than average to be sampled from, pick all
		temp.Balance.Attribute <- as.character(sorted.count.table[i,toBalance])
		temp.filtered.cell.ids <- rownames(seurat.object@meta.data[seurat.object@meta.data[,toBalance]== temp.Balance.Attribute,])[
				rownames(seurat.object@meta.data[seurat.object@meta.data[,toBalance]== temp.Balance.Attribute,]) %in% filtered.cell.ids]
		if ( sorted.count.table$n[i]   <= number.of.cells.left.to.distribute/(nrow(sorted.count.table)-i+1) ) {
			temp.output <- temp.filtered.cell.ids
			number.of.cells.left.to.distribute <- number.of.cells.left.to.distribute - length(temp.output)
			}
		# else randomly select the number of cells left to distribute	
		else{
			set.seed(7)
			temp.output <- sample(x=temp.filtered.cell.ids, size=floor(number.of.cells.left.to.distribute/(nrow(sorted.count.table)-i+1)))
			number.of.cells.left.to.distribute <- number.of.cells.left.to.distribute - length(temp.output)
			}
		output.vector[[j]] <- c(output.vector[[j]], temp.output)
		}
	}
output.vector <- do.call(c, output.vector)
return(output.vector)
}

mouse_cells <- myBalancedRandomSamplelingBig(numberToPullperLevel=500, toSampleFrom="cell.type.ident", toBalance="orig.ident", seurat.object=mouse)
length(mouse_cells) #13461

mouse_sub <- subset(mouse, cells=mouse_cells)
table(mouse_sub$orig.ident, mouse_sub$cell.type.ident)

human_cells <- myBalancedRandomSamplelingBig(numberToPullperLevel=500, toSampleFrom="cell.type.ident", toBalance="subject.ident", seurat.object=FullDataset_Object)
length(human_cells) #17149

human_sub <- subset(FullDataset_Object, cells=human_cells)
table(human_sub$subject.ident, human_sub$cell.type.ident)

# function, see above
mouse_markers_all <- myMultithreadFindMarkersDOR(seurat.object=mouse_sub, group.name="cell.type.ident", logFC.min.filter=0.2, ncores=8)
human_markers_all <- myMultithreadFindMarkersDOR(seurat.object=human_sub, group.name="cell.type.ident", logFC.min.filter=0.2, ncores=8)


###########
# translate human gene to mouse and vice versa
colnames(mouse_markers_all)[colnames(mouse_markers_all)=="gene"] <- "Gene.name.mouse"
colnames(human_markers_all)[colnames(human_markers_all)=="gene"] <- "Gene.name.human"

mouse_markers_all <- mouse_markers_all %>% plyr::join(x=., y=translator, by="Gene.name.mouse") 
human_markers_all <- human_markers_all %>% plyr::join(x=., y=translator, by="Gene.name.human") 

mouse_markers_all$origin <- "mouse"
human_markers_all$origin <- "human"

markers_all <- bind_rows(human_markers_all, mouse_markers_all)
markers_all$Gene.name.merged <- paste(markers_all$Gene.name.human , markers_all$Gene.name.mouse, sep="_")

# select specific marker genes per cell type
markers_all <- markers_all %>% filter(p_val_adj < 0.05)

Lymph_genes_human <- markers_all %>% filter(cell.type.ident == "EC_lymph" & origin=="human") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
Lymph_genes_mouse <- markers_all %>% filter(cell.type.ident == "EC_lymph" & origin=="mouse") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged) 
common_lymph_genes <- intersect(Lymph_genes_human,Lymph_genes_mouse)

human_venous <- markers_all %>% filter(origin=="human" & cell.type.ident=="EC_venous") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
human_capillary_A <- markers_all %>% filter(origin=="human" & cell.type.ident=="EC_capillary_A") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
human_capillary_B <- markers_all %>% filter(origin=="human" & cell.type.ident=="EC_capillary_B") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
human_arterial <- markers_all %>% filter(origin=="human" & cell.type.ident=="EC_arterial") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
human_lymph <- markers_all %>% filter(origin=="human" & cell.type.ident=="EC_lymph") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)

mouse_venous <- markers_all %>% filter(origin=="mouse" & cell.type.ident=="EC_venous") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
mouse_capillary_A <- markers_all %>% filter(origin=="mouse" & cell.type.ident=="EC_capillary_A") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
mouse_capillary_B <- markers_all %>% filter(origin=="mouse" & cell.type.ident=="EC_capillary_B") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
mouse_arterial <- markers_all %>% filter(origin=="mouse" & cell.type.ident=="EC_arterial") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)
mouse_lymph <- markers_all %>% filter(origin=="mouse" & cell.type.ident=="EC_lymph") %>% drop_na(Gene.name.human,Gene.name.mouse) %>% pull(Gene.name.merged)

human_EC_general_markers <- Reduce(intersect, list(human_venous, human_capillary_A, human_capillary_B, human_arterial, human_lymph))
mouse_EC_general_markers <- Reduce(intersect, list(mouse_venous, mouse_capillary_A, mouse_capillary_B, mouse_arterial, mouse_lymph))
common_EC_general_genes <- intersect(human_EC_general_markers,mouse_EC_general_markers)
		
human_EC_vasc_markers <- Reduce(intersect, list(human_venous, human_capillary_A, human_capillary_B, human_arterial))
human_EC_vasc_markers <- human_EC_vasc_markers[! human_EC_vasc_markers %in% Lymph_genes_human]
mouse_EC_vasc_markers <- Reduce(intersect, list(mouse_venous, mouse_capillary_A, mouse_capillary_B, mouse_arterial))
mouse_EC_vasc_markers <- mouse_EC_vasc_markers[! mouse_EC_vasc_markers %in% Lymph_genes_mouse]
common_EC_vasc_genes <- intersect(mouse_EC_vasc_markers,human_EC_vasc_markers)

common_lymph_genes <- common_lymph_genes[! common_lymph_genes %in% c(human_venous, human_capillary_A, human_capillary_B, human_arterial,
	mouse_venous, mouse_capillary_A, mouse_capillary_B, mouse_arterial)]

genes_to_plot <- c(common_lymph_genes, common_EC_general_genes, common_EC_vasc_genes)

# subset stuff to plot, then scale 
avg.exp.mouse.df.bigEC <- avg.exp.mouse.df %>% tibble::column_to_rownames("Gene.name.merged") %>% dplyr::select(.,-"EC_bronchial") %>% .[genes_to_plot, ] 
avg.exp.human.df.bigEC <- avg.exp.human.df %>% tibble::column_to_rownames("Gene.name.merged") %>% dplyr::select(.,-"EC_bronchial") %>% .[genes_to_plot, ]

avg.exp.human.df.bigEC <- avg.exp.human.df.bigEC[, !colnames(avg.exp.human.df.bigEC) %in% c("Gene.name.human","Gene.name.mouse")]
avg.exp.mouse.df.bigEC <- avg.exp.mouse.df.bigEC[, !colnames(avg.exp.mouse.df.bigEC) %in% c("Gene.name.human","Gene.name.mouse")]

#### heatmaps unity normalized:
avg.exp.human.df.bigEC.normalized <- apply(avg.exp.human.df.bigEC , MARGIN=1, FUN= myUnityNormalize) 
avg.exp.mouse.df.bigEC.normalized <- apply(avg.exp.mouse.df.bigEC , MARGIN=1, FUN= myUnityNormalize) 
rownames(avg.exp.human.df.bigEC.normalized) <- paste(rownames(avg.exp.human.df.bigEC.normalized), "human", sep="_")
rownames(avg.exp.mouse.df.bigEC.normalized) <- paste(rownames(avg.exp.mouse.df.bigEC.normalized), "mouse", sep="_")
merged.bigEC.normalized <- rbind(avg.exp.human.df.bigEC.normalized, avg.exp.mouse.df.bigEC.normalized)	

# sort by average unity normalized scpre, groupbed by cell type
sort_scores_ECbig <- c(
	sort(apply(merged.bigEC.normalized[c("Stromal_human", "Stromal_mouse", "Epithelial_human", "Epithelial_mouse", "Myeloid_human", "Myeloid_mouse", 
		"Lymphoid_human","Lymphoid_mouse"),common_EC_general_genes], 2, function(x) median(x)), decreasing=FALSE),
	sort(apply(merged.bigEC.normalized[c("Stromal_human", "Stromal_mouse", "Epithelial_human", "Epithelial_mouse", "Myeloid_human", "Myeloid_mouse", 
		"Lymphoid_human","Lymphoid_mouse"),common_EC_vasc_genes], 2, function(x) median(x)), decreasing=FALSE),
	sort(apply(merged.bigEC.normalized[c("Stromal_human", "Stromal_mouse", "Epithelial_human", "Epithelial_mouse", "Myeloid_human", "Myeloid_mouse", 
		"Lymphoid_human","Lymphoid_mouse"),common_lymph_genes], 2, function(x) median(x)), decreasing=FALSE))

merged.bigEC.normalized <- merged.bigEC.normalized[ c("EC_arterial_human", "EC_arterial_mouse","EC_capillary_A_human", "EC_capillary_A_mouse",
	"EC_capillary_B_human", "EC_capillary_B_mouse","EC_venous_human", "EC_venous_mouse", "EC_lymph_human", "EC_lymph_mouse",
	"Stromal_human", "Stromal_mouse", "Epithelial_human", "Epithelial_mouse", "Myeloid_human", "Myeloid_mouse", "Lymphoid_human","Lymphoid_mouse"),
	names(sort_scores_ECbig)] #   

# change orientation of heatmap
merged.bigEC.normalized <- t(merged.bigEC.normalized)

genes_to_show <- c("CALCRL_Calcrl","ERG_Erg","PECAM1_Pecam1", "CLDN5_Cldn5", "CDH5_Cdh5", "TIE1_Tie1", "CAV1_Cav1",
	"CAVIN2_Cavin2", "RAMP2_Ramp2", "SASH1_Sash1", #EC_general
	"FLT1_Flt1",  "CLEC14A_Clec14a", "SLCO2A1_Slco2a1", "BMPR2_Bmpr2", "ACVRL1_Acvrl1", "ESAM_Esam","TMEM100_Tmem100",
	"GATA2_Gata2", "EPAS1_Epas1", #EC_vasc
	"SEMA3D_Sema3d","PROX1_Prox1", "PDPN_Pdpn",  "RELN_Reln","TBX1_Tbx1","KLHL4_Klhl4", "SCN3A_Scn3a")    #EC_lymph
genes_to_show <- names(sort_scores_ECbig)[names(sort_scores_ECbig) %in% genes_to_show ]	

row_annotation <- rowAnnotation(markers = anno_mark(at = which(names(sort_scores_ECbig) %in% genes_to_show ), labels=genes_to_show, 
	labels_gp = gpar(fontsize = 6)))

pdf("Species.comparison.merged.ECbig.black.UnityNomalized.08_06_20.pdf", width=4, height=6)
Heatmap(merged.bigEC.normalized , col = viridis(256), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    use_raster=FALSE, column_title = "unity normalized", right_annotation = row_annotation, column_names_gp = gpar(fontsize = 6))
dev.off()


colnames(markers_all)[colnames(markers_all)=="origin"] <- "species"
colnames(markers_ECsub)[colnames(markers_ECsub)=="origin"] <- "species"
markers_all$dataset <- "Full Lung Dataset"
markers_ECsub$dataset <- "Lung vascular EC_subset"

markers_to_save <- rbind(markers_all, markers_ECsub)

markers_to_save$conserved.lymphatic.marker <- "no"
markers_to_save$conserved.pan.endothelial.marker <- "no"
markers_to_save$conserved.pan.vascular.marker <- "no"
markers_to_save$conserved.arterial.EC.marker <- "no"
markers_to_save$conserved.venous.EC.marker <- "no"
markers_to_save$conserved.capillaryA.EC.marker <- "no"
markers_to_save$conserved.capillaryB.EC.marker <- "no"

markers_to_save$conserved.lymphatic.marker[markers_to_save$Gene.name.merged %in% common_lymph_genes] <- "yes"
markers_to_save$conserved.pan.endothelial.marker[markers_to_save$Gene.name.merged %in% common_EC_general_genes] <- "yes"
markers_to_save$conserved.pan.vascular.marker[markers_to_save$Gene.name.merged %in% common_EC_vasc_genes] <- "yes"
markers_to_save$conserved.arterial.EC.marker[markers_to_save$Gene.name.merged %in% common_art_genes] <- "yes"
markers_to_save$conserved.venous.EC.marker[markers_to_save$Gene.name.merged %in% common_ven_genes] <- "yes"
markers_to_save$conserved.capillaryA.EC.marker[markers_to_save$Gene.name.merged %in% common_capA_genes] <- "yes"
markers_to_save$conserved.capillaryB.EC.marker[markers_to_save$Gene.name.merged %in% common_capB_genes] <- "yes"

markers_to_save$nCells.in <- NULL
markers_to_save$nCells.out <- NULL

write.table(markers_to_save, file="balanced_testing.conserved_marker_genes.08_18_2020.txt", row.names = FALSE, sep="\t")



#######
# get basic stats 

unique(FullDataset_Object$cohort.ident)
FullDataset_Object
unique(FullDataset_Object$subject.ident)

unique(EC_Subset_Object$cohort.ident)
EC_Subset_Object
unique(EC_Subset_Object$subject.ident)

EC_cancer_Merged_Object
unique(EC_cancer_Merged_Object$cohort.ident)
EC_cancer_Merged_Object
unique(EC_cancer_Merged_Object$subject.ident)

write.table(table(FullDataset_Object$subject.ident, FullDataset_Object$cell.type.ident), file="FullDataset_Object.subject.by.celltype.txt")
write.table(table(EC_Subset_Object$subject.ident, EC_Subset_Object$cell.type.ident), file="EC_Subset_Object.subject.by.celltype.txt")
write.table(table(EC_cancer_Merged_Object$subject.ident, EC_cancer_Merged_Object$cell.type.ident.cancer), file="EC_cancer_Merged_Object.subject.by.celltype.txt")

write.table(prop.table(table(FullDataset_Object$subject.ident, FullDataset_Object$cell.type.ident), margin =1)*100, 
	file="FullDataset_Object.subject.by.celltype.Percentages.txt")
write.table(prop.table(table(EC_Subset_Object$subject.ident, EC_Subset_Object$cell.type.ident), margin =1)*100, 
	file="EC_Subset_Object.subject.by.celltype.Percentages.txt")

	
	
write.table(table(mouse$orig.ident, mouse$cell.type.ident), file="mouse.animal.by.celltype.txt")
write.table(table(mouse_EC$orig.ident, mouse_EC$cell.type.ident), file="mouse_EC.animal.by.celltype.txt")
	
write.table(prop.table(table(mouse$orig.ident, mouse$cell.type.ident), margin =1)*100, 
	file="mouse.animal.by.celltype.Percentages.txt")
write.table(prop.table(table(mouse_EC$orig.ident, mouse_EC$cell.type.ident), margin =1)*100, 
	file="mouse_EC.animal.by.celltype.Percentages.txt")
	

prop.table(table(FullDataset_Object$cell.type.ident))*100
                    AM                      M                  cMono
           22.52483420             9.76680256            10.13752117
                ncMono                   cDC1                   cDC2
            3.58373288             0.13816715             2.76944389
             DC_mature          DC_Langerhans                    pDC
            0.10084408             0.34559731             0.11735236
                  Mast               T_helper              T_cytotox
            1.08380466             5.88125520             5.97599839
                 T_reg                     NK                    ILC
            0.35277483             2.30936522             0.38866240
                     B                 Plasma                    AT1
            1.16598720             0.44823577             2.37432172
                   AT2                  Basal             Suprabasal
           13.56011886             2.79061755             0.54226120
                Goblet                   Club               Ciliated
            0.96537567             0.65817806             4.09692515
              Ionocyte                   PNEC Fibroblast_adventitial
            0.02296805             0.01614941             0.91872183
   Fibroblast_alveolar                    SMC               Pericyte
            0.53005943             0.23362809             0.12632425
              Mesothel            EC_arterial         EC_capillary_A
            0.06998076             0.67791622             0.55805174
        EC_capillary_B              EC_venous           EC_bronchial
            1.19649163             0.56522925             0.75938101
              EC_lymph             Cell_Cycle
            0.75148575             1.49543510

			
prop.table(table(EC_Subset_Object$cell.type.ident))*100
   # EC_arterial   EC_bronchial EC_capillary_A EC_capillary_B      EC_venous
      # 17.23682       15.13010       15.30181       37.00304       15.32823


#### make boxplot of ratio cap a to cap b
# human
cell_df <- as.data.frame.matrix(table(EC_Subset_Object$subject.ident, EC_Subset_Object$cell.type.ident))
cell_df$ratio_CapBtoCapA <- cell_df$EC_capillary_B / cell_df$EC_capillary_A
cell_df$ratio_CapBtoCapA[cell_df$EC_capillary_A==0 | cell_df$EC_capillary_B==0] <- rep(NA)

cell_df$ratio_CapAtoCapB <- cell_df$EC_capillary_A / cell_df$EC_capillary_B  
cell_df$ratio_CapAtoCapB[cell_df$EC_capillary_A==0 | cell_df$EC_capillary_B==0] <- rep(NA)

summary(cell_df[!is.na(cell_df$ratio_CapBtoCapA),]$ratio_CapBtoCapA)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.01439  0.85972  2.06410  3.14017  3.87500 15.33333

summary(cell_df[!is.na(cell_df$ratio_CapAtoCapB),]$ratio_CapAtoCapB)
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
 # 0.06522  0.25893  0.48494  2.45079  1.16360 69.50000

 
prop_cell_df <- as.data.frame.matrix(prop.table(table(EC_Subset_Object$subject.ident, EC_Subset_Object$cell.type.ident), margin=1)*100)
summary(prop_cell_df)
  EC_arterial      EC_bronchial     EC_capillary_A   EC_capillary_B
 Min.   : 0.000   Min.   :  0.000   Min.   : 0.000   Min.   : 0.000
 1st Qu.: 9.441   1st Qu.:  0.000   1st Qu.: 0.000   1st Qu.: 3.125
 Median :18.797   Median :  2.370   Median : 9.524   Median :20.508
 Mean   :21.506   Mean   : 13.296   Mean   :15.261   Mean   :25.298
 3rd Qu.:30.861   3rd Qu.:  9.774   3rd Qu.:20.701   3rd Qu.:41.543
 Max.   :60.000   Max.   :100.000   Max.   :96.528   Max.   :78.431
   EC_venous
 Min.   :  0.000
 1st Qu.:  6.215
 Median : 21.212
 Mean   : 24.638
 3rd Qu.: 37.931
 Max.   :100.000

# make tidy
prop_cell_df <- gather(prop_cell_df,"cell.type.ident", "percentage")
# make factor to order
prop_cell_df$cell.type.ident <- factor(prop_cell_df$cell.type.ident, 
	levels=c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))

pdf("Boxplots.Percentages.vascECs.human.pdf", useDingbats=FALSE)	
prop_cell_df %>% ggplot(aes(y=percentage, x=cell.type.ident, fill=cell.type.ident)) +  
	geom_boxplot(outlier.shape = NA ) + scale_fill_manual(values=celltype_colors) + theme_classic() +
	scale_y_continuous(limits = c(0, 85))
dev.off()


# make same, but now only if at least 4 out of 5 EC subtypes were profiled and from distal lung
prop_cell_df <- as.data.frame.matrix(prop.table(table(EC_Subset_Object$subject.ident[EC_Subset_Object$location.ident=="Distal_Lung"],
	EC_Subset_Object$cell.type.ident[EC_Subset_Object$location.ident=="Distal_Lung"]), margin=1)*100)
prop_cell_df <-prop_cell_df[rowSums(prop_cell_df!=0) >=4,]
# n=48
summary(prop_cell_df)
  EC_arterial       EC_bronchial    EC_capillary_A  EC_capillary_B
 Min.   : 0.6944   Min.   : 0.000   Min.   : 0.00   Min.   : 0.00
 1st Qu.:12.7259   1st Qu.: 0.000   1st Qu.: 5.37   1st Qu.:10.84
 Median :20.0000   Median : 3.214   Median :13.56   Median :30.68
 Mean   :22.4300   Mean   : 6.878   Mean   :17.61   Mean   :30.79
 3rd Qu.:29.9374   3rd Qu.: 7.945   3rd Qu.:27.95   3rd Qu.:44.11
 Max.   :55.3846   Max.   :57.746   Max.   :96.53   Max.   :78.43
   EC_venous
 Min.   : 1.389
 1st Qu.:10.951
 Median :21.320
 Mean   :22.288
 3rd Qu.:32.333
 Max.   :55.556

# make tidy
prop_cell_df <- gather(prop_cell_df,"cell.type.ident", "percentage")
# make factor to order
prop_cell_df$cell.type.ident <- factor(prop_cell_df$cell.type.ident, 
	levels=c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))

pdf("Boxplots.Percentages.vascECs.human.atleast4ECtypesProfiled.onlyDistalLung.pdf", useDingbats=FALSE)	
prop_cell_df %>% ggplot(aes(y=percentage, x=cell.type.ident, fill=cell.type.ident)) +  
	geom_boxplot(outlier.shape = NA ) + scale_fill_manual(values=celltype_colors) + theme_classic() +
	scale_y_continuous(limits = c(0, 85))
dev.off()


 # make same, but now only if all 5 EC subtypes were profiled and from distal lung
prop_cell_df <- as.data.frame.matrix(prop.table(table(EC_Subset_Object$subject.ident[EC_Subset_Object$location.ident=="Distal_Lung"],
	EC_Subset_Object$cell.type.ident[EC_Subset_Object$location.ident=="Distal_Lung"]), margin=1)*100)
prop_cell_df <-prop_cell_df[rowSums(prop_cell_df!=0) ==5,]
# n=29
summary(prop_cell_df)
 EC_arterial     EC_bronchial     EC_capillary_A   EC_capillary_B
 Min.   : 2.26   Min.   : 0.2841   Min.   : 2.817   Min.   : 1.986
 1st Qu.:12.72   1st Qu.: 2.3697   1st Qu.: 6.061   1st Qu.:18.919
 Median :21.48   Median : 6.6667   Median :15.819   Median :32.143
 Mean   :21.86   Mean   : 8.0826   Mean   :17.514   Mean   :32.201
 3rd Qu.:26.67   3rd Qu.: 9.3863   3rd Qu.:27.273   3rd Qu.:43.333
 Max.   :47.75   Max.   :57.7465   Max.   :39.286   Max.   :70.641
   EC_venous
 Min.   : 4.804
 1st Qu.:10.000
 Median :16.742
 Mean   :20.341
 3rd Qu.:30.806
 Max.   :48.739


# make tidy
prop_cell_df <- gather(prop_cell_df,"cell.type.ident", "percentage")
# make factor to order
prop_cell_df$cell.type.ident <- factor(prop_cell_df$cell.type.ident, 
	levels=c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))

pdf("Boxplots.Percentages.vascECs.human.all5ECtypesProfiled.onlyDistalLung.pdf", useDingbats=FALSE)	
prop_cell_df %>% ggplot(aes(y=percentage, x=cell.type.ident, fill=cell.type.ident)) +  
	geom_boxplot(outlier.shape = NA ) + scale_fill_manual(values=celltype_colors) + theme_classic() +
	scale_y_continuous(limits = c(0, 85))
dev.off()


 
# mouse
cell_df_mouse <- as.data.frame.matrix(table( mouse_EC$orig.ident, mouse_EC$cell.type.ident)) 
cell_df_mouse$ratio_CapBtoCapA <- cell_df_mouse$EC_capillary_B / cell_df_mouse$EC_capillary_A
cell_df_mouse$ratio_CapBtoCapA[cell_df_mouse$EC_capillary_A==0 | cell_df_mouse$EC_capillary_B==0] <- rep(NA)

cell_df_mouse$ratio_CapAtoCapB <- cell_df_mouse$EC_capillary_A / cell_df_mouse$EC_capillary_B 
cell_df_mouse$ratio_CapAtoCapB[cell_df_mouse$EC_capillary_A==0 | cell_df_mouse$EC_capillary_B==0] <- rep(NA)

summary(cell_df_mouse[!is.na(cell_df_mouse$ratio_CapBtoCapA),]$ratio_CapBtoCapA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.218   2.543   4.176   4.870   6.848  10.887

 summary(cell_df_mouse[!is.na(cell_df_mouse$ratio_CapAtoCapB),]$ratio_CapAtoCapB)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.09185 0.14603 0.23956 0.28519 0.39350 0.82108

prop_cell_df_mouse <- as.data.frame.matrix(prop.table(table(mouse_EC$orig.ident, mouse_EC$cell.type.ident), margin=1)*100)
prop_cell_df_mouse <- prop_cell_df_mouse[,c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial")]
summary(prop_cell_df_mouse)
  EC_arterial    EC_capillary_A   EC_capillary_B    EC_venous
 Min.   :1.237   Min.   : 7.794   Min.   :45.34   Min.   : 3.645
 1st Qu.:1.634   1st Qu.:11.301   1st Qu.:61.34   1st Qu.: 7.706
 Median :2.305   Median :16.305   Median :68.79   Median :10.195
 Mean   :3.354   Mean   :17.901   Mean   :68.33   Mean   : 9.901
 3rd Qu.:4.489   3rd Qu.:23.006   3rd Qu.:77.06   3rd Qu.:12.783
 Max.   :9.639   Max.   :37.229   Max.   :84.85   Max.   :14.458
  EC_bronchial
 Min.   :0.0000
 1st Qu.:0.2934
 Median :0.5349
 Mean   :0.5179
 3rd Qu.:0.6581
 Max.   :1.3575

# make tidy
prop_cell_df_mouse <- gather(prop_cell_df_mouse,"cell.type.ident", "percentage")
# make factor to order
prop_cell_df_mouse$cell.type.ident <- factor(prop_cell_df_mouse$cell.type.ident, 
	levels=c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))

pdf("Boxplots.Percentages.vascECs.mouse.pdf", useDingbats=FALSE)	
prop_cell_df_mouse %>% ggplot(aes(y=percentage, x=cell.type.ident, fill=cell.type.ident)) +  
	geom_boxplot(outlier.shape = NA ) + scale_fill_manual(values=celltype_colors) + theme_classic() +
	scale_y_continuous(limits = c(0, 85))
dev.off()



cell_df$species <-"human"
cell_df_mouse$species <-"mouse"

plot_df <- rbind(cell_df[!is.na(cell_df$ratio_CapBtoCapA), c("species", "ratio_CapBtoCapA", "ratio_CapAtoCapB")], 
	cell_df_mouse[!is.na(cell_df_mouse$ratio_CapBtoCapA), c("species", "ratio_CapBtoCapA", "ratio_CapAtoCapB")])

pdf("Boxplots.Ratio.CapBtoCapA.pdf", useDingbats=FALSE)	
plot_df %>% ggplot(aes(y=ratio_CapBtoCapA, x=species, fill=species)) +  
	geom_boxplot(outlier.shape = NA ) + scale_fill_manual(values=species_colors) +
	scale_y_continuous(limits = quantile(plot_df$ratio_CapBtoCapA, c(0.1, 0.9)))  + theme_classic()
dev.off()

pdf("Boxplots.Ratio.CapAtoCapB.pdf", useDingbats=FALSE)	
plot_df %>% ggplot(aes(y=ratio_CapAtoCapB, x=species, fill=species)) +  
	geom_boxplot(outlier.shape = NA ) + scale_fill_manual(values=species_colors) +
	scale_y_continuous(limits = quantile(plot_df$ratio_CapAtoCapB, c(0.1, 0.9)))  + theme_classic()
dev.off()


wilcox.test(x=plot_df$ratio_CapBtoCapA[plot_df$species=="human"], y=plot_df$ratio_CapBtoCapA[plot_df$species=="mouse"])
# p-value = 0.004881

table(plot_df$species)
# human mouse
# 46    18



####################
# make basic stats for cohorts
# human
human_char <- read.table("Supp_Table_S1_SampleOverview.txt", sep="\t", header=T)

# NEC50 was excluded cause of cancer cells
age_df <- human_char %>% distinct(Subject, .keep_all = TRUE) %>% filter(Subject != "EC50") %>% group_by(Cohort) %>% summarize(
	n=length(Subject),
	mean=mean(Age),
	min=min(Age),
	max=max(Age),
	median=median(Age),
	q25= quantile(Age, probs = c(0.25)), 
	q75= quantile(Age, probs = c(0.75)))
write.table(age_df,file="BasicStats_human_byCohort_Age.txt")

summary(human_char %>% distinct(Subject, .keep_all = TRUE) %>% filter(Subject != "EC50") %>% .[,"Age"])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 1.00   32.25   54.00   49.33   63.50   88.00
  
write.table(human_char %>% distinct(Subject, .keep_all = TRUE) %>% filter(Subject != "EC50") %>% group_by(Cohort, Sex) %>% summarise(Sex.n = n()) %>% mutate(Sex.Per = Sex.n / sum(Sex.n)*100),
	file="BasicStats_human_byCohort_Sex.txt")
	
write.table(human_char %>% distinct(Subject, .keep_all = TRUE) %>% filter(Subject != "EC50") %>% group_by(Cohort, Location) %>% summarise(Location.n = n()) %>% mutate(Location.Per = Location.n / sum(Location.n)*100),
	file="BasicStats_human_byCohort_Location.txt")
	
write.table(human_char %>% distinct(Subject, .keep_all = TRUE) %>% filter(Subject != "EC50") %>% drop_na(Ever_Smoker) %>% group_by(Cohort, Ever_Smoker) %>% summarise(Ever_Smoker.n = n()) %>% 
		mutate(Ever_Smoker.Per = Ever_Smoker.n / sum(Ever_Smoker.n)*100),
	file="BasicStats_human_byCohort_Ever_Smoker.txt")	

# mouse
mouse_char <- read.table("Supp_Table_S6_MouseDatasets.txt", sep="\t", header=T)

age_df <- mouse_char %>% drop_na(AvergageAgeWeeks) %>% group_by(Authors) %>% summarize(
	n=length(animal.ident),
	mean=mean(AvergageAgeWeeks),
	min=min(AvergageAgeWeeks),
	max=max(AvergageAgeWeeks),
	median=median(AvergageAgeWeeks),
	q25= quantile(AvergageAgeWeeks, probs = c(0.25)), 
	q75= quantile(AvergageAgeWeeks, probs = c(0.75)))
age_df <- rbind(age_df, c("Raredon_et_al",rep(NA,7)))

write.table(age_df,file="BasicStats_mouse_byCohort_Age.txt")
	
summary(mouse_char %>% .[,"AvergageAgeWeeks"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
   # 8.00   10.75   12.50   30.81   34.00   91.00       2

  
write.table(mouse_char %>% group_by(Authors, Sex) %>% summarise(Sex.n = n()) %>% mutate(Sex.Per = Sex.n / sum(Sex.n)*100),
	file="BasicStats_mouse_byCohort_Sex.txt")

	
pdf("Boxplot.FullDataset_Object.Age.pdf")	
human_char %>% distinct(Subject, .keep_all = TRUE) %>% filter(Subject != "EC50") %>% ggplot(aes(x=Cohort,y=Age,fill=Cohort)) + geom_boxplot() + 
	scale_fill_manual(values=c(cohort_colors,  Yale_Baylor="#0A71B6")) +
	theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
	
pdf("Boxplot.MouseIntegrated.Age.pdf")	
mouse_char %>% ggplot(aes(x=Authors,y=AvergageAgeWeeks,fill=Authors)) + geom_boxplot() + 
	theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
		
library(reshape2) 
pdf("Histogram.FullDataset_Object.Age.pdf")	
human_char %>% distinct(Subject, .keep_all = TRUE) %>% filter(Subject != "EC50") %>% ggplot(aes(x=Age,fill=Cohort)) + geom_histogram(position = "stack", binwidth=5) + 
	scale_fill_manual(values=c(cohort_colors,  Yale_Baylor="#0A71B6")) +
	theme_classic() 
dev.off()



###########################################
# Stuff for revision
##########################################
#### compare with Lafyatis PH dataset
# get PH data
SC207 <- Read10X_h5("SC207_IPAH_36yearFemale_V2_raw_feature_bc_matrix.h5")
SC161 <- Read10X_h5("SC161_IPAH_50yearFemale_V2_raw_feature_bc_matrix.h5")
SC17 <- Read10X_h5("SC17_IPAH_21yearMale_V1_raw_feature_bc_matrix.h5")
SC14 <- Read10X_h5("SC14_CON_76yearMale_V1_raw_feature_bc_matrix.h5")
SC31 <- Read10X_h5("SC31_CON_56yearMale_V1_raw_feature_bc_matrix.h5")
SC31D <- Read10X_h5("SC31D_biologicalDuplicate_CON_56yearMale_V1raw_feature_bc_matrix.h5")
SC45 <- Read10X_h5("SC45_CON_55yearMale_V2_raw_feature_bc_matrix.h5")
SC56 <- Read10X_h5("SC56_CON_57yearFemale_V2_raw_feature_bc_matrix.h5")
SC59 <- Read10X_h5("SC59_CON_18yearMale_V2_raw_feature_bc_matrix.h5")
SC155 <- Read10X_h5("SC155_CON_lowerLobe_23yearFemale_V2_raw_feature_bc_matrix.h5")
SC155D <- Read10X_h5("SC156_CON_upperLobe_23yearFemale_samepatientSC155_V2_raw_feature_bc_matrix.h5")


# the used a number of genes of 200 filter
SC207_obj <- CreateSeuratObject(SC207, min.features = 200, names.field = 0)
SC207_obj$orig.ident <- "SC207"
SC161_obj <- CreateSeuratObject(SC161, min.features = 200, names.field = 0)
SC161_obj$orig.ident <- "SC161"
SC17_obj <- CreateSeuratObject(SC17, min.features = 200, names.field = 0)
SC17_obj$orig.ident <- "SC17"

SC14_obj <- CreateSeuratObject(SC14, min.features = 200, names.field = 0)
SC14_obj$orig.ident <- "SC14"
SC31_obj <- CreateSeuratObject(SC31, min.features = 200, names.field = 0)
SC31_obj$orig.ident <- "SC31"
SC31D_obj <- CreateSeuratObject(SC31D, min.features = 200, names.field = 0)
SC31D_obj$orig.ident <- "SC31D"
SC45_obj <- CreateSeuratObject(SC45, min.features = 200, names.field = 0)
SC45_obj$orig.ident <- "SC45"
SC56_obj <- CreateSeuratObject(SC56, min.features = 200, names.field = 0)
SC56_obj$orig.ident <- "SC56"
SC59_obj <- CreateSeuratObject(SC59, min.features = 200, names.field = 0)
SC59_obj$orig.ident <- "SC59"
SC155_obj <- CreateSeuratObject(SC155, min.features = 200, names.field = 0)
SC155_obj$orig.ident <- "SC155"
SC155D_obj <- CreateSeuratObject(SC155D, min.features = 200, names.field = 0)
SC155D_obj$orig.ident <- "SC155D"


Lafyatis_PH <- merge(x=SC207_obj, y=list( SC161_obj, SC17_obj) , add.cell.ids = c("SC207", "SC161", "SC17"), merge.data = FALSE)
Lafyatis_PH$disease.ident <- "PH"

Lafyatis_CTRL <- merge(x=SC14_obj, y=list(SC31_obj, SC31D_obj, SC45_obj, SC56_obj, SC59_obj, SC155_obj, SC155D_obj) , 
	add.cell.ids = c("SC14", "SC31_obj", "SC31D", "SC45", "SC56", "SC59", "SC155", "SC155D"), merge.data = FALSE)
Lafyatis_CTRL$disease.ident <- "CTRL"

Lafyatis <- merge(Lafyatis_CTRL, Lafyatis_PH)
Lafyatis <- subset(Lafyatis, subset = nCount_RNA > 1000)

Lafyatis[["percent.mito"]] <- PercentageFeatureSet(Lafyatis, pattern = "^MT-")
Lafyatis <- subset(Lafyatis, subset = percent.mito <= 20)

Lafyatis <- NormalizeData(Lafyatis)
Lafyatis <- FindVariableFeatures(Lafyatis, selection.method = "vst", nfeatures = 3000) 
Lafyatis <- ScaleData(Lafyatis, vars.to.regress=c("percent.mito", "nCount_RNA" ))
Lafyatis <- RunPCA(Lafyatis, npcs = 50)
Lafyatis <- RunUMAP(Lafyatis, reduction = "pca", dims = 1:30)
Lafyatis <- FindNeighbors(Lafyatis, reduction = "pca", dims = 1:30)
Lafyatis <- FindClusters(Lafyatis, resolution = 0.15)

Lafyatis_EC <- subset(Lafyatis, subset=seurat_clusters %in% c(3))

DefaultAssay(Lafyatis_EC) <- "RNA"
LafyatisEC.list <- SplitObject(Lafyatis_EC, split.by = "disease.ident")
for (i in 1:length(LafyatisEC.list )) {
    LafyatisEC.list[[i]] <- NormalizeData(LafyatisEC.list[[i]])
    LafyatisEC.list[[i]] <- FindVariableFeatures(LafyatisEC.list[[i]], selection.method = "vst", nfeatures = 3000) 
}

LafyatisEC.anchors <- FindIntegrationAnchors(object.list = LafyatisEC.list, dims = 1:50, anchor.features = 3000) 
LafyatisEC.integrated <- IntegrateData(anchorset = LafyatisEC.anchors, dims = 1:50)

DefaultAssay(LafyatisEC.integrated) <- "integrated"
LafyatisEC.integrated <- ScaleData(LafyatisEC.integrated)
LafyatisEC.integrated <- RunPCA(LafyatisEC.integrated, npcs = 50)
LafyatisEC.integrated <- RunUMAP(LafyatisEC.integrated, reduction = "pca", dims = 1:20,  n.neighbors = 50)
LafyatisEC.integrated <- FindNeighbors(LafyatisEC.integrated, reduction = "pca", dims = 1:20, k.param = 50)
LafyatisEC.integrated <- FindClusters(LafyatisEC.integrated, resolution = 0.3)

LafyatisEC.integrated$subject.ident <- LafyatisEC.integrated$orig.ident
LafyatisEC.integrated$subject.ident[LafyatisEC.integrated$orig.ident == "SC155D"] <- "SC155"
LafyatisEC.integrated$subject.ident[LafyatisEC.integrated$orig.ident == "SC31D"] <- "SC31"

DefaultAssay(LafyatisEC.integrated) <- "RNA"
LafyatisEC.markers <- FindAllMarkers(LafyatisEC.integrated, only.pos=T)
write.table(LafyatisEC.markers, "LafyatisEC.markers.temp.12_06_20.txt", row.names = FALSE, sep="\t")


Idents(LafyatisEC.integrated) <- "seurat_clusters"
LafyatisEC.integrated <- RenameIdents(LafyatisEC.integrated, "0" = "EC_bronchial", "1" = "EC_capillary_B", "2" = "EC_venous", 
    "3" = "EC_arterial", "4" = "EC_capillary_A", "5" = "Multiplet", "6" = "Multiplet", "7" = "Multiplet")
LafyatisEC.integrated$cell.type.ident <- Idents(LafyatisEC.integrated)
# DimPlot(LafyatisEC.integrated, group.by="cell.type.ident", label=T)

# rerun without multiplets
Lafyatis_EC <- subset(LafyatisEC.integrated, subset=cell.type.ident != "Multiplet")

DefaultAssay(Lafyatis_EC) <- "RNA"
LafyatisEC.list <- SplitObject(Lafyatis_EC, split.by = "disease.ident")
for (i in 1:length(LafyatisEC.list )) {
    LafyatisEC.list[[i]] <- NormalizeData(LafyatisEC.list[[i]])
    LafyatisEC.list[[i]] <- FindVariableFeatures(LafyatisEC.list[[i]], selection.method = "vst", nfeatures = 3000) 
}

LafyatisEC.anchors <- FindIntegrationAnchors(object.list = LafyatisEC.list, dims = 1:50, anchor.features = 3000) 
LafyatisEC.integrated <- IntegrateData(anchorset = LafyatisEC.anchors, dims = 1:50)

DefaultAssay(LafyatisEC.integrated) <- "integrated"
LafyatisEC.integrated <- ScaleData(LafyatisEC.integrated)
LafyatisEC.integrated <- RunPCA(LafyatisEC.integrated, npcs = 50)
LafyatisEC.integrated <- RunUMAP(LafyatisEC.integrated, reduction = "pca", dims = 1:20,  n.neighbors = 50)

LafyatisEC.integrated <- FindNeighbors(LafyatisEC.integrated, reduction = "pca", dims = 1:20, k.param = 50)
LafyatisEC.integrated <- FindClusters(LafyatisEC.integrated, resolution = 0.25)

Idents(LafyatisEC.integrated) <- "seurat_clusters"
LafyatisEC.integrated <- RenameIdents(LafyatisEC.integrated, "0" = "EC_bronchial", "1" = "EC_capillary_B", "2" = "EC_venous", 
    "3" = "EC_arterial", "4" = "EC_capillary_A")
LafyatisEC.integrated$cell.type.ident <- Idents(LafyatisEC.integrated)

DefaultAssay(LafyatisEC.integrated) <- "RNA"
Idents(LafyatisEC.integrated) <- "cell.type.ident"
LafyatisEC.markers <- FindAllMarkers(LafyatisEC.integrated, only.pos=T)
write.table(LafyatisEC.markers, "LafyatisEC.markers.12_06_20.txt", row.names = FALSE, sep="\t")

save(LafyatisEC.integrated, file="LafyatisEC.integrated")


######### UMAPs
# UMAP of EC object by cell type, cohort, subject
pdf(file = "LafyatisEC.integrated.UMAP.disease.ident.pdf", useDingbats=FALSE)
DimPlot(LafyatisEC.integrated, reduction = "umap", group.by="disease.ident", cols= c("CTRL"="#63B698","PH"="#C0589E"), 
	cells = sample(colnames(LafyatisEC.integrated)))
dev.off()
pdf(file = "LafyatisEC.integrated.UMAP.cell.type.ident.pdf", useDingbats=FALSE)
DimPlot(LafyatisEC.integrated, reduction = "umap", group.by="cell.type.ident", cols=celltype_colors, 
	cells = sample(colnames(LafyatisEC.integrated)))
dev.off()
pdf(file = "LafyatisEC.integrated.UMAP.subject.ident.pdf", useDingbats=FALSE)
DimPlot(LafyatisEC.integrated, reduction = "umap", group.by="subject.ident", 
	cells = sample(colnames(LafyatisEC.integrated))) 
dev.off()

# make differentials CTRL vs PH
celltypes <- unique(LafyatisEC.integrated$cell.type.ident)
Idents(LafyatisEC.integrated) <- "disease.ident"

results.list <- list()
for(i in 1:length(celltypes)){
results.list[[celltypes[i]]] <- FindMarkers(subset(LafyatisEC.integrated, subset=cell.type.ident==celltypes[i]), 
	ident.1="PH", ident.2="CTRL")
results.list[[celltypes[i]]]$cell.type.ident <- celltypes[i]
results.list[[celltypes[i]]] <- tibble::rownames_to_column(results.list[[celltypes[i]]], "gene")
cat("Done for", celltypes[i],"!\n",sep=" ") 
}

results.df <- do.call(rbind, results.list)
write.table(results.df, "LafyatisEC.PHvsCTRL.12_06_20.txt", row.names = FALSE, sep="\t")

#### Make Cell Type Marker Heatmap
celltypes_to_plot <- as.factor(c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))

genes_for_heatmap <- c("DKK2", "IGFBP3","FBLN5", "SERPINE2", "CLDN10", "GJA5", "CXCL12","BMX", "LTBP4", "HEY1","SOX5","SEMA3G",  # arterial
	"SOSTDC1", "EDNRB", "HPGD", "CYP3A5", "PRKG1", "TBX2", "RCSD1", "EDA", "B3GALNT1", "EXPH5", "NCALD","S100A4", # cap A
	"CA4", "AFF3", "ADGRL2","BTNL9","RGCC", "ADGRF5", "KIAA1217",#common capA and capB
	"FCN3", "IL7R", "CD36",  "NRXN3","SLC6A4", "GPIHBP1","ARHGAP18", "IL18R1", # cap B
	"CPE", "CLU", "C7","PTGS1","EFEMP1","MMRN1","PKHD1L1", "PDZRN4", "DKK3","PLAT", "CDH11","HDAC9", # venous
	"ACKR1","IGFBP7","MCTP1", "VWF", #common bronch and venous
	"COL15A1",  "ZNF385D", "EBF1", "TSHZ2","FLRT2","OLFM1", "CPXM2","PLVAP", "TPD52L1","PDE7B","VWA1", "SPRY1") # peribronchial
### same, but with multithreading, first make all possible combinations
LafyatisEC.integrated$cellBarcode <- colnames(LafyatisEC.integrated) 
cellTypes <- levels(as.factor(LafyatisEC.integrated$cell.type.ident)) 
cellTypes <- cellTypes[cellTypes %in% celltypes_to_plot]
meta.data.sub <- LafyatisEC.integrated@meta.data[,c("disease.ident", "cell.type.ident", "subject.ident", "cellBarcode")]

get.CT.DS.subj.vector <- function(cellTypes){
	tmp.meta.data <- meta.data.sub %>% filter(cell.type.ident== cellTypes)
	cohorts <- unique(tmp.meta.data$disease.ident)
	subjects <- unique(tmp.meta.data$subject.ident)
	tmp.CT.DS.subj.vector <- vector()
	for(j in 1:length(cohorts)){ 
		for(k in 1:length(subjects)){
		temp.cells <- tmp.meta.data %>% filter(disease.ident==cohorts[j] & subject.ident==subjects[k]) %>% pull(cellBarcode)
		if ( length(temp.cells) >1 ) { # >=
			tmp.CT.DS.subj.vector <- c(tmp.CT.DS.subj.vector, paste(cellTypes, cohorts[j], subjects[k], sep="__"))
			}
		}
	}
	cat("Completed for ", cellTypes, ".\n", sep="")
	return(tmp.CT.DS.subj.vector)
}
celltype_cohort_subject.list <- parallel::mclapply(cellTypes, get.CT.DS.subj.vector, mc.cores=ncores)
celltype_cohort_subject <- unlist(celltype_cohort_subject.list)

DefaultAssay(LafyatisEC.integrated) <- "RNA"
get.SubjectcohortCellTypeAvg <- function(celltype_cohort_subject){
	temp.cell.type <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][1]
	temp.cohort <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][2]
	temp.subject <- strsplit(as.character(celltype_cohort_subject),"__")[[1]][3]
	temp.meta.data <- LafyatisEC.integrated@meta.data[,c("disease.ident", "cell.type.ident", "subject.ident", "cellBarcode")]
	temp.cells <- temp.meta.data %>% filter(cell.type.ident==temp.cell.type & disease.ident==temp.cohort & 
		subject.ident==temp.subject) %>% pull(cellBarcode) 
	if (length(temp.cells) > 1) { 
		tmp.df <- as.data.frame(rowMeans(GetAssayData(LafyatisEC.integrated)[,temp.cells]))	   
     } else { 
		tmp.df <- as.data.frame(GetAssayData(LafyatisEC.integrated)[,temp.cells])
		cat("Subject",temp.subject,"only has 1",temp.cell.type,"cell, using singlet for",temp.cohort,"representation...\n",sep=" ") 
    } 
	colnames(tmp.df) <- paste(celltype_cohort_subject)
	return(tmp.df)
}

collapsed.mtx.list <- parallel::mclapply(celltype_cohort_subject, get.SubjectcohortCellTypeAvg, mc.cores=ncores)
collapsed.SubjectcohortCellTypeAvg.mtx <- Matrix(as.matrix(do.call(cbind, collapsed.mtx.list)), sparse = TRUE) 
dim(collapsed.SubjectcohortCellTypeAvg.mtx) 

heatmap_metadata <- as.data.frame(cbind(colnames(collapsed.SubjectcohortCellTypeAvg.mtx), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 1), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 2), 
	sapply(strsplit(as.character(colnames(collapsed.SubjectcohortCellTypeAvg.mtx)),"__"), `[`, 3)))
colnames(heatmap_metadata) <- c("cell.ident","cell.type.ident", "disease.ident", "subject.ident")

heatmap_metadata$cell.type.ident <- factor(heatmap_metadata$cell.type.ident, levels=celltypes_to_plot)
heatmap_metadata$disease.ident <- factor(heatmap_metadata$disease.ident, levels=c("CTRL","PH"))
	
cell_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, disease.ident, subject.ident ) %>% pull(cell.ident)
cohort_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, disease.ident, subject.ident ) %>% pull(disease.ident)
celltype_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, disease.ident, subject.ident ) %>% pull(cell.type.ident)
subject_order <- heatmap_metadata %>% filter(cell.type.ident %in% celltypes_to_plot) %>% 
	arrange(cell.type.ident, disease.ident, subject.ident ) %>% pull(subject.ident )

heatmap_df <-  as.matrix(collapsed.SubjectcohortCellTypeAvg.mtx[genes_for_heatmap,cell_order])
heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))


heatmap_cohort_annotation <- HeatmapAnnotation(cell_type=celltype_order, disease=cohort_order, subject=subject_order, 
    col = list(disease=PH_colors, cell_type=celltype_colors)) #, subject=subject_colors; show_legend = c("subject" = FALSE)
	
pdf("LafyatisEC.Heatmap.Markers.AveragePerSubject.12_08_20.pdf", height=10, width=16)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_split=celltype_order, column_title = NULL, use_raster=FALSE)
dev.off()


#############
# Make average per celltype plots for important GWAS hits
#### downloaded manually from OMIM
# Primary pulmonary arterial hypertension, incl teleangiectasia
# also from https://erj.ersjournals.com/content/53/1/1801899
# and https://www.nature.com/articles/s41569-019-0242-x/tables/2
PH_OMIM_genes <- c("BMPR2", "EIF2AK4","SMAD9", "CAV1",  "ENG", "ACVRL1" ,   "GDF2", "ABCC8", "SOX17", "ATP13A3", "AQP1") 
	#also BMP10, SMAD1, SAMD4; these dont work: "TBX4","KCNK3",
PVOD_OMIM_genes <- c("EIF2AK4", "BMPR2") # overlap
# ALVEOLAR CAPILLARY DYSPLASIA WITH MISALIGNMENT OF PULMONARY VEINS; ACDMPV
ACDMPV_OMIM_genes <- c("FOXF1")
# PULMONARY HYPOPLASIA, PRIMARY
pulm_hypoplasia_OMIM_genes <- c("ZFPM2") # these dont work: "TBX4",  "FGF10", "TBX5"

DefaultAssay(Big_Merged_Object) <- "RNA"
avg.exp.for.OMIM.df <- data.frame(
	"EC_arterial"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident=="EC_arterial")]),
	"EC_capillary_A"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident=="EC_capillary_A")]),
	"EC_capillary_B"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident=="EC_capillary_B")]),
	"EC_venous"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident=="EC_venous")]),
	"EC_bronchial"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident=="EC_bronchial")]),
	"EC_lymph"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident=="EC_lymph")]),
	"Epithelial"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident %in% 
		c("AT1", "AT2", "Basal", "Ciliated", "Club" , "Goblet", "Ionocyte", "PNEC", "Suprabasal"))]),
	"Stromal"=rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident %in% 
		c("Fibroblast_adventitial", "Fibroblast_alveolar", "Mesothel","Myofibroblast", "Pericyte", "SMC"))]),
	"Myloid" = rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident %in% 
		c("AM", "cMono", "ncMono", "DC_Langerhans", "DC_mature", "M", "Mast","pDC", "cDC1", "cDC2"))]),
	"Lymphoid" =rowMeans(GetAssayData(Big_Merged_Object, slot="data")[,WhichCells(Big_Merged_Object,expression=cell.type.ident %in% 
		c("T_helper", "B", "ILC","NK", "Plasma", "T_cytotox", "T_reg"))])
)

genes_to_plot <- c(PH_OMIM_genes, ACDMPV_OMIM_genes, pulm_hypoplasia_OMIM_genes)
heatmap_df <-  as.matrix(avg.exp.for.OMIM.df[genes_to_plot,])
heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))

pdf("EC_Merged_Object.Heatmap.OMIM.Disease.genes.05_20_20.pdf")
Heatmap(heatmap_df_normalized, col = viridis(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
    column_title = NULL, use_raster=FALSE)
dev.off()


#####################
# EC vitro data 
load("ECinVitro.Robj")

## Add percent.mito
ECinVitro[["percent.mito"]] <- PercentageFeatureSet(ECinVitro, pattern = "^MT-")
## Add percent.heatshock
heatshock.genes <- grep(pattern = "^HSP", x = rownames(ECinVitro), value = TRUE)
heatshock.genes <- c(heatshock.genes, grep(pattern = "^DNAJ", x = rownames(ECinVitro), value = TRUE))
ECinVitro[["percent.heatshock"]]<- Matrix::colSums(GetAssayData(ECinVitro, slot = "counts")[heatshock.genes, ])/
	Matrix::colSums(GetAssayData(ECinVitro, slot = "counts"))*100

ECinVitro <- NormalizeData(ECinVitro)
ECinVitro <- FindVariableFeatures(ECinVitro, selection.method = "vst", nfeatures = 1000) #3000
ECinVitro <- ScaleData(ECinVitro, vars.to.regress = c("nCount_RNA","percent.mito")) # 
ECinVitro <- RunPCA(ECinVitro, npcs = 50)
ECinVitro <- RunUMAP(ECinVitro, reduction = "pca", dims = 1:10) 
ECinVitro <- FindNeighbors(ECinVitro, reduction = "pca", dims = 1:10)
ECinVitro <- FindClusters(ECinVitro, resolution = 0.1, algorithm=4) 
ECinVitro$seurat_clusters <- as.character(ECinVitro$seurat_clusters)
ECinVitro$seurat_clusters[ECinVitro$seurat_clusters %in% c(3)] <- rep(2)

Idents(ECinVitro) <- "seurat_clusters"
markers.table <- FindAllMarkers(ECinVitro, only.pos=T)
write.table(markers.table, file="ECinVitro.markers.txt", sep="\t")

p1 <- DimPlot(ECinVitro, reduction = "umap", group.by="orig.ident", label=FALSE, cols=ECinVitro_colors, cells = sample(colnames(ECinVitro)), pt.size=1) 
png(file = "ECinVitro.DimPlot.orig.ident.png", width = 10, height = 10,units = 'in', res = 300)
p1 + NoLegend()
dev.off()

pdf("ECinVitro.UMAP.Legends.pdf", useDingbats=FALSE)
grid.draw(cowplot::get_legend(p1))
grid.newpage()
dev.off()

#### Make Cell Type Marker Heatmap
celltypes_to_plot <- as.factor(c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))
genes_for_heatmap <- c("PECAM1", "CDH5", "ERG", "CLDN5", "TIE1", # generau EC markers
	"DKK2", "IGFBP3","FBLN5", "SERPINE2", "CLDN10", "GJA5", "CXCL12","BMX", "LTBP4", "HEY1","SOX5","SEMA3G",  # arterial
	"SOSTDC1", "EDNRB", "HPGD", "CYP3A5", "PRKG1", "TBX2", "RCSD1", "EDA", "B3GALNT1", "EXPH5", "NCALD","S100A4", # cap A
	"CA4", "AFF3", "ADGRL2","BTNL9","RGCC", "ADGRF5", "KIAA1217",#common capA and capB
	"FCN3", "IL7R", "CD36",  "NRXN3","SLC6A4", "GPIHBP1","ARHGAP18", "IL18R1", # cap B
	"CPE", "CLU", "C7","PTGS1","EFEMP1","MMRN1","PKHD1L1", "PDZRN4", "DKK3","PLAT", "CDH11","HDAC9", # venous
	"ACKR1","IGFBP7","MCTP1", "VWF", #common bronch and venous
	"COL15A1",  "ZNF385D", "EBF1", "TSHZ2","FLRT2","OLFM1", "CPXM2","PLVAP", "TPD52L1","PDE7B","VWA1", "SPRY1") # peribronchial
### same, but with multithreading, first make all possible combinations
ECinVitro$cell.ident <- colnames(ECinVitro) 
heatmap_metadata <- ECinVitro@meta.data[,c( "orig.ident", "cell.ident")]

cell_order <- heatmap_metadata %>% 	arrange(orig.ident ) %>% pull(cell.ident)
subject_order <- heatmap_metadata %>% arrange(orig.ident ) %>% pull(orig.ident )
	
heatmap_df <-  GetAssayData(ECinVitro, slot="data")[genes_for_heatmap,cell_order]
heatmap_df_normalized <- t(apply(heatmap_df, MARGIN=1, FUN=myUnityNormalize))

heatmap_cohort_annotation <- HeatmapAnnotation( Sample=subject_order, col = list(Sample=ECinVitro_colors))
	
pdf("ECinVitro.Heatmap.Markers.12_10_20.pdf", height=10, width=16, useDingbats=FALSE)
Heatmap(heatmap_df_normalized, col = inferno(256), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    top_annotation=heatmap_cohort_annotation, column_title = NULL, use_raster=TRUE, 
	raster_device = "CairoPNG") 
dev.off()

####
# show, a) cell cycle, b) endoMT, c) lymphatic, d) tip like features; 4 per row; + general EC
png(file = "ECinVitro.FeaturePlot.MarkerGenes.png", width = 16, height = 21,units = 'in', res = 450)
FeaturePlot(ECinVitro, c("PECAM1", "CDH5","ERG", "TIE1", #general EC
						"MKI67", "TOP2A", "CENPF", "ASPM", #cell cylce
						"CXCR4","DLL4", "UNC5B", "CD34", # Tip
						"PDPN","RELN","PROX1", "TBX1", # lymph
						"CDH2","VIM", "FN1", "LTBP1", # EndoMT
						"ZEB1", "ZEB2", "TWIST2", "SNAI1"), order=T) & DarkTheme() & scale_color_viridis()
dev.off()


########################################
# Compare technial details

temp_data <- EC_Merged_Object@meta.data[,c("subject.ident", "nCount_RNA", "nFeature_RNA", "cell.type.ident", "cohort.ident")]
temp_data$highlight <- "no"
temp_data$highlight[temp_data$subject.ident %in%  c("115C", "1931C", "026C", "NDRI3")] <- "yes"

set.seed(13)
p1 <- temp_data %>%
    group_by(subject.ident, cohort.ident, highlight) %>%
	summarise(nGene = median(nFeature_RNA)) %>%
	ggplot(aes(y=nGene, x=cohort.ident, fill=cohort.ident)) +  
	geom_boxplot(outlier.shape = NA) + theme_classic() + scale_fill_manual(values=cohort_colors) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_jitter(width = 0.2,  aes(color = highlight) ) + scale_color_manual(values=c(no="black", yes="#39ff14")) +
	scale_y_continuous(limits = c(0, 4000))
pdf("BoxPlot.nFeature_RNA.ByCohort.pdf", useDingbats = FALSE)
print(p1)
dev.off()	

set.seed(11)
p2 <- temp_data %>%
    group_by(subject.ident, cohort.ident, highlight) %>%
	summarise(nUMI = median(nCount_RNA)) %>%
	ggplot(aes(y=nUMI, x=cohort.ident, fill=cohort.ident)) +  
	geom_boxplot(outlier.shape = NA) + theme_classic() + scale_fill_manual(values=cohort_colors) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_jitter(width = 0.2,  aes(color = highlight) ) + scale_color_manual(values=c(no="black", yes="#39ff14")) +
	scale_y_continuous(limits = c(0, 8400))
pdf("BoxPlot.nCount_RNA.ByCohort.pdf", useDingbats = FALSE)
print(p2)
dev.off()	

temp_data %>%
    group_by(subject.ident, cohort.ident) %>%
	summarise(nUMI = median(nCount_RNA)) %>%  kruskal.test(x=.$nUMI, g=.$cohort.ident)
# p-value = 0.008342

temp_data %>%
    group_by(subject.ident, cohort.ident) %>%
	summarise(nGene = median(nFeature_RNA)) %>%  kruskal.test(x=.$nGene, g=.$cohort.ident)
#p-value = 0.006933

temp_data %>%
    group_by(subject.ident, cohort.ident) %>%
	summarise(nUMI = median(nCount_RNA)) %>% 
    split(.$cohort.ident) %>% purrr::map(summary)


temp_data$cell.type.ident <- factor(temp_data$cell.type.ident, levels=
    c("EC_arterial","EC_capillary_A", "EC_capillary_B",  "EC_venous","EC_bronchial"))
    
p1 <- temp_data %>%
    group_by(subject.ident, cell.type.ident, highlight) %>%
	summarise(nUMI = median(nCount_RNA)) %>% 
    arrange(highlight) %>%
	ggplot(aes(y=nUMI, x=cell.type.ident, fill=cell.type.ident)) +  
	geom_boxplot(outlier.shape = NA) + theme_classic() + scale_fill_manual(values=celltype_colors) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_jitter(width = 0.2,  aes(color = highlight)) + scale_color_manual(values=c(no="black", yes="#39ff14")) +
    scale_y_continuous(limits = c(0, 23000))
pdf("BoxPlot.nCount_RNA.ByCellType.pdf", useDingbats = FALSE)
print(p1)
dev.off()	

p2 <- temp_data %>%
    group_by(subject.ident, cell.type.ident, highlight) %>%
	summarise(nGene = median(nFeature_RNA)) %>% 
    arrange(highlight) %>%
	ggplot(aes(y=nGene, x=cell.type.ident, fill=cell.type.ident)) +  
	geom_boxplot(outlier.shape = NA) + theme_classic() + scale_fill_manual(values=celltype_colors) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_jitter(width = 0.2,  aes(color = highlight)) + scale_color_manual(values=c(no="black", yes="#39ff14")) +
    scale_y_continuous(limits = c(0, 5200))
pdf("BoxPlot.nFeature_RNA.ByCellType.pdf", useDingbats = FALSE)
print(p2)
dev.off()	








