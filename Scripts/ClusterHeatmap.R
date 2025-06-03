#### Cluster Heatmap
#Continuing work from UMAP PCLS R script
DE_cluster<- FindAllMarkers(seurat_object2, assay = 'RNA', only.pos = T)

head(DE_cluster)

#                                 ensembl - symbol translation  ####
# load the biomaRt package
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")}
library(biomaRt)

# Set up BioMart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Assuming DE_tops$gene contains the Ensembl IDs
ensembl_ids <- DE_cluster$gene

# Step 3: Retrieve gene symbols
genes_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                    filters = 'ensembl_gene_id', 
                    values = ensembl_ids, 
                    mart = mart)

# merge this data back
DE_tops_with_symbols <- merge(DE_cluster, genes_info, by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE, no.dups=T)

str(DE_tops_with_symbols)

# Now DE_tops_with_symbols contains your original data plus a column 'external_gene_name' for the gene symbols
topheatmap<- DE_tops_with_symbols %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 50, wt = desc(avg_log2FC))

str(topheatmap)

# topheatmap$gene will be the list of genes for the final heatmap

list_of_genes_heatmap_ensemble<- topheatmap$gene
list_of_genes_heatmap_symbol<- topheatmap$external_gene_name
pheno


#                                 Cluster info into phenoData   ####

head(seurat_object2@meta.data)
clusters_seurat_heatmap <-  seurat_object2@meta.data
str(clusters_seurat_heatmap)
clusters_seurat_heatmap <- clusters_seurat_heatmap[,15:17]
phenodata_march24<- left_join(pheno_data, clusters_seurat_heatmap)
head(phenodata_march24)
tail(phenodata_march24, 100)
phenodata_march24<- phenodata_march24[order(phenodata_march24$seurat_clusters),]

table(phenodata_march24$seurat_clusters, phenodata_march24$cryo)

#                                 DEseq2 objects                ####

head(phenodata_march24)


# countdata converting ensembl to symbols
dim(master_table_numeric)
names_keys<- rownames(master_table_numeric)

library(annotables)
library(DESeq2)
names_translated<- grch38 %>% dplyr::filter(ensgene %in% names_keys)
head(names_translated)
dim(names_translated)
head(master_table_numeric)

seurat_clusters

head(master_table_numeric)
head(phenodata_march24)
row_names_order1 <- rownames(phenodata_march24)
count_data_reordered <- master_table_numeric[, match(row_names_order1, colnames(master_table_numeric))]
all(colnames(count_data_reordered) == row_names_order1)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data_reordered,
                              colData = phenodata_march24,
                              design = ~ seurat_clusters)

# Run the DESeq pipeline
dds <- DESeq(dds)

BiocManager::install("biobroom")
library(broom)
# Comparing Condition1 vs Condition2
res1_vs_2 <- results(dds, contrast=c("seurat_clusters", "0", "1"))
res1_vs_2<- as.data.frame(res1_vs_2)
res1_vs_2$gene <- rownames(res1_vs_2)
head(res1_vs_2)
# Comparing Condition1 vs Condition3
res1_vs_3 <- results(dds, contrast=c("seurat_clusters", "0", "2"))
res1_vs_3<- as.data.frame(res1_vs_3)
res1_vs_3$gene <- rownames(res1_vs_3)
head(res1_vs_3)
# Comparing Condition2 vs Condition3
res2_vs_3 <- results(dds, contrast=c("seurat_clusters", "1", "2"))
res2_vs_3<- as.data.frame(res2_vs_3)
res2_vs_3$gene <- rownames(res2_vs_3)
head(res2_vs_3)
#adding the annotable info to DGE tables
res1_vs_2<- res1_vs_2 %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("gene" = "ensgene")) %>% 
  dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, padj, symbol, description, chr)


res1_vs_3<- res1_vs_3 %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("gene" = "ensgene")) %>% 
  dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, padj, symbol, description, chr)


res2_vs_3<- res2_vs_3 %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("gene" = "ensgene")) %>% 
  dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, padj, symbol, description, chr) 
 

# Example: Selecting genes with padj < 0.05 in any comparison
# Assuming the gene column is named 'gene', adjust if it's named differently
sig_genes1_vs_2 <- res1_vs_2 %>%
  filter(padj < 0.01, baseMean > 50, abs(log2FoldChange) > 2) %>%
  pull(gene)

sig_genes1_vs_3 <- res1_vs_3 %>%
  filter(padj < 0.01, baseMean > 50, abs(log2FoldChange) > 2) %>%
  pull(gene)

sig_genes2_vs_3 <- res2_vs_3 %>%
  filter(padj < 0.01, baseMean > 50, abs(log2FoldChange) > 2) %>%
  pull(gene)

sig_genes <- unique(c(sig_genes1_vs_2, sig_genes1_vs_3, sig_genes2_vs_3))

vst_data <- vst(dds, blind=FALSE)
vst_data1 <- assay(vst_data)
vst_data1 <- as.data.frame(vst_data1)
vst_sig <- vst_data1[rownames(vst_data1) %in% sig_genes,]
heat <- t(scale(t(vst_sig)))

head(grch38)
heat_w_symbols<- as.data.frame(heat) %>%
  mutate(gene = rownames(heat)) %>%
  dplyr::inner_join(grch38, by = c("gene" = "ensgene"))

heat_w_symbols$symbol

#                                           Heatmap generation  #####
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(viridis)

Heatmap(heat,
        name = "Expression",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        #top_annotation = HeatmapAnnotation(df = data.frame(SampleType = rep(c("Type1", "Type2"), each = ncol(expression_data)/2))),
        row_title = "Top  Genes",
        column_title = "Samples")


phenodata_march24_2 <- phenodata_march24[,2:8]
phenodata_march24_2 <- phenodata_march24_2[,-5]
phenodata_march24_2 <- phenodata_march24_2[,-5]
head(phenodata_march24_2)
rownames(phenodata_march24_2)<- phenodata_march24$sampleName
magma(20)
# Define discrete and continuous color palettes
categorical_colors <- list(
  "cryo" = c("Cryo" = "dodgerblue3", "Fresh" = ""),
  "Test_Substance" = c("LPS" = "gray", "TGF_Beta" = "black", "Vehicle_Control"="red"),
  "Exposure_Time" = c("0"="#2F1163FF", "6"="#9F2F7FFF","24"="#E65163FF", "48"="#F9785DFF", "168"="#FEA36FFF", "336"="#FEC98DFF"), 
  "seurat_clusters" = c("0" = "darkorchid", "1" = "lightseagreen", "2" = "gold1"))


# Create the annotation
ha <- HeatmapAnnotation(df = phenodata_march24_2,
                        col = categorical_colors,
                        which = "column", # Position the annotation above the heatmap
                        simple_anno_size  = unit(3, "mm"), 
                        border = T) # Adjust heights as needed

plot(ha)

#Row Annotation
heat_w_symbols
row_genes <- c("IL24", "PI3", "CCL2", "LUM", "IL6", "IL1B", "PRELP", "KRT17", "HAS2", "MMP7", "COL1A1", "IL1RN", "GDF15", "SFTA3", "PDGFB", "KYNU")
row_anno1<- heat_w_symbols[heat_w_symbols$symbol %in% row_genes,]
row_anno1
gene_indices <- match(row_anno1$gene, rownames(heat))
head(row_anno1)

ra = rowAnnotation(foo = anno_mark(at = gene_indices, 
                                   labels = row_genes))
ra
split_1<- phenodata_march24_2$seurat_clusters
heatmap1<- Heatmap(heat,
                   name = "score",
                   #column_order = row_names_order,
                   cluster_rows = T,
                   row_km = 4 ,
                   column_split = split_1,
                   top_annotation = ha,
                   show_row_names = F,
                   show_column_names = F,
                   column_labels = F,
                   right_annotation = ra,
                   row_gap = unit(1, "mm"), border = T,
                   ) # Custom color gradient
heatmap1
save.image("PCLS_Feb_2024.RData")

#Extrsacting cluster info
#Here's a breakdown of how this works and how you can apply it to your heatmap1 object
#1. Extract Row Dendrogram: First, you extract the row dendrogram from the heatmap. This dendrogram represents the hierarchical clustering of rows.
heatmap1
r.dend <- row_dend(heatmap1)

#2. Extract clusters
rcl.list <- row_order(heatmap1)

#3.Check cluster size 
lapply(rcl.list, function(x) length(x))

#$`4`
#[1] 31

#$`3`
#[1] 25

#$`2`
#[1] 21

#$`1`
#[1] 38

#4. Extract genes from each cluster
# Assuming 'scaled_data_heatmap' is the original matrix used for the heatmap
for (i in 1:length(rcl.list)) {
  if (i == 1) {
    clu <- row.names(heat)[rcl.list[[i]]]
    out <- cbind(clu, paste("cluster", i, sep=""))
  } else {
    clu <- row.names(heat)[rcl.list[[i]]]
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
out
str(out)
colnames(out) <- c("gene", "cluster")
out2<- as.data.frame(out)  %>%
  dplyr::inner_join(grch38, by = c("gene" = "ensgene"))
out2

#Second Round Row Annotation



# Define discrete and continuous color palettes
categorical_colors <- list(
  "cryo" = c("Cryo" = "dodgerblue3", "Fresh" = ""),
  "Test_Substance" = c("LPS" = "gray", "TGF_Beta" = "black", "Vehicle_Control"="red"),
  "Exposure_Time" = c("0"="#2F1163FF", "6"="#9F2F7FFF","24"="#E65163FF", "48"="#F9785DFF", "168"="#FEA36FFF", "336"="#FEC98DFF"), 
  "seurat_clusters" = c("0" = "darkorchid", "1" = "lightseagreen", "2" = "gold1"))


# Create the annotation
ha <- HeatmapAnnotation(df = phenodata_march24_2,
                        col = categorical_colors,
                        which = "column", # Position the annotation above the heatmap
                        simple_anno_size  = unit(2, "mm"), 
                        border = T, 
                        annotation_legend_param = list(
                          cryo = list(direction = "horizontal"),
                          Test_Substance = list(direction = "horizontal"),
                          Exposure_Time = list(direction = "horizontal"),
                          seurat_clusters = list(direction = "horizontal"),
                          Concentration = list(direction = "horizontal")
                          ), annotation_name_gp = gpar( fontsize = 6)) 
pdf(file = "heatmap_fig1.pdf")
set.seed(123)
heatmap1 <- Heatmap(heat,
                    name = "score",
                    cluster_rows = TRUE,
                    #cluster_columns  = FALSE,
                    cluster_column_slices = TRUE,
                    row_km = 4,
                    column_split = split_1,
                    top_annotation = ha,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    column_labels = FALSE,
                    #right_annotation = ra,
                    row_gap = unit(1, "mm"), 
                    border = TRUE, heatmap_legend_param = list(direction = "horizontal"))

dev.off()

save.image("PCLS_Feb_2024.RData")

set.seed(123)
heatmap1 =
  draw(heatmap1, heatmap_legend_side = "bottom", annotation_legend_side="bottom"); row_order(heatmap1)

clusters_anno<- row_order(heatmap1)
clusters_anno
names(clusters_anno) <- paste0("Cluster", 1:4)

# Map indices to gene names for each cluster
gene_names_clusters <- lapply(clusters_anno, function(indices) {
  rownames(heat)[indices]
})

# Now gene_names_clusters will contain the gene names for each cluster
print(gene_names_clusters)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)

gene_names_clusters$'1'
EnsDb.Hsapiens.v86

# ensembl to symbols
gene_info_cluster1 <-  AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = gene_names_clusters$'1', columns = c("GENEID", "SYMBOL", "ENTREZID"), keytype = "GENEID")
gene_info_cluster2 <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = gene_names_clusters$'2', columns = c("GENEID", "SYMBOL","ENTREZID"), keytype = "GENEID")
gene_info_cluster3 <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = gene_names_clusters$'3', columns = c("GENEID", "SYMBOL","ENTREZID"), keytype = "GENEID")
gene_info_cluster4 <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = gene_names_clusters$'4', columns = c("GENEID", "SYMBOL","ENTREZID"), keytype = "GENEID")

library(clusterProfiler)
library(enrichplot)
enrichResult1 <- enrichGO(gene = gene_info_cluster1$ENTREZID, keyType= "ENTREZID", ont= "BP", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", qvalueCutoff = 0.01)
enrichResult2 <- enrichGO(gene = gene_info_cluster2$ENTREZID, keyType= "ENTREZID", ont= "BP", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", qvalueCutoff = 0.01)
enrichResult3 <- enrichGO(gene = gene_info_cluster3$ENTREZID, keyType= "ENTREZID", ont= "BP", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", qvalueCutoff = 0.01)
enrichResult4 <- enrichGO(gene = gene_info_cluster4$ENTREZID, keyType= "ENTREZID", ont= "BP", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", qvalueCutoff = 0.01)
significant_terms1 <- subset(enrichResult1@result, qvalue <= 0.25)
significant_terms2 <- subset(enrichResult2@result, qvalue <= 0.01)
significant_terms3 <- subset(enrichResult3@result, qvalue <= 0.01)
significant_terms4 <- subset(enrichResult4@result, qvalue <= 0.01)

# Combine the results into a list of tables
combinedResults <- list(
  Cluster1 = significant_terms1$ID,
  Cluster2 = significant_terms2$ID,
  Cluster3 = significant_terms3$ID,
  Cluster4 = significant_terms4$ID
)
combinedResults_SS <- list(
  Cluster1 = enrichResult1@result$ID,
  Cluster2 = enrichResult2@result$ID,
  Cluster3 = enrichResult3@result$ID,
  Cluster4 = enrichResult4@result$ID
)
install.packages("gridtext")
library(gridtext)
simplifyGOFromMultipleLists(combinedResults, padj_cutoff = 0.001)

#Saved as multiplelist.png


clusters_annosub<- clusters_anno

enrichResult3@result$ID

set.seed(123)
library(circlize)
library(simplifyEnrichment)
draw(heatmap1, heatmap_legend_side = "bottom", annotation_legend_side="bottom") + rowAnnotation(go = anno_word_cloud_from_GO(align_to  = clusters_anno,term =  combinedResults))
set.seed(123)
Heatmap(heat,
                    name = "score",
                    cluster_rows = TRUE,
                    #cluster_columns  = FALSE,
                    cluster_column_slices = TRUE,
                    row_km = 4,
                    column_split = split_1,
                    top_annotation = ha,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    column_labels = FALSE,
                    #right_annotation = ra,
                    row_gap = unit(1, "mm"), 
                    border = TRUE, heatmap_legend_param = list(direction = "horizontal")) + rowAnnotation(go = anno_word_cloud_from_GO(align_to  = clusters_anno,term =  combinedResults))
#draw(heatmap1)

dim(heat)

str(combinedResults_SS)
31+25+21+38
# Call count_words with test terms

# Inspect the frequencies
print(test_freq)






row_genes <- c("TRAF1","SAA1","MMP9","CCL7","KYNU","FGF2","GNLY","COL1A2", "CEACAM6", "IL24", "PI3", "CCL2", "LUM", "IL6", "IL1B", "PRELP", "KRT17", "HAS2", "MMP7", "COL1A1", "IL1RN", "GDF15", "SFTA3", "PDGFB")
row_anno1<- heat_w_symbols[heat_w_symbols$symbol %in% row_genes,]
row_anno1
gene_indices <- match(row_anno1$gene, rownames(heat))
head(row_anno1)

ra = rowAnnotation(genes = anno_mark(at = gene_indices, 
                                     labels = row_genes, labels_gp = gpar(fontsize = 6)))
ra
draw(heatmap1, heatmap_legend_side = "bottom", annotation_legend_side="bottom")




# Initialize a vector of zeros with length equal to the number of rows in your matrix
cluster_indices <- rep(0, length = 115)

# Assign cluster numbers to each row based on your clusters_anno list
for (cluster_name in names(clusters_anno)) {
  cluster_number <- as.integer(sub("Cluster", "", cluster_name))  # Extract the number from the cluster name
  row_indices <- clusters_anno[[cluster_name]]
  cluster_indices[row_indices] <- cluster_number
}

# Output the vector
print(cluster_indices)
combinedResults <- list(
  "1" = significant_terms1$ID,
  "2" = significant_terms2$ID,
  "3" = significant_terms3$ID,
  "4" = significant_terms4$ID
)
combinedResults <- list(
  "1" = significant_terms1,
  "2" = significant_terms2$ID,
  "3" = significant_terms3$ID,
  "4" = significant_terms4$ID
)
combinedResults_SS <- list(
  "1" = enrichResult1@result$ID,
  "2" = enrichResult2@result$ID,
  "3" = enrichResult3@result$ID,
  "4" = enrichResult4@result$ID
)
set.seed(123)
pdf("heatmap_with_anno.pdf")
heatmap1<- Heatmap(heat,
        name = "score",
        cluster_rows = F,
        #cluster_columns  = FALSE,
        cluster_column_slices = TRUE,
        #row_km = 4,
        column_split = split_1,
        top_annotation = ha,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_labels = FALSE,
        #right_annotation = ra,
        row_gap = unit(1, "mm"),
        row_split  = cluster_indices,
        border = TRUE, heatmap_legend_param = list(direction = "horizontal")) + rowAnnotation(GO= anno_word_cloud_from_GO(cluster_indices,combinedResults, fontsize_range = c(2,10)))
heatmap1 = draw(heatmap1)
row_order(heatmap1)
dev.off()
#draw(heatmap1)

draw(heatmap1, heatmap_legend_side = "bottom", annotation_legend_side="bottom") + rowAnnotation(GO= anno_word_cloud_from_GO(cluster_indices,combinedResults, max_words = 30))
dim(heat)


set.seed(123)
heatmap2 =
  draw(heatmap1, heatmap_legend_side = "bottom", annotation_legend_side="bottom") + rowAnnotation(GO= anno_word_cloud_from_GO(cluster_indices,combinedResults, max_words = 30))



heatmap3<- Heatmap(heat,
        name = "score",
        cluster_rows = F,
        column_split = split_1,
        top_annotation = ha,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_labels = FALSE,
        #right_annotation = ra,
        row_gap = unit(1, "mm"),
        row_split  = cluster_indices,
        border = TRUE,row_title = NULL, column_title = NULL,
        show_row_dend = FALSE, show_column_dend = FALSE) + rowAnnotation(GO= anno_word_cloud_from_GO(cluster_indices,combinedResults, exclude_words = "regulation"))

