#TGF-B analysis

# From previous step at ClusterHeatmap. Countdata converting ensembl to symbols
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

#Final Files with all sampels are as below (included LPS), so need to substract LPS from it and analyze rest

# Create a DESeqDataSet object
#dds <- DESeqDataSetFromMatrix(countData = count_data_reordered,
                              #colData = phenodata_march24,
                              #design = ~ seurat_clusters)
library(DESeq2)
#Remove LPS
count_data_reordered
phenodata_march24_nolps<- phenodata_march24[phenodata_march24$Test_Substance != 'LPS', ]
count_data_reordered_woLPS<- count_data_reordered[,colnames(count_data_reordered) %in% phenodata_march24_nolps$S_ID]
str(count_data_reordered_woLPS)
str(phenodata_march24_nolps)
phenodata_march24_nolps$Exposure_Time <- factor(phenodata_march24_nolps$Exposure_Time)

phenodata_march24_nolps$Test_Substance <- factor(phenodata_march24_nolps$Test_Substance)
phenodata_march24_nolps$cryo <- factor(phenodata_march24_nolps$cryo)

#table(phenodata_march24_nolps$Exposure_Time, phenodata_march24_nolps$Test_Substance)
#TGF_Beta Vehicle_Control
#6          0               6
#24        18               6
#48        18               6
#168       18               6
#336       18               6

#Need to eliminate time point 6h to be able to have a full matrix rank

phenodata_filtered <- phenodata_march24_nolps[phenodata_march24_nolps$Exposure_Time != 6,]
phenodata_25 <- phenodata_filtered[phenodata_filtered$Concentration == 25 | phenodata_filtered$Concentration == 0, ]
counts_25 <- count_data_reordered_woLPS[, rownames(phenodata_25)]
str(counts_25)
str(phenodata_25)
phenodata_25$Exposure_Time <- as.numeric(as.character(phenodata_25$Exposure_Time))
# Center and scale Exposure_Time
phenodata_25$Exposure_Time <- scale(phenodata_25$Exposure_Time, center = TRUE, scale = FALSE)

dds_25 <- DESeqDataSetFromMatrix(countData = counts_25,
                                       colData = phenodata_25,
                                       design = ~ Test_Substance + cryo + Exposure_Time)

table(phenodata_25)
dds_2b5 <- DESeq(dds_25b)
# Basic result without shrinkage for comparison
res_25b <- results(dds_2b5, contrast=c("Test_Substance", "TGF_Beta", "Vehicle_Control"))
plotMA(res_25b, main="MA Plot without Shrinkage", ylim=c(-6,6))

head(res_25b)
res_25b_df <- as.data.frame(res_25b)
res_25b_sorted <- res_25b_df[order(res_25b_df$log2FoldChange, decreasing = TRUE), ]
head(res_25b_sorted)
#addsymbols
res_25TGFB_sorted_with_symbols <- addGeneSymbols(res_25b_sorted, EnsDb.Hsapiens.v86)
head(res_25TGFB_sorted_with_symbols)

library(ggplot2)
library(ggrepel)
res_25TGFB_sorted_with_symbols$geneCategory <- with(res_25TGFB_sorted_with_symbols, ifelse(padj < 0.05 & log2FoldChange > 1, "Significantly Up",
                                                                 ifelse(padj < 0.05 & log2FoldChange < -1, "Significantly Down", "Not Significant")))
table(res_25TGFB_sorted_with_symbols$geneCategory)


library(ggplot2)
library(ggrepel)


res_25TGFB_sorted_with_symbols$Significance <- with(res_25TGFB_sorted_with_symbols, ifelse(padj < 0.05 & log2FoldChange > 1, "Significantly Up",
                                                                                           ifelse(padj < 0.05 & log2FoldChange < -1, "Significantly Down", "Not Significant")))
# Define color mapping for significance
color_values <- c("Significantly Up" = "indianred1","Significantly Down"="blue", "Not Significant" = "gray")

# Subset the top 15 up-regulated genes
top15_up <- res_25TGFB_sorted_with_symbols[res_25TGFB_sorted_with_symbols$Significance == "Significantly Up",]
top15_up <- top15_up[order(-top15_up$log2FoldChange),][1:20,]

# Subset the top 15 down-regulated genes
top15_down <- res_25TGFB_sorted_with_symbols[res_25TGFB_sorted_with_symbols$Significance == "Significantly Down",]
top15_down <- top15_down[order(top15_down$log2FoldChange),][1:10,]

# Combine the two subsets
top_genes <- rbind(top15_up, top15_down)

volcano_TGF<-ggplot(res_25TGFB_sorted_with_symbols, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(fill = Significance, color = Significance), shape = 21, stroke = 0.5) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  geom_label_repel(
    data = subset(top_genes, Significance != "Not Significant"),
    aes(label = SYMBOL, fill = Significance), # Labels with gene symbols
    fontface = 'bold',
    color = 'black', # Text color
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    max.overlaps = 70
  ) +
  scale_color_manual(values = c("black", "black", "black")) + # Border color for points
  labs(x = "Log2 Fold Change", y = "-Log10(adjusted p-value)", title = "Volcano Plot: Fresh vs Cryo across all samples") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "none") # Hide legend if not needed

ggsave("volcano_TGF-B_vs_VEhicle_global.pdf", plot = volcano_TGF, width = 11, height = 8.5, unit = "in")


install.packages("TmixClust")

#Before using the package, need to prepare the data, columns are timepoints and rows are genes. 
# TMixClust Data handling ####
# Assuming 'dds' is your DESeq2 object
# Step 1: Extract normalized counts
normalized_counts <- counts(dds_2b5, normalized = TRUE)
# Get the condition and time information from the colData
condition <- colData(dds)$Test_Substance
time_points <- colData(dds)$Exposure_Time

# Convert time_points to character to avoid issues if it's a factor
time_points <- as.character(time_points)

# Combine condition and time into a single factor indicating the grouping for each sample
condition_time <- factor(paste(condition, time_points, sep = "_"))


# Step 1: Ensure normalized_counts is a matrix with rownames
normalized_counts_matrix <- matrix(normalized_counts, nrow = nrow(normalized_counts), ncol = ncol(normalized_counts))
rownames(normalized_counts_matrix) <- rownames(normalized_counts)
colnames(normalized_counts_matrix) <- colnames(normalized_counts)

# Step 2: Adjust the split operation to maintain a matrix structure for each subset
split_counts <- lapply(split(colnames(normalized_counts_matrix), f = condition_time), function(cols) {
  normalized_counts_matrix[, cols, drop = FALSE]
})
# Step 3: Now calculate the average counts for each gene across samples within each condition-time combination
avg_counts_list <- lapply(split_counts, rowMeans)
# Convert the list back to a matrix for clustering analysis
avg_counts_matrix <- do.call(cbind, avg_counts_list)
# Ensure the rownames are preserved
rownames(avg_counts_matrix) <- rownames(normalized_counts)

desired_order <- c(24, 48, 168, 336)
# Generate the column names in the desired order for each condition
ordered_colnames <- c(paste("TGF_Beta", desired_order, sep = "_"), 
                      paste("Vehicle_Control", desired_order, sep = "_"))

# Reorder the columns of avg_counts_matrix according to the specified order
avg_counts_matrix_ordered <- avg_counts_matrix[, ordered_colnames]
# Check the structure of the reordered matrix to confirm
str(avg_counts_matrix_ordered)
TGF_matrix<- avg_counts_matrix_ordered[,1:4]
Control_matrix <- avg_counts_matrix_ordered[,5:8]

#TMixClust ####
library(TMixClust)
TGF_matrix<- as.data.frame(TGF_matrix)
Control_matrix<- as.data.frame(Control_matrix)

# Calculate the base mean expression for each gene across all conditions and time points
base_mean_expression <- rowMeans(avg_counts_matrix_ordered)

ggplot(res_25b, aes(x = baseMean)) +
  geom_histogram(binwidth = 1, fill = "dodgerblue", color = "black") +
  scale_x_log10() + # Log-transform the x-axis to better visualize the distribution
  labs(title = "Distribution of Base Mean Expressions", x = "Base Mean Expression (log scale)", y = "Frequency") +
  theme_minimal()

percentile_75 <- quantile(res_25b$baseMean, probs = 0.90, na.rm = TRUE)
histogram_plot <- ggplot(res_25b, aes(x = baseMean)) +
  geom_histogram(binwidth = 1, fill = "dodgerblue", color = "black") +
  geom_vline(xintercept = percentile_75, color = "red", linetype = "dashed", size = 1) +
  scale_x_log10() +
  labs(title = "Distribution of Base Mean Expressions",
       x = "Base Mean Expression (log scale)",
       y = "Frequency") +
  theme_minimal()

print(histogram_plot)

percentile_90 <- quantile(res_25b$baseMean, probs = 0.90, na.rm = TRUE)
#29
# Filter the results to keep only genes above the 90th percentile
filtered_genes <- res_25b[res_25b$baseMean >= percentile_90, ]
filtered_genes <- filtered_genes[!is.na(filtered_genes$pvalue) & filtered_genes$pvalue < 0.10, ]


# Filter genes based on this threshold
filtered_matrix <- as.data.frame(avg_counts_matrix_ordered[rownames(filtered_genes), ])


# Calculate the mean and standard deviation for each gene across time points
gene_means <- apply(filtered_matrix, 1, mean)
gene_sds <- apply(filtered_matrix, 1, sd)

# Subtract the mean and divide by the standard deviation for each gene
z_scores <- sweep(filtered_matrix, 1, gene_means, FUN="-")
z_scores <- sweep(z_scores, 1, gene_sds, FUN="/")

TGF_matrix<- z_scores[,1:4]
Control_matrix <- z_scores[,5:8]

plot_time_series_df(TGF_matrix)
Control_matrix
#TGF_Mix ####
TGF_TMix = TMixClust(TGF_matrix, nb_clusters = 3)

plot_time_series_df(TGF_matrix[TGF_TMix$em_cluster_assignment == 1,])
plot_time_series_df(TGF_matrix[TGF_TMix$em_cluster_assignment == 2,])
plot_time_series_df(TGF_matrix[TGF_TMix$em_cluster_assignment == 3,])

cluster1_tgf<- TGF_matrix[TGF_TMix$em_cluster_assignment == 1,]
cluster2_tgf<- TGF_matrix[TGF_TMix$em_cluster_assignment == 2,]
cluster3_tgf<- TGF_matrix[TGF_TMix$em_cluster_assignment == 3,]



# Now combine all three into a single data frame with an additional 'Cluster' column
all_clusters_tgf <- rbind(
  cbind(cluster1_tgf, Cluster = 1),
  cbind(cluster2_tgf, Cluster = 2),
  cbind(cluster3_tgf, Cluster = 3)
)
test1<- addGeneSymbols(all_clusters_tgf, EnsDb.Hsapiens.v86)
dim(test1)
head(test1)
head(phenodata_25)
pheno_tg<- phenodata_25[phenodata_25$Test_Substance == "TGF_Beta",]
count_tg<- counts_25[,colnames(counts_25) %in% pheno_tg$sampleName]
dds_25 <- DESeqDataSetFromMatrix(countData = count_tg,
                                 colData = pheno_tg,
                                 design = ~  cryo + Exposure_Time)
vst_data <- vst(dds_25, blind=FALSE)
vst_data1 <- assay(vst_data)
vst_data1 <- as.data.frame(vst_data1)
vst_sig <- vst_data1[rownames(vst_data1) %in% test1$gene_ids,]
heat <- t(scale(t(vst_sig)))
dim(vst_sig)
split_1<- pheno_tg$Exposure_Time
split_test <- list(Exposure_Time = pheno_tg$Exposure_Time, cryo = pheno_tg$cryo)


heatmap1<- Heatmap(heat,
                   name = "score",
                   #column_order = row_names_order,
                   #cluster_rows = T,
                   row_split = test1$Cluster, 
                   #cluster_columns = T,
                   #row_km = 4 ,
                   column_split = split_test,
                   cluster_row_slices = F,
                   cluster_column_slices = F,
                   #top_annotation = ha,
                   show_row_names = F,
                   show_column_names = F,
                   column_labels = F,
                   #right_annotation = ra,
                   row_gap = unit(1, "mm"), border = T,
) # Custom color gradient
heatmap1

devtools::install_github("Ylefol/TimeSeriesAnalysis@master")


#Control_Mix ####
Control_TMix = TMixClust(Control_matrix, nb_clusters = 3)
Control_TMix
plot_time_series_df(Control_matrix[Control_TMix$em_cluster_assignment == 1,])
plot_time_series_df(Control_matrix[Control_TMix$em_cluster_assignment == 2,])
plot_time_series_df(Control_matrix[Control_TMix$em_cluster_assignment == 3,])

cluster1_control<- Control_matrix[Control_TMix$em_cluster_assignment == 1,]
cluster2_control<- Control_matrix[Control_TMix$em_cluster_assignment == 2,]
cluster3_control<- Control_matrix[Control_TMix$em_cluster_assignment == 3,]



# Now combine all three into a single data frame with an additional 'Cluster' column
all_clusters_control <- rbind(
  cbind(cluster1_control, Cluster = 1),
  cbind(cluster2_control, Cluster = 2),
  cbind(cluster3_control, Cluster = 3)
)
test2<- addGeneSymbols(all_clusters_control, EnsDb.Hsapiens.v86)
dim(test1)
dim(test2)

table(rownames(test1) == rownames(test2))
res_25TGFB_sorted_with_symbols

deseq2_time_tgf<- res_25TGFB_sorted_with_symbols[res_25TGFB_sorted_with_symbols$gene_ids %in% test1$gene_ids, ]
table(deseq2_time_tgf$Significance)
head(test1)
head(test2)
head(phenodata_25)
pheno_tg<- phenodata_25[phenodata_25$Test_Substance == "TGF_Beta",]
count_tg<- counts_25[,colnames(counts_25) %in% pheno_tg$sampleName]
dds_25 <- DESeqDataSetFromMatrix(countData = count_tg,
                                 colData = pheno_tg,
                                 design = ~  cryo + Exposure_Time)
vst_data <- vst(dds_25, blind=FALSE)
vst_data1 <- assay(vst_data)
vst_data1 <- as.data.frame(vst_data1)
vst_sig <- vst_data1[rownames(vst_data1) %in% test1$gene_ids,]
heat <- t(scale(t(vst_sig)))
dim(vst_sig)
split_1<- pheno_tg$Exposure_Time
heatmap1<- Heatmap(heat,
                   name = "score",
                   #column_order = row_names_order,
                   cluster_rows = T,
                   #row_split = test1$Cluster, 
                   cluster_columns = T,
                   #row_km = 4 ,
                   column_split = split_test,
                   #cluster_row_slices = F,
                   #cluster_column_slices = F,
                   #top_annotation = ha,
                   show_row_names = F,
                   show_column_names = F,
                   column_labels = F,
                   #right_annotation = ra,
                   row_gap = unit(1, "mm"), border = T,
) # Custom color gradient
heatmap1

devtools::install_github("Ylefol/TimeSeriesAnalysis@master")

###TMix did not work lets do the Deseq to 336 and to 24 and 48  
devtools::install_github("Ylefol/TimeSeriesAnalysis@master")



phenodata_filtered <- phenodata_march24_nolps[phenodata_march24_nolps$Exposure_Time != 6,]
phenodata_25 <- phenodata_filtered[phenodata_filtered$Concentration == 25 | phenodata_filtered$Concentration == 0, ]
#phenodata_25s <- phenodata_25[phenodata_25$cryo == 'Fresh',]
phenodata_25s
counts_25s <- count_data_reordered_woLPS[, rownames(phenodata_25s)]
str(counts_25)
str(phenodata_25)
#phenodata_25$Exposure_Time <- as.numeric(as.character(phenodata_25$Exposure_Time))
# Center and scale Exposure_Time
#phenodata_25$Exposure_Time <- scale(phenodata_25$Exposure_Time, center = TRUE, scale = FALSE)

dds_25s <- DESeqDataSetFromMatrix(countData = counts_25s,
                                 colData = phenodata_25s,
                                 design = ~ Test_Substance + Exposure_Time)
ds<- DESeq(dds_25s)
res_fin<- results(ds, contrast=c("Exposure_Time", "336", "24"))
summary(res_fin)
res_fin2<- results(ds, contrast=c("Exposure_Time", "168", "48"))
summary(res_fin2)
res_fin3<- results(ds, contrast=c("Exposure_Time", "336", "48"))
res_fin4<- results(ds, contrast=c("Exposure_Time", "336", "168"))
res_fin5<- results(ds, contrast=c("Exposure_Time", "168", "24"))
res_fin6<- results(ds, contrast=c("Exposure_Time", "48", "24"))
summary(res_fin3)
summary(res_fin4)
summary(res_fin5)
summary(res_fin6)

significant_genes <- lapply(list(res_fin, res_fin2, res_fin3, res_fin4, res_fin5, res_fin6), function(x) {
  sig <- subset(x, padj < 0.05)  # Adjust this threshold as necessary
  rownames(sig)
})

# Step 2: Combine the lists to get a unique list of all DEGs
# For union (all unique DEGs in any list)
all_degs <- unique(unlist(significant_genes))

vst_sig <- vst_data1[rownames(vst_data1) %in% all_degs,]
heat <- t(scale(t(vst_sig)))
dim(vst_sig)
split_1<- pheno_tg$Exposure_Time
# Define discrete and continuous color palettes
categorical_colors <- list(
  "cryo" = c("Cryo" = "dodgerblue3", "Fresh" = "firebrick1"),
  "Test_Substance" = c("LPS" = "gray", "TGF_Beta" = "black", "Vehicle_Control"="white"),
  "Exposure_Time" = c("0"="#2F1163FF", "6"="#9F2F7FFF","24"="#E65163FF", "48"="#F9785DFF", "168"="#FEA36FFF", "336"="#FEC98DFF"), 
  "Concentration" = c("0" = "lightseagreen", "25" = "darkorchid"))

head(pheno_tg)
pheno_for_annotations<- pheno_tg[,2:5]

# Create the annotation
ham <- HeatmapAnnotation(df = pheno_for_annotations,
                         col = categorical_colors,
                         which = "column", # Position the annotation above the heatmap
                         simple_anno_size  = unit(3, "mm"), 
                         border = T,
                         show_annotation_name = F) # Adjust heights as needed
heatmap1<- Heatmap(heat,
                   name = "TGF-B",
                   #column_order = row_names_order,
                   cluster_rows = T,
                   row_km = 3,
                   #row_split = test1$Cluster, 
                   cluster_columns = F,
                   #row_km = 4 ,
                   column_split = split_test,
                   #cluster_row_slices = F,
                   #cluster_column_slices = F,
                   top_annotation = ham,
                   show_row_names = F,
                   show_column_names = F,
                   column_labels = F,
                   #right_annotation = ra,
                   row_gap = unit(0.2, "mm"), border = T,
) # Custom color gradient
heatmap1


save.image("PCLS_April24.RData")
library(dplyr)
res_fin<- as.data.frame(res_fin)
res_fin_s<- addGeneSymbols(res_fin, EnsDb.Hsapiens.v86)
res_finss<- res_fin_s[res_fin_s$padj<0.05,]
res_fin_ss<- top_n(res_finss, 40, log2FoldChange )  
res_fin_ss$SYMBOL



plot(ha)

pheno_control<- phenodata_25[phenodata_25$Test_Substance == "Vehicle_Control",]
count_control<- counts_25[,colnames(counts_25) %in% pheno_control$sampleName]
dds_control <- DESeqDataSetFromMatrix(countData = count_control,
                                 colData = pheno_control,
                                 design = ~  cryo + Exposure_Time)
vst_data_control <- vst(dds_control, blind=FALSE)
vst_data2 <- assay(vst_data_control)
vst_data2 <- as.data.frame(vst_data2)
vst_sig2 <- vst_data2[rownames(vst_data2) %in% all_degs,]
heat2 <- t(scale(t(vst_sig2)))
dim(vst_sig2)

split_2 <- factor(c(24, 24, 24,24, 24, 24, 48, 48, 48, 48, 48, 48, 168, 168, 168, 168, 168, 168, 336, 336, 336, 336, 336, 336), levels = c(24, 48, 168, 336))
split_2 <- factor(split_2, levels = sort(unique(split_2)))
split_2
split_2test <- list(Exposure_Time = pheno_control$Exposure_Time, cryo = pheno_control$cryo)

# Define discrete and continuous color palettes
categorical_colors <- list(
  "cryo" = c("Cryo" = "dodgerblue3", "Fresh" = "firebrick1"),
  "Test_Substance" = c("LPS" = "gray", "TGF_Beta" = "black", "Vehicle_Control"="white"),
  "Exposure_Time" = c("0"="#2F1163FF", "6"="#9F2F7FFF","24"="#E65163FF", "48"="#F9785DFF", "168"="#FEA36FFF", "336"="#FEC98DFF"), 
  "Concentration" = c("0" = "lightseagreen", "25" = "darkorchid"))

head(pheno_control)
pheno_for_annotations2<- pheno_control[,2:5]
pheno_for_annotations2_sorted <- pheno_for_annotations2 %>%
  arrange(as.numeric(Exposure_Time))

# View the sorted data frame
print(pheno_for_annotations2_sorted)
# Create the annotation
ham2 <- HeatmapAnnotation(df = pheno_for_annotations2,
                         col = categorical_colors,
                         which = "column", # Position the annotation above the heatmap
                         simple_anno_size  = unit(3, "mm"), 
                         border = T,
                         annotation_label = c("Test_Substance", "Concentration", "Exposure_Time", "cryo")
                         ) # Adjust heights as needed
heatmap2<- Heatmap(heat2,
                   name = "Control",
                   #column_order = row_names_order,
                   #cluster_rows = T,
                   #row_km = 3,
                   #row_split = test1$Cluster, 
                   cluster_columns = F,
                   #row_km = 4 ,
                   column_split = split_2test,
                   #cluster_row_slices = F,
                   cluster_column_slices = F,
                   top_annotation = ham2,
                   show_row_names = F,
                   show_column_names = F,
                   column_labels = F,
                   #right_annotation = ra,
                   row_gap = unit(0.2, "mm"), border = T,
) # Custom color gradient
heatmap2
htlist = heatmap1 + heatmap2
draw(htlist)

pdf("Combined_temporal_heatmap.pdf", width = 20, height = 10) # Adjust width and height as needed
draw(htlist)
dev.off()

#Linegraphs
avg_counts_matrix_ordered

linegraph_matrix<- as.data.frame(avg_counts_matrix_ordered[rownames(avg_counts_matrix_ordered) %in% all_degs,])
linegraph_matrix<- addGeneSymbols(linegraph_matrix, EnsDb.Hsapiens.v86)
head(linegraph_matrix)

library(ggplot2)
library(tidyr)
library(dplyr)

linegraph_long <- linegraph_matrix %>%
  pivot_longer(cols = -c(gene_ids, SYMBOL), names_to = "Condition_Time", values_to = "Counts") %>%
  mutate(Condition = gsub("(.*)_(.*)", "\\1", Condition_Time), # Extract Condition
         Time = as.numeric(gsub("(.*)_(.*)", "\\2", Condition_Time)))

specific_genes <- c("COL1A1", "THY1", "APLN") # Replace with your actual gene symbols

filtered_data <- linegraph_long %>%
  filter(SYMBOL %in% specific_genes)

ggplot(filtered_data, aes(x = Time, y = Counts, group = interaction(SYMBOL, Condition))) +
  geom_line(aes(color = SYMBOL), size = 1) +  # Lines colored by gene
  geom_point(aes(shape = Condition, color = SYMBOL), size = 3) +  # Points colored by gene
  scale_color_manual(values = c("COL1A1" = "cyan3", "ADAMTS2" = "dodgerblue1", "ACOT8" = "darkorchid2")) +  # Custom colors for each gene
  scale_shape_manual(values = c("TGF_Beta" = 15, "Vehicle_Control" = 16)) +  # Square for TGF_Beta, Circle for Vehicle_Control
  theme_minimal() +
  labs(title = NULL,
       x = "Time (hours)",
       y = "Expression Counts",
       color = "Gene",
       shape = "Condition") +
  theme(legend.title = element_blank(),
        axis.title = element_text(face = "bold"), # Make axis titles bold
        axis.text = element_text(face = "bold"), # Make axis text (numbers) bold
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank())

significant_genes
library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to convert Ensembl IDs to gene symbols
convert_ids_to_symbols <- function(ensembl_ids) {
  genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                 filters = 'ensembl_gene_id', 
                 values = ensembl_ids, 
                 mart = ensembl)
  return(genes)
}

# Apply the function to each element in the list of significant genes
symbols_list <- lapply(significant_genes, convert_ids_to_symbols)

specific_genes <- c("CTHRC1", "POSTN", "GDF15")

# WINNERS: COL1A1, mcemp1, hck
# WINNERS 2: TFRC
filtered_data <- linegraph_long %>%
  filter(SYMBOL %in% specific_genes)

ggplot(filtered_data, aes(x = Time, y = Counts, group = interaction(SYMBOL, Condition))) +
  geom_line(aes(color = SYMBOL), size = 1) +  # Lines colored by gene
  geom_point(aes(shape = Condition, color = SYMBOL), size = 3) +  # Points HCK by gene
  scale_color_manual(values = c("CTHRC1" = "cyan3", "POSTN" = "dodgerblue1", "GDF15" = "darkorchid2")) +  # Custom colors for each gene
  scale_shape_manual(values = c("TGF_Beta" = 15, "Vehicle_Control" = 16)) +  # Square for TGF_Beta, Circle for Vehicle_Control
  theme_minimal() +
  labs(title = NULL,
       x = "Time (hours)",
       y = "Expression Counts",
       color = "Gene",
       shape = "Condition") +
  theme(legend.position = c(0.2, 0.82),
        legend.title = element_blank(),
        axis.title = element_text(face = "bold"), # Make axis titles bold
        axis.text = element_text(face = "bold"), # Make axis text (numbers) bold
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank())


heatmap1
specific_genes <- c("SFTPC", "PKM", "SLPI")
# WINNERS: COL1A1, mcemp1, hck
# WINNERS 2: TFRC
filtered_data <- linegraph_long %>%
  filter(SYMBOL %in% specific_genes)
ggplot(filtered_data, aes(x = Time, y = Counts, group = interaction(SYMBOL, Condition))) +
  geom_line(aes(color = SYMBOL), size = 1) +  # Lines colored by gene
  geom_point(aes(shape = Condition, color = SYMBOL), size = 3) +  # Points HCK by gene
  scale_color_manual(values = c("SFTPC" = "cyan3", "PKM" = "dodgerblue1", "SLPI" = "darkorchid2")) +  # Custom colors for each gene
  scale_shape_manual(values = c("TGF_Beta" = 15, "Vehicle_Control" = 16)) +  # Square for TGF_Beta, Circle for Vehicle_Control
  theme_minimal() +
  labs(title = NULL,
       x = "Time (hours)",
       y = "Expression Counts",
       color = "Gene",
       shape = "Condition") +
  theme(legend.position = c(0.2, 0.82),
        legend.title = element_blank(),
        axis.title = element_text(face = "bold"), # Make axis titles bold
        axis.text = element_text(face = "bold"), # Make axis text (numbers) bold
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank())


heatmap1
r.dend <- row_dend(heatmap1)

#2. Extract clusters
rcl.list <- row_order(heatmap1)

#3.Check cluster size 
lapply(rcl.list, function(x) length(x))

#$`2`
#[1] 212

#$`1`
#[1] 177

#$`3`
#[1] 220


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
library(annotables)
out2<- as.data.frame(out)  %>%
  dplyr::inner_join(grch38, by = c("gene" = "ensgene"))
out2
out2[out2$cluster == "cluster2",]
