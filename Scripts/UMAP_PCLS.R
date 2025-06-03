#UMAP OF THE MASTER_COUNT_TABLE_2, this table was produce form the former script Count_table_Cleaning.R

# Assuming 'master_table_subset2' is already loaded and contains the full dataset with 'ensemble' and 'hgnc_symbol' columns
library(umap)
library(ggplot2)
library(tidyverse)
library(Seurat) # For variable features selection, optional

#Assuming 'master_table_subset2' is already loaded and contains the full dataset
# Prepare the data: Select only numeric columns for UMAP and transpose
# Ensure gene symbols are unique before setting as row names
master_table_numeric <- as.matrix(master_table_subset2[, -c(1, 2)]) # Removing 'ensemble' and 'hgnc_symbol'
unique_gene_symbols <- (master_table_subset2$ensemble)
unique_gene_symbols
rownames(master_table_numeric) <- unique_gene_symbols
# Transpose the table to have samples as rows and genes as columns
samples_data_transposed <- t(master_table_numeric)
str(samples_data_transposed)
head(samples_data_transposed)
# Optional: Feature selection based on variance (using Seurat package for demonstration)
seurat_object <- CreateSeuratObject(master_table_numeric)

# Extracting the names of high variance genes
str(seurat_object)
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = 'dispersion', nfeatures = 4000)

# Perform dimensionality reduction, e.g., PCA followed by UMAP
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- RunUMAP(seurat_object, dims = 1:10) # adjust dimensions as appropriate

# Plot UMAP
DimPlot(seurat_object, reduction = "umap")


# Set the S_ID column from pheno_Clean2 as row names
pheno_for_seurat<- data.frame(pheno_Clean2)
rownames(pheno_for_seurat) <- pheno_for_seurat$S_ID

# Subset the pheno_Clean2 data frame to match the sample names in your Seurat object
# This step ensures that the phenotype data aligns with the samples in your Seurat object
matched_pheno_data <- pheno_for_seurat[colnames(seurat_object), ]

# Add the matched phenotype data as metadata to the Seurat object
seurat_object <- AddMetaData(seurat_object, metadata = matched_pheno_data)

# Now, your Seurat object has the phenotype data included as metadata
# You can check the metadata by accessing the meta.data slot of your Seurat object
head(seurat_object@meta.data)
# If you want to see the plot without running this code due to package dependencies, you can install and load required packages as needed.
# Generate your UMAP plot with DimPlot
p <- DimPlot(seurat_object, reduction = "umap", group.by  = c('Test_Substance', 'cryo', 'Exposure_Time', 'Concentration')) +
  theme(text = element_text(family = "Avenir"), # Set global text elements to Calibri
        axis.title = element_text(family = "Avenir"), # Set axis titles to Calibri
        axis.text = element_text(family = "Avenir"), # Set axis text to Calibri
        legend.title = element_text(family = "Avenir"), # Set legend title to Calibri
        legend.text = element_text(family = "Avenir")) # Set legend text to Calibri
p

seurat_object2 <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object2 <- FindClusters(seurat_object2, resolution = 0.5)
DimPlot(seurat_object2)
head(seurat_object2@meta.data)
p2 <- DimPlot(seurat_object2, reduction = "umap", group.by  = c('Test_Substance', 'cryo', 'Exposure_Time', 'Concentration', 'seurat_clusters')) +
  theme(text = element_text(family = "Avenir"), # Set global text elements to Calibri
        axis.title = element_text(family = "Avenir"), # Set axis titles to Calibri
        axis.text = element_text(family = "Avenir"), # Set axis text to Calibri
        legend.title = element_text(family = "Avenir"), # Set legend title to Calibri
        legend.text = element_text(family = "Avenir")) # Set legend text to Calibri
p2
#Extract umap coordinates
umap_data <- Embeddings(seurat_object2, "umap")
cluster_assignments <- seurat_object2$seurat_clusters
umap_df <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], Cluster = as.factor(cluster_assignments))

library(ggplot2)

# Updated ggplot code with after_stat() for ggplot2 version 3.4.0 and above
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Cluster), alpha = 0.5) + # Scatter plot of UMAP points
  geom_density_2d_filled(aes(fill = Cluster, alpha = after_stat(level)), show.legend = FALSE) + # Density overlay with updated notation
  scale_fill_viridis_d() + # A visually pleasing color scale
  guides(color = guide_legend(title = "Cluster"), alpha = 'none') + # Customize legend title and remove alpha legend
  theme_minimal() +
  theme(text = element_text(family = "Avenir"),
        axis.title = element_text(family = "Avenir"),
        axis.text = element_text(family = "Avenir"),
        legend.title = element_text(family = "Avenir"),
        legend.text = element_text(family = "Avenir")) +
  labs(title = "UMAP with Density Overlay",
       x = "UMAP 1",
       y = "UMAP 2")

library(dplyr)

library(dplyr)
library(ggplot2)

# Initialize an empty data frame for storing hull coordinates
hulls_df <- data.frame(Cluster = character(), UMAP1 = numeric(), UMAP2 = numeric(), stringsAsFactors = FALSE)

# Loop through each cluster to calculate convex hulls and compile into hulls_df
unique_clusters <- unique(umap_df$Cluster)

for (cluster in unique_clusters) {
  # Subset points for the current cluster
  cluster_points <- subset(umap_df, Cluster == cluster)
  
  # Calculate convex hull indices
  hull_indices <- chull(cluster_points$UMAP1, cluster_points$UMAP2)
  
  # Add hull coordinates to hulls_df
  hulls_df <- rbind(hulls_df, data.frame(Cluster = cluster, 
                                         UMAP1 = cluster_points$UMAP1[hull_indices], 
                                         UMAP2 = cluster_points$UMAP2[hull_indices]))
}

# Ensure 'Cluster' is a factor for consistent coloring
hulls_df$Cluster <- factor(hulls_df$Cluster)

# Plot UMAP Results with Convex Hulls
ggplot() +
  geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster), alpha = 0.5) +
  geom_polygon(data = hulls_df, aes(x = UMAP1, y = UMAP2, fill = Cluster, group = Cluster), color = NA, alpha = 0.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "UMAP with Cluster Convex Hulls", x = "UMAP 1", y = "UMAP 2")

#Bring other metadata to the umap
# Add additional metadata fields from Seurat object's meta data to umap_df
umap_df$Test_Substance <- seurat_object2@meta.data$Test_Substance[match(row.names(umap_df), row.names(seurat_object2@meta.data))]
umap_df$Cryo <- seurat_object2@meta.data$cryo[match(row.names(umap_df), row.names(seurat_object2@meta.data))]
umap_df$Exposure_Time <- seurat_object2@meta.data$Exposure_Time[match(row.names(umap_df), row.names(seurat_object2@meta.data))]
umap_df$Concentration <- seurat_object2@meta.data$Concentration[match(row.names(umap_df), row.names(seurat_object2@meta.data))]

# Plot UMAP Results with Points Colored by "Test_Substance"
# Plot UMAP Results with Custom Colors for Points Based on "Test_Substance"
ggplot() +
  geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2, color = Test_Substance), alpha = 0.5) +
  geom_polygon(data = hulls_df, aes(x = UMAP1, y = UMAP2, fill = Cluster, group = Cluster), color = NA, alpha = 0.2) +
  scale_color_manual(values = c("LPS" = "orange1", 
                                "Vehicle_Control" = "royalblue2", 
                                "TGF_Beta" = "red1")) +
  scale_fill_viridis_d() + # This is for the cluster fill
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(family = "Avenir", size = 14), # Set legend text font
        legend.title = element_text(family = "Avenir", size = 14), # Set legend title font
        axis.line = element_line(color = "black", size = 0.5)) +
  labs(x = "UMAP 1", y = "UMAP 2")

install.packages("extrafont")
library(extrafont)
font_import()  # You only need to run this once to import fonts into R
loadfonts(device = "pdf")# For Windows, use loadfonts(device = "win")

final_f1<- ggplot() +
  geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2, color = Test_Substance), alpha = 0.5) +
  geom_polygon(data = hulls_df, aes(x = UMAP1, y = UMAP2, fill = Cluster, group = Cluster), color = NA, alpha = 0.2) +
  scale_color_manual(values = c("LPS" = "orange1", 
                                "Vehicle_Control" = "royalblue2", 
                                "TGF_Beta" = "red1")) +
  scale_fill_viridis_d() + # This is for the cluster fill
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
        legend.text = element_text(family = "Avenir", size = 10), # Set legend text font
        legend.title = element_text(family = "Avenir", size = 10), # Set legend title font
        axis.line = element_line(color = "black", size = 0.5)) +
  labs(title="Treatment by Cluster", x = "UMAP 1", y = "UMAP 2")
final_f1


######Combined UMAPS

install.packages("patchwork")
library(ggplot2)
library(patchwork)
head(umap_df)

# UMAP colored by 'cryo'
p_cryo <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cryo)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Fresh" = "goldenrod1", "Cryo" = "purple1")) +
  theme_minimal() +
  labs(title = "Preservation", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none", axis.line = element_line(color='black', size=0.5), plot.title = element_text(hjust = 0.5, face = "bold")) # Center and bold title)

# UMAP colored by 'Concentration'
p_concentration <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = factor(Concentration))) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(palette = "Spectral") +
  theme_minimal() +
  labs(title = "Concentration", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none", axis.line = element_line(color='black', size=0.5), plot.title = element_text(hjust = 0.5, face = "bold")) # Center and bold title)

# UMAP colored by 'Exposure_Time'
p_exposure_time <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = factor(Exposure_Time))) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(title = "Exposure Time", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none", axis.line = element_line(color='black', size=0.5), plot.title = element_text(hjust = 0.5, face = "bold")) # Center and bold title)

# Combine and arrange the plots
p_combined <- (final_f1 | (p_cryo / p_concentration / p_exposure_time)) + 
  plot_layout(ncol = 2, widths = c(3, 1)) &
  theme(legend.position = "bottom")

# Print the combined plot
# Applying theme settings to include a gray box around the legends
# Removing legend titles from the combined plot
p_final <- p_combined + plot_annotation(tag_levels = 'A') & 
  plot_layout(guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    legend.background = element_rect(fill = "gray90", color = NA, size = 0.1), # Gray background with black border
    legend.title = element_blank() # Remove legend titles
  )


p_final
ggsave(plot=p_final, filename = "figure1AtoD_8_10.svg", device=svg, width = 8, height = 10)
ggsave(plot=p_final, filename = "figure1AtoD_12_14.svg", device=svg, width = 10, height = 14)
ggsave(plot=p_final, filename = "figure1AtoD_10_12.svg", device=svg, width = 10, height = 12)



