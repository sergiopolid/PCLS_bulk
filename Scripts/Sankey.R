#Sankey/Parallel Plot  
#Draw a parallel/alluvial plot right next to the UMAP
#As per Edy, we got to nail down on the message of the UMAP

setwd("~/Dropbox (Partners HealthCare)/Kim Lab Miscellaneous/PCLS/PCLS_R_Studio_Project/PCLS_Project/Figure 1/Sankey")
ls()

#Seems the file to work with is: seurat_object2@meta.data
sankey<- seurat_object2@meta.data
install.packages("ggplot2")
install.packages("ggalluvial")
library(ggalluvial)

# Sample data creation (Replace this with your dataset)
# df <- data.frame(
#   Cluster = sample(1:3, 100, replace = TRUE),
#   Condition = sample(c("Condition1", "Condition2"), 100, replace = TRUE),
#   Cryo_Fresh = sample(c("Cryo", "Fresh"), 100, replace = TRUE),
#   Timepoint = sample(c("T1", "T2", "T3"), 100, replace = TRUE)
# )

# Convert factors to ensure correct ordering
sankey$seurat_clusters <- factor(sankey$seurat_cluster)
sankey$Test_Substance <- factor(sankey$Test_Substance, levels = unique(sankey$Test_Substance))
sankey$cryo <- factor(sankey$cryo, levels = unique(sankey$cryo))
sankey$Exposure_Time <- factor(sankey$Exposure_Time, levels = unique(sankey$Exposure_Time))
sankey$Concentration <- factor(sankey$Concentration, levels = unique(sankey$Concentration))
# Define custom colors for each category
# Adjust these colors to fit your preferences
colors <- c("seurat_clusters" = "#1f77b4", "Test_Substance" = "#ff7f0e", "cryo" = "#2ca02c", "Exposure_Time" = "#d62728")

# Creating the alluvial plot
ggplot(data = sankey,
       aes(axis1 = cryo, axis2 = Test_Substance, axis3 = seurat_clusters, axis4 = Exposure_Time, axis5 = Concentration)) +
  geom_alluvium(aes(fill = factor(seurat_clusters)), width = 0, color="black", size=0.1) + # Use fill to color by seurat_clusters
  geom_stratum(width = 1/6, alpha=0.3) + # Add strata to the sides
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), 
             fill = "white", color = "black", size = 3, fontface = "bold") + # Add labels to the strata
  scale_fill_manual(values = c("0" = "darkorchid", "1" = "lightseagreen", "2" = "gold1"),
                    name = "Clusters") +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom", # Move legend to the bottom
        legend.background = element_rect(fill = "grey90", color = "grey50"), # Grey background for legend
        legend.box.background = element_rect(color = "grey50", linetype = "solid"), # Optional: Border for legend box
        legend.box.margin = margin(6, 6, 6, 6)) # Optional: Adjust legend box margin

ggsave("alluvial_plot.svg", device = "svg", width = 8, height = 6, bg = "white")
ggsave("alluvial_plot.pdf", device = "pdf", width = 8, height = 6)


#Option2
ggplot(data = sankey,
       aes(axis1 = seurat_clusters, axis2 = cryo, axis3 = Test_Substance, axis4 = Exposure_Time)) +
  geom_alluvium(aes(fill = factor(seurat_clusters)), width = 0, color="black", size=0.1) + # Use fill to color by seurat_clusters
  geom_stratum(width = 1/6, alpha=0.3) + # Add strata to the sides
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), 
             fill = "white", color = "black", size = 3, fontface = "bold") + # Add labels to the strata
  scale_fill_manual(values = c("0" = "darkorchid", "1" = "lightseagreen", "2" = "gold1"),
                    name = "Clusters") +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom", # Move legend to the bottom
        legend.background = element_rect(fill = "grey90", color = "grey50"), # Grey background for legend
        legend.box.background = element_rect(color = "grey50", linetype = "solid"), # Optional: Border for legend box
        legend.box.margin = margin(6, 6, 6, 6)) # Optional: Adjust legend box margin

ggsave("alluvial_plot2.svg", device = "svg", width = 8, height = 6, bg = "white")
ggsave("alluvial_plot2.pdf", device = "pdf", width = 8, height = 6)


#Option3
ggplot(data = sankey,
       aes(axis1 = seurat_clusters, axis2 = Test_Substance, axis3 = cryo, axis4 = Exposure_Time)) +
  geom_alluvium(aes(fill = factor(seurat_clusters)), width = 0, color="black", size=0.1) + # Use fill to color by seurat_clusters
  geom_stratum(width = 1/6, alpha=0.3) + # Add strata to the sides
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), 
             fill = "white", color = "black", size = 3, fontface = "bold") + # Add labels to the strata
  scale_fill_manual(values = c("0" = "darkorchid", "1" = "lightseagreen", "2" = "gold1"),
                    name = "Clusters") +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom", # Move legend to the bottom
        legend.background = element_rect(fill = "grey90", color = "grey50"), # Grey background for legend
        legend.box.background = element_rect(color = "grey50", linetype = "solid"), # Optional: Border for legend box
        legend.box.margin = margin(6, 6, 6, 6)) # Optional: Adjust legend box margin

ggsave("alluvial_plot3.svg", device = "svg", width = 8, height = 6, bg = "white")
ggsave("alluvial_plot3.pdf", device = "pdf", width = 8, height = 6)

