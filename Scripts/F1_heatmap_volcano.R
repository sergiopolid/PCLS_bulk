#                                       DEseq2: Fresh vs. Cryo (Figure 1). ####

# Initialization 
# Install DESeq2 if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
# Load DESeq2
library(DESeq2)

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("Bioconductor/biomaRt", ref = "3.15") # replace X_Y with the Bioconductor release that matches your R version

# Load the biomaRt library
library(biomaRt)
library(dplyr)

#Load count table and pheno table

str(matched_pheno_data)
#Clean the pheno table
matched_pheno_data2<- matched_pheno_data[,-(1:2)]
str(matched_pheno_data2)
matched_pheno_data2<- matched_pheno_data2[,-(1:4)]
str(matched_pheno_data2)
#138 obs of 6 variables

str (master_table_numeric)
#[1:65988, 1:138]

#Changing the name of genes to symbols

# Select the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listMarts()
# Your list of Ensembl IDs
ensembl_ids <- row.names(master_table_numeric)
str(ensembl_ids)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# Get the corresponding HGNC symbols
hgnc_symbols <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                      filters = 'ensembl_gene_id',
                      values = ensembl_ids,
                      mart = ensembl)
# Check the results
str(hgnc_symbols)
duplicated_genes <- master_table_numeric[duplicated(rownames(master_table_numeric)) | duplicated(rownames(master_table_numeric), fromLast = TRUE), ]
#none

# Ensure the row names of master_table_numeric are character if not already
rownames(master_table_numeric) <- as.character(rownames(master_table_numeric))

# Create a named vector from hgnc_symbols for mapping
# This vector uses ensembl_gene_id as names and hgnc_symbol as values
symbols_vector <- setNames(hgnc_symbols$hgnc_symbol, hgnc_symbols$ensembl_gene_id)

# Replace the row names in master_table_numeric with the mapped HGNC symbols
# Unmatched Ensembl IDs will have NA as their row names, which you might want to handle
new_row_names <- symbols_vector[rownames(master_table_numeric)]
missing_symbols <- is.na(new_row_names) # Identify which rows did not have a matching symbol

# Optionally, handle rows with no matching symbols
# For simplicity, we'll keep the original Ensembl ID as the row name where no symbol is found
new_row_names[missing_symbols] <- rownames(master_table_numeric)[missing_symbols]

# Apply the new row names back to the master_table_numeric
rownames(master_table_numeric) <- new_row_names

# Check for any missing symbols if needed
sum(missing_symbols)

table(missing_symbols)

str(master_table_numeric)
head(master_table_numeric)


# Assuming 'count_data' is your matrix with gene counts
# and 'pheno_data' is a dataframe extracted from your provided phenotypic data
count_data<- master_table_numeric
head(count_data)

# The 'pheno_data' dataframe should look something like this:
pheno_data <- matched_pheno_data2
head(pheno_data)
# Ensure that the row names in 'pheno_data' match the column names in 'count_data'
#rownames(pheno_data) = pheno_data$sampleName

#                                               DESeq2 ####

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = pheno_data,
                              design = ~ cryo)

anyDuplicated(rownames(master_table_numeric))
# Run the DESeq pipeline
dds <- DESeq(dds)

# Extract results for the 'cryo' condition
results_cryo <- results(dds, contrast = c("cryo", "Fresh", "Cryo"))

# Order results by p-value
results_cryo <- results_cryo[order(results_cryo$pvalue), ]

str(results_cryo)
# Subset for significant results
sig_results_cryo <- subset(results_cryo, baseMean > 10)
dim(sig_results_cryo)
#10475 6

#sig_results_cryo <- subset(sig_results_cryo, padj < 0.10)
#dim(sig_results_cryo)

# View the top results
head(sig_results_cryo, 50)

#                                                 Volcano                      ####

# Plot a volcano plot or other suitable visualization
library(dplyr)
# Add a column for adjusted p-values or use p-value directly for simplicity
sig_results_cryo$logpadj <- -log10(sig_results_cryo$padj)
sig_results_cryo
sig_results_cryo$hgnc_symbol<- row.names(sig_results_cryo)
sig_results_cryo<- as.data.frame(sig_results_cryo)
# Flag the top 20 upregulated and top 20 downregulated genes based on logFC
top_upregulated <- sig_results_cryo[order(sig_results_cryo$log2FoldChange, decreasing = TRUE), ][1:20, ]
# Find top 20 downregulated genes
top_downregulated <- sig_results_cryo[order(sig_results_cryo$log2FoldChange), ][1:20, ]
# Combine the top upregulated and downregulated genes and mark them
top_genes <- rbind(
  transform(top_upregulated, topGene = "Top 20 Up"),
  transform(top_downregulated, topGene = "Top 20 Down")
)
head(sig_results_cryo)

head(top_genes)
library(ggplot2)
library(ggrepel)

# Create the volcano plot
head(sig_results_cryo)
head(top_genes)

#Joining tables to highlight the top 20 genes I want
volcanoplot_table <- left_join(sig_results_cryo, top_genes, by='hgnc_symbol')
volcanoplot_table <- volcanoplot_table[,-(9:15)]
head(volcanoplot_table)
volcanoplot_table <- volcanoplot_table %>%
  mutate(topGene = replace_na(topGene, "Other"))

#Remove three first genes, because of way to high pvals that offset the heatmap

# Step 1: Create a new column for gene classification based on your criteria
volcanoplot_table$geneCategory <- with(volcanoplot_table, ifelse(padj.x < 0.05 & log2FoldChange.x > 1.5, "Significantly Up",
                                                                 ifelse(padj.x < 0.05 & log2FoldChange.x < -1.5, "Significantly Down", "Not Significant")))
head(volcanoplot_table)
# Step 2: Plot with ggplot2, coloring dots based on the geneCategory
plot1 <- ggplot(volcanoplot_table, aes(x = log2FoldChange.x, y = logpadj.x)) +
  geom_point(aes(fill = geneCategory), color = "black", alpha = 0.5, shape = 21, stroke = 0.5) + # Points with conditional fill
  scale_fill_manual(values = c("Significantly Up" = "red", "Significantly Down" = "blue", "Not Significant" = "gray")) +
  geom_label_repel(data = subset(volcanoplot_table, geneCategory %in% c("Significantly Down", "Significantly Up")),
                   aes(label = hgnc_symbol), size = 3, box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"), max.overlaps = 70,
                   fill = "white", label.size = 0.1, segment.color = 'grey50') +
  labs(x = "Log2 Fold Change", y = "-Log10(adjusted p-value)", title = "Fresh vs Cryo across all samples") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(color = "black", size = 14), # X-axis title
        axis.title.y = element_text(color = "black", size = 14), # Y-axis title
        axis.text.x = element_text(color = "black", size = 12), # X-axis text
        axis.text.y = element_text(color = "black", size = 12)) # Y-axis text

ggsave("my_volcano_plot.svg", plot = plot1, width = 16, height = 12)
ggsave("my_volcano_plot.pdf", plot = plot1, width = 16, height = 12)

# Export results to CSV if needed
write.csv(as.data.frame(sig_cryo_fresh_symbol), "DE_cryo_vs_fresh.csv")
str(master_table_numeric)
master_table_numeric_symbol<- replace_ensembl_with_hgnc(master_table_numeric, hgnc_symbols)
str(master_table_numeric_symbol)


#                                           Heatmap generation ####

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

#Listing DE genes
table(volcanoplot_table$geneCategory)
subset_for_heatmap<- as.data.frame(volcanoplot_table)
subset_for_heatmap<- subset_for_heatmap[subset_for_heatmap$geneCategory != "Not Significant",]

#NAs issues
table(is.na(count_data))
# Check if all DE genes are present in count_data
all(subset_for_heatmap$hgnc_symbol %in% rownames(count_data))
# Summarize NA presence
summary(rowSums(is.na(count_data_subset_heatmap)))
summary(colSums(is.na(count_data_subset_heatmap)))

# Identify DE genes not in count_data
missing_genes <- subset_for_heatmap$hgnc_symbol[!subset_for_heatmap$hgnc_symbol %in% rownames(count_data)]
print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
#[1] "Missing genes: AGER.6, .12335"
subset_for_heatmap2<- subset_for_heatmap[subset_for_heatmap$hgnc_symbol != 'AGER.6',]
subset_for_heatmap2<- subset_for_heatmap2[subset_for_heatmap2$hgnc_symbol != '.12335',]
length(subset_for_heatmap2$hgnc_symbol)
length(rownames(count_data))
#subsetting the count matrix for the DE genes
str(count_data)
head(count_data)
count_data_subset_heatmap<- count_data[rownames(count_data) %in% subset_for_heatmap2$hgnc_symbol,]
#By doing these, tehre are 4 genes duplicated so will sum the counts
index1<- rownames(count_data_subset_heatmap)
count_data_subset_heatmap2<- as.data.frame(count_data_subset_heatmap)
head(count_data_subset_heatmap2)
count_data_subset_heatmap2$geneSymbol<- index1
aggregated_data <- aggregate(. ~ geneSymbol, data = count_data_subset_heatmap2, sum)
head(aggregated_data)
aggregated_data$geneSymbol
rownames(aggregated_data)<- aggregated_data$geneSymbol
aggregated_data<- aggregated_data[,-1]
head(aggregated_data)
scaled_data_heatmap <- t(scale(t(aggregated_data)))
table(is.na(scaled_data_heatmap))

Heatmap(scaled_data_heatmap,
        name = "Expression",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        #top_annotation = HeatmapAnnotation(df = data.frame(SampleType = rep(c("Type1", "Type2"), each = ncol(expression_data)/2))),
        row_title = "Top  Genes",
        column_title = "Samples")


head(pheno_data)


magma(100)
# Define discrete and continuous color palettes
categorical_colors <- list(
  "cryo" = c("Cryo" = "dodgerblue3", "Fresh" = "deeppink4"),
  "Test_Substance" = c("LPS" = "lightseagreen", "TGF_Beta" = "lightyellow1", 
                       "Vehicle_Control"="lightskyblue3"),# Adjust according to your actual treatments
  "Exposure_Time" = c("0"="#2F1163FF", "6"="#9F2F7FFF","24"="#E65163FF", "48"="#F9785DFF", "168"="#FEA36FFF", "336"="#FCFDBFFF"))


# Create the annotation
ha <- HeatmapAnnotation(df = pheno_data[,(2:5)],
                        col = categorical_colors,
                        which = "column", # Position the annotation above the heatmap
                        annotation_height = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), border = T) # Adjust heights as needed

# Calculate the 1st and 99th percentiles of the data
low <- quantile(scaled_data_heatmap, probs = 0.01, na.rm = TRUE)
high <- quantile(scaled_data_heatmap, probs = 0.99, na.rm = TRUE)
n_colors <- 10

# Generate a color vector using the magma function from viridis package
magma_colors <- magma(n_colors)
magma_colors
# Define a color mapping function using these bounds with the magma palette
#col_fun <- colorRamp2(c(low, high), c(magma_colors[1], magma_colors[256]))
plot(col_fun)
# Define a color mapping function using these bounds
col_fun <- colorRamp2(c(low,1.48, high), c("#440154FF","#1F998AFF", "#FDE725FF"))
set.seed(342)

head(pheno_data)

split_cryo <- pheno_data$cryo
heatmap1<- Heatmap(scaled_data_heatmap,
        name = "Z-score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        top_annotation = ha, # Use the annotation created above
        row_title = "Top Genes",
        column_title = "Samples",
        col = col_fun, row_km = 3, row_gap = unit(1, "mm"), border = T) # Custom color gradient

save.image(heatmap1, file="heatmap1.svg")

install.packages("viridis")
library(viridis)


save.image("PCLS_Feb_2024.RData")

#extract dendogram and clsutering genes from the heatmap
#The approach you found leverages the row_dend and row_order functions to extract and work with the clustering results from a heatmap created with the ComplexHeatmap package. This method does indeed provide a way to access the cluster information indirectly by using the row order of dendrogram leaves, which corresponds to the clustered rows in the heatmap.
#Here's a breakdown of how this works and how you can apply it to your heatmap1 object
#1. Extract Row Dendrogram: First, you extract the row dendrogram from the heatmap. This dendrogram represents the hierarchical clustering of rows.

r.dend <- row_dend(heatmap1)

#2. Extract clusters
rcl.list <- row_order(heatmap1)

#3.Check cluster size 
lapply(rcl.list, function(x) length(x))

#4. Extract genes from each cluster
# Assuming 'scaled_data_heatmap' is the original matrix used for the heatmap
for (i in 1:length(rcl.list)) {
  if (i == 1) {
    clu <- row.names(scaled_data_heatmap)[rcl.list[[i]]]
    out <- cbind(clu, paste("cluster", i, sep=""))
  } else {
    clu <- row.names(scaled_data_heatmap)[rcl.list[[i]]]
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
out
str(out)

out1<- out[1:14]
out2<- out[15:56]
out3<- out[57:90]
out1
out2
out3


#Crete GO from clusters
setwd("~/Dropbox (Partners HealthCare)/Kim Lab Miscellaneous/PCLS/pcls_data/pcls_data-2/Figure1/Enrichments")
list.files()

t1 <- read_excel("HuBMAP_ASCTplusB_augmented_2022_table.xlsx")
t2 <- read_excel("HuBMAP_ASCTplusB_augmented_2022_table2.xlsx", 
                 +     col_types = c("text", "text", "text", 
                                     +         "text", "numeric", "numeric", "text", 
                                     +         "numeric", "text"))

t3 <- read_excel("Reactome_2022_table1.xlsx", 
                    col_types = c("text", "text", "text", 
                                "text", "numeric", "numeric", "text", 
                               "numeric", "text"))
t4 <- read_excel("Reactome_2022_table2.xlsx",
                 col_types = c("text", "text", "text", 
                               "text", "numeric", "numeric", "text", 
                               "numeric", "text"))

str(t1)
str(t2)
#t1
t1$`Combined Score`
top5_df <- t1[order(-t1$`Combined Score`),]
head(t1)
top5_df <- top5_df[1:10,]
top5_df

# Create the bar graph
a1<- ggplot(top5_df, aes(x = `Combined Score`, y = reorder(Term, `Combined Score`))) +
  geom_bar(stat = "identity", fill = "deepskyblue1", color='black') +
  geom_text(aes(label = Term, x = `Combined Score` / 2), hjust = 0.5, color = "white", fontface='bold', size = 4) +  # Adjust x position for centering
  labs(x = "Combined Score", y = "", title = "Top 5 Gene Ontologies by Combined Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 12, face = "bold"),  # X axis numbers bold and bigger
    axis.title.x = element_text(size = 14, face = "bold"),  # X axis title bigger
    axis.text.y = element_blank(),  # Remove Y axis labels
    axis.ticks.y = element_blank(),  # Remove Y axis ticks
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),  # Black line for X axis
    axis.line.y = element_line(color = "black", size = 0.5)   # Black line for Y axis
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +  # Adjust expansion to ensure consistent axis appearance
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))  # Ensure X axis starts at 0 and has a slight expansion to the right

# Calculate the range for the X axis based on your data
x_max <- ceiling(max(top5_df$`Combined Score`)/500) * 500

# Add vertical dashed lines at every 500 units
if(x_max > 0) {
  for(i in seq(500, x_max, by = 500)) {
    a1 <- a1 + geom_vline(xintercept = i, linetype = "dashed", color = "azure3")
  }
}

a1
a1 <- ggplot() +
  # Add vertical dashed lines first to ensure they're in the background
  geom_vline(data = data.frame(x = seq(500, x_max, by = 500)), aes(xintercept = x), linetype = "dashed", color = "gray") +
  # Then add the bars, so they're drawn over the lines
  geom_bar(data = top5_df, aes(x = `Combined Score`, y = reorder(Term, `Combined Score`)), stat = "identity", fill = "deepskyblue1", color = "black") +
  geom_text(data = top5_df, aes(label = Term, x = `Combined Score` / 2, y = reorder(Term, `Combined Score`)), hjust = 0.5, color = "white", fontface='bold', size = 4) +
  labs(x = "Combined Score", y = "", title = "Top 5 Gene Ontologies by Combined Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))
a1
ggsave('enrichment_reactome1.svg', a1, device='svg', width=10, height=10)

t2
#t2 or enrichement by cell types in upregulated 
t2$`Combined Score`
top5_df <- t2[order(-t2$`Combined Score`),]
head(t2)
top5_df <- top5_df[1:10,]
top5_df

# Create the bar graph
a1<- 0
a1 <- ggplot() +
  # Add vertical dashed lines first to ensure they're in the background
  geom_vline(data = data.frame(x = seq(500, x_max, by = 500)), aes(xintercept = x), linetype = "dashed", color = "gray") +
  # Then add the bars, so they're drawn over the lines
  geom_bar(data = top5_df, aes(x = `Combined Score`, y = reorder(Term, `Combined Score`)), stat = "identity", fill = "darkviolet", color = "black") +
  geom_text(data = top5_df, aes(label = Term, x = `Combined Score` / 2, y = reorder(Term, `Combined Score`)), hjust = 0.5, color = "white", fontface='bold', size = 4) +
  labs(x = "Combined Score", y = "", title = "Top 5 Gene Ontologies by Combined Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 2500))
a1
ggsave('enrichment_celltype2.svg', a1, device='svg', width=10, height=10)



t3
#t2 or enrichement by cell types in upregulated 
t3$`Combined Score`
top5_df <- t3[order(-t3$`Combined Score`),]
head(t3)
top5_df <- top5_df[1:10,]
top5_df

# Create the bar graph
a1<- 0
a1 <- ggplot() +
  # Add vertical dashed lines first to ensure they're in the background
  geom_vline(data = data.frame(x = seq(500, x_max, by = 500)), aes(xintercept = x), linetype = "dashed", color = "gray") +
  # Then add the bars, so they're drawn over the lines
  geom_bar(data = top5_df, aes(x = `Combined Score`, y = reorder(Term, `Combined Score`)), stat = "identity", fill = "gold", color = "black") +
  geom_text(data = top5_df, aes(label = Term, x = 1000 , y = reorder(Term, `Combined Score`)), hjust = 0.5, color = "black", fontface='bold', size = 4) +
  labs(x = "Combined Score", y = "", title = "Top 5 Gene Ontologies by Combined Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 2000))
a1
ggsave('enrichment_reactome1.svg', a1, device='svg', width=10, height=10)


t4
#t2 or enrichement by cell types in upregulated 
t4$`Combined Score`
top5_df <- t4[order(-t4$`Combined Score`),]
head(t4)
top5_df <- top5_df[1:10,]
top5_df

# Create the bar graph
a1<- 0
a1 <- ggplot() +
  # Add vertical dashed lines first to ensure they're in the background
  geom_vline(data = data.frame(x = seq(500, x_max, by = 500)), aes(xintercept = x), linetype = "dashed", color = "gray") +
  # Then add the bars, so they're drawn over the lines
  geom_bar(data = top5_df, aes(x = `Combined Score`, y = reorder(Term, `Combined Score`)), stat = "identity", fill = "cyan4", color = "black") +
  geom_text(data = top5_df, aes(label = Term, x = 600  , y = reorder(Term, `Combined Score`)), hjust = 0.5, color = "black", fontface='bold', size = 4) +
  labs(x = "Combined Score", y = "", title = "Top 5 Gene Ontologies by Combined Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 1500))
a1
ggsave('enrichment_reactome2.svg', a1, device='svg', width=10, height=10)


