#MatrixPlot of Pathways

# Install and load necessary packages ####

setwd("~/Dropbox (Partners HealthCare)/Kim Lab Miscellaneous/PCLS/PCLS_R_Studio_Project/PCLS_Proejct/Figure 1/Matrix Plot")

if (!requireNamespace("enrichR", quietly = TRUE)) 
install.packages("enrichR")
library(enrichR)

# Set the path to the directory containing your CSV files

folder_path <- "~/Dropbox (Partners HealthCare)/Kim Lab Miscellaneous/PCLS/pcls_data/pcls_data-2/DE_gene_lists/overlaps_scatterplots"


# Read all CSV files in the folder ####
csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)

# Initialize a list to store enrichment results
enrichment_results <- list()

# Loop through each file
for (file_path in csv_files) {
  # Read the genes from the CSV file
  genes <- read.csv(file_path, stringsAsFactors = FALSE)$hgnc_symbol_Cryo # Assuming 'Gene' is the column name
  
  # Set the Enrichr libraries you wish to query
  # This is an example, you may want to use other libraries depending on your analysis
  databases <- c("Reactome_2022")
  
  # Connect to Enrichr
  enrichr_dbs <- enrichR::listEnrichrDbs()
  
  # Check if the databases are available and proceed if true
  if (all(databases %in% enrichr_dbs$libraryName)) {
    # Perform enrichment analysis
    enrich_results <- enrichR::enrichr(genes, databases)
    
    # No longer subsetting to top 3, keep all results
    # enrich_results is a list of data frames, one for each database
    enrichment_results[[basename(file_path)]] <- enrich_results
    
  } else {
    warning("Some databases are not available in Enrichr.")
  }
}

# Print or save the top 3 enrichment results for each gene list
print(enrichment_results)

head(enrichment_results$Category_high_low_24hrs.csv)
head(enrichment_results$Category_low_high_24hrs.csv)

str(enrichment_results)
# Optionally, write the results to a file
write(enrichment_results, file = "enrichment_analysis_results_reactome.txt")

#Flattening enrichment results intoa. simple table ####
# Initialize an empty data frame for the flattened data
flattened_data <- data.frame(
  Condition = character(),
  Source = character(),
  Term = character(),
  Combined.Score = numeric(),
  Genes = character(),
  Adjusted.P.value = numeric(),
  Overlap = character(),
  stringsAsFactors = FALSE # to keep text columns as characters
)

# Iterate over the list of enrichment results
for (condition_name in names(enrichment_results)) {
  condition_data <- enrichment_results[[condition_name]]
  
  # Iterate over each source within a condition (e.g., Reactome_2022, BioPlanet_2019)
  for (source_name in names(condition_data)) {
    source_data <- condition_data[[source_name]]
    
    # Check if source_data is not empty
    if(nrow(source_data) > 0) {
      # Prepare a temporary data frame with the data for this condition and source
      temp_data <- data.frame(
        Condition = condition_name,
        Source = source_name,
        Term = source_data$Term,
        Combined.Score = source_data$Combined.Score,
        Adjusted.P.value = source_data$Adjusted.P.value,
        Genes = source_data$Genes,
        Overlap = source_data$Overlap,
        stringsAsFactors = FALSE
      )
      
      # Append the temporary data to the flattened data frame
      flattened_data <- rbind(flattened_data, temp_data)
    }
  }
}

# The flattened_data data frame now contains the flattened information

flattened_data
table(flattened_data$Condition)

# Define the pairs based on your specification ####
condition_pairs <- list(
  TGF_B_168hrs = c("Category_high_low_168hrs.csv", "Category_low_high_168hrs.csv"),
  TGF_B_24hrs = c("Category_high_low_24hrs.csv", "Category_low_high_24hrs.csv"),
  TGF_B_48hrs = c("Category_high_low_48hrs.csv", "Category_low_high_48hrs.csv"),
  TGF_B_336hrs = c("Category_high_low_336hrs.csv", "Category_low_high_336hrs.csv"),
  LPS_6hrs = c("LPS_Category_high_low_6hrs.csv", "LPS_Category_low_high_6hrs.csv"),
  LPS_24hrs = c("LPS_Category_high_low_24hrs.csv", "LPS_Category_low_high_24hrs.csv"),
  LPS_48hrs = c("LPS_Category_high_low_48hrs.csv", "LPS_Category_low_high_48hrs.csv")
)
library(dplyr)

# Assuming you have `flattened_data` and `condition_pairs` defined as shown previously

library(dplyr)

# Assuming you have `flattened_data` and `condition_pairs` defined as shown previously

# Modified function to extract and merge data with defaults for missing pathways
extract_and_merge_data_with_defaults <- function(flattened_data, pair_names, default_combined_score = 1, default_p_value = 1) {
  df1 <- flattened_data %>% filter(Condition == pair_names[1])
  df2 <- flattened_data %>% filter(Condition == pair_names[2])
  
  # Ensure placeholder data frames match original df structures
  columns_needed <- c("Condition","Source", "Term", "Combined.Score", "Adjusted.P.value", "Genes", "Overlap") # Add other columns as necessary
  
  # Find unique terms in both dataframes
  unique_terms_df1 <- unique(df1$Term)
  unique_terms_df2 <- unique(df2$Term)
  
  # Find terms not present in df2 and create a placeholder dataframe with matching columns
  missing_in_df2 <- setdiff(unique_terms_df1, unique_terms_df2)
  if (length(missing_in_df2) > 0) {
    placeholder_df2 <- data.frame(Condition = pair_names[2],
                                  Source = "Reactome_2022",
                                  Term = missing_in_df2, 
                                  Combined.Score = rep(default_combined_score, length(missing_in_df2)), 
                                  Adjusted.P.value = rep(default_p_value, length(missing_in_df2)),
                                  Genes = rep(NA, length(missing_in_df2)),
                                  Overlap = rep(NA, length(missing_in_df2)))
    placeholder_df2 <- placeholder_df2[columns_needed] # Ensure the columns match
    df2 <- rbind(df2, placeholder_df2)
  }
  
  # Find terms not present in df1 and create a placeholder dataframe with matching columns
  missing_in_df1 <- setdiff(unique_terms_df2, unique_terms_df1)
  if (length(missing_in_df1) > 0) {
    placeholder_df1 <- data.frame(Condition = rep(pair_names[1],length(missing_in_df1)),
                                  Source = rep("Reactome_2022",length(missing_in_df1)),
                                  Term = missing_in_df1, 
                                  Combined.Score = rep(default_combined_score, length(missing_in_df1)), 
                                  Adjusted.P.value = rep(default_p_value, length(missing_in_df1)),
                                  Genes = rep(NA, length(missing_in_df1)),
                                  Overlap = rep(NA, length(missing_in_df1)))
    placeholder_df1 <- placeholder_df1[columns_needed] # Ensure the columns match
    df1 <- rbind(df1, placeholder_df1)
  }
  
  # Now merge the adjusted dataframes
  aggregated_data <- merge(df1, df2, by = "Term", suffixes = c(".df1", ".df2"))
  return(aggregated_data)
}

# Loop through each pair to aggregate data with defaults for missing pathways
all_pairs_aggregated_data_with_defaults <- list()

for (pair_name in names(condition_pairs)) {
  pair <- condition_pairs[[pair_name]]
  aggregated_data <- extract_and_merge_data_with_defaults(flattened_data, pair)
  aggregated_data$Pair <- pair_name  # Add an identifier for the pair
  all_pairs_aggregated_data_with_defaults[[pair_name]] <- aggregated_data
}

str(all_pairs_aggregated_data_with_defaults)


TGFB_combined <- all_pairs_aggregated_data_with_defaults[c(1,2,3,4)]
LPS_combined <- all_pairs_aggregated_data_with_defaults[c(5,6,7)]

# Combine and filter for TGF-B
TGFB_data_filtered <- do.call(rbind, lapply(names(TGFB_combined), function(x) {
  TGFB_combined[[x]] %>%
    mutate(Timepoint = x) %>%
    filter(Adjusted.P.value.df1 < 0.01 | Adjusted.P.value.df2 < 0.01)
})) %>%
  bind_rows()

# Combine and filter for LPS
LPS_data_filtered <- do.call(rbind, lapply(names(LPS_combined), function(x) {
  LPS_combined[[x]] %>%
    mutate(Timepoint = x) %>%
    filter(Adjusted.P.value.df1 < 0.01 | Adjusted.P.value.df2 < 0.01)
})) %>%
  bind_rows()


timepoint_shapes <- c('6hrs' = 16, '24hrs' = 15, '48hrs' = 17, '336hrs' = 18, '168hrs' = 18)


plot_data_filtered <- function(data, title) {
  p <- ggplot(data, aes(x = Combined.Score.df2, y = Combined.Score.df1)) +
    geom_point(aes(shape = factor(Timepoint)), color = "black", size = 4) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(title = title, x = "Combined Score Cryo (log10)", y = "Combined Score Fresh (log10)") +
    theme_minimal() +
    theme( axis.line = element_line(color = "black"),
           plot.title = element_blank(),
           legend.title = element_blank(), legend.position = "bottom") +
    labs(shape = "Timepoint") 
  
  return(p)
}

# Generate plots again with the updated function
TGFB_plot_filtered <- plot_data_filtered(TGFB_data_filtered, "TGF-B Conditions (Filtered)")
LPS_plot_filtered <- plot_data_filtered(LPS_data_filtered, "LPS Conditions (Filtered)")

# Display plots
print(TGFB_plot_filtered)
print(LPS_plot_filtered)

#Plotly ####
plot_data_filtered <- function(data, title) {
  # Define the hover text
  hover_text <- paste("Term:", data$Term)
  
  # Create the plotly object
  plotly_obj <- plot_ly(data, x = ~Combined.Score.df2, y = ~Combined.Score.df1, 
                        type = 'scatter', mode = 'markers',
                        hoverinfo = 'text', text = ~hover_text,
                        marker = list(size = 10, color = 'black')) %>%
    layout(title = title,
           xaxis = list(type = "log", title = "Combined Score Cryo (log10)"),
           yaxis = list(type = "log", title = "Combined Score Fresh (log10)"),
           hoverlabel = list(bgcolor = "white"),
           shapes = list(
             list(type = "line",
                  x0 = 1, y0 = 1, x1 = 10000, y1 = 10000, # Adjust x1, y1 as per your axis limits
                  line = list(color = "black", width = 2, dash = "dash"))
           )
    )
  
  return(plotly_obj)
  
}
# Replace 'your_data_frame' with your actual data frame and 'your_title' with your desired plot title
plotly_plot <- plot_data_filtered(TGFB_data_filtered, 'your_title')
plotly_plot

plotly_plot <- plot_data_filtered(LPS_data_filtered, 'your_title')
plotly_plot


