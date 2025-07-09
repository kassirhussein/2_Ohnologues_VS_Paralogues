# Load required libraries
library(tidyverse)        # Includes dplyr, ggplot2, tidyr, etc.
library(UpSetR)           # For creating UpSet plots
library(grid)             # For advanced graphical layout
library(readxl)           # To read Excel files
library(writexl)          # To write Excel files
library(VennDiagram)      # To draw Venn diagrams

# --------------------------------------------------------------
# Step 1: Load protein-coding gene list
# --------------------------------------------------------------

# Read protein-coding genes from text file
protein_coding_genes <- read_delim("C:/Users/HP-ssd/Desktop/Short term project/protein coding genes/gene_with_protein_product.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

# Extract gene symbols into a list
protein_coding_genes_list <- protein_coding_genes$symbol

# --------------------------------------------------------------
# Step 2: Load FUSIL dataset
# --------------------------------------------------------------

# Read in the FUSIL dataset containing gene symbols
fusil_m_gene <-  read_delim("C:/Users/HP-ssd/Desktop/Short term project2/fusil.csv")

# Check the number of unique gene symbols
length(unique(fusil_m_gene$gene_symbol))

# --------------------------------------------------------------
# Step 3: Load human gene paralogues from BioMart
# --------------------------------------------------------------

# Read in paralogues from CSV
human_gene_paralogues <- read.csv("C:/Users/HP-ssd/Desktop/Short term project2/paralogues/human_gene_paralogues.csv")

# Clean and rename relevant columns
human_gene_paralogues <- human_gene_paralogues %>%
  dplyr::select(-1,-2,-4) %>%                              # Drop unwanted columns
  rename(gene_symbol = external_gene_name)                # Rename for consistency

# Check number of unique genes
length(unique(human_gene_paralogues$gene_symbol))

# --------------------------------------------------------------
# Step 4: Merge FUSIL and paralogues datasets
# --------------------------------------------------------------

paralogue_fusil <- human_gene_paralogues %>%
  left_join(fusil_m_gene) %>%                             # Merge FUSIL info
  dplyr::select(-6, -7) %>%                               # Remove unnecessary columns
  mutate(hsapiens_paralog_associated_gene_name = na_if(hsapiens_paralog_associated_gene_name,"")) %>%  # Replace empty strings with NA
  distinct()                                              # Remove duplicate rows

# --------------------------------------------------------------
# Step 5: Load ohnologs (relaxed criteria)
# --------------------------------------------------------------

ohnologs_relaxed <- read_delim("C:/Users/HP-ssd/Desktop/Short term project/ohnologs/hsapiens.Pairs.Relaxed.2R.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

# --------------------------------------------------------------
# Step 6: Start analysis loop for different similarity thresholds
# --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(VennDiagram)
library(ggplot2)
library(grid)

# Define thresholds for paralogue similarity
thresholds <- c(30, 50, 70)

# Loop through each similarity threshold
for (threshold in thresholds) {
  
  # --------------------------------------------------------------
  # Step 1: Filter ohnologs relevant to FUSIL genes
  # --------------------------------------------------------------
  ohnologs_filtered <- ohnologs_relaxed %>%
    filter(Symbol1 %in% fusil_m_gene$gene_symbol) %>%     # Keep only those in FUSIL
    dplyr::select(Symbol1, Symbol2) %>%                   # Select gene and its paralogue
    rename(gene_symbol = Symbol1, gene_paralogue = Symbol2)  # Rename for consistency
  
  # --------------------------------------------------------------
  # Step 2: Filter BioMart paralogues based on threshold
  # --------------------------------------------------------------
  biomart_filtered <- human_gene_paralogues %>%
    filter(hsapiens_paralog_perc_id > threshold) %>%      # Apply similarity threshold
    filter(gene_symbol %in% protein_coding_genes$symbol) %>%  # Keep only protein-coding
    na.omit() %>%
    dplyr::select(gene_symbol, hsapiens_paralog_associated_gene_name) %>%  # Keep necessary columns
    rename(gene_paralogue = hsapiens_paralog_associated_gene_name)
  
  # --------------------------------------------------------------
  # Step 3: Count paralogues per gene for each source
  # --------------------------------------------------------------
  count_biomart <- biomart_filtered %>%
    group_by(gene_symbol) %>%
    tally(name = "n_biomart")                             # Count from BioMart
  
  count_ohnologue <- ohnologs_filtered %>%
    group_by(gene_symbol) %>%
    tally(name = "n_ohnologue")                           # Count from Ohnologues
  
  # --------------------------------------------------------------
  # Step 4: Merge paralogue counts and identify overlap
  # --------------------------------------------------------------
  merged_counts <- full_join(count_biomart, count_ohnologue, by = "gene_symbol") %>%
    replace_na(list(n_biomart = 0, n_ohnologue = 0)) %>%  # Fill missing values with 0
    mutate(
      in_biomart = n_biomart > 0,
      in_ohnologue = n_ohnologue > 0
    )
  
  # --------------------------------------------------------------
  # Step 5: Generate Venn Diagram of overlapping genes
  # --------------------------------------------------------------
  biomart_genes <- merged_counts$gene_symbol[merged_counts$in_biomart]
  ohnologue_genes <- merged_counts$gene_symbol[merged_counts$in_ohnologue]
  
  grid.newpage()
  draw.pairwise.venn(
    area1 = length(biomart_genes),
    area2 = length(ohnologue_genes),
    cross.area = length(intersect(biomart_genes, ohnologue_genes)),
    category = c("BioMart", "Ohnologue"),
    fill = c("skyblue", "orange"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.pos = c(-20, 20)
  )
  
  # Add plot title
  grid.text(
    paste("Gene Overlap Between BioMart (>", threshold, "% Similarity) and Ohnologue Sets"),
    x = 0.5, y = 0.95, gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  # --------------------------------------------------------------
  # Step 6: Create scatter plot comparing paralogue counts
  # --------------------------------------------------------------
  plot_data <- merged_counts %>%
    group_by(n_biomart, n_ohnologue) %>%
    summarise(gene_count = n(), .groups = "drop")
  
  print(
    ggplot(plot_data, aes(x = n_biomart, y = n_ohnologue, size = gene_count)) +
      geom_point(alpha = 0.7, color = "steelblue") +                   # Dot size = number of genes
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Reference line
      scale_size_continuous(name = "Gene Count") +
      labs(
        title = paste("Paralogue Count Comparison: nbiomart (> ", threshold, "%) vs nohnologue"),
        x = "nbiomart paralogue count",
        y = "nohnologue paralogue count"
      ) +
      theme_minimal()
  )
}



# ----------------------------------------------------------------------
# Comparing the 2 orthologues datasets with %similarity >30
# ----------------------------------------------------------------------

# Filter ohnologs to include only FUSIL genes
ohnologs_relaxed_filtered <- ohnologs_relaxed %>%
  filter(Symbol1 %in% fusil_m_gene$gene_symbol) %>%
  dplyr::select(3, 4) %>%
  rename("gene_symbol" = "Symbol1") %>%
  rename("gene_paralogue" = "Symbol2")

# Filter BioMart paralogues with %identity >30 and keep protein-coding genes
biomart_paralogue_relaxed_over30 <- human_gene_paralogues %>%
  filter(hsapiens_paralog_perc_id > 30) %>%
  filter(gene_symbol %in% protein_coding_genes$symbol) %>%
  na.omit() %>%
  dplyr::select(1, 2) %>%
  rename("gene_paralogue" = "hsapiens_paralog_associated_gene_name")

# Get overlapping gene pairs between the two sources
common_df <- merge(biomart_paralogue_relaxed_over30, ohnologs_relaxed_filtered, by = c("gene_symbol", "gene_paralogue"))

# Check if all ohnolog genes are in the BioMart dataset
all(unique(ohnologs_relaxed_filtered$gene_symbol) %in% paralogue_fusil$gene_symbol)

# Identify unique gene-paralogue pairs in each dataset
only_in_biomart <- setdiff(biomart_paralogue_relaxed_over30, ohnologs_relaxed_filtered)
only_in_biomart_df <- as.data.frame(only_in_biomart)

only_in_ohnologues <- setdiff(ohnologs_relaxed_filtered, biomart_paralogue_relaxed_over30)
only_in_ohnologues_df <- as.data.frame(only_in_ohnologues)

# Count paralogues per gene in each dataset
count_only_in_biomart <- biomart_paralogue_relaxed_over30 %>%
  group_by(gene_symbol) %>%
  tally() %>%
  rename("n_biomart" = "n")

count_only_in_ohnologues <- ohnologs_relaxed_filtered %>%
  group_by(gene_symbol) %>%
  tally() %>%
  rename("n_ohnologue" = "n")

# Check for unmatched gene symbols
all(count_only_in_ohnologues$gene_symbol %in% count_only_in_biomart$gene_symbol)
setdiff(count_only_in_ohnologues$gene_symbol, count_only_in_biomart$gene_symbol)

# Merge counts into a single dataframe
merged_counts_30 <- full_join(count_only_in_biomart, count_only_in_ohnologues, by = "gene_symbol") %>%
  mutate_all(~replace_na(., 0))

# Add logical presence flags
merged_counts_30$in_biomart <- as.logical(as.integer(merged_counts_30$n_biomart > 0))
merged_counts_30$in_ohnologue <- as.logical(as.integer(merged_counts_30$n_ohnologue > 0))

# Save to CSV
write.csv(merged_counts_30, "C:/Users/HP-ssd/Desktop/merged_counts_30.csv")

# Create Venn diagram
biomart_genes <- merged_counts_30$gene_symbol[merged_counts_30$in_biomart == 1]
ohnologue_genes <- merged_counts_30$gene_symbol[merged_counts_30$in_ohnologue == 1]

venn_plot_30 <- draw.pairwise.venn(
  area1 = length(biomart_genes),
  area2 = length(ohnologue_genes),
  cross.area = length(intersect(biomart_genes, ohnologue_genes)),
  category = c("BioMart", "Ohnologue"),
  fill = c("skyblue", "orange"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-20, 20)
)

grid.text("Gene Overlap Between BioMart (30% Similarity) and Ohnologue Sets", 
          x = 0.5, y = 0.95, gp = gpar(fontsize = 12, fontface = "bold"))

# Create scatter plot of paralogue count comparison
plot_data <- merged_counts_30 %>%
  group_by(n_biomart, n_ohnologue) %>%
  summarise(gene_count = n(), .groups = "drop")

ggplot(plot_data, aes(x = n_biomart, y = n_ohnologue, size = gene_count)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Paralogue Count Comparison: nbiomart (>30% Similarity) vs nohnologue",
    x = "nbiomart paralogue count",
    y = "nohnologue paralogue count",
    size = "Overlapping genes"
  ) +
  theme_minimal()


# -------------------------------------------------------------------------
# Comparing the 2 orthologue datasets with % similarity >50
# -------------------------------------------------------------------------

# Load ohnolog dataset
ohnologs_relaxed <- read_delim("C:/Users/HP-ssd/Desktop/Short term project/ohnologs/hsapiens.Pairs.Relaxed.2R.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

# Filter ohnologs for genes present in the FUSIL list and select relevant columns
ohnologs_relaxed_filtered <- ohnologs_relaxed %>%
  filter(Symbol1 %in% fusil_m_gene$gene_symbol) %>%
  dplyr::select(3, 4) %>%
  rename("gene_symbol" = "Symbol1") %>%
  rename("gene_paralogue" = "Symbol2")

# Filter BioMart paralogues with >50% identity and protein-coding genes
biomart_paralogue_relaxed_over50 <- human_gene_paralogues %>%
  filter(hsapiens_paralog_perc_id > 50) %>%
  filter(gene_symbol %in% protein_coding_genes$symbol) %>%
  na.omit() %>%
  dplyr::select(1, 2) %>%
  rename("gene_paralogue" = "hsapiens_paralog_associated_gene_name")

# Find common gene-paralogue pairs across both datasets
common_df <- merge(biomart_paralogue_relaxed_over50, ohnologs_relaxed_filtered,
                   by = c("gene_symbol", "gene_paralogue"))

# Check whether all ohnolog gene symbols exist in the BioMart + FUSIL merged dataset
all(unique(ohnologs_relaxed_filtered$gene_symbol) %in% paralogue_fusil$gene_symbol)

# Identify gene-paralogue combinations that are exclusive to each dataset
only_in_biomart <- setdiff(biomart_paralogue_relaxed_over50, ohnologs_relaxed_filtered)
only_in_biomart_df <- as.data.frame(only_in_biomart)

only_in_ohnologues <- setdiff(ohnologs_relaxed_filtered, biomart_paralogue_relaxed_over50)
only_in_ohnologues_df <- as.data.frame(only_in_ohnologues)

# Count paralogues per gene in BioMart and Ohnologue datasets
count_only_in_biomart <- biomart_paralogue_relaxed_over50 %>%
  group_by(gene_symbol) %>%
  tally() %>%
  rename("n_biomart" = "n")

count_only_in_ohnologues <- ohnologs_relaxed_filtered %>%
  group_by(gene_symbol) %>%
  tally() %>%
  rename("n_ohnologue" = "n")

# Check if all gene symbols in ohnologs are also present in biomart count
all(count_only_in_ohnologues$gene_symbol %in% count_only_in_biomart$gene_symbol)
setdiff(count_only_in_ohnologues$gene_symbol, count_only_in_biomart$gene_symbol)

# Merge paralogue counts into one data frame and replace NA with 0
merged_counts_50 <- full_join(count_only_in_biomart, count_only_in_ohnologues, by = "gene_symbol") %>%
  mutate_all(~replace_na(., 0))

# Create logical flags for gene presence in either dataset
merged_counts_50$in_biomart <- as.logical(as.integer(merged_counts_50$n_biomart > 0))
merged_counts_50$in_ohnologue <- as.logical(as.integer(merged_counts_50$n_ohnologue > 0))

# Export the merged results to CSV
write.csv(merged_counts_50, "C:/Users/HP-ssd/Desktop/merged_counts_50.csv")

# -------------------------------------------------------------------------
# Generate Venn Diagram of overlapping genes
# -------------------------------------------------------------------------

# Extract genes present in each dataset
biomart_genes <- merged_counts_50$gene_symbol[merged_counts_50$in_biomart == 1]
ohnologue_genes <- merged_counts_50$gene_symbol[merged_counts_50$in_ohnologue == 1]

# Draw the Venn diagram
venn_plot_50 <- draw.pairwise.venn(
  area1 = length(biomart_genes),
  area2 = length(ohnologue_genes),
  cross.area = length(intersect(biomart_genes, ohnologue_genes)),
  category = c("BioMart", "Ohnologue"),
  fill = c("skyblue", "orange"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-20, 20)
)

# Add title to Venn diagram
grid.text("Gene Overlap Between BioMart (50% Similarity) and Ohnologue Sets", 
          x = 0.5, y = 0.95, gp = gpar(fontsize = 12, fontface = "bold"))

# -------------------------------------------------------------------------
# Generate scatter plot comparing paralogue counts between the datasets
# -------------------------------------------------------------------------

# Summarize number of genes by their paralogue count in both datasets
plot_data <- merged_counts_50 %>%
  group_by(n_biomart, n_ohnologue) %>%
  summarise(gene_count = n(), .groups = "drop")

# Plot paralogue count comparison as a scatter plot
ggplot(plot_data, aes(x = n_biomart, y = n_ohnologue, size = gene_count)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Paralogue Count Comparison: nbiomart (>50% Similarity) vs nohnologue",
    x = "nbiomart paralogue count",
    y = "nohnologue paralogue count",
    size = "Overlapping genes"
  ) +
  theme_minimal()


# ------------------------------------------------------------------------------
# Comparing the 2 orthologue datasets with % similarity >70
# ------------------------------------------------------------------------------

# Load ohnolog dataset (relaxed)
ohnologs_relaxed <- read_delim("C:/Users/HP-ssd/Desktop/Short term project/ohnologs/hsapiens.Pairs.Relaxed.2R.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

# Filter ohnologs to retain only FUSIL genes and relevant columns
ohnologs_relaxed_filtered <- ohnologs_relaxed %>%
  filter(Symbol1 %in% fusil_m_gene$gene_symbol) %>%
  dplyr::select(3, 4) %>%
  rename("gene_symbol" = "Symbol1") %>%
  rename("gene_paralogue" = "Symbol2")

# Filter BioMart paralogues with >70% identity, keep only protein-coding genes
biomart_paralogue_relaxed_over70 <- human_gene_paralogues %>%
  filter(hsapiens_paralog_perc_id > 70) %>%
  filter(gene_symbol %in% protein_coding_genes$symbol) %>%
  na.omit() %>%
  dplyr::select(1, 2) %>%
  rename("gene_paralogue" = "hsapiens_paralog_associated_gene_name")

# Find common gene-paralogue pairs in both datasets
common_df <- merge(biomart_paralogue_relaxed_over70, ohnologs_relaxed_filtered,
                   by = c("gene_symbol", "gene_paralogue"))

# Check whether all ohnolog gene symbols are in the BioMart-FUSIL set
all(unique(ohnologs_relaxed_filtered$gene_symbol) %in% paralogue_fusil$gene_symbol)

# Identify gene-paralogue pairs unique to each dataset
only_in_biomart <- setdiff(biomart_paralogue_relaxed_over70, ohnologs_relaxed_filtered)
only_in_biomart_df <- as.data.frame(only_in_biomart)

only_in_ohnologues <- setdiff(ohnologs_relaxed_filtered, biomart_paralogue_relaxed_over70)
only_in_ohnologues_df <- as.data.frame(only_in_ohnologues)

# Count paralogues per gene in each dataset
count_only_in_biomart <- biomart_paralogue_relaxed_over70 %>%
  group_by(gene_symbol) %>%
  tally() %>%
  rename("n_biomart" = "n")

count_only_in_ohnologues <- ohnologs_relaxed_filtered %>%
  group_by(gene_symbol) %>%
  tally() %>%
  rename("n_ohnologue" = "n")

# Check whether all ohnolog genes are in BioMart count
all(count_only_in_ohnologues$gene_symbol %in% count_only_in_biomart$gene_symbol)
setdiff(count_only_in_ohnologues$gene_symbol, count_only_in_biomart$gene_symbol)

# Merge counts and fill missing values with 0
merged_counts_70 <- full_join(count_only_in_biomart, count_only_in_ohnologues, by = "gene_symbol") %>%
  mutate_all(~replace_na(., 0))

# Create logical flags indicating gene presence in each dataset
merged_counts_70$in_biomart <- as.logical(as.integer(merged_counts_70$n_biomart > 0))
merged_counts_70$in_ohnologue <- as.logical(as.integer(merged_counts_70$n_ohnologue > 0))

# Export merged dataset
write.csv(merged_counts_70, "C:/Users/HP-ssd/Desktop/merged_counts_70.csv")

# ------------------------------------------------------------------------------
# Venn Diagram: overlap of genes between BioMart and Ohnologues at >70%
# ------------------------------------------------------------------------------

# Extract gene sets for Venn diagram
biomart_genes <- merged_counts_70$gene_symbol[merged_counts_70$in_biomart == 1]
ohnologue_genes <- merged_counts_70$gene_symbol[merged_counts_70$in_ohnologue == 1]

# Draw Venn diagram
venn_plot_70 <- draw.pairwise.venn(
  area1 = length(biomart_genes),
  area2 = length(ohnologue_genes),
  cross.area = length(intersect(biomart_genes, ohnologue_genes)),
  category = c("BioMart", "Ohnologue"),
  fill = c("skyblue", "orange"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-20, 20)
)

# Add title above the Venn diagram
grid.text("Gene Overlap Between BioMart (70% Similarity) and Ohnologue Sets", 
          x = 0.5, y = 0.95, gp = gpar(fontsize = 12, fontface = "bold"))

# ------------------------------------------------------------------------------
# Scatter plot: compare number of paralogues per gene between sources
# ------------------------------------------------------------------------------

# Summarize frequency of gene counts per (n_biomart, n_ohnologue) combination
plot_data <- merged_counts_70 %>%
  group_by(n_biomart, n_ohnologue) %>%
  summarise(gene_count = n(), .groups = "drop")

# Plot comparison as scatter plot
ggplot(plot_data, aes(x = n_biomart, y = n_ohnologue, size = gene_count)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Paralogue Count Comparison: nbiomart (>70% Similarity) vs nohnologue",
    x = "nbiomart paralogue count",
    y = "nohnologue paralogue count",
    size = "Overlapping genes"
  ) +
  theme_minimal()


# ---------------------------- DRAFT ANALYSIS -----------------------------------
# Compare orthologue datasets using full BioMart and Ohnologue data (not threshold-based)
# ------------------------------------------------------------------------------

# Load ohnologs dataset (relaxed version)
ohnologs_relaxed <- read_delim("C:/Users/HP-ssd/Desktop/Short term project/ohnologs/hsapiens.Pairs.Relaxed.2R.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

# Filter ohnologs to only include FUSIL genes (TODO: Confirm necessity of this filter)
ohnologs_relaxed_filtered <- ohnologs_relaxed %>%
  filter(Symbol1 %in% fusil_m_gene$gene_symbol) %>%
  dplyr::select(3, 4) %>%
  rename("gene_symbol" = "Symbol1") %>%
  rename("gene_paralogue" = "Symbol2")

# Prepare BioMart paralogue dataset with % identity included
biomart_paralogue <- human_gene_paralogues %>%
  filter(gene_symbol %in% protein_coding_genes$symbol) %>%
  na.omit() %>%
  dplyr::select(1, 2, 3) %>%
  rename("gene_paralogue" = "hsapiens_paralog_associated_gene_name")

# Optional: write to CSV if needed
# write.csv(biomart_paralogue, "C:/Users/HP-ssd/Desktop/biomart_paralogue.csv")

# Count paralogues per gene from BioMart dataset, grouped by % identity
count_only_in_biomart <- biomart_paralogue %>%
  group_by(gene_symbol, hsapiens_paralog_perc_id) %>%
  tally() %>%
  rename("n_biomart" = "n")

# Count paralogues per gene from ohnolog dataset
count_only_in_ohnologues <- ohnologs_relaxed_filtered %>%
  group_by(gene_symbol) %>%
  tally() %>%
  rename("n_ohnologue" = "n")

# ------------------------------------------------------------------------------
# Join BioMart and Ohnologue paralogue counts
# ------------------------------------------------------------------------------

joined_paralogue <- count_only_in_biomart %>%
  full_join(count_only_in_ohnologues, by = "gene_symbol") %>%
  mutate_all(~replace_na(., 0))  # Fill NAs with 0

# Add binary flags for presence in each dataset
joined_paralogue$in_biomart <- as.numeric(as.integer(joined_paralogue$n_biomart > 0))
joined_paralogue$in_ohnologue <- as.numeric(as.integer(joined_paralogue$n_ohnologue > 0))

# Optional: write to CSV if needed
# write.csv(joined_paralogue, "C:/Users/HP-ssd/Desktop/joined_paralogue.csv")

# ------------------------------------------------------------------------------
# Prepare data for density plot
# ------------------------------------------------------------------------------

# Convert binary presence flags to factor labels
joined_paralogue$in_biomart <- factor(joined_paralogue$in_biomart, 
                                      levels = c(0, 1), 
                                      labels = c("Not in Biomart", "In Biomart"))

joined_paralogue$in_ohnologue <- factor(joined_paralogue$in_ohnologue, 
                                        levels = c(0, 1), 
                                        labels = c("Not in Ohnologues", "In Ohnologues"))

# Define combined presence group
joined_paralogue$group <- case_when(
  joined_paralogue$in_biomart == "In Biomart" & joined_paralogue$in_ohnologue == "In Ohnologues" ~ "In Both",
  joined_paralogue$in_biomart == "In Biomart" & joined_paralogue$in_ohnologue == "Not in Ohnologues" ~ "Biomart Only",
  joined_paralogue$in_biomart == "Not in Biomart" & joined_paralogue$in_ohnologue == "In Ohnologues" ~ "Ohnologue Only",
  TRUE ~ "In Neither"
)

# Reorder factor levels for plotting
joined_paralogue$group <- factor(joined_paralogue$group,
                                 levels = c("In Neither", "Biomart Only", "Ohnologue Only", "In Both"))

# ------------------------------------------------------------------------------
# Density plot: Paralog % Identity across gene groups
# ------------------------------------------------------------------------------

ggplot(joined_paralogue, aes(x = hsapiens_paralog_perc_id)) +
  geom_density(aes(color = group), linewidth = 1) +
  labs(
    title = "Density Plot of Paralog % Identity",
    subtitle = "Solid = BioMart, Dashed = Ohnologue",
    x = "Paralog % Identity",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

# ------------------------------------------------------------------------------
# Density plot after removing 0% identity entries
# ------------------------------------------------------------------------------

joined_paralogue_clean <- joined_paralogue %>%
  filter(hsapiens_paralog_perc_id > 0)

ggplot(joined_paralogue_clean, aes(x = hsapiens_paralog_perc_id)) +
  geom_density(aes(color = group), linewidth = 1) +
  labs(
    title = "Density Plot of Paralog % Identity",
    subtitle = "Solid = BioMart, Dashed = Ohnologue",
    x = "Paralog % Identity",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )
