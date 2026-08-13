# What this script does:
# Uses clonotyped sequence files to calculate relative and absolute V and J gene frequencies;
# Processes all files containing "YClon_clonotyped.tsv" in their names in a loop;
# Calculates the mean V and J annotation identity for each gene;
# Exports the results to separate worksheets in an Excel file for each sample.

library(data.table)
library(tidyverse)
library(dplyr)
library(openxlsx)

# Directory containing the clonotyped TSV files
input_dir <- "path/to/clonotyped_files"

# Find all clonotyped files
clonotyped_files <- sort(list.files(
  path = input_dir,
  recursive = TRUE,
  all.files = TRUE,
  pattern = "YClon_clonotyped\\.tsv$",
  full.names = FALSE
))

# Process each clonotyped file
for (i in seq_along(clonotyped_files)) {
  
  print(clonotyped_files[i])
  
  clonotype_data <- fread(
    file.path(input_dir, clonotyped_files[i]),
    header = TRUE,
    sep = "\t"
  )
  
  clonotype_data <- clonotype_data %>%
    select(
      sequence_id,
      v_call,
      v_identity,
      j_call,
      j_identity
    )
  
  
  # ============================================================
  # V gene frequency
  # ============================================================
  
  v_absolute_frequency <- as.data.frame(
    table(clonotype_data$v_call, useNA = "ifany")
  )
  
  colnames(v_absolute_frequency) <- c(
    "v_call",
    "absolute_frequency"
  )
  
  v_relative_frequency <- as.data.frame(
    prop.table(
      table(clonotype_data$v_call, useNA = "ifany")
    )
  )
  
  colnames(v_relative_frequency) <- c(
    "v_call",
    "relative_frequency"
  )
  
  v_gene_frequency <- merge(
    v_absolute_frequency,
    v_relative_frequency,
    by = "v_call"
  )
  
  
  # ============================================================
  # J gene frequency
  # ============================================================
  
  j_absolute_frequency <- as.data.frame(
    table(clonotype_data$j_call, useNA = "ifany")
  )
  
  colnames(j_absolute_frequency) <- c(
    "j_call",
    "absolute_frequency"
  )
  
  j_relative_frequency <- as.data.frame(
    prop.table(
      table(clonotype_data$j_call, useNA = "ifany")
    )
  )
  
  colnames(j_relative_frequency) <- c(
    "j_call",
    "relative_frequency"
  )
  
  j_gene_frequency <- merge(
    j_absolute_frequency,
    j_relative_frequency,
    by = "j_call"
  )
  
  
  # ============================================================
  # V gene identity
  # ============================================================
  
  v_identity_summary <- clonotype_data %>%
    filter(
      !is.na(v_call),
      !is.na(v_identity)
    ) %>%
    group_by(v_call) %>%
    summarise(
      mean_v_identity = mean(v_identity, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_v_identity))
  
  
  # ============================================================
  # J gene identity
  # ============================================================
  
  j_identity_summary <- clonotype_data %>%
    filter(
      !is.na(j_call),
      !is.na(j_identity)
    ) %>%
    group_by(j_call) %>%
    summarise(
      mean_j_identity = mean(j_identity, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_j_identity))
  
  
  # ============================================================
  # Export results to Excel
  # ============================================================
  
  workbook <- createWorkbook()
  
  addWorksheet(workbook, "V_Gene_Frequency")
  writeData(
    workbook,
    "V_Gene_Frequency",
    v_gene_frequency
  )
  
  addWorksheet(workbook, "J_Gene_Frequency")
  writeData(
    workbook,
    "J_Gene_Frequency",
    j_gene_frequency
  )
  
  addWorksheet(workbook, "V_Gene_Identity")
  writeData(
    workbook,
    "V_Gene_Identity",
    v_identity_summary
  )
  
  addWorksheet(workbook, "J_Gene_Identity")
  writeData(
    workbook,
    "J_Gene_Identity",
    j_identity_summary
  )
  
  
  # Output file name
  output_filename <- str_replace(
    basename(clonotyped_files[i]),
    "\\.tsv$",
    "_gene_frequency_identity_analysis.xlsx"
  )
  
  
  # Save Excel workbook
  saveWorkbook(
    workbook,
    file.path(input_dir, output_filename),
    overwrite = TRUE
  )
}


# List all generated Excel files
list.files(
  input_dir,
  pattern = "xlsx$",
  full.names = TRUE
)