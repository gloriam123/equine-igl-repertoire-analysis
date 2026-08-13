#script para calcular a freq genica dos clones e fornecer isso em diferentes arquivos do excel

install.packages("tidyverse")

library(data.table)
library(tidyverse)
library(dplyr)
library(openxlsx)

meus_arquivos <- "/Users/Gloria/Desktop/Tese/IgBLAST/anotados_IgBLAST_T2T/clones_T2T"

arquivos_clonotipados <- sort(list.files(
  path = meus_arquivos,
  recursive = TRUE,
  all.files = TRUE,
  pattern = "YClon_clonotyped\\.tsv$",
  full.names = FALSE
))

for (y in seq_along(arquivos_clonotipados)) {
  print(arquivos_clonotipados[y])
  
  CLONOTIPAGEM <- fread(
    file.path(meus_arquivos, arquivos_clonotipados[y]),
    header = TRUE,
    sep = "\t"
  )
  
  CLONOTIPAGEM <- CLONOTIPAGEM %>%
    select(sequence_id, v_call, v_identity, j_call, j_identity)
  
  # Frequência V
  ValorVabsoluto <- as.data.frame(table(CLONOTIPAGEM$v_call, useNA = "ifany"))
  colnames(ValorVabsoluto) <- c("v_call", "Freq_Abs")
  
  ValorVrelativo <- as.data.frame(prop.table(table(CLONOTIPAGEM$v_call, useNA = "ifany")))
  colnames(ValorVrelativo) <- c("v_call", "Freq_Rel")
  
  TabelaGeneV <- merge(ValorVabsoluto, ValorVrelativo, by = "v_call")
  
  # Frequência J
  ValorJabsoluto <- as.data.frame(table(CLONOTIPAGEM$j_call, useNA = "ifany"))
  colnames(ValorJabsoluto) <- c("j_call", "Freq_Abs")
  
  ValorJrelativo <- as.data.frame(prop.table(table(CLONOTIPAGEM$j_call, useNA = "ifany")))
  colnames(ValorJrelativo) <- c("j_call", "Freq_Rel")
  
  TabelaGeneJ <- merge(ValorJabsoluto, ValorJrelativo, by = "j_call")
  
  # Identidade V
  IdentidadeV <- CLONOTIPAGEM %>%
    filter(!is.na(v_call), !is.na(v_identity)) %>%
    group_by(v_call) %>%
    summarise(
      media_v_identity = mean(v_identity, na.rm = TRUE),
      n = n()
    ) %>%
    arrange(desc(media_v_identity))
  
  # Identidade J
  IdentidadeJ <- CLONOTIPAGEM %>%
    filter(!is.na(j_call), !is.na(j_identity)) %>%
    group_by(j_call) %>%
    summarise(
      media_j_identity = mean(j_identity, na.rm = TRUE),
      n = n()
    ) %>%
    arrange(desc(media_j_identity))
  
  workbook <- createWorkbook()
  addWorksheet(workbook, "Freq_V")
  writeData(workbook, "Freq_V", TabelaGeneV)
  addWorksheet(workbook, "Freq_J")
  writeData(workbook, "Freq_J", TabelaGeneJ)
  addWorksheet(workbook, "Identidade_V")
  writeData(workbook, "Identidade_V", IdentidadeV)
  addWorksheet(workbook, "Identidade_J")
  writeData(workbook, "Identidade_J", IdentidadeJ)
  
  nome_saida <- str_replace(
    basename(arquivos_clonotipados[y]),
    "\\.tsv$",
    "_analise_frequencia_identidade.xlsx"
  )
  
  saveWorkbook(
    workbook,
    file.path(meus_arquivos, nome_saida),
    overwrite = TRUE
  )
}


list.files("/Users/Gloria/Desktop/Tese/IgBLAST/anotados_IgBLAST_T2T/clones_T2T", pattern = "xlsx$", full.names = TRUE)