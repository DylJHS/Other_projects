########################################################################################
### IGNORE THIS SECTION
#######################################################################################
#Load the packages

library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(tidyr)
library(tools)
library(readr)
#######################################################################################
#######################################################################################




## --- LOOK HERE ----
#######################################################################################
### SET THE INPUT BELOW
#######################################################################################

# Set the data file name
data_name <- "SummaryData_Volcanoplot__CBX1.T_vs_NLS.T_noRibosome.txt"

# The analysis will look at the fold change of cond_1 over cond_2 (log2(cond_1/cond_2))
cond_1 <- "CBX1"
cond_2 <- "NLS"

# Define the genes of interest
goi_cellulo <- c(
  "PRDM2", "ZNF584", "C5orf24", "SP140L", "SCMH1", "CBX4"
)
goi_vitro <- c(
  "ZNF384", "TEAD3", "ZNF219", "SMARCAL1", "PRDM10"
)

# Define the controls
cntrls <- c("POGZ", "MPHOSPH8", "ZNF644", "CBX1")

# Define the known interactors
interactors <- c(
  "MIS12", "TRIM28", "EHMT2", "SETDB1", "ADNP2", "LBR"
)

# Set the pvalue cutoff
p_val_cutoff <- 0.05

# Set the log 2 fold change cutoff
fc_cutoff <- 4

#######################################################################################
#######################################################################################




#######################################################################################
### CAN IGNORE THE REST
#######################################################################################
# Set the working directory to the location of this script
this_file <- normalizePath(sys.frame(1)$ofile)
setwd(dirname(this_file))

# Load the data
# Check that the data file exists
if (!file.exists(data_name)) {
  stop("Data file does not exist.")
}

# get the extension type
suffix <- tools::file_ext(data_name)
# Load the data based on the file type
if (suffix == "csv"){
  ori_data <- read.csv(data_name, header = TRUE)
} else if (suffix == "txt"){
  ori_data <- read_delim(data_name, delim = "\t", col_names = TRUE, show_col_types = FALSE)
}

# Configure the data
all_interest <- c(goi_cellulo, goi_vitro, interactors, cntrls)

# Set the columns to keep
keep_columns <- c(
  c("Gene.names", "p_value", "log2FC"), 
  names(select(ori_data, starts_with(cond_1))), 
  names(select(ori_data, starts_with(cond_2)))
)

# Reduce data based on the cols to keep
data <- ori_data %>%
  mutate(p_value = 10^(-ori_data$`pvalue.-log10`)) %>% 
  rowwise() %>% 
  mutate(
    log2FC = FoldChange
  )%>% 
  ungroup() %>%
  select(
    all_of(keep_columns)
  ) 

# Set all columns but gene names to numeric
data[2:ncol(data)] <- lapply(data[2:ncol(data)], function(x) as.numeric(as.character(x)))

# order the data 
gene_flags <- data$Gene.names %in% goi_cellulo | 
  data$Gene.names %in% cntrls |
  data$Gene.names %in% interactors |
  data$Gene.names %in% goi_vitro
data <- rbind(
  data[!(gene_flags), ],
  data[gene_flags, ]
)

# Identify all significant genes
significant_genes <- data %>%
  filter(
    p_value <= p_val_cutoff & abs(log2FC) >= fc_cutoff
  ) %>% 
  pull(Gene.names)

# Map colours
keyvals.colour <- ifelse(data$Gene.names %in% goi_cellulo, '#9C0824',
                         ifelse(data$Gene.names %in% goi_vitro, '#FA9081',
                                ifelse(data$Gene.names %in% interactors, '#FF7F00',
                                       ifelse(data$Gene.names %in% cntrls, '#0000FF',
                                              ifelse(data$Gene.names %in% significant_genes, '#09622A', 'grey')
                                       )
                                )
                         )
)

names(keyvals.colour) <- ifelse(data$Gene.names %in% goi_cellulo, 'Gene of Interest (In cellulo)',
                                ifelse(data$Gene.names %in% goi_vitro, 'Gene of Interest (In vitro)',
                                       ifelse(data$Gene.names %in% interactors, 'Known Interactor',
                                              ifelse(data$Gene.names %in% cntrls, 'Control',
                                                     ifelse(data$Gene.names %in% significant_genes, 'Significant Gene', 'Other')
                                              )
                                       )
                                )
)



# Map sizes
keyvals.size <- ifelse(
  data$Gene.names %in% all_interest, 3, 1.5
)


# Volcano plot
xmin <- min(data$log2FC, na.rm = TRUE) - 0.5
xmax <- max(data$log2FC, na.rm = TRUE) + 0.5
ymax <- max(-log10(data$p_value), na.rm = TRUE) + 1

vol_plt <- EnhancedVolcano(
  data,
  lab = ifelse(data$Gene.names %in% all_interest, data$Gene.names, ""),
  selectLab = all_interest,
  x = 'log2FC',
  y = 'p_value',
  ylim =  c(0, ymax),
  xlim = c(xmin, xmax),
  subtitle = "",
  pCutoff = p_val_cutoff,
  FCcutoff = fc_cutoff,
  title = paste0("Differential Protein Abundance: ", cond_1, " vs ", cond_2),
  colCustom = keyvals.colour,
  pointSize = keyvals.size,
  caption = paste0("total = ", nrow(data), " genes"),
  legendPosition = "right",
  legendLabSize = 12,
  drawConnectors = TRUE,
  widthConnectors = 0.45,
  colConnectors = 'black',
  arrowheads = FALSE,
)

# Label the sides of the plot
vol_plt <- vol_plt +
  annotate(
    geom = "text",
    x = c(xmin + 0.2, xmax - 0.2),
    y = ymax - 0.5,
    label = c(cond_2, cond_1),
    size = 8,
    fontface = "bold",
    color = c("#A4365E", "#4A3D8D"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "italic"),
    plot.subtitle = element_text(hjust = 0.5, size = 16)
  )

# Save as pdf
ggsave(
  filename = paste(cond_1, cond_2, "volcano_plot.pdf", sep = "_"),
  plot = vol_plt,
  device = "pdf",
  width = 14,
  height = 8,
  dpi = 300
)

rm(list = ls())
