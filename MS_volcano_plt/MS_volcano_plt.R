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
goi <- c(
  "PRDM2", "ZNF584", "C5orf24", "SP140L", "SCMH1", "CBX4", 
  "SMARCAL1", "TEAD3", "PRDM10", "ZNF219", "ZNF644", "ZNF384"
)

# Define the controls
cntrls <- c("POGZ", "MPHOSPH8")

# Set the pvalue cutoff
p_val_cutoff <- 0.05

# Set the log 2 fold change cutoff
fc_cutoff <- 0.25

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
all_interest <- c(goi, cntrls)

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
    log2FC = log2(mean(c_across(starts_with(cond_1))) / mean(c_across(starts_with(cond_2))))
  ) %>% 
  ungroup() %>%
  select(
    all_of(keep_columns)
  ) 

# Set all columns but gene names to numeric
data[2:ncol(data)] <- lapply(data[2:ncol(data)], function(x) as.numeric(as.character(x)))

# order the data 
gene_flags <- data$Gene.names %in% goi | 
  data$Gene.names %in% cntrls
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

# Map colours: red for genes of interest, grey otherwise
keyvals.colour <- ifelse(data$Gene.names %in% goi, '#920000',
                         ifelse(data$Gene.names %in% cntrls, '#0000FF',
                                ifelse(data$Gene.names %in% significant_genes, '#09622A', 'grey')
                         )
)

names(keyvals.colour) <- ifelse(data$Gene.names %in% goi, 'Gene of Interest',
                                ifelse(data$Gene.names %in% cntrls, 'Control',
                                       ifelse(data$Gene.names %in% significant_genes, 'Significant', 'Not Significant')
                                )
)

# Map sizes: 4 for genes of interest, 1 for others
keyvals.size <- ifelse(
  data$Gene.names %in% all_interest, 3, 1.5
)

names(keyvals.size) <- ifelse(
  data$Gene.names %in% goi, 'Gene of Interest', 
  ifelse(data$Gene.names %in% cntrls, 'Control', 'Other')
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
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16)
  )

# Save as pdf
ggsave(
  filename = paste(cond_1, cond_2, "volcano_plot.pdf", sep = "_"),
  plot = vol_plt,
  device = "pdf",
  width = 12,
  height = 8,
  dpi = 300
)