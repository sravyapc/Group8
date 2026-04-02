# install packages (only needed once)
install.packages("readxl")
install.packages("pheatmap")

# load packages (run every session)
library(readxl)
library(pheatmap)
library(dplyr)

# ---------------------------
# LOAD AND CLEAN DATA
# ---------------------------
data <- read_excel("C:/Users/sravy_zbboaem/Downloads/GSE61271_gene_names_and_GSM_headers.xlsx")
data <- as.data.frame(data)


# fix missing gene names
data$Gene <- as.character(data$Gene)
data$Gene[is.na(data$Gene) | data$Gene == ""] <- paste0("Unknown_", which(is.na(data$Gene) | data$Gene == ""))

data <- data[!grepl("^ENSSSCG", data$Gene, ignore.case = TRUE), ]

# set row names
rownames(data) <- make.unique(data$Gene)

# remove non-numeric columns
data$Gene <- NULL
data$Gene_name <- NULL

# convert to numeric
data[] <- lapply(data, function(x) as.numeric(as.character(x)))

# remove bad rows
data <- data[complete.cases(data), ]






# remove zero-variance genes
row_var <- apply(data, 1, var, na.rm = TRUE)
data <- data[row_var > 0, ]

# ---------------------------
# SELECT TOP VARIABLE GENES
# ---------------------------
gene_variance <- apply(data, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:50]
heatmap_data <- data[top_genes, ]

# ---------------------------
# LOAD METADATA
# ---------------------------
metadata <- read.csv("C:/Users/sravy_zbboaem/Downloads/GSE61271_sample_metadata.csv")

# match metadata to samples
metadata2 <- metadata[match(colnames(heatmap_data), metadata$geo_accession), ]

# create annotation
annotation_col <- data.frame(
  BiologicalState = metadata2$biological_state,
  Sex = metadata2$Sex
)

rownames(annotation_col) <- colnames(heatmap_data)

# ---------------------------
# PLOT HEATMAP
# ---------------------------
pheatmap(
  as.matrix(heatmap_data),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_colnames = TRUE,
  angle_col = 45,
  fontsize_col = 7,
  main = "Top 50 Most Variable Genes"
)

