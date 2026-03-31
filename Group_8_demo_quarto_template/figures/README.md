# Figures

````{r}`

install.packages("readxl")
install.packages("pheatmap")

library(readxl)
library(pheatmap)
library(dplyr)


data <- read_excel("C:/Users/sravy_zbboaem/Downloads/GSE61271_gene_names_and_GSM_headers.xlsx")
data <- as.data.frame(data)


data$Gene <- as.character(data$Gene)
data$Gene[is.na(data$Gene) | data$Gene == ""] <- paste0("Unknown_", which(is.na(data$Gene) | data$Gene == ""))

data <- data[!grepl("^ENSSSCG", data$Gene, ignore.case = TRUE), ]

rownames(data) <- make.unique(data$Gene)

data$Gene <- NULL
data$Gene_name <- NULL

data[] <- lapply(data, function(x) as.numeric(as.character(x)))

data <- data[complete.cases(data), ]

row_var <- apply(data, 1, var, na.rm = TRUE)
data <- data[row_var > 0, ]

gene_variance <- apply(data, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:50]
heatmap_data <- data[top_genes, ]

metadata <- read.csv("C:/Users/sravy_zbboaem/Downloads/GSE61271_sample_metadata.csv")

metadata2 <- metadata[match(colnames(heatmap_data), metadata$geo_accession), ]

annotation_col <- data.frame(
  BiologicalState = metadata2$biological_state,
  Sex = metadata2$Sex
)

rownames(annotation_col) <- colnames(heatmap_data)


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

```` `

Export final figures here (png/pdf/svg).  
Use consistent names, e.g.:
- `fig01_flowchart.png`
- `fig02_eda_histograms.png`
- `fig03_model_roc.png`
