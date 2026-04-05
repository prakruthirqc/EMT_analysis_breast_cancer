# ==========================================
# Transcriptomic Analysis of TGF-beta Signaling and EMT
# Breast Cancer Cell Line Dataset (Heiser 2012)
# Author: Prakruthi R
# ==========================================

# Set working directory
setwd("C:/bioinfo_project/data")

# Load data
expr <- read.delim("gene_expression.txt", row.names = 1, check.names = FALSE)
pheno <- read.delim("phenotype.txt", check.names = FALSE)

# Check data
head(expr[,1:5])
head(pheno)
dim(expr)
dim(pheno)

# Genes of interest
genes <- c("TGFB1","TGFBR1","SMAD2","SMAD3","SMAD4",
           "CDH1","VIM","SNAI1","TWIST1","ZEB1","MKI67","BCL2")

# Extract gene expression
gene_expr <- expr[rownames(expr) %in% genes, ]
gene_expr <- t(gene_expr)
gene_expr <- as.data.frame(gene_expr)

head(gene_expr)
colnames(pheno)

# Merge phenotype data
gene_expr$CellLine <- rownames(gene_expr)

gene_expr <- merge(gene_expr,
                   pheno[, c("sampleName_ICR", "Sub_type_ICR")],
                   by.x = "CellLine",
                   by.y = "sampleName_ICR")

head(gene_expr)

# Load libraries
library(ggplot2)
library(pheatmap)

# ==========================================
# Boxplots for Gene Expression
# ==========================================

ggplot(gene_expr, aes(x=Sub_type_ICR, y=TGFB1, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("TGFB1 Expression Across Breast Cancer Subtypes")
ggsave("TGFB1_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=TGFBR1, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("TGFBR1 Expression Across Breast Cancer Subtypes")
ggsave("TGFBR1_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=SMAD2, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("SMAD2 Expression Across Breast Cancer Subtypes")
ggsave("SMAD2_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=SMAD3, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("SMAD3 Expression Across Breast Cancer Subtypes")
ggsave("SMAD3_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=SMAD4, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("SMAD4 Expression Across Breast Cancer Subtypes")
ggsave("SMAD4_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=CDH1, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("CDH1 Expression Across Breast Cancer Subtypes")
ggsave("CDH1_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=VIM, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("VIM Expression Across Breast Cancer Subtypes")
ggsave("VIM_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=SNAI1, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("SNAI1 Expression Across Breast Cancer Subtypes")
ggsave("SNAI1_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=TWIST1, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("TWIST1 Expression Across Breast Cancer Subtypes")
ggsave("TWIST1_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=ZEB1, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("ZEB1 Expression Across Breast Cancer Subtypes")
ggsave("ZEB1_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=MKI67, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("MKI67 Expression Across Breast Cancer Subtypes")
ggsave("MKI67_boxplot.png")

ggplot(gene_expr, aes(x=Sub_type_ICR, y=BCL2, fill=Sub_type_ICR)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("BCL2 Expression Across Breast Cancer Subtypes")
ggsave("BCL2_boxplot.png")

# ==========================================
# EMT Heatmap
# ==========================================

emt <- gene_expr[, c("CDH1","VIM","SNAI1","TWIST1","ZEB1")]

pheatmap(emt,
         scale = "row",
         filename = "EMT_heatmap.png")

# ==========================================
# Volcano Plot (Basal vs Luminal)
# ==========================================

subset_data <- gene_expr[gene_expr$Sub_type_ICR %in% c("BASAL","Luminal"), ]

subset_data$Group <- ifelse(subset_data$Sub_type_ICR=="BASAL","Basal","Luminal")

volcano_genes <- c("TGFB1","SMAD2","SMAD3","SMAD4","CDH1","VIM","SNAI1","TWIST1","ZEB1","MKI67","BCL2")

results <- data.frame(Gene=character(), logFC=numeric(), pvalue=numeric())

for (g in volcano_genes) {
  basal <- subset_data[subset_data$Group=="Basal", g]
  luminal <- subset_data[subset_data$Group=="Luminal", g]
  
  basal <- as.numeric(basal)
  luminal <- as.numeric(luminal)
  
  basal <- basal[!is.na(basal)]
  luminal <- luminal[!is.na(luminal)]
  
  if(length(basal) >= 2 & length(luminal) >= 2){
    ttest <- t.test(basal, luminal)
    logFC <- mean(basal) - mean(luminal)
    
    results <- rbind(results, data.frame(Gene=g, logFC=logFC, pvalue=ttest$p.value))
  }
}

print(results)

png("volcano_plot.png")
plot(results$logFC, -log10(results$pvalue),
     pch=19,
     xlab="Log Fold Change (Basal - Luminal)",
     ylab="-Log10 P-value",
     main="Volcano Plot: Basal vs Luminal")

abline(h=-log10(0.05), col="red")
text(results$logFC, -log10(results$pvalue), labels=results$Gene, pos=3, cex=0.8)

dev.off()