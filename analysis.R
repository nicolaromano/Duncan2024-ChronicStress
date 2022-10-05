# This is the analysis pipeline for the data in Duncan et al. 2022
# This is largely based on Love et al 2014 - Beginner’s guide to using the DESeq2 package

library(Rsamtools)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(GenomicFeatures)
library(biomaRt)
library(genefilter)

library(DESeq2)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pheatmap)
library(sva)

############ VERY SLOW - NEED TO DO THIS ONLY ONCE TO GENERATE OVERLAPS.RDS ###################
# summarizeOverlaps works in batches, of size yieldSize
bamfilenames <- dir("../Exp1/Aligned/GRCm39r107/", "bam$", full.names = TRUE)

bfl <- BamFileList(bamfilenames, yieldSize = 5e6)

# Exons by gene
annot <- makeTxDbFromGFF("~/Storage/Genomes/Mus_musculus.GRCm39.107.gtf", format="gtf")
exonsByGene <- exonsBy(annot, by="gene");
overlaps <- summarizeOverlaps(exonsByGene, bfl, singleEnd = FALSE,
                              ignore.strand = TRUE, fragments = TRUE)
saveRDS(overlaps, "outs/overlaps_exp1.rds")



bamfilenames <- dir("../Exp2/Aligned/", "bam$", full.names = TRUE)

bfl <- BamFileList(bamfilenames, yieldSize = 5e6)

# Exons by gene
annot <- makeTxDbFromGFF("~/Storage/Genomes/Mus_musculus.GRCm39.107.gtf", format="gtf")
exonsByGene <- exonsBy(annot, by="gene");
overlaps <- summarizeOverlaps(exonsByGene, bfl, singleEnd = FALSE,
                              ignore.strand = TRUE, fragments = TRUE)
saveRDS(overlaps, "outs/overlaps_exp2.rds")

###############################################################################################

# Save ENSEMBL/MGI gene names and descriptions into a data frame so that we don't need 
# to query ENSEMBL each time we want to convert from one to the other

# Convert Ensemble ids to gene names
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
grep("mgi", listAttributes(mart)$name, value = T)
genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
               values = rownames(overlaps), mart = mart)

genes <- genes[-which(duplicated(genes$ensembl_gene_id)),]
saveRDS(genes, "outs/genenames.rds")

###############################################################################################

# Or directly read the rds files
overlaps_exp1 <- readRDS("outs/overlaps_exp1.rds")
overlaps_exp2 <- readRDS("outs/overlaps_exp2.rds")
genes <- readRDS("outs/genenames.rds")

# Convenience functions to convert between ENSEMBL and MGI IDs

get_ensembl_from_mgi <- function(mgi) {
  genes$ensembl_gene_id[match(mgi, genes$mgi_symbol)]
}

get_mgi_from_ensembl <- function(ensembl_id) {
  res <- genes$mgi_symbol[match(ensembl_id, genes$ensembl_gene_id)]
  
  # Sometimes the match is empty or NA
  res[is.na(res)] <- ensembl_id[is.na(res)]
  res[nchar(res) == 0] <- ensembl_id[nchar(res) == 0]
  
  res
}

# Get expression matrices
mtx1 <- assay(overlaps_exp1)

# Rename columns for clarity
# e.g. from CM-1Aligned.sortedByCoord.out.bam to CM-1-exp_1
sapply(colnames(mtx1), function(x){
  group <- strsplit(x, split = "-")[[1]][1]
  num <- sub(".*-([0-9])Aligned.*", "\\1", x)
  
  paste(group, num, "exp_1", sep="-")
}) -> colnames(mtx1)


mtx2 <- assay(overlaps_exp2)

# Rename columns for clarity
group <- substr(colnames(mtx2), 1, 1)
group <- ifelse(group == "A", "CM", ifelse(group == "B", "CM_aged", "CSR_12wk"))
colnames(mtx2) <- paste(group, substr(colnames(mtx2), 2, 2), "exp_2", sep="-")

head(mtx1)
head(mtx2)

# Discard any gene that is not common
common_genes_e1e2 <- intersect(rownames(mtx1), rownames(mtx2))
counts <- cbind(mtx1[common_genes_e1e2,], mtx2[common_genes_e1e2,])

head(counts)
batch <- c(rep(1, ncol(mtx1)), rep(2, ncol(mtx2)))

data.frame(col = colnames(counts)) %>% 
  mutate(Group = factor(rep(c("CTRL", "CS", "CS-R4", 
                              "CTRL", "Aged", "CS-R12"), 
                            c(3, 3, 3, 
                              3, 4, 3)))) %>% 
  select(-col) -> metadata

metadata$Group <- factor(metadata$Group, levels = c("CTRL", "Aged",
                                                    "CS", "CS-R4",  
                                                    "CS-R12"))

##### ComBat-seq batch-correction #####
adjusted_counts <- ComBat_seq(counts, batch = factor(batch), covar_mod = metadata)

dds <- DESeqDataSetFromMatrix(adjusted_counts,
                              colData = cbind(metadata, 
                                              Batch=factor(batch),
                                              Replicate=factor(sapply(colnames(counts), function(x){
                                                strsplit(x, "-")[[1]][2]}))),
                              design = ~ Batch + Replicate + Group)

nrow(dds) # ~57k genes
# Only keep genes with a count of 10 or higher in at least 2 samples
dds <- dds[rowSums(counts(dds) >= 10) >= 2,]
nrow(dds) # ~23k genes left

# Remove pseudogenes (16291 genes left)
data.frame(ensembl = rownames(dds),
           mgi = get_mgi_from_ensembl(rownames(dds))) %>% 
  filter(!str_detect(mgi, "Rik$")) %>% 
  filter(!str_detect(mgi, "^Gm")) %>% 
  filter(!str_detect(mgi, "Rik$")) %>% 
  filter(!str_detect(mgi, "-ps[0-9]+$")) %>% 
  select(ensembl) -> filtered_genes

dds <- dds[filtered_genes$ensembl,]

# Apply variance-stabilising transform (for PCA calculation)
dds.vst <- vst(dds, blind = FALSE)

# Principal component analysis 
png("PCA_by_batch.png", res=150, width = 1200, height = 1200)
plotPCA(dds.vst, intgroup = "Batch") +
  ylim(-20, 40) +
  xlim(-20, 40) +
  labs(col = "Batch") +
  ggtitle("Principal component analysis") +
  theme(title = element_text(size=16),
        axis.title = element_text(size=14),
        axis.text = element_text(size=13),
        legend.text = element_text(size=14),
        legend.title = element_text(size=13))
dev.off()

png("PCA_by_group.png", res=150, width = 1200, height = 1200)
plotPCA(dds.vst, intgroup = "Group") +
  ylim(-20, 40) +
  xlim(-20, 40) +
  labs(col = "Group") +
  ggtitle("Principal component analysis") +
  theme(title = element_text(size=16),
        axis.title = element_text(size=14),
        axis.text = element_text(size=13),
        legend.text = element_text(size=14),
        legend.title = element_text(size=13))
dev.off()

# Perform DE analysis by
# - estimation of size factors
# - estimation of dispersion
# - Negative Binomial GLM fitting and Wald statistics
dds <- DESeq(dds)
contr <- resultsNames(dds)

contr <- contr[c(2, 6:9)] # Batch 2 vs 1 + all groups vs CTRL

# Get pairwise comparisons with log2FC shrinkage
DE_res <- lapply(contr, function(ctr){
  res <- results(dds, name = ctr, alpha = .05)
  res <- lfcShrink(dds, res = res, 
                   coef = ctr, type = "apeglm")
  res
})

# Set the list names to the contrasts for ease of use
names(DE_res) <- contr

# MA plots
MA_plots <- lapply(names(DE_res), function(x){
  res <- DE_res[[x]] %>% 
    as.data.frame() %>% 
    arrange(padj < 0.05)
  
  tb <- table(res$padj < 0.05, 
              sign(res$log2FoldChange))
  
  g <- res %>% 
    as.data.frame() %>% 
    rownames_to_column("ENSEMBL_id") %>% 
    mutate(MGI_name = get_mgi_from_ensembl(ENSEMBL_id)) %>% 
    mutate(padj = replace_na(padj, 1)) %>% 
    ggplot(aes(x = log10(baseMean), y = log2FoldChange)) +
    geom_point(aes(col = ifelse(padj < 0.05,
                                ifelse(log2FoldChange < 0, 
                                       "sigdown", "sigup"), "nonsig"))) +
    scale_color_manual(values = c("lightgray",  # Not changed
                                  "#ed2f2f", # Downregulated
                                  "#80c123") # Upregulated
                       ) +
    xlab(expression(log[10]~"(mean expression)")) +
    ylab(expression(log[2]~"fold change")) +
    annotate(geom = "text", x = Inf, y = Inf, 
             label = paste(tb[2, 2], "genes up\n",
                           tb[2, 1], "genes down"),
             hjust = 1, vjust = 1) +
    ggtitle(gsub("(_|Group_)", " ", x)) +
    ylim(-15, 15) +
    theme(legend.position = "none")
})

do.call("grid.arrange", c(MA_plots, ncol=3))

png("MA_plots.png", res = 250, width = 2400, height = 1600)
do.call("grid.arrange", c(MA_plots, ncol=3))
dev.off()


## Clustering of expression pattern
expr <- as.data.frame(assay(dds))
rownames(expr) <- make.unique(get_mgi_from_ensembl(rownames(expr)))

# Get only DE genes (in any comparison)
sigCS <- DE_res$Group_CS_vs_CTRL %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames() %>% 
  get_mgi_from_ensembl()

sigCSR4 <- DE_res$Group_CS.R4_vs_CTRL %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames() %>% 
  get_mgi_from_ensembl()

sigCSR12 <- DE_res$Group_CS.R12_vs_CTRL %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames() %>% 
  get_mgi_from_ensembl()

# All DE genes (1242 genes) in all comparisons
allDEgenes <- union(union(sigCS, sigCSR4), sigCSR12)

# z-score expression
z_score <- function(x) {
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

expr_z_score <- t(apply(expr[allDEgenes,], 1, z_score))

# Calculate mean expression in the different groups
meanExpr <- data.frame(CTRL = expr_z_score %>% 
                         as.data.frame() %>% 
                         select(starts_with("CM-")) %>% 
                         rowMeans(na.rm = TRUE),
                       CS = expr_z_score %>% 
                         as.data.frame() %>% 
                         select(starts_with("CS-")) %>% 
                         rowMeans(na.rm = TRUE),
                       CSR4 = expr_z_score %>% 
                         as.data.frame() %>% 
                         select(starts_with("CSR-")) %>% 
                         rowMeans(na.rm = TRUE),
                       CSR12 = expr_z_score %>% 
                         as.data.frame() %>% 
                         select(starts_with("CSR_12wk")) %>% 
                         rowMeans(na.rm = TRUE))

n_cutree <- 4
leg_br <- seq(-1.5, 1.5, 0.5)

png("DE_heatmap.png", height = 1500, width = 1500, res = 150)
ph <- pheatmap(meanExpr, cluster_cols = FALSE, show_rownames = FALSE,
               main = "All DE genes",
               clustering_method = "ward.D2",
               cutree_rows = n_cutree, treeheight_row = 100,
               legend_breaks = leg_br)
dev.off()

ordered_genes <- rownames(expr_z_score[ph$tree_row$order,])

ct <- cutree(ph$tree_row, n_cutree)

plot_profile <- function(group)
{
  meanExpr[names(ct[ct==group]),] %>% 
    as.data.frame() %>% 
    rownames_to_column("Gene") %>% 
    pivot_longer(cols=-Gene, names_to = "Group", values_to = "ZScore") %>% 
    mutate(Group = factor(Group, levels = c("CTRL", "CS", "CSR4", "CSR12"))) -> expr_long
  
  g <- ggplot(expr_long, aes(Group, ZScore, group=Gene)) +
    geom_line(col = "lightgray", alpha = 0.5) +
    geom_point(size = 0.1, alpha = 0.3) +
    stat_summary(fun = mean, geom = "line", col ="red", 
                 size = 1, group = 1) +
    ggtitle(paste(length(unique(expr_long$Gene)), "genes")) +
    theme_minimal()+
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=12))
  
  g
}

png("expression_profiles.png", height = 1500, width = 1500, res = 150)
do.call("grid.arrange", c(lapply(1:4, plot_profile), ncol=2))
dev.off()

