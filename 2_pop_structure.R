
#========================================
# 0. Setup paths
#========================================
data_dir    <- file.path(getwd(), "Data")
results_dir <- file.path(getwd(), "Results")
fig_dir     <- results_dir
if (!dir.exists(results_dir)) dir.create(results_dir)
if (!dir.exists(fig_dir))   dir.create(fig_dir)

#========================================
# 1. Load libraries
#========================================
library(pophelper)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(igraph)
library(ggraph)
library(reshape2)
library(RColorBrewer)
library(ape)
library(phangorn)

#========================================
# 2. PCA plot
#========================================
# Read group‐to‐sample mapping
groupinfo <- read.csv(
  file.path(data_dir, "ID_match.csv"),
  stringsAsFactors = FALSE
)

# Read PCA results
eigenvec <- read.table(
  file.path(data_dir, "filtered_PCA.eigenvec"),
  header = FALSE, stringsAsFactors = FALSE
)
eigenval <- read.table(
  file.path(data_dir, "filtered_PCA.eigenval"),
  header = FALSE, stringsAsFactors = FALSE
)

# Name PCA columns: FamilyID, IndividualID, then PCs
colnames(eigenvec) <- c(
  "FamilyID", "IndividualID",
  paste0("PC", 1:(ncol(eigenvec) - 2))
)

# Compute percentage variance per PC
variance_explained <- eigenval$V1 / sum(eigenval$V1) * 100
variance_df <- data.frame(
  PC       = paste0("PC", seq_along(variance_explained)),
  Variance = variance_explained
)

# Bar plot: variance explained
p_var <- ggplot(variance_df, aes(PC, Variance)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Variance Explained by PCs",
    x     = "Principal Component",
    y     = "Percent Variance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(fig_dir, "PCA_variance.png"),
  plot   = p_var,
  width  = 5, height = 4, dpi = 300
)

# Scatter: PC1 vs PC2 (uncolored)
pca_basic <- ggplot(eigenvec, aes(PC1, PC2)) +
  geom_point(color = "grey20", size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA: PC1 vs PC2",
    x     = paste0("PC1 (", round(variance_explained[1],2), "%)"),
    y     = paste0("PC2 (", round(variance_explained[2],2), "%)")
  ) +
  theme(panel.border = element_rect(fill = NA, color = "black"))
ggsave(
  file.path(fig_dir, "PCA_basic.png"),
  plot   = pca_basic,
  width  = 5, height = 5, dpi = 300
)

# PCA colored by group
eigenvec_grouped <- left_join(
  groupinfo, eigenvec,
  by = c("sample" = "IndividualID")
)
pca_grouped <- ggplot(eigenvec_grouped, aes(PC1, PC2, color = kin)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA by Group",
    x     = paste0("PC1 (", round(variance_explained[1],2), "%)"),
    y     = paste0("PC2 (", round(variance_explained[2],2), "%)"),
    color = "Group"
  ) +
  theme(panel.border = element_rect(fill = NA, color = "black"))
ggsave(
  file.path(fig_dir, "PCA_by_group.png"),
  plot   = pca_grouped,
  width  = 5, height = 5, dpi = 300
)

# PCA with non‐overlapping labels
pca_labeled <- pca_grouped +
  geom_text_repel(aes(label = sample), size = 3)
ggsave(
  file.path(fig_dir, "PCA_labeled.png"),
  plot   = pca_labeled,
  width  = 6, height = 6, dpi = 300
)

#========================================
# 3. IBS analysis
#========================================
# Read lower‐triangular IBS matrix
ibs_lines <- readLines(file.path(data_dir, "filtered_IBS.mibs"))
ibs_rows  <- strsplit(ibs_lines, "\\s+")
maxc      <- max(lengths(ibs_rows))
ibs_low   <- t(sapply(ibs_rows, function(r) c(as.numeric(r), rep(NA, maxc - length(r)))))
ibs_mat   <- ibs_low
ibs_mat[upper.tri(ibs_mat)] <- t(ibs_low)[upper.tri(ibs_low)]
ids       <- read.table(file.path(data_dir, "filtered_IBS.mibs.id"), stringsAsFactors = FALSE)$V1
rownames(ibs_mat) <- colnames(ibs_mat) <- ids

# Dendrogram
ibs_dist <- as.dist(1 - ibs_mat)
hc       <- hclust(ibs_dist, method = "complete")
png(
  file.path(fig_dir, "IBS_dendrogram.png"),
  width = 6, height = 5, units = "in", res = 300
)
plot(hc, main = "IBS Clustering", xlab = "", sub = "")
dev.off()

# Network at threshold = 0.8
thr <- 0.8
adj <- ibs_mat > thr
diag(adj) <- FALSE
g   <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
comps <- components(g)
V(g)$Group <- factor(comps$membership)
pal        <- brewer.pal(max(comps$membership), "Set3")

p_net <- ggraph(g, layout = "fr") +
  geom_edge_link(color = "grey50") +
  geom_node_point(aes(color = Group), size = 4, alpha = 0.8) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = pal) +
  theme_void() +
  ggtitle(paste("IBS Network >", thr))
ggsave(
  file.path(fig_dir, "IBS_network.png"),
  plot   = p_net,
  width  = 6, height = 6, dpi = 300
)

#========================================
# 4. Phylogenetic tree
#========================================
dist_mat <- as.matrix(read.table(file.path(data_dir, "filtered_distance.dist"), header = FALSE))
tip_ids  <- sub("\\t.*", "", readLines(file.path(data_dir, "filtered_distance.dist.id")))
rownames(dist_mat) <- colnames(dist_mat) <- tip_ids

tree   <- nj(as.dist(dist_mat))
rooted <- root(tree, outgroup = tip_ids[1], resolve.root = TRUE)

# Color tips by group
gi_filt <- filter(groupinfo, sample %in% rooted$tip.label)
pal2    <- setNames(rainbow(length(unique(gi_filt$kin))), unique(gi_filt$kin))
tipcols <- sapply(rooted$tip.label, function(l) pal2[gi_filt$kin[gi_filt$sample == l]])

png(
  file.path(fig_dir, "phylo_tree.png"),
  width = 6, height = 6, units = "in", res = 300
)
plot.phylo(rooted, tip.color = tipcols,
           cex = 0.6, edge.width = 1,
           main = "Phylogenetic Tree")
legend("topright", legend = names(pal2), col = pal2, pch = 19, title = "Group")
dev.off()

#========================================
# 5. Admixture CV errors
#========================================
admixture_files <- list.files(
  path    = data_dir,
  pattern = "admixture_K\\d+\\.out",
  full.names = TRUE
)
cv_errors <- data.frame(K = integer(), CV = numeric())
for (f in admixture_files) {
  txt <- readLines(f)
  K   <- as.numeric(sub("admixture_K(\\d+)\\.out", "\\1", basename(f)))
  line<- grep("CV error", txt, value = TRUE)
  if (length(line)) {
    cv <- as.numeric(sub(".*: (.*)", "\\1", line))
    cv_errors <- rbind(cv_errors, data.frame(K = K, CV = cv))
  }
}
cv_errors <- cv_errors[order(cv_errors$K), ]

p_cv <- ggplot(cv_errors, aes(K, CV)) +
  geom_line(color = "blue") + geom_point() +
  theme_minimal() +
  labs(title = "Admixture CV Error", x = "K", y = "CV Error") +
  theme(panel.border = element_rect(fill = NA, color = "black"))
ggsave(
  file.path(fig_dir, "admixture_cv.png"),
  plot   = p_cv,
  width  = 6, height = 4, dpi = 300
)

#========================================
# 6. Population structure (.Q files)
#========================================
fam    <- read.table(file.path(data_dir, "plink_filtered_selected.fam"), stringsAsFactors = FALSE)
labels <- fam$V2
qfiles <- list.files(path = data_dir, pattern = "\\.Q$", full.names = TRUE)
qfiles <- qfiles[order(nchar(qfiles), qfiles)]
qlist  <- readQ(files = qfiles)
qlist  <- lapply(qlist, setNames, labels)
# after qlist <- readQ(files=...)
names(qlist) <- paste0("K", sapply(qlist, ncol))

for (k in seq_along(qlist)) {
  thisK <- names(qlist)[k]
  plotQ(
    qlist[k],                 # single‐element named list
    imgoutput  = "sep",
    clustercol = brewer.pal(lengths(qlist[k]), "Set1"),
    splab      = thisK,
    exportplot = TRUE,
    returnplot = FALSE,
    exportpath = fig_dir,
    imgtype    = "png",
    width      = 20,
    height     = 5
  )
}

