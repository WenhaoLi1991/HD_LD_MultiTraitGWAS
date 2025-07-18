#========================================
# 0. Load libraries and set directories
#========================================
library(qqman)
library(dplyr)
library(purrr)
library(openxlsx)

# Define data and results directories
data_dir    <- file.path(getwd(), "Data")
results_dir <- file.path(getwd(), "Results")
if (!dir.exists(results_dir)) dir.create(results_dir)

#========================================
# 1. Parameters
#========================================
expend.bp        <- 100000    # flanking window around each SNP
gwas.threshold   <- 1e-4      # p-value cutoff
sel.trait        <- "YPP"     # trait to analyze
output.dir       <- file.path(results_dir, sel.trait)
if (!dir.exists(output.dir)) dir.create(output.dir)

#========================================
# 2. Load GWAS results for selected trait
#========================================
load(file.path(
  results_dir,
  paste0("GWASDropsxE_", sel.trait, ".Rdata")
))

# Extract the GWAS table for this trait
gwas_data <- GWASDropsxE$GWAResult$phenoDat %>%
  filter(trait == sel.trait) %>%
  select(snp, chr, ps = ps, pValue = pValue)

# Replace any zero p‑values with the smallest nonzero p
min_pval <- min(gwas_data$pValue[gwas_data$pValue > 0], na.rm = TRUE)
gwas_data$pValue[gwas_data$pValue == 0] <- min_pval

# Remove missing or infinite p‑values
gwas_data <- gwas_data %>%
  filter(!is.na(pValue) & is.finite(pValue))

# Prefix chromosome with "chr"
gwas_data$chr <- paste0("chr", gwas_data$chr)

#========================================
# 3. Identify significant SNPs
#========================================
gwas_sig <- gwas_data %>% filter(pValue < gwas.threshold)

#========================================
# 4. Group overlapping significant ranges
#========================================
group_overlapping_ranges <- function(df) {
  df <- df %>%
    arrange(ps) %>%
    mutate(
      range_start = ps - expend.bp,
      range_end   = ps + expend.bp
    )
  grp      <- 0
  current_end <- -Inf
  df$group <- NA_integer_
  for (i in seq_len(nrow(df))) {
    if (df$range_start[i] > current_end) {
      grp <- grp + 1
      current_end <- df$range_end[i]
    } else {
      current_end <- max(current_end, df$range_end[i])
    }
    df$group[i] <- grp
    df$range_start[i] <- min(df$range_start[i], df$range_start[df$group == grp])
    df$range_end[i]   <- max(df$range_end[i], df$range_end[df$group == grp])
  }
  df
}

grouped_data <- gwas_sig %>%
  group_by(chr) %>%
  group_split() %>%
  map_df(group_overlapping_ranges) %>%
  mutate(group = paste0(chr, "_", group))

# Write full grouped table
write.table(
  grouped_data,
  file = file.path(output.dir, "grouped_data.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

#========================================
# 5. Export per‑group SNP lists and ranges
#========================================
unique_groups <- unique(grouped_data$group)

for (grp in unique_groups) {
  sub <- grouped_data %>% filter(group == grp)
  chr <- unique(sub$chr)
  rs  <- min(sub$range_start)
  re  <- max(sub$range_end)
  
  all_chr <- gwas_data %>% filter(chr == chr)
  sel_idx <- which.min(abs(all_chr$ps - rs)):which.min(abs(all_chr$ps - re))
  result_current <- all_chr[sel_idx, ]
  
  write.table(
    result_current %>% select(chr, ps, pValue),
    file = file.path(output.dir, paste0("gwas_pvalue_", grp, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  write.table(
    data.frame(chr = chr, range_start = rs, range_end = re),
    file = file.path(output.dir, paste0("gwas_range_", grp, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
}
