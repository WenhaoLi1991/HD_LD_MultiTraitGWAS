# Load required packages
library(ggplot2)
library(statgenGWAS)
library(statgenGxE)
library(statgenSTA)
library(statgenQTLxT)
library(dplyr)

# Clear workspace
rm(list = ls())

# Define directories
data_dir    <- file.path(getwd(), "Data")
results_dir <- file.path(getwd(), "Results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Load precomputed GWAS/Q×E data
load(file = file.path(data_dir, "alldropdata.Rdata"))

# Inspect phenotype objects
str(dropsPheno_simQTL)
str(dropsPheno)
# (Optionally use dropsPheno_simQTL in place of dropsPheno)
# dropsPheno <- dropsPheno_simQTL

# Split phenotypes by experiment
dropsPhenoList <- split(dropsPheno, dropsPheno$Experiment)

# Remove the Experiment column from each trial’s data
dropsPhenoList <- lapply(dropsPhenoList, function(df) df[ , setdiff(names(df), "Experiment")])

# Names of trials
names(dropsPhenoList)

# Build GData object: genotypes, map, and phenotypes
gDataDrops <- createGData(
  geno  = dropsMarkers,
  map   = dropsMap,
  pheno = dropsPhenoList
)

# Impute and recode markers, removing duplicates
set.seed(1234)
gDataDropsDedup <- codeMarkers(
  gDataDrops,
  impute  = TRUE,
  verbose = TRUE
)

# Define traits to analyze
alltraits <- c(
  "DtSS","DtTS","GP","YPH","YPP",
  "KSR_logit","EL","ED","ARN",
  "AKNR","HKW","TBL","EH"
)

#--------------------------------------
# High-density trial: multi-trait GWAS
#--------------------------------------
GWASDrops_high_dens <- runMultiTraitGwas(
  gData         = gDataDropsDedup,
  traits        = alltraits,
  kinshipMethod = "vanRaden",
  trials        = "high_dens",
  covModel      = "fa",
  thrType       = "bonf"
)

# Save the high-density results
save(
  GWASDrops_high_dens,
  file = file.path(results_dir, "GWASDrops_high_dens_fa_bonf.Rdata")
)

# Generate and style QTL effect plot
g_high <- plot(
  GWASDrops_high_dens,
  plotType  = "qtl",
  yThr      = 6,
  normalize = TRUE
)
effectPlot_high_dens <- g_high +
  scale_color_manual(
    name     = "Allelic effect",
    labels   = c("negative", "positive"),
    values   = c("green4", "tomato"),
    na.translate = FALSE
  ) +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(
      fill     = NA,
      color    = "black",
      size     = 0.5,
      linetype = "solid"
    )
  )

# Display and save the plot
print(effectPlot_high_dens)
ggsave(
  filename = file.path(results_dir, "effectPlot_high_dens.png"),
  plot     = effectPlot_high_dens,
  width    = 8,
  height   = 4,
  units    = "in"
)

#--------------------------------------
# Low-density trial: multi-trait GWAS
#--------------------------------------
GWASDrops_low_dens <- runMultiTraitGwas(
  gData         = gDataDropsDedup,
  traits        = alltraits,
  kinshipMethod = "vanRaden",
  trials        = "low_dens",
  covModel      = "fa",
  thrType       = "bonf"
)

# Save the low-density results
save(
  GWASDrops_low_dens,
  file = file.path(results_dir, "GWASDrops_low_dens_fa_bonf.Rdata")
)

# Generate and style QTL effect plot for low density
g_low <- plot(
  GWASDrops_low_dens,
  plotType  = "qtl",
  yThr      = 6,
  normalize = TRUE
)
effectPlot_low_dens <- g_low +
  scale_color_manual(
    name     = "Allelic effect",
    labels   = c("negative", "positive"),
    values   = c("green4", "tomato"),
    na.translate = FALSE
  ) +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(
      fill     = NA,
      color    = "black",
      size     = 0.5,
      linetype = "solid"
    )
  )

print(effectPlot_low_dens)
ggsave(
  filename = file.path(results_dir, "effectPlot_low_dens.png"),
  plot     = effectPlot_low_dens,
  width    = 8,
  height   = 4,
  units    = "in"
)

#--------------------------------------
# Extract and reshape p‑values: low density
#--------------------------------------
GWASres_low <- GWASDrops_low_dens$GWAResult$low_dens

# Build wide table of p‑values by trait
GWASres_list_low <- lapply(
  unique(GWASres_low$trait),
  function(tr) {
    df <- subset(GWASres_low, trait == tr)
    setNames(
      data.frame(pValue = df$pValue),
      paste0(tr, "_pValue")
    )
  }
)
GWASres.df.low_dens <- do.call(
  cbind,
  c(
    list(
      SNP = GWASres_low$snp[GWASres_low$trait == unique(GWASres_low$trait)[1]],
      chr = GWASres_low$chr[GWASres_low$trait == unique(GWASres_low$trait)[1]],
      pos = GWASres_low$pos[GWASres_low$trait == unique(GWASres_low$trait)[1]]
    ),
    GWASres_list_low
  )
)

# Generate Manhattan‐style plot for low density
library(CMplot)
CMplot(
  GWASres.df.low_dens[, 1:4],
  plot.type  = "m",
  r          = 0.4,
  col        = c("grey20", "grey70"),
  chr.labels = paste0("Chr", 1:10),
  threshold  = c(1e-7, 1e-6),
  cir.chr.h  = 1.5,
  amplify    = TRUE,
  threshold.lty = c(1, 2),
  threshold.col = c("red", "blue"),
  signal.line   = 1,
  signal.cex    = 0.8,
  signal.col    = c("red", "blue"),
  cex           = 0.8,
  chr.den.col   = c("darkgreen", "yellow", "tomato"),
  bin.size      = 1e6,
  outward       = FALSE,
  file          = "png",
  dpi           = 300,
  file.output   = TRUE,
  verbose       = TRUE,
  file.name     = "low_dens",
  width         = 30,
  threshold.lwd = 10
)

#--------------------------------------
# Extract and reshape p‑values: high density
#--------------------------------------
GWASres_high <- GWASDrops_high_dens$GWAResult$high_dens

GWASres_list_high <- lapply(
  unique(GWASres_high$trait),
  function(tr) {
    df <- subset(GWASres_high, trait == tr)
    setNames(
      data.frame(pValue = df$pValue),
      paste0(tr, "_pValue")
    )
  }
)
GWASres.df.high_dens <- do.call(
  cbind,
  c(
    list(
      SNP = GWASres_high$snp[GWASres_high$trait == unique(GWASres_high$trait)[1]],
      chr = GWASres_high$chr[GWASres_high$trait == unique(GWASres_high$trait)[1]],
      pos = GWASres_high$pos[GWASres_high$trait == unique(GWASres_high$trait)[1]]
    ),
    GWASres_list_high
  )
)

# Generate Manhattan‐style plot for high density
CMplot(
  GWASres.df.high_dens[, 1:4],
  plot.type  = "m",
  r          = 0.4,
  col        = c("grey20", "grey70"),
  chr.labels = paste0("Chr", 1:10),
  threshold  = c(1e-7, 1e-6),
  cir.chr.h  = 1.5,
  amplify    = TRUE,
  threshold.lty = c(1, 2),
  threshold.col = c("red", "blue"),
  signal.line   = 1,
  signal.cex    = 0.8,
  signal.col    = c("red", "blue"),
  cex           = 0.8,
  chr.den.col   = c("darkgreen", "yellow", "tomato"),
  bin.size      = 1e6,
  outward       = FALSE,
  file          = "png",
  dpi           = 300,
  file.output   = TRUE,
  verbose       = TRUE,
  file.name     = "high_dens",
  width         = 30,
  threshold.lwd = 10
)
