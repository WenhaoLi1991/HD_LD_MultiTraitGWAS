
#――――――――――――――――――――――――――――――――
# 0. Setup
#――――――――――――――――――――――――――――――――
getwd()
library(statgenGWAS)
library(statgenGxE)
library(statgenSTA)
library(statgenQTLxT)
library(dplyr)

#――――――――――――――――――――――――――――――――
# 1. Load raw data
#――――――――――――――――――――――――――――――――
data_dir           <- file.path(getwd(), "Data")
dropsMarkers_raw_df<- read.delim2(file.path(data_dir, "dropsMarkers_raw.txt"))
IDmatch            <- read.csv(   file.path(data_dir, "ID_match.csv"))
pheno_raw          <- read.csv(   file.path(data_dir, "pheno_raw.csv"))

#――――――――――――――――――――――――――――――――
# 2. Process genotype calls (dropsMarkers)
#――――――――――――――――――――――――――――――――
# 2.1 Strip off metadata columns, keep only marker calls
dropsMarkers_GT <- dropsMarkers_raw_df[ , -c(1:9)]

# 2.2 Convert “0/0”, “0/1”, “1/1” → 0,1,2
dropsMarkers_converted <- dropsMarkers_GT %>%
  mutate(across(everything(), ~ case_when(
    . == "0/0"            ~ 0,
    . %in% c("0/1","1/0") ~ 1,
    . == "1/1"            ~ 2,
    TRUE                  ~ NA_real_
  )))

# 2.3 Transpose and re‑label rows/cols
dropsMarkers <- t(dropsMarkers_converted)
sample.names   <- colnames(dropsMarkers_converted)
rownames(dropsMarkers) <- sub("_.*", "", sample.names)
colnames(dropsMarkers)  <- paste0("X", dropsMarkers_raw_df$ID)

#――――――――――――――――――――――――――――――――
# 3. Synchronize markers & phenotypes
#――――――――――――――――――――――――――――――――
# 3.1 Prepare pheno + ID match, ensure column names align
colnames(pheno_raw)[1] <- "ID"
colnames(IDmatch)[1]  <- "sample"

pheno_raw2 <- pheno_raw %>%
  left_join(IDmatch, by = "ID")  # adds column `sample`

# 3.2 Find common samples, filter & order both tables
common_samples <- intersect(pheno_raw2$sample, rownames(dropsMarkers))

pheno_raw2_filtered <- pheno_raw2 %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

dropsMarkers_filtered <- dropsMarkers[common_samples, , drop = FALSE]

# 3.3 Reassign back to original names
pheno_raw2_ordered <- pheno_raw2_filtered
dropsMarkers       <- dropsMarkers_filtered

#――――――――――――――――――――――――――――――――
# 4. Load & tidy map (dropsMap)
#――――――――――――――――――――――――――――――――
dropsMap <- read.csv(file.path(data_dir, "dropsMap.csv"))
colnames(dropsMap)[1:3] <- c("SNP.names", "chr", "pos")
rownames(dropsMap) <- dropsMap$SNP.names

#――――――――――――――――――――――――――――――――
# 5. Build dropsPheno (for G×E)
#――――――――――――――――――――――――――――――――
# 5.1 Convert SS/TS into days since reference
reference_date <- as.Date("2023-04-26")
pheno_raw2_ordered <- pheno_raw2_ordered %>%
  mutate(
    SS_date = as.Date(SS, format = "%Y/%m/%d"),
    TS_date = as.Date(TS, format = "%Y/%m/%d"),
    DtSS     = as.integer(SS_date - reference_date),
    DtTS     = as.integer(TS_date - reference_date)
  )

# 5.2 Select & reorder the final phenotype columns
selcol <- c(
  "sample","ID","DtSS","DtTS","GP",
  "YPH","YPP","KSR","EL","ED","ARN","AKNR","HKW","TBL",
  "PH","EH","LR","BSR","MCH",
  "DensityIndex","Experiment"
)

dropsPheno <- pheno_raw2_ordered %>%
  select(any_of(selcol)) %>%
  rename(genotype = sample) %>%
  rename(Experiment = last_col()) %>%
  select(
    Experiment, genotype,
    DtSS, DtTS, GP,
    YPH, YPP, KSR, EL, ED, ARN, AKNR, HKW, TBL,
    PH, EH, LR, BSR, MCH
  )

# 5.3 Final object for G×E
dropsPheno_GxE <- dropsPheno


######################################################### all drops ############################################################################################

# List of variables to keep
variables_to_keep <- c("dropsMap", "dropsMarkers", "dropsPheno")

# Remove all other variables
rm(list = setdiff(ls(), variables_to_keep))

# Free up memory
gc()

# Confirm what's left in the environment
print(ls())

colnames(dropsMap)
head(dropsMap)
head(dropsMarkers)[1:6,1:6]
head(dropsPheno)[1:6,1:10]

###############################################################################################################################################################
# sim 3 replicates 
###############################################################################################################################################################
alltraits <- c('DtSS','DtTS','GP',
               'YPH','YPP','KSR_logit','EL','ED','ARN','AKNR','HKW','TBL',
               'PH','EH')


source("dropPheno_simQTL.R")
# dropsPheno <- dropsPheno[alltraits]

# Set random seed for reproducibility
set.seed(2025)

# List of traits to simulate (exclude Experiment and genotype)
trait_cols <- colnames(dropsPheno)[3:ncol(dropsPheno)]

# Create an empty list to store results
dropsPheno_list <- list()

# Loop over each row (i.e. each genotype × experiment row)
for (i in 1:nrow(dropsPheno)) {
  
  # Extract row
  row <- dropsPheno[i, ]
  
  # Simulate replicates for each trait
  for (rep in 1:3) {  # two replicates
    
    sim_row <- row
    
    # For each trait
    for (trait in trait_cols) {
      
      M <- row[[trait]]  # mean value
      
      # Simulate replicate value ~ Normal(mean, small sd)
      # You can adjust the SD_factor (e.g. 0.05 = 5% of mean)
      SD_factor <- 0.05
      simulated_value <- rnorm(1, mean = M, sd = abs(M) * SD_factor + 1e-6)  # avoid sd=0
      
      sim_row[[trait]] <- simulated_value
    }
    
    # Add replicate info
    sim_row$Rep <- rep
    
    # Store the row
    dropsPheno_list[[length(dropsPheno_list) + 1]] <- sim_row
  }
}

# Combine into a dataframe
dropsPheno <- do.call(rbind, dropsPheno_list)




# data transform
data <- dropsPheno$LR
data_corrected <- ifelse(data == 0, 0.001, ifelse(data == 100, 99.999, data))
data_proportions <- data_corrected / 100
logit_trans <- log(data_proportions / (1 - data_proportions + 0.001))
dropsPheno$LR_logit <- logit_trans

data <- dropsPheno$BSR
data_corrected <- ifelse(data == 0, 0.001, ifelse(data == 100, 99.999, data))
data_proportions <- data_corrected / 100
logit_trans <- log(data_proportions / (1 - data_proportions + 0.001))
dropsPheno$BSR_logit <- logit_trans

data <- dropsPheno$KSR
data_corrected <- ifelse(data == 0, 0.001, ifelse(data == 100, 99.999, data))
data_proportions <- data_corrected / 100
logit_trans <- log(data_proportions / (1 - data_proportions + 0.001))
dropsPheno$KSR_logit  <- logit_trans

###################################
data <- dropsPheno$LR
data_corrected <- ifelse(data == 0, 0.001, ifelse(data == 100, 99.999, data))
data_proportions <- data_corrected / 100
logit_trans <- log(data_proportions / (1 - data_proportions + 0.001))
dropsPheno$LR_logit <- logit_trans

data <- dropsPheno$BSR
data_corrected <- ifelse(data == 0, 0.001, ifelse(data == 100, 99.999, data))
data_proportions <- data_corrected / 100
logit_trans <- log(data_proportions / (1 - data_proportions + 0.001))
dropsPheno$BSR_logit <- logit_trans

data <- dropsPheno$KSR
data_corrected <- ifelse(data == 0, 0.001, ifelse(data == 100, 99.999, data))
data_proportions <- data_corrected / 100
logit_trans <- log(data_proportions / (1 - data_proportions + 0.001))
dropsPheno$KSR_logit  <- logit_trans

# Check
str(dropsPheno)
head(dropsPheno)

# dropsPheno <- dropsPheno
###############################################################################################################################################################
###############################################################################################################################################################
# Step 1: Rename levels in "Experiment" column
# dropsPheno$Experiment <- ifelse(
#   dropsPheno$Experiment == "低", "low_dens",
#   ifelse(dropsPheno$Experiment == "高", "high_dens", dropsPheno$Experiment))

# Step 2: Define the columns for which you want to calculate the mean
# selected_trait <- c("生育期.天","倒伏率",'倒伏率_logit', "空杆率", '空杆率_logit',"出籽率", '出籽率_logit',"公顷产量.kg", "单株产量.kg", 
#                         "穗长.cm", "穗粗.cm", "穗行数均值", "行粒数均值", 
#                         "百粒重.g", "秃尖.cm", "株高.cm", "穗位高.cm", "收获时含水量")

# selected_trait <- c("GP","LR","LR_logit","BSR","BSR_logit","YPH","YPP","KSR","KSR_logit","EL","ED","ARN", "AKNR","HKW","TBL","PH","EH","MCH")




# selected_trait <- c('DtSS','DtTS','GP',
#                         'YPH','YPP',"KSR",'KSR_logit','EL','ED','ARN','AKNR','HKW','TBL',
#                         'PH','EH',
#                         "LR",'LR_logit',"BSR",'BSR_logit',
#                         'MCH')
data_dir <- file.path(getwd(), "Data")
load(file=file.path(data_dir,'alldropdata.Rdata'))


selected_trait <-  c('DtSS','DtTS','GP','YPP','EL','ED','ARN','AKNR','HKW','TBL','EH')


# str(dropsPheno)

# Step 3: Group by genotype and calculate the mean of low_dens and high_dens
library(dplyr)
mean_rows <- dropsPheno %>%
  filter(Experiment %in% c("low_dens", "high_dens")) %>%  # Keep only low_dens and high_dens
  group_by(genotype) %>%  # Group by genotype
  summarise(across(all_of(selected_trait), mean, na.rm = TRUE)) %>%  # Calculate the mean for each column
  mutate(Experiment = "mean")  # Add "mean" to the Experiment column

# Step 4: Add the mean rows to the original data frame
dropsPheno_with_means <- bind_rows(dropsPheno, mean_rows)

# Step 5: View the final data frame
print(head(dropsPheno_with_means))

order(dropsPheno_with_means$genotype)

###############################################################################################################################################################
###############################################################################################################################################################
# sel.trait <- "穗长.cm"


options("statgen.trialColors" = c("#FF9700", "#00C9D1"))

GxEres.ls <- list()

# dropsPheno_real <- dropsPheno
# dropsPheno <- dropsPheno

for (i in 1:length(selected_trait)) {
  sel.trait <- selected_trait[i]
  ###############################################################################################################################################################
  ###############################################################################################################################################################
  dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "genotype", trial = "Experiment")
  
  # dropsVarComp <- gxeVarComp(TD = dropsTD, trait = sel.trait, trials = names(dropsTD))
  # dropsVarComp <- gxeVarComp(TD = dropsTD, trait = sel.trait, nestingFactor = "Experiment")
  dropsVarComp <- gxeVarComp(TD = dropsTD, trait = sel.trait)
  
  
  summary(dropsVarComp)
  
  diagnostics(dropsVarComp)
  
  vc(dropsVarComp)
  
  herit(dropsVarComp)
  
  predGeno <- predict(dropsVarComp)
  
  # head(predGeno)
  
  predGenoTrial <- predict(dropsVarComp, predictLevel = "trial")
  
  # head(predGenoTrial)
  
  GxEres.ls[[sel.trait]][["dropsTD"]] <- dropsTD
  GxEres.ls[[sel.trait]][["dropsVarComp"]] <- dropsVarComp
  
  ###############################################################################################################################################################
  ###############################################################################################################################################################
  
  dropsTD_addMean <- statgenSTA::createTD(data = dropsPheno_with_means, genotype = "genotype", trial = "Experiment")
  
  dropsFW <- gxeFw(TD = dropsTD_addMean, trait = sel.trait)

  summary(dropsFW)
  
  dropsGGE <- gxeGGE(TD = dropsTD_addMean, trait = sel.trait)
  
  summary(dropsGGE) 
  
  GxEres.ls[[sel.trait]][["dropsTD"]] <- dropsTD
  GxEres.ls[[sel.trait]][["dropsVarComp"]] <- dropsVarComp
  GxEres.ls[[sel.trait]][["dropsTD_addMean"]] <- dropsTD_addMean
  GxEres.ls[[sel.trait]][["dropsFW"]] <- dropsFW
  GxEres.ls[[sel.trait]][["dropsGGE"]] <- dropsGGE

}

###############################################################################################################################################################
###############################################################################################################################################################

library(gridExtra)
library(ggplot2)
library(grid)

# Create an empty list to store all the plots
all_plots <- list()
all_GGEbiplot <- list()
# Loop over each trait and generate 5 plots for each trait
for (i in 1:length(selected_trait)) {
  
  sel.trait <- selected_trait[i]
  
  # Extract necessary data for each trait
  dropsTD <- GxEres.ls[[sel.trait]][["dropsTD"]]
  dropsVarComp <- GxEres.ls[[sel.trait]][["dropsVarComp"]]
  dropsTD_addMean <- GxEres.ls[[sel.trait]][["dropsTD_addMean"]]
  dropsFW <- GxEres.ls[[sel.trait]][["dropsFW"]]
  dropsGGE <- GxEres.ls[[sel.trait]][["dropsGGE"]]
  
  # Convert base R or incompatible plots to grobs
  plot1 <- grid.grabExpr(plot(dropsVarComp))  # Wrap base R plot as grob
  
  plot2 <- grid.grabExpr(plot(dropsTD, plotType = "box", traits = sel.trait, colorTrialBy = "trial", orderBy = "descending"))
  plot3 <- grid.grabExpr(plot(dropsTD, plotType = "scatter", traits = sel.trait, colorTrialBy = "trial"))
  
  # Use ggplot-compatible plots directly
  plot4 <- plot(dropsGGE, plotType = "GGE", scale = 0.5) +  
    # ggtitle(paste0("GGE for ", sel.trait, ' (h^2=', round(herit(dropsVarComp), 2), ')')) +
    ggtitle(sel.trait) +
    guides(shape = guide_legend(nrow = 3, byrow = FALSE)) +
    ggplot2::theme(
      legend.position = "top",
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA,
                                           color = "black",
                                           size = 0.5,
                                           linetype = "solid"))
  
  plot5 <- plot(dropsFW, plotType = "line") +  
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))
  
  # Add the 5 plots for this trait into the all_plots list
  all_plots <- c(all_plots, list(plot1, plot2, plot3, plot4, plot5))
  all_GGEbiplot <- c(all_GGEbiplot,list(plot4))
  # all_GGEbiplot <- c(all_GGEbiplot,list(plot3))
}

png(file.path(getwd(),'Results','GxE_plots.png'), width = 4000, height = 6000, res = 300)  # Adjust dimensions and resolution

grid.arrange(
  grobs = all_plots,  # All collected plots
  ncol = 5  # Number of columns in the panel
)
dev.off()  # Close the device



png(file.path(getwd(),'Results','GGE_plots.png'), width = 2000, height = 2000, res = 250)  # Adjust dimensions and resolution

grid.arrange(
  grobs = all_GGEbiplot,  # All collected plots
  ncol = 4  # Number of columns in the panel
)
dev.off()  # Close the device


# names(all_GGEbiplot)

for (trait in selected_trait) {
  
  
  print(trait)
  
  
  png(file.path(getwd(),'Results',paste0("GGEbiplot_",trait,'.png')), width = 700, height = 700, res = 300)  # Adjust dimensions and resolution
  
  print(all_GGEbiplot[[which(selected_trait == trait)]])
  # grid.arrange(
  #   grobs = all_GGEbiplot,  # All collected plots
  #   ncol = 4  # Number of columns in the panel
  # )
  dev.off()  # Close the device
  
  
  
}


# make plot one by plot type
#vc plot###############################################################################################################

vc.ls <- list()
for (i in 1:length(selected_trait)) {

sel.trait <- selected_trait[i]

dropsVarComp <- GxEres.ls[[sel.trait]][["dropsVarComp"]]
vcovPerc <- dropsVarComp$fullRandVC$vcovPerc

vc.ls[[sel.trait]] <- data.frame(residuals = vcovPerc[4],
                                 genotype_by_trial = vcovPerc[3],
                                genotype = vcovPerc[2],
                                trial = vcovPerc[1])

}

vc.df <- do.call('rbind', vc.ls)
vc.df$trait <- rownames(vc.df)

# Load ggplot2
library(ggplot2)
library(reshape2)

# Reshape the data from wide to long format for ggplot
vc_long <- melt(vc.df, id.vars = "trait", 
                variable.name = "component", 
                value.name = "proportion")

vc_long$trait <- factor(vc_long$trait, levels = selected_trait)
vc_long$component <- factor(vc_long$component, levels = c(
                                                          "genotype",
                                                          "trial",
                                                          "genotype_by_trial",
                                                          "residuals" ))
vc_long$proportion <- round(vc_long$proportion*100,2)
vcplot <- ggplot(vc_long, aes(x = trait, y = proportion, fill = component)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  # scale_fill_manual(values = c("residuals" = "#CC79A8", 
  #                              "genotype" = "#0EA079",
  #                              "trial" = "#5AB4E5")) +

  scale_fill_manual(values = c("residuals" = "lightsteelblue", 
                               "genotype_by_trial" = "#00C9D1",
                               "genotype" = "#FF9700",
                               "trial" = "#FF4830"),) +


  geom_text(aes(label = round(proportion, 2)),  # Add labels with rounded proportions
            position = position_stack(vjust = 0.5),  # Position labels in the center of each segment
            size = 3,  # Text size
            color = "black") +  # Label text color
  # theme_minimal() +
  labs(title = "",
       x = "Trait",
       y = "Explained Proportion",
       fill = "Variance components:") +
  # theme(
  #   axis.text.x = element_text(angle = 45, hjust = 1),
  #   plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  #   # legend.position = 'top'
  # )+
  theme(
    legend.position = 'top',
    axis.text.x = element_text(angle = 45, hjust = 1),
    # legend.title = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank()
    # strip.text.x = element_blank()
    # panel.border = element_rect(fill = NA, color = "black",size = 0.5, linetype = "solid"
    # )
  )

print(vcplot)



png(file.path(getwd(),'Results','vc_sim_plot.png'), width = 2000, height = 2000, res = 300)  # Adjust dimensions and resolution
print(vcplot)
dev.off()  # Close the device

#correlation plot#################################################################################################################
# Load required libraries
library(GGally)
library(ggplot2)
library(dplyr)

# Ensure `Experiment` is a factor
dropsPheno$Experiment <- as.factor(dropsPheno$Experiment)

# Select only numeric columns for ggpairs
traits <- selected_trait

# Verify numeric data types
dropsPheno_numeric <- dropsPheno %>%
  select(all_of(selected_trait)) %>%
  mutate_if(is.character, as.numeric)

# Combine with `Experiment` for coloring
dropsPheno_filtered <- dropsPheno %>%
  select(Experiment, all_of(selected_trait))

library(GGally)
library(ggplot2)

# Define custom colors for the groups
# custom_colors <- c("high_dens" = "#B42B22",  # Red for high_dens
#                    "low_dens" = "#315A89")  # Blue for low_dens
dropsPheno_filtered$Experiment <- factor(dropsPheno_filtered$Experiment, levels = c("low_dens", 'high_dens'))

Experimentcolor <- ifelse(dropsPheno_filtered$Experiment  == "high_dens","orange2","cyan3")


pairplot <- ggpairs(dropsPheno_filtered, columns=2:ncol(dropsPheno_filtered), aes(colour=Experimentcolor),
        lower=list(continuous="smooth"), 
        diag=list(continuous="densityDiag")
        # upper = list(continuous = "blank")
        )+
  scale_color_manual(values = unique(Experimentcolor))+
  scale_fill_manual(values = unique(Experimentcolor))+


  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank()
  )

# Print the pair plot
# print(pairplot)
# scale_fill_manual(values = c( "Microbe_1"="#D3D1D2", "Microbe_2"="#CC79A8", "Microbe_3"="#7AC0EA", "Microbe_4"="#0EA079","Microbe_5"="#68C49D",
#                               "Microbe_6"="#EACEDF", "Microbe_7"="#5AB4E5","Microbe_8"="#D99BBB","Microbe_9"="#98D1F0", "Microbe_10"="#7F7F7F")) +  

png(file.path(getwd(),'Results','pairplot_new.png'), width = 4000, height = 6000, res = 300)  # Adjust dimensions and resolution
print(pairplot)
dev.off()  # Close the device




#boxplot#########################################################################################################################
library(gghalves)

plot_list <- list()
for (i in 1:length(selected_trait)) {
trait <- selected_trait[i]
pheno_long <- dropsPheno %>% select('Experiment', trait)
colnames(pheno_long)[2] <- 'value'
custom_colors <- c('high_dens' = "#FF9700", 'low_dens' = "#00C9D1")


pheno_long$Experiment <- factor(pheno_long$Experiment, levels = c("low_dens", 'high_dens'))

p <- ggplot(pheno_long, aes(x = Experiment, y = value)) +
  # Add half violin plots (left for H, right for L)
  geom_half_violin(data = pheno_long %>% filter(Experiment == 'low_dens'), aes(fill = Experiment), alpha = 0.8, width = 0.7, side = "l") +  # Right side for L group
  geom_half_violin(data = pheno_long %>% filter(Experiment == 'high_dens'), aes(fill = Experiment), alpha = 0.8, width = 0.7, side = "r") +  # Left side for H group
 
  # Add boxplot
  geom_boxplot(aes(fill = Experiment), width = 0.2) +  # Boxplot with no outliers
  
  # Add points
  # geom_point(aes(color = Experiment), size = 2, alpha = 0.6) +  # Points for each data point
  
  # Add lines connecting paired data points
  geom_line( linewidth = 0.2) +  # Connecting lines based on Difference
  
  # Manually set colors for violin and points using scale_fill_manual
  scale_fill_manual(values = custom_colors) +  # Custom fill colors for H and L groups
  scale_color_manual(values = unique(custom_colors))+
  # Manually set colors for lines indicating upward/downward trends
  # scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "blue")) +  # Red for upward, Blue for downward
  
  # Labels and title
  labs(title = paste(trait), x = "", y = paste( "Value")) +
  
  # Theme settings
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'none',
    # legend.title = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    # axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "black",size = 0.5, linetype = "solid"
    # )
  ))

# Store the plot in the list
plot_list[[trait]] <- p
}

# Combine all the plots using cowplot::plot_grid
combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)

# Print the combined plot
print(combined_plot)

ggsave(filename = file.path(getwd(),'Results','density_boxplot.png'), plot = combined_plot, width = 8, height =8, units = "in")
