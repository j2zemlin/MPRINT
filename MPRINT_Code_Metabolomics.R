

library(dplyr)
library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(limma) 
library(patchwork)
library(rstatix)
library(Spectra) 
library(MsBackendMgf) 
library(homologueDiscoverer)

###########################################################################################

# Read in data
feature_table <- read_csv("mprint2_quant.csv") 
colnames(feature_table)[3] <- "RT"

metadata <- read.csv("updated_metadata_FBMN_new.csv") 

# Clean filename column but KEEP metadata as a data frame
metadata$filename <- gsub("\\.mzML$", "", metadata$filename)
metadata_clean <- metadata   # full metadata table
#sample_order <- read.csv("metabolomics/sequence.csv") %>% dplyr::select(-1)

annotations <- read.csv("annotations_mprint.csv")
annotations$Scan <- as.character(annotations$Scan)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "Scan")) %>% 
  dplyr::select(1:4,18)

# Data table
data <- feature_table %>% 
  column_to_rownames("row ID") %>% 
  dplyr::select(contains("Peak")) %>% t() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% arrange(SampleID) %>% 
  distinct(SampleID, .keep_all = TRUE)

data$SampleID <- gsub(".mzML Peak area", "", data$SampleID)

#merges

# Make sure all columns are numeric before summing
data <- data %>%
  mutate(across(c(`5081`, `5477`, `5483`, `5637`, `6604`), as.numeric)) %>%
  mutate(`5081_5477_5483_5637_6604` = `5081` + `5477` + `5483` + `5637` + `6604`)


# Merge bile acid features
features_to_merge <- c("5081", "5477", "5430", "5127", "5456", "4715", "4966", 
                       "5017", "5491", "4946", "5483", "5477", "5470", "5051", 
                       "5141", "5454", "5632", "4724", "6158", "4988", "4848", 
                       "5385", "5483")

# Ensure all are numeric and create the merged column
data <- data %>%
  mutate(across(all_of(features_to_merge), as.numeric)) %>%
  rowwise() %>%
  mutate(merged_cluster_A = sum(c_across(all_of(features_to_merge)), na.rm = TRUE)) %>%
  ungroup()

# Clean data of Blank & QC files and filter for files only within metadata 
data_clean <- data %>% 
  dplyr::filter(!str_detect(SampleID, "BLANK|QC")) %>%       # drop BLANK/QC
  dplyr::select_if(~ any(. != 0)) %>%                        # drop all-zero features
  dplyr::filter(SampleID %in% metadata_clean$filename)       # keep only samples in metadata


# RCLR transformation of data
data_final_clr <- decostand(data_clean %>% column_to_rownames("SampleID"), method = "rclr")


# PCA with whole scores
PCA_whole <- mixOmics::pca(data_final_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("SampleID") %>%
  left_join(metadata_clean, by = c("SampleID" = "filename"))

# Define Grouping
i <- "ATTRIBUTE_abtreatment"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Mprint - Whole"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) + coord_fixed()

PCA_plot
# Print the plot to the Viewer (RStudio)
print(PCA_plot)

# PERMANOVA
dist_metabolites <- vegdist(data_final_clr, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, PCA_whole_scores$ATTRIBUTE_Treatment2)
anova(disper_donor)
permanova <- adonis2(dist_metabolites ~ ATTRIBUTE_Treatment2, PCA_whole_scores, na.action = na.omit, by = "terms")
permanova
# Include results from permanova in the unsupervised PCA (R2, F-value, and P-value)

# PLSDA - Cohort
PLSDA_cohort <- mixOmics::plsda(data_final_clr%>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                PCA_whole_scores$ATTRIBUTE_abtreatment, ncomp = 2, scale = TRUE)

# Check the structure of the PLSDA model
str(PLSDA_cohort)

# PLSDA scores and metadata merge
PLSDA_cohort_scores <- data.frame(PLSDA_cohort$variates$X) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metadata_clean, by = c("SampleID" = "filename"))

#PLSDA Plot

custom_colors <- c(
  "Amp" = "#F37D79",
  "Augmentin" = "#ED3426",
  "Mock" = "#2C3685"
)

#custom shapes
custom_shapes <- c(
  "Augmentin" = 15,   # solid square
  "Amp" = 17,  # solid triangle
  "Mock" = 16         # solid circle
)
PLSDA_cohort_plot <- PLSDA_cohort_scores %>%
  ggscatter(x = "comp1", y = "comp2", 
            color = "ATTRIBUTE_abtreatment", 
            shape = "ATTRIBUTE_abtreatment",    
            alpha = 0.7, 
            title = "PLSDA - Cohort",
            xlab = paste("Component 1 (", round(PLSDA_cohort$prop_expl_var$X[1]*100, digits = 1), "%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_cohort$prop_expl_var$X[2]*100, digits = 1), "%)", sep = ""),
            legend.title = "Group", 
            ggtheme = theme_classic()) +
stat_ellipse(aes(color = ATTRIBUTE_abtreatment), type = "t", linetype = 2) +
  geom_point(data = PLSDA_cohort_scores %>%
               group_by(ATTRIBUTE_abtreatment) %>%
               summarise(across(matches("comp"), mean)), 
             aes(comp1, comp2, color = ATTRIBUTE_abtreatment),
             size = 3, shape = 8, show.legend = FALSE) +  
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme(plot.title = element_text(size = 10), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)) + 
  coord_fixed()

# Print the plot to the Viewer (RStudio)
print(PLSDA_cohort_plot)

#Test the efficacy of PLSDA model performance
# Evaluate model performance with repeated M-fold cross-validation
set.seed(123)  # for reproducibility
PLSDA_cohort_perf <- perf(
  PLSDA_cohort,
  validation = "Mfold",
  folds = 4,
  nrepeat = 100,
  progressBar = TRUE,
  auc = TRUE
)
# Plot performance
plot(PLSDA_cohort_perf, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
# Print AUC scores
auc_score <- PLSDA_cohort_perf$auc
print(auc_score)


# Extract loadings from PLSDA model for component 1
Loadings_cohort <- plotLoadings(PLSDA_cohort, comp = 1, method = "mean", contrib = "max", plot = FALSE)

# Extract component 1 loadings as a proper data frame
Loadings_cohort <- as.data.frame(PLSDA_cohort$loadings$X[, 1, drop = FALSE])  # keep as dataframe

# Rename and move rownames to a column
Loadings_cohort <- Loadings_cohort %>%
  rownames_to_column(var = "ID") %>%
  dplyr::rename(Loading = `comp1`)  # or use "V1" if unnamed

# Extract VIPs and filter based on component 1
VIPs_cohort <- as.data.frame(mixOmics::vip(PLSDA_cohort))
VIPs_cohort$ID <- rownames(VIPs_cohort)
VIPs_cohort_filter <- VIPs_cohort %>% dplyr::filter(comp1 > 1)  # keep features with VIP > 1

# Combine VIPs and Loadings
VIPs_cohort_Load <- VIPs_cohort_filter %>% 
  dplyr::select(ID, comp1) %>% 
  left_join(Loadings_cohort, by = "ID") %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% 
  arrange(desc(comp1))

# Write final table to CSV
write_csv(x = VIPs_cohort_Load, file = "VIPs_cohort_Load_comp1.csv")

########################################################################################################

#VIP dotplot and heatmap code
top_feature_ids <- as.character(c(
  3793, 683, 7100, 7258, 1272, 2828, 427, 1426, 5108, 380, 6508, 325, 528, 4175, 2531, 320,
  6179, 521, 3913, 7009, 296, 643, 7117, 7089, 1295, 6718, 4531, 4500, 7506, 7426, 1365, 5111,
  7427, 7428, 1743, 7272, 7227, 327, 5051, 2468, 731, 7094, 555, 218, 1792, 1415, 769, 4118, 7096, 5491
))

# Assume metadata_clean links SampleID to ATTRIBUTE_abtreatment
vip_abundance <- data_clean %>%
  dplyr::select(SampleID, all_of(top_feature_ids)) %>%
  left_join(metadata_clean, by = c("SampleID" = "filename"))

# Get mean peak area for each Feature × Treatment group
group_means <- vip_abundance %>%
  pivot_longer(cols = all_of(top_feature_ids),
               names_to = "Feature", values_to = "peak_area") %>%
  group_by(ATTRIBUTE_abtreatment, Feature) %>%
  summarise(mean_abundance = mean(peak_area, na.rm = TRUE), .groups = "drop")

# Bin mean abundance into Low / Medium / High per treatment
abundance_binned <- group_means %>%
  pivot_wider(names_from = ATTRIBUTE_abtreatment, values_from = mean_abundance) %>%
  mutate(across(-Feature, ~cut(.x,
                               breaks = quantile(.x, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE),
                               labels = c("Low", "Medium", "High"),
                               include.lowest = TRUE)))

# VIP info (from your PLSDA)
vip_info <- vip_df %>%
  filter(ID %in% top_feature_ids) %>%
  mutate(Feature = as.character(ID)) %>%
  dplyr::select(Feature, comp1, Compound_Name)

# Merge abundance bins with VIP info
abundance_categorized_vip <- abundance_binned %>%
  left_join(vip_info, by = "Feature")

# Rank by descending order
abundance_categorized_vip <- abundance_categorized_vip %>%
  arrange(comp1) %>%
  mutate(Feature = factor(Feature, levels = unique(Feature)))

# Heatmap data (long format + numeric encoding of Low / Medium / High)
abundance_tile_data <- abundance_categorized_vip %>%
  pivot_longer(cols = c("Amp", "Augmentin", "Mock"),
               names_to = "Group", values_to = "Level") %>%
  mutate(Numeric_Level = case_when(
    Level == "Low" ~ 0,
    Level == "Medium" ~ 0.5,
    Level == "High" ~ 1,
    TRUE ~ NA_real_
  )) %>%
  # Ensure same Feature ordering as the VIP table
  mutate(Feature = factor(Feature, levels = levels(abundance_categorized_vip$Feature)))

# VIP dotplot (ordered by VIP score) 
vip_dot_plot <- ggplot(abundance_categorized_vip,
                       aes(x = comp1, y = forcats::fct_rev(Feature))) +
  geom_point(color = "#5c7ac4", size = 3) +
  labs(x = "VIP Score (comp1)", y = "Feature ID", title = "Top VIP Features") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(0, 5, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Heatmap plot (aligned ordering)
abundance_tile_plot <- ggplot(abundance_tile_data,
                              aes(x = Group,
                                  y = forcats::fct_rev(Feature),
                                  fill = Numeric_Level)) +
  geom_tile(color = "white", width = 0.95, height = 0.95) +
  scale_fill_gradientn(
    colours = c("#8357A0", "#F5AECA", "#F56B9E"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("Low", "Medium", "High"),
    name = "Level"
  ) +
  coord_fixed(ratio = 1.2) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_blank(),  # y labels come from left panel
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# Combine both heatmap and dotplot 
combined_vip_plot <- vip_dot_plot + abundance_tile_plot +
  patchwork::plot_layout(widths = c(2, 1))

# Show Plot
combined_vip_plot


##################################
#Boxplots for ampicillin and augmetin detection
# Define feature groups to merge
group_1_features <- as.character(c(1401, 1387, 756, 1384, 755))  # Amoxicillin
group_2_features <- as.character(c(3342, 3377, 3002, 5064, 4701, 3897))  # Ampicillin

# Merge intensities and keep SampleID
data_merged <- data_clean %>%
  mutate(
    Merged_Amoxicillin = rowSums(across(all_of(group_1_features)), na.rm = TRUE),
    Merged_Ampicillin = rowSums(across(all_of(group_2_features)), na.rm = TRUE)
  ) %>%
  dplyr::select(SampleID, Merged_Ampicillin, Merged_Amoxicillin)

# Reshape to long and merge metadata
merged_long <- data_merged %>%
  pivot_longer(cols = starts_with("Merged_"), names_to = "feature_id", values_to = "peak_area") %>%
  left_join(metadata_clean, by = c("SampleID" = "filename")) %>%
  filter(!is.na(ATTRIBUTE_abtreatment), ATTRIBUTE_Age %in% c("mom", "infant")) %>%
  mutate(
    ATTRIBUTE_abtreatment = factor(ATTRIBUTE_abtreatment, levels = c("Mock", "Amp", "Augmentin")),
    ATTRIBUTE_Age = factor(ATTRIBUTE_Age, levels = c("mom", "infant"))  # for consistent ordering
  )

# Custom colors by treatment
custom_colors <- c(
  "Amp" = "#F37D79",
  "Augmentin" = "#ED3426",
  "Mock" = "#2C3685"
)

# Plot with faceting by feature and age
merged_plot <- ggboxplot(merged_long,
                         x = "ATTRIBUTE_abtreatment",
                         y = "peak_area",
                         fill = "ATTRIBUTE_abtreatment",
                         add = "jitter",
                         add.params = list(color = "ATTRIBUTE_abtreatment", alpha = 0.7),
                         width = 0.2,
                         lwd = 0.5,
                         ylab = "Peak area",
                         palette = custom_colors) +
  facet_grid(ATTRIBUTE_Age ~ feature_id, scales = "free_y") +
  labs(title = "Summed Intensities of Ampicillin and Amoxicillin Features in Mothers vs Infants",
       x = "Antibiotic Treatment") +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank())

# Display plot
print(merged_plot)

# Boxplots for ampicillin and augmentin detection
# Define feature groups to merge
group_1_features <- as.character(c(1401, 1387, 756, 1384, 755))  # Amoxicillin
group_2_features <- as.character(c(3342, 3377, 3002, 5064, 4701, 3897))  # Ampicillin

# Merge intensities and keep SampleID + log10 transform
data_merged <- data_clean %>%
  mutate(
    Merged_Amoxicillin = rowSums(across(all_of(group_1_features)), na.rm = TRUE),
    Merged_Ampicillin = rowSums(across(all_of(group_2_features)), na.rm = TRUE)
  ) %>%
  mutate(
    log_Merged_Amoxicillin = log10(Merged_Amoxicillin + 1),
    log_Merged_Ampicillin  = log10(Merged_Ampicillin + 1)
  ) %>%
  dplyr::select(SampleID, log_Merged_Ampicillin, log_Merged_Amoxicillin)

# Reshape to long and merge metadata
merged_long <- data_merged %>%
  pivot_longer(cols = starts_with("log_"), names_to = "feature_id", values_to = "peak_area_log") %>%
  left_join(metadata_clean, by = c("SampleID" = "filename")) %>%
  filter(!is.na(ATTRIBUTE_abtreatment), ATTRIBUTE_Age %in% c("mom", "infant")) %>%
  mutate(
    ATTRIBUTE_abtreatment = factor(ATTRIBUTE_abtreatment, levels = c("Mock", "Amp", "Augmentin")),
    ATTRIBUTE_Age = factor(ATTRIBUTE_Age, levels = c("mom", "infant"))
  )

# Colors by treatment
custom_colors <- c(
  "Amp" = "#F37D79",
  "Augmentin" = "#ED3426",
  "Mock" = "#2C3685"
)

# Plot with faceting
merged_plot <- ggboxplot(merged_long,
                         x = "ATTRIBUTE_abtreatment",
                         y = "peak_area_log",
                         fill = "ATTRIBUTE_abtreatment",
                         add = "jitter",
                         add.params = list(
                           shape = 21,      # circle with fill and outline
                           fill = "white",  # white interior of point
                           color = "black", # black outline
                           alpha = 0.9,
                           size = 1.5       # point size
                         ),
                         width = 0.2,
                         lwd = 0.5,
                         ylab = "log10(Peak area + 1)",
                         palette = custom_colors) +
  facet_grid(ATTRIBUTE_Age ~ feature_id, scales = "free_y") +
  labs(title = "Log10 Summed Intensities of Ampicillin and Amoxicillin Features\nin Mothers vs Infants",
       x = "Antibiotic Treatment") +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank())

print(merged_plot)




######################################################################################################



###Clean metabolomics feature table for Bile acid boxplots
# Read feature table and clean column names
feature_table <- fread("mprint2_quant.csv") %>%
  rename_all(~gsub("\\ Peak area", "", .))

# Check cleaned column names
colnames(feature_table)

# Transpose and reformat data using updated column names
data_transpose <- feature_table %>%
  column_to_rownames("row ID") %>%
  dplyr::select(contains(".mzML")) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("filename")

# Extract blank samples
data_blank <- data_transpose %>%
  filter(str_detect(filename, "Blank"))

# Calculate mean, SD, and CV for blank features
blank_feature_info <- data.frame(
  Feature = colnames(data_blank)[-1],
  Mean_blank = data_blank %>% column_to_rownames("filename") %>% colMeans(),
  SD_blank = data_blank %>% column_to_rownames("filename") %>% apply(2, sd)
) %>%
  mutate(CV_blank = SD_blank / Mean_blank) %>%
  filter(Mean_blank > 0) %>%
  arrange(desc(Mean_blank))

# Extract sample data (excluding pools)
data_sample <- data_transpose %>%
  filter(str_detect(filename, "ME")) %>%
  filter(!str_detect(filename, "Pool"))

# Calculate mean, SD, and CV for sample features
sample_feature_info <- data.frame(
  Feature = colnames(data_sample)[-1],
  Mean_sample = data_sample %>% column_to_rownames("filename") %>% colMeans(),
  SD_sample = data_sample %>% column_to_rownames("filename") %>% apply(2, sd)
) %>%
  mutate(CV_sample = SD_sample / Mean_sample) %>%
  filter(Mean_sample > 0) %>%
  arrange(desc(Mean_sample))

# Manually specified features to remove
specified_features <- c("5661", "4075", "3861", "5291", "4484")

# Determine features to remove based on blank/sample ratios and manual list
feature_to_remove <- blank_feature_info %>%
  left_join(sample_feature_info, by = "Feature") %>%
  filter(Mean_blank > 0) %>%
  mutate(Sample_Blank = Mean_sample / Mean_blank) %>%
  filter(Sample_Blank < 5 | is.na(Sample_Blank)) %>%
  bind_rows(blank_feature_info %>% filter(Feature %in% specified_features)) %>%
  distinct(Feature, .keep_all = TRUE)

# Load metadata
metadata <- read_csv("metabolomics_metadata.csv")
data_transpose$filename <- gsub("\\.mzML$", "", data_transpose$filename)

data_transpose <- as_tibble(data_transpose)

# Clean and merge data
data_clean <- data_transpose %>%
  dplyr::select(-dplyr::any_of(feature_to_remove$Feature)) %>%
  dplyr::filter(!str_detect(filename, "6mix|Pool|Blank")) %>%
  dplyr::left_join(
    metadata_clean %>%
      dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, ATTRIBUTE_timepoint, ATTRIBUTE_Day),
    by = "filename"
  )

# Replace negative values with 0
data_clean[data_clean < 0] <- 0

# Select specific features to keep
values_to_keep <- c("7547", "6178", "6062", "5380", "6057", "5484", "7581", "6217", "6604",
                    "6704", "5081", "5477", "5430", "5127", "5456", "4715", "4966", "5017",
                    "5491", "4946", "5483", "5470", "5051", "5141", "5454", "5632",
                    "4724", "6158", "4988", "4848", "5385", "5092", "5202",
                    "5519", "5637", "6813", "6424", "6423", "6419", "7427")

# Load list of bile acid-related features
bile_list <- read_tsv("bile_list_output.tsv")

scan_numbers <- bile_list$X.Scan.
scan_numbers <- as.character(scan_numbers)
data_clean_bile <- data_clean %>%
  dplyr::select(filename, dplyr::any_of(scan_numbers))

#filtering for mom only
merged_df <- data_clean_bile %>%
  left_join(metadata_clean, by = "filename")
filtered_df <- merged_df %>%
  filter(ATTRIBUTE_Age == 'mom')
filtered_df <- filtered_df %>%
  mutate(
    ATTRIBUTE_timepoint = as.factor(ATTRIBUTE_timepoint),
    ATTRIBUTE_abtreatment = as.factor(ATTRIBUTE_abtreatment)
  )
df_long <- filtered_df %>%
  pivot_longer(cols = all_of(scan_numbers),
               names_to = 'ScanNumber',
               values_to = 'Intensity')
df_long <- df_long %>%
  mutate(Log_Intensity = log(Intensity + 1))
df_long <- df_long %>%
  mutate(ScanNumber = as.character(ScanNumber)) %>%
  left_join(
    bile_list %>% mutate(X.Scan. = as.character(X.Scan.)),
    by = c('ScanNumber' = 'X.Scan.')
  )
ba_order <- bile_list %>%
  distinct(X.Scan., BA) %>%
  arrange(BA) %>%
  pull(X.Scan.)
df_long <- df_long %>%
  mutate(ScanNumber = factor(ScanNumber, levels = ba_order))
mean_intensity <- df_long %>%
  group_by(ScanNumber) %>%
  summarise(mean_intensity = mean(Log_Intensity, na.rm = TRUE))
df_long <- df_long %>%
  left_join(mean_intensity, by = 'ScanNumber') %>%
  mutate(Intensity_Deviation = Log_Intensity - mean_intensity)
df_long <- df_long %>%
  complete(ATTRIBUTE_timepoint, ATTRIBUTE_abtreatment, ScanNumber, fill = list(Intensity_Deviation = NA))
df_wide <- df_long %>%
  pivot_wider(names_from = c(ATTRIBUTE_timepoint, ATTRIBUTE_abtreatment),
              values_from = Intensity_Deviation) %>%
  replace_na(list(Intensity_Deviation = 0))


#--Bile acid box plots

tauro_summary <- df_long %>%
  filter(ScanNumber %in% tauro_features) %>%
  group_by(filename) %>%
  summarise(Mean_Tauro_Intensity = mean(Log_Intensity, na.rm = TRUE))


tauro_summary <- tauro_summary %>%
  left_join(metadata_clean, by = "filename")

tauro_summary <- tauro_summary %>%
  filter(ATTRIBUTE_Age == "mom") %>%
  mutate(Time_Cluster = case_when(
    ATTRIBUTE_timepoint %in% c(1, 2) ~ "Early",
    ATTRIBUTE_timepoint %in% c(3, 4, 5) ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  group_by(filename, Time_Cluster) %>%
  summarise(Mean_Tauro_Intensity = mean(Mean_Tauro_Intensity, na.rm = TRUE), .groups = "drop")

#-- Unconjugated bile acid box plots
tauro_features_unconjugated <- c("6424", "6419", "6423", "7427")

custom_colors <- c(
  "Amp" = "#F37D79",
  "Augmentin" = "#ED3426",
  "Mock" = "#2C3685"
)

# Filter for mom samples only and group timepoints into Early and Late
df_long_mom <- df_long %>%
  filter(ATTRIBUTE_Age == "mom") %>%
  mutate(Time_Cluster = case_when(
    ATTRIBUTE_timepoint %in% c(1, 2) ~ "Early",
    ATTRIBUTE_timepoint %in% c(3, 4, 5) ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Time_Cluster))


tauro_summary <- df_long %>%
  filter(ScanNumber %in% tauro_features_unconjugated, ATTRIBUTE_Age == "mom") %>%
  group_by(filename) %>%
  summarise(Mean_Tauro_Intensity = mean(Log_Intensity, na.rm = TRUE)) %>%
  left_join(metadata_clean, by = "filename") %>%
  mutate(Time_Cluster = case_when(
    ATTRIBUTE_timepoint %in% c(1, 2) ~ "Early",
    ATTRIBUTE_timepoint %in% c(3, 4, 5) ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Time_Cluster))


ggplot(tauro_summary, aes(x = interaction(Time_Cluster, ATTRIBUTE_abtreatment), 
                          y = Mean_Tauro_Intensity, 
                          fill = ATTRIBUTE_abtreatment)) +
  stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 1.5) +
  geom_boxplot(
    width = 0.4,
    outlier.shape = NA,
    linewidth = 1.3,
    fatten = 1
  ) +
  geom_jitter(
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8,
    size = 3.5,
    width = 0.2,
    alpha = 0.9
  ) +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Early.Amp", "Late.Amp"),
      c("Early.Augmentin", "Late.Augmentin"),
      c("Early.Mock", "Late.Mock")
    ),
    p.adjust.method = "BH",
    label = "p.format"
  ) +
  labs(
    title = "UnConjugated Bile Acid Intensity by Treatment and Time Cluster",
    y = "Mean Log-Intensity of Tauro-Conjugated Bile Acids",
    x = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none"
  )


#---------AMINE CONJUGATED
amine_conjugated <- c("7547", "6604", "5637", "6813", "7581", "5519", "6217", "6704")

custom_colors <- c(
  "Amp" = "#F37D79",
  "Augmentin" = "#ED3426",
  "Mock" = "#2C3685"
)

# Filter for mom samples only and group timepoints into Early and Late
df_long_mom <- df_long %>%
  filter(ATTRIBUTE_Age == "mom") %>%
  mutate(Time_Cluster = case_when(
    ATTRIBUTE_timepoint %in% c(1, 2) ~ "Early",
    ATTRIBUTE_timepoint %in% c(3, 4, 5) ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Time_Cluster))


tauro_summary <- df_long %>%
  filter(ScanNumber %in% amine_conjugated, ATTRIBUTE_Age == "mom") %>%
  group_by(filename) %>%
  summarise(Mean_Tauro_Intensity = mean(Log_Intensity, na.rm = TRUE)) %>%
  left_join(metadata_clean, by = "filename") %>%
  mutate(Time_Cluster = case_when(
    ATTRIBUTE_timepoint %in% c(1, 2) ~ "Early",
    ATTRIBUTE_timepoint %in% c(3, 4, 5) ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Time_Cluster))


ggplot(tauro_summary, aes(x = interaction(Time_Cluster, ATTRIBUTE_abtreatment), 
                          y = Mean_Tauro_Intensity, 
                          fill = ATTRIBUTE_abtreatment)) +
  stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 1.5) +
  geom_boxplot(
    width = 0.4,
    outlier.shape = NA,
    linewidth = 1.3,
    fatten = 1
  ) +
  geom_jitter(
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8,
    size = 3.5,
    width = 0.2,
    alpha = 0.9
  ) +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Early.Amp", "Late.Amp"),
      c("Early.Augmentin", "Late.Augmentin"),
      c("Early.Mock", "Late.Mock")
    ),
    p.adjust.method = "BH",
    label = "p.format"
  ) +
  labs(
    title = "Amine-Conjugated Bile Acid Intensity by Treatment and Time Cluster",
    y = "Mean Log-Intensity of Tauro-Conjugated Bile Acids",
    x = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none"
  )

#---

#-- Taurine conjugated features
tauro_features_conjugated <- c("5483", "5491", "5477", "5141", "5454", "5470", "5484", "5092",
                               "5051", "4966", "5017", "5081", "4848", "6062", "4715", "6057",
                               "5456", "6178", "5127", "5430", "5385", "5380", "4988", "6158",
                               "5202", "4724", "4946")

custom_colors <- c(
  "Amp" = "#F37D79",
  "Augmentin" = "#ED3426",
  "Mock" = "#2C3685"
)

# Filter for mom samples only and group timepoints into Early and Late
df_long_mom <- df_long %>%
  filter(ATTRIBUTE_Age == "mom") %>%
  mutate(Time_Cluster = case_when(
    ATTRIBUTE_timepoint %in% c(1, 2) ~ "Early",
    ATTRIBUTE_timepoint %in% c(3, 4, 5) ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Time_Cluster))


tauro_summary <- df_long %>%
  filter(ScanNumber %in% tauro_features_conjugated, ATTRIBUTE_Age == "mom") %>%
  group_by(filename) %>%
  summarise(Mean_Tauro_Intensity = mean(Log_Intensity, na.rm = TRUE)) %>%
  left_join(metadata_clean, by = "filename") %>%
  mutate(Time_Cluster = case_when(
    ATTRIBUTE_timepoint %in% c(1, 2) ~ "Early",
    ATTRIBUTE_timepoint %in% c(3, 4, 5) ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Time_Cluster))


ggplot(tauro_summary, aes(x = interaction(Time_Cluster, ATTRIBUTE_abtreatment), 
                          y = Mean_Tauro_Intensity, 
                          fill = ATTRIBUTE_abtreatment)) +
  stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 1.5) +
  geom_boxplot(
    width = 0.4,
    outlier.shape = NA,
    linewidth = 1.3,
    fatten = 1
  ) +
  geom_jitter(
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8,
    size = 3.5,
    width = 0.2,
    alpha = 0.9
  ) +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Early.Amp", "Late.Amp"),
      c("Early.Augmentin", "Late.Augmentin"),
      c("Early.Mock", "Late.Mock")
    ),
    p.adjust.method = "BH",
    label = "p.format"
  ) +
  labs(
    title = "Tauro-Conjugated Bile Acid Intensity by Treatment and Time Cluster",
    y = "Mean Log-Intensity of Tauro-Conjugated Bile Acids",
    x = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none"
  )



