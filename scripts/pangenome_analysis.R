####-------ROARY, PIGGY, SCOARY PROCESSING---------####
####-------ROARY---------####
# Read in Roary gene_presence_absence
roary_raw_meta <- read_csv("raw_data/roary/ST58.roary.meta.csv", 
                      col_names = TRUE,
                      col_types = cols(
                        .default = col_character(),
                        `Non-unique Gene name` = col_character(),
                        `No. isolates` = col_double(),
                        `No. sequences` = col_double(),
                        `Avg sequences per isolate` = col_character(),
                        `Genome Fragment` = col_character(),
                        `Order within Fragment` = col_character(),
                        `Accessory Fragment` = col_character(),
                        `Accessory Order with Fragment` = col_character(),
                        `Min group size nuc` = col_double(),
                        `Max group size nuc` = col_double(),
                        `Avg group size nuc` = col_double()
                      ), skip_empty_rows = FALSE, guess_max = Inf)

# Read in binary gene presence absence
roary_tab <- read.table("raw_data/roary/ST58.roary.gene_presence_absence.Rtab", header =TRUE, check.names = FALSE)

# Join metadata and binary table
roary_filter <- left_join(roary_raw_meta, roary_tab)

# Create percentage carriage column to filter on
roary_filter <- roary_filter %>% mutate(Carriage = `No. isolates`/753*100, 
                          Gene_Unique = paste(`Non-unique Gene name`,
                                              Gene, sep = "_")) %>% 
  select(Gene_Unique, 3:4, Carriage, everything())

# Get rid of NAs
roary_filter$Gene_Unique <- gsub("NA_", "", roary_filter$Gene_Unique)

# Define the core: genes present in >= 99% of strains
roary_core <- 
  roary_filter %>% filter(Carriage >= 99)

# Define the conservative accessory: genes present in 15% > strains <= 90%
roary_accessory <-
  roary_filter %>% filter(Carriage < 99) %>% filter(Carriage > 15)

####-------ROARY-SCOARY---------####
# Read in Scoary results
bap_1.scoary <-
  read_csv("raw_data/scoary/roary/BAP1_roary.scoary.csv")

bap_2.scoary <-
  read_csv("raw_data/scoary/roary/BAP2_roary.scoary.csv")

bap_3.scoary <-
  read_csv("raw_data/scoary/roary/BAP3_roary.scoary.csv")

bap_4.scoary <-
  read_csv("raw_data/scoary/roary/BAP4_roary.scoary.csv")

bap_5.scoary <-
  read_csv("raw_data/scoary/roary/BAP5_roary.scoary.csv")

bap_6.scoary <-
  read_csv("raw_data/scoary/roary/BAP6_roary.scoary.csv")

# Combine into list
bap_scoary_list <- list(bap_1 = bap_1.scoary,
                        bap_2 = bap_2.scoary,
                        bap_3 = bap_3.scoary,
                        bap_4 = bap_4.scoary,
                        bap_5 = bap_5.scoary,
                        bap_6 = bap_6.scoary)

# Processing loop for each file
for (f in 1:length(bap_scoary_list)){
  data <- bap_scoary_list[[f]]
  name <- names(bap_scoary_list[f])
  # Filter hypotheticals and combine gene and non-unique gene name to create a unique name that splits paralogs
  data <- data %>% 
    mutate(BAP = rep(name, nrow(data)), 
           Gene_Unique = paste(`Non-unique Gene name`, Gene, sep = "_")) %>% 
    filter(!grepl("hypothetical", Annotation), Benjamini_H_p < 1E-50) %>% 
    select(Gene_Unique, everything())
  
  # Get rid of NAs
  data$Gene_Unique <- gsub("NA_", "", data$Gene_Unique)
  
  # Split into over and under-represented ColV clade genes based on the Odds Ratio
  over_rep <- data %>% filter(Odds_ratio > 1) %>%
    # slice(0:20) %>%
    dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
    dplyr::rename(
      Pos_present = Number_pos_present_in,
      Neg_present = Number_neg_present_in,
      Pos_absent = Number_pos_not_present_in,
      Neg_absent = Number_neg_not_present_in
    ) 
  
  under_rep <- data %>% filter(Odds_ratio < 1) %>%
    # slice(0:20) %>%
    dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
    dplyr::rename(
      Pos_present = Number_pos_present_in,
      Neg_present = Number_neg_present_in,
      Pos_absent = Number_pos_not_present_in,
      Neg_absent = Number_neg_not_present_in
    ) 
  assign(paste(names(bap_scoary_list[f]), "over", sep = "_"), over_rep)
  assign(paste(names(bap_scoary_list[f]), "under", sep = "_"), under_rep)
}

# ONLY BAP2 AND BAP6 HAVE GENES MEETING THE 1E-50 CUTOFF SO JUST COMBINE THESE FOR SUPP TABLE 3

bap_data <- rbind(bap_2_over, bap_2_under, bap_6_over, bap_6_under) %>% select(BAP, everything())

bap_data <- bap_data %>% mutate(BAP = toupper(BAP), BAP = gsub("_", "", BAP))

write_csv(bap_data, "outputs/data/TableS3.roary_scoary.csv")