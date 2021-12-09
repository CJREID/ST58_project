####------SETUP-------####
# Source the data processing and metadata processing scripts
source("scripts/enterobase_abricateR.R")
source("scripts/enterobase_metadata_processing.R")

####------READ AND PROCESS DATA-------####
# Define the genome names
genome_list <- "raw_data/enterobase_colv/enterobase_genome_list.txt"

# Run abricateR to process ColV data
EB_abricateR(genomelist = genome_list,
                      abricate_data = "raw_data/enterobase_colv/colV.tab",
                      perc_ident = 90,
                      perc_length = 95,
                      name = "enterobase")

# Assign output of abricateR function
geno <- ColV.enterobase.N90.L95

# Read in serotype and fimH data
entero_OH <- read_tsv("raw_data/enterobase_colv/enterobase_OH.txt") %>% select(Barcode, `O Antigen`, `H Antigen`)
entero_fimH <- read_tsv("raw_data/enterobase_colv/enterobase_fimH.txt") %>% select(Barcode, fimH = `fimH (fimTyper)`)

# Read in metadata sheet. This has been manually curated. 
meta <- read_delim(
        "raw_data/enterobase_colv/enterobase_Manualprocessed_metadata.txt",
        "\t",
        escape_double = FALSE,
        trim_ws = TRUE
)

# Select rows and sources of interest
meta <- meta %>% select(Name, Barcode, Assembly_barcode, Source = Final_Source, Length, Coverage) %>%
  filter(Source == "Bovine" | 
           Source == "Human_other" | 
           Source == "Porcine" | 
           Source == "Poultry" | 
           Source == "Human_ExPEC" | 
           Source == "Companion_animal")

# Read in ST data
STs <- read_delim("raw_data/enterobase_colv/enterobase_STs.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(ST, Barcode)

# Join genotype and metadata based on co-occuring "Assembly_barcodes"
geno_meta <- inner_join(geno, meta)

# Exclude overly large, small or fragmented genomes
geno_meta <- geno_meta %>% mutate(Coverage = as.numeric(Coverage), Length =as.numeric(Length)) %>% 
  filter(Length >= 4500000,Length <= 6500000,Coverage >= 20)


# Join the STs to geno_meta
geno_meta <- left_join(geno_meta, STs, by = "Barcode") 

# Get rid of Enterobase 'negative STs'
geno_meta <- geno_meta %>% filter(ST >= 1)

# Convert to factor class
geno_meta$ST.x <- as.factor(geno_meta$ST)

# Create another ST column by adding "ST" to ST numbers
geno_meta <- geno_meta %>% mutate(ST = paste0("ST", ST.x))

# Calculate ColV carriage rate in full collection
total_pos <- as.integer(geno_meta %>% select(ColV) %>% summarise(sum(ColV)))
total_seqs <- nrow(geno_meta)
ColV_carriage_all = total_pos/total_seqs*100

# Recode Source names
geno_meta$Source <- recode(geno_meta$Source, 
                           Human_ExPEC = "ExPEC", 
                           Human_other = "Human Other",
                           Companion_animal = "Companion Animal")

# Join serotype and fimH data
geno_meta <- left_join(geno_meta, entero_OH, by = "Barcode")
geno_meta <- left_join(geno_meta, entero_fimH, by = "Barcode") %>% 
  select(Name, Assembly_barcode, Barcode, Source, Length, Coverage, ST, ST.x, `O Antigen`, `H Antigen`, fimH, everything())

# Define most abundant STs as those with >100 assemblies
topSTs <- 
  geno_meta %>% group_by(ST) %>% dplyr::summarise(Total_count = n()) %>% filter(Total_count >100)

# Summarise total ColV numbers for every ST by source
ST <- geno_meta %>% group_by(ST, Source, ColV) %>% dplyr::summarise(ColV_source_count = n()) %>% filter(ColV == 1) 

# Summarise ColV numbers by ST
ST2 <- geno_meta %>% group_by(ST, ColV) %>% filter(ColV == 1) %>% summarise(ColV_count = n())

#Join the above two dataframes
intermediate1 <- left_join(ST, ST2)

# Filter out STs with less than 100 genomes and less than 10% ColV Carriage
colV_summary <- left_join(intermediate1, topSTs) %>%
  filter(!is.na(Total_count)) %>% 
  dplyr::mutate(`ColV Source Carriage` = ColV_source_count/Total_count*100,  `ColV ST Carriage` = ColV_count/Total_count*100) %>%
  filter(`ColV ST Carriage` >=10) %>% arrange(desc(`ColV ST Carriage`))

# Extract a list of common STs with >10% ColV carriage
top_colV_STs <- unique(colV_summary$ST)

# Create table 
geno_ST <- left_join(geno_meta %>% filter(ST %in% top_colV_STs), colV_summary %>% select(`ColV ST Carriage`, ST))

# Process Source vs ColV carriage
entero_source <- left_join(geno_meta %>% group_by(Source) %>% dplyr::summarise(`Total Source` = n()),
                     geno_meta %>% group_by(Source, ColV) %>% dplyr::summarise(Total_ColV = n())) %>%
  mutate(`Source Carriage` = `Total_ColV`/`Total Source`*100)

# Write Supplementary Data 4
write_csv(geno_meta, "outputs/data/SuppData4.enterobase_colv.csv")
