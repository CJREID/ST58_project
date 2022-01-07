####------SETUP -------####
# Source ABRicate processing script
source("scripts/abricateR.R")
source("scripts/ST58_metadata_curation.R")

# Define sample names - final collection of 752 sequences
sample_names <-
  read_table("raw_data/meta/sample_names.txt", col_names = FALSE)

# Read in metadata
metadata <- read_csv("raw_data/meta/ST58.metadata.csv")

# Create subdirectory to store figures and other outputs generated in these scripts
if (dir.exists("outputs")){
} else {
  dir.create("outputs")
  dir.create("outputs/data")
  dir.create("outputs/figures")
}

####------ASSEMBLY-STATS-------####
# Read in assembly statistics
assembly_stats <- read_tsv("raw_data/assembly-stats/ST58.assemblystats.tsv", skip =1,
                           col_names = c("Name", 
                                         "Assembly_length", 
                                         "n_contigs", 
                                         "Mean_contig_length",
                                         "Longest_contig", 
                                         "Shortest_contig",
                                         "Ns",
                                         "Gaps",
                                         "N50",
                                         "N50n",
                                         "N70",
                                         "N70n",
                                         "N90",
                                         "N90n")) %>% 
  mutate(Name = gsub("\\.fasta", "",Name))

# Check maximum and minimum genome sizes are normal
min(assembly_stats$Assembly_length)  # Min = 4,602,717
max(assembly_stats$Assembly_length)  # Max = 5,890,098 - both normal

# Check contigs < 1000
max(assembly_stats$n_contigs) # Maximum 874

####------POINTFINDER------####
# Read in mutation data
fqr_mut <-
  read_tsv("raw_data/pointfinder/ST58.FQR.mutations.tsv", 
           col_names = c("Mutation",
                         "Nucleotide change",
                         "Amino acid change",
                         "Resistance",
                         "PMID",
                         "Name"))

# Read in phenotypic prediction data
fqr_pheno <-
  read_tsv("raw_data/pointfinder/ST58.FQR.phenotype.prediction.tsv",
           skip = 1,
           col_names = c("Name",
                         "Nalidixic Acid",
                         "Ciprofloxacin",
                         "Unknown")) %>%
  select(-Unknown)

# Join mutations and phenotypes and filter out strains with no predicted resistance
fqr_full <- left_join(fqr_pheno, fqr_mut, by = "Name") %>% filter(`Nalidixic Acid` == 1)

# All strains remaining are predicted to be Nalidixic Acid and Ciprofloxacin resistant so make a new column summarising this as FQR
fqr_full$FQR <- fqr_full$`Nalidixic Acid`

# Select the name and FQR columns
fqr_pos <- unique(fqr_full %>% select(Name, FQR))

# Convert "1" to "Yes"
fqr_pos$FQR <- fqr_pos$FQR %>% recode(`1` = "Yes")

# Get remaining FQR-negative sample names
fqr_neg <- sample_names %>% filter(X1 %notin% fqr_pos$Name) %>% select(Name = X1)

# Give them an FQR column of "No"
fqr_neg$FQR <- rep("No", nrow(fqr_neg))

# Stick them back together for joining to the metadata below
fqr_all <- rbind(fqr_pos, fqr_neg)
####------ARIBA-------####
# Read ARIBA gene screening data
ariba_filt <- read_csv("raw_data/ariba/ST58.ARIBA.csv")

# Filter ARIBA data to the strains that remain in 'metadata'
ariba_filt <- ariba_filt %>% filter(name %in% sample_names$X1)

# Extract ordered names to assign as rownames for heatmap matrices later
strain_names <- ariba_filt$name

# Reorder columns
ariba_filt <-
  ariba_filt %>% dplyr::select(name, O_type, H_type, OH_type, everything())

# Splits the gene hit types (virulence, mobile genetic elements, resistance, plasmid genes) to separate dfs 
vir <- ariba_filt %>% dplyr::select(name, matches("^v_"))
mge <- ariba_filt %>% dplyr::select(starts_with("i_"))
arg <- ariba_filt %>% dplyr::select(starts_with("r_"))
plas <- ariba_filt %>% dplyr::select(starts_with("p_"))

####------VIRULENCE -------####
# Sum the total gene carriage by strain
vir_totals <- rowSums(dplyr::select(vir,-name))

# Order the genes alphabetically
vir <- vir %>% dplyr::select(order(starts_with("v_")),-name)

# Trim gene prefixes
colnames(vir) <-
  gsub(pattern = "^[v]_", "", colnames(vir), perl = TRUE)

# Clean up colnames for messy genes
vir <- vir %>% dplyr::rename(
  papGII_1 = papGII,
  papGII_2 = papGII_CP018957.1,
  papGIII = papGIII_CP010151.1,
  bmaE = bmaE_M15677.1,
  irp2 = irp2__gi,
  sitA = sitA__gi,
  kpsMT_II_1 = kpsMT_II_CFT073,
  kpsMT_II_2 = kpsMT_II_i14,
  kpsMT_III = kpsMT_III_AF007777.1,
  malX = malX_AF081286.1
)

# Extract column names
vir_names <- colnames(vir)

# Create a dataframe 'vir_heat' for use with ggtree/gheatmap
vir_heat <- vir

# Convert the columns to character class
vir_heat[vir_names] <- lapply(vir_heat[vir_names], as.character)

# Replace zeroes with 'Absent' and any non-zero value with 'Present'
vir_heat[vir_heat >= 1] <- "Present"
vir_heat[vir_heat == 0] <- "Absent"

# Assign rownames back to df
rownames(vir_heat) <- strain_names

# Convert all columns to integers
vir <- vir %>% mutate_all(as.integer)

# Sum VAGs for each strain
vir <- vir %>% mutate(`Total VAGs` = rowSums(.))

# Re-assign rownames
vir$Name <- strain_names

# Sum VAGs for the collection
vir_sums <- colSums(vir %>% select(-Name))

####------MOBILE GENETIC ELEMENTS-------####
# Clean the column names and reassign row names for each gene type dataframe
colnames(mge) <- gsub("^i_", "", colnames(mge), perl = TRUE)
rownames(mge) <- strain_names
mge <- lapply(mge, as.numeric)

# Recreate Name column for joining
mge$Name <- strain_names
mge <- as.data.frame(mge)

# Convert to a matrix for heatmap use
mge_heat <- as.matrix(mge %>% select(-Name))

# Replace zeroes with 'Absent' and any non-zero value with 'Present'
mge_heat[mge_heat >= 1] <- "Present"
mge_heat[mge_heat == 0] <- "Absent"

rownames(mge_heat) <- strain_names

####------AMR GENES------####
# Clean up column names
colnames(arg) <- gsub("^[r,i]_", "", colnames(arg), perl = TRUE)
colnames(arg) <- gsub("___", "_", colnames(arg), perl = TRUE)
colnames(arg) <- gsub("__", "_", colnames(arg), perl = TRUE)
colnames(arg) <- gsub("_$", "", colnames(arg), perl = TRUE)

# Re-assign rownames
rownames(arg) <- strain_names

# Remove mdfA, which is not a legitimate horizontally acquired antimicrobial resistance gene
arg <- arg %>% dplyr::select(-mdf_A)

# Convert columns to integers
arg <- arg %>% mutate_all(as.integer)

# Sum up the ARGs per strain
arg <- arg %>% mutate(`Total ARGs` = rowSums(.))

# Re-create Name column for joining
arg$Name <- strain_names
arg_sums <- colSums(arg %>% select(-Name))

# Add predicted FQR for the heatmap
fqr <- fqr_all %>% mutate(FQR = gsub("Yes", "1", FQR), FQR = gsub("No", 0, FQR))
arg <- left_join(arg, fqr)

# Create arg_heat for use with ggtree/gheatmap
arg_heat <- arg %>% select(-`Total ARGs`,-Name)

# Convert the columns to character class
arg_heat <- arg_heat %>% mutate_all(as.character)

# Replace zeroes with 'Absent' and any non-zero value with 'Present'
arg_heat[arg_heat >= 1] <- "Present"
arg_heat[arg_heat == 0] <- "Absent"

# Re-assign rownames
rownames(arg_heat) <- strain_names

####------PLASMID REPLICON GENES------####
# Clean up column names
colnames(plas) <- gsub("^p_", "", colnames(plas), perl = TRUE)
colnames(plas) <- gsub("_$", "", colnames(plas), perl = TRUE)
colnames(plas) <- gsub("__", "_", colnames(plas), perl = TRUE)
colnames(plas) <-
  gsub("_CIT|_I_gamma|_CAV1320", "", colnames(plas), perl = TRUE)
colnames(plas) <- gsub(".G", "_G", colnames(plas), perl = TRUE)

# Re-assign rownames
rownames(plas) <- strain_names
plas$Name <- strain_names

# Create plas_heat Process for use with ggtree/gheatmap
plas_heat <- plas %>% select(-Name)

# Convert the columns to character class
plas_heat <- plas_heat %>% mutate_all(as.character)

# Replace zeroes with 'Absent' and any non-zero value with 'Present'
plas_heat[plas_heat >= 1] <- "Present"
plas_heat[plas_heat == 0] <- "Absent"

# Re-assign rownames
rownames(plas_heat) <- strain_names

####------SEROTYPES------####
# OH for metadata
# Select appropriate columns from ARIBA data - view this object for comprehensive O and H types before it gets filtered to major and minor types
OH <-
  ariba_filt %>% dplyr::select(name, OH_type, O_type, H_type) %>% dplyr::rename(Name = name)

# Define major types as those with >= 10 representatives and minor types as those with < 10 for OH, O and H separately
major_OH <-
  OH %>% dplyr::select(Name, OH_type) %>% group_by(OH_type) %>% filter(n() >= 10)
minor_OH <-
  OH %>% dplyr::select(Name, OH_type) %>% group_by(OH_type) %>% filter(n() < 10)
major_O <-
  OH %>% dplyr::select(Name, O_type) %>% group_by(O_type) %>% filter(n() >= 10)
minor_O <-
  OH %>% dplyr::select(Name, O_type) %>% group_by(O_type) %>% filter(n() < 10)
major_H <-
  OH %>% dplyr::select(Name, H_type) %>% group_by(H_type) %>% filter(n() >= 10)
minor_H <-
  OH %>% dplyr::select(Name, H_type) %>% group_by(H_type) %>% filter(n() < 10)

# Change the value of the types with less than 10 representatives to "Other" for visualisation on tree
minor_OH$OH_type <- "Other"
minor_H$H_type <- "Other"
minor_O$O_type <- "Other"

# Stick them back together row-wise
all_OH_types <-
  rbind(major_OH, minor_OH) %>% dplyr::rename(ID = Name)
all_O_types <- 
  rbind(major_O, minor_O) %>% dplyr::rename(ID = Name)
all_H_types <- 
  rbind(major_H, minor_H) %>% dplyr::rename(ID = Name)

# Join into one table by Name
final_OH <- left_join(all_O_types, all_H_types)
final_OH <- left_join(final_OH, all_OH_types)
final_OH <- final_OH %>% dplyr::rename(Name = ID)

####------F PLASMID pMLST------####
# Read data in
raw <- read.delim("raw_data/ariba/ST58.ariba.F.pMLST.tsv", stringsAsFactors = FALSE)

# Select relevant columns and generate a coverage column by dividing ref_base_assembled by ref_len
raw1 <- raw %>% 
  select(filename, 
         ref_name, 
         ref_len, 
         ref_base_assembled, 
         pc_ident) %>% mutate(coverage = (ref_base_assembled/ref_len*100))

# Filter out non F-type hits
f1 <- raw1 %>% filter(grepl("FII|FIA|FIB",raw1$ref_name))

# Fix reference gene names
f1$ref_name <- gsub("^.*__.*__F", "F", f1$ref_name)
f1$ref_name <- gsub("__.*", "", f1$ref_name)

# Filter based on coverage and identity
f2 <- f1 %>% filter(coverage >=100, pc_ident >=90)

# Cast the data into long format
f2_cast <- reshape2::dcast(f2, filename ~ ref_name, fun.aggregate = max, value.var = "pc_ident")
f2_cast[f2_cast == -Inf] <- 0

# Simplify colnames
colnames(f2_cast) <- gsub(x = colnames(f2_cast), "FIB_", "B")
colnames(f2_cast) <- gsub(x = colnames(f2_cast), "FIA_", "A")
colnames(f2_cast) <- gsub(x = colnames(f2_cast), "FII_", "F")

# Split into separate alleles
FA <- f2_cast %>% select(filename, contains("A"))
FB <- f2_cast %>% select(filename, contains("B"))
FF <- f2_cast %>% select(filename, contains("F"))

# MELT SUB-DFs
# A
FA_melt <- reshape2::melt(FA, id.vars = "filename", variable.name = "A-type")
FA_melt_filt <- FA_melt[FA_melt$value >= 90,]
FA_melt_filt$`A-type` <- as.character(FA_melt_filt$`A-type`) 

# Check that samples don't have more than one allele
if (length(unique(FA_melt_filt$filename)) < ncol(FA_melt_filt)){
  print("Some strains may have more than one allele!")
} else{
  print("All strains contain one hit for this replicon")
}

# B
FB_melt <- reshape2::melt(FB, id.vars = "filename", variable.name = "B-type")
FB_melt_filt <- FB_melt[FB_melt$value >= 90,]
FB_melt_filt$`B-type` <- as.character(FB_melt_filt$`B-type`)

# Check that samples don't have more than one allele
if (length(unique(FB_melt_filt$filename)) < ncol(FB_melt_filt)){
  print("Some strains may have more than one allele!")
} else{
  print("All strains contain one hit for this replicon")
}

# F
FF_melt <- reshape2::melt(FF, id.vars = "filename", variable.name = "F-type")
FF_melt_filt <- FF_melt[FF_melt$value >= 90,]
FF_melt_filt$`F-type` <- as.character(FF_melt_filt$`F-type`) 

# Check that samples don't have more than one allele
if (length(unique(FF_melt_filt$filename)) < ncol(FF_melt_filt)){
  print("Some strains may have more than one allele!")
} else{
  print("All strains contain one hit for this replicon")
}

# Duplicate sorting
F_dup <- FF_melt_filt[duplicated(FF_melt_filt$filename),]
F_duplicates <- FF_melt_filt %>% filter(filename %in% F_dup$filename)
F_non_duplicates <- FF_melt_filt %>% filter(filename %notin% F_dup$filename)

# Get the sample names that have duplicate entries
dup_names <- as.character(unique(F_duplicates$filename))

# This loop splits the F_duplicate df into individual dfs for each sample
datalist <- list()
for (i in seq_along(dup_names)){
  datalist[[i]] <- F_duplicates %>% filter(grepl(dup_names[i], F_duplicates$filename))
}

# This loop then selects the best hit for each sample, but will still result in dfs that have two hits with the same % identity
filtered_datalist <- list()
for (i in seq_along(datalist)){
  filtered_datalist[[i]] <- datalist[[i]] %>% filter(value == max(value))
}

#Initialise two empty dataframes 
split1 <- data.frame()
split2 <- data.frame()

# Row bind the double hits in filtered_datalist to split1 and the singles to split2
for (i in seq_along(filtered_datalist)){
  if (nrow(filtered_datalist[[i]]) > 1){
    split1 <- rbind(split1, filtered_datalist[[i]])
  }
  else{
    split2 <- rbind(split2, filtered_datalist[[i]])
  }
}

# Create a vector that designates duplicate hits as 'a' or 'b' for each sample in split1
split_vect <- rep(c("A", "B"), nrow(split1)/2)

# Add as a column
split1$F_split <- split_vect

# Cast the data into long format
F_dup_solve <- reshape2::dcast(data = split1, filename ~ F_split, value.var = "F-type", drop = FALSE)

# Create a new column with the two hits pasted together
F_dup_solve <- F_dup_solve %>% mutate(`F-type` = paste(A, B, sep = "_")) %>% select(filename, `F-type`)

# Bring back the non-duplicates
F_non_duplicates <- F_non_duplicates %>% select(filename, `F-type`)

# Select appropriate columns from split2
split2 <- split2 %>% select(filename, `F-type`)

# Stick it back together row-wise
FF_processed_dups <- rbind(F_non_duplicates, F_dup_solve, split2)

# Stick EVERYTHING back together
FFA <- left_join(FF_processed_dups, FA_melt_filt %>% select(-value), by = "filename")
FFAB <- left_join(FFA, FB_melt_filt %>% select(-value), by = "filename") %>% select(filename, `F-type`, `A-type`, `B-type`) 
FFAB <- apply(FFAB, 2, as.character)
FFAB <- as.data.frame(FFAB, stringsAsFactors = FALSE)
FFAB <- FFAB %>% replace_na(list(filename = "-", `F-type` = "-", `A-type` = "-", `B-type` = "-"))

# Resolve it 
F_type <- FFAB %>% mutate(`F Plasmid` = paste(`F-type`, `A-type`, `B-type`, sep = ":"))

# Split into types with 10 or more, and less than 10 representatives
major_F <- F_type %>% group_by(`F Plasmid`) %>% filter(n() >= 10)
minor_F <- F_type %>% group_by(`F Plasmid`) %>% filter(n() < 10)

# Designate types with less than 10 as "Other"
minor_F$`F Plasmid` <- "Other"

# Stick them back together
F_types <- rbind(major_F, minor_F)

# Select appropriate columns
all_F_types <- F_types %>% select(Name = filename, `F Plasmid`)

####------I1 PLASMID pMLST------####
I1_pMLST_raw <- read_tsv("raw_data/mlst/I1_pMLST_output.tsv",
         col_names = c("Name", "DB", "ST", "repI1", "ardA",
                       "trbA", "sogS", "pilL"))

I1_pMLST <- I1_pMLST_raw %>% mutate(ST = gsub("-", NA, ST)) %>% select(Name, I1_pMLST = ST)

# Split into types with 10 or more, and less than 10 representatives
majorI1 <- I1_pMLST %>% group_by(I1_pMLST) %>% filter(n() >= 5)
minorI1 <- I1_pMLST %>% group_by(I1_pMLST) %>% filter(n() < 5)

# Designate types with less than 10 as "Other"
minorI1$I1_pMLST <- "Other"

# Stick them back together
I1_simple <- rbind(majorI1, minorI1) %>% rename(I1_simple = I1_pMLST)

####------COLV------------####
# abricateR function takes concatenated Abricate output for any gene screening data (ARGs, VAGs etc) as well as Abricate data for ColV gene screening
# and returns gene co-occurence information as well as applying the Liu criteria to infer which strains most likely carry a ColV plasmid
# The default setting is nucleotide ID 90% and hit length 95%
abricateR(file = "raw_data/abricate/ST58.abricate.contig_cocarriage.tab",
          output = "ST58", 
          ColV_Liu_data = "raw_data/abricate/ST58.abricate.colv.summary.tsv")

# This object has all the ColV gene data 
colV <- ST58.colv.N90.L95

# This one has gene co-occurence data
gene_cooccurence <- `ST58_co-occurence_N90L95`

# Create a new column of yes/no ColV presence based on Liu criteria of at lease one gene from four or more groups
colV <- colV %>% rename(ColV_sum = ColV, Name = name)
# Substitute 0 for 'No'
colV$ColV <- gsub("0", "No", colV$ColV_sum)
# Substitute 1 for 'Yes'
colV$ColV <- gsub("1", "Yes", colV$ColV)
# Convert ColV to integer class
colV$ColV_sum <- as.integer(colV$ColV_sum)

####------pCERC4 MAPPING-------####
# This section involves converting ABRicate hits across the pCERC4 plasmid backbone into a heatmap of nucleotide identity across 100bp bins
# Paths to input files
path_to_tree <- "raw_data/trees/ST58.CGA.ogroot.contree"
path_to_abricate <- "raw_data/abricate/ST58.pCERC4.abricate.tab"
plasrefname <- "pCERC4"

# Minimum hit thresholds
# Minimum hit length (as a percentage [i.e. 0.5 = 0.5%])
min_hit_length <- 0.5
# Minimum nucleotide ID (also as a percentage [i.e. 90 = 90%])
min_hit_id <- 90

# Read in the abricate genotype data sheet
#(small number of rows for colname reassignment)
#This is to reduce memory requirements
abricate_hits <-
  read_delim(
    path_to_abricate,
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    n_max = 10
  )

# Colname reassignment
colnames(abricate_hits)[c(1, 10:11)] <-
  c("name", "perc_coverage", "perc_identity")

# Extract column names for later reassignment
abricate_hits_colnames <- colnames(abricate_hits)

# Read in full abricate data sheet
abricate_hits <-
  read_delim(
    path_to_abricate,
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    col_names = FALSE,
    skip = 1
  )

#Remove cases where there are multiple headers from concatenation of abricate reports
abricate_hits <- abricate_hits %>% filter(X2 != "SEQUENCE")

#Colname reassignment
colnames(abricate_hits) <- abricate_hits_colnames

#Convert percent coverage and identity to numeric type to allow filtering
abricate_hits$perc_coverage <- as.numeric(abricate_hits$perc_coverage)
abricate_hits$perc_identity <- as.numeric(abricate_hits$perc_identity)

#Filter to perc_identity > 90%
abricate_hits <- abricate_hits %>% filter(perc_identity > min_hit_id)
abricate_hits <- abricate_hits %>% filter(perc_coverage > min_hit_length)

#Trim excess characters in the assembly names and reassign this to rownames
abricate_hits$name <- gsub("\\..*", "", abricate_hits$name)

#Read in the tree file
tree <-
  read.tree(file = path_to_tree)

#Subset the hits to strains within the tree (this saves memory)
abricate_hits <- abricate_hits %>% filter(name %in% tree$tip.label)

#Extract from coverage column the coverage range and the length of the reference
abricate_hits$coverage_range <- gsub("\\/.*","", abricate_hits$COVERAGE)
abricate_hits$ref_length <- gsub(".*\\/","", abricate_hits$COVERAGE)

#Save length of plasmid reference as a variable for later use
ref_length <- as.numeric(unique(abricate_hits$ref_length))

#Replace the '-' with a ':' in the coverage
abricate_hits$coverage_range <- gsub("-",":", abricate_hits$coverage_range)

#Create a column for start and end coordinates of hits
abricate_hits$end <- gsub("[0-9]+:","", abricate_hits$coverage_range)
abricate_hits$start <- gsub(":[0-9]+","", abricate_hits$coverage_range)

#Select columns of interest
abricate_hits <- abricate_hits %>%
  select(name, gene = GENE, ref_length, start, end, percentage = perc_coverage)

#Convert start and end coordinates to numeric
abricate_hits$start <- as.numeric(abricate_hits$start)
abricate_hits$end <- as.numeric(abricate_hits$end)

#Create an empty matrix equal to length of ref plasmid
empty_plasrow <- rep(0, times = unique(abricate_hits$ref_length))

#Create an empty matrix with n rows (n = sample size) with ncol == length(ref plasmid)
empty_plasmatrix <- matrix(rep(empty_plasrow,
                               times = length(unique(abricate_hits$name))),
                           nrow = length(unique(abricate_hits$name)))

#Create a list of levels for sample names, a list of start coords and a list of end coords
#and bind these in a list of lists
start_ends <- list(as.list(as.integer(factor(abricate_hits$name, levels = unique(abricate_hits$name)))),
                   as.list(as.integer(abricate_hits$start)),
                   as.list(as.integer(abricate_hits$end)))

name_fact_as_int <- as.integer(factor(abricate_hits$name, levels = unique(abricate_hits$name)))

name_fact <- abricate_hits$name

check_df <- as.data.frame(cbind(name_fact, name_fact_as_int))

check <- unique(check_df$name_fact_as_int)

if (sum(c(1:length(unique(abricate_hits$name))) != check) != 0) {
  print("Error - have a look at the check_df dataframe in the global environment. Names may have been sorted incorrectly.")
  assign("check_df", check_df, envir=globalenv())
  break
}

#Create a counter
counter <- 0

# Map the BLAST hits to our matrix of bp coordinates
for (i in 1:nrow(abricate_hits)){
  sample <- start_ends[[1]][[i]]
  start_coord <- start_ends[[2]][[i]]
  end_coord <- start_ends[[3]][[i]]
  empty_plasmatrix[sample, start_coord:end_coord] <- 1
  counter <- counter + 1
  if(counter %% 1000 == 0){
    message(paste(counter, "out of ", nrow(abricate_hits), "hits processed"))}
}

#Rename matrix
base_matrix <- empty_plasmatrix

#Remove old matrix
rm(empty_plasmatrix)

# Convert matrix to a dataframe
base_matrix <- as.data.frame(base_matrix, stringsAsFactors = FALSE)

# Assign sample names to rows
rownames(base_matrix) <- unique(abricate_hits$name)

# Get the length of the reference sequence
df_length <- length(abricate_hits)

# Generate indices that cover blocks of 100 columns in base_matrix
bin_ranges <- c(seq(from = 1, to = ref_length, by = 100))

bin_ranges2 <- c(seq(from = 100, to = ref_length, by = 100), ref_length)

# Split the ranges into two lists of  [1] start and [2] end indices 
bin_splits <- list(bin_ranges, bin_ranges2)

# Initialise empty vector for loop below
binned_hits <- vector()

#Create a counter
counter <- 0

# Binning loop
for (i in 1:length(bin_ranges2)){
  # If the last bin is only 1 base long then the as.matrix line won't work,
  #so we have to include the if statement below:
  if(i == length(bin_ranges) & bin_ranges[length(bin_ranges)] == bin_ranges2[length(bin_ranges2)]){
    row_sum <- 1
  }else{
    # Generate row sums (i.e. Number of matching bases) for 100 column chunks of the base_matrix
    row_sum <- as.matrix(rowSums(base_matrix[,bin_splits[[1]][i]:bin_splits[[2]][i]]))
  }
  # Bind them together in a new vector
  binned_hits <- cbind(binned_hits, row_sum)
  counter <- counter + 1
  if(counter %% 50 == 0){
    message(paste(counter, " samples out of ", length(bin_splits[[1]]), " processed" ))}
}

#Save rownames
nems <- rownames(binned_hits)

#Convert binned_hits to a data frame
binned_hits <- as.data.frame(binned_hits)

hits_df <- as.data.frame(rowSums(binned_hits))

hits_df$working_names <- rownames(hits_df)

colnames(hits_df) <- c("Percent_hit","working_name")

hits_df$Percent_hit <- round((hits_df$Percent_hit/ref_length) * 100)

# Write the DF
#write.csv(x = binned_hits, file = paste0("outputs/data/plasmid_mapping/",plasrefname, "_plasmid_coverage.csv"), row.names = TRUE)

#write.csv(x = hits_df, file = paste0("outputs/data/plasmid_mapping/",plasrefname, "_plasmid_coverage_percentage.csv"), row.names = TRUE)

####------FIMH----------####
# Read in data and tidy up column names
fimH <- read_csv("raw_data/meta/ST58.fimH.enterobase.csv") %>% 
  dplyr::select(Barcode, `fimH (fimTyper)`) %>% 
  dplyr::rename(fimH = `fimH (fimTyper)`) %>% 
  filter(Barcode %in% metadata$Barcode)

# Define major and minor types by those with more/less than 10 representatives, respectively
major_fimH <-
  fimH %>% dplyr::select(Barcode, fimH) %>% group_by(fimH) %>% filter(n() > 10)
minor_fimH <-
  fimH %>% dplyr::select(Barcode, fimH) %>% group_by(fimH) %>% filter(n() <= 10)
minor_fimH$fimH <- 
  "Other"
fimH_edit <-
  rbind(major_fimH, minor_fimH)

####------dfrA5-IS26 SIGNATURE----------####
# Process the dfrA5-IS26 signature data
abricateR("raw_data/abricate/ST58.IS26_signature.tab", output = "IS26_signature",identity = 99, length = 90)
signature_pos <- `IS26_signature_simple_summary_N99L90` %>% rename(IS26_signature = `848_gb|CP014489.1|:130176-132039`, Name = name)
signature_neg <- sample_names %>% filter(X1 %notin% signature_pos$Name) %>% rename(Name = X1)
signature_neg$IS26_signature <- rep(0, nrow(signature_neg))
IS26_signature <- rbind(signature_pos, signature_neg)

# Convert binary to text
IS26_signature[IS26_signature == 1] <- "Yes"
IS26_signature[IS26_signature == 0] <- "No"
  
####------FASTBAPS------####
# Read in BAPS clustering data
bap_gene <- read_csv("raw_data/fastBAPS/ST58.CGA.fastbaps_clusters.csv") %>% select(Name = Isolates, `Level 1`) %>% arrange(Name)

# Tree
gene_tree <- tree

# Extract names
strain_names <- bap_gene$Name

# Re-assign names to rownames
rownames(bap_gene) <- strain_names

# Generate data frame for plotting panels for BAP clustering of each alignment
plot.df.gene <- data.frame(id = rownames(bap_gene), 
                           fastbaps1 = bap_gene$`Level 1`, 
                           stringsAsFactors = FALSE)

# Create a character column of the values to make plotting easier
plot.df.gene$gene_bap <- as.character(plot.df.gene$fastbaps1)

# Create ggtrtee object with branches coloured by ColV presence/absence
g_gene <- ggtree(gene_tree, branch.length = "none")

# Create facet plots showing fastbaps grouping for each phylogeny
facet_gene <- facet_plot(g_gene, 
                         panel = "fastbaps",
                         data = plot.df.gene,
                         geom = geom_raster, 
                         aes(x = fastbaps1,
                             fill = gene_bap))


# Create value column of 1's for casting purposes
bap_count <- rep(1, nrow(bap_gene))

# Bind it to the data
gene_cast <- cbind(bap_gene, bap_count)

# Cast each level to make a column representing each cluster within a level
bap.df.gene <- gene_cast %>% select(Name, `Level 1`, bap_count) %>% reshape2::dcast(Name ~ `Level 1`)

# Loops to generate column names for each BAP group i.e. BAP1, BAP2 etc 
# Core gene 
bap.gene.cols <- c("Name")
for (f in (seq_along(names(bap.df.gene)))[1:ncol(bap.df.gene)-1]){
  add <- paste("BAP", f, sep = "")
  bap.gene.cols <- c(bap.gene.cols, add)
}

# Assign colnames back to respective data
colnames(bap.df.gene) <- bap.gene.cols

# Recode NAs as 0
bap.df.gene[is.na(bap.df.gene)] <- 0

# Get the sample order from the tree for each
d.gene <- fortify(gene_tree)
dd.gene <- subset(d.gene, isTip)
colorder.gene <- dd.gene$label[order(dd.gene$y, decreasing=TRUE)]

# Re-order dataframes based on their repsective trees
bap.gene.ordered <- bap.df.gene[match(colorder.gene, bap.df.gene$Name),]

# Create table of Names and BAP groups to add to metadata
meta_bap <- reshape2::melt(bap.df.gene) %>% filter(value == 1) %>% select(Name, BAP = variable)

# Write traits file for use with Scoary
write_csv(bap.gene.ordered, "raw_data/scoary/scoary_allBAPS_traits.csv")

####-----COMBINE WITH METADATA------------####
# Join ColV plasmid carriage data
metadata <- left_join(metadata, colV) %>% dplyr::select(colnames(metadata), ColV, ColV_sum)

# Replace NAs with No and 0 for the ColV and ColV_sum columns respectively
metadata$ColV <- replace_na(metadata$ColV, "No")
metadata$ColV_sum <- replace_na(metadata$ColV_sum, 0)

# Join OH data
metadata <- left_join(metadata, final_OH) %>% dplyr::select(colnames(metadata), O_type, H_type, OH_type)

# Join fimH data
metadata <- left_join(metadata, fimH_edit)

# Join F-pMLST data to metadata
metadata <- left_join(metadata, all_F_types)

# Join IS26 signature data
metadata <- left_join(metadata, IS26_signature)

# Join pointfinder data
metadata <- left_join(metadata, fqr_all)

# Join I1 pMLST data
metadata <- left_join(metadata, I1_pMLST) %>% mutate(I1_pMLST = as.factor(I1_pMLST))
metadata <- left_join(metadata, I1_simple)

# Join BAPS data
metadata <- left_join(metadata, meta_bap)

# Read in and join SRA accession numbers for in-house sequences
SRA <- read_tsv("raw_data/meta/SRA_accessions.tsv")
SRA <- SRA %>% select(Name = sample_name, SRA_accession = accession)

# Read in and join BioSample accession numbers for in-house sequences
biosample <- read_tsv("raw_data/meta/biosample_accessions.tsv")
biosample <- biosample %>% select(Name = sample_name, Biosample =  accession)

# Join Biosample and SRA data
our_accessions <- left_join(biosample, SRA)

# Read in accessions from in-houe sequences that are part of other Bioprojects
our_previous_accessions <- read_tsv("raw_data/meta/previous_biosamples.tsv")

# Join them all together including 'public' accessions (object from ST58_metadata_curation.R)
accessions <- rbind(our_accessions, public, our_previous_accessions)

# Join all Biosample and SRA accessions to metadata
metadata <- left_join(metadata, accessions, by = "Name")

# Join assembly stats
metadata <- left_join(metadata, assembly_stats, by = "Name")

# Order columns
metadata <- metadata %>% select(Name, Barcode, Biosample, SRA_accession, 
                                `Data Origin`, PMID, Niche, Source, `Collection Year`, 
                                Continent, Country, BAP, everything(), -ColV_sum)

meta_cols <- names(metadata)

# Create combined metadata and gene screening tables for data exploration
meta_arg <- left_join(metadata, arg, by = "Name") %>% select(-FQR.y)
meta_vir <- left_join(metadata, vir, by = "Name")
meta_plas <- left_join(metadata, plas, by = "Name")
meta_mge <- left_join(metadata, mge, by = "Name")
meta_colv <- left_join(metadata, colV, by = "Name") %>% select(-ColV.y) %>% rename(ColV = ColV.x)


# Combine all of them for Supplementary Table 2
meta_genes <- left_join(meta_arg, meta_vir %>% select(Name, bmaE:ncol(meta_vir)), by = "Name") %>% rename(`fimH allele` = fimH, fimH = fimH.y)
meta_genes <- left_join(meta_genes, meta_plas %>% select(Name, Col156:ncol(meta_plas)), by = "Name")
meta_genes <- left_join(meta_genes, meta_mge %>% select(Name, IS1:ncol(meta_mge)), by = "Name")
meta_genes <- left_join(meta_genes, meta_colv %>% select(Name, cvaA:ncol(meta_colv)), by = "Name") %>%
  select(-contains(".x")) %>% rename_with(~ gsub("\\.y", "", .x)) %>% mutate(across(where(is.numeric), ~replace_na(.x, 0)))

# Write Supplementary Table 1
write_csv(metadata, "outputs/data/SuppData1.full_metadata.csv")

# Write Supplementary Table 2
write_csv(meta_genes, "outputs/data/SuppData2.gene_screening.csv")
