####-----PROCESSING-----####
# Read in raw Enterobase metadata sheet
raw_data <-
  read_delim(
    "raw_data/meta/ST58.raw_enterobase_metadata.txt",
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

# Read in spreadsheet with Enterobase barcodes and working names
samples <-
  read_csv(
    "raw_data/meta/names_barcodes.csv",
    col_types = "c",
    col_names = c('Working_name', 'Barcode')
  )

#Define column names and merge barcodes into main spreadsheet
meta <- merge(raw_data, samples, by = "Barcode")

#Remove whitespace from colnames
colnames(meta) <- gsub(" ", "_", colnames(meta))
#Set row names to working names
rows <- meta$Working_name

#Subset to only include sources and names
sources <- meta[, c(1, 4, 49, 3, 17:19)]
ours <- filter(sources, grepl('ithree', Lab_Contact))
NDs <- sources %>% filter(Source_Details == "ND")
sources <-
  filter(sources,
         !grepl('ithree', Lab_Contact) & sources$Source_Details != "ND")

# Extract SRA and Biosamples for public genomes
public <- meta %>% filter(!grepl('ithree', Lab_Contact)) %>% select(Name = Working_name, Biosample = Sample_ID, SRA_accession = Working_name)

#Split source info into separate dataframes/tibbles or whatever they are
niches <- unique(sources$`Source Niche`)
types <- unique(sources$`Source Type`)
details <- unique(sources$`Source Details`)

##This section involves splitting off dataframes with similar metadata that needs to be edited.
ours <-
  ours %>% mutate(New_Niche = paste(ours$Source_Niche, "__", sep = ""))
ours <-
  ours %>% mutate(New_Type = paste(ours$Source_Type, "__", sep = ""))
our_livestock <- ours %>% filter(Source_Niche == "Livestock")
our_livestock$New_Type <-
  gsub("Avian", "Poultry", our_livestock$New_Type)
our_other <-
  ours %>% filter(
    Source_Niche != "Food" &
      Source_Niche != "Environment" &
      Source_Niche != "Human" & Source_Niche != "Livestock"
  )
our_food <- ours %>% filter(Source_Niche == "Food")
our_food$New_Type <- gsub("Avian", "Poultry", our_food$New_Type)
our_water <- ours %>% filter(Source_Niche == "Environment")
our_human <- ours %>% filter(Source_Niche == "Human")

# Add New Niche and Type columns by adding underscores to the existing Niche and Type

# Changing food samples from livestock to fall under their ;ivestock categories whilst keeping Plant Foods under Foods Niche
our_food$New_Niche <-
  ifelse(!grepl("Plant", our_food$Source_Type),
         "Livestock__",
         our_food$New_Niche)

# Spltting water into 'Wastewater' and 'Water'
our_water$New_Type <-
  ifelse(
    grepl(".*[Ww]astewater.*", our_water$Source_Details),
    "Wastewater__",
    our_water$New_Type
  )
# Humans
our_human$New_Type <- "ExPEC__"

# Edit ND samples
NDs <-
  NDs %>% mutate(New_Niche = paste(Source_Niche, "__", sep = ""))
NDs <- 
  NDs %>% mutate(New_Type = paste(Source_Type, "__", sep = ""))
NDs$New_Type <-
  ifelse(grepl("Human", NDs$Source_Type), "Human Other__", NDs$New_Type)

# Remove anything with Animal Feed as the Source_Type
sources <- sources %>% filter(Source_Type != "Animal Feed")

# This section is all the regex substitutions used to curate the metadata before manual editing

# Change anything with the Bos taurus, beef, cattle, cow, bovine or bovinae into Bovine
sources$New_Niche <-
  gsub(
    ".*[Bb]os [Tt]aurus.*|.*[Bb]os [Pp]rimigenius.*|.*[Bb]eef.*|.*[Vv]eal.*|.*[Cc]attle.*|.*[Cc]ow.*|.*[Bb]ovine.*|[Bb]ovinae.*",
    "Livestock__",
    sources$Source_Details
  )

sources$New_Type <-
  gsub(
    ".*[Bb]os [Tt]aurus.*|.*[Bb]os [Pp]rimigenius.*|.*[Bb]eef.*|.*[Vv]eal.*|.*[Cc]attle.*|.*[Cc]ow.*|.*[Bb]ovine.*|[Bb]ovinae.*",
    "Bovine__",
    sources$Source_Details
  )

# Remove semi-colons from new columns
sources$New_Niche <- gsub(";", "", sources$New_Niche)
sources$New_Type <- gsub(";", "", sources$New_Type)

# Change anything with Porcine, pork, pig, sus scrofa or sus domesti into Porcine
sources$New_Niche <-
  gsub(
    ".*[Ss]wine.*|.*[Pp]orcine.*|.*PORK.*|.*pork.*|.*Pork.*|.*[Pp]ig.*|.*[Ss]us [Ss]crofa.*|.*[Ss]us [Dd]omesti.*",
    "Livestock__",
    sources$New_Niche
  )

sources$New_Type <-
  gsub(
    ".*[Ss]wine.*|Swine|.*[Pp]orcine.*|.*PORK.*|.*pork.*|.*Pork.*|.*[Pp]ig.*|.*[Ss]us [Ss]crofa.*|.*[Ss]us [Dd]omesti.*",
    "Porcine__",
    sources$New_Type
  )

# Change anything to do with poultry, gallus, egg, turkey, chicken or duck into poultry
sources$New_Niche <-
  gsub(
    ".*[Bb]roiler.*|.*[Pp]oultry.*|.*[Gg]allus.*|.*[Ee]gg.*|.*[Tt]urkey.*|.*[Cc]hicken.*|.*[Dd]uck.*",
    "Livestock__",
    sources$New_Niche
  )
sources$New_Type <-
  gsub(
    ".*[Bb]roiler.*|.*[Pp]oultry.*|.*[Gg]allus.*|.*[Ee]gg.*|.*[Tt]urkey.*|.*[Cc]hicken.*|.*[Dd]uck.*",
    "Poultry__",
    sources$New_Type
  )

# Change Canis/Canine/Dog/Canius (some people have incorrect spelling) into Companion animal
sources$New_Niche <-
  gsub(
    ".*[Cc]anis.*|.*[Cc]anis].*|.*[Cc]anine.*|^[^Prairie.*[Dd]og.*|.*[Cc]anius.*|dog",
    "Companion Animal__",
    sources$New_Niche
  )

sources$New_Type <-
  gsub(
    ".*[Cc]anis.*|.*[Cc]anis].*|.*[Cc]anine.*|^[^Prairie.*[Dd]og.*|.*[Dd]og.*|.*[Cc]anius.*|dog",
    "Canine__",
    sources$New_Type
  )

# As above but for feline, cat (avoiding 'catheter') and felis
sources$New_Niche <-
  gsub(
    ".*[Ff]eline].*|.*[Cc]at[^heter].*|.*[Ff]elis.*",
    "Companion Animal__",
    sources$New_Niche
  )

sources$New_Type <-
  gsub(".*[Ff]eline].*|.*[Cc]at[^heter].*|.*[Ff]elis.*",
       "Feline__",
       sources$New_Type)

# As above but for horses
sources$New_Niche <-
  gsub(".*[Ee]quine.*|.*[Hh]orse.*",
       "Companion Animal__",
       sources$New_Niche)

sources$New_Type <-
  gsub(".*[Ee]quine.*|.*[Hh]orse.*", "Equine__", sources$New_Type)

# Find all urinary/urine associated strains, ExPEC and UPEC and label them as Human_ExPEC. Avoid words like "during" with the '[ary|e]'
sources$New_Niche <-
  gsub(".*UTI.*|.*[Uu]rin[ary|e].*|.*ExPEC.*|.*UPEC.*",
       "Human__",
       sources$New_Niche)

sources$New_Type <-
  gsub(".*UTI.*|.*[Uu]rin[ary|e].*|.*ExPEC.*|.*UPEC.*",
       "ExPEC__",
       sources$New_Type)

# Avoiding "bloody diarrheoa", find all with blood, or the previous format and assign to Human_ExPEC
sources$New_Niche <-
  gsub(".*[Bb]lood.*",
       "Human__",
       sources$New_Niche)
sources$New_Type <-
  gsub(".*[Bb]lood.*",
       "ExPEC__",
       sources$New_Type)

# Assign bloody diarrheoa designations as Human_other
sources$New_Niche <-
  gsub(".*[Bb]loody diar.*", "Human__", sources$New_Niche)

sources$New_Type <-
  gsub(".*[Bb]loody diar.*", "IPEC__", sources$New_Type)

# A common pattern from uploaders is a certain format beginning with Homo sapiens and ending with respiratory or wound. Assign these as Human ExPEC
sources$New_Niche <-
  gsub(
    ".*Homo sapiens.*[Ff]aeces|.*Homo sapiens.*[Ff]aecal|.*Homo sapiens.*[Ff]eces|.*Homo sapiens.*[Ff]ecal.*",
    "Human__",
    sources$New_Niche
  )

sources$New_Type <-
  gsub(
    ".*Homo sapiens.*[Ff]aeces|.*Homo sapiens.*[Ff]aecal|.*Homo sapiens.*[Ff]eces|.*Homo sapiens.*[Ff]ecal.*",
    "Human Other__",
    sources$New_Type
  )

# Assign human strains not ending in sputum or fluid to Other
sources$New_Niche <-
  gsub(
    ".*Homo sapiens.*[^putum]$|.*Homo sapiens.*[^fluid]$",
    "Human__",
    sources$New_Niche
  )

sources$New_Type <-
  gsub(".*Homo sapiens.*[^putum]$|.*Homo sapiens.*[^fluid]$",
       "Human Other__",
       sources$New_Type)

# Remainder ending in sapiens; convert to Other
sources$New_Niche <-
  gsub(".*Homo sapiens$",
       "Human__",
       sources$New_Niche)

sources$New_Type <-
  gsub(".*Homo sapiens$",
       "Human Other__",
       sources$New_Type)

# Assign Wastewater to Niche Environment and Type Wastewater
sources$New_Niche <-
  gsub("^Wastewater.*|.*[Ww]astewater",
       "Environment__",
       sources$New_Niche)

sources$New_Type <-
  gsub("^Wastewater.*|.*[Ww]astewater",
       "Wastewater__",
       sources$New_Type)

# Sheep to Livestock/Ovine
sources$New_Niche <-
  gsub(".*aries.*|.*sheep.*", "Livestock__", sources$New_Niche)

sources$New_Type <-
  gsub(".*aries.*|.*sheep.*", "Ovine__", sources$New_Type)

# Steller's Eider is a type of duck so will be Wild Animal, Avian...not Food, Plant as described by Enterobase
sources$New_Niche <-
  gsub('.*Eider.*', "Wild Animal__", sources$New_Niche)

sources$New_Type <-
  gsub('.*Eider.*', "Avian__", sources$New_Type)

# Arugula is Food, Plant
sources$New_Niche <- gsub("Arugula", "Food__", sources$New_Niche)
sources$New_Type <- gsub("Arugula", "Plant__", sources$New_Type)

# Eagle and gull to Wild Animal/Avian
sources$New_Niche <-
  gsub(".*Eagle.*|.*Gull.*", "Wild Animal__", sources$New_Niche)

sources$New_Type <-
  gsub(".*Eagle.*|.*Gull.*", "Avian__", sources$New_Type)

# Caecal samples all belonged to one of our poultry collections so will recode as Livestock, Poultry
sources$New_Niche <-
  gsub("Caecal", "Livestock__", sources$New_Niche)

sources$New_Type <- gsub("Caecal", "Poultry__", sources$New_Type)

# Bison ground meat is Livestock, Bovine
sources$New_Niche <-
  gsub(".*Bison.*", "Livestock__", sources$New_Niche)

sources$New_Type <- gsub(".*Bison.*", "Bovine__", sources$New_Type)

# Oyster
sources$New_Niche <-
  gsub(".*Oyster.*", "Livestock__", sources$New_Niche)

sources$New_Type <-
  gsub(".*Oyster.*", "Livestock Other__", sources$New_Type)

# Bind sources and our spreadsheets
sources <-
  rbind(sources,
        NDs,
        our_other,
        our_food,
        our_water,
        our_human,
        our_livestock)

# Check the counts
easy_source_count <-
  sources %>% group_by(New_Niche, New_Type) %>% summarise(counts = n())
easy_source_count <- arrange(easy_source_count, New_Niche, New_Type)

# Split out spreadsheet for manual curation by selecting samples for which New Niche hasn't yet been edited
manual_edits <- sources %>% filter(!grepl("__$", New_Niche))
auto_edits <- sources %>% filter(grepl("__$", New_Niche))

# Write the rows that require manual editing to a spreadsheet
write_csv(manual_edits, "raw_data/meta/ST58_manual_edits_raw.csv")

## The above file had New_Niche and New_Type entered manually, based on source details and NCBI data
## Twenty-one samples that could not be classified properly were given New_Niche and New_Type "EXCLUDE"

#Read in the processed manual edit sheet
manual_edits_processed <-
  read.csv("raw_data/meta/ST58_manual_edits_processed.csv")

# Remove excluded samples
manual_edits_processed <-
  manual_edits_processed %>% filter(New_Niche != "EXCLUDE")

#Bind manual edits to auto edits to produce final source metadata sheet
sources <- rbind(auto_edits, manual_edits_processed)

#Correct mistake where Type was Swine instead of Porcine
sources$New_Type <- gsub("Swine", "Porcine", sources$New_Type)

#Remove three sequences that didn't have cgMLST information ERR3333589, ERR3333590 and SRR9852831
nocgMLST <- c("ERR3333589", "ERR3333590", "SRR9852831")
sources <- sources %>% filter(!Working_name %in% nocgMLST)

#Remove sequence SRR3098930 where assembly is >6MB
sources <- sources %>% filter(!Working_name %in% c("SRR3098930"))

#Use sources to subset the original Enterobase metadata to just the strains we're still working with
filt_meta <- raw_data %>% filter(Barcode %in% sources$Barcode)

#Remove the underscores from the New_ columns and give them new names
sources$`Source Niche` <- gsub("__$", "", sources$New_Niche)
sources$`Source Type` <- gsub("__$", "", sources$New_Type)

#Select the columns we want from sources and the ones we want from Enterobase
clean_sources <-
  sources %>% select(Barcode, Working_name, `Source Niche`, `Source Type`)

clean_filt_meta <-
  filt_meta %>% select(Barcode, `Collection Year`, Continent, Country)

#Bind them together
clean_metadata <- left_join(clean_sources, clean_filt_meta)

#Rename the 'Working_name' column to just 'Name'
clean_metadata <-
  clean_metadata %>% select(Name = Working_name, everything())

# Eight sequences were excluded during Gubbins analysis due to missing data threshold
# Remove these from the final metadata sheet
out <- c("SRR1544291", "SRR8502430", "ESBL_Ecoli_6.20", "Vlees_ESBL_2013_0550", "ERR2206024", "SRR8502394", "OT_ESBL_0408", "SafeFoodERA_ESBL_119")

clean_metadata <- 
  clean_metadata %>% 
  filter(Name %notin% out) %>% 
  select(Name, Barcode, Niche = `Source Niche`,
         Source = `Source Type`, `Collection Year`, 
         Continent, Country)

# Read in column of data origin (In-house, collaborator or Public) and PMID so readers can find methods associated with in-house strains 
origins_PMID <-read_csv("raw_data/meta/ST58.dataorigin.PMID.csv")

clean_metadata <- left_join(clean_metadata, origins_PMID)

# Write the csv 
write_csv(clean_metadata, "raw_data/meta/ST58.metadata.csv")
