####-----PROCESSING-----####
# Read in metadata sheet
meta_eb <- read_delim("raw_data/enterobase_colv/enterobase_raw_metadata.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

# Replace illegal chars with underscores
colnames(meta_eb) <- gsub(" ", "_", colnames(meta_eb))
colnames(meta_eb) <- gsub("/", "_", colnames(meta_eb))

# Filter out strains wher the source niche is unknown
meta_eb <- meta_eb %>% filter(Source_Niche != "ND")

# Filter out strains where the continent is unresolved
meta_eb <- meta_eb %>% filter(Continent != "Unresolved")

# Filter out strains where the Source_Niche is laboratory
meta_eb <- meta_eb %>% filter(Source_Niche != "Laboratory")

# Filter out strains where the species is non-E. coli
meta_eb <- meta_eb %>% filter(Species == "Escherichia coli")

# Note that the replacements below have been each been curated to ensure that 'rat blood' doesnt get misassigned to 'human_ExPEC'.
# Note that the below changes are case insensitive for the first letter of a given word.
# Assign anything with rat as other

meta_eb$Revised_Source_Details <- gsub(".*rat;", "Other", meta_eb$Source_Details)

# Assign goat related strains as other
meta_eb$Revised_Source_Details <- gsub(".*[Cc]aprine.*|.*[Gg]oat.*", "Other__", meta_eb$Revised_Source_Details)

# Assign sheep related strains as other
meta_eb$Revised_Source_Details <- gsub(".*[^Bb][Oo]vis].*|.*[Oo]vine.*|.*[Ss]heep.*", "Other__", meta_eb$Revised_Source_Details)

# Change anything with the Bos taurus, beef, cattle, cow, bovine or bovinae into Bovine
meta_eb$Revised_Source_Details <- gsub(".*[Bb]os [Tt]aurus.*|.*[Bb]eef.*|.*[Vv]eal.*|.*[Cc]attle.*|.*[Cc]ow.*|.*[Bb]ovine.*|[Bb]ovinae", "Bovine__", meta_eb$Revised_Source_Details)

# Change anything with Porcine, pork, pig, sus scrofa or sus domesti into Porcine
meta_eb$Revised_Source_Details <- gsub(".*[Ss]wine.*|.*[Pp]orcine.*|.*PORK.*|.*pork.*|.*Pork.*|.*[Pp]ig.*|.*[Ss]us [Ss]crofa.*|.*[Ss]us [Dd]omesti.*", "Porcine__", meta_eb$Revised_Source_Details)

# Change anything to do with poultry, gallus, egg, turkey, chicken or duck into poultry
meta_eb$Revised_Source_Details <- gsub(".*[Bb]roiler.*|.*[Pp]oultry.*|.*[Gg]allus.*|.*[Ee]gg.*|.*[Tt]urkey.*|.*[Cc]hicken.*|.*[Dd]uck.*", "Poultry__", meta_eb$Revised_Source_Details)

# Change Canis/Canine/Dog/Canius (some people have incorrect spelling) into Companion animal
meta_eb$Revised_Source_Details <- gsub(".*[Cc]anis.*|.*[Cc]anis].*|.*[Cc]anine.*|.*[Dd]og.*|.*[Cc]anius.*", "Companion_animal__", meta_eb$Revised_Source_Details)

# As above but for feline, cat (avoiding 'catheter') and felis
meta_eb$Revised_Source_Details <- gsub(".*[Ff]eline].*|.*[Cc]at[^heter].*|.*[Ff]elis.*", "Companion_animal__", meta_eb$Revised_Source_Details)

# As above but for horses
meta_eb$Revised_Source_Details <- gsub(".*[Ee]quine.*|.*[Hh]orse.*", "Companion_animal__", meta_eb$Revised_Source_Details)

# Find all urinary/urine associated strains, ExPEC and UPEC and label them as Human_ExPEC. Avoid words like "during" with the '[ary|e]'
meta_eb$Revised_Source_Details <- gsub(".*UTI.*|.*[Uu]rin[ary|e].*|.*ExPEC.*|.*UPEC.*", "Human_ExPEC__", meta_eb$Revised_Source_Details)

# A common pattern from uploaders is a certain format beginning with Homo sapiens and ending with respiratory or wound. Assign these as Human ExPEC
meta_eb$Revised_Source_Details <- gsub(".*Homo sapiens.*[respiratory|wound].*", "Human_ExPEC__", meta_eb$Revised_Source_Details)

# Avoiding "bloody diarrheoa", find all with blood, or the previous format and assign to Human_ExPEC
meta_eb$Revised_Source_Details <- gsub(".*[Bb]lood[^y].*|.*Homo sapiens.*Age.*male.*[Bb]lood$.*", "Human_ExPEC__", meta_eb$Revised_Source_Details)

# Assign bloody diarrheoa designations as Human_other
meta_eb$Revised_Source_Details <- gsub(".*[Bb]loody diar.*", "Human_other__", meta_eb$Revised_Source_Details)

# Assign more human feces strains to human_other
meta_eb$Revised_Source_Details <- gsub(".*Homo sapiens.*feces.*Age.*", "Human_other__", meta_eb$Revised_Source_Details)

# Special case, curated
meta_eb$Revised_Source_Details <- gsub(".*Human Clinical [^S].*", "Human_ExPEC__", meta_eb$Revised_Source_Details)

# Special case, curated
meta_eb$Revised_Source_Details <- gsub(".*Human Clinical S.*", "Human_other__", meta_eb$Revised_Source_Details)

# Assign more rat and mouse associated strains to Other
meta_eb$Revised_Source_Details <- gsub(".*[Mm]us mu.*|.*[Mm]ouse.*|.*[Rr]at.*|.*[Rr]odent.*", "Other__", meta_eb$Revised_Source_Details)

# Special case, curated
meta_eb$Revised_Source_Details <- gsub(".*homo.*", "Human_other__", meta_eb$Revised_Source_Details)

# Human strains yet to be categorised to Human_other
meta_eb$Revised_Source_Details <- gsub("human; Homo sapiens|^Homo sapiens$|^Human; Homo sapiens$", "Human_other__", meta_eb$Revised_Source_Details)

# Write the automatically processed metadata sheet so it can be manually curated
write_delim(meta_eb, file = "raw_data/enterobase_colv/enterobase_Rprocessed_metadata.txt", delim = "\t")

# Please note the manually curated dataset "enterobase_Manualprocessed_metadata.txt" is present in the raw_data/enterobase_colv folder and each entry that has been manually curated 
# has a '1' in the column 'Manually_curated"