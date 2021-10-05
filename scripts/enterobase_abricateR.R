EB_abricateR <- function(genomelist, abricate_data, perc_ident = 90, perc_length = 95, name = "output"){
        
        require(magrittr)
        require(dplyr)
        require(ggplot2)
        require(reshape2)
  
        #Read in the full list of genomes. Note it is written in as a dataframe to allow joining with other dataframes later on
        #(redundant step, overwritten later).
        genome_list <- read_delim(
                genomelist,
                "\t",
                escape_double = FALSE,
                col_names = FALSE,
                trim_ws = TRUE,
                skip = 1
        )
        
        genome_list$X1 <- gsub(".* ESC","ESC",genome_list$X1)
        genome_list$X1 <- gsub(".fasta","",genome_list$X1)
        
        #change column name of the previously read in list.
        colnames(genome_list) <- "Assembly_barcode"
        
        #read in the abricate output
        genodata <- read_delim(abricate_data,
                               "\t",
                               escape_double = FALSE,
                               trim_ws = TRUE)
        
        #change column names to remove the "%" and "#" signs in the column names. These interfere later on when we try to call these columns
        colnames(genodata) <- gsub("%", "PERC_", colnames(genodata))
        colnames(genodata) <- gsub("#FILE", "Assembly_barcode", colnames(genodata))
        
        genodata$PERC_COVERAGE <- as.numeric(genodata$PERC_COVERAGE)
        genodata$PERC_IDENTITY <- as.numeric(genodata$PERC_IDENTITY)
        
        #Remove extra headers
        genodata <- genodata %>% filter(Assembly_barcode != "#FILE")
        
        #apply our filter to only allow genes which have a percent coverage of >=95 and a percent identity of >=90
        genodata <-
                genodata %>% filter(PERC_COVERAGE >= perc_length) %>% filter(PERC_IDENTITY >= perc_ident)
        
        #trim the path off of the rear of the assembly filenames
        genodata$Assembly_barcode <-
                gsub(".fasta", "", genodata$Assembly_barcode)
        
        # appends a 1 to the final column of each row (used later by dcast)
        rep(x = 1, times = nrow(genodata)) -> genodata$gene_present
        
        #cast the data to generate a table showing which genes are present and which genes are absent
        #i.e. make the list wide rather than long
        simple_summary <- reshape2::dcast(
                data = genodata,
                Assembly_barcode ~ GENE,
                value.var = 'gene_present',
                drop = FALSE
        )
        
        simple_summary <- left_join(genome_list, simple_summary)
        
        #create a new column for each group of genes and assign to each cell in this column the sum of the number of genes from each group present in a given sample
        simple_summary <-
                simple_summary %>% mutate(grp1_cvaABC_cvi = cvaA + cvaB + cvaC + cvi)
        simple_summary <-
                simple_summary %>% mutate(grp2_iroBCDEN = iroB + iroC + iroD + iroE + iroN)
        simple_summary <-
                simple_summary %>% mutate(grp3_iucABCD_iutA = iucA + iucB + iucC + iucD + iutA)
        simple_summary <-
                simple_summary %>% mutate(grp4_etsABC = etsA + etsB + etsC)
        simple_summary <-
                simple_summary %>% mutate(grp5_ompT_hlyF = ompT + hlyF)
        simple_summary <-
                simple_summary %>% mutate(grp6_sitABCD = sitA + sitB + sitC + sitD)
        
        #Save the row names to a different object so we can later reassign them after editing the table
        assembly_name <- simple_summary$Assembly_barcode
        
        #Take all rows and all columns except for the first column and reassign the dataframe "simple summary" as such
        simple_summary <-
                simple_summary[1:nrow(simple_summary), 2:ncol(simple_summary)]
        
        #Replace gene hits such that if any gene was detected more than once, or if more than one gene in a given group were detected, a corresponding gene/gene group is not counted twice
        simple_summary[simple_summary > 1] <- 1
        
        #reassign the Assembly_barcode column
        simple_summary$Assembly_barcode <- assembly_name
        
        #change the column order so that Assembly_barcode comes first.
        simple_summary <-
                simple_summary %>% select(Assembly_barcode, everything())
        
        #Add an extra column which sums the number of gene groups detected
        simple_summary <-
                simple_summary %>% mutate(
                        grp_sum = grp1_cvaABC_cvi + grp2_iroBCDEN + grp3_iucABCD_iutA + grp4_etsABC +
                                grp5_ompT_hlyF + grp6_sitABCD
                )
        
        #If 3 or less gene groups were detected within a given sample, replace this value with No (ColV plasmid not considered present)
        simple_summary$ColV <- gsub("1|2|3", "No", simple_summary$grp_sum)
        
        #If 4 or more gene groups were detected within a given sample, replace this value with Yes (ColV plasmid considered present)
        simple_summary$ColV <- gsub("4|5|6", "Yes", simple_summary$ColV)
        
        simple_summary$ColV <- gsub("No", "0", simple_summary$ColV)
        simple_summary$ColV <- gsub("Yes", "1", simple_summary$ColV)
        
        #Replace all NAs with zeroes, as these strains did not carry and genes from the ColV abricate screen at all and instead come up as NA
        #CHECK THIS
        simple_summary <- type.convert(x = simple_summary)
        simple_summary[is.na(simple_summary)] <- 0
                
        assign(paste("ColV", name, paste("N", perc_ident, sep = ""), paste("L", perc_length, sep =""), sep = "."), simple_summary, envir=globalenv())
        assign(paste("ColV", name, "genodata", sep = "."), genodata, envir=globalenv())

}
