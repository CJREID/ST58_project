####-----PACKAGES-----####
# Load required packages
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(ggtree)
library(RColorBrewer)

####-----SETUP-----####
# USEFUL FUNCTION FOR SUBSETTING
'%notin%' <- Negate('%in%')

# SET WORKING DIRECTORY - ADD THE PATH TO THE `ST58_project` REPOSITORY ON YOUR COMPUTER
# SCRIPT WILL NOT WORK IF YOU DO NOT CHANGE THIS
wrkdir <- "/Volumes/126451/WORK/Projects/ST58/revision/ST58_project"
setwd(wrkdir)

# SOURCE DATA PROCESSING SCRIPTS
source("scripts/gene_analysis.R")
source("scripts/pangenome_analysis.R")
source("scripts/enterobase_analysis.R")

# Metadata edit for ease of use with ggplot
meta2 <- metadata %>% 
  mutate(O_type = gsub("Other", "Other O Type", O_type),
         H_type = gsub("Other", "Other H Type", H_type),
         `F Plasmid` = gsub("Other", "Other F Plasmid", `F Plasmid`),
         `F Plasmid` = replace_na(`F Plasmid`, "No F Plasmid"),
         fimH = gsub("Other", "Other fimH", fimH),
         OH_type = gsub("\\*/O100\\*", "", OH_type),
         OH_type = gsub("Other", "Other Serotype", OH_type)) %>% 
  group_by_all() %>% dplyr::summarise(`Count` = n())
  
# Edited metadata for use with tree and some plots that don't require a count
meta3 <- meta2 %>% ungroup()

####-----DEFINE COLOURS-----####
#Define Source Type colours based off colours used for the tree
type_cols <- c(
  "Avian" = "#a2d1cd",
  "Bovine" = "#5d39a8",
  "Canine" = "#b2df8a",
  "Equine" = "#33a02c",
  "ExPEC" = "#cb56c7",
  "Feline" = "#4d4040",
  "Human Other" = "#f09424",
  "Livestock Other" = "#9083cb",
  "Other Mammal" = "#c28b4c",
  "Ovine" = "#cd9dae",
  "Plant" = "#c64a34",
  "Porcine" = "#55868c",
  "Poultry" = "#cccb51",
  "Wastewater" = "#b2436e",
  "Water" = "#1f78b4"
)

F_cols <- c("F18:-:B1" = "#8DD3C7",
  "F18:-:B16" = "#FFFFB3",
  "F18:-:B8" = "#BEBADA",
  "F18:A6:B1" = "#FB8072",
  "F2:-:B1" = "#80B1D3",
  "F24:-:B1" = "#FDB462",
  "F36:A4:B1" = "#B3DE69",
  "F57:A18:B23" = "#FCCDE5",
  "F70:-:-" = "#D9D9D9",
  "F89:A6:B16" = "#BC80BD", 
  "Other F Plasmid" = "#CCEBC5",
  "No F Plasmid" = "#ADACAC")

# Colours for pheatmap of Liu ColV genes
pheat_cols = list(`F Plasmid` = c("F18:-:B1" = "#8DD3C7",
                    "F18:-:B16" = "#FFFFB3",
                    "F18:-:B8" = "#BEBADA",
                    "F18:A6:B1" = "#FB8072",
                    "F2:-:B1" = "#80B1D3",
                    "F24:-:B1" = "#FDB462",
                    "F36:A4:B1" = "#B3DE69",                       
                    "F89:A6:B16" = "#BC80BD", 
                    "Other F Plasmid" = "#CCEBC5",
                    "No F Plasmid" = "white"),
                  `Source` = type_cols)

# Colours for Source designations used in the Enterobase ColV analysis
type_cols_colv <- c(
  "Bovine" = "#5d39a8",
  "Companion Animal" = "#4d4040",
  "ExPEC" = "#cb56c7",
  "Human Other" = "#f09424",
  "Porcine" = "#55868c",
  "Poultry" = "#cccb51"
)

# Colours for BAP groupings
bap_cols <- c("BAP1" = "#99ccff",
              "BAP2" = "#e3e386",
              "BAP3" = "#ffb59e",
              "BAP4" = "#e0a500",
              "BAP5" = "#7ed67e",
              "BAP6" = "#bcacdd",
              "0" = "black")

bap_cols_graph <- c("BAP1" = "#99ccff",
                    "BAP2" = "#e3e386",
                    "BAP3" = "#ffb59e",
                    "BAP4" = "#e0a500",
                    "BAP5" = "#7ed67e",
                    "BAP6" = "#bcacdd")

bap_cols_tree <- c("BAP1" = "#99ccff",
              "BAP2" = "#e3e386",
              "BAP3" = "#ffb59e",
              "BAP4" = "#e0a500",
              "BAP5" = "#7ed67e",
              "BAP6" = "#bcacdd",
              "BAP7" = NULL,
              "BAP8" = NULL,
              "0" = NULL)

# Annotation colours for the low SNP distance heatmap
anno_colors <- list(Source = c(
  "Avian" = "#a2d1cd",
  "Bovine" = "#5d39a8",
  "Canine" = "#b2df8a",
  "Equine" = "#33a02c",
  "Ovine" = "#cd9dae",
  "Porcine" = "#55868c",
  "Poultry" = "#cccb51",
  "Water" = "#1f78b4",
  "ExPEC" = "#cb56c7",
  "Human Other" = "#f09424"),
  ColV = c("ColV-Pos" = "#8dd3c7",
           "ColV-Neg" = "#c97461"))

# H_type colours
OH_vars <- unique(meta3$OH_type)
OH_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(OH_vars))
names(OH_cols) <- sort(OH_vars)

# Colours for the functional categories
func_cols <- colorRampPalette(brewer.pal(12, "Set3"))(20)

# O_type colours
O_vars <- unique(meta3$O_type)
O_cols <- colorRampPalette(brewer.pal(9, "Set1"))(9)
names(O_cols) <- sort(O_vars)

# H_type colours
H_vars <- unique(meta3$H_type)
H_cols <- colorRampPalette(brewer.pal(8, "Set1"))(12)
names(H_cols) <- sort(H_vars)

# fimH colours
fimH_vars <- unique(meta3$fimH)
fimH_cols <- colorRampPalette(brewer.pal(8, "Set1"))(10)
names(fimH_cols) <- sort(fimH_vars)

# ColV colours
colv_cols <- c( "Yes" = "#8dd3c7", "No" =  "white")

# Combined colours for Figure 2
heat_cols <- c(bap_cols_tree, type_cols, O_cols, H_cols, fimH_cols, F_cols, colv_cols)

####------FIGURE 1 - METADATA SUMMARY-------####
# Niche by Source
nicheplot <- ggplot(meta2, aes(fct_infreq(factor(Niche)), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  scale_x_discrete(name = "Niche")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 600), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 8, vjust = , face = "plain"), legend.position = "right") 

# Continent by Source
continentplot <- ggplot(meta2, aes(fct_infreq(factor(Continent)), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_x_discrete(name = "Continent")+
  scale_fill_manual(values = type_cols) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 600), n.breaks = 6)

# Years 2000+ by Source
year2kplot <- ggplot(subset(meta2, `Collection Year` %in% c(2000:2020)),
                     aes(as.character(`Collection Year`), `Count`)) + 
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 125), n.breaks = 6) +
  scale_x_discrete(name = "Collection Year")

# List plots for figure 1
type_plots <- list(nicheplot, continentplot)

# Arrange Year, Niche and Continent plots on one page
figure1 <- ggarrange(year2kplot, ggarrange(plotlist = type_plots,
                                           ncol = 2,
                                           labels =c("b)","c)"), legend = "none"),
                     nrow =2, labels = "a)",
                     common.legend = TRUE,
                     legend = "bottom")

# Save the plot as a figure
ggsave("Figure1_metadata.png",
       figure1, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 2 - PHYLOGENETIC TREE------####
# Define BAP group names for use with the groupOTU function
bap_groups <- list(BAP1 = unlist(c(metadata %>% filter(BAP == "BAP1") %>% select(Name))),
                   BAP2 = unlist(c(metadata %>% filter(BAP == "BAP2") %>% select(Name))),
                   BAP3 = unlist(c(metadata %>% filter(BAP == "BAP3") %>% select(Name))),
                   BAP4 = unlist(c(metadata %>% filter(BAP == "BAP4") %>% select(Name))),
                   BAP5 = unlist(c(metadata %>% filter(BAP == "BAP5") %>% select(Name))),
                   BAP6 = unlist(c(metadata %>% filter(BAP == "BAP6") %>% select(Name))))

# Group the tree by BAP clusters
grouped_tree <- groupOTU(gene_tree, bap_groups)

# Draw circular tree bound to metadata
CGA_tree <- ggtree(gene_tree, layout = "fan", open.angle = 7.5, size =.2) %<+% meta3

# Split tree data into separate dataframe
CGA_filt <- CGA_tree$data

# Extract sample names for each BAP group in order they appear in the tree
b1 <- CGA_filt %>% filter(BAP == "BAP1") %>% pull(label)
b2 <- CGA_filt %>% filter(BAP == "BAP2") %>% pull(label)
b3 <- CGA_filt %>% filter(BAP == "BAP3") %>% pull(label)
b4 <- CGA_filt %>% filter(BAP == "BAP4") %>% pull(label)
b5 <- CGA_filt %>% filter(BAP == "BAP5") %>% pull(label)
b6 <- CGA_filt %>% filter(BAP == "BAP6") %>% pull(label)

# Use names defined above to extract ancestral node numbers for each BAP group
bap1_node <- MRCA(CGA_tree, c(b1[1], b1[length(b1)]))
bap2_node <- MRCA(CGA_tree, c(b2[1], b2[length(b2)]))
bap3_node <- MRCA(CGA_tree, c(b3[1], b3[length(b3)]))
bap4_node <- MRCA(CGA_tree, c(b4[1], b4[length(b4)]))
bap5_node <- MRCA(CGA_tree, c(b5[1], b5[length(b5)]))
bap6_node <- MRCA(CGA_tree, c(b6[1], b6[length(b6)]))

# Create dataframe for BAP nodes
d <- data.frame(node=c(bap1_node, bap2_node, bap3_node, bap4_node, bap5_node, bap6_node), 
                type=c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"))

# Draw circular tree with BAP groups highlighted
CGA_tree <- rotate_tree(CGA_tree, 150) + 
  geom_highlight(data=d, alpha=.3, 
                 aes(node=node, fill=type),
                 extendto = 0.00288) +
  scale_fill_manual(values=bap_cols_tree)

# Extract names from meta3
meta_names <- meta3$Name

# Define metadata for tree
meta_heat <- meta3 %>% select(Source, O_type, H_type, fimH, `F Plasmid`, ColV)

# Set names as rownames for gheatmap
rownames(meta_heat) <- meta_names

# Draw Figure 2 ST58 phylogenetic tree with metadata
figure2 <- gheatmap(CGA_tree,
                    meta_heat,
                    font.size = 1.6,
                    hjust = .5,
                    colnames =TRUE,
                    colnames_position = "top",
                    colnames_angle = 58,
                    colnames_offset_y = 7.4,
                    width = .4,
                    offset = -.00007,
                    color = rgb(0, 0, 0, alpha = .1)) +
  scale_fill_manual(values = heat_cols, na.value = "white", name = "Data") + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

# NOTE: The legend of this figure has been manually rearranged for presentation
# Save figure 2
ggsave("Figure2_ST58_tree.png", 
       figure2, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 3 - BAP GROUPING STATS------####
# Plot of Source distributions within BAP groups
bapsource <- ggplot(meta3 %>% filter(BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9"), aes(BAP)) + geom_bar(aes(fill = Source)) +
  scale_fill_manual(values = type_cols) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 400), n.breaks = 6) +
  scale_x_discrete(name = NULL, labels = c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

# Plot of F plasmid type distributions within BAP groups
bapF <- ggplot(meta3 %>% filter(!is.na(`F Plasmid`), BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9"), aes(BAP)) + geom_bar(aes(fill = `F Plasmid`)) +
  scale_fill_manual(name = "F Plasmid", values = F_cols) +
  scale_y_continuous(name = NULL, expand = c(0, 0), limits = c(0, 400), n.breaks = 6) +
  scale_x_discrete(name = NULL, labels = c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

# Plot of serotype distributions within BAP groups
bapOH <- ggplot(meta3 %>% filter(BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9"), aes(BAP)) + geom_bar(aes(fill = OH_type)) +
  scale_fill_manual(values = OH_cols, name = "Serotype") +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 400), n.breaks = 6) +
  scale_x_discrete(name = NULL, labels = c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

# Plot of fimH distributions within BAP groups
bapfimH <- ggplot(meta3 %>% filter(BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9"), aes(BAP)) + geom_bar(aes(fill = fimH)) +
  scale_fill_manual(values = fimH_cols) +
  scale_y_continuous(name = NULL, expand = c(0, 0), limits = c(0, 400), n.breaks = 6) +
  scale_x_discrete(name = NULL, labels = c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

bap_colv <- ggplot(meta3 %>% filter(BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9"), aes(BAP)) + geom_bar(aes(fill = ColV)) +
  scale_fill_manual(values = c("Yes" = "#00BFC4FF", "No" = "#F8766DFF")) +
  scale_y_continuous(name = NULL, expand = c(0, 0), limits = c(0, 400), n.breaks = 6) +
  scale_x_discrete(name = NULL, labels = c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

# Plot both on one page
figure3 <- ggarrange(bapsource, bapF, bap_colv, bapOH, bapfimH,
                     ncol = 3, nrow = 2, 
                     labels =c("a)","b)", "c)", "d)", "e)"), 
                     hjust = c(-0.5, .8, .8, -0.5, .8),
                     legend = "right")

# Save the plot as a figure
ggsave("Figure3_BAP_slices.png", 
       figure3, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 4 - pCERC4 HEATMAP-------####
# NOTE THIS FIGURE IS A COMPOSITE OF THREE TREE-HEATMAPS ALIGNED TOGETHER
# Draw the tree
tree2 <- ggtree(grouped_tree, aes(color=group), branch.length = "none") %<+% metadata +
  scale_color_manual(values = bap_cols)

# Create vector of F types with strain names as rownames
F_types_ggtree <- metadata %>% select(`F Plasmid`)
rownames(F_types_ggtree) <- metadata$Name

# As above except with ColV presence/absence
colV_ggtree <- metadata %>% select(ColV)
rownames(colV_ggtree) <- metadata$Name

# Tree with F type bands
figure4_pt1 <- gheatmap(tree2,
                     data = F_types_ggtree,
                     font.size = 2,
                     hjust = 0,
                     colnames =FALSE,
                     width = 7,
                     offset = -153.5,
                     color = NULL) + 
  scale_fill_manual(values = F_cols, na.value = "white") +
  theme(legend.position = "none")

# Tree with ColV +/- bands
figure4_pt2 <- gheatmap(tree2,
                        data = colV_ggtree,
                        font.size = 2,
                        hjust = 0,
                        colnames =FALSE,
                        width = 3,
                        offset = -65,
                        color = NULL) + 
                        scale_fill_manual(values = c("Yes" = "#df03fc", "No" = "white")) +
                        theme(legend.position = "none")

# Tree with pCERC4 binned BLAST hits heatmaps
figure4_pt3 <- gheatmap(tree2, 
                        data = binned_hits,
                        font.size = 2,
                        hjust = 0,
                        colnames =FALSE,
                        width = 20,
                        offset = 10,
                        color = NULL) + 
      scale_fill_gradient(low = "white", high = "#8dd3c7", na.value = "white") +
      theme(legend.position = "none")

# Save images - PLEASE NOTE THIS FIGURE IN THE MANUSCRIPT HAS HAD ITS LEGENDS RE-DRAWN AND A SCHEMATIC MAP OF THE PCERC4 PLASMID ADDED MANUALLY
ggsave("Figure4_pt1_unedited_pCERC4_map.png",
       figure4_pt1, 
       path = "outputs/figures/", 
       device = "png",
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

ggsave("Figure4_pt2_unedited_pCERC4_map.png",
       figure4_pt2, 
       path = "outputs/figures/", 
       device = "png",
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

ggsave("Figure4_pt3_unedited_pCERC4_map.png",
       figure4_pt3, 
       path = "outputs/figures/", 
       device = "png",
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 5 - ST58 SOURCE AND ST DISTRIBUTION-------####
# ColV by Source in primary ST58 collection
colvplot <- ggplot(meta2, aes(fct_rev(fct_infreq(factor(ColV))), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 450)) +
  scale_x_discrete(name = "ST58 - ColV", labels = c("Positive", "Negative"))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "bottom") 

bap_source_plot <- ggplot(meta2 %>% filter(BAP == "BAP2" | BAP == "BAP6"), aes(BAP, `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 450), name = "Proportion") +
  scale_x_discrete(name = "ST58 - ColV")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "bottom") 

# Save as picture
ggsave("Figure5_SourceColV.png",
       colvplot, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 6 - SCOARY GENES WITH PHYLOGENY-------####
# Get gene presence absence from full Roary table
roary_filter2 <- roary_filter %>% select(Gene_Unique, `103764`:ncol(roary_filter))

# Select BAP2 significant genes from fastBAPS data
bap2map_A <- bap_data %>% filter(BAP == "BAP2") %>% select(Gene_Unique, Annotation, Odds_ratio, Benjamini_H_p)

# Pull out the presence/absence of the significant genes
bap2map <- left_join(bap2map_A, roary_filter2) %>% select(Gene = Gene_Unique, everything(), -Annotation, -Benjamini_H_p, -Odds_ratio)

# Get the gene names for rownames later
bap2genes <- c(bap2map$Gene)

# Transpose so genes are columns and de-select sample names
bap2map <- t(bap2map) %>% as.data.frame %>% slice(2:n())

# Replace Vx colnames with gene names
names(bap2map) <-bap2genes

# Select BAP6 significant genes from fastBAPS data
bap6map_A <- bap_data %>% filter(BAP == "BAP6") %>% select(Gene_Unique, Annotation, Odds_ratio, Benjamini_H_p)

# Pull out the presence/absence of the significant genes
bap6map <- left_join(bap6map_A, roary_filter2) %>% select(Gene = Gene_Unique, everything(), -Annotation, -Benjamini_H_p, -Odds_ratio)

# Get the gene names for rownames later
bap6genes <- c(bap6map$Gene)

# Transpose so genes are columns and de-select sample names
bap6map <- t(bap6map) %>% as.data.frame %>% slice(2:n())

# Replace Vx colnames with gene names
names(bap6map) <-bap6genes

# Create a column of simple gene names
find_paralogs <- bap2map_A %>% 
  mutate(Gene = str_extract(Gene_Unique, "[A-z]{3,4}")) %>%
  mutate(Gene = gsub("_", "", Gene))

# Split out over-represented 
paralogs_pos <- find_paralogs %>% filter(Odds_ratio >= 1) %>% select(-Odds_ratio)

# Create column designating this positive association
paralogs_pos$Association <- rep("Pos", nrow(paralogs_pos))

# Split out under-represented
paralogs_neg <- find_paralogs %>% filter(Odds_ratio < 1) %>% select(-Odds_ratio)

# Create column designating negative association
paralogs_neg$Association <- rep("Neg", nrow(paralogs_neg))

# Match gene names from both over- and under-represented groups
paralogs <- paralogs_neg %>% filter(paralogs_neg$Gene %in% paralogs_pos$Gene)

# Arrange in a meaningful format
paralogs <- paralogs %>% ungroup %>% group_by(Annotation, Gene, Association) %>% summarise(n()) %>% arrange(Gene) %>% select(Annotation, Gene)

# Create another tree object for plotting
tree3 <- ggtree(grouped_tree, aes(color=group), branch.length = "none") %<+% metadata +
  scale_color_manual(values = bap_cols)

# Replace the long 'group_nnnn' names with 'b' to denote alternative allele of a gene
names(bap2map) <- gsub("_group.*", "_b", names(bap2map))
names(bap6map) <- gsub("_group.*", "_b", names(bap6map))

# Make another version of the bap2 genes ordered the same as the tree
# Extract tip order from phylotree
# d <- fortify(gene_tree)
# dd <- subset(d, isTip)
# colorder <- dd$label[order(dd$y, decreasing=TRUE)]
# colordnum1 <-match(colorder,rownames(bap2map))
# bap2mapB <- bap2map[colordnum1,]
# 
# which(rownames(bap2mapB) == "SRR7707874")

# Plot gene tree with BAP3 over- and under-represented genes
bap2treemap <- gheatmap(tree3,
         data = bap2map,
         font.size = 2,
         hjust = 1,
         colnames = TRUE,
         colnames_position = "bottom",
         colnames_angle = 90,
         colnames_offset_y = -3,
         colnames_offset_x = ,
         width = 12,
         offset = ,
         color = NULL) + 
  scale_fill_manual(name = "Genes", values = c("1" = "#e3e386", "0" = "#D35FB7"), labels = c("Present", "Absent")) +
  theme(legend.position = "right") +ylim(-30, NA)

# Plot gene tree with BAP5 over- and under-represented genes
bap6treemap <- gheatmap(tree3,
         data = bap6map,
         font.size = 2,
         hjust = 1,
         colnames =TRUE,
         colnames_position = "bottom",
         colnames_angle = 90,
         colnames_offset_y = -3,
         colnames_offset_x = ,
         width = 5,
         offset = ,
         color = NULL) + 
  scale_fill_manual(name = "Genes", values = c("1" = "#5d39a8", "0" = "#E66100"), labels = c("Present", "Absent"))+
  theme(legend.position = "right") +ylim(-50, NA)

# Visualise together for fun
# ggarrange(bap2treemap, bap6treemap, 
#          ncol = 2, nrow = 1, labels =c("a)","b)"), common.legend = TRUE, legend = "none")

# Save BAP3 heatmap
ggsave("Figure6_scoary_BAP2.png",
       bap2treemap, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 7 - COMPARE ARGs and VAGs IN BAPS and ColV+/--------####
# Filter out single member BAP groups
meta_arg_bapfilt <- meta_arg %>% filter(BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9")
meta_vir_bapfilt <- meta_vir %>% filter(BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9")

# Generate boxplot of Total ARGs by BAP group
bap_argplot <- ggboxplot(meta_arg_bapfilt,
                         x = "BAP", y = "Total ARGs",
                         color = "BAP", palette = bap_cols_graph,
                         order = c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"),
                         ylab = "ARGs", xlab = "BAP Group")

# Generate boxplot of Total VAGs by BAP group
bap_virplot <-  ggboxplot(meta_vir_bapfilt, 
                          x = "BAP", y = "Total VAGs",
                          color = "BAP", palette = bap_cols_graph,
                          order = c("BAP1", "BAP2", "BAP3", "BAP4", "BAP5", "BAP6"),
                          ylab = "VAGs", xlab = "BAP Group")

# Generate boxplot of Total ARGs by ColV +/-
colv_argplot <-  ggboxplot(meta_arg, x = "ColV", y = "Total ARGs",
                           color = "ColV", palette = c("#92f7e6", "#f52020"),
                           order = c("Yes", "No"),
                           ylab = "ARGs", xlab = "ColV") +
  scale_x_discrete(label = c("Positive", "Negative"))

# Generate boxplot of Total VAGs by ColV +/-
colv_virplot <-  ggboxplot(meta_vir, x = "ColV", y = "Total VAGs",
                           color = "ColV", palette = c("#92f7e6", "#f52020"),
                           order = c("Yes", "No"),
                           ylab = "VAGs", xlab = "ColV")+
  scale_x_discrete(label = c("Positive", "Negative"))

# Kruskall Wallis tests indicate that there is a difference between BAP groups and ColV+/ColV- in terms of ARG and VAG carriage rates
kruskal.test(`Total ARGs` ~ BAP, data = meta_arg_bapfilt)
kruskal.test(`Total VAGs` ~ BAP, data = meta_vir_bapfilt)
kruskal.test(`Total ARGs` ~ ColV, data = meta_arg)
kruskal.test(`Total VAGs` ~ ColV, data = meta_vir)

# Pairwise Wilcoxon tests to determine the pairwise differences between BAP groups with BH adjustment for multiple testing
pairwise.wilcox.test(meta_arg_bapfilt$`Total ARGs`, meta_arg_bapfilt$BAP,
                     p.adjust.method = "BH")
pairwise.wilcox.test(meta_vir_bapfilt$`Total VAGs`, meta_vir_bapfilt$BAP,
                     p.adjust.method = "BH")

# Wilcoxon Rank-Sum test for difference between two groups (ColV+ and ColV-) for ARG and VAG carriage 
wilcox.test(x = as.numeric(unlist(meta_arg %>% filter(ColV == "Yes") %>% select(`Total ARGs`))), 
            y = as.numeric(unlist(meta_arg %>% filter(ColV == "No") %>% select(`Total ARGs`))),
            conf.int = TRUE)

wilcox.test(x = as.numeric(unlist(meta_vir %>% filter(ColV == "Yes") %>% select(`Total VAGs`))), 
            y = as.numeric(unlist(meta_vir %>% filter(ColV == "No") %>% select(`Total VAGs`))),
            conf.int = TRUE)

# Create a list of significantly different group-pairs as determined by Wilcoxon tests above
arg_comparisons <- list(c("BAP2", "BAP3"),
                        c("BAP2", "BAP4"),
                        c("BAP2", "BAP5"),
                        c("BAP2", "BAP6"))

vir_comparisons <- list(c("BAP1", "BAP2"),
                        c("BAP2", "BAP3"),
                        c("BAP2", "BAP4"),
                        c("BAP2", "BAP5"),
                        c("BAP2", "BAP6"))

# Plot with signficant differences annotated
bap_argplot_p <- bap_argplot + stat_compare_means(comparisons = arg_comparisons, method = "wilcox.test")    
bap_virplot_p <- bap_virplot + stat_compare_means(comparisons = vir_comparisons, method = "wilcox.test")    

colv_argplot_p <- colv_argplot + stat_compare_means(method = "t.test", label.x = 1.25)    
colv_virplot_p <- colv_virplot + stat_compare_means(method = "t.test", label.x = 1.25) 

# Plot them all together
figure7 <- ggarrange(bap_argplot_p, bap_virplot_p, colv_argplot_p, colv_virplot_p, 
                     ncol = 2, nrow = 2, labels =c("a)","b)","c)","d)"), common.legend = TRUE, legend = "none")

# Save figure
ggsave("Figure7_ARGsVAGs.png",
       figure7, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 8 - SOURCE DISTRIBUTION COLV IN ENTEROBASE-------####
# Plot Source vs ColV carriage, coloured by Source
source_v_colvperc <- ggplot(entero_source %>% filter(ColV == "1"), aes(x = reorder(Source, -`Source Carriage`), y =`Source Carriage`))+
  geom_col(aes(fill = Source)) +
  theme_classic() +
  scale_fill_manual(name = "Source", 
                    labels = c("Bovine", "Companion Animal", "ExPEC", "Human Other", "Porcine", "Poultry"),
                    values = type_cols_colv) +
  theme(axis.text.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = .5, size = 16, face = "plain"),
        legend.position = "none") + 
  scale_y_continuous(name = "ColV Carriage Rate (%)", expand = c(0, 0), limits = c(0, 100), breaks = c(seq(0,100,20)))+
  scale_x_discrete(name = "Source",
                   labels = c("Poultry",
                              "Porcine",
                              "ExPEC",
                              str_wrap("Companion Animal", 10),
                              str_wrap("Human Other", 5),
                              "Bovine"))

# ColV by Source in Enterobase genome collection
colventero_plot <- ggplot(geno_meta, aes(x = fct_rev(fct_infreq(factor(recode(ColV, "0" = "Negative", "1" = "Positive"))))))+ 
  geom_bar(aes(fill = `Source`))+
  scale_fill_manual(values = type_cols_colv) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 30000)) +
  scale_x_discrete(name = "Enterobase - ColV")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# Plot ST vs ColV Carriage, coloured by Source
ST_v_colvperc <- ggplot(colV_summary, aes(x = reorder(ST, -`ColV ST Carriage`), y =`ColV Source Carriage`))+
  geom_col(aes(fill = Source))+
  theme_classic()+
  scale_fill_manual(name = "Source", 
                    labels = c("Bovine", "Companion Animal", "ExPEC", "Human Other", "Porcine", "Poultry"),
                    values = type_cols_colv) +
  theme(axis.text.x = element_text(color = "grey20", size = 7, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 7, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = .5, size = 16, face = "plain"),
        legend.position = "right") + 
  scale_y_continuous(name = "ColV Carriage Rate (%)", expand = c(0, 0),limits = c(0, 100), breaks = c(seq(0,100,20)))+
  scale_x_discrete(name = "Sequence Type")

# Plot absolute counts of the top ColV carrying STs and fill by Source
ST_count <- ggplot(geno_ST, aes(x = fct_infreq(factor(ST))))+
  geom_bar(aes(fill = Source))+
  theme_classic()+
  scale_fill_manual(name = "Source", 
                    labels = c("Bovine", "Companion Animal", "ExPEC", "Human Other", "Porcine", "Poultry"),
                    values = type_cols_colv) +
  theme(axis.text.x = element_text(color = "grey20", size = 7, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 7, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = .5, size = 16, face = "plain"),
        legend.position = "right") + 
  scale_y_continuous(name = "Count", expand = c(0, 0),limits = c(0, 600))+
  scale_x_discrete(name = "Sequence Type")

# Make a list of plots to be drawn together
colv_plot_list <- list(source_v_colvperc, colventero_plot, ST_v_colvperc, ST_count)

# Draw Figure 6
figure8 <- ggarrange(plotlist = colv_plot_list, 
                     ncol = 2,
                     nrow=2,
                     common.legend = TRUE,
                     labels =c("a)","b)", "c)", "d)"),
                     legend = "bottom")

ggsave("Figure8_EnterobaseColV.png",
       figure8, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------FIGURE 9 - LOW SNP DISTANCE HEATMAP-------####
# Read SNP data in - assumes you have a blank top left cell
core_snps.A <- read_csv(file = "raw_data/snp-dists/ST58_core_SNP_dists.csv") %>% rename(Name = ...1)

# Extract names
rows <- core_snps.A$Name

# Process core.A SNPs
melt.core.A <- reshape2::melt(core_snps.A)

# Combine metadata for LHS strains and rename columns
melt.core.A1 <- left_join(melt.core.A, metadata) %>% 
  select(Name, variable, value, Source, Country, `Collection Year`, ColV) %>% 
  rename(Name1 = Name,
         Name = variable,
         Source1 = Source,
         Country1 = Country,
         ColV1 = ColV,
         `Collection Year1` = `Collection Year`)

# Combine metadata for RHS strains and rename columns, then filter for Human in Source1 and non-human in Source2
melt.core.A2 <- left_join(melt.core.A1, metadata) %>%
  select(Name1, Name, Source1, Source, Country1, Country, `Collection Year1`, `Collection Year`, ColV1, ColV, value)%>% 
  rename(Name2 = Name, Source2 = Source, Country2 = Country, `Collection Year2` = `Collection Year`, ColV2 = ColV, `Core SNPs` = value)%>%
  filter(Source1 %in% c("ExPEC", "Human Other")) %>%
  filter(Source2 %notin% c("ExPEC", "Human Other"))

# Recode ColV presencer/absence for readability
melt.core.A2$ColV1 <- melt.core.A2$ColV1 %>% dplyr::recode(Yes = "ColV-Pos", No = "ColV-Neg")
melt.core.A2$ColV2 <- melt.core.A2$ColV2 %>% dplyr::recode(Yes = "ColV-Pos", No = "ColV-Neg")

# Re-order columns
melt.core.A2 <- melt.core.A2 %>% select(Name1, Source1, Country1, `Collection Year1`, ColV1,
                            `Core SNPs`,
                            Name2, Source2, Country2, `Collection Year2`, ColV2)

# Mutate the strain column to contain condensed info
SNP_pairs <- melt.core.A2 %>% dplyr::mutate(`Strain 1` = paste(Name1, Source1, Country1, `Collection Year1`, ColV1, sep = " - "),
                                      `Strain 2` = paste(Name2, Source2, Country2, `Collection Year2`, ColV2, sep = " - "))


# Organise the columns
SNP_pairs <- SNP_pairs %>% select(`Strain 1`, `Core SNPs`, `Strain 2`)

# Filter Core SNPs equal to or less than 20 and remove self-comparisons
SNP_int <- SNP_pairs %>%
  filter(`Core SNPs` <= 20, `Strain 1` != `Strain 2`)

# Define the annotation rows for heatmap
anno_rows <- melt.core.A2 %>% 
  filter(`Core SNPs` <= 20) %>% 
  select(Name1, Source = Source1, 
         ColV = ColV1) %>% distinct()

# Convert Name column to rownames
rownem <- anno_rows$Name1
anno_rows <- anno_rows %>% select(-Name1)
rownames(anno_rows) <- rownem

# As above except for columns
anno_cols_snp <- melt.core.A2 %>% 
  filter(`Core SNPs` <= 20) %>%
  select(Name2, 
         Source = Source2, 
         ColV = ColV2) %>% distinct()

# Convert name to rownames
colnem <-anno_cols_snp$Name2
anno_cols_snp <- anno_cols_snp %>% select(-Name2)
rownames(anno_cols_snp) <- colnem

# Create DF of paired names and SNP counts
paircast <- melt.core.A2 %>% select(Name1, `Core SNPs`, Name2) %>% filter(`Core SNPs` <= 20)

# Cast to SNP matrix
paircast2 <- reshape2::dcast(paircast, Name1 ~ Name2, value.var = "Core SNPs")

# Move Name column to rownames and convert to matrics
SNP_rows <- paircast2$Name1
paircast_mat <- as.matrix(paircast2 %>% select(-Name1))
rownames(paircast_mat) <- SNP_rows

# Set NAs outside the SNP range so they can be ignored
paircast_mat[is.na(paircast_mat)] <- 100

# Make quick heatmap with clustering to group the most closely related strains
out <- pheatmap(mat = paircast_mat,
                cluster_cols = TRUE,
                cluster_rows = TRUE)

# Extract row and column orders from the above heatmap
roword <- rownames(paircast_mat[out$tree_row[["order"]],])
colord <- colnames(paircast_mat[,out$tree_col[["order"]]])

# Reorder matrix
paircast_mat2 <- paircast_mat[roword,]
paircast_mat2 <- paircast_mat2[,colord]

# Reset 100 values as NAs
paircast_mat2[paircast_mat2 == 100] <- NA

# Create Figure 9
figure9 <- pheatmap(mat = paircast_mat2,
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    annotation_col = anno_cols_snp,
                    annotation_row = anno_rows,
                    color = c(brewer.pal(n = 9, "BuPu")),
                    breaks = c(seq(from = 0, to = 20, by = 2)),
                    fontsize_col = 6,
                    fontsize_row =6,
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    annotation_colors = anno_colors,
                    annotation_names_col = TRUE,
                    annotation_names_row = TRUE,
                    legend = TRUE,
                    annotation_legend = TRUE,
                    border_color = "black",
                    na.cols = "grey")

# Save figure - PLEASE NOTE THIS FIGURE IN THE MANUSCRIPT HAS HAD ITS LEGENDS RE-DRAWN MANUALLY FOR CONSISTENCY OF PRESENTATION
ggsave("Figure9.lowSNPheatmap.png", 
       figure9, 
       path = "outputs/figures/", 
       device = "png", 
       width = 250,
       height = 180,
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 1 - EXPANDED COLLECTION YEARS-------####
# Year by Source
yearplot <- ggplot(meta2, aes(`Collection Year`, `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 125), n.breaks = 6)+
  scale_x_continuous(name = "Collection Year", n.breaks = 30)+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# Save figure
ggsave("FigureS1_allyears.png",
       yearplot, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 150, 
       unit ="mm", 
       dpi = 600) 

####------SUPP FIGURE 2 - OH TYPES BY SOURCE AND COLV-------####
# OH by Source
oh_sourceplot <- ggplot(meta2, aes(fct_infreq(factor(OH_type)), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300), n.breaks =6) +
  scale_x_discrete(name = "Serotype") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 12, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# OH by ColV
oh_colvplot <- ggplot(meta2, aes(fct_infreq(factor(OH_type)), `Count`)) +
  geom_bar(aes(fill = `ColV`), stat="identity") +
  scale_fill_manual(values = c("No" = "#F8766DFF", "Yes" = "#00BFC4FF"))+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300), n.breaks =6) +
  scale_x_discrete(name = "Serotype") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 12, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# Plot together
figureS2 <- ggarrange(oh_sourceplot, oh_colvplot, ncol = 2, labels =c("a)","b)"), legend = "bottom")

# Save figure
ggsave("FigureS2_OH.png", 
       figureS2, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 150, 
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 3 - FIMH BY SOURCE AND COLV-------####
# fimH by Source
fimH_sourceplot <- ggplot(meta2, aes(fct_infreq(factor(fimH)), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_x_discrete(name="fimH")+
  scale_fill_manual(values = type_cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350), n.breaks =6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# fimH by ColV
fimH_colvplot <- ggplot(meta2, aes(fct_infreq(factor(fimH)), `Count`)) +
  geom_bar(aes(fill = `ColV`), stat="identity") +
  scale_fill_manual(values = c("No" = "#F8766DFF", "Yes" = "#00BFC4FF")) +
  scale_x_discrete(name="fimH")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 350), n.breaks = 5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# Plot together
figureS3 <- ggarrange(fimH_sourceplot, fimH_colvplot, ncol = 2, labels =c("a)","b)"), legend = "bottom")

# Save figure
ggsave("FigureS3_fimH.png", 
       figureS3, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 150, 
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 4 - F-pMLST BY SOURCE AND COLV-------####
# F-plasmid pMLST by ColV proportional
ftype_colv_absplot <- ggplot(meta2, aes(fct_infreq(factor(`F Plasmid`)), `Count`)) +
  geom_bar(aes(fill = `ColV`), position="stack", stat="identity") + 
  scale_fill_manual(values = c("No" = "#F8766DFF", "Yes" = "#00BFC4FF")) +
  scale_x_discrete(name = "F plasmid pMLST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 300), n.breaks = 6) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),
        legend.position = "right") 

# F-plasmid pMLST by ColV absolute
ftype_colv_percplot <- ggplot(meta2, aes(fct_infreq(factor(`F Plasmid`)), `Count`)) +
  geom_bar(aes(fill = `ColV`), position="fill", stat="identity") + 
  scale_x_discrete(name = "F plasmid pMLST")+
  scale_fill_manual(values = c("No" = "#F8766DFF", "Yes" = "#00BFC4FF")) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 1), n.breaks = 5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),
        legend.position = "right") 

# F-plasmid pMLST by Source proportional
ftype_source_percplot <- ggplot(meta2, aes(fct_infreq(factor(`F Plasmid`)), `Count`)) +
  geom_bar(aes(fill = `Source`), position="fill", stat="identity") + 
  scale_fill_manual(values = type_cols) +
  scale_x_discrete(name = "F plasmid pMLST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 1), n.breaks = 5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),
        legend.position = "right") 

# F-plasmid pMLST by Source absolute
ftype_source_absplot <- ggplot(meta2, aes(fct_infreq(factor(`F Plasmid`)), `Count`)) +
  geom_bar(aes(fill = `Source`), position="stack", stat="identity") + 
  scale_fill_manual(values = type_cols) +
  scale_x_discrete(name = "F plasmid pMLST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 250), n.breaks = 5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),
        legend.position = "right") 

# Plot together
figureS4 <- ggarrange(ftype_source_absplot, ftype_colv_absplot, ncol = 2, labels =c("a)","b)"), legend = "bottom")

# Save figure
ggsave("FigureS4_Ftype.png", 
       figureS4, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 150, 
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 5 - LIU COLV GENE PRESENCE/ABSENCE W/ COLBM MARKERS-------####
# Select ColV positive strains
meta_colv_pos <- meta_colv %>%
  filter(ColV_sum == 1) %>%
  mutate(Total_genes = rowSums(across(cvaA:sitD)),
         `F Plasmid` = gsub("Other", "Other F Plasmid", `F Plasmid`),
         `F Plasmid` = replace_na(`F Plasmid`, "No F Plasmid"))

colbm <- meta_vir %>% select(Name, cba, cma) %>% filter(Name %in% meta_colv_pos$Name)

meta_colv_pos <- left_join(meta_colv_pos, colbm)

# Extract names
r_nam <- meta_colv_pos$Name

# Select genes
genes_only <- meta_colv_pos %>% select(cvaA:sitD, cba, cma)

genes_n_meta <- meta_colv_pos %>% select(cvaA:sitD, cba, cma)

# Select sources for annotation
source <- meta_colv_pos %>% dplyr::select(`Source`)

# Select F types for annotation
f_type <- meta_colv_pos %>% dplyr::select(`F Plasmid`)

# Assign rownames
rownames(source) <- r_nam
rownames(genes_only) <- r_nam
rownames(f_type) <- r_nam

# Combine Source and F-types
anno_cols = cbind(`Source` = source, `F Plasmid` = f_type)
rownames(anno_cols)<-r_nam

# Transpose gene presence/absence matrix
flip_genes <- t(genes_only)

# Plot heatmap
figureS5 <- pheatmap(flip_genes,
                     cluster_cols = TRUE,
                     cluster_rows = FALSE,
                     annotation_col = anno_cols,
                     color = c("white", "#8dd3c7"),
                     breaks = c(0,0.5,1),
                     fontsize_col = 8,
                     show_rownames = TRUE,
                     show_colnames = FALSE,
                     annotation_names_col = TRUE,
                     annotation_colors = pheat_cols,
                     legend = FALSE,
                     annotation_legend = TRUE,
                     border_color = "black")

# Save figure - PLEASE NOTE THIS FIGURE IN THE MANUSCRIPT HAS HAD ITS LEGENDS RE-DRAWN MANUALLY
ggsave("FigureS5_colV_genes_raw.png",
       figureS5, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)


####------SUPP FIGURE 6 - COLV CARRIAGE BY SOURCE-------####
# Proportional Source by ColV
sourcecolv_percplot <- ggplot(meta2, aes(fct_infreq(factor(Source)), `Count`)) +
  geom_bar(aes(fill = `ColV`), position="fill", stat="identity") + 
  scale_fill_manual(values = c("No" = "#F8766DFF", "Yes" = "#00BFC4FF")) +
  scale_x_discrete(name = "Source")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), limits = c(0, 1), n.breaks = 5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# Absolute Source by ColV
sourcecolv_absplot <- ggplot(meta2, aes(fct_infreq(factor(`Source`)), `Count`)) +
  geom_bar(aes(fill = `ColV`), position="stack", stat="identity") + 
  scale_fill_manual(values = c("No" = "#F8766DFF", "Yes" = "#00BFC4FF")) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 300), n.breaks = 6) +
  scale_x_discrete(name = "Source")+
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right") 

# Plot together
figureS6 <- ggarrange(sourcecolv_absplot, sourcecolv_percplot, 
                      ncol = 2, labels =c("a)","b)"), 
                      legend = "bottom", common.legend = TRUE)

# Save figure
ggsave("FigureS6_SourceColV.png", 
       figureS6, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 150, 
       unit ="mm", 
       dpi = 600)


####------SUPP FIGURE 7 - BAP6 SCOARY HEATMAP-------####
# Save BAP6 heatmap - bap6treemap object originates in FIGURE 6 - SCOARY GENES WITH PHYLOGENY section of this script
ggsave("FigureS7_scoary_BAP5.png",
       bap6treemap, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 8 - ARG HEATMAP-------####
# Define colours for presence and absence
heat_names <- 
  factor(c("Present", "Absent"), levels = c("Present", "Absent"))
heat_cols <- 
  factor(c("#8dd3c7", "#ededed"), levels = c("#8dd3c7", "#ededed"))

# Draw tree with ColV tips
heat_t <- 
  ggtree(gene_tree, branch.length="none") %<+% metadata  +
  geom_tippoint(aes(color = ColV, hjust=0, offset = 2), size =.7) +
  scale_color_manual(values = c("Yes" = "#df03fc", "No" = "white"), 
                     guide = guide_legend(reverse = TRUE, title = "ColV-Pos"))
meta_tree <- metadata %>% filter(BAP != "BAP7" & BAP != "BAP8" & BAP != "BAP9")

# Draw tree with BAP tips
heat_t_bap <- ggtree(gene_tree, branch.length="none") %<+% meta_tree +
  geom_tippoint(aes(color = BAP, hjust=0, offset = 2), size =.7) +
  scale_color_manual(values = bap_cols, 
                     guide = guide_legend(title = "BAP Group"))

# Draw ARG heatmap
figureS8 <- 
  gheatmap(heat_t_bap, arg_heat, width = 40, colnames = T, colnames_angle = 65, colnames_offset_y = 1.5, 
         colnames_offset_x = 0, offset = 5, colnames_position = "top", font.size = 2, hjust=0, color = NULL) +
          ylim(NA, 800) + 
  scale_fill_manual(values = c("#8dd3c7", "black"), 
                    breaks = heat_names, guide = guide_legend(title = "Genes"))

# Save figure
ggsave("FigureS8_arg_heatmap.png",
       figureS8, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 9 - VIRULENCE HEATMAP-------####
# Draw virulence heatmap
figureS9 <- 
  gheatmap(heat_t_bap, vir_heat, width = 40, colnames = T, colnames_angle = 65, colnames_offset_y = 1.5, 
         colnames_offset_x = 0, offset = 5, colnames_position = "top", font.size = 2.3, hjust=0, color = NULL) +
  ylim(NA, 800) + scale_fill_manual(values = c("#8dd3c7", "black"), breaks = heat_names, guide = guide_legend(title = "Genes"))

# Save figure
ggsave("FigureS9_vir_heatmap.png",
       figureS9, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 10 - PLASMID HEATMAP-------####
# Draw plasmid replicon heatmap
figureS10 <- 
  gheatmap(heat_t_bap, plas_heat, width = 20, colnames = T, colnames_angle = 65, colnames_offset_y = 1.5, 
         colnames_offset_x = 0, offset = -5, colnames_position = "top", font.size = 2, hjust=0, color = NULL) +
  ylim(NA, 800) + scale_fill_manual(values = c("#8dd3c7", "black"), breaks = heat_names, guide = guide_legend(title = "Genes"))

# Save figure
ggsave("FigureS10_plas_heatmap.png",
       figureS10, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####------SUPP FIGURE 11 - PAIRWISE SNP DISTANCE MATRIX-------####
# Melt all SNP values into one column using core_snps.A from Figure 9 section
core_snp_distro <- reshape2::melt(core_snps.A)

# Calculate mean
mean(core_snp_distro$value)

# Plot distribution of core SNP values
# ggplot(core_snp_distro, aes(x=value)) + geom_density()

# Extract tip order from phylotree
d <- fortify(gene_tree)
dd <- subset(d, isTip)
colorder <- dd$label[order(dd$y, decreasing=TRUE)]

# Fix problem strain names
core_snps.A$Name <- gsub("ESBL_E.coli_2.64", "ESBL_Ecoli_2.64", core_snps.A$Name)
rows <- core_snps.A$Name

# Back to snps processing
core_snps.A <- core_snps.A %>% select(-Name)
rownames(core_snps.A) <- rows
colordnum<-match(colorder,names(core_snps.A))
core_snps.A <- core_snps.A[,colordnum]

# Use this to quickly compare pairwise SNP between two sequences
# core_snps.A["ERR2019093", "SS_2_G2"]

# Write reordered SNPs to file
# write.csv(core_snps.A, "outputs/reorder.snp-dists.csv")

# Extract SNP distance ranges
mo1k<-core_snps.A[] >1000
mo250<-core_snps.A[] >250 & core_snps.A[] <=1000
mo100<-core_snps.A[] > 100 & core_snps.A[] <=250
les100<-core_snps.A[] <= 100 & core_snps.A[] > 75
les75<-core_snps.A[] <= 75 & core_snps.A[] > 50
les50<-core_snps.A[] <= 50 & core_snps.A[] > 25
les25<-core_snps.A[] <= 25 & core_snps.A[] > 15
les15<-core_snps.A[] <= 15

# Convert all columns to character class
column_names <- names(core_snps.A)
core_snps.A[column_names] <- lapply(core_snps.A[column_names], as.character)

# Assign character values to each range
core_snps.A[mo1k] <- "1000+"
core_snps.A[mo250] <- "250-1000"
core_snps.A[mo100] <- "100-250"
core_snps.A[les100] <- "76-100"
core_snps.A[les75] <- "51-75"
core_snps.A[les50] <- "26-50"
core_snps.A[les25] <- "16-25"
core_snps.A[les15] <- "0-15"

# Restore rownames
rownames(core_snps.A) <- rows

# SNP range breaks
SNP_brks<-c("0-15",
        "16-25",
        "26-50",
        "51-75",
        "76-100",
        "100-250",
        "250-1000",
        "1000+") 

# SNP range colours
SNP_cols<-c(brewer.pal(n = 8, name = "BuPu"))

# Generate tree with tips coloured by Source Type
tr <- ggtree(gene_tree, branch.length="none") %<+% metadata  +
  geom_tippoint(aes(color = BAP, hjust=0,offset = 2), size =.7) +
  scale_color_manual(values = bap_cols_tree)

# Draw heatmap
figureS11 <- gheatmap(tr, core_snps.A, width = 50, colnames = T, colnames_angle = 75, colnames_offset_y = 0, 
         colnames_offset_x = 0, offset = 10, colnames_position = "top", font.size = .2, hjust=0, color = NULL) +
  ylim(NA, 800) + scale_fill_manual(values=SNP_cols, breaks = rev(SNP_brks), name = "SNP count")

# Save figure
ggsave("FigureS11_SNP_heatmap.png", 
       figureS11, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 1200)

####------EXTRA PLOTS-------####
# H by Source
hplot <- ggplot(meta2, aes(fct_infreq(factor(`H_type`)), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  scale_x_discrete(name = "H Type") +
  scale_y_continuous(name = "H Type", expand = c(0, 0), limits = c(0, 300), n.breaks = 6) +
  theme_classic()

# O by Source
oplot <- ggplot(meta2, aes(fct_infreq(factor(O_type)), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  scale_x_discrete(name = "O Type") +
  scale_y_continuous(name = "O Type",expand = c(0, 0), limits = c(0, 250), n.breaks = 5) +
  theme_classic()

# Source
sourceplot <- ggplot(meta2, aes(fct_infreq(factor(Source)), `Count`)) +
  geom_bar(aes(fill = `Source`), stat="identity") +
  scale_fill_manual(values = type_cols) +
  scale_x_discrete(name = "Source") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 275)) +
  theme_classic()

# IncI1
source_g15 <- meta_plas %>% group_by(Source) %>% summarise(n()) %>% filter(`n()` >15)

meta_plas$IncI1 <- meta_plas$IncI1 %>% recode("1" = "Yes", "0" = "No")

I1_bap_perc <- meta_plas %>% filter(!is.na(IncI1)) %>% 
  ggplot(aes(x=BAP, fill = IncI1)) +
  geom_bar(position = "fill")

I1_bap_ST <- meta_plas %>% filter(!is.na(I1_simple)) %>% 
  ggplot(aes(x=BAP, fill = I1_simple)) +
  geom_bar(position = "fill")

I1_bap_abs <- meta_plas %>% filter(!is.na(I1_simple)) %>% 
  ggplot(aes(x=BAP, fill = I1_simple)) +
  geom_bar(position = "stack")
  


