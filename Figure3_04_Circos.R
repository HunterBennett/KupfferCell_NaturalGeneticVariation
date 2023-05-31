# Hunter Bennett
# Glass Lab | Kupffer strains project | April 11 2023
# Visualize signaling networks in the liver identified by NicheNet using Circos plots

# clean environment (if desired)
rm(list = ls())

# change working directory
setwd('~/strains_github/results/Figure3/')
getwd()

# load packages
library(nichenetr)
library(tidyverse)
library(circlize)

#### USER INPUT REQUIRED. 1. set global analysis variables ####
out_prefix = 'control_kupffer'
out_directory = './Circos/'
dir.create(out_directory)
target_membership_path = './nichenet_circos_diff_gene_membership.txt'
c57_active_target_links_path = './C57BL6J/c57bl6_kupffer_target_spec_circos_ligand_target_links_df.txt'
balb_active_target_links_path = './BALBcJ/balb_kupffer_target_spec_circos_ligand_target_links_df.txt'
aj_active_target_links_path = './AJ/aj_kupffer_target_spec_circos_ligand_target_links_df.txt'

# import linkage data frame for circos plot
c57_active_ligand_target_links_df = as_tibble(read.table(c57_active_target_links_path,
                                                     sep='\t',
                                                     header=TRUE))
balb_active_ligand_target_links_df = as_tibble(read.table(balb_active_target_links_path,
                                                         sep='\t',
                                                         header=TRUE))
aj_active_ligand_target_links_df = as_tibble(read.table(aj_active_target_links_path,
                                                         sep='\t',
                                                         header=TRUE))
# combine data frames
active_ligand_target_links_df = rbind(c57_active_ligand_target_links_df, balb_active_ligand_target_links_df, aj_active_ligand_target_links_df)

# read in group membership
target_membership = as_tibble(read.table(target_membership_path,
                                         sep='\t',
                                         header=TRUE))

target_membership$target <- toupper(target_membership$target)

# merge group membership
active_ligand_target_links_df = active_ligand_target_links_df %>% left_join(target_membership)

# filter to select ligands of interest
ligands_of_interest = c('LEP', 'APP', 'PLG', 'APOE', 'BMP2', 'ADAM17')
active_ligand_target_links_df = active_ligand_target_links_df %>% filter(ligand %in% ligands_of_interest)

# filter low affinity interactions
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.2) # select cutoff
# cutoff_include_all_ligands = 0

# remove ligands and targets
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

# prepare circos visualization #
active_ligand_target_links_df
active_ligand_target_links_df$target_group %>% unique()

#### USER INPUT REQUIRED. 2. set color dict variables. should only have to change for different cell types ####
# set color scheme
grid_col_ligand =c("hep-specific" = "#8c510a",
                   "lsec-specific" = "#dfc27d",
                   "stellate-specific" = "#35978f",
                   "hep-lsec" = "#bf812d",
                   "hep-stel" = "#fdae6b",
                   "lsec-stel" = "#80cdc1",
                   "general-ligands" = "#bdbdbd")

grid_col_target =c(
  "aj_union" = "#e41a1c",
  "aj_c57_shared" = "#dd1c77",
  "aj_balb_shared" = "#756bb1",
  "balb_union" = "#377eb8",
  "balb_c57_shared" = "#43a2ca",
  "c57_union" = "#31a354")

# extract the lists into tibbles
grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_group = grid_col_target %>% names(), color_target_type = grid_col_target)

# add coloring to circos links
circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand, target, weight)

# pull colors to list
ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# extract ligands from link data frame
unique(active_ligand_target_links_df$ligand_type)

hep_specific_ligands = as.character(unique(active_ligand_target_links_df$ligand[active_ligand_target_links_df$ligand_type=='hep-specific']))
lsec_specific_ligands = as.character(unique(active_ligand_target_links_df$ligand[active_ligand_target_links_df$ligand_type=='lsec-specific']))
stellate_specific_ligands = as.character(unique(active_ligand_target_links_df$ligand[active_ligand_target_links_df$ligand_type=='stellate-specific']))
hep_lsec_ligands = as.character(unique(active_ligand_target_links_df$ligand[active_ligand_target_links_df$ligand_type=='hep-lsec']))
hep_stellate_ligands = as.character(unique(active_ligand_target_links_df$ligand[active_ligand_target_links_df$ligand_type=='hep-stel']))
lsec_stellate_ligands = as.character(unique(active_ligand_target_links_df$ligand[active_ligand_target_links_df$ligand_type=='lsec-stel']))
general_ligands = as.character(unique(active_ligand_target_links_df$ligand[active_ligand_target_links_df$ligand_type=='general-ligands']))

# extract targets from link data frame
unique(active_ligand_target_links_df$target_group)

aj_union_targets = as.character(unique(active_ligand_target_links_df$target[active_ligand_target_links_df$target_group=='aj_union']))
aj_balb_shared_targets = as.character(unique(active_ligand_target_links_df$target[active_ligand_target_links_df$target_group=='aj_balb_shared']))
aj_c57_shared_targets = as.character(unique(active_ligand_target_links_df$target[active_ligand_target_links_df$target_group=='aj_c57_shared']))
balb_union_targets = as.character(unique(active_ligand_target_links_df$target[active_ligand_target_links_df$target_group=='balb_union']))
balb_c57_shared_targets = as.character(unique(active_ligand_target_links_df$target[active_ligand_target_links_df$target_group=='balb_c57_shared']))
c57_union_targets = as.character(unique(active_ligand_target_links_df$target[active_ligand_target_links_df$target_group=='c57_union']))


# order ligands and targets by group - remember that the paste step will keep pasting spaces each time its rerun.
# so if you run it repeatedly you will add more trailing white space and mess up the intersect for ligand_order
target_order = c(aj_union_targets,
                 aj_balb_shared_targets,
                 balb_union_targets,
                 balb_c57_shared_targets,
                 c57_union_targets,
                 aj_c57_shared_targets)

ligand_order = c(hep_specific_ligands,
                 lsec_specific_ligands,
                 stellate_specific_ligands,
                 hep_lsec_ligands,
                 hep_stellate_ligands,
                 lsec_stellate_ligands,
                 general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

#### USER INPUT REQUIRED. 3. Set gap sizes in circos plot, set which ligand/target types are present ####
# define gap sizes below, and comment out ligand and target types not present in the data.

# define segment gaps
width_same_cell_same_ligand_type = 0.5
width_different_cell = 5
width_ligand_target = 0.5
width_same_cell_same_target_type = 0.5
width_same_cell_different_target_type = 0.5

circos_links$ligand_type %>% unique() # view what ligand types are in links
circos_links$target_group %>% unique()

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "hep-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "lsec-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  # rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "stellate-specific") %>% distinct(ligand) %>% nrow() -1)),
  # width_different_cell,
  # rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "hep-lsec") %>% distinct(ligand) %>% nrow() -1)), # no ligands fit this category
  # width_different_cell,
  # rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "hep-stel") %>% distinct(ligand) %>% nrow() -1)), # no ligands fit this category
  # width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "lsec-stel") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "general-ligands") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_group == "aj_union") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_group == "aj_balb_shared") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_group == "balb_union") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_group == "balb_c57_shared") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_group == "c57_union") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_group == "aj_c57_shared") %>% distinct(target) %>% nrow() -1)),
  width_different_cell
)

# make the plot
circos.clear()
circos.par(gap.degree = gaps)
chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = FALSE,
             grid.col = grid_col,
             transparency = transparency,
             diffHeight = 0.005,
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))

# we go back to the first track and customize sector labels
circos.track(track.index = 1,panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)

circos.clear()

# save to a pdf file
transparent_circos_path = paste(out_directory,
                                out_prefix,
                                '_ligand_target_circos_transparent.pdf',
                                sep="")
pdf(transparent_circos_path, width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = FALSE,
             grid.col = grid_col,
             transparency = transparency,
             diffHeight = 0.005,
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()