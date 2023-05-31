# Hunter Bennett
# Glass Lab | Kupffer strains project | April 11 2023
# BALB/cJ Kupffer cell analysis using nichenet
# Notebook to run NicheNet regression on BALB/cJ specific genes (union set)

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
geneset = "./nichenet_balbUnion.tsv" # pick set of target genes
tpm = 10 # TPM to define an 'expressed' gene
out_prefix = 'balb_kupffer'
out_directory = './BALBcJ/'
dir.create(out_directory)
expression_matrix_path = './all_liver_cells_rna_tpm.txt'
n_links = 500 # number of ligand-target links to extract for analysis

# load ligand-target model developed by NicheNet authors
# Weighted matrix of the interaction affinity between ligands and TARGET GENES (skips the receptor)
# need to re-read the nichenet paper to figure out where this comes from
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5, 1:5]

# load expression data matrix with gene names
# in rows and the samples in columns
expression = read.csv(expression_matrix_path,sep='\t')

# rename first column so can be easily accessed
names(expression)[names(expression)=="X"] <- "gene"
expression$gene <- toupper(expression$gene)
expression[1:5, 1:5]
names(expression)

# extract cell columns
aj_hep = grep('aj_hepatocyte', names(expression))
balb_hep = grep('balbc_hepatocyte', names(expression))
c57_hep = grep('c57_hepatocyte', names(expression))

aj_lsec = grep('aj_lsec', names(expression))
balb_lsec = grep('balbc_lsec', names(expression))
c57_lsec = grep('c57_lsec', names(expression))

aj_stel = grep('aj_stel', names(expression))
balb_stel = grep('balbc_stel', names(expression))
c57_stel = grep('c57_stel', names(expression))

aj_kup = grep('aj_kup', names(expression))
balb_kup = grep('balbc_kup', names(expression))
c57_kup = grep('c57_kup', names(expression))

# extract cell gene expression values
# this requires setting an arbitrary value to the levels required for a gene to be
# "expressed". Also need to convert string to all upper case for string matching
# since all the nichenet rows and columns are all upper case

# extract hepatocyte genes
expressed_genes_aj_hep = expression[rowMeans(expression[, aj_hep])>tpm, "gene"]
expressed_genes_balb_hep = expression[rowMeans(expression[, balb_hep])>tpm, "gene"]
expressed_genes_c57_hep = expression[rowMeans(expression[, c57_hep])>tpm, "gene"]
expressed_genes_hep = union(union(expressed_genes_aj_hep, expressed_genes_balb_hep), expressed_genes_c57_hep)

# extract lsec genes
expressed_genes_aj_lsec = expression[rowMeans(expression[, aj_lsec])>tpm, "gene"]
expressed_genes_balb_lsec = expression[rowMeans(expression[, balb_lsec])>tpm, "gene"]
expressed_genes_c57_lsec = expression[rowMeans(expression[, c57_lsec])>tpm, "gene"]
expressed_genes_lsec = union(union(expressed_genes_aj_lsec, expressed_genes_balb_lsec), expressed_genes_c57_lsec)

# extract stellate genes
expressed_genes_aj_stel = expression[rowMeans(expression[, aj_stel])>tpm, "gene"]
expressed_genes_balb_stel = expression[rowMeans(expression[, balb_stel])>tpm, "gene"]
expressed_genes_c57_stel = expression[rowMeans(expression[, c57_stel])>tpm, "gene"]
expressed_genes_stel = union(union(expressed_genes_aj_stel, expressed_genes_balb_stel), expressed_genes_c57_stel)

# extract kupffer genes
expressed_genes_aj_kup = expression[rowMeans(expression[, aj_kup])>tpm, "gene"]
expressed_genes_balb_kup = expression[rowMeans(expression[, balb_kup])>tpm, "gene"]
expressed_genes_c57_kup = expression[rowMeans(expression[, c57_kup])>tpm, "gene"]
expressed_genes_kup = union(union(expressed_genes_aj_kup, expressed_genes_balb_kup), expressed_genes_c57_kup)

# extract all niche genes. we now have lists of expressed genes for each cell type of interest
expressed_genes_aj_niche = union(union(expressed_genes_aj_hep, expressed_genes_aj_lsec), expressed_genes_aj_stel)
expressed_genes_balb_niche = union(union(expressed_genes_balb_hep, expressed_genes_balb_lsec), expressed_genes_balb_stel)
expressed_genes_c57_niche = union(union(expressed_genes_c57_hep, expressed_genes_c57_lsec), expressed_genes_c57_stel)

# read in target genes
spec = read.table(geneset, sep='\t', header=TRUE)
head(spec)

names(spec)[1] <- 'gene'
spec$gene <- toupper(spec$gene)
spec[1:5, 'gene']


#### USER INPUT REQUIRED 2. Set Target Gene Set and Niche Gene Set ####
# foreground set is the set of target genes whose expression we are
# trying to predict using ligands in nichenet.
# niche set is the set of potential ligands in the niche driving
# that cells gene expression profile
foreground_set = expressed_genes_balb_kup
niche_set = expressed_genes_balb_niche

# filter based on tpm cutoff
geneset_oi <- spec[spec$gene %in% foreground_set, "gene"]
print(length(geneset_oi))
length(spec$gene)

# filter genes in ligand/target matrix
# this is the 'foreground' gene set of genes specific
geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
geneset_oi <- as.vector(geneset_oi)

# background is all KC expressed genes that are contained in the nichenet ligand target matrix
background_expressed_genes = foreground_set[foreground_set %in% rownames(ligand_target_matrix)]

# check the size of the vectors - background should be substantially larger than foreground
print(length(geneset_oi))
print(length(background_expressed_genes))

# read in niche-net's defined interaction matrix
# which records the affinity of interaction between ligands "from" column
# with their cognate receptors "to" column
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

# If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions
# and only keep ligand-receptor interactions that are described in curated databases
# To do this: uncomment following line of code:
lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

# find ligands expressed in target cell of interest, in this case, bulk 130d organoid
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,niche_set)
# optional: append receptors - Leptin
expressed_ligands <- append(expressed_ligands, c("LEP", "ADIPOQ", "GHRL", "FGF19", "INS", "GCG"))

# find receptors expressed in target cell of interest, in this case, organoid microglia
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,foreground_set)

# trim network so that it only includes ligand/receptor pairs that are contained in our dataset
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

# extract possible ligands (since each ligand can activate multiple receptors this is a redundant list that needs to be cleaned)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)


# Perform Niche Net ligand analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                              background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix,
                                              potential_ligands = potential_ligands)


# ligand_activities is a table with the pearson correlation stats for each ligand
# as well as the auroc and aupr predictive stats.
ligand_activities %>% arrange(-pearson) 

# save ligand activities
ligand_activities_string = paste(out_directory,
                                out_prefix,
                                '_circos_ligand_activity.txt',
                                sep="")
ligand_activities %>% arrange(-pearson) %>% write_tsv(file=ligand_activities_string)

#### USER INPUT REQUIRED. 3. Select number of ligands to pass to circos for plotting and change strain for rowmeans ####
# update: want to display many ligands since the output matrices are useful for examining
# high scoring ligands, in this case we will try and export all ligands with correlation 
# above 0.04
n_ligands = 10 # ligands for circos
n_ligands_df = 100 # ligands for df

# write out top ligand-targets to df
best_upstream_ligands_extended = ligand_activities %>% top_n(n_ligands_df, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

# Infer target genes of top ranked ligands
active_ligand_target_links_df_towrite = best_upstream_ligands_extended %>% lapply(get_weighted_ligand_target_links,
                                                                                  geneset = geneset_oi,
                                                                                  ligand_target_matrix = ligand_target_matrix,
                                                                                  n = n_links) %>% bind_rows()

# write to df
extended_ligand_target_path = paste(out_directory,
                                    out_prefix,
                                    '_circos_ligand_target_links_df_extended.txt',
                                    sep="")

write.table(active_ligand_target_links_df_towrite,
            file=extended_ligand_target_path,
            sep='\t',
            row.names=FALSE)

# write out top ligand-receptors to DF

# get the ligand-receptor network of the top-ranked ligands
lr_network_top_extended = lr_network %>% filter(from %in% best_upstream_ligands_extended & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors_extended = lr_network_top_extended %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df_extended = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands_extended & to %in% best_upstream_receptors_extended)

# write to df
extended_ligand_receptor_path = paste(out_directory,
                                      out_prefix,
                                      '_circos_ligand_receptor_links_df_extended.txt',
                                      sep="")

write.table(lr_network_top_df_extended,
            file=extended_ligand_receptor_path,
            sep='\t',
            row.names=FALSE)


# extract bets upstream ligands
best_upstream_ligands = ligand_activities %>% top_n(n_ligands, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

# fix mg_expression and bulk_expresssion data frames
rownames(expression) <- expression$gene
# expression <- expression[, -1]

# create tibble with average expression values
ligand_expression_tbl = tibble(
  ligand = best_upstream_ligands, 
  hep = rowMeans(expression[best_upstream_ligands, balb_hep]),
  lsec = rowMeans(expression[best_upstream_ligands, balb_lsec]),
  stel = rowMeans(expression[best_upstream_ligands, balb_stel]))

#### USER INPUT REQUIRED. 4. manually curate strain specificity ####
# expression tibble to decide which ligands fit which category
ligand_expression_tbl

hep_spec_ligands = c('FGF1', 'APOE', 'VEGFA', 'VTN')
lsec_spec_ligands = c('CXCL9')
stel_spec_ligands = c('CXCL12', 'SPP1')
hep_lsec_ligands = c()
hep_stel_ligands = c()
lsec_stel_ligands = c('TGFB1', 'ITGB1')
general_ligands = c('LEP')

# Save version without specifying which strain they came from.
ligand_type_indication_df = tibble(
  ligand_type = c(rep("hep-specific", times = hep_spec_ligands %>% length()),
                  rep("lsec-specific", times = lsec_spec_ligands %>% length()),
                  rep("stellate-specific", times = stel_spec_ligands %>% length()),
                  rep("hep-lsec", times = hep_lsec_ligands %>% length()),
                  rep("hep-stel", times = hep_stel_ligands %>% length()),
                  rep("lsec-stel", times = lsec_stel_ligands %>% length()),
                  rep("general-ligands", times = general_ligands %>% length())),
  ligand = c(hep_spec_ligands, 
             lsec_spec_ligands,
             stel_spec_ligands,
             hep_lsec_ligands,
             hep_stel_ligands,
             lsec_stel_ligands,
             general_ligands))

ligand_type_indication_df

# Infer target genes of top ranked ligands
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
                                                                 geneset = geneset_oi,
                                                                 ligand_target_matrix = ligand_target_matrix,
                                                                 n = n_links) %>% bind_rows()

# if you want to make circos plots for multiple gene sets
# combine the different data frames and differentiate which target belongs to which gene set via the target type
active_ligand_target_links_df = active_ligand_target_links_df %>% inner_join(ligand_type_indication_df) 
active_ligand_target_links_df

# export here for circos plots
ligand_target_df_string = paste(out_directory,
                                out_prefix,
                                '_circos_ligand_target_links_df.txt',
                                sep="")

write.table(active_ligand_target_links_df,
            file=ligand_target_df_string,
            sep='\t',
            row.names=FALSE)

# Save version with strain specificity for combining circos plots
ligand_type_indication_df = tibble(
  ligand_type = c(rep("hep-specific", times = hep_spec_ligands %>% length()),
                  rep("lsec-specific", times = lsec_spec_ligands %>% length()),
                  rep("stellate-specific", times = stel_spec_ligands %>% length()),
                  rep("hep-lsec", times = hep_lsec_ligands %>% length()),
                  rep("hep-stel", times = hep_stel_ligands %>% length()),
                  rep("lsec-stel", times = lsec_stel_ligands %>% length()),
                  rep("general-ligands", times = general_ligands %>% length())),
  ligand = c(hep_spec_ligands, 
             lsec_spec_ligands,
             stel_spec_ligands,
             hep_lsec_ligands,
             hep_stel_ligands,
             lsec_stel_ligands,
             general_ligands))

# Infer target genes of top ranked ligands
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
                                                                 geneset = geneset_oi,
                                                                 ligand_target_matrix = ligand_target_matrix,
                                                                 n = n_links) %>% bind_rows()
# if you want to make circos plots for multiple gene sets
# combine the different data frames and differentiate which target belongs to which gene set via the target type
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = out_prefix) %>% inner_join(ligand_type_indication_df) 

# export here for combined circos plots
target_spec_ligand_target_df_string = paste(out_directory,
                                            out_prefix,
                                            "_target_spec_circos_ligand_target_links_df.txt",
                                            sep="")
write.table(active_ligand_target_links_df,
            file=target_spec_ligand_target_df_string,
            sep='\t',
            row.names=FALSE)



