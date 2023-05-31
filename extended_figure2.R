# set user specific  working directory
setwd('./strains_github/results/ExtendedData_Figure2/')
getwd()


source("fn_clear_session.R")
source("pheatmap2.R")
invisible(lapply(
  c(
    "tidyverse",
    "ggrepel",
    "cowplot",
    "scales",
    "ggrastr",
    "grid",
    "gridExtra",
    "viridis",
    "knitr",
    "skimr",
    "ggplotify",
    "RColorBrewer",
    "patchwork",
    "ggpubr"
  ),
  suppressPackageStartupMessages(library),
  character.only = TRUE,
))

tpm <- read_tsv("./data_figure1_tpm.txt") %>% rename(Gene = 1)
diff <- read_tsv("./data_figure1_deseq2.tsv")

################################################################################
# composite heatmaps

goi_all <- diff %>%
  filter(str_detect(contrast, "BMDM.+.BMDM.+|KC.+.KC.+")) %>%
  filter(str_detect(contrast, "CB6F1J", negate = T)) %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>%
  pull(Gene) %>% unique()

all_data <-
  tpm %>%
  filter(Gene %in% goi_all) %>%
  select(-contains("CB6F1J")) %>%
  filter(if_any(-1, ~ .x > 8)) %>%
  pivot_longer(-Gene) %>%
  mutate(value = log2(value + 1)) %>%
  pivot_wider(names_from = name, values_from = value)

all_annotation <-
  data.frame(CellType = factor(
    c(
      "BMDM",
      "BMDM",
      "BMDM",
      "BMDM",
      "BMDM",
      "BMDM",
      "KC",
      "KC",
      "KC",
      "KC",
      "KC",
      "KC"
    )
  ),
  Strain = factor(
    c(
      "AJ",
      "AJ",
      "BALBcJ",
      "BALBcJ",
      "C57BL6J",
      "C57BL6J",
      "AJ",
      "AJ",
      "BALBcJ",
      "BALBcJ",
      "C57BL6J",
      "C57BL6J"
    )
  ))
rownames(all_annotation) <- colnames(all_data[-1])

gene_list_all <-
  all_data %>%
  column_to_rownames(var = "Gene") %>%
  pheatmap2(
    .,
    scale = "row",
    clustering_method = "ward.D",
    cutree_rows = 6
  )

gene_list_all$clusterGene |> as.data.frame() |> tibble() |> write_csv("./ed_fig2_gene_list_all.csv")

hm_all <-
  all_data %>%
  column_to_rownames(var = "Gene") %>%
  pheatmap3(
    .,
    scale = "row",
    clustering_method = "ward.D",
    cutree_rows = 6,
    annotation_names_row = F,
    treeheight_row = 5,
    treeheight_col = 5,
    fontsize = 6,
    fontsize_row = 5,
    show_colnames = F,
    show_rownames = F,
    color = colorRampPalette(rev(brewer.pal(
      n = 11, name =
        "RdBu"
    )))(17),
    legend = T,
    annotation_legend = T,
    annotation_col = all_annotation,
    annotation_colors =
      list(
        CellType = c(BMDM = "#261433", KC = "#F7C5A5"),
        Strain = c(
          AJ = "#C42503",
          BALBcJ = "#455BCD",
          C57BL6J = "#A2FC3C"
        ),
        Clusters = c(
          C1 = "#26172A",
          C2 = "#3F366F",
          C3 = "#38639D",
          C4 = "#3491A8",
          C5 = "#47BFAD",
          C6 = "#A8E1BC"
        )
      )
  ) %>% as.ggplot()

meta_all <- read_csv("./metascape_all_HeatmapSelectedGO.csv")

ann_plot_meta_all <- data.frame(
  Clusters = factor(c("C1", "C2", "C3", "C4", "C5", "C6"))
)
row.names(ann_plot_meta_all) <- ann_plot_meta_all$Clusters

hm_meta_all <- meta_all |> pivot_longer(cols = -c(1:2)) |> mutate(name = str_remove_all(name, "_LogP_")) |>
  unite("term", c("GO", "Description"), sep = "; ") |>
  mutate(term = str_wrap(term, width = 60, indent = 0, exdent = 0, whitespace_only = TRUE)) |>
  pivot_wider(names_from = name, values_from = value) |>
  column_to_rownames(var = "term") %>%
  pheatmap(
    .,
    scale = "none",
    fontsize = 6,
    fontsize_row = 5,
    treeheight_row = 5,
    treeheight_col = 5,
    show_colnames = F,
    color = viridis_pal(option = "rocket")(17),
    annotation_col = ann_plot_meta_all,
    annotation_colors =
      list(
        Clusters = c(
          C1 = "#26172A",
          C2 = "#3F366F",
          C3 = "#38639D",
          C4 = "#3491A8",
          C5 = "#47BFAD",
          C6 = "#A8E1BC"
        )
      )
  ) %>%
  as.ggplot()

################################################################################
# BMDM heatmaps

goi_bmdm <- diff %>%
  filter(str_detect(contrast, "BMDM.+.BMDM.+")) %>%
  filter(str_detect(contrast, "CB6F1J", negate = T)) %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>%
  pull(Gene) %>% unique()

bmdm_data <-
  tpm %>% 
  filter(Gene %in% goi_bmdm) %>%
  # slice_sample(n = 100) %>%
  select(-contains("CB6F1J")) %>%
  select(-contains("KC")) %>%
  filter(if_any(-1, ~ .x > 8)) %>%
  pivot_longer(-Gene) %>%
  mutate(value = log2(value + 1)) %>%
  pivot_wider(names_from = name, values_from = value)

bmdm_annotation <-
  data.frame(CellType = factor(
    c(
      "BMDM",
      "BMDM",
      "BMDM",
      "BMDM",
      "BMDM",
      "BMDM"
    )
  ),
  Strain = factor(
    c(
      "AJ",
      "AJ",
      "BALBcJ",
      "BALBcJ",
      "C57BL6J",
      "C57BL6J"
    )
  ))

rownames(bmdm_annotation) <- colnames(bmdm_data[-1])

gene_list_bmdm <-
  bmdm_data %>% 
  column_to_rownames(var = "Gene") %>%
  pheatmap2(
    .,
    scale = "row",
    clustering_method = "ward.D",
    cutree_rows = 8
  )

gene_list_bmdm$clusterGene |> as.data.frame() |> tibble() |> write_csv("./ed_fig2_gene_list_bmdm.csv")

viridis_pal(option = "mako")(8) 

hm_bmdm <-
  bmdm_data %>%
  column_to_rownames(var = "Gene") %>%
  pheatmap3(
    .,
    scale = "row",
    clustering_method = "ward.D",
    cutree_rows = 8,
    annotation_names_row = F,
    treeheight_row = 5,
    treeheight_col = 5,
    fontsize = 6,
    fontsize_row = 5,
    show_colnames = F,
    show_rownames = F,
    color = colorRampPalette(rev(brewer.pal(
      n = 11, name =
        "RdBu"
    )))(17),
    legend = T,
    annotation_legend = T,
    annotation_col = all_annotation,
    annotation_colors =
      list(
        CellType = c(BMDM = "#261433", KC = "#F7C5A5"),
        Strain = c(
          AJ = "#C42503",
          BALBcJ = "#455BCD",
          C57BL6J = "#A2FC3C"
        ),
        Clusters = c(
          C1 = "#0B0405FF",
          C2 = "#2E1E3CFF",
          C3 = "#413D7BFF",
          C4 = "#37659EFF",
          C5 = "#348FA7FF",
          C6 = "#40B7ADFF",
          C7 = "#8AD9B1FF",
          C8 = "#DEF5E5FF"
        )
      )
  ) %>% as.ggplot()

meta_bmdm <- read_csv("./metascape_bmdm_HeatmapSelectedGOParent.csv")

ann_plot_meta_bmdm <- data.frame(
  Clusters = factor(c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"))
)
row.names(ann_plot_meta_bmdm) <- ann_plot_meta_bmdm$Clusters

hm_meta_bmdm <- meta_bmdm |> 
  mutate(Description = str_replace_all(Description, "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
                                  "adaptive imm. response based on somatic rec. of imm. receptors built  Ig domains")) |>
  pivot_longer(cols = -c(1:2)) |> mutate(name = str_remove_all(name, "_LogP_")) |>
  unite("term", c("GO", "Description"), sep = "; ") |>
  mutate(term = str_wrap(term, width = 60, indent = 0, exdent = 0, whitespace_only = TRUE)) |>
  pivot_wider(names_from = name, values_from = value) |>
  column_to_rownames(var = "term") %>%
  pheatmap(
    .,
    scale = "none",
    fontsize = 6,
    fontsize_row = 5,
    treeheight_row = 5,
    treeheight_col = 5,
    show_colnames = F,
    color = viridis_pal(option = "rocket")(17),
    annotation_col = ann_plot_meta_bmdm,
    annotation_colors =
      list(
        Clusters = c(
          C1 = "#0B0405FF",
          C2 = "#2E1E3CFF",
          C3 = "#413D7BFF",
          C4 = "#37659EFF",
          C5 = "#348FA7FF",
          C6 = "#40B7ADFF",
          C7 = "#8AD9B1FF",
          C8 = "#DEF5E5FF"
        )
      )
  ) %>%
  as.ggplot()

################################################################################
# KC heatmaps
goi_kc <- diff %>%
  filter(str_detect(contrast, "KC.+.KC.+")) %>%
  filter(str_detect(contrast, "CB6F1J", negate = T)) %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>%
  pull(Gene) %>% unique()

kc_data <-
  tpm %>% 
  filter(Gene %in% goi_kc) %>%
  select(-contains("CB6F1J")) %>%
  select(-contains("BMDM")) %>%
  filter(if_any(-1, ~ .x > 8)) %>%
  pivot_longer(-Gene) %>%
  mutate(value = log2(value + 1)) %>%
  pivot_wider(names_from = name, values_from = value)

kc_annotation <-
  data.frame(CellType = factor(
    c(
      "KC",
      "KC",
      "KC",
      "KC",
      "KC",
      "KC"
    )
  ),
  Strain = factor(
    c(
      "AJ",
      "AJ",
      "BALBcJ",
      "BALBcJ",
      "C57BL6J",
      "C57BL6J"
    )
  ))
rownames(kc_annotation) <- colnames(kc_data[-1])

gene_list_kc <-
  kc_data %>% 
  column_to_rownames(var = "Gene") %>%
  pheatmap2(
    .,
    scale = "row",
    clustering_method = "ward.D",
    cutree_rows = 7
  )

gene_list_kc$clusterGene |> as.data.frame() |> tibble() |> write_csv("./ed_fig2_gene_list_kc.csv")

viridis_pal(option = "mako")(7) 

hm_kc <-
  kc_data %>%
  column_to_rownames(var = "Gene") %>%
  pheatmap3(
    .,
    scale = "row",
    clustering_method = "ward.D",
    cutree_rows = 7,
    annotation_names_row = F,
    treeheight_row = 5,
    treeheight_col = 5,
    fontsize = 6,
    fontsize_row = 5,
    show_colnames = F,
    show_rownames = F,
    color = colorRampPalette(rev(brewer.pal(
      n = 11, name =
        "RdBu"
    )))(17),
    legend = T,
    annotation_legend = T,
    annotation_col = all_annotation,
    annotation_colors =
      list(
        CellType = c(BMDM = "#261433", KC = "#F7C5A5"),
        Strain = c(
          AJ = "#C42503",
          BALBcJ = "#455BCD",
          C57BL6J = "#A2FC3C"
        ),
        Clusters = c(
          C1 = "#0B0405FF",
          C2 = "#342346FF",
          C3 = "#40498EFF",
          C4 = "#357BA2FF",
          C5 = "#38AAACFF",
          C6 = "#78D6AEFF",
          C7 = "#DEF5E5FF"
        )
      )
  ) %>% as.ggplot()

meta_kc <- read_csv("./metascape_kc_HeatmapSelectedGO.csv")

ann_plot_meta_kc <- data.frame(
  Clusters = factor(c("C1", "C2", "C3", "C4", "C5", "C6", "C7"))
)
row.names(ann_plot_meta_kc) <- ann_plot_meta_kc$Clusters

hm_meta_kc <- meta_kc |> pivot_longer(cols = -c(1:2)) |> mutate(name = str_remove_all(name, "_LogP_")) |>
  unite("term", c("GO", "Description"), sep = "; ") |>
  mutate(term = str_wrap(term, width = 60, indent = 0, exdent = 0, whitespace_only = TRUE)) |>
  pivot_wider(names_from = name, values_from = value) |>
  column_to_rownames(var = "term") %>%
  pheatmap(
    .,
    scale = "none",
    fontsize = 6,
    fontsize_row = 5,
    treeheight_row = 5,
    treeheight_col = 5,
    show_colnames = F,
    color = viridis_pal(option = "rocket")(17),
    annotation_col = ann_plot_meta_kc,
    annotation_colors =
      list(
        Clusters = c(
          C1 = "#0B0405FF",
          C2 = "#342346FF",
          C3 = "#40498EFF",
          C4 = "#357BA2FF",
          C5 = "#38AAACFF",
          C6 = "#78D6AEFF",
          C7 = "#DEF5E5FF"
        )
      )
  ) %>%
  as.ggplot()

################################################################################
# Figure assembly
design <- "
112222
334444
556666
"
group1 <- hm_all + hm_meta_all + hm_bmdm + hm_meta_bmdm + hm_kc + hm_meta_kc +
  plot_layout(design = design, guides = 'collect') +
  plot_annotation(tag_levels = c("a")) & 
  theme(plot.tag = element_text(size = 7, color = "black", face = "bold"))
group1
ggsave(filename = "./extended_figure2.pdf", height = 10, width = 7.5, units = "in")

df_list <- list(all_data, bmdm_data, kc_data) |> set_names(c("BMDM & KC", "BMDM alone", "KC alone"))

tmp <- map_dfr(.x = df_list, ~ pull(., Gene) |> length()) |> 
  pivot_longer(1:3) |> rename(" " = 1, "Total DEG #" = 2)

tb1 <- ggtexttable(
  tmp,
  rows = NULL,
  theme = ttheme(
    base_style = "light",
    colnames.style = colnames_style(
      color = "black",
      face = "plain",
      size = 6,
      fill = NA
    ),
    tbody.style = tbody_style(
      color = "black",
      face = "plain",
      size = 6,
      linewidth = 0.25,
      linecolor = "gray",
      fill = NA
    )
  )
)
tb1

df_list2 <- list(bmdm_data, kc_data)
map(.x = df_list2, ~ pull(., Gene)) |> flatten() |> length()
map(.x = df_list2, ~ pull(., Gene)) |> reduce(union) |> length()
map(.x = df_list2, ~ pull(., Gene)) |> reduce(symdiff) |> length()
map(.x = df_list2, ~ pull(., Gene)) |> reduce(intersect) |> length()
