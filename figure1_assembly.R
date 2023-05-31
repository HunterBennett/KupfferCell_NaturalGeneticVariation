# set user specific working directory
setwd('./strains_github/results/Figure1/')
getwd()

source("fn_clear_session.R")
source("fn_figure1_assembly.R")

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
    "ggfortify",
    "tidytext"
  ),
  suppressPackageStartupMessages(library),
  character.only = TRUE,
))

.lineWidth <- 0.5
lineWidth <- .lineWidth
.pointSize <- 7
smallPointSize <- 6

pca_matrix <- read_tsv("./data_figure1_tpm.txt") %>% rename(gene = 1) %>%
  select(-contains("CB6")) %>%
  column_to_rownames("gene") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

sample_info <- read_tsv("./data_figure1_tpm.txt") %>% select(-1) %>%
  select(-contains("CB6")) %>%
  colnames() %>% tibble() %>% rename(sample = 1) %>%
  separate(
    sample,
    into = c("cell", "strain", "treatment", "replicate"),
    sep =  "_",
    remove = F
  )

sample_pca <- prcomp(pca_matrix)
pc.plot <-
  autoplot(sample_pca, data = sample_info, colour = "strain", shape = "cell") +
  scale_colour_manual(values = c("#C42503", "#455BCD", "#A2FC3C")) +
  # theme_classic() +
  theme(
    text = element_text(size = smallPointSize, colour = "black"),
    plot.title  = element_text(size = smallPointSize, colour = "black"),
    axis.ticks.x = element_line(linewidth = lineWidth * 0.4668623, colour = "black"),
    axis.ticks.y = element_line(linewidth = lineWidth * 0.4668623, colour = "black"),
    axis.title  = element_text(size = smallPointSize, colour = "black"),
    axis.text.x  = element_text(size = smallPointSize, colour = "black"),
    axis.text.y  = element_text(size = smallPointSize, colour = "black"),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.1, 'cm'),
    legend.text = element_text(size = smallPointSize, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    # panel.border = element_blank(),
    axis.line = element_line(linewidth = lineWidth * 0.4668623, colour = "black"),
    plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm")
  ) 

tpm <- read_tsv("./data_figure1_tpm.txt") %>% rename(Gene = 1)
diff <- read_tsv("./data_figure1_deseq2.tsv")
vars1 <- list(
  "AJ vs. BALBcJ" =
    list(
      contrast = c(
        "KC_AJ_Untreated-KC_BALBcJ_Untreated",
        "BMDM_AJ_Untreated-BMDM_BALBcJ_Untreated"
      ),
      group = c("KC", "BMDM"),
      labels = c("KC only", "BMDM only", "BMDM & KC", ""),
      lowerMean = "KC_BALBcJ_Untreated",
      upperMean = "KC_AJ_Untreated",
      upperLabel = "AJ",
      lowerLabel = "BALBcJ"
    ),
  "AJ vs. C57BL6J" =
    list(
      contrast = c(
        "KC_AJ_Untreated-KC_C57BL6J_Untreated",
        "BMDM_AJ_Untreated-BMDM_C57BL6J_Untreated"
      ),
      group = c("KC", "BMDM"),
      labels = c("KC only", "BMDM only", "BMDM & KC", ""),
      lowerMean = "KC_C57BL6J_Untreated",
      upperMean = "KC_AJ_Untreated",
      upperLabel = "AJ",
      lowerLabel = "C57BL6J"
    ),
  "BALBcJ vs. C57BL6J" =
    list(
      contrast = c(
        "KC_BALBcJ_Untreated-KC_C57BL6J_Untreated",
        "BMDM_BALBcJ_Untreated-BMDM_C57BL6J_Untreated"
      ),
      group = c("KC", "BMDM"),
      labels = c("KC only", "BMDM only", "BMDM & KC", ""),
      lowerMean = "KC_C57BL6J_Untreated",
      upperMean = "KC_BALBcJ_Untreated",
      upperLabel = "BALBcJ",
      lowerLabel = "C57BL6J"
    )
)

pa <- 
  vars1 %>% map(
    ., ~ functionMergeProcessMAPlot(variableList = .) +
      xlab("Mean TPM") +
      ylab(paste0("Log2FC(",.x[[6]][[1]], "/", .x[[7]][[1]], ")"))
  )
pa <- wrap_plots(pa[[1]], pa[[2]], pa[[3]]) + plot_layout(guides = "collect")

tpm <- read_tsv("./data_figure1_tpm.txt", show_col_types = FALSE) %>% 
  rename(Gene = 1) %>%
  pivot_longer(
    cols = -1,
    names_to = c("Cell", "Strain", "Treatment", "Replicate"),
    names_sep = "_", values_to = "TPM"
  ) %>%
  group_by(Gene, Cell, Strain, Treatment)

strainList <- list("AJ", "BALBcJ", "C57BL6J")
.diff <- filter(diff, str_detect(contrast, c("KC|CB6F1J"), negate = T))

goiList <-
  list(aj = list("Ccl2", "Ch25h", "Csf1", "Rgs1"),
       balb = list("Abcg1", "Acss2", "Arg2", "Ifi203"),
       c57 = list("Aoah", "C1qb", "Marco", "Sdc1")
  )

.tpm <- tpm %>% filter(Strain != "CB6F1J")
plotList_bmdm <- goiList %>% map(barChartFunction)

goiList <-
  list(aj = list("Atf5", "Cd300e", "Dnase1l3", "H2-DMb2"),
       balb = list("Adam33", "Cxcl14", "Pf4", "Trem2"),
       c57 = list("Cd40", "Cxcl13", "Irak3", "Mefv")
  )

plotList_kc <- goiList %>% map(barChartFunction)

metascape <- read_tsv("./data_figure1_metascape.txt")
pm <-
  metascape |> select(-c(1:3)) |> rename("C57BL6/J" = 2, "BALBc/J" = 3, "A/J" = 4) |>
  pivot_longer(-1, names_to = "strain", values_to = "log10p") |>
  mutate(log10p = abs(log10p)) |>
  filter(log10p > 0) |>
  mutate(Description = str_to_sentence(Description)) |>
  mutate(strain = as.factor(strain),
         Description = as.factor(Description),
         Description = reorder_within(Description, log10p, strain)) |>
  ggplot(aes(x = log10p, y = Description, fill = strain)) +
  # ggplot(aes(reorder_within(x = `-log10(p)`, y = Description, fill = strain), `-log10(p)`)) +
  geom_col() +
  # coord_flip() + 
  # scale_x_reordered() +
  # scale_x_reverse() +
  facet_wrap(~ strain, scales = "free", drop = T) +
  scale_y_reordered(labels = function(x) str_wrap(x, width = 41)) +
  
  scale_fill_manual(values = c("#C42503", "#455BCD", "#A2FC3C")) +
  scale_colour_manual(values = c("#C42503", "#455BCD", "#A2FC3C")) +
  ylab(NULL) + xlab("-log10(p)") +
  theme(
    text = element_text(size = smallPointSize, colour = "black"),
    line = element_line(size = lineWidth * 0.4668623, colour = "black"),
    rect = element_rect(size = lineWidth * 0.4668623),
    axis.ticks.x = element_line(size = lineWidth * 0.4668623, colour = "black"),
    axis.ticks.y = element_line(size = lineWidth * 0.4668623, colour = "black"),
    axis.text.x  = element_text(size = smallPointSize, colour = "black"),
    axis.text.y  = element_text(size = smallPointSize, colour = "black"),
    axis.title = element_text(size = smallPointSize, colour = "black"),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "none",
    legend.title = element_text(size = smallPointSize, colour = "black"),
    legend.key = element_blank(),
    legend.key.size = unit(0.1, 'cm'),
    legend.background = element_blank(),
    legend.text = element_text(size = smallPointSize, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = lineWidth * 0.4668623, colour = "black"),
    plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm"),
    panel.spacing = unit(0.1, "cm"),
    strip.text = element_text(
      margin = margin(0,0,0,0),
      face = "plain",
      size = smallPointSize
    ),
    strip.placement = "outside",
    strip.clip = "off",
    strip.background = element_blank()
  )
pm

lvl1 <- wrap_plots(pc.plot, pa, widths = c(0.25, 0.75))
lvl2 <- wrap_plots(plotList_bmdm, guides = "collect")
lvl3 <- wrap_plots(plotList_kc, guides = "collect")
lvl2 <- wrap_plots(wrap_plots(plotList_bmdm), wrap_plots(plotList_kc), nrow = 2, guides = "collect", tag_level = "new")
lvl4 <- wrap_plots(pm, nrow = 1)

fig1 <-
  plot_grid(
  lvl1,
  lvl2,
  lvl4,
  nrow = 3,
  rel_heights = c(37.5, 65, 45), #37.5
  labels = "auto",
  label_size = .pointSize,
  label_fontface = "plain",
  axis = "lr"
)

ggsave(plot = fig1, filename = "figure1.pdf", height = 152.5, width = 180, units = "mm")
