# set user specific working directory
setwd('./strains_github/results/Figure6/')
getwd()

# Load functions, libraries, and setup environment
source("fn_clear_session.R")
source("fn_process_plot.R")
source("fn_process_peaks_strict.R")

lapply(X = c("tidyverse", "ggrepel", "cowplot", "grid", 'ggrastr', 'knitr',
             "scales", "gridExtra", "ggpubr", "patchwork", "ggridges"), 
       FUN = library, character.only = TRUE)

lineWidth <- 0.5
pointSize <- 0.25
textSize <- 7
smallTextSize <- 5

# Load data
tpm <- read_tsv("./data_figure6_tpm_means.tsv")
diffRNA <- read_tsv("./data_figure6_deseq2.tsv")

diffPeaks <- file.path(
  "./data_atac_f0_f1/",
  list.files(path = "./data_atac_f0_f1/", pattern = ".txt")
) %>% map_df( ~ read_plus(.)) %>% rename(Position = 1)

# analysis and plot for panel 6a
vars <-
  list(
    contrast = c(
      "KC_BALBcJ_F0_LPS-KC_BALBcJ_F0_Untreated",
      "KC_C57BL6J_F0_LPS-KC_C57BL6J_F0_Untreated"
    ),
    group = c(" BALBcJ", " C57BL6J"),
    labels = c("BALBcJ", "C57BL6J", "Shared", "Same"),
    labelKey = c("First", "Second", "Both", "NotDiff")
  )

fig6a <- mergeProcessRatioFunction() + ggtitle("RNA-seq") + 
  xlab("Log2(LPS/untreated)\nC57BL6J") + ylab("BALBcJ\nLog2(LPS/untreated)")

# analysis for panel 6b
diffPeaks %>% select(contrast) %>% unique()

vars <-
  list(
    contrast = c(
      "BALBCJ_F0_LPS-BALBCJ_F0_Control",
      "C57BL6J_F0_LPS-C57BL6J_F0_Control"
    ),
    group = c(".balb", ".c57"),
    labels = c("BALB", "C57", "Both", "Same"),
    labelKey = c("First", "Second", "Both", "NotDiff")
  )

counts <- inner_join(
  diffPeaks %>% filter(contrast == "BALBCJ_F0_LPS-BALBCJ_F0_Control"),
  diffPeaks %>% filter(contrast == "C57BL6J_F0_LPS-C57BL6J_F0_Control"), by = "Position",
  suffix = c(vars$group[1], vars$group[2])) %>%
  select(contains("Position") | contains("normTag") & -contains("F1")) %>%
  select(where(~ !all(is.na(.x))))

merged <- functionProcessPeaks(full_join, diffPeaks) %>%
  inner_join(., counts) %>% filter(if_any(.cols = 7:10, ~ .x > 4))

fig6b <-
  functionRatioPlot(merged,
                    textSize = textSize,
                    smallTextSize = smallTextSize,
                    lineWidth = lineWidth,
                    pointSize = pointSize) +
  ggtitle("ATAC-seq") + xlab("Log2(LPS/untreated)\nC57BL6J") + ylab("BALBcJ\nLog2(LPS/untreated)")

# analysis for panel 6e
coordinates <- read_tsv("./data_atac_f0_f1_normalized.txt") %>% rename("Position" = 1)

vars <-
  list(
    contrast = c(
      "BALBCJ_F0_LPS-BALBCJ_F0_Control",
      "C57BL6J_F0_LPS-C57BL6J_F0_Control"
    ),
    group = c(".balb", ".c57"),
    labels = c("BALB", "C57", "Both", "Same"),
    labelKey = c("First", "Second", "Both", "NotDiff")
  )

counts <- inner_join(
  diffPeaks %>% filter(contrast == "BALBCJ_F0_LPS-BALBCJ_F0_Control"),
  diffPeaks %>% filter(contrast == "C57BL6J_F0_LPS-C57BL6J_F0_Control"), by = "Position",
  suffix = c(vars$group[1], vars$group[2])) %>%
  select(contains("Position") | contains("normTag") & -contains("F1")) %>%
  select(where(~ !all(is.na(.x))))

merged <- functionProcessPeaks(full_join, diffPeaks) %>%
  inner_join(., counts) %>% filter(if_any(.cols = 7:10, ~ .x > 4))

filter(diffPeaks, contrast == "BALBCJ_F0_Control-C57BL6J_F0_Control") %>% colnames()

.list <-
  list(
    f0Peaks_balb_LPS_up =
      merged %>%
      filter(DEG == "BALB" & log2FoldChange.balb > 1) %>%
      inner_join(coordinates, .),
    f0Peaks_balb_LPS_down =
      merged %>%
      filter(DEG == "BALB" & log2FoldChange.balb < -1) %>%
      inner_join(coordinates, .),
    f0Peaks_c57_LPS_up =
      merged %>%
      filter(DEG == "C57" & log2FoldChange.c57 > 1) %>%
      inner_join(coordinates, .),
    f0Peaks_c57_LPS_down =
      merged %>%
      filter(DEG == "C57" & log2FoldChange.c57 < -1) %>%
      inner_join(coordinates, .),
    f0Peaks_f1_same =
      merged %>%
      filter(DEG == "Same") %>%
      inner_join(coordinates, .)
  )

c57Prepped <- .list %>% pluck("f0Peaks_c57_LPS_up") %>% select(1) %>%
  inner_join(.,
             filter(diffPeaks, contrast == "C57BL6J_F0_Control-BALBCJ_F0_Control"),
             by = "Position") %>%
  inner_join(
    .,
    filter(diffPeaks, contrast == "C57BL6J_F0_LPS-BALBCJ_F0_LPS"),
    by = "Position",
    suffix = c(".control", ".lps")
  ) %>%
  mutate("Basal accessibility" = factor(
    ifelse(
      log2FoldChange.control > 1 & padj.control < 0.05,
      "Low basal",
      ifelse(log2FoldChange.control < -1 &
               padj.control < 0.05, "High basal", "Equal basal")
    )
  )) %>% 
  select(-contains("F1")) %>% 
  select(1, 9, 12, 22, 25, 26) %>% drop_na() %>%
  rename("balb_ctrl" = 2, "c57_ctrl" = 3, "balb_lps"= 4, "c57_lps" = 5) %>%
  pivot_longer(2:5, names_to = "Treatment", values_to = "Counts") %>%
  separate(Treatment, into = c("Strain", "Treatment")) %>% 
  mutate(Strain = str_replace_all(Strain, c("balb" =  "BALBcJ", "c57" = "C57BL6J"))) %>%
  mutate(Treatment = str_replace_all(Treatment, c(ctrl = "PBS", lps = "LPS")))  %>%
  mutate(Treatment = fct_relevel(Treatment, c("PBS", "LPS"))) %>%
  mutate(`Basal accessibility` = fct_relevel(`Basal accessibility`, c("Low basal", "Equal basal", "High basal")))

balbPrepped <- .list %>% pluck("f0Peaks_balb_LPS_up") %>% select(1) %>%
  inner_join(.,
             filter(diffPeaks, contrast == "BALBCJ_F0_Control-C57BL6J_F0_Control"),
             by = "Position") %>%
  inner_join(
    .,
    filter(diffPeaks, contrast == "BALBCJ_F0_LPS-C57BL6J_F0_LPS"),
    by = "Position",
    suffix = c(".control", ".lps")
  ) %>%
  mutate("Basal accessibility" = factor(
    ifelse(
      log2FoldChange.control > 1 & padj.control < 0.05,
      "Low basal",
      ifelse(log2FoldChange.control < -1 &
               padj.control < 0.05, "High basal", "Equal basal")
    )
  )) %>%
  select(-contains("F1")) %>% 
  select(1, 9, 12, 22, 25, 26) %>% drop_na() %>%
  rename("balb_ctrl" = 2, "c57_ctrl" = 3, "balb_lps"= 4, "c57_lps" = 5) %>%
  pivot_longer(2:5, names_to = "Treatment", values_to = "Counts") %>%
  separate(Treatment, into = c("Strain", "Treatment")) %>% 
  mutate(Strain = str_replace_all(Strain, c("balb" =  "BALBcJ", "c57" = "C57BL6J"))) %>%
  mutate(Treatment = str_replace_all(Treatment, c(ctrl = "PBS", lps = "LPS")))  %>%
  mutate(Treatment = fct_relevel(Treatment, c("PBS", "LPS"))) %>%
  mutate(`Basal accessibility` = fct_relevel(`Basal accessibility`, c("Low basal", "Equal basal", "High basal")))

a <- 
  balbPrepped |>
  mutate(Strain = str_replace_all(Strain, pattern = c("BALBcJ" = "Responsive", "C57BL6J" = "Non-responsive")))
b <- 
  c57Prepped |>
  mutate(Strain = str_replace_all(Strain, pattern = c("C57BL6J" = "Responsive", "BALBcJ" = "Non-responsive")))

fig6e <- full_join(a, b) |> 
  mutate(`Basal accessibility` = str_remove_all(
    `Basal accessibility`,
    " basal"
  )) |>
  mutate(`Basal accessibility` = fct_relevel(`Basal accessibility`, c("Low", "Equal", "High"))) |>
  ggplot(aes(
    x = Counts, y = Treatment, fill = Strain
  )) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = 4,
    alpha = 0.8,
    size = lineWidth * 0.4668623,
    color = "black"
  ) +
  facet_wrap(
    `Basal accessibility` ~ Treatment,
    strip.position = "left", 
    ncol = 1, scales = "free_y"
  ) +
  scale_fill_viridis_d(option = "rocket", end = 0.6) +
  xlab("Log2(ATAC-seq +1)") + ylab("Peak density") +
  labs(fill = "Lipopolysaccharide\nresponsiveness") +
  theme(
    text = element_text(size = textSize, colour = "black"),
    line = element_line(linewidth = lineWidth * 0.4668623, colour = "black"),
    rect = element_rect(linewidth = lineWidth * 0.4668623),
    plot.title  = element_text(size = textSize, colour = "black", hjust = 0),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title  = element_text(size = smallTextSize, colour = "black"),
    axis.text.x  = element_text(size = smallTextSize, colour = "black"),
    axis.text.y  = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.line.y = element_blank(),
    legend.key.size = unit(0.25, 'cm'),
    legend.position = "none",
    legend.text = element_text(size = smallTextSize, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = lineWidth * 0.4668623, colour = "black"),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(-0.5, "cm"),
    panel.ontop = TRUE,
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    strip.text = element_text(face = "plain", size = smallTextSize),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.clip = "off"
  )

# data for panel 6f pie chart
fig6f <- full_join(a, b) |> 
  mutate(`Basal accessibility` = str_remove_all(
    `Basal accessibility`,
    " basal"
  )) |>
  mutate(`Basal accessibility` = fct_relevel(`Basal accessibility`, c("Low", "Equal", "High"))) |>
  distinct(Position, .keep_all = T) |>
  group_by(`Basal accessibility`) %>%
  summarize(Count = n()) %>% ungroup() %>% 
  kable(format = "pandoc")

fig6f %>% clipr::write_clip()

# analysis for panel 6h table
vars <-
  list(
    contrast = c(
      "BALBCJ_F1_LPS-C57BL6J_F1_LPS",
      "BALBCJ_F0_LPS-C57BL6J_F0_LPS"
    ),
    group = c(".f1", ".f0"),
    labels = c("Mixed", "Trans", "Cis", "Same"),
    labelKey = c("First", "Second", "Both", "NotDiff")
  )

counts_lps <- inner_join(
  diffPeaks %>% filter(contrast == "BALBCJ_F1_LPS-C57BL6J_F1_LPS"),
  diffPeaks %>% filter(contrast == "BALBCJ_F0_LPS-C57BL6J_F0_LPS"), by = "Position", 
  suffix = c(vars$group[1], vars$group[2])) %>%
  select(contains("Position") | contains("normTag") & -contains("Control")) %>%
  select(where(~ !all(is.na(.x)))) %>%
  filter(if_any(.cols = 2:3, ~ .x > 4))

fig6h <- full_join(
  functionProcessPeaks(full_join, diffPeaks) %>%
    inner_join(., counts_lps) %>%
    inner_join(.,
               balbPrepped %>% pivot_wider(
                 id_cols = 1:2,
                 names_from = 3:4,
                 values_from = 5
               )) %>%
    select(DEG, `Basal accessibility`) %>%
    mutate("LPS_Induced" = "BALBcJ"),
  functionProcessPeaks(full_join, diffPeaks) %>%
    inner_join(., counts_lps) %>%
    inner_join(.,
               c57Prepped %>% pivot_wider(
                 id_cols = 1:2,
                 names_from = 3:4,
                 values_from = 5
               )) %>% 
    select(DEG, `Basal accessibility`) %>%
    mutate("LPS_Induced" = "C57BL6J")) %>%
  rename("Category" = 1, "Basal accessibility" = 2) %>% group_by(Category, `Basal accessibility`) %>%
  summarize(Count = n()) %>% ungroup() %>% arrange(`Basal accessibility`) %>%
  pivot_wider(names_from = Category, values_from = "Count") %>%
  select(`Basal accessibility`, Same, Cis, Trans, Mixed) %>%
  knitr::kable(format = "simple") 

fig6h %>% clipr::write_clip()

# export panels 6a, 6b, and 6e as a pdf file
fig6 <- (fig6a + fig6b) / fig6e
ggsave(filename = "figure6.pdf", fig6, height = 6, width = 6, units = "in")

# c57Prepped %>% pivot_wider(id_cols = 1:2, names_from = 3:4, values_from = 5) %>%
#   inner_join(coordinates, .) %>%
#   write_csv("peaks_lps_induced_c57bl6j.csv")
# 
# balbPrepped %>% pivot_wider(id_cols = 1:2, names_from = 3:4, values_from = 5) %>%
#   inner_join(coordinates, .) %>%
#   write_csv("peaks_lps_induced_balbcj.csv")
