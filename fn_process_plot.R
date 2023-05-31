read_plus <- function(flnm) {
  read_tsv(flnm, col_names = TRUE, col_types = "cddddddcdd") %>%
    mutate(sourceFile = flnm)
}

functionProcessData <- function(joinType, tpmData, diffData, variableList) {
  
  
  vars <- variableList
  
  tmpVar <- list(stat = c("log2FoldChange", "padj"))
  vars2 <-
    cross2(pluck(tmpVar, "stat"), pluck(vars, "group")) %>% map(lift(str_c)) %>% as_vector()
  fc <- keep(vars2, str_detect(vars2, "log2"))
  padj <- keep(vars2, str_detect(vars2, "padj"))
  labels <- vars$labels
  
  merged <-
    joinType(
      diffData %>% filter(., contrast == vars$contrast[1]),
      diffData %>% filter(., contrast == vars$contrast[2]),
      by = "Gene",
      suffix = c(vars$group[1], vars$group[2])
    ) %>%
    select(Gene | contains(tmpVar$stat)) %>%
    mutate(across(all_of(fc), ~ replace_na(., 0))) %>%
    mutate(across(all_of(padj), ~ replace_na(., 1))) %>%
    mutate(., DEG = factor(ifelse(
      abs(.data[[fc[1]]]) > 1 & .data[[padj[1]]] < 0.05 &
        (abs(.data[[fc[2]]]) < 1 | .data[[padj[2]]] > 0.05),
      labels[1],
      ifelse(
        (abs(.data[[fc[1]]]) < 1 | .data[[padj[1]]] > 0.05) &
          abs(.data[[fc[2]]]) > 1 & .data[[padj[2]]] < 0.05,
        labels[2],
        ifelse(
          abs(.data[[fc[1]]]) > 1 & .data[[padj[1]]] < 0.05 &
            abs(.data[[fc[2]]]) > 1 & .data[[padj[2]]] < 0.05,
          labels[3],
          labels[4]
        )
      )
    ))) %>%
    mutate(DEG = fct_relevel(DEG, c(labels[4], labels[3], labels[2], labels[1]))) %>%
    inner_join(., tpmData)
}

functionProcessPeaks <- function(joinType, diffData) {
  
  
  tmpVar <- list(stat = c("log2FoldChange", "padj"))
  vars2 <-
    cross2(pluck(tmpVar, "stat"), pluck(vars, "group")) %>% map(lift(str_c)) %>% as_vector()
  fc <- keep(vars2, str_detect(vars2, "log2"))
  padj <- keep(vars2, str_detect(vars2, "padj"))
  labels <- vars$labels
  
  merged <-
    joinType(
      diffData %>% filter(., contrast == vars$contrast[1]),
      diffData %>% filter(., contrast == vars$contrast[2]),
      by = "Position",
      suffix = c(vars$group[1], vars$group[2])
    ) %>%
    select(Position | contains(tmpVar$stat)) %>%
    mutate(across(all_of(fc), ~ replace_na(., 0))) %>%
    mutate(across(all_of(padj), ~ replace_na(., 1))) %>%
    mutate(., DEG = factor(ifelse(
      abs(.data[[fc[1]]]) > 1 & .data[[padj[1]]] < 0.05 &
        abs(.data[[fc[2]]]) < 1 & .data[[padj[2]]] > 0.05,
      labels[1],
      ifelse(
        abs(.data[[fc[1]]]) < 1 & .data[[padj[1]]] > 0.05 &
          abs(.data[[fc[2]]]) > 1 & .data[[padj[2]]] < 0.05,
        labels[2],
        ifelse(
          abs(.data[[fc[1]]]) > 1 & .data[[padj[1]]] < 0.05 &
            abs(.data[[fc[2]]]) > 1 & .data[[padj[2]]] < 0.05,
          labels[3],
          labels[4]
        )
      )
    ))) %>%
    mutate(DEG = fct_relevel(DEG, c(labels[4], labels[3], labels[2], labels[1])))
}

functionRatioPlot <- function(processedData, textSize, smallTextSize, lineWidth, pointSize) {
  
  textSize <- textSize
  smallTextSize <- smallTextSize
  lineWidth <- lineWidth * 0.4668623
  pointSize <- pointSize
  
  tmpVar <- list(stat = c("log2FoldChange", "padj"))
  vars2 <-
    cross2(pluck(tmpVar, "stat"), pluck(vars, "group")) %>% map(lift(str_c)) %>% as_vector()
  fc <- keep(vars2, str_detect(vars2, "log2"))
  padj <- keep(vars2, str_detect(vars2, "padj"))
  labels <- vars$labels
  
  pal_p3 <- c(viridis_pal(option = "viridis")(4))
  grob_p3 <- grobTree(
    grid.points(
      x = 0.05,
      y = 0.95,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[4]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[4], ": ",
            select(processedData, DEG) %>%
              filter(DEG == labels[4]) %>%
              tally(),
            sep = ""
      ),
      x = 0.05,
      y = 0.95,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = smallTextSize)
    ),
    grid.points(
      x = 0.05,
      y = 0.87,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[1]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[3], ": ",
            select(processedData, DEG) %>%
              filter(DEG == labels[3]) %>%
              tally(),
            sep = ""
      ),
      x = 0.05,
      y = 0.87,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = smallTextSize)
    ),
    grid.points(
      x = 0.05,
      y = 0.8,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[3]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[2], ": ",
            select(processedData, DEG) %>%
              filter(DEG == labels[2]) %>%
              tally(),
            sep = ""
      ),
      x = 0.05,
      y = 0.8,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = smallTextSize)
    ),
    grid.points(
      x = 0.05,
      y = 0.73,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[2]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[1], ": ",
            select(processedData, DEG) %>%
              filter(DEG == labels[1]) %>%
              tally(),
            sep = ""
      ),
      x = 0.05,
      y = 0.73,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = smallTextSize)
    )
  )
  
  .x <- sym(fc[2])
  .y <- sym(fc[[1]])
  
  processedData %>% 
    arrange(match(DEG, rev(labels))) %>%
    ggplot(.,
           aes(
             x = !!(.x),
             y = !!(.y),
             colour = DEG
           )) +
    ggrastr::rasterise(geom_point(size = pointSize, shape = "bullet"), dpi = 600) +
    geom_abline(slope = 0, size = lineWidth) +
    geom_abline(slope = 1, size = lineWidth) +
    geom_vline(xintercept = 0, size = lineWidth) +
    theme_bw() +
    theme(
      text = element_text(size = smallTextSize, colour = "black"),
      line = element_line(size = lineWidth, colour = "black"),
      rect = element_rect(size = lineWidth, colour = "black"),
      plot.title  = element_text(size = textSize, colour = "black"),
      axis.title  = element_text(size = smallTextSize, colour = "black"),
      axis.ticks = element_line(size = lineWidth, colour = "black"),
      axis.text  = element_text(size = smallTextSize, colour = "black"),
      axis.line = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = lineWidth, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm"),
      plot.background = element_blank()
    ) +
    scale_color_manual(values = c("#FDE725FF", "#440154FF", "#35B779FF", "#31688EFF")) +
    annotation_custom(grob_p3) 
}

mergeProcessRatioFunction <- function() {
  functionProcessData(full_join, tpm, diffRNA, vars) %>%
    functionRatioPlot(.,
                      textSize = textSize,
                      smallTextSize = smallTextSize,
                      lineWidth = lineWidth,
                      pointSize = pointSize
    )
}