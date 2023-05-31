#Clear plots
clearSession <- function() {
  if (!is.null(dev.list()))
    dev.off()
  # Clear console
  cat("\014")
  # Clean workspace
  rm(list = ls())
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
    inner_join(., tpmData) %>%
    mutate(meanTPM = log2(rowMeans(select(
      ., contains(vars$lowerMean) | contains(vars$upperMean)
    ))
    + 1)) %>%
    arrange(., desc(DEG))
}

functionMAPlot <- function(processedData, variableList, textSize, lineWidth, pointSize) {
  vars <- variableList
  processedData %>% 
    mutate(DEG = fct_relevel(DEG, vars$labels[4], vars$labels[3], vars$labels[2], vars$labels[1])) %>%
    arrange(DEG) %>%
    ggplot(aes(
      x = meanTPM,
      y = log2FoldChangeKC,
      color = DEG
    )) +
    # geom_point(size = 1) +
    scale_size_area() +
    scale_color_viridis_d(option = "rocket", begin = 0.85, end = 0) +
    # annotation_custom(upper1) +
    # annotation_custom(lower1) +
    scale_y_continuous(expand = c(.25, .05)) +
    theme(
      text = element_text(size = smallPointSize, colour = "black"),
      plot.title  = element_text(size = smallPointSize, colour = "black"),
      axis.ticks.x = element_line(size = lineWidth * 0.4668623, colour = "black"),
      axis.ticks.y = element_line(size = lineWidth * 0.4668623, colour = "black"),
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
      panel.border = element_blank(),
      axis.line = element_line(size = lineWidth * 0.4668623, colour = "black"),
      plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm")
    ) +
    geom_abline(slope = 0, size = lineWidth * 0.4668623) +
    ggrastr::rasterise(geom_point(size = 0.1, shape = "bullet"), dpi = 600) #+
    #geom_text_repel(
    #  data = subset(processedData, DEG ==vars$labels[1]),
    #  aes(x = meanTPM,
    #      y = log2FoldChangeKC, label = Gene),
    #  min.segment.length = 0,
    #  box.padding = 0.3,
    #  segment.size  = 0.25 * 0.4668623,
    #  fontface = "italic",
    #  segment.linetype = 1,
    #  color = "gray21",
    #  max.overlaps = 11,
    #  seed = 42,
    #  force = 1,
    #  max.time = 310,
    #  max.iter = 1000000,
    #  size = smallPointSize * 0.4668623
    #)
}

functionMergeProcessMAPlot <- function(variableList) {
  functionProcessData(left_join, tpm, diff, variableList) %>%
    functionMAPlot(
      processedData = .,
      lineWidth = lineWidth,
      pointSize = smallPointSize,
      variableList
    )
  }

barChartFunction <- function(varList) {
  tpmSummarized <-
    .tpm %>%
    filter(Gene %in% varList) %>%
    summarize(mean = mean(TPM),
              sd = sd(TPM),
              n = n())
  tpmNonSummarized <- .tpm %>%
    filter(Gene %in% varList)
  ggplot() +
    geom_bar(
      data = tpmSummarized,
      aes(y = mean, x = Cell, fill = Strain),
      stat = "identity",
      position = position_dodge()
    ) +
    geom_point(
      data = tpmNonSummarized,
      aes(y = TPM, x = Cell, group = Strain),
      alpha = 1,
      size = 0.1,
      shape = 21,
      fill = "black",
      colour = "black",
      position = position_dodge(width = 1)
    ) + coord_flip() +
    # geom_errorbar(aes(xmax = mean + sd, xmin = mean - sd),
    #               width = 0,
    #               size = 0.25, position = position_dodge(0.9)) +
    facet_wrap(~ Gene,
               nrow = 1,
               scales = "free_x") +
    scale_fill_manual(values = c("#C42503", "#455BCD", "#A2FC3C")) +
    scale_colour_manual(values = c("#C42503", "#455BCD", "#A2FC3C")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    xlab("") + ylab("Mean TPM") +
    theme(
      text = element_text(size = smallPointSize, colour = "black"),
      line = element_line(size = lineWidth * 0.4668623, colour = "black"),
      rect = element_rect(size = lineWidth * 0.4668623),
      axis.ticks.x = element_line(size = lineWidth * 0.4668623, colour = "black"),
      axis.ticks.y = element_line(size = lineWidth * 0.4668623, colour = "black"),
      axis.text.x  = element_text(
        hjust = 1,
        vjust = 1,
        angle = 45,
        size = smallPointSize,
        colour = "black"
      ),
      axis.text.y  = element_text(size = smallPointSize, colour = "black"),
      axis.title = element_text(size = smallPointSize, colour = "black"),
      # plot.title = element_text(hjust = 0, face = "italic"),
      plot.subtitle = element_text(hjust = 0),
      legend.position = "right",
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
      # strip.text.x = element_text(
      #   face = "italic",
      #   size = smallPointSize,
      #   hjust = 0
      # ),
      strip.text = element_text(
        # vjust = -2.0,
        margin = margin(0,0,0,0),
        face = "italic",
        size = smallPointSize
      ),
      # panel.margin.y = unit(-2, "lines"),
      # strip.text.y = element_text(
      #   face = "plain",
      #   size = smallPointSize,
      #   hjust = 0.5,
      #   margin = unit(0, 0, 0, 0, "mm")
      # ),
      strip.placement = "outside",
      strip.clip = "off",
      strip.background = element_blank()
    )
}