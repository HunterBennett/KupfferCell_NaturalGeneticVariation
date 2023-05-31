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