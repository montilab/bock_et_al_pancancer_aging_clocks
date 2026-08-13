reactable_pretty <- function(results) {
  results %>%
    dplyr::mutate(across(where(~ is.numeric(.) && any(. %% 1 != 0, na.rm = TRUE)), ~ signif(., digits = 3))) %>%
    reactable::reactable(
      searchable = TRUE, resizable = TRUE, striped = TRUE, highlight = TRUE, showPageSizeOptions = TRUE,
      pageSizeOptions = c(10, 25, 50, 100)
    )
}

kable_pretty <- function(df, font_size = 10, full_width = TRUE, ...) {
  df %>%
    knitr::kable(format = "html", escape = FALSE, ...) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "condensed"),
      full_width = full_width,
      font_size = font_size,
      position = "left"
    )
}
