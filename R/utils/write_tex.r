# rounding helper
mround <- function(x, digits) sprintf(paste0("%.", digits, "f"), round(x, digits))

# LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    readr::write_lines("output/values/species_values.tex", append = append)
}