library(haven)

d <- read_dta("examples/penn_sample.dta")
cat("Attrs of id column:\n")
print(attributes(d$id))
cat("\nAll column attrs:\n")
for (col in names(d)) {
  attrs <- attributes(d[[col]])
  if (!is.null(attrs$labels) || !is.null(attrs$label)) {
    cat(col, ":", attrs$label, "\n")
    if (!is.null(attrs$labels)) print(attrs$labels)
  }
}

# Check jasa2
cat("\n\n--- jasa2.dta ---\n")
j <- read_dta("examples/jasa2.dta")
cat("Cols:", paste(names(j), collapse=", "), "\n")
cat("Dims:", nrow(j), "x", ncol(j), "\n")
# Look for string/country columns
for (col in names(j)) {
  if (is.character(j[[col]])) {
    cat("String col:", col, "- sample:", paste(head(unique(j[[col]]), 5), collapse=", "), "\n")
  }
}
print(head(j, 3))
