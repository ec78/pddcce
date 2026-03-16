library(haven)

d <- read_dta("examples/penn_sample.dta")
mapping <- read.csv("validation/id_country_map_gdponly.csv", stringsAsFactors=FALSE)

# Use short country names for GAUSS examples (avoid very long official names)
name_overrides <- c(
  "Bolivia (Plurinational State of)" = "Bolivia",
  "Congo, Democratic Republic"       = "DR Congo",
  "China, Hong Kong SAR"             = "Hong Kong",
  "Iran (Islamic Republic of)"       = "Iran",
  "Republic of Korea"                = "Korea",
  "Syrian Arab Republic"             = "Syria",
  "U.R. of Tanzania: Mainland"       = "Tanzania",
  "United States of America"         = "USA",
  "Venezuela (Bolivarian Republic of)" = "Venezuela"
)

mapping$name <- mapping$country
for (long in names(name_overrides)) {
  mapping$name[mapping$country == long] <- name_overrides[long]
}

cat("Country mapping (id -> name):\n")
print(mapping[order(mapping$id), c("id","isocode","name")])

# Merge country name into data
d2 <- merge(d, mapping[, c("id","isocode","name")], by="id", all.x=TRUE)
cat("\nAny unmapped IDs:", sum(is.na(d2$name)), "\n")

# Replace numeric id with string country name
# Keep original id as backup column
d2$id_orig <- d2$id
d2$id <- d2$name

# Reorder columns: id (string), year, then data columns
d2 <- d2[, c("id", "year", "log_rgdpo", "log_hc", "log_ck", "log_ngd")]
d2 <- d2[order(d2$id, d2$year), ]

cat("\nNew data head:\n")
print(head(d2))
cat("\nDims:", nrow(d2), "x", ncol(d2), "\n")
cat("Unique country names:", length(unique(d2$id)), "\n")

# Save as Stata .dta file
write_dta(d2, "examples/penn_world.dta")
cat("\nSaved to examples/penn_world.dta\n")
cat("id column is now string country name\n")
cat("Sort order: by country name then year\n")
