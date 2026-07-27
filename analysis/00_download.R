## analysis/00_download.R -- fetch the Kossakowski et al. (2017) ESM dataset
##
## Kossakowski, J. J., Groot, P. C., Haslbeck, J. M. B., Borsboom, D., &
## Wichers, M. (2017). Data from 'Critical Slowing Down as a Personalized Early
## Warning Signal for Depression'. Journal of Open Psychology Data, 5: 1.
## https://doi.org/10.5334/jopd.29    OSF: https://osf.io/j4fg8    CC-BY 4.0
##
## The raw archive is NOT committed. It is 796 KB of someone's clinical diary,
## it is one command to fetch, and the checksum below pins exactly which bytes
## every downstream result was computed from. /data-raw is gitignored and
## nothing ever writes to it except this script.
##
## `osfr` is not used: it is a dependency with a large tree for what is one
## authenticated-free HTTP GET against a documented public API.

source("R/setup.R")

RAW  <- file.path(PROJ, "data-raw")
ZIP  <- file.path(RAW, "ESMdata.zip")
SHA  <- "465ae9862f6d8a1d10ff322dfc1a88a87b422988501555d9cbedd85fdb535c5f"
URL  <- "https://osf.io/download/c6xt4/"     # osfstorage/ESMdata.zip on node j4fg8

dir.create(RAW, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(ZIP)) {
  cat("downloading from OSF...\n")
  utils::download.file(URL, ZIP, mode = "wb", quiet = TRUE)
}

got <- tools::sha256sum(ZIP)[[1]]
if (!identical(unname(got), SHA)) {
  stop("checksum mismatch.\n  expected ", SHA, "\n  got      ", got,
       "\nThe OSF file changed, or the download is corrupt. Do not proceed.")
}
cat("checksum verified:", substr(SHA, 1, 16), "...\n")

utils::unzip(ZIP, exdir = RAW, overwrite = TRUE)
csv <- file.path(RAW, "ESMdata", "ESMdata.csv")
stopifnot("ESMdata.csv not found after unzip" = file.exists(csv))

d <- read.csv(csv, stringsAsFactors = FALSE)
cat(sprintf("ESMdata.csv: %d rows x %d columns\n", nrow(d), ncol(d)))
cat(sprintf("dates %s to %s (%d calendar days, %d with observations)\n",
            min(d$date), max(d$date),
            as.numeric(diff(range(as.Date(d$date, "%d/%m/%y")))) + 1,
            length(unique(d$date))))
cat("phases:", paste(names(table(d$phase)), table(d$phase), sep = "=", collapse = "  "), "\n")
