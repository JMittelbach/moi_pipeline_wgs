#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

allow_fallback <- identical(Sys.getenv("TACTCV_ALLOW_BINOMMIX_FALLBACK"), "1")
if (allow_fallback && !requireNamespace("moimix", quietly = TRUE)) {
  cat("moimix=unavailable; using explicit flexmix-equivalent pilot fallback\n")
  quit(save = "no", status = 0L)
}

# moimix imports these Bioconductor packages.  They are intentionally installed
# through BiocManager rather than pinned as Conda packages: the older Conda
# builds are not available for every supported platform (notably osx-arm64),
# while BiocManager selects the Bioconductor release matching this R runtime.
bioc_packages <- c(
  "SeqArray", "SeqVarTools", "IRanges", "S4Vectors", "GenomicRanges",
  "BiocParallel"
)
missing_bioc <- bioc_packages[!vapply(
  bioc_packages,
  requireNamespace,
  FUN.VALUE = logical(1),
  quietly = TRUE
)]
if (length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
  }
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE,
    quiet = TRUE
  )
}

if (!requireNamespace("moimix", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("R package 'remotes' is missing; rerun ./setup.sh to recreate/update the Conda environment")
  }
  # Gaia can reach GitHub over git, but its GitHub API/PAT lookup is not
  # reliable.  Install the pinned source revision through git instead of
  # remotes::install_github(), which requires the API endpoint.
  moimix_ref <- "802eaf1fab653690b1b1f1475c879b5189ee40ae"
  remotes::install_git(
    "https://github.com/bahlolab/moimix.git",
    ref = moimix_ref,
    upgrade = "never",
    dependencies = FALSE
  )
}

if (!requireNamespace("moimix", quietly = TRUE)) {
  stop("moimix installation did not produce an importable package")
}

cat(sprintf("moimix=%s\n", as.character(packageVersion("moimix"))))
