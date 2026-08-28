#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: 07_binommix.R COUNTS.tsv SAMPLE_ID K_COMMA_LIST COVERAGE_THRESHOLD NITER OUTPUT.tsv")
}

counts_path <- args[[1L]]
sample_id <- args[[2L]]
k_values <- as.integer(strsplit(args[[3L]], ",", fixed = TRUE)[[1L]])
coverage_threshold <- as.integer(args[[4L]])
niter <- as.integer(args[[5L]])
output_path <- args[[6L]]

has_moimix <- requireNamespace("moimix", quietly = TRUE)
if (!has_moimix && !requireNamespace("flexmix", quietly = TRUE)) {
  stop("neither R package 'moimix' nor 'flexmix' is installed; run ./setup.sh")
}
if (length(k_values) < 1L || any(!is.finite(k_values)) || any(k_values < 1L | k_values > 5L)) {
  stop("moimix k values must be integers in 1..5")
}
if (coverage_threshold < 0L || niter < 1L) stop("invalid coverage threshold or iteration count")

required <- c("chrom", "pos", "ref", "alt", "AD", "DP")
first_line <- readLines(counts_path, n = 1L, warn = FALSE)
has_header <- length(first_line) == 1L &&
  identical(strsplit(first_line[[1L]], "\t", fixed = TRUE)[[1L]], required)
input <- read.delim(counts_path, header = has_header, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "")
if (!has_header) names(input) <- required
if (!identical(names(input), required)) stop("counts file has an unexpected header")

split_ad <- strsplit(as.character(input$AD), ",", fixed = TRUE)
valid_ad <- vapply(split_ad, function(x) length(x) == 2L && all(grepl("^[0-9]+$", x)), logical(1))
valid_dp <- grepl("^[0-9]+$", as.character(input$DP))
input <- input[valid_ad & valid_dp, , drop = FALSE]
if (nrow(input) == 0L) stop("no valid allele-count rows")
ad_ref <- as.integer(vapply(strsplit(as.character(input$AD), ",", fixed = TRUE), `[`, "", 1L))
ad_alt <- as.integer(vapply(strsplit(as.character(input$AD), ",", fixed = TRUE), `[`, "", 2L))
dp <- as.integer(input$DP)
keep <- dp >= coverage_threshold & (ad_ref + ad_alt) > 0L
if (!any(keep)) stop("no sites pass coverage threshold")
ad_ref <- ad_ref[keep]
ad_alt <- ad_alt[keep]
site_ids <- paste(input$chrom[keep], input$pos[keep], sep = ":")

# moimix::binommix expects an alleleCounts object with one row per sample and
# one column per site.  The object is deliberately constructed from the fixed
# REF/ALT count table rather than from a sample-derived VCF.
ref <- matrix(ad_ref, nrow = 1L, dimnames = list(sample = sample_id, site = site_ids))
alt <- matrix(ad_alt, nrow = 1L, dimnames = list(sample = sample_id, site = site_ids))
counts <- structure(
  list(
    ref = ref,
    alt = alt,
    # moimix uses dosage==1 to retain heterozygous/informative sites.  The
    # fixed REF/ALT panel is already restricted to biallelic markers, so the
    # explicit dosage matrix is one for every retained site.
    dosage = matrix(1L, nrow = 1L, ncol = length(site_ids),
                    dimnames = list(sample = sample_id, site = site_ids))
  ),
  class = "alleleCounts"
)

result <- tryCatch({
  if (has_moimix) {
    fit <- moimix::binommix(
      counts,
      sample.id = sample_id,
      k = k_values,
      coverage_threshold = 0L,
      niter = niter
    )
  } else {
    # Exact binomial-mixture model used by moimix::binommix, kept explicit for
    # macOS ARM smoke tests where gdsfmt/SeqArray cannot currently compile.
    # The count matrices are sample-by-site; flexmix expects one observation
    # (and a two-column success/failure response) per site.
    y_obs <- cbind(as.numeric(alt), as.numeric(ref))
    fit <- list(
      fits = flexmix::initFlexmix(
        y_obs ~ 1,
        k = k_values,
        model = flexmix::FLXMRglm(y_obs ~ 1, family = "binomial"),
        control = list(iter.max = niter, minprior = 0),
        nrep = 5
      ),
      fallback = TRUE
    )
  }

  fits <- fit$fits
  if (is.null(fits)) stop("binommix returned no candidate fits")
  if (inherits(fits, "stepFlexmix")) {
    # Both moimix and the explicit fallback return a stepFlexmix object.
    # BIC() is vector-valued for this S4 class; iterating over the object
    # itself fails on current flexmix versions.
    bic <- as.numeric(BIC(fits))
    best <- which.min(ifelse(is.finite(bic), bic, Inf))
    best_k <- as.integer(fits@k[[best]])
    best_fit <- flexmix::getModel(fits, which(fits@k == best_k)[[1L]])
  } else {
    bic <- vapply(fits, BIC, numeric(1))
    best <- which.min(ifelse(is.finite(bic), bic, Inf))
    best_k <- as.integer(k_values[[best]])
    best_fit <- fits[[best]]
  }
  if (!any(is.finite(bic))) stop("all candidate BinomMix BIC values are non-finite")
  theta <- if (has_moimix) tryCatch(
    moimix::getTheta(fits, criterion = "BIC"),
    error = function(e) moimix::getTheta(best_fit)
  ) else list()
  if (has_moimix) {
    pi_hat <- if (!is.null(theta$pi.hat)) paste(format(theta$pi.hat, digits = 12), collapse = ";") else ""
    mu_hat <- if (!is.null(theta$mu.hat)) paste(format(theta$mu.hat, digits = 12), collapse = ";") else ""
  } else {
    pi_hat <- paste(format(flexmix::prior(best_fit), digits = 12), collapse = ";")
    mu_hat <- paste(format(as.numeric(flexmix::parameters(best_fit)), digits = 12), collapse = ";")
  }
  # moimix does not expose a formal MOI confidence interval.  Keep the
  # relative BIC support so the summary can show an honest model-selection
  # confidence proxy without calling it a calibrated probability.
  finite_bic <- bic[is.finite(bic)]
  bic_delta <- if (length(finite_bic) > 1L) sort(finite_bic)[[2L]] - min(finite_bic) else NA_real_
  bic_weights <- if (length(finite_bic)) {
    weights <- exp(-0.5 * (bic - min(finite_bic)))
    weights[!is.finite(weights)] <- 0
    weights / sum(weights)
  } else {
    rep(NA_real_, length(bic))
  }
  data.frame(
    sample_id = sample_id,
    status = "estimated",
    reason = if (has_moimix) "" else "moimix_unavailable_flexmix_equivalent",
    model_k = best_k,
    callable_sites = ncol(ref),
    bic = as.numeric(bic[[best]]),
    bic_delta = bic_delta,
    bic_weight = as.numeric(bic_weights[[best]]),
    pi_hat = pi_hat,
    mu_hat = mu_hat,
    stringsAsFactors = FALSE
  )
}, error = function(error) {
  data.frame(
    sample_id = sample_id,
    status = "abstain_model_failure",
    reason = conditionMessage(error),
    model_k = NA_integer_,
    callable_sites = ncol(ref),
    bic = NA_real_,
    bic_delta = NA_real_,
    bic_weight = NA_real_,
    pi_hat = "",
    mu_hat = "",
    stringsAsFactors = FALSE
  )
})

write.table(result, output_path, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "", col.names = TRUE)
