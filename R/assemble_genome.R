# ============================================================================
# Real-data ingestion: read per-chromosome monomer fastas into the standard
# genome data contract (one row per monomer: chrom, num, bponly, dir, aligned)
# -- the same shape sim_genome() produces, so downstream analyses can score
# simulated and observed genomes identically.
# ============================================================================

#' Assemble an observed single-haplotype genome from per-chromosome fastas
#'
#' Reads per-chromosome monomer fastas for ONE taxon / ONE haplotype and
#' returns the standard genome data contract shared with \code{sim_genome()}:
#' one row per monomer with integer chromosome rank, within-chromosome
#' position, gapless sequence, orientation, plus the gapped sequence (column
#' \code{aligned}) for divergence computation.
#'
#' @param files Character vector of per-chromosome fasta paths.
#' @param chrom_field,num_field,dir_field 1-based header fields (split on
#'   \code{sep}) holding the chromosome id, monomer position, and strand.
#' @param sep Header field separator.
#' @return data.table with columns \code{chrom}, \code{num}, \code{bponly},
#'   \code{dir}, \code{aligned}.
#' @export
assemble_genome_dt <- function(files, chrom_field = 3L, num_field = 4L,
                               dir_field = 7L, sep = "_") {
  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("need Biostrings")
  parts <- lapply(files, function(f) {
    dna <- Biostrings::readDNAStringSet(f); raw <- as.character(dna)
    nm  <- names(dna); hd <- strsplit(nm, sep, fixed = TRUE)
    chr <- vapply(hd, `[`, "", chrom_field)
    data.table(chrom_id = chr,
               num     = as.numeric(vapply(hd, `[`, "", num_field)),
               dir     = vapply(hd, `[`, "", dir_field),
               bponly  = gsub("-", "", raw, fixed = TRUE),
               aligned = raw)
  })
  dt <- rbindlist(parts)
  # integer chrom rank (gapless 1..K per contract), preserving Chr-id order
  dt[, chrom := as.integer(factor(chrom_id, levels = sort(unique(chrom_id))))]
  dt[, chrom_id := NULL]
  setorder(dt, chrom, num)
  dt[]
}
