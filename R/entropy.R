#' Shannon entropy of a character vector
#'
#' Computes the Shannon entropy (in bits, log base 2) of the frequency
#' distribution of elements in \code{x}. Higher values indicate greater
#' diversity among monomer sequences.
#'
#' @param x Character vector of monomer sequences.
#'
#' @return A single non-negative numeric value (entropy in bits).
#' @export
#'
#' @examples
#' shannon_entropy(c("AAA", "AAA", "AAA"))  # 0 — all identical
#' shannon_entropy(c("AAA", "CCC", "GGG"))  # log2(3) — all distinct
shannon_entropy <- function(x) {
  p <- table(x)
  p <- p[p > 0] / sum(p)
  -sum(p * log2(p))
}

#' Summarise entropy for a simulation result
#'
#' Extracts the Shannon entropy, total array size, and unique-monomer
#' count from a \code{\link{sim_array}} result.
#'
#' @param sim Output from \code{\link{sim_array}}.
#'
#' @return A one-row \code{data.table} with columns:
#'   \describe{
#'     \item{H}{Shannon entropy of the final array (bits).}
#'     \item{total}{Total number of monomers in the final array.}
#'     \item{uniqueN}{Number of distinct monomer sequences.}
#'     \item{meanN}{Mean copy number per unique sequence (\code{total / uniqueN}).}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- sim_array(max_generations = 200, mu_total = 3e-5)
#' get_sim_entropies(sim)
#' }
get_sim_entropies <- function(sim) {
  mono <- sim$monomers
  dt <- data.table(
    H       = shannon_entropy(mono$bponly),
    total   = nrow(mono),
    uniqueN = uniqueN(mono$bponly)
  )
  dt[, meanN := total / uniqueN]
  dt
}
