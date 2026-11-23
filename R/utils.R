# Internal ggplot2 axis helper: approximately-n pretty integer breaks
# (used by plot_repeat_fingerprint).
int_breaks_pretty <- function(n = 5) {
  function(x) {
    b <- pretty(range(x, na.rm = TRUE), n)
    b[b == round(b)]
  }
}
