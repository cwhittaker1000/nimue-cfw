## Helper: compute mapping from canonical 17 age-bands to user groups
##
## Canonical bands: [0,5), [5,10), ..., [75,80), [80,100)
compute_age_mapping <- function(age_breaks = NULL) {

  # Canonical 17-group boundaries used by squire
  base_breaks <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100)
  base_lower  <- base_breaks[-length(base_breaks)]  # 17 lower bounds

  if (is.null(age_breaks)) {
    # Default: keep the 17 canonical groups
    age_breaks <- base_breaks
  }
  if (!is.numeric(age_breaks)) {
    stop("`age_breaks` must be numeric.")
  }
  if (length(age_breaks) < 2L) {
    stop("`age_breaks` must have length >= 2.")
  }
  if (age_breaks[1] != 0 || tail(age_breaks, 1) != 100) {
    stop("`age_breaks` must start at 0 and end at 100.")
  }
  if (is.unsorted(age_breaks, strictly = TRUE)) {
    stop("`age_breaks` must be strictly increasing.")
  }
  if (!all(age_breaks %in% base_breaks)) {
    stop("`age_breaks` must be multiples of 5 and align with 0, 5, ..., 80, 100.")
  }

  # Map each canonical 5-year band to a new group index
  mapping <- findInterval(base_lower, age_breaks,
                          rightmost.closed = FALSE,
                          all.inside = TRUE)

  list(
    mapping    = mapping,         # length 17, values in 1:(length(age_breaks)-1)
    base_breaks = base_breaks,
    age_breaks = age_breaks,
    N_age      = length(age_breaks) - 1L
  )
}

## Helper: aggregate a 17-length age vector to new groups, population-weighted
aggregate_age_vector <- function(x, mapping, N_age_new, weights, name = "vector") {

  if (length(x) != length(mapping)) {
    stop("`", name, "` must have length 17 (canonical age groups).")
  }
  if (length(weights) != length(mapping)) {
    stop("`weights` must have length 17 (canonical age groups).")
  }

  out <- numeric(N_age_new)
  for (g in seq_len(N_age_new)) {
    idx <- which(mapping == g)
    w   <- weights[idx]
    if (sum(w) == 0) {
      # Fallback: if somehow no population in this group, just average
      out[g] <- mean(x[idx])
    } else {
      out[g] <- sum(x[idx] * w) / sum(w)
    }
  }
  out
}


## Helper: aggregate population + 17x17 contact matrix to new groups
aggregate_population_and_matrix <- function(pop17, C17, mapping, N_age_new) {

  if (length(pop17) != length(mapping)) {
    stop("Population length must match number of canonical age groups (17).")
  }
  if (!all(dim(C17) == c(length(mapping), length(mapping)))) {
    stop("Contact matrix must be 17 x 17 after processing.")
  }

  # Contact numbers at fine resolution: M_ab = C_ab * N_a
  M_old <- C17 * pop17  # row-wise scaling

  pop_new <- numeric(N_age_new)
  M_new   <- matrix(0, nrow = N_age_new, ncol = N_age_new)

  for (g in seq_len(N_age_new)) {
    idx_g <- which(mapping == g)
    pop_new[g] <- sum(pop17[idx_g])
  }

  for (g in seq_len(N_age_new)) {
    idx_g <- which(mapping == g)
    for (h in seq_len(N_age_new)) {
      idx_h <- which(mapping == h)
      M_new[g, h] <- sum(M_old[idx_g, idx_h, drop = FALSE])
    }
  }

  # Convert back to per-capita rates: C*_gh = M*_gh / N_g
  C_new <- M_new
  for (g in seq_len(N_age_new)) {
    if (pop_new[g] > 0) {
      C_new[g, ] <- C_new[g, ] / pop_new[g]
    } else {
      C_new[g, ] <- 0
    }
  }

  list(
    population = pop_new,
    contact    = C_new
  )
}
