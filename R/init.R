
#' Check and set up initial values for vaccine model
#'
#' @inheritParams run
#'
#' @return Checked initial values data.frame
init <- function(population_list, seeding_cases, seeding_age_order = NULL, init = NULL){

  if (is.null(init)) {

    ## ------------------------------------------------------------------------
    ## 1. Convert population_list to a matrix: N_age x N_locations
    ## ------------------------------------------------------------------------
    if (is.list(population_list)) {
      N_locations <- length(population_list)
      pop_matrix  <- do.call(cbind, population_list)  # N_age x N_locations
    } else if (is.matrix(population_list)) {
      pop_matrix  <- population_list
      N_locations <- ncol(pop_matrix)
    } else {
      # Single location (backward compatibility)
      N_locations <- 1
      pop_matrix  <- matrix(population_list, ncol = 1)
    }

    N_age <- nrow(pop_matrix)
    if (N_age < 1) {
      stop("Population must have at least one age group.")
    }

    ## ------------------------------------------------------------------------
    ## 2. Check / expand seeding_cases
    ## ------------------------------------------------------------------------
    if (!is.numeric(seeding_cases)) {
      stop("`seeding_cases` must be numeric.")
    }
    if (any(seeding_cases < 0)) {
      stop("`seeding_cases` must be >= 0 for all locations.")
    }

    if (length(seeding_cases) == 1L) {
      seeding_cases <- rep(seeding_cases, N_locations)
    } else if (length(seeding_cases) != N_locations) {
      stop("`seeding_cases` must be length 1 or length N_locations (", N_locations, ").")
    }
    ## ------------------------------------------------------------------------
    ## 3. Check / expand seeding_age_order
    ## ------------------------------------------------------------------------
    if (!is.null(seeding_age_order)) {
      if (!is.list(seeding_age_order)) {
        # Assume same order for all locations
        seeding_age_order <- rep(list(seeding_age_order), N_locations)
      } else if (length(seeding_age_order) != N_locations) {
        stop("`seeding_age_order` must be NULL, a vector, or a list of length N_locations")
      }

      # Optional: sanity check that indices are within [1, N_age]
      for (loc in seq_len(N_locations)) {
        idx <- seeding_age_order[[loc]]
        if (any(idx < 1 | idx > N_age)) {
          stop("All indices in seeding_age_order[[", loc,
               "]] must be between 1 and N_age = ", N_age)
        }
      }
    }

    ## ------------------------------------------------------------------------
    ## 4. Allocate state arrays: N_age x N_vaccine x N_locations
    ## ------------------------------------------------------------------------
    N_vaccine <- 4L  # consistent with parameters() where N_vaccine = 4

    S_0     <- array(0, dim = c(N_age, N_vaccine, N_locations))
    E_0     <- array(0, dim = c(N_age, N_vaccine, N_locations))
    D_0     <- array(0, dim = c(N_age, N_vaccine, N_locations))
    R_0     <- array(0, dim = c(N_age, N_vaccine, N_locations))
    ICase_0 <- array(0, dim = c(N_age, N_vaccine, N_locations))
    IMild_0 <- array(0, dim = c(N_age, N_vaccine, N_locations))
    IHosp_0 <- array(0, dim = c(N_age, N_vaccine, N_locations))

    ## Default seeding age groups (when seeding_age_order is NULL):
    ## - If N_age == 17: keep original behaviour (indices 8–11).
    ## - Otherwise: choose a "middle" block of up to 4 age groups.
    if (N_age == 17) {
      default_seed_indices <- c(8, 9, 10, 11)
    } else {
      num_mid <- min(4L, N_age)
      start   <- floor((N_age - num_mid) / 2) + 1L
      default_seed_indices <- seq.int(start, length.out = num_mid)
    }

    ## ------------------------------------------------------------------------
    ## 5. Loop over locations and seed infections
    ## ------------------------------------------------------------------------
    for (loc in seq_len(N_locations)) {

      population_loc     <- pop_matrix[, loc]
      seeding_cases_loc  <- seeding_cases[loc]

      # Distribute seeds over age groups
      if (is.null(seeding_age_order)) {
        # Default: seed into a "middle-aged" block of age groups
        raw_seeding_cases <- rep(0, N_age)
        n_seed_groups     <- length(default_seed_indices)
        if (seeding_cases_loc > 0 && n_seed_groups > 0) {
          seed_draw <- stats::rmultinom(
            n    = 1,
            size = seeding_cases_loc,
            prob = rep(1 / n_seed_groups, n_seed_groups)
          )
          raw_seeding_cases[default_seed_indices] <- as.vector(seed_draw)
        }
      } else {
        # User-specified seeding order for this location
        seeding_age_order_loc <- seeding_age_order[[loc]]
        raw_seeding_cases     <- rep(floor(seeding_cases_loc / N_age), N_age)
        remainder             <- seeding_cases_loc %% N_age
        if (remainder > 0) {
          for (s in seq_len(remainder)) {
            raw_seeding_cases[seeding_age_order_loc[s]] <- raw_seeding_cases[seeding_age_order_loc[s]] + 1
          }
        }
      }

      # Set initial S and E for this location (unvaccinated stratum)
      S <- population_loc - raw_seeding_cases
      E <- raw_seeding_cases

      S_0[, 1, loc] <- S
      E_0[, 1, loc] <- E

      # Check population conservation
      if (!all((rowSums(S_0[,, loc]) + rowSums(E_0[,, loc])) == population_loc)) {
        stop(paste0("Row sums of init should be identical to population for location ", loc))
      }
    }

    ## ------------------------------------------------------------------------
    ## 6. Wrap up
    ## ------------------------------------------------------------------------
    init <- list(
      S_0     = S_0,
      E_0     = E_0,
      D_0     = D_0,
      R_0     = R_0,
      ICase_0 = ICase_0,
      IMild_0 = IMild_0,
      IHosp_0 = IHosp_0
    )
  }

  return(init)
}

