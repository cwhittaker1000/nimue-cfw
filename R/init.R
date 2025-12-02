
#' Check and set up initial values for vaccine model
#'
#' @inheritParams run
#'
#' @return Checked initial values data.frame
init <- function(population_list, seeding_cases, seeding_age_order = NULL, init = NULL){

  if (is.null(init)) {

    # Determine number of locations and convert population to matrix format
    if (is.list(population_list)) {
      N_locations <- length(population_list)
      pop_matrix <- do.call(cbind, population_list)  # 17 x N_locations
    } else if (is.matrix(population_list)) {
      pop_matrix <- population_list  # Assume already 17 x N_locations
      N_locations <- ncol(pop_matrix)
    } else {
      # Single location (backward compatibility)
      N_locations <- 1
      pop_matrix <- matrix(population_list, ncol = 1)
    }

    # Check dimensions
    if (nrow(pop_matrix) != 17) {
      stop("population must be divided up into 17x 5-year age bands spanning 0 to 80+")
    }

    # Check seeding_cases is correct length
    if (length(seeding_cases) == 1) {
      seeding_cases <- rep(seeding_cases, N_locations)
    } else if (length(seeding_cases) != N_locations) {
      stop("seeding_cases must be length 1 or length N_locations")
    }

    # Check seeding_age_order
    if (!is.null(seeding_age_order)) {
      if (!is.list(seeding_age_order)) {
        # Assume same order for all locations
        seeding_age_order <- rep(list(seeding_age_order), N_locations)
      } else if (length(seeding_age_order) != N_locations) {
        stop("seeding_age_order must be NULL, a vector, or a list of length N_locations")
      }
    }

    # Initialize output arrays: 17 x 6 x N_locations
    S_0 <- array(0, dim = c(17, 6, N_locations))
    E_0 <- array(0, dim = c(17, 6, N_locations))
    D_0 <- array(0, dim = c(17, 6, N_locations))
    R_0 <- array(0, dim = c(17, 6, N_locations))
    ICase_0 <- array(0, dim = c(17, 6, N_locations))
    IMild_0 <- array(0, dim = c(17, 6, N_locations))
    IHosp_0 <- array(0, dim = c(17, 6, N_locations))
    age_group_indices <- c(8, 9, 10, 11) # age_group indices corresponding to middle-aged travellers

    # Loop over locations
    for (loc in seq_len(N_locations)) {
      population_loc <- pop_matrix[, loc]
      seeding_cases_loc <- seeding_cases[loc]

      # Distribute seeds
      if (is.null(seeding_age_order)) {
        raw_seeding_cases <- rep(0, 17)
        raw_seeding_cases[age_group_indices] <- as.vector(
          stats::rmultinom(1, size = seeding_cases_loc, prob = rep(0.25, 4))
        )
      } else {
        seeding_age_order_loc <- seeding_age_order[[loc]]
        raw_seeding_cases <- rep(floor(seeding_cases_loc / 17), 17)
        for (s in seq_len(seeding_cases_loc %% 17)) {
          raw_seeding_cases[seeding_age_order_loc[s]] <-
            raw_seeding_cases[seeding_age_order_loc[s]] + 1
        }
      }

      # Set initial conditions for this location
      S <- population_loc - raw_seeding_cases
      E <- raw_seeding_cases

      # First vaccine group (unvaccinated) gets the initial infections
      S_0[, 1, loc] <- S
      E_0[, 1, loc] <- E

      # Check population conservation
      if (!all((rowSums(S_0[,,loc]) + rowSums(E_0[,,loc])) == population_loc)) {
        stop(paste0("Row sums of init should be identical to population for location ", loc))
      }

      init <- list(
        S_0 = S_0,
        E_0 = E_0,
        D_0 = D_0,
        R_0 = R_0,
        ICase_0 = ICase_0,
        IMild_0 = IMild_0,
        IHosp_0 = IHosp_0
      )

    }
  }
  return(init)
}
