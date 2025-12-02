
#' Check and set up initial values for vaccine model
#'
#' @inheritParams run
#'
#' @return Checked initial values data.frame
init <- function(population, seeding_cases, seeding_age_order = NULL,
                 init = NULL){

  if (is.null(init)) {

  if (length(population) != 17) {
    stop("population must be divided up into 17x 5-year age bands spanning 0 to 80+")
  }
  assert_int(seeding_cases)
  age_group_indices <- c(8, 9, 10, 11) # age_group indices corresponding to middle-aged travellers

  empty <- matrix(0, nrow = 17, ncol = 6)

  # distribute seeds either based on middle aged travellers randomly or ordinally by seeding_age_order
  if(is.null(seeding_age_order)) {
    raw_seeding_cases <- rep(0, length(population))
    raw_seeding_cases[age_group_indices] <- as.vector(stats::rmultinom(1, size = seeding_cases, prob = rep(0.25, 4)))
  } else {
    raw_seeding_cases <- rep(floor(seeding_cases/17), 17)
    for(s in seq_len(seeding_cases %% 17)) {
      raw_seeding_cases[seeding_age_order[s]] <- raw_seeding_cases[seeding_age_order[s]] + 1
    }
  }

  S = population - raw_seeding_cases
  S_0 = matrix(c(S, rep(0, 17*5)), nrow = 17, ncol = 6)
  E = raw_seeding_cases
  E_0 = matrix(c(E, rep(0, 17*5)), nrow = 17, ncol = 6)


  if(!all((rowSums(S_0) + rowSums(E_0)) == population)){
    stop("Row sums of init should be identical to population")
  }

  empty_inits <- c("D_0", "R_0", "ICase_0", "IMild_0", "IHosp_0")
  empties <- lapply(empty_inits, function(x){empty})
  names(empties) <- empty_inits

  init <- c(
    list(
      S_0 = S_0,
      E_0 = E_0
    ), empties)

  }

  return(init)
}
