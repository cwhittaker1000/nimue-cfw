#' Return the default probabilities for modelling defined in \code{squire}
#' For more info see \href{squire parameters vignette}{https://mrc-ide.github.io/squire/articles/parameters.html}
#' @return list of default probabilities
default_probs <- function() {
  prob_hosp <- c(
    0.000840764, 0.001182411, 0.001662887, 0.002338607, 0.003288907,
    0.004625365, 0.006504897, 0.009148183, 0.012865577, 0.018093546,
    0.025445917, 0.035785947, 0.050327683, 0.0707785, 0.099539573,
    0.1399878, 0.233470395)
  list(
    prob_hosp = prob_hosp,
    prob_death_hosp = c(rep(0.95, length(prob_hosp))), ### NOTE - ALWAYS BE CAREFUL TO SPECIFY YOUR OWN PROB_DEATH_HOSP
    p_dist = rep(1, length(prob_hosp)),
    rel_infectiousness = rep(1, 17),
    rel_infectiousness_vaccinated = rep(1,17)
  )
}
probs <- default_probs()

#' Return the default hospital durations for modelling defined in \code{squire}
#' For more info see \href{squire parameters vignette}{https://mrc-ide.github.io/squire/articles/parameters.html}
#' @return list of default durations
default_durations <- function() {
  squire_durs <- list(
    dur_IHosp = 8,
    dur_R = Inf,
    dur_E  = 4.6,
    dur_IMild = 2.1,
    dur_ICase = 4.5
  )
  return(squire_durs)
}

durs <- default_durations()

#' Return the default vaccine parameters for modelling
#' @return list of default vaccine parameters
default_vaccine_pars <- function() {
  list(dur_V = 365,
       vaccine_efficacy_infection = rep(0.95, 17),
       vaccine_efficacy_disease = rep(0.95, 17),
       max_vaccine = 1000,
       tt_vaccine = 0,
       dur_vaccine_delay = 14,
       vaccine_coverage_mat = matrix(0.8, ncol = 17, nrow = 1))
}

vaccine_pars <- default_vaccine_pars()

#' Vaccine parameters
#'
#' @details All durations are in days.
#'
#' @inheritParams run
## need to build in time-varying Rt
#
countries <- c("France", "United Kingdom")
dur_E <- durs$dur_E
dur_R <- durs$dur_R
dur_IHosp <- durs$dur_IHosp
dur_ICase <- durs$dur_ICase
dur_IMild <- durs$dur_IMild
dur_V <- Inf

prob_hosp = probs$prob_hosp
prob_death_hosp = probs$prob_death_hosp
p_dist = probs$p_dist
rel_infectiousness = probs$rel_infectiousness
rel_infectiousness_vaccinated = probs$rel_infectiousness_vaccinated

dur_V <- vaccine_pars$dur_V
vaccine_efficacy_infection <- vaccine_pars$vaccine_efficacy_infection
vaccine_efficacy_disease <- vaccine_pars$vaccine_efficacy_disease
dur_vaccine_delay <- vaccine_pars$dur_vaccine_delay
max_vaccine <- vaccine_pars$max_vaccine
tt_vaccine <- vaccine_pars$tt_vaccine
vaccine_coverage_mat <- vaccine_pars$vaccine_coverage_mat
q <- NULL
pi_travel <- NULL

parameters <- function(

  # Demography
  countries = NULL,
  age_breaks = NULL,

  # Transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  ## Mobility parameters
  q_flight = NULL,              # length(countries) vector
  pi_travel_flight = NULL,      # n_countries x n_countries matrix (home x destination)
  q_non_flight = NULL,              # length(countries) vector
  pi_travel_non_flight = NULL,      # n_countries x n_countries matrix (home x destination)

  # Initial state, duration, reps
  time_period = 365,
  seeding_cases,
  seeding_age_order = NULL,
  init = NULL,

  # Parameters
  # Probabilities
  prob_hosp = probs$prob_hosp,
  prob_death_hosp = probs$prob_death_hosp,
  p_dist = probs$p_dist,
  rel_infectiousness = probs$rel_infectiousness,
  rel_infectiousness_vaccinated = probs$rel_infectiousness_vaccinated,

  # Durations
  dur_E,
  dur_IMild,
  dur_ICase,
  dur_IHosp,

  # Vaccine
  dur_V,
  vaccine_efficacy_infection,
  vaccine_efficacy_disease,
  max_vaccine,
  max_vaccine_set = NULL,
  tt_vaccine,
  dur_vaccine_delay,
  vaccine_coverage_mat) {

  ## ---------------------------------------------------------------------------
  ## 1. Load raw population and contact matrices (canonical 17 age groups)
  ## ---------------------------------------------------------------------------
  n_countries <- length(countries)
  population_raw_list <- vector(mode = "list", length = n_countries)
  contact_matrix_list <- vector(mode = "list", length = n_countries)

  for (i in seq_len(n_countries)) {
    cpm <- parse_country_population_mixing_matrix(country = countries[i],
                                                  population = NULL,
                                                  contact_matrix_set = NULL)
    population_raw_list[[i]] <- cpm$population  # length 17
    contact_matrix_list[[i]] <- cpm$contact_matrix_set  # 16x16 contact matrix
  }

  ## ---------------------------------------------------------------------------
  ## 2. Age grouping: mapping + aggregation
  ## ---------------------------------------------------------------------------
  age_info <- compute_age_mapping(age_breaks)
  mapping      <- age_info$mapping       # length 17
  age_breaks   <- age_info$age_breaks    # possibly defaulted
  N_age        <- age_info$N_age         # new number of age groups

  # Aggregate population & contact matrices to chosen age groups
  population_list       <- vector("list", n_countries)
  baseline_matrix_list  <- vector("list", n_countries)  # for Rt->beta
  mixing_matrix_list    <- vector("list", n_countries)  # per-capita mixing
  for (i in seq_len(n_countries)) {

    # Get baseline population with 17 age-groups
    pop17 <- population_raw_list[[i]]

    # First: process original contact matrix to 17x17
    C17 <- process_contact_matrix_scaled_age(contact_matrix_list[[i]], pop17)

    # Aggregate to new age groups
    agg <- aggregate_population_and_matrix(pop17, C17, mapping, N_age)
    population_list[[i]]      <- agg$population  # length N_age
    baseline_matrix_list[[i]] <- agg$contact     # N_age x N_age

    # Mixing matrix used inside the model: divide by population (as before)
    mixing_matrix_list[[i]] <- div_pop(agg$contact, agg$population)
  }

  ## ---------------------------------------------------------------------------
  ## 3. Initialise state with aggregated populations
  ## ---------------------------------------------------------------------------
  mod_init <- init(population_list, seeding_cases, seeding_age_order, init)

  ## Build 3D mixing array: [age, age, location]
  matrices_set_array <- array(
    unlist(mixing_matrix_list, use.names = FALSE),
    dim = c(
      N_age,
      N_age,
      n_countries
    )
  )

  ## ---------------------------------------------------------------------------
  ## 4. Durations and time-varying Rt
  ## ---------------------------------------------------------------------------
  gamma_E        <- 1 / dur_E
  gamma_IMild    <- 1 / dur_IMild
  gamma_ICase    <- 1 / dur_ICase
  gamma_IHosp    <- 1 / dur_IHosp
  gamma_V        <- 1 * 1 / dur_V
  gamma_vaccine_delay <- 1 * 1 / dur_vaccine_delay

  n_tt <- length(tt_R0)
  if (n_tt < 1L) {
    stop("`tt_R0` must have length >= 1")
  }

  ## q: length must equal number of countries; default = no travel
  if (is.null(q)) {
    q <- rep(0, n_countries)
  } else {
    if (!is.numeric(q)) {
      stop("`q` must be numeric.")
    }
    if (length(q) != n_countries) {
      stop("`q` must have length equal to `length(countries)`.")
    }
  }

  ## ---------------------------------------------------------------------------
  ## 4a. Mobility: q_flight, q_non_flight, pi_travel_flight, pi_travel_non_flight
  ## ---------------------------------------------------------------------------

  ## q_flight: length must equal number of countries; default = no flight travel
  if (is.null(q_flight)) {
    q_flight <- rep(0, n_countries)
  } else {
    if (!is.numeric(q_flight)) {
      stop("`q_flight` must be numeric.")
    }
    if (length(q_flight) != n_countries) {
      stop("`q_flight` must have length equal to `length(countries)`.")
    }
  }

  ## q_non_flight: length must equal number of countries; default = no non-flight travel
  if (is.null(q_non_flight)) {
    q_non_flight <- rep(0, n_countries)
  } else {
    if (!is.numeric(q_non_flight)) {
      stop("`q_non_flight` must be numeric.")
    }
    if (length(q_non_flight) != n_countries) {
      stop("`q_non_flight` must have length equal to `length(countries)`.")
    }
  }

  ## pi_travel_flight: n_countries x n_countries matrix (home x destination)
  ## default = all zeros (no between-location mixing via flights)
  if (is.null(pi_travel_flight)) {
    pi_travel_flight <- matrix(0, nrow = n_countries, ncol = n_countries)
  } else {
    if (!is.matrix(pi_travel_flight)) {
      pi_travel_flight <- as.matrix(pi_travel_flight)
    }
    expected_dim <- c(n_countries, n_countries)
    if (!identical(dim(pi_travel_flight), expected_dim)) {
      stop("`pi_travel_flight` must be a ", n_countries, " x ", n_countries,
           " matrix (rows = home, cols = destination).")
    }
  }

  ## pi_travel_non_flight: n_countries x n_countries matrix (home x destination)
  ## default = all zeros (no between-location mixing via non-flight travel)
  if (is.null(pi_travel_non_flight)) {
    pi_travel_non_flight <- matrix(0, nrow = n_countries, ncol = n_countries)
  } else {
    if (!is.matrix(pi_travel_non_flight)) {
      pi_travel_non_flight <- as.matrix(pi_travel_non_flight)
    }
    expected_dim <- c(n_countries, n_countries)
    if (!identical(dim(pi_travel_non_flight), expected_dim)) {
      stop("`pi_travel_non_flight` must be a ", n_countries, " x ", n_countries,
           " matrix (rows = home, cols = destination).")
    }
  }

  ## Setting up the time-varying Rt for each country via the Rt matrix
  if (is.matrix(R0)) {
    Rt_mat <- R0
    expected_dim <- c(n_tt, n_countries)
    if (!identical(dim(Rt_mat), expected_dim)) {
      stop("When `R0` is a matrix, it must have nrow = length(tt_R0) and ncol = length(countries).")
    }

  } else if (is.list(R0)) {
    if (length(R0) != n_countries) {
      stop("When `R0` is a list, it must have length equal to `length(countries)`.")
    }
    Rt_mat <- matrix(NA_real_, nrow = n_tt, ncol = n_countries)

    for (i in seq_len(n_countries)) {
      R0_i <- R0[[i]]
      if (length(R0_i) != n_tt) {
        stop("Each element of `R0` list must have length equal to `length(tt_R0)`.")
      }
      Rt_mat[, i] <- R0_i
    }

  } else {

    ## 1) Single country: R0 and tt_R0 both allowed, must have the same length.
    if (n_countries == 1L) {
      if (length(R0) != n_tt) {
        stop("For a single country, `length(R0)` must equal `length(tt_R0)`.")
      }
      Rt_mat <- matrix(R0, ncol = 1L, nrow = n_tt)
    }  else {

    ## 2) Multiple countries: (i) one scalar per country (no time variation) -> length(R0) = n_countries, length(tt_R0) = 1
    ##                        (ii) one shared time-varying vector across countries -> length(R0) = length(tt_R0)
    ##                        (iii) full time x country matrix supplied as a numeric vector
      if (length(R0) == n_countries && n_tt == 1L) { # (i) different scalar per country, single timepoint
        Rt_mat <- matrix(R0, ncol = n_countries, nrow = 1L)
      } else if (length(R0) == n_tt) {  # (ii) same trajectory across countries
        Rt_mat <- matrix(rep(R0, each = n_countries), ncol = n_countries, nrow = n_tt, byrow = TRUE)
      } else if (length(R0) == n_tt * n_countries) {
        Rt_mat <- matrix(R0, nrow = n_tt, ncol = n_countries, byrow = TRUE)
      } else {
        stop(
          "For multiple countries, `R0` must either:\n", "  * have length `length(countries)` with `length(tt_R0) == 1`, or\n",
          "  * have length `length(tt_R0)` (shared trajectory across countries), or\n",
          "  * be supplied as a matrix (preferred) or vector with length `length(tt_R0) * length(countries)`."
        )
      }
    }
  }

  ## ---------------------------------------------------------------------------
  ## 4b. Vaccine capacity: build max_vaccine_set (time x location)
  ## ---------------------------------------------------------------------------
  n_tt_vac <- length(tt_vaccine)
  if (n_tt_vac < 1L) {
    stop("`tt_vaccine` must have length >= 1")
  }

  if (!is.null(max_vaccine_set)) {
    # User supplied full schedule (preferred for location-specific capacity)
    max_vaccine_set <- as.matrix(max_vaccine_set)
    expected_dim <- c(n_tt_vac, n_countries)
    if (!identical(dim(max_vaccine_set), expected_dim)) {
      stop("`max_vaccine_set` must be a matrix with nrow = length(tt_vaccine) and ncol = length(countries).")
    }
  } else {
    # Construct max_vaccine_set from simpler max_vaccine inputs
    if (!is.numeric(max_vaccine)) {
      stop("`max_vaccine` must be numeric when `max_vaccine_set` is NULL.")
    }

    if (length(max_vaccine) == 1L) {
      # Same capacity for all locations and all times
      max_vaccine_set <- matrix(max_vaccine, nrow = n_tt_vac, ncol = n_countries)
    } else if (length(max_vaccine) == n_tt_vac) {
      # Same time-varying schedule for all locations
      max_vaccine_set <- matrix(max_vaccine, nrow = n_tt_vac, ncol = n_countries)
    } else if (length(max_vaccine) == n_countries) {
      # Different constant capacity per location over time
      max_vaccine_set <- matrix(max_vaccine, nrow = n_tt_vac, ncol = n_countries, byrow = TRUE)
    } else {
      stop("When `max_vaccine_set` is NULL, `max_vaccine` must have length 1, length(tt_vaccine), or length(countries).")
    }
  }

  ## ---------------------------------------------------------------------------
  ## 5. Aggregate age-specific parameters to new age groups (if needed)
  ## -------------------------------------------------------------------

  ## Global population by canonical 17 age bands (sum across countries)
  ## Note: currently this adjustment for prob_hosp etc is being done on global pop for reweighting and calculation
  ##       of probs - maybe at some point we might want to make this country-specific
  pop_weights17 <- Reduce(`+`, population_raw_list)  # each element is length 17

  # prob_hosp
  if (length(prob_hosp) == length(mapping)) {
    prob_hosp <- aggregate_age_vector(prob_hosp, mapping, N_age,
                                      weights = pop_weights17,
                                      name = "prob_hosp")
  } else if (length(prob_hosp) != N_age) {
    stop("`prob_hosp` must have length 17 or length equal to the number of age groups implied by `age_breaks`.")
  }

  # prob_death_hosp
  if (length(prob_death_hosp) == length(mapping)) {
    prob_death_hosp <- aggregate_age_vector(prob_death_hosp, mapping, N_age,
                                            weights = pop_weights17,
                                            name = "prob_death_hosp")
  } else if (length(prob_death_hosp) != N_age) {
    stop("`prob_death_hosp` must have length 17 or length equal to the number of age groups implied by `age_breaks`.")
  }

  # rel_infectiousness
  if (length(rel_infectiousness) == length(mapping)) {
    rel_infectiousness <- aggregate_age_vector(rel_infectiousness, mapping, N_age,
                                               weights = pop_weights17,
                                               name = "rel_infectiousness")
  } else if (length(rel_infectiousness) != N_age) {
    stop("`rel_infectiousness` must have length 17 or length equal to the number of age groups implied by `age_breaks`.")
  }

  # rel_infectiousness
  if (length(rel_infectiousness_vaccinated) == length(mapping)) {
    rel_infectiousness_vaccinated <- aggregate_age_vector(rel_infectiousness_vaccinated, mapping, N_age,
                                               weights = pop_weights17,
                                               name = "rel_infectiousness_vaccinated")
  } else if (length(rel_infectiousness_vaccinated) != N_age) {
    stop("`rel_infectiousness_vaccinated` must have length 17 or length equal to the number of age groups implied by `age_breaks`.")
  }

  # p_dist
  if (length(p_dist) == length(mapping)) {
    p_dist <- aggregate_age_vector(p_dist, mapping, N_age,
                                   weights = pop_weights17,
                                   name = "p_dist")
  } else if (length(p_dist) != N_age) {
    stop("`p_dist` must have length 17 or length equal to the number of age groups implied by `age_breaks`.")
  }
  # vaccine_efficacy_infection
  if (length(vaccine_efficacy_infection) == length(mapping)) {
    vaccine_efficacy_infection <- aggregate_age_vector(vaccine_efficacy_infection, mapping, N_age,
                                                       weights = pop_weights17,
                                                       name = "vaccine_efficacy_infection")
  } else if (length(vaccine_efficacy_infection) != 1L && length(vaccine_efficacy_infection) != N_age) {
    stop("`vaccine_efficacy_infection` must have length 1, 17, or length equal to the number of age groups implied by `age_breaks`.")
  }

  # vaccine_efficacy_disease
  if (length(vaccine_efficacy_disease) == length(mapping)) {
    vaccine_efficacy_disease <- aggregate_age_vector(vaccine_efficacy_disease, mapping, N_age,
                                                     weights = pop_weights17,
                                                     name = "vaccine_efficacy_disease")
  } else if (length(vaccine_efficacy_disease) != 1L && length(vaccine_efficacy_disease) != N_age) {
    stop("`vaccine_efficacy_disease` must have length 1, 17, or length equal to the number of age groups implied by `age_breaks`.")
  }
  ## 5b. Check vaccine_coverage_mat has N_age columns
  if (is.null(dim(vaccine_coverage_mat))) {
    vaccine_coverage_mat <- matrix(vaccine_coverage_mat, nrow = 1)
  }
  if (ncol(vaccine_coverage_mat) != N_age) {
    stop("`vaccine_coverage_mat` must have N_age = ", N_age, " columns (one per age group). Got ", ncol(vaccine_coverage_mat), " columns instead.")
  }

  ## ---------------------------------------------------------------------------
  ## 6. Convert Rt to betas for each country using aggregated matrices
  ## ---------------------------------------------------------------------------
  beta_set <- matrix(NA_real_, nrow = n_tt, ncol = n_countries)
  for (i in seq_len(n_countries)) {

    # Country-specific baseline mixing matrix (already aggregated)
    baseline_matrix <- baseline_matrix_list[[i]]

    # Time-varying beta for this country
    beta_set[, i] <- vapply(
      seq_len(n_tt),
      function(j) {
        beta_est_infectiousness(
          dur_IMild = dur_IMild, dur_ICase = dur_ICase,
          prob_hosp = prob_hosp,
          mixing_matrix = baseline_matrix,
          rel_infectiousness = rel_infectiousness,
          R0 = Rt_mat[j, i]
        )
      },
      numeric(1)
    )
  }

  ## ---------------------------------------------------------------------------
  ## 7. p_dist matrix, vaccine parameters, VE formatting
  ## ---------------------------------------------------------------------------
  N_vaccine_states <- 4
  p_dist_mat <- matrix(rep(p_dist, N_vaccine_states),
                       nrow = N_age,
                       ncol = N_vaccine_states)
  p_dist_mat <- p_dist_mat / mean(p_dist_mat)

  p_dist <- p_dist_mat

  # Format vaccine-specific parameters
  gamma_vaccine <- c(0, gamma_vaccine_delay, gamma_V, 0)
  rel_infectiousness_vaccinated <- format_rel_inf_vacc_for_odin(rel_inf_vacc = rel_infectiousness_vaccinated, N_age = N_age)

  # Vaccine efficacy setup for the model
  vaccine_efficacy_infection_odin_array <- format_ve_i_for_odin(vaccine_efficacy_infection = vaccine_efficacy_infection, N_age = N_age)
  prob_hosp_odin_array <- format_ve_d_for_odin(vaccine_efficacy_disease = vaccine_efficacy_disease, prob_hosp = prob_hosp, N_age = N_age)
  prob_death_hosp_odin_array <- format_ve_d_for_odin(vaccine_efficacy_disease = 0, prob_hosp = prob_death_hosp, N_age = N_age)

  ## ---------------------------------------------------------------------------
  ## 8. Collate parameters
  ## ---------------------------------------------------------------------------

  # Collate Parameters Into List
  pars <- c(mod_init,
            list(N_age = N_age,
                 N_locations = n_countries,
                 age_breaks = age_breaks,
                 q_flight = q_flight,
                 pi_travel_flight = pi_travel_flight,
                 q_non_flight = q_non_flight,
                 pi_travel_non_flight = pi_travel_non_flight,
                 gamma_E = gamma_E,
                 gamma_IMild = gamma_IMild,
                 gamma_ICase = gamma_ICase,
                 gamma_IHosp = gamma_IHosp,
                 prob_hosp = prob_hosp_odin_array,
                 prob_death_hosp = prob_death_hosp_odin_array,
                 rel_infectiousness = rel_infectiousness,
                 rel_infectiousness_vaccinated = rel_infectiousness_vaccinated,
                 p_dist = p_dist,
                 mix_mat_set = matrices_set_array,
                 tt_beta = tt_R0,
                 beta_set = beta_set,
                 population_list = population_list,
                 max_vaccine       = max_vaccine,
                 max_vaccine_set   = max_vaccine_set,
                 tt_vaccine        = tt_vaccine,
                 vaccine_efficacy_infection = vaccine_efficacy_infection_odin_array,
                 vaccine_coverage_mat = vaccine_coverage_mat,
                 N_vaccine = 4,
                 N_prioritisation_steps = nrow(vaccine_coverage_mat),
                 gamma_vaccine = gamma_vaccine))

  class(pars) <- c("vaccine_parameters", "nimue_parameters")

  return(pars)
}

#' Estimate beta parameter for explicit model
#'
#' @param dur_IMild Duration of mild infectiousness (days)
#' @param dur_ICase Delay between symptom onset and requiring hospitalisation (days)
#' @param prob_hosp Probability of hospitilisation by ages
#' @param rel_infectiousness Relative infectiousness of age categories relative
#'   to maximum infectiousness age category
#' @param mixing_matrix Mixing matrix
#' @param R0 Basic reproduction number
#'
#' @return Beta parameter
#' @export
#'
# #' @examples
beta_est_infectiousness <- function(dur_IMild,
                                    dur_ICase,
                                    prob_hosp,
                                    rel_infectiousness,
                                    mixing_matrix,
                                    R0) {

  # assertions
  assert_single_pos(dur_ICase, zero_allowed = FALSE)
  assert_single_pos(dur_IMild, zero_allowed = FALSE)
  assert_numeric(prob_hosp)
  assert_numeric(rel_infectiousness)
  assert_same_length(prob_hosp, rel_infectiousness)
  assert_numeric(mixing_matrix)
  assert_square_matrix(mixing_matrix)
  assert_same_length(mixing_matrix[,1], prob_hosp)
  assert_pos(R0, zero_allowed = FALSE)

  if(sum(is.na(prob_hosp)) > 0) {
    stop("prob_hosp must not contain NAs")
  }

  if(sum(is.na(rel_infectiousness)) > 0) {
    stop("rel_infectiousness must not contain NAs")
  }

  if(sum(is.na(mixing_matrix)) > 0) {
    stop("mixing_matrix must not contain NAs")
  }

  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  adjusted_eigen <- Re(eigen(mixing_matrix*relative_R0_by_age*rel_infectiousness)$values[1])
  R0 / adjusted_eigen

}

#' @noRd
format_rel_inf_vacc_for_odin <- function(rel_inf_vacc, N_age) {

  # rel_inf_vacc must be length 1 or N_age
  if (length(rel_inf_vacc) == 1L) {
    rel_inf_vacc <- rep(rel_inf_vacc, N_age)
  }

  if (length(rel_inf_vacc) != N_age) {
    stop("`rel_inf_vacc` must be length 1 or length N_age = ", N_age)
  }

  # age x vaccine-class matrix (N_age x 4)
  # Classes: [1:2] = 1,
  #          [3]   = rel_inf_vacc,
  #          [4]   = 1
  mat <- matrix(
    c(
      rep(1, N_age * 2L),
      rel_inf_vacc,
      rep(1, N_age)
    ),
    nrow = N_age,
    ncol = 4
  )

  return(mat)
}


#' @noRd
format_ve_i_for_odin <- function(vaccine_efficacy_infection, N_age) {

  # vaccine_efficacy_infection must be length 1 or N_age
  if (length(vaccine_efficacy_infection) == 1L) {
    vaccine_efficacy_infection <- rep(vaccine_efficacy_infection, N_age)
  }

  if (length(vaccine_efficacy_infection) != N_age) {
    stop("Parameter `vaccine_efficacy_infection` must be length 1 or length N_age = ", N_age)
  }

  # Convert to relative susceptibility (1 - VE)
  ve_i <- 1 - vaccine_efficacy_infection

  # age x vaccine-class matrix (N_age x 4)
  # Classes: [1:2] = 1,
  #          [3]   = ve_i,
  #          [4]   = 1 (waned)
  vaccine_efficacy_infection_mat <- matrix(
    c(
      rep(1, N_age * 2L),
      ve_i,
      rep(1, N_age)
    ),
    nrow = N_age,
    ncol = 4
  )

  return(vaccine_efficacy_infection_mat)
}

#' @noRd
format_ve_d_for_odin <- function(vaccine_efficacy_disease,
                                 prob_hosp,
                                 N_age) {

  # vaccine_efficacy_disease must be length 1 or N_age
  if (length(vaccine_efficacy_disease) == 1L) {
    vaccine_efficacy_disease <- rep(vaccine_efficacy_disease, N_age)
  }

  if (length(vaccine_efficacy_disease) != N_age) {
    stop("Parameter `vaccine_efficacy_disease` must be length 1 or length N_age = ", N_age)
  }

  if (length(prob_hosp) != N_age) {
    stop("Parameter `prob_hosp` must have length N_age = ", N_age)
  }

  # Per-age hospitalisation prob (or death prob) for vaccinated classes
  prob_hosp_vaccine <- (1 - vaccine_efficacy_disease) * prob_hosp

  # age x vaccine-class matrix (N_age x 4)
  # Classes:
  #   [1:2] = prob_hosp (unvaccinated and vaccinated but no protection),
  #   [3]   = prob_hosp_vaccine (vaccinated and protected),
  #   [4]   = prob_hosp (waned)
  prob_hosp_mat <- matrix(
    c(
      prob_hosp,
      prob_hosp,
      prob_hosp_vaccine,
      prob_hosp
    ),
    nrow = N_age,
    ncol = 4
  )

  return(prob_hosp_mat)
}


