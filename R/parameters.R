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

parameters <- function(

  # Demography
  countries = NULL,

  # Transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  ## Mobility parameters
  q = NULL,              # length(countries) vector
  pi_travel = NULL,      # n_countries x n_countries matrix (home x destination)

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
  tt_vaccine,
  dur_vaccine_delay,
  vaccine_coverage_mat) {

  # Handle country population args
  population_list <- vector(mode = "list", length = length(countries))
  contact_matrix_list <- vector(mode = "list", length = length(countries))
  for (i in 1:length(countries)) {
    cpm <- parse_country_population_mixing_matrix(country = countries[i],
                                                  population = NULL,
                                                  contact_matrix_set = NULL)
    population_list[[i]] <- cpm$population
    contact_matrix_list[[i]] <- cpm$contact_matrix_set  # Keep as list
  }

  # Initial state and matrix formatting
  # ----------------------------------------------------------------------------

  # Initialise initial conditions
  mod_init <- init(population_list, seeding_cases, seeding_age_order, init)

  # Convert contact matrices to input matrices
  stopifnot(length(population_list) == length(contact_matrix_list))
  matrices_set_list <- Map(
    function(contact_set, pop) {
      matrix_set_explicit(list(contact_set), pop)
    },
    contact_matrix_list,
    population_list
  )
  matrices_set_array <- array(
    unlist(matrices_set_list, use.names = FALSE),
    dim = c(
      nrow(matrices_set_list[[1]]),
      ncol(matrices_set_list[[1]]),
      length(matrices_set_list)
    )
  )

  # If a vector is put in for matrix targeting
  ### NOTE: COME BACK TO THIS!!!
  # if(is.vector(vaccine_coverage_mat)){
  #   vaccine_coverage_mat <- matrix(vaccine_coverage_mat, ncol = 17)
  # }

  # Convert and Generate Parameters As Required
  # ----------------------------------------------------------------------------

  # durations
  gamma_E = 1/dur_E
  gamma_IMild = 1/dur_IMild
  gamma_ICase = 1/dur_ICase
  gamma_IHosp = 1/dur_IHosp
  gamma_V <- 2 * 1/dur_V
  gamma_vaccine_delay <- 2 * 1 / dur_vaccine_delay

  ## Setting up time-varying Rt
  n_countries <- length(countries)
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

  ## pi_travel: n_countries x n_countries matrix (home x destination)
  ## default = all zeros (no between-location mixing while away)
  if (is.null(pi_travel)) {
    pi_travel <- matrix(0, nrow = n_countries, ncol = n_countries)
  } else {
    if (!is.matrix(pi_travel)) {
      pi_travel <- as.matrix(pi_travel)
    }
    expected_dim <- c(n_countries, n_countries)
    if (!identical(dim(pi_travel), expected_dim)) {
      stop("`pi_travel` must be a ", n_countries, " x ", n_countries,
           " matrix (rows = home, cols = destination).")
    }
  }

  ## Setting up the time-varying Rt for each country via the Rt matrix
  if (is.list(R0)) {
    if (length(R0) != n_countries) {
      stop("When `R0` is a list, it must have length equal to `length(countries)`.")
    }
    Rt_mat <- matrix(NA_real_, nrow = n_countries, ncol = n_tt)

    for (i in seq_len(n_countries)) {
      R0_i <- R0[[i]]
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
      if (length(R0) == n_countries && n_tt == 1L) { # (i) different scalar per country, single timepoint
        Rt_mat <- matrix(R0, ncol = n_countries, nrow = 1L)
      } else if (length(R0) == n_tt) {  # (ii) same trajectory across countries
        Rt_mat <- matrix(rep(R0, each = n_countries), ncol = n_countries, nrow = n_tt, byrow = TRUE)
      } else {
        stop(
          "For multiple countries, `R0` must either:\n", "  * have length `length(countries)` with `length(tt_R0) == 1`, or\n",
          "  * have length `length(tt_R0)` (shared trajectory across countries)."
        )
      }
    }
  }

  ## Converting Rt to betas for each country
  # beta_set is i timepoints and j countries
  beta_set <- matrix(NA_real_, nrow = n_tt, ncol = n_countries)
  for (i in seq_len(n_countries)) {

    # Country-specific baseline mixing matrix
    baseline_matrix <- process_contact_matrix_scaled_age(contact_matrix_list[[i]], population_list[[i]] )

    # Time-varying beta for this country
    beta_set[, i] <- vapply(seq_len(n_tt), function(j) { beta_est_infectiousness(
      dur_IMild = dur_IMild, dur_ICase = dur_ICase,
      prob_hosp = prob_hosp,
      mixing_matrix = baseline_matrix,
      rel_infectiousness = rel_infectiousness,
      R0 = Rt_mat[j, i])
      }, numeric(1))
  }

  # normalise to sum to 1
  p_dist <- matrix(rep(p_dist, 6), nrow = 17, ncol = 6)
  p_dist <- p_dist/mean(p_dist)

  # Format vaccine-specific parameters
  gamma_vaccine <- c(0, gamma_vaccine_delay, gamma_vaccine_delay, gamma_V, gamma_V, 0)
  rel_infectiousness_vaccinated <- format_rel_inf_vacc_for_odin(rel_infectiousness_vaccinated)

  # Vaccine efficacy setup for the model
  # First the vaccine efficacy infection
  vaccine_efficacy_infection_odin_array <- format_ve_i_for_odin(vaccine_efficacy_infection = vaccine_efficacy_infection)
  # Second the vaccine efficacy disease affecting prob_hosp
  prob_hosp_odin_array <- format_ve_d_for_odin(vaccine_efficacy_disease = vaccine_efficacy_disease, prob_hosp = prob_hosp)
  prob_death_hosp_odin_array <- format_ve_d_for_odin(vaccine_efficacy_disease = 0, prob_hosp = prob_death_hosp)

  # Collate Parameters Into List
  pars <- c(mod_init,
            list(N_age = length(population_list[[1]]),
                 N_locations = length(population_list),
                 q = q,
                 pi_travel = pi_travel,
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
                 max_vaccine = max_vaccine,
                 tt_vaccine = tt_vaccine,
                 vaccine_efficacy_infection = vaccine_efficacy_infection_odin_array,
                 vaccine_coverage_mat = vaccine_coverage_mat,
                 N_vaccine = 6,
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
format_rel_inf_vacc_for_odin <- function(rel_inf_vacc) {

  if(length(rel_inf_vacc) == 1){
    rel_inf_vacc <- rep(rel_inf_vacc, 17)
  }

  return(matrix(c(rep(1, 17 * 3),
                  rel_inf_vacc, rel_inf_vacc,
                  rep(1, 17)), nrow = 17, ncol = 6))

}

#' @noRd
format_ve_i_for_odin <- function(vaccine_efficacy_infection) {

  # vaccine_efficacy_infection must be length 1 or 17
  if (length(vaccine_efficacy_infection) == 1L) {
    vaccine_efficacy_infection <- rep(vaccine_efficacy_infection, 17)
  }

  if (length(vaccine_efficacy_infection) != 17L) {
    stop("Parameter `vaccine_efficacy_infection` must be length 1 or length 17")
  }

  # Convert to relative susceptibility (1 - VE)
  ve_i <- 1 - vaccine_efficacy_infection

  # age x vaccine-class matrix (17 x 6)
  # Classes: [1:3] = 1, [4:5] = ve_i, [6] = 1 (as in original nimue logic)
  vaccine_efficacy_infection_mat <- matrix(
    c(
      rep(1, 17 * 3),
      ve_i, ve_i,
      rep(1, 17)
    ),
    nrow = 17, ncol = 6
  )

  return(vaccine_efficacy_infection_mat)
}

#' @noRd
format_ve_d_for_odin <- function(vaccine_efficacy_disease,
                                 prob_hosp) {

  # vaccine_efficacy_disease must be length 1 or 17
  if (length(vaccine_efficacy_disease) == 1L) {
    vaccine_efficacy_disease <- rep(vaccine_efficacy_disease, 17)
  }

  if (length(vaccine_efficacy_disease) != 17L) {
    stop("Parameter `vaccine_efficacy_disease` must be length 1 or length 17")
  }

  if (length(prob_hosp) != 17L) {
    stop("Parameter `prob_hosp` must be length 17")
  }

  # Per-age hospitalisation prob for vaccinated classes
  prob_hosp_vaccine <- (1 - vaccine_efficacy_disease) * prob_hosp

  # age x vaccine-class matrix (17 x 6)
  # Classes: [1:3] = prob_hosp, [4:5] = prob_hosp_vaccine, [6] = prob_hosp
  prob_hosp_mat <- matrix(
    c(
      prob_hosp, prob_hosp, prob_hosp,
      prob_hosp_vaccine, prob_hosp_vaccine,
      prob_hosp
    ),
    nrow = 17, ncol = 6
  )

  return(prob_hosp_mat)
}

