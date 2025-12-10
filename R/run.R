#' Run the vaccine model
#'
#' @param countries Character vector of locations to include in the simulation.
#'   Populations and baseline contact matrices are sourced for each country.
#' @param age_breaks Numeric vector specifying the upper bound of each age
#'   group. If \code{NULL}, the default 17-group structure is used.
#' @param R0 Basic reproduction number for each time period specified in
#'   \code{tt_R0}. Default is 3.
#' @param tt_R0 Time points (in days) at which \code{R0} changes. Defaults to 0
#'   (no change).
#' @param beta_set Optional alternative parameterisation via beta rather than
#'   \code{R0}. If \code{NULL}, beta is derived from \code{R0}.
#' @param q_flight Daily probability of travelling by flight for each country.
#' @param pi_travel_flight Matrix of destination probabilities for flight travel
#'   (home country in rows, destination in columns).
#' @param q_non_flight Daily probability of travelling by non-flight transport
#'   for each country.
#' @param pi_travel_non_flight Matrix of destination probabilities for
#'   non-flight travel (home country in rows, destination in columns).
#' @param time_period Length of the simulation in days. Default is 365.
#' @param replicates Number of stochastic replicates. Currently fixed to 1 in
#'   the underlying model.
#' @param seed Random seed used for simulations. Default draws from
#'   \code{runif(1, 0, 1e8)}.
#' @param prob_hosp Probability of hospitalisation by age group.
#' @param prob_death_hosp Probability of death among hospitalised cases by age
#'   group.
#' @param rel_infectiousness Relative infectiousness per age category.
#' @param rel_infectiousness_vaccinated Relative infectiousness for vaccinated
#'   individuals compared to unvaccinated individuals by age group.
#' @param dur_E Mean duration of incubation period (days).
#' @param dur_IMild Mean duration of mild infection (days).
#' @param dur_ICase Mean duration from symptom onset to hospital admission
#'   (days).
#' @param dur_IHosp Mean duration of hospitalisation (days).
#' @param dur_V Mean duration of vaccine-derived immunity (days).
#' @param vaccine_efficacy_infection Vaccine efficacy against infection, either
#'   a single value or a vector by age group.
#' @param vaccine_efficacy_disease Vaccine efficacy against severe disease,
#'   either a single value or a vector by age group.
#' @param max_vaccine Maximum number of vaccine doses available per day.
#' @param max_vaccine_set Optional matrix specifying vaccine capacity over time
#'   and by location (rows are time, columns are countries).
#' @param tt_vaccine Time points (in days) when vaccine capacity changes.
#' @param dur_vaccine_delay Delay (in days) before vaccine-derived immunity
#'   becomes effective.
#' @param vaccine_coverage_mat Matrix of target vaccine coverage by time (rows)
#'   and age group (columns).
#' @param seeding_cases Initial number of cases used to seed the epidemic.
#' @param seeding_age_order Optional integer vector controlling the order in
#'   which seeds are allocated to age groups. If \code{NULL}, seeds are
#'   distributed randomly within working ages.
#' @param init Optional list of initial conditions for the simulation.
#' @param use_dde Logical indicating whether to use the \code{dde}
#'   implementation (default) or \code{deSolve}.
#' @param ... Additional arguments passed to the underlying odin model
#'   \code{run} method.
#'
#' @return Simulation output
#' @export
run <- function(

  # demography
  countries = NULL,
  age_breaks = NULL,   # user-specified age group breaks [0, a), [a, b), ..., [e, 100)

  # transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  ## Mobility parameters
  q_flight = NULL,              # length(countries) vector
  pi_travel_flight = NULL,      # n_countries x n_countries matrix (home x destination)
  q_non_flight = NULL,              # length(countries) vector
  pi_travel_non_flight = NULL,      # n_countries x n_countries matrix (home x destination)

  # initial state, duration, reps
  time_period = 365,
  replicates = 1,
  seed = stats::runif(1, 0, 100000000),

  # parameters
  # probabilities
  prob_hosp = probs$prob_hosp,
  prob_death_hosp = probs$prob_death_hosp,

  # onward infectiousness
  rel_infectiousness = probs$rel_infectiousness,
  rel_infectiousness_vaccinated = probs$rel_infectiousness_vaccinated,

  # durations
  dur_E  = durs$dur_E,
  dur_IMild = durs$dur_IMild,
  dur_ICase = durs$dur_ICase,
  dur_IHosp = durs$dur_IHosp,

  # vaccine
  dur_V = vaccine_pars$dur_V,
  vaccine_efficacy_infection = vaccine_pars$vaccine_efficacy_infection,
  vaccine_efficacy_disease = vaccine_pars$vaccine_efficacy_disease,
  max_vaccine = vaccine_pars$max_vaccine,
  max_vaccine_set = NULL,   # time x location matrix
  tt_vaccine = vaccine_pars$tt_vaccine,
  dur_vaccine_delay = vaccine_pars$dur_vaccine_delay,
  vaccine_coverage_mat = vaccine_pars$vaccine_coverage_mat,

  # misc
  seeding_cases = 20,
  seeding_age_order = NULL,
  init = NULL,
  use_dde = TRUE,
  ...
) {

  # Grab function arguments
  args <- as.list(environment())
  set.seed(seed)


  ## Warning
  if (is.null(prob_death_hosp)) {
    warning(
      "`prob_death_hosp` has not been specified; using an internal, slightly weird default instead. ",
      "Please set `prob_death_hosp` explicitly for any serious analysis.",
      call. = FALSE
    )
  }

  # create parameter list
  pars <- parameters(countries = countries,
                     age_breaks = age_breaks,
                     R0 = R0,
                     tt_R0 = tt_R0 ,
                     beta_set = beta_set,
                     q_flight = q_flight,
                     pi_travel_flight = pi_travel_flight,
                     q_non_flight = q_non_flight,
                     pi_travel_non_flight = pi_travel_non_flight,
                     time_period = time_period,
                     seeding_cases = seeding_cases,
                     seeding_age_order = seeding_age_order,
                     prob_hosp = prob_hosp,
                     prob_death_hosp = prob_death_hosp,
                     rel_infectiousness = rel_infectiousness,
                     rel_infectiousness_vaccinated = rel_infectiousness_vaccinated,
                     dur_E = dur_E,
                     dur_IMild = dur_IMild,
                     dur_ICase = dur_ICase,
                     dur_IHosp = dur_IHosp,
                     dur_V = dur_V,
                     vaccine_efficacy_infection = vaccine_efficacy_infection,
                     vaccine_efficacy_disease = vaccine_efficacy_disease,
                     max_vaccine = max_vaccine,
                     max_vaccine_set  = max_vaccine_set,
                     tt_vaccine = tt_vaccine,
                     dur_vaccine_delay = dur_vaccine_delay,
                     vaccine_coverage_mat = vaccine_coverage_mat,
                     init = init)

  # Set model type
  replicates <- 1
  mod_gen = vaccine_simplified_multiloc

  # Running the Model
  mod <- mod_gen$new(user = pars, unused_user_action = "ignore", use_dde = use_dde)

  # Daily output by default
  t <- round(seq(from = 1, to = time_period))

  results <- mod$run(t, ...)

  # coerce to array
  results <- array(results, dim = c(dim(results), 1), dimnames = dimnames(results))

  # Summarise inputs
  parameters <- args
  parameters$age_breaks <- age_breaks

  out <- list(output = results, parameters = parameters, model = mod, odin_parameters = pars)
  out <- structure(out, class = "nimue_simulation")
  return(out)

}

