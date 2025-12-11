# devtools::load_all()  # from the package root, or library(yourpkg)

countries <- c("France", "United Kingdom")
library(nimue)

squire:::get_lmic_countries()

tt_vaccine <- c(0, 30, 60)
max_vaccine_set <- rbind(
  c(1000,  500),
  c(2000,  800),
  c(3000, 20)
)
squire::income_group$country
# simple, no mobility (defaults q_* = 0, pi_* = 0)
sim <- run(
  countries        = countries,
  R0               = c(3, 3),
  tt_R0            = 0,
  time_period      = 60,
  seeding_cases    = 20,
  max_vaccine      = 1000,
  max_vaccine_set  = max_vaccine_set,
  tt_vaccine       = tt_vaccine,
  use_dde          = TRUE
)

str(sim$output)
sim$odin_parameters$N_locations      # should be 2
sim$odin_parameters$max_vaccine_set  # should be 3 x 2 matrix

# trivial symmetric kernels & small away fractions
q_flight      <- c(0.1, 0.2)
q_non_flight  <- c(0.1, 0.1)
pi_travel_flight     <- matrix(c(0,1, 1,0), nrow = 2, byrow = TRUE)
pi_travel_non_flight <- matrix(c(0,1, 1,0), nrow = 2, byrow = TRUE)

sim2 <- run(
  countries            = countries,
  R0                   = c(3, 3),
  tt_R0                = 0,
  time_period          = 60,
  seeding_cases        = 20,
  max_vaccine_set      = max_vaccine_set,
  tt_vaccine           = tt_vaccine,
  q_flight             = q_flight,
  q_non_flight         = q_non_flight,
  pi_travel_flight     = pi_travel_flight,
  pi_travel_non_flight = pi_travel_non_flight,
  use_dde              = TRUE
)

tt_vaccine <- c(0, 30, 60)
max_vaccine_set <- rbind(
  c(1000,  500),
  c(2000,  800),
  c(3000, 20)
)

# simple, no mobility (defaults q_* = 0, pi_* = 0)
sim <- run(
  countries        = squire::income_group$country,
  R0               = 3,
  tt_R0            = 0,
  time_period      = 60,
  seeding_cases    = 20,
  max_vaccine      = 0,
  max_vaccine_set  = NULL,
  tt_vaccine       = 0,
  use_dde          = TRUE
)
