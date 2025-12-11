library(squire)   # for income_group
library(tictoc)
library(dplyr)
library(ggplot2)
# library(nimuecfw)  # or whatever your package name is

# All countries in the income_group table
countries <- squire::income_group$country
n_loc <- length(countries)  # should be 201

# Away fractions: 2% via flights, 2% via non-flight (gravity)
q_flight     <- rep(0.02, n_loc)
q_non_flight <- rep(0.02, n_loc)
# Total away fraction per location = 0.04, so 96% time at home

# Flight travel kernel: uniform across all *other* locations
pi_travel_flight <- matrix(1 / (n_loc - 1), nrow = n_loc, ncol = n_loc)
diag(pi_travel_flight) <- 0  # no self-travel; each row still sums to 1

# Non-flight (gravity-style) travel kernel: here, same structure as flights,
# but you can later replace this with a distance-based matrix
pi_travel_non_flight <- matrix(1 / (n_loc - 1), nrow = n_loc, ncol = n_loc)
diag(pi_travel_non_flight) <- 0

# Simple test run: no vaccination (capacity = 0)
tic()
age_breaks <- c(0, 20, 40, 60, 80, 100)
sim <- run(
  countries              = countries,
  age_breaks             = age_breaks,            # use default 17 age groups
  R0                     = 3,
  tt_R0                  = 0,
  time_period            = 60,
  seeding_cases          = 20,
  max_vaccine            = 0,              # no vaccines at any time/location
  max_vaccine_set        = NULL,           # let parameters() build a 0 matrix
  tt_vaccine             = 0,              # single timepoint is fine
  q_flight               = q_flight,
  pi_travel_flight       = pi_travel_flight,
  q_non_flight           = q_non_flight,
  pi_travel_non_flight   = pi_travel_non_flight,
  use_dde                = TRUE,
  vaccine_coverage_mat = NULL,
)
toc()
output1 <- format_multiloc(sim) %>%
  filter(compartment == "IMild")

ggplot(output1, aes(x = t, y = value, group = location_index))  +
  geom_line() +
  ylab("Infections") +
  xlab("Time") +
  theme_bw()



# Quick sanity checks
str(sim$output)
sim$odin_parameters$N_locations         # should be 201
sim$odin_parameters$q_flight[1:5]
rowSums(sim$odin_parameters$pi_travel_flight)[1:5]
