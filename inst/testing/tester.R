library(dplyr)
library(ggplot2)
library(nimue)

## check final size is right when specifying particular R0s - check vs squire
q_flight      <- c(0, 0)
q_non_flight  <- c(0, 0)
pi_travel_flight     <- matrix(c(0,0, 0,0), nrow = 2, byrow = TRUE)
pi_travel_non_flight <- matrix(c(0,0, 0,0), nrow = 2, byrow = TRUE)

x <- run(countries = c("Nigeria", "United Kingdom"),
         R0 = list(c(2, 2), c(1.2, 1.2)),
         q_flight             = q_flight,
         q_non_flight         = q_non_flight,
         pi_travel_flight     = pi_travel_flight,
         pi_travel_non_flight = pi_travel_non_flight,
         tt_R0 = c(0, 220),
         seeding_cases = c(10, 5),
         time_period = 1000)

squire:::get_population("France")

output1 <- format_multiloc(x)

ggplot(output1, aes(x = t, y = value, col = compartment))  +
  geom_line(size = 1) +
  ylab("Infections") +
  facet_wrap(location_index ~ ., scales = "free_y") +
  xlab("Time") +
  theme_bw()

output1 %>%
  filter(compartment == "D") %>%
  group_by(location_index) %>%
  summarise(D = max(value))

output1 %>%
  filter(compartment == "S") %>%
  group_by(location_index) %>%
  summarise(N = max(value))


1636420 / 206139567
1702891 / 67885984


library(nimue)   # or your package name if this is a fork

# Define 5 age groups:
# [0,20), [20,40), [40,60), [60,80), [80,100)
age_breaks <- c(0, 20, 40, 60, 80, 100)

sim <- run(
  countries   = c("France", "Nigeria"),
  age_breaks  = age_breaks,
  R0          = 3,
  tt_R0       = 0,
  time_period = 365,
  seeding_cases = 20,
  use_dde     = TRUE,
  vaccine_coverage_mat = matrix(0.8, ncol = 5, nrow = 1)
)

# Quick sanity checks
sim$odin_parameters$N_age        # should be 5
sim$odin_parameters$age_breaks   # should be c(0, 20, 40, 60, 80, 100)
dim(sim$odin_parameters$mix_mat_set)  # c(5, 5, 2) for 2 locations

# For example, inspect the 5x5 mixing matrix for France:
sim$odin_parameters$mix_mat_set[ , , 1]

output1 <- format_multiloc(sim)

ggplot(output1, aes(x = t, y = value, col = compartment))  +
  geom_line(size = 1) +
  ylab("Infections") +
  facet_wrap(location_index ~ ., scales = "free_y") +
  xlab("Time") +
  theme_bw()

### testing vacciantion
library(nimue)   # or your package name if this is a fork

# Define 5 age groups:
# [0,20), [20,40), [40,60), [60,80), [80,100)
age_breaks <- c(0, 20, 40, 60, 80, 100)

sim <- run(
  countries   = c("United Kingdom"),
  age_breaks  = age_breaks,
  R0          = 1.5,
  max_vaccine = seq(0, 2500000, length.out = 100),
  tt_vaccine = seq(0, 10, length.out = 100),
  vaccine_efficacy_disease = rep(0.95, 17),
  vaccine_efficacy_infection = rep(0, 17),
  tt_R0       = 0,
  time_period = 365,
  seeding_cases = 20,
  dur_vaccine_delay = 2,
  dur_V = 1e5,
  use_dde     = TRUE,
  vaccine_coverage_mat = matrix(0.95, ncol = 5, nrow = 1)
)

output1 <- format_multiloc(sim) %>%
  filter(compartment == "D")

ggplot(output1, aes(x = t, y = value, col = compartment))  +
  geom_line(size = 1) +
  ylab("Infections") +
  facet_wrap(location_index ~ ., scales = "free_y") +
  xlab("Time") +
  theme_bw()

output1 %>%
  filter(compartment == "D") %>%
  group_by(location_index) %>%
  summarise(D = max(value))

78434
830670

1 - (78434 / 830667)

##


# Vaccine capacity time points
tt_vaccine <- c(0, 30, 60)  # days

# Location-specific max_vaccine per time:
# rows = times in tt_vaccine, cols = locations (France, UK)
max_vaccine_set <- rbind(
  c(10^7,  10^6),  # at t = 0
  c(10^7,  10^6),  # at t = 30
  c(10^7, 10^6)   # at t = 60
)

sim <- run(
  countries   = c("United Kingdom", "Nigeria"),
  age_breaks  = age_breaks,
  R0          = 3,
  vaccine_efficacy_disease = rep(0.95, 17),
  vaccine_efficacy_infection = rep(0, 17),
  tt_R0       = 0,
  time_period = 365,
  seeding_cases = 20,
  dur_vaccine_delay = 2,
  dur_V = 1e5,
  max_vaccine_set        = max_vaccine_set,
  tt_vaccine             = tt_vaccine,
  use_dde     = TRUE,
  vaccine_coverage_mat = matrix(0.95, ncol = 5, nrow = 1)
)

output1 <- format_multiloc(sim)
ggplot(output1, aes(x = t, y = value, col = compartment))  +
  geom_line(size = 1) +
  ylab("Infections") +
  facet_wrap(location_index ~ ., scales = "free_y") +
  xlab("Time") +
  theme_bw()

output1 %>%
  filter(compartment == "D") %>%
  group_by(location_index) %>%
  summarise(D = max(value))

167435 / 1741582
1730267 / 1730267

0.95 * 0.95
