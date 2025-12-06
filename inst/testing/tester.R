library(dplyr)
library(ggplot2)
library(nimue)

## check final size is right when specifying particular R0s - check vs squire

x <- run(countries = c("Nigeria", "United Kingdom"),
         R0 = list(c(2, 2), c(1.2, 1.2)),
         q = c(0, 0.9),
         pi_travel = matrix(data = c(0, 1, 1, 0), nrow = 2, ncol = 2, byrow = TRUE),
         tt_R0 = c(0, 220),
         seeding_cases = c(10, 10),
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
  countries   = c("France", "Germany"),
  age_breaks  = age_breaks,
  R0          = 3,
  tt_R0       = 0,
  time_period = 365,
  seeding_cases = 20,
  use_dde     = TRUE
)

# Quick sanity checks
sim$odin_parameters$N_age        # should be 5
sim$odin_parameters$age_breaks   # should be c(0, 20, 40, 60, 80, 100)
dim(sim$odin_parameters$mix_mat_set)  # c(5, 5, 2) for 2 locations

# For example, inspect the 5x5 mixing matrix for France:
sim$odin_parameters$mix_mat_set[ , , 1]

