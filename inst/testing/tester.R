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
