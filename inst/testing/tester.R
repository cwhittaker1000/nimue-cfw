library(dplyr)
library(ggplot2)
library(nimue)

x <- run(countries = c("Nigeria", "United Kingdom"),
         seeding_cases = c(10, 0))

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
