#!/usr/bin/env Rscript
library(tidyverse)

plot <- ggplot(mpg, aes(displ, hwy, colour = class)) +
    geom_point()
mtcars |> write_tsv("cars.tsv")
ggsave("cars.png", plot = plot)
