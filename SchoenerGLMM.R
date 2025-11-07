#This script analyses and plot the niche overlaps (Schoener index values) within flowering community and between flowering community and pollinators
 
#As the response values are between 0 and 1, a beta regression is used (with 0 and 1 inflation if needed) 

library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(purrr)      # for looping over models
library(tibble)     # for tibble/data frame formatting

#Schoener brms

#All the flowering species
AllFlowers <- read.csv(file = file.path("./Data", "Schoener.csv"))
M3.1 <- brm(Schoener ~ Selevation + (1|Meadow/Site),
            data = AllFlowers,
            family = 'zero_one_inflated_beta',
            iter = 1000,
            chains = 2,
            cores = 2)

#All flowering species and insects total
Insects <- read.csv(file = file.path("./Data", "Schoener Flower Insect.csv"))
M3.2 <- brm(Schoener ~ Elevation + (1|Meadow/Site),
            data = Insects,
            family = 'beta',
            iter = 1000,
            chains = 2,
            cores = 2)

#All flowering species and Coleoptera
Coleoptera <- read.csv(file = file.path("./Data", "Schoener Flower Coleoptera.csv"))
M3.3 <- brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = Coleoptera,
             family = 'zero_inflated_beta',
             iter = 1000,
             chains = 2,
             cores = 2)

#All flowering species and Diptera
Diptera <- read.csv(file = file.path("./Data", "Schoener Flower Diptera.csv"))
M3.4 <- brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = Diptera,
             family = 'zero_inflated_beta',
             iter = 1000,
             chains = 2,
             cores = 2)

#All flowering species and Hymenoptera
Hymenoptera <- read.csv(file = file.path("./Data", "Schoener Flower Hymenoptera.csv"))
M3.5 <- brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = Hymenoptera,
             family = 'beta',
             iter = 1000,
             chains = 2,
             cores = 2)

#All flowering species and Anthophila
Anthophila <- read.csv(file = file.path("./Data", "Schoener Flower Anthophila.csv"))
M3.6 <- brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = Anthophila,
             family = 'zero_inflated_beta',
             iter = 1000,
             chains = 2,
             cores = 2)

#All flowering species and Muscidae
Muscidae <- read.csv(file = file.path("./Data", "Schoener Flower Muscidae.csv"))
M3.7 <- brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = Muscidae,
             family = 'zero_inflated_beta',
             iter = 1000,
             chains = 2,
             cores = 2)

#All flowering species and Syrphidae
Syrphidae <- read.csv(file = file.path("./Data", "Schoener Flower Syrphidae.csv"))
M3.8 <- brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = Syrphidae,
             family = 'zero_inflated_beta',
             iter = 1000,
             chains = 2,
             cores = 2)



#All flowering species with their respective phenologies 
FlowersPhenology <- read.csv(file = file.path("./Data", "Schoener Flower Phenology.csv"))
M3.9 <- brm(Schoener ~ SpeciesAveragePhenology * Elevation + (1|Meadow/Site),
             data = FlowersPhenology,
             family = 'zero_one_inflated_beta',
             inits = 0,
             iter = 1000,
             chains = 2,
             cores = 2)

#Alternative with Standardised values the stabilise the model fitting
M3.10 <-brm(Schoener ~ SSpeciesAveragePhenology * SElevation + (1|Meadow/Site),
             data = FlowersPhenology,
             family = 'zero_one_inflated_beta',
             iter = 1000,
             chains = 2,
             cores = 2)



# Extracting coefficiences from brms####
summary(M6.7)
fixef(M3.1,
      summary = TRUE,
      robust = FALSE,
      probs = c(0.025, 0.975),
      pars = NULL
)


# Put your models into a named list
models <- list(
  M3.1 = M3.1,
  M3.2 = M3.2,
  M3.3 = M3.3,
  M3.4 = M3.4,
  M3.5 = M3.5,
  M3.6 = M3.6,
  M3.7 = M3.7,
  M3.8 = M3.8,
  #M3.9 = M3.9,
  M3.10 = M3.10
)

# Function to extract fixed effects and calculate ORs
extract_model_summary <- function(model, model_name) {
  fe <- fixef(model, summary = TRUE, probs = c(0.025, 0.975))
  
  fe_df <- as.data.frame(fe) %>%
    rownames_to_column("Variable") %>%
    mutate(
      Model = model_name,
      OR = exp(Estimate),
      OR_CI_low = exp(Q2.5),
      OR_CI_high = exp(Q97.5)
    )
  
  return(fe_df)
}

# Apply across all models
model_summaries <- purrr::imap_dfr(models, extract_model_summary)

# Clean up column names
model_summaries <- model_summaries %>%
  rename(
    Estimate = Estimate,
    SE = Est.Error,
    CI_low = Q2.5,
    CI_high = Q97.5
  ) %>%
  select(Model, Variable, Estimate, SE, CI_low, CI_high, OR, OR_CI_low, OR_CI_high)

# View combined summary
print(model_summaries)

#Plots####

#All in one####
make_fits <- function(model, data, label) {
  new_data <- data %>%
    distinct(Meadow, Site) %>%
    crossing(Elevation = seq(min(data$Elevation),
                             max(data$Elevation),
                             length.out = 100))
  
  fits <- fitted(model,
                 newdata = new_data,
                 re_formula = NA,
                 probs = c(0.025, 0.975)) %>%
    as.data.frame() %>%
    bind_cols(new_data) %>%
    mutate(Model = label)
  
  return(fits)
}

#Fits
fits_insects     <- make_fits(M3.2, Insects,     "All insects")
fits_hymenoptera <- make_fits(M3.5, Hymenoptera, "Hymenoptera")
fits_coleoptera  <- make_fits(M3.3, Coleoptera,  "Coleoptera")
fits_diptera     <- make_fits(M3.4, Diptera,     "Diptera")

# Combine all fits
fits_all <- bind_rows(fits_insects, fits_coleoptera, fits_diptera, fits_hymenoptera)

# Plot
p <- ggplot(fits_all, aes(x = Elevation, y = Estimate, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  labs(title = "Schoener overlap vs Elevation",
       y = "Schoener overlap",
       x = "Elevation") +
  theme_minimal() +
  theme(legend.position = "top")

p


#2+2 plots####
make_fits <- function(model, data, label) {
  new_data <- data %>%
    distinct(Meadow, Site) %>%
    crossing(Elevation = seq(min(data$Elevation),
                             max(data$Elevation),
                             length.out = 100))
  
  fits <- fitted(model,
                 newdata = new_data,
                 re_formula = NA,
                 probs = c(0.025, 0.975)) %>%
    as.data.frame() %>%
    bind_cols(new_data) %>%
    mutate(Model = label)
  
  return(fits)
}

#Fits
fits_insects     <- make_fits(M3.2, Insects,     "All insects")
fits_coleoptera  <- make_fits(M3.3, Coleoptera,  "Coleoptera")
fits_diptera     <- make_fits(M3.4, Diptera,     "Diptera")
fits_hymenoptera <- make_fits(M3.5, Hymenoptera, "Hymenoptera")

# Combine all fits and add grouping variable for facets
fits_all <- bind_rows(fits_insects, fits_coleoptera, fits_diptera, fits_hymenoptera) %>%
  mutate(Panel = case_when(
    Model %in% c("All insects", "Hymenoptera") ~ "Insects + Hymenoptera",
    Model %in% c("Coleoptera", "Diptera") ~ "Coleoptera + Diptera"
  ))

# Plot
p <- ggplot(fits_all, aes(x = Elevation, y = Estimate, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Panel, nrow = 1) +   # panels side by side
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Schoener overlap vs Elevation",
       y = "Schoener overlap",
       x = "Elevation") +
  theme_minimal() +
  theme(legend.position = "top")

p


# Plot 1: All insects + Hymenoptera
fits_group1 <- bind_rows(fits_insects, fits_hymenoptera)

p1 <- ggplot(fits_group1, aes(x = Elevation, y = Estimate, 
                              color = Model, fill = Model, linetype = Model)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_linetype_manual(values = c("All insects" = "solid",
                                   "Hymenoptera" = "dashed")) +
  scale_color_manual(values = c("All insects" = "blue",
                                "Hymenoptera" = "azure4")) +
  scale_fill_manual(values = c("All insects" = "blue",
                               "Hymenoptera" = "azure4")) +
  labs(
    y = "Schoener overlap",
    x = "Elevation") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    axis.line = element_line(color = "black")
  )


# Plot 2: Coleoptera + Diptera
fits_group2 <- bind_rows(fits_coleoptera, fits_diptera)

p2 <- ggplot(fits_group2, aes(x = Elevation, y = Estimate, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Coleoptera" = "forestgreen",
                                "Diptera"    = "gold2")) +
  scale_fill_manual(values = c("Coleoptera" = "forestgreen",
                               "Diptera"    = "gold2")) +
  labs(
    y = "Schoener overlap",
    x = "Elevation") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    axis.line = element_line(color = "black")
  )

# Print separately
p1
p2



#Muscidae, Syrphidae, Anthophila####
fits_anthophila <- make_fits(M3.6, Anthophila, "Anthophila")
fits_muscidae   <- make_fits(M3.7, Muscidae,   "Muscidae")
fits_syrphidae  <- make_fits(M3.8, Syrphidae,  "Syrphidae")

# Combine all fits
fits_all <- bind_rows(fits_anthophila, fits_muscidae, fits_syrphidae)


# Plot

p <- ggplot(fits_all, aes(x = Elevation, y = Estimate, color = Model, fill = Model, linetype = Model)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Anthophila" = "firebrick",
                                "Muscidae"   = "black",
                                "Syrphidae"  = "royalblue3")) +
  scale_fill_manual(values = c("Anthophila" = "firebrick",
                               "Muscidae"   = "black",
                               "Syrphidae"  = "royalblue3")) +
  scale_linetype_manual(values = c("Anthophila" = "dashed",
                                   "Muscidae"   = "dashed",
                                   "Syrphidae"  = "dashed")) +
  labs(#title = "Schoener overlap vs Elevation",
       y = "Schoener overlap",
       x = "Elevation") +
  theme_minimal(base_size = 20) + 
  theme(
    legend.position = "top",
    text = element_text(color = "black"),        # all text black
    axis.text = element_text(color = "black"),   # axis numbers black
    axis.title = element_text(color = "black"),  # axis labels black
    plot.title = element_text(color = "black"),  # title black
    axis.line = element_line(color = "black")   # add black axis lines
    #,panel.grid = element_blank()                 # optional: remove gray grid
  )

p
