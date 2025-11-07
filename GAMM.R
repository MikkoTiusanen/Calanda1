# GAMMs with interactions####
library(mgcv)
library(ggplot2)
library(geometry)
library(dplyr)
library(purrr)
library(broom)
library(tibble)

data <- read.csv(file = file.path("./Data", "FloweringInsectDensity.csv"))
data[9:176] <- lapply(data[9:176], as.numeric)

#Run GAMMs, first with elevation smoother, then without
#Total flowers
M1.1 <- gamm(
  Totalflowers ~ s(DOY) + s(MinOfelevation) + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)

summary(M1.1$gam)

M1.1 <- gamm(
  Totalflowers ~ s(DOY) + MinOfelevation + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)

summary(M1.1$gam)


#Flowering species
M1.2 <- gamm(
  Species ~ s(DOY) + s(MinOfelevation) + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)

summary(M1.2$gam)

M1.2 <- gamm(
  Species ~ s(DOY) + MinOfelevation + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)
summary(M1.2$gam)


#Insect density
M1.3 <- gamm(
  InsectDensity ~ s(DOY) + s(MinOfelevation) + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)
summary(M1.3$gam)

M1.3 <- gamm(
  InsectDensity ~ s(DOY) + MinOfelevation + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)
summary(M1.3$gam)

#Insects/ flower
M1.4 <- gamm(
  InsectsPerFlower ~ s(DOY) + s(MinOfelevation) + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)
summary(M1.4$gam)

M1.4 <- gamm(
  InsectsPerFlower ~ s(DOY) + MinOfelevation + ti(DOY, MinOfelevation),
  random = list(Meadow = ~1, Site = ~1 | Meadow),
  data = data
)
summary(M1.4$gam)


#Extracting model coefficients####
models <- list(
  TotalFlowers     = M1.1,
  Species          = M1.2,
  InsectDensity    = M1.3,
  InsectsPerFlower = M1.4
)

# Get the gam object
get_gam <- function(model) {
  if (inherits(model, "gam")) return(model)
  if (is.list(model) && !is.null(model$gam)) return(model$gam)
  stop("Model does not appear to be a 'gam' or 'gamm' object (no $gam found).")
}

# Find first matching column index
first_match_idx <- function(names_vec, patterns) {
  for (pat in patterns) {
    idx <- grep(pat, names_vec, ignore.case = TRUE)
    if (length(idx) > 0) return(idx[1])
  }
  integer(0)
}

# Extract smoothers
extract_smooth <- function(model, model_name) {
  gam_obj <- tryCatch(get_gam(model), error = function(e) return(NULL))
  if (is.null(gam_obj)) return(NULL)
  ssum <- summary(gam_obj)
  sm <- ssum$s.table
  if (is.null(sm)) return(NULL)
  sm_df <- as.data.frame(sm, stringsAsFactors = FALSE)
  sm_df$term <- rownames(sm_df)
  sm_df$model <- model_name
  rownames(sm_df) <- NULL
  
  # normalize names
  nm <- make.names(colnames(sm_df))
  
  # few different names ways to find the wanted parameters 
  edf_idx <- first_match_idx(nm, c("^edf$"))
  refdf_idx <- first_match_idx(nm, c("ref\\.?df", "ref_df", "ref\\.df"))
  F_idx <- first_match_idx(nm, c("^F$", "^F\\.value$", "F_value", "F\\.stat", "F\\."))
  p_idx <- first_match_idx(nm, c("p\\.?value", "^p$", "^Pr", "p\\.value", "Pr\\."))
  
  # choose a numeric column (other than edf/ref) as F if none found
  numeric_cols <- which(sapply(sm_df, is.numeric))
  if (length(F_idx) == 0) {
    exclude <- c(edf_idx, refdf_idx)
    cand <- setdiff(numeric_cols, exclude)
    if (length(cand) > 0) F_idx <- cand[1]
  }
  
  # assemble a tibble with consistent column names (use NA when missing)
  res <- tibble(
    model  = sm_df$model,
    term   = sm_df$term,
    edf    = if (length(edf_idx)) as.numeric(sm_df[[edf_idx]]) else NA_real_,
    ref_df = if (length(refdf_idx)) as.numeric(sm_df[[refdf_idx]]) else NA_real_,
    F      = if (length(F_idx)) as.numeric(sm_df[[F_idx]]) else NA_real_,
    p_value = if (length(p_idx)) as.numeric(sm_df[[p_idx]]) else NA_real_
  )
  return(res)
}

# Extract model parameters
extract_parametric <- function(model, model_name) {
  gam_obj <- tryCatch(get_gam(model), error = function(e) return(NULL))
  if (is.null(gam_obj)) return(NULL)
  ssum <- summary(gam_obj)
  pm <- ssum$p.table
  if (is.null(pm)) return(NULL)
  pm_df <- as.data.frame(pm, stringsAsFactors = FALSE)
  pm_df$term <- rownames(pm_df)
  pm_df$model <- model_name
  rownames(pm_df) <- NULL
  
  # normalize names
  nm <- make.names(colnames(pm_df))
  
  # detect columns
  est_idx  <- first_match_idx(nm, c("^Estimate$", "estimate"))
  se_idx   <- first_match_idx(nm, c("Std..Error", "Std.Error", "Std_Error", "Std\\.Error", "^se$"))
  t_idx    <- first_match_idx(nm, c("^t.value$", "t.value", "t_value", "^t$"))
  p_idx    <- first_match_idx(nm, c("Pr\\.\\.", "Pr", "p\\.?value", "^p$"))
  
  # back up if some not found
  numeric_cols <- which(sapply(pm_df, is.numeric))
  if (length(est_idx) == 0 && length(numeric_cols) > 0) est_idx <- numeric_cols[1]
  if (length(se_idx) == 0 && length(numeric_cols) > 1) se_idx <- numeric_cols[min(2, length(numeric_cols))]
  if (length(t_idx) == 0) {
    # try any column with 't' in name
    t_idx <- first_match_idx(nm, c("t\\.value", "t_value", "\\bt\\b"))
  }
  if (length(p_idx) == 0) {
    p_idx <- first_match_idx(nm, c("Pr", "p.value", "p_value"))
  }
  
  # build tibble
  res <- tibble(
    model    = pm_df$model,
    term     = pm_df$term,
    estimate = if (length(est_idx)) suppressWarnings(as.numeric(pm_df[[est_idx]])) else NA_real_,
    std.error= if (length(se_idx)) suppressWarnings(as.numeric(pm_df[[se_idx]])) else NA_real_,
    t.value  = if (length(t_idx)) suppressWarnings(as.numeric(pm_df[[t_idx]])) else NA_real_,
    p_value  = if (length(p_idx)) suppressWarnings(as.numeric(pm_df[[p_idx]])) else NA_real_
  )
  return(res)
}

# Apply to all models
safe_map2_dfr <- function(.x, .y, fn) {
  purrr::map2_dfr(.x, .y, function(m, n) {
    tryCatch(fn(m, n), error = function(e) {
      warning(sprintf("Failed to extract from model '%s': %s", n, e$message))
      NULL
    })
  })
}

smooth_terms     <- safe_map2_dfr(models, names(models), extract_smooth)
parametric_terms <- safe_map2_dfr(models, names(models), extract_parametric)

# Combine to a table
all_terms <- bind_rows(
  if (nrow(smooth_terms) > 0)     mutate(smooth_terms,     type = "smooth")     else tibble(),
  if (nrow(parametric_terms) > 0) mutate(parametric_terms, type = "parametric") else tibble()
)

# Show the results
message("Smooth terms (edf / Ref.df / F / p):")
print(smooth_terms)
message("Parametric terms (Estimate / Std. Error / t / p):")
print(parametric_terms)
message("All terms combined:")
print(all_terms)


#Plots####
# Flower abundance
# Extract the GAM component from gamm
gam_model <- M1.1$gam

# Build a prediction grid
DOY_seq <- seq(min(data$DOY, na.rm = TRUE),
               max(data$DOY, na.rm = TRUE), length.out = 1000)
Elev_seq <- seq(min(data$MinOfelevation, na.rm = TRUE),
                max(data$MinOfelevation, na.rm = TRUE), length.out = 1000)

newdat <- expand.grid(DOY = DOY_seq,
                      MinOfelevation = Elev_seq)

# Predict (type="response" gives fitted values)
newdat$fit <- predict(gam_model, newdata = newdat, type = "response")


# Make a matrix of observed predictor space
obs <- cbind(data$DOY, data$MinOfelevation)
grid <- cbind(newdat$DOY, newdat$MinOfelevation)

# Keep only points inside convex hull of observed data
inside <- inhulln(convhulln(obs), grid)
newdat$fit[!inside] <- NA

ggplot(newdat, aes(x = DOY, y = MinOfelevation, fill = fit)) +
  geom_tile(na.rm = TRUE) +
  scale_fill_gradientn(colors = c("blue4",
                                  "mediumorchid",
                                  "red",
                                  "yellow"), 
                       na.value = "white") +
  labs(fill = "Flower abundance",
       x = "Day of Year",
       y = "Elevation") +
  annotate("text", x = min(newdat$DOY), y = max(newdat$MinOfelevation), 
           label = "A", hjust = -0.5, vjust = 1.5, size = 14, fontface = "bold") +
  theme_minimal() +
  theme(
    # Axis settings
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    
    # Remove grid
    panel.grid = element_blank(),
    
    # Legend settings
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14)
  )

# Species number
# Extract the GAM component from gamm
gam_model <- M1.2$gam

# Build a prediction grid
DOY_seq <- seq(min(data$DOY, na.rm = TRUE),
               max(data$DOY, na.rm = TRUE), length.out = 1000)
Elev_seq <- seq(min(data$MinOfelevation, na.rm = TRUE),
                max(data$MinOfelevation, na.rm = TRUE), length.out = 1000)

newdat <- expand.grid(DOY = DOY_seq,
                      MinOfelevation = Elev_seq)

# Predict (type="response" gives fitted values)
newdat$fit <- predict(gam_model, newdata = newdat, type = "response")


# Make a matrix of observed predictor space
obs <- cbind(data$DOY, data$MinOfelevation)
grid <- cbind(newdat$DOY, newdat$MinOfelevation)

# Keep only points inside convex hull of observed data
inside <- inhulln(convhulln(obs), grid)
newdat$fit[!inside] <- NA

ggplot(newdat, aes(x = DOY, y = MinOfelevation, fill = fit)) +
  geom_tile(na.rm = TRUE) +
  scale_fill_gradientn(colors = c("blue4",
                                  "mediumorchid",
                                  "red",
                                  "yellow"), 
                       na.value = "white") +
  labs(fill = "Species number",
       x = "Day of Year",
       y = "Elevation") +
  annotate("text", x = min(newdat$DOY), y = max(newdat$MinOfelevation), 
           label = "B", hjust = -0.5, vjust = 1.5, size = 14, fontface = "bold") +
  theme_minimal() +
  theme(
    # Axis settings
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    
    # Remove grid
    panel.grid = element_blank(),
    
    # Legend settings
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14)
  )

# Insect density
# Extract the GAM component from gamm
gam_model <- M1.3$gam

# Build a prediction grid
DOY_seq <- seq(min(data$DOY, na.rm = TRUE),
               max(data$DOY, na.rm = TRUE), length.out = 1000)
Elev_seq <- seq(min(data$MinOfelevation, na.rm = TRUE),
                max(data$MinOfelevation, na.rm = TRUE), length.out = 1000)

newdat <- expand.grid(DOY = DOY_seq,
                      MinOfelevation = Elev_seq)

# Predict (type="response" gives fitted values)
newdat$fit <- predict(gam_model, newdata = newdat, type = "response")


# Make a matrix of observed predictor space
obs <- cbind(data$DOY, data$MinOfelevation)
grid <- cbind(newdat$DOY, newdat$MinOfelevation)

# Keep only points inside convex hull of observed data
inside <- inhulln(convhulln(obs), grid)
newdat$fit[!inside] <- NA

ggplot(newdat, aes(x = DOY, y = MinOfelevation, fill = fit)) +
  geom_tile(na.rm = TRUE) +
  scale_fill_gradientn(colors = c("blue4",
                                  "mediumorchid",
                                  "red",
                                  "yellow"), 
                       na.value = "white") +
  labs(fill = "Insect density",
       x = "Day of Year",
       y = "Elevation") +
  annotate("text", x = min(newdat$DOY), y = max(newdat$MinOfelevation), 
           label = "C", hjust = -0.5, vjust = 1.5, size = 14, fontface = "bold") +
  theme_minimal() +
  theme(
    # Axis settings
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    
    # Remove grid
    panel.grid = element_blank(),
    
    # Legend settings
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14)
  )

# Extract the GAM component from gamm
gam_model <- M1.4$gam

# Build a prediction grid
DOY_seq <- seq(min(data$DOY, na.rm = TRUE),
               max(data$DOY, na.rm = TRUE), length.out = 1000)
Elev_seq <- seq(min(data$MinOfelevation, na.rm = TRUE),
                max(data$MinOfelevation, na.rm = TRUE), length.out = 1000)

newdat <- expand.grid(DOY = DOY_seq,
                      MinOfelevation = Elev_seq)

# Predict (type="response" gives fitted values)
newdat$fit <- predict(gam_model, newdata = newdat, type = "response")


# Make a matrix of observed predictor space
obs <- cbind(data$DOY, data$MinOfelevation)
grid <- cbind(newdat$DOY, newdat$MinOfelevation)

# Keep only points inside convex hull of observed data
inside <- inhulln(convhulln(obs), grid)
newdat$fit[!inside] <- NA

ggplot(newdat, aes(x = DOY, y = MinOfelevation, fill = fit)) +
  geom_tile(na.rm = TRUE) +
  scale_fill_gradientn(
    colors = c("blue4", "mediumorchid", "red", "yellow"), 
    na.value = "white",
    breaks = c(0.0003, 0.0002, 0.0001, 0),
    labels = c("0.0003", "0.0002", "0.0001", "0")
  ) +
  labs(fill = "Insect/flowers",
       x = "Day of Year",
       y = "Elevation") +
  annotate("text", x = min(newdat$DOY), y = max(newdat$MinOfelevation), 
           label = "D", hjust = -0.5, vjust = 1.5, size = 14, fontface = "bold") +
  theme_minimal() +
  theme(
    # Axis settings
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    
    # Remove grid
    panel.grid = element_blank(),
    
    # Legend settings
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14)
  )
