###############################################################
##   One-Stage DLNM by Region (Categorical Exposure)          ##
###############################################################

rm(list = ls())  # clear environment to avoid conflicts

# ------------------------------
# 0. Libraries & Setup
# ------------------------------
library(dlnm)       # distributed lag non-linear models
library(splines)    # spline functions
library(gnm)        # generalized nonlinear models (handles eliminate=)
library(ggplot2)    # plotting
library(dplyr)      # data wrangling
library(qs)         # fast save/load
library(readr)      # csv I/O
library(MASS)       # for mvrnorm
library(geobr)      # Brazil municipality shapefiles
library(sf)         # geometry dropping

# Create output folder if missing
dir.create("modeling_outputs_categorical_region", showWarnings = FALSE, recursive = TRUE)

# ------------------------------
# 1. Load & Stack Data (≥100 events set)
# ------------------------------
data_list <- readRDS("data_list_100_ischemic.rds")

stacked <- bind_rows(lapply(seq_along(data_list), function(i) {
  df <- data_list[[i]]
  df$id <- names(data_list)[i]
  df
}))

# ------------------------------
# 2. Add Region Metadata
# ------------------------------
muni_meta <- geobr::read_municipality(code_muni = "all", year = 2020) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(abbrev_state, name_muni, name_region) %>%
  dplyr::mutate(
    name_region = dplyr::recode(name_region,
                                "Norte"        = "North",
                                "Nordeste"     = "Northeast",
                                "Sudeste"      = "Southeast",
                                "Sul"          = "South",
                                "Centro Oeste" = "Central-West"
    )
  )

stacked <- stacked %>%
  left_join(muni_meta, by = c("abbrev_state","name_muni")) %>%
  mutate(
    id      = factor(id),
    stratum = factor(stratum),
    block   = interaction(id, stratum, drop = TRUE),
    region  = factor(name_region)
  )

# ------------------------------
# 2. Define categorical exposures (by municipality: 1%, 10%, 90%, 99%)
# ------------------------------
stacked <- stacked %>%
  group_by(id) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    # Municipality-specific thresholds
    p01 = quantile(app_temp, 0.01, na.rm = TRUE),
    p10 = quantile(app_temp, 0.10, na.rm = TRUE),
    p90 = quantile(app_temp, 0.90, na.rm = TRUE),
    p99 = quantile(app_temp, 0.99, na.rm = TRUE),
    
    # Heat indicators
    EHE_90 = as.integer(app_temp >= p90),
    EHE_99 = as.integer(app_temp >= p99),
    
    # Cold indicators
    ECE_10 = as.integer(app_temp <= p10),
    ECE_01 = as.integer(app_temp <= p01),
    
    # Multi-day persistence
    EHE_90_2 = as.integer(EHE_90 == 1 & lead(EHE_90, 1, default = 0) == 1),
    EHE_90_3 = as.integer(EHE_90 == 1 & lead(EHE_90, 1, default = 0) == 1 &
                            lead(EHE_90, 2, default = 0) == 1),
    EHE_99_2 = as.integer(EHE_99 == 1 & lead(EHE_99, 1, default = 0) == 1),
    EHE_99_3 = as.integer(EHE_99 == 1 & lead(EHE_99, 1, default = 0) == 1 &
                            lead(EHE_99, 2, default = 0) == 1),
    
    ECE_10_2 = as.integer(ECE_10 == 1 & lead(ECE_10, 1, default = 0) == 1),
    ECE_10_3 = as.integer(ECE_10 == 1 & lead(ECE_10, 1, default = 0) == 1 &
                            lead(ECE_10, 2, default = 0) == 1),
    ECE_01_2 = as.integer(ECE_01 == 1 & lead(ECE_01, 1, default = 0) == 1),
    ECE_01_3 = as.integer(ECE_01 == 1 & lead(ECE_01, 1, default = 0) == 1 &
                            lead(ECE_01, 2, default = 0) == 1)
  ) %>%
  ungroup()

#Inspect thresholds by municipality
stacked %>%
  group_by(id) %>%
  summarise(
    p01 = unique(p01),
    p10 = unique(p10),
    p90 = unique(p90),
    p99 = unique(p99)
  ) %>%
  print(n = 20)

# Define which indicators to run
indicators <- c("EHE_90","EHE_99","EHE_90_2","EHE_90_3",
                "EHE_99_2","EHE_99_3",
                "ECE_10","ECE_01","ECE_10_2","ECE_10_3",
                "ECE_01_2","ECE_01_3")


indicator_region_summary <- stacked %>%
  group_by(region) %>%
  summarise(
    # Total flagged days
    across(all_of(indicators), ~ sum(. == 1, na.rm = TRUE), .names = "days_{.col}"),
    # Total deaths on flagged days
    across(all_of(indicators), ~ sum(events[. == 1], na.rm = TRUE), .names = "deaths_{.col}"),
    # Mean temperature on flagged days
    across(all_of(indicators), ~ mean(app_temp[. == 1], na.rm = TRUE), .names = "meanT_{.col}"),
    # Median temperature on flagged days
    across(all_of(indicators), ~ median(app_temp[. == 1], na.rm = TRUE), .names = "medianT_{.col}"),
    # Average days per municipality
    across(all_of(indicators), ~ sum(. == 1, na.rm = TRUE) / n_distinct(id[. == 1]), .names = "avg_days_per_muni_{.col}"),
    .groups = "drop"
  )

print(indicator_region_summary)

# Save to CSV
readr::write_csv(indicator_region_summary, "modeling_outputs_categorical_region/indicator_region_summary.csv")

# ------------------------------
# 4. Region-level model fitting
# ------------------------------

lag_used <- 7
or_table <- tibble()

# Loop over regions
for (reg in unique(stacked$region)) {
  cat("=== Region:", reg, "===\n")
  
  stacked_region <- stacked %>% filter(region == reg)
  
  # Loop over indicators
  for (var in indicators) {
    cat("Fitting model for:", var, "in region:", reg, "\n")
    
    x <- stacked_region[[var]]
    
    # Skip if indicator never occurs or too rare
    if (is.null(x) || sum(x, na.rm = TRUE) < 3) {
      cat("  -> Skipping", var, "in", reg, "(insufficient exposure days)\n")
      next
    }
    
    # Build crossbasis for binary exposure (lagged effect up to lag_used)
    cb <- crossbasis(x, lag = lag_used,
                     argvar = list(fun = "lin"),   # linear for binary
                     arglag = list(fun = "ns", knots = logknots(lag_used, 2)))
    
    # Only keep strata with events
    ind <- tapply(stacked_region$events, stacked_region$block, sum)[stacked_region$block]
    
    # Fit conditional Poisson model
    fit <- tryCatch(
      gnm(events ~ cb + pm25 + so2 + co + no2 + go3,
          family    = quasipoisson,
          eliminate = block,
          subset    = ind > 0,
          data      = stacked_region),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      cat("  -> Model failed for", var, "in", reg, "\n")
      next
    }
    
    # Skip if cb terms were dropped
    if (!any(grepl("^cb", names(coef(fit))))) {
      cat("  -> Model dropped cb terms for", var, "in", reg, "\n")
      next
    }
    
    # Predict cumulative RR for exposed vs not exposed
    cp <- tryCatch(
      crosspred(cb, fit, cen = 0, at = 1, cumul = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(cp)) {
      cat("  -> crosspred failed for", var, "in", reg, "\n")
      next
    }
    
    rr  <- cp$cumRRfit ["1", "lag7"]
    lcl <- cp$cumRRlow ["1", "lag7"]
    ucl <- cp$cumRRhigh["1", "lag7"]
    
    or_table <- bind_rows(or_table, tibble(
      region   = reg,
      variable = var,
      OR  = rr,
      LCL = lcl,
      UCL = ucl
    ))
    
    # Cleanup to save memory
    rm(cb, fit, cp); gc()
  }
}

print(or_table)
readr::write_csv(or_table, "modeling_outputs_categorical_region/or_table_by_region.csv")

# ------------------------------
# 5. Regional Attributable Burden (AF/AN) via Simulation
# ------------------------------
library(Matrix)

nsim <- 100
set.seed(123)

# helper function to compute AF/AN
af_an_from_beta <- function(beta_l, x, y, trunc_excess = TRUE, event_mask = NULL) {
  n <- length(x)
  L <- length(beta_l) - 1
  eta <- numeric(n)
  
  # convolve lag kernel with exposure history
  for (l in 0:L) {
    idx <- (l + 1):n
    if (length(idx) > 0) {
      eta[idx] <- eta[idx] + beta_l[l + 1] * x[seq_len(n - l)]
    }
  }
  
  RR <- exp(eta)
  RR <- pmax(RR, 1e-12)  # guard against underflow
  
  # attributable numbers per day
  if (trunc_excess) {
    AN <- y * pmax(0, 1 - 1 / RR)
  } else {
    AN <- y * (1 - 1 / RR)
  }
  
  AN_tot <- sum(AN, na.rm = TRUE)
  AF_tot <- if (sum(y, na.rm = TRUE) > 0) AN_tot / sum(y, na.rm = TRUE) else NA_real_
  
  # --- NEW: per-event-day metric ---
  if (is.null(event_mask)) {
    event_mask <- (x == 1)  # assume binary event coding
  }
  idx_ev <- which(event_mask)
  AN_per_event_day <- if (length(idx_ev) > 0) mean(AN[idx_ev], na.rm = TRUE) else NA_real_
  
  list(AF = AF_tot, AN = AN_tot, AN_per_event_day = AN_per_event_day)
}

afan_rows <- list()

for (reg in unique(stacked$region)) {
  cat("=== Region:", reg, "AF/AN ===\n")
  stacked_region <- stacked %>% filter(region == reg)
  
  for (var in indicators) {
    cat("Running AF/AN for:", var, "in region:", reg, "\n")
    
    x <- stacked_region[[var]]
    y <- stacked_region$events
    
    # skip if no exposure days
    if (sum(x == 1, na.rm = TRUE) == 0 || sum(y, na.rm = TRUE) == 0) {
      afan_rows[[length(afan_rows) + 1]] <- data.frame(
        region = reg, variable = var, days = sum(x == 1, na.rm = TRUE), 
        events = sum(y[x == 1], na.rm = TRUE),
        AF_med = NA_real_, AF_low = NA_real_, AF_high = NA_real_,
        AN_med = NA_real_, AN_low = NA_real_, AN_high = NA_real_
      )
      next
    }
    
    # rebuild model for this region + indicator
    cb <- crossbasis(x, lag = lag_used,
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "ns", knots = logknots(lag_used, 2)))
    ind <- tapply(stacked_region$events, stacked_region$block, sum)[stacked_region$block]
    
    fit <- tryCatch(
      gnm(events ~ cb + pm25 + so2 + co + no2 + go3,
          family = quasipoisson,
          eliminate = block,
          subset = ind > 0,
          data = stacked_region),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    coef_all <- coef(fit)
    vcov_all <- vcov(fit)
    ix <- grep("^cb", names(coef_all))
    if (length(ix) == 0) next
    
    b <- coef_all[ix]
    V <- vcov_all[ix, ix, drop = FALSE]
    
    # lag basis
    L <- max(attr(cb, "lag"))   # gives 7
    arglag <- attr(cb, "arglag")
    B <- do.call(dlnm::onebasis, c(list(x = 0:L), arglag))
    
    # ensure covariance PD
    V_pd <- tryCatch({ chol(V); V }, error = function(e) as.matrix(nearPD(V, corr = FALSE)$mat))
    
    # simulate
    b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V_pd)
    beta_draws <- t(B %*% t(b_draws))
    
    AFv <- numeric(nsim); ANv <- numeric(nsim); ANpev <- numeric(nsim)
    for (i in seq_len(nsim)) {
      out <- af_an_from_beta(beta_draws[i, ], x, y, trunc_excess = TRUE)
      AFv[i] <- out$AF
      ANv[i] <- out$AN
      ANpev[i] <- out$AN_per_event_day
    }
    
    afan_rows[[length(afan_rows) + 1]] <- data.frame(
      region = reg,
      variable = var,
      days = sum(x == 1, na.rm = TRUE),
      events = sum(y[x == 1], na.rm = TRUE),
      AF_med = median(AFv, na.rm = TRUE),
      AF_low = quantile(AFv, 0.025, na.rm = TRUE),
      AF_high= quantile(AFv, 0.975, na.rm = TRUE),
      AN_med = median(ANv, na.rm = TRUE),
      AN_low = quantile(ANv, 0.025, na.rm = TRUE),
      AN_high= quantile(ANv, 0.975, na.rm = TRUE),
      AN_per_event_med  = median(ANpev, na.rm = TRUE),
      AN_per_event_low  = quantile(ANpev, 0.025, na.rm = TRUE),
      AN_per_event_high = quantile(ANpev, 0.975, na.rm = TRUE)
    )
    
    rm(cb, fit); gc()
  }
}

afan_table_region <- bind_rows(afan_rows)
print(afan_table_region)

readr::write_csv(afan_table_region, "modeling_outputs_categorical_region/afan_table_by_region.csv")
