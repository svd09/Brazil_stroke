###############################################################
##   One-Stage DLNM (Categorical Exposure, 4 cutpoints)      ##
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

# Create output folder if missing
dir.create("modeling_outputs_categorical", showWarnings = FALSE, recursive = TRUE)

# ------------------------------
# 1. Load & Stack Data (≥100 events set)
# ------------------------------
data_list <- readRDS("data_list_100_all.rds")

stacked <- bind_rows(lapply(seq_along(data_list), function(i) {
  df <- data_list[[i]]
  df$id <- names(data_list)[i]
  df
}))

stacked <- stacked %>%
  mutate(
    id      = factor(id),
    stratum = factor(stratum),
    block   = interaction(id, stratum, drop = TRUE)
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


# Build summary table (overall)
indicator_summary <- lapply(indicators, function(var) {
  # Subset the rows where the indicator is 1
  subset_data <- stacked[stacked[[var]] == 1, ]
  
  data.frame(
    variable      = var,
    n_days        = nrow(subset_data),                           # total days
    deaths        = sum(subset_data$events, na.rm = TRUE),      # total deaths
    mean_temp     = mean(subset_data$app_temp, na.rm = TRUE),   # mean temp
    median_temp   = median(subset_data$app_temp, na.rm = TRUE), # median temp
    avg_days_per_muni = nrow(subset_data) / length(unique(subset_data$id)) # avg days per municipality
  )
})

indicator_summary <- do.call(rbind, indicator_summary)

# Print or save
print(indicator_summary)
readr::write_csv(indicator_summary, "modeling_outputs_categorical/indicator_summary.csv")

# ------------------------------
# 3. Fit One-Stage Models (per categorical indicator)
# ------------------------------

lag_used <- 7

# Create an empty results table
or_table <- tibble()

# will hold cb + model per indicator for AF/AN sims
results_models <- vector("list", length(indicators))
names(results_models) <- indicators


# Loop over each indicator
for (var in indicators) {
  cat("Fitting model for:", var, "\n")
  
  x <- stacked[[var]]
  
  # Build crossbasis for binary exposure (lagged effect up to lag_used)
  cb <- crossbasis(x, lag = lag_used,
                   argvar = list(fun = "lin"),   # linear because binary
                   arglag = list(fun = "ns", knots = logknots(lag_used, 2)))
  
  # Only keep strata with events
  ind <- tapply(stacked$events, stacked$block, sum)[stacked$block]
  
  # Fit conditional Poisson model
  fit <- tryCatch(
    gnm(events ~ cb + pm25 + so2 + co + no2 + go3,
        family    = quasipoisson,
        eliminate = block,
        subset    = ind > 0,
        data      = stacked),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    # Predict cumulative RR for exposed vs not exposed
    cp <- crosspred(cb, fit, cen = 0, at = 1, cumul = TRUE)
    
    # Extract summary results
    #rr   <- cp$allRRfit[which(cp$predvar == 1)]
    #lcl  <- cp$allRRlow[which(cp$predvar == 1)]
    #ucl  <- cp$allRRhigh[which(cp$predvar == 1)]
    
    rr  <- cp$cumRRfit ["1", "lag7"]
    lcl <- cp$cumRRlow ["1", "lag7"]
    ucl <- cp$cumRRhigh["1", "lag7"]
    
    # Append to results table
    or_table <- bind_rows(or_table, tibble(
      variable = var,
      OR  = rr,
      LCL = lcl,
      UCL = ucl
    ))
    
    # save objects for AF/AN simulations
    results_models[[var]] <- list(cb = cb, model = fit)
    
  }
  
  # Cleanup to save memory
  rm(cb, fit, cp); gc()
}

print(or_table)
write_csv(or_table, "modeling_outputs_categorical/or_table.csv")

# ------------------------------
# 5. Attributable burden (AF/AN) via simulation
# ------------------------------
library(MASS)     # mvrnorm
library(Matrix)   # nearPD for numerical safety

nsim <- 100
set.seed(123)

# compute AF & AN from a lag-kernel (beta_l), exposure x, and counts y
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

for (ind in indicators) {
  res <- results_models[[ind]]
  if (is.null(res)) next  # model failed / was skipped
  
  cb <- res$cb
  m  <- res$model
  x  <- stacked[[ind]]
  y  <- stacked$events
  
  # if no exposed days or no events, skip cleanly
  if (sum(x == 1, na.rm = TRUE) == 0 || sum(y, na.rm = TRUE) == 0) {
    afan_rows[[length(afan_rows) + 1]] <- data.frame(
      variable = ind, days = sum(x == 1, na.rm = TRUE), events = sum(y[x == 1], na.rm = TRUE),
      AF_med = NA_real_, AF_low = NA_real_, AF_high = NA_real_,
      AN_med = NA_real_, AN_low = NA_real_, AN_high = NA_real_
    )
    next
  }
  
  # pull cb coefficients & vcov
  coef_all <- coef(m)
  vcov_all <- vcov(m)
  ix <- grep("^cb", names(coef_all))
  if (length(ix) == 0) {
    # nothing to simulate if cb terms not found
    afan_rows[[length(afan_rows) + 1]] <- data.frame(
      variable = ind, days = sum(x == 1, na.rm = TRUE), events = sum(y[x == 1], na.rm = TRUE),
      AF_med = NA_real_, AF_low = NA_real_, AF_high = NA_real_,
      AN_med = NA_real_, AN_low = NA_real_, AN_high = NA_real_
    )
    next
  }
  b <- coef_all[ix]
  V <- vcov_all[ix, ix, drop = FALSE]
  
  # rebuild lag basis to map cb coefs -> lag-specific effects
  L      <- max(attr(cb, "lag"))   # gives 7
  arglag <- attr(cb, "arglag")
  B <- do.call(dlnm::onebasis, c(list(x = 0:L), arglag))
  
  # ensure covariance is positive definite
  V_pd <- tryCatch({ chol(V); V }, error = function(e) as.matrix(nearPD(V, corr = FALSE)$mat))
  
  # simulate coefficient draws and map to lag kernels
  b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V_pd)
  beta_draws <- t(B %*% t(b_draws))  # each row is lag-specific log-RR (length L+1)
  
  AFv <- numeric(nsim)
  ANv <- numeric(nsim)
  ANpev <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    out <- af_an_from_beta(beta_draws[i, ], x, y, trunc_excess = TRUE)
    AFv[i] <- out$AF
    ANv[i] <- out$AN
    ANpev[i] <- out$AN_per_event_day
  }
  
  afan_rows[[length(afan_rows) + 1]] <- data.frame(
    variable = ind,
    days   = sum(x == 1, na.rm = TRUE),
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
}

afan_table <- bind_rows(afan_rows)
print(afan_table)
readr::write_csv(afan_table, "modeling_outputs_categorical/afan_table.csv")

# ------------------------------
# 6. Save results
# ------------------------------
#qs::qsave(results_models, "modeling_outputs_categorical/categorical_models.qs")
