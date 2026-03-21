# ============================================================================
# utils.R — Shared Utility Functions
# ============================================================================
#
# Contains numerical integration, control-profile expansion, and data helpers
# used by both the ODE solver and the network optimization pipeline.
#
# WHAT'S NEW vs. the original notebook:
#   - Functions previously defined inline in notebook cells are now centralized.
#   - get_ode_initial_guess() is NEW: downsamples the ODE optimal control to K
#     piecewise-constant segments for optimizer warm-start (Recommendation #8).
# ============================================================================

# --------------------------------------------------------------------------
# Simpson's rule integration of the SEITR objective functional
# UNCHANGED from the original.
# --------------------------------------------------------------------------
calculate_objective_functional <- function(results_df) {
  n <- nrow(results_df)
  if (n %% 2 == 0) n <- n - 1
  times        <- results_df$time[1:n]
  E            <- results_df$E[1:n]
  I_vals       <- results_df$I[1:n]
  control_cost <- results_df$control_cost[1:n]
  weights      <- rep(2, n)
  weights[1]   <- 1
  weights[n]   <- 1
  weights[seq(2, n - 1, by = 2)] <- 4
  h       <- mean(diff(times))
  J_E     <- (h / 3) * sum(weights * E)
  J_I     <- (h / 3) * sum(weights * I_vals)
  J_W     <- (h / 3) * sum(weights * control_cost)
  J_total <- J_E + J_I + J_W
  return(list(J_E = J_E, J_I = J_I, J_W = J_W, J_total = J_total))
}

# --------------------------------------------------------------------------
# Expand K piecewise-constant segment parameters to a full time-grid profile
# UNCHANGED from the original.
# --------------------------------------------------------------------------
expand_u1 <- function(params_seg, total_length, interval_length, zeta) {
  prof <- rep(params_seg, each = interval_length)[1:total_length]
  pmin(zeta, pmax(0, prof))
}

# --------------------------------------------------------------------------
# NEW: Extract an informed initial guess from the ODE optimal control.
#
# FIX: The original used (k-1)*interval_length+1 as idx_start, which
# overruns total_steps when K does not divide total_steps evenly (e.g.,
# K=20, total_steps=101 → interval_length=6, segments 18-20 start at
# indices 103,109,115 which are past the 101-element vector → NaN).
#
# The fix uses the SAME segment assignment as expand_u1: build a vector
# of segment IDs via rep(1:K, each=interval_length)[1:total_steps], then
# average the ODE values within each segment. This guarantees every
# segment maps to at least one time point and the segmentation is
# consistent with how expand_u1 will later reconstruct the profile.
# --------------------------------------------------------------------------
get_ode_initial_guess <- function(ode_u1, ode_time, t_max, K, dt = 1) {
  net_times       <- seq(0, t_max, by = dt)
  u1_interp       <- approx(x = ode_time, y = ode_u1, xout = net_times, rule = 2)$y
  total_steps     <- length(net_times)
  interval_length <- ceiling(total_steps / K)

  # Segment assignment vector — identical logic to expand_u1
  seg_ids <- rep(1:K, each = interval_length)[1:total_steps]

  seg_values <- numeric(K)
  for (k in 1:K) {
    idx <- which(seg_ids == k)
    if (length(idx) > 0) {
      seg_values[k] <- mean(u1_interp[idx])
    } else {
      # Fallback: segment has no assigned points (should not happen now)
      seg_values[k] <- seg_values[max(1L, k - 1L)]
    }
  }
  return(seg_values)
}
