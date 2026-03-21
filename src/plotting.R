# ============================================================================
# plotting.R — Visualization Functions for SEITR Experiments
# ============================================================================
#
# Contains all ggplot-based visualization functions used in the analysis.
#
# WHAT'S NEW vs. the original notebook:
#   - Extracted from inline notebook cells into reusable, documented functions.
#   - plot_experiment_diagnostics() is NEW: combines optimization trace +
#     control profile + compartment dynamics into a single call, replacing
#     the manually assembled grid.arrange blocks.
#   - Deprecated ggplot2 size= parameter replaced with linewidth= to suppress
#     warnings in ggplot2 >= 3.4.
#   - All functions accept the structured output from run_full_experiment(),
#     eliminating the need to manually extract df_avg, df_min, df_max, etc.
# ============================================================================

library(ggplot2)
library(gridExtra)

# --------------------------------------------------------------------------
# Plot a single SEITR compartment with min/max shading and no-control dashed
#
# UNCHANGED logic from the original plot_compartment_band(), but:
#   - Uses linewidth= instead of deprecated size= for line geoms
#   - Accepts data directly rather than requiring global variable names
# --------------------------------------------------------------------------
plot_compartment_band <- function(df_avg, df_min, df_max, df_avg_no_control,
                                  comp, color, title) {
  df <- data.frame(
    time           = df_avg$time,
    avg            = df_avg[[comp]],
    min_val        = df_min[[comp]],
    max_val        = df_max[[comp]],
    avg_no_control = df_avg_no_control[[comp]]
  )
  ggplot(df, aes(x = time)) +
    geom_ribbon(aes(ymin = min_val, ymax = max_val), fill = color, alpha = 0.2) +
    geom_line(aes(y = avg), color = color, linewidth = 1) +
    geom_line(aes(y = avg_no_control), color = color, linetype = "dashed", linewidth = 1) +
    labs(title = paste0(title, " (Network Opt)"), x = "Time", y = "Count") +
    theme_minimal()
}

# --------------------------------------------------------------------------
# Plot compartment with ODE overlay (controlled + uncontrolled ODE curves)
#
# UNCHANGED logic from original plot_compartment_band_overlay(), but:
#   - Uses linewidth= instead of deprecated size=
#   - Accepts ODE data.frame with standardized column names
# --------------------------------------------------------------------------
plot_compartment_band_overlay <- function(df_avg, df_min, df_max,
                                          df_avg_no_control, comp, color, title,
                                          ode_results, after, before) {
  df <- data.frame(
    time           = df_avg$time,
    avg            = df_avg[[comp]],
    min_val        = df_min[[comp]],
    max_val        = df_max[[comp]],
    avg_no_control = df_avg_no_control[[comp]]
  )
  ggplot(df, aes(x = time)) +
    geom_ribbon(aes(ymin = min_val, ymax = max_val), fill = color, alpha = 0.2) +
    geom_line(aes(y = avg), color = color, linewidth = 1.2) +
    geom_line(aes(y = avg_no_control), color = color, linetype = "dashed", linewidth = 1.2) +
    geom_line(data = ode_results, aes(x = time, y = !!as.name(after)),
              color = "black", linewidth = 0.75, linetype = "dashed") +
    geom_line(data = ode_results, aes(x = time, y = !!as.name(before)),
              color = "gray20", linewidth = 0.75, linetype = "dotted") +
    labs(title = title, x = "Time", y = "Count") +
    theme_minimal()
}

# --------------------------------------------------------------------------
# Plot ODE optimal control results (convergence + control + compartments)
#
# UNCHANGED visualization logic, but packaged as a single function call.
# --------------------------------------------------------------------------
plot_ode_results <- function(ode_sol) {
  # Convergence trace
  df_obj <- data.frame(Iteration = ode_sol$iterations,
                       Objective = ode_sol$objective_history)
  p_obj <- ggplot(df_obj, aes(x = Iteration, y = Objective)) +
    geom_line(color = "black", linewidth = 1) +
    labs(title = "Objective Functional", x = "Iterations", y = "J(u)") +
    theme_minimal()

  # Optimal control profile
  df_u1 <- data.frame(Time = ode_sol$time, u1 = ode_sol$u1)
  p_u1 <- ggplot(df_u1, aes(x = Time, y = u1)) +
    geom_line(color = "red", linewidth = 1) +
    labs(title = "Optimal Control u1", x = "Time (Days)", y = "u1(t)") +
    theme_minimal()

  gridExtra::grid.arrange(p_obj, p_u1, ncol = 2)

  # Compartment trajectories (before/after optimization)
  df_states <- data.frame(
    Time = ode_sol$time,
    S_before = ode_sol$S_uncontrolled, S_after = ode_sol$S,
    E_before = ode_sol$E_uncontrolled, E_after = ode_sol$E,
    I_before = ode_sol$I_uncontrolled, I_after = ode_sol$I,
    T_before = ode_sol$T_uncontrolled, T_after = ode_sol$T_state,
    R_before = ode_sol$R_uncontrolled, R_after = ode_sol$R
  )

  plot_comp <- function(df, before, after, title, color) {
    ggplot(df, aes(x = Time)) +
      geom_line(aes(y = !!as.name(before)), color = "black", linewidth = 1) +
      geom_line(aes(y = !!as.name(after)), color = color, linewidth = 1) +
      labs(title = title, x = "Time (Days)", y = title) +
      theme_minimal()
  }

  p_S <- plot_comp(df_states, "S_before", "S_after", "Susceptibles", "blue")
  p_E <- plot_comp(df_states, "E_before", "E_after", "Exposed", "orange")
  p_I <- plot_comp(df_states, "I_before", "I_after", "Infected", "red")
  p_T <- plot_comp(df_states, "T_before", "T_after", "Treatment", "purple")
  p_R <- plot_comp(df_states, "R_before", "R_after", "Recovered", "green")

  gridExtra::grid.arrange(p_S, p_E, p_I, p_T, p_R, ncol = 2)
}

# --------------------------------------------------------------------------
# NEW: Plot full experiment diagnostics from run_full_experiment() output
#
# Produces two panels:
#   Panel 1: Objective functional (network mean + min/max band) + ODE overlay,
#            and the optimal control profile (network vs ODE).
#   Panel 2: Five compartment plots with min/max bands, no-control baseline,
#            and ODE trajectory overlays.
# --------------------------------------------------------------------------
plot_experiment_diagnostics <- function(exp_result, ode_sol, w1) {

  df_avg   <- exp_result$results_opt$avg
  df_min   <- exp_result$results_opt$min
  df_max   <- exp_result$results_opt$max
  df_noctl <- exp_result$results_noctl$avg  # May be NULL if run_no_control=FALSE
  u1_opt   <- exp_result$u1_profile

  # Ensure u1_profile length matches data
  if (length(u1_opt) < nrow(df_avg)) {
    u1_opt <- c(u1_opt, rep(tail(u1_opt, 1), nrow(df_avg) - length(u1_opt)))
  } else if (length(u1_opt) > nrow(df_avg)) {
    u1_opt <- u1_opt[seq_len(nrow(df_avg))]
  }

  # --- Panel 1: Objective + Control ---

  # ODE instantaneous objective interpolated to network grid
  n_ode <- length(ode_sol$time)
  if (n_ode %% 2 == 0) n_ode <- n_ode - 1
  ode_times <- ode_sol$time[1:n_ode]
  J_ode_inst <- ode_sol$E[1:n_ode] + ode_sol$I[1:n_ode] + w1 * ode_sol$u1[1:n_ode]^2
  J_ode_interp <- approx(x = ode_times, y = J_ode_inst,
                          xout = df_avg$time, rule = 2)$y

  df_obj_mm <- data.frame(
    time = df_avg$time,
    min  = df_min$E + df_min$I + df_min$control_cost,
    max  = df_max$E + df_max$I + df_max$control_cost,
    avg  = df_avg$E + df_avg$I + df_avg$control_cost,
    ode  = J_ode_interp
  )

  p_obj <- ggplot(df_obj_mm, aes(x = time)) +
    geom_ribbon(aes(ymin = min, ymax = max), fill = "gray70", alpha = 0.3) +
    geom_line(aes(y = avg), color = "black", linewidth = 1) +
    geom_line(aes(y = ode), color = "blue", linewidth = 0.7,
              linetype = "dashed", alpha = 0.8) +
    labs(title = "Obj Func: Net (black), ODE (blue)",
         x = "Time", y = "Instantaneous Objective") +
    theme_minimal()

  # Control profile
  df_u1_net <- data.frame(time = df_avg$time, u1 = u1_opt)
  df_u1_ode <- data.frame(time = ode_sol$time, u1 = ode_sol$u1)

  p_u1 <- ggplot(df_u1_net, aes(x = time, y = u1)) +
    geom_line(color = "red", linewidth = 1) +
    geom_line(data = df_u1_ode, aes(x = time, y = u1),
              color = "blue", linewidth = 0.5, linetype = "dashed") +
    labs(title = "Control: Net (red), ODE (blue)", x = "Time", y = "u1(t)") +
    theme_minimal()

  gridExtra::grid.arrange(p_obj, p_u1, ncol = 2)

  # --- Panel 2: Compartment dynamics with ODE overlay ---
  ode_df <- data.frame(
    time     = ode_sol$time,
    S_after  = ode_sol$S,          E_after  = ode_sol$E,
    I_after  = ode_sol$I,          T_after  = ode_sol$T_state,
    R_after  = ode_sol$R,
    S_before = ode_sol$S_uncontrolled, E_before = ode_sol$E_uncontrolled,
    I_before = ode_sol$I_uncontrolled, T_before = ode_sol$T_uncontrolled,
    R_before = ode_sol$R_uncontrolled
  )

  # Use no-control baseline if available; otherwise use avg as placeholder
  # (no-control dashed line will overlap the controlled line in that case)
  noctl_df <- if (!is.null(df_noctl)) df_noctl else df_avg

  p_S <- plot_compartment_band_overlay(df_avg, df_min, df_max, noctl_df,
           "S", "blue", "Susceptible", ode_df, "S_after", "S_before")
  p_E <- plot_compartment_band_overlay(df_avg, df_min, df_max, noctl_df,
           "E", "orange", "Exposed", ode_df, "E_after", "E_before")
  p_I <- plot_compartment_band_overlay(df_avg, df_min, df_max, noctl_df,
           "I", "red", "Infected", ode_df, "I_after", "I_before")
  p_T <- plot_compartment_band_overlay(df_avg, df_min, df_max, noctl_df,
           "T", "purple", "Treatment", ode_df, "T_after", "T_before")
  p_R <- plot_compartment_band_overlay(df_avg, df_min, df_max, noctl_df,
           "R", "green", "Recovered", ode_df, "R_after", "R_before")

  gridExtra::grid.arrange(p_S, p_E, p_I, p_T, p_R, ncol = 2)

  # Print summary
  if (!is.null(exp_result$optim_J)) {
    cat("Optimizer-reported J:", round(exp_result$optim_J, 4), "\n")
  }
  cat("Forward check J:    ", round(exp_result$J_opt$J_total, 4), "\n")
  if (!is.null(exp_result$J_noctl)) {
    cat("No-control J:       ", round(exp_result$J_noctl$J_total, 4), "\n")
  }
}
