# ============================================================================
# experiment.R — Flexible Experiment Runner
# ============================================================================
#
# Provides run_experiment() with a modular design: the user can run
# optimization, forward simulation, or both. This replaces the monolithic
# run_full_experiment() that forced all stages to run together.
#
# USAGE PATTERNS:
#
#   1. Optimization + forward check (full pipeline):
#      run_experiment(..., K = 5, n_starts = 10)
#
#   2. Forward simulation only (with a known control profile):
#      run_experiment(..., u1_profile = my_u1)
#
#   3. No-control baseline only:
#      run_experiment(..., u1_profile = rep(0, total_steps))
#
#   4. Skip no-control baseline:
#      run_experiment(..., K = 5, run_no_control = FALSE)
# ============================================================================

library(parallel)
library(optimParallel)

# --------------------------------------------------------------------------
# run_experiment: Flexible experiment runner
#
#   If u1_profile is NULL AND K is provided, runs optimization first.
#   If u1_profile is provided, skips optimization entirely.
#   If run_no_control = TRUE, also runs a u1=0 baseline for comparison.
#
#   Returns a list with all results and metadata.
# --------------------------------------------------------------------------
run_experiment <- function(
  # Network specification
  network_type, n, n_par1, n_par2,
  # Epidemiological parameters
  Lambda, beta1, beta2, beta3, alpha1, alpha2, delta_I, delta_T, mu,
  # Control parameters
  w1, zeta,
  # Simulation setup
  t_max, dt = 1,
  init_S, init_E, init_I, init_T, init_R,
  # --- Mode control ---
  u1_profile    = NULL,      # If provided: skip optimization, use this profile
  K             = NULL,      # If provided (and u1_profile is NULL): run optimization
  n_starts      = 10,        # Multi-start count (optimization only)
  ode_init_guess = NULL,     # ODE warm-start segments (optimization only)
  run_no_control = TRUE,     # Also run u1=0 baseline?
  # Simulation replicates
  num_exp_forward = 20L,     # Replicates for forward simulation
  # Parallel setup
  n_cores = max(1, parallel::detectCores() - 2),
  cpp_path = "src/seitr_kernel.cpp"
) {

  total_steps <- length(seq(0, t_max, by = dt))

  # Determine mode
  do_optimize <- is.null(u1_profile) && !is.null(K)

  if (is.null(u1_profile) && is.null(K)) {
    stop("Either u1_profile (for forward sim) or K (for optimization) must be provided.")
  }

  cat("=== Experiment:", network_type, "n=", n, "par1=", n_par1)
  if (do_optimize) cat(" K=", K)
  cat(" ===\n")

  # ===========================================================================
  # OPTIMIZATION (only if u1_profile not provided)
  # ===========================================================================
  optim_output <- NULL

  if (do_optimize) {
    interval_length <- ceiling(total_steps / K)

    cat("Setting up parallel cluster with", n_cores, "cores...\n")
    cl <- makeCluster(n_cores)
    setDefaultCluster(cl)

    # Compile C++ kernel on each worker
    cpp_path_abs <- normalizePath(cpp_path)
    clusterExport(cl, "cpp_path_abs", envir = environment())
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_path_abs)
    })

    # Reproducible parallel RNG
    RNGkind("L'Ecuyer-CMRG")
    clusterSetRNGStream(cl, 12345)

    # Export parameters as individual globals for worker access
    .exp_net_type   <- network_type
    .exp_n          <- n
    .exp_n_par1     <- n_par1
    .exp_n_par2     <- n_par2
    .exp_Lambda     <- Lambda
    .exp_beta1      <- beta1
    .exp_beta2      <- beta2
    .exp_beta3      <- beta3
    .exp_alpha1     <- alpha1
    .exp_alpha2     <- alpha2
    .exp_delta_I    <- delta_I
    .exp_delta_T    <- delta_T
    .exp_mu         <- mu
    .exp_w1         <- w1
    .exp_zeta       <- zeta
    .exp_t_max      <- t_max
    .exp_dt         <- dt
    .exp_init_S     <- init_S
    .exp_init_E     <- init_E
    .exp_init_I     <- init_I
    .exp_init_T     <- init_T
    .exp_init_R     <- init_R
    .exp_total_steps     <- total_steps
    .exp_interval_length <- interval_length

    clusterExport(cl, c(
      ".exp_net_type", ".exp_n", ".exp_n_par1", ".exp_n_par2",
      ".exp_Lambda", ".exp_beta1", ".exp_beta2", ".exp_beta3",
      ".exp_alpha1", ".exp_alpha2", ".exp_delta_I", ".exp_delta_T", ".exp_mu",
      ".exp_w1", ".exp_zeta", ".exp_t_max", ".exp_dt",
      ".exp_init_S", ".exp_init_E", ".exp_init_I", ".exp_init_T", ".exp_init_R",
      ".exp_total_steps", ".exp_interval_length"
    ), envir = environment())

    clusterExport(cl, c("expand_u1", "calculate_objective_functional"),
                  envir = globalenv())

    # Objective functions for early/late stages
    .obj_early <- function(params_seg) {
      set.seed(12345)
      u1p <- expand_u1(params_seg, .exp_total_steps, .exp_interval_length, .exp_zeta)
      res <- run_seitr_simulation_cpp(
        network_type_str = .exp_net_type, n = .exp_n,
        n_par1 = .exp_n_par1, n_par2 = .exp_n_par2,
        Lambda = .exp_Lambda, beta1 = .exp_beta1, beta2 = .exp_beta2, beta3 = .exp_beta3,
        alpha1 = .exp_alpha1, alpha2 = .exp_alpha2,
        delta_I = .exp_delta_I, delta_T = .exp_delta_T, mu = .exp_mu,
        init_S = .exp_init_S, init_E = .exp_init_E, init_I = .exp_init_I,
        init_T = .exp_init_T, init_R = .exp_init_R,
        t_max = .exp_t_max, u1_profile = u1p, num_exp = 5L, w1 = .exp_w1
      )
      calculate_objective_functional(res$avg)$J_total
    }

    .obj_late <- function(params_seg) {
      set.seed(12345)
      u1p <- expand_u1(params_seg, .exp_total_steps, .exp_interval_length, .exp_zeta)
      res <- run_seitr_simulation_cpp(
        network_type_str = .exp_net_type, n = .exp_n,
        n_par1 = .exp_n_par1, n_par2 = .exp_n_par2,
        Lambda = .exp_Lambda, beta1 = .exp_beta1, beta2 = .exp_beta2, beta3 = .exp_beta3,
        alpha1 = .exp_alpha1, alpha2 = .exp_alpha2,
        delta_I = .exp_delta_I, delta_T = .exp_delta_T, mu = .exp_mu,
        init_S = .exp_init_S, init_E = .exp_init_E, init_I = .exp_init_I,
        init_T = .exp_init_T, init_R = .exp_init_R,
        t_max = .exp_t_max, u1_profile = u1p, num_exp = 20L, w1 = .exp_w1
      )
      calculate_objective_functional(res$avg)$J_total
    }

    clusterExport(cl, c(".obj_early", ".obj_late"), envir = environment())

    # Multi-start optimization
    lower_bounds <- rep(0, K)
    upper_bounds <- rep(zeta, K)
    optim_results <- vector("list", n_starts)

    for (i in seq_len(n_starts)) {
      if (i == 1 && !is.null(ode_init_guess)) {
        init_seg <- ode_init_guess
        cat("  Start", i, "/", n_starts, ": ODE warm-start\n")
      } else {
        init_seg <- runif(K, 0, zeta)
        cat("  Start", i, "/", n_starts, ": random init\n")
      }

      # Stage 1: cheap exploration
      res1 <- tryCatch(
        optimParallel(
          par = init_seg, fn = .obj_early, method = "L-BFGS-B",
          lower = lower_bounds, upper = upper_bounds,
          control = list(maxit = 500, factr = 1e8),
          parallel = list(cl = cl, forward = FALSE, loginfo = FALSE)
        ),
        error = function(e) {
          cat("    Stage 1 error:", e$message, "\n")
          list(par = init_seg, value = .obj_early(init_seg))
        }
      )

      # Stage 2: refinement
      optim_results[[i]] <- tryCatch(
        optimParallel(
          par = res1$par, fn = .obj_late, method = "L-BFGS-B",
          lower = lower_bounds, upper = upper_bounds,
          control = list(maxit = 1000, factr = 1e7),
          parallel = list(cl = cl, forward = FALSE, loginfo = FALSE)
        ),
        error = function(e) {
          cat("    Stage 2 error:", e$message, "\n")
          res1
        }
      )
      cat("    -> J =", round(optim_results[[i]]$value, 4), "\n")
    }

    stopCluster(cl)

    # Select best
    best_idx <- which.min(sapply(optim_results, function(x) {
      if (is.null(x$value)) Inf else x$value
    }))
    best <- optim_results[[best_idx]]
    cat("Best start:", best_idx, "with J =", round(best$value, 4), "\n")

    u1_segments <- best$par
    u1_profile  <- expand_u1(u1_segments, total_steps, interval_length, zeta)

    optim_output <- list(
      best        = best,
      all_results = optim_results,
      best_idx    = best_idx,
      u1_segments = u1_segments,
      K           = K,
      interval_length = interval_length
    )
  }

  # ===========================================================================
  # FORWARD SIMULATION (always runs)
  # ===========================================================================
  cat("Running forward simulation (", num_exp_forward, "replicates)...\n")
  results_opt <- run_seitr_simulation_cpp(
    network_type_str = network_type, n = n, n_par1 = n_par1, n_par2 = n_par2,
    Lambda = Lambda, beta1 = beta1, beta2 = beta2, beta3 = beta3,
    alpha1 = alpha1, alpha2 = alpha2, delta_I = delta_I, delta_T = delta_T, mu = mu,
    init_S = init_S, init_E = init_E, init_I = init_I,
    init_T = init_T, init_R = init_R,
    t_max = t_max, u1_profile = u1_profile, num_exp = num_exp_forward, w1 = w1
  )
  J_opt <- calculate_objective_functional(results_opt$avg)

  # ===========================================================================
  # NO-CONTROL BASELINE (optional)
  # ===========================================================================
  results_noctl <- NULL
  J_noctl       <- NULL

  if (run_no_control) {
    cat("Running no-control baseline (", num_exp_forward, "replicates)...\n")
    u1_zero <- rep(0, total_steps)
    results_noctl <- run_seitr_simulation_cpp(
      network_type_str = network_type, n = n, n_par1 = n_par1, n_par2 = n_par2,
      Lambda = Lambda, beta1 = beta1, beta2 = beta2, beta3 = beta3,
      alpha1 = alpha1, alpha2 = alpha2, delta_I = delta_I, delta_T = delta_T, mu = mu,
      init_S = init_S, init_E = init_E, init_I = init_I,
      init_T = init_T, init_R = init_R,
      t_max = t_max, u1_profile = u1_zero, num_exp = num_exp_forward, w1 = w1
    )
    J_noctl <- calculate_objective_functional(results_noctl$avg)
  }

  # ===========================================================================
  # REPORT
  # ===========================================================================
  if (!is.null(optim_output)) {
    cat("  Optimizer J  =", round(optim_output$best$value, 4), "\n")
  }
  cat("  Forward J    =", round(J_opt$J_total, 4), "\n")
  if (!is.null(J_noctl)) {
    cat("  No-control J =", round(J_noctl$J_total, 4), "\n")
  }
  cat("=== Done ===\n\n")

  return(list(
    # Optimization output (NULL if u1_profile was supplied directly)
    optim_output  = optim_output,
    optim_J       = if (!is.null(optim_output)) optim_output$best$value else NULL,
    # Control profile used
    u1_profile    = u1_profile,
    # Forward simulation results
    results_opt   = results_opt,
    J_opt         = J_opt,
    # No-control baseline (NULL if run_no_control = FALSE)
    results_noctl = results_noctl,
    J_noctl       = J_noctl,
    # Experiment metadata
    params = list(
      network_type = network_type, n = n, n_par1 = n_par1, n_par2 = n_par2,
      K = K, total_steps = total_steps
    )
  ))
}
