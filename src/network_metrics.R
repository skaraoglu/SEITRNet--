# ============================================================================
# network_metrics.R — Network Topology Metrics for Cross-Chapter Analysis
# ============================================================================
#
# Computes network-level summary statistics on contact networks generated
# by igraph. These are the SAME metrics used in Chapter 4's brain network
# analysis (H4: clustering C, path length L, small-world sigma, modularity Q),
# enabling direct cross-chapter comparison.
#
# The C++ simulation kernel uses an internal flat adjacency matrix and does
# not expose network structure to R. This file generates fresh igraph
# networks (with the same parameters and generators) and computes metrics
# on those networks. Since metrics are measured at t=0 (before demographic
# events), this is statistically equivalent to measuring inside the kernel.
#
# Chapter 1's Algorithm 1 already records (P_k, C, l_avg, |C_max|) per
# timestep during simulation. Here we replicate and extend those metrics
# in a standalone function for the post-hoc Experiment 2 analysis.
# ============================================================================

library(igraph)

# --------------------------------------------------------------------------
# Generate an igraph network (matching the C++ kernel's generators)
# --------------------------------------------------------------------------
generate_igraph_network <- function(network_type, n, n_par1, n_par2 = 10L) {
  if (network_type == "ER") {
    igraph::sample_gnp(n = n, p = n_par1, directed = FALSE)
  } else if (network_type == "BA") {
    m <- max(1L, round(n_par1 * n))
    igraph::sample_pa(n = n, power = 1, m = m, directed = FALSE)
  } else if (network_type == "WS") {
    igraph::sample_smallworld(dim = 1, size = n, nei = n_par2, p = n_par1)
  } else {
    stop(paste("Unknown network type:", network_type))
  }
}

# --------------------------------------------------------------------------
# Compute topology metrics on a single igraph object
#
# Returns a named list matching Chapter 4's H4 metrics plus degree stats.
# --------------------------------------------------------------------------
compute_topology_metrics <- function(g) {
  n_nodes <- igraph::vcount(g)
  n_edges <- igraph::ecount(g)
  deg     <- igraph::degree(g)

  # Degree distribution moments
  deg_mean     <- mean(deg)
  deg_var      <- var(deg)
  deg_skew     <- if (sd(deg) > 0) mean(((deg - mean(deg)) / sd(deg))^3) else 0
  deg_max      <- max(deg)

  # Clustering coefficient (global transitivity) — same as Chapter 4
  C <- igraph::transitivity(g, type = "global")
  if (is.nan(C)) C <- 0

  # Average shortest path length (largest connected component) — same as Chapter 4
  comps    <- igraph::components(g)
  lcc_ids  <- which(comps$membership == which.max(comps$csize))
  lcc_frac <- length(lcc_ids) / n_nodes

  if (length(lcc_ids) > 1) {
    g_lcc <- igraph::induced_subgraph(g, lcc_ids)
    L     <- igraph::mean_distance(g_lcc, directed = FALSE)
  } else {
    L <- Inf
  }

  # Small-world index sigma = (C/C_rand) / (L/L_rand) — same formula as Chapter 4
  # C_rand and L_rand for an ER random graph with same n and mean degree
  p_equiv <- deg_mean / (n_nodes - 1)
  C_rand  <- p_equiv  # Expected clustering for ER G(n,p)
  L_rand  <- if (p_equiv > 0 && p_equiv < 1) {
    log(n_nodes) / log(deg_mean)  # Approximation for connected ER
  } else {
    1
  }
  sigma <- if (C_rand > 0 && L_rand > 0 && L > 0 && is.finite(L)) {
    (C / max(C_rand, 1e-10)) / (L / max(L_rand, 1e-10))
  } else {
    NA_real_
  }

  # Modularity (Louvain community detection) — same as Chapter 4
  if (n_edges > 0) {
    comm <- igraph::cluster_louvain(igraph::as.undirected(g))
    Q    <- igraph::modularity(comm)
  } else {
    Q <- 0
  }

  # Degree assortativity
  r_assort <- tryCatch(
    igraph::assortativity_degree(g, directed = FALSE),
    error = function(e) NA_real_
  )

  list(
    n_nodes   = n_nodes,
    n_edges   = n_edges,
    deg_mean  = deg_mean,
    deg_var   = deg_var,
    deg_skew  = deg_skew,
    deg_max   = deg_max,
    C         = C,
    L         = L,
    sigma     = sigma,
    Q         = Q,
    r_assort  = r_assort,
    lcc_frac  = lcc_frac
  )
}

# --------------------------------------------------------------------------
# Compute averaged metrics across multiple network realisations
#
# Generates n_reps networks and averages the metrics. This accounts for
# the randomness in network generation (especially for ER and WS).
# --------------------------------------------------------------------------
compute_avg_topology_metrics <- function(network_type, n, n_par1, n_par2 = 10L,
                                          n_reps = 20L) {
  all_metrics <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    g <- generate_igraph_network(network_type, n, n_par1, n_par2)
    all_metrics[[r]] <- compute_topology_metrics(g)
  }

  # Average across replicates
  metric_names <- names(all_metrics[[1]])
  avg <- list()
  for (nm in metric_names) {
    vals <- sapply(all_metrics, `[[`, nm)
    avg[[nm]]           <- mean(vals, na.rm = TRUE)
    avg[[paste0(nm, "_sd")]] <- sd(vals, na.rm = TRUE)
  }
  avg$network_type <- network_type
  avg$n            <- n
  avg$n_par1       <- n_par1
  avg$n_par2       <- n_par2
  avg$n_reps       <- n_reps

  avg
}

# --------------------------------------------------------------------------
# Compute control profile shape descriptors
#
# Quantifies how the optimised u1(t) differs from the ODE reference.
# These are the descriptors from Experiment 4 in the roadmap.
# --------------------------------------------------------------------------
compute_control_shape <- function(u1_profile, ode_u1_interp = NULL, dt = 1) {
  T_max <- (length(u1_profile) - 1) * dt
  u     <- u1_profile

  # Total effort (L1 norm, integrated)
  total_effort <- sum(u) * dt

  # Peak intensity
  peak <- max(u)

  # Effective duration: fraction of time where u > 10% of peak
  if (peak > 0) {
    eff_duration <- mean(u > 0.1 * peak)
  } else {
    eff_duration <- 0
  }

  # Front-loading index: effort in first half / total effort
  half_idx <- ceiling(length(u) / 2)
  front_loading <- if (total_effort > 0) {
    sum(u[1:half_idx]) / sum(u)
  } else {
    NA_real_
  }

  # Temporal variability: CV of u across time
  cv <- if (mean(u) > 0) sd(u) / mean(u) else NA_real_

  # ODE correlation (if ODE reference provided)
  ode_cor <- if (!is.null(ode_u1_interp) && length(ode_u1_interp) == length(u)) {
    cor(u, ode_u1_interp, use = "complete.obs")
  } else {
    NA_real_
  }

  list(
    total_effort  = total_effort,
    peak          = peak,
    eff_duration  = eff_duration,
    front_loading = front_loading,
    cv            = cv,
    ode_cor       = ode_cor
  )
}
