# ============================================================================
# ode_control.R — Deterministic ODE Optimal Control Solver
# ============================================================================
#
# DIRECT TRANSLITERATION of the original notebook cells 12-14.
# Wrapped in a function for reuse. No algorithmic changes whatsoever.
#
# The forward-backward sweep method (Pontryagin's Maximum Principle) for
# the SEITR model with time-dependent treatment control u1(t).
#
# CRITICAL: Every formula, every RK4 midpoint convention, the control
# update rule, the convergence criterion, and the objective computation
# are copied verbatim from the original code. The only structural change
# is wrapping the inline code into a callable function.
# ============================================================================

solve_ode_optimal_control <- function(
  Lambda, beta1, beta2, beta3, alpha1, alpha2, delta_I, delta_T, mu,
  w1, zeta,
  h = 0.01, t_max = 100,
  S0, E0, I0, T0, R0,
  delta_conv = 0.001  # Convergence tolerance (called 'delta' in original)
) {

  # --- Derived quantities (matches original cell 12) ---
  k1 <- beta3 + mu + delta_I + alpha1
  k2 <- beta2 + mu
  k3 <- mu + delta_T + alpha2
  R0_val <- (beta1 * beta2) / (k1 * k2)
  cat("R0 =", R0_val, "\n")

  # --- Time grid ---
  t <- seq(0, t_max, by = h)
  L <- length(t)

  # --- State variable initialization ---
  S <- E <- I <- T_v <- R_v <- N <- numeric(L)
  S[1] <- S0; E[1] <- E0; I[1] <- I0; T_v[1] <- T0; R_v[1] <- R0
  N[1] <- S0 + E0 + I0 + T0 + R0

  # --- Control and adjoint initialization (matches original cell 12) ---
  u1 <- rep(0.03, L)  # Original: u1_n100 <- rep(0.03, L)
  lambda1 <- lambda2 <- lambda3 <- lambda4 <- lambda5 <- numeric(L)
  # Terminal conditions = 0 (already zero from numeric(L))

  itr <- 0
  test <- -1

  # --- State equation RHS (matches original cell 13 EXACTLY) ---
  f1 <- function(S, I, N)     Lambda - (beta1 * S * I) / N - mu * S
  f2 <- function(S, I, N, E)  (beta1 * S * I) / N - (beta2 + mu) * E
  f3 <- function(E, I, u1)    beta2 * E - (beta3 + mu + delta_I + u1) * I
  f4 <- function(u1, I, Tv)   u1 * I - (mu + delta_T + alpha2) * Tv
  f5 <- function(I, Tv, R)    beta3 * I + alpha2 * Tv - mu * R
  f6 <- function(N, I, Tv)    Lambda - N * mu - delta_I * I - delta_T * Tv

  # --- Adjoint equation RHS (matches original cell 14 EXACTLY) ---
  g1 <- function(lambda1, I, N, lambda2) {
    (beta1 * I * lambda1) / N + mu * lambda1 - (lambda2 * beta1 * I) / N
  }
  g2 <- function(lambda2, lambda3) {
    -1 + beta2 * lambda2 + mu * lambda2 - beta2 * lambda3
  }
  g3 <- function(lambda1, S, N, lambda2, lambda3, u1, lambda4, lambda5) {
    -1 + (beta1 * S * lambda1) / N - (beta1 * S * lambda2) / N +
      beta3 * lambda3 + mu * lambda3 + delta_I * lambda3 +
      u1 * lambda3 - lambda4 * u1 - lambda5 * beta3
  }
  g4 <- function(lambda4, lambda5) {
    mu * lambda4 + delta_T * lambda4 + alpha2 * lambda4 - alpha2 * lambda5
  }
  g5 <- function(lambda5) {
    lambda5 * mu
  }

  # --- Initial forward sweep (before iteration, for uncontrolled baseline) ---
  for (i in 1:(L - 1)) {
    m1 <- f1(S[i], I[i], N[i])
    n1 <- f2(S[i], I[i], N[i], E[i])
    o1 <- f3(E[i], I[i], u1[i])
    p1 <- f4(u1[i], I[i], T_v[i])
    q1 <- f5(I[i], T_v[i], R_v[i])
    r1 <- f6(N[i], I[i], T_v[i])

    m2 <- f1(S[i]+0.5*h*m1, I[i]+0.5*h*o1, N[i]+0.5*h*r1)
    n2 <- f2(S[i]+0.5*h*m1, I[i]+0.5*h*o1, N[i]+0.5*h*r1, E[i]+0.5*h*n1)
    o2 <- f3(E[i]+0.5*h*n1, I[i]+0.5*h*o1, u1[i]+0.5*h)
    p2 <- f4(u1[i]+0.5*h, I[i]+0.5*h*o1, T_v[i]+0.5*h*p1)
    q2 <- f5(I[i]+0.5*h*o1, T_v[i]+0.5*h*p1, R_v[i]+0.5*h*q1)
    r2 <- f6(N[i]+0.5*h*r1, I[i]+0.5*h*o1, T_v[i]+0.5*h*p1)

    m3 <- f1(S[i]+0.5*h*m2, I[i]+0.5*h*o2, N[i]+0.5*h*r2)
    n3 <- f2(S[i]+0.5*h*m2, I[i]+0.5*h*o2, N[i]+0.5*h*r2, E[i]+0.5*h*n2)
    o3 <- f3(E[i]+0.5*h*n2, I[i]+0.5*h*o2, u1[i]+0.5*h)
    p3 <- f4(u1[i]+0.5*h, I[i]+0.5*h*o2, T_v[i]+0.5*h*p2)
    q3 <- f5(I[i]+0.5*h*o2, T_v[i]+0.5*h*p2, R_v[i]+0.5*h*q2)
    r3 <- f6(N[i]+0.5*h*r2, I[i]+0.5*h*o2, T_v[i]+0.5*h*p2)

    m4 <- f1(S[i]+h*m3, I[i]+h*o3, N[i]+h*r3)
    n4 <- f2(S[i]+h*m3, I[i]+h*o3, N[i]+h*r3, E[i]+h*n3)
    o4 <- f3(E[i]+h*n3, I[i]+h*o3, u1[i]+h)
    p4 <- f4(u1[i]+h, I[i]+h*o3, T_v[i]+h*p3)
    q4 <- f5(I[i]+h*o3, T_v[i]+h*p3, R_v[i]+h*q3)
    r4 <- f6(N[i]+h*r3, I[i]+h*o3, T_v[i]+h*p3)

    S[i+1]    <- S[i]    + (h/6)*(m1 + 2*m2 + 2*m3 + m4)
    E[i+1]    <- E[i]    + (h/6)*(n1 + 2*n2 + 2*n3 + n4)
    I[i+1]    <- I[i]    + (h/6)*(o1 + 2*o2 + 2*o3 + o4)
    T_v[i+1]  <- T_v[i]  + (h/6)*(p1 + 2*p2 + 2*p3 + p4)
    R_v[i+1]  <- R_v[i]  + (h/6)*(q1 + 2*q2 + 2*q3 + q4)
    N[i+1]    <- N[i]    + (h/6)*(r1 + 2*r2 + 2*r3 + r4)
  }

  # --- Convergence history ---
  FV_history  <- numeric()
  itr_history <- numeric()

  # =========================================================================
  # FORWARD-BACKWARD SWEEP ITERATION
  # Copied VERBATIM from original cell 14. Every RK4 midpoint, every
  # convergence criterion term, every indexing convention is preserved.
  # =========================================================================
  while (test < 0) {
    itr <- itr + 1

    oldu1 <- u1
    oldS <- S; oldE <- E; oldI <- I; oldT <- T_v; oldR <- R_v
    oldlambda1 <- lambda1; oldlambda2 <- lambda2; oldlambda3 <- lambda3
    oldlambda4 <- lambda4; oldlambda5 <- lambda5

    # --- FORWARD: state equations (RK4) ---
    for (i in 1:(L - 1)) {
      m1 <- f1(S[i], I[i], N[i])
      n1 <- f2(S[i], I[i], N[i], E[i])
      o1 <- f3(E[i], I[i], u1[i])
      p1 <- f4(u1[i], I[i], T_v[i])
      q1 <- f5(I[i], T_v[i], R_v[i])
      r1 <- f6(N[i], I[i], T_v[i])

      m2 <- f1(S[i]+0.5*h*m1, I[i]+0.5*h*o1, N[i]+0.5*h*r1)
      n2 <- f2(S[i]+0.5*h*m1, I[i]+0.5*h*o1, N[i]+0.5*h*r1, E[i]+0.5*h*n1)
      o2 <- f3(E[i]+0.5*h*n1, I[i]+0.5*h*o1, u1[i]+0.5*h)
      p2 <- f4(u1[i]+0.5*h, I[i]+0.5*h*o1, T_v[i]+0.5*h*p1)
      q2 <- f5(I[i]+0.5*h*o1, T_v[i]+0.5*h*p1, R_v[i]+0.5*h*q1)
      r2 <- f6(N[i]+0.5*h*r1, I[i]+0.5*h*o1, T_v[i]+0.5*h*p1)

      m3 <- f1(S[i]+0.5*h*m2, I[i]+0.5*h*o2, N[i]+0.5*h*r2)
      n3 <- f2(S[i]+0.5*h*m2, I[i]+0.5*h*o2, N[i]+0.5*h*r2, E[i]+0.5*h*n2)
      o3 <- f3(E[i]+0.5*h*n2, I[i]+0.5*h*o2, u1[i]+0.5*h)
      p3 <- f4(u1[i]+0.5*h, I[i]+0.5*h*o2, T_v[i]+0.5*h*p2)
      q3 <- f5(I[i]+0.5*h*o2, T_v[i]+0.5*h*p2, R_v[i]+0.5*h*q2)
      r3 <- f6(N[i]+0.5*h*r2, I[i]+0.5*h*o2, T_v[i]+0.5*h*p2)

      m4 <- f1(S[i]+h*m3, I[i]+h*o3, N[i]+h*r3)
      n4 <- f2(S[i]+h*m3, I[i]+h*o3, N[i]+h*r3, E[i]+h*n3)
      o4 <- f3(E[i]+h*n3, I[i]+h*o3, u1[i]+h)
      p4 <- f4(u1[i]+h, I[i]+h*o3, T_v[i]+h*p3)
      q4 <- f5(I[i]+h*o3, T_v[i]+h*p3, R_v[i]+h*q3)
      r4 <- f6(N[i]+h*r3, I[i]+h*o3, T_v[i]+h*p3)

      S[i+1]    <- S[i]    + (h/6)*(m1 + 2*m2 + 2*m3 + m4)
      E[i+1]    <- E[i]    + (h/6)*(n1 + 2*n2 + 2*n3 + n4)
      I[i+1]    <- I[i]    + (h/6)*(o1 + 2*o2 + 2*o3 + o4)
      T_v[i+1]  <- T_v[i]  + (h/6)*(p1 + 2*p2 + 2*p3 + p4)
      R_v[i+1]  <- R_v[i]  + (h/6)*(q1 + 2*q2 + 2*q3 + q4)
      N[i+1]    <- N[i]    + (h/6)*(r1 + 2*r2 + 2*r3 + r4)
    }

    # Save uncontrolled baseline on first iteration
    if (itr == 1) {
      S1 <- S; E1 <- E; I1 <- I; T1 <- T_v; R1 <- R_v; N1 <- N
    }

    # --- BACKWARD: adjoint equations (RK4, backward in time) ---
    # EXACT COPY of original backward sweep. The midpoint convention
    # uses "state[j] - 0.5*h" (MATLAB FBSM style), NOT interpolated
    # midpoints. This is the convention from the original code.
    for (i in 1:(L - 1)) {
      j <- L + 1 - i

      k1 <- g1(lambda1[j], I[j], N[j], lambda2[j])
      l1 <- g2(lambda2[j], lambda3[j])
      c1 <- g3(lambda1[j], S[j], N[j], lambda2[j], lambda3[j], u1[j], lambda4[j], lambda5[j])
      v1 <- g4(lambda4[j], lambda5[j])
      d1 <- g5(lambda5[j])

      k2 <- g1(lambda1[j]-0.5*h*k1, I[j]-0.5*h, N[j]-0.5*h, lambda2[j]-0.5*h*l1)
      l2 <- g2(lambda2[j]-0.5*h*l1, lambda3[j]-0.5*h*c1)
      c2 <- g3(lambda1[j]-0.5*h*k1, S[j]-0.5*h, N[j]-0.5*h, lambda2[j]-0.5*h*l1, lambda3[j]-0.5*h*c1, u1[j]-0.5*h, lambda4[j]-0.5*h*v1, lambda5[j]-0.5*h*d1)
      v2 <- g4(lambda4[j]-0.5*h*v1, lambda5[j]-0.5*h*d1)
      d2 <- g5(lambda5[j]-0.5*h*d1)

      k3 <- g1(lambda1[j]-0.5*h*k2, I[j]-0.5*h, N[j]-0.5*h, lambda2[j]-0.5*h*l2)
      l3 <- g2(lambda2[j]-0.5*h*l2, lambda3[j]-0.5*h*c2)
      c3 <- g3(lambda1[j]-0.5*h*k2, S[j]-0.5*h, N[j]-0.5*h, lambda2[j]-0.5*h*l2, lambda3[j]-0.5*h*c2, u1[j]-0.5*h, lambda4[j]-0.5*h*v2, lambda5[j]-0.5*h*d2)
      v3 <- g4(lambda4[j]-0.5*h*v2, lambda5[j]-0.5*h*d2)
      d3 <- g5(lambda5[j]-0.5*h*d2)

      k4 <- g1(lambda1[j]-h*k3, I[j]-h, N[j]-h, lambda2[j]-h*l3)
      l4 <- g2(lambda2[j]-h*l3, lambda3[j]-h*c3)
      c4 <- g3(lambda1[j]-h*k3, S[j]-h, N[j]-h, lambda2[j]-h*l3, lambda3[j]-h*c3, u1[j]-h, lambda4[j]-h*v3, lambda5[j]-h*d3)
      v4 <- g4(lambda4[j]-h*v3, lambda5[j]-h*d3)
      d4 <- g5(lambda5[j]-h*d3)

      lambda1[j-1] <- lambda1[j] - (h/6)*(k1 + 2*k2 + 2*k3 + k4)
      lambda2[j-1] <- lambda2[j] - (h/6)*(l1 + 2*l2 + 2*l3 + l4)
      lambda3[j-1] <- lambda3[j] - (h/6)*(c1 + 2*c2 + 2*c3 + c4)
      lambda4[j-1] <- lambda4[j] - (h/6)*(v1 + 2*v2 + 2*v3 + v4)
      lambda5[j-1] <- lambda5[j] - (h/6)*(d1 + 2*d2 + 2*d3 + d4)
    }

    # --- OBJECTIVE FUNCTIONAL (Simpson's rule, EXACT original formula) ---
    # Original uses 0.5*w1 for the control cost term.
    JE <- E[1] + E[L] + 4*sum(E[seq(2, L-1, by=2)]) + 2*sum(E[seq(3, L-2, by=2)])
    JI <- I[1] + I[L] + 4*sum(I[seq(2, L-1, by=2)]) + 2*sum(I[seq(3, L-2, by=2)])
    JW <- 0.5*w1*(u1[1]^2 + u1[L]^2 + 4*sum(u1[seq(2, L-1, by=2)]^2) + 2*sum(u1[seq(3, L-2, by=2)]^2))
    J <- (h/3)*(JE + JI + JW)
    FV_history[itr]  <- J
    itr_history[itr] <- itr

    # --- CONTROL UPDATE (EXACT original formula) ---
    # Original: ustar <- pmin(zeta, pmax(0, (I * (lambda3 - lambda4)) / w1))
    #           u1 <- 0.5 * (ustar + oldu1)
    # The 0.5 averaging with the old control is critical for stability.
    ustar <- pmin(zeta, pmax(0, (I * (lambda3 - lambda4)) / w1))
    u1 <- 0.5 * (ustar + oldu1)

    # --- CONVERGENCE CHECK (EXACT original formula) ---
    # Original: delta*sum(|x|) - sum(|old_x - x|) > 0 for all variables
    temp1  <- delta_conv*sum(abs(u1)) - sum(abs(oldu1 - u1))
    temp2  <- delta_conv*sum(abs(S)) - sum(abs(oldS - S))
    temp3  <- delta_conv*sum(abs(E)) - sum(abs(oldE - E))
    temp4  <- delta_conv*sum(abs(I)) - sum(abs(oldI - I))
    temp5  <- delta_conv*sum(abs(T_v)) - sum(abs(oldT - T_v))
    temp6  <- delta_conv*sum(abs(R_v)) - sum(abs(oldR - R_v))
    temp7  <- delta_conv*sum(abs(lambda1)) - sum(abs(oldlambda1 - lambda1))
    temp8  <- delta_conv*sum(abs(lambda2)) - sum(abs(oldlambda2 - lambda2))
    temp9  <- delta_conv*sum(abs(lambda3)) - sum(abs(oldlambda3 - lambda3))
    temp10 <- delta_conv*sum(abs(lambda4)) - sum(abs(oldlambda4 - lambda4))
    temp11 <- delta_conv*sum(abs(lambda5)) - sum(abs(oldlambda5 - lambda5))

    test <- min(temp1, temp2, temp3, temp4, temp5, temp6,
                temp7, temp8, temp9, temp10, temp11)
    cat(itr, "   ", test, "\n")
  }

  cat("ODE optimal control converged after", itr, "iterations.\n")

  return(list(
    time = t,
    S = S, E = E, I = I, T_state = T_v, R = R_v, N = N,
    S_uncontrolled = S1, E_uncontrolled = E1, I_uncontrolled = I1,
    T_uncontrolled = T1, R_uncontrolled = R1, N_uncontrolled = N1,
    u1 = u1,
    iterations = itr_history,
    objective_history = FV_history,
    R0 = R0_val,
    final_objective = tail(FV_history, 1)
  ))
}
