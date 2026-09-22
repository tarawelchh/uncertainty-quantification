library(fields) # for rdist()
source("scripts/setup.R")
set.seed(1)

history_match_2D <- function(var1, var2, index1, index2, z) {
  ### Define wave 1 run locations: set random seed above ###
  xD_w1 <- lhd_maximin(nl = 14)
  ### Define run locations ###
  n_run <- 14
  xD <- matrix(rep(0.5, n_run * 5), nrow = n_run, ncol = 5)
  colnames(xD) <- input_names
  lhd_points <- lhd_maximin(n_run)
  xD[, var1] <- lhd_points[, 1]
  xD[, var2] <- lhd_points[, 2]
  xD_1_2_w1 <- cbind(xD[, var1], xD[, var2])
  D <- friedman_matrix(xD)

  ### Define 50x50 grid of prediction points xP for emulator evaluation ###
  x_vals <- seq(0, 1, length.out = 50)
  grid_2D <- setNames(expand.grid(x_vals, x_vals), c(var1, var2))
  xP <- matrix(rep(0.5, 2500 * 2), nrow = 2500, ncol = 2)
  colnames(xP) <- c(var1, var2)
  xP[, var1] <- grid_2D[, var1]
  xP[, var2] <- grid_2D[, var2]
  xP_full <- matrix(rep(0.5, 2500 * 5), nrow = 2500, ncol = 5)
  xP_full[, index1] <- xP[, var1]
  xP_full[, index2] <- xP[, var2]

  sigma <- sd(D)
  E_f <- mean(D)

  ### Evaluate emulator over 201 prediction points xP ###
  sigma_e <- 0.1
  sigma_epsilon <- 0.01

  em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD_1_2_w1, D = D, theta = 0.4, sigma = sd(D), E_f = mean(D)))
  E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
  Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))

  ### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
  Imp_mat <- sqrt((E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2))
  Imp_mat[Imp_mat > 40] <- 40
  ### Define colours and levels for implausibility plots ###
  imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 18, 2), 20, 30, 41)

  ### plot wave 1 implausibility and wave 1 runs only ###
  emul_fill_cont_V2(
    cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = xD, x_grid = x_vals,
    xD_col = "purple", color.palette = imp_cols, main = "Implausibility I(x): Wave 1"
  )

  ### plot wave 1 emulator expectation ###
  emul_fill_cont_V2(
    cont_mat = E_D_fx_mat, cont_levs = NULL, xD = xD, x_grid = x_vals,
    color.palette = exp_cols, main = "Emulator Adjusted Expectation E_D[f(x)]"
  )
  ### plot wave 1 emulator variance ###
  emul_fill_cont_V2(
    cont_mat = Var_D_fx_mat, cont_levs = NULL, xD = xD, x_grid = x_vals,
    color.palette = var_cols, main = "Emulator Adjusted Variance Var_D[f(x)]"
  )

  # --- Step 1: Create Candidates from Implausibility Matrix ---
  grid_points <- setNames(expand.grid(x_vals, x_vals), c(var1, var2))
  grid_points$Imp <- as.vector(Imp_mat)

  # Keep only the "Green" points (Non-Implausible)
  candidates <- grid_points[grid_points$Imp <= 3, c(var1, var2)]

  # --- Step 2: The "Exclusion Zone" (New Logic) ---
  # Calculate distance matrix between Candidates (rows) and Wave 1 points (cols)
  dists <- rdist(candidates, xD_1_2_w1)

  # Find the minimum distance from each candidate to the CLOSEST Wave 1 point
  min_dists <- apply(dists, 1, min)

  # Define a "Safety Radius"
  # (e.g., 5% of the domain width to ensure we don't stack points)
  radius <- 0.05

  # Filter: Keep only candidates that are far away from Wave 1
  survivors <- candidates[min_dists > radius, ]

  # --- Step 3: Run K-Means on Survivors ---
  n_wave2 <- 10

  if (nrow(survivors) < n_wave2) {
    print("Warning: Exclusion radius too large! Taking all survivors.")
    xD_wave2 <- survivors
  } else {
    # Now K-means finds centers that are spaced out from EACH OTHER
    # And Step 2 ensured they are spaced out from WAVE 1
    km_result <- kmeans(survivors, centers = n_wave2, nstart = 30)
    xD_wave2 <- km_result$centers
  }

  # Format output
  xD_1_2_w2 <- as.matrix(xD_wave2)
  ### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
  emul_fill_cont_V2(
    cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = rbind(xD_1_2_w1, xD_1_2_w2),
    x_grid = x_vals, xD_col = rep(c("purple", "pink"), c(nrow(xD_1_2_w1), nrow(xD_1_2_w2))),
    color.palette = imp_cols, main = "Implausibility I(x): Wave 1"
  )


  ### loop over adding the wave 2 runs: add k runs ###
  for (k in 0:10) { # k=0: wave 1, k>0 add k wave 2 runs sequentially

    xD_1_2 <- xD_1_2_w1
    if (k > 0) xD_1_2 <- rbind(xD_1_2, xD_1_2_w2[1:k, ]) # k=0: wave 1, k>0 add wave 2 runs sequentially

    ### Perform 14 + k runs of model and store as D (would take days for realistic example!) ###
    xD <- matrix(rep(0.5, (n_run + k) * 5), nrow = (n_run + k), ncol = 5)
    colnames(xD) <- input_names
    xD[, c(var1, var2)] <- xD_1_2
    D <- friedman_matrix(xD)
    sigma <- sd(D)
    E_f <- mean(D)

    ### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
    em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD_1_2, D = D, theta = 0.4, sigma, E_f))
    E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
    Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))

    ### Evaluate true function and store in matrix for diagnostic comparisons ###
    fxP_mat <- matrix(friedman_matrix(xP_full), nrow = length(x_vals), ncol = length(x_vals))

    ### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
    Imp_mat <- sqrt((E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2))
    Imp_mat[Imp_mat > 40] <- 40
    ### Calculate Imp Measure for True f(x) Over All 50x50 = 2500 input points in xP ###
    Imp_true_mat <- sqrt((fxP_mat - z)^2 / (sigma_e^2 + sigma_epsilon^2))
    Imp_true_mat[Imp_true_mat > 40] <- 40

    ### Define colours and levels for implausibility plots ###
    imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 18, 2), 20, 30, 41)

    ### if k=0 plot wave 1 runs only ###
    if (k == 0) {
      emul_fill_cont_V2(
        cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = xD_1_2,
        x_grid = x_vals, xD_col = "purple", color.palette = imp_cols,
        main = "Implausibility I(x): Wave 1"
      )
    }
    ### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
    emul_fill_cont_V2(
      cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = rbind(xD_1_2, xD_1_2_w2),
      x_grid = x_vals, xD_col = rep(c("purple", "pink"), c(nrow(xD_1_2) + k, nrow(xD_1_2_w2) - k)), # cover unevaluated w2 points in pink
      color.palette = imp_cols, main = "Implausibility I(x): Wave 2"
    )

    ### once last run done so k=8, plot implausibility for true function f(x) to compare ###
    if (k == nrow(xD_1_2_w2)) {
      emul_fill_cont_V2(
        cont_mat = Imp_true_mat, cont_levs = imp_levs,
        cont_levs_lines = 3, xD = xD_1_2, x_grid = x_vals, xD_col = "purple", plot_xD = FALSE,
        color.palette = imp_cols, main = "Implausibility I(x) using True f(x)"
      )
    }
  }

  emul_fill_cont(
    cont_mat = fxP_mat, cont_levs = seq(7.5, 17.5, 0.5), xD = xD_1_2, x_grid = x_vals,
    color.palette = exp_cols, plot_xD = FALSE, main = "True Computer Model Function f(x)"
  )


  em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD_1_2, D = D, theta = 0.4, sigma = sd(D), E_f = mean(D)))
  E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
  Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))


  emul_fill_cont_V2(
    cont_mat = E_D_fx_mat, cont_levs = NULL, xD = xD, x_grid = x_vals,
    color.palette = exp_cols, main = "Emulator Adjusted Expectation E_D[f(x)]"
  )
  ### plot wave 1 emulator variance ###
  emul_fill_cont_V2(
    cont_mat = Var_D_fx_mat, cont_levs = NULL, xD = xD, x_grid = x_vals,
    color.palette = var_cols, main = "Emulator Adjusted Variance Var_D[f(x)]"
  )
}

history_match_2D("x1", "x2", 1, 2, 14.5)
history_match_2D("x1", "x3", 1, 3, 14.5)
