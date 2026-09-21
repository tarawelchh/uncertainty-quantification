source("scripts/setup.R")
set.seed(1441)
library(pdist)
library(fields)

n_run <- 12
xD <- matrix(rep(0, n_run * 2), nrow = n_run, ncol = 2)
xD_full <- matrix(rep(1, n_run * 5), nrow = n_run, ncol = 5)
colnames(xD_full) <- input_names
lhd_points <- lhd_maximin(n_run)
xD[, 1] <- lhd_points[, 1]
xD[, 2] <- lhd_points[, 2]
xD_full[, "x1"] <- lhd_points[, 1]
xD_full[, "x2"] <- lhd_points[, 2]
xD_w1 <- xD
D <- friedman_matrix(xD_full)

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_vals <- seq(0, 1, length.out = 50)
grid_2D <- expand.grid(x1 = x_vals, x2 = x_vals) # This creates 2500 rows
xP <- matrix(rep(1, 2500 * 2), nrow = 2500, ncol = 2)
colnames(xP) <- c("x1", "x2")
xP[, "x1"] <- grid_2D[, "x1"]
xP[, "x2"] <- grid_2D[, "x2"]
xP_full <- matrix(rep(1, 2500 * 5), nrow = 2500, ncol = 5)
xP_full[, 1:2] <- xP

sigma <- 0.1
E_f <- mean(D)


### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- simple_BL_emulator_v3(xP = xP, xD = xD, D = D, theta = 0.45, sigma, mean(D))
E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))

### Defined Best Run output so far = f_plus and model discrepancy ###
f_plus <- max(D)
sigma_epsilon <- 0.2

### Evaluate true function f(x) and store in matrix for (possible) diag. comparisons ###
fxP_mat <- matrix(friedman_matrix(xP_full), nrow = length(x_vals), ncol = length(x_vals))
### find true max (we pretend we don't know this!) ###
f_plus_true <- max(fxP_mat)

### plot emulator expectation and true function ###
emul_fill_cont(
  cont_mat = E_D_fx_mat, cont_levs = NULL, xD = xD, x_grid = x_vals,
  color.palette = exp_cols, main = "Emulator Adjusted Expectation E_D[f(x)]"
)

emul_fill_cont(
  cont_mat = fxP_mat, cont_levs = NULL, xD = xD, x_grid = x_vals,
  color.palette = exp_cols, plot_xD = FALSE, main = "True Computer Model Function f(x)"
)

### Calculate new Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat <- (f_plus - E_D_fx_mat) / sqrt(Var_D_fx_mat + sigma_epsilon^2)
Imp_mat[Imp_mat < 0] <- 0.0001 # replace negatives with small value for plotting purposes
Imp_mat[Imp_mat > 40] <- 40

### Calculate true Implausibility Measure Over All 201 input points in xP ###
Imp_true_mat <- (f_plus_true - fxP_mat) / sqrt(sigma_epsilon^2)
Imp_true_mat[Imp_true_mat > 40] <- 40 # limit high implausibilities for plotting purposes

imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 18, 2), 20, 30, 41)

### Plot New Optimisation Implausibility ###
emul_fill_cont_V2(
  cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3,
  xD = xD, x_grid = x_vals, xD_col = "purple", color.palette = imp_cols,
  main = "Max Implausibility I_M(x)"
)

### Plot New Optimisation Implausibility for true function with true max ###
emul_fill_cont_V2(
  cont_mat = Imp_true_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = xD, x_grid = x_vals,
  xD_col = "purple", plot_xD = FALSE, color.palette = imp_cols,
  main = "Implausibility I(x) using True f(x)"
)


# create candidates from implausibility matrix
grid_points <- expand.grid(x1 = x_vals, x2 = x_vals)
grid_points$Imp <- as.vector(Imp_mat)
# keep non-implausible points
candidates <- grid_points[grid_points$Imp <= 3, c("x1", "x2")]
# calculate distance matrix between candidates (rows) and wave 1 points (cols)
dists <- rdist(candidates, xD)
# find the min dist from each candidate to the closest Wave 1 point
min_dists <- apply(dists, 1, min)


radius <- 0.05
survivors <- candidates[min_dists > radius, ] # dont want points on top

# run kmeans on survivors
n_wave2 <- 7

if (nrow(survivors) < n_wave2) {
  print("radius of exclusion too large - no survivors")
  xD_wave2 <- survivors
} else {
  # kmeans finds centres that are spaced out from EACH OTHER
  # above radius ensured they are spaced out from WAVE 1
  km_result <- kmeans(survivors, centers = n_wave2, nstart = 30)
  xD_wave2 <- km_result$centers
}

xD_wave2 <- as.matrix(xD_wave2)

### loop over adding the wave 2 runs: add k runs ###
for (k in 0:7) { # k=0: wave 1, k>0 add k wave 2 runs sequentially
  xD <- xD_w1
  if (k > 0) xD <- rbind(xD, xD_wave2[1:k, ]) # k=0: wave 1, k>0 add wave 2 runs sequentially

  ### Perform 14 + k runs of model and store as D (would take days for realistic example!) ###
  xD_full <- matrix(rep(1, (n_run + k) * 5), nrow = (n_run + k), ncol = 5)
  xD_full[, 1:2] <- xD
  D <- friedman_matrix(xD_full)
  sigma <- sd(D)
  E_f <- mean(D)

  ### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
  em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD, D = D, theta = 0.45, sigma, E_f))
  E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
  Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))

  ### Evaluate true function and store in matrix for diagnostic comparisons ###
  fxP_mat <- matrix(friedman_matrix(xP_full), nrow = length(x_vals), ncol = length(x_vals))

  ### Calculate new Implausibility Measure Over All 50x50 = 2500 input points in xP ###
  Imp_mat <- (f_plus - E_D_fx_mat) / sqrt(Var_D_fx_mat + sigma_epsilon^2)
  Imp_mat[Imp_mat < 0] <- 0.0001 # replace negatives with small value for plotting purposes
  Imp_mat[Imp_mat > 40] <- 40

  ### Calculate true Implausibility Measure Over All 201 input points in xP ###
  Imp_true_mat <- (f_plus_true - fxP_mat) / sqrt(sigma_epsilon^2)
  Imp_true_mat[Imp_true_mat > 40] <- 40 # limit high implausibilities for plotting purposes


  ### Define colours and levels for implausibility plots ###
  imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 18, 2), 20, 30, 41)

  ### if k=0 plot wave 1 runs only ###
  if (k == 0) {
    emul_fill_cont_V2(
      cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = xD,
      x_grid = x_vals, xD_col = "purple", color.palette = imp_cols,
      main = "Implausibility I(x): Wave 1"
    )
  }
  ### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
  emul_fill_cont_V2(
    cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = rbind(xD, xD_wave2),
    x_grid = x_vals, xD_col = rep(c("purple", "pink"), c(nrow(xD) + k, nrow(xD_wave2) - k)), # cover unevaluated w2 points in pink
    color.palette = imp_cols, main = "Implausibility I(x): Wave 2"
  )
  ### once last run done so k=8, plot implausibility for true function f(x) to compare ###
  if (k == nrow(xD_wave2)) {
    emul_fill_cont_V2(
      cont_mat = Imp_true_mat, cont_levs = imp_levs,
      cont_levs_lines = 3, xD = xD, x_grid = x_vals, xD_col = "purple", plot_xD = FALSE,
      color.palette = imp_cols, main = "Implausibility I(x) using True f(x)"
    )
  }
}
