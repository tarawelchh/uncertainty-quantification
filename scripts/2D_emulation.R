source("scripts/setup.R")
set.seed(123)
n_run <- 20
xD <- matrix(rep(0.5, n_run * 5), nrow = n_run, ncol = 5)
colnames(xD) <- input_names
lhd_points <- lhd_maximin(n_run)
xD[, "x1"] <- lhd_points[, 1]
xD[, "x2"] <- lhd_points[, 2]
xD_1_2 <- cbind(xD[, "x1"], xD[, "x2"])
D <- friedman_matrix(xD)
### Perform 20 runs of model and store as D (this would takes days for realistic example!) ###

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_vals <- seq(0, 1, length.out = 50)
grid_2D <- expand.grid(x1 = x_vals, x2 = x_vals) # This creates 2500 rows
xP <- matrix(rep(0.5, 2500 * 2), nrow = 2500, ncol = 2)
colnames(xP) <- c("x1", "x2")
xP[, "x1"] <- grid_2D[, "x1"]
xP[, "x2"] <- grid_2D[, "x2"]

sigma <- sd(D)
E_f <- mean(D)

em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD_1_2, D = D, theta = 0.3, sigma, E_f))
head(em_out)

### store emulator output as matrices to aid plotting ###
E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))

emul_fill_cont(E_D_fx_mat, xD = xD_1_2, x_grid = x_vals)

theta_seq <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.65)

### loop over vector of theta values ###
for (i in 1:length(theta_seq)) {
  ### Evaluate emulator over 201 prediction points xP and store in matrices ###
  em_out <- t(apply(xP[, 1:2], 1, simple_BL_emulator_v2, xD = xD_1_2, D = D, theta = theta_seq[i], sigma, E_f))
  E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
  Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))

  ### Plot filled contour plot of emulator expectation ###
  emul_fill_cont(
    cont_mat = E_D_fx_mat, xD = xD_1_2, x_grid = x_vals,
    color.palette = magma, main = paste("Emul. Adj. Expectation E_D[f(x)], theta =", theta_seq[i])
  )
}

for (i in 1:length(theta_seq)) {
  em_out <- t(apply(xP[, 1:2], 1, simple_BL_emulator_v2, xD = xD_1_2, D = D, theta = theta_seq[i], sigma, E_f))
  Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_vals), ncol = length(x_vals))
  emul_fill_cont(
    cont_mat = Var_D_fx_mat, xD = xD_1_2, x_grid = x_vals,
    color.palette = magma, main = paste("Emul. Adj. Var Var_D[f(x)], theta =", theta_seq[i])
  )
}
