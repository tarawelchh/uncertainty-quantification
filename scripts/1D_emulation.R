source("scripts/setup.R")
### Define run locations ###
xD <- matrix(rep(0.5, 30), nrow = 6, ncol = 5)
colnames(xD) <- input_names
xD[, "x3"] <- t(seq(0, 1, length.out = 6))

### Perform 5 runs of model and store as D (this would takes days for realistic example!) ###
D <- friedman_matrix(xD)
xD_x3 <- xD[, "x3"]
E_f <- mean(D)
sigma <- sd(D)

x <- 0.1
em_out_1pt <- simple_BL_emulator_v1(x, xD_x3, D, 0.25, sigma, E_f)

### Prediction Interval at x=0.1 ###
low1 <- em_out_1pt["ExpD_f(x)"] - 3 * sqrt(em_out_1pt["VarD_f(x)"])
up1 <- em_out_1pt["ExpD_f(x)"] + 3 * sqrt(em_out_1pt["VarD_f(x)"])
c(low1, up1, use.names = FALSE)

### Evaluate emulator over 201 prediction points xP ###
xP <- seq(0.001, 0.999, len = 201)
nP <- length(xP) # number of prediction points
em_out <- t(sapply(xP, simple_BL_emulator_v1, xD = xD_x3, D = D, theta = 0.25, sigma, E_f)) # t(): gives transpose of the matrix
head(em_out)

plot_BL_emulator_V1(em_out, xP, xD_x3, D)

### make sequence of sigma values to use for emulation ###
theta_seq <- c(0.01, 0.2, 0.3, 0.4, 0.35, 0.45, 0.65, 0.9)

### for loop over different sigma values in sigma_seq ###
for (i in 1:length(theta_seq)) {
  ### Evaluate emulator over 201 prediction points xP with sigma=sigma_seq[i] ###
  em_out <- t(sapply(xP, simple_BL_emulator_v1, xD = xD_x3, D = D, theta = theta_seq[i], sigma = sigma, E_f = E_f))

  ### Plot emulator output in each case, note use of "paste" for plot title ###
  plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD_x3, D = D, maintitle = paste("Theta =", theta_seq[i]))
}

sigma_seq <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5)

### for loop over different sigma values in sigma_seq ###
for (i in 1:length(sigma_seq)) {
  ### Evaluate emulator over 201 prediction points xP with sigma=sigma_seq[i] ###
  em_out <- t(sapply(xP, simple_BL_emulator_v1, xD = xD_x3, D = D, theta = 0.25, sigma = sigma_seq[i], E_f = E_f))

  ### Plot emulator output in each case, note use of "paste" for plot title ###
  plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD_x3, D = D, maintitle = paste("Sigma =", sigma_seq[i]))
}

em_out <- t(sapply(xP, simple_BL_emulator_v1, xD = xD_x3, D = D, theta = 0.9, sigma = sigma, E_f = E_f))
plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD_x3, D = D, maintitle = paste("Theta =", 0.65))
