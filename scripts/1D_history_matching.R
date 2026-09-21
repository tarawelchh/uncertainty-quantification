source("scripts/setup.R")
plot_1D_Implausibility_V1 <- function(
  Implaus, # the vector of implausibilities to plot
  xP, # the vector of input points where the implausibilities were evaluated
  imp_cutoff = 3, # implausibility cutoff
  maintitle = NULL # title of the plot, blank if NULL
) {
  imp_breaks <- c(seq(0, 2.8, 0.2), seq(3, 9.5, 0.5), seq(10, 80, 10)) # colour steps for imp

  imp_cut_index <- cut(Implaus, imp_breaks, labels = FALSE) # which imp intervals each pt lies in
  cols <- turbo(length(imp_breaks) - 1, begin = 0.15, end = 1) # assign reduced "turbo" colours

  plot(xP, Implaus,
    ylim = c(0, 20), ty = "l", lwd = 2, # plot imp with optional main title
    xlab = "Input parameter x", ylab = "Implausibility I(x)", main = maintitle
  )
  abline(h = imp_cutoff, col = "green", lwd = 2) # draw horizontal implausibility cutoff
  rect_wid <- 1 / (length(xP) - 1) # set up simple rectangles to colour x-axis
  for (i in 1:length(xP)) rect(xP[i] - rect_wid / 2, -1, xP[i] + rect_wid / 2, 0, col = cols[imp_cut_index][i], border = NA)
  abline(h = 0, col = 1, lwd = 1) # additional black line for clarity

  legend("topright",
    legend = c("Implausibility I(x)", "Implausibility Cutoff c"),
    lty = c(1, 1), col = c(1, "green"), lwd = 2
  ) # legend
}

history_match_1D_wave1 <- function(xD_xi, variable_name, z) {
  ### Evaluate emulator over 201 prediction points xP ###
  xP <- seq(0.001, 0.999, len = 201)
  ### Define run locations ###
  xD <- matrix(rep(0.5, 5 * length(xD_xi)), nrow = length(xD_xi), ncol = 5)
  colnames(xD) <- input_names
  xD[, variable_name] <- xD_xi

  ### Perform 6 runs of model and store as D (this would takes days for realistic example!) ###
  D <- friedman_matrix(xD)
  ### Evaluate emulator over 201 prediction points xP ###
  em_out <- t(sapply(xP, simple_BL_emulator_v1, xD = xD_xi, D = D, theta = 0.25, sigma = sd(D), E_f = mean(D)))

  ### Defined Extra Objects for HM and Implausibility ###
  sigma_e <- 0.1
  sigma_epsilon <- 0.02

  ### Plot emulator output with observation errors ###
  plot_BL_emulator_V2(
    em_out = em_out, xP = xP, xD = xD_xi, D = D, maintitle = "Wave 1 Emulator Output: 6 runs",
    z = z, sigma_e = sigma_e, sigma_epsilon = sigma_epsilon
  )

  ### Calculate Implausibility Measure Over All 201 input points in xP ###
  Implaus <- sqrt((em_out[, "ExpD_f(x)"] - z)^2 / (em_out[, "VarD_f(x)"] + sigma_e^2 + sigma_epsilon^2))
  ### Plot Implausibility function ###
  ### Plot the implausbility ###
  Imp_mat[Imp_mat > 40] <- 40
  plot_1D_Implausibility_V1(Implaus = Implaus, xP = xP)
}

##############################################################################################
### 1D HM example wave 2: add three points sequentially ###
history_match_1D_wave2 <- function(xD_xi, xD_xi_w2, variable_name, z) {
  for (k in c(1:2)) { # k=1: plot the emulator, k=2 plot the implausibility
    for (j in 0:3) { # if j==0 do wave 1, if j>0 add j wave 2 points
      ### Define run locations ###
      if (j == 0) {
        current_xD_xi <- xD_xi
      } else {
        current_xD_xi <- c(xD_xi, xD_xi_w2[1:j])
      }

      size <- length(xD_xi) + j
      xD <- matrix(rep(0.5, 5 * size), nrow = size, ncol = 5)
      colnames(xD) <- input_names
      print(xD[, "x1"])
      xD[, variable_name] <- current_xD_xi
      ### Perform 6 runs and store as D (this would takes days for realistic example!) ###
      D <- friedman_matrix(xD)

      ### Evaluate emulator over 201 prediction points xP ###
      em_out <- t(sapply(xP, simple_BL_emulator_v1, xD = current_xD_xi, D = D, theta = 0.25, sigma = sd(D), E_f = mean(D)))

      ### Plot emulator output with observation errors ###
      if (k == 1) {
        plot_BL_emulator_V2(
          em_out = em_out, xP = xP, xD = current_xD_xi, D = D,
          maintitle = paste("Wave", c(1, 2, 2, 2)[j + 1], "Emulator Output:", length(current_xD_xi), "runs"),
          z = z, sigma_e = sigma_e, sigma_epsilon = sigma_epsilon
        )
      }

      ### Calculate Implausibility Measure Over All 201 input points in xP ###
      if (k == 2) {
        Implaus <- sqrt((em_out[, "ExpD_f(x)"] - z)^2 /
          (em_out[, "VarD_f(x)"] + sigma_e^2 + sigma_epsilon^2))
      }

      ### Plot the implausbility ###
      if (k == 2) {
        plot_1D_Implausibility_V1(
          Implaus = Implaus, xP = xP,
          maintitle = paste("Wave", c(1, 2, 2, 2)[j + 1], "Implausibility:", length(current_xD_xi), "runs")
        )
      }
    }
  }
}


xD_x1 <- c(0, 0.2, 0.4, 0.6, 1) # the wave 1 inputs
xD_x1_w2 <- c(0.3, 0.45, 0.9) # the wave 2 inputs
z_x1 <- 14.5

history_match_1D_wave1(xD_x1, "x1", z_x1)
history_match_1D_wave2(xD_x1, xD_x1_w2, "x1", z_x1)

xD_x3 <- c(0, 0.4, 0.6, 0.86, 1)
z_x3 <- 17.021
xD_x3_w2 <- c(0.1, 0.5, 0.8, 0.85)
history_match_1D_wave1(xD_x3, "x3", z_x3)
history_match_1D_wave2(xD_x3, xD_x1_w2, "x3", z_x3)
