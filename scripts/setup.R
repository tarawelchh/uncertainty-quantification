library(viridis)
library(viridisLite)
friedman <- function(xx) {
    ##########################################################################
    #
    # FRIEDMAN FUNCTION
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUT:
    #
    # xx = c(x1, x2, x3, x4, x5)
    #
    ##########################################################################

    x1 <- xx[1]
    x2 <- xx[2]
    x3 <- xx[3]
    x4 <- xx[4]
    x5 <- xx[5]

    term1 <- 10 * sin(pi * x1 * x2)
    term2 <- 20 * (x3 - 0.5)^2
    term3 <- 10 * x4
    term4 <- 5 * x5

    y <- term1 + term2 + term3 + term4
    return(y)
}

friedman_matrix <- function(X) {
    return(apply(X, 1, friedman))
}

input_names <- c("x1", "x2", "x3", "x4", "x5")


##############################################################################################
### Define simple Bayes Linear emulator for single input ###
simple_BL_emulator_v1 <- function(x, # the emulator prediction point
                                  xD, # the run input locations xD
                                  D, # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1, # the correlation lengths
                                  sigma = 1, # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0 # prior expectation of f: E(f(x)) = 0
) {
    # store length of runs D
    n <- length(D)

    ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
    Cov_fx_fxdash <- function(x, xdash) sigma^2 * exp(-(x - xdash)^2 / theta^2)


    ### Define 5 objects needed for BL adjustment ###
    # Create E[D] vector
    E_D <- rep(E_f, n)

    # Create Var_D matrix:
    Var_D <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) for (j in 1:n) Var_D[i, j] <- Cov_fx_fxdash(xD[i], xD[j])

    # Create E[f(x)]
    E_fx <- E_f

    # Create Var_f(x)
    Var_fx <- sigma^2

    # Create Cov_fx_D row vector
    Cov_fx_D <- matrix(0, nrow = 1, ncol = n)
    for (j in 1:n) Cov_fx_D[1, j] <- Cov_fx_fxdash(x, xD[j])


    ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
    ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)
    VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)

    ### return emulator expectation and variance ###
    return(c("ExpD_f(x)" = ED_fx, "VarD_f(x)" = VarD_fx))
}
### End simple Bayes Linear emulator for single input ###
##############################################################################################

##############################################################################################
### Define simple Bayes Linear emulator for single input in 2D ###

simple_BL_emulator_v2 <- function(x, # the emulator prediction point
                                  xD, # the run input locations xD
                                  D, # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1, # the correlation lengths
                                  sigma = 1, # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0 # prior expectation of f: E(f(x)) = 0
) {
    # store length of runs D
    n <- length(D)

    ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
    # Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2)    # XXX Old 1D version
    Cov_fx_fxdash <- function(x, xdash) sigma^2 * exp(-sum((x - xdash)^2) / theta^2) # XXX New 2D version


    ### Define 5 objects needed for BL adjustment ###
    # Create E[D] vector
    E_D <- rep(E_f, n)

    # Create Var_D matrix:
    Var_D <- matrix(0, nrow = n, ncol = n)
    # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])  # XXX Old 1D version
    for (i in 1:n) for (j in 1:n) Var_D[i, j] <- Cov_fx_fxdash(xD[i, ], xD[j, ]) # XXX New 2D version

    # Create E[f(x)]
    E_fx <- E_f

    # Create Var_f(x)
    Var_fx <- sigma^2

    # Create Cov_fx_D row vector
    Cov_fx_D <- matrix(0, nrow = 1, ncol = n)
    # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])    # XXX Old 1D version
    for (j in 1:n) Cov_fx_D[1, j] <- Cov_fx_fxdash(x, xD[j, ]) # XXX New 2D version


    ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
    ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)
    VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)

    ### return emulator expectation and variance ###
    return(c("ExpD_f(x)" = ED_fx, "VarD_f(x)" = VarD_fx))
}
### End simple Bayes Linear emulator for single input in 2D ###
##############################################################################################

##############################################################################################
### Define more efficient Bayes Linear emulator for Multiple inputs in 2D ###

simple_BL_emulator_v3 <- function(xP, # the set of emulator prediction points
                                  xD, # the run input locations xD
                                  D, # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1, # the correlation lengths (can be a vector)
                                  sigma = 1, # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0, # prior expectation of f: E(f(x)) = 0
                                  using_pdist = 1 # if you have installed pdist package
) {
    # store length of runs D and number of prediction points xP
    n <- length(D)
    nP <- nrow(xP) # XXX New V3

    # # Rescale each input by theta. Works for different theta for each input and for same theta
    xP <- t(t(xP) / theta) # XXX New V3: solution to Exercise 8.3, CP 3.
    xD <- t(t(xD) / theta) # XXX New V3: solution to Exercise 8.3, CP 3.

    ### Define Cov structure of f(x): Cov[f(x),f(xdash)], now to act on matrix of distances ###
    # Cov_fx_fxdash <- function(x,xdash) sig^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX Old V2
    Cov_fx_fxdash <- function(dist_matrix) sigma^2 * exp(-(dist_matrix)^2) # XXX New dist V3


    ### Define 5 objects needed for BL adjustment ###
    # Create E[D] vector
    E_D <- rep(E_f, n)

    # Create Var_D matrix:
    # Var_D <- matrix(0,nrow=n,ncol=n)                                        # XXX Old V2
    # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX Old V2
    Var_D <- Cov_fx_fxdash(as.matrix(dist(xD))) # XXX New dist V3

    # Create E[f(x)]
    E_fx <- rep(E_f, nP)

    # Create Var_f(x)
    Var_fx <- rep(sigma^2, nP)

    # Create Cov_fx_D row vector now using pdist() function if available, if not use dist()
    # Cov_fx_D <- matrix(0,nrow=1,ncol=n)                       # XXX Old V2
    # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX Old V2
    if (using_pdist) Cov_fx_D <- Cov_fx_fxdash(as.matrix(pdist(xP, xD))) # XXX New V3
    if (!using_pdist) {
        Cov_fx_D <- Cov_fx_fxdash(as.matrix(dist(rbind(xP, xD)))[1:nP, (nP + 1):(nP + n)])
    } # XXX NewV3

    # find inverse of Var_D using Cholesky decomposition (check Wikipedia if interested!)
    Var_D_inv <- chol2inv(chol(Var_D)) # more efficient way to calculate inverse of Cov mat

    ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
    cov_fx_D_Var_D_inv <- Cov_fx_D %*% Var_D_inv # Need this twice so pre-calculate here
    ED_fx <- E_fx + cov_fx_D_Var_D_inv %*% (D - E_D) # adj. expectation of ALL xP pts at once
    # VarD_fx <-  Var_fx - cov_fx_D_Var_D_inv %*% t(Cov_fx_D)       # XXX Old V2
    VarD_fx <- Var_fx - apply(cov_fx_D_Var_D_inv * Cov_fx_D, 1, sum) # fast way to get diagonals
    # and hence all nP variances (Note: does not do full np x np covariance matrix)

    ### return emulator expectation and variance ###
    return(cbind("ExpD_f(x)" = c(ED_fx), "VarD_f(x)" = VarD_fx))
}
### End simple Bayes Linear emulator for single input in 2D ###
##############################################################################################


##############################################################################################
### Function to plot simple emulator output
plot_BL_emulator_V1 <- function(
  em_out, # a nP x 2 matrix of emulator outputs
  xP, # vector of nP inputs where emulator evaluated
  xD, # the run input locations
  D, # the run outputs D = (f(x^1),...,f(x^n))
  maintitle = NULL # title of the plot, blank if NULL
) {
    ### plot emulator output ###
    plot(xP,
        em_out[, "ExpD_f(x)"],
        ylim = c(13.5, 21),
        ty = "l", col = "blue", lwd = 2.5,
        xlab = expression("x3"),
        ylab = "Emulator Output",
        main = maintitle
    ) # main argument: plot title
    lines(xP, em_out[, "ExpD_f(x)"] + 3 * sqrt(em_out[, "VarD_f(x)"]), col = "red", lwd = 2.5)
    lines(xP, em_out[, "ExpD_f(x)"] - 3 * sqrt(em_out[, "VarD_f(x)"]), col = "red", lwd = 2.5)

    ### plot true function: we would not normally be able to do this! ###
    xP_full <- matrix(rep(0.5, 1005), nrow = 201, ncol = 5)
    colnames(xP_full) <- input_names
    xP_full[, "x3"] <- xP
    lines(xP, friedman_matrix(xP_full), lwd = 2, lty = 1)

    ### Plot the runs ###
    points(xD, D, pch = 21, col = 1, bg = "green", cex = 1.5)
    legend("top",
        legend = c(
            "Emulator Expectation", "Emulator Prediction Interval",
            "True function f(x)", "Model Evaluations"
        ),
        lty = c(1, 1, 1, NA, 2), pch = c(NA, NA, NA, 16, NA), col = c("blue", "red", 1, "green", 1),
        lwd = c(2.5, 2.5, 2.5, 2.5, 1), pt.cex = 1.3
    )
}
### Function to plot simple emulator output
##############################################################################################

##############################################################################################
### Function to plot simple emulator output
plot_BL_emulator_V2 <- function(
  em_out, # a nP x 2 matrix of emulator outputs
  xP, # vector of nP inputs where emulator evaluated
  xD, # the run input locations
  D, # the run outputs D = (f(x^1),...,f(x^n))
  maintitle = NULL, # title of the plot, blank if NULL
  z = NULL, # the observed data z
  sigma_e = NULL, # the observation errors SD s.t. Var[e] = sigma_e^2
  sigma_epsilon = NULL, # the model discrepancy SD s.t. Var[epsilon] = sigma_epsilon^2
  plot_true = FALSE # don't plot true function unless this is TRUE
) {
    ### plot emulator output ###
    plot(xP, em_out[, "ExpD_f(x)"],
        ylim = c(7, 30), ty = "l", col = "blue", lwd = 2.5,
        xlab = "Input parameter x", ylab = "Output f(x)", main = maintitle
    ) # main argument: plot title
    lines(xP, em_out[, "ExpD_f(x)"] + 3 * sqrt(em_out[, "VarD_f(x)"]), col = "red", lwd = 2.5)
    lines(xP, em_out[, "ExpD_f(x)"] - 3 * sqrt(em_out[, "VarD_f(x)"]), col = "red", lwd = 2.5)

    ### plot true function: we would not normally be able to do this! ###
    if (plot_true) lines(xP, friedman_matrix(xP), lwd = 2, lty = 1)

    ### Plot the runs ###
    points(xD, D, pch = 21, col = 1, bg = "green", cex = 1.5)

    ### Plot the Observed data plus errors due to obs and MD ###
    if (!is.null(z)) {
        abline(h = z, lwd = 1.4)
        abline(h = z + 3 * sqrt(sigma_e^2 + sigma_epsilon^2), lty = 2, lwd = 1.2)
        abline(h = z - 3 * sqrt(sigma_e^2 + sigma_epsilon^2), lty = 2, lwd = 1.2)
    }

    if (plot_true) {
        legend("top",
            legend = c(
                "Emulator Expectation",
                "Emulator Prediction Interval",
                "True function f(x)", "Model Evaluations"
            ),
            lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 16), col = c("blue", "red", 1, "green"), lwd = 2.5, pt.cex = 1.3
        )
    }
    if (!is.null(z)) {
        legend("top",
            legend = c(
                "Emulator Expectation",
                "Emulator Prediction Interval",
                "Model Evaluations", "Observation z",
                "3 sigma interval"
            ),
            lty = c(1, 1, NA, 1, 2), pch = c(NA, NA, 16, NA, NA), col = c("blue", "red", "green", 1, 1),
            lwd = c(2.5, 2.5, NA, 1.4, 1.2), pt.cex = 1.3
        )
    }
}
### Function to plot simple emulator output
##############################################################################################

## MAKE LATIN HYPERCUBE##
lhd_maximin <- function(nl = 16) { # nl = number of points in LHD

    x_lhd <- cbind("x1" = sample(0:(nl - 1)), "x2" = sample(0:(nl - 1))) / nl + 0.5 / nl # create LHD

    ### Maximin loop: performs swaps on 1st of two closest points with another random point
    for (i in 1:1000) {
        mat <- as.matrix(dist(x_lhd)) + diag(10, nl) # creates matrix of distances between points
        # note the inflated diagonal
        closest_runs <- which(mat == min(mat), arr.ind = TRUE) # finds pairs of closest runs
        ind <- closest_runs[sample(nrow(closest_runs), 1), 1] # chooses one of close runs at random
        swap_ind <- sample(setdiff(1:nl, ind), 1) # randomly selects another run to swap with
        x_lhd2 <- x_lhd # creates second version of LHD
        x_lhd2[ind[1], 1] <- x_lhd[swap_ind, 1] # swaps x_1 vals between 1st close run & other run
        x_lhd2[swap_ind, 1] <- x_lhd[ind[1], 1] # swaps x_1 vals between 1st close run & other run
        if (min(dist(x_lhd2)) >= min(dist(x_lhd)) - 0.00001) { # if min distance between points is same or better
            x_lhd <- x_lhd2 # we replace LHD with new LHD with the swap
            # cat("min dist =",min(dist(x_lhd)),"Iteration = ",i,"\n") # write out min dist
        }
    }
    ### plot maximin LHD ###
    plot(x_lhd,
        xlim = c(0, 1), ylim = c(0, 1), pch = 16, xaxs = "i", yaxs = "i", col = "blue",
        xlab = "x1", ylab = "x2", cex = 1.4
    )
    abline(h = (0:nl) / nl, col = "grey60")
    abline(v = (0:nl) / nl, col = "grey60")
    return(x_lhd)
}


### define filled contour plot function for emulator output ###
emul_fill_cont <- function(
  cont_mat, # matrix of values we want contour plot of
  cont_levs = NULL, # contour levels (NULL: automatic selection)
  nlev = 20, # approx no. of contour levels for auto select
  plot_xD = TRUE, # plot the design runs TRUE or FALSE
  xD = NULL, # the design points if needed
  xD_col = "green", # colour of design runs
  x_grid, # grid edge locations that define xP
  ... # extra arguments passed to filled.contour
) {
    ### Define contour levels if necessary ###
    if (is.null(cont_levs)) cont_levs <- pretty(cont_mat, n = nlev)

    ### create the filled contour plot ###
    filled.contour(x_grid, x_grid, cont_mat,
        levels = cont_levs, xlab = "x1", ylab = "x2", ...,
        plot.axes = {
            axis(1)
            axis(2) # sets up plotting in contour box
            contour(x_grid, x_grid, cont_mat, add = TRUE, levels = cont_levs, lwd = 0.8) # plot contour lines
            if (plot_xD) points(xD, pch = 21, col = 1, bg = xD_col, cex = 1.5)
        }
    ) # plot design points
}


### define filled contour plot function for emulator output ###
emul_fill_cont_V2 <- function(
  cont_mat, # matrix of values we want contour plot of
  cont_levs = NULL, # contour levels (NULL: automatic selection)
  cont_levs_lines = NULL, # contour levels for lines (NULL: automatic selection)
  nlev = 20, # approx no. of contour levels for auto select
  plot_xD = TRUE, # plot the design runs TRUE or FALSE
  xD = NULL, # the design points if needed
  xD_col = "green", # colour of design runs
  x_grid, # grid edge locations that define xP
  ... # extra arguments passed to filled.contour
) {
    ### Define contour levels if necessary ###
    if (is.null(cont_levs)) cont_levs <- pretty(cont_mat, n = nlev)

    ### create the filled contour plot ###
    filled.contour(x_grid, x_grid, cont_mat,
        levels = cont_levs, xlab = "x1", ylab = "x2", ...,
        plot.axes = {
            axis(1)
            axis(2) # sets up plotting in contour box
            if (is.null(cont_levs_lines)) contour(x_grid, x_grid, cont_mat, add = TRUE, levels = cont_levs, lwd = 0.8) # plot usual contour lines
            if (!is.null(cont_levs_lines)) {
                contour(x_grid, x_grid, cont_mat, add = TRUE, levels = cont_levs, lwd = 0.4, labels = "") # plot thin contour lines
                contour(x_grid, x_grid, cont_mat, add = TRUE, levels = cont_levs_lines, lwd = 2) # plot thick contour lines
            }
            if (plot_xD) points(xD, pch = 21, col = 1, bg = xD_col, cex = 1.5)
        }
    ) # plot design points
}
### end define filled contour plot function for emulator output ###
##############################################################################################

exp_cols <- magma
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo

imp_cols <- function(n) turbo(n, begin = 0.15, end = 1)
