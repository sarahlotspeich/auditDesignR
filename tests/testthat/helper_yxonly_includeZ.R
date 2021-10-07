library(testthat)
library(auditDesignR)

# Set sample sizes ----------------------------------------
N <- 10000 # Phase-I = N
n <- 400 # Phase-II/audit size = n

# Misclassification rates ---------------------------------
fpr_Ystar <- 0.25; tpr_Ystar <- 0.75
fpr_Xstar <- 0.25; tpr_Xstar <- 0.75

eta <- vector() # Nuisance parameters eta

# Parameter values for P(Z) -----------------------------------
pZ <- 0.25

# True parameter values for P(Y|X,Z) --------------------------
pY <- 0.3
eta[1] <- log(pY / (1 - pY))
beta <- 0.3
eta[2] <- - 0.25

# True parameter value for P(X|Z) -----------------------------
pX <- 0.1
eta[12] <- log(pX / (1 - pX))
eta[13] <- 0.5

# Parameters for error model P(X*|X,Y,Z) ----------------------
eta[8] <- - log((1 - fpr_Xstar) / fpr_Xstar)
eta[10] <- - log((1 - tpr_Xstar) / tpr_Xstar) - eta[8]
eta[9] <- (beta + 0.15)
eta[11] <- 1

# Parameters for error model P(Y*|X*,Y,X) -------------------
eta[3] <- - log((1 - fpr_Ystar) / fpr_Ystar)
eta[5] <- - log((1 - tpr_Ystar) / tpr_Ystar) - eta[3]
eta[4] <- eta[6] <- (beta + 0.25) / 2
eta[7] <- 1

set.seed(918)

complete_data <- expand.grid(Y = c(0, 1), X = c(0, 1), Ystar = c(0, 1), Xstar = c(0, 1), Z = c(0, 1), V = c(0, 1))
s <- score(comp_dat = complete_data, Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", addl_covar = "Z", Validated = "V", beta = beta, eta = eta)

# Generate Phase I data -------------------------------------
Z <- rbinom(n = N, size = 1, prob = sigmoid(log(pZ / (1 - pZ))))
X <- rbinom(n = N, size = 1, prob = sigmoid(eta[12] + eta[13] * Z))
Y <- rbinom(n = N, size = 1, prob = sigmoid(eta[1] + eta[2] * Z + beta * X))
Xstar <- rbinom(n = N, size = 1, prob = sigmoid(eta[8] + eta[9] * Y + eta[10] * X + eta[11] * Z))
Ystar <- rbinom(n = N, size = 1, prob = sigmoid(eta[3] + eta[4] * Xstar + eta[5] * Y + eta[6] * X * eta[7] * Z))
sim_dat <- data.frame(Y, X, Ystar, Xstar, Z)

# Cross-tabulate Phase I strata -----------------------------
N000 <- with(sim_dat, sum(Ystar == 0 & Xstar == 0 & Z == 0))
N010 <- with(sim_dat, sum(Ystar == 0 & Xstar == 1 & Z == 0))
N100 <- with(sim_dat, sum(Ystar == 1 & Xstar == 0 & Z == 0))
N110 <- with(sim_dat, sum(Ystar == 1 & Xstar == 1 & Z == 0))
N001 <- with(sim_dat, sum(Ystar == 0 & Xstar == 0 & Z == 1))
N011 <- with(sim_dat, sum(Ystar == 0 & Xstar == 1 & Z == 1))
N101 <- with(sim_dat, sum(Ystar == 1 & Xstar == 0 & Z == 1))
N111 <- with(sim_dat, sum(Ystar == 1 & Xstar == 1 & Z == 1))
stratN <- list(N000 = N000, N010 = N010, N100 = N100, N110 = N110, N001 = N001, N011 = N011, N101 = N101, N111 = N111)
