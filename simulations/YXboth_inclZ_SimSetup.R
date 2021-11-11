###########################################################
# Outcome and exposure misclassification ##################
# with additional error-free covariate ####################
###########################################################

#Run once: devtools::install_github("sarahlotspeich/auditDesignR", ref = "main")
# library(auditDesignR)
R.utils::sourceDirectory("R/")
tic("inclZ")
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

# Design 1: SRS
tic("design 1")
V_srs <- sample_srs(phI = N, phII = n)
mle_srs <- twophase_mle(dat = cbind(V = V_srs, sim_dat), Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", addl_covar = "Z", Validated = "V")
beta_srs <- mle_srs$mod_Y_val$Est[3]
toc()
# design 1: 176.1 sec elapsed

# Design 2: CC*
tic("design 2")
V_cc <- sample_cc(dat = sim_dat, phI = N, phII = n, sample_on = "Ystar")
mle_cc <- twophase_mle(dat = cbind(V = V_cc, sim_dat), Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", addl_covar = "Z", Validated = "V")
beta_cc <- mle_cc$mod_Y_val$Est[3]
toc()
# design 2: 171.67 sec elapsed


# Design 3: BCC*
tic("design 3")
V_bcc <- sample_bcc(dat = sim_dat, phI = N, phII = n, sample_on = c("Ystar", "Xstar", "Z"))
mle_bcc <- twophase_mle(dat = cbind(V = V_bcc, sim_dat), Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", addl_covar = "Z", Validated = "V")
beta_bcc <- mle_bcc$mod_Y_val$Est[3]
toc()
# design 3: 200.55 sec elapsed


# Design 4: optMLE
tic("design 4")
grid_search <- optMLE_grid(phI = N, phII = n, phI_strat = stratN, min_n = 10, sample_on = c("Ystar", "Xstar", "Z"), indiv_score = s)
if (grid_search$findOptimal) {
  opt_des <- grid_search$min_var_design
  V_optMLE <- sample_optMLE(dat = sim_dat, sample_on = c("Ystar", "Xstar", "Z"), des = opt_des)
  mle_optMLE <- twophase_mle(dat = cbind(V = V_optMLE, sim_dat), Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", addl_covar = "Z", Validated = "V")
  beta_optMLE <- mle_optMLE$mod_Y_val$Est[3]
}
toc()
# design 4: 341.81 sec elapsed


# Design 5: optMLE-2
tic("design 5")
V_wave1 <- sample_bcc(dat = sim_dat, phI = N, phII = (n / 2), sample_on = c("Ystar", "Xstar", "Z"))

wave1_strat <- data.frame(n000 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 0 & Z == 0)),
                          n010 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 1 & Z == 0)),
                          n100 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 0 & Z == 0)),
                          n110 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 1 & Z == 0)),
                          n001 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 0 & Z == 1)),
                          n011 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 1 & Z == 1)),
                          n101 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 0 & Z == 1)),
                          n111 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 1 & Z == 1)))
mle_wave1 <- twophase_mle(dat = cbind(V = V_wave1, sim_dat), Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", addl_covar = "Z", Validated = "V")
beta_hat <- mle_wave1$mod_Y_val$Est[3]
eta_hat <- with(mle_wave1, c(mod_Y_val$Est[1:2], mod_Y_unval$Est, mod_X_unval$Est, mod_X_val$Est))
s_hat <- score(comp_dat = complete_data, Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", Validated = "V", addl_covar = "Z", beta = beta_hat, eta = eta_hat)
grid_search <- optMLE_grid(phI = N, phII = (n / 2), phI_strat = stratN, phIIa_strat = wave1_strat, min_n = 0, sample_on = c("Ystar", "Xstar", "Z"), indiv_score = s_hat)
if (grid_search$findOptimal) {
  opt_des2 <- grid_search$min_var_design
  opt_des2[, c("n000", "n010", "n100", "n110", "n001", "n011", "n101", "n111")] <- opt_des2[, c("n000", "n010", "n100", "n110", "n001", "n011", "n101", "n111")] - with(wave1_strat, c(n000, n010, n100, n110, n001, n011, n101, n111))
  V_optMLE2 <- pmax(V_wave1, sample_optMLE(dat = cbind(V = V_wave1, sim_dat), sample_on = c("Ystar", "Xstar", "Z"), des = opt_des2, wave1_Validated = "V"))
  mle_optMLE2 <- twophase_mle(dat = cbind(V = V_optMLE2, sim_dat), Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", addl_covar = "Z", Validated = "V")
  beta_optMLE2 <- mle_optMLE2$mod_Y_val$Est[3]
}
toc()
# 628.58 sec elapsed

toc()

# 1527.81

