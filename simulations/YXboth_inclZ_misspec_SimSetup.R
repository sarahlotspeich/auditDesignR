###########################################################
# Outcome and exposure misclassification ##################
# with additional error-free covariate ####################
# possible model misspecification at design stage #########
###########################################################

#Run once: devtools::install_github("sarahlotspeich/auditDesignR", ref = "main")
library(auditDesignR)
tic("misspec")
# Set sample sizes ----------------------------------------
N <- 10000 # Phase-I = N
n <- 400 # Phase-II/audit size = n

# Coefficient on interaction between X and Z --------------
delta <- 0.25

# Misclassification rates ---------------------------------
fpr_Ystar <- 0.25; tpr_Ystar <- 0.75
fpr_Xstar <- 0.25; tpr_Xstar <- 0.75

eta <- vector() # Nuisance parameters eta
# e <- vector() # NEEDED TO ADD TO MAKE WORK
# Parameter values for P(Z) -----------------------------------
pZ <- 0.25

# True parameter values for P(Y|X,Z) --------------------------
pY <- 0.3
eta[1] <- log(pY / (1 - pY))
beta <- 0.3
eta[2] <- 0.25

# Parameters for error model P(Y*|X*,Y,X,Z,XZ) ----------------
e[3] <- - log((1 - fpr_Ystar) / fpr_Ystar)
e[5] <- - log((1 - tpr_Ystar) / tpr_Ystar) - e[3]
e[4] <- e[6] <- (beta + 0.25) / 2
e[7] <- 1
e[8] <- delta

# Parameters for error model P(X*|Y,X,Z) ----------------------
e[9] <- - log((1 - p0_Xstar) / p0_Xstar)
e[11] <- - log((1 - p1_Xstar) / p1_Xstar) - e[9]
e[10] <- (beta + 0.15)
e[12] <- 1
e[13] <- delta

# True parameter value for P(X|Z) -----------------------------
pX <- 0.1
e[14] <- log(pX / (1 - pX))
e[15] <- 0.5

set.seed(918)

# Nuisance parameters for the optMLE and optMLE-2 designs ---
eta_int <- e
# Nuisance parameters for the optMLE* and optMLE-2* designs -
eta_me <- e[-c(8, 13)]

# Generate Phase I data -------------------------------------
Z <- rbinom(n = N, size = 1, prob = sigmoid(log(0.25 / (1 - 0.25))))
X <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[12] + eta_me[13] * Z))
XZ <- X * Z
Y <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[1] + eta_me[2] * Z + beta * X))
misspec <- "Ystar+Xstar" # Interactions in both error models
#misspec <- "Ystar" # Interactions in only Y* error model
#misspec <- "Xstar" # Interactions in only X* error model
if (misspec == "Ystar+Xstar") {
  Xstar <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[8] + eta_me[9] * Y + eta_me[10] * X + eta_me[11] * Z + delta * XZ))
  Ystar <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[3] + eta_me[4] * Xstar + eta_me[5] * Y + eta_me[6] * X * eta_me[7] * Z + delta * XZ))
} else if (misspec == "Ystar") {
  Xstar <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[8] + eta_me[9] * Y + eta_me[10] * X + eta_me[11] * Z))
  Ystar <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[3] + eta_me[4] * Xstar + eta_me[5] * Y + eta_me[6] * X * eta_me[7] * Z + delta * XZ))
} else if (misspec == "Xstar") {
  Xstar <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[8] + eta_me[9] * Y + eta_me[10] * X + eta_me[11] * Z  + delta * XZ))
  Ystar <- rbinom(n = N, size = 1, prob = sigmoid(eta_me[3] + eta_me[4] * Xstar + eta_me[5] * Y + eta_me[6] * X * eta_me[7] * Z))
}
sim_dat <- data.frame(Y, X, Ystar, Xstar, Z, XZ)

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
V_srs <- sample_srs(phI = N, phII = n)
mle_srs <- twophase_mle_int(dat = cbind(V = V_srs, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                            X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                            Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                            int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
beta_srs <- mle_srs$mod_Y_val$Est[3]

# Design 2: CC*
V_cc <- sample_cc(dat = sim_dat, phI = N, phII = n, sample_on = "Ystar")
mle_cc <- twophase_mle_int(dat = cbind(V = V_cc, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                           X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                           Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                           int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
beta_cc <- mle_cc$mod_Y_val$Est[3]

# Design 3: BCC*
V_bcc <- sample_bcc(dat = sim_dat, phI = N, phII = n, sample_on = c("Ystar", "Xstar", "Z"))
mle_bcc <- twophase_mle_int(dat = cbind(V = V_bcc, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                            X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                            Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                            int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
beta_bcc <- mle_bcc$mod_Y_val$Est[3]

# Design 4: optMLE
complete_data <- expand.grid(Y = c(0, 1), X = c(0, 1), Ystar = c(0, 1), Xstar = c(0, 1), Z = c(0, 1), V = c(0, 1))
s <- score_w_int(comp_dat = complete_data, Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar",
                 addl_covar = "Z", Validated = "V", beta = beta, eta = eta_int,
                 int_Y_unval = grepl("Ystar", misspec), int_X_unval = grepl("Xstar", misspec))
grid_search <- optMLE_grid(phI = N, phII = n, phI_strat = stratN, min_n = 10, sample_on = c("Ystar", "Xstar", "Z"), indiv_score = s)
if (grid_search$findOptimal) {
  opt_des <- grid_search$min_var_design
  V_optMLE <- sample_optMLE(dat = sim_dat, sample_on = c("Ystar", "Xstar", "Z"), des = opt_des)
  mle_optMLE <- twophase_mle_int(dat = cbind(V = V_optMLE, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                                 X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                                 Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                                 int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
  beta_optMLE <- mle_optMLE$mod_Y_val$Est[3]
}

# Design 5: optMLE*
s <- score(comp_dat = complete_data, Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar",
           addl_covar = "Z", Validated = "V", beta = beta, eta = eta_me)
grid_search <- optMLE_grid(phI = N, phII = n, phI_strat = stratN, min_n = 10, sample_on = c("Ystar", "Xstar", "Z"), indiv_score = s)
if (grid_search$findOptimal) {
  opt_des <- grid_search$min_var_design
  V_optMLEstar <- sample_optMLE(dat = sim_dat, sample_on = c("Ystar", "Xstar", "Z"), des = opt_des)
  mle_optMLEstar <- twophase_mle_int(dat = cbind(V = V_optMLEstar, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                                     X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                                     Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                                     int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
  beta_optMLEstar <- mle_optMLE$mod_Y_val$Est[3]
}

# Design 6: optMLE-2
V_wave1 <- sample_bcc(dat = sim_dat, phI = N, phII = (n / 2), sample_on = c("Ystar", "Xstar", "Z"))
wave1_strat <- data.frame(n000 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 0 & Z == 0)),
                          n010 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 1 & Z == 0)),
                          n100 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 0 & Z == 0)),
                          n110 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 1 & Z == 0)),
                          n001 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 0 & Z == 1)),
                          n011 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 0 & Xstar == 1 & Z == 1)),
                          n101 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 0 & Z == 1)),
                          n111 = with(sim_dat[V_wave1 == 1, ], sum(Ystar == 1 & Xstar == 1 & Z == 1)))
mle_wave1 <- twophase_mle_int(dat = cbind(V = V_wave1, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                              X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                              Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                              int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
beta_hat <- mle_wave1$mod_Y_val$Est[3]
eta_hat <- with(mle_wave1, c(mod_Y_val$Est[1:2], mod_Y_unval$Est, mod_X_unval$Est, mod_X_val$Est))
s_hat <- score_w_int(comp_dat = complete_data, Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar",
                     addl_covar = "Z", Validated = "V", beta = beta_hat, eta = eta_hat,
                     int_Y_unval = grepl("Ystar", misspec), int_X_unval = grepl("Xstar", misspec))
grid_search <- optMLE_grid(phI = N, phII = (n / 2), phI_strat = stratN, phIIa_strat = wave1_strat, min_n = 0, sample_on = c("Ystar", "Xstar", "Z"), indiv_score = s_hat)
if (grid_search$findOptimal) {
  opt_des2 <- grid_search$min_var_design
  opt_des2[, c("n000", "n010", "n100", "n110", "n001", "n011", "n101", "n111")] <- opt_des2[, c("n000", "n010", "n100", "n110", "n001", "n011", "n101", "n111")] - with(wave1_strat, c(n000, n010, n100, n110, n001, n011, n101, n111))
  V_optMLE2 <- pmax(V_wave1, sample_optMLE(dat = cbind(V = V_wave1, sim_dat), sample_on = c("Ystar", "Xstar", "Z"), des = opt_des2, wave1_Validated = "V"))
  mle_optMLE2 <- twophase_mle_int(dat = cbind(V = V_optMLE2, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                                  X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                                  Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                                  int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
  beta_optMLE2 <- mle_optMLE2$mod_Y_val$Est[3]
}

# Design 7: optMLE-2*
beta_hat <- mle_wave1$mod_Y_val$Est[3]
eta_hat <- with(mle_wave1, c(mod_Y_val$Est[1:2], mod_Y_unval$Est, mod_X_unval$Est, mod_X_val$Est))
s_hat_me <- score_w_int(comp_dat = complete_data, Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar",
                        addl_covar = "Z", Validated = "V", beta = beta_hat, eta = eta_hat,
                        int_Y_unval = FALSE, int_X_unval = FALSE)
grid_search <- optMLE_grid(phI = N, phII = (n / 2), phI_strat = stratN, phIIa_strat = wave1_strat, min_n = 0, sample_on = c("Ystar", "Xstar", "Z"), indiv_score = s_hat)
if (grid_search$findOptimal) {
  opt_des2 <- grid_search$min_var_design
  opt_des2[, c("n000", "n010", "n100", "n110", "n001", "n011", "n101", "n111")] <- opt_des2[, c("n000", "n010", "n100", "n110", "n001", "n011", "n101", "n111")] - with(wave1_strat, c(n000, n010, n100, n110, n001, n011, n101, n111))
  V_optMLE2star <- pmax(V_wave1, sample_optMLE(dat = cbind(V = V_wave1, sim_dat), sample_on = c("Ystar", "Xstar", "Z"), des = opt_des2, wave1_Validated = "V"))
  mle_optMLE2star <- twophase_mle_int(dat = cbind(V = V_optMLE2star, sim_dat), Y_val = "Y", Y_unval = "Ystar",
                                  X_val = "X", X_unval = "Xstar", addl_covar = "Z", interact = "XZ",
                                  Validated = "V", int_Y_unval = grepl(pattern = "Ystar", misspec),
                                  int_X_unval = grepl(pattern = "Xstar", misspec), noSE = FALSE)
  beta_optMLE2star <- mle_optMLE2star$mod_Y_val$Est[3]
}

toc()
