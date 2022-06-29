###########################################################
# Outcome and exposure misclassification ##################
###########################################################

#Run once: devtools::install_github("sarahlotspeich/auditDesignR", ref = "main")
library(auditDesignR)

# Set sample sizes ----------------------------------------
N <- 10000 # Phase-I = N
n <- 400 # Phase-II/audit size = n

# Misclassification rates ---------------------------------
fpr_Ystar <- 0.1; tpr_Ystar <- 0.9
fpr_Xstar <- 0.1; tpr_Xstar <- 0.9

eta <- vector() # Nuisance parameters eta

# True parameter values for P(Y|X) --------------------------
pY <- 0.3
eta[1] <- log(pY / (1 - pY)) 
beta <- 0.3

# True parameter value for P(X) -----------------------------
pX <- 0.1 
eta[9] <- log(pX / (1 - pX)) 

# Parameters for error model P(X*|X,Y) ----------------------
eta[6] <- - log((1 - fpr_Xstar) / fpr_Xstar) 
eta[8] <- - log((1 - tpr_Xstar) / tpr_Xstar) - eta[6] 
eta[7] <- (beta + 0.15) 

# Parameters for error model P(Y*|X*,Y,X) -------------------
eta[2] <- - log((1 - fpr_Ystar) / fpr_Ystar)
eta[4] <- - log((1 - tpr_Ystar) / tpr_Ystar) - eta[2] 
eta[3] <- eta[5] <- (beta + 0.25) / 2

set.seed(918)

# Generate Phase I data -------------------------------------
X <- rbinom(n = N, size = 1, prob = sigmoid(eta[9]))
Y <- rbinom(n = N, size = 1, prob = sigmoid(eta[1] + beta * X))
Xstar <- rbinom(n = N, size = 1, prob = sigmoid(eta[6] + eta[7] * Y + eta[8] * X))
Ystar <- rbinom(n = N, size = 1, prob = sigmoid(eta[2] + eta[3] * Xstar + eta[4] * Y + eta[5] * X))
sim_dat <- data.frame(Y, X, Ystar, Xstar)

# Cross-tabulate Phase I strata -----------------------------
N00 <- table(Ystar, Xstar)[1,1]
N01 <- table(Ystar, Xstar)[1,2]
N10 <- table(Ystar, Xstar)[2,1]
N11 <- table(Ystar, Xstar)[2,2]
stratN <- list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)

# Design 1: SRS ---------------------------------------------
V_srs <- sample_srs(phI = N, 
                    phII = n)
mle_srs <- twophase_mle(dat = cbind(V = V_srs, sim_dat), 
                        Y_val = "Y", 
                        Y_unval = "Ystar", 
                        X_val = "X", 
                        X_unval = "Xstar", 
                        Validated = "V")
beta_srs <- mle_srs$mod_Y_val$Est[2]

# Design 2: CC* 
V_cc <- sample_cc(dat = sim_dat, 
                  phI = N, 
                  phII = n, 
                  sample_on = "Ystar")
mle_cc <- twophase_mle(dat = cbind(V = V_cc, sim_dat), 
                       Y_val = "Y", 
                       Y_unval = "Ystar", 
                       X_val = "X", 
                       X_unval = "Xstar", 
                       Validated = "V")
beta_cc <- mle_cc$mod_Y_val$Est[2]

# Design 3: BCC* 
V_bcc <- sample_bcc(dat = sim_dat, 
                    phI = N, 
                    phII = n, 
                    sample_on = c("Ystar", "Xstar"))
mle_bcc <- twophase_mle(dat = cbind(V = V_bcc, sim_dat), 
                        Y_val = "Y", 
                        Y_unval = "Ystar", 
                        X_val = "X", 
                        X_unval = "Xstar", 
                        Validated = "V")
beta_bcc <- mle_bcc$mod_Y_val$Est[2]

# Design 4: optMLE
## Construct complete data and calculate score vectors for each 
complete_data <- expand.grid(Y = c(0, 1), 
                             X = c(0, 1), 
                             Ystar = c(0, 1), 
                             Xstar = c(0, 1), 
                             V = c(0, 1))
s <- score(comp_dat = complete_data, 
           Y_val = "Y", 
           Y_unval = "Ystar", 
           X_val = "X", 
           X_unval = "Xstar", 
           Validated = "V", 
           beta = beta, 
           eta = eta)
## Use adaptive grid search to find optMLE 
grid_search <- optMLE_grid(phI = N, 
                           phII = n, 
                           phI_strat = stratN, 
                           min_n = 10, 
                           sample_on = c("Ystar", "Xstar"), 
                           indiv_score = s, 
                           return_full_grid = FALSE)
if (grid_search$findOptimal) {
  opt_des <- grid_search$min_var_design
  ### Apply optMLE to select subjects 
  V_optMLE <- sample_optMLE(dat = sim_dat, 
                            sample_on = c("Ystar", "Xstar"), 
                            des = opt_des)
  mle_optMLE <- twophase_mle(dat = cbind(V = V_optMLE, sim_dat), 
                             Y_val = "Y", 
                             Y_unval = "Ystar", 
                             X_val = "X", 
                             X_unval = "Xstar", 
                             Validated = "V")
  beta_optMLE <- mle_optMLE$mod_Y_val$Est[2]
} 

# Design 5: optMLE-2
## Wave 1 (na = n / 2)
V_wave1 <- sample_bcc(dat = sim_dat, 
                      phI = N, 
                      phII = (n / 2), 
                      sample_on = c("Ystar", "Xstar"))
wave1_strat <- data.frame(n00 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[1, 1],
                          n01 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[1, 2],
                          n10 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[2, 1],
                          n11 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[2, 2])
### Estimate model parameters using Wave 1 data 
mle_wave1 <- twophase_mle(dat = cbind(V = V_wave1, sim_dat), 
                          Y_val = "Y", 
                          Y_unval = "Ystar", 
                          X_val = "X", 
                          X_unval = "Xstar", 
                          Validated = "V")
beta_wave1 <- mle_wave1$mod_Y_val$Est[2]
eta_wave1 <- with(mle_wave1, c(mod_Y_val$Est[1], mod_Y_unval$Est, mod_X_unval$Est, mod_X_val$Est))
### Design next wave of audits using Wave 1 parameters
s_hat <- score(comp_dat = complete_data, Y_val = "Y", Y_unval = "Ystar", X_val = "X", X_unval = "Xstar", Validated = "V", beta = beta_hat, eta = eta_hat)
grid_search <- optMLE_grid(phI = N, phII = (n / 2), phI_strat = stratN, phIIa_strat = wave1_strat, min_n = 0, sample_on = c("Ystar", "Xstar"), indiv_score = s_hat)
if (grid_search$findOptimal) {
  opt_des2 <- grid_search$min_var_design
  ### Subtract subjects already sampled in Wave 1
  opt_des2[, c("n00", "n01", "n10", "n11")] <- opt_des2[, c("n00", "n01", "n10", "n11")] - with(wave1_strat, c(n00, n01, n10, n11))
  V_wave1 <- sample_optMLE(dat = cbind(V = V_wave1, sim_dat), 
                           sample_on = c("Ystar", "Xstar"), 
                           des = opt_des2, 
                           wave1_Validated = "V")
  V <- pmax(V_wave1, V_wave2)
  ### Estimate model parameters using Waves 1 and 2 data 
  mle_wave2 <- twophase_mle(dat = cbind(V = V, sim_dat), 
                            Y_val = "Y", 
                            Y_unval = "Ystar", 
                            X_val = "X", 
                            X_unval = "Xstar", 
                            Validated = "V")
  ### Final parameter estimate for \beta
  beta_optMLE2 <- mle_wave2$mod_Y_val$Est[2]
} 

# Design 6: optMLE-3
## Wave 1 (na = n / 2)
V_wave1 <- sample_bcc(dat = sim_dat, 
                      phI = N, 
                      phII = (n / 2), 
                      sample_on = c("Ystar", "Xstar"))
wave1_strat <- data.frame(n00 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[1, 1],
                          n01 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[1, 2],
                          n10 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[2, 1],
                          n11 = with(sim_dat[V_wave1 == 1, ], 
                                     table(Ystar, Xstar))[2, 2])
### Estimate model parameters using Wave 1 data 
mle_wave1 <- twophase_mle(dat = cbind(V = V_wave1, sim_dat), 
                          Y_val = "Y", 
                          Y_unval = "Ystar", 
                          X_val = "X", 
                          X_unval = "Xstar", 
                          Validated = "V")
beta_wave1 <- mle_wave1$mod_Y_val$Est[2]
eta_wave1 <- with(mle_wave1, 
                  c(mod_Y_val$Est[1], 
                    mod_Y_unval$Est, 
                    mod_X_unval$Est, 
                    mod_X_val$Est))
### Design next wave of audits using Wave 1 parameters
s_wave1 <- score(comp_dat = complete_data, 
                 Y_val = "Y", 
                 Y_unval = "Ystar", 
                 X_val = "X", 
                 X_unval = "Xstar", 
                 Validated = "V", 
                 beta = beta_wave1, 
                 eta = eta_wave1)
grid_search <- optMLE_grid(phI = N, 
                           phII = (n / 4), 
                           phI_strat = stratN, 
                           phIIa_strat = wave1_strat, 
                           min_n = 0, 
                           sample_on = c("Ystar", "Xstar"), 
                           indiv_score = s_wave1)
## Wave 2 (nb = n / 4)
if (grid_search$findOptimal) {
  opt_des2 <- grid_search$min_var_design
  ### Subtract subjects already sampled in Wave 1
  opt_des2[, c("n00", "n01", "n10", "n11")] <- opt_des2[, c("n00", "n01", "n10", "n11")] - with(wave1_strat, c(n00, n01, n10, n11))
  V_wave2 <- sample_optMLE(dat = cbind(V = V_wave1, sim_dat), 
                           sample_on = c("Ystar", "Xstar"), 
                           des = opt_des2, 
                           wave1_Validated = "V")
  wave2_strat <- data.frame(n00 = with(sim_dat[V_wave2 == 1, ], 
                                       table(Ystar, Xstar))[1, 1],
                            n01 = with(sim_dat[V_wave2 == 1, ], 
                                       table(Ystar, Xstar))[1, 2],
                            n10 = with(sim_dat[V_wave2 == 1, ], 
                                       table(Ystar, Xstar))[2, 1],
                            n11 = with(sim_dat[V_wave2 == 1, ], 
                                       table(Ystar, Xstar))[2, 2])
  V <- pmax(V_wave1, V_wave2)
  ### Estimate model parameters using Waves 1 and 2 data 
  mle_wave2 <- twophase_mle(dat = cbind(V = V, sim_dat), 
                            Y_val = "Y", Y_unval = "Ystar", 
                            X_val = "X", X_unval = "Xstar", 
                            Validated = "V")
  beta_wave2 <- mle_wave2$mod_Y_val$Est[2]
  eta_wave2 <- with(mle_wave2, c(mod_Y_val$Est[1], mod_Y_unval$Est, mod_X_unval$Est, mod_X_val$Est))
  ### Design next wave of audits using Wave 1 and 2 parameters
  s_wave2 <- score(comp_dat = complete_data, 
                   Y_val = "Y", 
                   Y_unval = "Ystar", 
                   X_val = "X", 
                   X_unval = "Xstar", 
                   Validated = "V", 
                   beta = beta_wave2, 
                   eta = eta_wave2)
  grid_search <- optMLE_grid(phI = N, 
                             phII = (n / 4), 
                             phI_strat = stratN, 
                             phIIa_strat = wave1_strat + wave2_strat, 
                             min_n = 0, 
                             sample_on = c("Ystar", "Xstar"), 
                             indiv_score = s_wave2)
  
  ## Wave 3 (nc = n / 4)
  if (grid_search$findOptimal) {
    opt_des3 <- grid_search$min_var_design
    ### Subtract subjects already sampled in Waves 1 and 2
    opt_des3[, c("n00", "n01", "n10", "n11")] <- opt_des3[, c("n00", "n01", "n10", "n11")] - with(wave1_strat, c(n00, n01, n10, n11))
    opt_des3[, c("n00", "n01", "n10", "n11")] <- opt_des3[, c("n00", "n01", "n10", "n11")] - with(wave2_strat, c(n00, n01, n10, n11))
    V_wave3 <- sample_optMLE(dat = cbind(V = V, sim_dat), 
                             sample_on = c("Ystar", "Xstar"), 
                             des = opt_des2, 
                             wave1_Validated = "V")
    V <- pmax(V_wave1, V_wave2, V_wave3)
    ### Estimate model parameters using Waves 1, 2, and 3 data 
    mle_wave3 <- twophase_mle(dat = cbind(V = V, sim_dat), 
                              Y_val = "Y", Y_unval = "Ystar", 
                              X_val = "X", X_unval = "Xstar", 
                              Validated = "V")
    ### Final parameter estimate for \beta
    beta_optMLE3 <- mle_wave3$mod_Y_val$Est[2]
  } 
} 
