# Optimal Multi-Wave Validation for Secondary Use Data with Outcome and Exposure Misclassification
## Lotspeich, Amorim, Shaw, Tao, and Shepherd
The complete R package `auditDesignR` and code for the simulation settings included in the paper. 

### Install
To install the package, run the following in your `R` console: `devtools::install_github("sarahlotspeich/auditDesignR", ref = "main")`.

### Validation Study Designs

  - *Simple random sampling (SRS):* All subjects in Phase I have equal probability of inclusion in Phase II.
  - *Unvalidated case-control sampling (CC*):* Subjects are stratified on Y* and separate random samples of size n/2 are drawn from each stratum.
  - *Unvalidated balanced case-control sampling (BCC*):* Subjects are jointly stratified on (Y*, X*) and separate random samples of size n/4 subjects drawn from each stratum. 
  - *Optimal design (optMLE):* Subjects are jointly stratified on (Y*,X*), and stratum sizes are chosen following Section 2.2. The optMLE is included as a "gold standard" design since it requires knowing the parameters θ.
  - *Two-wave approximate optimal design (optMLE-2):* Subjects are jointly stratified on (Y*,X*). In the first wave, n/2 subjects are selected using BCC*, and in the second wave the remaining subjects are chosen following the design in Section 2.4.

### Simulation settings 
Inside the `simulations` subdirectory, you will find the following: 

  - `YXboth_SimSetup.R`: simulations with outcome and exposure misclassification, defining sampling strata on error-prone outcome and exposure, Y* and X* (Table 1 in the Manuscript)
  - `YXboth_inclZ_SimSetup.R`: simulations with outcome and exposure misclassification, defining sampling strata on error-prone outcome and exposure, Y* and X*, and error-free binary covariate Z (Table S2 in the Supplemental Materials for the Manuscript).
  - `YXboth_inclZ_misspec_SimSetup.R`: 
  - `Yonly_SimSetup.R`: simulations with outcome misclassification alone, defining sampling strata on error-prone outcome Y* and error-free exposure X (top of Table S3 in the Supplemental Materials for the Manuscript).
  - `Xonly_SimSetup.R`: simulations with exposure misclassification alone, defining sampling strata on error-free outcome Y and error-prone exposure X* (top of Table S3 in the Supplemental Materials for the Manuscript).
