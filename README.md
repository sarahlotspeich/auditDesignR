# Optimal Multi-Wave Validation for Secondary Use Data with Outcome and Exposure Misclassification
## Lotspeich, Amorim, Shaw, Tao, and Shepherd
The complete R package `auditDesignR` and code for the simulation settings included in the paper. 

### Install
To install the package, run the following in your `R` console: `devtools::install_github("sarahlotspeich/auditDesignR", ref = "main")`.

### Simulation settings 
Inside the `simulations` subdirectory, you will find the following: 

#### Table 1: Simulation results under outcome and exposure misclassification

Data were generated for a Phase I sample of N = 10000 subjects according to equation (1). True X and Y were generated from Bernoulli distributions with $p_{x} = P(X=1)$ and $P(Y=1|X) = [1 + \exp\{-(\beta_0 + 0.3 X)\}]^{-1}$. The approximate outcome prevalence $p_{y} = P(Y=1|X=0)$ was used to define $\beta_0 = \log\{p_{y} / (1 - p_{y})\}$. Error-prone $Y^*$ and $X^*$ were generated from Bernoulli distributions with $P(X^*=1|Y,X) = [1 + \exp\{-(\gamma_0 + 0.45 Y + \gamma_1 X )\}]^{-1}$ and $P(Y^*=1|X^*,Y,X) = [1 + \exp\{-(\alpha_0 + 0.275X^* + \alpha_1 Y + 0.275 X)\}]^{-1}$, where $(\gamma_0, \gamma_1)$ and $(\alpha_0, \alpha_1)$ control the strength of the relationship between error-prone and error-free values. We define the "baseline" false positive and true positive rates for $X^*$, denoted $FPR_{0}(X^*)$ and $TPR_{0}(X^*)$, respectively, as the false positive and true positive rates of $X^*$ when $Y=0$. Similarly, $FPR_{0}(Y^*)$ and $TPR_{0}(Y^*)$ are the false positive and true positive rates for $Y^*$ when $X=X^*=0$.  With these definitions, we have $\alpha_{0} = -\log\left\{\frac{1 - FPR_{0}(Y^*)}{FPR_{0}(Y^*)}\right\}$, $\alpha_1 =  -\log\left\{\frac{1 - TPR_{0}(Y^*)}{TPR_{0}(Y^*)}\right\} - \alpha_0$, $\gamma_{0} = -\log\left\{\frac{1 - FPR_{0}(X^*)}{FPR_{0}(X^*)}\right\}$, and $\gamma_1 = -\log\left\{\frac{1 - TPR_{0}(X^*)}{TPR_{0}(X^*)}\right\} - \gamma_0$. Using the designs in Section 3.1, $n = 400$ subjects were selected in Phase II. 

File: `Table1_SimSetup.R`

#### Table 2: Simulation results under outcome and exposure misclassification with available error-free covariate information

Following equation (1), data were generated for a Phase I sample of N = 10000 subjects. Error-free binary covariate Z was generated from a Bernoulli distribution with P(Z = 1) = 0.25, 0.5. True X and Y were generated from Bernoulli distributions with P(X=1|Z) = 1/[1 + exp{-(-2.2 + 0.5Z)}] and P(Y = 1|X, Z)=1/[1 + exp{-(-0.85 + 0.3X + bz Z)}], for bz = -0.25, 0, 0.25. Baseline misclassification rates were fixed at FPR0 = 0.25 and TPR0 = 0.75 such that X* and Y^* were generated from Bernoulli distributions with P(X*=1|Y,X,Z) = 1/[1 + \exp\{-(-1.1 + 0.45 Y + 2.2 X + λZ)\}] and P(Y*=1|X*,Y,X,Z) = 1/[1 + \exp\{-(-1.1 + 0.275X* + 2.2 Y + 0.275 X + λZ)\}], where λ = -1, 0, 1. In Phase II, n = 400 subjects were selected by extensions of Section 3.1 to sample on (Y*, X*, Z).
