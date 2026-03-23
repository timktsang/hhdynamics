# Statistical Methodology

## Overview

This vignette describes the statistical model, likelihood, prior
distributions, and MCMC algorithm implemented in **hhdynamics**. For a
practical workflow guide, see
[`vignette("hhdynamics-intro")`](https://timktsang.github.io/hhdynamics/articles/hhdynamics-intro.md).

The model is a discrete-time Bayesian household transmission model that
separates community and within-household infection sources. It was
originally described in Tsang et al. (2014) and has been extended here
to support covariate effects, missing data imputation, and optional
serial interval estimation.

## Model

### Study design

The data come from a case-ascertained household study. Each household is
enrolled through an index case (the first identified case), and
household contacts are followed prospectively for a defined period. The
observed data for each contact include whether they became infected and,
if so, the day of symptom onset. Covariates (age, sex, etc.) may be
recorded for each individual.

### Infection hazard

For household contact $i$ in household $h$ on day $t$, the hazard of
infection is:

$$\lambda_{i}(t) = \left\lbrack \beta_{c} + \sum\limits_{\substack{j \in \mathcal{I}_{h} \\ j \neq i}}\frac{\beta_{h} \cdot w\left( t - t_{j} \right) \cdot \exp\left( \mathbf{x}_{j}^{\top}{\mathbf{α}}_{\text{inf}} + \varepsilon_{j} \right)}{\left( n_{h} - 1 \right)^{\gamma}} \right\rbrack \cdot \exp\left( \mathbf{z}_{i}^{\top}{\mathbf{α}}_{\text{sus}} \right)$$

where:

- $\beta_{c} > 0$ is the daily rate of infection from the community,
  assumed constant over the follow-up period.
- $\beta_{h} > 0$ is the baseline rate of person-to-person transmission
  within households.
- $w( \cdot )$ is the serial interval distribution: a discrete
  probability mass function where $w(d)$ gives the probability of
  transmission $d$ days after the infector’s symptom onset. By default,
  a 14-day influenza serial interval from Tsang et al. (2014) is used.
- $t_{j}$ is the symptom onset day of infected household member $j$.
- $\mathcal{I}_{h}$ is the set of all infected members in household $h$
  (including the index case and any tertiary cases).
- $\mathbf{x}_{j}$ is the vector of infectivity covariates for person
  $j$, and ${\mathbf{α}}_{\text{inf}}$ are the corresponding log-scale
  effects.
- $\varepsilon_{j} \sim \text{Normal}\left( 0,\sigma_{\varepsilon}^{2} \right)$
  is an optional individual-level random effect on infectivity
  ($\sigma_{\varepsilon}$ is estimated; disabled by default).
- $n_{h}$ is the household size and $\gamma \geq 0$ is a scaling
  parameter (currently fixed at 0, so there is no frequency-dependent
  adjustment).
- $\mathbf{z}_{i}$ is the vector of susceptibility covariates for person
  $i$, and ${\mathbf{α}}_{\text{sus}}$ are the corresponding log-scale
  effects.

The daily probability that contact $i$ becomes infected on day $t$,
given they are still susceptible, is:

$$P\left( {\text{infected on day}\mspace{6mu}}t \mid {\text{susceptible at}\mspace{6mu}}t \right) = 1 - \exp\left( - \lambda_{i}(t) \right)$$

### Covariates

Covariate effects enter the hazard multiplicatively on the log scale:

- **Infectivity covariates** ($\mathbf{x}_{j}$): modify the transmission
  rate of the infector. A positive $\alpha_{\text{inf},k}$ means that
  level increases the hazard of infecting others. The exponentiated
  coefficient $\exp\left( \alpha_{\text{inf},k} \right)$ is the relative
  infectiousness compared to the reference level.
- **Susceptibility covariates** ($\mathbf{z}_{i}$): modify the
  susceptibility of the contact. A positive $\alpha_{\text{sus},k}$
  means that level is more susceptible. The exponentiated coefficient
  $\exp\left( \alpha_{\text{sus},k} \right)$ is the relative
  susceptibility compared to the reference level.

Factor covariates are encoded as dummy variables using treatment
contrasts (the first level is the reference).

## Likelihood

### Individual likelihood

For contact $i$ in household $h$, let $T_{i}$ be the infection day (if
infected) or the end of follow-up plus one (if not infected), and let
$s_{h}$ be the start of the observation period (the index case’s onset
day). Let $\delta_{i} = 1$ if infected, $0$ otherwise.

The log-likelihood contribution of contact $i$ is:

$$\ell_{i} = - \sum\limits_{t = s_{h}}^{T_{i} - 2}\lambda_{i}(t) + \delta_{i} \cdot \log\left( 1 - \exp\left( - \lambda_{i}\left( T_{i} - 1 \right) \right) \right)$$

The first term captures the survival probability (not being infected on
days before $T_{i}$), and the second term captures the probability of
infection on day $T_{i}$ (only included if the contact was actually
infected).

### Full likelihood

The full data log-likelihood is the sum over all non-index contacts in
all households:

$$\ell({\mathbf{θ}}) = \sum\limits_{h = 1}^{H}\sum\limits_{i = 1}^{n_{h} - 1}\ell_{i}$$

where
${\mathbf{θ}} = \left( \beta_{c},\beta_{h},\gamma,{\mathbf{α}}_{\text{inf}},{\mathbf{α}}_{\text{sus}},\sigma_{\varepsilon} \right)$
is the full parameter vector.

When random effects are included, the likelihood also includes the
random effects density:

$$\ell_{\text{RE}} = \sum\limits_{h = 1}^{H}\sum\limits_{j \in \mathcal{I}_{h}}\log\phi\left( \varepsilon_{j};0,\sigma_{\varepsilon}^{2} \right)$$

where $\phi\left( \cdot ;0,\sigma^{2} \right)$ is the normal density.

## Prior distributions

| Parameter               | Prior                      | Description                                                         |
|-------------------------|----------------------------|---------------------------------------------------------------------|
| $\sigma_{\varepsilon}$  | Uniform(0.009, 5)          | RE standard deviation (fixed at 1 when random effects are disabled) |
| $\beta_{c}$             | Uniform($10^{- 18}$, 9.99) | Community infection rate                                            |
| $\beta_{h}$             | Uniform($10^{- 18}$, 9.99) | Household transmission rate                                         |
| $\gamma$                | Uniform(0, 1)              | Household size parameter (currently fixed at 0)                     |
| $\alpha_{\text{inf},k}$ | Normal(0, 3)               | Each infectivity covariate effect                                   |
| $\alpha_{\text{sus},k}$ | Normal(0, 3)               | Each susceptibility covariate effect                                |
| Weibull shape           | Uniform(0.1, 10)           | SI shape (only when `estimate_SI = TRUE`)                           |
| Weibull scale           | Uniform(0.1, 20)           | SI scale (only when `estimate_SI = TRUE`)                           |

The uniform priors on $\beta_{c}$ and $\beta_{h}$ are deliberately wide
and effectively non-informative for realistic transmission rates. The
Normal(0, 3) prior on covariate effects is weakly informative — it
assigns most probability mass to relative risks between
$\exp( - 6) \approx 0.002$ and $\exp(6) \approx 400$, which is very
permissive for epidemiological covariate effects.

## MCMC algorithm

### Metropolis-Hastings sampler

Parameters are updated one at a time via Metropolis-Hastings with
Gaussian random-walk proposals. At iteration $b$:

1.  For each parameter $\theta_{k}$ with a free indicator
    (`move[k] = 1`):
    1.  Propose $\theta_{k}^{*} = \theta_{k}^{(b - 1)} + \epsilon$,
        where
        $\epsilon \sim \text{Normal}\left( 0,\sigma_{k}^{2} \right)$.
    2.  Compute the log acceptance ratio:
        $$r = \left\lbrack \ell\left( {\mathbf{θ}}^{*} \right) + \log\pi\left( {\mathbf{θ}}^{*} \right) \right\rbrack - \left\lbrack \ell\left( {\mathbf{θ}}^{(b - 1)} \right) + \log\pi\left( {\mathbf{θ}}^{(b - 1)} \right) \right\rbrack$$
    3.  Accept with probability $\min\left( 1,\exp(r) \right)$.

If the proposed parameter falls outside the prior support, the prior
log-density is $- \infty$ and the proposal is automatically rejected.

### Adaptive proposal tuning

After the first 500 iterations, proposal standard deviations
$\sigma_{k}$ are set to the empirical posterior standard deviation of
the chain so far, with multiplicative adjustments based on acceptance
rates:

| Acceptance rate | Adjustment               |
|-----------------|--------------------------|
| \< 10%          | $\sigma_{k} \times 0.5$  |
| 10–15%          | $\sigma_{k} \times 0.8$  |
| 15–20%          | $\sigma_{k} \times 0.95$ |
| 20–30%          | No change (target range) |
| 30–40%          | $\sigma_{k} \times 1.05$ |
| 40–90%          | $\sigma_{k} \times 1.2$  |
| \> 90%          | $\sigma_{k} \times 2.0$  |

This adaptive scheme typically achieves stable acceptance rates within
1000–2000 iterations.

### Data augmentation

At each MCMC iteration, after updating all parameters, the following
latent variables are updated for each household member via a joint
Metropolis-Hastings step:

1.  **Infection times**: For infected non-index contacts, a new onset
    time is proposed uniformly over the follow-up window. The proposal
    is accepted or rejected based on the full household likelihood
    ratio.

2.  **Random effects** (when enabled): For each infected individual, a
    new random effect $\varepsilon_{j}^{*}$ is proposed from
    $\varepsilon_{j} + \text{Normal}\left( 0,\sigma_{\varepsilon} \right)$.

3.  **Missing covariates**: For individuals with missing factor
    covariates, a new level is proposed uniformly from all levels of the
    factor. This is a symmetric proposal, so no correction term is
    needed.

All three augmentation steps are proposed jointly for each household
member and accepted or rejected together. This joint update is more
efficient than separate updates because the onset time, random effect,
and covariate value are correlated in the posterior. The augmentation
step is parallelized across households using RcppParallel.

### Serial interval estimation

When `estimate_SI = TRUE`, two additional parameters (Weibull shape and
scale) are added to the MCMC. At each iteration, the serial interval PMF
is recomputed from the current Weibull parameters:

$$w(d) = F_{W}(d + 1;k,\lambda) - F_{W}(d;k,\lambda),\quad d = 1,2,\ldots,14$$

where $F_{W}( \cdot ;k,\lambda)$ is the Weibull CDF with shape $k$ and
scale $\lambda$. The Weibull parameters are updated via the same
Metropolis-Hastings random-walk scheme as the other parameters.

## Convergence diagnostics

Because the sampler uses a single chain with adaptive tuning, standard
multi-chain diagnostics (such as $\widehat{R}$) are not directly
applicable. Instead, we recommend:

1.  **Trace plots**: `plot_diagnostics(fit)` shows trace plots and
    posterior densities for all parameters. Look for stationarity after
    burn-in and good mixing (no long excursions or sticky regions).

2.  **Effective sample size (ESS)**:
    `plot_diagnostics(fit, show_ess = TRUE)` or
    `table_parameters(fit, show_ess = TRUE)` report ESS using the
    initial positive sequence estimator (Geyer, 1992). ESS values below
    100 suggest poor mixing; consider increasing `n_iteration` or
    investigating identifiability.

3.  **Acceptance rates**: `table_parameters(fit)` shows acceptance
    rates. Rates outside 15–40% indicate poor tuning, which may occur if
    the burn-in was insufficient for the adaptive scheme to converge. In
    practice, the default burn-in of 5000 iterations is usually
    sufficient.

4.  **Log-likelihood trace**: `fit$log_likelihood` stores the full
    log-likelihood trace (before thinning) for visual inspection.

### Recommended MCMC settings

For most datasets:

- `n_iteration = 50000`: sufficient for convergence with
  moderate-dimensional models.
- `burnin = 10000`: allows the adaptive tuning to stabilize.
- `thinning = 10`: reduces autocorrelation in stored samples.

For datasets with many covariates, missing data, or
`estimate_SI = TRUE`, consider longer chains (100,000+ iterations) and
examine ESS for each parameter.

## Interpretation of parameters

### Base parameters

- `community`: Reported as a daily probability on the probability scale
  via $1 - \exp\left( - \beta_{c} \right)$. This is the daily
  probability of infection from community sources, independent of
  household transmission.
- `household`: Reported as a daily probability via
  $1 - \exp\left( - \beta_{h} \right)$. This is the baseline daily
  probability of person-to-person transmission within a household from a
  single infected member at peak serial interval.

### Covariate effects

Covariate effects are estimated on the log scale. The
[`summary()`](https://rdrr.io/r/base/summary.html) method reports both
the raw (log-scale) estimates and the exponentiated estimates:

- $\exp\left( \alpha_{k} \right) > 1$: This level has higher
  infectivity/susceptibility than the reference.
- $\exp\left( \alpha_{k} \right) < 1$: This level has lower
  infectivity/susceptibility than the reference.
- $\exp\left( \alpha_{k} \right) = 1$: No difference from the reference
  level.

For example, if the susceptibility effect for age group “6-17” is
$\widehat{\alpha} = - 0.3$ with 95% CrI $( - 0.8,0.2)$, then
$\exp\left( \widehat{\alpha} \right) = 0.74$ with 95% CrI $(0.45,1.22)$.
This means children aged 6-17 have an estimated 26% lower susceptibility
than the reference group, but the credible interval includes 1 (no
difference).

### Serial interval parameters

When `estimate_SI = TRUE`, the output includes `si_shape` and `si_scale`
— the Weibull shape ($k$) and scale ($\lambda$) parameters. The mean
serial interval is $\lambda \cdot \Gamma(1 + 1/k)$ and the variance is
$\lambda^{2}\left\lbrack \Gamma(1 + 2/k) - \Gamma(1 + 1/k)^{2} \right\rbrack$.

## Computational details

- The likelihood and data augmentation are parallelized across
  households using **RcppParallel** (Intel TBB backend), providing
  significant speedups for large datasets.
- Forward simulation
  ([`simulate_data()`](https://timktsang.github.io/hhdynamics/reference/simulate_data.md))
  is forced to single-threaded mode because the R random number
  generator is not thread-safe.
- The C++ backend uses **RcppArmadillo** for linear algebra operations.

## References

Tsang TK, Lau EHY, Cauchemez S, Cowling BJ. Association between antibody
titers and protection against influenza virus infection within
households. *Journal of Infectious Diseases*. 2014;210(5):684–692.
<https://doi.org/10.1093/infdis/jiu186>

Geyer CJ. Practical Markov chain Monte Carlo. *Statistical Science*.
1992;7(4):473–483.
