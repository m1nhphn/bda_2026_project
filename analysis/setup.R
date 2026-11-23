options(repos = c(CRAN = "https://cloud.r-project.org"))
if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)
}
if(!require(tidyr)){
    install.packages("tidyr")
    library(tidyr)
}
if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
}
if(!require(patchwork)){
    install.packages("patchwork")
    library(patchwork)
}
if(!require(posterior)){
    install.packages("posterior")
    library(posterior)
}
if(!require(bayesplot)){
    install.packages("bayesplot")
    library(bayesplot)
}
if(!require(ggdist)){
    install.packages("ggdist")
    library(ggdist)
}
if(!require(loo)){
    install.packages("loo")
    library(loo)
}
if(!require(priorsense)){
    install.packages("priorsense")
    library(priorsense)
}
if(!require(knitr)){
    install.packages("knitr")
    library(knitr)
}
if(!require(brms)){
    install.packages("brms")
    library(brms)
}
if(!require(cmdstanr)){
    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    library(cmdstanr)
}
cmdstan_installed <- function(){
  res <- try(out <- cmdstanr::cmdstan_path(), silent = TRUE)
  !inherits(res, "try-error")
}
if(!cmdstan_installed()){
    install_cmdstan()
}

ggplot2::theme_set(theme_minimal(base_size = 14))
bayesplot::bayesplot_theme_set(theme_minimal(base_size = 14))
options(brms.backend = "cmdstanr")
SEED <- 2026
set.seed(SEED)
dir.create("fits", showWarnings = FALSE)

cigar <- read.csv("Cigar.csv") |>
  select(-rownames) |>
  mutate(
    state = factor(state),
    year_full = 1900 + year,
    rprice = price / cpi * 100,
    rincome = ndi / cpi * 100,
    rpimin = pimin / cpi * 100,
    log_sales = log(sales),
    log_rprice = log(rprice),
    log_rincome = log(rincome),
    log_rpimin = log(rpimin),
    log_rprice_c = log_rprice - mean(log_rprice),
    log_rincome_c = log_rincome - mean(log_rincome),
    log_rpimin_c = log_rpimin - mean(log_rpimin),
    year_c = (year_full - mean(year_full)) / 10
  )

specs <- tibble(
  model = c("M1a", "M1b", "M2a", "M2b"),
  controls = c(FALSE, TRUE, FALSE, TRUE),
  hierarchical = c(FALSE, FALSE, TRUE, TRUE)
)

make_formula <- function(controls, hierarchical) {
  rhs <- c("log_rprice_c",
           if (controls) c("log_rpimin_c", "log_rincome_c", "year_c"),
           if (hierarchical) "(1 + log_rprice_c | state)")
  bf(as.formula(paste("log_sales ~", paste(rhs, collapse = " + "))), family = gaussian())
}
formulas <- setNames(Map(make_formula, specs$controls, specs$hierarchical), specs$model)

make_priors <- function(controls, hierarchical, f = 1) {
  p <- c(
    set_prior(sprintf("normal(4.8, %s)", 1 * f), class = "Intercept"),
    set_prior(sprintf("normal(0, %s)", 1 * f), class = "b"),
    set_prior(sprintf("normal(0, %s)", 0.5 * f), class = "sigma")
  )
  if (controls) {
    p <- c(p, set_prior(sprintf("normal(0, %s)", 0.5 * f), class = "b", coef = "year_c"))
  }
  if (hierarchical) {
    p <- c(p,
           set_prior(sprintf("normal(0, %s)", 0.5 * f), class = "sd", group = "state"),
           set_prior("lkj(2)", class = "cor"))
  }
  p
}
priors <- setNames(Map(make_priors, specs$controls, specs$hierarchical), specs$model)

key_pars <- function(controls, hierarchical) {
  c("b_log_rprice_c",
    if (controls) c("b_log_rpimin_c", "b_log_rincome_c", "b_year_c"),
    "sigma",
    if (hierarchical) c("sd_state__Intercept", "sd_state__log_rprice_c"))
}
pars <- setNames(Map(key_pars, specs$controls, specs$hierarchical), specs$model)

fit_prior <- function(m) {
  brm(formulas[[m]], data = cigar, prior = priors[[m]],
      sample_prior = "only", chains = 4, iter = 2000, warmup = 1000,
      seed = SEED, refresh = 0,
      file = file.path("fits", paste0("wf_", tolower(m), "_prior")))
}

fit_model <- function(m, f = 1, tag = "") {
  i <- match(m, specs$model)
  brm(formulas[[m]], data = cigar,
      prior = make_priors(specs$controls[i], specs$hierarchical[i], f),
      chains = 4, iter = if (specs$hierarchical[i]) 4000 else 2000, warmup = 1000,
      seed = SEED, refresh = 0, save_pars = save_pars(all = TRUE),
      file = file.path("fits", paste0("wf_", tolower(m), tag)))
}

get_fits <- function() setNames(lapply(specs$model, fit_model), specs$model)
