# BDA 2026 Project — Instructions for Claude Code

This repository is our project for **Bayesian Data Analysis (CS-E5710), Aalto 2026**.
The goal is a full Bayesian workflow analysis of one dataset, written up as a report
of **at most 20 pages** (appendix excluded) and turned into a 15-minute presentation.

- **Report deadline: Sun 29 Nov 2026, 23:59** (FeedbackFruits, peer-graded with the course rubric)
- Presentation: 7–11 Dec 2026 (on campus, 15 min + 5–10 min questions)
- Course page: https://avehtari.github.io/BDA_course_Aalto/project.html
- Rubric: https://avehtari.github.io/BDA_course_Aalto/project_rubric.html

> **Fill in before starting**
> - Dataset: Cigarette consumption panel for 46 US states, 1963–1992 (Baltagi & Levin 1992;
>   also in the R package `plm` as `Cigar`). File: `Cigar.csv` (to be moved to `data/`).
>   1,380 rows, no missing values. Earlier analyses to cite: Baltagi & Levin (1986, 1992),
>   Baltagi, *Econometric Analysis of Panel Data*.
> - Research question: How strongly is cigarette consumption per person associated with real
>   cigarette prices in US states (1963–1992), and does a hierarchical model with state-specific
>   price elasticities predict better than a pooled linear regression?
> - Outcome variable and type: `sales` (packs per person per year), continuous and positive;
>   modelled as `log(sales)`. Main covariate: `log(price / cpi)`, the log real price. Other
>   covariates: `log(ndi / cpi)` (real income), `log(pimin / cpi)` (lowest real price in
>   neighbouring states), and a centred year trend.
> - Grouping variable(s) for the hierarchical model: `state` (46 states × 30 years each).
>   Planned models: M1 pooled; M2 with a varying intercept and price slope by state, `(1 + log_rprice | state)`.
> - Group members: `<names>`

---

## 1. Ground rules for Claude Code

1. **We are the authors.** Course policy says AI is allowed but most of the work must be
   done by the students, and AI use must be reported. So:
   - Explain every modelling choice in comments and in plain language to us before
     moving on, so we can defend it in the oral exam.
   - Keep a running log in `AI_USAGE.md`: date, what you were asked, what you produced,
     what we changed. We'll use it to write the AI-use statement in the report.
   - Don't write the final report prose on your own. Give us drafts or bullet points, and we'll write it.
2. **Never invent results.** Every number, figure and claim in the report must come from code
   that has actually been run in this repo. If something fails, say so.
3. **Don't overclaim causality** unless the design supports it.
4. **Ask before** changing the dataset, the research question or the set of models.

## 2. Tooling

- Language: **R** (≥ 4.x), using **brms** with the **cmdstanr** backend (`backend = "cmdstanr"`).
- Use: `brms`, `cmdstanr`, `posterior`, `bayesplot`, `loo`, `priorsense`, `ggplot2`, `dplyr`.
- Report: **Quarto** (`report/report.qmd`) rendered to PDF (fallback: HTML → print to PDF).
- Seeds: set `seed =` in every `brm()` call and `set.seed()` for any other randomness.
- Cache fitted models with `file = "fits/<name>"` in `brm()` so re-rendering is fast.

## 3. Repository layout

```
data/            raw data (read-only) + data/processed/
R/01_data.R      load, clean, explore, subsample if > ~1e6 rows
R/02_models.R    model formulas, priors, fits
R/03_checks.R    convergence diagnostics + posterior predictive checks
R/04_sensitivity.R  prior sensitivity (priorsense + alternative priors)
R/05_compare.R   PSIS-LOO comparison
fits/            cached brms fits (git-ignored)
figures/         saved figures for report and slides
report/report.qmd
slides/          presentation
AI_USAGE.md
```

## 4. Hard requirements (from the course page — all must be satisfied)

### Data
- Not allowed: titanic, mtcars, iris, Palmer penguins, the baseball batting data,
  any dataset used in course demos (e.g. bodyfat, diabetes), or datasets used in other CS courses.
- Ideally 100 < n < 1,000,000 (subsample if larger) and has group structure for a hierarchical model.
- If the data has been analysed online before: cite it and state clearly how our analysis differs.

### Models
- **At least two models**, e.g. pooled vs. hierarchical, linear vs. non-linear, or two
  observation models (Poisson vs. negative binomial, etc.).
- At least some models must include **covariates**.
- Hierarchical = a parameter's prior depends on other estimated parameters
  (e.g. `b_j ~ normal(mu, sigma)`, `mu ~ normal(0,1)`, `sigma ~ exponential(1)`).
- **Every parameter needs an explicit proper prior.** brms defaults are flat (improper)
  for `b` coefficients, so always set them with `set_prior()`. Check with `default_prior()`/
  `prior_summary()` and **report any default proper priors we left unchanged.**
- No hard bounds / `uniform(a, b)` unless values outside are truly impossible.
- For group-level SDs prefer half-normal or half-t priors (course notes recommend these over half-Cauchy).
- Do **prior predictive checks** (`sample_prior = "only"`) to show priors give plausible data.

### Inference
- Report the exact `brm()` call(s): iterations, warmup, chains, `adapt_delta` etc., **with a
  reason for each choice** (default: 4 chains, 2000 iter, 1000 warmup unless there's cause to change).

### Convergence diagnostics — for ALL models, with interpretation
- R-hat (target < 1.01), bulk and tail ESS, number of divergences, max treedepth hits.
- Summarise in text when everything is fine. Full tables go in the appendix.
- If problems: describe the fix (reparameterise, adapt_delta, better priors) and refit.

### Posterior predictive checks — for ALL models, with interpretation
Choose by outcome type:
- Continuous → `pp_check(type = "dens_overlay")` plus stat checks
- Discrete counts → rootogram (`type = "rootogram"`), bars, or stat checks (KDE only if many distinct values)
- Binary / ordinal / categorical → **calibration plots** (e.g. reliability diagrams, `bayesplot` PAV-adjusted calibration), `bars`
- Also use LOO-PIT (`pp_check(type = "loo_pit_overlay")` or `ppc_loo_pit_qq`) where appropriate.
- Say what misfit you see and whether/how the model was improved.

### Prior sensitivity — for ALL models, with interpretation
- `priorsense::powerscale_sensitivity()` and `powerscale_plot_dens()`.
- Plus at least one refit with alternative priors (wider/narrower), comparing the key estimates.
- State whether the estimates of interest changed and what that means.

### Model comparison
- `loo()` for each model, check Pareto k̂ (report any k̂ > 0.7 and handle via `moment_match`/`reloo`).
- `loo_compare()` with **elpd_diff and its SE**, interpreted (difference small if |diff| < ~4 or < 2·SE).
- Add a practical metric if relevant (e.g. classification accuracy, RMSE on the outcome scale).
- Justify the preferred model by predictive performance, interpretability, simplicity and the goal.

### Numbers
- Report a sensible number of digits (follow the course "Digits" notebook: don't print more
  digits than the Monte Carlo error supports). Use `posterior::summarise_draws()` + rounding.

## 5. Report structure (`report/report.qmd`)

1. **Introduction:** motivation, problem, main modelling idea, one illustrative figure
2. **Data:** source, description, preprocessing, earlier analyses and how ours differs
3. **Models:** mathematical notation for each model (observation model + priors)
4. **Priors:** justification (weakly vs. informative, ideally backed by domain sources), prior predictive check
5. **Model code:** brms formulas and Stan code if custom (long code → appendix, referenced in text)
6. **Inference:** sampling command and settings with reasons
7. **Convergence diagnostics:** for all models, interpreted
8. **Posterior predictive checks:** for all models, interpreted
9. **Prior sensitivity analysis:** for all models, interpreted
10. **Model comparison:** LOO with uncertainty, choice of model, justified
11. **Results:** posterior estimates for the research question, on an interpretable scale
12. **Discussion:** limitations, problems, possible improvements
13. **Conclusion:** main finding answering the research question
14. **Self-reflection:** what the group learned
15. **AI-use statement:** whether and how AI was used (from `AI_USAGE.md`)
16. **References** · **Appendix** (full code, full diagnostic tables)

Style: tell the story of what we did and why, not a list of outputs. Keep code output short in the
main text (`echo: false`/`output: false` where appropriate). Number and caption all figures/tables
(`#fig-`, `#tbl-`, `#eq-` labels) and reference them in the text. Fonts in figures must be readable.

## 6. Pre-submission checklist (mirrors the peer-review rubric)

- [ ] Introduction describes both the data and the analysis, with motivation
- [ ] ≥ 2 models clearly described, reasonable for the problem
- [ ] All priors listed, proper, justified; unchanged defaults reported; no flat priors on `b`
- [ ] Model code included and readable; sampling command visible
- [ ] R-hat, ESS, divergences/treedepth reported **for all models** and **discussed**
- [ ] PPCs **for all models**, appropriate to outcome type, **discussed**
- [ ] Prior sensitivity **for all models**, **discussed**
- [ ] LOO comparison with SE, Pareto k̂ checked, **discussed**
- [ ] Problems and improvements discussed
- [ ] Clear conclusion tied to the research question
- [ ] Self-reflection section present
- [ ] Statistical terms used correctly (credible interval, not confidence interval; "draws", ESS ≠ "effective samples")
- [ ] Reasonable digits; reasonable amount of code output
- [ ] Main body ≤ 20 pages (excl. title page, TOC, references, appendix)
- [ ] AI-use statement included
- [ ] Report renders to PDF without errors, and author names are removed if anonymity is required

## 7. Presentation (after the report)

- 15 min hard limit. Everyone speaks. First slide: title + all names.
- Explain and justify the observation model and priors; show that every workflow step was done.
- Last slide = take-home messages + contact/link (**not** just "Thank you"/"Questions?").
- Prepare for oral-exam questions: why these priors, what R-hat/ESS/divergences mean, how PSIS-LOO works,
  what the PPCs show, how to interpret elpd_diff.
