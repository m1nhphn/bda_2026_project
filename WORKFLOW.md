# Workflow — Cigarette prices and consumption (BDA 2026 project)

How we go from `Cigar.csv` to the report. Each stage lists what we do, why, and the
script that does it.

**Research question.** How strongly are cigarette sales per person associated with the
real cigarette price in US states (1963–1992)? Does a hierarchical model with
state-specific price effects predict better than a pooled model?


---

## Pipeline at a glance

```
1. Data and preprocessing        Cigar.csv -> real prices, logs, centring
          |
2. Choose predictors             price + income + neighbour price + year
          |                      (collinearity checked, NOT used for selection)
          |
3. Models (2 x 2 design)
                     PRICE ONLY        PRICE + CONTROLS
     POOLED          M1a               M1b
     HIERARCHICAL    M2a               M2b
          |
4. Priors + prior predictive check
          |
5. Fit with brms                 one joint posterior per model
          |
6. Check every model             convergence, PPC, prior sensitivity
          |
7. Compare with LOO              elpd_diff +/- SE, Pareto k
          |
8. Results and conclusion        price elasticity, limitations
```

---

## How to run

One Quarto file per stage, all in the `analysis/` folder together with `Cigar.csv`. Each one loads `setup.R`
(packages, data, model specs, priors), so each file runs on its own.

**RStudio:** open a stage file in `analysis/` (e.g. `analysis/01_data.qmd`) and click **Run All** or **Render**.
Run them in order the first time: `04_fit.qmd` fits the models and caches them in
`fits/`, and the later stages load them from there.

**Terminal** (needs Quarto installed):

```bash
cd "/Users/hieungo/Desktop/report bda /analysis"
quarto render 01_data.qmd
```

**The folder name ends with a space.** If you see `cannot open file 'setup.R'`,
R is not running in the `analysis/` folder.

| File | Stage | First run |
|---|---|---|
| `01_data.qmd` | 1 Data and figures | seconds |
| `02_predictors.qmd` | 2 Collinearity | seconds |
| `03_models_priors.qmd` | 3–4 Models, priors, prior predictive | a few minutes (4 prior fits) |
| `04_fit.qmd` | 5 Fit M1a, M1b, M2a, M2b | a few minutes |
| `05_checks.qmd` | 6 Convergence, PPC, residuals | ~1 minute |
| `06_sensitivity.qmd` | 6 Prior sensitivity | ~10 minutes (8 refits) |
| `07_compare.qmd` | 7 LOO comparison | ~1 minute |
| `08_results.qmd` | 8 Elasticities | seconds |

After the first run, fits load from `fits/` and every file takes seconds.

---

## Stage 1 — Data and preprocessing ✅ done

**Data.** `Cigar.csv` has 46 US states × 30 years (1963–1992), so 1,380 rows with no
missing values. It's the Baltagi & Levin panel, also in the R package `plm` as `Cigar`.
Earlier analyses must be cited: Baltagi & Levin (1986, 1992) and Baltagi,
*Econometric Analysis of Panel Data*.

| Raw variable | Meaning |
|---|---|
| `state`, `year` | state code, year 63–92 |
| `price` | average price per pack (cents, nominal) |
| `pimin` | minimum price in neighbouring states (cents, nominal) |
| `ndi` | disposable income per person (nominal) |
| `cpi` | consumer price index (1983 = 100) |
| `sales` | packs per person (aged 14+) per year |
| `pop`, `pop16` | population, population aged 16+ |

**Preprocessing:**

| New variable | Formula | Why |
|---|---|---|
| `rprice`, `rpimin`, `rincome` | value / `cpi` × 100 | remove inflation (1983 prices) |
| `log_sales`, `log_rprice`, … | `log()` | slope = elasticity (% change per % change) |
| `log_rprice_c`, … | minus the mean | the intercept means an average state-year |
| `year_c` | (year − 1977.5) / 10 | year trend per decade, centred |

**Visualisation:** `01_data.qmd`

| Figure | Report section | Shows |
|---|---|---|
| `fig_intro_price_sales` | Introduction | pooled line vs. one line per state |
| `fig_time_trends` | Data | real price and sales over time |
| `fig_outcome_dist` | Data | sales vs. log sales distribution |
| `fig_covariates` | Data | log sales vs. each control |
| `fig_state_slopes` | Data | separate price slope per state (exploratory only) |

---

## Stage 2 — Choose predictors

We choose predictors for subject-matter reasons, **before** fitting.

| Predictor | Role | Reason |
|---|---|---|
| `log_rprice_c` | **main** | the research question; its coefficient is the price elasticity |
| `log_rpimin_c` | control | cross-border buying when neighbours are cheaper |
| `log_rincome_c` | control | income affects demand |
| `year_c` | control | smoking changed over time for non-price reasons |

**Not used:** `cpi` (already used for deflating), `pop` (sales are already per person),
`pop16` share (overlaps with year and income, r ≈ 0.75), lagged sales (future work).

**Collinearity check (only a check, not used for selection):**
- log real price ↔ log neighbour price: r = 0.81, so expect uncertain, trading-off
  coefficients
- income ↔ year: r = 0.66

---

## Stage 3 — Models (2 × 2 design)

Each comparison changes **one thing at a time**.

| | Price only | Price + controls |
|---|---|---|
| **Pooled** | **M1a**: `log_sales ~ log_rprice_c` | **M1b**: `log_sales ~ log_rprice_c + log_rpimin_c + log_rincome_c + year_c` |
| **Hierarchical** | **M2a**: M1a + `(1 + log_rprice_c \| state)` | **M2b**: M1b + `(1 + log_rprice_c \| state)` |

- **Rows (M1 vs. M2):** does letting states differ help? This is the main question.
- **Columns (a vs. b):** do the controls change the price effect? This is the robustness check.
- **If time is short:** M1b vs. M2b is the core comparison.

Likelihood: Gaussian on `log_sales`. If the posterior predictive checks show heavy
tails, try Student-t.

---

## Stage 4 — Priors + prior predictive check

All priors must be explicit and proper (the brms default for `b` is flat).

| Parameter | Proposed prior | Status |
|---|---|---|
| Intercept | normal(4.8, 1) | to justify (log of ~120 packs) |
| Slopes `b` | normal(0, 1) | to justify (published elasticities ≈ −0.3 to −0.5) |
| `year_c` | normal(0, 0.5) | to justify |
| `sigma` | half-normal(0, 0.5) | to justify |
| State SDs (M2) | half-normal(0, 0.5) | to justify |
| Correlation (M2) | lkj(2) | to justify |

Prior predictive check: `brm(..., sample_prior = "only")` + `pp_check()` for each model.

---

## Stage 5 — Fit with brms

- brms + cmdstanr, 4 chains, cache fits in `fits/`.
- Default 2000 iterations / 1000 warmup. Change it only with a reason, and report it.
- A hierarchical model gives **one joint posterior**: population means, between-state
  SDs, and every state's intercept and slope, estimated together (partial pooling).
  States are **not** fitted separately and then combined.

---

## Stage 6 — Check every model

1. **Convergence:** R-hat < 1.01, bulk/tail ESS > 400, 0 divergences, no max tree depth hits
2. **Posterior predictive checks:** density overlay, stats (sd, min/max), LOO-PIT
3. **Prior sensitivity:** `priorsense::powerscale_sensitivity()` + refit with wider and
   narrower priors, then compare the key estimates
4. **Residuals over time:** the years within a state are correlated, which is a
   limitation to discuss

---

## Stage 7 — Compare with LOO

- `loo()` for all models. Check Pareto k̂; if k̂ > 0.7, use `moment_match = TRUE`
  (needs `save_pars(all = TRUE)` when fitting).
- `loo_compare()`: report elpd_diff **and** its SE. A difference is small if |diff| < 4 or
  < 2·SE.
- Optional: leave-states-out (grouped K-fold by state) to test predicting a new state.

---

## Stage 8 — Results and conclusion

- Price elasticity (posterior mean + 90% credible interval). Example wording:
  "a 10% higher real price is associated with ~X% fewer packs per person".
- State-specific elasticities from the hierarchical model.
- Limitations: association, not causation; correlated years within a state;
  cross-border buying; omitted habit dynamics.

---

## Files

```
README.md              project rules and course requirements
WORKFLOW.md            this file
AI_USAGE.md            log of AI use (for the AI-use statement)
_old/                  earlier scripts, figures and report draft (not used)
analysis/
  Cigar.csv            raw data (read-only)
  setup.R              packages, data preprocessing, model specs, priors, fit helpers
  01_data.qmd … 08_results.qmd   one file per stage (see "How to run")
  fits/                cached brms fits (wf_m1a.rds, wf_m2b_prior.rds, …)
```

---

## Progress

- [x] Stage 1: data, preprocessing, visualisation
- [ ] Stage 2: confirm the predictor choice as a group
- [ ] Stage 3: decide 2 × 2 vs. M1b/M2b only
- [ ] Stage 4: justify the priors, prior predictive check
- [ ] Stage 5: fit all models
- [ ] Stage 6: diagnostics, PPC, sensitivity for **all** models
- [ ] Stage 7: LOO comparison
- [ ] Stage 8: results, discussion, conclusion
- [ ] Report: self-reflection, AI-use statement, references, appendix
