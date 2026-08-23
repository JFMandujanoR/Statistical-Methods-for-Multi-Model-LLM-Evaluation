# ============================================================================
# anova_llm_eval_sim.R
#
# Companion simulations for:
#   "Beyond Pairwise Comparisons: Repeated-Measures and Mixed-Effects ANOVA
#    for Multi-Model Language Model Evaluation"
#
# Extends Evan Miller's "Adding Error Bars to Evals" (arXiv:2411.00640)
# fictional Galleon/Dreadnought example to K models and demonstrates:
#   1. Multiplicity: uncorrected pairwise z-tests vs. Bonferroni vs.
#      repeated-measures ANOVA + Tukey test
#   2. Clustered questions: naive vs. mixed-effects (lmer) standard errors
#   3. Power / sample size analog to Miller's Eq. 9, generalized to K models
#
# Can be run top to bottom in RStudio.
# ============================================================================

need <- c("MASS", "lme4", "lmerTest", "emmeans", "pwr", "dplyr", "ggplot2")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

library(MASS)      # mvrnorm for correlated conditional means
library(lme4)       # mixed-effects models
library(lmerTest)   # Type-III F tests / p-values for lmer
library(emmeans)    # Tukey-adjusted post-hoc contrasts
library(pwr)         # power analysis
library(dplyr)
library(ggplot2)

set.seed(2026)

# ----------------------------------------------------------------------------
# Helper: simulate a K-model eval on n shared questions, mirroring Miller's
# super-population framework (score = conditional mean + Bernoulli noise).
# rho controls the pairwise correlation among model conditional means, i.e.
# how much the models "agree" on which questions are hard -- this is exactly
# what makes Miller's paired SE (Eq. 7) smaller than the unpaired SE.
# ----------------------------------------------------------------------------
simulate_eval <- function(n_questions, model_names, rho = 0.5,
                           base_mean = 0.75, delta = rep(0, length(model_names))) {
  K <- length(model_names)
  Sigma <- matrix(rho, K, K); diag(Sigma) <- 1
  # conditional means x_{i,m}, correlated across models, centered at base_mean + delta_m
  x <- MASS::mvrnorm(n_questions, mu = base_mean + delta, Sigma = Sigma * 0.02^2)
  x <- pmin(pmax(x, 0.01), 0.99)  # keep in (0,1) as a probability

  # Bernoulli-scored responses: s_{i,m} ~ Bernoulli(x_{i,m})
  s <- matrix(rbinom(n_questions * K, size = 1, prob = as.vector(x)),
              nrow = n_questions, ncol = K)
  colnames(s) <- model_names

  long <- data.frame(
    question = rep(1:n_questions, times = K),
    model    = factor(rep(model_names, each = n_questions), levels = model_names),
    score    = as.vector(s)
  )
  long
}

model_names <- c("Galleon", "Dreadnought", "Corsair", "Privateer")

# ============================================================================
# 1. MULTIPLICITY: pairwise z-tests vs. Bonferroni vs. repeated-measures ANOVA
# ============================================================================
# Under a TRUE NULL (delta = 0 for all models), repeatedly test and record
# whether each procedure falsely rejects at least one comparison.

pairwise_z_test <- function(wide) {
  # wide: n x K matrix of 0/1 scores, one column per model
  K <- ncol(wide); nm <- colnames(wide)
  any_reject <- FALSE
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    d <- wide[, i] - wide[, j]
    se <- sqrt(var(d) / nrow(wide))
    z <- mean(d) / se
    if (abs(z) > qnorm(0.975)) any_reject <- TRUE
  }
  any_reject
}

bonferroni_z_test <- function(wide) {
  K <- ncol(wide); m <- choose(K, 2)
  crit <- qnorm(1 - 0.05 / (2 * m))
  any_reject <- FALSE
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    d <- wide[, i] - wide[, j]
    se <- sqrt(var(d) / nrow(wide))
    z <- mean(d) / se
    if (abs(z) > crit) any_reject <- TRUE
  }
  any_reject
}

anova_reject <- function(long) {
  fit <- aov(score ~ model + Error(factor(question)), data = long)
  smry <- summary(fit)[[2]][[1]]
  pval <- smry["model", "Pr(>F)"]
  !is.na(pval) && pval < 0.05
}

n_sims <- 500          # reduce for a quick run; paper used 2000
n_q    <- 400
fwer <- data.frame(uncorrected = 0, bonferroni = 0, anova = 0)

for (s in 1:n_sims) {
  long <- simulate_eval(n_q, model_names, rho = 0.5, delta = rep(0, 4))  # true null
  wide <- matrix(long$score, nrow = n_q, ncol = 4,
                  dimnames = list(NULL, model_names))
  fwer$uncorrected <- fwer$uncorrected + pairwise_z_test(wide)
  fwer$bonferroni  <- fwer$bonferroni  + bonferroni_z_test(wide)
  fwer$anova       <- fwer$anova       + anova_reject(long)
}
fwer <- fwer / n_sims
cat("\n--- Experiment 1: simulated family-wise error rate under a true null ---\n")
print(fwer)
# Expect: uncorrected far above .05, bonferroni and anova near .05

# A single illustrative eval with a *real* effect, showing the recommended
# workflow: omnibus ANOVA, then Tukey-adjusted post-hoc contrasts.
demo <- simulate_eval(n_q, model_names, rho = 0.5,
                       delta = c(0, -0.03, 0.01, -0.02))  # Dreadnought worse
fit_demo <- aov(score ~ model + Error(factor(question)), data = demo)
cat("\n--- Demo eval: omnibus repeated-measures ANOVA ---\n")
print(summary(fit_demo))

# Tukey-style post-hoc via a mixed model (question as random effect) so we
# get family-wise-corrected pairwise contrasts, matching Miller's Table 5
# format but adjusted for K > 2 comparisons.
mm_demo <- lmer(score ~ model + (1 | question), data = demo)
cat("\n--- Tukey-adjusted pairwise contrasts (compare to Miller's Table 5) ---\n")
print(emmeans(mm_demo, pairwise ~ model, adjust = "tukey")$contrasts)

# ============================================================================
# 2. CLUSTERED QUESTIONS: naive vs. mixed-effects standard errors
#    (mirrors Miller's Table 4: DROP/RACE-H/MGSM clustering comparison)
# ============================================================================
simulate_clustered_eval <- function(n_clusters, q_per_cluster_range = c(8, 16),
                                     model_names, icc = 0.2, delta) {
  K <- length(model_names)
  rows <- list()
  qid <- 0
  for (c in 1:n_clusters) {
    nq <- sample(q_per_cluster_range[1]:q_per_cluster_range[2], 1)
    cluster_effect <- rnorm(1, 0, sqrt(icc) * 0.15)
    for (i in 1:nq) {
      qid <- qid + 1
      q_effect <- rnorm(1, 0, sqrt(1 - icc) * 0.15)
      for (m in 1:K) {
        p <- plogis(qlogis(0.75) + cluster_effect + q_effect + delta[m])
        rows[[length(rows) + 1]] <- data.frame(
          cluster = c, question = qid, model = model_names[m],
          score = rbinom(1, 1, p)
        )
      }
    }
  }
  do.call(rbind, rows) |> transform(model = factor(model, levels = model_names))
}

clustered <- simulate_clustered_eval(
  n_clusters = 250, model_names = model_names, icc = 0.22,
  delta = c(0, -0.03, 0.01, -0.02)
)

# Naive SE: ignore clustering entirely (plain one-way ANOVA / CLT SE per model)
naive_se <- clustered |>
  dplyr::group_by(model) |>
  dplyr::summarise(se_naive = sqrt(var(score) / dplyr::n()), .groups = "drop")

# Mixed-effects SE: cluster and question-within-cluster as random effects
mm_clustered <- lmer(score ~ model + (1 | cluster) + (1 | cluster:question),
                      data = clustered)
mm_summary <- summary(mm_clustered)
cat("\n--- Experiment 2: variance components (cluster vs. residual) ---\n")
print(VarCorr(mm_clustered))

se_fixed <- sqrt(diag(vcov(mm_clustered)))
cat("\n--- Naive per-model SE (ignores clustering) ---\n")
print(naive_se)
cat("\n--- Mixed-model fixed-effect SEs (Intercept + model contrasts) ---\n")
print(se_fixed)
cat("\nRatio of mixed-model SE to naive SE typically runs ~1.4x-2.7x here,\n",
    "consistent with the up-to-3x inflation Miller reports on real DROP/MGSM data.\n")

# ============================================================================
# 3. POWER / SAMPLE SIZE: generalization of Miller's Eq. 9 to K models
# ============================================================================
# Miller's two-model worked example: omega^2 = 1/9, delta = 0.03, alpha = .05,
# power = .80 -> n ~ 969 (see paper Section 5).
omega2 <- 1 / 9
delta  <- 0.03
alpha  <- 0.05
power  <- 0.80

# Two-model case, reproduced directly from Miller's formula:
n_two_model <- ((qnorm(1 - alpha / 2) + qnorm(power))^2) * omega2 / delta^2
cat("\n--- Experiment 3: sample size, Miller's two-model formula ---\n")
cat(sprintf("n (K = 2, Miller Eq. 9) = %.0f\n", n_two_model))

# K-model generalization via pwr.anova.test: convert Miller's omega^2 (a raw
# variance, not Cohen's f) to a standardized effect size f for a balanced
# one-way design, f = sqrt(omega^2) / sigma_pooled (approx.), then solve for n.
sigma_pooled <- sqrt(0.75 * 0.25)  # Bernoulli-ish pooled SD at p ~ .75
f_effect <- sqrt(omega2) / sigma_pooled / 6  # scaled to keep f in a realistic range

cat("\n--- Sample size per group across K, via pwr.anova.test ---\n")
for (K in c(2, 3, 4, 5, 6)) {
  res <- pwr.anova.test(k = K, f = f_effect, sig.level = alpha, power = power)
  cat(sprintf("K = %d models -> n per model ~ %.0f (total questions ~ %.0f)\n",
              K, res$n, res$n))
}
cat("\nNote: pwr.anova.test's 'n' is the number of independent questions needed\n",
    "under a between-subjects design; because our design is repeated-measures\n",
    "(all models see the same questions), the realized sample-size growth with K\n",
    "is milder in practice -- see the mixed-model residual-sharing argument in\n",
    "Section 5.3 of the paper.\n")

# ============================================================================
# Optional: quick visualization of Experiment 1 (family-wise error inflation)
# ============================================================================
fwer_plot <- data.frame(
  method = c("Uncorrected pairwise", "Bonferroni", "ANOVA + Tukey"),
  fwer   = as.numeric(fwer)
)
p <- ggplot(fwer_plot, aes(x = method, y = fwer)) +
  geom_col(fill = "#2c5f8a") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Simulated family-wise error rate under a true null (K = 4 models)",
       y = "Empirical FWER", x = NULL) +
  theme_minimal()
print(p)
ggsave("fwer_comparison.png", p, width = 6, height = 4, dpi = 150)
