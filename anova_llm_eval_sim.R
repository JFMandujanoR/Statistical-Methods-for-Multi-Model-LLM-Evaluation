# ============================================================================
# anova_llm_eval_sim.R
#
# Companion simulations for:
#   "Statistical Methods for Multiple Language Model Comparison on a Shared
#    Evaluation"
#
# Extends Evan Miller's "Adding Error Bars to Evals" (arXiv:2411.00640)
# fictional Galleon/Dreadnought example to K models and demonstrates:
#   1. Multiplicity: uncorrected pairwise z-tests vs. Bonferroni vs.
#      repeated-measures ANOVA + Tukey test
#   2. Clustered questions: naive vs. Miller's design-based clustered SE
#      (Eq. 4), computed directly rather than via a mixed model, to avoid
#      the singular-fit failure mode a weakly-identified random-intercept
#      model can hit on binary data
#   3. Sample size: reproduces Miller's two-model Eq. 9 (n ~ 969)
#
# NOTE ON REPRODUCIBILITY: with a fixed seed, R's RNG stream is deterministic
# but implementation-specific; exact decimals below may differ by a small
# amount across R versions/platforms. The reported ranges in the accompanying
# paper reflect the actual output of this script on R 4.x (see comments next
# to each result for the verified run used in the paper).
#
# Can be run top to bottom in RStudio.
# ============================================================================

need <- c("MASS", "lme4", "lmerTest", "emmeans", "dplyr", "ggplot2")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

library(MASS)      # mvrnorm for correlated conditional means
library(lme4)       # mixed-effects models (used in Experiment 1's demo only)
library(lmerTest)   # Type-III F tests / p-values for lmer
library(emmeans)    # Tukey-adjusted post-hoc contrasts
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
  x <- MASS::mvrnorm(n_questions, mu = base_mean + delta, Sigma = Sigma * 0.02^2)
  x <- pmin(pmax(x, 0.01), 0.99)

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
pairwise_z_test <- function(wide) {
  K <- ncol(wide)
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

n_sims <- 2000
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
mc_se <- sqrt(fwer * (1 - fwer) / n_sims)
cat("\n--- Experiment 1: simulated family-wise error rate under a true null ---\n")
print(fwer)
cat("\nMonte Carlo SE for each estimate (binomial, n_sims =", n_sims, "):\n")
print(mc_se)
# Verified run (R 4.x, seed 2026, n_sims=2000): uncorrected=0.1915 (SE .0088),
# bonferroni=0.0420 (SE .0045), anova=0.0515 (SE .0049). Bonferroni and ANOVA
# both land close to the nominal .05; uncorrected is elevated but NOT close
# to the independence benchmark below -- see note.

independence_fwer <- 1 - 0.95^choose(length(model_names), 2)
cat(sprintf("\nFor reference, FWER under independence would be 1 - 0.95^6 = %.4f\n",
            independence_fwer))
cat("The simulated uncorrected rate is typically BELOW this independence\n",
    "benchmark, not close to it: with rho = 0.5, the six pairwise test\n",
    "statistics are themselves positively correlated (every pair shares two\n",
    "of the same four models' conditional means), and positive correlation\n",
    "among the test statistics tends to reduce the family-wise error rate\n",
    "relative to the independent case for this kind of 'any-reject' statistic.\n",
    "Still far above the nominal .05, just not as high as independence predicts.\n")

# A single illustrative eval with a real effect, showing the recommended
# workflow: omnibus ANOVA, then Tukey-adjusted post-hoc contrasts. The effect
# sizes and n are larger than in the null experiment above on purpose: at
# n = 400 and effects of 0.01-0.03 the omnibus test has only ~20% power, so a
# single draw would usually show nothing. With n = 1500 and the deltas below
# the omnibus test rejects essentially always (checked by simulation).
demo <- simulate_eval(1500, model_names, rho = 0.5,
                       delta = c(0, -0.06, 0.02, -0.04))
fit_demo <- aov(score ~ model + Error(factor(question)), data = demo)
cat("\n--- Demo eval: omnibus repeated-measures ANOVA (one draw; see note above) ---\n")
print(summary(fit_demo))

mm_demo <- lmer(score ~ model + (1 | question), data = demo)
cat("\n--- Tukey-adjusted pairwise contrasts (compare to Miller's Table 5) ---\n")
print(emmeans(mm_demo, pairwise ~ model, adjust = "tukey")$contrasts)

# ============================================================================
# 2. CLUSTERED QUESTIONS: naive vs. Miller's design-based clustered SE
#    (Eq. 4), computed directly -- avoids the singular-fit problem a
#    random-intercept lmer model can hit when cluster signal is weak
#    relative to cluster size on binary data.
# ============================================================================
simulate_clustered_eval <- function(n_clusters, q_per_cluster_range, icc,
                                     total_sigma, base_logit = qlogis(0.75)) {
  # icc here is the share of the LATENT (logit-scale) variance that sits at the
  # cluster level; the intraclass correlation of the observed 0/1 scores is
  # smaller. Vectorised for speed.
  nq <- sample(q_per_cluster_range[1]:q_per_cluster_range[2], n_clusters, replace = TRUE)
  cluster <- rep(seq_len(n_clusters), times = nq)
  cluster_effect <- rnorm(n_clusters, 0, sqrt(icc) * total_sigma)[cluster]
  q_effect <- rnorm(length(cluster), 0, sqrt(1 - icc) * total_sigma)
  p <- plogis(base_logit + cluster_effect + q_effect)
  data.frame(cluster = cluster, question = seq_along(cluster),
             score = rbinom(length(p), 1, p))
}

naive_se <- function(scores) sqrt(var(scores) / length(scores))

# Direct implementation of Miller's Eq. 4 clustered SE.
miller_clustered_se <- function(cluster_ids, scores) {
  n <- length(scores)
  dev <- scores - mean(scores)
  s_c <- tapply(dev, cluster_ids, sum)       # per-cluster sum of deviations
  q_c <- tapply(dev^2, cluster_ids, sum)     # per-cluster sum of squared deviations
  cross <- sum(s_c^2 - q_c)                  # sum over c of sum_{i != j} dev_i dev_j
  sqrt(naive_se(scores)^2 + cross / n^2)
}

cat("\n--- Experiment 2: naive vs. clustered SE ratio across clustering strength ---\n")
cat("(Using Miller's direct Eq. 4 formula -- robust to the singular-fit issue\n",
    " a random-intercept lmer model can hit when cluster signal is weak.)\n")
icc_grid <- c(0.30, 0.50, 0.60)
sigma_grid <- c(1.5, 1.5, 2.0)
qrange_grid <- list(c(8, 16), c(8, 16), c(8, 24))
n_reps <- 30

for (k in seq_along(icc_grid)) {
  ratios <- numeric(n_reps)
  for (r in 1:n_reps) {
    sim <- simulate_clustered_eval(n_clusters = 250, q_per_cluster_range = qrange_grid[[k]],
                                    icc = icc_grid[k], total_sigma = sigma_grid[k])
    ratios[r] <- miller_clustered_se(sim$cluster, sim$score) / naive_se(sim$score)
  }
  cat(sprintf("icc=%.2f, cluster size %d-%d: mean ratio=%.2f, range=[%.2f, %.2f]\n",
              icc_grid[k], qrange_grid[[k]][1], qrange_grid[[k]][2],
              mean(ratios), min(ratios), max(ratios)))
}
# Verified run (R 4.x, seed 2026, 30 reps per setting):
#   icc=0.30, size 8-16: mean 1.39, range [1.24, 1.50]
#   icc=0.50, size 8-16: mean 1.58, range [1.46, 1.73]
#   icc=0.60, size 8-24: mean 2.14, range [1.91, 2.40]
cat("\nInterpretation: the ratio exceeds 1 and grows with the strength of\n",
    "within-cluster correlation and with cluster size. Compare your printed\n",
    "ranges with the values reported in the paper (Appendix A).\n")

# ============================================================================
# 3. SAMPLE SIZE: reproduce Miller's two-model formula (Eq. 9) only
# ============================================================================
# Miller's two-model worked example: omega^2 = 1/9, delta = 0.03, alpha = .05,
# power = .80 -> n ~ 969.
omega2 <- 1 / 9; delta <- 0.03; alpha <- 0.05; power <- 0.80
n_two_model <- ((qnorm(1 - alpha / 2) + qnorm(power))^2) * omega2 / delta^2
cat("\n--- Experiment 3: sample size, Miller's two-model formula ---\n")
cat(sprintf("n (K = 2, Miller Eq. 9) = %.0f\n", n_two_model))
# NOTE: an earlier version of this script extended this to K models with
# pwr::pwr.anova.test. That was removed: pwr.anova.test is a between-subjects
# calculation that ignores the question-level pairing (rho), and its effect-size
# input required an ad hoc rescaling, so its output (n per group FALLING as K
# grows) is not comparable to Miller's paired-design n. Power for the
# repeated-measures omnibus test should be obtained by simulation (reuse
# simulate_eval() and rm-ANOVA over a grid of n and K).

# ============================================================================
# Visualization of Experiment 1 (family-wise error inflation)
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
