# ============================================================================
# anova_llm_eval_mmlupro_real.R
#
# Reproduces the real-data application (paper Section 5.2): six openly
# available ~7B models, scored on n = 1,497 real, shared MMLU-Pro questions
# (a category-stratified subsample of the official cached predictions
# released with the MMLU-Pro benchmark).
#
# Data source (real, not simulated):
#   Wang, Y., Ma, X., Zhang, G., et al. "MMLU-Pro: A More Robust and
#   Challenging Multi-Task Language Understanding Benchmark." NeurIPS 2024
#   Track on Datasets and Benchmarks. arXiv:2406.01574.
#   Cached model predictions: github.com/TIGER-AI-Lab/MMLU-Pro/eval_results
#
# The accompanying file `mmlu_pro_real_data.csv` (shipped alongside this
# script) contains one row per (question_id, model) pair with columns:
#   question_id, model, correct (0/1), category (one of 14 MMLU-Pro subjects)
# It is the exact subsample used to produce Table 1 and the numbers in
# Section 5.2 of the paper -- running this script reproduces them exactly
# (no internet access required; the CSV is self-contained).
#
# Run top to bottom in RStudio, with mmlu_pro_real_data.csv in the working
# directory (or edit the path below).
# ============================================================================

need <- c("dplyr", "tidyr", "lme4", "lmerTest", "emmeans")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)

df <- read.csv("mmlu_pro_real_data.csv", stringsAsFactors = FALSE)
df$model <- factor(df$model)
df$category <- factor(df$category)
df$question_id <- factor(df$question_id)

n_q <- length(unique(df$question_id))
K <- length(levels(df$model))
cat(sprintf("K = %d models, n = %d shared questions, %d total observations\n",
            K, n_q, nrow(df)))

# ----------------------------------------------------------------------------
# Table 1: per-model mean + naive SE
# ----------------------------------------------------------------------------
desc <- df %>%
  group_by(model) %>%
  summarise(mean_acc = mean(correct), se_naive = sd(correct) / sqrt(n()), .groups = "drop") %>%
  arrange(desc(mean_acc))
cat("\n--- Table 1: per-model mean accuracy and naive SE ---\n")
print(desc)

wide <- df %>%
  dplyr::select(question_id, model, correct) %>%
  pivot_wider(names_from = model, values_from = correct) %>%
  arrange(question_id)
wide_mat <- as.matrix(wide[, -1])
rownames(wide_mat) <- wide$question_id

# ----------------------------------------------------------------------------
# 1-2. Unpaired vs. paired pairwise z-tests
# ----------------------------------------------------------------------------
models <- colnames(wide_mat)
pairs <- combn(models, 2, simplify = FALSE)

unpaired <- do.call(rbind, lapply(pairs, function(p) {
  a <- wide_mat[, p[1]]; b <- wide_mat[, p[2]]
  se <- sqrt(var(a) / n_q + var(b) / n_q)
  z <- (mean(a) - mean(b)) / se
  data.frame(model_a = p[1], model_b = p[2], diff = mean(a) - mean(b),
             se = se, z = z, p = 2 * (1 - pnorm(abs(z))))
}))
cat("\n--- Unpaired pairwise z-tests ---\n")
print(unpaired[order(unpaired$p), ])
cat("Uncorrected significant (.05):", sum(unpaired$p < .05), "/", nrow(unpaired), "\n")

paired <- do.call(rbind, lapply(pairs, function(p) {
  d <- wide_mat[, p[1]] - wide_mat[, p[2]]
  se <- sd(d) / sqrt(n_q)
  z <- mean(d) / se
  data.frame(model_a = p[1], model_b = p[2], diff = mean(d),
             se = se, z = z, p = 2 * (1 - pnorm(abs(z))))
}))
paired <- paired[order(paired$p), ]
m <- nrow(paired)
paired$bonferroni_sig <- paired$p < 0.05 / m
paired$holm_sig <- sapply(seq_len(m), function(i) paired$p[i] < 0.05 / (m - i + 1))
cat("\n--- Paired pairwise z-tests + Bonferroni/Holm ---\n")
print(paired)
cat("Uncorrected significant:", sum(paired$p < .05), "/", m, "\n")
cat("Bonferroni significant:", sum(paired$bonferroni_sig), "/", m, "\n")
cat("Holm significant:", sum(paired$holm_sig), "/", m, "\n")

# ----------------------------------------------------------------------------
# 3. Omnibus repeated-measures (randomized-block) ANOVA vs. naive one-way
# ----------------------------------------------------------------------------
fit_naive <- aov(correct ~ model, data = df)
cat("\n--- Naive one-way ANOVA (ignores question blocking) ---\n")
print(summary(fit_naive))

fit_block <- aov(correct ~ model + question_id, data = df)
cat("\n--- Repeated-measures / randomized-block ANOVA (model + question) ---\n")
print(anova_tbl <- summary(fit_block))

# ----------------------------------------------------------------------------
# 4. Mixed-effects model with subject category as a random effect
# ----------------------------------------------------------------------------
mm <- lmer(correct ~ model + (1 | category), data = df)
cat("\n--- Mixed-effects model: correct ~ model + (1|category) ---\n")
print(anova(mm))
cat("\n--- Variance components ---\n")
print(VarCorr(mm))

cat("\n--- Estimated-marginal-means pairwise contrasts, Tukey-type adjustment ---\n")
cat("(Approximation: FWER control is not exact here since the design is\n",
    " unbalanced across categories -- see the paper's caveat in Section 4.)\n")
print(emmeans(mm, pairwise ~ model, adjust = "tukey")$contrasts)
