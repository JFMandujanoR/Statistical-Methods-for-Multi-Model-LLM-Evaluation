# ============================================================================
# anova_llm_eval_mmlupro_real.R
#
# Reproduces the real-data application (paper Section 5.2 / Appendix B): six
# openly available ~7B models, scored on n = 1,497 real, shared MMLU-Pro
# questions (a category-stratified subsample of the official cached
# predictions released with the MMLU-Pro benchmark).
#
# Data source (real, not simulated):
#   Wang, Y., Ma, X., Zhang, G., et al. "MMLU-Pro: A More Robust and
#   Challenging Multi-Task Language Understanding Benchmark." NeurIPS 2024
#   Track on Datasets and Benchmarks. arXiv:2406.01574.
#   Cached model predictions: github.com/TIGER-AI-Lab/MMLU-Pro/eval_results
#
# The accompanying file `mmlu_pro_real_data.csv` contains one row per
# (question_id, model) pair with columns:
#   question_id, model, correct (0/1), category (one of 14 MMLU-Pro subjects)
# Subsample: fixed seed 2026, category-stratified, n = 1,497 of 12,032
# questions (~12%), all 14 categories retained proportionally.
#
# This version fits BOTH the under-specified mixed model (category random
# intercept only) AND the correctly specified one from the paper's Eq. 3,
# (1|category) + (1|category:question), to reproduce the diagnostic
# discussed in Appendix B: a cluster-only random intercept is crossed with
# the model fixed effect and so contributes little to between-model
# contrasts, understating precision relative to the paired tests. Only the
# nested question-within-category term recovers the item-level pairing.
#
# BUGFIX (this version): the Bonferroni flag on the fully specified mixed
# model's contrasts must be computed from RAW (unadjusted) p-values derived
# from each contrast's z-ratio, NOT from the already Tukey-adjusted p-value
# emmeans returns. Applying a second correction on top of an already-
# corrected p-value is wrong and was too conservative in an earlier version
# of this script (it found only 1 Bonferroni-significant pair instead of 2,
# with a third landing just outside the threshold -- see printed output).
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
# 1-2. Unpaired vs. paired pairwise z-tests (Table 3)
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
cat("\n--- Table 3: Paired pairwise z-tests + Bonferroni/Holm ---\n")
print(paired)
cat("Uncorrected significant:", sum(paired$p < .05), "/", m, "\n")
cat("Bonferroni significant:", sum(paired$bonferroni_sig), "/", m, "\n")
cat("Holm significant:", sum(paired$holm_sig), "/", m, "\n")

# ----------------------------------------------------------------------------
# 3. Omnibus repeated-measures (randomized-block) ANOVA vs. naive one-way,
#    plus the Friedman test as a non-parametric robustness check appropriate
#    for the binary outcome (Table 5 / Section 5.2)
# ----------------------------------------------------------------------------
fit_naive <- aov(correct ~ model, data = df)
cat("\n--- Naive one-way ANOVA (ignores question blocking) ---\n")
print(summary(fit_naive))

fit_block <- aov(correct ~ model + question_id, data = df)
cat("\n--- Repeated-measures / randomized-block ANOVA (model + question) ---\n")
print(anova_tbl <- summary(fit_block))

friedman_result <- friedman.test(wide_mat)
cat("\n--- Friedman test (non-parametric omnibus, Demsar 2006) ---\n")
print(friedman_result)

# ----------------------------------------------------------------------------
# 4a. UNDER-SPECIFIED mixed-effects model: category random intercept only
#     (kept here deliberately as the diagnostic example discussed in
#     Appendix B -- this is what NOT to stop at)
# ----------------------------------------------------------------------------
mm_undersp <- lmer(correct ~ model + (1 | category), data = df)
cat("\n--- [DIAGNOSTIC] Under-specified mixed model: correct ~ model + (1|category) ---\n")
print(anova(mm_undersp))
cat("\n--- Variance components (under-specified) ---\n")
print(VarCorr(mm_undersp))

cat("\n--- Tukey-adjusted contrasts, under-specified model (compare to Table 4) ---\n")
cat("(Note the contrast SEs here are close to the naive, unpaired SE --\n",
    " this is the symptom described in Appendix B: a category-only random\n",
    " intercept is crossed with model and contributes little to between-\n",
    " model contrasts.)\n")
print(emmeans(mm_undersp, pairwise ~ model, adjust = "tukey")$contrasts)

# ----------------------------------------------------------------------------
# 4b. CORRECTLY SPECIFIED mixed-effects model: category random intercept +
#     question nested within category, matching Eq. 3 in the paper. This is
#     the model that should be reported as the headline result (Table 4).
# ----------------------------------------------------------------------------
mm <- lmer(correct ~ model + (1 | category) + (1 | category:question_id), data = df)
cat("\n--- [CORRECT] Fully specified mixed model: correct ~ model + (1|category) + (1|category:question_id) ---\n")
print(anova(mm))
cat("\n--- Variance components (fully specified) ---\n")
print(VarCorr(mm))

cat("\n--- Table 4: Tukey-adjusted pairwise contrasts, fully specified model ---\n")
cat("(Approximation: Tukey's exact FWER guarantee assumes a balanced,\n",
    " non-clustered design -- see Section 4's caveat. For an exact\n",
    " guarantee under this clustered design, use Holm/Bonferroni on the\n",
    " RAW p-values below, not the Tukey-adjusted ones.)\n")
contrasts_full <- emmeans(mm, pairwise ~ model, adjust = "tukey")$contrasts
print(contrasts_full)

contrasts_df <- as.data.frame(contrasts_full)
# BUGFIX: compute Bonferroni significance from the RAW p-value implied by
# each contrast's z-ratio, not from contrasts_df$p.value (which is already
# Tukey-adjusted). Stacking Bonferroni on top of a Tukey adjustment is
# over-conservative and was wrong in an earlier version of this script.
contrasts_df$p_raw <- 2 * (1 - pnorm(abs(contrasts_df$z.ratio)))
contrasts_df$bonferroni_sig <- contrasts_df$p_raw < 0.05 / nrow(contrasts_df)
cat("\n--- Same contrasts with CORRECTED Bonferroni flag (raw p, exact FWER guarantee) ---\n")
print(contrasts_df[, c("contrast", "estimate", "SE", "z.ratio", "p_raw", "bonferroni_sig")])
cat("\nBonferroni significant (correct):", sum(contrasts_df$bonferroni_sig), "/", nrow(contrasts_df), "\n")

cat("\n--- Summary: does the fully specified mixed model recover the paired-test precision? ---\n")
cat(sprintf("Paired-test SE for M-Inst-v0.1 vs M-v0.1: %.4f\n",
            paired$se[paired$model_a == "Mistral-7B-Instruct-v0.1" & paired$model_b == "Mistral-7B-v0.1" |
                      paired$model_b == "Mistral-7B-Instruct-v0.1" & paired$model_a == "Mistral-7B-v0.1"]))
cat("Mixed-model (fully specified) contrast SE: see 'SE' column above -- should be close.\n")
