# Statistical-Methods-for-Multi-Model-LLM-Evaluation

Statistical methods for comparing more than two language models on a shared evaluation, extending Miller (2024)'s paired and clustered standard-error framework to *K > 2* models via repeated-measures ANOVA and mixed-effects models.

## Files

- `anova_llm_eval_sim.R` — simulation study (family-wise error rate, clustered SE, power analysis)
- `anova_llm_eval_mmlupro_real.R` — real-data application, run on `mmlu_pro_real_data.csv`
- `mmlu_pro_real_data.csv` — real, item-level correctness data for 6 models on 1,497 shared MMLU-Pro questions
- `fwerComparison.png` — figure comparing family-wise error rates across correction methods

## Data source

`mmlu_pro_real_data.csv` is a subsample built from the official cached model predictions released with:

> Yubo Wang et al. "MMLU-Pro: A More Robust and Challenging Multi-Task Language Understanding Benchmark." NeurIPS 2024. arXiv:2406.01574.
> Predictions from [TIGER-AI-Lab/MMLU-Pro](https://github.com/TIGER-AI-Lab/MMLU-Pro) (`eval_results/`), Apache License 2.0.

## Citation
Mandujano Reyes, Juan Francisco. "Statistical Methods for Multiple Language Model Comparison on a Shared Evaluation." arXiv preprint arXiv:2608.22659 (2026).
