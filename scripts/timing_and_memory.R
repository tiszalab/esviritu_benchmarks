# timing_and_memory.R -------------------------------------------------------
# Analyze and plot runtime + peak memory usage across tools, using
# data/time_and_mem/timing.mock_and_neg.1.tsv (per-step /usr/bin/time logs)
# and data/time_and_mem/sr_esviritu_readtotals.tsv (input read depth per
# sample, for scaling plots).

library(tidyverse)
library(here)

paths <- list(
  timing  = here::here("data/time_and_mem/timing.mock_and_neg.1.tsv"),
  reads   = here::here("data/time_and_mem/sr_esviritu_readtotals.tsv"),
  out_fig = here::here("charts"),
  out_tab = here::here("stats")
)

# ---- Tool naming & aesthetics (matches scripts/benchmark_stats.Rmd) -------
tool_levels  <- c("esviritu", "centrifuger", "kraken2", "metabuli", "sylph")
tool_labels  <- c(esviritu = "EsViritu", centrifuger = "Centrifuger",
                  kraken2 = "Kraken2", metabuli = "Metabuli", sylph = "Sylph")
tool_palette <- c(esviritu = "#90D4CC", centrifuger = "#FBA72A",
                  kraken2 = "#D3D4D8", metabuli = "#CB7A5C", sylph = "#5785C1")

theme_bench <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

save_fig <- function(plot, name, w = 8, h = 5) {
  ggsave(file.path(paths$out_fig, name), plot, width = w, height = h)
  message("wrote ", file.path(paths$out_fig, name))
}

save_tab <- function(df, name) {
  readr::write_csv(df, file.path(paths$out_tab, name))
  message("wrote ", file.path(paths$out_tab, name))
}

# ---- Load -------------------------------------------------------------
raw <- read_tsv(paths$timing, show_col_types = FALSE)
read_totals <- read_tsv(paths$reads, show_col_types = FALSE)

# Pull a "key=value" out of the comma-separated `params` column.
extract_kv <- function(x, key) str_match(x, paste0(key, "=([^,]+)"))[, 2]

raw <- raw |>
  mutate(
    mode = extract_kv(params, "mode"),
    mhl  = extract_kv(params, "mhl"),
    ms   = extract_kv(params, "ms"),
    conf = extract_kv(params, "conf"),
    mhg  = extract_kv(params, "mhg"),
    prec = extract_kv(params, "prec"),
    ptax = extract_kv(params, "ptax"),
    mnk  = extract_kv(params, "mnk")
  )

# ---- Collapse multi-step pipelines into one total runtime + peak memory --
# Centrifuger: classify (depends on mhl only) -> quant (depends on mhl, ms).
# Wall time is additive; peak memory is the max across the two steps.
cent_classify <- raw |>
  filter(tool == "centrifuger", step == "classify") |>
  select(benchmark, sample, mhl, classify_sec = real_sec, classify_mem = max_mem_kb)

cent_quant <- raw |>
  filter(tool == "centrifuger", step == "quant") |>
  select(benchmark, sample, mhl, ms, quant_sec = real_sec, quant_mem = max_mem_kb)

cent_total <- cent_quant |>
  left_join(cent_classify, by = c("benchmark", "sample", "mhl")) |>
  mutate(tool = "centrifuger",
         param_combo = paste0(mhl, "_ms", ms),
         total_sec   = classify_sec + quant_sec,
         peak_mem_kb = pmax(classify_mem, quant_mem)) |>
  select(benchmark, tool, sample, param_combo, total_sec, peak_mem_kb)

# Sylph: sketch (no params) -> profile (depends on mnk). One sketch per
# sample is reused across all mnk settings.
syl_sketch <- raw |>
  filter(tool == "sylph", step == "sketch") |>
  select(benchmark, sample, sketch_sec = real_sec, sketch_mem = max_mem_kb)

syl_profile <- raw |>
  filter(tool == "sylph", step == "profile") |>
  select(benchmark, sample, mnk, profile_sec = real_sec, profile_mem = max_mem_kb)

syl_total <- syl_profile |>
  left_join(syl_sketch, by = c("benchmark", "sample")) |>
  mutate(tool = "sylph",
         param_combo = paste0("mnk", mnk),
         total_sec   = sketch_sec + profile_sec,
         peak_mem_kb = pmax(sketch_mem, profile_mem)) |>
  select(benchmark, tool, sample, param_combo, total_sec, peak_mem_kb)

# EsViritu, Kraken2, Metabuli each run as a single step per parameter set.
single_total <- raw |>
  filter(tool %in% c("esviritu", "kraken2", "metabuli")) |>
  mutate(param_combo = case_when(
           tool == "esviritu" ~ mode,
           tool == "kraken2"  ~ paste0("conf", conf, "_mhg", mhg),
           tool == "metabuli" ~ paste0("prec", prec, "_ms", ms, "_ptax", ptax)
         ),
         total_sec   = real_sec,
         peak_mem_kb = max_mem_kb) |>
  select(benchmark, tool, sample, param_combo, total_sec, peak_mem_kb)

totals <- bind_rows(single_total, cent_total, syl_total) |>
  left_join(read_totals, by = "sample") |>
  mutate(
    base_tool   = factor(tool, levels = tool_levels, labels = tool_labels[tool_levels]),
    total_min   = total_sec / 60,
    peak_mem_gb = peak_mem_kb / 1e6,
    # Distinguish EsViritu's two run modes on plots that show all tools.
    tool_disp   = if_else(tool == "esviritu",
                          paste0(tool_labels[["esviritu"]], " (", param_combo, ")"),
                          as.character(base_tool)),
    benchmark   = factor(benchmark, levels = c("mock", "negative1"),
                        labels = c("Mock community", "Negative control"))
  )

save_tab(totals, "timing_memory_all_runs.csv")

# ---- Summary table: median/mean runtime + memory per tool ----------------
summary_tbl <- totals |>
  group_by(benchmark, base_tool) |>
  summarise(
    n_runs        = n(),
    median_min    = median(total_min),
    mean_min      = mean(total_min),
    median_mem_gb = median(peak_mem_gb),
    mean_mem_gb   = mean(peak_mem_gb),
    .groups = "drop"
  ) |>
  arrange(benchmark, base_tool)
save_tab(summary_tbl, "timing_memory_summary_by_tool.csv")
print(summary_tbl, n = Inf)

# ---- Figure 1: total runtime by tool (across all samples & param sets) ---
p_time <- totals |>
  ggplot(aes(base_tool, total_min, colour = tolower(as.character(base_tool)))) +
  geom_boxplot(outlier.shape = NA, colour = "grey60") +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.6) +
  scale_y_log10() +
  scale_colour_manual(values = tool_palette, guide = "none") +
  facet_wrap(~ benchmark) +
  labs(x = NULL, y = "Wall-clock runtime, minutes (log scale)",
       title = "Total runtime by tool",
       subtitle = "Multi-step pipelines summed; spread reflects all parameter sets tested") +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_time, "timing_runtime_by_tool.pdf", w = 8, h = 5)

# ---- Figure 2: peak memory by tool ---------------------------------------
p_mem <- totals |>
  ggplot(aes(base_tool, peak_mem_gb, colour = tolower(as.character(base_tool)))) +
  geom_boxplot(outlier.shape = NA, colour = "grey60") +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.6) +
  scale_y_log10() +
  scale_colour_manual(values = tool_palette, guide = "none") +
  facet_wrap(~ benchmark) +
  labs(x = NULL, y = "Peak memory, GB (log scale)",
       title = "Peak memory usage by tool",
       subtitle = "Multi-step pipelines: max across steps; spread reflects all parameter sets tested") +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_mem, "timing_memory_by_tool.pdf", w = 8, h = 5)

# ---- Figure 3: runtime scaling with input read depth ---------------------
p_scale <- totals |>
  filter(!is.na(read_total)) |>
  ggplot(aes(read_total, total_min, colour = tolower(as.character(base_tool)))) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(aes(group = base_tool), method = "lm", se = FALSE, linewidth = 0.7) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10() +
  scale_colour_manual(values = tool_palette, name = "Tool",
                      labels = function(x) tool_labels[x]) +
  labs(x = "Input read pairs (log scale)", y = "Wall-clock runtime, minutes (log scale)",
       title = "Runtime scaling with input depth") +
  theme_bench()
save_fig(p_scale, "timing_runtime_vs_depth.pdf", w = 7.5, h = 5)

# ---- Figure 4: parameter sensitivity of runtime (multi-setting tools) ----
p_param_time <- totals |>
  filter(tool != "esviritu") |>
  ggplot(aes(param_combo, total_min, colour = tolower(as.character(base_tool)))) +
  geom_boxplot(outlier.shape = NA, colour = "grey60") +
  geom_jitter(width = 0.15, height = 0, size = 1.4, alpha = 0.6) +
  scale_y_log10() +
  scale_colour_manual(values = tool_palette, guide = "none") +
  facet_wrap(~ base_tool, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Wall-clock runtime, minutes (log scale)",
       title = "Runtime sensitivity to parameter choice") +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_fig(p_param_time, "timing_runtime_by_param.pdf", w = 10, h = 5)

# ---- Figure 5: parameter sensitivity of peak memory -----------------------
p_param_mem <- totals |>
  filter(tool != "esviritu") |>
  ggplot(aes(param_combo, peak_mem_gb, colour = tolower(as.character(base_tool)))) +
  geom_boxplot(outlier.shape = NA, colour = "grey60") +
  geom_jitter(width = 0.15, height = 0, size = 1.4, alpha = 0.6) +
  scale_y_log10() +
  scale_colour_manual(values = tool_palette, guide = "none") +
  facet_wrap(~ base_tool, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Peak memory, GB (log scale)",
       title = "Memory sensitivity to parameter choice") +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_fig(p_param_mem, "timing_memory_by_param.pdf", w = 10, h = 5)

message("Done: figures written to ", paths$out_fig, "; tables written to ", paths$out_tab)
