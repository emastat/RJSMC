# Test markdown generation for official_simulation.R
# Run this to verify the report is compact and matrices are correct (no full simulation).

# Source only the parts we need: define build function and create mock results
devtools::load_all()

# Mock results_summary (4 runs, minimal matrices)
n_runs <- 4L
base_settings <- list(
  n_particle = 5000L,
  n_ite = 50000L,
  burn_in = 20000L,
  thinning = 5L,
  method = "turcotte"
)
results_summary <- list()
for (run_idx in 1:n_runs) {
  cm_bp <- matrix(c(100L, 2L, 0L, 1L, 10L, 0L, 0L, 0L, 0L), 3L, 3L,
    dimnames = list(c("True_0", "True_1", "True_2"), c("Est_0", "Est_1", "Est_2")))
  cm_v <- matrix(0L, 5L, 5L)
  diag(cm_v) <- c(50L, 30L, 20L, 40L, 60L)
  dimnames(cm_v) <- list(
    paste0("True_V_", 1:5),
    paste0("Est_V_", 1:5)
  )
  results_summary[[run_idx]] <- list(
    run_number = run_idx,
    length_UI = c(0.5, 1, 2, 5)[run_idx],
    elapsed_time_seconds = 100 + run_idx * 20,
    n_UI = 50L - run_idx * 10,
    n_particles = 5000L,
    n_intervals = 1000L - run_idx * 50,
    n_observations = 4683L,
    breakpoint_confusion_matrix = cm_bp,
    V_state_confusion_matrix = cm_v,
    avg_error_breakpoints = 0.01 * run_idx,
    avg_error_V_state = 0.01 * run_idx,
    breakpoint_correct = sum(diag(cm_bp)),
    breakpoint_total = sum(cm_bp),
    breakpoint_accuracy = sum(diag(cm_bp)) / sum(cm_bp),
    V_state_correct = sum(diag(cm_v)),
    V_state_total = sum(cm_v),
    V_state_accuracy = sum(diag(cm_v)) / sum(cm_v)
  )
}

# Build markdown using the same function as official_simulation.R
# (copy of the function so this test is self-contained and can run without running SMC)
build_official_simulation_md <- function(results_summary, base_settings, n_obs, n_Bvec, U) {
  n_runs <- length(results_summary)
  md <- character(0)
  md[1] <- "# Official Simulation Results\n\n"
  md[1] <- paste0(md[1], "This document summarizes the results of running the SMC algorithm with different Update Interval lengths.\n\n")
  md[1] <- paste0(md[1], "## Simulation Setup\n\n")
  md[1] <- paste0(md[1], "- **Number of particles**: ", base_settings$n_particle, "\n")
  md[1] <- paste0(md[1], "- **Number of iterations**: ", base_settings$n_ite, "\n")
  md[1] <- paste0(md[1], "- **Burn-in**: ", base_settings$burn_in, "\n")
  md[1] <- paste0(md[1], "- **Thinning**: ", base_settings$thinning, "\n")
  md[1] <- paste0(md[1], "- **Method**: ", base_settings$method, "\n")
  md[1] <- paste0(md[1], "- **Total observations**: ", n_obs, "\n")
  md[1] <- paste0(md[1], "- **Number of true breakpoints**: ", n_Bvec, "\n")
  md[1] <- paste0(md[1], "- **Number of V states**: ", U, "\n\n")
  md[1] <- paste0(md[1], "## Results Summary\n\n### Performance Summary Table\n\n")
  md[1] <- paste0(md[1], "| length_UI | Run Time (s) | Breakpoint Accuracy | V State Accuracy | Avg Error (Breakpoints) | Avg Error (V State) |\n")
  md[1] <- paste0(md[1], "|-----------|---------------|---------------------|------------------|-------------------------|---------------------|\n")
  for (run_idx in seq_len(n_runs)) {
    res <- results_summary[[run_idx]]
    bp_acc <- sprintf("%.4f", res$breakpoint_accuracy)
    v_acc <- if (!is.na(res$V_state_accuracy)) sprintf("%.4f", res$V_state_accuracy) else "N/A"
    bp_err <- sprintf("%.4f", res$avg_error_breakpoints)
    v_err <- if (!is.na(res$avg_error_V_state)) sprintf("%.4f", res$avg_error_V_state) else "N/A"
    md[1] <- paste0(md[1], sprintf("| %.1f | %.2f | %s | %s | %s | %s |\n",
      res$length_UI, res$elapsed_time_seconds, bp_acc, v_acc, bp_err, v_err))
  }
  md[1] <- paste0(md[1], "\n")
  for (run_idx in seq_len(n_runs)) {
    res <- results_summary[[run_idx]]
    md[1] <- paste0(md[1], sprintf("## Run %d: length_UI = %.1f\n\n### Run Configuration\n\n", run_idx, res$length_UI))
    md[1] <- paste0(md[1], sprintf("- **Update Interval length**: %s\n- **Number of Update Intervals**: %d\n- **Number of particles**: %d\n- **Number of evaluation intervals**: %d\n- **Number of observations**: %d\n- **Computation time**: %.2f seconds\n\n",
      res$length_UI, res$n_UI, res$n_particles, res$n_intervals, res$n_observations, res$elapsed_time_seconds))
    md[1] <- paste0(md[1], "### Breakpoint Detection Results\n\n**Confusion Matrix** (rows = true, columns = estimated):\n\n")
    bp_cm_lines <- capture.output(print(res$breakpoint_confusion_matrix))
    md[1] <- paste0(md[1], "```\n", paste(bp_cm_lines, collapse = "\n"), "\n```\n\n")
    md[1] <- paste0(md[1], sprintf("- **Total intervals evaluated**: %d\n- **Correct predictions**: %d\n- **Accuracy**: %.4f (%.2f%%)\n- **Average error**: %.4f\n\n",
      res$breakpoint_total, res$breakpoint_correct, res$breakpoint_accuracy, res$breakpoint_accuracy * 100, res$avg_error_breakpoints))
    md[1] <- paste0(md[1], "**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).\n\n")
    if (!is.null(res$V_state_confusion_matrix)) {
      md[1] <- paste0(md[1], "### V State Classification Results\n\n**Confusion Matrix** (rows = true, columns = estimated):\n\n")
      v_cm_lines <- capture.output(print(res$V_state_confusion_matrix))
      md[1] <- paste0(md[1], "```\n", paste(v_cm_lines, collapse = "\n"), "\n```\n\n")
      md[1] <- paste0(md[1], sprintf("- **Total observations evaluated**: %d\n- **Correct predictions**: %d\n- **Accuracy**: %.4f (%.2f%%)\n- **Average error**: %.4f\n\n",
        res$V_state_total, res$V_state_correct, res$V_state_accuracy, res$V_state_accuracy * 100, res$avg_error_V_state))
      md[1] <- paste0(md[1], sprintf("**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to %d).\n\n", U))
    } else {
      md[1] <- paste0(md[1], "### V State Classification Results\n\n*V state evaluation was not performed for this run.*\n\n")
    }
    md[1] <- paste0(md[1], "---\n\n")
  }
  md[1] <- paste0(md[1], "## Comparative Analysis\n\n### Breakpoint Detection Performance\n\n")
  md[1] <- paste0(md[1], "| length_UI | Accuracy | Avg Error |\n|-----------|----------|-----------|\n")
  for (run_idx in seq_len(n_runs)) {
    res <- results_summary[[run_idx]]
    md[1] <- paste0(md[1], sprintf("| %.1f | %.4f | %.4f |\n", res$length_UI, res$breakpoint_accuracy, res$avg_error_breakpoints))
  }
  md[1] <- paste0(md[1], "\n### V State Classification Performance\n\n| length_UI | Accuracy | Avg Error |\n|-----------|----------|-----------|\n")
  for (run_idx in seq_len(n_runs)) {
    res <- results_summary[[run_idx]]
    if (!is.na(res$V_state_accuracy)) md[1] <- paste0(md[1], sprintf("| %.1f | %.4f | %.4f |\n", res$length_UI, res$V_state_accuracy, res$avg_error_V_state))
  }
  md[1] <- paste0(md[1], "\n### Computational Performance\n\n| length_UI | Run Time (seconds) | Update Intervals |\n|-----------|---------------------|------------------|\n")
  for (run_idx in seq_len(n_runs)) {
    res <- results_summary[[run_idx]]
    md[1] <- paste0(md[1], sprintf("| %.1f | %.2f | %d |\n", res$length_UI, res$elapsed_time_seconds, res$n_UI))
  }
  md[1] <- paste0(md[1], "\n## Notes\n\n")
  md[1] <- paste0(md[1], "- **Breakpoint Confusion Matrix**: 3x3 matrix (rows = true, columns = estimated breakpoint count 0/1/2).\n")
  md[1] <- paste0(md[1], sprintf("- **V State Confusion Matrix**: %dx%d matrix (rows = true, columns = estimated V state).\n", U, U))
  md[1] <- paste0(md[1], "- **Average Error**: Proportion of incorrect predictions.\n\n")
  md[1]
}

md_content <- build_official_simulation_md(
  results_summary,
  base_settings,
  n_obs = 4683L,
  n_Bvec = 52L,
  U = 5L
)

# Assertions
stopifnot("Markdown should be a single string" = length(md_content) == 1L)
nc <- nchar(md_content)
stopifnot("Markdown should be compact (under 25k chars for 4 runs)" = nc < 25000L)
title_count <- lengths(regmatches(md_content, gregexpr("# Official Simulation Results", md_content, fixed = TRUE)))[1]
stopifnot("Title should appear exactly once" = title_count == 1L)
# First breakpoint code block: should contain matrix rows (True_0, True_1, True_2) and no duplicate title inside
first_bp_block <- sub(".*?```\\n(.*?)```.*", "\\1", md_content)
stopifnot("First code block should contain True_0" = grepl("True_0", first_bp_block, fixed = TRUE))
stopifnot("First code block should contain Est_0" = grepl("Est_0", first_bp_block, fixed = TRUE))
stopifnot("Code block must not contain repeated full document" = !grepl("# Official Simulation Results", first_bp_block, fixed = TRUE))

cat("All markdown generation tests passed.\n")
cat("  - Single string output\n")
cat("  - Length:", nc, "characters (compact)\n")
cat("  - Title appears once\n")
cat("  - Breakpoint matrix block is correct (True_0, Est_0, no duplication)\n")
# Write sample to file for visual check
if (!dir.exists("results_explanation")) dir.create("results_explanation", recursive = TRUE)
writeLines(md_content, "results_explanation/test_markdown_sample.md")
cat("Sample report written to results_explanation/test_markdown_sample.md\n")
cat("\nYou can run the full simulation with: source(\"tests/official_simulation.R\")\n")
