# Official Simulation (LOW RATE): SMC with different Update Interval lengths
#
# Same design as official_simulation_light.R but using the english_words_low_rate
# dataset (lower message rates: 30 and 1000 messages/hour).
#
# Configuration: same as LIGHT (1000 particles, 10000 iterations, 4000 burn-in, thinning 5).
# Report: results_explanation/official_simulation_low_rate.md
#
# Prerequisite: run tests/english_words_low_rate_sim_data.R from package root so that
# english_words_low_rate is available in data/.

# Clear environment
rm(list=ls())
set.seed(32453)

# Load required libraries
devtools::load_all()
library(mclust)
library(fitdistrplus)

# Load the "english_words_low_rate" data
dat <- english_words_low_rate
Bvec <- dat$Bvec
Vvec <- dat$Vvec
Zvec <- dat$Zvec
Qvec <- dat$Qvec
Fvec <- dat$Fvec

Tvec <- dat$Tvec
Yvec <- dat$Yvec

# Prepare time series data
ts_data <- list()
ts_data$Yvec <- Yvec
ts_data$Tvec <- Tvec

# Prepare parameters (same for all runs)
parameters <- list()
parameters$U <- dat$U
parameters$W <- dat$W
parameters$K <- dat$K
parameters$lambdamat <- dat$lambdamat
parameters$keyvec <- dat$keyvec
parameters$etavec <- dat$etavec
parameters$key0vec <- dat$key0vec
parameters$eta0vec <- dat$eta0vec
parameters$alphavec <- dat$alphavec
parameters$muvec <- dat$muvec
parameters$probvec_V <- dat$probvec_V
parameters$probvec_Z <- dat$probvec_Z
parameters$probvec_Q <- dat$probvec_Q
parameters$probvec_F <- dat$probvec_F
parameters$P0 <- dat$P0
parameters$minimum_n <- dat$minimum_n

# Prepare base settings (LIGHT config; length_UI will vary)
base_settings <- list()
base_settings$num_logs <- dat$num_logs
base_settings$Jss1 <- 1/3
base_settings$Js1s <- 1/3
base_settings$Smax <- 150
base_settings$thinning <- 5
base_settings$method <- "turcotte"

# Define the different length_UI values to test
length_UI_values <- c(0.5, 1, 2, 4)

# Per-run SMC configuration (indexed to match length_UI_values)
n_ite_values      <- c(1000, 4000, 16000, 64000)
n_particle_values <- c(100,   400,  1000,  1000)
burn_in_values    <- c(400,  1600,  6400, 26000)

n_runs <- length(length_UI_values)

# Initialize storage for results
results_summary <- list()

cat("========================================\n")
cat("OFFICIAL SIMULATION (LOW RATE): Starting runs\n")
cat("========================================\n\n")

# Run simulation for each length_UI value
for(run_idx in 1:n_runs) {
  length_UI <- length_UI_values[run_idx]

  cat(sprintf("Run %d/%d: length_UI = %.1f\n", run_idx, n_runs, length_UI))
  cat("----------------------------------------\n")

  # Set length_UI and per-run SMC parameters
  settings <- base_settings
  settings$length_UI  <- length_UI
  settings$n_particle <- n_particle_values[run_idx]
  settings$n_ite      <- n_ite_values[run_idx]
  settings$burn_in    <- burn_in_values[run_idx]

  # Run SMC (C++) then post-processing into RJSMC
  cat("  Running SMC algorithm...\n")
  start_time <- Sys.time()
  out_cpp <- RJSMC::SMC(ts_data = ts_data,
                        parameters = parameters,
                        settings = settings)
  end_time <- Sys.time()
  elapsed_cpp <- as.numeric(difftime(end_time, start_time, units = "secs"))
  out_SMC <- RJSMC::smc_post_processing(out_cpp,
                                        parameters = parameters,
                                        settings = settings,
                                        interval_length = 0.01,
                                        elapsed_time = elapsed_cpp)
  elapsed_time <- out_SMC@elapsed_time
  cat(sprintf("  SMC completed in %.2f seconds\n", elapsed_time))

  # Run breakpoint_cx evaluation
  cat("  Running breakpoint_cx evaluation...\n")
  evaluation_results <- breakpoint_cx(
    obj = out_SMC,
    true_Bvec = dat$Bvec,
    Tvec = Tvec,
    U = dat$U,
    storage_weight = out_SMC@storage_weight,
    UI_bounds = out_SMC@UI_bounds,
    Vvec = dat$Vvec
  )

  # Store key results
  results_summary[[run_idx]] <- list(
    run_number = run_idx,
    length_UI = length_UI,
    elapsed_time_seconds = elapsed_time,
    n_UI = length(out_SMC@UI_bounds) - 1,
    n_particles = settings$n_particle,
    n_intervals = nrow(evaluation_results$intervals),
    n_observations = length(Tvec),

    breakpoint_confusion_matrix = evaluation_results$confusion_matrix,
    V_state_confusion_matrix = evaluation_results$V_confusion_matrix,
    avg_error_breakpoints = evaluation_results$avg_error_breakpoints,
    avg_error_V_state = evaluation_results$avg_error_V_state,

    breakpoint_correct = sum(diag(evaluation_results$confusion_matrix)),
    breakpoint_total = sum(evaluation_results$confusion_matrix),
    breakpoint_accuracy = sum(diag(evaluation_results$confusion_matrix)) / sum(evaluation_results$confusion_matrix),

    V_state_correct = if(!is.null(evaluation_results$V_confusion_matrix)) sum(diag(evaluation_results$V_confusion_matrix)) else NA,
    V_state_total = if(!is.null(evaluation_results$V_confusion_matrix)) sum(evaluation_results$V_confusion_matrix) else NA,
    V_state_accuracy = if(!is.null(evaluation_results$V_confusion_matrix)) sum(diag(evaluation_results$V_confusion_matrix)) / sum(evaluation_results$V_confusion_matrix) else NA
  )

  cat(sprintf("  Run %d completed successfully\n\n", run_idx))
}

cat("========================================\n")
cat("All runs completed!\n")
cat("========================================\n\n")

# Generate markdown report
cat("Generating markdown report...\n")

if(dir.exists("tests")) {
  results_dir <- "results_explanation"
} else {
  results_dir <- "../results_explanation"
}

if(!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Build markdown report (same structure as official_simulation_light)
build_official_simulation_md <- function(results_summary, base_settings, n_obs, n_Bvec, U) {
  n_runs <- length(results_summary)
  md <- character(0)
  md[1] <- "# Official Simulation Results (Low Rate)\n\n"
  md[1] <- paste0(md[1], "This document summarizes the results of running the SMC algorithm with the **english_words_low_rate** dataset (30 and 1000 messages/hour) and different Update Interval lengths.\n\n")
  md[1] <- paste0(md[1], "## Simulation Setup\n\n")
  md[1] <- paste0(md[1], "- **Dataset**: english_words_low_rate\n")
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
  md[1] <- paste0(md[1], "| length_UI | Accuracy | Avg Error |\n|-----------|----------|----------|\n")
  for (run_idx in seq_len(n_runs)) {
    res <- results_summary[[run_idx]]
    md[1] <- paste0(md[1], sprintf("| %.1f | %.4f | %.4f |\n", res$length_UI, res$breakpoint_accuracy, res$avg_error_breakpoints))
  }
  md[1] <- paste0(md[1], "\n### V State Classification Performance\n\n| length_UI | Accuracy | Avg Error |\n|-----------|----------|----------|\n")
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
  md[1] <- paste0(md[1], "- **Dataset**: english_words_low_rate (30 and 1000 messages/hour).\n")
  md[1] <- paste0(md[1], "- **Breakpoint Confusion Matrix**: 3x3 matrix (rows = true, columns = estimated breakpoint count 0/1/2).\n")
  md[1] <- paste0(md[1], sprintf("- **V State Confusion Matrix**: %dx%d matrix (rows = true, columns = estimated V state).\n", U, U))
  md[1] <- paste0(md[1], "- **Average Error**: Proportion of incorrect predictions.\n\n")
  md[1]
}

md_content <- build_official_simulation_md(
  results_summary,
  base_settings,
  n_obs = length(Tvec),
  n_Bvec = length(dat$Bvec),
  U = dat$U
)

# Write markdown file (low-rate run)
md_file_path <- file.path(results_dir, "official_simulation_low_rate.md")
writeLines(md_content, md_file_path)

cat(sprintf("Markdown report generated: %s\n", md_file_path))
cat("\nSimulation (LOW RATE) completed successfully!\n")
