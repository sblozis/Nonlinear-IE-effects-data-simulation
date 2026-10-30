# ============================================================
# Simulation 2: Fixed Exponential-Decay MELS Model
# The initial-elevation (IE) effect decays exponentially but is fixed
# level-2 random-effects block is 2 x 2 (u_b0, u_sigma).
# Correct model: fixed exponential decay, random intercept, random scale
# Misspecified model: decay replaced by a single-day indicator at time = 1
# ============================================================

RESET_CACHE <- FALSE

checkpoint_dirs <- c(
  "sim_expDecay_checkpoints_correct_logsd",
  "sim_expDecay_checkpoints_wrong_fixedIE_logsd"
)

if (RESET_CACHE) {
  unlink(checkpoint_dirs, recursive = TRUE)
}

tmp_root <- file.path(tempdir(), "fixed_decay_cmdstan_tmp")
dir.create(tmp_root, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_root)

n_chains       <- 4L
n_warmup       <- 2000L
n_iter         <- 8000L  
n_warmup_refit <- 2000L
n_iter_refit   <- 10000L 
n_reps         <- 250L

divergence_refit_threshold <- 5L

fit_stages <- list(
  list(iter = n_iter,       warmup = n_warmup,       adapt_delta = 0.98,  max_treedepth = 12L, seed_base = 123L),
  list(iter = n_iter_refit, warmup = n_warmup_refit, adapt_delta = 0.99,  max_treedepth = 15L, seed_base = 10000L),
  list(iter = 16000L,       warmup = 4000L,          adapt_delta = 0.995, max_treedepth = 15L, seed_base = 20000L),
  list(iter = 24000L,       warmup = 6000L,          adapt_delta = 0.995, max_treedepth = 15L, seed_base = 30000L)
)
max_stage <- length(fit_stages) - 1L

library(brms)
library(cmdstanr)
library(dplyr)
library(purrr)
library(furrr)
library(future)
library(tibble)
library(posterior)

set.seed(123)

# True parameter values 
# Model:
#   y_ij         = b0 + u_b0,i + b1 * exp(-lambda * time_ij) + e_ij,
#                  e_ij ~ N(0, sigma_i^2)
#   log(sigma_i) = log_sigma0 + u_sigma,i         
#   (u_b0, u_sigma)' ~ MVN(0, Sigma),  Sigma = D R D
#     D = diag(sd_int, sd_scale),  R = [[1, cor_int_scale], [cor_int_scale, 1]]
true_vals <- c(
  b0            = 2.5,
  b1            = 0.5,    
  lambda        = 1.15,   
  log_sigma0    = -0.75,
  sd_int        = sqrt(0.7),
  sd_scale      = sqrt(0.125),
  cor_int_scale = -0.20
)

re_terms <- c("int", "scale")

build_R <- function(tv) {
  R <- diag(2L)
  dimnames(R) <- list(re_terms, re_terms)
  R["int", "scale"] <- R["scale", "int"] <- tv[["cor_int_scale"]]
  R
}

check_corr_matrix <- function(R, tol = 1e-8) {
  stopifnot(is.matrix(R), nrow(R) == ncol(R), isSymmetric(unname(R)),
            all(abs(diag(R) - 1) < 1e-12), all(abs(R[upper.tri(R)]) < 1))
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= tol) {
    stop(sprintf("Correlation matrix is not positive definite (min eigenvalue %.3g).", min(ev)))
  }
  invisible(ev)
}

build_Sigma <- function(tv) {
  sds <- c(tv[["sd_int"]], tv[["sd_scale"]])
  if (any(!is.finite(sds)) || any(sds <= 0)) {
    stop("All random-effect SDs must be finite and strictly positive.")
  }
  R <- build_R(tv)
  check_corr_matrix(R)
  D <- diag(sds)
  dimnames(D) <- list(re_terms, re_terms)
  D %*% R %*% D
}

local({
  Sigma <- build_Sigma(true_vals)
  sds   <- c(true_vals[["sd_int"]], true_vals[["sd_scale"]])
  stopifnot(isTRUE(all.equal(unname(sqrt(diag(Sigma))), sds, tolerance = 1e-12)))
  stopifnot(isTRUE(all.equal(unname(stats::cov2cor(Sigma)),
                             unname(build_R(true_vals)), tolerance = 1e-12)))
  R_chk        <- matrix(c(1, -0.20, -0.20, 1), nrow = 2L)
  D_logvar     <- diag(c(sqrt(0.7), sqrt(0.5)))  
  T_half       <- diag(c(1, 1 / 2))
  Sigma_logvar <- D_logvar %*% R_chk %*% D_logvar
  stopifnot(isTRUE(all.equal(unname(T_half %*% Sigma_logvar %*% T_half),
                             unname(Sigma), tolerance = 1e-12)))
  stopifnot(isTRUE(all.equal(true_vals[["log_sigma0"]], -1.5 / 2)))
})

true_Sigma_vals <- local({
  S    <- build_Sigma(true_vals)
  keep <- which(upper.tri(S, diag = TRUE), arr.ind = TRUE)
  stats::setNames(S[keep], paste0("Sigma_", re_terms[keep[, 1]], "_", re_terms[keep[, 2]]))
})

n_clusters    <- 300L
n_per_cluster <- 10L

simulate_rep <- function(rep_id, n_clusters, n_per_cluster, tv,
                         seed_offset = 0L) {
  set.seed(123 + seed_offset + rep_id)

  Sigma <- build_Sigma(tv)  
  u <- MASS::mvrnorm(n = n_clusters, mu = c(0, 0), Sigma = Sigma)

  cluster <- rep(seq_len(n_clusters), each = n_per_cluster)
  time    <- rep(0:(n_per_cluster - 1L), times = n_clusters)
  day1    <- as.numeric(time == 0L)   

  b0_i    <- tv[["b0"]] + u[cluster, 1L]              
  sigma_i <- exp(tv[["log_sigma0"]] + u[cluster, 2L])    
  mu_i    <- b0_i + tv[["b1"]] * exp(-tv[["lambda"]] * time)

  y <- stats::rnorm(n_clusters * n_per_cluster, mean = mu_i, sd = sigma_i)

  data.frame(y = y, cluster = cluster, time = time, day1 = day1)
}

safe_max <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

get_diag_vals <- function(fit_i) {
  n_divergent <- tryCatch(
    {
      np <- brms::nuts_params(fit_i)
      sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
    },
    error = function(e) NA_integer_
  )

  draws <- posterior::as_draws_df(fit_i, variable = "^(b_|sd_|cor_)", regex = TRUE)
  population_names <- posterior::variables(draws)

  if (length(population_names) == 0L) {
    return(list(n_divergent = n_divergent, max_rhat = NA_real_,
                min_ess_bulk = NA_real_, min_ess_tail = NA_real_))
  }

  summaries <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")

  low_parameters <- summaries |>
    dplyr::filter(ess_bulk < 400 | ess_tail < 400) |>
    dplyr::pull(variable)
  if (length(low_parameters) > 0L) {
    message("  [Warning] ESS below 400 for: ", paste(low_parameters, collapse = ", "))
  }

  list(
    n_divergent  = n_divergent,
    max_rhat     = safe_max(summaries$rhat),
    min_ess_bulk = safe_min(summaries$ess_bulk),
    min_ess_tail = safe_min(summaries$ess_tail)
  )
}

needs_refit <- function(diagnostics) {
  anyNA(unlist(diagnostics, use.names = FALSE)) ||
    diagnostics$n_divergent > divergence_refit_threshold ||
    diagnostics$max_rhat > 1.01 ||
    diagnostics$min_ess_bulk < 400 ||
    diagnostics$min_ess_tail < 400
}

get_correlation_draw <- function(draws, term1, term2) {
  candidates <- c(
    paste0("cor_cluster__", term1, "__", term2),
    paste0("cor_cluster__", term2, "__", term1)
  )
  found <- intersect(candidates, names(draws))
  if (length(found) != 1L) {
    stop("Could not uniquely identify the brms correlation parameter for ",
         term1, " and ", term2, ". Available correlation parameters: ",
         paste(grep("^cor_", names(draws), value = TRUE), collapse = ", "))
  }
  draws[[found]]
}

primary_parameters <- names(true_vals)

model_specs <- list(
  correct = list(
    fixed     = c(b0 = "b_b0_Intercept", b1 = "b_b1_Intercept",
                  lambda = "b_lambda_Intercept", log_sigma0 = "b_sigma_Intercept"),
    brms_term = c(int = "b0_Intercept", scale = "sigma_Intercept"),
    notes     = character(0)
  ),
  wrong_fixed = list(
    fixed     = c(b0 = "b_Intercept", b1 = "b_day1", log_sigma0 = "b_sigma_Intercept"),
    brms_term = c(int = "Intercept", scale = "sigma_Intercept"),
    notes     = c(
      b1     = "estimated as the day-1 indicator coefficient (not the decay amplitude)",
      lambda = "not estimated: the day-1 model has no decay rate"
    )
  )
)

expected_variables <- function(model_type) {
  spec <- model_specs[[model_type]]
  keep <- which(upper.tri(diag(length(re_terms)), diag = TRUE), arr.ind = TRUE)
  c(
    intersect(primary_parameters,
              c(names(spec$fixed), "sd_int", "sd_scale", "cor_int_scale")),
    paste0("Sigma_", re_terms[keep[, 1]], "_", re_terms[keep[, 2]])
  )
}

summarise_draw_vector <- function(values, rep_id, parameter, type) {
  tibble::tibble(rep = rep_id, type = type, variable = parameter,
                 mean = mean(values), median = stats::median(values),
                 sd = stats::sd(values),
                 lower_2.5  = unname(stats::quantile(values, 0.025)),
                 upper_97.5 = unname(stats::quantile(values, 0.975)))
}

extract_draws <- function(fit, rep_id, model_type) {
  spec <- model_specs[[model_type]]
  if (is.null(spec)) stop("Unknown model_type: ", model_type)

  draws <- as.data.frame(posterior::as_draws_df(
    fit, variable = "^(b_|sd_cluster__|cor_cluster__)", regex = TRUE))

  required <- c(unname(spec$fixed), paste0("sd_cluster__", spec$brms_term[re_terms]))
  missing <- setdiff(required, names(draws))
  if (length(missing) > 0L) {
    stop("Required posterior draw(s) not found: ", paste(missing, collapse = ", "))
  }

  for (p in names(spec$fixed)) draws[[p]] <- draws[[spec$fixed[[p]]]]
  for (term in re_terms) {
    draws[[paste0("sd_", term)]] <- draws[[paste0("sd_cluster__", spec$brms_term[[term]])]]
  }
  draws$cor_int_scale <- get_correlation_draw(
    draws, spec$brms_term[["int"]], spec$brms_term[["scale"]])
  parameters <- intersect(primary_parameters, expected_variables(model_type))
  k       <- length(re_terms)
  sd_mat  <- as.matrix(draws[paste0("sd_", re_terms)])
  cor_vec <- draws$cor_int_scale
  keep    <- which(upper.tri(diag(k), diag = TRUE), arr.ind = TRUE)
  sigma_names <- paste0("Sigma_", re_terms[keep[, 1]], "_", re_terms[keep[, 2]])
  sigma_draws <- t(vapply(seq_len(nrow(draws)), function(d) {
    R <- matrix(c(1, cor_vec[d], cor_vec[d], 1), nrow = 2L)
    D <- diag(sd_mat[d, ], nrow = k)
    (D %*% R %*% D)[keep]
  }, numeric(nrow(keep))))
  colnames(sigma_draws) <- sigma_names

  dplyr::bind_rows(
    purrr::map_dfr(parameters, function(p) summarise_draw_vector(draws[[p]], rep_id, p, "primary")),
    purrr::map_dfr(sigma_names, function(p) summarise_draw_vector(sigma_draws[, p], rep_id, p, "derived_Sigma"))
  )
}

summarise_sim <- function(all_parameters, model_label, model_type) {
  notes <- model_specs[[model_type]]$notes
  truth <- dplyr::bind_rows(
    tibble::tibble(type = "primary", variable = names(true_vals), true_value = unname(true_vals)),
    tibble::tibble(type = "derived_Sigma", variable = names(true_Sigma_vals),
                   true_value = unname(true_Sigma_vals)))
  unmatched <- dplyr::anti_join(dplyr::distinct(all_parameters, type, variable), truth,
                                by = c("type", "variable"))
  if (nrow(unmatched) > 0L) stop("Estimates without a true value: ",
                                 paste(unmatched$variable, collapse = ", "))
  estimates <- all_parameters |>
    dplyr::inner_join(truth, by = c("type", "variable")) |>
    dplyr::group_by(type, variable) |>
    dplyr::summarise(
      n_reps        = dplyr::n(),
      avg_est       = mean(rep_mean, na.rm = TRUE),
      bias          = avg_est - true_value[1L],
      rmse          = sqrt(mean((rep_mean - true_value)^2, na.rm = TRUE)),
      emp_sd        = stats::sd(rep_mean, na.rm = TRUE),
      avg_post_sd   = mean(rep_sd, na.rm = TRUE),
      coverage_95   = mean(true_value >= lower_2.5 & true_value <= upper_97.5, na.rm = TRUE),
      avg_median    = mean(rep_median, na.rm = TRUE),
      bias_median   = avg_median - true_value[1L],
      rmse_median   = sqrt(mean((rep_median - true_value)^2, na.rm = TRUE)),
      mcse_bias     = emp_sd / sqrt(n_reps),
      mcse_coverage = sqrt(coverage_95 * (1 - coverage_95) / n_reps),
      .groups = "drop")
  truth |>
    dplyr::left_join(estimates, by = c("type", "variable")) |>
    dplyr::filter(type == "primary" | !is.na(n_reps)) |>
    dplyr::mutate(
      model  = model_label,
      n_reps = dplyr::coalesce(n_reps, 0L),
      status = dplyr::case_when(
        variable %in% names(notes) ~ unname(notes[variable]),
        n_reps > 0L ~ "estimated",
        TRUE ~ "not estimated"
      )
    ) |>
    dplyr::select(model, type, variable, status, true_value, n_reps, avg_est, bias, rmse,
                  emp_sd, avg_post_sd, coverage_95, avg_median, bias_median, rmse_median,
                  mcse_bias, mcse_coverage)
}

bform_correct <- brms::bf(
  y ~ b0 + b1 * exp(-lambda * time),
  b0 ~ 1 + (1 |q| cluster),
  b1 ~ 1,
  lambda ~ 1,
  sigma ~ 1 + (1 |q| cluster),
  nl = TRUE
)

priors_correct <- c(
  prior(normal(0, 5), nlpar = "b0"),
  prior(normal(0, 2), nlpar = "b1"),
  prior(normal(1, 1), nlpar = "lambda", lb = 0),
  prior(normal(0, 1.5), class = "Intercept", dpar = "sigma"),
  prior(student_t(3, 0, 1), class = "sd", nlpar = "b0"),
  prior(student_t(3, 0, 1), class = "sd", dpar = "sigma"),
  prior(lkj(1), class = "cor")
)

bform_wrong_fixedIE <- brms::bf(
  y ~ day1 + (1 |q| cluster),
  sigma ~ 1 + (1 |q| cluster)
)

priors_wrong_fixedIE <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b", coef = "day1"),
  prior(normal(0, 1.5), class = "Intercept", dpar = "sigma"),
  prior(student_t(3, 0, 1), class = "sd"),
  prior(student_t(3, 0, 1), class = "sd", dpar = "sigma"),
  prior(lkj(1), class = "cor")
)

run_full_pipeline_parallel <- function(model_label, checkpoint_dir,
                                       bform, priors, model_type) {
  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  config_id <- rlang::hash(list(
    true_vals = true_vals, design = c(n_clusters, n_per_cluster),
    mcmc = c(n_chains, n_warmup, n_iter, n_warmup_refit, n_iter_refit, divergence_refit_threshold),
    formula = utils::capture.output(print(bform)), priors = as.data.frame(priors),
    model_type = model_type,
    output = "sd-cor-median-Sigma-v2"))
  expected_vars <- expected_variables(model_type)

  message("=== ", model_label, ": compiling Stan model once ===")
  compilation_data <- simulate_rep(
    rep_id = 1L, n_clusters = n_clusters,
    n_per_cluster = n_per_cluster, tv = true_vals
  )

  compiled_fit <- brms::brm(
    formula = bform, data = compilation_data, prior = priors,
    chains = 1L, iter = 20L, warmup = 10L, cores = 1L,
    backend = "cmdstanr", seed = 123, refresh = 0, init = 0,
    control = list(adapt_delta = 0.98, max_treedepth = 12)
  )

  fit_at_stage <- function(dat, i, stage) {
    s <- fit_stages[[stage + 1L]]
    tryCatch(
      update(
        compiled_fit, newdata = dat, recompile = FALSE,
        chains = n_chains, iter = s$iter, warmup = s$warmup,
        cores = 1L, seed = s$seed_base + i, refresh = 0, init = 0,
        control = list(adapt_delta = s$adapt_delta, max_treedepth = s$max_treedepth)
      ),
      error = function(e) {
        message("  Replication ", i, ", stage ", stage, " failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  diag_from_cache <- function(cached) {
    list(n_divergent  = cached$n_divergent[1L],  max_rhat     = cached$max_rhat[1L],
         min_ess_bulk = cached$min_ess_bulk[1L], min_ess_tail = cached$min_ess_tail[1L])
  }

  runner <- function(i) {
    on.exit(unlink(list.files(tempdir(), pattern = "\\.csv$", full.names = TRUE)),
            add = TRUE)

    checkpoint_file <- file.path(checkpoint_dir, sprintf("rep%03d.rds", i))
    cached <- NULL
    start_stage <- 0L

    if (file.exists(checkpoint_file)) {
      cached <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
      required_cache_columns <- c(
        "rep", "type", "variable", "mean", "median", "sd",
        "lower_2.5", "upper_97.5", "n_divergent", "max_rhat",
        "min_ess_bulk", "min_ess_tail", "was_refit", "config_id"
      )

      cache_valid <- !is.null(cached) &&
        nrow(cached) > 0L &&
        all(required_cache_columns %in% names(cached)) &&
        all(cached$rep == i) &&
        all(cached$config_id == config_id) &&
        setequal(cached$variable, expected_vars) &&
        !all(is.na(cached$max_rhat))

      if (cache_valid) {
        if (!"fit_stage" %in% names(cached)) {
          cached$fit_stage <- if (isTRUE(cached$was_refit[1L])) 1L else 0L
        }
        cached$converged <- !needs_refit(diag_from_cache(cached))
        if (cached$converged[1L] || cached$fit_stage[1L] >= max_stage) {
          return(cached)
        }
        start_stage <- cached$fit_stage[1L] + 1L
        message(sprintf(
          "  Replication %d: cached fit (stage %d, max Rhat=%s) fails diagnostics; escalating to stage %d",
          i, cached$fit_stage[1L], format(cached$max_rhat[1L], digits = 5), start_stage
        ))
      } else {
        cached <- NULL
      }
    }

    dat <- simulate_rep(
      rep_id = i, n_clusters = n_clusters,
      n_per_cluster = n_per_cluster, tv = true_vals
    )

    fit_i <- NULL
    diagnostics <- NULL
    final_stage <- NA_integer_

    for (stage in seq.int(start_stage, max_stage)) {
      candidate <- fit_at_stage(dat, i, stage)
      if (is.null(candidate)) next
      fit_i <- candidate
      diagnostics <- get_diag_vals(fit_i)
      final_stage <- stage
      if (!needs_refit(diagnostics)) break
      if (stage < max_stage) {
        message(sprintf(
          paste0("  Replication %d, stage %d fails diagnostics: divergences=%s, ",
                 "max Rhat=%s, bulk ESS=%s, tail ESS=%s; escalating to stage %d"),
          i, stage, diagnostics$n_divergent,
          format(diagnostics$max_rhat, digits = 5),
          format(diagnostics$min_ess_bulk, digits = 6),
          format(diagnostics$min_ess_tail, digits = 6),
          stage + 1L
        ))
      }
    }

    if (is.null(fit_i)) return(cached)

    result <- tryCatch(
      extract_draws(fit_i, i, model_type = model_type),
      error = function(e) {
        message("  Draw extraction failed for replication ", i, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(result)) return(cached)

    result$n_divergent  <- diagnostics$n_divergent
    result$max_rhat     <- diagnostics$max_rhat
    result$min_ess_bulk <- diagnostics$min_ess_bulk
    result$min_ess_tail <- diagnostics$min_ess_tail
    result$was_refit    <- final_stage > 0L
    result$fit_stage    <- final_stage
    result$converged    <- !needs_refit(diagnostics)
    result$config_id    <- config_id

    saveRDS(result, checkpoint_file)
    result
  }

  results_list <- furrr::future_map(
    seq_len(n_reps),
    runner,
    .options = furrr::furrr_options(seed = TRUE),
    .progress = TRUE
  )

  successful_reps <- sum(!vapply(results_list, is.null, logical(1)))
  failed_reps <- n_reps - successful_reps

  all_parameters <- dplyr::bind_rows(results_list) |>
    dplyr::rename(rep_mean = mean, rep_median = median, rep_sd = sd)
  if (nrow(all_parameters) == 0L) {
    stop("No replications completed successfully for model: ", model_label)
  }

  summary_all   <- summarise_sim(all_parameters, model_label, model_type)
  summary_table <- dplyr::filter(summary_all, type == "primary")
  derived_table <- dplyr::filter(summary_all, type == "derived_Sigma")

  replication_diagnostics <- all_parameters |>
    dplyr::distinct(rep, n_divergent, max_rhat, min_ess_bulk, min_ess_tail,
                    was_refit, fit_stage, converged)

  diagnostic_summary <- replication_diagnostics |>
    dplyr::summarise(
      requested_reps         = n_reps,
      successful_reps        = successful_reps,
      failed_reps            = failed_reps,
      reps_refit             = sum(was_refit, na.rm = TRUE),
      reps_with_divergences  = sum(n_divergent > 0, na.rm = TRUE),
      max_divergences        = safe_max(n_divergent),
      reps_rhat_over_1.01    = sum(max_rhat > 1.01, na.rm = TRUE),
      worst_rhat             = safe_max(max_rhat),
      reps_low_bulk_ess      = sum(min_ess_bulk < 400, na.rm = TRUE),
      reps_low_tail_ess      = sum(min_ess_tail < 400, na.rm = TRUE),
      reps_not_converged     = sum(!converged, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(model = model_label, .before = 1)

  stage_counts <- vapply(
    0:max_stage,
    function(s) sum(replication_diagnostics$fit_stage == s, na.rm = TRUE),
    integer(1)
  )
  names(stage_counts) <- paste0("reps_final_stage_", 0:max_stage)
  diagnostic_summary <- dplyr::bind_cols(diagnostic_summary, tibble::as_tibble_row(stage_counts))

  converged_table <- if (any(all_parameters$converged)) {
    summarise_sim(dplyr::filter(all_parameters, converged), model_label, model_type) |>
      dplyr::filter(type == "primary")
  } else {
    NULL
  }

  list(
    results = summary_table,
    results_converged = converged_table,
    derived = derived_table,
    diagnostics = diagnostic_summary,
    replication_diagnostics = replication_diagnostics
  )
}

n_workers <- min(future::availableCores(), n_reps)
message(sprintf("Running pipeline using %d worker core(s)...", n_workers))
future::plan(future::multisession, workers = n_workers)

correct_exp_out <- run_full_pipeline_parallel(
  model_label = "correct_expDecay",
  checkpoint_dir = checkpoint_dirs[1L],
  bform = bform_correct,
  priors = priors_correct,
  model_type = "correct"
)

wrong_fixed_out <- run_full_pipeline_parallel(
  model_label = "wrong_fixedIE",
  checkpoint_dir = checkpoint_dirs[2L],
  bform = bform_wrong_fixedIE,
  priors = priors_wrong_fixedIE,
  model_type = "wrong_fixed"
)

future::plan(future::sequential)

combined_results <- dplyr::bind_rows(correct_exp_out$results, wrong_fixed_out$results)
combined_results_converged <- dplyr::bind_rows(correct_exp_out$results_converged,
                                               wrong_fixed_out$results_converged)
combined_derived <- dplyr::bind_rows(correct_exp_out$derived, wrong_fixed_out$derived)
combined_diagnostics <- dplyr::bind_rows(correct_exp_out$diagnostics, wrong_fixed_out$diagnostics)
replication_diagnostics <- dplyr::bind_rows(
  dplyr::mutate(correct_exp_out$replication_diagnostics, model = "correct_expDecay", .before = 1),
  dplyr::mutate(wrong_fixed_out$replication_diagnostics, model = "wrong_fixedIE", .before = 1)
)

utils::write.csv(combined_results, "sim_expDecay_summary_results.csv", row.names = FALSE)
utils::write.csv(combined_results_converged, "sim_expDecay_summary_results_converged_only.csv",
                 row.names = FALSE)
utils::write.csv(combined_derived, "sim_expDecay_derived_Sigma.csv", row.names = FALSE)
utils::write.csv(combined_diagnostics, "sim_expDecay_summary_diagnostics.csv", row.names = FALSE)
utils::write.csv(replication_diagnostics, "sim_expDecay_replication_diagnostics.csv", row.names = FALSE)
