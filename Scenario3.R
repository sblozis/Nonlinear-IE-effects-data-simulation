# ============================================================
# Simulation 3: Nonlinear Exponential-Decay MELS Model
# IE effect at time = 1 is random across subjects
# decay rate parameter is fixed across subjects
# Correct model: random asymptote, random initial-elevation effect,
#                and random residual scale
# Misspecified model: fixed initial-elevation effect
# ============================================================

RESET_CACHE <- FALSE

checkpoint_dirs <- c(
  "sim2_checkpoints_correct_logsd",
  "sim2_checkpoints_wrong_logsd"
)

if (RESET_CACHE) {
  unlink(checkpoint_dirs, recursive = TRUE)
}

tmp_root <- file.path(tempdir(), "sim2_cmdstan_tmp")
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
#   log_sigma0    = fixed intercept of the residual SD on the log-SD scale,
#                   i.e., log(sigma_i) = log_sigma0 + u_sigma
#   sd_int        = SD of the random intercept (asymptote) u_b0
#   sd_ie         = SD of the random initial-elevation (IE) effect u_b1
#   sd_scale      = SD of the random scale u_sigma on the log-SD scale
#   cor_int_ie    = Corr(u_b0, u_b1)
#   cor_int_scale = Corr(u_b0, u_sigma)
#   cor_ie_scale  = Corr(u_b1, u_sigma)
true_vals2 <- c(
  b0            = 2.5,
  b1            = 0.5,
  lambda        = 1.15,
  log_sigma0    = -0.75,
  sd_int        = sqrt(0.7),
  sd_ie         = sqrt(0.25),
  sd_scale      = sqrt(0.125),
  cor_int_ie    = 0.15,
  cor_int_scale = -0.20,
  cor_ie_scale  = -0.10
)

re_terms <- c("int", "ie", "scale")   

build_R2 <- function(tv) {
  R <- diag(3L)
  dimnames(R) <- list(re_terms, re_terms)
  R["int", "ie"]    <- R["ie", "int"]    <- tv[["cor_int_ie"]]
  R["int", "scale"] <- R["scale", "int"] <- tv[["cor_int_scale"]]
  R["ie", "scale"]  <- R["scale", "ie"]  <- tv[["cor_ie_scale"]]
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

build_Sigma2 <- function(tv) {
  sds <- c(tv[["sd_int"]], tv[["sd_ie"]], tv[["sd_scale"]])
  if (any(!is.finite(sds)) || any(sds <= 0)) {
    stop("All random-effect SDs must be finite and strictly positive.")
  }
  R <- build_R2(tv)
  check_corr_matrix(R)
  D <- diag(sds)
  dimnames(D) <- list(re_terms, re_terms)
  D %*% R %*% D
}

local({
  Sigma <- build_Sigma2(true_vals2)
  sds   <- c(true_vals2[["sd_int"]], true_vals2[["sd_ie"]], true_vals2[["sd_scale"]])
  stopifnot(isTRUE(all.equal(unname(sqrt(diag(Sigma))), sds, tolerance = 1e-12)))
  stopifnot(isTRUE(all.equal(unname(stats::cov2cor(Sigma)),
                             unname(build_R2(true_vals2)), tolerance = 1e-12)))
  R_old <- matrix(c( 1.00,  0.15, -0.20,
                     0.15,  1.00, -0.10,
                    -0.20, -0.10,  1.00), nrow = 3L, byrow = TRUE)
  D_logvar     <- diag(c(sqrt(0.7), sqrt(0.25), sqrt(0.5)))  
  Sigma_logvar <- D_logvar %*% R_old %*% D_logvar
  T_half       <- diag(c(1, 1, 1 / 2))                        
  stopifnot(isTRUE(all.equal(unname(T_half %*% Sigma_logvar %*% T_half),
                             unname(Sigma), tolerance = 1e-12)))
  stopifnot(isTRUE(all.equal(true_vals2[["log_sigma0"]], -1.5 / 2)))
})

true_Sigma_vals2 <- local({
  S    <- build_Sigma2(true_vals2)
  keep <- which(upper.tri(S, diag = TRUE), arr.ind = TRUE)
  stats::setNames(S[keep], paste0("Sigma_", re_terms[keep[, 1]], "_", re_terms[keep[, 2]]))
})

n_clusters2    <- 300L
n_per_cluster2 <- 10L

simulate_rep2 <- function(rep_id, n_clusters, n_per_cluster, tv,
                           seed_offset = 0L) {
  set.seed(123 + seed_offset + rep_id)

  Sigma2 <- build_Sigma2(tv)

  u <- MASS::mvrnorm(
    n = n_clusters,
    mu = c(0, 0, 0),
    Sigma = Sigma2
  )

  cluster <- rep(seq_len(n_clusters), each = n_per_cluster)
  time    <- rep(0:(n_per_cluster - 1L), times = n_clusters)

  b0_i    <- tv[["b0"]] + u[cluster, 1L]             
  b1_i    <- tv[["b1"]] + u[cluster, 2L]                 
  sigma_i <- exp(tv[["log_sigma0"]] + u[cluster, 3L])    
  mu_i    <- b0_i + b1_i * exp(-tv[["lambda"]] * time)

  y <- stats::rnorm(
    n_clusters * n_per_cluster,
    mean = mu_i,
    sd = sigma_i
  )

  data.frame(y = y, cluster = cluster, time = time)
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
  draw_names <- posterior::variables(draws)
  population_names <- draw_names[
    grepl("^(b_|sd_|cor_)", draw_names) & !grepl("^r_", draw_names)
  ]

  if (length(population_names) == 0L) {
    return(list(
      n_divergent = n_divergent,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  summaries <- posterior::summarise_draws(
    posterior::subset_draws(draws, variable = population_names),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  low_parameters <- summaries |>
    dplyr::filter(ess_bulk < 400 | ess_tail < 400) |>
    dplyr::pull(variable)

  if (length(low_parameters) > 0L) {
    message(
      "  [Warning] ESS below 400 for: ",
      paste(low_parameters, collapse = ", ")
    )
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
    stop(
      "Could not uniquely identify the brms correlation parameter for ",
      term1,
      " and ",
      term2,
      ". Available correlation parameters: ",
      paste(grep("^cor_", names(draws), value = TRUE), collapse = ", ")
    )
  }

  draws[[found]]
}

primary_parameters <- names(true_vals2)
model_terms <- function(random_b1) if (random_b1) re_terms else c("int", "scale")
brms_term <- c(int = "b0_Intercept", ie = "b1_Intercept", scale = "sigma_Intercept")
cor_name  <- function(pair) paste0("cor_", pair[[1]], "_", pair[[2]])

expected_variables <- function(random_b1) {
  terms <- model_terms(random_b1)
  keep  <- which(upper.tri(diag(length(terms)), diag = TRUE), arr.ind = TRUE)
  c(intersect(primary_parameters,
              c("b0", "b1", "lambda", "log_sigma0", paste0("sd_", terms),
                vapply(utils::combn(terms, 2L, simplify = FALSE), cor_name, ""))),
    paste0("Sigma_", terms[keep[, 1]], "_", terms[keep[, 2]]))
}

summarise_draw_vector <- function(values, rep_id, parameter, type) {
  tibble::tibble(rep = rep_id, type = type, variable = parameter,
                 mean = mean(values), median = stats::median(values),
                 sd = stats::sd(values),
                 lower_2.5  = unname(stats::quantile(values, 0.025)),
                 upper_97.5 = unname(stats::quantile(values, 0.975)))
}

extract_draws_sim2 <- function(fit, rep_id, random_b1 = TRUE) {
  draws <- as.data.frame(posterior::as_draws_df(
    fit, variable = "^(b_|sd_cluster__|cor_cluster__)", regex = TRUE))
  terms      <- model_terms(random_b1)
  term_pairs <- utils::combn(terms, 2L, simplify = FALSE)

  required <- c("b_b0_Intercept", "b_b1_Intercept", "b_lambda_Intercept",
                "b_sigma_Intercept", paste0("sd_cluster__", brms_term[terms]))
  missing <- setdiff(required, names(draws))
  if (length(missing) > 0L) {
    stop("Required posterior draw(s) not found: ", paste(missing, collapse = ", "))
  }

  draws$b0         <- draws$b_b0_Intercept
  draws$b1         <- draws$b_b1_Intercept
  draws$lambda     <- draws$b_lambda_Intercept
  draws$log_sigma0 <- draws$b_sigma_Intercept
  for (term in terms) {
    draws[[paste0("sd_", term)]] <- draws[[paste0("sd_cluster__", brms_term[[term]])]]
  }
  for (pair in term_pairs) {
    draws[[cor_name(pair)]] <- get_correlation_draw(
      draws, brms_term[[pair[[1]]]], brms_term[[pair[[2]]]])
  }
  parameters <- intersect(primary_parameters, expected_variables(random_b1))

  k        <- length(terms)
  sd_mat   <- as.matrix(draws[paste0("sd_", terms)])
  cor_mat  <- as.matrix(draws[vapply(term_pairs, cor_name, "")])
  pair_idx <- utils::combn(k, 2L)
  keep     <- which(upper.tri(diag(k), diag = TRUE), arr.ind = TRUE)
  sigma_names <- paste0("Sigma_", terms[keep[, 1]], "_", terms[keep[, 2]])
  sigma_draws <- t(vapply(seq_len(nrow(draws)), function(d) {
    R <- diag(k)
    R[t(pair_idx)]        <- cor_mat[d, ]
    R[t(pair_idx[2:1, ])] <- cor_mat[d, ]
    D <- diag(sd_mat[d, ], nrow = k)
    (D %*% R %*% D)[keep]
  }, numeric(nrow(keep))))
  colnames(sigma_draws) <- sigma_names

  dplyr::bind_rows(
    purrr::map_dfr(parameters, function(p) summarise_draw_vector(draws[[p]], rep_id, p, "primary")),
    purrr::map_dfr(sigma_names, function(p) summarise_draw_vector(sigma_draws[, p], rep_id, p, "derived_Sigma"))
  )
}

summarise_sim2 <- function(all_parameters, model_label) {
  truth <- dplyr::bind_rows(
    tibble::tibble(type = "primary", variable = names(true_vals2), true_value = unname(true_vals2)),
    tibble::tibble(type = "derived_Sigma", variable = names(true_Sigma_vals2),
                   true_value = unname(true_Sigma_vals2)))
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
    dplyr::mutate(model = model_label, n_reps = dplyr::coalesce(n_reps, 0L),
                  status = dplyr::if_else(n_reps > 0L, "estimated",
                                          "not estimated: no random IE effect in this model")) |>
    dplyr::select(model, type, variable, status, true_value, n_reps, avg_est, bias, rmse,
                  emp_sd, avg_post_sd, coverage_95, avg_median, bias_median, rmse_median,
                  mcse_bias, mcse_coverage)
}

bform2_correct <- brms::bf(
  y ~ b0 + b1 * exp(-lambda * time),
  b0 ~ 1 + (1 |q| cluster),
  b1 ~ 1 + (1 |q| cluster),
  lambda ~ 1,
  sigma ~ 1 + (1 |q| cluster),
  nl = TRUE
)

priors2_correct <- c(
  prior(normal(0, 5), nlpar = "b0"),
  prior(normal(0, 2), nlpar = "b1"),
  prior(normal(1, 1), nlpar = "lambda", lb = 0),
  prior(normal(0, 1.5), class = "Intercept", dpar = "sigma"),
  prior(student_t(3, 0, 1), class = "sd", nlpar = "b0"),
  prior(student_t(3, 0, 1), class = "sd", nlpar = "b1"),
  prior(student_t(3, 0, 1), class = "sd", dpar = "sigma"),
  prior(lkj(1), class = "cor")
)

bform2_wrong <- brms::bf(
  y ~ b0 + b1 * exp(-lambda * time),
  b0 ~ 1 + (1 |q| cluster),
  b1 ~ 1,
  lambda ~ 1,
  sigma ~ 1 + (1 |q| cluster),
  nl = TRUE
)

priors2_wrong <- c(
  prior(normal(0, 5), nlpar = "b0"),
  prior(normal(0, 2), nlpar = "b1"),
  prior(normal(1, 1), nlpar = "lambda", lb = 0),
  prior(normal(0, 1.5), class = "Intercept", dpar = "sigma"),
  prior(student_t(3, 0, 1), class = "sd", nlpar = "b0"),
  prior(student_t(3, 0, 1), class = "sd", dpar = "sigma"),
  prior(lkj(1), class = "cor")
)

run_full_pipeline2_parallel <- function(model_label, checkpoint_dir,
                                         bform, priors, random_b1) {
  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  
    config_id <- rlang::hash(list(
    true_vals = true_vals2, design = c(n_clusters2, n_per_cluster2),
    mcmc = c(n_chains, n_warmup, n_iter, n_warmup_refit, n_iter_refit, divergence_refit_threshold),
    formula = utils::capture.output(print(bform)), priors = as.data.frame(priors),
    output = "sd-cor-median-Sigma-v2"))
  expected_vars <- expected_variables(random_b1)

  message("=== ", model_label, ": compiling Stan model once ===")
  compilation_data <- simulate_rep2(
    rep_id = 1L,
    n_clusters = n_clusters2,
    n_per_cluster = n_per_cluster2,
    tv = true_vals2
  )

  compiled_fit <- brms::brm(
    formula = bform,
    data = compilation_data,
    prior = priors,
    chains = 1L,
    iter = 20L,
    warmup = 10L,
    cores = 1L,
    backend = "cmdstanr",
    seed = 123,
    refresh = 0,
    init = 0,
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

    dat <- simulate_rep2(
      rep_id = i, n_clusters = n_clusters2,
      n_per_cluster = n_per_cluster2, tv = true_vals2
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
      extract_draws_sim2(fit_i, i, random_b1 = random_b1),
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
    summary_all   <- summarise_sim2(all_parameters, model_label)
  summary_table <- dplyr::filter(summary_all, type == "primary")
  derived_table <- dplyr::filter(summary_all, type == "derived_Sigma")

  replication_diagnostics <- all_parameters |>
    dplyr::distinct(
      rep,
      n_divergent,
      max_rhat,
      min_ess_bulk,
      min_ess_tail,
      was_refit,
      fit_stage,
      converged
    )

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
    summarise_sim2(dplyr::filter(all_parameters, converged), model_label) |>
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
future::plan(future::multisession, workers = n_workers)

correct_exp_out <- run_full_pipeline2_parallel(
  model_label = "correct_expDecay",
  checkpoint_dir = checkpoint_dirs[1L],
  bform = bform2_correct,
  priors = priors2_correct,
  random_b1 = TRUE
)

wrong_exp_out <- run_full_pipeline2_parallel(
  model_label = "misspecified_fixed_effect",
  checkpoint_dir = checkpoint_dirs[2L],
  bform = bform2_wrong,
  priors = priors2_wrong,
  random_b1 = FALSE
)

future::plan(future::sequential)

combined_results <- dplyr::bind_rows(
  correct_exp_out$results,
  wrong_exp_out$results
)

combined_results_converged <- dplyr::bind_rows(
  correct_exp_out$results_converged,
  wrong_exp_out$results_converged
)

combined_derived <- dplyr::bind_rows(
  correct_exp_out$derived,
  wrong_exp_out$derived
)

combined_diagnostics <- dplyr::bind_rows(
  correct_exp_out$diagnostics,
  wrong_exp_out$diagnostics
)

replication_diagnostics <- dplyr::bind_rows(
  dplyr::mutate(
    correct_exp_out$replication_diagnostics,
    model = "correct_expDecay",
    .before = 1
  ),
  dplyr::mutate(
    wrong_exp_out$replication_diagnostics,
    model = "misspecified_fixed_effect",
    .before = 1
  )
)

utils::write.csv(
  combined_results,
  "sim2_expDecay_summary_results.csv",
  row.names = FALSE
)
utils::write.csv(
  combined_results_converged,
  "sim2_expDecay_summary_results_converged_only.csv",
  row.names = FALSE
)
utils::write.csv(
  combined_derived,
  "sim2_expDecay_derived_Sigma.csv",
  row.names = FALSE
)
utils::write.csv(
  combined_diagnostics,
  "sim2_expDecay_summary_diagnostics.csv",
  row.names = FALSE
)
utils::write.csv(
  replication_diagnostics,
  "sim2_expDecay_replication_diagnostics.csv",
  row.names = FALSE
)
