# Simulation code for 2S-PA Growth Modeling
# Original script: https://github.com/marklhc/awc-growth-supp/blob/main/runsim_awcgrowth.R

# Global Setup ------------------------------------------------------------


library(SimDesign)
library(RPushbullet)

DESIGNFACTOR <- createDesign(
  n = c(100, 250, 1000),
  r_ni = c(0, .25, .55),
  kappa2 = c(0, 0.25), 
  misspecify = c(FALSE, TRUE)
)

# Function to generate covariance matrix for minor factors
get_ucov <- function(p, scale = sqrt(0.1), seed = 123, n = 5) {
  set.seed(seed)
  W <- matrix(rnorm(p * n), nrow = n)
  WtW <- crossprod(W)
  D <- diag(1 / sqrt(diag(WtW))) * scale
  D %*% WtW %*% D
}
# Global covariance matrix for minor factors
UCOV <- get_ucov(20)

# Function for Simulations ------------------------------------------------


generate <- function(condition, fixed_objects = NULL) {
  n_waves <- fixed_objects$n_waves  # number of waves
  n_items <- fixed_objects$n_items  # number of items
  # Loadings and error terms for growth model
  Lambda_growth <- cbind(1, seq_len(n_waves) - 1)
  Theta_growth <- diag(fixed_objects$theta_growth, n_waves)
  # Means and variances for growth factors
  Alpha_growth <- c(fixed_objects$kappa1, condition$kappa2)
  Psi_growth <- fixed_objects$Psi_growth
  if (condition$misspecify) {
    # Add minor quadratic trend
    Lambda_growth <- cbind(Lambda_growth, 
                           (seq_len(n_waves) - (n_waves + 1) / 2)^2)
    Alpha_growth <- c(Alpha_growth, fixed_objects$kappa3)
    Psi_growth <- rbind(
      cbind(Psi_growth, 0), 
      c(0, 0, fixed_objects$phi33)
    )
  }
  # Variances/Covariances and means of latent outcomes
  Psi_cfa <- Lambda_growth %*% Psi_growth %*% t(Lambda_growth) +
    Theta_growth
  Alpha_cfa <- Lambda_growth %*% Alpha_growth
  # Measurement parameters (loadings and intercepts)
  lambda_vec <- fixed_objects$lambda_vec
  lambda_mat <- matrix(
    lambda_vec,
    nrow = n_items, ncol = n_waves
  )
  nu_vec <- fixed_objects$nu_vec
  nu_mat <- matrix(
    nu_vec, 
    nrow = n_items, ncol = n_waves
  )
  # Add misspecification to loadings and intercepts
  r_ni_char <- as.character(condition$r_ni)
  if (condition$r_ni > 0) {
    lambda_mat <- lambda_mat + fixed_objects$lambda_deviation[[r_ni_char]]
    nu_mat <- nu_mat + fixed_objects$nu_deviation[[r_ni_char]]
  }
  # Convert to block loading matrix
  Lambda_cfa <- as.matrix(
    Matrix::bdiag(
      lapply(seq_len(ncol(lambda_mat)),
             function(i, x = lambda_mat) x[ , i])
    )
  )
  Nu_cfa <- as.vector(nu_mat)
  ar_rho <- fixed_objects$ar_rho
  if (condition$misspecify) {
    Theta_cfa <- ar_rho^abs(outer(seq_len(n_waves), seq_len(n_waves), "-")) %x%
      diag(1 - lambda_vec^2 - 0.1)
  } else {
    Theta_cfa <- ar_rho^abs(outer(seq_len(n_waves), seq_len(n_waves), "-")) %x%
      diag(1 - lambda_vec^2)
  }
  # Variances/Covariances and means of observed indicators
  Sigma <- Lambda_cfa %*% Psi_cfa %*% t(Lambda_cfa) + Theta_cfa
  Mu <- Lambda_cfa %*% Alpha_cfa + Nu_cfa
  # Add misspecification due to minor factors
  if (condition$misspecify) {
    Sigma <- Sigma + fixed_objects$ucov
  }
  dat <- mvnfast::rmvn(condition$n, mu = Mu, sigma = Sigma)
  colnames(dat) <- paste0("y",
                          outer(seq_len(n_items), seq_len(n_waves), paste0))
  as.data.frame(dat)
}
# Test:
# dat <- generate(DESIGNFACTOR[1, ], fixed_objects = list(
#   ucov = UCOV,
#   n_waves = 4,
#   n_items = 5,
#   theta_growth = .50,
#   kappa1 = 0,
#   Psi_growth = matrix(c(.5, .089, .089, .1), nrow = 2),
#   kappa3 = -0.01,
#   phi33 = 0.004,
#   lambda_vec = c(.8, .5, .7, .65, .7),
#   nu_vec = c(0, .5, -.25, .25, -.5),
#   lambda_deviation =
#     list(`0.25` = rbind(0, 0, 0, 0, seq(0, 0.3, length.out = 4)),
#          `0.55` = rbind(0, 0,
#                         c(0, 0.2, 0.3, 0.1),
#                         c(0, -0.05, 0, 0.05),
#                         seq(0, 0.3, length.out = 4))),
#   nu_deviation =
#     list(`0.25` = rbind(0, 0, 0,
#                         seq(0, 0.75, length.out = 4),
#                         c(0, 0.125, -0.125, 0)),
#          `0.55` = rbind(c(0, 0.75, 0.5, 0.25),
#                         c(0, 0.25, 0, 0.50),
#                         0,
#                         seq(0, 0.75, length.out = 4),
#                         c(0, 0.125, -0.125, 0))),
#   ar_rho = .20  # autoregressive coefficient
# ))

fit_config <- function(dat, n_waves, n_items) {
  ynames <- outer(seq_len(n_items), seq_len(n_waves), paste0)
  load_lines <- apply(ynames, 2,
                      function(x) paste0("y", x, collapse = " + "))
  cfa_model <- paste0(
    c(paste0("eta", seq_len(n_waves), " =~ ", load_lines,
             collapse = "\n  "),
      paste0(
        apply(ynames, 1, function(x) {
          combn(paste0("y", x), 2, FUN = paste0, collapse = " ~~ ")
        }),
        collapse = "\n  "
      )
    ),
    collapse = "\n  "
  )
  lavaan::cfa(cfa_model, data = dat, std.lv = TRUE,
              meanstructure = TRUE)
}
# Test:
# fit_config(dat, n_waves = 4, n_items = 5)

align_lambda_nu <- function(config_pars, ...) {
  config_lambda <- config_pars$lambda
  n_waves <- ncol(config_lambda)
  n_items <- nrow(config_lambda) %/% n_waves
  config_lambda <- crossprod(config_lambda,
                             rep(1, n_waves) %x% diag(n_items))
  config_nu <- matrix(config_pars$nu,
                      nrow = n_waves,
                      ncol = n_items,
                      byrow = TRUE)
  # Assume equal sample size
  aligned_pars <- sirt::invariance.alignment(config_lambda, config_nu,
                                             ...)
  aligned_pars[c("lambda.aligned", "nu.aligned", "pars")]
}

fit_2nd_growth_w_labels <- function(dat, n_waves, n_items,
                                    lambda_prefix = NULL,
                                    nu_prefix = NULL) {
  yindices <- outer(seq_len(n_items), seq_len(n_waves), paste0)
  ynames <- t(yindices)
  ynames[] <- paste0("y", ynames)
  load_lines <- apply(`[<-`(ynames,
                            paste0(lambda_prefix, " * ", ynames)),
                      1,
                      paste0, collapse = " + ")
  growth_model <- paste0(
    c(paste0("eta", seq_len(n_waves), " =~ .8 * ", ynames[ , 1], " + ",
             load_lines,
             collapse = "\n  "),
      paste0("i =~ ", paste0("1 * eta", seq_len(n_waves),
                             collapse = " + ")),
      paste0("s =~ ", paste0(seq_len(n_waves) - 1, " * eta",
                             seq_len(n_waves), collapse = " + ")),
      paste0(ynames[ , 1], " ~ 0 * 1 + ",
             nu_prefix[ , 1], " * 1", collapse = "\n  "),
      paste0(ynames[ , -1], " ~ ",
             nu_prefix[ , -1],
             " * 1", collapse = "\n  "),
      paste0(
        apply(ynames, 2, function(x) {
          combn(x, 2, FUN = paste0, collapse = " ~~ ")
        }),
        collapse = "\n  "
      ),
      paste0("eta", seq_len(n_waves),
             " ~~ psi * eta", seq_len(n_waves), collapse = "\n  "),
      paste0("eta", seq_len(n_waves), " ~ 0 * 1", collapse = "\n  "),
      "# Sum of loadings and intercepts",
      paste0("sum_load := .8 + ",
             paste0(lambda_prefix[1, -1], collapse = " + ")),
      paste0("sum_int := ",
             paste0(nu_prefix[1, -1], collapse = " + "))
    ),
    collapse = "\n  "
  )
  lavaan::growth(growth_model, data = dat)
}

fit_2nd_growth_w_starts <- function(dat, n_waves, n_items,
                                    lambda_start = NULL,
                                    nu_start = NULL) {
  yindices <- outer(seq_len(n_items), seq_len(n_waves), paste0)
  ynames <- t(yindices)
  ynames[] <- paste0("y", ynames)
  lambda_labels <- t(yindices)
  lambda_labels[] <- paste0("l", lambda_labels)
  nu_labels <- t(yindices)
  nu_labels[] <- paste0("n", nu_labels)
  lambda_start_org <- lambda_start
  lambda_start[] <- paste0(lambda_labels, " * ", ynames,
                           " + start(", lambda_start, ")")
  nu_start_org <- nu_start
  nu_start[] <- paste0(nu_labels, " * 1 + start(", nu_start, ")")
  load_lines <- apply(`[<-`(ynames,
                            paste0(lambda_start, " * ", ynames)),
                      1,
                      paste0, collapse = " + ")
  growth_model <- paste0(
    c(paste0("eta", seq_len(n_waves), " =~ ",
             load_lines,
             collapse = "\n  "),
      paste0("i =~ ", paste0("1 * eta", seq_len(n_waves),
                             collapse = " + ")),
      paste0("s =~ ", paste0(seq_len(n_waves) - 1, " * eta",
                             seq_len(n_waves), collapse = " + ")),
      paste0(ynames[ , 1], " ~ ", nu_start_org[ , 1],
             " * 1", collapse = "\n  "),
      paste0(ynames[ , -1], " ~ ",
             nu_start[ , -1],
             " * 1", collapse = "\n  "),
      paste0(
        apply(ynames, 2, function(x) {
          combn(x, 2, FUN = paste0, collapse = " ~~ ")
        }),
        collapse = "\n  "
      ),
      paste0("eta", seq_len(n_waves),
             " ~~ psi * eta", seq_len(n_waves), collapse = "\n  "),
      paste0("eta", seq_len(n_waves), " ~ 0 * 1", collapse = "\n  "),
      "# Sum of loadings and intercepts",
      paste0("sum_load := ", lambda_start_org[1, 1], " + ",
             paste0(lambda_labels[1, -1], collapse = " + ")),
      paste0("sum_int := ", nu_start_org[1, 1], " + ",
             paste0(nu_labels[1, -1], collapse = " + "))
    ),
    collapse = "\n  "
  )
  lavaan::growth(growth_model, data = dat)
}

# Get factor scores
get_fs <- function(config_fit, aligned_pars, n_waves, n_items) {
  pars_cf <- lavInspect(config_fit, what = "est")
  ey <- config_fit@implied$mean[[1]]
  yc <- sweep(config_fit@Data@X[[1]], MARGIN = 2, STATS = ey)
  lam <- Matrix::bdiag(split(aligned_pars$lambda.aligned, 
                             rep(seq_len(n_waves), n_items)))
  lam <- as.matrix(lam)
  th <- pars_cf$theta
  thinv_lam <- solve(th, lam)
  fs_mat <- solve(crossprod(thinv_lam, lam), t(thinv_lam))
  alpha <- aligned_pars$pars$alpha0
  fs <- t(fs_mat %*% t(yc) + c(alpha))
  # error variances of factor scores
  ev <- solve(crossprod(thinv_lam, lam))
  colnames(fs) <- paste0("fs_", paste0("eta", seq_len(n_waves)))
  list(fs = fs, ev = ev)
}

fit_2spa_growth <- function(fs, ev_fs, n_waves) {
  growth_model <- paste0(c(
    paste0("eta", seq_len(n_waves), " =~ 1 * fs_", "eta", seq_len(n_waves)), 
    "# Constrain error variances", 
    paste0("fs_eta", seq_len(n_waves), " ~~ ", 
           diag(ev_fs), " * fs_eta", seq_len(n_waves)), 
    "# Constrain error covariances", 
    paste0("fs_eta", col(lower.tri(ev_fs))[lower.tri(ev_fs)], 
           " ~~ ", ev_fs[lower.tri(ev_fs)], 
           " * fs_eta", row(lower.tri(ev_fs))[lower.tri(ev_fs)]), 
    "# Constrain intercepts", 
    paste0("fs_eta", seq_len(n_waves), " ~ 0 * 1"), 
    "# Linear Growth Model", 
    paste0("i =~ ", paste0((paste0("1 * eta", seq_len(n_waves))), 
                           collapse = " + ")), 
    paste0("s =~ ", paste0(paste0(seq_len(n_waves) - 1, " * eta", seq_len(n_waves)), 
                           collapse = " + ")), 
    "# Variance-covariances of intercepts and slopes", 
    "i ~~ i", "s ~~ s", "i ~~ s", 
    "# Means of level and slope", 
    "i ~ 1", "s ~ 1", 
    "# Fixed disturbances of latent outcomes to zero", 
    paste0("eta", seq_len(n_waves), " ~ 0 * 1")
  ), collapse = "\n")
  lavaan::sem(growth_model, data = fs)
}

# Helper function to extract output
extract_res <- function(object, pars = c("i~1", "s~1", "i~~i", "s~~s"),
                        true_load = .67 * 5, true_int = 0, tspa = FALSE) {
  coefs <- lavaan::coef(object, type = "user")
  scaled_ses <- sqrt(diag(lavaan::vcov(object)[pars, pars]))
  scaled_ests <- coefs[pars]
  if (!tspa) {
    sum_load <- coefs["sum_load"]
    sum_int <- coefs["sum_int"]
    scaled_ests[1] <- (sum_int - true_int + sum_load * scaled_ests[1]) /
      true_load
    scaled_ests[2] <- scaled_ests[2] * sum_load / true_load
    scaled_ests[3:4] <- (sum_load / true_load)^2 * scaled_ests[3:4]
    scaled_ses[1:2] <- (sum_load / true_load) * scaled_ses[1:2]
    scaled_ses[3:4] <- (sum_load / true_load) * scaled_ses[3:4]
  }
  c(scaled_ests, scaled_ses)
}

analyse <- function(condition, dat, fixed_objects = NULL) {
  n_waves <- fixed_objects$n_waves
  n_items <- fixed_objects$n_items
  # Initialize output
  methods <- c("TSPA")
  pars <- c("meani", "means", "vari", "vars")
  stats <- c("est", "ase")
  out <- vector("list", length(methods))
  names(out) <- methods
  # Configural model
  config_fit <- fit_config(dat, n_waves = n_waves, n_items = n_items)
  # Alignment
  config_pars <- lavaan::lavInspect(config_fit, what = "est")
  aligned_pars <- align_lambda_nu(config_pars, fixed = TRUE)
  # Growth model (2SPA)
  fs_ev <- get_fs(config_fit, aligned_pars, n_waves, n_items)
  growth_fit_2spa <- fit_2spa_growth(fs_ev$fs, fs_ev$ev, n_waves)
  out$TSPA <- extract_res(growth_fit_2spa, tspa = TRUE)
  # Reorder output
  out <- unlist(out)[c(aperm(array(seq_len(length(methods) * length(stats) *
                                             length(pars)),
                                   dim = c(length(pars),
                                           length(stats),
                                           length(methods))), c(3, 1, 2)))]
  # Add names
  names_out <- outer(
    outer(stats, methods, FUN = paste, sep = "_"), pars,
    FUN = paste, sep = "_")
  names(out) <- aperm(names_out, perm = c(2, 3, 1))
  out
}
# Test:
# analyse(DESIGNFACTOR[1, ], dat,
#         fixed_objects = list(n_waves = 4, n_items = 5))

evaluate <- function(condition, results, fixed_objects = NULL) {
  # Specify population value
  methods <- c("TSPA")
  pars <- c("meani", "means", "vari", "vars")
  pop <- rep(c(fixed_objects$kappa1, condition$kappa2, 
               diag(fixed_objects$Psi_growth)), each = length(methods))
  results_est <- results[ , grep("^est_", colnames(results))]
  results_ase <- results[ , grep("^ase_", colnames(results))]
  results_ci <- matrix(
    as.matrix(
      rbind(results_est - qnorm(.975) * results_ase,
            results_est + qnorm(.975) * results_ase)
    ),
    nrow = nrow(results)
  )
  ret_ci <- SimDesign::ECR(
    sweep(results_ci, 2, rep(pop, each = 2), "-"),
    parameter = 0
  )
  names(ret_ci) <- gsub("est_", "ci_", colnames(results_est))
  ret <- c(
    bias = SimDesign::bias(
      results_est,
      parameter = pop),
    std_bias = SimDesign::bias(
      results_est,
      parameter = pop,
      type = "standardized"),
    rmse = SimDesign::RMSE(
      results_est,
      parameter = pop
    ),
    coverage = ret_ci
  )
  ret
}


#-------------------------------------------------------------------

res <-
  runSimulation(
    design = DESIGNFACTOR,
    replications = 2500,
    generate = generate,
    analyse = analyse,
    summarise = evaluate,
    packages = c("Matrix", "lavaan", "mvnfast", "sirt"),
    fixed_objects = list(
      ucov = UCOV,
      n_waves = 4,
      n_items = 5,
      theta_growth = .50,
      kappa1 = 0,
      Psi_growth = matrix(c(.5, .089, .089, .1), nrow = 2),
      kappa3 = -0.01,
      phi33 = 0.004,
      lambda_vec = c(.8, .5, .7, .65, .7),
      nu_vec = c(0, .5, -.25, .25, -.5),
      lambda_deviation =
        list(`0.25` = rbind(0, 0, 0, 0, seq(0, 0.3, length.out = 4)),
             `0.55` = rbind(0, 0,
                            c(0, 0.2, 0.3, 0.1),
                            c(0, -0.05, 0, 0.05),
                            seq(0, 0.3, length.out = 4))),
      nu_deviation =
        list(`0.25` = rbind(0, 0, 0,
                            seq(0, 0.75, length.out = 4),
                            c(0, 0.125, -0.125, 0)),
             `0.55` = rbind(c(0, 0.75, 0.5, 0.25),
                            c(0, 0.25, 0, 0.50),
                            0,
                            seq(0, 0.75, length.out = 4),
                            c(0, 0.125, -0.125, 0))),
      ar_rho = .20  # autoregressive coefficient
    ), 
    seed = rep(670084, nrow(DESIGNFACTOR)),
    save = TRUE,
    save_results = TRUE,
    filename = "simulation/simres_2spa_only",
    save_details = list(
      save_results_dirname = "simulation/simdetails_2spa_only"
    ),
    # allow_na = TRUE,
    parallel = TRUE,
    ncores = min(4L, parallel::detectCores() - 1), 
    notification = "condition"
  )
