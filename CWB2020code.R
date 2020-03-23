suppressPackageStartupMessages(library(tidyverse)) # Required library

## Part 1: Archaeology; Monte Carlo and Bayes ####

# logit: Compute a logit transformation
# Input:
#   x : scalar or vector of values between 0 and 1
# Output:
#   The logit transformation of x

logit <- function(x) {
  log(x) - log(1 - x)
}

# ilogit: Compute the inverse logit transformation
# Input:
#   x : scalar or vector
# Output:
#   The inverse logit transformation of x; values between 0 and 1

ilogit <- function(x) {
  1 / (1 + exp(-x))
}


# arch_param_to_theta: Convert from (N,phi) to (theta1,theta2)
# Input:
#   param : A vector with two elements, N and phi
#   y : data vector
# Output:
#   A length 2 vector with converted parameter values

arch_param_to_theta <- function(param, y) {
  c(
    param[1],
    logit(param[2])
  )
}

# arch_param_from_theta: Convert from (theta1,theta2) to (N,phi)
# Input:
#   theta : A vector with two elements, theta1 and theta2
#   y : data vector
# Output:
#   A length 2 vector with converted parameter values

arch_param_from_theta <- function(theta, y) {
  c(
    theta[1],
    ilogit(theta[2])
  )
}


# arch_loglike: Compute the log-likelihood for multiple Binom(N,phi) observations
# Input:
#   param : Either a vector with two elements, N and phi, or
#           a data.frame with two columns, named N and phi
#   y : data vector
# Output:
#   The log-likelihood for each N,phi combination. The implementation
#   internally uses the log-gamma function to provide a differentiable
#   function of N, which allows optim() to treat N as a continuous variable.
#   If N < max(y), the log-likelihood is defined to be -Inf.

arch_loglike <- function(param, y) {
  if (is.data.frame(param)) {
    ok <- param[["N"]] >= max(y)
    ll <- numeric(nrow(param))
    ll[!ok] <- -Inf
    for (yy in y) {
      ll[ok] <- ll[ok] +
        lchoose(param[["N"]][ok], yy) +
        yy * log(param[["phi"]][ok]) +
        (param[["N"]][ok] - yy) * log(1 - param[["phi"]][ok])
    }
    ll
  } else if (is.matrix(param)) {
    arch_loglike(data.frame(N = param[, 1],
                            phi = param[, 2]),
                 y)
  } else {
    sum(lchoose(param[1], y) +
      y * log(param[2]) +
      (param[1] - y) * log(1 - param[2]))
  }
}



# arch_loglike_from_theta: Compute the log-likelihood for multiple Binom(N,phi) observations
# Input:
#   theta : A vector with two elements, theta1 and theta2
#   y : data vector
# Output:
#   The log-likelihood for the (N,phi) pair obtained by arch_param_from_theta().
#   The implementation internally uses the log-gamma function to provide a
#   differentiable function of N, which allows optim() to treat N as a
#   continuous variable.
#   If N < max(y), the log-likelihood is defined to be -Inf.

arch_loglike_from_theta <- function(theta, y) {
  param <- arch_param_from_theta(theta, y)
  arch_loglike(param, y)
}


# arch_likelihood_estimation: ML estimation for the archaeology model
# Input:
#   y : A vector of observations from Bin(N,phi)
# Output:
#   A data.frame with a single row, and 6 columns:
#     N : The ML estimate of N (treating N as a continuous parameter)
#     N_lower, N_upper : A 95% confidence interval attempt for N
#     phi : The ML estimate of phi
#     phi_lower, phi_upper : A 95% confidence interval attempt for phi

arch_likelihood_estimation <- function(y) {
  opt <- optim(arch_param_to_theta(c(2 * max(y), 1/2), y),
    arch_loglike_from_theta,
    y = y,
    method = "Nelder-Mead",
    control = list(fnscale = -1,
                   maxit = 5000),
    hessian = TRUE
  )
  std_dev <- diag(solve(-opt$hessian))^0.5
  lower <- arch_param_from_theta(opt$par - qnorm(0.975) * std_dev, y)
  upper <- arch_param_from_theta(opt$par - qnorm(0.025) * std_dev, y)
  param <- arch_param_from_theta(opt$par, y)
  data.frame(
    N = param[1],
    N_lower = lower[1],
    N_upper = upper[1],
    phi = param[2],
    phi_lower = lower[2],
    phi_upper = upper[2]
  )
}






## Part 2: 3D printer filament data revisited; ####
# Leave-one-out cross validation and randomisation tests

# score_se: Compute the Squared Error score
# Input:
#   pred : data.frame with (at least) a column "mu"
#   y : data vector
# Output:
#   a vector of scores

score_se <- function(pred, y) {
  (y - pred$mu)^2
}

# score_ds: Compute the Dawid-Sebastiani score
# Input:
#   pred : data.frame with (at least) columns "mu" and "sigma"
#   y : data vector
# Output:
#   a vector of scores

score_ds <- function(pred, y) {
  ((y - pred$mu) / pred$sigma)^2 + 2 * log(pred$sigma)
}

# Input:
#   pred: data.frame with (at least) columns 'lwr' and 'upr'
#   y: data vector
#   alpha: The target probability mass to fall outside the intervals
# Output:
#   A vector of Interval score values

score_interval <- function(pred, y, alpha = 0.1) {
  L <- pred$lwr
  U <- pred$upr
  U - L + 2 / alpha * (pmax(0, L - y) + pmax(0, y - U))
}



# model_Z: Compute model matrices for use by calc_EV
#   Same input/output format as in CWA

model_Z <- function(formulas, data) {
  ZE <- model.matrix(formulas$E, data = data)
  ZV <- model.matrix(formulas$V, data = data)
  list(ZE = cbind(ZE, ZV * 0), ZV = cbind(ZE * 0, ZV))
}

# calc_EV: Evaluate the expectation and variance for model A and B
# Input:
#   theta : A vector with theta1 and theta2
#   data : A data.frame containing the required predictors, including CAD_Weight
#   model : A string denoting which model to use, either "A" or "B".
#   Sigma_theta : If not NULL, an estimate of the covariance matrix for
#                 the uncertainty of estimated theta1 and theta2
# Output:
#   A list with four elements:
#     E : E(y|theta,x)
#     V : Var(y|theta,x)
#     VE : Var(E(y|theta,x)|x) or NULL
#     EV : E(Var(y|theta,x)|x) or NULL

calc_EV <- function(theta, data, model = c("A", "B"),
                    Sigma_theta = NULL) {
  model <- match.arg(model)
  if (model == "A") {
    Z <- model_Z(
      list(
        E = ~ 1 + CAD_Weight,
        V = ~ 1 + CAD_Weight
      ),
      data
    )
    VE <- EV <- NULL
    if (!is.null(Sigma_theta)) {
      # E(Var(y|theta,x)|x)
      EV <- exp(Z$ZV %*% theta + rowSums(Z$ZV * (Z$ZV %*% Sigma_theta)) / 2)
      # Var(E(y|theta,x)|x)
      VE <- rowSums(Z$ZE * (Z$ZE %*% Sigma_theta))
    }
    out <- list(
      E = Z$ZE %*% theta,
      V = exp(Z$ZV %*% theta),
      VE = VE,
      EV = EV
    )
  } else {
    Z <- model_Z(
      list(
        E = ~ 1 + CAD_Weight,
        V = ~ 1 + I(CAD_Weight^2)
      ),
      data
    )
    VE <- EV <- NULL
    if (!is.null(Sigma_theta)) {
      # E(Var(y|theta,x)|x)
      # (pmin: Ignore large Sigma_theta values)
      EV <- Z$ZV %*% exp(theta + pmin(0.5^2, diag(Sigma_theta)) / 2)
      # Var(E(y|theta,x)|x)
      VE <- rowSums(Z$ZE * (Z$ZE %*% Sigma_theta))
    }
    out <- list(
      E = Z$ZE %*% theta,
      V = Z$ZV %*% exp(theta),
      VE = VE,
      EV = EV
    )
  }
  out
}


# neg_log_like: Evaluate the negated log-likelihood for model A and B
# Input:
#   theta : A vector with theta1 and theta2
#   data : A data.frame containing the required predictors, including CAD_Weight
#   model : A string denoting which model to use, either "A" or "B".

neg_log_like <- function(theta, data, model = c("A", "B")) {
  model <- match.arg(model)
  EV <- calc_EV(theta, data, model)
  - sum(dnorm(data[["Actual_Weight"]],
              mean = EV$E,
              sd = EV$V^0.5,
              log = TRUE
  ))
}


# estimate_and_predict: Estimate model and compute predictions
#
# Input:
#   data : A data.frame with data used for parameter estimation
#   newdata : A data.frame with data used for prediction
#   model : A string determining which model to use, either "A" or "B"
#   alpha: The target error probability for the prediction intervals
#   df: Degrees of freedom for t-quantiles for the interval construction.
#       Inf = Normal distribution.
# Output:
#   A data.frame with columns 'mu', 'sigma', 'lwr', 'upr', also see
#   estimate_and_predict_output_template

estimate_and_predict <- function(data, newdata, model = c("A", "B"),
                                 alpha = 0.1, df = Inf) {
  model <- match.arg(model)
  if (model == "A") {
    theta_start <- c(-0.1, 1.07, -2, 0.05)
  } else {
    theta_start <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(theta_start,
    neg_log_like,
    data = data,
    model = model,
    hessian = TRUE,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )
  pred_EV <- calc_EV(opt$par,
    data = newdata, model = model,
    Sigma_theta = solve(opt$hessian)
  )
  pred_sd <- (pred_EV$EV + pred_EV$VE)^0.5

  q <- qt(1 - alpha / 2, df = df)
  lwr <- pred_EV$E - q * pred_sd
  upr <- pred_EV$E + q * pred_sd
  data.frame(mu = pred_EV$E, sigma = pred_sd,
             lwr = lwr, upr = upr)
}

# estimate_and_predict_output_template:
#   Pre-initialise data storage for storing output from estimate_and_predict
# Input:
#   n : The number of output rows to pre-initialise
# Output:
#   A data.frame with columns 'mu', 'sigma', 'lwr', 'upr', with
#   n rows filled with zeros.

estimate_and_predict_output_template <- function(n) {
  data.frame(mu = numeric(n),
             sigma = numeric(n),
             lwr = numeric(n),
             upr = numeric(n))
}


# calc_scores: Evaluate the SE, DS, and Interval scores
# Input:
#   pred : prediction information from estimate_and_predict
#   y : A vector of observations, one for each row of pred
#   alpha : The error probability to use in the Interval score. Default = 0.1

calc_scores <- function(pred, y, alpha = 0.1) {
  data.frame(
    SE = score_se(pred, y),
    DS = score_ds(pred, y),
    Interval = score_interval(pred, y, alpha = alpha)
  )
}
