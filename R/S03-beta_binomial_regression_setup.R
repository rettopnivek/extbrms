# Functions for beta-binomial custom family
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2021-10-01

# Table of contents
# 1) beta_binomial_regression_setup

#### 1) beta_binomial_regression_setup ####
#' Title
#'
#' Description.
#'
#' @param global Logical; if \code{TRUE} uses
#'   the \code{<<-} assignment to add all
#'   functions and terms to the global environment.
#'
#' @details Forthcoming.
#'
#' @return Output.
#'
#' @examples
#' # Examples
#'
#' @export

beta_binomial_regression_setup <- function( global = FALSE ) {

  # Create custom family for beta-binomial
  # distribution with mean/precision parameterization
  beta_binomial_m_prec <- custom_family(
    "beta_binomial_m_prec", dpars = c("mu", "phi"),
    links = c("logit", "log"), lb = c(NA, 0),
    type = "int", vars = "vint1[n]"
  )

  # Define the corresponding Stan density function
  beta_binomial_m_prec_stan_functions <- "
    real beta_binomial_m_prec_lpmf(int y, real mu, real phi, int N) {
      return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
    }
    int beta_binomial_m_prec_rng(real mu, real phi, int N) {
      return beta_binomial_rng(N, mu * phi, (1 - mu) * phi);
    }
  "

  #
  beta_binomial_m_prec_stanvar <- stanvar(
    scode = beta_binomial_m_prec_stan_functions,
    block = "functions"
  )

  beta_binomial_m_prec_lpmf <- function( y, mu, phi, N ) {

    a.dnm <- mu * phi
    b.dnm <- (1-mu)*phi

    a.num <- y + a.dnm
    b.num <- N - y + b.dnm

    p1 = beta( a.num, b.num ) / beta( a.dnm, b.dnm )
    p2 = choose( N, y )

    return( log( p2 * p1 ) )
  }

  beta_binomial_m_prec_rng <- function( mu, phi, N ) {

    a.beta <- mu * phi
    b.beta <- (1-mu)*phi

    nobs <- max( length( a.beta ), length( b.beta ), length( N ) )

    prob <- rbeta( nobs, a.beta, b.beta )
    out <- rbinom( nobs, N, prob )

    return( out )
  }

  log_lik_beta_binomial_m_prec <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    phi <- brms::get_dpar(prep, "phi", i = i)
    trials <- prep$data$vint1[i]
    y <- prep$data$Y[i]
    beta_binomial_m_prec_lpmf(y, mu, phi, trials)
  }

  posterior_predict_beta_binomial_m_prec <- function(i, prep, ...) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    phi <- brms::get_dpar(prep, "phi", i = i)
    trials <- prep$data$vint1[i]
    beta_binomial_m_prec_rng(mu, phi, trials)
  }

  if ( global ) {

    beta_binomial_m_prec <<- beta_binomial_m_prec
    beta_binomial_m_prec_stan_functions <<-
      beta_binomial_m_prec_stan_functions
    beta_binomial_m_prec_stanvar <<-
      beta_binomial_m_prec_stanvar
    beta_binomial_m_prec_lpmf <<-
      beta_binomial_m_prec_lpmf
    beta_binomial_m_prec_rng <<-
      beta_binomial_m_prec_rng

    log_lik_beta_binomial_m_prec <<-
      log_lik_beta_binomial_m_prec

    posterior_predict_beta_binomial_m_prec <<-
      posterior_predict_beta_binomial_m_prec

  } else {

    out <- list(
      beta_binomial_m_prec = beta_binomial_m_prec,
      beta_binomial_m_prec_stan_functions =
        beta_binomial_m_prec_stan_functions,
      beta_binomial_m_prec_stanvar =
        beta_binomial_m_prec_stanvar,
      beta_binomial_m_prec_lpmf =
        beta_binomial_m_prec_lpmf,
      beta_binomial_m_prec_rng =
        beta_binomial_m_prec_rng,
      log_lik_beta_binomial_m_prec =
        log_lik_beta_binomial_m_prec,
      posterior_predict_beta_binomial_m_prec =
        posterior_predict_beta_binomial_m_prec
    )

    return( out )
  }

}
