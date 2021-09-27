#' Specify Algorithm Settings for 'brm'
#'
#' A convenience function used to specify settings
#' when sampling from the joint posterior of a model
#' fit using \code{\link[brms]{brm}}. Can be used to
#' set the number of samples from the posterior to draw,
#' the warm-up period used pre-sampling, the number
#' of chains and cores, and additional control options.
#'
#' @param desired_iter The total number of samples to draw
#'   from the posterior (the sum of individual chain).
#' @param rng_seed RNG seed used to ensure reproducibility.
#' @param warmup The number of iterations to use for the
#'   warm-up period before drawing samples from the posterior.
#' @param adapt_delta Proportion governing the gradient descent
#'   algorithm; higher values result in slower but more
#'   stable estimates (i.e., lower risk of divergent estimates
#'   where gradient descent overshoots).
#' @param max_treedepth Integer governing steps used by
#'   gradient descent; higher values result in slower but
#'   more stable estimation.
#' @param cores Number of cores to use for parallel processing.
#' @param chains Number of chains to estimate (higher number
#'   provides better estimate of convergence diagnostics).
#' @param alpha Cut-off for significance (will be adjusted for
#'   Monte Carlo).
#' @param quick_fit Logical; if \code{TRUE} adjusts settings
#'   for a quicker, less stable fit (good for initial model
#'   specification).
#'
#' @return A list with the algorithm settings.
#' \enumerate{
#'   \item total: The total number of posterior samples to draw.
#'   \item warmup: The number of iterations to use for the warm-up
#'     period before sampling from the posterior.
#'   \item iter: The number of iterations to use per chain, the sum
#'     of the warm-up iterations and the division of the total samples
#'     by the number of chains.
#'   \item adapth_delta: A proportion governing the gradient descent
#'     settings.
#'   \item max_treedepth: An integer governing the gradient descent
#'     settings.
#'   \item cut_off: The cut-off for statistical significant after
#'     correcting for Monte Carlo error (by taking the lower
#'     boundary of the 95% confidence interval around the
#'     specified alpha).
#' }
#'
#' @examples
#' # Default
#' alg <- algorithm_settings()
#'
#' # Settings for faster, less stable model fitting
#' alg <- algorithm_settings( desired_iter = 2000,
#'                            adapt_delta = .8,
#'                            max_treedepth = 12 )
#'
#' @export

# 5.1)
algorithm_settings <- function( desired_iter = 10000,
                                rng_seed = NULL,
                                warmup = 500,
                                adapt_delta = .9,
                                max_treedepth = 15,
                                cores = 4,
                                chains = 4,
                                alpha = .05,
                                quick_fit = FALSE ) {

  # Initialize list with desired settings
  algorithm <- list(
    total = desired_iter,
    warmup = warmup,
    iter = round( desired_iter/chains + warmup ),
    cores = cores,
    chains = chains,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    rng_seed = NA
  )
  # If no seed specified, create random seed
  if ( is.null( rng_seed ) ) {
    rng_seed <- round( 1e5 * runif( 1 ) )
  }
  algorithm$rng_seed = rng_seed

  if ( quick_fit ) {
    algorithm$total <- 2000
    algorithm$warmup <- 250
    algorithm$iter <- 250 + 2000 / chains
    algorithm$adapt_delta <- .8
    algorithm$max_treedepth <- 12
  }

  # Lower boundary of 95% confidence interval for MCMC error
  # around significance cut-off
  algorithm$cut_off <-
    qbinom( .025, algorithm$total, alpha ) / algorithm$total

  return( algorithm )
}

