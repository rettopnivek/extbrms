# General-purpose functions for estimation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-01-22

# Table of contents
# 1) algorithm_settings
# 2) Summary functions for Monte Carlo samples
#   2.1) pv
#   2.2) ci
# 3) %s%
#   3.1) Determine class and extract names
#   3.2) Define functions
#   3.3) Convert to R expression
# 4) pull_param
#   4.1) Determine class and extract samples
#   4.2) User-defined subset
#   4.3) Pre-defined parameter subsets
# 5) fit_with_brm
#   5.1) Setup
#   5.2) Fit model
#   5.3) Model predictions
#   5.4) Additional processing

#### 1) algorithm_settings ####
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
#' alg <- algorithm_settings( quick_fit = TRUE )
#'
#' @export

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

#### 2) Summary functions for Monte Carlo samples ####

#### 2.1) pv ####
#' Compute P-values From Monte Carlo Samples
#'
#' A general-purpose function that computes the
#' proportion of values from a sample generated
#' via Monte Carlo methods that are above or below
#' a comparison point. Can be used, for example,
#' to estimate a posterior p-value from MCMC
#' samples.
#'
#' @param x A vector of values.
#' @param alternative A character string specifying the
#'   alternative hypothesis, either \code{'two-sided'}
#'   (default), \code{'greater'}, or \code{'less'}.
#'   You can specify just the initial letter.
#' @param comparison A cut-off value; the function
#'   computes the proportion of values in \code{x} that fall
#'   above or below \code{comparison}.
#'
#' @return A proportion between 0 and 1.
#'
#' @examples
#' set.seed( 831 ) # For reproducibility
#' x <- rnorm( 1000, mean = 1 )
#' # Two-sided comparison against 0
#' pv( x )
#' # One-sided comparison against 0
#' pv( x, alternative = 'greater' )
#' pv( x, alternative = 'less' )
#' # Non-zero comparison
#' pv( x, comparison = 1 )
#' pv( x, comparison = 1, alternative = 'greater' )
#' pv( x, comparison = 1, alternative = 'less' )
#'
#' @export

pv <- function( x,
                alternative = 'two-sided',
                comparison = 0 ) {

  check <- TRUE

  # Two-sided test
  if (alternative %in% c( "t", "two-sided" )) {
    check <- FALSE

    if (median(x) > comparison) {
      out <- mean(x < comparison)
    } else {
      out <- mean(x > comparison)
    }
    out <- out * 2
  }

  # Test if greater than comparison
  if (alternative %in% c( "g", "greater" ) ) {
    check <- FALSE
    out <- mean(x < comparison)
  }

  # Test if less than comparison
  if (alternative %in% c( "l", "less" ) ) {
    check <- FALSE
    out <- mean(x > comparison)
  }

  # Informative error message if
  # 'alternative' misspecified
  if (check) {
    err_msg <- paste(
      "Please specify 'alternative' as",
      "either 'two-sided', 'greater', or",
      "'less'."
    )
    stop(err_msg)
  }

  return( out )
}

#### 2.2) ci ####
#' Compute Credible Intervals From Monte Carlo Samples
#'
#' Function to compute credible intervals of a desired
#' width from a vector of Monte Carlo samples via
#' the base R \code{\link[stats]{quantile}} function.
#'
#' @param x A vector of values, Monte Carlo samples for
#'   a prior or posterior distribution.
#' @param width The desired width of the credible interval.
#' @param index The index for the bounds to return, either
#'   \code{1} for the lower bound, \code{2} for the upper
#'   bound, for the vector \code{1:2} for both (default).
#'
#' @return A vector giving the lower and/or upper bound
#'   for the credible interval of the specified width.
#'
#' @examples
#' set.seed( 831 ) # For reproducibility
#' x <- rnorm( 1000, mean = 1 )
#' # 95% CI
#' ci( x )
#' # Approx. one standard deviation from
#' # the mean (interval is close to [0,1])
#' ci( x, width = .68 )
#' # Lower or upper bound
#' ci( x, .68, 1 )
#' ci( x, .68, 2 )
#'
#' @export

ci <- function( x, width = .95, index = 1:2 ) {

  bnds <- .5 + c( -.5, .5 )*width

  q <- quantile( x, probs = bnds )

  return( q[index] )
}

#### 3) %s% ####
#' Operator to Select Parameters
#'
#' An operator that takes an input
#' (such as a \code{\link[brms]{brmsfit}}
#' object, matrix, or output from the
#' \code{\link{fit_with_brm}} function)
#' and extracts a character vector with
#' the parameter names that meet a
#' specified set of criteria.
#'
#' @param x Either a \code{\link[brms]{brmsfit}}
#'   object,a matrix of posterior samples, a
#'   character vector, or a list with an
#'   element \code{fit} that is a
#'   \code{\link[brms]{brmsfit}} object.
#' @param y A character string, to be
#'   converted into an R expression used
#'   to match and select parameter names.
#'
#' @details The operator converts the character
#' string on the right-hand side into an
#' expression making calls to
#' \code{\link[base]{grepl}} and
#' \code{&}, \code{|}, and \code{!}.
#'
#' @return A character vector.
#'
#' @examples
#' # Example vector of parameter names
#' x <- c(
#'   'b_Intercept', 'b_speed', 'sigma',
#'   'prior_Intercept', 'prior_sigma', 'lp__'
#' )
#'
#' # Extract all parameters
#' x %s% ''
#'
#' # Extract parameters with 'b_'
#' x %s% 'b_'
#'
#' # Extract parameters with both 'b_' and 'speed'
#' x %s% 'b_ & speed'
#'
#' # Extract parameters with 'b_' but not 'speed'
#' x %s% 'b_ & - speed'
#'
#' # Extract parameters with 'b_' or 'sigma'
#' x %s% 'b_ | sigma'
#'
#' # Combinations of '&' and '|'
#' x %s% '( b_ | sigma ) & - prior'
#'
#' @export

`%s%` <- function( x, y) {

  #### 3.1) Determine class and extract names ####

  prm <- NULL

  if ( class( x )[1] == 'brmsfit' ) {
    prm <- variables( x )
  }

  if ( class( x )[1] == 'matrix' ) {
    prm <- colnames( x )
  }

  if ( class( x )[1] == 'character' ) {
    prm <- x
  }

  if ( class( x )[1] == 'list' ) {
    if ( !is.null( x$fit ) ) {
      prm <- variables( x$fit )
    }
  }

  if ( is.null( prm ) ) {

    err_msg <- paste0(
      "First argument must be an object of class ",
      "'brmsfit', 'matrix', 'character', or 'list'"
    )
    stop( err_msg )

  }

  #### 3.2) Define functions ####

  # Shorthand for 'grepl' function
  qg <- function(x, y) {
    entries <- grepl( x, y, fixed = TRUE )
    return( entries )
  }

  # Shorthand for 'strsplit' function
  qs <- function(x, y) {
    out <- strsplit( x, split = y, fixed = TRUE )[[1]]
    return( out )
  }

  #### 3.3) Convert to R expression ####

  if ( y == '' ) {
    return( prm )
  }

  chr_vec <- qs( y, ' ' )
  mod <- chr_vec %in% c( '-', '|', '&', '(', ')' )
  chr_vec[ !mod ] <-
    paste0( "qg('", chr_vec[ !mod ], "',prm) " )

  chr_vec[ chr_vec == '-' ] <- '!'

  chr_vec[ chr_vec == '|' ] <- '| '
  chr_vec[ chr_vec == '&' ] <- '& '

  chr_vec <- paste( chr_vec, collapse = '' )

  return( prm[ eval( parse( text = chr_vec ) ) ] )
}

#### 4) pull_param ####
#' Get Posterior Samples for Specified Parameters
#'
#' Function to extract the matrix of posterior
#' samples from a \code{\link[brms]{brmsfit}}
#' object for a specified subset of parameters.
#'
#' @param x Either a \code{\link[brms]{brmsfit}}
#'   object,a matrix of posterior samples, or a
#'   list with an element \code{fit} that is a
#'   \code{\link[brms]{brmsfit}} object.
#' @param expr An optional character string,
#'   to be converted into an R expression used
#'   to match and select parameter names.
#' @param predefined A character string, matched
#'   against a set of pre-defined labels for
#'   different types of parameters.
#'
#' @return A matrix of posterior samples.
#'
#' @examples
#'
#' \dontrun{
#' data( cars )
#' # Setting for algorithm for quicker estimation
#' alg <- algorithm_settings( quick_fit = T, rng_seed = 111 )
#' # Fit model using 'brm'
#' est <- fit_with_brm( cars, dist ~ speed, algorithm = alg, status = T )
#' # Extract regression coefficients
#' PS <- pull_param( est )
#' # Extract residual standard deviation
#' PS <- pull_param( est, predefined = 'special' )
#' }
#'
#' @export

pull_param <- function( x, expr = NULL,
                        predefined = 'b_' ) {

  #### 4.1) Determine class and extract samples ####

  ps <- NULL

  if ( class( x )[1] == 'brmsfit' ) {
    ps <- as.matrix( x )
  }

  if ( class( x )[1] == 'matrix' ) {
    ps <- x
  }

  if ( class( x )[1] == 'list' ) {
    if ( !is.null( x$fit ) ) {
      ps <- as.matrix( x$fit )
    }
  }

  if ( is.null( ps ) ) {

    err_msg <- paste0(
      "First argument must be an object of class ",
      "'brmsfit', 'matrix', or 'list'"
    )
    stop( err_msg )

  }

  #### 4.2) User-defined subset ####

  if ( !is.null( expr ) ) {

    prm <- x %s% expr

    return( ps[,prm] )

  }

  #### 4.3) Pre-defined parameter subsets ####

  # Possible inputs to specify subset of
  # parameters to returns
  types <- list(
    all = c( 'All', 'all', 'a' ),
    population = c('Population', 'population',
                   'Pop', 'pop', 'p',
                   '2nd', 'Second', 'second',
                   'b_' ),
    cluster = c( 'Group', 'group',
                 'Cluster', 'cluster',
                 'Group-level', 'group-level',
                 '1st', 'First', 'first',
                 'r_' ),
    intercept = c( 'Intercept', 'intercept',
                   'B0', 'b0', 'i' ),
    variance = c( 'Variance', 'variance',
                  'Var', 'var', 'v', 'sd_' ),
    covariance = c( 'Covariance', 'covariance',
                    'Correlation', 'correlation',
                    'Cov', 'cov', 'Cor', 'cor',
                    'cor_', 'c' ),
    special = c( 'Special', 'special',
                 's' ),
    prior = c( 'Priors', 'priors',
               'Prior', 'prior',
               'pr' )
  )

  prm <- NULL

  if ( predefined %in% types$population ) {
    prm <- x %s% 'b_ & - prior_ & - [ & - __'
  }

  if ( predefined %in% types$cluster ) {
    prm <- x %s% 'r_ & - prior_'
  }

  if ( predefined %in% types$intercept ) {
    prm <- x %s% 'b_Intercept & - prior_ & - [ & - __'
  }

  if ( predefined %in% types$variance ) {
    prm <- x %s% 'sd_ & - prior_'
  }

  if ( predefined %in% types$covariance ) {
    prm <- x %s% 'cor_ & - prior_ & - ['
  }

  if ( predefined %in% types$prior ) {
    prm <- x %s% 'prior_'
  }

  if ( predefined %in% types$special ) {

    special <- c(
      'sigma',
      'phi',
      'nu',
      'shape',
      'ndt',
      'alpha',
      'xi',
      'beta',
      'kappa'
    )

    special <- paste(
      special, collapse = ' | '
    )

    prm <- x %s% paste0( '( ', special, ' ) & - prior_' )

  }

  if ( is.null( prm ) ) {
    stop( 'Input does not match any pre-defined options' )
  }

  return( ps[, prm ] )
}

#### 5) fit_with_brm ####
#' Fit Model to Data Using 'brm'
#'
#' A convenience function that makes a call to the
#' \code{\link[brms]{brm}} function with the
#' specified user settings. Arguments are designed
#' to improve compatibility with packages like
#' 'dplyr', and to streamline specification
#' of algorithm settings.
#'
#' @param data An object of class \code{data.frame} (or one that can
#'   be coerced to that class) containing data of all variables used
#'   in the model.
#' @param formula An object of class \code{\link[stats]{formula}},
#'   \code{\link[brms]{brmsformula}}, or \code{\link[brms]{mvbrmsformula}}
#'   (or one that can be coerced to that classes):
#'   A symbolic description of the model to be fitted. The
#'   details of model specification are explained in
#'   \code{\link[brms]{brmsformula}}.
#' @param prior One or more \code{brmsprior} objects created by
#'   \code{\link[brms]{set_prior}} or related functions and
#'   combined using the \code{c} method or the \code{+} operator.
#'   See also \code{\link[brms]{get_prior}} for more help.
#' @param algorithm A list specifying the number of iterations to run,
#'   the warm-up period, the number of chains and cores to use,
#'   the settings for the \code{adapt_delta} and \code{max_treedepth}
#'   control options, and the seed to use for reproducibility. Output
#'   from the \code{\link[extbrms]{algorithm_settings}} function.
#' @param sample_prior A character string, either code{'no'},
#'   (no sampling from the prior distributions), \code{'yes'}
#'   (sample from both the posterior and the prior distributions),
#'   or \code{'only'} (sample only from the prior distribution).
#' @param output An optional list to which the output should be
#'   added.
#' @param track_time Logical; if \code{TRUE}, tracks the run time
#'   of the function.
#' @param predict_probs A vector of probabilities, input passed
#'   on to the \code{probs} argument in \code{\link[brms]{predict.brmsfit}}.
#' @param status Logical; if \code{TRUE} prints brief
#'   messages to the console to track progress for debugging
#'   purposes.
#' @param ... Additional arguments to the \code{\link[brms]{brm}}
#'   function.
#'
#' @details Check the help page for \code{\link[brms]{brm}} for
#' more details.
#'
#' @return A list including the elements...
#' \itemize{
#'   \item data: The data frame used as input to the
#'     \code{brm} function;
#'   \item priors: The priors extracted via the
#'     \code{prior_summary} function;
#'   \item fit: The output from the \code{brm} function;
#'   \item predicted: A matrix simulated data based on
#'     the sets posterior samples (rows) per each observation
#'     in the data set (columns);
#'   \item algorithm: The list with the algorithm settings;
#'   \item parameters: A list with subsets of parameter names;
#'   \item time: How long it took to fit the model.
#' }
#'
#' @examples
#' # Example for simple linear regression
#' \dontrun{
#' data( cars )
#' # Setting for algorithm for quicker estimation
#' alg <- algorithm_settings( quick_fit = T, rng_seed = 111 )
#' # Fit model using 'brm'
#' est <- fit_with_brm( cars, dist ~ speed, algorithm = alg )
#' }
#'
#' # Example for logistic regression
#' \dontrun{
#' data( mtcars )
#' # Standardize predictors
#' dtf <- mtcars[,c('vs','wt','disp')]
#' dtf$wt <- scale( dtf$wt )[,1]
#' dtf$disp <- scale( dtf$disp )[,1]
#' # Setting for algorithm for quicker estimation
#' alg <- algorithm_settings( quick_fit = T, rng_seed = 112 )
#' # Specify priors
#' prs <- get_prior(
#'   vs ~ wt + disp,
#'   family = bernoulli( link = 'logit' ),
#'   data = mtcars
#' )
#' # Empirical Bayes prior on intercept
#' p <- mean( dtf$vs )
#' lo <- log( p / ( 1 - p ) )
#' prs$prior[ prs$class == 'Intercept' ] <-
#'   paste0( 'normal( ', round( lo, 2 ), ', 2.5 )' )
#' # Weakly regularizing priors on coefficients
#' prs$prior[ prs$class == 'b' & prs$coef != '' ] <-
#'   'normal( 0.0, 1.0 )'
#' est <- fit_with_brm(
#'   dtf, vs ~ wt + disp,
#'   family = bernoulli( link = 'logit' ),
#'   prior = prs,
#'   algorithm = alg
#' )
#' }
#'
#' @export

fit_with_brm <- function( data, formula,
                          prior = NULL, algorithm = NULL,
                          sample_prior = 'yes',
                          output = NULL,
                          track_time = T,
                          predict_probs = NULL,
                          status = FALSE,
                          ... ) {

  #### 5.1) Setup ####

  if ( status ) message( 'Start estimation' )

  if ( is.null( algorithm ) ) {
    if ( status )
      message( '- Settings for estimation algorithm' )

    algorithm <- extbrms::algorithm_settings()
  }

  if ( is.null( predict_probs ) ) {
    if ( status )
      message( '- Default prediction intervals' )
    predict_probs <- c(
      .025, # Lower boundary - 95%
      .05, # Lower boundary - 90%
      round( pnorm( -1 ), 3 ), # Lower boundary - 1 SD for normal
      .25, # Lower boundary - 50%
      .75, # Upper boundary - 50%
      round( pnorm( 1 ), 3 ), # Upper boundary - 1 SD for normal
      .95, # Upper boundary - 90%
      .975 # Upper boundary - 95%
    )
  }

  if ( is.null( output ) ) {
    if ( status )
      message( '- Initialize output' )
    output <- list(
      data = data,
      priors = NULL,
      fit = NULL,
      predicted = NULL,
      algorithm = algorithm,
      parameters = list()
    )

  } else {
    output$data <- data
    output$algorithm <- algorithm
  }

  #### 5.2) Fit model ####

  if ( is.null( output$time ) & track_time ) {
    if ( status )
      message( '- Track estimation time' )
    output$time <- Sys.time()
  }

  if ( status )
    message( '- Fit model to data' )

  output$fit <- brm(
    formula,
    data = data,
    prior = prior,
    iter = algorithm$iter,
    warmup = algorithm$warmup,
    chains = algorithm$chains,
    cores = algorithm$cores,
    sample_prior = sample_prior,
    control = list(
      adapt_delta = algorithm$adapt_delta,
      max_treedepth = algorithm$max_treedepth
    ),
    seed = algorithm$rng_seed,
    ...
  )

  #### 5.3) Model predictions ####

  if ( status )
    message( '- Generate model predictions' )

  # Simulated data for each set of posterior samples
  output$predicted <- posterior_predict( output$fit )

  n_sim <- ncol( output$predicted )
  n_obs <- nrow( output$data )

  if ( n_sim != n_obs ) {

    wrn_msg <- paste0(
      "Number of predicted observations does not match ",
      "number of rows in data passed to 'brm' function - ",
      "check for missing data in outcomes and predictors."
    )
    warning( wrn_msg )

  } else {

    # Check if predictions are for categorical/ordinal data

    # Predict only one row of data and use
    # only 100 draws from posterior for speed
    nd <- output$data[1,]
    check <- predict( output$fit, newdata = nd, draw_ids = 1:100 )

    # Categorical or ordinal data
    if ( !all( c( 'Estimate', 'Q2.5', 'Q97.5' ) %in%
               colnames( check ) ) ) {

      if ( status )
        message( '- Predictions for categorical or ordinal data' )

      # Close 'Categorical or ordinal data'
    } else {

      if ( status )
        message( '- Predictions for standard data' )

      # Add predictions to saved data
      output$data$Y.hat.Mean <-
        apply( output$predicted, 2, mean )
      output$data$Y.hat.Median <-
        apply( output$predicted, 2, median )
      output$data$Y.hat.SD <-
        apply( output$predicted, 2, sd )

      K <- length( predict_probs )
      if ( K > 0 ) {
        prd <- matrix( NA, n_obs, K )
        for ( k in 1:K ) {
          prd[,k] <-
            apply( output$predicted, 2, quantile,
                   prob = predict_probs[k] )
        }
      }

      nms <- rep( 'Y.hat.', K )
      i <- predict_probs < .5
      nms[i] <- paste0( nms[i], '.LB.', predict_probs[i] )
      i <- predict_probs >.5
      nms[i] <- paste0( nms[i], '.UB.', predict_probs[i] )

      colnames( prd ) <- nms
      output$data <- cbind( output$data, prd )

      # Close else for 'Categorical or ordinal data'
    }

  }

  #### 5.4) Additional processing ####

  if ( status )
    message( '- Extract priors for model' )

  output$priors <- prior_summary( output$fit )

  if ( status )
    message( '- Additional processing' )

  # Create subsets of parameter names for
  # easier extraction

  all_param <- variables( output$fit )

  cls <- unique( output$priors$class )

  cls <- cls[ cls != '' ]

  lst <- lapply( cls, function(s)
    all_param %s% paste0( s, ' & - ( prior_ | [ )' ) )
  names( lst ) <- cls

  # Remove standard deviations for group terms
  # from 'Intercept' category
  if ( 'Intercept' %in% names( lst ) ) {
    lst$Intercept <-
      lst$Intercept %s% '- sd_'
  }

  # Save all parameters
  lst$all <- all_param

  # Save priors
  if ( any( grepl( 'prior_', all_param, fixed = TRUE ) ) ) {
    lst$priors <- all_param %s% 'prior_'
  }

  output$parameters <- lst

  if ( track_time ) {
    output$time <- Sys.time() - output$time
  }

  return( output )
}
