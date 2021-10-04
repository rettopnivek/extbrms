# General-purpose functions for estimation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2021-10-01

# Table of contents
# 1) algorithm_settings
# 2) Summary functions for Monte Carlo samples
#   2.1) p_value.monte_carlo
#   2.2) ci
# 3) pull_param
#   3.1) Custom inclusion/exclusion criteria
#   3.2) Pre-defined parameter subsets
#     3.2.1) Return all parameters
#     3.2.2) Return population-level parameters
#     3.2.3) Return cluster-level parameters
#     3.2.4) Return intercept parameter
#     3.2.5) Return variance parameters
#     3.2.6) Return covariance parameters
#     3.2.7) Return samples drawn from priors
#     3.2.8) Return unique parameter subsets
# 4) fit_with_brm

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

#### 2.1) p_value.monte_carlo ####
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
#' p_value.monte_carlo( x )
#' # One-sided comparison against 0
#' p_value.monte_carlo( x, alternative = 'greater' )
#' p_value.monte_carlo( x, alternative = 'less' )
#' # Non-zero comparison
#' p_value.monte_carlo( x, comparison = 1 )
#' p_value.monte_carlo( x, comparison = 1, alternative = 'greater' )
#' p_value.monte_carlo( x, comparison = 1, alternative = 'less' )
#'
#' @export

p_value.monte_carlo <- function( x,
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

#### 2.3) summary.monte_carlo ####
#' Title
#'
#' Description.
#'
#' @param x ...
#' @param form ...
#' @param tidy ...
#'
#' @details
#'
#' @return Output.
#'
#' @examples
#' # Examples
#'
#' @export

summary.monte_carlo <- function( x,
                                 form = NULL,
                                 tidy = TRUE ) {

  default_labels <- NULL

  if ( is.null( form ) ) {
    form <- '[[M]] ([[SD]])||[[CI.025]] to [[CI.975]]||[[P|3]]'
    default_labels <- c( 'Estimate_M_SD',
                         'CI_95',
                         'p_value' )
  }

  # Shorthand for 'grepl' function
  qg <- function(string, term) {
    return( grepl( term, string, fixed = TRUE ) )
  }

  # Shorthand for 'strsplit' function
  qs <- function(string, term) {
    strsplit( string, split = term, fixed = TRUE )[[1]]
  }

  # Function to determine where input should be
  # split up into separate columns based on the
  # symbol '||'
  find_columns <- function( form ) {

    if ( any( qg( form, '||' ) ) ) {
      return( qs( form, '||' ) )
    } else {
      return( form )
    }

  }

  # Function to split up input into separate
  # components for terms bracketed by '[[' and ']]'
  parse_input <- function( form ) {

    comp <- qs( form, ']]' )

    keep <- qg( comp, '[[' )
    comp <- comp[ keep ]

    comp <- sapply( comp, function(x) qs( x, '[[' )[2] )

    return( comp )
  }

  # ...
  transform_var <- function( x, type ) {

    types <- list(
      identity = c(
        "Identity", "identity",
        "id", "I", "i", "1"
      ),
      exponent = c(
        "Exponent", "exponent",
        "Exp", "exp",
        "E", "e", "2"
      ),
      logistic = c(
        "Logistic", "logistic",
        "L", "l", "3"
      )
    )

    if ( type %in% types$identity ) {
      return( x )
    }

    if ( type %in% types$exponent ) {
      return( exp( x ) )
    }

    if ( type %in% types$logistic ) {
      return( 1/( 1 + exp( -x ) ) )
    }

  }

  stat_from_symbol <- function( x, symb ) {

    # ...

    ns <- length( x )
    qnt <- qbinom( c( .025, .975 ), ns, .05 )/ns

    num <- as.character( diff( qnt ) )
    num <- qs( num, "0." )[2]
    num <- qs( num, "" )

    n_zeros <- min( which( num != "0" ) ) - 1

    if ( round( diff( qnt ), n_zeros ) != 0 ) {
      n_zeros <- n_zeros - 1
    }

    digits <- n_zeros

    # ...
    parts <- c( "", digits, "i" )

    # ...
    if ( qg( symb, '|' ) ) {

      val <- qs( symb, '|' )

      parts[1] <- val[1]

      if ( qg( parts[2], '!' ) ) {

        val[2] <- qs( val[2], '!' )[1]

      }

      parts[2] <- val[2]
    }

    if ( qg( '!', symb ) ) {

      val <- qs( symb, '!' )

      if ( parts[1] == "" ) {
        parts[1] <- val[1]
      }

      parts[3] <- val[2]

    }

    if ( parts[1] == "" ) {
      parts[1] <- symb
    }

    # ...

    if ( parts[1] == "M" ) {

      out <- mean( transform_var( x, parts[3] ) )

    }

    if ( parts[1] %in% c( "Md" ) ) {

      out <- median( transform_var( x, parts[3] ) )

    }

    if ( parts[1] == "SD" ) {

      out <- sd( transform_var( x, parts[3] ) )

    }

    if ( parts[1] %in% c( "P" ) ) {

      out <- p_value.monte_carlo( transform_var( x, parts[3] ) )

    }

    if ( qg( parts[1], "CI." ) ) {

      p <- as.numeric( qs( parts[1], "CI" )[2] )

      out <- quantile( transform_var( x, parts[3] ), probs = p )

    }

    out <- format(
      out,
      digits = as.numeric( parts[2] ),
      nsmall = as.numeric( parts[2] )
    )

    return( out )
  }

  # Split up into separate columns
  all_forms <- find_columns( form )

  # Number of rows (one per each parameter)
  J <- ncol( x )

  # Number of of columns
  K <- length( all_forms )

  # Initialize matrix for results
  out <- matrix( '', J, K )

  # Loop over rows
  for ( j in 1:J ) {

    # Loop over columns
    for ( k in 1:K ) {

      output <- all_forms[k]

      all_symbs <- parse_input( output )

      for ( i in 1:length( all_symbs ) ) {

        val <- stat_from_symbol( x[,j], all_symbs[i] )
        output <- gsub( paste0( '[[', all_symbs[i], ']]' ),
                        val, output, fixed = TRUE )

      }

      out[j,k] <- output

      # Close 'Loop over columns'
    }

    # Close 'Loop over rows'
  }

  if ( is.null( default_labels ) ) {
    default_labels <- paste0( 'C_', 1:K )
  }

  out <- cbind( colnames( x ),
                out )
  colnames( out ) <- c(
    'Variable',
    default_labels
  )

  out <- data.frame( out, stringsAsFactors = F )

  if ( tidy ) {

    out$Variable <- gsub( 'b_', '', out$Variable )
    out$Variable <- gsub( 'sd_', '', out$Variable )
    out$Variable <- gsub( 'cor_', '', out$Variable )
    out$Variable <- gsub( '_', ' ', out$Variable )
    out$Variable <- gsub( '.', ' ', out$Variable )
    out$Variable <- gsub( 'IV.', '', out$Variable )
    out$Variable <- gsub( 'DC.', '', out$Variable )
    out$Variable <- gsub( 'EC.', '', out$Variable )
    out$Variable <- gsub( 'ZS.', '', out$Variable )
    out$Variable <- gsub( 'RW.', '', out$Variable )

  }

  return( out )
}

#### 3) pull_param ####
#' Title
#'
#' Description.
#'
#' @param prm ...
#' @param type ...
#' @param include ...
#' @param exclude ...
#' @param special ...
#' @param and ...
#'
#' @details
#'
#' @return Output.
#'
#' @examples
#' # Examples
#'
#' @export

pull_param <- function( x,
                        type = 'b_',
                        include = NULL,
                        exclude = NULL,
                        special = NULL,
                        and = c( TRUE, FALSE ) ) {

  prm <- NULL
  ps <- NULL

  if ( class( x )[1] == 'brmsfit' ) {
    prm <- variables( x )
    ps <- as.matrix( x )
  }

  if ( class( x )[1] == 'matrix' ) {
    prm <- colnames( x )
    ps <- x
  }


  if ( class( x )[1] == 'character' ) {
    prm <- x
  }

  if ( is.null( prm ) ) {

    err_msg <- paste0(
      "First argument must be an object of class ",
      "'brmsfit', 'matrix', or 'character'"
    )
    stop( err_msg )

  }

  # Shorthand for 'grepl' function
  qg <- function(x, y) {
    entries <- grepl( x, y, fixed = TRUE )
    return( entries )
  }

  #### 3.1) Custom inclusion/exclusion criteria ####

  # If inclusion/exclusion criteria were specified
  if ( !is.null( include ) | !is.null( exclude ) ) {

    # Number of parameter names
    n = length( prm )

    # Initialize output
    log_vec <- rep( FALSE, length( prm ) )

    # If inclusion criteria were specified
    if ( !is.null( include ) ) {

      # Number of criteria
      ni <- length( include )

      # Initialize logical matrix
      log_mat <- matrix( FALSE, n, ni )

      # Loop over criteria
      for ( k in 1:ni ) {

        # Check if criteria met
        log_mat[,k] = qg( include[k], prm )

        # Output for and/or
        if ( and[1] ) {
          log_vec <- rowSums( log_mat ) == ni
        } else {
          log_vec <- rowSums( log_mat ) > 0
        }

        # Close 'Loop over criteria'
      }

      # Close 'If inclusion criteria were specified'
    }

    # If exclusion criteria were specified
    if ( !is.null( exclude ) ) {

      # Number of criteria
      ne <- length( exclude )

      # Initialize logical matrix
      log_mat <- matrix( FALSE, n, ne )

      # Loop over criteria
      for ( k in 1:ne ) {

        # Check if criteria met
        log_mat[,k] = qg( exclude[k], prm )

        # Output for and/or
        if ( and[2] ) {
          log_vec <-
            log_vec &
            !( rowSums( log_mat ) == ne )
        } else {
          log_vec <-
            log_vec &
            !( rowSums( log_mat ) > 0 )
        }

        # Close 'Loop over criteria'
      }

      # Close 'If exclusion criteria were specified'
    }

    # If matrix of posterior samples was provided
    if ( !is.null( ps ) ) {
      return( ps[, prm[ log_vec ] ] )
    }

    # Return subset of names
    return( prm[ log_vec ] )

    # Close 'If inclusion/exclusion criteria were specified'
  }

  #### 3.2) Pre-defined parameter subsets ####

  # Possible inputs to specify subset of
  # parameters to returns
  types <- list(
    all = c( 'All', 'all', 'a' ),
    population = c('Population', 'population',
                   'Pop', 'pop', 'p',
                   '2nd', 'Second', 'second',
                   'b_' ),
    cluster = c( 'Cluster', 'cluster',
                 'Group', 'group',
                 'Group-level', 'group-level',
                 '1st', 'First', 'first',
                 'r_' ),
    intercept = c( 'Intercept', 'intercept',
                   'B0', 'b0', 'i' ),
    variance = c( 'Variance', 'variance',
                  'Var', 'var', 'v', 'sd_', 'v' ),
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

  # Initialize output
  out <- NULL

  #### 3.2.1) Return all parameters ####

  # Note: Excludes priors and log-probability
  if ( type %in% types$all ) {

    entries <-
      !qg( 'prior', prm ) &
      !qg( 'lp__', prm )

    if ( any( entries ) ) {
      out <- prm[ entries ]
    }

    # Close 'Return all parameters'
  }

  #### 3.2.2) Return population-level parameters ####

  if ( type %in% types$population ) {

    entries <-
      qg( 'b_', prm ) &
      !qg( 'prior', prm ) &
      !qg( 'lp__', prm )

    if ( any( entries ) ) {

      out <- prm[ entries ]

      # Make sure intercept term is always first
      if ( any( out == 'b_Intercept' ) ) {

        out <- c(
          out[ out == 'b_Intercept' ],
          out[ out != 'b_Intercept' ]
        )

      }

    }

    # Close 'Return population-level parameters'
  }

  #### 3.2.3) Return cluster-level parameters ####

  if ( type %in% types$cluster ) {

    entries <-
      qg( 'r_', prm ) &
      !qg( 'prior', prm ) &
      !qg( 'lp__', prm )

    if ( any( entries ) ) {
      out <- prm[ entries ]
    }

    # Close 'Return cluster-level parameters'
  }

  #### 3.2.4) Return intercept parameter ####

  if ( type %in% types$intercept ) {

    entries <-
      qg( 'b_Intercept', prm ) &
      !qg( 'prior', prm ) &
      !qg( 'lp__', prm )

    if ( any( entries ) ) {
      out <- prm[ entries ]
    }

    # Close 'Return intercept parameter'
  }

  #### 3.2.5) Return variance parameters ####

  if ( type %in% types$variance ) {

    entries <-
      qg( 'sd_', prm ) &
      !qg( 'prior', prm ) &
      !qg( 'lp__', prm )

    if ( any( entries ) ) {
      out <- prm[ entries ]
    }

    # Close 'Return variance parameters'
  }

  #### 3.2.6) Return covariance parameters ####

  if ( type %in% types$covariance ) {

    entries <-
      qg( 'cor_', prm ) &
      !qg( 'prior', prm ) &
      !qg( 'lp__', prm )

    if ( any( entries ) ) {
      out <- prm[ entries ]
    }

    # Close 'Return covariance parameters'
  }

  #### 3.2.7) Return samples drawn from priors ####

  if ( type %in% types$prior ) {

    entries <-
      qg( 'prior', prm ) &
      !qg( 'lp__', prm )

    if ( any( entries ) ) {
      out <- prm[ entries ]
    }

    # Close 'Return samples drawn from priors'
  }


  #### 3.2.8) Return unique parameter subsets ####

  if ( type %in% types$special ) {

    if ( is.null( special ) ) {
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
    }

    log_vec <- rep( FALSE, length( prm ) )
    for ( k in 1:length( special ) ) {

      entries <-
        qg( special[k], prm ) &
        qg( 'prior', prm ) &
        !qg( 'lp__', prm )
      log_vec[ entries ] <- TRUE

    }
    entries <- log_vec

    if ( any( entries ) ) {
      out <- prm[ entries ]
    }

    # Close 'Return unique parameter subsets'
  }

  if ( !is.null( out ) ) {

    # If matrix of posterior samples was provided
    if ( !is.null( ps ) ) {
      return( ps[, out ] )
    }

    return( out )
  }

  stop( 'Type of parameter not found' )
}

#### 4) fit_with_brm ####
#' Fit 'brm' Model to Data
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
#' @param predict.probs A vector of probabilities, input passed
#'   on to the \code{probs} argument in
#'   \code{\link[brms]{predict.brmsfit}}.
#' @param ... Additional arguments to the \code{\link[brms]{brm}}
#'   function.
#'
#' @details Check the help page for \code{\link[brms]{brm}} for
#' more details.
#'
#' @return A list including the elements...
#' \itemize{
#'   \item data: The data frame used as input to the \code{brm} function.
#'   \item priors: The priors extracted via the \code{prior_summary} function.
#'   \item fit: The output from the \code{brm} function.
#'   \item algorithm: The list with the algorithm settings.
#' }
#'
#' @examples
#' # Example for simple linear regression
#' data( cars )
#' # Setting for algorithm for quicker estimation
#' alg <- algorithm_settings( quick_fit = T, rng_seed = 111 )
#' # Fit model using 'brm'
#' est <- fit_with_brm( cars, dist ~ speed, algorithm = alg )
#'
#' # Example for logistic regression
#' # Simulate data
#' X <- cbind( 1, rep( 0:1, each = 100 ) )
#' p <- 1 / ( 1 + exp( -( X %*% c( -1, .5 ) ) ) )
#' df <- data.frame(
#'   Y = rbinom( nrow( X ), 1, p ),
#'   X = X[,2]
#' )
#'
#' # Empirical Bayes prior on intercept
#' obs <- mean( df$Y[ df$X == 0 ] )
#' obs <- round( log( obs / (1 - obs) ) , 2 )
#' prs <- c(
#'   prior_string( paste0( 'normal( ', obs, ', 2.5 )' ), class = 'Intercept' ),
#'   prior_string( 'normal( 0.0, 2.5 )', class = 'b' )
#' )
#'
#' # Setting for algorithm for quicker estimation
#' alg <- algorithm_settings( quick_fit = T, rng_seed = 112 )
#' est <- fit_with_brm(
#'   df, Y ~ X,
#'   family = bernoulli( link = 'logit' ),
#'   prior = prs,
#'   algorithm = alg
#' )
#'
#' @export

fit_with_brm <- function( data, formula,
                          prior = NULL, algorithm = NULL,
                          sample_prior = 'yes',
                          output = NULL,
                          track_time = T,
                          predict.probs = NULL,
                          ... ) {

  if ( is.null( algorithm ) ) {
    algorithm <- extbrms::algorithm_settings()
  }

  if ( is.null( predict.probs ) ) {
    predict.probs <- c(
      .025, # Lower boundary - 95%
      .05, # Lower boundary - 90%
      pnorm( -1 ), # Lower boundary - 1 SD for normal
      .25, # Lower boundary - 50%
      .5, # Median
      .75, # Upper boundary - 50%
      pnorm( 1 ), # Upper boundary - 1 SD for normal
      .95, # Upper boundary - 90%
      .975 # Upper boundary - 95%
    )
  }

  if ( is.null( output ) ) {

    output <- list(
      data = data,
      priors = NULL,
      fit = NULL,
      predicted = NULL,
      algorithm = algorithm,
      parameters = list(
        all = '',
        population = ''
      )
    )

  } else {
    output$data <- data
    output$algorithm <- algorithm
  }

  if ( is.null( output$time ) & track_time ) {
    output$time <- Sys.time()
  }

  output$fit <- brm(
    formula,
    data = data,
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

  predicted_y <- predict( output$fit, probs = predict.probs )
  prd_lbl <- paste0( 'Y.hat.LB.', round( predict.probs, 3 ) )
  prd_lbl[ predict.probs > .5 ] <-
    gsub( '.LB.', '.UB.', prd_lbl[ predict.probs > .5 ], fixed = TRUE )

  colnames( predicted_y ) <- c(
    'Y.hat.Mean',
    'Y.hat.SD',
    prd_lbl
  )

  if ( nrow( predicted_y ) == nrow( output$data ) ) {
    output$data <- cbind( output$data, predicted_y )
  } else {
    wrn_msg <- paste0(
      "Number of predicted observations does not match ",
      "number of rows in data passed to 'brm' function - ",
      "check for missing data in outcomes and predictors."
    )
    warning( wrn_msg )
  }

  output$priors <- prior_summary( output$fit )

  if ( track_time ) {
    output$time <- Sys.time() - output$time
  }

  return( output )
}




