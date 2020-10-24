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
#'   \item predicted: The output from the \code{predict.brmsfit} function.
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
                          ... ) {

  if ( is.null( algorithm ) ) {
    algorithm <- algorithm_settings()
  }

  if ( is.null( output ) ) {

    output <- list(
      data = data,
      priors = NULL,
      fit = NULL,
      predicted = NULL,
      algorithm = algorithm
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
  output$predicted <- predict( output$fit )
  output$priors <- prior_summary( output$fit )

  if ( track_time ) {
    output$time <- Sys.time() - output$time
  }

  return( output )
}




