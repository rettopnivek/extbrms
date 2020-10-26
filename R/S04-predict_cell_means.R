#' Predict Cell Means Based on Grouping Factors
#'
#' A convenience function to generate model predictions
#' over specified grouping factors for population-level
#' effects.
#'
#' @param est Output from the \code{\link{fit_with_brm}}
#'   function.
#' @param grp A vector with the column names for the
#'   grouping factors.
#' @param dm An optional data frame with the grouping
#'   factors and the weights for the population-level
#'   effects needed to compute the cell means.
#' @param lwr A list with 1) the column names for the
#'   grouping terms used in a mixed effects model, and
#'   2) any additional variables not part of the
#'   population-level effects.
#' @param ... Additional parameters for
#'   \code{\link[brms]{predict.brmsfit}}.
#'
#' @return Output from the \code{\link[brms]{predict.brmsfit}}
#'   function.
#'
#' @examples
#' # Load in example for linear mixed effects model
#' data("sleepstudy", package = "lme4")
#' df <- sleepstudy
#' df$LC.ZS.Days <- scale( df$Days )[,1]
#' # Fit model with linear contrast for days
#' # with subject-level intercepts and
#' # slopes
#' alg <- algorithm_settings( desired_iter = 5000,
#'                            cores = 1,
#'                            adapt_delta = .8,
#'                            rng_seed = 111 )
#' est <- fit_with_brm(
#'   df, Reaction ~ LC.ZS.Days + (1 + LC.ZS.Days|Subject),
#'   algorithm = alg
#' )
#'
#' # Predict averages per day (no subject-level
#' # uncertainty)
#' tst <- predict_cell_means(
#'   est, grp = 'Days'
#' )
#' print( round( tst, 1 ) )
#'
#' # Predict averages per day (with subject-level
#' # uncertainty)
#' tst <- predict_cell_means(
#'   est, grp = 'Days', lwr = list( 'ID', 'LC.ZS.Days' )
#' )
#' print( round( tst, 1 ) )
#'
#' @export

predict_cell_means = function( est,
                               grp = NULL, dm = NULL,
                               lwr = NULL,
                               ... ) {

  # Extract data used to when fitting the model
  df <- est$data

  # Extract posterior samples for
  # population-level effects
  P <- as.matrix( est$fit )

  # Column names for population-level effects
  iv <- colnames( P )
  iv <- iv[ grepl( 'b_', iv, fixed = T ) ]
  P <- P[,iv]

  # Number of posterior samples
  S <- nrow( P )

  # Match labels with those in data frame
  sel <- gsub( 'b_', '', iv, fixed = T )
  sel <- sel[ sel != 'Intercept' ]

  # If grouping terms (lower levels of
  # hierarchy) were specified
  if ( !is.null( lwr ) ) {
    # Add variables used in lower level
    if ( !is.null( lwr[[2]] ) ) {
      sel = unique( c( sel, lwr[[2]] ) )
    }
  }

  # Specify weights for computing
  # cell means
  if ( is.null( dm ) ) {

    if ( is.null( grp ) ) {
      stop( paste0(
        "Please provide vector 'grp' with ",
        "column names for grouping factors"
      ), call. = F )
    }

    if ( length( grp ) > 1 ) {
      dm <- aggregate(
        df[,sel], df[,grp],
        function(x) mean(x)
      )
    } else {
      dm <- aggregate(
        df[,sel], list( df[[grp]] ),
        function(x) mean(x)
      )
    }

    colnames( dm )[ 1:length( grp ) ] = grp
    colnames( dm )[ -( 1:length( grp ) ) ] = sel

  }

  # Number of cells
  nc = nrow( dm )

  # Initialize data frame for new data
  nd = df[ 1:nrow( dm ), ]

  # Loop over cells
  for ( i in 1:nc ) {

    # Update grouping factors
    for ( j in 1:length( grp ) ) {
      nd[[ grp[j] ]][i] = dm[[ grp[j] ]][i]
    }

    # Update predictors
    for ( j in 1:length( sel ) ) {
      nd[[ sel[j] ]][i] = dm[[ sel[j] ]][i]
    }

  }

  if ( is.null( lwr ) ) {
    # Exclude uncertainty from lower levels

    out = predict( est$fit, newdata = nd, re_formula = NA,
                   ... )

  } else {
    # Include uncertainty from lower levels

    for ( j in 1:length( lwr[[1]] ) ) {

      if ( is.factor( nd[[ lwr[[1]][j] ]] ) ) {
        nd[[ lwr[[1]][j] ]] =
          factor( nd[[ lwr[[1]][j] ]],
                  levels =
                    c( levels( nd[[ lwr[[1]][j] ]] ), 'NEW' )
          )
        nd[[ lwr[[1]][j] ]] = 'NEW'

      } else {
        nd[[ lwr[[1]][j] ]] = 'NEW'
      }

    }

    out = predict( est$fit, newdata = nd, allow_new_levels = T,
                   ... )

  }

  return( out )
}



