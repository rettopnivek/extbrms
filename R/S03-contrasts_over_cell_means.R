#' Compute Omnibus Tests for Cell Means
#'
#' Convenience function that computes
#' an omnibus test over a set of
#' contrasts. First, posterior samples
#' over cell means (based on grouping
#' factors in the fitted data) are
#' extracted. Next, user-specified
#' contrasts are computed over the cell
#' means. Finally, an omnibus test checking
#' whether all contrasts are equal to
#' zero is computed.
#'
#' @param est Output from the \code{\link{fit_with_brm}}
#'   function.
#' @param contrast_list A list of list, each internal
#'   list specifying a contrast to compute over the
#'   cell means, where \code{I} gives the column indices
#'   to include from the cell means, and \code{W} gives
#'   the weights.
#' @param grp A vector with the column names for the
#'   grouping factors.
#' @param dm An optional data frame with the grouping
#'   factors and the weights for the population-level
#'   effects needed to compute the cell means.
#' @param digits The number of digits to round to when
#'   reporting posterior p-values.
#'
#' @return A list with the estimate, estimated
#'   posterior p-value, a character string with
#'   a nicely formatted p-value, and a matrix with
#'   the posterior samples for the estimated cell means.
#'
#' @examples
#' # Example data set with 2 x 3 design
#' data( "ToothGrowth" )
#' df <- ToothGrowth
#' # Outcome is length of tooth
#' df$Y <- df$len
#' # Define terms for main effects of delivery type
#' # and dose and interaction between them
#' # and standardize all predictors
#' df$PC.ZS.Type <- scale( df$supp == 'OJ' )
#' df$PC.ZS.Dose <- scale( df$dose )
#' df$PC.ZS.Type_x_dose <- df$PC.ZS.Type * df$PC.ZS.Dose
#'
#' # Fit linear model using 'brms'
#' alg <- algorithm_settings( quick_fit = T, rng_seed = 111 )
#' prs <- c(
#'   prior_string( 'normal( 18.81, 3 )', class = 'Intercept' ),
#'   prior_string( 'normal( 0.0, 0.8 )', class = 'b' ),
#'   prior_string( 'student_t( 3, 0.0, 7.6 )', class = 'sigma' )
#' )
#' est <- fit_with_brm( df, Y ~ PC.ZS.Type + PC.ZS.Dose + PC.ZS.Type_x_dose,
#'                      prior = prs, algorithm = alg )
#'
#' # Grouping factors
#' grp <- c( 'supp', 'dose' )
#' # Design matrix over grouping factors
#' dm <- aggregate( df[,5:7], df[,grp], mean )
#' # Data frame for results
#' res <- data.frame(
#'   Effect = c( 'Type', 'Dose', 'Type x Dose' ),
#'   p_value = '',
#'   stringsAsFactors = F
#' )
#'
#' # Create list specifying contrasts to
#' # compute main effect for delivery type
#' cl <- list( C1 = list( I = 1:6,
#'                        W = c(-1,1)[ (dm$supp == 'VC') + 1 ] ) )
#' tst <- omnibus_tests_for_cell_means(
#'   est, cl, grp = c( 'supp', 'dose' )
#' )
#' res$p_value[1] <- tst$char
#'
#' # Create list specifying contrasts to
#' # compute main effect for dose
#' cl <- list(
#'   C1 = list( I = 1:4,
#'              W = c(-1,-1,1,1) ),
#'   C2 = list( I = c(1:2,5:6),
#'              W = c(-1,-1,1,1) )
#' )
#' tst <- omnibus_tests_for_cell_means(
#'   est, cl, grp = c( 'supp', 'dose' )
#' )
#' res$p_value[2] <- tst$char
#'
#' # Create list specifying contrasts to
#' # compute interaction between type
#' # and dose
#' cl <- list(
#'   C1 = list( I = 1:4,
#'              W = c(1,-1,-1,1) ),
#'   C2 = list( I = c(1:2,5:6),
#'              W = c(1,-1,-1,1) )
#' )
#' tst <- omnibus_tests_for_cell_means(
#'   est, cl, grp = c( 'supp', 'dose' )
#' )
#' res$p_value[3] <- tst$char
#'
#' # Display results
#' print( res )
#'
#' # Classical approach
#' print( anova( lm( Y ~ supp*dose, data = df ) ) )
#'
#' @export

omnibus_tests_for_cell_means <- function( est,
                                       contrast_list,
                                       grp = NULL,
                                       dm = NULL,
                                       digits = 3 ) {

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

  # Create matrix for posterior samples for cell means
  CM <- matrix( NA, S, nrow( dm ) )

  # Specify design matrix
  DM <- as.matrix( dm[,sel] )
  # Add column if needed for intercept
  if ( any( grepl( 'b_Intercept', iv ) ) ) {
    DM <- cbind( 1, DM )
  }
  # Loop over cell means
  for ( i in 1:ncol( CM ) ) {

    # Loop over samples
    for ( j in 1:S ) {
      # Compute cell means
      CM[j,i] <- P[j,] %*% DM[i,]
    }

  }

  nc <- length( contrast_list )

  M <- matrix( NA, S, nc )
  for ( i in 1:nc ) {
    M[,i] <- t( colMeans( t( CM[ , contrast_list[[i]]$I ] ) *
                            contrast_list[[i]]$W ) )
  }

  out <- list(
    estimate = mean( rowMeans( M ) ),
    p_value = mean( rowSums( M > 0 ) == ncol( M ) ),
    char = '',
    cell_means = CM
  )
  if ( out$p_value > .5 ) out$p_value = 1 - out$p_value
  out$p_value = out$p_value * 2

  out$char = paste0( 'p = ', round( out$p_value, digits ) )
  if ( out$p_value == 0 )
    out$char = paste0( 'p < 0.',
                       paste(
                         c( rep( '0', digits - 1 ),'1' ),
                         collapse = '' ) )

  return( out )
}

