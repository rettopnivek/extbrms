% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S01-General_purpose_functions_for_estimation.R
\name{fit_with_brm}
\alias{fit_with_brm}
\title{Fit Model to Data Using 'brm'}
\usage{
fit_with_brm(
  data,
  formula,
  prior = NULL,
  algorithm = NULL,
  sample_prior = "yes",
  output = NULL,
  track_time = T,
  predict_probs = NULL,
  status = FALSE,
  ...
)
}
\arguments{
\item{data}{An object of class \code{data.frame} (or one that can
be coerced to that class) containing data of all variables used
in the model.}

\item{formula}{An object of class \code{\link[stats]{formula}},
\code{\link[brms]{brmsformula}}, or \code{\link[brms]{mvbrmsformula}}
(or one that can be coerced to that classes):
A symbolic description of the model to be fitted. The
details of model specification are explained in
\code{\link[brms]{brmsformula}}.}

\item{prior}{One or more \code{brmsprior} objects created by
\code{\link[brms]{set_prior}} or related functions and
combined using the \code{c} method or the \code{+} operator.
See also \code{\link[brms]{get_prior}} for more help.}

\item{algorithm}{A list specifying the number of iterations to run,
the warm-up period, the number of chains and cores to use,
the settings for the \code{adapt_delta} and \code{max_treedepth}
control options, and the seed to use for reproducibility. Output
from the \code{\link[extbrms]{algorithm_settings}} function.}

\item{sample_prior}{A character string, either code{'no'},
(no sampling from the prior distributions), \code{'yes'}
(sample from both the posterior and the prior distributions),
or \code{'only'} (sample only from the prior distribution).}

\item{output}{An optional list to which the output should be
added.}

\item{track_time}{Logical; if \code{TRUE}, tracks the run time
of the function.}

\item{predict_probs}{A vector of probabilities, input passed
on to the \code{probs} argument in \code{\link[brms]{predict.brmsfit}}.}

\item{status}{Logical; if \code{TRUE} prints brief
messages to the console to track progress for debugging
purposes.}

\item{...}{Additional arguments to the \code{\link[brms]{brm}}
function.}
}
\value{
A list including the elements...
\itemize{
\item data: The data frame used as input to the
\code{brm} function;
\item priors: The priors extracted via the
\code{prior_summary} function;
\item fit: The output from the \code{brm} function;
\item predicted: A matrix simulated data based on
the sets posterior samples (rows) per each observation
in the data set (columns);
\item algorithm: The list with the algorithm settings;
\item parameters: A list with subsets of parameter names;
\item time: How long it took to fit the model.
}
}
\description{
A convenience function that makes a call to the
\code{\link[brms]{brm}} function with the
specified user settings. Arguments are designed
to improve compatibility with packages like
'dplyr', and to streamline specification
of algorithm settings.
}
\details{
Check the help page for \code{\link[brms]{brm}} for
more details.
}
\examples{
# Example for simple linear regression
\dontrun{
data( cars )
# Setting for algorithm for quicker estimation
alg <- algorithm_settings( quick_fit = T, rng_seed = 111 )
# Fit model using 'brm'
est <- fit_with_brm( cars, dist ~ speed, algorithm = alg )
}

# Example for logistic regression
\dontrun{
data( mtcars )
# Standardize predictors
dtf <- mtcars[,c('vs','wt','disp')]
dtf$wt <- scale( dtf$wt )[,1]
dtf$disp <- scale( dtf$disp )[,1]
# Setting for algorithm for quicker estimation
alg <- algorithm_settings( quick_fit = T, rng_seed = 112 )
# Specify priors
prs <- get_prior(
  vs ~ wt + disp,
  family = bernoulli( link = 'logit' ),
  data = mtcars
)
# Empirical Bayes prior on intercept
p <- mean( dtf$vs )
lo <- log( p / ( 1 - p ) )
prs$prior[ prs$class == 'Intercept' ] <-
  paste0( 'normal( ', round( lo, 2 ), ', 2.5 )' )
# Weakly regularizing priors on coefficients
prs$prior[ prs$class == 'b' & prs$coef != '' ] <-
  'normal( 0.0, 1.0 )'
est <- fit_with_brm(
  dtf, vs ~ wt + disp,
  family = bernoulli( link = 'logit' ),
  prior = prs,
  algorithm = alg
)
}

}
