% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.rpair.R
\name{plot.rpair}
\alias{plot.rpair}
\title{Plot coefficients from a "rpair" object}
\usage{
\method{plot}{rpair}(x, xvar = c("norm", "lambda", "dev"), legend = FALSE)
}
\arguments{
\item{x}{An "rpair" object.}

\item{xvar}{Value on the x-axis. "norm" plots against the L1-norm of the coefficients, "lambda"
against the log-lambda sequence, and "dev" against the percent deviance explained.}

\item{legend}{Whether to include a legend. Default: FALSE.}

\item{label}{Whether to include variable labels. Default: FALSE.}
}
\value{
ggplot object
}
\description{
Produces a coefficient profile plot of the coefficient paths for a fitted "rpair" object. Supports "norm",
"lambda" and "dev" plot types. Plots are created as ggplot objects.
}
\examples{
efit = rpair_gloss(rpair::ds1_x, rpair::ds1_y, standardize = F, pmax = 50, loss_type = "exp")
plot(efit, xvar="norm")
plot(efit, xvar="lambda", label=T)

}
\author{
KC
}
