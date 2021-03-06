% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combFun.R
\name{summary.CURobustML}
\alias{summary.CURobustML}
\title{Summarize the estimates from the \code{CURobustML} routine}
\usage{
\method{summary}{CURobustML}(object, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{object}{object of class \code{CURobustML}.}

\item{digits}{integer indicating the number of decimal places}

\item{...}{additional arguments ...}
}
\description{
coefficient table containing estimates, standard errors, z-value, and p-value. If exists, bootstrap standard errors are used by default.
}
