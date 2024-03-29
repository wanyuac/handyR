% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumuCurve.R
\name{cumuCurve}
\alias{cumuCurve}
\title{Calculating cumulative proportions of objects in a sample population}
\usage{
cumuCurve(obj, pop)
}
\arguments{
\item{obj}{The object whose cumulative curve is of our interest.}

\item{pop}{The population from which individual objects are sampled.
Initially, this function is designed for processing counts in the results of the function hist(...). In this scenario,
both obj and pop are counts.}
}
\description{
This function calculates cumulative proportions of objects in a sample population.
There are two input numeric vectors, and both vectors are binned using the same width (for example, binned by 1000).
Therefore, all inputs are counts.
}
\author{
Yu Wan (\email{wanyuac@gmail.com})
}
