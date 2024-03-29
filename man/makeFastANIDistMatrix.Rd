% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylo__makeFastANIDistMatrix.R
\name{makeFastANIDistMatrix}
\alias{makeFastANIDistMatrix}
\title{Make a distance matrix from FastANI output (tab-delimited file)}
\usage{
makeFastANIDistMatrix(f, keep_asym = FALSE, frac = FALSE, suffix = "fasta")
}
\arguments{
\item{f}{Path to the tab-delimited output file of FastANI}

\item{keep_asym}{A logical flag specifying whether to keep the original asymmetric distance matrix.}

\item{frac}{A logical flag specifying whether to convert percentages to decimal fractions. This option
does not affect the tree topology as the change is proportional.}

\item{suffix}{Filename extension of input FASTA files for FastANI, such as fna and fasta.}
}
\value{
One or two n-by-n distance matrices (depending on keep_asm), where n denotes the number of genomes.
}
\description{
This function coverts a FastANI output into a distance matrix.
}
\author{
Yu Wan (\email{wanyuac@126.com})
}
