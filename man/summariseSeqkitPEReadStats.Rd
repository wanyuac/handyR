% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioinf__summariseSeqkitPEReadStats.R
\name{summariseSeqkitPEReadStats}
\alias{summariseSeqkitPEReadStats}
\title{Summarising statistics of paired-end read sets}
\usage{
summariseSeqkitPEReadStats(
  tsv,
  header = c("file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len",
    "max_len", "Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num", "Q20_perc", "Q30_perc",
    "AvgQual", "GC", "sum_N"),
  ref_len = 5e+06,
  ext = ".fastq.gz",
  suf_R = FALSE,
  sort_by_name = "increasing",
  sort_by_depth = NULL
)
}
\arguments{
\item{tsv}{Tab-delimited output (a TSV file) from command 'seqkit stats'}

\item{header}{Column names matching seqkit stats's output}

\item{ref_len}{Length of the reference genome in base pairs. It could be the mean or median genome size of a species.}

\item{ext}{Filename extension of FASTQ files in the input TSV file. For example, '.fastq.gz' or '.fq.gz'.}

\item{suf_R}{A logical value indicating whether 'R' is used in the filename suffices. For instance, suf_R = TRUE when read files are ended with '_R1.fastq.gz' and '_R2.fastq.gz'.}

\item{sort_by_name}{A string with values "decreasing", "increasing (default), or NULL indicating whether the output data frame will be sorted by isolate names in a specific order. This argument overrides "sort_by_depth" if the former is not NULL.}

\item{sort_by_depth}{A string with values "decreasing", "increasing", or NULL (no sorting) indicating whether the output data frame will be sorted in a specific order for sequencing depths and isolate names.}
}
\value{
A data frame with summary statistics including sequencing depths.
}
\description{
This function summarises statistics from seqkit for paired-end read sets. To use this function, the input TSV file must
be generated using command 'seqkit stats --all'.
}
\author{
Yu Wan <wanyuac@gmail.com>
}
