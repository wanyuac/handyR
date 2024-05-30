#' @title Summarising statistics of paired-end read sets
#' @description This function summarises statistics from seqkit for paired-end read sets. To use this function, the input TSV file must
#' be generated using command 'seqkit stats --all'.
#' @param tsv Tab-delimited output (a TSV file) from command 'seqkit stats'
#' @param header Column names matching seqkit stats's output
#' @param ref_len Length of the reference genome in base pairs
#' @param ext Filename extension of FASTQ files in the input TSV file. For example, '.fastq.gz' or '.fq.gz'.
#' @param suf_R A logical value indicating whether 'R' is used in the filename suffices. For instance, suf_R = TRUE when read files are ended with '_R1.fastq.gz' and '_R2.fastq.gz'.
#' @param sort_by_name A string with values "decreasing", "increasing (default), or NULL indicating whether the output data frame will be sorted by isolate names in a specific order. This argument overrides "sort_by_depth" if the former is not NULL.
#' @param sort_by_depth A string with values "decreasing", "increasing", or NULL (no sorting) indicating whether the output data frame will be sorted in a specific order for sequencing depths and isolate names.
#' @return A data frame with summary statistics including sequencing depths.
#' @author Yu Wan <wanyuac@@gmail.com>
#' @export
# (C) Copyright 2023 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Release: 4 Jan 2023; last update: 30 May 2024.

summariseSeqkitPEReadStats <- function(tsv,
                                       header = c("file", "format", "type", "num_seqs", "sum_len", "min_len",
                                                  "avg_len", "max_len", "Q1", "Q2", "Q3", "sum_gap", "N50",
                                                  "N50_num", "Q20_perc", "Q30_perc", "AvgQual", "GC"),
                                       ref_len = 5e6,
                                       ext = ".fastq.gz",
                                       suf_R = FALSE,
                                       sort_by_name = "increasing",
                                       sort_by_depth = NULL) {
    rs <- read.delim(file = tsv, stringsAsFactors = FALSE)
    header[1] <- "file"  # Fix the name of the first column for the ease of coding
    names(rs) <- header
    rs$file <- gsub(pattern = ext, replacement = "", x = basename(rs$file), fixed = TRUE)  # Remove filename extensions (e.g., ".fastq.gz")
    trim_len <- ifelse(suf_R, 3, 2)  # "_R[1,2]", then set suf_R = TRUE; otherwise ("_[1,2]"), set suf_R = FALSE (default)
    rs <- rs[order(rs$file, decreasing = FALSE), ]  # Ensure _1 always goes before _2 for the same isolate
    rs$file <- sapply(rs$file, function(x) substr(x, start = 1, stop = nchar(x) - trim_len), simplify = TRUE, USE.NAMES = FALSE)  # Remove "_1", "_2" or "_R1", "_R2" from column "file"
    isolates <- unique(rs$file)  # Isolate names
    n <- length(isolates)  # Number of isolates
    out <- do.call(rbind.data.frame, lapply(isolates, .perIsolateStats, rs, ref_len))
    m <- sum(is.na(out$Sum_len))
    if (m > 0) {
        print(paste("Warning: read sets of", m, "out of", n, "isolates could not be summarised.", sep = " "))
    } else {
        print(paste("All read sets of", n, "isolates have been summarised.", sep = " "))
    }
    if (! is.null(sort_by_name)) {
        out <- out[order(out$Isolate, decreasing = (sort_by_name == "decreasing")), ]
    } else if (! is.null(sort_by_depth)) {
        out <- out[order(out$Depth, out$Isolate, decreasing = (sort_by_depth == "decreasing")), ]
    }
    rownames(out) <- NULL
    return(out)
}

.perIsolateStats <- function(i, rs, ref_len) {
    r <- subset(rs, file == i)  # Returns two rows for paired-end reads
    if (nrow(r) == 2) {
        total_bases <- sum(r$sum_len)
        total_reads <- sum(r$num_seqs)
        # Note on 9 Apr 2024: removed 'N50 = paste(r$N50, collapse = ";")' and 'Q3_len = paste(r$Q3, collapse = ";")' as they have become optional in new seqkit versions.
        out <- data.frame(Isolate = i,
                          Read_num = total_reads,
                          Sum_len = total_bases,
                          Max_len = paste(r$max_len, collapse = ";"),
                          Q3_len = paste(r$Q3, collapse = ";"),
                          Avg_len = round(total_bases / total_reads, digits = 0),
                          Min_len = paste(r$min_len, collapse = ";"),
                          Depth = round(total_bases / ref_len, digits = 2),
                          Mean_GC = round((r$sum_len[1] * r$GC[1] + r$sum_len[2] * r$GC[2]) / sum(r$sum_len), digits = 2),
                          stringsAsFactors = FALSE)
    } else {
        out <- data.frame(Isolate = i,
                          Read_num = NA,
                          Sum_len = NA,
                          Max_len = NA,
                          Q3_len = NA,
                          Avg_len = NA,
                          Min_len = NA,
                          Depth = NA,
                          Mean_GC = NA,
                          stringsAsFactors = FALSE)
    }
    return(out)
}
