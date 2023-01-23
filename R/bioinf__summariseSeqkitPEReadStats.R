#' @title Summarising statistics of paired-end read sets
#' @description This function summarises statistics from seqkit for paired-end read sets.
#' @param tsv Tab-delimited output (a TSV file) from command 'seqkit stats'
#' @param ref_len Length of the reference genome in base pairs
#' @param ext Filename extension of FASTQ files in the input TSV file. For example, '.fastq.gz' or '.fq.gz'.
#' @param suf_R A logical value indicating whether 'R' is used in the filename suffices. For instance, suf_R = TRUE when read files are ended with '_R1.fastq.gz' and '_R2.fastq.gz'.
#' @param srt A string with values "decreasing", "increasing", or NULL (no sorting) indicating whether the output data frame will be sorted in a specific order for sequencing depths and isolate names.
#' @return A data frame with summary statistics including sequencing depths.
#' @author Yu Wan <wanyuac@@126.com>
#' @export
# (C) Copyright 2023 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 4 Jan 2023; last update: 22 Jan 2023.

summariseSeqkitPEReadStats <- function(tsv, ref_len = 5e6, ext = ".fastq.gz", suf_R = FALSE, srt = "increasing") {
    rs <- read.delim(file = tsv, stringsAsFactors = FALSE)
    names(rs)[1] <- "File"
    rs$File <- gsub(pattern = ext, replacement = "", x = basename(rs$File), fixed = TRUE)
    trim_len <- ifelse(suf_R, 3, 2)  # suf_R = TRUE: "_R[1,2]"; otherwise, "_[1,2]"
    isolates <- unique(sapply(rs$File, function(x) substr(x, start = 1, stop = nchar(x) - trim_len), simplify = TRUE, USE.NAMES = FALSE))
    n <- length(isolates)
    out <- do.call(rbind.data.frame, lapply(isolates, .perIsolateStats, rs, ref_len, suf_R))
    m <- sum(is.na(out$Sum_len))
    if (m > 0) {
        print(paste("Warning: read sets of", m, "out of", n, "isolates could not be summarised.", sep = " "))
    } else {
        print(paste("All read sets of", n, "isolates have been summarised.", sep = " "))
    }
    if (! is.null(srt)) {
        out <- out[order(out$Depth, out$Isolate, decreasing = (srt == "decreasing")), ]
    }
    rownames(out) <- NULL
    return(out)
}

.perIsolateStats <- function(i, rs, ref_len, suf_R) {
    suf <- ifelse(suf_R, "_R", "_")
    r <- subset(rs, File %in% c(paste0(i, suf, "1"), paste0(i, suf, "2")))
    if (nrow(r) == 2) {
        r <- r[order(r$File, decreasing = FALSE), ]
        total_bases <- sum(r$sum_len)
        total_reads <- sum(r$num_seqs)
        out <- data.frame(Isolate = i, Sum_len = total_bases, Read_num = total_reads,
                          Depth = round(total_bases / ref_len, digits = 2),
                          Max_len = paste(r$max_len, collapse = ";"), N50 = paste(r$N50, collapse = ";"),
                          Avg_len = round(total_bases / total_reads, digits = 0),
                          Min_len = paste(r$min_len, collapse = ";"), Q3_len = paste(r$Q3, collapse = ";"),
                          stringsAsFactors = FALSE)
    } else {
        out <- data.frame(Isolate = i, Sum_len = NA, Read_num = NA, Depth = NA, Max_len = NA, N50 = NA, Avg_len = NA,
                          Min_len = NA, Q3_len = NA, stringsAsFactors = FALSE)
    }
    return(out)
}
