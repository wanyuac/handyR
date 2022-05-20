#' @title Calculating sequencing depths from the output of 'seqkit stats'
#' @description This function calculates the sequencing depth of each genome using the output of
#' command 'seqkit stats' (https://github.com/shenwei356/seqkit)
#' @param reads A data frame imported from the output file of 'seqkit stats'. The first column of
#' input file names must not contain paths. For example, sample_1.fastq.gz is valid whereas
#' reads/filtered/sample_1.fastq.gz is invalid. Each FASTQ file should follow the pattern *_{1, 2}.fastq.gz.
#' @param len Length (bp) of the reference genome
#' @param srt Whether to sort rows by sequencing depths (default: TRUE)
#' @returns A data frame of isolate names and sequencing dapths
#' @export
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 20 May 2022

seqDepth <- function(reads, len, srt = TRUE) {
    # d: a data frame of three columns: isolate/genome, num_seqs, sum_len
    # len: whole-genome length of the reference strain (including all its plasmids, if they exist)
    # srt: whether sort rows in a descending order of sequencing depths (FALSE: sort by genome names)
    # Return: a data frame of total read lengths and fold coverages across isolates
    # Example: estimateSequencingDepth(reads[, c("file", "num_seqs", "sum_len", "min_len", "max_len")], len = 2.5e6, srt = TRUE)
    reads <- reads[, c("file", "num_seqs", "sum_len", "min_len", "max_len")]
    names(reads) <- c("Genome", "Num_reads", "Sum_len", "Min_len", "Max_len")
    reads$Genome <- sapply(reads$Genome, function(x) gsub(pattern = "_1.fastq.gz", replacement = "", x = x, fixed = TRUE), simplify = TRUE, USE.NAMES = FALSE)
    reads$Genome <- sapply(reads$Genome, function(x) gsub(pattern = "_2.fastq.gz", replacement = "", x = x, fixed = TRUE), simplify = TRUE, USE.NAMES = FALSE)
    genomes <- unique(reads$Genome)
    fc <- data.frame(Genome = character(0), Num_reads = integer(0), Sum_bp = integer(0), Min_R1 = integer(0),
                     Min_R2 = integer(0), Max_R1 = integer(0), Max_R2 = integer(0), stringsAsFactors = FALSE)
    for (g in genomes) {
        rs <- subset(reads, Genome == g)  # Two rows
        fc <- rbind.data.frame(fc, data.frame(Genome = g, Num_reads = sum(rs$Num_reads), Sum_bp = sum(rs$Sum_len),
                                              Min_R1 = rs$Min_len[1], Min_R2 = rs$Min_len[2], Max_R1 = rs$Max_len[1],
                                              Max_R2 = rs$Max_len[2], stringsAsFactors = FALSE))
    }
    fc$Avg_len <- round(fc$Sum_bp / fc$Num_reads, digits = 2)  # Average read length
    fc$Depth <- round(fc$Sum_bp / len, digits = 2)
    if (srt) {
        fc <- fc[order(fc$Depth, fc$Avg_len, decreasing = TRUE), ]
    } else {
        fc <- fc[order(fc$Genome, decreasing = FALSE), ]
    }
    return(fc[, c("Genome", "Depth", "Sum_bp", "Num_reads", "Avg_len", "Min_R1", "Min_R2", "Max_R1", "Max_R2")])
}
