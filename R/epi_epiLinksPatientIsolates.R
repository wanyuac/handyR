#' @title Identifying epidemiological links with isolates from a specific patient
#' @description Finding out epidemiologically linked isolates within the same each genetic cluster c as isolates of a specific patient.
#' @param p A patient ID that can be found in L
#' @param c Genetic cluster of interest
#' @param L A data frame of four columns: "Isolate", "Cluster", "Patient", "Lab"; "Cluster" stands for genetic clusters, and "Lab" may represent a hospital or healthcare trust.
#' @param D A symmetric matrix of maximum days between possible collection dates of isolates
#' @param G A symmetric matrix of genetic distances
#' @param d_max maximum days to define an epidemiological link
#' @param e_max maximum days to define an episode of infection or carriage
#' @export epiLinksPatientIsolates
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 28 August 2022; last update: 29 August 2022

epiLinksPatientIsolates <- function(p, c, L, D, G, d_max = 90, e_max = 90) {
    require(igraph)

    # Stage 1: edges representing isolates from the same patient in the same episode of infection or carriage, which may collect different labs ###############
    L_pc <- subset(L, (Patient == p) & (Cluster == c))
    output <- .getEpiLinks(c, L_pc, D, G, e_max)  # Edges between isolates of genetic cluster c and from patient p regardless of labs

    # Stage 2: isolates from other patients that are epidemiologically linked to isolates from patient p ###############
    output_tmp <- data.frame(Isolate_1 = character(), Isolate_2 = character(), Cluster = character(), Divergence = integer(), Days = integer(),
                             Patient_1 = character(), Patient_2 = character(), Lab_1 = character(), Lab_2 = character(), stringsAsFactors = FALSE)
    for (l in unique(L_pc$Lab)) {  # Labs where isolates of cluster c were collected from patient p
        output_tmp <- rbind.data.frame(output_tmp, .getEpiLinks(c, subset(L, (Lab == l) & (Cluster == c)), D, G, d_max))
    }

    # Stage 3: Identify connected components that contain isolates from patient p (regardless of labs) ###############
    E <- output_tmp[, c("Isolate_1", "Isolate_2")]
    names(E) <- c("from", "to")
    net <- igraph::graph_from_data_frame(d = E, directed = FALSE)
    cc <- igraph::components(net, mode = "strong")
    for (i in 1 : cc$no) {  # Search isolates of patient p against each connected component
        isolates <- names(which(cc$membership == i))  # Isolates in the i-th component
        if (any(L_pc$Isolate %in% isolates)) {
            output <- rbind.data.frame(output, subset(output_tmp, (Isolate_1 %in% isolates) & (Isolate_2 %in% isolates)), stringsAsFactors = FALSE)
        }
    }
    output <- output[!duplicated(output[, c(1, 2)]), ]  # Duplicated rows are created when mulitiple epidemiologically linked isolates from the same patient were obtained at the same lab.

    return(output)
}

.getEpiLinks <- function(c, L, D, G, t_max = 90) {
    # This is a subordinate function of epiLinksPatientIsolates
    output <- data.frame(Isolate_1 = character(), Isolate_2 = character(), Cluster = character(), Divergence = integer(), Days = integer(),
                         Patient_1 = character(), Patient_2 = character(), Lab_1 = character(), Lab_2 = character(), stringsAsFactors = FALSE)
    n <- nrow(L)
    if (n > 1) {  # At least two isolates
        for (i in 1 : (n - 1)) {
            for (j in (i + 1) : n) {
                r_i <- L[i, ]
                r_j <- L[j, ]
                s_i <- r_i$Isolate
                s_j <- r_j$Isolate
                t <- D[s_i, s_j]  # Days between collection dates
                if (t <= t_max) {  # Same episode or transmission chain: an epidemiological link is established between these two isolates
                    output <- rbind.data.frame(output, data.frame(Isolate_1 = s_i, Isolate_2 = s_j, Cluster = c, Divergence = G[s_i, s_j], Days = t,
                                                                  Patient_1 = r_i$Patient, Patient_2 = r_j$Patient, Lab_1 = r_i$Lab, Lab_2 = r_j$Lab, stringsAsFactors = FALSE),
                                               stringsAsFactors = FALSE)
                }
            }
        }
    }  # Otherwise, no epidemiological link is possible.

    return(output)
}
