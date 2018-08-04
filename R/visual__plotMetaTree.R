#' @title Drawing a metatree (i.e. tree + data + metadata) to either the screen or a file.
#' @description This script is originally developed by Dr Kathryn Holt in the University of Melbourne.
#' The original script (plotTree.R) was acquired from https://github.com/katholt/plotTree/tree/master on 15 Sept 2015.
#' There are two important changes in my version:
#' (1) The function name was changed from plotTree to plotMetaTree to avoid the conflict with the plotTree function in the phytools package.
#' (2) Options "w" and "h" were renamed to "img.width" and "img.height" to avoid the problem:
#'     Error in switch(units, in = res, cm = res/2.54, mm = res/25.4, px = 1) *  :
#'         non-numeric argument to binary operator
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Editions: 16 Sept, 13 Oct 2015, 27 Feb 2016, 14 Nov 2016, 12 Apr 2017, 5 Aug 2018

plotMetaTree <- function(tree = NULL, heatmapData = NULL, barData = NULL,
                         infoFile = NULL, blockFile = NULL, snpFile = NULL,
                         encodedSNPMatrix = NULL, snp.pos = NULL, snp.pch = "|", bar.cex = 0.25,
                         ladderise = NULL, ancestral.reconstruction = FALSE,
                         heatmap.colours = rev(gray(seq(from = 0, to = 1, by = 0.1))),
                         heatmapBreaks = NULL, heatmapDecimalPlaces = 1, heatmap.blocks = NULL,
                         cluster = FALSE, cluster.method = "ward.D", dist.method = "euclidean",
                         treeWidth = 20, infoWidth = 10, dataWidth = 50,
                         barDataWidth = 10, blockPlotWidth = 10, edgeWidth = 1,
                         labelHeight = 10, mainHeight = 100, barDataCol = 2,
                         show.tipLabels = FALSE, tipLabelSize = 1, offset = 0,
                         tip.colour.cex = 0.5, tipColours = NULL,
                         colourNodesBy = NULL, infoCols = NULL, blwd = 5,
                         block.colour = "black", snp.colour = "red", gapChar = "?",
                         genome.size = 5E6, genome.offset = 0,
                         img.fmt = "png", img.name = "img", img.width = 640, img.height = 480, img.units = "px", img.res = 72,
                         legend = TRUE, legend.pos = "bottomleft", pie.cex = 0.5,
                         axis = FALSE, axisPos = 3, edge.color = "black",
                         infoCex = 0.8, colLabelCex = 0.8,
                         vlines.heatmap = NULL, vlines.heatmap.col = "gray75", vlines.heatmap.lty = 1,
                         hlines.heatmap = NULL, hlines.heatmap.col = "gray75", hlines.heatmap.lty = 1) {

    require(ape)

    # Prepare the tree, which may be either a tree file name or a variable name.
    if (is.character(tree)) {
        t <- read.tree(tree)
    } else {
        t <- tree
    }

    # Choose whether to ladderise the tree or not. ladderize is a function in the package ape.
    if (is.null(ladderise)) {
        tl <- t
    }
    else if (ladderise == "descending") {
        tl <-ladderize(phy = t, right = TRUE)  # specifies whether the smallest clade is on the right-hand side
    }
    else if (ladderise == "ascending") {
        tl <- ladderize(phy = t, right = FALSE)
    }
    else {
        print("Ladderisation option should be exactly 'ascending' or 'descending'.\nAny other command will raise this error.\nLeave this option empty to order branches as per input tree.")
    }

    # Get the ladderised tip order
    tips <- tl$edge[, 2]  # An edge is defined as a two-element vector from x to y. Here, "tips" actually comprises internal nodes.

    # Extract the order of tips from the ladderised tree, excluding internal nodes
    tip.order <- tips[tips <= length(tl$tip.label)]  # retain only labels that are not more than the sample size - they represent samples (tips) on the tree
    tip.label.order <- tl$tip.label[tip.order] # for ordering data. note that for tiplabel(), the order is the same as in t$tip (= tl$tip)

    # Prepare the heat map
    if (!is.null(heatmapData)) {
        # read a heatmap data set and convert it into to data frame
        x <- .readMatrix(heatmapData)

        # match rows of the heatmap matrix with tips of the tree
        # Variables x and the tree must have the same number of taxa.
        y.ordered <- x[tip.label.order, ]  # Match taxa on the heat map with those on the tree. y.ordered is a matrix ordered by rows.

        # whether to perform clustering on columns or not?
        if (cluster) {
            if (cluster.method == "square" & ncol(y.ordered) == nrow(y.ordered)) {
                # reorder columns to follow the row order in the matrix for the heat map
                original.order <- 1 : nrow(x)
                names(original.order) <- rownames(x)
                reordered <- original.order[tip.label.order]
                y.ordered <- y.ordered[, rev(as.numeric(reordered))]
            } else {
                h <- hclust(d = dist(x = t(na.omit(y.ordered)), method = dist.method), method = cluster.method)
                y.ordered <- y.ordered[, h$order]  # recorder columns (genes) of the heatmap matrix
                # The order of rows (tips) remain the unchanged.
            }
        }
    }

    # Prepare the bar plot
    if (!is.null(barData)) {
        b <- .readMatrix(barData)
        barData <- b[, 1]
        names(barData) <- rownames(b)
    }

    # Prepare sample information
    if (!is.null(infoFile)) {
        info <- .readMatrix(infoFile)
        info.ordered <- info[rev(tip.label.order),]
    } else {
        info.ordered <- NULL
    }

    # Colouring tips by sample information
    # Specify a categorical trait in the info matrix to colour nodes and infer ancestral states
    ancestral <- NULL
    nodeColourSuccess <- NULL
    if (!is.null(colourNodesBy) & !is.null(infoFile)) {
        if (colourNodesBy %in% colnames(info.ordered)) {
            nodeColourSuccess <- TRUE
            loc1 <- info.ordered[, which(colnames(info.ordered) == colourNodesBy)]

            # assign values
            tipLabelSet <- character(length(loc1))
            names(tipLabelSet) <- rownames(info.ordered)
            groups <- table(loc1, exclude = "")  # find out categories in the subject column of the info matrix
            n <- length(groups)  # number of categories
            groupNames<-names(groups)

            # set colours
            colours <- ifelse(is.null(tipColours), rainbow(n), tipColours)

            # assign colours based on values
            for (i in 1 : n) {
                g <- groupNames[i]
                tipLabelSet[loc1 == g] <- colours[i]
            }
            tipLabelSet <- tipLabelSet[tl$tip]

            # ancestral reconstruction
            if (ancestral.reconstruction) {
                ancestral <- ace(loc1, tl, type = "discrete")  # a function in the ape package for ancestral character estimation
            }
        }
    }

    # Initialise an external device for plotting
    if (img.fmt == "png") {
        png(filename = paste(img.name, img.fmt, sep = "."), width = img.width, height = img.height,
            units = img.units, res = img.res, type = "cairo-png")
    } else {
        img.fmt <- "pdf"
        pdf(filename = paste(img.name, img.fmt, sep = "."), width = img.width, height = img.height,
            units = img.units, res = img.res)
    }

    # Set up a plotting layout
    doBlocks <- !is.null(blockFile) | !is.null(snpFile) | !is.null(encodedSNPMatrix)

    l <- .getLayout(infoFile = infoFile, heatmapData = heatmapData, barData = barData, doBlocks = doBlocks,
                    treeWidth = treeWidth, infoWidth = infoWidth, dataWidth = dataWidth, edgeWidth = edgeWidth,
                    labelHeight = labelHeight, mainHeight = mainHeight, barDataWidth = barDataWidth,
                    blockPlotWidth = blockPlotWidth)
    layout(l$m, widths = l$w, heights = l$h)  # specifying complex plot arrangements for the current image device
    print(l$m)

    # Draw the tree
    par(mar = c(0, 0, 0, 0))  # no margin around the tree
    tlp <- plot.phylo(tl, no.margin = TRUE, show.tip.label = show.tipLabels, label.offset = offset,
                      edge.width = edgeWidth, edge.color = edge.color, xaxs = "i", yaxs = "i",
                      y.lim = c(0, length(tl$tip) + 0.5), cex = tipLabelSize)  # The distance between two neighbouring taxa equals 1.

    # colour tips by the selected categorical trait
    if (!is.null(nodeColourSuccess)) {
        tiplabels(col= tipLabelSet, pch = 16, cex = tip.colour.cex)

        if (ancestral.reconstruction) {
            nodelabels(pie = ancestral$lik.anc, cex = pie.cex, piecol = colours)
        }

        if (legend) {
            legend(legend.pos, legend = groupNames, fill = colours)
        }
    }

    if (axis) {
        axisPhylo(axisPos)
    }

    # Plot sample information
    if (!is.null(infoFile)) { # if info is provided
        par(mar = rep(0, 4))

        if (!is.null(infoCols)) {
            infoColNumbers <- which(colnames(info.ordered) %in% infoCols)
        } else {
            infoColNumbers <- 1:ncol(info.ordered)
        }

        plot(NA, axes = FALSE, pch = "",
             xlim = c(0, length(infoColNumbers) + 1.5), ylim = c(0.5, length(tl$tip) + 0.5),
             xaxs = "i", yaxs = "i")

        # plot all info columns
        for (i in 1 : length(infoColNumbers)) {
            j <- infoColNumbers[i]
            text(x = rep(i + 1, nrow(info.ordered) + 1), y = c((nrow(info.ordered)) : 1),
                 info.ordered[, j], cex = infoCex)
        }
    }

    # Draw the heat map
    if (!is.null(heatmapData)) {
        # set levels to of values to which colours will be designated
        # Colours scale across the whole range of values in the y.ordered matrix, not across single columns.
        if (is.null(heatmapBreaks)) {
            heatmapBreaks <- seq(min(y.ordered, na.rm = TRUE), max(y.ordered, na.rm = TRUE),
                                 length.out = length(heatmap.colours) + 1)  # the minimum and maximum of the y.ordered matrix
        }  # heatmapBreaks specifies boundaries between colour segments.

        # plot heatmap
        par(mar = rep(0, 4), xpd = TRUE)
        # nrow(t(y.ordered)) = ncol(y.ordered), which equals the number of variables for every taxon.
        # Hence, both x and y specify midpoints between coloured cells. Namely, colour boxes are X-centred and Y -centred at
        # middle points. Each box has a width and height of 1.
        image(x = (1 : ncol(y.ordered)) - 0.5, y = (1 : nrow(y.ordered)) - 0.5, z = as.matrix(t(y.ordered)),
              col = heatmap.colours, breaks = heatmapBreaks, axes = FALSE, xaxs = "i", yaxs = "i", xlab = "", ylab = "")  # heatmapBreaks works for z values
        print("Image has drawn")
        # draw vertical lines between columns on the heat map
        # Note that the left-most of the heatmap is v = 0 and the right-most of it is v = ncol(matrix).
        if (!is.null(vlines.heatmap)) {
            for (line.v in vlines.heatmap) {
                abline(v = line.v, col = vlines.heatmap.col, lty = vlines.heatmap.lty)
            }
        }

        # draw horizontal lines defining rows over heatmap
        # Note that the bottom of the heatmap is h = 0 and the top of it is h = nrow(matrix).
        if (!is.null(hlines.heatmap)) {
            for (line.h in hlines.heatmap) {
                abline(h = line.h, col = hlines.heatmap.col, lty = hlines.heatmap.lty)
            }
        }

        # overlay blocks on the heatmap
        if (!is.null(heatmap.blocks)) {
            for (coords in heatmap.blocks) {
                rect(xleft = coords[1], 0, coords[2], ncol(y.ordered), col = vlines.heatmap.col, border = NA)
            }
        }

        # draw labels on the heat map
        par(mar = rep(0, 4))
        plot(NA, axes = FALSE, xaxs = "i", yaxs = "i", ylim = c(0, 2), xlim = c(0.5, ncol(y.ordered) + 0.5))
        text(1 : ncol(y.ordered) - 0.5, rep(0, ncol(x)), colnames(y.ordered), srt = 90, cex = colLabelCex, pos = 4)

        # Finally, add a bar to illustrate the colour scale of the heatmap at the top-left corner
        if (img.height < 5000) {
            par(mar = c(5, 0, 5, 2))  # otherwise, an error occurs: "figure margins too large"
        } else {
            print("Set margins")
            par(mar = c(10, 0, 30, 2))
        }

        # Function as.matrix(heatmapBreaks) returns an n-by-1 matrix, where n is the number of breaks.
        image(z = as.matrix(heatmapBreaks), breaks = heatmapBreaks, col = heatmap.colours,
              yaxt = "n", axes = FALSE)  # By default, boundaries range from 0 to 1 on the X axis.
        axis(1, at = seq(0, 1, length.out = length(heatmapBreaks)),
             labels = round(heatmapBreaks, digits = heatmapDecimalPlaces))
    }

    # bar plot
    if (!is.null(barData)) {
        par(mar = rep(0, 4))
        barplot(barData[tip.label.order], horiz = TRUE, axes = FALSE, xaxs = "i",
                yaxs = "i", xlab = "", ylab = "", ylim = c(0.25, length(barData) + 0.25),
                xlim = c((-1) * max(barData, na.rm = TRUE) / 20, max(barData, na.rm = TRUE)),
                col = barDataCol, border = 0, width = 0.5, space = 1, names.arg = NA)

        # scale for barData plot
        par(mar = c(2, 0, 0, 0))
        plot(NA, yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "",
             ylim = c(0,2), xlim = c((-1) * max(barData, na.rm = TRUE) / 20, max(barData, na.rm = TRUE)),
             frame.plot = FALSE)
    }

    # draw SNP and recombination blocks
    if (doBlocks) {
        # According to xlim, the range of data points on the X axis is the genome size + 1.5 when
        # genome.offset = 0. SNPs are plotted at their genomic positions and hence reveal their
        # actual distribution on the reference genome.

        par(mar = rep(0, 4))
        plot(NA, axes = FALSE, pch = "",
             xlim = c(genome.offset, genome.offset + genome.size + 1.5),
             ylim = c(0.5, length(tl$tip) + 0.5), xaxs = "i", yaxs = "i")  # create an empty plotting area

        # Plot snps as vertical bars
        # Expected format of the CSV file:
        #   First column: SNP positions (integers). This column will become row names of the data frame "snps".
        #   Excluding the 1st column, namely, contents of "snps".
        #       First row: strain names, where the first one is the reference strain.
        #       Everyone of rest of rows: allele (nucleotides or gaps) of a SNP
        # Assuming L SNPs found in N strains, the dimension of "snps" should be (L+1)-by-N.
        if (!is.null(snpFile)) {
            encodedSNPMatrix <- NULL  # The option snpFile overrides encodedSNPMatrix to avoid plotting SNP information twice.

            # read SNP genotypes
            snps <- read.csv(snpFile, header = FALSE , row.names = 1, stringsAsFactors = FALSE)  # in case colnames start with numbers or contain dashes, which R does not like as column headers
            snps.strainCols <- snps[1, ]  # values of the first row, or column names in the CSV file: strain names
            snps <- snps[-1, ]  # rest of rows: a matrix of every SNP
            ref.alleles <- snps[, 1]  # alleles across L loci in the reference genome

            for (strain in tip.label.order){
                # print SNPs compared to reference (ancestral) alleles in column 1
                query.alleles <- snps[, which(snps.strainCols == strain)]  # alleles of L loci in the query genome

                # identify SNPs which are neither all the same to the reference nor are gaps
                s <- rownames(snps)[ref.alleles != query.alleles & query.alleles != gapChar & ref.alleles != gapChar]
                y <- which(tip.label.order == strain)

                if (length(s) > 0) {  # when there are different alleles present in the query genome comparing to the reference genome
                    for (x in s) {
                        points(x, y, pch = "|", col = snp.colour, cex = bar.cex)
                    }
                }
            }
        }

        # Draw a SNP matrix of encoded genotypes (e.g., minor allele = 1 and major allele = 0)
        # This feature only displays positions of SNPs, regardless how many alleles are present at each SNP (They will all be coloured the same).
        # Structure of the encodedSNPMatrix: n(sample)-by-n(SNPs). It must be a matrix.
        #   rownames: sample names; colnames: SNP names.
        # snp.pos is a named vector of integers (by SNP names) that must be supplied.
        if (is.matrix(encodedSNPMatrix)) {
            if (is.null(snp.pos)) {
                print("Error: A named integer vector of SNP positions must be provided.")
                dev.off()
                stop()
            } else {
                strains <- rownames(encodedSNPMatrix)
                for (snp in colnames(encodedSNPMatrix)) {
                    strains.2plot <- strains[encodedSNPMatrix[, snp] != 0]  # strains having the current minor allele
                    pos <- snp.pos[[snp]]
                    for (s in strains.2plot) {
                        points(x = pos, y = which(tip.label.order == s), pch = snp.pch, col = snp.colour, cex = bar.cex)  # "s" must present on the tree
                    }
                }
            }
        }

        # plot blocks
        if (!is.null(blockFile)){
            blocks <- read.delim(blockFile, header = FALSE)
            for (i in 1 : nrow(blocks)) {
                if (as.character(blocks[i, 1]) %in% tip.label.order) {
                    y <- which(tip.label.order == as.character(blocks[i, 1]))
                    x1 <- blocks[i, 2]
                    x2 <- blocks[i, 3]
                    lines(c(x1, x2), c(y, y), lwd = blwd, lend = 2, col = block.colour)
                }
            }
        }

    }

    dev.off()  # close the plotting device

    # Return ordered info and ancestral reconstruction object
    if (!is.null(heatmapData)){
        mat <- as.matrix(y.ordered)
        mat <- mat[nrow(mat) : 1, ]  # reverse the order of rows to match the phylogeny as illustrated in the plot
    }
    else {
        mat <- NULL
    }

    return(list(OTUs = tip.label.order, mat = mat, info = info.ordered, anc = ancestral))
}

.readMatrix <- function(m) {
    # Read a CSV file when a matrix or data frame is not provided.
    # The CSV file should have at least a column for row names, namely, taxon names.
    # m: a name of a matrix/data frame, or a name of a CSV file

    if (is.matrix(m)) {
        x <- data.frame(m)
    }
    else if (is.data.frame(m)) {
        x <- m
    }
    else {
        x <- read.csv(m, row.names = 1)  # take the first column as row names
    }

    return(x)
}

.getLayout <- function(infoFile, heatmapData, barData, doBlocks,
                      treeWidth = 10, infoWidth = 10, dataWidth = 30,
                      edgeWidth = 1, barDataWidth = 10, blockPlotWidth = 10,
                      labelHeight = 10, mainHeight = 100) {

    # tree
    w <- c(edgeWidth, treeWidth)
    m <- cbind(c(0, 0, 0), c(0, 1, 0)) # first two columns, edge + tree
    x <- 1

    # info
    if (!is.null(infoFile)) {  # when info is provided
        x <- x + 1
        m <- cbind(m, c(0, x, 0))
        w <- c(w, infoWidth)
    }

    # heatmap
    if (!is.null(heatmapData)) {
        x <- x + 1
        m <- cbind(m, c(x + 1, x, 0)) # add heatmap & labels
        x <- x + 2
        m[1, 2] <- x # add heatmap scale above tree
        w <- c(w, dataWidth)
    }

    # barplot
    if (!is.null(barData)) {
        x <- x + 1
        m <- cbind(m, c(0, x, x + 1)) # barplot and scale bar
        x <- x + 1
        w <- c(w, barDataWidth)
    }

    if (doBlocks) {
        x <- x + 1
        m <- cbind(m, c(0, x, 0)) # recomb blocks
        w <- c(w, blockPlotWidth)
    }

    # empty edge column
    m <- cbind(m, c(0, 0, 0))
    w <- c(w, edgeWidth)

    if (!is.null(heatmapData) | !is.null(barData)) {
        h <- c(labelHeight, mainHeight, labelHeight)
    }
    else {
        h <- c(edgeWidth, mainHeight, edgeWidth)
    }

    return(list(m = as.matrix(m), w = w, h = h))
}
