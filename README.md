# handyR  
This repository provides R functions for statistical analysis. It was called R_toolkit and built in 2015. I converted it into a package on 4 Aug 2018 to simplify its use. I also renamed it to make it easier to type and remember.

## A list of scripts
* [makeNamedList.R](#makeNamedList)
* [df2matrix.R](#df2matrix)
* [matrix2df.R](#matrix2df)
* [maskMatrix.R](#maskMatrix)
* [matrixCutoff.R](#matrixCutoff)
* [maxRange.R](#maxRange)
* [mkContingencyTable.R](#mkContingencyTable)
* [cumuCurve.R](#cumuCurve)
* [subsetDfCon.R](#subsetDfCon)

Data visualisation
* [plotMetaTree.R](#plotMetaTree)
* [scatterPlot.R](#scatterPlot)
* [Fishers\_exact\_test.R](#FishersExactTest)
* [buildNetwork.R](#buildNetwork)
* [labelGenerator.R](#labelGenerator)

## Manual
### <a name="makeNamedList"></a>makeNamedList.R
```R
# definition
makeNamedList(names.1 = NA, names.2 = NA)

# example usage
x <- makeNamedList(names.1 = c("car", "horse"), names.2 = c("price", "speed"))
y <- makeNamedList(names.1 = c("car", "horse"))
z <- makeNamedList(names.2 = c("car", "horse"))
```
This function returns a list of one or two levels with names for each level. For example, the variable x in the command above becomes a two-level named list where an element of it can be referred to as x[["car"]][["price"]] or x[["horse"]][["speed"]] and so forth. By contrast, a normal named list is returned if either names.1 or names.2 is left blank. Therefore, the list y equals to z in the aforementioned example, and they are both a one-level list (y[["car"]] and y[["horse"]], etc).  

### <a name="df2matrix"></a>df2matrix.R
```R
df2matrix(df, diag = 1, replace.na = 1)
```
This function coverts a data frame into a symmetric matrix. The data frame consists of three columns: name1, name2, value. The union of name1 and name2 will become row names and column names in the matrix, with the value goes to cells (name1, name2) and (name2, name1). Note that currently this function does not work for situations where name1 = name2. Instead, the diagnal values will be set with the argument "diag". NA values in the data frame will be replaced with the argument "replace.na".

### <a name="matrix2df"></a>matrix2df.R
```R
matrix2df(m, diag = FALSE, small = TRUE)
```
This is an inverse of the function [df2matrix](#df2matrix). It converts a symmetric matrix into a data frame of three columns: V1, V2 and value. Values of V1 and V2 are row names and column names of the input matrix. The argument "diag" determines whether the diagonal values should be transfered to the data frame (in this case, V1 = V2).

A list with two vectors will be returned if the option for small-matrix mode (small) is turned off (set to FALSE). In the large-matrix mode, which is a memory-efficient way for processing an extremely large matrix. Values in the matrix are transferred to a single vector, where row names or column names are discarded. The second vector is included in the list in order to help users to access row names and column names. This vector records the index of a value in the first vector, where a new row in the matrix starts being read.

### <a name="maskMatrix"></a>maskMatrix.R
```R
maskMatrix(m, df, default = 0, keep.diag = TRUE) 
```
This function generates a matrix comprised of only 0 or 1 in accordance with an input data frame so that you can filter values in a target matrix (m) using matrix multiplexing. The data frame records pairs of item names that you want to keep. By default, values in the input matrix that are not recorded in the data frame will be replaced by naughts. However, you can change the default factor at the argument "default". The arugment "keep.diag" determines whether diagonal values in the matrix should be kept intact.

### <a name="matrixCutoff"></a>matrixCutoff.R
```R
matrixCutoff(m, cutoff, comp = ">", new.val = 0)
```
This function resets a value to a new value (new.val) in an input matrix (m) if it is smaller/greater than the threshold. The direction of comparison is termined by the argument "comp".

### <a name="maxRange"></a>maxRange.R
```R
maxRange(v1, v2)
```
This function returns the maximum range covering both numeric vectors v1 and v2, even though they are non-overlapping.

### <a name="mkContingencyTable"></a>mkContingencyTable.R
```R
contingencyTable(a, b, c, d)
```
This programme constructs a two-by-two contingency table. a: the number of joint presence; d: the number of joint absence.

### <a name="cumuCurve"></a>cumuCurve.R
This function calculates cumulative proportions of objects in a sample population. There are two input numeric vectors, and both vectors are binned using the same width (for example, binned by 1000):  
   1. obj: the object whose cumulative curve is of our interest.
   2. pop: the population from which individual objects are sampled.  

Initially, this function is designed for processing counts in the results of the function hist(...). In this scenario, # both obj and pop are counts.

### <a name="subsetDfCon"></a>subsetDfCon.R
This function take a subset of a data frame in accordance with one or more conditions provided in the second data frame.  
Inputs:  
   1. subject: the data frame to be subset
   2. conditions: a data frame whose columns provide conditions for selecting rows in the subject.
   3. col.subject: column names or indices in the subject data frame to be searched for
   4. col.con: column names or indices in the condition data frame to be used as conditions for searching rows in the subject data frame

Note that col.subject and col.con must be matched in terms of the order of comparisons.  
Example usage:
```R
# x$a1 will be searched for values in y$a1, and so forth:
d <- subsetDfCon(subject = x, df.con = y, col.subject = c("a1", "a2), col.con = c("a1", "b1"))
# x[, 2] will be searched for values in y[, 3], for example:
d <- subsetDfCon(subject = x, df.con = y, col.subject = c(1, 2, 3), col.con = c(1, 3, 5))
```

### <a name="plotMetaTree"></a>plotMetaTree.R

```R
x <- plotMetaTree(tree, ladderise = NULL, heatmapData = NULL, barData = NULL, infoFile = NULL,
                   blockFile = NULL, snpFile = NULL, gapChar = "?", genome_size = 5E6,
                   blwd = 5, block_colour = "black", snp_colour="red", genome_offset = 0,
                   colourNodesBy = NULL, infoCols = NULL,
                   outputPNG = NULL, outputPDF = NULL, img.width = 640, img.height = 480,
                   heatmap.colours = rev(gray(seq(0,1,0.1))),
                   tip.labels = FALSE, tipLabelSize = 1, offset = 0, tip.colour.cex = 0.5, tipColours = NULL, lwd = 1.5,
                   legend = TRUE, legend.pos = "bottomleft", ancestral.reconstruction = FALSE,
                   cluster = FALSE, cluster.method = "ward.D", dist.method = "euclidean",
                   axis = FALSE, axisPos = 3, edge.color = "black", infoCex = 0.8, colLabelCex=0.8,
                   treeWidth = 10, infoWidth = 10, dataWidth = 30, edgeWidth = 1,
                   labelHeight = 10, mainHeight = 100, barDataWidth = 10, blockPlotWidth = 10,
                   barDataCol = 2, heatmapBreaks = NULL, heatmapDecimalPlaces = 1,
                   heatmap.blocks = NULL, pie.cex = 0.5,
                   vlines.heatmap = NULL, vlines.heatmap.col = "gray50", vlines.heatmap.lty = 1,
                   hlines.heatmap = NULL, hlines.heatmap.col = "gray50", hlines.heatmap.lty = 1)
```

This function draws a phylogenetic tree and metadata into a single diagram. To be more specific, it combines a phylogenetic tree, a heat map, a bar plot and a textual annotations in parallel into a single plot. This script is a modification of the script [plotTree.R](https://github.com/katholt/plotTree) developed by Dr Kathryn Holt from the University of Melbourne. However, there are two important changes in my script:
* The function name was changed from plotTree to plotMetaTree in order to avoid the conflict with the plotTree function in the phytools package.
* Options "w" and "h" were renamed to "img.width" and "img.height" to avoid the problem that: *Error in switch(units, \`in\` = res, cm = res/2.54, mm = res/25.4, px = 1) \*: non-numeric argument to binary operator\*  

Note that it would be better to use "binary" as the dist.method for binary data.  

The function returns a list of four elements:
* OTUs: tip names from the phylogenetic tree as ordered in the plot
* mat: the matrix of data drawn besides the tree
* info: information of tips following the same order as shown in the plot
* anc: ancestral nodes from the tree

### <a name="scatterPlot"></a>scatterPlot.R

```R
scatterPlot(x, y, z, x.bg = NA, y.bg = NA, z.bg = NA, log = "both", log.base = 10, replace.zero = 1e-5,
            breaks.min = 0, breaks.max = 1, breaks.n = 51, col.grad = c("blue", "green", "red"),
            log.z = FALSE, log.z.base = 10,
            output = "plot.png", w = 800, h = 640,
            title = "Y~X", x.lab = "X", y.lab = "Y", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
            x.extension = c(0, 0), y.extension = c(0, 0), margin = c(4.5, 4.5, 3, 2),
            joint.col = FALSE, col.bg = "grey", pch.bg = 1, cex.bg = 1, 
            pch = 19, cex = 1,
            show.range = FALSE, x.range = NA, x.bg.range = NA, range.pch = ".", cex.range = 1,
            add.line = FALSE, line.a = NULL, line.b = NULL, line.h = NULL, line.v = NULL, line.style = 3, line.col = "grey",
            turnoff.dev = TRUE)
```

This function draws a scatter plot of x against y, and colour coded by z.

### <a name="FishersExactTest"></a>Fishers\_exact\_test.R

```R
plotFishersP(r1, c1, n, alpha, panels = 3, prefix = "FishersProbability_", font.size = 1.5, w = 1200, h = 1000)

# example usage
p <- plotFishersP(r1 = 80, c1 = 200, n = 300, alpha = 0.05, panels = 3, prefix = "FishersProbability_3panels_", font.size = 2)

```
This function simulates and plots probabilities of a two-sided Fisher's exact test. It returns a list that consists of a vector of exact p-values and a vector of accumulative p-values over all test inputs (from 0 to n).

### <a name="buildNetwork"></a>buildNetwork.R

```R
network <- buildNetwork(df, w)

# Example usage

network <- buildNetwork((df[df$val <= max.val, ])[, c("node1", "node2", "val")])
forceNetwork(Links = network[["links"]], Nodes = network[["nodes"]], Source = "node1", Target = "node2",
                   Value = "val", NodeID = "name", Group = "group", linkWidth = 1, linkColour = "grey",
                   zoom = TRUE)  # takes as input the three ordered columns of the original data frame
```

This function prepares a list for construction a network using the function forceNetwork from the networkD3 package. The argument w is the weight multiplied with the link values for data transformation.

### <a name="labelGenerator"></a>labelGenerator.R
This function generates labels for an axis at a given interval so that ticks can be denser than labels.

```R
labelGenerator(ticks = seq(0, 10, by = 1), interval = 5, labels = c("0", "5k", "10k"))
```
