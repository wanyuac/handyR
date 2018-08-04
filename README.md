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
