### HARMAN ### 

# =================================================================================== #
#
# You need to provide two csv files. 1) data 2) sample information
# Follow the data format of the examples closely. Sample names in data.csv must be first 
# row and EXACTLY match the first column in sample_info.csv. Also ensure that sample order 
# is the same order for both .CSV to fix issues with sample names later.
# 
# !!! I would return all the commands one by one. Not to be run without supervision !!!
#
# =================================================================================== #


# =================================================================================== #
# Set up packages and directory
# =================================================================================== #

# if you haven't installed these packages...
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("Biobase","Harman","msmsEDA","knitr"))

library(Biobase)
library(Harman)
library(msmsEDA)
library(knitr)

# sets working directory to the same location as this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# or hash the above command and just set it to whatever folder 
setwd('E:\\R\\XCMS\\170522_leggat\\Harman') # don't forget the apostrophes around the path


# =================================================================================== #
# Harman
# =================================================================================== #

# Import data and make ExpressionSet ================================================ #
# replace "data.csv" and/or "sample_info.csv" with "yourfilename.csv" if you need
data_mat <- as.matrix(read.table("data.csv", header = TRUE,sep = ","))
sample_info <- read.table("sample_info.csv", sep = ",", header = TRUE, row.names = 1)

# grouping info
pd_t <-data.frame(sample_info_timepts)
# metabolite info -> not necessary to name them, just generates list of numbers from
# 1 to end of metabolite list
fd <-data.frame(featurenames = 1:dim(data_mat)[1])

# create ExpressionSet
ALLset <-MSnSet(exprs = data_mat,
                pData = pd,
                fData = fd)

## if error message shows: "Error in validObject(.Object):invalid class "MSnSet" object:sampleNames 
## differ between assayData and phenoData"
## -> mismatch between colnames(data_mat) and rownames(pd): can be fixed by unhash the below, 
## ensure that BOTH both files have the same sample order!! 
# rownames(pd) <- colnames(data_mat) 


# DATA MATRIX CHECK ================================================================= #

# Any NA's in the set? TRUE = yes, FALSE = no
table(is.na(data_mat))
# Are any rows full of zeroes? TRUE = yes, FALSE = no
table(apply(data_mat, 1, function(x) sum(x == 0) == nrow(data_mat)))

## If you answered YES to any of the above questions, then unhash the next command AND
## ALLset.pp <- pp.msms.data(ALLset)  # ANDDDDD
## !!!!!!!!!!!! change all instances below of 'ALLset' to 'ALLset.pp' !!!!!!!!!!!!


## Pre-Harman Prep =================================================================== #

## combine variables (not counting batch number) if no. variables > 1
## Harman only takes one experimental variable

## EXAMPLE: experiment with 3 variables - combined to make 18 groups 
## then unhash next line (after adjusting for number variables and column names: best to stick with var1, var2, ..
# expt <-paste(pData(ALLset)$var1, pData(ALLset)$var2, pData(ALLset)$var3, sep = "")
## var1, var2, var3 columns are essentailly concatenated to form one column. 
## (in place of var1, var2, var3, these will be your column names from sample_info.csv)
## '$' allows you to choose a certain column
## I have also specified that the variables be separated by no space (sep = "")
## expt is all variables EXCEPT batch number

## IF you only have ONE variable (e.g. var1) AND batch number, then run from here:
## replace var1 with variable name if different
expt <- pData(ALLset)$var1

# batch information
batch <-pData(ALLset)$batch

# visualise groupings + batch
table(expt, batch)

# data is logged to ^2
exprs <- exprs(ALLset)
log_exprs <- log(exprs(ALLset) + 1, 2)

# use next line if you would like to plot 2 plots in a row. 
# adjust plot layout using by chaning: c( number of rows, plots per row) if requiring to adjust plot layout
par(mfrow=c(1,2))
# Have a look at the distribution after logging
plot(density(exprs))
plot(density(log_exprs))

 
# Harman ============================================================================ # 

# saves Harman result as an object 'hm' 
# confidence limit can be changed. below it is set at 0.95. 
# increase from 0.95 and it will be more conservative
hm <-harman(log_exprs, expt = expt, batch = batch, 0.95)

# Plotting ========================================================================== #

# summary in text of groupings and PC scores
summary(hm)

# print the summary into separate csv file. possibly other better ways, but difficulties with using summary and write.csv circumvented
sink("hm_timepts.csv")
print(summary(hm))
sink()


# Before and after PC SCORES (not of data)
plot(hm) # both original and corrected
pcaPlot(hm, this = "original") # separate
pcaPlot(hm, this = "corrected") # separate
plot(hm, 1, 7) # can plot selected PCs, this plots PC1 v PC4

# arrowPlot
arrowPlot(hm)

# Post-Harman ======================================================================= #

# rebuild data, de-log
log_exprs_hm <- reconstructData(hm)
exprs_hm <- 2^log_exprs_hm - 1

# save corrected data as csv file
write.csv(exprs_hm, 'corrected_data_timepts.csv')
# extract corrections and confidence from Harman and save as .csv
write.csv(hm$stats$confidence, 'confidence_timepts.csv')
write.csv(hm$stats$correction, 'correction_timepts.csv')


# More Plotting ===================================================================== #
# PC plots of uncorrected logged DATA
prcompPlot(log_exprs, colFactor=batch, pchFactor=expt, 
	main='Original (colour by batch), PC1 v PC2')
prcompPlot(log_exprs, colFactor=expt, pchFactor=batch, 
	main='Original (colour by expt), PC1 v PC2')

# PC plots of hm corrected logged DATA
## if adrian wants to understand this better, read the fucking email trail
prcompPlot(log_exprs_hm, colFactor=batch, pchFactor=expt, 
	main='Corrected (colour by batch), PC1 v PC2')
prcompPlot(log_exprs_hm, colFactor=expt, pchFactor=batch, 
	main='Corrected (colour by expt), PC1 v PC2')

# PC plots of specific axes by batch: just change the pc_x and pc_y numbers
prcompPlot(log_exprs, pc_x = 15, pc_y = 16, colFactor=batch, pchFactor=expt,
           main='Original (colour by batch), PC15 v PC16', palette = "rainbow")
prcompPlot(log_exprs_hm, pc_x = 15, pc_y = 16, colFactor=batch, pchFactor=expt,
           main='Corrected (colour by batch), PC15 v PC16')

# PC plots of specific axes by experiment: just change the pc_x and pc_y numbers
prcompPlot(log_exprs, pc_x = 15, pc_y = 16, colFactor=expt, pchFactor=batch,
           main='Original (colour by expt), PC15 v PC16', legend = FALSE)
prcompPlot(log_exprs_hm, pc_x = 15, pc_y = 16, colFactor=expt, pchFactor=batch,
           main='Corrected (colour by expt), PC15 v PC16', legend = FALSE)


# Custom changes to legend: paste following to define new plot function. use adrian_prcompPlot like prcompPlot.
# use 'legend_position' and 'elsize' to alter position and size of legend. -> see 'legend' function for more
# info/options to alter legend. 


adrian_prcompPlot <- function (object, pc_x = 1, pc_y = 2, scale = FALSE, colFactor = NULL, 
          pchFactor = NULL, palette = "rainbow", legend = TRUE, legend_position = "topright", elsize = 0.7, offset = 0, iffset = 0, ...) 
{
  if (class(object) == "data.frame") {
    object <- as.matrix(object)
  }
  if (class(object) == "matrix" & typeof(object) %in% c("double", 
                                                        "integer")) {
    object <- stats::prcomp(t(object), retx = TRUE, center = TRUE, 
                            scale. = scale)
  }
  if (class(object) != "prcomp") {
    stop("Require an instance of 'prcomp', a matrix of type 'double' or\n 'integer', or a data.frame coercible to such a matrix.")
  }
  if (is.null(colFactor)) {
    legend <- FALSE
    colFactor <- rep("", ncol(object$x))
  }
  if (is.null(pchFactor)) {
    pchFactor <- rep(1, ncol(object$x))
  }
  if (!is.factor(colFactor)) {
    colFactor <- factor(colFactor)
  }
  if (!is.factor(pchFactor)) {
    pchFactor <- factor(pchFactor)
  }
  if (length(colFactor) != ncol(object$x)) {
    stop("The length of colFactor and object do not match.")
  }
  if (length(pchFactor) != ncol(object$x)) {
    stop("The length of pchFactor and object do not match.")
  }
  mypchs <- (1:length(levels(pchFactor)))[pchFactor]
  factor_names <- levels(colFactor)
  num_levels <- length(factor_names)
  mypalette <- match.fun(palette)(num_levels)
  mycols <- mypalette[colFactor]
  lposition <- legend_position
  ofs <- offset
  ifs <- iffset
  lsize <- elsize
  graphics::plot(object$x[, pc_x], object$x[, pc_y], xlab = paste("PC", 
                                                                  pc_x, sep = ""), ylab = paste("PC", pc_y, sep = ""), 
                 col = mycols, pch = mypchs, ...)
  if (legend == TRUE) {
    graphics::legend(x = lposition, inset = c(ofs,ifs), legend = factor_names, 
                     fill = mypalette, cex = lsize, bg = "transparent", xpd = NA)
  }
}


adrian_prcompPlot(log_exprs, pc_x = 1, pc_y = 2, colFactor=expt, pchFactor=batch,
           main='Original (colour by expt), PC1 v PC2', legend = TRUE, legend_position = "bottomright", elsize = 0.4, offset = -0.2, iffset = 0)
adrian_prcompPlot(log_exprs_hm, pc_x = 1, pc_y = 2, colFactor=expt, pchFactor=batch,
           main='Corrected (colour by expt), PC1 v PC2', legend = TRUE, legend_position = "bottomright", elsize = 0.4, offset = -0.4, iffset = 0)
 





 

# =================================================================================== #
# gPCA
# =================================================================================== #

# not completely necessary to do this bit. but if you want go ahead

# if you haven't installed this package...
# install.packages('gPCA')

library(gPCA)

# gPCA on raw data before hm (no log)
pre_hm_gpca <- gPCA.batchdetect(t(exprs), batch = batch, nperm=1e3, center = F)
pre_log_hm_gpca <- gPCA.batchdetect(t(log_exprs), batch = batch, nperm=1e3, center = F)

gDist(pre_hm_gpca)
gDist(pre_log_hm_gpca)


# gPCA delta ; gPCA p value
print(c('delta:',pre_hm_gpca$delta)); print(c('pval: ',pre_hm_gpca$p.val))
# gPCA % variation explained by batch
((pre_hm_gpca$varPCg1-pre_hm_gpca$varPCu1)/pre_hm_gpca$varPCg1)*100


# gPCA on data after hm (no log)
hm_gpca <-gPCA.batchdetect(t(exprs_hm), batch = batch, center = F)
gDist(hm_gpca)
# gPCA delta ; gPCA p value
print(c('delta:',hm_gpca$delta)); print(c('pval: ',hm_gpca$p.val))
# gPCA % variation explained by batch
((hm_gpca$varPCg1-hm_gpca$varPCu1)/hm_gpca$varPCg1)*100

