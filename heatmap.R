## DGT Manuscript - 08.01.2017 --------------------------------
## ------------------------------------------------------------
## This script conducts a correlation analysis on an soil 
## dataset and generates a heatmap 
## ------------------------------------------------------------

# import entire data set, which contains both experiments
soil.dat <- read.csv("soil_data.csv") 
str(soil.dat)

# convert ojbect from integer to factor
soil.dat$object <- as.factor(soil.dat$object) 
head(soil.dat)

# remove slag data so treatments are consistent across sites and experiments
analysis.dat <- soil.dat[!(soil.dat$source == "slag") & 
                           !(soil.dat$source == "incorp"),] 

# check to make sure subsetting worked properly
str(analysis.dat)
analysis.dat

## -------------------------------------------------------- ##
## Testing assumptions and running Pearson correlation test
## -------------------------------------------------------- ##

# first need to check whether each variable's distribution is 
# approximately normal
result <- analysis.dat[,c(8:20)]
lshap <- lapply(result, shapiro.test)
lres <- sapply(lshap, `[`, c("statistic","p.value"))
df_shap <- t(lres)
df_shap

# log transforming the whole data.frame, since
# the majority of the columns failed shapiro.wilk
result.trans <- log(result)
head(result.trans)

# Correlation matrix with p-values. 

cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}

# Use this to dump the cor.prob output to a 4 column matrix
# with row/column indices, correlation, and p-value.
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}

# correlation matrix with p-values
cor.prob(result.trans)

# "flatten" that table
flattenSquareMatrix(cor.prob(result.trans)) ##### it's not computing correlations for variables with missing data
saveRDS(corr_matrix, file = "~/Google Drive/Personal/M.Sc. thesis/Data evaluation/R stats/corr_table.txt")

# create a correlation heatmap
# install.packages("GGally")
library(GGally)
ggcorr(result.trans[, -2], 
       method = c("pairwise.complete.obs", "pearson"),
       label_alpha = TRUE, label = TRUE)

