# From example in manual:
# https://cran.r-project.org/web/packages/SmartSVA/SmartSVA.pdf
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3808-1

# Also see some tutorials and Andy's paper:
# http://genomicsclass.github.io/book/pages/pca_svd.html
# https://www.biorxiv.org/content/early/2017/03/26/120899

# Also see CMS paper and tool:
# https://github.com/haschard/CMS/blob/master/CMS_v1.0.R
# https://www.nature.com/articles/ng.3975.pdf
library(SmartSVA)
## Methylation M values (CpG by Sample)
Y <- matrix(rnorm(20*1000), 1000, 20)
df <- data.frame(pred=gl(2, 10))
## Determine the number of SVs
Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
## Add one extra dimension to compensate potential loss of 1 degree of freedom
##  in confounded scenarios (very important)
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
mod <- model.matrix( ~ pred, df)
sv.obj <- smartsva.cpp(Y, mod, mod0=NULL, n.sv=n.sv)
## Speed comparison to traditional SVA
## Not run:
## Methylation M values (CpG by Sample, 27K by 1,000)
Y <- matrix(rnorm(1000*27000), 27000, 1000)
df <- data.frame(pred=gl(2, 500))
## Determine the number of SVs
Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
n.sv <- 50
mod <- model.matrix( ~ pred, df)
time_smart <- system.time(sv.obj1 <- smartsva.cpp(Y, mod, mod0=NULL, B=5, alpha = 1, VERBOSE=TRUE, n.sv=n.sv))
time_sva <- system.time(sv.obj2 <- sva(Y, mod, mod0=NULL, B=5,  n.sv=n.sv))
## Check if the solutions are the same
head(sv.obj1$sv)
head(sv.obj2$sv)
# Not the same
# Opposite sign, needs rounding, otherwise looks like the same.
identical(abs(round(sv.obj1$sv[50, 50], 4)),
          abs(round(sv.obj2$sv[50, 50], 4)))
identical(abs(round(sv.obj1$sv[1:50, 1:50], 4)),
          abs(round(sv.obj2$sv[1:50, 1:50], 4)))
identical(abs(round(sv.obj1$sv[, 1], 4)),
          abs(round(sv.obj2$sv[, 1], 4)))
identical(abs(round(sv.obj1$sv[1, ], 4)),
          abs(round(sv.obj2$sv[1, ], 4)))

# Explore
dim(sv.obj1$sv)
dim(sv.obj2$sv)

str(sv.obj1) # smart
str(sv.obj2) # normal

