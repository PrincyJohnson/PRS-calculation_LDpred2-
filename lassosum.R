install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
library(devtools)
install_github("tshmak/lassosum")
library(lassosum)
library(methods)
library(magrittr)
library(parallel)

##################33
cl <- makeCluster(2)
cl
sum_stat <- "Height.QC.gz"
sum_stat
bfile <- "EUR.QC"
bfile

# Read in and process the covariates
covariate <- fread("EUR.cov")
pcs <- fread("EUR.eigenvec") %>%
  setnames(., colnames(.), c("FID","IID", paste0("PC",1:6)))
# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- merge(covariate, pcs)
cov
# We will need the EUR.hg19 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
ld.file <- "EUR.hg19"
ld.file
# output prefix
prefix <- "EUR"
# Read in the target phenotype file
target.pheno <- fread("EUR.height")[,c("FID", "IID", "Height")]
# Read in the summary statistics
ss <- fread(sum_stat)
ss
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]
ss
# Transform the P-values into correlation
cor <- p2cor(p = ss$P,
             n = ss$N,
             sign = log(ss$OR)
)
cor
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]

# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline(
  cor = cor,
  chr = ss$CHR,
  pos = ss$BP,
  A1 = ss$A1,
  A2 = ss$A2,
  ref.bfile = bfile,
  test.bfile = bfile,
  LDblocks = ld.file, 
  cluster=cl
)
out
# Store the R2 results
target.res <- validate(out, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov))
target.res
# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
r2

