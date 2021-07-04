getwd();
install.packages("remotes")
library(remotes)
remotes::install_github("https://github.com/privefl/bigsnpr.git")
library(bigsnpr)
setwd("C:/Users/rraag/Desktop/McGill PhD/PING analysis/Data/GWAS")
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

####################Reading phenoypes and covariate files######################

library(data.table)
library(magrittr)
pheno_file <- fread("EUR.height")
pheno_file
cov_file <- fread("EUR.cov")
cov_file
pc_file <- fread("EUR.eigenvec")
pc_file
##Rename colums for Pc file##############3
colnames(pc_file) <- c("FID","IID",paste0("PC",1:6))
pc_file
###########Merging the files to form a required file ################
pheno <- merge(pheno_file,cov_file) %>% merge(.,pc_file)
pheno



##########################Obtaining SNPs#########################


info <- readRDS(url("https://ndownloader.figshare.com/files/25503788"))

##################3
# Read in the summary statistic file
sumstats <- bigreadr::fread2("Height.QC.gz") 
sumstats
# LDpred 2 require the header to follow the exact naming
names(sumstats) <-
  c("chr",
    "pos",
    "rsid",
    "a1",
    "a0",
    "n_eff",
    "beta_se",
    "p",
    "OR",
    "INFO",
    "MAF")
# Transform the OR into log(OR)
sumstats$beta <- log(sumstats$OR)
# Filter out hapmap SNPs
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
sumstats

######
# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
# preprocess the bed file (only need to do once for each data set)
snp_readBed("EUR.QC.bed")
# now attach the genotype object
obj.bigSNP <- snp_attach("EUR.QC.rds")
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats, map)
info_snp
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
View(POS2)
# calculate LD
for (chr in 1:22) {
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = NCORES,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))
#########################

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
y <- pheno[fam.order, on = c("FID", "IID")]
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
  paste0("Height~Sex+", .) %>%
  as.formula %>%
  lm(., data = y) %>%
  summary
null.r2 <- null.model$r.squared
null.r2
# Prepare data for grid model
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid.param <-
  expand.grid(p = p_seq,
              h2 = h2_seq,
              sparse = c(FALSE, TRUE))
# Get adjusted beta from grid model
beta_grid <-
  snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
beta_grid

###############333
if(is.null(obj.bigSNP)){
  obj.bigSNP <- snp_attach("EUR.QC.rds")
}
genotype <- obj.bigSNP$genotypes
# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
pred_grid <- big_prodMat(   genotype, 
                            beta_grid, 
                            ind.col = info_snp$`_NUM_ID_`)
View(pred_grid)

###############################3
reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
  paste0("Height~PRS+Sex+", .) %>%
  as.formula
reg.dat <- y
max.r2 <- 0
for(i in 1:ncol(pred_grid)){
  reg.dat$PRS <- pred_grid[,i]
  grid.model <- lm(reg.formula, dat=reg.dat) %>%
    summary  
  if(max.r2 < grid.model$r.squared){
    max.r2 <- grid.model$r.squared
  }
}
(result <- data.table(
  grid = max.r2 - null.r2,
  null = null.r2
))
reg.dat
reg.formula
