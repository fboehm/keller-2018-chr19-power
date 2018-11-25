##First read in the arguments listed at the command line

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
print(args)
##args is now a list of character vectors
print(args$argname)
proc_num <- as.numeric(args$argname)
print(proc_num)
trait_indic <- proc_num %/% 25 + 1 # define trait_indic 
run_num <- as.numeric(args$run_num)
print(run_num)
#(nsnp <- as.numeric(args$nsnp))
#(s1 <- as.numeric(args$s1))
start_snp <- 400 # full scan is 1000 by 1000 and starts at index 400 in probs object
###############

library(dplyr)
# load expression traits
readRDS("data-to-condor/chr19_expr_10mb.rds") -> locals
trait_id <- colnames(locals)[trait_indic]
readRDS("data-to-condor/expr_asah2.rds") -> asah2
"ENSMUSG00000024887" -> asah2_id
# define phenames
phenames <- c(trait_id, asah2_id)

# load chr2 allele probabilities
readRDS("genoprobs_chr19-keller.rds") -> geno # genoprobs_chr19.rds is on SQUID

# load kinship matrix (LOCO, ie, for chromosome 2, ie, doesn't use chr2 data-to-condor)
readRDS("data-to-condor/kinship-chr19-keller.rds") -> kinship

# load covariates
readRDS("data-to-condor/addcovar-keller.rds") -> covar

# load annot
annot <- readRDS("data-to-condor/annot.rds")

# create matrix of two expression traits
pheno <- as.matrix(cbind(locals[ , trait_indic, drop = FALSE], asah2))
rownames(pheno) <- rownames(locals)
# get only shared subjects with no missing pheno or covariates

# remove subjects with missing data

id2keep <- rownames(locals)
gg <- geno
gg2 <- gg[rownames(gg) %in% id2keep, , , drop = FALSE]
kk <- kinship
kk2 <- kk[rownames(kk) %in% id2keep, colnames(kk) %in% id2keep, drop = FALSE]
cc2 <- covar[rownames(covar) %in% id2keep, , drop = FALSE]


# define s1 and nsnp
asah2_index <- annot %>%
  filter(gene_id == asah2_id) %>%
  select(marker_index) %>%
  unlist()
trait_index <- annot %>%
  filter(gene_id == trait_id) %>%
  select(marker_index) %>%
  unlist()
nsnp <- 1000
#nsnp <- 5 # for trial runs only!

# two-dimensional scan
library(qtl2pleio)
probs = gg2
kinship = kk2
addcovar = cc2[ , -5] # need to remove column 5 because we have no mice from wave 5
                  #covariates = cc2[ , -5], # need to remove column 5 because we have no mice from wave 5
#start_snp = s1
n_snp = nsnp
max_iter = 10000
# code below is from scan_pvl's code
#########
id2keep <- make_id2keep(probs = probs,
                        pheno = pheno,
                        addcovar = addcovar,
                        kinship = kinship
)
# remove - from id2keep vector - subjects with a missing phenotype or covariate
pheno <- subset_input(input = pheno, id2keep = id2keep)
subjects_phe <- check_missingness(pheno)
id2keep <- intersect(id2keep, subjects_phe)
if (!is.null(addcovar)) {
  addcovar <- subset_input(input = addcovar, id2keep = id2keep)
  subjects_cov <- check_missingness(addcovar)
  id2keep <- intersect(id2keep, subjects_cov)
}
# Send messages if there are two or fewer subjects
if (length(id2keep) == 0){stop("no individuals common to all inputs")}
if (length(id2keep) <= 2){
  stop(paste0("only ", length(id2keep),
              " common individual(s): ",
              paste(id2keep, collapse = ": ")))
}
# subset inputs to get all without missingness
probs <- subset_input(input = probs, id2keep = id2keep)
pheno <- subset_input(input = pheno, id2keep = id2keep)
if (!is.null(kinship)) {
  kinship <- subset_kinship(kinship = kinship, id2keep = id2keep)
}
if (!is.null(addcovar)) {
  addcovar <- subset_input(input = addcovar, id2keep = id2keep)
}
if (!is.null(kinship)){
  # covariance matrix estimation
  message(paste0("starting covariance matrices estimation with data from ", length(id2keep), " subjects."))
  # first, run gemma2::MphEM(), by way of calc_covs(), to get Vg and Ve
  cc_out <- calc_covs(pheno, kinship, max_iter = 10000, max_prec = 1 / 100000000, covariates = addcovar)
  Vg <- cc_out$Vg
  Ve <- cc_out$Ve
  message("covariance matrices estimation completed.")
  # define Sigma
  Sigma <- calc_Sigma(Vg, Ve, kinship)
}

# define Sigma_inv
Sigma_inv <- solve(Sigma)
# prepare table of marker indices for each call of scan_pvl
mytab <- prep_mytab(d_size = 2, n_snp = n_snp)
# subset mytab based on indices for each job, where 25 jobs constitute a single 1000 by 1000 marker scan
s_index <- (proc_num %% 25) * 40000 # each 2d scan - 1000 by 1000 markers - is broken into 25 jobs; each job is 40000 model fits, 
# ie, 200 by 200 grid - NO! this is wrong.
# each job is 40k model fits, but it is not a 200 by 200 grid, due
# to the way that we make mytab.
# First index of mytab is rep(1:1000, times = 40), ie, it runs 1 to 1000 
# and repeats 40 times
#s_index tells us where, in mytab, to start the subsetting to get mytab_sub
mytab_sub <- mytab[(s_index + 1):(s_index + 40000), ]
# set up parallel analysis
list_result <- parallel::mclapply(
                                  X = as.data.frame(t(mytab_sub)),
                                  FUN = fit1_pvl,
                                  addcovar = addcovar,
                                  probs = probs,
                                  inv_S = Sigma_inv,
                                  S = Sigma,
                                  start_snp = start_snp, # set start_snp = 1 and use mytab_sub 
                                  # to indicate which markers to put into design matrix
                                  pheno = pheno,
                                  mc.cores = 1 # use only one core per job; ie, no parallelization at this level!
)
mytab_sub$loglik <- unlist(list_result)
# trick is to get the correct marker ids here!

# first, isolate the 1000 marker ids
mm <- dimnames(probs)[[3]][start_snp:(start_snp + n_snp - 1)] # markers for the entire scan
m1 <- rep(mm, times = 40)
## m2 is the tricky part, because it depends on job number
# job 0: 1 to 40 with trait 1
# job 1: 41 to 80 (tr1)
# job 24: 961 to 1000 (with trait 1)
# job 25: 1 to 40 with trait 2
m2_ind <- 40 * (proc_num %/% 25) 
m2_pre <- mm[(m2_ind + 1):(m2_ind + 40)]
m2 <- rep(m2_pre, each = 1000)
mytab2 <- tibble::as_tibble(m1, m2)
mytab2$loglik <- mytab_sub$loglik
## save results
fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
write.table(mytab2, fn_out, quote = FALSE)
devtools::session_info()
q("no")


