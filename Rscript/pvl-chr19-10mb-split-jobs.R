##First read in the arguments listed at the command line

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
print(args)
##args is now a list of character vectors
print(args$argname)
proc_num <- as.numeric(args$argname)
print(proc_num)
trait_indic <- proc_num + 1 # define hot_indic 
run_num <- as.numeric(args$run_num)
print(run_num)
#(nsnp <- as.numeric(args$nsnp))
#(s1 <- as.numeric(args$s1))

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
<<<<<<< HEAD
start_snp = s1
=======
#start_snp = s1
>>>>>>> cfee50398e6545f8e9d358ae2ff0ea232f62308f
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
if (!is.null(addcovar)) {
  addcovar <- drop_depcols(addcovar)
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
  cc_out <- calc_covs(pheno, kinship, max_iter = max_iter, max_prec = max_prec, covariates = addcovar)
  Vg <- cc_out$Vg
  Ve <- cc_out$Ve
  message("covariance matrices estimation completed.")
  # define Sigma
  Sigma <- calc_Sigma(Vg, Ve, kinship)
}
if (is.null(kinship)){
  # get Sigma for Haley Knott regression without random effect
  Ve <- var(pheno) # get d by d covar matrix
  Sigma <- calc_Sigma(Vg = NULL, Ve = Ve, n_mouse = nrow(pheno))
}

# define Sigma_inv
Sigma_inv <- solve(Sigma)
# prepare table of marker indices for each call of scan_pvl
mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)
# subset mytab based on indices for each job, where 25 jobs constitute a single 1000 by 1000 marker scan
<<<<<<< HEAD
s_index <- ((proc_num + 1) %% 25) * 40000 # each 2d scan is broken into 25 jobs; each job is 40000 model fits, ie, 200 by 200 grid
=======
s_index <- (proc_num %% 25) * 40000 # each 2d scan - 1000 by 1000 markers - is broken into 25 jobs; each job is 40000 model fits, ie, 200 by 200 grid
>>>>>>> cfee50398e6545f8e9d358ae2ff0ea232f62308f
mytab_sub <- mytab[(s_index + 1):(s_index + 40000), ]
# set up parallel analysis
list_result <- parallel::mclapply(
                                  X = as.data.frame(t(mytab_sub)),
                                  FUN = fit1_pvl,
                                  addcovar = addcovar,
                                  probs = probs,
                                  inv_S = Sigma_inv,
                                  S = Sigma,
                                  start_snp = 1, # set start_snp = 1 and use mytab_sub 
                                  # to indicate which markers to put into design matrix
                                  pheno = pheno,
                                  mc.cores = 1 # use only one core per job; ie, no parallelization at this level!
)
mytab$loglik <- unlist(list_result)
marker_id <- dimnames(probs)[[3]][start_snp:(start_snp + n_snp - 1)]
mytab2 <- tibble::as_tibble(apply(FUN = function(x) marker_id[x], X = mytab[, -ncol(mytab)], MARGIN = 2))
mytab2$loglik <- mytab$loglik
#########


fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
write.table(s_out, fn_out, quote = FALSE)
devtools::session_info()
q("no")
