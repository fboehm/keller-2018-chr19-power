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
pheno <- cbind(locals[ , trait_indic, drop = FALSE], asah2)
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
s1 <- min(asah2_index, trait_index) - 15
#nsnp <- abs(asah2_index - trait_index) + 1 + 15 + 15
nsnp <- 5 # for trial runs only!

# two-dimensional scan
library(qtl2pleio)
s_out <- scan_pvl(probs = gg2,
                  pheno = pheno,
                  kinship = kk2,
#                  addcovar = cc2[ , -5], # need to remove column 5 because we have no mice from wave 5
                  covariates = cc2[ , -5], # need to remove column 5 because we have no mice from wave 5
                  start_snp = s1,
                  n_snp = nsnp, max_iter = 10000
)


fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
write.table(s_out, fn_out, quote = FALSE)
devtools::session_info()
q("no")
